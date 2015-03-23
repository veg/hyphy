/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon@cfenet.ubc.ca)
  Steven Weaver (sweaver@ucsd.edu)
  
Module Developers:
	Lance Hepler (nlhepler@gmail.com)
	Martin Smith (martin.audacis@gmail.com)

Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdf22@cam.ac.uk)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "calcnode.h"
#include "likefunc.h"
#include "scfg.h"
#include <math.h>
#include <float.h>

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


//#define _UBER_VERBOSE_DUMP_MATRICES
//#define _UBER_VERBOSE_DUMP 27

extern  long likeFuncEvalCallCount,
        matrixExpCount;

#ifdef MDSOCL
int launchmdsocl(long siteCount,
                 long nodeCount,
                 long alphabetDimension,
                 _SimpleList& updateNodes,
                 _SimpleList& flatParents,
                 _SimpleList& flatNodes,
                 _SimpleList& flatCLeaves,
                 _SimpleList& flatLeaves,
                 _SimpleList& flatTree,
                 _Parameter* iNodeCache,
                 long* lNodeFlags,
                 _SimpleList taggedInternals,
                 _GrowingVector* lNodeResolutions);
#endif

#ifdef  _SLKP_LFENGINE_REWRITE_


_Parameter          _lfScalerUpwards          = pow(2.,200.),
                    _lfScalingFactorThreshold = 1./_lfScalerUpwards,
                    _logLFScaler            = 200. *log(2.);

_GrowingVector      _scalerMultipliers,
                    _scalerDividers;

//_Parameter          eval_buffer [600];


/*----------------------------------------------------------------------------------------------------------*/

inline void _handle4x4_pruning_case (double* childVector, double* tMatrix, double* parentConditionals) {
#ifdef _SLKP_USE_SSE_INTRINSICS
        double tv [4]     __attribute__ ((aligned (16))) = {childVector[0],
                                                            childVector[1],
                                                            childVector[2],
                                                            childVector[3]};
                                                    
        __m128d buffer0 = _mm_loadu_pd (tv),
                buffer1 = _mm_loadu_pd (tv+2),
                matrix01 = _mm_loadu_pd (tMatrix),
                matrix12 = _mm_loadu_pd (tMatrix+2),
                matrix34 = _mm_loadu_pd (tMatrix+4),
                matrix56 = _mm_loadu_pd (tMatrix+6),
                reg_storage  = _mm_mul_pd (buffer0, matrix01),
                reg_storage2 = _mm_mul_pd (buffer0, matrix34);
               
        matrix34     = _mm_mul_pd(buffer1, matrix12);
        matrix56     = _mm_mul_pd(buffer1, matrix56);
        reg_storage  = _mm_add_pd (reg_storage, matrix34);
        reg_storage2 = _mm_add_pd (reg_storage2, matrix56);
        reg_storage  = _mm_hadd_pd (reg_storage,reg_storage2);
        matrix01 = _mm_loadu_pd (parentConditionals);
        matrix01 = _mm_mul_pd (reg_storage, matrix01);
        _mm_storeu_pd (parentConditionals, matrix01);
        
                
        
        matrix01 = _mm_loadu_pd (tMatrix+8);
        matrix12 = _mm_loadu_pd (tMatrix+10);
        matrix34 = _mm_loadu_pd (tMatrix+12);
        matrix56 = _mm_loadu_pd (tMatrix+14);
        reg_storage  = _mm_mul_pd (buffer0, matrix01);
        reg_storage2 = _mm_mul_pd (buffer0, matrix34);                
              
        matrix34     = _mm_mul_pd(buffer1, matrix12);
        matrix56     = _mm_mul_pd(buffer1, matrix56);
        reg_storage  = _mm_add_pd (reg_storage, matrix34);
        reg_storage2 = _mm_add_pd (reg_storage2, matrix56);
        reg_storage  = _mm_hadd_pd (reg_storage,reg_storage2);
        
        matrix01 = _mm_loadu_pd (parentConditionals+2);
        matrix01 = _mm_mul_pd (reg_storage, matrix01);
        _mm_storeu_pd (parentConditionals+2, matrix01);

               
#else  
        _Parameter t1 = childVector[0] - childVector[3],
                   t2 = childVector[1] - childVector[3],
                   t3 = childVector[2] - childVector[3],
                   t4 = childVector[3];
                   
        parentConditionals [0] *= tMatrix[0]  * t1 + tMatrix[1] * t2 + tMatrix[2] * t3 + t4;
        parentConditionals [1] *= tMatrix[4]  * t1 + tMatrix[5] * t2 + tMatrix[6] * t3 + t4;
        parentConditionals [2] *= tMatrix[8]  * t1 + tMatrix[9] * t2 + tMatrix[10] * t3 + t4;
        parentConditionals [3] *= tMatrix[12] * t1 + tMatrix[13] * t2 + tMatrix[14] * t3 + t4;
#endif

}

/*----------------------------------------------------------------------------------------------------------*/

_Parameter  acquireScalerMultiplier (long s)
{
    if (s>0) {
        if (s >= _scalerMultipliers.used)
            for (long k = _scalerMultipliers.used; k <= s; k++) {
                _scalerMultipliers.Store (exp (-_logLFScaler * k));
            }
        return _scalerMultipliers.theData[s];
    }
    s = -s;
    if (s >= _scalerDividers.used)
        for (long k = _scalerDividers.used; k <= s; k++) {
            _scalerDividers.Store (exp (_logLFScaler * k));
        }
    return _scalerDividers.theData[s];
}

/*----------------------------------------------------------------------------------------------------------*/
void        _TheTree::ExponentiateMatrices  (_List& expNodes, long tc, long catID)
{
    _List           matrixQueue,
                    nodesToDo;
                    
    _SimpleList     isExplicitForm;
    bool            hasExpForm = false;

    for (unsigned long nodeID = 0; nodeID < expNodes.lLength; nodeID++) {
        long didIncrease = matrixQueue.lLength;
        _CalcNode* thisNode = (_CalcNode*) expNodes(nodeID);
        if (thisNode->RecomputeMatrix (catID, categoryCount, nil, &matrixQueue,&isExplicitForm)) {
             hasExpForm = true;
        }
        #ifdef _UBER_VERBOSE_DUMP
          if (likeFuncEvalCallCount == _UBER_VERBOSE_DUMP)
            printf ("NodeID %d (%s). Old length %ld, new length %ld\n", nodeID, thisNode->GetName()->sData, didIncrease,matrixQueue.lLength); 
        #endif
        if ((didIncrease = (matrixQueue.lLength - didIncrease))) {
            for (long copies = 0; copies < didIncrease; copies++) {
                nodesToDo << thisNode;
            }
        }
    }
    
    //printf ("%ld %d\n", nodesToDo.lLength, hasExpForm);

    unsigned long matrixID;
    
    _List * computedExponentials = hasExpForm? new _List (matrixQueue.lLength) : nil;
    
#ifdef _OPENMP
    unsigned long nt = cBase<20?1:(MIN(tc, matrixQueue.lLength / 3 + 1));
    matrixExpCount += matrixQueue.lLength;
#endif

    #pragma omp parallel for default(shared) private (matrixID) schedule(static) if (nt>1)  num_threads (nt)
    for  (matrixID = 0; matrixID < matrixQueue.lLength; matrixID++) {
        if (isExplicitForm.lData[matrixID] == 0) { // normal matrix to exponentiate
            ((_CalcNode*) nodesToDo(matrixID))->SetCompExp (((_Matrix*)matrixQueue(matrixID))->Exponentiate(), catID);
        } else {
             (*computedExponentials) [matrixID] = ((_Matrix*)matrixQueue(matrixID))->Exponentiate();
        }
    }
    
    
  


    if (computedExponentials) {
        _CalcNode * current_node         = nil;
        _List       buffered_exponentials;
        
        for (unsigned long mx_index = 0; mx_index < nodesToDo.lLength; mx_index++) {
            if (isExplicitForm.lData[mx_index]) {
                _CalcNode *next_node = (_CalcNode*) nodesToDo (mx_index);
                //printf ("%x %x\n", current_node, next_node);
                if (next_node != current_node) {
                    if (current_node) {
                        current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
                    }
                    current_node = next_node;
                    buffered_exponentials.Clear(true);
                    buffered_exponentials.AppendNewInstance((*computedExponentials)(mx_index));
                 }
                else {
                    buffered_exponentials.AppendNewInstance((*computedExponentials)(mx_index));
                }
            } else {
                if (current_node) {
                    current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
                }
                current_node = nil;
            }
        }
        if (current_node) {
            current_node->RecomputeMatrix (catID, categoryCount, nil, nil, nil, &buffered_exponentials);
        }
        DeleteObject(computedExponentials);
    #ifdef _UBER_VERBOSE_DUMP_MATRICES
      if (likeFuncEvalCallCount == _UBER_VERBOSE_DUMP) {
        fprintf (stderr, "\n T_MATRIX = {"); 
        for (unsigned long nodeID = 0; nodeID < flatLeaves.lLength + flatTree.lLength - 1; nodeID++) {
            bool    isLeaf     = nodeID < flatLeaves.lLength;

            _CalcNode * current_node = isLeaf? (((_CalcNode**) flatCLeaves.lData)[nodeID]):
                              (((_CalcNode**) flatTree.lData)  [nodeID - flatLeaves.lLength]);
            if (nodeID) {
              fprintf (stderr, ",");
            }
            fprintf (stderr, "\n\"%s\":%s", current_node->GetName()->sData, _String((_String*)current_node->GetCompExp()->toStr()).sData); 
            
          }
        fprintf (stderr, "\n};\n"); 
      }
    #endif


    }
}

/*----------------------------------------------------------------------------------------------------------*/

long        _TheTree::DetermineNodesForUpdate   (_SimpleList& updateNodes, _List* expNodes, long catID, long addOne, bool canClear)
{
    nodesToUpdate.Populate (flatLeaves.lLength + flatTree.lLength - 1, 0, 0);
    _CalcNode       *currentTreeNode;
    long            lastNodeID = -1;

    // look for nodes with model changes and mark the path up to the root as needing an update


    if (addOne >= 0) {
        nodesToUpdate.lData[addOne] = 1;
    }

    if (forceRecalculationOnTheseBranches.lLength) {
        for (unsigned long markedNode = 0; markedNode < forceRecalculationOnTheseBranches.lLength; markedNode++) {
            nodesToUpdate.lData[forceRecalculationOnTheseBranches.lData[markedNode]] = 1;
        }

        if (canClear) {
            forceRecalculationOnTheseBranches.Clear();
        }
    }

    for (unsigned long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++) {
        bool    isLeaf     = nodeID < flatLeaves.lLength;

        currentTreeNode = isLeaf? (((_CalcNode**) flatCLeaves.lData)[nodeID]):
                          (((_CalcNode**) flatTree.lData)  [nodeID - flatLeaves.lLength]);

        if (currentTreeNode->NeedToExponentiate (catID)) {
            if (expNodes) {
                (*expNodes) << currentTreeNode;
                //printf ("EXP>%s\n", currentTreeNode->GetName()->sData);
                lastNodeID = nodeID;
            } else {
                currentTreeNode->RecomputeMatrix (catID, categoryCount, nil);
            }

            nodesToUpdate.lData[nodeID] = 1;
        }

        if (nodesToUpdate.lData[nodeID]) {
            nodesToUpdate.lData[flatParents.lData[nodeID]+flatLeaves.lLength] = 1;
        }
    }


    // one more pass to pick up all descendants of changed internal nodes

    for (unsigned long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
        if (nodesToUpdate.lData[flatLeaves.lLength+flatParents.lData[nodeID]] && nodesToUpdate.lData[nodeID] == 0) {
            nodesToUpdate.lData[nodeID] = 1;
        }

    // write out all changed nodes

    for (unsigned long nodeID = 0; nodeID < nodesToUpdate.lLength; nodeID++)
        if (nodesToUpdate.lData[nodeID]) {
            updateNodes << nodeID;
        }

    if (expNodes && expNodes->lLength == 1) {
        return lastNodeID;
    }

    return -1;
}

/*----------------------------------------------------------------------------------------------------------*/

void        _TheTree::FillInConditionals        (_DataSetFilter*        theFilter, _Parameter*  iNodeCache,  _SimpleList*   tcc)
// this utility function will simply fill in all the conditional probability vectors for internal nodes,
// including those that were skipped due to column sorting optimization
// this is useful to avoid code duplication for other functions (e.g. ancestral sampling) that
// make use of conditional probability vectors, but would not benefit from subtree caching
{
    if (!tcc) {
        return;
    }

    long            alphabetDimension     =         theFilter->GetDimension(),
                    siteCount           =         theFilter->NumberDistinctSites();

    for  (long nodeID = 0; nodeID < flatTree.lLength; nodeID++) {
        _Parameter * conditionals       = iNodeCache +(nodeID  * siteCount) * alphabetDimension;
        long        currentTCCIndex     = siteCount * nodeID,
                    currentTCCBit        = currentTCCIndex % _HY_BITMASK_WIDTH_;

        currentTCCIndex /= _HY_BITMASK_WIDTH_;
        for (long siteID = 0; siteID < siteCount; siteID++, conditionals += alphabetDimension) {
            if (siteID  && (tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0) {
                for (long k = 0; k < alphabetDimension; k++) {
                    conditionals[k] = conditionals[k-alphabetDimension];
                }
            }
            if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                currentTCCBit   = 0;
                currentTCCIndex ++;
            }

        }
    }
}

/*----------------------------------------------------------------------------------------------------------*/

#ifdef MDSOCL
_Parameter _TheTree::OCLLikelihoodEvaluator (			_SimpleList&		     updateNodes, 
														_DataSetFilter*		 theFilter,
                                                        _Parameter*			 iNodeCache,
                                                        long	   *			 lNodeFlags,
                                                        _GrowingVector*		 lNodeResolutions,
														_OCLEvaluator& OCLEval)
{

	_SimpleList		taggedInternals					(flatNodes.lLength, 0, 0);
	//printf("Launching a tree in OpenCL...\n");
	return ((_Parameter) OCLEval.launchmdsocl(			updateNodes,
														flatParents,
														flatNodes,
														flatCLeaves,
														flatLeaves,
														flatTree,
														theProbs,
														theFilter->theFrequencies,
														lNodeFlags,
														taggedInternals,
														lNodeResolutions));
}

#endif

/*----------------------------------------------------------------------------------------------------------*/

_Parameter  _TheTree::VerySimpleLikelihoodEvaluator   (_SimpleList&          updateNodes,
        _DataSetFilter*      theFilter,
        _Parameter*          iNodeCache,
        long       *             lNodeFlags,
        _GrowingVector*      lNodeResolutions)
{
// set some useful global variables
    _SimpleList     taggedInternals                 (flatNodes.lLength, 0, 0);
    // as we are traversing the tree, we are ensuring that each node which changes
    // forces the parent to change and preps the parent cache for computations
    // by zeroing everything out

    // flatNodes is the array of indices (post-order traveral) of _internal nodes_
    // for example ((A,C)N1,D,(B,E)N2)Root will have
    // N1 (index 0), N2 (index 1), Root (index 2) in this array

    long            alphabetDimension     =         theFilter->GetDimension(),
                    // the number of characters (and the dimension of transition matrices)

                    siteCount           =         theFilter->NumberDistinctSites();
    // how many unique sites are there

    for  (long nodeID = 0; nodeID < updateNodes.lLength; nodeID++) {
        long    nodeCode   = updateNodes.lData [nodeID],
                parentCode = flatParents.lData [nodeCode];
        // the INDEX of the parent node for the current node;
        // this list (a member of _TheTree is computed when the tree is created
        // for ((A,C)N1,D,(B,E)N2)Root this array will store
        // 0 (A), 0(C), 2(D), 1(B), 1(E), 2(N1), 2(N3), -1 (root)

        bool    isLeaf     = nodeCode < flatLeaves.lLength;

        if (!isLeaf) {
            nodeCode -=  flatLeaves.lLength;
        }

        _Parameter * parentConditionals = iNodeCache +  (parentCode  * siteCount) * alphabetDimension;
        // this is a convenience pointer into the global cache
        // each node will have siteCount*alphabetDimension contiguous doubles



        if (taggedInternals.lData[parentCode] == 0)
            // mark the parent for update and clear its conditionals if needed
        {
            taggedInternals.lData[parentCode]     = 1; // only do this once
            for (long k = 0, k3 = 0; k < siteCount; k++)
                for (long k2 = 0; k2 < alphabetDimension; k2++) {
                    parentConditionals [k3++] = 1.0;
                }
        }



        _Parameter  *       tMatrix = (isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
                                       ((_CalcNode*) flatTree    (nodeCode)))->GetCompExp(0)->theData;
        // in the host code the transition matrix is retrieved from the _CalcNode object
        // in the devide code there will probably be another array to grab it from

        _Parameter  *       childVector; // the vector of conditional probabilities for
        // THIS node (nodeCode)

        if (!isLeaf) {
            childVector = iNodeCache + (nodeCode * siteCount) * alphabetDimension;
        }
        // if THIS node is internal, simply look up the conditional vector in the same cache
        // as the parent (just a different offset)

        for (long siteID = 0; siteID < siteCount; siteID++,
                parentConditionals += alphabetDimension) {
            _Parameter  sum      = 0.0;

            if (isLeaf)
                // the leaves do NOT have a conditional probability vector because most of them are fully resolved
                // For example the codon AAG is going to map to (AAA)0,(AAC)0,(AAG)1,.....,(TTT)0, so we can
                // just store this as index 2, which says to put a '1' in the 2-nd index of the array (and assume the rest are 0)
                // this is stored in lNodeFlags
                // For ambiguous codons, e.g. AAR = {AAA,AAG}, the host code will resolve this to 1,0,1,0,...0 at prep time (because this will never change)
                // and stores it in the lNodeResolutions (say at index K) and puts -K-1 into lNodeFlags, so that we now have to look it up here
            {
                long siteState = lNodeFlags[nodeCode*siteCount + siteID] ;

                if (siteState >= 0)
                    // a single character state; sweep down the appropriate column
                    // note that one of the loops (over child state) drops out, since there is only one state
                {
                    long matrixIndex =  siteState;
                    for (long k = 0; k < alphabetDimension; k++, matrixIndex += alphabetDimension) {
                        parentConditionals[k] *= tMatrix[matrixIndex];
                    }
                    continue; // nothing else to do, move to the next site
                } else {
                    childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                }
                // look up the resolution for the ambugious node -- this will have to be on the device as well
                // but can be in constant memory
            }

            _Parameter *matrixPointer = tMatrix;

            for (long p = 0; p < alphabetDimension; p++) {
                _Parameter      accumulator = 0.0;

                for (long c = 0; c < alphabetDimension; c++) {
                    accumulator +=  matrixPointer[c]   * childVector[c];
                }

                matrixPointer                 += alphabetDimension;

                sum += (parentConditionals[p] *= accumulator);
            }

            // if (sum < small_number) -- handle underflow

            childVector    += alphabetDimension;
            // shift the position in childvector to the next site
        }
    }

    // now just process the root and return the likelihood

    _Parameter  * rootConditionals = iNodeCache + alphabetDimension * ((flatTree.lLength-1)  * siteCount),
                  // the root is always the LAST internal node in all lists
                  result = 0.0;


    for (long siteID = 0; siteID < siteCount; siteID++) {
        _Parameter accumulator = 0.;
        for (long p = 0; p < alphabetDimension; p++,rootConditionals++) {
            accumulator += *rootConditionals * theProbs[p];
        }
        /* theProbs is a member variable of the tree, which basically determines
           what probability there is to observe a given character at the root
           in simple cases it is fixed for the duration of optimization, but for more
           complex models it may change from iteration to iteration
         */
        result += log(accumulator) * theFilter->theFrequencies [siteID];
        // correct for the fact that identical alignment columns may appear more than once
    }


    return result;
}


/*----------------------------------------------------------------------------------------------------------*/

_Parameter      _TheTree::ComputeTreeBlockByBranch  (                   _SimpleList&        siteOrdering,
        _SimpleList&        updateNodes,
        _SimpleList*        tcc,
        _DataSetFilter*     theFilter,
        _Parameter*         iNodeCache,
        long      *         lNodeFlags,
        _Parameter*         scalingAdjustments,
        _GrowingVector*     lNodeResolutions,
        long&               overallScaler,
        long                siteFrom,
        long                siteTo,
        long                catID,
        _Parameter*         storageVec,
        long*               siteCorrectionCounts,
        long                setBranch,
        long*               setBranchTo
                                                    )
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
    // process the leaves first

    _SimpleList     taggedInternals                 (flatNodes.lLength, 0, 0);
    long            alphabetDimension     =         theFilter->GetDimension(),
                    siteCount           =         theFilter->NumberDistinctSites(),
                    alphabetDimensionmod4  =      alphabetDimension-alphabetDimension%4;

    _CalcNode       *currentTreeNode;
    long            localScalerChange     =         0;

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    for  (unsigned long nodeID = 0; nodeID < updateNodes.lLength; nodeID++) {
        long    nodeCode   = updateNodes.lData [nodeID],
                parentCode = flatParents.lData [nodeCode];

        bool    isLeaf     = nodeCode < flatLeaves.lLength;

        if (!isLeaf) {
            nodeCode -=  flatLeaves.lLength;
        }

        _Parameter * parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedInternals.lData[parentCode] == 0)
            // mark the parent for update and clear its conditionals if needed
        {
            taggedInternals.lData[parentCode]     = 1;
            _Parameter      _hprestrict_ *localScalingFactor      = scalingAdjustments + parentCode*siteCount;

            bool    matchSet   = (parentCode == setBranch);

            if (alphabetDimension == 4) {
                long k3     = 0;
                if (matchSet)
                    for (long k = siteFrom; k < siteTo; k++, k3+=4) {
                        parentConditionals [k3]   = 0.;
                        parentConditionals [k3+1] = 0.;
                        parentConditionals [k3+2] = 0.;
                        parentConditionals [k3+3] = 0.;
                        parentConditionals [k3+setBranchTo[siteOrdering.lData[k]]] = localScalingFactor[k];
                    }
                else
                    for (long k = siteFrom; k < siteTo; k++, k3+=4) {
                        _Parameter scaler = localScalingFactor[k];
                        parentConditionals [k3]   = scaler;
                        parentConditionals [k3+1] = scaler;
                        parentConditionals [k3+2] = scaler;
                        parentConditionals [k3+3] = scaler;
                    }
            } else {
                long k3     = 0;
                if (matchSet) {
                    for (long k = siteFrom; k < siteTo; k++) {
                        for (long k2 = 0; k2 < alphabetDimension; k2++) {
                            parentConditionals [k3+k2] = 0.;
                        }

                        parentConditionals[k3 + setBranchTo[siteOrdering.lData[k]]] = localScalingFactor[k];
                        k3                 +=   alphabetDimension;
                    }
                } else {
                    for (long k = siteFrom; k < siteTo; k++) {
                        _Parameter scaler = localScalingFactor[k];
                        for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
                            parentConditionals [k3] = scaler;
                        }
                    }
                }
            }
        }

        currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
                          ((_CalcNode*) flatTree    (nodeCode));

        _Parameter  *       _hprestrict_ transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
        _Parameter  *       childVector,
                    *     lastUpdatedSite;

        if (!isLeaf) {
            childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
        }

        long currentTCCIndex        ,
             currentTCCBit            ,
             parentTCCIIndex      ,
             parentTCCIBit            ;

        if (tcc) {
            parentTCCIIndex = siteCount * parentCode + siteFrom;
            parentTCCIBit   = parentTCCIIndex % _HY_BITMASK_WIDTH_;
            parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
            if (! isLeaf) {
                currentTCCIndex = siteCount * nodeCode + siteFrom;
                currentTCCBit   = currentTCCIndex % _HY_BITMASK_WIDTH_;
                currentTCCIndex /= _HY_BITMASK_WIDTH_;
            }
        }

        //long successiveSkips = 0;

        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += alphabetDimension) {
            if (tcc) {
                if (parentTCCIBit == _HY_BITMASK_WIDTH_) {
                    parentTCCIBit   = 0;
                    parentTCCIIndex ++;
                }

                if (siteID > siteFrom && (tcc->lData[parentTCCIIndex] & bitMaskArray.masks[parentTCCIBit]) > 0) {
                    if (!isLeaf) {
                        childVector     += alphabetDimension;
                        if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                            currentTCCBit   = 0;
                            currentTCCIndex ++;
                        }
                    }
                    parentTCCIBit++;
                    /*
                    if (likeFuncEvalCallCount == 4 || likeFuncEvalCallCount == 5) {
                        printf ("\nSKIPPED site %ld @ eval %ld\n", siteID, likeFuncEvalCallCount);
                    }*/   
                    continue;
                }
                parentTCCIBit++;
            }

            _Parameter  *tMatrix = transitionMatrix,
                         sum         = 0.0;

            char        didScale = 0;

            if (isLeaf) {
                long siteState;

                if (setBranch == nodeCode + flatTree.lLength) {
                    siteState = setBranchTo[siteOrdering.lData[siteID]] ;
                } else {
                    siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.lData[siteID]] ;
                }
                if (siteState >= 0)
                    // a single character state; sweep down the appropriate column
                {
                    tMatrix  +=  siteState;
                    if (alphabetDimension == 4) {
                        parentConditionals[0] *= tMatrix[0];
                        parentConditionals[1] *= tMatrix[4];
                        parentConditionals[2] *= tMatrix[8];
                        parentConditionals[3] *= tMatrix[12];
                    } else {
                        long k = 0;
                        for (; k < alphabetDimensionmod4; k+=4, tMatrix += alphabetDimension+alphabetDimension+alphabetDimension+alphabetDimension) {
                            parentConditionals[k]   *= tMatrix[0];
                            parentConditionals[k+1] *= tMatrix[alphabetDimension];
                            parentConditionals[k+2] *= tMatrix[alphabetDimension+alphabetDimension];
                            parentConditionals[k+3] *= tMatrix[alphabetDimension+alphabetDimension+alphabetDimension];
                        }
                        for (; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
                            parentConditionals[k] *= *tMatrix;
                        }
                    }
                    continue;
                } else {
                    childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                }
            } else {
                if (tcc) {
                    if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom) {
                        for (long k = 0; k < alphabetDimension; k++) {
                            childVector[k] = lastUpdatedSite[k];
                        }
                    }
                    if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                        currentTCCBit   = 0;
                        currentTCCIndex ++;
                    }
                    lastUpdatedSite = childVector;
                }
            }

            if (alphabetDimension == 4) { // special case for nuc data

                _handle4x4_pruning_case (childVector, tMatrix, parentConditionals);

                // handle scaling if necessary
                // the check for sum > 0.0 is necessary for 'degenerate' log-L functions (-infinity)

                sum     = parentConditionals [0] + parentConditionals [1] + parentConditionals [2] + parentConditionals [3];
                if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                    _Parameter tryScale                                 = scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards;
                    if (tryScale < HUGE_VAL) {
                        parentConditionals [0]                             *= _lfScalerUpwards;
                        parentConditionals [1]                             *= _lfScalerUpwards;
                        parentConditionals [2]                             *= _lfScalerUpwards;
                        parentConditionals [3]                             *= _lfScalerUpwards;
                        localScalerChange                                  += theFilter->theFrequencies [siteOrdering.lData[siteID]];
                        scalingAdjustments [parentCode*siteCount + siteID]  = tryScale;
                        didScale                                            = 1;
                    }
                } else {
                    if (sum > _lfScalerUpwards) {
                        scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalingFactorThreshold;
                        parentConditionals [0]                             *= _lfScalingFactorThreshold;
                        parentConditionals [1]                             *= _lfScalingFactorThreshold;
                        parentConditionals [2]                             *= _lfScalingFactorThreshold;
                        parentConditionals [3]                             *= _lfScalingFactorThreshold;
                        localScalerChange                                  -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
                        didScale                                            = -1;
                    }
                }

                childVector += 4;
            } else {
                _Parameter sum = 0.0;

#ifndef _SLKP_SSE_VECTORIZATION_

                if (alphabetDimension > alphabetDimensionmod4){
                    
                    /*_Parameter maxParentP = _lfScalingFactorThreshold;
                    
                    for (long p = 0; p < alphabetDimension; p++)
                        if (parentConditionals[p] > maxParentP) {
                            maxParentP = parentConditionals[p];
                        } 
                        
                    maxParentP *= DBL_EPSILON;  */      
                    
#ifdef _SLKP_USE_SSE_INTRINSICS
                    double buffer[2] __attribute__ ((aligned (16)));
#endif
#ifdef _SLKP_USE_AVX_INTRINSICS
                    double buffer[4] __attribute__ ((aligned (32)));
#endif
                    for (long p = 0; p < alphabetDimension; p++) {
                        /*if (parentConditionals[p] < maxParentP) {
                            tMatrix               += alphabetDimension;
                            continue;
                        }*/
                        
                        _Parameter      accumulator = 0.0;

#ifdef _SLKP_USE_SSE_INTRINSICS

                        __m128d buffer1,
                                buffer2,
                                buffer3 = _mm_setzero_pd(),
                                buffer4 = _mm_setzero_pd(),
                                load1, 
                                load2,
                                load3,
                                load4;
                                
                        
                        if (((long int)tMatrix & 0x1111b) == 0 && ((long int)childVector & 0x1111b) == 0){ 
                           for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                                load1 = _mm_load_pd (tMatrix+c);
                                load2 = _mm_load_pd (tMatrix+c+2);
                                load3 = _mm_load_pd (childVector+c);
                                load4 = _mm_load_pd (childVector+c+2);
                                buffer1 = _mm_mul_pd (load1, load3);
                                buffer2 = _mm_mul_pd (load2, load4);
                                buffer3 = _mm_add_pd (buffer1,buffer3);
                                buffer4 = _mm_add_pd (buffer2,buffer4);
                            }      
                        } else {
                           for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                                load1 = _mm_loadu_pd (tMatrix+c);
                                load2 = _mm_loadu_pd (tMatrix+c+2);
                                load3 = _mm_loadu_pd (childVector+c);
                                load4 = _mm_loadu_pd (childVector+c+2);
                                buffer1 = _mm_mul_pd (load1, load3);
                                buffer2 = _mm_mul_pd (load2, load4);
                                buffer3 = _mm_add_pd (buffer1,buffer3);
                                buffer4 = _mm_add_pd (buffer2,buffer4);
                            }      
                        
                        }
                        
                        buffer3 = _mm_add_pd (buffer3, buffer4);    
                        _mm_store_pd (buffer, buffer3);
                        accumulator = buffer[0] + buffer[1];
#elif defined _SLKP_USE_AVX_INTRINSICS
                        
                        __m256d buffer1,
                                buffer3 = _mm256_setzero_pd(),
                                load1,
                                load3;
                        
                        
                        if (((long int)tMatrix & 0x11111b) == 0 && ((long int)childVector & 0x11111b) == 0){
                            for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                                load1 = _mm256_load_pd (tMatrix+c);
                                load3 = _mm256_load_pd (childVector+c);
                                buffer1 = _mm256_mul_pd (load1, load3);
                                buffer3 = _mm256_add_pd (buffer1,buffer3);
                             }
                        } else {
                            for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                                load1 = _mm256_loadu_pd (tMatrix+c);
                                load3 = _mm256_loadu_pd (childVector+c);
                                buffer1 = _mm256_mul_pd (load1, load3);
                                buffer3 = _mm256_add_pd (buffer1,buffer3);
                            }
                            
                        }
                        
                        _mm256_store_pd (buffer, buffer3);
                        accumulator = buffer[0] + buffer[1] + buffer[2] + buffer[3];
#else
                        for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                        // 4 - unroll the loop
                             _Parameter pr1 =    tMatrix[c]   * childVector[c],
                                        pr2 =    tMatrix[c+1] * childVector[c+1],
                                        pr3 =    tMatrix[c+2] * childVector[c+2],
                                        pr4 =    tMatrix[c+3] * childVector[c+3];
                            pr1 += pr2;
                            pr3 += pr4;
                            accumulator += pr1+pr3;
                        }
#endif

                        for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
                            accumulator +=  tMatrix[c] * childVector[c];
                        }

                        tMatrix               += alphabetDimension;
                        sum += (parentConditionals[p] *= accumulator);
                    }
                }
                else
                    for (long p = 0; p < alphabetDimension; p++) {
                        _Parameter      accumulator = 0.0;

                        for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
                            accumulator +=  tMatrix[c]   * childVector[c] +
                                            tMatrix[c+1] * childVector[c+1] +
                                            tMatrix[c+2] * childVector[c+2] +
                                            tMatrix[c+3] * childVector[c+3];

                        tMatrix               += alphabetDimension;
                        sum += (parentConditionals[p] *= accumulator);
                    }
#else
                for (long p = 0; p < alphabetDimension; p++) {
                    _Parameter      accumulator = 0.0;

                    for (long c = 0; c < alphabetDimension; c++) {
                        accumulator +=  tMatrix[c]   * childVector[c];
                    }

                    tMatrix               += alphabetDimension;
                    sum += (parentConditionals[p] *= accumulator);
                }

#endif


                if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                    _Parameter tryScale                                 = scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards;
                    if (tryScale < HUGE_VAL) {
                        scalingAdjustments [parentCode*siteCount + siteID] = tryScale;
                        for (long c = 0; c < alphabetDimension; c++) {
                            parentConditionals [c] *= _lfScalerUpwards;
                        }

                        localScalerChange                                      += theFilter->theFrequencies [siteOrdering.lData[siteID]];
                        didScale                                                = 1;
                    }
                } else {
                    if (sum > _lfScalerUpwards) {
                        scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalingFactorThreshold;
                        for (long c = 0; c < alphabetDimension; c++) {
                            parentConditionals [c] *= _lfScalingFactorThreshold;
                        }
                        localScalerChange                                  -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
                        didScale                                            = -1;
                    }
                }
                childVector    += alphabetDimension;
            }

            if (didScale) {
                if (siteCorrectionCounts) {
                    siteCorrectionCounts [siteOrdering.lData[siteID]] += didScale;
                }

                //printf ("NS: site %d node %d \n", siteOrdering.lData[siteID], parentCode);

                if (tcc) {
                    // look ahead to see if we need to correct for downstream cached nodes
                    long cparentTCCIIndex   =   parentTCCIIndex,
                         cparentTCCIBit   =   parentTCCIBit;

                    _Parameter              scM;
                    if (didScale < 0) {
                        scM = _lfScalingFactorThreshold;
                    } else {
                        scM = _lfScalerUpwards;
                    }


                    for (long sid = siteID + 1; sid < siteTo; sid++,cparentTCCIBit++) {
                        if (cparentTCCIBit == _HY_BITMASK_WIDTH_) {
                            cparentTCCIBit   = 0;
                            cparentTCCIIndex ++;
                        }

                        if ((tcc->lData[cparentTCCIIndex] & bitMaskArray.masks[cparentTCCIBit]) > 0) {
                            if (siteCorrectionCounts) {
                                siteCorrectionCounts [siteOrdering.lData[sid]] += didScale;
                            }
                            scalingAdjustments   [parentCode*siteCount + sid] *= scM;
                            localScalerChange                               += didScale * theFilter->theFrequencies [siteOrdering.lData[sid]];
                        } else {
                            break;
                        }
                    }
                }
            }
        }
    }

    // assemble the entire likelihood

    _Parameter _hprestrict_ * rootConditionals = iNodeCache + alphabetDimension * (siteFrom + (flatTree.lLength-1)  * siteCount);
    _Parameter                result = 0.0,
                              correction = 0.0;


    for (long siteID = siteFrom, rootIndex = 0L; siteID < siteTo; siteID++) {
        _Parameter accumulator = 0.;
      
        if (setBranch == flatTree.lLength-1) {
            long                rootState = setBranchTo[siteOrdering.lData[siteID]];
            accumulator         = rootConditionals[rootIndex + rootState] * theProbs[rootState];
            rootIndex           += alphabetDimension;
        } else
            for (long p = 0; p < alphabetDimension; p++,rootIndex++) {
                accumulator += rootConditionals[rootIndex] * theProbs[p];
            }

        /*#pragma omp critical
        {
            if (likeFuncEvalCallCount == 0) {
                eval_buffer [siteID] = accumulator;
            } else {
                if (!CheckEqual (eval_buffer [siteID], accumulator)) {
                    printf ("EVAL %ld, Site %ld/%ld, LogL %g, %g (%d)\n", likeFuncEvalCallCount, siteFrom, siteID, accumulator, eval_buffer [siteID],localScalerChange);
                }
                eval_buffer [siteID] = accumulator;
            }
        }*/

        if (storageVec) {
            storageVec [siteOrdering.lData[siteID]] = accumulator;
        } else {
            if (accumulator <= 0.0) {
                result = -A_LARGE_NUMBER;
                #pragma omp critical
                {
                    ReportWarning (_String("Site ") & (1L+siteOrdering.lData[siteID]) & " evaluated to a 0 probability in ComputeTreeBlockByBranch");
                }
                break;
            }
          
            _Parameter term,
                       temp_sum;
          
            if (theFilter->theFrequencies [siteOrdering.lData[siteID]] > 1) {
                term = log(accumulator) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
            } else {
                term = log(accumulator);
            }
          // Kahan sum
            term     -= correction;
            temp_sum = result + term;
            correction = (temp_sum - result) - term;
            result = temp_sum;
          
        }
    }

    if (!storageVec && localScalerChange) {
        #pragma omp atomic
        overallScaler += localScalerChange;
    }

    return result;
}

#ifdef _SLKP_DEBUG_

/*----------------------------------------------------------------------------------------------------------*/
void echoNodeList (_SimpleList& theNodes, _SimpleList& leaves, _SimpleList& iNodes)
{
    for (long n = 0; n < theNodes.lLength; n++) {
        node<long>* nd = (node<long>*)(theNodes(n)<leaves.lLength?leaves(theNodes(n)):iNodes(theNodes(n)-leaves.lLength));
        printf ("%d %d %s\n", n, theNodes(n), LocateVar(nd->in_object)->GetName()->sData);
    }
}

#endif

/*----------------------------------------------------------------------------------------------------------*/

void            _TheTree::ComputeBranchCache    (
    _SimpleList&            siteOrdering,
    long                    brID,
    _Parameter*         cache,
    _Parameter*         iNodeCache,
    _DataSetFilter*     theFilter,
    long           *        lNodeFlags,
    _Parameter*         scalingAdjustments,
    long        *           siteCorrectionCounts,
    _GrowingVector*     lNodeResolutions,
    long&                   overallScaler,
    long                    siteFrom,
    long                    siteTo,
    long                    catID,
    _SimpleList*            tcc,
    _Parameter*         siteRes
)
{

    //printf ("ComputeBranchCache\n");

    _SimpleList taggedNodes (flatLeaves.lLength + flatNodes.lLength, 0, 0),
                nodesToProcess,
                rootPath;

    long        myParent               = brID       -flatLeaves.lLength,
                alphabetDimension     =            theFilter->GetDimension(),
                alphabetDimensionmod4  =         alphabetDimension - alphabetDimension % 4,
                siteCount               =            theFilter->NumberDistinctSites();

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    do {
        taggedNodes.lData[myParent+flatLeaves.lLength] = 1;
        myParent = flatParents.lData[myParent+flatLeaves.lLength];
    } while (myParent >= 0);


    for (unsigned long k = 0; k <  flatLeaves.lLength+flatNodes.lLength; k++) {
        myParent = flatParents.lData[k];
        if (taggedNodes.lData[myParent+flatLeaves.lLength] == 1 && taggedNodes.lData[k] == 0) {
            if (myParent != brID - flatLeaves.lLength) {
                nodesToProcess << k;
            }
        }
        if (taggedNodes.lData[k]) {
            rootPath << k;
        }
    }

    //printf ("ComputeBranchCache at branch %d; siteOdering %s\n",
    //      brID, _String((_String*)siteOrdering.toStr()).sData);

    //echoNodeList (rootPath,flatLeaves,flatNodes );
    //echoNodeList (nodesToProcess,flatLeaves,flatNodes);

    _Parameter * state = cache + alphabetDimension * siteFrom,
                 * childVector;

    long        localScalerChange = 0;

    // first populate the downward looking vector of conditionals

    if (brID < flatLeaves.lLength) { // a leaf
        for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension) {
            long siteState = lNodeFlags[brID*siteCount + siteOrdering.lData[siteID]] ;
            if (siteState >= 0)
                // a single character state; sweep down the appropriate column
            {
                for (long s = 0; s < alphabetDimension; s++) {
                    state[s] = 0.;
                }
                state[siteState] = 1.;
            } else {
                childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                for (long s = 0; s < alphabetDimension; s++) {
                    state[s] = childVector[s];
                }
            }
        }
    } else { // an internal branch
        long        nodeCode = brID - flatLeaves.lLength;
        _Parameter *lastUpdated = iNodeCache + (nodeCode * siteCount + siteFrom) * alphabetDimension;

        long currentTCCIndex        ,
             currentTCCBit            ;

        if (tcc) {
            currentTCCIndex = siteCount * nodeCode + siteFrom;
            currentTCCBit   = currentTCCIndex % _HY_BITMASK_WIDTH_;
            currentTCCIndex /= _HY_BITMASK_WIDTH_;
        }

        for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension) {
            if (tcc) {
                if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) == 0) {
                    lastUpdated = iNodeCache + (nodeCode * siteCount + siteID) * alphabetDimension;
                }
            }

            for (long s = 0; s < alphabetDimension; s++) {
                state[s] = lastUpdated[s];
            }

            if (tcc) {
                if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                    currentTCCBit   = 0;
                    currentTCCIndex ++;
                }
            } else {
                lastUpdated += alphabetDimension;
            }
        }
    }

    state = cache + alphabetDimension * siteCount;

    taggedNodes.Populate (flatTree.lLength, 0, 0);
    rootPath.Flip ();

    for  (long nodeID = 0; nodeID < nodesToProcess.lLength + rootPath.lLength - 2; nodeID++) {
        bool    notPassedRoot = nodeID<nodesToProcess.lLength;

        long    nodeCode   = notPassedRoot?nodesToProcess.lData [nodeID]:rootPath.lData[nodeID-nodesToProcess.lLength],
                parentCode = notPassedRoot?flatParents.lData [nodeCode]:(rootPath.lData[nodeID-nodesToProcess.lLength+1] - flatLeaves.lLength);

        bool    isLeaf     = nodeCode < flatLeaves.lLength;

        if (!isLeaf) {
            nodeCode -=  flatLeaves.lLength;
        }

        _Parameter * parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedNodes.lData[parentCode] == 0)
            // mark the parent for update and clear its conditionals if needed
        {
            taggedNodes.lData[parentCode]     = 1;
            _Parameter      _hprestrict_ *localScalingFactor      = scalingAdjustments + parentCode*siteCount;
            if (alphabetDimension == 4) {
                long k3     = 0;
                for (long k = siteFrom; k < siteTo; k++, k3+=4) {
                    _Parameter scaler = localScalingFactor[k];
                    parentConditionals [k3]   = scaler;
                    parentConditionals [k3+1] = scaler;
                    parentConditionals [k3+2] = scaler;
                    parentConditionals [k3+3] = scaler;
                }
            } else {
                long k3     = 0;
                for (long k = siteFrom; k < siteTo; k++) {
                    _Parameter scaler = localScalingFactor[k];
                    for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
                        parentConditionals [k3] = scaler;
                    }
                }
            }
        }

        _CalcNode    * currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
                                         ((_CalcNode*) flatTree    (notPassedRoot?nodeCode:parentCode));

        _Parameter  *       _hprestrict_ transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;

        _Parameter  *       childVector,
                    *     lastUpdatedSite;

        if (!isLeaf) {
            lastUpdatedSite = childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
        }

        long currentTCCIndex        ,
             currentTCCBit            ,
             parentTCCIIndex      ,
             parentTCCIBit            ;

        if (tcc) {
            parentTCCIIndex = siteCount * parentCode + siteFrom;
            parentTCCIBit   = parentTCCIIndex % _HY_BITMASK_WIDTH_;
            parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
            if (! isLeaf) {
                currentTCCIndex = siteCount * nodeCode + siteFrom;
                currentTCCBit   = currentTCCIndex % _HY_BITMASK_WIDTH_;
                currentTCCIndex /= _HY_BITMASK_WIDTH_;
            }
        }

        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += alphabetDimension) {
            _Parameter  *tMatrix = transitionMatrix;

            char canScale = !notPassedRoot;

            if (isLeaf) {
                long siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.lData[siteID]] ;
                if (siteState >= 0)
                    // a single character state; sweep down the appropriate column
                {
                    if (alphabetDimension == 4) { // special case for nuc data
                        parentConditionals[0] *= tMatrix[siteState];
                        parentConditionals[1] *= tMatrix[siteState+4];
                        parentConditionals[2] *= tMatrix[siteState+8];
                        parentConditionals[3] *= tMatrix[siteState+12];
                    } else {
                        tMatrix  +=  siteState;
                        for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
                            parentConditionals[k] *= *tMatrix;
                        }

                    }
                    continue;
                } else {
                    childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                }
                canScale = 0;
            } else {
                if (tcc&&notPassedRoot) {
                    if ((tcc->lData[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom)
                        // the value of this conditional vector needs to be copied from a previously stored site
                        // subtree duplication
                        for (long k = 0; k < alphabetDimension; k++) {
                            childVector[k] = lastUpdatedSite[k];
                        }
                    else {
                        lastUpdatedSite = childVector;
                    }

                    if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                        currentTCCBit   = 0;
                        currentTCCIndex ++;
                    }
                    if (++parentTCCIBit == _HY_BITMASK_WIDTH_) {
                        parentTCCIBit   = 0;
                        parentTCCIIndex ++;
                    }
                }
            }

            _Parameter sum      = .0;
            char       didScale = 0;

            if (alphabetDimension == 4) { // special case for nuc data
                _handle4x4_pruning_case (childVector, tMatrix, parentConditionals);

                if (canScale) {
                    sum     = parentConditionals [0] + parentConditionals [1] + parentConditionals [2] + parentConditionals [3];
                    if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                        _Parameter tryScale                                 = scalingAdjustments [nodeCode*siteCount + siteID] * _lfScalerUpwards;
                        if (tryScale < HUGE_VAL) {
                            parentConditionals[0]                             *= _lfScalerUpwards;
                            parentConditionals[1]                             *= _lfScalerUpwards;
                            parentConditionals[2]                             *= _lfScalerUpwards;
                            parentConditionals[3]                             *= _lfScalerUpwards;

                            localScalerChange                                  += theFilter->theFrequencies [siteOrdering.lData[siteID]];
                            didScale                                            = 1;
                        }
                    } else {
                        if (sum > _lfScalerUpwards) {
                            parentConditionals [0]                             *= _lfScalingFactorThreshold;
                            parentConditionals [1]                             *= _lfScalingFactorThreshold;
                            parentConditionals [2]                             *= _lfScalingFactorThreshold;
                            parentConditionals [3]                             *= _lfScalingFactorThreshold;

                            localScalerChange                                  -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
                            didScale                                            = -1;
                        }
                    }
                }
                childVector += 4;
            } else {
#ifdef _SLKP_USE_SSE_INTRINSICS
                    double buffer[2] __attribute__ ((aligned (16)));
#endif
#ifndef _SLKP_SSE_VECTORIZATION_
                for (long p = 0; p < alphabetDimension; p++) {
                    _Parameter      accumulator = 0.0;
#ifdef _SLKP_USE_SSE_INTRINSICS

                        __m128d buffer1,
                                buffer2,
                                buffer3 = _mm_setzero_pd(),
                                buffer4 = _mm_setzero_pd(),
                                load1, 
                                load2,
                                load3,
                                load4;
                                
                        
                        if (((long int)tMatrix & 0x1111b) == 0 && ((long int)childVector & 0x1111b) == 0){ 
                           for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                                load1 = _mm_load_pd (tMatrix+c);
                                load2 = _mm_load_pd (tMatrix+c+2);
                                load3 = _mm_load_pd (childVector+c);
                                load4 = _mm_load_pd (childVector+c+2);
                                buffer1 = _mm_mul_pd (load1, load3);
                                buffer2 = _mm_mul_pd (load2, load4);
                                buffer3 = _mm_add_pd (buffer1,buffer3);
                                buffer4 = _mm_add_pd (buffer2,buffer4);
                            }      
                        } else {
                           for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                                load1 = _mm_loadu_pd (tMatrix+c);
                                load2 = _mm_loadu_pd (tMatrix+c+2);
                                load3 = _mm_loadu_pd (childVector+c);
                                load4 = _mm_loadu_pd (childVector+c+2);
                                buffer1 = _mm_mul_pd (load1, load3);
                                buffer2 = _mm_mul_pd (load2, load4);
                                buffer3 = _mm_add_pd (buffer1,buffer3);
                                buffer4 = _mm_add_pd (buffer2,buffer4);
                            }      
                        
                        }
                        
                        buffer3 = _mm_add_pd (buffer3, buffer4);    
                        _mm_store_pd (buffer, buffer3);
                        accumulator = buffer[0] + buffer[1];
                    
#else
                   for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
                    {
                        _Parameter  pr1 =    tMatrix[c]   * childVector[c],
                                    pr2 =    tMatrix[c+1] * childVector[c+1],
                                    pr3 =    tMatrix[c+2] * childVector[c+2],
                                    pr4 =    tMatrix[c+3] * childVector[c+3];
                         pr1 += pr2;
                         pr3 += pr4;
                         accumulator += pr1+pr3;
                    }
#endif

                    for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
                        accumulator +=  tMatrix[c] * childVector[c];
                    }

                    tMatrix               += alphabetDimension;
                    sum += (parentConditionals[p] *= accumulator);
                }
#else
                for (long p = 0; p < alphabetDimension; p++) {
                    _Parameter      accumulator = 0.0;

                    for (long c = 0; c < alphabetDimension; c++) { // 4 - unroll the loop
                        accumulator +=  tMatrix[c]   * childVector[c];
                    }


                    tMatrix               += alphabetDimension;
                    sum += (parentConditionals[p] *= accumulator);
                }
#endif
                if (canScale) {
                    if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                        _Parameter tryScale                                 = scalingAdjustments [nodeCode*siteCount + siteID] * _lfScalerUpwards;
                        if (tryScale < HUGE_VAL) {
                            for (long c = 0; c < alphabetDimension; c++) {
                                parentConditionals [c] *= _lfScalerUpwards;
                            }

                            localScalerChange                                      += theFilter->theFrequencies [siteOrdering.lData[siteID]];
                            didScale                                                = 1;
                        }
                    } else {
                        if (sum > _lfScalerUpwards) {
                            for (long c = 0; c < alphabetDimension; c++) {
                                parentConditionals [c] *= _lfScalingFactorThreshold;
                            }

                            localScalerChange                                  -= theFilter->theFrequencies [siteOrdering.lData[siteID]];
                            didScale                                            = -1;
                        }
                    }
                }
                childVector    += alphabetDimension;
            }

            if (didScale&&siteCorrectionCounts) {
                siteCorrectionCounts [siteOrdering.lData[siteID]] += didScale;
                siteRes[siteOrdering.lData[siteID]] *= didScale<0?_lfScalingFactorThreshold:_lfScalerUpwards;
            }
        }
    }


    _Parameter _hprestrict_ *rootConditionals   = iNodeCache +  (rootPath.lData[rootPath.lLength-2] - flatLeaves.lLength)  * siteCount * alphabetDimension;

    for (long ii = siteFrom * alphabetDimension; ii < alphabetDimension*siteTo; ii++) {
        state[ii] = rootConditionals[ii];
    }

    if (!siteCorrectionCounts && localScalerChange) {
        #pragma omp atomic
        overallScaler += localScalerChange;
    }
}

/*----------------------------------------------------------------------------------------------------------*/

_Parameter          _TheTree::ComputeLLWithBranchCache (
    _SimpleList&            siteOrdering,
    long                    brID,
    _Parameter*         cache,
    _DataSetFilter*     theFilter,
    long                    siteFrom,
    long                    siteTo,
    long                    catID,
    _Parameter*         storageVec
)
{
    long        alphabetDimension      = theFilter->GetDimension(),
                //alphabetDimensionmod4  = alphabetDimension - alphabetDimension % 4,
                siteCount            =  theFilter->NumberDistinctSites();

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    _Parameter _hprestrict_ *branchConditionals = cache              + siteFrom * alphabetDimension;
    _Parameter _hprestrict_ *rootConditionals   = branchConditionals + siteCount * alphabetDimension;
    _Parameter  result = 0.0,
                correction = 0.0;


    //printf ("ComputeLLWithBranchCache @ %d catID = %d branchID = %d\n", likeFuncEvalCallCount, catID, brID);

    _CalcNode *givenTreeNode = brID < flatLeaves.lLength ? (((_CalcNode**) flatCLeaves.lData)[brID]):
                               (((_CalcNode**) flatTree.lData)   [brID - flatLeaves.lLength]);

    _Parameter  *   _hprestrict_ transitionMatrix = givenTreeNode->GetCompExp(catID)->theData;

    for (long siteID = siteFrom; siteID < siteTo; siteID++) {
        _Parameter accumulator = 0.;
        


        if (alphabetDimension == 4) {
            accumulator =    rootConditionals[0] * theProbs[0] *
                             (branchConditionals[0] *  transitionMatrix[0] + branchConditionals[1] *  transitionMatrix[1] + branchConditionals[2] *  transitionMatrix[2] + branchConditionals[3] *  transitionMatrix[3]) +
                             rootConditionals[1] * theProbs[1] *
                             (branchConditionals[0] *  transitionMatrix[4] + branchConditionals[1] *  transitionMatrix[5] + branchConditionals[2] *  transitionMatrix[6] + branchConditionals[3] *  transitionMatrix[7]) +
                             rootConditionals[2] * theProbs[2] *
                             (branchConditionals[0] *  transitionMatrix[8] + branchConditionals[1] *  transitionMatrix[9] + branchConditionals[2] *  transitionMatrix[10] + branchConditionals[3] *  transitionMatrix[11]) +
                             rootConditionals[3] * theProbs[3] *
                             (branchConditionals[0] *  transitionMatrix[12] + branchConditionals[1] *  transitionMatrix[13] + branchConditionals[2] *  transitionMatrix[14] + branchConditionals[3] *  transitionMatrix[15]);
            rootConditionals += 4;
        } else {
            long       rmx = 0;
            for (long p = 0; p < alphabetDimension; p++,rootConditionals++) {
                _Parameter     r2 = 0.;
                for (long c = 0; c < alphabetDimension; c++, rmx++) {
                    r2 += branchConditionals[c] *  transitionMatrix[rmx];
                }
                accumulator += *rootConditionals * theProbs[p] * r2;
            }

        }

        branchConditionals += alphabetDimension;

        if (storageVec) {
            storageVec [siteOrdering.lData[siteID]] = accumulator;
        } else {
            if (accumulator <= 0.0) {
                result = -A_LARGE_NUMBER;
                #pragma omp critical
                {
                    ReportWarning (_String("Site ") & (1L+siteOrdering.lData[siteID]) & " evaluated to a 0 probability in ComputeLLWithBranchCache");
                }
                break;
            }
            _Parameter term;
            if (theFilter->theFrequencies [siteOrdering.lData[siteID]] > 1) {
                term = log(accumulator) * theFilter->theFrequencies [siteOrdering.lData[siteID]] - correction;
            } else {
                term = log(accumulator) - correction;
            }
            _Parameter temp_sum = result + term;
            correction = (temp_sum - result) - term;
            result = temp_sum;
          
            //result += log(accumulator) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
        }
    }
    return result;
}

/*----------------------------------------------------------------------------------------------------------*/

_Parameter      _TheTree::ComputeTwoSequenceLikelihood
(
    _SimpleList   & siteOrdering,
    _DataSetFilter* theFilter,
    long      *         lNodeFlags,
    _GrowingVector* lNodeResolutions,
    long                siteFrom,
    long                siteTo,
    long                catID,
    _Parameter*     storageVec
)
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
    // process the leaves first

    long            alphabetDimension      =            theFilter->GetDimension(),
                    siteCount            =            theFilter->NumberDistinctSites(),
                    alphabetDimensionmod4  =          alphabetDimension-alphabetDimension%4;

    _CalcNode       *theNode               =            ((_CalcNode*) flatCLeaves (0));
    _Parameter  *   _hprestrict_ transitionMatrix
    =           theNode->GetCompExp(catID)->theData,
    result                 =            0.;

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    for (long siteID = siteFrom; siteID < siteTo; siteID++) {
        _Parameter  *tMatrix = transitionMatrix,
                     sum     = 0.;

        long siteState1 = lNodeFlags[siteOrdering.lData[siteID]],
             siteState2 = lNodeFlags[siteCount + siteOrdering.lData[siteID]];
        
        if (siteState1 >= 0)
            // a single character state; sweep down the appropriate column
        {
            if (siteState2 >= 0) { // both completely resolved;
                sum = tMatrix[siteState1*alphabetDimension + siteState2];
            } else { // first resolved, second is not
                _Parameter* childVector = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;
                tMatrix   +=  siteState1*alphabetDimension;
                if (alphabetDimension == 4) { // special case for nuc data
                    sum = tMatrix[0] * childVector[0] +
                          tMatrix[1] * childVector[1] +
                          tMatrix[2] * childVector[2] +
                          tMatrix[3] * childVector[3];
                } else {
                    for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
                        sum +=  tMatrix[c]   * childVector[c] +
                                tMatrix[c+1] * childVector[c+1] +
                                tMatrix[c+2] * childVector[c+2] +
                                tMatrix[c+3] * childVector[c+3];

                    for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
                        sum +=  tMatrix[c] * childVector[c];
                    }
                }
            }
            sum *= theProbs[siteState1];
        } else {
            if (siteState2 >=0 ) { // second resolved, but not the first
               _Parameter* childVector = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension;
                tMatrix                +=  siteState2;
                if (alphabetDimension == 4) { // special case for nuc data
                    sum = tMatrix[0] * childVector[0]  * theProbs[0]+
                          tMatrix[4] * childVector[1]  * theProbs[1]+
                          tMatrix[8] * childVector[2]  * theProbs[2]+
                          tMatrix[12] * childVector[3] * theProbs[3];

                } else {
                    for (long c = 0; c < alphabetDimensionmod4; c+=4, tMatrix += 4*alphabetDimension) // 4 - unroll the loop
                        sum +=  tMatrix[0]   * childVector[c] * theProbs[c]+
                                tMatrix[alphabetDimension] * childVector[c+1] * theProbs[c+1]+
                                tMatrix[alphabetDimension+alphabetDimension] * childVector[c+2] * theProbs[c+2]+
                                tMatrix[alphabetDimension+alphabetDimension+alphabetDimension] * childVector[c+3] * theProbs[c+3];

                    for (long c = alphabetDimensionmod4; c < alphabetDimension; c++, tMatrix += alphabetDimension) {
                        sum +=  tMatrix[0] * childVector[c] * theProbs[c];
                    }
                }
            } else
                // both unresolved
            {
                _Parameter *childVector1 = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension,
                            *childVector2 = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;

                if (alphabetDimension == 4) { // special case for nuc data
                    sum = (tMatrix[0] * childVector2[0] + tMatrix[1] * childVector2[1] + tMatrix[2] * childVector2[2] + tMatrix[3] * childVector2[3])     * childVector1[0] * theProbs[0]+
                          (tMatrix[4] * childVector2[0] + tMatrix[5] * childVector2[1] + tMatrix[6] * childVector2[2] + tMatrix[7] * childVector2[3])     * childVector1[1] * theProbs[1]+
                          (tMatrix[8] * childVector2[0] + tMatrix[9] * childVector2[1] + tMatrix[10] * childVector2[2] + tMatrix[11] * childVector2[3])   * childVector1[2] * theProbs[2] +
                          (tMatrix[12] * childVector2[0] + tMatrix[13] * childVector2[1] + tMatrix[14] * childVector2[2] + tMatrix[15] * childVector2[3]) * childVector1[3] * theProbs[3];

                } else {
                    for (long r = 0; r < alphabetDimension; r++) { // 4 - unroll the loop
                        if (childVector1[r] > 0.0) {
                            _Parameter sum2 = 0.0;
                            for (long c = 0; c < alphabetDimensionmod4; c+=4, tMatrix += 4) // 4 - unroll the loop
                                sum2   +=  tMatrix[0] * childVector2[c] +
                                           tMatrix[1] * childVector2[c+1]+
                                           tMatrix[2] * childVector2[c+2]+
                                           tMatrix[3] * childVector2[c+3];

                            for (long c = alphabetDimensionmod4; c < alphabetDimension; c++, tMatrix ++) {
                                sum2 +=  tMatrix[0] * childVector2[c];
                            }

                            sum += sum2 * childVector1[r] * theProbs[r];
                        } else {
                            tMatrix += alphabetDimension;
                        }
                    }
                }
            }
        }
        if (storageVec) {
            storageVec [siteOrdering.lData[siteID]] = sum;
        } else {
            if (sum <= 0.0) {
                return -A_LARGE_NUMBER;
            } else {
                //printf ("%d: %g\n", siteID, sum);
                result += log(sum) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
            }
        }
    }

    return result;
}

//_______________________________________________________________________________________________

void     _TheTree::SampleAncestorsBySequence (_DataSetFilter* dsf, _SimpleList& siteOrdering, node<long>* currentNode, _AVLListX* nodeToIndex, _Parameter* iNodeCache,
        _List& result, _SimpleList* parentStates, _List& expandedSiteMap, _Parameter* catAssignments, long catCount)

// must be called initially with the root node


// dsf:                         the filter to sample from
// siteOrdering:                the map from cache ordering to actual pattern ordering
// currentNode:                 the node index to sample for
// nodeToIndex:                 an AVL that maps the address of an internal node pointed to by node<long> to its order in the tree postorder traversal
// iNodeCache:                  internal node likelihood caches
// results:                     the list that will store sampled strings
// parentStates:                sampled states for the parent of the current node
// expandedSiteMap:             a list of simple lists giving site indices for each unique column pattern in the alignment
// catAssignments:              a vector assigning a (partition specific) rate category to each site (nil if no rate variation)
// catCount:                    the number of rate classes

// this needs to be updated to deal with traversal caches!
{
    long                      childrenCount     = currentNode->get_num_nodes();

    if (childrenCount) {
        long            siteCount                       = dsf->NumberDistinctSites  (),
                        alphabetDimension              = dsf->GetDimension         (),
                        nodeIndex                       = nodeToIndex->GetXtra (nodeToIndex->Find ((BaseRef)currentNode)),
                        unitLength                     = dsf->GetUnitLength(),
                        catBlockShifter                    = catAssignments?(dsf->NumberDistinctSites()*GetINodeCount()):0;


        _CalcNode *     currentTreeNode = ((_CalcNode*) flatTree (nodeIndex));
        _SimpleList     sampledStates     (dsf->GetSiteCount (), 0, 0);

        _Parameter  *       _hprestrict_ transitionMatrix = (catAssignments|| !parentStates)?nil:currentTreeNode->GetCompExp()->theData;
        _Parameter  *       _hprestrict_ conditionals     = catAssignments?nil:(iNodeCache + nodeIndex  * siteCount * alphabetDimension);
        _Parameter  *       _hprestrict_ cache            = new _Parameter [alphabetDimension];

        for (long           pattern = 0; pattern < siteCount; pattern++) {
            _SimpleList*    patternMap = (_SimpleList*) expandedSiteMap (siteOrdering.lData[pattern]);
            if (catAssignments) {
                long localCatID = catAssignments[siteOrdering.lData[pattern]];
                if (parentStates) {
                    transitionMatrix = currentTreeNode->GetCompExp(localCatID)->theData;
                }

                conditionals     = iNodeCache + localCatID*alphabetDimension*catBlockShifter + (pattern + nodeIndex  * siteCount) * alphabetDimension;
            }

            for (long site = 0; site < patternMap->lLength; site++) {
                long        siteID =   patternMap->lData[site];

                _Parameter  randVal  = genrand_real2(),
                            totalSum = 0.;

                _Parameter  *       _hprestrict_  matrixRow;

                if  (parentStates == nil) {
                    matrixRow = theProbs;
                } else {
                    matrixRow = transitionMatrix + parentStates->lData[siteID] * alphabetDimension;
                }

                for (long i = 0; i<alphabetDimension; i++) {
                    totalSum += (cache[i] = matrixRow[i]*conditionals[i]);
                }

                randVal *= totalSum;
                totalSum    = 0.0;
                long        sampledChar = -1;
                while       (totalSum < randVal) {
                    sampledChar ++;
                    totalSum += cache[sampledChar];
                }

                sampledStates.lData[siteID] = sampledChar;
            }

            if (catAssignments == nil) {
                conditionals += alphabetDimension;
            }
        }

        delete [] cache;

        _SimpleList  conversion;
        _AVLListXL   conversionAVL (&conversion);

        _String * sampledSequence = new _String (siteCount*unitLength, true);
        _String  letterValue (unitLength,false);
        for (long charIndexer = 0; charIndexer < sampledStates.lLength; charIndexer++) {
            dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(sampledStates.lData[charIndexer]), unitLength, letterValue.sData, &conversionAVL);
            (*sampledSequence) << letterValue;
        }
        sampledSequence->Finalize();
        result.AppendNewInstance(sampledSequence);
        //printf ("%d: %s\n", nodeIndex, sampledSequence->sData);

        for (long child = 1; child <= childrenCount; child ++) {
            SampleAncestorsBySequence (dsf, siteOrdering, currentNode->go_down(child), nodeToIndex, iNodeCache, result, &sampledStates, expandedSiteMap, catAssignments, catCount);
        }
    }
}


//_______________________________________________________________________________________________

_List*   _TheTree::RecoverAncestralSequences (_DataSetFilter* dsf,
        _SimpleList& siteOrdering,
        _List& expandedSiteMap,
        _Parameter* iNodeCache,
        _Parameter* catAssignments,
        long catCount,
        long* lNodeFlags,
        _GrowingVector* lNodeResolutions,
        bool              alsoDoLeaves
                                             )


// dsf:                         the filter to sample from
// siteOrdering:                the map from cache ordering to actual pattern ordering
// expandedSiteMap:             a list of simple lists giving site indices for each unique column pattern in the alignment
// iNodeCache:                  internal node likelihood caches
// catAssignments:              a vector assigning a (partition specific) rate category to each site
// catCount:                    the number of rate classes
// alsoDoLeaves:                if true, also return ML reconstruction of observed (or partially observed) sequences

{
    long            patternCount                    = dsf->NumberDistinctSites  (),
                    alphabetDimension                = dsf->GetDimension         (),
                    unitLength                        = dsf->GetUnitLength        (),
                    iNodeCount                        = GetINodeCount             (),
                    leafCount                     = GetLeafCount              (),
                    siteCount                        = dsf->GetSiteCount         (),
                    allNodeCount                    = 0,
                    stateCacheDim                    = (alsoDoLeaves? (iNodeCount + leafCount): (iNodeCount)),
                    *stateCache                        = new long [patternCount*(iNodeCount-1)*alphabetDimension],
    *leafBuffer                     = new long [(alsoDoLeaves?leafCount*patternCount:1)*alphabetDimension];

    // a Patterns x Int-Nodes x CharStates integer table
    // with the best character assignment for node i given that its parent state is j for a given site

    _Parameter          *buffer                         = new _Parameter [alphabetDimension];
    // iNodeCache will be OVERWRITTEN with conditional pair (i,j) conditional likelihoods

    checkPointer    (stateCache);
    checkPointer    (leafBuffer);

    _SimpleList     taggedInternals (iNodeCount, 0, 0),
                    postToIn;

    MapPostOrderToInOderTraversal (postToIn);
    // all nodes except the root

    allNodeCount = iNodeCount + leafCount - 1;

    for  (long nodeID = 0; nodeID < allNodeCount; nodeID++) {
        long    parentCode = flatParents.lData [nodeID],
                nodeCode   = nodeID;

        bool    isLeaf     = nodeID < flatLeaves.lLength;


        if (!isLeaf) {
            nodeCode -=  flatLeaves.lLength;
            AddBranchToForcedRecomputeList (nodeCode);
        }

        _Parameter * parentConditionals = iNodeCache + parentCode * alphabetDimension * patternCount;

        if (taggedInternals.lData[parentCode] == 0)
            // mark the parent for update and clear its conditionals if needed
        {
            taggedInternals.lData[parentCode]     = 1;
            long k3     = 0;
            for (long k = 0; k < patternCount; k++) {
                for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
                    parentConditionals [k3] = 1.;
                }
            }
        }

        _CalcNode *          currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):((_CalcNode*) flatTree    (nodeCode));
        _Parameter  *       _hprestrict_ transitionMatrix = catAssignments?nil:currentTreeNode->GetCompExp()->theData;
        // this will need to be toggled on a per site basis
        _Parameter  *       childVector;

        if (!isLeaf) {
            childVector = iNodeCache + (nodeCode * patternCount) * alphabetDimension;
        }

        for (long siteID = 0; siteID < patternCount; siteID++, parentConditionals += alphabetDimension) {
            if (catAssignments) {
                transitionMatrix = currentTreeNode->GetCompExp(catAssignments[siteOrdering.lData[siteID]])->theData;
            }

            _Parameter  _hprestrict_ *tMatrix = transitionMatrix;
            if (isLeaf) {
                long siteState = lNodeFlags[nodeCode*patternCount + siteOrdering.lData[siteID]] ;
                if (siteState >= 0)
                    // a fully resolved leaf
                {
                    tMatrix  +=  siteState;
                    for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
                        parentConditionals[k] *= *tMatrix;
                    }
                    if (alsoDoLeaves) {
                        for (long k = 0; k < alphabetDimension; k++) {
                            leafBuffer[k] = siteState;
                        }
                        leafBuffer += alphabetDimension;
                    }

                    continue;
                } else
                    // an ambiguous leaf
                {
                    childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                }

            }

            // now repopulate this vector as necessary -- if we are here this means
            // that the subtree below has been completely processed,
            // the i-th cell of childVector contains the likelihood of the _optimal_
            // assignment in the subtree below given that the character at the current
            // node is i.

            // hence, given parent state 'p', we optimize
            // max_i pr (p->i) childVector [i] and store it in the p cell of vector childVector

            _Parameter overallMax                     = 0.0;

            long       *stateBuffer                   = isLeaf?leafBuffer:stateCache;

            // check for degeneracy

            long howManyOnes = 0;
            for (long k = 0; k < alphabetDimension; k++) {
                howManyOnes += childVector[k]==1.;
            }

            if (howManyOnes == alphabetDimension) {
                for (long k = 0; k < alphabetDimension; k++) {
                    stateBuffer[k] = -1;
                }
            } else {
                for (long p = 0; p < alphabetDimension; p++) {
                    _Parameter max_lik = 0.;
                    long       max_idx = 0;

                    for (long c = 0; c < alphabetDimension; c++) {
                        _Parameter thisV = tMatrix[c] * childVector[c];
                        if (thisV > max_lik) {
                            max_lik = thisV;
                            max_idx = c;
                        }
                    }

                    stateBuffer [p] = max_idx;
                    buffer [p]      = max_lik;

                    if (max_lik > overallMax) {
                        overallMax = max_lik;
                    }

                    tMatrix += alphabetDimension;
                }

                if (overallMax > 0.0 && overallMax < _lfScalingFactorThreshold) {
                    for (long k = 0; k < alphabetDimension; k++) {
                        buffer[k] *= _lfScalerUpwards;
                    }
                }

                // buffer[p] now contains the maximum likelihood of the tree
                // from this point forward given that parent state is p
                // and stateBuffer[p] stores the maximimizing assignment
                // for this node

                for (long k = 0; k < alphabetDimension; k++) {
                    long stateValue = stateBuffer[k];
                    if (stateValue >= 0) {
                        parentConditionals[k] *= buffer[k];
                    }
                }
            }

            if (isLeaf) {
                if (alsoDoLeaves) {
                    leafBuffer += alphabetDimension;
                }
            } else {
                stateCache += alphabetDimension;
            }

            childVector += alphabetDimension;
        }
    }

    _List      *result = new _List;
    for (long k = 0; k < stateCacheDim; k++) {
        result->AppendNewInstance (new _String(siteCount*unitLength,false));
    }

    _Parameter   _hprestrict_ * rootConditionals = iNodeCache + alphabetDimension * ((iNodeCount-1)  * patternCount);
    _SimpleList  parentStates (stateCacheDim,0,0),
                 conversion;

    stateCache -= patternCount*(iNodeCount-1)*alphabetDimension;
    if (alsoDoLeaves) {
        leafBuffer -= patternCount*leafCount*alphabetDimension;
    }

    _AVLListXL    conversionAVL (&conversion);
    _String       codeBuffer    (unitLength, false);

    for (long siteID = 0; siteID < patternCount; siteID++, rootConditionals += alphabetDimension) {
        _Parameter max_lik = 0.;
        long       max_idx = 0;

        long howManyOnes = 0;
        for (long k = 0; k < alphabetDimension; k++) {
            howManyOnes += rootConditionals[k]==1.;
        }

        _SimpleList*    patternMap = (_SimpleList*) expandedSiteMap (siteOrdering.lData[siteID]);

        if (howManyOnes != alphabetDimension) {
            for (long c = 0; c < alphabetDimension; c++) {
                _Parameter thisV = theProbs[c] * rootConditionals[c];
                if (thisV > max_lik) {
                    max_lik = thisV;
                    max_idx = c;
                }
            }

            parentStates.lData[iNodeCount-1] = max_idx;
            for  (long nodeID = iNodeCount-2; nodeID >=0 ; nodeID--) {
                long parentState = parentStates.lData[flatParents.lData [nodeID+flatLeaves.lLength]];
                if (parentState == -1) {
                    parentStates.lData[nodeID] = -1;
                } else {
                    parentStates.lData[nodeID] = stateCache[(patternCount*nodeID+siteID)*alphabetDimension + parentState];
                }
            }
            if (alsoDoLeaves)
                for  (long nodeID = 0; nodeID <leafCount ; nodeID++) {
                    long parentState = parentStates.lData[flatParents.lData [nodeID]];
                    if (parentState == -1) {
                        parentStates.lData[nodeID+iNodeCount] = -1;
                    } else {
                        parentStates.lData[nodeID+iNodeCount] = leafBuffer[(patternCount*nodeID+siteID)*alphabetDimension + parentState];
                    }
                }
        } else {
            parentStates.Populate(stateCacheDim,-1,0);
        }

        for  (long nodeID = 0; nodeID < stateCacheDim ; nodeID++) {
            dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(parentStates.lData[nodeID]), unitLength, codeBuffer.sData, &conversionAVL);
            _String  *sequence   = (_String*) (*result)(nodeID<iNodeCount?postToIn.lData[nodeID]:nodeID);

            for (long site = 0; site < patternMap->lLength; site++) {
                char* storeHere = sequence->sData + patternMap->lData[site]*unitLength;
                for (long charS = 0; charS < unitLength; charS ++) {
                    storeHere[charS] = codeBuffer.sData[charS];
                }
            }
        }
    }

    delete [] stateCache;
    delete [] leafBuffer;
    delete [] buffer;

    return result;
}
//_______________________________________________________________________________________________

void     _TheTree::SetupCategoryMapsForNodes (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers)
{
    _CalcNode* iterator = DepthWiseTraversal (true);
    while (iterator) {
        iterator->SetupCategoryMap (containerVariables,classCounter,multipliers);
        iterator = DepthWiseTraversal(false);
    }
}
//_______________________________________________________________________________________________

void     _CalcNode::SetupCategoryMap (_List& containerVariables, _SimpleList& classCounter, _SimpleList& multipliers)
{

    long    totalCategories = classCounter.Element(-1),
            globalCatCount  = containerVariables.lLength-1,
            localCategories = 1,
            catCount        = categoryVariables.lLength-1,
            entriesPerCat   = 2+catCount;

    //for (long k = 0; k<categoryVariables.lLength;k++)
    //    printf ("%ld\n", categoryVariables(k));//, ((_Variable*)categoryVariables(k))->GetName()->sData);
  
    if (catCount<0) {
        remapMyCategories.Clear();
    } else {

        remapMyCategories.Populate (totalCategories*entriesPerCat,0,0);

        _SimpleList     remappedIDs,
                        rateMultiplers (categoryVariables.lLength,1,0),
                        categoryValues (globalCatCount+1,0,0);

        for (long myCatID = 0; myCatID <= catCount; myCatID++) {
            long coordinate = containerVariables.FindPointer(LocateVar(categoryVariables.lData[myCatID]));
            if (coordinate < 0) {
                WarnError ("Internal error in SetupCategoryMap. Please report to spond@ucsd.edu");
            }
            localCategories *= classCounter.lData[coordinate];
            remappedIDs << coordinate;
        }

        for (long myCatID = catCount-1; myCatID >= 0; myCatID--) {
            rateMultiplers.lData[myCatID] = rateMultiplers.lData[myCatID+1]*classCounter.lData[remappedIDs.lData[myCatID+1]];
        }

        for (long currentRateCombo  = 0; currentRateCombo < totalCategories; currentRateCombo++) {
            long copyRateCombo = currentRateCombo;
            for (long variableID = 0; variableID <= globalCatCount; variableID++) {
                categoryValues.lData[variableID] = copyRateCombo / multipliers.lData[variableID];
                copyRateCombo = copyRateCombo%multipliers.lData[variableID];
                //printf ("%d %d %d %d\n", currentRateCombo, variableID, multipliers.lData[variableID], categoryValues.lData[variableID]);
            }

            long localCatID = 0;

            for  (long localVariableID = 0; localVariableID<=catCount; localVariableID++) {
                localCatID += rateMultiplers.lData[localVariableID] * categoryValues.lData[remappedIDs.lData[localVariableID]];
            }

            long offset = currentRateCombo * entriesPerCat;
            remapMyCategories.lData[offset] = localCatID;
           //printf ("[%ld] = %ld (%ld)\n", offset, localCatID, );

            offset++;
            for  (long localVariableID = 0; localVariableID<=catCount; localVariableID++) {
                remapMyCategories[offset++] = categoryValues.lData[remappedIDs.lData[localVariableID]];
            }

        }
    }

  //printf ("Node remap at %s yielded %s\n", GetName()->sData, _String((_String*)remapMyCategories.toStr()).sData);

}

//_______________________________________________________________________________________________

_Parameter   _TheTree::Process3TaxonNumericFilter (_DataSetFilterNumeric* dsf, long catID)
{

    _Parameter *l0 =  dsf->probabilityVectors.theData +
                      dsf->categoryShifter * catID + dsf->theNodeMap.lData[0]*dsf->shifter,
                      *l1 = dsf->probabilityVectors.theData +
                            dsf->categoryShifter * catID + dsf->theNodeMap.lData[1]*dsf->shifter,
                            *l2 = dsf->probabilityVectors.theData +
                                  dsf->categoryShifter * catID + dsf->theNodeMap.lData[2]*dsf->shifter,
                                  * matrix0 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[0]->in_object)))->GetCompExp(catID)->theData,
                                    * matrix1 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[1]->in_object)))->GetCompExp(catID)->theData,
                                      * matrix2 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[2]->in_object)))->GetCompExp(catID)->theData,
                                        overallResult = 0.;

    long        patternCount =  dsf->NumberDistinctSites();

    _Parameter  currentAccumulator = 1.;

    for (long patternIndex = 0; patternIndex < patternCount; patternIndex ++, l0+=4, l1+=4, l2+=4) {
        _Parameter rp0 = l0[0] * matrix0[0]+ l0[1]  * matrix0[1]  + l0[2] * matrix0[2]  + l0[3] * matrix0[3];
        _Parameter rp1 = l0[0] * matrix0[4]+ l0[1]  * matrix0[5]  + l0[2] * matrix0[6]  + l0[3] * matrix0[7];
        _Parameter rp2 = l0[0] * matrix0[8]+ l0[1]  * matrix0[9]  + l0[2] * matrix0[10] + l0[3] * matrix0[11];
        _Parameter rp3 = l0[0] * matrix0[12]+ l0[1] * matrix0[13] + l0[2] * matrix0[14] + l0[3] * matrix0[15];

        rp0 *= l1[0] * matrix1[0] + l1[1] * matrix1[1]  + l1[2] * matrix1[2]  + l1[3] * matrix1[3];
        rp1 *= l1[0] * matrix1[4] + l1[1] * matrix1[5]  + l1[2] * matrix1[6]  + l1[3] * matrix1[7];
        rp2 *= l1[0] * matrix1[8] + l1[1] * matrix1[9]  + l1[2] * matrix1[10] + l1[3] * matrix1[11];
        rp3 *= l1[0] * matrix1[12]+ l1[1] * matrix1[13] + l1[2] * matrix1[14] + l1[3] * matrix1[15];

        rp0 *= l2[0] * matrix2[0] + l2[1] * matrix2[1]  + l2[2] * matrix2[2]  + l2[3] * matrix2[3];
        rp1 *= l2[0] * matrix2[4] + l2[1] * matrix2[5]  + l2[2] * matrix2[6]  + l2[3] * matrix2[7];
        rp2 *= l2[0] * matrix2[8] + l2[1] * matrix2[9]  + l2[2] * matrix2[10] + l2[3] * matrix2[11];
        rp3 *= l2[0] * matrix2[12]+ l2[1] * matrix2[13] + l2[2] * matrix2[14] + l2[3] * matrix2[15];

        _Parameter  result = theProbs[0]*rp0+
                             theProbs[1]*rp1+
                             theProbs[2]*rp2+
                             theProbs[3]*rp3;


        if (result<=0.0) {
            return -A_LARGE_NUMBER;
        }

        long patternFreq = dsf->theFrequencies[patternIndex];
        for  (long freqIterator = 0; freqIterator < patternFreq; freqIterator++) {
            _Parameter tryMultiplication = currentAccumulator*result;
            if (tryMultiplication > 1.e-300) {
                currentAccumulator = tryMultiplication;
            } else {
                overallResult += myLog (currentAccumulator);
                currentAccumulator = result;
            }
        }
    }
    return overallResult + myLog (currentAccumulator);
}

//__________________________________________________________________________________
// Tree cluster comparison functions
//__________________________________________________________________________________

void        _TreeTopology::ComputeClusterTable (_SimpleList& result, _SimpleList& pswRepresentation)
{
    long            leafCount = pswRepresentation.Element(-2),
                    leafCode  = 0,
                    L,R;

    result.Clear    ();
    result.Populate (3*leafCount,-1,0);

    for (long k = 0; k < pswRepresentation.lLength-2; k+=2) {
        if (pswRepresentation.lData[k] < leafCount) { // is a leaf
            R = leafCode++;
        } else {
            long row;
            L = pswRepresentation.lData[k-2*pswRepresentation.lData[k+1]];
            if (k == pswRepresentation.lLength-4 /* root */
                    || pswRepresentation.lData[k+3] == 0) {
                row = R;
            } else {
                row = L;
            }

            result.lData[row*3] = L;
            result.lData[row*3+1] = R;
        }
    }
}
//__________________________________________________________________________________
_String*        _TreeTopology::ConvertFromPSW                       (_AVLListX& nodeMap,_SimpleList& pswRepresentation)
{
    _String * result = new _String (128L, true);
    if (pswRepresentation.lLength > 4) {
        long leafCount = pswRepresentation.Element (-2);
        // traverse backwards
        bool lastLeaf = false;
        _SimpleList     bounds;

        for (long k = pswRepresentation.lLength-4; k>=0; k-=2) {
            if (lastLeaf) {
                (*result) << ',';
            }
            if (pswRepresentation.lData[k] >= leafCount) { //
                (*result) << ')';
                lastLeaf = false;
                bounds   << k-2*pswRepresentation.lData[k+1];
            } else {
                _String nodeName (*(_String*)nodeMap.Retrieve(pswRepresentation.lData[k]));
                nodeName.Flip();
                (*result) << nodeName;
                while (bounds.Element(-1) == k && bounds.lLength) {
                    (*result) << '(';
                    bounds.Pop();
                }
                lastLeaf = true;
            }
        }
    }
    result->Finalize();
    result->Flip();
    return result;
}

//__________________________________________________________________________________
bool        _TreeTopology::ConvertToPSW (_AVLListX& nodeMap, _List* inames, _SimpleList& pswRepresentation, bool reference)
{
    if (reference == false) {
        nodeMap.Clear();
    }

    pswRepresentation.Clear();

    long    leafIndex  = 0,
            iNodeCount = -1,
            level      = 0;

    _String nodeName;

    DepthWiseTLevel (level,&GetRoot());
    _SimpleList levelBuffer;

    while (currentNode) {
        GetNodeName (currentNode,nodeName,false);


        while (levelBuffer.countitems() <= level) {
            levelBuffer << 0;
        }

        if (IsCurrentNodeATip()) {

            pswRepresentation << leafIndex;
            pswRepresentation << 0;
            if (reference) {
                long remapped = nodeMap.Find(&nodeName);
                if (remapped < 0) {
                    return false;
                } else {
                    remapped = nodeMap.GetXtra (remapped);
                    if (remapped >= 0) {
                        pswRepresentation << remapped;
                    } else {
                        return false;
                    }
                }

                leafIndex++;
            } else {
                nodeMap.Insert(nodeName.makeDynamic(), leafIndex++, false);
            }
        } else {
            pswRepresentation << iNodeCount;
            pswRepresentation << levelBuffer.lData[level];
            if (reference) {
                pswRepresentation << 0;
            } else {
                (*inames) && &nodeName;
            }

            iNodeCount--;
        }
        if (level) {
            levelBuffer.lData[level-1] += levelBuffer.lData[level]+1;
        }
        levelBuffer.lData[level]   = 0;
        DepthWiseTLevel  (level,nil);
    }

    for (long k = 0; k < pswRepresentation.lLength; k+=(reference?3:2))
        if (pswRepresentation.lData[k] < 0) {
            pswRepresentation.lData[k] = leafIndex-pswRepresentation.lData[k]-1;
        }

    pswRepresentation << leafIndex;
    pswRepresentation << (-iNodeCount-1);

    return true;
}


//__________________________________________________________________________________

_AssociativeList *   _TreeTopology::SplitsIdentity (_PMathObj p)
// compare tree topologies
{
    _Matrix     * result = (_Matrix*) checkPointer(new _Matrix (2,1,false,true)),
                  * result2 = nil;

    _FString    * treeR  = new _FString();

    _Constant * bc = (_Constant*) BranchCount ();
    result->theData[0] = bc->Value();
    result->theData[1] = -1;


    if (p && (p->ObjectClass() == TOPOLOGY || p->ObjectClass() == TREE)) {
        _List           avlSupport,
                        iNames;

        _AVLListX       nameMap  (&avlSupport);

        _SimpleList     psw, psw2, clusters, inodeList;


        ConvertToPSW            (nameMap, &iNames, psw);
        ComputeClusterTable     (clusters, psw);
        if (((_TreeTopology*)p)->ConvertToPSW           (nameMap, nil, psw2, true)) {
            _SimpleList           workSpace;
            long leafCount      = psw.Element (-2);

            for (unsigned long k = 0; k < psw2.lLength-3; k+=3) {

                if (psw2.lData[k] < leafCount) {
                    workSpace << 1;
                    workSpace << 1;
                    workSpace << psw2.lData[k+2];
                    workSpace << psw2.lData[k+2];
                } else {
                    _SimpleList quad;
                    quad << leafCount+1;
                    quad << 0;
                    quad << 0;
                    quad << 1;

                    long  w = psw2.lData[k+1];
                    while (w > 0) {
                        _SimpleList quad2;
                        quad2 << workSpace.Pop();
                        quad2 << workSpace.Pop();
                        quad2 << workSpace.Pop();
                        quad2 << workSpace.Pop();
                        w -= quad2.lData[3];
                        quad.lData[0] = MIN(quad2.lData[0],quad.lData[0]);
                        quad.lData[1] = MAX(quad2.lData[1],quad.lData[1]);
                        quad.lData[2] += quad2.lData[2];
                        quad.lData[3] += quad2.lData[3];
                    }

                    if (quad.lData[2] == quad.lData[1] - quad.lData[0] + 1) {
                        if (clusters.lData[3*quad.lData[0]] == quad.lData[0] && clusters.lData[3*quad.lData[0]+1] == quad.lData[1]) {
                            clusters.lData[3*quad.lData[0]+2] = 1;
                        } else if (clusters.lData[3*quad.lData[1]] == quad.lData[0] && clusters.lData[3*quad.lData[1]+1] == quad.lData[1]) {
                            clusters.lData[3*quad.lData[1]+2] = 1;
                        }
                    }
                    quad.Flip();
                    workSpace << quad;
                }
            }

            psw2.Clear();
            long matchCount = 0,
                 iNodeCount = 0;

            long L, R;

            _SimpleList leafSpans (leafCount,0,0),
                        iNodesTouched;

            for (unsigned long k = 0; k < psw.lLength-2; k+=2) {
                if (psw.lData[k] < leafCount) {
                    R = psw.lData[k];
                    psw2 << R;
                    psw2 << 0;
                    leafSpans.lData[R] = (psw2.lLength>>1);
                } else {
                    long ll = k-2*psw.lData[k+1];
                    L       = psw.lData[ll];
                    if ((clusters.lData[3*L] == L && clusters.lData[3*L+1] == R && clusters.lData[3*L+2] > 0)
                            || (clusters.lData[3*R] == L && clusters.lData[3*R+1] == R && clusters.lData[3*R+2] > 0)) {
                        L = (psw2.lLength>>1) - leafSpans.lData[L] + 1;
                        psw2 << leafCount+iNodeCount++;
                        psw2 << L;

                        iNodesTouched << psw.lData[k];
                    }
                }
            }

            for (unsigned long k = 0; k < psw2.lLength; k+=2)
                if (psw2.lData[k] < leafCount) {
                    psw2.lData[k+1] = 0;
                } else {
                    matchCount++;
                }

            psw2 << leafCount;
            psw2 << iNodeCount;

            result->theData[0] = psw.Element (-1);
            result->theData[1] = matchCount;


            *treeR->theString  = ConvertFromPSW (nameMap, psw2);

            _List sharedNames;
            for (long k = 0; k < iNodesTouched.lLength; k++) {
                sharedNames << iNames (iNodesTouched(k)-leafCount);
            }

            result2 = new _Matrix (sharedNames);
        }

    }

    DeleteObject (bc);

    _AssociativeList * resultList = new _AssociativeList;
    resultList->MStore ("CLUSTERS", result, false);
    if (result2) {
        resultList->MStore ("NODES",    result2, false);
    }
    resultList->MStore ("CONSENSUS", treeR, false);
    return resultList;
}

//_______________________________________________________________________________________________

_VariableContainer*     _CalcNode::ParentTree(void) {
    _String parentTree = ParentObjectName();
    _VariableContainer * theParent = (_VariableContainer *)FetchVar(LocateVarByName(parentTree));
    if (theParent && theParent->ObjectClass () == TREE) {
        return theParent;
    }
    return nil;
}


#endif
