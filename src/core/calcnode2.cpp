/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
  Art FY Poon    (apoon42@uwo.ca)
  Steven Weaver (sweaver@temple.edu)

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

#include <math.h>
#include <float.h>

#include "calcnode.h"
#include "likefunc.h"
#include "scfg.h"
#include "function_templates.h"
#include "global_things.h"


//#define _UBER_VERBOSE_DUMP_MATRICES
//#define _UBER_VERBOSE_DUMP 27

extern  long likeFuncEvalCallCount,
        matrix_exp_count;

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
                 hyFloat* iNodeCache,
                 long* lNodeFlags,
                 _SimpleList taggedInternals,
                 _Vector* lNodeResolutions);
#endif

hyFloat          _lfScalerUpwards          = pow(2.,200.),
                    _lfScalingFactorThreshold = 1./_lfScalerUpwards,
                    _logLFScaler            = 200. *log(2.);

_Vector      _scalerMultipliers,
                    _scalerDividers;

//hyFloat          eval_buffer [600];


/*----------------------------------------------------------------------------------------------------------*/

inline void _handle4x4_pruning_case (double const* childVector, double const* tMatrix, double* parentConditionals, void* transposed_mx) {
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


  /*
   A1*B1 + A2*B2 + A3*B3 + A4*B4, where A4 = 1-A1-A2-A3 can be done with three multiplications
   and 3 extra additions, like

   A1*(B1-B4) + A2*(B2-B4) + A3*(B3-B4) + B4

   */

#elif defined _SLKP_USE_AVX_INTRINSICS


      __m256d c3     = _mm256_set1_pd(childVector[3]),
              c0     = _mm256_sub_pd(_mm256_set1_pd(childVector[0]),c3),
              c1     = _mm256_sub_pd(_mm256_set1_pd(childVector[1]),c3),
              c2     = _mm256_sub_pd(_mm256_set1_pd(childVector[2]),c3),
              t0,t1,t2;

          if (transposed_mx) {
             t0    = ((__m256d*)transposed_mx)[0];
             t1    = ((__m256d*)transposed_mx)[1];
             t2    = ((__m256d*)transposed_mx)[2];

          } else {
            t0     = (__m256d) {tMatrix[0],tMatrix[4],tMatrix[8],tMatrix[12]};
            t1     = (__m256d) {tMatrix[1],tMatrix[5],tMatrix[9],tMatrix[13]};
            t2     = (__m256d) {tMatrix[2],tMatrix[6],tMatrix[10],tMatrix[14]};
          }

        // load transition matrix by column

        __m256d sum01 = _mm256_add_pd (_mm256_mul_pd(c0,t0),_mm256_mul_pd(c1,t1)),
                sum23 = _mm256_add_pd (_mm256_mul_pd(c2,t2), c3);

        _mm256_storeu_pd(parentConditionals, _mm256_mul_pd (_mm256_loadu_pd (parentConditionals), _mm256_add_pd (sum01, sum23)));



#else
        // 12 multiplications, 16 additions, 3 subtractions

        hyFloat t1 = childVector[0] - childVector[3],
                   t2 = childVector[1] - childVector[3],
                   t3 = childVector[2] - childVector[3],
                   t4 = childVector[3];

        parentConditionals [0] *= (tMatrix[0]  * t1 + tMatrix[1] * t2) + (tMatrix[2] * t3 + t4);
        parentConditionals [1] *= (tMatrix[4]  * t1 + tMatrix[5] * t2) + (tMatrix[6] * t3 + t4);
        parentConditionals [2] *= (tMatrix[8]  * t1 + tMatrix[9] * t2) + (tMatrix[10] * t3 + t4);
        parentConditionals [3] *= (tMatrix[12] * t1 + tMatrix[13] * t2) + (tMatrix[14] * t3 + t4);
#endif

}


/*----------------------------------------------------------------------------------------------------------*/

hyFloat  acquireScalerMultiplier (long s)
{
    if (s>0) {
        if (s >= _scalerMultipliers.get_used())
            for (long k = _scalerMultipliers.get_used(); k <= s; k++) {
                _scalerMultipliers.Store (exp (-_logLFScaler * k));
            }
        return _scalerMultipliers.theData[s];
    }
    s = -s;
    if (s >= _scalerDividers.get_used())
        for (long k = _scalerDividers.get_used(); k <= s; k++) {
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
            printf ("NodeID %ld (%s). Old length %ld, new length %ld (%ld)\n", nodeID, thisNode->GetName()->sData, didIncrease,matrixQueue.lLength, isExplicitForm.lLength);
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
    matrix_exp_count += matrixQueue.lLength;
#endif

    #pragma omp parallel for default(shared) private (matrixID) schedule(static) if (nt>1)  num_threads (nt)
    for  (matrixID = 0; matrixID < matrixQueue.lLength; matrixID++) {
        if (isExplicitForm.lData[matrixID] == 0 || !hasExpForm) { // normal matrix to exponentiate
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

        if (currentTreeNode->NeedNewCategoryExponential (catID)) {
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

void        _TheTree::FillInConditionals        (_DataSetFilter const*        theFilter, hyFloat*  iNodeCache,  _SimpleList*   tcc)
// this utility function will simply fill in all the conditional probability vectors for internal nodes,
// including those that were skipped due to column sorting optimization
// this is useful to avoid code duplication for other functions (e.g. ancestral sampling) that
// make use of conditional probability vectors, but would not benefit from subtree caching
{
    if (!tcc) {
        return;
    }

    long            alphabetDimension     =         theFilter->GetDimension(),
                    siteCount           =         theFilter->GetPatternCount();

    for  (long nodeID = 0; nodeID < flatTree.lLength; nodeID++) {
        hyFloat * conditionals       = iNodeCache +(nodeID  * siteCount) * alphabetDimension;
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
hyFloat _TheTree::OCLLikelihoodEvaluator (			_SimpleList&		     updateNodes,
														_DataSetFilter*		 theFilter,
                                                        hyFloat*			 iNodeCache,
                                                        long	   *			 lNodeFlags,
                                                        _Vector*		 lNodeResolutions,
														_OCLEvaluator& OCLEval)
{

	_SimpleList		taggedInternals					(flatNodes.lLength, 0, 0);
	//printf("Launching a tree in OpenCL...\n");
	return ((hyFloat) OCLEval.launchmdsocl(			updateNodes,
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

hyFloat  _TheTree::VerySimpleLikelihoodEvaluator   (_SimpleList&          updateNodes,
        _DataSetFilter*      theFilter,
        hyFloat*          iNodeCache,
        long       *             lNodeFlags,
        _Vector*      lNodeResolutions)
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

                    siteCount           =         theFilter->GetPatternCount();
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

        hyFloat * parentConditionals = iNodeCache +  (parentCode  * siteCount) * alphabetDimension;
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



        hyFloat  *       tMatrix = (isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
                                       ((_CalcNode*) flatTree    (nodeCode)))->GetCompExp(0)->theData;
        // in the host code the transition matrix is retrieved from the _CalcNode object
        // in the devide code there will probably be another array to grab it from

        hyFloat  *       childVector; // the vector of conditional probabilities for
        // THIS node (nodeCode)

        if (!isLeaf) {
            childVector = iNodeCache + (nodeCode * siteCount) * alphabetDimension;
        }
        // if THIS node is internal, simply look up the conditional vector in the same cache
        // as the parent (just a different offset)

        for (long siteID = 0; siteID < siteCount; siteID++,
                parentConditionals += alphabetDimension) {
            hyFloat  sum      = 0.0;

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
                // look up the resolution for the ambiguous node -- this will have to be on the device as well
                // but can be in constant memory
            }

            hyFloat *matrixPointer = tMatrix;

            for (long p = 0; p < alphabetDimension; p++) {
                hyFloat      accumulator = 0.0;

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

    hyFloat  * rootConditionals = iNodeCache + alphabetDimension * ((flatTree.lLength-1)  * siteCount),
                  // the root is always the LAST internal node in all lists
                  result = 0.0;


    for (long siteID = 0; siteID < siteCount; siteID++) {
        hyFloat accumulator = 0.;
        for (long p = 0; p < alphabetDimension; p++,rootConditionals++) {
            accumulator += *rootConditionals * theProbs[p];
        }
        /* theProbs is a member variable of the tree, which basically determines
           what probability there is to observe a given character at the root
           in simple cases it is fixed for the duration of optimization, but for more
           complex models it may change from iteration to iteration
         */
        result += log(accumulator) * theFilter->theFrequencies.get (siteID);
        // correct for the fact that identical alignment columns may appear more than once
    }


    return result;
}


/*----------------------------------------------------------------------------------------------------------*/

hyFloat      _TheTree::ComputeTreeBlockByBranch  (                   _SimpleList&        siteOrdering,
        _SimpleList&        updateNodes,
        _SimpleList*        tcc,
        _DataSetFilter const*     theFilter,
        hyFloat*         iNodeCache,
        long      *         lNodeFlags,
        hyFloat*         scalingAdjustments,
        _Vector*     lNodeResolutions,
        long&               overallScaler,
        long                siteFrom,
        long                siteTo,
        long                catID,
        hyFloat*         storageVec,
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
                    siteCount           =         theFilter->GetPatternCount(),
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

        hyFloat * parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedInternals.lData[parentCode] == 0)
            // mark the parent for update and clear its conditionals if needed
        {
            taggedInternals.lData[parentCode]     = 1;
            hyFloat      _hprestrict_ *localScalingFactor      = scalingAdjustments + parentCode*siteCount;

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
                        hyFloat scaler = localScalingFactor[k];
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
                        hyFloat scaler = localScalingFactor[k];
                        for (long k2 = 0; k2 < alphabetDimension; k2++, k3++) {
                            parentConditionals [k3] = scaler;
                        }
                    }
                }
            }
        }

        currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):
                          ((_CalcNode*) flatTree    (nodeCode));

        hyFloat  const * transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
        hyFloat  *       childVector,
                    *       lastUpdatedSite;

#ifdef _SLKP_USE_AVX_INTRINSICS
      __m256d tmatrix_transpose [4] = {
        (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
        (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
        (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
        (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
      };
#endif

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

            hyFloat  const *tMatrix = transitionMatrix;
            hyFloat  sum         = 0.0;

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

#ifdef _SLKP_USE_AVX_INTRINSICS
              _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, tmatrix_transpose);
#else
              _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, nil);
#endif
                // handle scaling if necessary
                // the check for sum > 0.0 is necessary for 'inadmissible' log-L functions (-infinity)
                // for example if a change must happen on a zero-branch length

                sum     = (parentConditionals [0] + parentConditionals [1]) + (parentConditionals [2] + parentConditionals [3]);
                if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                    hyFloat tryScale                                 = scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards;
                    if (tryScale < HUGE_VAL) {
                        parentConditionals [0]                             *= _lfScalerUpwards;
                        parentConditionals [1]                             *= _lfScalerUpwards;
                        parentConditionals [2]                             *= _lfScalerUpwards;
                        parentConditionals [3]                             *= _lfScalerUpwards;
                        localScalerChange                                  += theFilter->theFrequencies.get (siteOrdering.lData[siteID]);
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
                        localScalerChange                                  -= theFilter->theFrequencies.get (siteOrdering.lData[siteID]);
                        didScale                                            = -1;
                    }
                }

                childVector += 4;
            } else {
                hyFloat sum = 0.0;

                if (alphabetDimension > alphabetDimensionmod4){


                    for (long p = 0L; p < alphabetDimension; p++) {
                         hyFloat      accumulator = 0.0;

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
                        double buffer[2] __attribute__ ((aligned (16)));
                        _mm_store_pd (buffer, buffer3);
                        accumulator = buffer[0] + buffer[1];

#elif defined _SLKP_USE_AVX_INTRINSICS // end _SLKP_USE_SSE_INTRINSICS

                        __m256d sum256 = _mm256_setzero_pd();

                        /*if (((long int)tMatrix & 0x11111b) == 0 && ((long int)childVector & 0x11111b) == 0){
                            for (long c = 0L; c < alphabetDimensionmod4; c+=4L) {
                                sum = _mm256_add_pd (sum,_mm256_mul_pd (_mm256_load_pd (tMatrix+c), _mm256_load_pd (childVector+c)));
                             }
                        } else {*/
                            for (long c = 0; c < alphabetDimensionmod4; c+=4L) {
                              __m256d matrix_quad = _mm256_loadu_pd (tMatrix+c),
                                      child_quad = _mm256_loadu_pd (childVector+c),
                                      prod = _mm256_mul_pd (matrix_quad, child_quad);

                                sum256 = _mm256_add_pd (sum256,prod);
                            }

                        //}
                        //double buffer[4] __attribute__ ((aligned (32)));
                        //_mm256_store_pd(buffer, sum256);

                        //accumulator = (buffer[0] + buffer[1]) + (buffer[2] + buffer[3]);

                        accumulator = _avx_sum_4(sum256);
                        //NOT sure why copy to doubles and add is faster
                        // that AVX instructions
#else // _SLKP_USE_AVX_INTRINSICS
                        for (long c = 0; c < alphabetDimensionmod4; c+=4) {
                        // 4 - unroll the loop
                             hyFloat pr1 =    tMatrix[c]   * childVector[c],
                                        pr2 =    tMatrix[c+1] * childVector[c+1],
                                        pr3 =    tMatrix[c+2] * childVector[c+2],
                                        pr4 =    tMatrix[c+3] * childVector[c+3];
                            pr1 += pr2;
                            pr3 += pr4;
                            accumulator += pr1+pr3;
                        }
#endif // regular code

                        for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
                            accumulator +=  tMatrix[c] * childVector[c];
                        }

                        tMatrix               += alphabetDimension;
                        sum += (parentConditionals[p] *= accumulator);
                    }
                }
                else {
                    #ifdef _SLKP_USE_AVX_INTRINSICS
                    if (alphabetDimension == 20UL) {
                      for (long p = 0; p < alphabetDimension; p++) {
                              __m256d t_matrix[5] = {_mm256_loadu_pd(tMatrix),
                                _mm256_loadu_pd(tMatrix+4UL),
                                _mm256_loadu_pd(tMatrix+8UL),
                                _mm256_loadu_pd(tMatrix+12UL),
                                _mm256_loadu_pd(tMatrix+16UL)},

                              c_vector[5] = {_mm256_loadu_pd(childVector),
                                _mm256_loadu_pd(childVector+4UL),
                                _mm256_loadu_pd(childVector+8UL),
                                _mm256_loadu_pd(childVector+12UL),
                                _mm256_loadu_pd(childVector+16UL)};

                              t_matrix[0] = _mm256_mul_pd(t_matrix[0], c_vector[0]);
                              t_matrix[1] = _mm256_mul_pd(t_matrix[1], c_vector[1]);
                              t_matrix[2] = _mm256_mul_pd(t_matrix[2], c_vector[2]);
                              t_matrix[3] = _mm256_mul_pd(t_matrix[3], c_vector[3]);
                              t_matrix[4] = _mm256_mul_pd(t_matrix[4], c_vector[4]);

                              t_matrix[0] = _mm256_add_pd (t_matrix[0],t_matrix[1]);
                              t_matrix[2] = _mm256_add_pd (t_matrix[2],t_matrix[3]);
                              t_matrix[0] = _mm256_add_pd (t_matrix[0],t_matrix[2]);

                              tMatrix               += 20UL;
                              sum += (parentConditionals[p] *= _avx_sum_4(_mm256_add_pd (t_matrix[0],t_matrix[4])));
                           }
                      } else
                     #endif // _SLKP_USE_AVX_INTRINSICS

                    for (long p = 0; p < alphabetDimension; p++) {
                        hyFloat      accumulator = 0.0;

                        for (long c = 0; c < alphabetDimensionmod4; c+=4) // 4 - unroll the loop
                            accumulator +=  tMatrix[c]   * childVector[c] +
                                            tMatrix[c+1] * childVector[c+1] +
                                            tMatrix[c+2] * childVector[c+2] +
                                            tMatrix[c+3] * childVector[c+3];

                        tMatrix               += alphabetDimension;
                        sum += (parentConditionals[p] *= accumulator);
                    }
                }



                if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                    hyFloat tryScale                                 = scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards;
                    if (tryScale < HUGE_VAL) {
                        scalingAdjustments [parentCode*siteCount + siteID] = tryScale;
                        for (long c = 0; c < alphabetDimension; c++) {
                            parentConditionals [c] *= _lfScalerUpwards;
                        }

                        localScalerChange                                      += theFilter->theFrequencies.get(siteOrdering.lData[siteID]);
                        didScale                                                = 1;
                    }
                } else {
                    if (sum > _lfScalerUpwards) {
                        scalingAdjustments [parentCode*siteCount + siteID] *= _lfScalingFactorThreshold;
                        for (long c = 0; c < alphabetDimension; c++) {
                            parentConditionals [c] *= _lfScalingFactorThreshold;
                        }
                        localScalerChange                                  -= theFilter->theFrequencies.get (siteOrdering.lData[siteID]);
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

                    hyFloat              scM;
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
                            localScalerChange                               += didScale * theFilter->theFrequencies (siteOrdering.lData[sid]);
                        } else {
                            break;
                        }
                    }
                }
            }
        }
    }

    // assemble the entire likelihood

    hyFloat _hprestrict_ * rootConditionals = iNodeCache + alphabetDimension * (siteFrom + (flatTree.lLength-1)  * siteCount);
    hyFloat                result = 0.0,
                              correction = 0.0;


    for (long siteID = siteFrom, rootIndex = 0L; siteID < siteTo; siteID++) {
        hyFloat accumulator = 0.;

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
                    hy_global::ReportWarning (_String("Site ") & (1L+siteOrdering.lData[siteID]) & " evaluated to a 0 probability in ComputeTreeBlockByBranch");
                }
                break;
            }

            hyFloat term = log(accumulator),
                       temp_sum;

            long       site_frequency = theFilter->theFrequencies (siteOrdering.lData[siteID]);

            if (site_frequency > 1L) {
                term *= site_frequency;
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
        printf ("%ld %ld %s\n", n, theNodes(n), LocateVar(nd->in_object)->GetName()->sData);
    }
}

#endif

/*----------------------------------------------------------------------------------------------------------*/

void            _TheTree::ComputeBranchCache    (
    _SimpleList&            siteOrdering,
    long                    brID,
    hyFloat*         cache,
    hyFloat*         iNodeCache,
    _DataSetFilter const*     theFilter,
    long           *        lNodeFlags,
    hyFloat*         scalingAdjustments,
    long        *           siteCorrectionCounts,
    _Vector*     lNodeResolutions,
    long&                   overallScaler,
    long                    siteFrom,
    long                    siteTo,
    long                    catID,
    _SimpleList*            tcc,
    hyFloat*         siteRes
)
{


    //printf ("ComputeBranchCache\n");

    _SimpleList taggedNodes (flatLeaves.lLength + flatNodes.lLength, 0, 0),
                nodesToProcess,
                rootPath;

    long        myParent               = brID       -flatLeaves.lLength,
                alphabetDimension     =            theFilter->GetDimension(),
                alphabetDimensionmod4  =         alphabetDimension - alphabetDimension % 4,
                siteCount               =            theFilter->GetPatternCount();

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    do {
        taggedNodes.lData[myParent+flatLeaves.lLength] = 1;
        myParent = flatParents.lData[myParent+flatLeaves.lLength];
    } while (myParent >= 0);


    for (unsigned long k = 0UL; k <  flatLeaves.lLength+flatNodes.lLength; k++) {
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

    /*printf ("\n\nComputeBranchCache at branch %ld; siteOdering %s\n",
            brID, _String((_String*)siteOrdering.toStr()).sData);

    echoNodeList (rootPath,flatLeaves,flatNodes );
    printf ("\n");
    echoNodeList (nodesToProcess,flatLeaves,flatNodes);
    */

    hyFloat * state = cache + alphabetDimension * siteFrom,
                 * childVector;

    long        localScalerChange = 0;

    // first populate the downward looking vector of conditionals

    if (brID < flatLeaves.lLength) { // a leaf
        for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension) {
            long siteState = lNodeFlags[brID*siteCount + siteOrdering.lData[siteID]] ;
            if (siteState >= 0) {
                // a single character state; sweep down the appropriate column
                for (unsigned long s = 0UL; s < alphabetDimension; s++) {
                    state[s] = 0.;
                }
                state[siteState] = 1.;
            } else {
                childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                for (unsigned long s = 0UL; s < alphabetDimension; s++) {
                    state[s] = childVector[s];
                }
            }
        }
    } else { // an internal branch
        long        nodeCode = brID - flatLeaves.lLength;
        hyFloat *lastUpdated = iNodeCache + (nodeCode * siteCount + siteFrom) * alphabetDimension;

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

        hyFloat * parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedNodes.lData[parentCode] == 0)
            // mark the parent for update and clear its conditionals if needed
        {
            //printf ("Resetting parentCode = %ld\n", parentCode);
            taggedNodes.lData[parentCode]     = 1;
            hyFloat     const *localScalingFactor      = scalingAdjustments + parentCode*siteCount;
            if (alphabetDimension == 4L) {
                long k3     = 0;
                for (long k = siteFrom; k < siteTo; k++, k3+=4) {
                    hyFloat scaler = localScalingFactor[k];
                    parentConditionals [k3]   = scaler;
                    parentConditionals [k3+1] = scaler;
                    parentConditionals [k3+2] = scaler;
                    parentConditionals [k3+3] = scaler;
                }
            } else {
                unsigned long k3     = 0UL;
                for (long k = siteFrom; k < siteTo; k++) {
                    hyFloat scaler = localScalingFactor[k];
                    for (unsigned long k2 = 0UL; k2 < alphabetDimension; k2++, k3++) {
                        parentConditionals [k3] = scaler;
                    }
                }
            }
        }

        _CalcNode    * currentTreeNode = (_CalcNode*) (isLeaf?  flatCLeaves (nodeCode):
                                                       flatTree    (notPassedRoot?nodeCode:parentCode));

        //printf ("isLeaf = %d, nodeCode = %ld, parentCode = %ld, matrix from %s, parent name %s\n", isLeaf, nodeCode, parentCode, currentTreeNode->GetName()->sData, ((_CalcNode    *)flatTree(parentCode))->GetName()->sData);

        hyFloat  const * _hprestrict_ transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;

        #ifdef _SLKP_USE_AVX_INTRINSICS
              __m256d tmatrix_transpose [4] = {
                (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
              };
        #endif

       hyFloat  *       childVector,
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
            hyFloat  const *tMatrix = transitionMatrix;

            char canScale = !notPassedRoot;

            if (isLeaf) {
                long siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.lData[siteID]] ;
                if (siteState >= 0) {
                    // a single character state; sweep down the appropriate column
                     if (alphabetDimension == 4) { // special case for nuc data
                        parentConditionals[0] *= tMatrix[siteState];
                        parentConditionals[1] *= tMatrix[siteState+4];
                        parentConditionals[2] *= tMatrix[siteState+8];
                        parentConditionals[3] *= tMatrix[siteState+12];
                    } else {
                        tMatrix  +=  siteState;
                        for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
                          parentConditionals[k] *= *tMatrix;
                            //printf ("Leaf %ld %g %g\n", k, parentConditionals[k], *tMatrix);
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

            hyFloat sum      = .0;
            char       didScale = 0;

            if (alphabetDimension == 4) { // special case for nuc data
                #ifdef _SLKP_USE_AVX_INTRINSICS
                      _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, tmatrix_transpose);
                #else
                      _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, nil);
                #endif

                if (canScale) {
                    sum     = (parentConditionals [0] + parentConditionals [1]) + (parentConditionals [2] + parentConditionals [3]);
                    if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                        hyFloat tryScale                                 = scalingAdjustments [nodeCode*siteCount + siteID] * _lfScalerUpwards;
                        if (tryScale < HUGE_VAL) {
                            parentConditionals[0]                             *= _lfScalerUpwards;
                            parentConditionals[1]                             *= _lfScalerUpwards;
                            parentConditionals[2]                             *= _lfScalerUpwards;
                            parentConditionals[3]                             *= _lfScalerUpwards;

                            localScalerChange                                  += theFilter->theFrequencies (siteOrdering.lData[siteID]);
                            didScale                                            = 1;
                        }
                    } else {
                        if (sum > _lfScalerUpwards) {
                            parentConditionals [0]                             *= _lfScalingFactorThreshold;
                            parentConditionals [1]                             *= _lfScalingFactorThreshold;
                            parentConditionals [2]                             *= _lfScalingFactorThreshold;
                            parentConditionals [3]                             *= _lfScalingFactorThreshold;

                            localScalerChange                                  -= theFilter->theFrequencies (siteOrdering.lData[siteID]);
                            didScale                                            = -1;
                        }
                    }
                }
                childVector += 4;
            } else {
                for (long p = 0; p < alphabetDimension; p++) {

#ifdef _SLKP_USE_AVX_INTRINSICS
                  if (alphabetDimension == 20UL) {

                    __m256d t_matrix[5] = {_mm256_loadu_pd(tMatrix),
                                           _mm256_loadu_pd(tMatrix+4UL),
                                           _mm256_loadu_pd(tMatrix+8UL),
                                           _mm256_loadu_pd(tMatrix+12UL),
                                            _mm256_loadu_pd(tMatrix+16UL)},

                            c_vector[5] = {_mm256_loadu_pd(childVector),
                              _mm256_loadu_pd(childVector+4UL),
                              _mm256_loadu_pd(childVector+8UL),
                              _mm256_loadu_pd(childVector+12UL),
                              _mm256_loadu_pd(childVector+16UL)};

                    t_matrix[0] = _mm256_mul_pd(t_matrix[0], c_vector[0]);
                    t_matrix[1] = _mm256_mul_pd(t_matrix[1], c_vector[1]);
                    t_matrix[2] = _mm256_mul_pd(t_matrix[2], c_vector[2]);
                    t_matrix[3] = _mm256_mul_pd(t_matrix[3], c_vector[3]);
                    t_matrix[4] = _mm256_mul_pd(t_matrix[4], c_vector[4]);

                    t_matrix[0] = _mm256_add_pd (t_matrix[0],t_matrix[1]);
                    t_matrix[2] = _mm256_add_pd (t_matrix[2],t_matrix[3]);
                    t_matrix[0] = _mm256_add_pd (t_matrix[0],t_matrix[2]);

                    tMatrix               += 20UL;
                    sum += (parentConditionals[p] *= _avx_sum_4(_mm256_add_pd (t_matrix[0],t_matrix[4])));
                    continue;

                  }
#endif
                  hyFloat      accumulator = 0.0;

                  for (long c = 0; c < alphabetDimensionmod4; c+=4) { // 4 - unroll the loop
                      hyFloat  pr1 =    tMatrix[c]   * childVector[c],
                                  pr2 =    tMatrix[c+1] * childVector[c+1],
                                  pr3 =    tMatrix[c+2] * childVector[c+2],
                                  pr4 =    tMatrix[c+3] * childVector[c+3];
                       pr1 += pr2;
                       pr3 += pr4;
                       accumulator += pr1+pr3;
                  }

                  for (long c = alphabetDimensionmod4; c < alphabetDimension; c++) {
                      accumulator +=  tMatrix[c] * childVector[c];
                  }

                    tMatrix               += alphabetDimension;
                    //printf ("%ld %g %g\n", p, parentConditionals[p], accumulator);
                    sum += (parentConditionals[p] *= accumulator);
                }

                childVector    += alphabetDimension;

                if (canScale) {
                    if (sum < _lfScalingFactorThreshold && sum > 0.0) {
                        hyFloat tryScale                                 = scalingAdjustments [nodeCode*siteCount + siteID] * _lfScalerUpwards;
                        if (tryScale < HUGE_VAL) {
                            //printf ("tryScale < HUGE_VAL\n");
                            for (unsigned long c = 0UL; c < alphabetDimension; c++) {
                                parentConditionals [c] *= _lfScalerUpwards;
                            }

                            localScalerChange                                      += theFilter->theFrequencies.get (siteOrdering.lData[siteID]);
                            didScale                                                = 1;
                        }
                    } else {
                        if (sum > _lfScalerUpwards) {
                            //printf ("sum > _lfScalerUpwards\n");
                            for (unsigned long c = 0UL; c < alphabetDimension; c++) {
                              parentConditionals [c] *= _lfScalingFactorThreshold;
                            }

                            localScalerChange                                  -= theFilter->theFrequencies.get (siteOrdering.lData[siteID]);
                            didScale                                            = -1;
                        }
                    }
                }
            }

            if (didScale&&siteCorrectionCounts) {
                siteCorrectionCounts [siteOrdering.lData[siteID]] += didScale;
                siteRes[siteOrdering.lData[siteID]] *= didScale<0?_lfScalingFactorThreshold:_lfScalerUpwards;
            }
        }
    }



    //printf ("root name %s\n", ((_CalcNode    *)flatTree(rootPath.lData[rootPath.lLength-2] - flatLeaves.lLength))->GetName()->sData);

    hyFloat const _hprestrict_ *rootConditionals   = iNodeCache +  (rootPath.lData[rootPath.lLength-2] - flatLeaves.lLength)  * siteCount * alphabetDimension;

    state = cache + alphabetDimension * siteCount;
    for (unsigned long ii = siteFrom * alphabetDimension; ii < alphabetDimension*siteTo; ii++) {
        state[ii] = rootConditionals[ii];
        //printf ("Root conditional [%ld] = %g, node state [%ld] = %g\n", ii, state[ii], ii, cache[ii]);
    }

    if (!siteCorrectionCounts && localScalerChange) {
        #pragma omp atomic
        overallScaler += localScalerChange;

      //#pragma omp atomic
      // printf ("Rescale in ComputeBranchCache at branch %ld %ld\n", brID, localScalerChange);
    }
}

/*----------------------------------------------------------------------------------------------------------*/

const _CalcNode* _TheTree::GetNodeFromFlatIndex(long index) const {
  return index < flatLeaves.lLength ? (((_CalcNode**) flatCLeaves.lData)[index]):
  (((_CalcNode**) flatTree.lData)   [index - flatLeaves.lLength]);
}

/*----------------------------------------------------------------------------------------------------------*/

hyFloat          _TheTree::ComputeLLWithBranchCache (
    _SimpleList&            siteOrdering,
    long                    brID,
    hyFloat*         cache,
    _DataSetFilter const*     theFilter,
    long                    siteFrom,
    long                    siteTo,
    long                    catID,
    hyFloat*         storageVec
)
{
  auto bookkeeping =  [&siteOrdering, &storageVec, &theFilter] (const long siteID, const hyFloat accumulator, hyFloat& correction, hyFloat& result) ->  void {

    long direct_index = siteOrdering.lData[siteID];

    if (storageVec) {
      storageVec [direct_index] = accumulator;
    } else {
      if (accumulator <= 0.0) {
        throw (1L+direct_index);
      }

      hyFloat term;
      long       site_frequency = theFilter->theFrequencies.get(direct_index);
      if ( site_frequency > 1L) {
        term =  log(accumulator) * site_frequency - correction;
      } else {
          term = log(accumulator) - correction;
      }
      hyFloat temp_sum = result + term;
      correction = (temp_sum - result) - term;
      result = temp_sum;
        //result += log(accumulator) * theFilter->theFrequencies [siteOrdering.lData[siteID]];
    }
  };

  const unsigned long          alphabetDimension      = theFilter->GetDimension(),
                               siteCount              =  theFilter->GetPatternCount();

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    hyFloat const _hprestrict_ *branchConditionals = cache              + siteFrom * alphabetDimension;
    hyFloat const _hprestrict_ *rootConditionals   = branchConditionals + siteCount * alphabetDimension;
    hyFloat  result = 0.0,
                correction = 0.0;


    //printf ("ComputeLLWithBranchCache @ %d catID = %d branchID = %d\n", likeFuncEvalCallCount, catID, brID);

    _CalcNode const *givenTreeNode = GetNodeFromFlatIndex (brID);

    hyFloat  const *   _hprestrict_ transitionMatrix = givenTreeNode->GetCompExp(catID)->theData;


// cases by alphabet dimension

  try {
    switch (alphabetDimension) {
/****

  NUCLEOTIDES

 ****/
        case 4UL: {
          #ifdef _SLKP_USE_AVX_INTRINSICS
                  __m256d tmatrix_transpose [4] = {
                    (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                    (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                    (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                    (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                  };
          #endif

          for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
            hyFloat accumulator = 0.;
             #ifdef _SLKP_USE_AVX_INTRINSICS
               __m256d root_c = _mm256_loadu_pd (rootConditionals),
               probs  = _mm256_loadu_pd (theProbs),
               b_cond0 = _mm256_set1_pd(branchConditionals[0]),
               b_cond1 = _mm256_set1_pd(branchConditionals[1]),
               b_cond2 = _mm256_set1_pd(branchConditionals[2]),
               b_cond3 = _mm256_set1_pd(branchConditionals[3]),
               s01    = _mm256_add_pd ( _mm256_mul_pd (b_cond0, tmatrix_transpose[0]), _mm256_mul_pd (b_cond1, tmatrix_transpose[1])),
               s23    = _mm256_add_pd ( _mm256_mul_pd (b_cond2, tmatrix_transpose[2]), _mm256_mul_pd (b_cond3, tmatrix_transpose[3]));
               accumulator = _avx_sum_4(_mm256_mul_pd (_mm256_mul_pd (root_c, probs), _mm256_add_pd (s01,s23)));

              #else
                  accumulator =    rootConditionals[0] * theProbs[0] *
                  (branchConditionals[0] *  transitionMatrix[0] + branchConditionals[1] *  transitionMatrix[1] + branchConditionals[2] *  transitionMatrix[2] + branchConditionals[3] *  transitionMatrix[3]) +
                  rootConditionals[1] * theProbs[1] *
                  (branchConditionals[0] *  transitionMatrix[4] + branchConditionals[1] *  transitionMatrix[5] + branchConditionals[2] *  transitionMatrix[6] + branchConditionals[3] *  transitionMatrix[7]) +
                  rootConditionals[2] * theProbs[2] *
                  (branchConditionals[0] *  transitionMatrix[8] + branchConditionals[1] *  transitionMatrix[9] + branchConditionals[2] *  transitionMatrix[10] + branchConditionals[3] *  transitionMatrix[11]) +
                  rootConditionals[3] * theProbs[3] *
                  (branchConditionals[0] *  transitionMatrix[12] + branchConditionals[1] *  transitionMatrix[13] + branchConditionals[2] *  transitionMatrix[14] + branchConditionals[3] *  transitionMatrix[15]);
              #endif
              rootConditionals += 4UL;
              branchConditionals += 4UL;
              bookkeeping (siteID, accumulator, correction, result);
          } // siteID
        }
        break;
/****

 AMINOACIDS

 ****/
      case 20UL: {
        for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
          hyFloat accumulator = 0.;
          #ifdef _SLKP_USE_AVX_INTRINSICS
            __m256d bc_vector[5] = {_mm256_loadu_pd(branchConditionals),
              _mm256_loadu_pd(branchConditionals+4UL),
              _mm256_loadu_pd(branchConditionals+8UL),
              _mm256_loadu_pd(branchConditionals+12UL),
              _mm256_loadu_pd(branchConditionals+16UL)};


            hyFloat const * tm = transitionMatrix;

            for (unsigned long p = 0UL; p < 20UL; p++, rootConditionals++) {

              __m256d t_matrix[5] = { _mm256_loadu_pd(tm),
                                      _mm256_loadu_pd(tm+4UL),
                                      _mm256_loadu_pd(tm+8UL),
                                      _mm256_loadu_pd(tm+12UL),
                                      _mm256_loadu_pd(tm+16UL)};


              t_matrix[0] = _mm256_mul_pd(t_matrix[0], bc_vector[0]);
              t_matrix[1] = _mm256_mul_pd(t_matrix[1], bc_vector[1]);
              t_matrix[2] = _mm256_mul_pd(t_matrix[2], bc_vector[2]);
              t_matrix[3] = _mm256_mul_pd(t_matrix[3], bc_vector[3]);
              t_matrix[4] = _mm256_mul_pd(t_matrix[4], bc_vector[4]);

              t_matrix[0] = _mm256_add_pd (t_matrix[0],t_matrix[1]);
              t_matrix[1] = _mm256_add_pd (t_matrix[2],t_matrix[3]);
              t_matrix[3] = _mm256_add_pd (t_matrix[0],t_matrix[1]);

              tm += 20UL;

              accumulator += *rootConditionals * theProbs[p] * _avx_sum_4(_mm256_add_pd (t_matrix[3],t_matrix[4]));
            }
          #else // _SLKP_USE_AVX_INTRINSICS
            unsigned long rmx = 0UL;
            for (unsigned long p = 0UL; p < 20UL; p++,rootConditionals++) {
                hyFloat     r2 = 0.;

                  for (unsigned long c = 0UL; c < 20UL; c+=4UL, rmx +=4UL) {
                    r2 += (branchConditionals[c]   *  transitionMatrix[rmx] +
                          branchConditionals[c+1] *  transitionMatrix[rmx+1]) +
                          (branchConditionals[c+2] *  transitionMatrix[rmx+2] +
                          branchConditionals[c+3] *  transitionMatrix[rmx+3]);
                  }

                accumulator += *rootConditionals * theProbs[p] * r2;
            }
          #endif  // _SLKP_USE_AVX_INTRINSICS
          branchConditionals += 20UL;
          bookkeeping (siteID, accumulator, correction, result);

      } // siteID

      } // case 20
      break;
        /****

         CODONS

         ****/
      case 60UL:
      case 61UL:
      case 62UL:
      case 63UL: {
        for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
          hyFloat accumulator = 0.;

          unsigned long rmx = 0UL;
          for (unsigned long p = 0UL; p < alphabetDimension; p++,rootConditionals++) {
            hyFloat     r2 = 0.;
            unsigned       long c = 0UL;

           #ifdef _SLKP_USE_AVX_INTRINSICS

            __m256d sum256 = _mm256_setzero_pd ();

            for (; c < 60UL; c+=12UL, rmx +=12UL) {

              __m256d branches0 = _mm256_loadu_pd (branchConditionals+c),
                      branches1 = _mm256_loadu_pd (branchConditionals+c+4),
                      branches2 = _mm256_loadu_pd (branchConditionals+c+8),
                      matrix0   = _mm256_loadu_pd (transitionMatrix+rmx),
                      matrix1   = _mm256_loadu_pd (transitionMatrix+rmx+4),
                      matrix2   = _mm256_loadu_pd (transitionMatrix+rmx+8);

                branches0 = _mm256_mul_pd(branches0, matrix0);
                branches1 = _mm256_mul_pd(branches1, matrix1);
                branches2 = _mm256_mul_pd(branches2, matrix2);

                branches0 = _mm256_add_pd (branches0,branches2);
                sum256 = _mm256_add_pd (branches0,_mm256_add_pd (sum256, branches1));
             }

              r2 = _avx_sum_4(sum256);

           #else // _SLKP_USE_AVX_INTRINSICS
             for (; c < 60UL; c+=4UL, rmx +=4UL) {
                r2 += (branchConditionals[c]   *  transitionMatrix[rmx] +
                       branchConditionals[c+1] *  transitionMatrix[rmx+1]) +
                (branchConditionals[c+2] *  transitionMatrix[rmx+2] +
                 branchConditionals[c+3] *  transitionMatrix[rmx+3]);
              }
            #endif

            for (; c < alphabetDimension; c++, rmx ++) {
              r2 += branchConditionals[c]   *  transitionMatrix[rmx];
            }

            accumulator += *rootConditionals * theProbs[p] * r2;
          }

          branchConditionals += alphabetDimension;
          bookkeeping (siteID, accumulator, correction, result);

        }
      } // cases 60-63
      break;
      default: { // valid alphabetDimension >= 2

        if (alphabetDimension % 2) { // odd
          unsigned long alphabetDimension_minus1 = alphabetDimension-1;
          for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
            hyFloat accumulator = 0.;

            unsigned long rmx = 0UL;
            for (unsigned long p = 0UL; p < alphabetDimension; p++,rootConditionals++) {
              hyFloat     r2 = 0.;

              for (unsigned long c = 0UL; c < alphabetDimension_minus1; c+=2UL, rmx +=2UL) {
                r2 +=  branchConditionals[c]   *  transitionMatrix[rmx] +
                       branchConditionals[c+1] *  transitionMatrix[rmx+1];
              }

              r2 += branchConditionals[alphabetDimension_minus1]   *  transitionMatrix[rmx++];

              accumulator += *rootConditionals * theProbs[p] * r2;
            }

            branchConditionals += alphabetDimension;
            bookkeeping (siteID, accumulator, correction, result);

          }
        } else {
          for (unsigned long siteID = siteFrom; siteID < siteTo; siteID++) {
            hyFloat accumulator = 0.;

            unsigned long rmx = 0UL;
            for (unsigned long p = 0UL; p < alphabetDimension; p++,rootConditionals++) {
              hyFloat     r2 = 0.;

              for (unsigned long c = 0UL; c < alphabetDimension; c+=2UL, rmx +=2UL) {
                r2 +=  branchConditionals[c]   *  transitionMatrix[rmx] +
                       branchConditionals[c+1] *  transitionMatrix[rmx+1];
              }

              accumulator += *rootConditionals * theProbs[p] * r2;
            }

            branchConditionals += alphabetDimension;
            bookkeeping (siteID, accumulator, correction, result);

          }
        }

      } // default
    } // switch (alphabetDimension)
  } catch (long site) {
    #pragma omp critical
    {
      hy_global::ReportWarning (_String("Site ") & _String(site) & " evaluated to a 0 probability in ComputeLLWithBranchCache");
    }
    return -A_LARGE_NUMBER;
  }
  return result;
}

/*----------------------------------------------------------------------------------------------------------*/

hyFloat      _TheTree::ComputeTwoSequenceLikelihood
(
    _SimpleList   & siteOrdering,
    _DataSetFilter const* theFilter,
    long      *         lNodeFlags,
    _Vector* lNodeResolutions,
    long                siteFrom,
    long                siteTo,
    long                catID,
    hyFloat*     storageVec
)
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
    // process the leaves first

    long            alphabetDimension      =            theFilter->GetDimension(),
                    siteCount            =            theFilter->GetPatternCount(),
                    alphabetDimensionmod4  =          alphabetDimension-alphabetDimension%4;

    _CalcNode       *theNode               =            ((_CalcNode*) flatCLeaves (0));
    hyFloat  *   _hprestrict_ transitionMatrix
    =           theNode->GetCompExp(catID)->theData,
    result                 =            0.;

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }

    for (long siteID = siteFrom; siteID < siteTo; siteID++) {
        hyFloat  *tMatrix = transitionMatrix,
                     sum     = 0.;

        long siteState1 = lNodeFlags[siteOrdering.lData[siteID]],
             siteState2 = lNodeFlags[siteCount + siteOrdering.lData[siteID]];

        if (siteState1 >= 0)
            // a single character state; sweep down the appropriate column
        {
            if (siteState2 >= 0) { // both completely resolved;
                sum = tMatrix[siteState1*alphabetDimension + siteState2];
            } else { // first resolved, second is not
                hyFloat* childVector = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;
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
               hyFloat* childVector = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension;
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
                hyFloat *childVector1 = lNodeResolutions->theData + (-siteState1-1) * alphabetDimension,
                            *childVector2 = lNodeResolutions->theData + (-siteState2-1) * alphabetDimension;

                if (alphabetDimension == 4) { // special case for nuc data
                    sum = (tMatrix[0] * childVector2[0] + tMatrix[1] * childVector2[1] + tMatrix[2] * childVector2[2] + tMatrix[3] * childVector2[3])     * childVector1[0] * theProbs[0]+
                          (tMatrix[4] * childVector2[0] + tMatrix[5] * childVector2[1] + tMatrix[6] * childVector2[2] + tMatrix[7] * childVector2[3])     * childVector1[1] * theProbs[1]+
                          (tMatrix[8] * childVector2[0] + tMatrix[9] * childVector2[1] + tMatrix[10] * childVector2[2] + tMatrix[11] * childVector2[3])   * childVector1[2] * theProbs[2] +
                          (tMatrix[12] * childVector2[0] + tMatrix[13] * childVector2[1] + tMatrix[14] * childVector2[2] + tMatrix[15] * childVector2[3]) * childVector1[3] * theProbs[3];

                } else {
                    for (long r = 0; r < alphabetDimension; r++) { // 4 - unroll the loop
                        if (childVector1[r] > 0.0) {
                            hyFloat sum2 = 0.0;
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
                result += log(sum) * theFilter->theFrequencies.get (siteOrdering.lData[siteID]);
            }
        }
    }

    return result;
}

//_______________________________________________________________________________________________

void     _TheTree::SampleAncestorsBySequence (_DataSetFilter const* dsf, _SimpleList const& siteOrdering, node<long>* currentNode, _AVLListX const* nodeToIndex, hyFloat const* iNodeCache,
        _List& result, _SimpleList* parentStates, _List& expandedSiteMap, hyFloat const* catAssignments, long catCount)

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
      long              siteCount                       = dsf->GetPatternCount  (),
                        alphabetDimension              = dsf->GetDimension         (),
                        nodeIndex                       = nodeToIndex->GetXtra (nodeToIndex->Find ((BaseRef)currentNode)),
                        unitLength                     = dsf->GetUnitLength(),
                        catBlockShifter                    = catAssignments?(dsf->GetPatternCount()*GetINodeCount()):0;


        _CalcNode *     currentTreeNode = ((_CalcNode*) flatTree (nodeIndex));
        _SimpleList     sampledStates     (dsf->GetSiteCountInUnits (), 0, 0);

        hyFloat  const *  transitionMatrix = (catAssignments|| !parentStates)?nil:currentTreeNode->GetCompExp()->theData;
        hyFloat  const *  conditionals     = catAssignments?nil:(iNodeCache + nodeIndex  * siteCount * alphabetDimension);
        hyFloat        *  cache            = new hyFloat [alphabetDimension];

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

                hyFloat  randVal  = genrand_real2(),
                            totalSum = 0.;

                hyFloat  const *   matrixRow;

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

        _StringBuffer * sampledSequence = new _StringBuffer (siteCount*unitLength);
        _String  letterValue ((unsigned long) unitLength);
        for (long charIndexer = 0; charIndexer < sampledStates.countitems(); charIndexer++) {
            dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(sampledStates.lData[charIndexer]), unitLength, letterValue, &conversionAVL);
            (*sampledSequence) << letterValue;
        }
        sampledSequence->TrimSpace();
        result.AppendNewInstance(sampledSequence);
        //printf ("%d: %s\n", nodeIndex, sampledSequence->sData);

        for (long child = 1L; child <= childrenCount; child ++) {
            SampleAncestorsBySequence (dsf, siteOrdering, currentNode->go_down(child), nodeToIndex, iNodeCache, result, &sampledStates, expandedSiteMap, catAssignments, catCount);
        }
    }
}


//_______________________________________________________________________________________________

_List*   _TheTree::RecoverAncestralSequences (_DataSetFilter const* dsf,
        _SimpleList const& siteOrdering,
        _List const& expandedSiteMap,
        hyFloat * iNodeCache,
        hyFloat const* catAssignments,
        long catCount,
        long* lNodeFlags,
        _Vector * lNodeResolutions,
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
    long            patternCount                     = dsf->GetPatternCount  (),
                    alphabetDimension                = dsf->GetDimension         (),
                    unitLength                       = dsf->GetUnitLength        (),
                    iNodeCount                       = GetINodeCount             (),
                    leafCount                        = GetLeafCount              (),
                    siteCount                        = dsf->GetSiteCountInUnits    (),
                    allNodeCount                     = 0,
                    stateCacheDim                    = (alsoDoLeaves? (iNodeCount + leafCount): (iNodeCount));

    long            *stateCache                     = new long [patternCount*(iNodeCount-1)*alphabetDimension],
                    *leafBuffer                     = new long [(alsoDoLeaves?leafCount*patternCount:1)*alphabetDimension];

    // a Patterns x Int-Nodes x CharStates integer table
    // with the best character assignment for node i given that its parent state is j for a given site

    hyFloat          *buffer                         = new hyFloat [alphabetDimension];
    // iNodeCache will be OVERWRITTEN with conditional pair (i,j) conditional likelihoods


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

        hyFloat * parentConditionals = iNodeCache + parentCode * alphabetDimension * patternCount;

        if (taggedInternals.lData[parentCode] == 0L) {
            // mark the parent for update and clear its conditionals if needed
            taggedInternals.lData[parentCode]     = 1L;
            InitializeArray(parentConditionals, patternCount*alphabetDimension, 1.);
        }

        _CalcNode *          currentTreeNode = isLeaf? ((_CalcNode*) flatCLeaves (nodeCode)):((_CalcNode*) flatTree    (nodeCode));
        hyFloat  const*        transitionMatrix = nil;

        if (!catAssignments) {
            _Matrix* comp_exp = currentTreeNode->GetCompExp();
            if (!comp_exp) {
                hy_global::HandleApplicationError(_String ("Internal error in ") & __PRETTY_FUNCTION__ & ". Transition matrix not computed for " & *currentTreeNode->GetName());
                return;
            }
            transitionMatrix = comp_exp->theData;
        }

        // this will need to be toggled on a per site basis
        hyFloat  *       childVector;

        if (!isLeaf) {
            childVector = iNodeCache + (nodeCode * patternCount) * alphabetDimension;
        }

        for (long siteID = 0; siteID < patternCount; siteID++, parentConditionals += alphabetDimension) {
            if (catAssignments) {
                transitionMatrix = currentTreeNode->GetCompExp(catAssignments[siteOrdering.lData[siteID]])->theData;
            }

            hyFloat  const *tMatrix = transitionMatrix;
            if (isLeaf) {
                long siteState = lNodeFlags[nodeCode*patternCount + siteOrdering.lData[siteID]] ;
                if (siteState >= 0L) { // a fully resolved leaf
                    tMatrix  +=  siteState;
                    for (long k = 0; k < alphabetDimension; k++, tMatrix += alphabetDimension) {
                        parentConditionals[k] *= *tMatrix;
                    }
                    if (alsoDoLeaves) {
                        InitializeArray (leafBuffer, alphabetDimension, (const long)siteState);
                        leafBuffer += alphabetDimension;
                    }

                    continue;
                } else {// an ambiguous leaf
                    childVector = lNodeResolutions->theData + (-siteState-1L) * alphabetDimension;
                }

            }

            // now repopulate this vector as necessary -- if we are here this means
            // that the subtree below has been completely processed,
            // the i-th cell of childVector contains the likelihood of the _optimal_
            // assignment in the subtree below given that the character at the current
            // node is i.

            // hence, given parent state 'p', we optimize
            // max_i pr (p->i) childVector [i] and store it in the p cell of vector childVector

            hyFloat overallMax                     = 0.0;

            long       *stateBuffer                   = isLeaf?leafBuffer:stateCache;

            // check for degeneracy

            bool completely_unresolved = ArrayAll (childVector, alphabetDimension, [] (hyFloat x, unsigned long) {return x == 1.;});

            if (completely_unresolved) {
                InitializeArray(stateBuffer, alphabetDimension, -1L);
            } else {
                for (long p = 0L; p < alphabetDimension; p++) {
                    hyFloat max_lik = 0.;
                    long       max_idx = 0L;

                    for (long c = 0L; c < alphabetDimension; c++) {
                        hyFloat thisV = tMatrix[c] * childVector[c];
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
                    for (long k = 0L; k < alphabetDimension; k++) {
                        buffer[k] *= _lfScalerUpwards;
                    }
                }

                // buffer[p] now contains the maximum likelihood of the tree
                // from this point forward given that parent state is p
                // and stateBuffer[p] stores the maximizing assignment
                // for this node

                for (long k = 0; k < alphabetDimension; k++) {
                    if (stateBuffer[k] >= 0L) {
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

    hyFloat   _hprestrict_ * rootConditionals = iNodeCache + alphabetDimension * ((iNodeCount-1)  * patternCount);
    _SimpleList  parentStates (stateCacheDim,0,0),
                 conversion;

    stateCache -= patternCount*(iNodeCount-1)*alphabetDimension;
    if (alsoDoLeaves) {
        leafBuffer -= patternCount*leafCount*alphabetDimension;
    }

    _AVLListXL    conversionAVL (&conversion);
    _String       codeBuffer    ((unsigned long)unitLength);

    for (long siteID = 0; siteID < patternCount; siteID++, rootConditionals += alphabetDimension) {
        hyFloat max_lik = 0.;
        long       max_idx = 0;

        long howManyOnes = 0;
        for (long k = 0; k < alphabetDimension; k++) {
            howManyOnes += rootConditionals[k]==1.;
        }

        _SimpleList const*    patternMap = (_SimpleList const*) expandedSiteMap.GetItem(siteOrdering.lData[siteID]);

        if (howManyOnes != alphabetDimension) {
            for (long c = 0; c < alphabetDimension; c++) {
                hyFloat thisV = theProbs[c] * rootConditionals[c];
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
            dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(parentStates.lData[nodeID]), unitLength, codeBuffer, &conversionAVL);
            _String  *sequence   = (_String*) (*result)(nodeID<iNodeCount?postToIn.lData[nodeID]:nodeID);

            for (long site = 0; site < patternMap->lLength; site++) {
                unsigned long offset = patternMap->lData[site]*unitLength;
                for (long charS = 0; charS < unitLength; charS ++) {
                    sequence->set_char (offset + charS, codeBuffer.char_at(charS));
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
    _TreeIterator ti (this, _HY_TREE_TRAVERSAL_POSTORDER);
    while (_CalcNode* iterator = ti.Next()) {
        iterator->SetupCategoryMap (containerVariables,classCounter,multipliers);
    }
}


//_______________________________________________________________________________________________

hyFloat   _TheTree::Process3TaxonNumericFilter (_DataSetFilterNumeric* dsf, long catID)
{

    hyFloat *l0 =  dsf->probabilityVectors.theData +
                      dsf->categoryShifter * catID + dsf->theNodeMap.lData[0]*dsf->shifter,
                      *l1 = dsf->probabilityVectors.theData +
                            dsf->categoryShifter * catID + dsf->theNodeMap.lData[1]*dsf->shifter,
                            *l2 = dsf->probabilityVectors.theData +
                                  dsf->categoryShifter * catID + dsf->theNodeMap.lData[2]*dsf->shifter,
                                  * matrix0 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[0]->in_object)))->GetCompExp(catID)->theData,
                                    * matrix1 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[1]->in_object)))->GetCompExp(catID)->theData,
                                      * matrix2 = ((_CalcNode*)(LocateVar(theRoot->nodes.data[2]->in_object)))->GetCompExp(catID)->theData,
                                        overallResult = 0.;

    long        patternCount =  dsf->GetPatternCount();

    hyFloat  currentAccumulator = 1.;

    for (long patternIndex = 0; patternIndex < patternCount; patternIndex ++, l0+=4, l1+=4, l2+=4) {
        hyFloat rp0 = l0[0] * matrix0[0]+ l0[1]  * matrix0[1]  + l0[2] * matrix0[2]  + l0[3] * matrix0[3];
        hyFloat rp1 = l0[0] * matrix0[4]+ l0[1]  * matrix0[5]  + l0[2] * matrix0[6]  + l0[3] * matrix0[7];
        hyFloat rp2 = l0[0] * matrix0[8]+ l0[1]  * matrix0[9]  + l0[2] * matrix0[10] + l0[3] * matrix0[11];
        hyFloat rp3 = l0[0] * matrix0[12]+ l0[1] * matrix0[13] + l0[2] * matrix0[14] + l0[3] * matrix0[15];

        rp0 *= l1[0] * matrix1[0] + l1[1] * matrix1[1]  + l1[2] * matrix1[2]  + l1[3] * matrix1[3];
        rp1 *= l1[0] * matrix1[4] + l1[1] * matrix1[5]  + l1[2] * matrix1[6]  + l1[3] * matrix1[7];
        rp2 *= l1[0] * matrix1[8] + l1[1] * matrix1[9]  + l1[2] * matrix1[10] + l1[3] * matrix1[11];
        rp3 *= l1[0] * matrix1[12]+ l1[1] * matrix1[13] + l1[2] * matrix1[14] + l1[3] * matrix1[15];

        rp0 *= l2[0] * matrix2[0] + l2[1] * matrix2[1]  + l2[2] * matrix2[2]  + l2[3] * matrix2[3];
        rp1 *= l2[0] * matrix2[4] + l2[1] * matrix2[5]  + l2[2] * matrix2[6]  + l2[3] * matrix2[7];
        rp2 *= l2[0] * matrix2[8] + l2[1] * matrix2[9]  + l2[2] * matrix2[10] + l2[3] * matrix2[11];
        rp3 *= l2[0] * matrix2[12]+ l2[1] * matrix2[13] + l2[2] * matrix2[14] + l2[3] * matrix2[15];

        hyFloat  result = theProbs[0]*rp0+
                             theProbs[1]*rp1+
                             theProbs[2]*rp2+
                             theProbs[3]*rp3;


        if (result<=0.0) {
            return -A_LARGE_NUMBER;
        }

        long patternFreq = dsf->theFrequencies[patternIndex];
        for  (long freqIterator = 0; freqIterator < patternFreq; freqIterator++) {
            hyFloat tryMultiplication = currentAccumulator*result;
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
    long            leafCount = pswRepresentation.Element(-2L),
                    leafCode  = 0L,
                    L = 0L,
                    R = 0L;

    result.Clear    ();
    result.Populate (3*leafCount,-1,0);

    for (long k = 0; k < pswRepresentation.lLength-2L; k+=2) {
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
    _StringBuffer * result = new _StringBuffer (128L);
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
    result->TrimSpace();
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
            iNodeCount = -1;

    _SimpleList levelBuffer;

    node_iterator<long> ni (theRoot, _HY_TREE_TRAVERSAL_POSTORDER);

    while (node<long> * currentNode = ni.Next (&levelBuffer)) {
        _String nodeName = GetNodeName (currentNode);

        while (levelBuffer.countitems() <= ni.Level()) {
            levelBuffer << 0;
        }

        if (currentNode->is_leaf()) {
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
            pswRepresentation << levelBuffer.lData[ni.Level()];
            if (reference) {
                pswRepresentation << 0;
            } else {
                (*inames) && &nodeName;
            }

            iNodeCount--;
        }
        if (ni.Level()) {
            levelBuffer.lData[ni.Level()-1] += levelBuffer.lData[ni.Level()]+1;
        }
        levelBuffer.lData[ni.Level()]   = 0;
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

_AssociativeList *   _TreeTopology::SplitsIdentity (_PMathObj p) {
// compare tree topologies
    _Matrix     * result = new _Matrix (2,1,false,true),
                  * result2 = nil;

    _String    * tree_repr  = nil;

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

            long L, R = -1L;

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


            tree_repr  = ConvertFromPSW (nameMap, psw2);

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
    resultList->MStore ("CONSENSUS", tree_repr ? new _FString (tree_repr) : new _FString(), false);
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

 
