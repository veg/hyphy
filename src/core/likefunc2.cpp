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

#include "likefunc.h"
#include <math.h>

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

#ifdef  _SLKP_LFENGINE_REWRITE_

_String _hyMarginalSupportMatrix ("marginal_support_matrix");

/*--------------------------------------------------------------------------------------------------*/

void    _LikelihoodFunction::DetermineLocalUpdatePolicy (void)
{
    for (unsigned long k = 0; k < theTrees.lLength; k ++) {
        unsigned long catCount = ((_TheTree*)LocateVar(theTrees(k)))->categoryCount;
        
        _List * lup = new _List,
              * mte = new _List;

        computedLocalUpdatePolicy.AppendNewInstance (new _SimpleList (catCount,0,0));
        
        for (unsigned long l = 0; l < catCount; l++) {
            lup->AppendNewInstance (new _SimpleList);
            mte->AppendNewInstance (new _List);
        }

        localUpdatePolicy.AppendNewInstance      (lup);
        matricesToExponentiate.AppendNewInstance (mte);
    }
}


/*--------------------------------------------------------------------------------------------------*/

void    _LikelihoodFunction::ComputeParameterPenalty (void){
  smoothingPenalty = 0.0;
  if (smoothingTerm > 0.0) {
      //printf ("\n_LikelihoodFunction::ComputeParameterPenalty\n");
      for (unsigned long k = 0; k < indexInd.lLength; k ++) {
        _Parameter lb = GetIthIndependentBound(k, true),
                   ub = GetIthIndependentBound(k, false),
                   mp = 0.5*(lb+ub),
                   span = ub-lb,
                   v  = GetIthIndependent(k);
                   
       _Parameter term = exp (50*log (2.*fabs (v-mp)/span));
       /*if (term > 0.0) {
        printf ("\n[_LikelihoodFunction::ComputeParameterPenalty %lu: %g %g %g %g]\n", k, lb, ub, v, term); 
       }*/
       smoothingPenalty += term;
        // (2.*(v-mp)/span)^50
      }
  }
}


/*--------------------------------------------------------------------------------------------------*/

void    _LikelihoodFunction::FlushLocalUpdatePolicy (void)
{
    computedLocalUpdatePolicy.Clear();
    localUpdatePolicy.Clear();
    matricesToExponentiate.Clear();
}

//_______________________________________________________________________________________
void            _LikelihoodFunction::PartitionCatVars     (_SimpleList& storage, long partIndex)
{
    if (partIndex < blockDependancies.lLength) {
        for (long bit = 0; bit < 32; bit++)
            if (CheckNthBit(blockDependancies.lData[partIndex], bit)) {
                storage << indexCat.lData[bit];
            }
    }
}

//_______________________________________________________________________________________
long            _LikelihoodFunction::TotalRateClassesForAPartition    (long partIndex, char mode)
{
    if (partIndex >= 0 && partIndex < categoryTraversalTemplate.lLength) {
        _List* myList = (_List*)categoryTraversalTemplate(partIndex);
        if (myList->lLength)
            if (mode == 0) {
                return ((_SimpleList*)((*myList)(1)))->Element(-1);
            } else {
                long hmmCats = 1;
                _SimpleList * catVars = (_SimpleList*)(*myList)(0);
                for (long id = 0; id < catVars->lLength; id++)
                    if (mode == 1) {
                        if (((_CategoryVariable*)catVars->lData[id])->IsHiddenMarkov()) {
                            hmmCats *= ((_SimpleList*)((*myList)(1)))->Element(id);
                        }
                    } else if (mode == 2) {
                        if (((_CategoryVariable*)catVars->lData[id])->IsConstantOnPartition()) {
                            hmmCats *= ((_SimpleList*)((*myList)(1)))->Element(id);
                        }
                    }
                return hmmCats;

            }
    } else if (partIndex < 0) {
        long catCount = 1;
        if (mode == 0)
            for (long k = 0; k < indexCat.lLength; k++) {
                catCount *= ((_CategoryVariable*)LocateVar (indexCat.lData[k]))->GetNumberOfIntervals();
            }
        else if (mode == 1) {
            for (long k = 0; k < categoryTraversalTemplate.lLength; k++) {
                long partHMMCount = TotalRateClassesForAPartition(k,1);
                catCount = MAX(partHMMCount,catCount);
            }
        }
        return catCount;
    }
    return 1;
}

//_______________________________________________________________________________________
void            _LikelihoodFunction::SetupCategoryCaches      (void)
{
    categoryTraversalTemplate.Clear();
    for (long partIndex = 0; partIndex < theDataFilters.lLength; partIndex++)
        if (blockDependancies.lData[partIndex] == 0) {
            _List * noCatVarList = new _List;
            noCatVarList->AppendNewInstance (new _List);
            noCatVarList->AppendNewInstance (new _SimpleList((long)1L));
            noCatVarList->AppendNewInstance (new _SimpleList((long)1L));
            noCatVarList->AppendNewInstance (new _SimpleList());
            noCatVarList->AppendNewInstance (new _SimpleList((long)0L));
            categoryTraversalTemplate.AppendNewInstance (noCatVarList);
        } else {
            _SimpleList       myCats;
            PartitionCatVars  (myCats, partIndex);
            _List*            catVarReferences = new _List,
            *             container        = new _List;

            _SimpleList *     catVarCounts     = new _SimpleList,
            *     catVarOffsets    = new _SimpleList (myCats.lLength,1,0),
            *     hmmAndCOP        = new _SimpleList (),
            *     varType          = new _SimpleList (myCats.lLength,1,0);

            long              totalCatCount    = 1L,
                              hmmCatCount      = 1L,
                              catVarFlags      = 0L,
                              varIndex;

            for ( varIndex = 0; varIndex < myCats.lLength; varIndex++) {
                _CategoryVariable * aCV = (_CategoryVariable *)LocateVar (myCats.lData[varIndex]);
                (*catVarReferences) << aCV;
                long                intervalCount = aCV->GetNumberOfIntervals();
                (*catVarCounts)     << intervalCount;

                if (aCV->IsHiddenMarkov() || aCV->IsConstantOnPartition()) {
                    if (aCV->IsConstantOnPartition()) {
                        if (catVarFlags & (_hyphyCategoryCOP|_hyphyCategoryHMM)) {
                            break;
                        }
                        varType->lData[varIndex] = _hyphyCategoryCOP;
                    } else {
                        if (catVarFlags & (_hyphyCategoryCOP|_hyphyCategoryHMM)) {
                            break;
                        }
                        varType->lData[varIndex] = _hyphyCategoryHMM;
                    }

                    (*hmmAndCOP) << intervalCount;
                    hmmCatCount *= intervalCount;
                } else {
                    varType->lData[varIndex] = _hyphyCategoryNormal;
                }

                catVarFlags |= varType->lData[varIndex];
                totalCatCount       *= intervalCount;
            }

            if (varIndex <  myCats.lLength) {
                WarnError ("Currently, HyPhy can support at most one HMM or Constant on Partition variable per partition");
                return;
            }

            (*catVarCounts) << totalCatCount;
            (*varType)      << catVarFlags;

            for (long varIndex = myCats.lLength-2; varIndex >= 0; varIndex--) {
                catVarOffsets->lData[varIndex] = catVarOffsets->lData[varIndex+1]*catVarCounts->lData[varIndex+1];
            }

            for (long varIndex = hmmAndCOP->lLength-2; varIndex >= 0; varIndex--) {
                hmmAndCOP->lData[varIndex] *= hmmAndCOP->lData[varIndex+1];
            }

            if (hmmAndCOP->lLength) {
                (*hmmAndCOP) << hmmCatCount;
            }

            container->AppendNewInstance (catVarReferences);
            container->AppendNewInstance (catVarCounts);
            container->AppendNewInstance (catVarOffsets);
            container->AppendNewInstance (hmmAndCOP);
            container->AppendNewInstance (varType);

            ((_TheTree*)LocateVar(theTrees(partIndex)))->SetupCategoryMapsForNodes(*catVarReferences,*catVarCounts,*catVarOffsets);

            categoryTraversalTemplate.AppendNewInstance(container);
        }

    if (indexCat.lLength) {
        if (siteResults) {
            DeleteObject (siteResults);
        }
        AllocateSiteResults();
    }
}

/*--------------------------------------------------------------------------------------------------*/

void    _LikelihoodFunction::RestoreScalingFactors (long index, long branchID, long patternCnt, long* scc, long *sccb)
{
    if (branchID >= 0) { // finished using an existing cache
        overallScalingFactors[index] = overallScalingFactorsBackup[index];
        if (sccb)
            for (long recoverIndex = 0; recoverIndex < patternCnt; recoverIndex++) {
                scc[recoverIndex] = sccb[recoverIndex];
            }
    }
}

/*--------------------------------------------------------------------------------------------------*/

bool    _LikelihoodFunction::ProcessPartitionList (_SimpleList& partsToDo, _Matrix* partitionList, _String caller)
{
    long    partCount = CountObjects(0);
    partsToDo.Populate (partCount, 0, 1);
    if (partitionList) {
        partitionList->ConvertToSimpleList (partsToDo);
        DeleteObject (partitionList);
        partsToDo.Sort();
        partsToDo.FilterRange (-1, partCount);
        if (partsToDo.lLength == 0) {
            WarnError (_String("An invalid partition specification in call to ") & caller);
            return nil;
        }
    }

    return true;
}


//_______________________________________________________________________________________

void    _LikelihoodFunction::ReconstructAncestors (_DataSet &target,_SimpleList& doTheseOnes, _String& baseResultID,  bool sample, bool doMarginal, bool doLeaves)
/*
    Reconstruct ancestors for a likelihood function using

-- target      :    the _DataSet object that will receive the results
-- doTheseOnes :    a _sorted_ array of partition indices to include in this operation; is assumed to contain valid indices (i.e. 0 -- number of partitions - 1)
-- baseResultID:    the HBL identifier of the dataset that will receive the result; used as a prefix for .marginal_support_matrix support matrix (when doMarginal = true)
-- sample      :    if true, an ancestral sample (weighted by likelihood) is drawn, otherwise an ML (or maginal) reconstruction is carried out
-- doMarginal  :    if sample == false, doMarginal determines how the ancestors are reconstructed; if true, the reconstruction is marginal (maximizes
                    the likelihood of each node while summing over the rest), otherwise it is joint.
-- doLeaves    :    if sample == false and doMarginal == false (for now) and doLeaves == true, then the procedure will also
                    reconstruct (joint ML) the best assignment of leaves

*/
{
    _DataSetFilter *dsf             = (_DataSetFilter*)dataSetFilterList (theDataFilters(doTheseOnes.lData[0]));
    _TheTree        *firstTree      = (_TheTree*)LocateVar(theTrees(doTheseOnes.lData[0]));

    target.SetTranslationTable      (dsf->GetData());
    target.ConvertRepresentations();

    computationalResults.ZeroUsed();
    PrepareToCompute();

    // check if we need to deal with rate variation
    _Matrix         *rateAssignments = nil;
    if  (!doMarginal && indexCat.lLength>0) {
        rateAssignments = (_Matrix*)checkPointer(ConstructCategoryMatrix(doTheseOnes,_hyphyLFConstructCategoryMatrixClasses,false));
    } else {
        Compute();    // need to do this to populate rate matrices
    }

    long siteOffset         = 0,
         patternOffset        = 0,
         sequenceCount       ;

    for (long i = 0; i<doTheseOnes.lLength; i++) {
        long       partIndex    = doTheseOnes.lData[i];
        _TheTree   *tree        = (_TheTree*)LocateVar(theTrees(partIndex));
        dsf = (_DataSetFilter*)dataSetFilterList (theDataFilters(partIndex));

        long    catCounter = 0;

        if (rateAssignments) {
            _SimpleList             pcats;
            PartitionCatVars        (pcats,partIndex);
            catCounter            = pcats.lLength;
        }

        if (i==0) {
            tree->AddNodeNamesToDS (&target,sample == false && doLeaves,!(doLeaves && doMarginal) ,2*(doMarginal == false && sample == false && doLeaves));
            // store internal or leaf node names in the dataset
            sequenceCount = target.GetNames().lLength;
        } else {
            if (!tree->Equal(firstTree)) { // incompatible likelihood function
                ReportWarning ((_String("Ancestor reconstruction had to ignore partition ")&_String(partIndex+1)&" of the likelihood function since it has a different tree topology than the first part."));
                continue;
            }
            _TranslationTable * mtt = target.GetTT()->MergeTables(dsf->GetData()->GetTT());
            if (mtt) {
                target.SetTranslationTable      (mtt);
                DeleteObject                    (mtt);
            } else {
                ReportWarning ((_String("Ancestor reconstruction had to ignore partition ")&_String(partIndex+1)&" of the likelihood function since it has a character alphabet incompatible with the first part."));
                continue;
            }
        }

        _List       * expandedMap   = dsf->ComputePatternToSiteMap(),
                      * thisSet;

        if (sample) {
            _AVLListX   * nodeMapper    = tree->ConstructNodeToIndexMap(true);
            thisSet                     = new _List;
            _SimpleList* tcc            = (_SimpleList*)treeTraversalMasks(partIndex);
            if (tcc) {
                long shifter = dsf->GetDimension()*dsf->NumberDistinctSites()*tree->GetINodeCount();
                for (long cc = 0; cc <= catCounter; cc++) {
                    tree->FillInConditionals(dsf, conditionalInternalNodeLikelihoodCaches[partIndex] + cc*shifter, tcc);
                }
            }
            tree->SampleAncestorsBySequence (dsf, *(_SimpleList*)optimalOrders.lData[partIndex],
                                             &tree->GetRoot(),
                                             nodeMapper,
                                             conditionalInternalNodeLikelihoodCaches[partIndex],
                                             *thisSet,
                                             nil,
                                             *expandedMap,
                                             catCounter?rateAssignments->theData+siteOffset:nil,
                                             catCounter);


            nodeMapper->DeleteAll(false);
            DeleteObject (nodeMapper);

        } else {
            if (doMarginal) {
                _Matrix  *marginals = new _Matrix;
                _String  supportMxID = baseResultID & '.' & _hyMarginalSupportMatrix;
                thisSet = RecoverAncestralSequencesMarginal (partIndex, *marginals, *expandedMap, doLeaves);
                CheckReceptacleAndStore(&supportMxID, "ReconstructAncestors", true, marginals, false);

            } else
                thisSet = tree->RecoverAncestralSequences (dsf,
                          *(_SimpleList*)optimalOrders.lData[partIndex],
                          *expandedMap,
                          conditionalInternalNodeLikelihoodCaches[partIndex],
                          catCounter?rateAssignments->theData+siteOffset:nil,
                          catCounter,
                          conditionalTerminalNodeStateFlag[partIndex],
                          (_GrowingVector*)conditionalTerminalNodeLikelihoodCaches(partIndex),
                          doLeaves
                                                          );

        }


        _String * sampledString = (_String*)(*thisSet)(0);

        for (long siteIdx = 0; siteIdx<sampledString->sLength; siteIdx++) {
            target.AddSite (sampledString->sData[siteIdx]);
        }

        for (long seqIdx = 1; seqIdx < sequenceCount; seqIdx++) {
            sampledString = (_String*)(*thisSet)(seqIdx);
            for (long siteIdx = 0; siteIdx<sampledString->sLength; siteIdx++) {
                target.Write2Site (siteOffset + siteIdx, sampledString->sData[siteIdx]);
            }
        }
        DeleteObject (thisSet);
        DeleteObject (expandedMap);
        siteOffset    += dsf->GetSiteCount();
        patternOffset += dsf->GetSiteCount();
    }


    target.Finalize();
    target.SetNoSpecies(target.GetNames().lLength);

    if (rateAssignments) {
        DeleteObject (rateAssignments);
    }

    DoneComputing ();

}

//_______________________________________________________________________________________________

void            _LikelihoodFunction::PopulateConditionalProbabilities   (long index, char runMode, _Parameter* buffer, _SimpleList& scalers, long branchIndex, _SimpleList* branchValues)
// this function computes site probabilties for each rate class (or something else that involves iterating over rate classes)
// see run options below

// run mode can be one of the following

// _hyphyLFConditionProbsRawMatrixMode : simply   populate an M (number of rate classes) x S (number of site patterns) matrix of conditional likelihoods
//   : expected minimum dimension of buffer is M*S
//   : scalers will have M*S entries laid out as S for rate class 0, S for rate class 1, .... S for rate class M-1

// _hyphyLFConditionProbsScaledMatrixMode : simply   populate an M (number of rate classes) x S (number of site patterns) and scale to the lowest multiplier
//   : expected minimum dimension of buffer is M*S
//   : scalers will have S entries

// _hyphyLFConditionProbsWeightedSum : compute  a sum for each site using weighted by the probability of a given category
//   : expected minimum dimension of buffer is 2*S
//   : scalers will have S entries
//   : **note that the behavior is different if there are HMM (or constant on partition) variables
//   : the size of the buffer is S*(N+1), where N is the cumulative number of categories in such variables for this partition
//   : the size of the scaler is also S*N
//   : the code will behave as _hyphyLFConditionProbsScaledMatrixMode with all other category variables
//   : summed CONDITIONED on the values of HMM/Constant on partition

// _hyphyLFConditionProbsMaxProbClass : compute the category index of maximum probability
//   : expected minimum dimension of buffer is 3*S -- the result goes into offset 0
//   : scalers will have S entries

// _hyphyLFConditionProbsClassWeights : compute the weight of each rate class index
//   : expected minimum dimension of buffer is M
//   : scalers will have no entries

// _hyphyLFConditionMPIIterate : compute conditional likelihoods of the partition using MPI
//   : run mode effectively the same as _hyphyLFConditionProbsWeightedSum
{
    _List               *traversalPattern       = (_List*)categoryTraversalTemplate(index),
                         *variables                = (_List*)((*traversalPattern)(0)),
                          *catWeigths               = nil;

    _SimpleList         *categoryCounts         = (_SimpleList*)((*traversalPattern)(1)),
                         *categoryOffsets     = (_SimpleList*)((*traversalPattern)(2)),
                          *hmmAndCOP                = (_SimpleList*)((*traversalPattern)(3)),
                           categoryValues         (categoryCounts->lLength,0,0);

    long                totalSteps              = categoryOffsets->lData[0] * categoryCounts->lData[0],
                        catCount              = variables->lLength-1,
                        blockLength               = BlockLength(index),
                        hmmCatSize               = hmmAndCOP->Element(-1),
                        hmmCatCount                = hmmAndCOP->lLength?(totalSteps/hmmCatSize):0,
                        currentHMMCat         = 1,
                        arrayDim              ;

    bool                isTrivial               = variables->lLength == 0,
                        switchingHMM         = false;

    _CategoryVariable   *catVariable;


    switch (runMode) {
    case _hyphyLFConditionProbsRawMatrixMode:
        arrayDim = catCount*blockLength;
        break;
    case _hyphyLFConditionProbsClassWeights:
        arrayDim = 0;
        break;
    default:
        arrayDim = hmmCatCount?blockLength*hmmCatSize:blockLength;
    }

    if (runMode == _hyphyLFConditionProbsWeightedSum || runMode == _hyphyLFConditionMPIIterate || runMode == _hyphyLFConditionProbsClassWeights) {
        if (runMode == _hyphyLFConditionProbsWeightedSum || runMode == _hyphyLFConditionMPIIterate) {
            long upperBound = hmmCatCount?hmmAndCOP->Element(-1)*blockLength:blockLength;
            for (long r = 0; r < upperBound; r++) {
                buffer[r] = 0.;
            }
        }
        catWeigths = new _List;
    } else if (runMode == _hyphyLFConditionProbsMaxProbClass)
        for (long r = 0, r2 = 2*blockLength; r < blockLength; r++, r2++) {
            buffer[r] = 0.0;
            buffer[r2] = 0.0;
        }

    for                 (long currentCat        = 0; currentCat <= catCount; currentCat++) {
        (catVariable = ((_CategoryVariable**)(variables->lData))[currentCat])->Refresh();
        catVariable->SetIntervalValue(0,true);
        if (runMode == _hyphyLFConditionProbsWeightedSum || runMode == _hyphyLFConditionMPIIterate || runMode == _hyphyLFConditionProbsClassWeights) {
            (*catWeigths) << catVariable->GetWeights();
        }
    }


    scalers.Populate    (arrayDim,0,0);

#ifdef __HYPHYMPI__
    _GrowingVector * computedWeights = nil;
    if (runMode == _hyphyLFConditionMPIIterate) {
        computedWeights = new _GrowingVector;
    }
    long                mpiTasksSent = 0;
#endif

    for                 (long pass = 0; pass < 1+(runMode == _hyphyLFConditionMPIIterate); pass++) {
        for                 (long currentRateCombo  = 0; currentRateCombo < totalSteps; currentRateCombo++) {

            // setting each category variable to its appropriate value

            _Parameter       currentRateWeight = 1.;
            if (pass == 0) {
                if (!isTrivial) {
                    long remainder = currentRateCombo % categoryCounts->lData[catCount];

                    if (hmmCatCount) {
                        currentHMMCat = currentRateCombo / hmmCatCount;
                        switchingHMM = (currentRateCombo % hmmCatCount) == 0;
                    }

                    if (currentRateCombo && remainder  == 0) {
                        categoryValues.lData[catCount] = 0;
                        (((_CategoryVariable**)(variables->lData))[catCount])->SetIntervalValue(0);
                        for (long uptick = catCount-1; uptick >= 0; uptick --) {
                            categoryValues.lData[uptick]++;
                            if (categoryValues.lData[uptick] == categoryCounts->lData[uptick]) {
                                categoryValues.lData[uptick] = 0;
                                (((_CategoryVariable**)(variables->lData))[uptick])->SetIntervalValue(0);
                            } else {
                                (((_CategoryVariable**)(variables->lData))[uptick])->SetIntervalValue(categoryValues.lData[uptick]);
                                break;
                            }
                        }
                    } else {
                        if (currentRateCombo) {
                            categoryValues.lData[catCount]++;
                            (((_CategoryVariable**)(variables->lData))[catCount])->SetIntervalValue(remainder);
                        }
                    }
                }

                if (runMode == _hyphyLFConditionProbsWeightedSum || runMode == _hyphyLFConditionProbsClassWeights || runMode == _hyphyLFConditionMPIIterate) {
                    for                 (long currentCat        = hmmCatCount; currentCat <= catCount; currentCat++) {
                        currentRateWeight *= ((_Matrix**)catWeigths->lData)[currentCat]->theData[categoryValues.lData[currentCat]];
                    }

#ifdef __HYPHYMPI__
                    if (runMode == _hyphyLFConditionMPIIterate && pass == 0) {
                        computedWeights->Store(currentRateWeight);
                    }
#endif
                    if (runMode == _hyphyLFConditionProbsClassWeights) {
                        buffer [currentRateCombo] = currentRateWeight;
                        continue;
                    } else if (currentRateWeight == 0.0) { // nothing to do, eh?
                        continue;
                    }
#ifdef __HYPHYMPI__
                    else {
                        if (runMode == _hyphyLFConditionMPIIterate) {
                            SendOffToMPI (currentRateCombo);
                            mpiTasksSent ++;
                            continue;
                        }
                    }
#endif
                }
            }

            long useThisPartitonIndex = currentRateCombo;

#ifdef __HYPHYMPI__
            if (runMode == _hyphyLFConditionMPIIterate) {
                MPI_Status     status;
                ReportMPIError(MPI_Recv (resTransferMatrix.theData, resTransferMatrix.GetSize(), MPI_DOUBLE, MPI_ANY_SOURCE , HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD,&status),true);
                useThisPartitonIndex = status.MPI_SOURCE-1;
                currentRateWeight    = computedWeights->theData[useThisPartitonIndex];
            }
#endif

            // now that the categories are set we can proceed with the computing step
            long             indexShifter                   = blockLength * useThisPartitonIndex;
            long             *siteCorrectors                = ((_SimpleList**)siteCorrections.lData)[index]->lLength?
                    (((_SimpleList**)siteCorrections.lData)[index]->lData) + indexShifter
                    :nil;


            if (runMode == _hyphyLFConditionProbsRawMatrixMode || runMode == _hyphyLFConditionProbsScaledMatrixMode)
                // populate the matrix of conditionals and scaling factors
            {
                _Parameter  _hprestrict_ *bufferForThisCategory = buffer + indexShifter;

                ComputeBlock    (index, bufferForThisCategory, useThisPartitonIndex, branchIndex, branchValues);
                if (usedCachedResults) {
                    bool saveFR = forceRecomputation;
                    forceRecomputation = true;
                    ComputeBlock    (index, bufferForThisCategory, useThisPartitonIndex, branchIndex, branchValues);
                    forceRecomputation = saveFR;
                }

                if (runMode == _hyphyLFConditionProbsRawMatrixMode)
                    for (long p = 0; p < blockLength; p++) {
                        scalers.lData[p+indexShifter] = siteCorrectors[p];
                    }
                else {
                    if (siteCorrectors) {
                        for (long r1 = 0; r1 < blockLength; r1++) {
                            long scv              = *siteCorrectors,
                                 scalerDifference = scv-scalers.lData[r1];

                            if (scalerDifference > 0)
                                // this class has a _bigger_ scaling factor than at least one other class
                                // hence it needs to be scaled down (unless it's the first class)
                            {
                                if (useThisPartitonIndex==0) { //(scalers.lData[r1] == -1)
                                    scalers.lData[r1] = scv;
                                } else {
                                    bufferForThisCategory[r1] *= acquireScalerMultiplier (scalerDifference);
                                }
                            } else {
                                if (scalerDifference < 0)
                                    // this class is a smaller scaling factor, i.e. its the biggest among all those
                                    // considered so far; all other classes need to be scaled down
                                {
                                    _Parameter scaled = acquireScalerMultiplier (-scalerDifference);
                                    for (long z = indexShifter+r1-blockLength; z >= 0; z-=blockLength) {
                                        buffer[z] *= scaled;
                                    }

                                    scalers.lData[r1] = scv;
                                }
                            }
                            siteCorrectors++;
                        }
                    }
                }
            } else {
                if (runMode == _hyphyLFConditionProbsWeightedSum || runMode == _hyphyLFConditionProbsMaxProbClass || runMode == _hyphyLFConditionMPIIterate) {
                    //if (branchIndex>=0)
                    //  ((_TheTree*)LocateVar(theTrees.lData[index]))->AddBranchToForcedRecomputeList (branchIndex+((_TheTree*)LocateVar(theTrees.lData[index]))->GetLeafCount());

#ifdef          __HYPHYMPI__
                    if (runMode == _hyphyLFConditionMPIIterate) {
                        long offset = resTransferMatrix.GetVDim();

                        for (long k = 0; k < blockLength; k++) {
                            buffer[blockLength+k] = resTransferMatrix.theData[k];
                            siteCorrectors[k]     = resTransferMatrix.theData[k+offset];
                        }
                    } else

#endif

                        ComputeBlock    (index, buffer + (hmmCatCount?hmmCatSize:1)*blockLength, useThisPartitonIndex, branchIndex, branchValues);

                    if (runMode != _hyphyLFConditionMPIIterate && usedCachedResults) {
                        bool saveFR = forceRecomputation;
                        forceRecomputation = true;
                        ComputeBlock    (index, buffer + (hmmCatCount?hmmCatSize:1)*blockLength, useThisPartitonIndex, branchIndex, branchValues);
                        forceRecomputation = saveFR;
                    }


                    if (runMode == _hyphyLFConditionProbsWeightedSum || runMode == _hyphyLFConditionMPIIterate) {
                        long lowerBound  = hmmCatCount?blockLength*currentHMMCat:0,
                             upperBound  = hmmCatCount?blockLength*(1+currentHMMCat):blockLength,
                             lowerBound2 = hmmCatCount?(hmmCatSize*blockLength):blockLength;


                        for (long r1 = lowerBound, r2 = lowerBound2; r1 < upperBound; r1++,r2++) {
                            if (siteCorrectors) {
                                long scv = *siteCorrectors;

                                if (scv < scalers.lData[r1]) { // this class has a _smaller_ scaling factor
                                    buffer[r1] = currentRateWeight * buffer[r2] + buffer[r1] * acquireScalerMultiplier (scalers.lData[r1] - scv);
                                    scalers.lData[r1] = scv;
                                } else {
                                    if (scv > scalers.lData[r1]) { // this is a _larger_ scaling factor
                                        buffer[r1] += currentRateWeight * buffer[r2] * acquireScalerMultiplier (scv - scalers.lData[r1]);
                                    } else { // same scaling factors
                                        buffer[r1] += currentRateWeight * buffer[r2];
                                    }
                                }

                                siteCorrectors++;
                            } else {
                                buffer[r1] += currentRateWeight * buffer[r2];
                            }

                        }
                    } else { // runMode = _hyphyLFConditionProbsMaxProbClass
                        for (long r1 = blockLength*2, r2 = blockLength, r3 = 0; r3 < blockLength; r1++,r2++,r3++) {
                            bool doChange = false;
                            if (siteCorrectors) {
                                long scv  = *siteCorrectors,
                                     diff = scv - scalers.lData[r3];

                                if (diff<0) { // this has a _smaller_ scaling factor
                                    _Parameter scaled = buffer[r1]*acquireScalerMultiplier (diff);
                                    if (buffer[r2] > scaled) {
                                        doChange = true;
                                    } else {
                                        buffer[r1] = scaled;
                                    }
                                    scalers.lData[r3] = scv;
                                } else {
                                    if (diff>0) { // this is a _larger_ scaling factor
                                        buffer[r2] *= acquireScalerMultiplier (-diff);
                                    }
                                    doChange = buffer[r2] > buffer[r1] && ! CheckEqual (buffer[r2],buffer[r1]);
                                }

                                siteCorrectors++;
                            } else {
                                doChange = buffer[r2] > buffer[r1] && ! CheckEqual (buffer[r2],buffer[r1]);
                            }

                            if (doChange) {
                                buffer[r1]         = buffer[r2];
                                buffer[r3]         = useThisPartitonIndex;
                            }
                        }
                    }
                }
            }
#ifdef __HYPHYMPI__
            if (--mpiTasksSent == 0) {
                break;
            }
#endif
        }
    }
#ifdef __HYPHYMPI__
    DeleteObject (computedWeights);
#endif
    DeleteObject (catWeigths);
}

//_______________________________________________________________________________________________

void            _LikelihoodFunction::ComputeSiteLikelihoodsForABlock    (long index, _Parameter* results, _SimpleList& scalers, long branchIndex, _SimpleList* branchValues, char mpiRunMode)
// assumes that results is at least blockLength slots long
{
    if (blockDependancies.lData[index]) {
        PopulateConditionalProbabilities(index, mpiRunMode == _hyphyLFMPIModeREL ?_hyphyLFConditionMPIIterate:_hyphyLFConditionProbsWeightedSum, results, scalers, branchIndex, branchValues);
    } else {
        ComputeBlock        (index, results, -1, branchIndex, branchValues);
        scalers.Clear       ();
        scalers.Duplicate   (siteCorrections(index));
    }
}

//_______________________________________________________________________________________________

_List*   _LikelihoodFunction::RecoverAncestralSequencesMarginal (long index, _Matrix & supportValues, _List& expandedSiteMap, bool doLeaves)
// index:           which part to process
// supportValues:   for each internal node and site stores alphabetDimension values for the
//              :   relative support of each residue at a given site
//              :   linearized 3D matrix
//              :   1st - node index (same order as flatTree)
//              :   2nd - site index (only unique patterns are stored)
//              :   3rd - the character

// doLeaves     :   compute support values leaves instead of internal nodes

{

    _DataSetFilter* dsf             = (_DataSetFilter*)dataSetFilterList (theDataFilters(index));
    _TheTree        *blockTree      = (_TheTree*)LocateVar(theTrees.lData[index]);

    long            patternCount                    = dsf->NumberDistinctSites  (),
                    alphabetDimension                = dsf->GetDimension         (),
                    unitLength                        = dsf->GetUnitLength        (),
                    iNodeCount                        = blockTree->GetINodeCount  (),
                    leafCount                     = blockTree->GetLeafCount   (),
                    matrixSize                       = doLeaves?leafCount:iNodeCount,
                    siteCount                        = dsf->GetSiteCount         (),
                    shiftForTheNode                 = patternCount * alphabetDimension;

    _Parameter      *siteLikelihoods                = new _Parameter [2*patternCount],
    *siteLikelihoodsSpecState       = new _Parameter [2*patternCount];

    _SimpleList     scalersBaseline,
                    scalersSpecState,
                    branchValues,
                    postToIn;

    blockTree->MapPostOrderToInOderTraversal (postToIn, doLeaves == false);
    supportValues.Clear                      ();
    CreateMatrix                             (&supportValues,matrixSize,shiftForTheNode,false,true,false);

    ComputeSiteLikelihoodsForABlock          (index, siteLikelihoods, scalersBaseline);
    // establish a baseline likelihood for each site

    if (doLeaves) {
        for                             (long currentChar = 0; currentChar < alphabetDimension; currentChar++) {
            branchValues.Populate           (patternCount,currentChar,0);
            for (long branchID = 0; branchID < leafCount; branchID ++) {
                blockTree->AddBranchToForcedRecomputeList (branchID);
                long mappedBranchID = postToIn.lData[branchID];
                ComputeSiteLikelihoodsForABlock (index, siteLikelihoodsSpecState, scalersSpecState,
                                                 branchID+iNodeCount, &branchValues);
                for (long siteID = 0; siteID < patternCount; siteID++) {
                    long scaleDiff = (scalersSpecState.lData[siteID]-scalersBaseline.lData[siteID]);
                    _Parameter ratio = siteLikelihoodsSpecState[siteID]/siteLikelihoods[siteID];

                    if (scaleDiff > 0) {
                        ratio *= acquireScalerMultiplier(scaleDiff);
                    }
                    supportValues.theData[mappedBranchID*shiftForTheNode + siteID*alphabetDimension + currentChar] = ratio;
                }
                blockTree->AddBranchToForcedRecomputeList (branchID);
            }
        }
    }

    else
        for                             (long currentChar = 0; currentChar < alphabetDimension-1; currentChar++)
            // the prob for the last char is  (1 - sum (probs other chars))
        {
            branchValues.Populate           (patternCount,currentChar,0);
            for (long branchID = 0; branchID < iNodeCount; branchID ++) {
                long mappedBranchID = postToIn.lData[branchID];
                ComputeSiteLikelihoodsForABlock (index, siteLikelihoodsSpecState, scalersSpecState, branchID, &branchValues);
                for (long siteID = 0; siteID < patternCount; siteID++) {
                    long scaleDiff = (scalersSpecState.lData[siteID]-scalersBaseline.lData[siteID]);
                    _Parameter ratio = siteLikelihoodsSpecState[siteID]/siteLikelihoods[siteID];
                    if (scaleDiff > 0) {
                        ratio *= acquireScalerMultiplier(scaleDiff);
                    }
                    supportValues.theData[mappedBranchID*shiftForTheNode + siteID*alphabetDimension + currentChar] = ratio;
                }
                blockTree->AddBranchToForcedRecomputeList (branchID+leafCount);
            }
        }

    _SimpleList  conversion;
    _AVLListXL   conversionAVL (&conversion);
    _String      codeBuffer    (unitLength, false);
    _List        *result       = new _List;

    for (long k = 0; k < matrixSize; k++) {
        result->AppendNewInstance (new _String(siteCount*unitLength,false));
    }

    for (long siteID = 0; siteID < patternCount; siteID++) {
        _SimpleList*    patternMap = (_SimpleList*) expandedSiteMap (siteID);

        for  (long nodeID = 0; nodeID < matrixSize ; nodeID++) {
            long            mappedNodeID = postToIn.lData[nodeID];
            _Parameter      max_lik     = 0.,
                            sum         = 0.,
                            *scores       = supportValues.theData + shiftForTheNode*mappedNodeID +  siteID*alphabetDimension;
            long            max_idx     = 0;

            for (long charID = 0; charID < alphabetDimension-(!doLeaves); charID ++) {
                sum+=scores[charID];
                if (scores[charID] > max_lik) {
                    max_idx = charID;
                    max_lik = scores[charID];

                }
            }

            //if (fabs(scores[alphabetDimension-1]+sum-1.) > 0.1)
            //  WarnError (_String("Bad monkey!") & scores[alphabetDimension-1] & ":" & (1.-sum) );

            if (doLeaves) {
                sum = 1./sum;
                for (long charID = 0; charID < alphabetDimension; charID ++) {
                    scores [charID] *= sum;
                    /*if (siteID == 16)
                        printf ("Site %ld Leaf %ld (%ld) Char %ld = %g\n", siteID, nodeID, mappedNodeID, charID,
                                supportValues.theData[mappedNodeID*shiftForTheNode + siteID*alphabetDimension + charID]);
                     */

                }
            } else {
                scores[alphabetDimension-1] = 1. - sum;

                if (scores[alphabetDimension-1] > max_lik) {
                    max_idx = alphabetDimension-1;
                }
            }

            dsf->ConvertCodeToLettersBuffered (dsf->CorrectCode(max_idx), unitLength, codeBuffer.sData, &conversionAVL);
            _String  *sequence   = (_String*) (*result)(mappedNodeID);

            for (long site = 0; site < patternMap->lLength; site++) {
                //if (patternMap->lData[site] == 119)
                //  printf ("%ld\n",
                //          siteID);
                char* storeHere = sequence->sData + patternMap->lData[site]*unitLength;
                for (long charS = 0; charS < unitLength; charS ++) {
                    storeHere[charS] = codeBuffer.sData[charS];
                }
            }

        }
    }
    delete [] siteLikelihoods;
    delete [] siteLikelihoodsSpecState;
    return result;
}

//__________________________________________________________________________________

_Parameter          _LikelihoodFunction::SumUpHiddenMarkov (const _Parameter * patternLikelihoods, _Matrix& hmm, _Matrix& hmf, _SimpleList * duplicateMap, const _SimpleList* scalers, long bl)
{
    long               ni           = hmm.GetHDim(),
                       mi           = duplicateMap?duplicateMap->lData[duplicateMap->lLength-1]:bl-1,
                       siteScaler    = duplicateMap?scalers->lData[mi]:((_SimpleList*)((_List*)scalers)->lData[0])->lData[mi];

    _Matrix            temp  (ni,1,false,true),
                       temp2 (ni,1,false,true);

    _Parameter         correctionFactor = 0; // correction factor

    for (long m=0, mi2 = mi; m<ni; m++,mi2 += bl) {
        long currentScaler = duplicateMap?scalers->lData[mi2]:((_SimpleList*)((_List*)scalers)->lData[m])->lData[mi];
        if (currentScaler < siteScaler) { // this class has a _smaller_ scaling factor
            _Parameter upby = acquireScalerMultiplier (siteScaler - currentScaler);
            for (long rescale = 0; rescale < m; rescale ++) {
                temp2.theData[rescale] *= upby;
            }

            temp2.theData[m] = patternLikelihoods[mi2];
            siteScaler = currentScaler;
        } else {
            if (currentScaler > siteScaler) { // this is a _larger_ scaling factor
                temp2.theData[m] = patternLikelihoods[mi2]*acquireScalerMultiplier (currentScaler - siteScaler);
            } else { // same scaling factors
                temp2.theData[m] = patternLikelihoods[mi2];
            }
        }
    }

    for (long i=duplicateMap?duplicateMap->lLength-2:bl-2; i>=0; i--) {
        _Parameter max        = 0.;
        long       siteScaler = duplicateMap?
                                scalers->lData[duplicateMap->lData[i]]:
                                ((_SimpleList*)((_List*)scalers)->lData[0])->lData[i];

        for (long k=0; k<ni; k++) {
            _Parameter scrap = 0.;

            mi    = duplicateMap?duplicateMap->lData[i]:i;

            for (long m=0; m<ni; m++,mi += bl) {
                long currentScaler = duplicateMap?
                                     scalers->lData[mi]:
                                     ((_SimpleList*)((_List*)scalers)->lData[m])->lData[i];

                if (currentScaler < siteScaler) { // this class has a _smaller_ scaling factor
                    _Parameter upby = acquireScalerMultiplier (siteScaler - currentScaler);
                    for (long rescale = 0; rescale < k; rescale ++) {
                        temp.theData[rescale] *= upby;
                    }

                    scrap = scrap * upby + hmm.theData[k*ni+m] * patternLikelihoods[mi] * temp2.theData[m];
                    siteScaler = currentScaler;
                } else {
                    if (currentScaler > siteScaler) // this is a _larger_ scaling factor
                        scrap += hmm.theData[k*ni+m] * patternLikelihoods[mi] * temp2.theData[m]
                                 *acquireScalerMultiplier (currentScaler - siteScaler);

                    else { // same scaling factors
                        scrap += hmm.theData[k*ni+m] * patternLikelihoods[mi] * temp2.theData[m];
                    }
                }
            }

            temp.theData[k] = scrap;

            if (scrap>max) {
                max = scrap;
            }
        }

        if (max <= 0.0) {
            return -A_LARGE_NUMBER;
        }

        correctionFactor -= log (max);
        if (siteScaler) {
            correctionFactor -= siteScaler * _logLFScaler;
        }

        max = 1./max;
        for (long k=0; k<ni; k++) {
            temp.theData[k] *= max;
        }

        _Parameter* swap = temp.theData;
        temp.theData     = temp2.theData;
        temp2.theData    = swap;
    }

    _Parameter scrap = 0.0;

    for (long k=0; k<ni; k++) {
        scrap += temp2.theData[k] * hmf.theData[k];
    }

    return myLog(scrap) - correctionFactor;
}

//__________________________________________________________________________________

void        _LikelihoodFunction::RunViterbi ( _Matrix & result,                 const _Parameter * patternLikelihoods,
        _Matrix & hmm,                  _Matrix& hmf,
        _SimpleList * duplicateMap,       const _SimpleList* scalers,
        long bl )
{
    long               ni           = hmm.GetHDim(),
                       siteCount    = duplicateMap?duplicateMap->lLength:bl;

    _Matrix            temp  (ni,1,false,true),
                       temp2 (ni,1,false,true);

    _SimpleList        pathRecovery (siteCount * ni, 0, 0);

    if  ((duplicateMap?duplicateMap->lLength:bl) > 1)
        // non-trivial case (more than one site)
    {
        for (long site = siteCount-1; site > 0; site --) {
            for (long parentState = 0; parentState < ni; parentState ++) {
                long            bestState     = 0,
                                mi            = duplicateMap?
                                                duplicateMap->lData[site]:
                                                site,

                                                currentScaler = duplicateMap?
                                                        scalers->lData[mi]:
                                                        ((_SimpleList*)((_List*)scalers)->lData[0])->lData[site];


                _Parameter      bestValue = log(patternLikelihoods[mi]*hmm.theData[parentState*ni]) + temp.theData[0];
                mi           += bl;

                if (currentScaler) {
                    bestValue -= currentScaler * _logLFScaler;
                }

                for (long currentState = 1; currentState < ni; currentState ++, mi += bl) {
                    currentScaler = duplicateMap?
                                    scalers->lData[mi]:
                                    ((_SimpleList*)((_List*)scalers)->lData[currentState])->lData[site];

                    _Parameter      currentValue = log(patternLikelihoods[mi]*hmm.theData[parentState*ni + currentState]) +
                                                   temp.theData[currentState];
                    if (currentScaler) {
                        currentValue -= currentScaler * _logLFScaler;
                    }

                    if (currentValue > bestValue) {
                        bestValue = currentValue;
                        bestState = currentState;
                    }
                }
                temp2.theData[parentState] = bestValue;
                pathRecovery.lData[site*ni + parentState] = bestState;
                //if (parentState != bestState && parentState == 1)
                //  printf ("%d %d -> %d\n", site, parentState, bestState);
            }
            _Parameter* swap = temp.theData;
            temp.theData     = temp2.theData;
            temp2.theData    = swap;
        }
    } else {
        for (long parentState = 0; parentState < ni; parentState ++) {
            temp.theData[parentState] = log(patternLikelihoods[parentState]) +
                                        (duplicateMap?scalers->lData[parentState]:((_SimpleList*)((_List*)scalers)->lData[parentState])->lData[0])*_logLFScaler;
        }
    }

    long            mi        = duplicateMap?duplicateMap->lData[0]:0,
                    bestState = 0;

    _Parameter      bestValue     = log(patternLikelihoods [mi]*hmf.theData[0]) + temp.theData[0] +
                                    (duplicateMap?scalers->lData[mi]:((_SimpleList*)((_List*)scalers)->lData[0])->lData[0])*_logLFScaler;

    mi+=bl;


    for (long initState = 1; initState < ni; initState ++, mi += bl) {
        long currentScaler = duplicateMap?scalers->lData[mi]:((_SimpleList*)((_List*)scalers)->lData[initState])->lData[0];

        _Parameter      currentValue = log(patternLikelihoods[mi]*hmf.theData[initState]) +
                                       temp.theData[initState];
        if (currentScaler) {
            currentValue -= currentScaler * _logLFScaler;
        }

        if (currentValue > bestValue) {
            bestValue = currentValue;
            bestState = initState;
        }
    }

    result.theData[0] = bestState;

    for (long site = 1; site < siteCount; site++) {
        result.theData[site] = pathRecovery.lData[site*ni + (long)result.theData[site-1]];
        //printf ("%d: (%g) 0-%ld 1-%ld\n",  site,result.theData[site-1], pathRecovery.lData[site*ni], pathRecovery.lData[site*ni+1]);
    }

}

//_______________________________________________________________________________________________


_Parameter mapParameterToInverval (_Parameter in, char type, bool inverse)
{
    switch (type) {
    case _hyphyIntervalMapExpit:
        if (inverse) {
            //return log (in / (1-in));
            return tan (M_PI * (in - 0.5));
        } else {
            //return 1. / (1. + exp(-in));
            return atan (in) * M_1_PI + 0.5;
        }
        break;
    case _hyphyIntervalMapSqueeze:
        if (inverse) {
            return in/(1.-in);
        } else {
            return in/(1.+in);
        }
        break;

    }
    return in;
}

//_______________________________________________________________________________________________

void _LikelihoodFunction::SetupParameterMapping (void)
{
    parameterTransformationFunction.Clear();
    parameterValuesAndRanges = new _Matrix (indexInd.lLength, 4, false, true);
    checkParameter(addLFSmoothing, smoothingTerm, 0.0);
    checkParameter(reduceLFSmoothing, smoothingReduction, 0.8);
    if (smoothingPenalty < 0.0) {
      smoothingPenalty = 0.0;
    }
    if (smoothingReduction <= 0.0 || smoothingReduction >= 1.0) {
      smoothingReduction = 0.8;
    }
    

    for (unsigned long pIndex = 0; pIndex < indexInd.lLength; pIndex++) {
        _Variable* cv        = GetIthIndependentVar(pIndex);
        _Parameter thisLB    = cv->GetLowerBound(),
                   thisUB    = cv->GetUpperBound(),
                   thisValue = cv->Compute()->Value();

        //parameterTransformationFunction << _hyphyIntervalMapID;
        if (thisLB >= 0.0 && thisUB <= 1.0) {
            parameterTransformationFunction << _hyphyIntervalMapID;
        } else if (thisLB >=0.0) {
            parameterTransformationFunction << _hyphyIntervalMapSqueeze;
        } else {
            parameterTransformationFunction << _hyphyIntervalMapExpit;
        }


        parameterValuesAndRanges->Store(pIndex,0,thisValue);
        parameterValuesAndRanges->Store(pIndex,1,mapParameterToInverval(thisValue,parameterTransformationFunction.Element(-1),false));
        parameterValuesAndRanges->Store(pIndex,2,mapParameterToInverval(thisLB,parameterTransformationFunction.Element(-1),false));
        parameterValuesAndRanges->Store(pIndex,3,mapParameterToInverval(thisUB,parameterTransformationFunction.Element(-1),false));
    }

}

//_______________________________________________________________________________________________

void _LikelihoodFunction::CleanupParameterMapping (void)
{
    smoothingPenalty = 0.0;
    smoothingTerm    = 0.0;
    DeleteObject (parameterValuesAndRanges);
    parameterValuesAndRanges = nil;
    parameterTransformationFunction.Clear();
}



//_______________________________________________________________________________________________

_Parameter _LikelihoodFunction::SumUpSiteLikelihoods (long index, const _Parameter * patternLikelihoods, const _SimpleList& patternScalers)
/*
 compute the likelihood of a partition (index), corrected for scaling,
 by summing pattern likelihoods from patternLikelihoods, weighted by pattern frequencies
 and corrected for scaling factors from patternScalers
*/
{

    _Parameter       logL             = 0.;
    _SimpleList      *catVarType      = (_SimpleList*)((*(_List*)categoryTraversalTemplate(index))(4));
    long             cumulativeScaler = 0,
                     categoryType     = catVarType->Element (-1);

    _SimpleList     * patternFrequencies = &((_DataSetFilter*)dataSetFilterList (theDataFilters(index)))->theFrequencies;

    // check to see if we need to handle HMM or COP variables
    if (categoryType & _hyphyCategoryHMM) {
        _CategoryVariable*hmmVar = (_CategoryVariable*)((*(_List*)(*(_List*)categoryTraversalTemplate(index))(0))(0));
        _Matrix          *hmm    = hmmVar->ComputeHiddenMarkov(),
                          *hmf    = hmmVar->ComputeHiddenMarkovFreqs();

        _SimpleList      *dmap   = &((_DataSetFilter*)dataSetFilterList (theDataFilters(index)))->duplicateMap;

        return           SumUpHiddenMarkov (patternLikelihoods,
                                            *hmm,
                                            *hmf,
                                            dmap,
                                            &patternScalers,
                                            patternFrequencies->lLength
                                           );
    } else {
        if (categoryType & _hyphyCategoryCOP) {
            WarnError ("Constant-on-partition categories are currently not supported by the evaluation engine");
        } else // simple sum clause

        {
            for              (long patternID = 0; patternID < patternFrequencies->lLength; patternID++) {
                long patternFrequency = patternFrequencies->lData[patternID];
                if (patternFrequency > 1) {
                    logL             += myLog(patternLikelihoods[patternID])*patternFrequency;
                    cumulativeScaler += patternScalers.lData[patternID]*patternFrequency;
                } else
                    // all this to avoid a double*long multiplication
                {
                    logL             += myLog(patternLikelihoods[patternID]);
                    cumulativeScaler += patternScalers.lData[patternID];
                }
            }
        }
    }

    return logL - cumulativeScaler * _logLFScaler;

}

//_______________________________________________________________________________________________
// return the AVL with parameters
// AVL will have the following entries
// "Categories"
// "Global Independent"
// "Global Constrained"
// "Local Independent"
// "Local Constrained"
// "Trees"
// "Models"
// "Base frequencies"
// "Datafilters"
// "Compute Template"

_AssociativeList* _LikelihoodFunction::CollectLFAttributes (void)
{
    _AssociativeList * resList = new _AssociativeList;

    _SimpleList         *vl,
                        list1;

    _List               modelList;

    InsertVarIDsInList (resList, "Categories", GetCategoryVars ());

    SplitVariableIDsIntoLocalAndGlobal (GetIndependentVars (), modelList);
    InsertVarIDsInList (resList, "Global Independent", *(_SimpleList*)modelList(0));
    InsertVarIDsInList (resList, "Local Independent",   *(_SimpleList*)modelList(1));

    SplitVariableIDsIntoLocalAndGlobal (GetDependentVars (), modelList);
    InsertVarIDsInList (resList, "Global Constrained", *(_SimpleList*)modelList(0));
    InsertVarIDsInList (resList, "Local Constrained",   *(_SimpleList*)modelList(1));


    list1.Clear();
    vl = &GetTheTrees ();
    modelList.Clear();

    for (long n=0; n<vl->lLength; n++) {
        list1 << vl->lData[n];
        _SimpleList partModels;
        ((_TheTree*)FetchVar (vl->lData[n]))->CompileListOfModels(partModels);
        if (partModels.lLength == 1) {
            modelList << modelNames (partModels.lData[0]);
        } else {
            modelList.AppendNewInstance(new _String ("__MULTIPLE__"));
        }
    }
    InsertVarIDsInList (resList, "Trees", list1);


    list1.Clear();
    vl = &GetTheFilters ();
    for (long p=0; p<vl->lLength; p++) {
        list1 << vl->lData[p];
    }

    InsertStringListIntoAVL (resList, "Datafilters", list1, dataSetFilterNamesList);
    InsertVarIDsInList (resList, "Base frequencies", GetBaseFreqs());
    {
        _SimpleList indexer         (modelList.lLength,0,1);
        InsertStringListIntoAVL     (resList, "Models", indexer, modelList);
    }

    _Formula        *computeT = HasComputingTemplate();
    resList->MStore (_String("Compute Template"), new _FString((_String*)(computeT?computeT->toStr():new _String)), false);

    return resList;
}

//_______________________________________________________________________________________________

void _LikelihoodFunction::UpdateBlockResult (long index, _Parameter new_value)
{
    if (computationalResults.GetUsed()>index) {
        computationalResults.theData[index] = new_value;
    } else {
        computationalResults.Store(new_value);
    }
}



//_______________________________________________________________________________________________
#ifdef __HYPHYMPI__
long RetrieveMPICount (char)
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

void MPISwitchNodesToMPIMode (long totalNodeCount)
{
    _String message = mpiLoopSwitchToOptimize & hyphyMPIOptimizerMode;

    // send a context switch signal
    for (long ni = 1; ni <= totalNodeCount; ni++) {
        MPISendString (message, ni);
    }
    // receive confirmation of successful switch
    for (long ni = 1; ni <= totalNodeCount; ni++) {
        long fromNode = ni;
        _String t (MPIRecvString (ni,fromNode));
        if (!t.Equal (&mpiLoopSwitchToOptimize)) {
            WarnError (_String("[MPI] Failed to confirm MPI mode switch at node ") & ni);
            return;
        } else {
            ReportWarning (_String("[MPI] Successful mode switch to mode ") & hyphyMPIOptimizerMode & " confirmed from node " & ni);
        }
    }
}

#endif
#endif
