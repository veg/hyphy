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

  //#define _UBER_VERBOSE_LF_DEBUG

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <math.h>


#include "likefunc.h"
#include "calcnode.h"
#include "site.h"
#include "batchlan.h"
#include "category.h"
#include "function_templates.h"
#include "global_object_lists.h"
#include "global_things.h"
#include "time_difference.h"
#include "scfg.h"
#include "tree_iterator.h"
#include "vector.h"

using namespace hyphy_global_objects;
using namespace hy_global;

  //#define _UBER_VERBOSE_LF_DEBUG 1

//#define    _COMPARATIVE_LF_DEBUG_DUMP
//#define    _COMPARATIVE_LF_DEBUG_CHECK

#if defined _COMPARATIVE_LF_DEBUG_CHECK
#include <signal.h>
    unsigned long  _comparative_lf_index = 0UL;
    extern _Matrix* _comparative_lf_debug_matrix;
#endif

#if defined _COMPARATIVE_LF_DEBUG_DUMP
    unsigned long  _comparative_lf_index = 0UL;
    extern _GrowingVector* _comparative_lf_debug_matrix;
#endif




void         BenchmarkThreads (_LikelihoodFunction*);

long         siteEvalCount  =   0;




#ifdef MDSOCL
	_OCLEvaluator *OCLEval;
#endif




#ifdef  __HYPHYMPI__

_List           parallelOptimizerTasks;

_Matrix         varTransferMatrix,
                resTransferMatrix;

long            MPICategoryCount,
                transferrableVars,
                hyphyMPIOptimizerMode = _hyphyLFMPIModeNone;


_String         mpiLoopSwitchToOptimize ("_CONTEXT_SWITCH_MPIPARTITIONS_"),
                mpiLoopSwitchToBGM      ("_BGM_SWITCH_");

#endif

#define     SQR(A) (A)*(A)
#define     GOLDEN_RATIO 1.618034
#define     GOLDEN_RATIO_R  0.61803399
#define     GOLDEN_RATIO_C  (1-GOLDEN_RATIO_R)

#define     _HY_SLOW_CONVERGENCE_RATIO 2.
#define     _HY_SLOW_CONVERGENCE_RATIO_INV 1./_HY_SLOW_CONVERGENCE_RATIO



#define   PERTURBATION_OF_ZERO    0.0

long      likeFuncEvalCallCount = 0L,
          lockedLFID         = -1L;

#define  STD_GRAD_STEP 1.0e-6




// some string constants

_String  const
                                kOptimizationHardLimit          ("OPTIMIZATION_TIME_HARD_LIMIT"),
                                kMaximumIterationsPerVariable   ("MAXIMUM_ITERATIONS_PER_VARIABLE"),
                                kUseLastResults                 ("USE_LAST_RESULTS"),
                                // the hard limit, in seconds, for the an optimizer to run
                                kGlobalStartingPoint            ("GLOBAL_STARTING_POINT"),
                                kRandomStartingPerturbations    ("RANDOM_STARTING_PERTURBATIONS"),
                                kStartingPrecision              ("STARTING_PRECISION"),
                                kAllowBoundary                  ("ALLOW_BOUNDARY"),
                                kBracketingPersistence          ("BRACKETING_PERSISTENCE"),
                                kIntermediatePrecision          ("INTERMEDIATE_PRECISION"),
                                keepOptimalOrder                ("KEEP_OPTIMAL_ORDER"),
                                optimizeSummationOrder          ("OPTIMIZE_SUMMATION_ORDER"),
                                optimizePartitionSize           ("OPTIMIZE_SUMMATION_ORDER_PARTITION"),
                                likefuncOutput                  ("LIKELIHOOD_FUNCTION_OUTPUT"),
                                categorySimulationMethod        ("CATEGORY_SIMULATION_METHOD"),
                                kUseInitialDistanceGuess        ("USE_DISTANCES"),
                                covariancePrecision             ("COVARIANCE_PRECISION"),
                                cacheSubtrees                   ("CACHE_SUBTREES"),
                                likeFuncCountVar                ("LF_CALL_COUNT"),
                                doShuffleOrder                  ("SHUFFLE_ORDER_OF_PARAMETERS"),
                                forceDistanceEstimates          ("FORCE_DISTANCE_ESTIMATES"),
                                useDuplicateMatrixCaching       ("USE_DUPLICATE_MATRIX_CACHING"),
                                useFullMST                      ("USE_MST_HEURISTIC"),
                                stateCountMatrix                ("STATE_COUNT_MATRIX"),
                                wStateCountMatrix               ("WSTATE_COUNT_MATRIX"),
                                allowSequenceMismatch           ("ALLOW_SEQUENCE_MISMATCH"),
                                shortMPIReturn                  ("SHORT_MPI_RETURN"),
                                mpiPrefixCommand                ("MPI_PREFIX_COMMAND"),
                                kSkipConjugateGradient          ("SKIP_CONJUGATE_GRADIENT"),
                                useIntervalMapping              ("USE_INTERVAL_MAPPING"),
                                intervalMappingMethod           ("INTERVAL_MAPPING_METHOD"),
                                useAdaptiveVariableStep         ("USE_ADAPTIVE_VARIABLE_STEP"),
                                storeRootSupportFlag            ("STORE_ROOT_SUPPORT"),
                                supportMatrixVariable           ("SUPPORT_MATRIX_LIST"),
                                optimizationStatusFile          ("SAVE_OPT_STATUS_TO"),
                                autoParalellizeLF               ("AUTO_PARALLELIZE_OPTIMIZE"),
                                lfExtraLFExportCode             ("LF_NEXUS_EXPORT_EXTRA"),
                                optimizationStringTemplate      ("OPTIMIZATION_PROGRESS_TEMPLATE"),
                                kMaximumIterations              ("MAXIMUM_OPTIMIZATION_ITERATIONS"),
                                // use
                                // $1 for status
                                // $2 for log L
                                // $3 for percent done
                                // $4 for time elapsed
                                // $5 for evals/second
                                // $6 for CPU load
                                optimizationStringStatus        ("OPTIMIZATION_PROGRESS_STATUS"),
                                optimizationStringQuantum       ("OPTIMIZATION_PROGRESS_QUANTUM"),
                                categoryMatrixScalers           (".site_scalers"),
                                categoryLogMultiplier           (".log_scale_multiplier"),
                                minimumSitesForAutoParallelize  ("MINIMUM_SITES_FOR_AUTO_PARALLELIZE"),
                                userSuppliedVariableGrouping    ("PARAMETER_GROUPING"),
                                kAddLFSmoothing                 ("LF_SMOOTHING_SCALER"),
                                kOptimizationPrecision          ("OPTIMIZATION_PRECISION"),
                                kOptimizationMethod             ("OPTIMIZATION_METHOD"),
                                kReduceLFSmoothing              ("LF_SMOOTHING_REDUCTION");




void        countingTraverse         (node<long>*, long&, long, long&, bool);
void        countingTraverseArbRoot  (node<long>*, node<long>*, long&, long, long&);
long        findAvailableSlot        (_SimpleList&, long&);
void        setComputingArrays       (node<long>*, node<long>*, _SimpleList&, _SimpleList&, _SimpleList &, _SimpleList&, _SimpleList&, long&);


_SimpleList Fibonacci;



hyFloat  precision,
            optimizationPrecMethod,
            maxItersPerVar,
            dFPrintFormat,
            dFDefaultWidth,
            assignedSeedVal = -1.0,
            categorySimMethod;

bool        forceRecomputation = false,
            isInOptimize       = false,
            usedCachedResults  = false;

long        bracketFCount   = 0,
            bracketCount   = 0,
            oneDFCount       = 0,
            oneDCount         = 0,
            categID       = 0,
            offsetCounter   = 1;

extern      _List
dataSetList,
likeFuncList;

_Variable*  siteWiseVar     = nil,
            *  blockWiseVar = nil;


_String  const *  progressFileString = nil;

node<long>* DepthWiseStepTraverserLevel  (long&, node<long>* root);
hyFloat  myLog (hyFloat);


//__________________________________________________________________________________


void         BenchmarkThreads (_LikelihoodFunction* lf) {
    long         alterIndex = 0;

    if      (lf->HasComputingTemplate()) {
        for (unsigned long k=0; k<lf->GetIndependentVars().lLength; k++)
            if (!LocateVar (lf->GetIndependentVars().lData[k])->IsGlobal()) {
                alterIndex = k;
                break;
            }
    }



#ifdef  _OPENMP
    lf->SetThreadCount (1);
#endif
    TimeDifference timer;
    lf->SetIthIndependent (alterIndex,lf->GetIthIndependent(alterIndex));
    lf->Compute           ();



hyFloat            tdiff = timer.TimeSinceStart();
#ifdef  _OPENMP
#ifdef  __HYPHYMPI__ 
    if ( hy_global::system_CPU_count > 1 && hy_mpi_node_rank == 0) {
#else
    if ( hy_global::system_CPU_count > 1) {
#endif

        hyFloat          minDiff = tdiff;
        long                bestTC  = 1;

        for (long k = 2; k <= hy_global::system_CPU_count; k++) {
            lf->SetThreadCount              (k);
            TimeDifference timer;
            lf->SetIthIndependent           (alterIndex,lf->GetIthIndependent(alterIndex));
            lf->Compute                     ();
            tdiff = timer.TimeSinceStart();
            if (tdiff < minDiff) {
                minDiff = tdiff;
                bestTC  = k;
            } else {
                break;
            }
        }
        lf->SetThreadCount (bestTC);
        ReportWarning       (_String("Auto-benchmarked an optimal number (") & bestTC & ") of threads.");
    } 
#endif
}


#if defined  __UNIX__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__

void        UpdateOptimizationStatus    (hyFloat, long, char, bool, _String * fileName = nil);
    
//_______________________________________________________________________________________

_CustomFunction::_CustomFunction (_String const & arg, _VariableContainer const * context ) {
    _String body        (arg),
            error ;

    _FormulaParsingContext fpc (&error, context);

    if (Parse (&myBody, body, fpc, nil) == HY_FORMULA_EXPRESSION) {
        _SimpleList myVars;
        _AVLList al (&myVars);
        myBody.ScanFForVariables(al,true,false,false);
        al.ReorderList();

        for (unsigned long k=0UL; k<myVars.lLength; k++) {
            if (LocateVar(myVars.lData[k])->IsIndependent()) {
                indexInd << myVars.lData[k];
            }
        }
    } else {
        throw (_String ("Error while parsing ") & arg.Enquote() & " in CustomFunction " & error);
    }
}
    
//__________________________________________________________________________________

class  _OptimiztionProgress {
private:
    hyFloat bigLastMax,
            sigma,
            s,
            sX,
            sXX,
            sY,
            sYY,
            sXY,
            nPercentDone,
            percentDone;
    
    long    loopCounter;
public:
    _OptimiztionProgress (void) {
        loopCounter = 0L;
        bigLastMax  = 0.;
        sigma = 0.;
        s = 0.;
        sX = 0.;
        sXX = 0.;
        sY = 0.;
        sYY = 0.;
        sXY = 0.;
        nPercentDone = 0.;
        percentDone = 0.;
    }
    
    hyFloat PushValue (hyFloat max) {
        loopCounter++;
        if (loopCounter == 1L) {
            bigLastMax = max;
            return get();
        }
        nPercentDone = max - bigLastMax;
        if (nPercentDone == 0.) {
            nPercentDone = kMachineEpsilon;
        }
        sigma = 2.0/nPercentDone;
        s+= sigma;
        sX+= sigma*(loopCounter);
        sXX += sigma*(loopCounter)*(loopCounter);
        sY+=nPercentDone*sigma;
        sXY+=(loopCounter)*nPercentDone*sigma;

        if (loopCounter>1.0) {
            nPercentDone = (sX*sY-s*sXY);
            if (nPercentDone != 0.0) {
                nPercentDone = (sXX*sY-sX*sXY)/nPercentDone;
            }
            if (nPercentDone != 0.0) {
                nPercentDone = loopCounter/nPercentDone;
            }
            nPercentDone -= .05;
            if ((nPercentDone>0.0)&&(nPercentDone<.9999)) {
                nPercentDone = pow (nPercentDone,4)*100.0;
                if (nPercentDone>percentDone) {
                    percentDone = nPercentDone;
                } else {
                    if (percentDone<50.0) {
                        percentDone+=1.0;
                    }
                }
            } else if (percentDone<50.0) {
                percentDone+=1.0;
            }
        }
        bigLastMax = max;
        return get();
    }
    
    hyFloat  get (void) const {
        return percentDone;
    }
};
    
    
//__________________________________________________________________________________

void        UpdateOptimizationStatus (hyFloat max, long pdone, char init, bool optimization, _String const* fileName) {
    static long     lCount;
    static long     lastDone;
    static double   elapsed_time = 0.0;
    static hyFloat
    update_quantum = 0.0;
    static _String  userReportString;
    static _String  userStatusString;

    static clock_t  userTimeStart;
    FILE           *outFile = fileName?doFileOpen (fileName->get_str(),"w"):nil;
    _FString*       t;

    static          TimeDifference timer;

    if (init==0) {
        lCount          = likeFuncEvalCallCount;
        timer.Start();
#ifndef _MINGW32_MEGA_
        setvbuf           (stdout,nil, _IONBF,1);
#endif
        lastDone        = 0;
        userTimeStart   = clock();
        checkParameter   (optimizationStringQuantum, update_quantum, 0.0);
        t = (_FString*)FetchObjectFromVariableByType (&optimizationStringTemplate,STRING);
        userReportString = t?t->get_str():kEmptyString;
        t = (_FString*)FetchObjectFromVariableByType (&optimizationStringStatus,STRING);
        userStatusString = t?t->get_str():kEmptyString;
        elapsed_time     = 0.0;
    } else if (init==1) {
        double timeDiff = timer.TimeSinceStart();

        //printf ("%g %g\n", timeDiff,elapsed_time);

        if (pdone<0) {
            pdone = lastDone;
        }
        lastDone = pdone;

        if (timeDiff == 0.0 || timeDiff < update_quantum) {
            return;
        } else {
            elapsed_time += timeDiff;
            timer.Start();
        }


        if (userReportString.nonempty()) {
            char buffer[255];

            _String reportString = userReportString.Replace ("$1",userStatusString, true);
            if (optimization) {
                snprintf (buffer, sizeof(buffer), "%15.10g", (double)max);
                reportString = reportString.Replace ("$2", buffer, true);
            } else {
                reportString = reportString.Replace ("$2", kEmptyString, true);
            }
            reportString = reportString.Replace ("$3", _String(pdone), true);
            _String       tStamp (_String::FormatTimeString(elapsed_time));
            reportString = reportString.Replace ("$4",tStamp, true);
            if (elapsed_time) {
                snprintf (buffer,sizeof(buffer),"%8.4g", (clock()-userTimeStart)/((hyFloat)CLOCKS_PER_SEC*(elapsed_time)));
                reportString = reportString.Replace ("$6", buffer, true);
                snprintf (buffer, sizeof(buffer), "%8.4g", (likeFuncEvalCallCount-lCount)/elapsed_time);
                reportString = reportString.Replace ("$5", buffer, true);
            } else {
                reportString = reportString.Replace ("$5", kEmptyString, true);
                reportString = reportString.Replace ("$6", kEmptyString, true);
            }

            if (outFile) {
                fprintf (outFile,"%s", reportString.get_str());
            } else
#ifndef _MINGW32_MEGA_
                printf ("\015%s", reportString.get_str());
#else
                SetStatusLine (reportString);
#endif
        } else {
            char buffer [1024];
            if (optimization) {
                if (outFile)
                    fprintf (outFile,"Likelihood function optimization status\nCurrent Maximum: %-14.8g (%ld %% done)\nLikelihood Function evaluations/second: %-8.4g", (double)max, pdone,
                             (likeFuncEvalCallCount-lCount)/elapsed_time);
                else {
                    long written = snprintf (buffer,1024,"Current Max: %-14.8g (%ld %% done) LF Evals/Sec: %-8.4g", (double)max, pdone, (likeFuncEvalCallCount-lCount)/elapsed_time);

                    if (elapsed_time) {
                        snprintf (buffer+written,1024-written, "CPU Load: %-8.4g", (clock()-userTimeStart)/((hyFloat)CLOCKS_PER_SEC*elapsed_time));
                    }
                }
            } else {
                snprintf (buffer, 1024, "Sites done: %g (%ld %% done)", (double)max, pdone);
            }

#ifndef _MINGW32_MEGA_
            printf ("\015%s", buffer);
#else
            SetStatusLine (_String(buffer));
#endif
        }


    } else {
        if (outFile) {
            fprintf (outFile,"DONE");
        } else {
#ifndef _MINGW32_MEGA_
            printf ("\033\015 ");
            setvbuf (stdout,nil,_IOLBF,1024);
#endif
        }
    }
    if (outFile) {
        fclose (outFile);
    }
}

#else


#endif

//__________________________________________________________________________________

hyFloat myLog (hyFloat arg)
{
    return (arg>0.0)?log(arg):-1000000.;
}


//_______________________________________________________________________________________

_LikelihoodFunction::_LikelihoodFunction (void)
{
    Init();
}

//_______________________________________________________________________________________

void _LikelihoodFunction::Init (void)
{
    siteResults         = nil;
    bySiteResults       = nil;
    hasBeenOptimized    = false;
    hasBeenSetUp        = 0;
    templateKind        = _hyphyLFComputationalTemplateNone;
    computingTemplate   = nil;
    mstCache            = nil;
    nonConstantDep      = nil;
    evalsSinceLastSetup = 0;
    siteArrayPopulated  = false;
    smoothingTerm       = 0.;
    smoothingPenalty    = 0.;

    conditionalInternalNodeLikelihoodCaches = nil;
    conditionalTerminalNodeStateFlag        = nil;
    siteScalingFactors                      = nil;
    branchCaches                            = nil;
    parameterValuesAndRanges                = nil;
    optimizatonHistory                      = nil;

#ifdef  _OPENMP
    lfThreadCount       = 1L;
#ifdef  __HYPHYMPI__
    if (hy_mpi_node_rank > 0)
        SetThreadCount      (1L);
    else
#endif

 SetThreadCount      (system_CPU_count);
#endif

}

//_______________________________________________________________________________________

_LikelihoodFunction::_LikelihoodFunction (_String& s, _VariableContainer* p)
// from triplets
//format: datasetfilter name, tree name, frequency matrix name; etc...
{
    Init     ();
    _List    tripletsRaw (&s,';'),
             tripletsSplit;
    for (unsigned long k = 0; k < tripletsRaw.lLength; k++) {
        _List thisTriplet (tripletsRaw(k),',');
        tripletsSplit << thisTriplet;
    }
    Construct(tripletsSplit,p);
}

//_______________________________________________________________________________________

_TheTree * _LikelihoodFunction::GetIthTree (long f) const {
  _Variable *tree_var = LocateVar (theTrees.lData[f]);
  if (tree_var && tree_var->ObjectClass() == TREE) {
    return (_TheTree*)tree_var;
  }
  return nil;
}

  //_______________________________________________________________________________________

_DataSetFilter const * _LikelihoodFunction::GetIthFilter (long f) const {
  return GetDataFilter (theDataFilters.lData[f]);
}

//_______________________________________________________________________________________

_String const* _LikelihoodFunction::GetIthFilterName (long f) const {
  return GetFilterName(theDataFilters.lData[f]);
}

//_______________________________________________________________________________________

_DataSetFilter * _LikelihoodFunction::GetIthFilterMutable (long f) const {
  return (_DataSetFilter *) GetDataFilter (theDataFilters.lData[f]);
}


//_______________________________________________________________________________________

_Matrix * _LikelihoodFunction::GetIthFrequencies (long f) const {
    return (_Matrix *)FetchObjectFromVariableByTypeIndex(theProbabilities.Element(f), MATRIX);
}

//_______________________________________________________________________________________

_String const * _LikelihoodFunction::GetIthFrequenciesName (long f) const {
  return LocateVar(theProbabilities.Element(f))->GetName ();
}

//_______________________________________________________________________________________


bool    _LikelihoodFunction::MapTreeTipsToData (long f, _String *errorMessage, bool leafScan) { // from triplets
  
    _TheTree*       t = GetIthTree(f);

    _TreeIterator   ti (t, _HY_TREE_TRAVERSAL_POSTORDER | fTreeIteratorTraversalSkipRoot);

    _DataSetFilter  * df = GetIthFilterMutable (f);
    long            dfDim = df->GetDimension(true);

    _List           tips;
    
    try {

        while (_CalcNode* iterator = ti.Next()) {
            if (ti.IsAtLeaf()) {
                tips.AppendNewInstance (new _String (iterator->ContextFreeName ()));
            }
            if (iterator->GetModelIndex () == HY_NO_MODEL) {
                throw _String ("Model is not associated with the node:") & iterator->ContextFreeName();
            } else if (iterator->GetModelDimension() != dfDim) {
                throw _String ("The dimension of the transition matrix at node ") & iterator->ContextFreeName ().Enquote() & " is not equal to the state count in the data filter associated with the tree.";
            }
        }

        // now that "tips" contains all the names of tree tips we can
        // scan thru the names in the datafilter and check whether there is a 1-1 match
        
        if ((t->IsDegenerate()?2:tips.lLength) != df->NumberSpecies()) {
            throw _String("The number of tree tips in ")&t->GetName()->Enquote () & _String((long)tips.lLength).Enquote('(', ')')
            & " is not equal to the number of sequences in the data filter associated with the tree " &
            '(' & _String((long)df->NumberSpecies()) & ").";
        }

        if (!t->IsDegenerate()) {
            long        j,
                        k;

            _SimpleList tipMatches;
            // produce a sorted list of sequence names

            j = df->FindSpeciesName (tips, tipMatches);

            if (j != tips.lLength) { // try numeric match
                if (hy_env::EnvVariableTrue(hy_env::try_numeric_sequence_match)) {

                    try {
                        tipMatches.Clear();


                        for (j=0L; j<tips.lLength; j++) {
                            _String *thisName = (_String*)tips(j);
                            long numeric_value = thisName->to_long();
                            if (numeric_value<=tips.lLength && numeric_value >= 0L && _String(numeric_value) == *thisName) {
                                tipMatches<<numeric_value;
                            } else {
                                throw (j);
                            }
                        }

                        if (j==tips.lLength) {
                            if (tipMatches.Find(0L) < 0) // map to indexing from 0
                                tipMatches.Offset (-1L);

                            _SimpleList const *dfMap = (_SimpleList const*)df->GetMap();

                            if (dfMap) {
                                for (unsigned long k = 0UL; k < tips.lLength; k++) {
                                    tipMatches.lData[k] = dfMap->lData[tipMatches.lData[k]];
                                }
                            }
                        }
                    } catch (unsigned long i) {
                        j = -1L;
                    }
                }
            }

            if (j==tips.lLength) { // all matched
                /*
                    20100913: SLKP need to check that reusing the datafilter will not mess up existing likelihood function dependendancies
                */

                _SimpleList * currentMap = (_SimpleList *)df->GetMap();
                if (! currentMap || ! currentMap->Equal (tipMatches)) {
                    for (unsigned long lfID = 0UL; lfID < likeFuncList.lLength; lfID++) {
                        _LikelihoodFunction* lfp = (_LikelihoodFunction*)likeFuncList(lfID);
                        if (lfp && lfp != this && lfp->DependOnDF (theDataFilters.lData[f])) {
                           throw _String ("Cannot reuse the filter '") & *GetObjectNameByType (HY_BL_DATASET_FILTER, theDataFilters.lData[f], false) &
                                       "' because it is already being used by likelihood function '" &
                                       *GetObjectNameByType (HY_BL_LIKELIHOOD_FUNCTION, lfID, false) & "', and the two likelihood functions impose different leaf-to-sequence mapping. " &
                                       "Create a copy the filter and pass it to the second likelihood function to resolve this issue.";
                        }
                    }
                    df->SetMap(tipMatches);
                }
                ReportWarning (_String ("The tips of the tree:") & *t->GetName() &" were matched with the species names from the data in the following numeric order (0-based) "& _String ((_String*)tipMatches.toStr()));
            } else {
                _String warnMsg = _String ("The leaf of the tree ") & *t->GetName() &" labeled " &(*(_String*)tips(j)).Enquote()
                                  &" had no match in the data filter. Please make sure that all leaf names correspond to a sequence name in the data file.";

                if (! hy_env::EnvVariableTrue(allowSequenceMismatch)) {
                    throw (warnMsg);
                }
                ReportWarning (warnMsg);
            }
        }
        if (leafScan) {
            ((_SimpleList*)leafSkips(f))->Clear();
            df->MatchStartNEnd(*(_SimpleList*)optimalOrders(f),*(_SimpleList*)leafSkips(f));
        }
    } catch (const _String err) {
        HandleOrStoreApplicationError (errorMessage, err);
        return false;
        
    }
    return true;

}

//_______________________________________________________________________________________

void     _LikelihoodFunction::Rebuild (bool rescan_parameters) {
  computationalResults.Clear();
  hasBeenSetUp     = 0;
  hasBeenOptimized = false;
  _String ignored_error;
  try {
    for (unsigned long k = 0UL; k < theDataFilters.lLength; k++) {
      if (! (GetIthFilter (k) && GetIthTree (k) && GetIthFrequencies(k) && CheckIthPartition(k, &ignored_error))) {
        throw (k);
      }
    }
  }
  catch (unsigned long code) {
    ReportWarning (_String ("Likelihood function cleared because partition index '") & (long) code & "' points to invalid components");
    Clear();
    return;
  }
  AllocateTemplateCaches();
  Setup(false);
  if (rescan_parameters) {
    RescanAllVariables();
  }
}

//_______________________________________________________________________________________

void     _LikelihoodFunction::Clear (void)
{
    DeleteCaches  ();

    //unsigned long partition_count = CountObjects(kLFCountPartitions);

    UnregisterListeners ();
    theTrees.Clear();

    theDataFilters.Clear();
    theProbabilities.Clear();
    indexInd.Clear();
    indexDep.Clear();
    indexCat.Clear();
    blockDependancies.Clear();
    computationalResults.Clear();
    partScalingCache.Clear();
    indVarsByPartition.Clear();
    depVarsByPartition.Clear();


    optimalOrders.Clear();
    leafSkips.Clear();
    hasBeenSetUp            = 0;
    hasBeenOptimized        = false;
    if (computingTemplate) {
        delete computingTemplate;
        computingTemplate = nil;
        templateKind      = _hyphyLFComputationalTemplateNone;
    }
    if (mstCache) {
        delete (mstCache);
        mstCache = nil;
    }

    if (optimizatonHistory) {
      DeleteObject(optimizatonHistory);
      optimizatonHistory = nil;
    }

    treeTraversalMasks.Clear();
    canUseReversibleSpeedups.Clear();
#ifdef  _OPENMP
#ifdef  __HYPHYMPI__
    if (hy_mpi_node_rank > 0)
        SetThreadCount      (1L);
    else
#endif

    SetThreadCount      (system_CPU_count);

#endif
}


//_______________________________________________________________________________________

void     _LikelihoodFunction::AllocateTemplateCaches (void) {
  partScalingCache.Clear();
  DeleteObject(bySiteResults);

  if (templateKind < 0 || templateKind == _hyphyLFComputationalTemplateBySite) {

    long max_filter_size = 0L;
    for (unsigned long f=0UL; f<theDataFilters.lLength; f++) {
      StoreIfGreater(max_filter_size, GetIthFilter (f)->GetSiteCountInUnits());
    }


#ifdef __HYPHYMPI__
    bySiteResults = new _Matrix    (theTrees.lLength+3,max_filter_size,false,true);
#else
    bySiteResults = new _Matrix    (theTrees.lLength+2,max_filter_size,false,true);
#endif
    for (unsigned long k = 0UL; k <= theTrees.lLength; k++) {
      partScalingCache < new _SimpleList(max_filter_size, 0,0);
    }
  } else {
    bySiteResults = nil;
  }
}
//_______________________________________________________________________________________

bool     _LikelihoodFunction::CheckIthPartition(unsigned long partition, _String * errorString, _String const * df, _String const * tree, _String const * efv) {
  _DataSetFilter const*     filter = GetIthFilter (partition);

  long filter_dimension          = filter->GetDimension(true),
       freq_dimension           = GetIthFrequencies(partition)->GetHDim ();

  if (freq_dimension  != filter_dimension) {

    if (df && efv) {
      HandleOrStoreApplicationError(errorString,_String("The dimension of the equilibrium frequencies vector ") &
               efv->Enquote() & " (" & freq_dimension & ") doesn't match the number of states in the dataset filter (" & filter_dimension & ") " & df->Enquote());
    }
    else {
      HandleOrStoreApplicationError (errorString, _String ("Incompatible dimensions between the filter (") & filter_dimension & ") and the frequencies matrix (" & freq_dimension & ")");
    }
    return false;
  }

  if (filter->IsNormalFilter() == false) { // do checks for the numeric filter
    if (filter->NumberSpecies() != 3UL || filter_dimension != 4L) {
          HandleOrStoreApplicationError(errorString,_String ("Datafilters with numerical probability vectors must contain exactly three sequences and contain nucleotide data. Had ") & (long) filter->NumberSpecies() & " sequences on alphabet of dimension " & (long) filter_dimension & '.');
      return false;
    }
  }

  return MapTreeTipsToData (partition, errorString);

}

//_______________________________________________________________________________________

bool     _LikelihoodFunction::Construct(_List& triplets, _VariableContainer* theP) {
/* SLKP v2.0 code cleanup 20090316 */

/* modified the code to take arguments as a pre-partitioned list,
   instead of the string;
   this will make building LFs from matrices of strings possible */

// from triplets
// format: datasetfilter name, tree name, frequency matrix name; etc...

    Clear ();
    long i = 0L;

    for (; i< (long)triplets.lLength-2L; i+=3L) {
        _String  object_name;
        long     objectID;

        // add datasetfilter
        object_name = AppendContainerName (*(_String*)triplets(i), theP);
        objectID   = FindDataFilter (object_name);

      //printf ("[_LikelihoodFunction::Construct] %s / %s\n", object_name.sData, GetFilterName(objectID)->sData);

        if (objectID < 0) {
            HandleApplicationError (_String("Could not locate a datafilter ")& object_name.Enquote());
            return false;
        } else {
          if (!RegisterChangeListenerForDataFilter(objectID,this)) {
            HandleApplicationError (_String("Could register likelihood function as listener for ")& object_name.Enquote());
            return false;
          }
          theDataFilters<<objectID;
        }

        // add the tree
        object_name = AppendContainerName (*(_String*)triplets(i+1), theP);
        _TheTree   * treeVar = (_TheTree*)FetchObjectFromVariableByType (&object_name, TREE);
        if (!treeVar) {
            HandleApplicationError (_String("Could not locate a tree variable named ")& object_name.Enquote());
            return false;
        } else {
            theTrees<<treeVar->get_index();
        }

        // add the matrix of probabilities
        object_name = AppendContainerName (*(_String*)triplets(i+2), theP);
        objectID   = LocateVarByName(object_name);
        _Matrix*   efv              = (_Matrix*)FetchObjectFromVariableByTypeIndex(objectID, MATRIX);
        if (!efv) {
            HandleApplicationError (_String("Could not locate a frequencies matrix named ") & object_name.Enquote());
            return false;
        } else {
            theProbabilities<<variableNames.GetXtra(objectID);
        }

      if (!CheckIthPartition (theTrees.lLength-1L, nil, (_String*)triplets.GetItem(i), (_String*)triplets.GetItem(i+1), (_String*)triplets.GetItem(i+2))) {
        return false;
      }

    }
    if (i && i == triplets.lLength-1) {
        _String templateFormulaString (ProcessLiteralArgument((_String*)triplets(i),theP));

        if (templateFormulaString.nonempty()) {
            siteWiseVar  = CheckReceptacle (&hy_env::sitewise_matrix,kEmptyString),
            blockWiseVar = CheckReceptacle (&hy_env::blockwise_matrix,kEmptyString);

            _String    copyString         (templateFormulaString);
            // do this because _Formula constructor will consume the string parameter
            _Formula   templateFormula    (templateFormulaString,theP);


            if (templateFormula.IsEmpty()|| terminate_execution) {
                HandleApplicationError (copyString.Enquote()  & " is not a valid formula specification in call to LikelihoodFunction");
                Clear     ();
                return    false;
            }

            bool  hasSiteMx         = templateFormula.DependsOnVariable(siteWiseVar->get_index()),
                  hasBlkMx            = templateFormula.DependsOnVariable(blockWiseVar->get_index());

            templateKind = _hyphyLFComputationalTemplateNone;
            long            templateFormulaOpCount = templateFormula.NumberOperations();

            if ( (hasBlkMx && hasSiteMx) || !(hasBlkMx||hasSiteMx) ) {
                if (hasBlkMx||hasSiteMx == false) {
                    if (templateFormulaOpCount==1)
                        // potentially an HMM
                    {
                        _Operation * firstOp = templateFormula.GetIthTerm (0);
                        if (firstOp->IsAVariable(false)) {
                            _Variable * hmmVar = LocateVar(firstOp->GetAVariable());
                            if (hmmVar->IsCategory() && ((_CategoryVariable*)hmmVar)->IsHiddenMarkov()) {
                                templateKind = -hmmVar->get_index()-1;
                                hasSiteMx    = true;
                            }
                        }
                    } else if (templateFormulaOpCount>=2) { // user function
                        long         nOps   =  templateFormula.GetIthTerm (templateFormulaOpCount-1)->UserFunctionID();
                        if (nOps>=0) { // user defined function
                            templateKind = _hyphyLFComputationalTemplateByPartition+1+nOps;
                            hasSiteMx    = true;
                        }
                    }
                }
                if (templateKind == _hyphyLFComputationalTemplateNone)
                    // error condition here
                {
                    HandleApplicationError ( copyString & " must depend either on " & hy_env::sitewise_matrix.Enquote() & " or " & hy_env::blockwise_matrix.Enquote() & " (but not on both). Alternatively, it could be a reference to a HM category variable or a user defined BL function." );
                    Clear     ();
                    return false;
                }
            }

            // determine the longest filter

            long             maxFilterSize = 0;

            if (hasSiteMx) {
                if (templateKind == _hyphyLFComputationalTemplateNone) {
                    templateKind = _hyphyLFComputationalTemplateBySite;
                }

             } else {
                templateKind = hasBlkMx?_hyphyLFComputationalTemplateByPartition : _hyphyLFComputationalTemplateNone;
            }

            // now test evaluate the formula

            if (templateKind==_hyphyLFComputationalTemplateBySite || templateKind==_hyphyLFComputationalTemplateByPartition) {
                _Matrix         testMx  (theTrees.lLength,1,false,true),
                                testMx2 (theTrees.lLength,1,false,true);

                for (long testCount = 0; testCount<25; testCount++) {
                    for (unsigned long di = 0; di < theTrees.lLength; di++) {
                        testMx.theData[di]  = genrand_real1 ();
                        testMx2.theData[di] = -genrand_real1 ();
                    }

                    siteWiseVar->SetValue  (&testMx);
                    blockWiseVar->SetValue (&testMx2);

                    HBLObjectRef    testResult = templateFormula.Compute();
                    _String     errMessage;
                    if (!testResult || terminate_execution) {
                        errMessage = _String ("Failed to evaluate computation template formula (")& copyString & ") in LikelihoodFunction constructor.";
                    } else {
                        if (testResult->ObjectClass() == NUMBER) {
                            if (testResult->Value()>0.0 ) {
                                errMessage = _String ("Computation template formula (")& copyString & ") in LikelihoodFunction constructor evaluated to a positive value = " & testResult->Value() & " (as a log-likelihood - must be non-positive).";
                            }
                        } else {
                            errMessage = _String ("Computation template formula (")& copyString & ") in LikelihoodFunction constructor evaluated to a non-scalar value.";
                        }
                    }

                    if (errMessage.nonempty()) {
                        HandleApplicationError (errMessage);
                        Clear ();
                        return false;
                    }
                }
            }

            computingTemplate = (_Formula*)templateFormula.makeDynamic();
            AllocateTemplateCaches();
        }
    }

    if (theTrees.lLength>1) {
        _SimpleList      checkDupTreeIDs (theTrees);
        checkDupTreeIDs.Sort ();
        for (long i=1; i<checkDupTreeIDs.lLength; i++)
            if (checkDupTreeIDs.lData[i] == checkDupTreeIDs.lData[i-1]) {
                HandleApplicationError (_String("The same tree - ")&*LocateVar(checkDupTreeIDs.lData[i])->GetName() &
                           " - can't be used in multiple partitions in the same likelihood function. You should create an independent tree for each partition, and constrain the parameters instead.");
                return false;
            }
    } else {
        if (theTrees.lLength == 0) {
            HandleApplicationError ("Too few arguments in call to _LikelihoodFunction::Construct");
            Clear ();
            return false;

        }
    }


    ScanAllVariables();
    Setup();
    return true;
}

//_______________________________________________________________________________________

_LikelihoodFunction::_LikelihoodFunction (_LikelihoodFunction& lf) // stack copy
{
    Clear();

    hasBeenOptimized    = lf.hasBeenOptimized;
    templateKind        = lf.templateKind;

    if (lf.computingTemplate) {
        computingTemplate   = (_Formula*)lf.computingTemplate->makeDynamic();
    } else {
        computingTemplate   = nil;
    }

    mstCache        = nil;
    nonConstantDep  = nil;

    Duplicate (&lf);
}

//_______________________________________________________________________________________

BaseRef _LikelihoodFunction::makeDynamic (void) const  // dynamic copy of this object
{
    _LikelihoodFunction * res = new _LikelihoodFunction;
    res->Duplicate(this);
    return res;
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::Duplicate (BaseRefConst obj) // duplicate an object into this one
{
    _LikelihoodFunction const* lf = (_LikelihoodFunction const*)obj;
    theTrees.Duplicate(&lf->theTrees);
    theProbabilities.Duplicate(&lf->theProbabilities);
    theDataFilters.Duplicate(&lf->theDataFilters);
    indexInd.Duplicate(&lf->indexInd);
    indexDep.Duplicate(&lf->indexDep);
    indexCat.Duplicate(&lf->indexCat);
    blockDependancies.Duplicate(&lf->blockDependancies);
    computationalResults.Duplicate(&lf->computationalResults);
    siteResults = nil;

    optimalOrders.Duplicate(&lf->optimalOrders);
    leafSkips.Duplicate (&lf->leafSkips);
    templateKind        = lf->templateKind;

    if (lf->optimizatonHistory) {
      optimizatonHistory = new _AssociativeList;
      optimizatonHistory->Duplicate (lf->optimizatonHistory);
    } else {
      optimizatonHistory  = nil;
    }

    if (lf->computingTemplate) {
        computingTemplate   = (_Formula*)lf->computingTemplate->makeDynamic();
    } else {
        computingTemplate   = nil;
    }

    if (lf->mstCache) {
        mstCache = new MSTCache;

        mstCache->computingOrder.Duplicate(&lf->mstCache->computingOrder);
        mstCache->storageOrder.Duplicate  (&lf->mstCache->storageOrder);
        mstCache->referenceOrder.Duplicate(&lf->mstCache->referenceOrder);
        mstCache->parentOrder.Duplicate (&lf->mstCache->parentOrder);
        mstCache->resultCache.Duplicate (&lf->mstCache->resultCache);
        mstCache->statesCache.Duplicate (&lf->mstCache->statesCache);
        mstCache->cacheSize.Duplicate(&lf->mstCache->cacheSize);
    }

    if (lf->bySiteResults) {
        bySiteResults = (_Matrix*)lf->bySiteResults->makeDynamic();
    } else {
        bySiteResults = nil;
    }

    if (lf->nonConstantDep) {
        nonConstantDep = (_SimpleList*)lf->nonConstantDep->makeDynamic();
    } else {
        nonConstantDep = nil;
    }
}


//_______________________________________________________________________________________
_SimpleList const&    _LikelihoodFunction::GetIndependentVars (void) const {
    return indexInd;
}

//_______________________________________________________________________________________
_SimpleList const&    _LikelihoodFunction::GetDependentVars (void) const {
    return indexDep;
}

//_______________________________________________________________________________________
_SimpleList const&    _LikelihoodFunction::GetCategoryVars (void) const {
    return indexCat;
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::GetGlobalVars (_SimpleList& rec) const {
    _Variable*      thisV;
    long            k;

    for (k=0; k<indexInd.lLength; k++) {
        thisV = LocateVar (indexInd.lData[k]);
        if (thisV->IsGlobal()) {
            rec << indexInd.lData[k];
        }
    }
    for (k=0; k<indexDep.lLength; k++) {
        thisV = LocateVar (indexDep.lData[k]);
        if (thisV->IsGlobal()) {
            rec << indexDep.lData[k];
        }
    }
}

//_______________________________________________________________________________________
hyFloat  _LikelihoodFunction::GetIthIndependent (long index, bool map) const {
    hyFloat return_value;

    if (parameterValuesAndRanges && map) {
        return_value = (*parameterValuesAndRanges)(index,1);
    } else {
        return_value = ((_Constant*) LocateVar (indexInd.lData[index])->Compute())->Value();
    }
    if (isnan(return_value)) {
        HandleApplicationError(*GetIthIndependentName (index) & " evaluated to a NaN; this can cause all kinds of odd behavior downstream, therefore it is safer to quit now");
    }
    return return_value;
}

  //_______________________________________________________________________________________
const _String*  _LikelihoodFunction::GetIthIndependentName  (long index) const {
  return LocateVar (indexInd.lData[index])->GetName();
}

//_______________________________________________________________________________________
const _String*  _LikelihoodFunction::GetIthDependentName  (long index) const {
    return LocateVar (indexDep.lData[index])->GetName();
}

//_______________________________________________________________________________________

hyFloat  _LikelihoodFunction::GetIthIndependentBound      (long index, bool isLower) const{
    if (parameterValuesAndRanges) {
        return (*parameterValuesAndRanges)(index,isLower?2:3);
    }
    if (isLower) {
        return GetIthIndependentVar(index)->GetLowerBound();
    }
    return GetIthIndependentVar(index)->GetUpperBound();

}
//_______________________________________________________________________________________
hyFloat  _LikelihoodFunction::GetIthDependent (long index) const {
    return ((_Constant*) LocateVar (indexDep.lData[index])->Compute())->Value();
}

//_______________________________________________________________________________________
_Variable*  _LikelihoodFunction::GetIthIndependentVar (long index) const {
    return LocateVar (indexInd.lData[index]);
}

  //_______________________________________________________________________________________
_CategoryVariable*  _LikelihoodFunction::GetIthCategoryVar (long index) const {
  return (_CategoryVariable*)LocateVar (indexCat.lData[index]);
}

//_______________________________________________________________________________________
_Variable*  _LikelihoodFunction::GetIthDependentVar (long index) const {
    return LocateVar (indexDep.lData[index]);
}
//_______________________________________________________________________________________
void    _LikelihoodFunction::SetIthIndependent (long index, hyFloat p) {
    if (parameterValuesAndRanges) {
        parameterValuesAndRanges->Store(index,1,p);
        p = mapParameterToInverval(p,parameterTransformationFunction.Element(index),true);
        parameterValuesAndRanges->Store(index,0,p);
    }
    //printf ("%10.10g\n", p);
    _Variable * v =(_Variable*) LocateVar (indexInd.lData[index]);
    v->SetValue (new _Constant (p), false);
    if (parameterValuesAndRanges) {
      hyFloat check_value = v->Value();
      if (p != check_value) {
        parameterValuesAndRanges->Store(index,0,check_value);
      }
    }
}

//_______________________________________________________________________________________
bool    _LikelihoodFunction::IsIthParameterGlobal (long index) const {
    _Variable * v =(_Variable*) LocateVar (indexInd.lData[index]);
    return v->IsGlobal();
}


//_______________________________________________________________________________________
bool    _LikelihoodFunction::CheckAndSetIthIndependent (long index, hyFloat p)
{
    _Variable * v =(_Variable*) LocateVar (indexInd.lData[index]);

    bool set;

    if (parameterValuesAndRanges) {
        parameterValuesAndRanges->Store(index,1,p);
        p = mapParameterToInverval(p,parameterTransformationFunction.Element(index),true);
        parameterValuesAndRanges->Store(index,0,p);
    }

    hyFloat oldValue = v->Value();

    if (p != INFINITY) {
      if (p!=0.0) {
          set = (fabs((oldValue-p)/p))>kMachineEpsilon;
      } else {
          set = fabs(oldValue-p)>kMachineEpsilon;
      }
    } else {
      set = true;
    }

    if (set) {
        v->SetValue (new _Constant (p), false);
        /**
         SLKP : because 'p' may be moved back into parameter bounds,
         if parameterValuesAndRanges is being used, we may need to update that the set
         */
        /*if (parameterValuesAndRanges) {
          hyFloat check_value = v->Value();
          if (p != check_value) {
            parameterValuesAndRanges->Store(index,0,check_value);
          }
        }*/

    }

    return set;
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::SetIthDependent (long index, hyFloat p)
{
    _Variable * v =(_Variable*) LocateVar (indexDep.lData[index]);
    v->SetValue (new _Constant (p),false);
}

//_______________________________________________________________________________________
long    _LikelihoodFunction::SetAllIndependent (_Matrix* v)
{
    unsigned long upto = v->GetSize();
    long     set_this_many = 0;
    upto = MIN (upto, indexInd.lLength);
    for (long k = 0; k < upto; k++) {
        set_this_many += CheckAndSetIthIndependent (k, v->theData[k]);
    }
    return set_this_many;
}

//_______________________________________________________________________________________
_Matrix*    _LikelihoodFunction::RemapMatrix(_Matrix* source, const _SimpleList& partsToDo) const {
    long hDim               =   source->GetHDim(),
         vDim                =   0,
         offsetInSource      =   0,
         offsetInTarget      =   0;

    for (unsigned long i=0; i<partsToDo.lLength; i++) {
        vDim += GetIthFilter (partsToDo.lData[i])->GetSiteCountInUnits();
    }

    _Matrix* res = new _Matrix (hDim,vDim,false,true);

    for (long aPart = 0; aPart<partsToDo.lLength; aPart++) {
        long partIndex = partsToDo.lData[aPart];
        _DataSetFilter  const * dsf = GetIthFilter (partIndex);
        long filterSize = dsf->GetSiteCountInUnits();

        if (HasHiddenMarkov(blockDependancies.lData[partIndex])>=0)
            // do nothing, just copy
        {
            for (long rowIndex = 0; rowIndex < hDim; rowIndex++)
                for (long columnIndex = 0; columnIndex < filterSize; columnIndex++) {
                    res->Store (rowIndex, columnIndex + offsetInTarget, (*source)(rowIndex, columnIndex + offsetInSource));
                }

            offsetInSource  +=  filterSize;
        } else {
            for (long rowIndex = 0; rowIndex < hDim; rowIndex++)
                for (long columnIndex = 0; columnIndex < filterSize; columnIndex++) {
                    res->Store (rowIndex, columnIndex + offsetInTarget, (*source)(rowIndex, dsf->duplicateMap.lData[columnIndex] + offsetInSource));
                }

            offsetInSource  +=  BlockLength(partIndex);
        }
        offsetInTarget  +=  filterSize;
    }
    res->AmISparse();
    return res;

}

//_______________________________________________________________________________________

#ifdef __HYPHYMPI__
void        _LikelihoodFunction::MPI_LF_Compute (long senderID, bool partMode)
#else
void        _LikelihoodFunction::MPI_LF_Compute (long, bool)
#endif
{
#ifdef __HYPHYMPI__

    if (!partMode) {
        AllocateSiteResults ();
    }

    _Matrix       variableStash (indexInd.lLength,1,false,true);
    hyFloat    siteLL = 0.;

    MPI_Status    status;
    //ReportWarning (_String ("Waiting on:") & (long) indexInd.lLength & " MPI_DOUBLES");
    ReportMPIError(MPI_Recv(variableStash.theData, indexInd.lLength, MPI_DOUBLE, senderID, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD,&status),false);

    //printf ("[ENTER NODE] %d\n", hy_mpi_node_rank);

    while (variableStash.theData[0]>=-1e100)
        // variableStash.theData[0] < -1e100 terminates the computation
    {
        //ReportWarning (_String("In at step  ") & loopie);
        bool    doSomething = false;
        for (long i = 0; i<indexInd.lLength; i++) {
            _Variable *anInd = LocateVar(indexInd.lData[i]);
            //ReportWarning (*anInd->GetName() & " = " & variableStash.theData[i]);
            if (anInd->HasChanged() || !CheckEqual(anInd->Value(), variableStash.theData[i])) {
                doSomething = true;
                SetIthIndependent (i,variableStash.theData[i]);
            }
        }
        if (doSomething) {
            if (partMode) {
                siteLL = Compute();
                /*if (hy_mpi_node_rank == 1) {
                    printf ("\033\015 Node %ld, value %g", siteLL);
                }*/
            } else {

                if (PreCompute()) {
                    //ReportWarning ("Computing...");
                    ComputeSiteLikelihoodsForABlock (0,siteResults->theData, siteScalerBuffer);
                    PostCompute();
                } else
                    // dependant condition failed
                {
                    siteResults->PopulateConstantMatrix (1e-100);
                    //printf ("[PWNED] %d\n", hy_mpi_node_rank);
                }

                for (long k = 0; k < siteScalerBuffer.lLength; k++) {
                    siteResults->theData[siteScalerBuffer.lLength+k] = siteScalerBuffer.lData[k];
                    //printf ("%d %d => %g %g\n", hy_mpi_node_rank, k, siteResults->theData[k], siteResults->theData[siteScalerBuffer.lLength+k]);
                }
            }

        }
        //else
        //  printf ("%d : NOTHING TO DO!\n", hy_mpi_node_rank);

        //printf ("%d [mode = %d] %d/%d\n", hy_mpi_node_rank, partMode, siteResults->GetSize(), siteScalerBuffer.lLength);

        if (partMode) {
            ReportMPIError(MPI_Send(&siteLL, 1, MPI_DOUBLE, senderID, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD),true);
        } else
            // need to send both
        {
            ReportMPIError(MPI_Send(siteResults->theData, siteResults->GetSize(), MPI_DOUBLE, senderID, HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD),true);
        }
        ReportMPIError(MPI_Recv(variableStash.theData, indexInd.lLength, MPI_DOUBLE, senderID, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD,&status),false);
    }
    //ReportWarning (_String("Exiting slave loop after step  ") & loopie);
#endif
}



//_______________________________________________________________________________________
_Matrix*    _LikelihoodFunction::ConstructCategoryMatrix (const _SimpleList& whichParts, unsigned runMode, bool remap, _String* storageID) {
    long                        hDim = 1,
                                vDim = 0,
                                currentOffset = 0;


    PrepareToCompute();
    if (runMode == _hyphyLFConstructCategoryMatrixConditionals || runMode == _hyphyLFConstructCategoryMatrixWeights)
        // just return the matrix with class weights
    {
        for (long whichPart=0; whichPart<whichParts.lLength; whichPart++) {
            long myCatCount = TotalRateClassesForAPartition (whichParts.lData[whichPart]);
            if (myCatCount > hDim) {
                hDim = myCatCount;
            }
        }
    }


    if (runMode == _hyphyLFConstructCategoryMatrixWeights) {
        _Matrix     * catWeights = new _Matrix (whichParts.lLength, hDim, false, true);
        _SimpleList scalers;
        for (long whichPart=0; whichPart<whichParts.lLength; whichPart++) {
            PopulateConditionalProbabilities (whichParts.lData[whichPart],
                                              _hyphyLFConditionProbsClassWeights,
                                              catWeights->theData + hDim*whichPart,
                                              scalers);

        }
        catWeights->Transpose();
        DoneComputing();
        return catWeights;
    }

    // compute the number of columns in the matrix

    if (templateKind < 0) {
        vDim    =   GetIthFilter(whichParts.lData[0])->GetSiteCountInUnits();
    } else
        for (long i=0; i<whichParts.lLength; i++)
            if (runMode != _hyphyLFConstructCategoryMatrixConditionals
                    && HasHiddenMarkov(blockDependancies.lData[whichParts.lData[i]])>=0) {
                vDim    +=  GetIthFilter(whichParts.lData[i])->GetSiteCountInUnits();
            }
    // all sites
            else {
                vDim    +=      BlockLength(i);
            }
    // only unique patterns for now



    if (runMode == _hyphyLFConstructCategoryMatrixClasses || runMode == _hyphyLFConstructCategoryMatrixSiteProbabilities)
        // just return the maximum conditional likelihood category
    {
        _Matrix      *result = new _Matrix (hDim,vDim,false,true),
                      *cache     = nil;
        _SimpleList  *scalerCache = nil;

        bool         done         = false;

        if (runMode == _hyphyLFConstructCategoryMatrixSiteProbabilities) {
            long bufferL = PartitionLengths (0, &whichParts);
            cache        = new _Matrix (bufferL,2,false,true);
            scalerCache  = new _SimpleList (bufferL,0,0);
        } else {
            if (templateKind < 0) { // HMM
                _CategoryVariable*hmmVar = (_CategoryVariable*)FetchVar (-templateKind-1);

                Compute                              ();

                RunViterbi                           (*result,bySiteResults->theData,
                                                      *hmmVar->ComputeHiddenMarkov(),
                                                      *hmmVar->ComputeHiddenMarkovFreqs(),
                                                      nil,
                                                      &partScalingCache,
                                                      vDim
                                                     );

                remap = false;
                done = true;

            }

        }

        // now proceed to compute each of the blocks
        if (!done)
            for (long whichPart=0; whichPart<whichParts.lLength; whichPart++) {
                long i = whichParts.lData[whichPart];

                if (runMode == _hyphyLFConstructCategoryMatrixSiteProbabilities) {
                    long partititonSpan             = BlockLength (i);
                    ComputeSiteLikelihoodsForABlock (i, cache->theData, *scalerCache);
                    for (long c = 0; c < partititonSpan; c++) {
                        result->theData[currentOffset+c] = log(cache->theData[c]);
                        if (scalerCache->lData[c]) {
                            result->theData[currentOffset+c] -= scalerCache->lData[c]*_logLFScaler;
                        }


                    }
                    currentOffset                  += partititonSpan;
                    continue;
                }

                if (blockDependancies.lData[i] > 0)
                    // if a partition that does not depend on category variables
                    // then the matrix is already populated with zeros
                {
                    _SimpleList      *catVarType             = (_SimpleList*)((*(_List*)categoryTraversalTemplate(i))(4));
                    _SimpleList   const * duplicate_site_map = (_SimpleList   const*)GetIthFilter (i) -> GetDuplicateSiteMap();

                    long             categoryType        = catVarType->Element (-1),
                                     blockSize          = BlockLength (i),
                                     siteCount          = duplicate_site_map->lLength;

                    // check to see if we need to handle HMM or COP variables
                    if (categoryType & _hyphyCategoryHMM) {
                        _CategoryVariable*hmmVar        = (_CategoryVariable*)((*(_List*)(*(_List*)categoryTraversalTemplate(i))(0))(0));

                        /*
                           run the Viterbi algorithm to reconstruct the most likely
                           HMM path
                        */
                        ComputeSiteLikelihoodsForABlock    (i, siteResults->theData, siteScalerBuffer);

                        RunViterbi (*result,
                                    siteResults->theData,
                                    *hmmVar->ComputeHiddenMarkov(),
                                    *hmmVar->ComputeHiddenMarkovFreqs(),
                                    duplicate_site_map,
                                    &siteScalerBuffer,
                                    blockSize);


                        currentOffset += siteCount;
                    } else {

                        long hasConstantOnPartition = HasHiddenMarkov(blockDependancies.lData[i],false);
                        if (hasConstantOnPartition<0) {
                            _SimpleList scalers (blockSize, 0, 0);
                            // allocate a vector of 3*blockSize

                            hyFloat  * likelihoodBuffer = new hyFloat[3*blockSize];
                            PopulateConditionalProbabilities (i,_hyphyLFConditionProbsMaxProbClass,likelihoodBuffer,scalers);

                            for (long i = 0; i < blockSize; i++) {
                                result->theData[currentOffset+i] = likelihoodBuffer[i];
                            }

                            delete       [] likelihoodBuffer;
                        } else {
                            HandleApplicationError ("This feature has not yet been implemented in the new LF engine framework");
                            return result;

                            long    mxDim = 1,
                                    bl = BlockLength (i),
                                    j;

                            for (j=0; j<=hasConstantOnPartition; j++)
                                if (CheckNthBit (blockDependancies.lData[i],j)) {
                                    _CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[j]);
                                    mxDim *= thisC->GetNumberOfIntervals();
                                }

                            _Matrix ccache (mxDim,1,false,true);

                            categID = 0;
                            RecurseConstantOnPartition    (i,0,blockDependancies.lData[j],hasConstantOnPartition,1,ccache);
                            categID = 0;

                            /*
                                FindMaxCategory(i,  HighestBit( blockDependancies.lData[i]), blockDependancies.lData[i], ff+1, currentOffset, result);
                            */

                            hyFloat maxValue = -1.e300;
                            long       mxIndex  = -1;

                            for (j=0; j<mxDim; j++)
                                if (ccache.theData[j] > maxValue) {
                                    maxValue = ccache.theData[j];
                                    mxIndex = j;
                                }

                            mxDim = 1;
                            for (j=HighestBit (blockDependancies.lData[i]); j>hasConstantOnPartition; j--)
                                if (CheckNthBit (blockDependancies.lData[i],j)) {
                                    _CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[j]);
                                    mxDim *= thisC->GetNumberOfIntervals();
                                }

                            mxIndex *= mxDim;

                            for   (long kk = 0; kk<bl; kk++) {
                                result->theData[currentOffset+kk] += mxIndex;
                            }

                        }
                        currentOffset+=blockSize;
                    }
                }
            }
        DoneComputing();
        DeleteObject (cache);
        DeleteObject (scalerCache);
        if (remap) {
            _Matrix * retObj = RemapMatrix(result, whichParts);
            DeleteObject (result);
            return retObj;
        }
        return result;
    } else {
        long maxPartSize = 0;
        for (long whichPart=0; whichPart<whichParts.lLength; whichPart++) {
            long myWidth        = BlockLength(whichParts.lData[whichPart]);
            maxPartSize         = MAX (maxPartSize,myWidth);
        }

        _Vector  allScalers (false);
        _SimpleList     scalers;
        // allocate a buffer big enough to store the matrix for each block

        _Matrix         *result = new _Matrix (hDim,vDim,false,true),
        *holder = new _Matrix (hDim, maxPartSize, false, true);

        maxPartSize = 0; // reused as a shift index

        for (long whichPart=0; whichPart<whichParts.lLength; whichPart++) {
            PopulateConditionalProbabilities (whichParts.lData[whichPart],
                                              _hyphyLFConditionProbsScaledMatrixMode,
                                              holder->theData,
                                              scalers);

            allScalers << scalers;
            long        thisBlockSiteCount = BlockLength(whichParts.lData[whichPart]);
            result->CopyABlock (holder, 0, maxPartSize, ((_SimpleList*)(*(_List*)categoryTraversalTemplate(whichParts.lData[whichPart]))(1))->Element(-1), thisBlockSiteCount);
            maxPartSize += thisBlockSiteCount;
        }
        DoneComputing   ();
        DeleteObject    (holder);

        if (remap) {
            if (storageID) {
                _Matrix * remappedCorrections = RemapMatrix(&allScalers, whichParts);
                _String   scalerID            = (*storageID) & categoryMatrixScalers;
                CheckReceptacleAndStore (&scalerID, kEmptyString, false, remappedCorrections, false);
                scalerID = (*storageID) & categoryLogMultiplier;
                CheckReceptacleAndStore (&scalerID, kEmptyString, false, new _Constant(_logLFScaler), false);
            }

            holder          = RemapMatrix(result, whichParts);
            DeleteObject    (result);
            return          holder;
        } else {
            return result;
        }
    }
    return nil;
}

//_______________________________________________________________________________________
long    _LikelihoodFunction::PartitionLengths       (char runMode,  _SimpleList const * filter)
{
    long maxDim = 0;

    for (long i=0; i<(filter?filter->lLength:theTrees.lLength); i++) {
        long bl = BlockLength (filter?filter->lData[i]:i);
        if (runMode == 0) {
            maxDim = MAX (maxDim, bl);
        } else {
            maxDim += bl;
        }
    }

    return maxDim;
}


//_______________________________________________________________________________________
void    _LikelihoodFunction::AllocateSiteResults        (void)
{
    long dim            = PartitionLengths(0),
         catSpan       = TotalRateClassesForAPartition(-1,1) + 1;

    siteResults         = new _Matrix (dim,catSpan,false,true);
    siteScalerBuffer.Populate (dim,0,0);
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::ZeroSiteResults        (void)
{
    if (siteResults) {
        long upperLimit = siteResults->GetSize();
        for (long k=0; k<upperLimit; k++) {
            siteResults->theData[k] = 0;
        }
        siteScalerBuffer.Populate (upperLimit,0,0);
    }
}


//_______________________________________________________________________________________

#ifndef __HYPHYMPI__
bool      _LikelihoodFunction::SendOffToMPI       (long)
{
    return false;
#else
bool      _LikelihoodFunction::SendOffToMPI       (long index) {
// dispatch an MPI task to node 'index+1'

/* 20170404 SLKP Need to check if the decision to recompute a partition is made correctly.
 In particular, need to confirm that changes to category variables are handled correctly (e.g. HaveParametersChanged, vs has changed */

    bool                sendToSlave = (computationalResults.GetSize() < parallelOptimizerTasks.lLength);
    _SimpleList     *   slaveParams = (_SimpleList*)parallelOptimizerTasks(index);

    for (unsigned long varID = 0UL; varID < slaveParams->lLength; varID++) {
        _Variable * aVar = LocateVar (slaveParams->lData[varID]);
        if (aVar->IsIndependent()) {
            varTransferMatrix.theData[varID] = aVar->Value();
        } else {
            varTransferMatrix.theData[varID] = aVar->Compute()->Value();
        }

        //printf ("%s => %g\n", aVar->GetName()->sData, varTransferMatrix.theData[varID]);
        sendToSlave = sendToSlave || aVar->HasChanged();
    }

    if (sendToSlave) {
        ReportMPIError(MPI_Send(varTransferMatrix.theData, slaveParams->lLength, MPI_DOUBLE, index+1 , HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD),true);
    }
    return sendToSlave;


#endif //__HYPHYMPI__

}


//_______________________________________________________________________________________

bool    _LikelihoodFunction::PreCompute         (void)
{

    useGlobalUpdateFlag = true;
    // mod 20060125 to only update large globals once
    unsigned long i = 0UL;

    _SimpleList * arrayToCheck = nonConstantDep?nonConstantDep:&indexDep;

    for (; i < arrayToCheck->lLength; i++) {
        _Variable* cornholio = LocateVar(arrayToCheck->lData[i]);
        hyFloat tp = cornholio->Compute()->Value();
        if (!cornholio->IsValueInBounds(tp)){
            ReportWarning (_String ("Failing bound checks on ") & *cornholio->GetName() & " = " & _String (tp, "%25.16g"));
            break;
        }
    }

    useGlobalUpdateFlag = false;
    // mod 20060125 to only update large globals once

    for (unsigned long j=0UL; j<arrayToCheck->lLength; j++) {
        _Variable* cornholio = LocateVar(arrayToCheck->lData[j]);
        if (cornholio->varFlags&HY_DEP_V_COMPUTED) {
            cornholio->varFlags -= HY_DEP_V_COMPUTED;
        }
    }

    return (i==arrayToCheck->lLength);
}


//_______________________________________________________________________________________

void    _LikelihoodFunction::PostCompute        (void) {
    _SimpleList * arrayToCheck = nonConstantDep?nonConstantDep:&indexDep;

    //useGlobalUpdateFlag = true;
    for (unsigned long i=0; i<arrayToCheck->lLength; i++) {
        LocateVar (arrayToCheck->lData[i])->Compute();
    }
    //useGlobalUpdateFlag = false;
    // mod 20060125 comment out the compute loop; seems redundant
    for (unsigned long i=0UL; i<indexInd.lLength; i++) {
        _Variable * this_p = GetIthIndependentVar(i);
        if (this_p->varFlags & HY_VARIABLE_CHANGED) {
          this_p->varFlags -= HY_VARIABLE_CHANGED;
        }
    }
}



//_______________________________________________________________________________________

void    _LikelihoodFunction::ComputeBlockForTemplate        (long i, bool force)
{
    long          blockWidth = bySiteResults->GetVDim();
    hyFloat   * resStore  = bySiteResults->theData+i*blockWidth;

    ComputeSiteLikelihoodsForABlock (i,bySiteResults->theData+i*blockWidth,*(_SimpleList*)partScalingCache(i));
    if (! usedCachedResults) {
        long    *   siteCorrectors  = ((_SimpleList**)siteCorrections.lData)[i]->lData,
                    upto           = ((_SimpleList**)siteCorrections.lData)[i]->lLength;
        for (long s = 0; s < upto; s++) {
            resStore[s] *= acquireScalerMultiplier(siteCorrectors[s]);
        }
    }

    if (force || !usedCachedResults)
        // remap compressed sites to full length
    {
        ComputeBlockForTemplate2 (i, resStore, resStore, blockWidth);
    }
}


//_______________________________________________________________________________________

void    _LikelihoodFunction::ComputeBlockForTemplate2       (long i, hyFloat* resTo, hyFloat* resFrom, long blockWidth)
{
    _DataSetFilter    const   *df = GetIthFilter(i);
    long*               dupMap = df->duplicateMap.lData,
                        dupL   = df->duplicateMap.lLength;

    if (resTo == resFrom) {
        _Matrix             temp (1,blockWidth,false,true);
        for (long k1 = 0; k1 < dupL; k1++) {
            temp.theData[k1] = resFrom[dupMap[k1]];
        }

        for (long k2 = 0; k2< dupL; k2++) {
            resTo[k2] = temp.theData[k2];
        }

        for (long k3 = dupL; k3 < blockWidth; k3++) {
            resTo[k3] = 1.;
        }
    } else {
        for (long k1 = 0; k1 < dupL; k1++) {
            resTo[k1] = resFrom[dupMap[k1]];
        }

        for (long k3 = dupL; k3 < blockWidth; k3++) {
            resTo[k3] = 1.;
        }
    }
}

//_______________________________________________________________________________________

hyFloat  _LikelihoodFunction::Compute        (void)
/*
    code cleanup SLKP: 20090317

    todo: need to ensure that partitions that are changed between calls
          are properly identified

          in prepare to compute calls, perhaps set up a map of
          independent variable ID -> partition/class pair
*/
{

    hyFloat result = 0.;

    if (!PreCompute()) {
        return -INFINITY;
    }

    /* GUI flag to verify whether MLEs have been altered
       after last optimization
     */

    if (!isInOptimize && hasBeenOptimized)
        for (unsigned long i=0; i<indexInd.lLength; i++)
            if (LocateVar (indexInd.lData[i])->HasChanged()) {
                hasBeenOptimized = false;
                break;
            }

    /*  compute modes:

            0. need only log-likelihood values for each partition; the likelihood function itself
              is computed by summing individual partition log-likelihoods

            1. need likelihoods by pattern; for each partition a vector of likelihoods and
              a vector of scaling parameters is returned

            2. need conditional likelihoods by pattern; for each partition a matrix of conditional
              likelihoods (site/rate class pairing) is returned along with a vector of scaling
              factors

            3. need log-likelihoods by block: for each partition, the log-likelihood value is
              returned

            4. MPI compute mode: independent partitions
    */

    char       computeMode = 0;
    if        (computingTemplate) {
        if (templateKind > _hyphyLFComputationalTemplateByPartition) {
            computeMode = 2;
        } else if (templateKind == _hyphyLFComputationalTemplateByPartition) {
            computeMode = 3;
        } else {
            computeMode = 1;
        }

        /*
         write out matrices of conditional likelihoods for each partition
         (into SITE_LIKELIHOODS_0, SITE_LIKELIHOODS_1 etc)
         as well as site scaling factors
         (into SITE_LIKELIHOODS_0.scalers, SITE_LIKELIHOODS_1.scalers ...)
         */

    }

#ifdef __HYPHYMPI__
    if ((hyphyMPIOptimizerMode ==  _hyphyLFMPIModePartitions || hyphyMPIOptimizerMode ==  _hyphyLFMPIModeAuto) &&   hy_mpi_node_rank == 0) {
        computeMode = 4;
    }
#endif

    bool done = false;
#ifdef _UBER_VERBOSE_LF_DEBUG
    fprintf (stderr, "\n*** Likelihood function evaluation %ld ***\n", likeFuncEvalCallCount+1);
    for (unsigned long i=0; i<indexInd.lLength; i++) {
        _Variable *v = LocateVar (indexInd.lData[i]);
        if (v->HasChanged()) {
          fprintf (stderr, "[CHANGED] ");
        }
        fprintf (stderr, "%s = %15.12g\n", v->GetName()->sData, v->theValue);
    }
#endif
    if (computeMode == 0 || computeMode == 3) {
        _Matrix     * blockMatrix = nil;
        if (computeMode == 3) {
            blockWiseVar->SetValue      (new _Matrix (theTrees.lLength,1,false,true), false);
            blockMatrix = (_Matrix*)blockWiseVar->GetValue();

        }
        for (unsigned long partID=0; partID<theTrees.lLength; partID++) {
            if (blockDependancies.lData[partID])
                // has category variables
            {
                if ( computationalResults.get_used()<=partID || HasBlockChanged(partID))
                    // first time computing or partition requires updating
                {
                    /* TODO: add HMM and constant on partition test
                       Roll into ComputeSiteLikelihoodsForABlock and SumUpSiteLikelihoods?
                    */

#ifdef __HYPHYMPI__
                    if (hy_mpi_node_rank == 0) {
                        ComputeSiteLikelihoodsForABlock    (partID, siteResults->theData, siteScalerBuffer, -1, nil, hyphyMPIOptimizerMode);
                    } else
#endif
                        ComputeSiteLikelihoodsForABlock    (partID, siteResults->theData, siteScalerBuffer);

#ifdef _UBER_VERBOSE_LF_DEBUG
                    fprintf (stderr, "Did compute %g\n", result);
#endif
                    hyFloat                       blockResult = SumUpSiteLikelihoods (partID, siteResults->theData, siteScalerBuffer);
                    UpdateBlockResult               (partID, blockResult);
                    if (blockMatrix) {
                        blockMatrix->theData[partID] = blockResult;
                    } else {
                        result += blockResult;
                    }
                } else {
                    if (blockMatrix) {
                        blockMatrix->theData[partID] =  computationalResults.theData[partID];
                    }  else {
                        result += computationalResults.theData[partID];
                    }
                }
            } else {
                hyFloat  blockResult =  ComputeBlock (partID);
                if (blockMatrix) {
                    blockMatrix->theData[partID] = blockResult;
                } else {
                    result += blockResult;
                }
                UpdateBlockResult       (partID, blockResult);
            }

        }
        if (blockMatrix) {
            result = computingTemplate->Compute()->Value();
        }
        done = true;
    } else if (computeMode == 1)
        // handle _hyphyLFComputationalTemplateBySite
    {
        long    blockWidth = bySiteResults->GetVDim();

#ifdef __HYPHYMPI__
        if (hyphyMPIOptimizerMode == _hyphyLFMPIModeSiteTemplate && hy_mpi_node_rank == 0) {
            long    totalSent = 0,
                    blockPatternSize = PartitionLengths (0);
            for (long blockID = 0; blockID < parallelOptimizerTasks.lLength; blockID ++) {
                bool sendToSlave = SendOffToMPI (blockID);
                if (sendToSlave) {
                    //printf ("[SEND]: %d\n", blockID, "\n");
                    totalSent++;
                }
                /*
                 else
                    for (long k = 0; k < blockPatternSize; k++)
                        printf ("[%d] %d %g %d\n", blockID, k, (bySiteResults->theData + blockID*theTrees.lLength)[k],  ((_SimpleList*)partScalingCache(blockID))->lData[k]);
                */
            }

            while (totalSent) {
                MPI_Status      status;

                ReportMPIError  (MPI_Recv (bySiteResults->theData + blockWidth*theTrees.lLength, 2*blockPatternSize, MPI_DOUBLE, MPI_ANY_SOURCE , HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD,&status),true);
                long partID = status.MPI_SOURCE-1;

                //printf ("[RECEIVE]: %d\n", partID, "\n");

                _DataSetFilter const * thisBlockFilter      = GetIthFilter(partID);
                thisBlockFilter->PatternToSiteMapper (bySiteResults->theData + blockWidth*theTrees.lLength, bySiteResults->theData + partID*blockWidth, blockWidth, 1.);
                _SimpleList* blockScalers = ((_SimpleList*)partScalingCache(partID));
                thisBlockFilter->PatternToSiteMapper (bySiteResults->theData + (blockWidth*theTrees.lLength+blockPatternSize), blockScalers->lData, blockWidth, 0L);
                totalSent--;
            }
        } else
#endif
            for (long partID=0; partID<theTrees.lLength; partID++) {
                ComputeSiteLikelihoodsForABlock      (partID, bySiteResults->theData + blockWidth*theTrees.lLength, *(_SimpleList*)partScalingCache(theTrees.lLength));
                if (usedCachedResults == false) {
                    _DataSetFilter const * thisBlockFilter      = GetIthFilter(partID);
                    thisBlockFilter->PatternToSiteMapper (bySiteResults->theData + blockWidth*theTrees.lLength, bySiteResults->theData + partID*blockWidth, blockWidth, 1.);
                    thisBlockFilter->PatternToSiteMapper (((_SimpleList*)partScalingCache(theTrees.lLength))->lData, ((_SimpleList*)partScalingCache(partID))->lData, blockWidth, 0L);
                }
            }


        if (templateKind < 0) { // HMM
            _CategoryVariable*hmmVar = (_CategoryVariable*)FetchVar (-templateKind-1);
            _Matrix          *hmm    = hmmVar->ComputeHiddenMarkov(),
                              *hmf    = hmmVar->ComputeHiddenMarkovFreqs();

            result           = SumUpHiddenMarkov (bySiteResults->theData,
                                                  *hmm,
                                                  *hmf,
                                                  nil,
                                                  &partScalingCache,
                                                  blockWidth
                                                 );

        } else {
            // no need to remap; just process directly based on partition indices

            siteArrayPopulated         = true;
            siteWiseVar->SetValue      (new _Matrix (theTrees.lLength,1,false,true), false);
            _SimpleList scalerCache    (theTrees.lLength,0,0);
            _Matrix     * siteMatrix = (_Matrix*)siteWiseVar->GetValue();

            for (long siteID = 0; siteID < blockWidth; siteID++) {
                // pass one to determine scaling factors
                long minScalingFactor = ((_SimpleList*)partScalingCache(0))->lData[siteID];
                scalerCache.lData [0] = minScalingFactor;
                for (long partID=1; partID<theTrees.lLength; partID++) {
                    scalerCache.lData [partID] = ((_SimpleList*)partScalingCache(partID))->lData[siteID];
                    if (minScalingFactor > scalerCache.lData [partID]) {
                        minScalingFactor = scalerCache.lData [partID];
                    }
                }
                for (long partID=0; partID<theTrees.lLength; partID++) {
                    siteMatrix->theData[partID] = bySiteResults->theData[partID*blockWidth+siteID];
                    long diff = scalerCache.lData[partID]-minScalingFactor;
                    if (diff) {
                        siteMatrix->theData[partID] *= acquireScalerMultiplier(diff);
                    }
                }



                result += computingTemplate->Compute()->Value();
                if (minScalingFactor) {
                    result-=_logLFScaler*minScalingFactor;
                }
            }
        }
        done = true;
    } else {
        if (computeMode == 4) {
#ifdef __HYPHYMPI__
            if (hy_mpi_node_rank == 0) {
                long    totalSent = 0;
                for (long blockID = 0; blockID < parallelOptimizerTasks.lLength; blockID ++) {
                    bool sendToSlave = SendOffToMPI (blockID);
                    if (sendToSlave) {
                        totalSent++;
                    } else {
                        result += computationalResults[blockID];
                    }
                }

                while (totalSent) {
                    MPI_Status      status;
                    hyFloat      blockRes;
                    ReportMPIError(MPI_Recv (&blockRes, 1, MPI_DOUBLE, MPI_ANY_SOURCE , HYPHY_MPI_DATA_TAG, MPI_COMM_WORLD,&status),true);
                    //printf ("Got %g from block %d \n", blockRes, status.MPI_SOURCE-1);

                    result            += blockRes;
                    /*if (status.MPI_SOURCE == 1 && computationalResults.GetUsed()) {
                        printf ("\033\015 COMPUTED / CACHED %g / %g             ", blockRes, computationalResults[0]);
                    }*/
                    UpdateBlockResult (status.MPI_SOURCE-1, blockRes);
                    totalSent--;
                }

                done = true;
            }
#endif
        }
    }

    if (done) {

        likeFuncEvalCallCount ++;
        evalsSinceLastSetup   ++;
        PostCompute ();
#ifdef _UBER_VERBOSE_LF_DEBUG
        fprintf (stderr, "%g\n", result);
#endif
        if (isnan (result)) {
            ReportWarning ("Likelihood function evaluation encountered a NaN (probably due to a parameterization error or a bug).");
            return -INFINITY;
        }
        
        if (result >= 0.) {
            if (result >= __DBL_EPSILON__ * 1.e4) {
                char buffer [2048];
                snprintf (buffer, 2047, "Internal error: Encountered a positive log-likelihood (%g) at evaluation %ld, mode %ld, template %ld. This is usually a consequence of 'infinite-like' parameter values", likeFuncEvalCallCount-1, result, computeMode, templateKind);
                _TerminateAndDump(buffer, true);
            } else {
                result = 0.;
            }
        }

        ComputeParameterPenalty ();

        hyFloat regularized_value = result - smoothingPenalty;
        return regularized_value;
    }

    HandleApplicationError ("Sorry; this likelihood feature has not yet been ported to the v2.0 LF engine in HyPhy");
    return -INFINITY;


}
//_______________________________________________________________________________________

long        _LikelihoodFunction::BlockLength(long index) const {
    return GetIthFilter (index)->GetPatternCount();
}
//_______________________________________________________________________________________

bool        _LikelihoodFunction::HasBlockChanged(long index) const {
    return GetIthTree (index)->HasChanged2();
}

//_______________________________________________________________________________________

void      _LikelihoodFunction::RecurseConstantOnPartition (long blockIndex, long index, long dependance, long highestIndex, hyFloat weight, _Matrix& cache)
{
    _CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[index]);

    if (index<highestIndex) {
        if ((!CheckNthBit(dependance,index))||thisC->IsHiddenMarkov()) {
            RecurseCategory (blockIndex, index+1, dependance,highestIndex,weight);
        } else {
            thisC->Refresh();
            long nI = thisC->GetNumberOfIntervals ();
            offsetCounter *= nI;
            for (long k = 0; k<nI; k++) {
                thisC->SetIntervalValue(k);
                RecurseConstantOnPartition(blockIndex,index+1,dependance, highestIndex,weight*thisC->GetIntervalWeight(k),cache);
                categID+=offsetCounter/nI;
            }
            offsetCounter/=nI;
            if (offsetCounter>1) {
                categID-=nI*offsetCounter;
            }
        }
    } else {
        long        category_count         = thisC->GetNumberOfIntervals (),
                    current_index_offset   = BlockLength(blockIndex),
                    high_bit               = HighestBit(blockDependancies.lData[blockIndex]);

        thisC->Refresh();

        hyFloat* site_results     = siteResults->fastIndex();
        _Matrix*    category_weights = thisC->GetWeights();

        _DataSetFilter const * data_filter = GetIthFilter(blockIndex);

        for (unsigned long category_index = 0UL;  category_index < category_count; category_index ++) {
            thisC->SetIntervalValue(category_index,!category_index);

            InitializeArray (site_results, current_index_offset, 0.0);

            //sR -= currentOffset;

            if (high_bit > index) {
                offsetCounter *= category_count;
                RecurseCategory (blockIndex,index+1,blockDependancies.lData[blockIndex],high_bit,1);
                offsetCounter /= category_count;
            } else {
                ComputeBlock (blockIndex,site_results);
            }

            hyFloat log_sum  = 0.0;
            for   (long kk = 0; kk<current_index_offset; kk++,site_results++) {
              log_sum +=  myLog (*site_results)*data_filter->GetFrequency(kk);
            }

            log_sum += myLog (category_weights->theData[category_index]*weight);
            cache.theData[categID] = log_sum;

            categID+=offsetCounter;
            site_results -= current_index_offset;
        }

        if (offsetCounter>1) {
            categID -= category_count*offsetCounter;
        }
    }
}

//_______________________________________________________________________________________

void      _LikelihoodFunction::RecurseCategory(long blockIndex, long index, long dependance, long highestIndex, hyFloat weight
#ifdef _SLKP_LFENGINE_REWRITE_
        ,_SimpleList* siteMultipliers, char runMode, hyFloat *runStorage,
        long branchIndex,              _SimpleList* branchValues
#endif
                                              )
{
    _CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[index]);
    if (index<highestIndex) {
        if ((!CheckNthBit(dependance,index))||thisC->IsHiddenMarkov())
            RecurseCategory (blockIndex, index+1, dependance,highestIndex,weight
#ifdef _SLKP_LFENGINE_REWRITE_
                             ,siteMultipliers,runMode,runStorage
#endif
                            );
        else {
            thisC->Refresh();
            long nI = thisC->GetNumberOfIntervals ();
            offsetCounter *= nI;
            for (long k = 0; k<nI; k++) {
                thisC->SetIntervalValue(k);
                RecurseCategory(blockIndex,index+1,dependance, highestIndex,weight*thisC->GetIntervalWeight(k)
#ifdef _SLKP_LFENGINE_REWRITE_
                                ,siteMultipliers,runMode,runStorage,branchIndex,branchValues
#endif
                               );
                categID+=offsetCounter/nI;
            }
            offsetCounter/=nI;
            if (offsetCounter>1) {
                categID-=nI*offsetCounter;
            }
        }
    } else {
        if (thisC->IsHiddenMarkov()) {
            if (offsetCounter == 1) { // this is the only categ for the block
                ComputeBlock (blockIndex,siteResults->fastIndex());
            }
        } else {
            long        hDim            = siteResults->GetHDim(),
                        nI              = thisC->GetNumberOfIntervals (),
                        currentOffset  = BlockLength(blockIndex);

            thisC->Refresh();

            hyFloat* sR  = siteResults->fastIndex();
            _Matrix*    cws = thisC->GetWeights();

#ifdef _SLKP_LFENGINE_REWRITE_
            long    *   siteCorrectors  = ((_SimpleList**)siteCorrections.lData)[blockIndex]->lLength?
                                          (((_SimpleList**)siteCorrections.lData)[blockIndex]->lData) + categID * currentOffset
                                          :nil;
#endif

            for (long k = 0; k<nI; k++) {
                thisC->SetIntervalValue     (k,!k);
                ComputeBlock                (blockIndex,sR+hDim);
                hyFloat                  localWeight = cws->theData[k]*weight ;
#ifdef _SLKP_LFENGINE_REWRITE_
                if (runMode == 1)
                    // decide what is the most likely category assignment
                {
                    for (long r1 = 0, r2 = hDim; r1 < currentOffset; r1++,r2++) {
                        bool doChange = false;
                        if (siteCorrectors) {
                            long scv = *siteCorrectors;

                            if (scv < siteMultipliers->lData[r1]) { // this has a _smaller_ scaling factor
                                hyFloat scaled = sR[r1]*acquireScalerMultiplier (siteMultipliers->lData[r1] - scv);
                                if (sR[r2] > scaled) {
                                    doChange = true;
                                } else {
                                    sR[r1] = scaled;
                                }
                                siteMultipliers->lData[r1] = scv;
                            } else {
                                if (scv > siteMultipliers->lData[r1]) { // this is a _larger_ scaling factor
                                    sR[r2] *= acquireScalerMultiplier (scv - siteMultipliers->lData[r1]);
                                }
                                doChange = sR[r2] > sR[r1] && ! CheckEqual (sR[r2],sR[r1]);
                            }

                            siteCorrectors++;
                        } else {
                            doChange = sR[r2] > sR[r1] && ! CheckEqual (sR[r2],sR[r1]);
                        }

                        if (doChange) {
                            sR[r1]         = sR[r2];
                            runStorage[r1] = categID;
                        }
                    }
                } else {
                    /*if (runMode == 2)
                    // write out conditional probs
                    {
                        for (long r1 = 0, r2 = hDim, r3 = categID*(hDim+columnShift); r1 < currentOffset; r1++,r2++,r3++)
                        {
                            if (siteCorrectors)
                            {
                                long scv = *siteCorrectors;
                                if (scv < siteMultipliers->lData[r1]) // this has a _smaller_ scaling factor
                                {
                                    hyFloat scaled = acquireScalerMultiplier (siteMultipliers->lData[r1] - scv);
                                    for (long z = r3-(hDim+columnShift); z >= 0; z-=(hDim+columnShift))
                                        runStorage[z] *= scaled;
                                    siteMultipliers->lData[r1] = scv;
                                }
                                else
                                {
                                    if (scv > siteMultipliers->lData[r1]) // this is a _larger_ scaling factor
                                        sR[r2] *= acquireScalerMultiplier (scv - siteMultipliers->lData[r1]);
                                }

                                siteCorrectors++;
                            }
                            runStorage[r3] = sR[r2];
                        }
                    }
                    else*/
#endif
                    for (long r1 = 0, r2 = hDim; r1 < currentOffset; r1++,r2++) {
#ifdef _SLKP_LFENGINE_REWRITE_
                        if (siteCorrectors) {
                            long scv = *siteCorrectors;

                            //if (likeFuncEvalCallCount == 43)
                            //{
                            //  printf ("catID %d site %d corrector %d (%d) value %g:%g\n", categID, r1, scv, siteMultipliers->lData[r1], sR[r2], sR[r2]*acquireScalerMultiplier (scv-siteMultipliers->lData[r1] ));
                            //}

                            if (scv < siteMultipliers->lData[r1]) { // this has a _smaller_ scaling factor
                                sR[r1] = localWeight * sR[r2] + sR[r1] * acquireScalerMultiplier (siteMultipliers->lData[r1] - scv);
                                siteMultipliers->lData[r1] = scv;
                            } else {
                                if (scv > siteMultipliers->lData[r1]) { // this is a _larger_ scaling factor
                                    sR[r1] += localWeight * sR[r2] * acquireScalerMultiplier (scv - siteMultipliers->lData[r1]);
                                } else { // same scaling factors
                                    sR[r1] += localWeight * sR[r2];
                                }
                            }

                            siteCorrectors++;
                        } else
#endif
                            sR[r1] += localWeight * sR[r2];

                    }
#ifdef _SLKP_LFENGINE_REWRITE_
                }
#endif
                categID+=offsetCounter;
#ifdef _SLKP_LFENGINE_REWRITE_
                if (offsetCounter > 1) {
                    siteCorrectors += (offsetCounter-1) * currentOffset;
                }
#endif
            }
            if (offsetCounter>1) {
                categID-=nI*offsetCounter;
            }
        }
    }
}

//_______________________________________________________________________________________
void        _LikelihoodFunction::RandomizeList (_SimpleList& orderList, long elements) {
    orderList.Clear();
    orderList.Populate (elements, 0, 1);
    orderList.Permute(0);
}

//_______________________________________________________________________________________
void        _LikelihoodFunction::CheckFibonacci (hyFloat shrinkFactor)
{
    long n=Fibonacci.lLength-1;
    if (n<0) {
        Fibonacci<<1;
        Fibonacci<<1;
        n+=2;
    }
    while (Fibonacci(n)<shrinkFactor) {
        Fibonacci<<Fibonacci(n)+Fibonacci(n-1);
        n++;
    }
}
//_______________________________________________________________________________________

void    _LikelihoodFunction::CheckDependentBounds (void) {
/*
   this function makes sure that a constrained optimization starts within the domain
   of allowed parameter values
*/

    unsigned long dep_var_count = indexDep.countitems();
    
    if  (dep_var_count == 0UL) { // nothing to do here
        return;
    }
    
    long index, i, j, badConstraint;
    
    _Matrix currentValues (dep_var_count,1,false,true),
    lowerBounds   (dep_var_count,1,false,true),
    upperBounds   (dep_var_count,1,false,true);
    
    bool    ohWell = false;
    
    _Variable*  cornholio;
    _SimpleList badIndices; //indices of dependent variables which are out of bounds
    
    if (nonConstantDep) {
        DeleteObject(nonConstantDep);
    }
    nonConstantDep = new _SimpleList;
    _SimpleList nonConstantIndices; // for error reporting
    
    for (index = 0; index< dep_var_count && !ohWell; index++) {
        // check whether any of the dependent variables are out of bounds
        cornholio                       =   GetIthDependentVar(index);
        
        currentValues.theData[index]    =   cornholio->Compute()->Value();
        lowerBounds.theData[index]      =   cornholio->GetLowerBound();
        upperBounds.theData[index]      =   cornholio->GetUpperBound();
        
        //fprintf (stderr, "_LikelihoodFunction::CheckDependentBounds variable %s (%d), current value %g, range %g to %g\n", cornholio->theName->sData, index, currentValues.theData[index], lowerBounds.theData[index], upperBounds.theData[index]);
        
        bool badApple = currentValues.theData[index]<lowerBounds.theData[index] || currentValues.theData[index]>upperBounds.theData[index];
        if (badApple) {
            badIndices<<index;
            //fprintf (stderr, "---> Constraint violated\n");
        }
        
        if (cornholio->IsConstant()) {
            badConstraint = indexDep.lData[index];
            ohWell = badApple;
            j      = index; // for error reporting at the bottom
        } else {
            (*nonConstantDep) << indexDep.lData[index];
            nonConstantIndices << index;
        }
    }
    
    if (badIndices.lLength && !ohWell) // one of the variables has left its prescribed bounds
                                       // build a table of dependancies
    {
        _Matrix     dependancies (MAX(3,indexDep.lLength),indexInd.lLength,true,true);
        
        // element (i,j) represents the dependance of i-th dep var on the j-th ind var
        // 0 - no dep,
        // 1 -> monotone increase,
        // -1 -> monotone decrease
        
        for (index = 0; index<indexInd.lLength; index++) {
            hyFloat          temp = GetIthIndependent(index);
            SetIthIndependent  (index,temp*1.000000000001);
            
            for (j=0; j < indexDep.lLength; j++) {
                hyFloat temp1 = GetIthDependent(j);
                if (temp1>currentValues[j]) {
                    dependancies.Store(j,index,1.0);
                } else if (temp1<currentValues[j]) {
                    dependancies.Store(j,index,-1.0);
                }
            }
            SetIthIndependent (index,temp);
        }
        
        //fprintf (stderr, "\n%s\n", _String((_String*)dependancies.toStr()).sData);
        
        // now we can go through the dependant variables which are out of bounds one at a time
        // and attempt to move them back in.
        for (index = badIndices.lLength-1; index>-1; index--) {
            long       badVarIndex      = badIndices.lData[index];
            hyFloat temp             = GetIthDependent(badVarIndex);
            
            currentValues[badVarIndex]  = temp;
            
            bool              tooLow = temp<lowerBounds[badVarIndex];
            
            // look for an independent which will increase (if tooLow) or decrease (too high) this variable
            // w/o decreasing the others
            // this is equivalent to searching for the matrix of dependancies for a column w/o negative
            // entries,such that row "index" has "1" in it
            
            hyFloat correlation = 0.0;
            
            for (j = 0; j < indexInd.lLength; j++) {
                if ((correlation = dependancies(badVarIndex,j)) != 0.) {
                    // fprintf (stderr, "## %d -> %g\n", j, correlation);
                    for (i=0; i<badIndices.lLength; i++) {
                        if (i != index) {
                            hyFloat depIJ = dependancies(badIndices.lData[i],j);
                            if (depIJ && depIJ != correlation) {
                                break;
                            }
                        }
                    }
                    if (i==badIndices.lLength) {
                        break;    //match found
                    }
                }
                
            }
            
            if (j == indexInd.lLength) { // no suitable variable found - try random permutations
                break;
            }
            
            //fprintf (stderr, "Will try to adjust indepdenent variable %s (%d)\n", GetIthIndependentName(j)->sData, j);
            
            // now nudge the independent variable (indexed by "j") upward (or downward),
            // until var badVarIndex is within limits again
            // try this by a trivial bisection
            
            hyFloat left         ,
            right      ,
            middle,
            lb         = tooLow?lowerBounds[badVarIndex]:upperBounds[badVarIndex];
            
            bool    decrement = (correlation < 0. && tooLow) || (correlation > 0. && !tooLow);
            
            if (decrement) { // need to decrement "j"
                left  = GetIthIndependentBound(j,true);
                right = GetIthIndependent(j);
            } else { // need to increment "j"
                right = GetIthIndependentBound(j,false);
                left  = GetIthIndependent(j);
            }
            
            
            temp=right-left>0.001?
            0.00001:
            (right-left)*0.0001;
            
            while (right-left>temp) {
                middle                  = (left+right)/2.0;
                SetIthIndependent         (j,middle);
                hyFloat                adjusted = GetIthDependent          (badVarIndex);
                if ((tooLow && adjusted > lb) || (!tooLow && adjusted < lb)) { // in range
                    if (decrement) {
                        left = middle;
                    } else {
                        right = middle;
                    }
                    
                } else {
                    if (decrement) {
                        right = middle;
                    } else {
                        left = middle;
                    }
                }
            }
            // take "right" as the new value for the ind
            SetIthIndependent(j,decrement?left:right);
        }
        
        // now we check to see whether all is well
        
        if (index==-1) {
            // verify that all dependent vars are now within allowed domains
            for (index = indexDep.lLength-1; index>-1; index--)
                // check whether any of the dependent variables are out of bounds
            {
                currentValues[index]=GetIthDependent(index);
                if (currentValues[index]<lowerBounds[index] || currentValues[index]>upperBounds[index]) {
                    break;
                }
            }
            if (index == -1) { // we did it!
                return;
            }
        }
        
        // if we get here, then some of the dependent variables couldn't be moved withing bound
        // we will try 10,000 random assignments of independent variables
        // hoping that one of them will do the trick
        // reuse the first two rows of "dependancies" to store
        // the lower and upper bounds for the dependent vars
        
        //fprintf (stderr, "Trying random permutations...\n");
        
        // compile the list of variables that are tied in with bad constraints
        _SimpleList _aux;
        _AVLList    tagged ( &_aux);
        
        for (long k = 0; k < badIndices.lLength; k++) {
            for (long k2 = 0; k2 < indexInd.lLength; k2++) {
                if (dependancies(badIndices.get(k), k2) != 0.0) {
                    tagged.Insert((BaseRef)k2);
                }
            }
        }
        
        tagged.ReorderList();
        
        // fprintf (stderr, "Tagged the following variables %s\n", _String((_String*)_aux.toStr()).sData);
        
        
        for (index = 0; index<indexInd.lLength; index++) {
            dependancies.Store (0,index,GetIthIndependentBound (index,true));
            dependancies.Store (1,index,(GetIthIndependentBound (index,false)>10?10:GetIthIndependentBound (index,true))-dependancies(0,index));
            dependancies.Store (2,index,GetIthIndependent (index));
        }
        
        // fprintf (stderr, "\n%s\n", _String((_String*)dependancies.toStr()).sData);
        
        
        for (i = 0L; i < 10000L; i++) {
            // choose random values for the variables that are involved with bad constraints
            for (long v = 0L; v < _aux.lLength; v++) {
                index = _aux.get(v);
                SetIthIndependent   (index,dependancies(0,index)+genrand_real2()*dependancies(1,index));
                //fprintf (stderr, "[%d] %s => %g\n", index, GetIthIndependentName(index)->sData, GetIthIndependent(index));
            }
            for (j = 0; j < nonConstantDep->lLength; j++) {
                // check whether any of the dependent variables are out of bounds
                long j_corrected = nonConstantIndices.get(j);
                currentValues.theData[j_corrected]    =   LocateVar(nonConstantDep->lData[j])->Compute()->Value();
                //fprintf (stderr, "[%d] %g (%g, %g)\n", j, j_corrected, currentValues.theData[j_corrected], lowerBounds.theData[j_corrected], upperBounds[j_corrected]);
                if (currentValues.theData[j_corrected]<lowerBounds.theData[j_corrected] || currentValues.theData[j_corrected]>upperBounds.theData[j_corrected]) {
                    //fprintf (stderr, "===| CHECK FAILED\n");
                    badConstraint = nonConstantDep->lData[j];
                    break;
                }
            }
            if (j == nonConstantDep->lLength) {
                break;
            }
        }
        
        ohWell = i==10000;
    }
    if (ohWell) {
        cornholio              = LocateVar(badConstraint);
        
      
        //fprintf (stderr, "%d\n", j);
        
        j = nonConstantIndices.get(j);
        
        _TerminateAndDump(_String("Constrained optimization failed, since a starting point within the domain specified for the variables couldn't be found.\nSet it by hand, or check your constraints for compatibility.\nFailed constraint:")
                          & (*cornholio->GetName()) & ":=" & _String((_String*)cornholio->GetFormulaString(kFormulaStringConversionSubstiteValues)) & " must be in [" & lowerBounds[j] & "," & upperBounds[j] &"]. Current value = " & currentValues[j] & ".");
        
        
    }
}
//_______________________________________________________________________________________
inline hyFloat sqr (hyFloat x)
{
    return x*x;
}


//_______________________________________________________________________________________
  void        _LikelihoodFunction::GetInitialValues (void) const {
      // compute pairwise distances for block index
    for (long partition_index=0UL; partition_index<theTrees.lLength; partition_index++) {
      _DataSetFilter const* filter       = GetIthFilter(partition_index);
      _TheTree      * tree         = GetIthTree  (partition_index);
      _Matrix       * eq_freqs     = GetIthFrequencies(partition_index);
      unsigned long tip_count      = filter->NumberSpecies();

      if (filter->GetData()->GetTT()->IsStandardNucleotide() && filter->IsNormalFilter() && tip_count<150 && eq_freqs->IsIndependent()) {
          // if not - use distance estimates

        if (tree->IsDegenerate()) {
          continue;
        }

        unsigned  long  inode_count = 0UL;

        _TreeIterator ti (tree, _HY_TREE_TRAVERSAL_PREORDER | fTreeIteratorTraversalSkipRoot);

        _SimpleList node_to_index_support ;
        _AVLListX node_to_index (&node_to_index_support);
        _List    root_paths;

        bool two = false;

        while (_CalcNode* iterator = ti.Next()) {
           if (ti.IsAtLeaf()) {
            _SimpleList * indexed_path = new _SimpleList;
            for (long node_index = ti.History().countitems() - 2; node_index >= 0; node_index--) {
              long lookup_index = node_to_index.Find ((BaseRef)ti.History().Element(node_index));
              if (lookup_index >= 0) {
                (*indexed_path) << node_to_index.GetXtra(lookup_index);
              }
            }
            indexed_path->Sort();
            _SimpleList nc;
            nc = ti.History();
            root_paths.AppendNewInstance(indexed_path);
          } else {
            node_to_index.Insert ((BaseRef)ti.GetNode(), inode_count++);
          }

          _SimpleList avl_storage;
          _AVLList node_variables (&avl_storage);
          iterator->ScanContainerForVariables(node_variables,node_variables);

          if (node_variables.countitems()>=2) {
            two = true;
          }
        }

          // first of all, construct the matrix of branch sums,
          // which will be topology dependent
          // the branches are associated with the node they terminate with

        _Matrix  theSystem  ((tip_count-1)*tip_count/2,inode_count+tip_count,true,true);

        unsigned long eq_index = 0UL;

        for (unsigned long row = 0UL; row < tip_count-1 ; row++) {
          for (unsigned long column = row + 1UL; column < tip_count; column++, eq_index++) {
            theSystem.Store(eq_index,row,1.0);
            theSystem.Store(eq_index,column,1.0);
            _SimpleList symmetric_diff;
            symmetric_diff.XOR(*(_SimpleList*)root_paths.GetItem(row), *(_SimpleList*)root_paths.GetItem(column));
            for (unsigned long unshared = 0UL; unshared < symmetric_diff.lLength; unshared++) {
              theSystem.Store(eq_index,tip_count + symmetric_diff.Element(unshared) ,1.0);
            }
          }
        }

         _Matrix transversions (theSystem.GetHDim(),1,false,true),
                 transitions   (theSystem.GetHDim(),1,false,true),
                 jc (theSystem.GetHDim(),1,false,true),
                 transpose(theSystem);

        transpose.Transpose();
        _Matrix AAst (transpose);
        AAst*=theSystem;
        _Matrix* LUSystem = (_Matrix*)AAst.LUDecompose();
        if (!LUSystem) {
          return;
        }
          // bingo - the system is now formulated
          // now we build the right-hand sides.
        AAst.Clear();
        _Matrix *freq;
        if (filter->GetUnitLength()==1) {
          freq = eq_freqs;
          eq_freqs->AddAReference();
        } else {
          freq = filter->HarvestFrequencies (1,1,false);
        }

        _Matrix diffs (1,7,false,true);

          // compute the universal frequency weights
        hyFloat tmp,t2,
        cR = (*freq)[0]+(*freq)[2],
        cY = (*freq)[1]+(*freq)[3],
        c1=2*(*freq)[0]*(*freq)[2]/cR,
        c2=2*(*freq)[1]*(*freq)[3]/cY,
        c3=2*((*freq)[0]*(*freq)[2]*cY/cR
              +(*freq)[1]*(*freq)[3]*cR/cY), comps, fP = 1.0-sqr((*freq)[0])-sqr((*freq)[1])-sqr((*freq)[2])-sqr((*freq)[3]);

        DeleteObject(freq);

        eq_index = 0UL;
        for (unsigned long row = 0UL; row <tip_count - 1; row++)
          for (unsigned long column = row +1UL; column<tip_count; column++, eq_index++) {
              // Tamura - Nei distance formula
            filter->ComputePairwiseDifferences (diffs,row,column) ;
            hyFloat Q = diffs.theData[1]+diffs.theData[3]+
            diffs.theData[4]+diffs.theData[6],P1=diffs.theData[2],
            P2=diffs.theData[5];
            comps = Q+P1+P2+diffs.theData[0];
            if (comps==0.0) {
              continue;
            }
            if (two) {
              Q/=comps;
              P1/=comps;
              P2/=comps;
              tmp = 1-Q/(2*cR*cY);
              if (tmp>0.0) {
                transversions[eq_index] = -2*cR*cY*log(tmp);
              } else {
                transversions[eq_index] = 200*cR*cY;
              }

              tmp = 1.0-1.0/c1*P1-.5/cR*Q;
              if (tmp>0.0) {
                t2 = -c1*log(tmp);
              } else {
                t2 = c1*100;
              }

              tmp = 1.0-1.0/c2*P2-.5/cY*Q;
              if (tmp>0.0) {
                t2 -= c2*log(tmp);
              } else {
                t2 += c2*100;
              }

              tmp = 1-.5/(cR*cY)*Q;
              if (tmp>0.0) {
                t2 += c3*log(tmp);
              } else {
                t2 += c3*100;
              }

              transitions[eq_index]= t2;
            }
            hyFloat P = (Q+P1+P2)/comps;
            P = 1.0-P/fP;
            if (P>0) {
              jc[eq_index]=-fP*log(P);
            } else {
              jc[eq_index]=20.0;
            }
          }

        _Matrix *trstEst=nil,
        *trvrEst=nil,
        *jcEst=nil;

        theSystem=transpose;
        theSystem*=jc;
        jc.Clear();
        jcEst = (_Matrix*)LUSystem->LUSolve(&theSystem);
        if (two) {
          theSystem=transpose;
          transpose*=transversions;
          theSystem*=transitions;
          transitions.Clear();
          transversions.Clear();
          trstEst = (_Matrix*)LUSystem->LUSolve(&theSystem);
          trvrEst = (_Matrix*)LUSystem->LUSolve(&transpose);
        }
        DeleteObject (LUSystem);
          // now that we have the estimates set the values of branch parameters
          // the numbering is as follows
          // the tips : 0 to # tips in the same order as in the tree string
          // the branches: # tips + 1 to # tips + #branches
          // within each branch:
          // if there is one ind variable present: set to sum of both
          // if two parameters present - set first one to transition and the 2nd one to transversion
          // for the third parameter and afterwards use the average of the two
          // produce a list of depthwise traversed nodes

        ti.Reset();
        ti.Next(); // skip the root

        unsigned long tip_index   = 0UL,
                      inode_index = 0UL;

        while (_CalcNode* iterator = ti.Next()) {
          _SimpleList   independent_vars_l,
                        dependent_vars_l;

          _AVLList independent_vars (&independent_vars_l),
                    dependent_vars  (&dependent_vars_l);

          iterator->ScanContainerForVariables(independent_vars,dependent_vars);

          independent_vars.ReorderList();
          dependent_vars.ReorderList();

          if (ti.IsAtLeaf()) {
            c3 = (*jcEst)[tip_index];
            if (two) {
              c1 = (*trstEst)[tip_index];
              c2 = (*trvrEst)[tip_index];
            }
            tip_index++;
          } else {
            c3 = (*jcEst)[tip_count + inode_index];
            if (two) {
              c1 = (*trstEst)[tip_count + inode_index];
              c2 = (*trvrEst)[tip_count + inode_index];
            }
            inode_index++;
          }
          if (independent_vars_l.countitems()==1) {
            _Variable * branch_parameter = LocateVar(independent_vars_l.GetElement(0));
            if (!branch_parameter->HasBeenInitialized () || !branch_parameter->HasChanged()) {
              ReportWarning(_String("Initial guess for ") & *branch_parameter->GetName() & " is " & (c3*.66667));
              branch_parameter->CheckAndSet (c3*.66667, true);
            }
          } else if (independent_vars_l.countitems() >=2) {
            for (unsigned long int p = 0UL; p < independent_vars_l.lLength; p++) {
              _Variable * branch_parameter = LocateVar(independent_vars_l.GetElement(p));
              if (! branch_parameter->HasBeenInitialized () || !branch_parameter->HasChanged()) {
                c3 = p == 0 ? c1 : ( p > 1 ? (c1+c2)*0.5 : c2);
                ReportWarning(_String("Initial guess for ") & *branch_parameter->GetName() & " is " & (c3));
                branch_parameter->CheckAndSet (c3, true);
              }
            }
          }
        }

        DeleteObject(trstEst);
        DeleteObject(trvrEst);
        DeleteObject(jcEst);
      } else {
        hyFloat      initValue = hy_env::EnvVariableGetNumber(kGlobalStartingPoint, 0.1);

        _SimpleList     indeps;
        _AVLList        iavl (&indeps);

        tree->ScanContainerForVariables (iavl, iavl);

        for (long vc = 0; vc < indeps.lLength; vc++) {
            //char buf[512];
          _Variable * localVar = LocateVar(indeps.lData[vc]);
          if (localVar->IsIndependent() && !localVar->HasChanged() && !localVar->IsGlobal()) {
            localVar->CheckAndSet (initValue);
              //      snprintf (buf, sizeof(buf),"[PRESET]%s = %g\n", localVar->GetName()->sData, localVar->Compute()->Value());
          }
            //else
            //  snprintf (buf, sizeof(buf),"[PRESET]%s = %g\n", localVar->GetName()->sData, localVar->Compute()->Value());
            //BufferToConsole(buf);
        }


      }
    }
  }
//_______________________________________________________________________________________

void        _LikelihoodFunction::SetReferenceNodes (void) {
    /*
    hyFloat          cv;

    checkParameter      (useDuplicateMatrixCaching, cv, 0.);

    if (cv>0.5) {
        _List               mappedNodes;
        _SimpleList         mappedTo,
                            canMap;

        for (unsigned i=0UL; i<theTrees.lLength; i++) {
            _TreeIterator ti (GetIthTree(i), _HY_TREE_TRAVERSAL_POSTORDER);

            while (_CalcNode* iterator = ti.Next()) {
                long rV = iterator->CheckForReferenceNode ();
                if (rV >= 0) {
                    mappedNodes << iterator;
                    mappedTo    << rV;
                } else {
                    canMap << iterator->get_index();
                }
            }
        }

        if (mappedNodes.lLength) {
            canMap.Sort();
            for (unsigned long i=0UL; i<mappedNodes.lLength; i++) {
                if (canMap.BinaryFind (mappedTo.lData[i])>=0) {
                    _CalcNode * travNode = (_CalcNode*)mappedNodes(i);
                    travNode->SetRefNode (mappedTo.lData[i]);
                    ((_CalcNode*)LocateVar(mappedTo.lData[i]))->AddRefNode();
                    ReportWarning (_String ("Matrix for node ") & *travNode->GetName() & " mapped to " &
                                   *LocateVar(mappedTo.lData[i])->GetName());
                }
            }
        }
    }*/
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::InitMPIOptimizer (void)
{
#ifdef __HYPHYMPI__
    parallelOptimizerTasks.Clear();
    transferrableVars  = 0;

    int                      totalNodeCount = RetrieveMPICount (0);
    hyFloat               aplf = 0.0;
    checkParameter           (autoParalellizeLF,aplf,0.0);
    hyphyMPIOptimizerMode  = round(aplf);

    if (hyphyMPIOptimizerMode == _hyphyLFMPIModeREL)
        // this is the branch to deal with multiple rate categories
    {

        long        cacheSize  = PartitionLengths (0),
                    categCount = TotalRateClassesForAPartition(-1);

        if (categCount == 0) {
            ReportWarning ("[MPI] Cannot use REL MPI optimization mode on nodes likelihood functions without rate categories");
            return;
        }

        if (totalNodeCount <= categCount) {
            ReportWarning ("[MPI] Cannot initialize REL MPI optimization because there must be at least one more node than there are rate categories");
            return;
        }

        if (theTrees.lLength != 1) {
            ReportWarning ("[MPI] Cannot initialize REL MPI because it does not support likelihood functions with multiple partitions");
            return;
        }

        if (computingTemplate != _hyphyLFComputationalTemplateNone)
            // have a custom computing template
        {
            ReportWarning ("[MPI] Cannot initialize REL MPI because it does not support likelihood functions with computational templates");
            return;
        }

        for (long i=0; i<indexCat.lLength; i++) {
            _CategoryVariable* thisCV = (_CategoryVariable*) LocateVar (indexCat.lData[i]);
            if (thisCV->IsHiddenMarkov()||thisCV->IsConstantOnPartition()) {
                ReportWarning ("[MPI] Cannot initialize REL MPI because it does not support 'Hidden Markov' or 'Constant on Partition' category variables");
                return;
            }
        }

        ReportWarning    (_String ("[MPI] with:") & categCount & " categories on " & (long)totalNodeCount & " MPI nodes");

        MPISwitchNodesToMPIMode (categCount);

        _StringBuffer          sLF (8192L);
        SerializeLF      (sLF,_hyphyLFSerializeModeCategoryAsGlobal);
        sLF.TrimSpace    ();
        long             senderID = 0;

        // sets up a de-categorized LF on each compute node
        for (long i = 1; i<=categCount; i++) {
            MPISendString (sLF,i);
        }

        for (long i = 1; i<=categCount; i++) {
            _String     *mapString  = MPIRecvString (i,senderID);

            // this will return the ';' separated string of partition variable names
            _List        varNames           = mapString->Tokenize (";"),
                         slaveNodeMapL;
            _AVLListX    slaveNodeMap       (&slaveNodeMapL);
            slaveNodeMap.PopulateFromList   (varNames);

            _SimpleList varMap,
                        map2;


            for (long vi = 0; vi < indexCat.lLength; vi++) {
                _Variable * catVar   = LocateVar(indexCat.lData[vi]);
                _String * searchTerm = catVar->GetName();
                long f = slaveNodeMap.Find (searchTerm);
                if (f<0) {
                    HandleApplicationError (_String ("[MPI] InitMPIOptimizer failed to map category var ") & *searchTerm & " for MPI node " & i &".\nHad the following variables\n" & _String((_String*)slaveNodeMap.toStr()));
                    return ;
                }
                varMap << slaveNodeMap.GetXtra(f);
                map2   << catVar->get_index();
                //printf ("%s->%d\n", searchTerm->sData, slaveNodeMap.GetXtra(f));
            }


            for (long vi = 0; vi < indexInd.lLength; vi++) {
                _Variable * indepVar  = LocateVar(indexInd.lData[vi]);
                _String  *  searchTerm = indepVar->GetName();
                long    f = slaveNodeMap.Find (searchTerm);
                if (f<0 && !indepVar->IsGlobal()) {
                    HandleApplicationError (_String ("[MPI] InitMPIOptimizer failed to map independent variable ") & *searchTerm & " for MPI node " & i &".\nHad the following variables\n" & _String((_String*)slaveNodeMap.toStr()));
                    return;
                } else if (f>0) {
                    f = slaveNodeMap.GetXtra(f);
                }

                if (f>=0) {
                    varMap << f;
                    map2   << indepVar->get_index();
                }
            }

            _SimpleList* mappedVariables = new _SimpleList (map2.lLength,0,0);

            for (long vi = 0; vi < map2.lLength; vi++) {
                mappedVariables->lData[varMap.lData[vi]] = map2.lData[vi];
            }
            //printf ("%d->%s\n", varMap.lData[vi], LocateVar( map2.lData[vi])->GetName()->sData);


            parallelOptimizerTasks.AppendNewInstance (mappedVariables);

            if (transferrableVars > 0) {
                if (transferrableVars != mappedVariables->lLength) {
                    HandleApplicationError (_String ("[MPI] InitMPIOptimizer failed: inconsistent variable counts between spawned instances."));
                }
            } else {
                transferrableVars = mappedVariables->lLength;
            }

            DeleteObject (mapString);

            //ReportWarning (((_String*)parallelOptimizerTasks.toStr())->getStr());
        }

        _Matrix::CreateMatrix (&varTransferMatrix, 1, transferrableVars, false, true, false);
        _Matrix::CreateMatrix (&resTransferMatrix, 2, cacheSize, false, true, false);
        ReportWarning(_String("[MPI] InitMPIOptimizer successful. ") & transferrableVars & " transferrable parameters, " & cacheSize & " sites to be cached.");
    } else {
        if (hyphyMPIOptimizerMode == _hyphyLFMPIModePartitions || hyphyMPIOptimizerMode == _hyphyLFMPIModeAuto || hyphyMPIOptimizerMode == _hyphyLFMPIModeSiteTemplate) {
            int     slaveNodes = totalNodeCount-1;

            if (slaveNodes < 2) {
                ReportWarning("[MPI] cannot initialize a partition (or auto) optimizer on fewer than 3 MPI nodes");
                return;
            }



            if (hyphyMPIOptimizerMode == _hyphyLFMPIModePartitions   && theDataFilters.lLength>1 ||
                    hyphyMPIOptimizerMode == _hyphyLFMPIModeAuto         && theDataFilters.lLength == 1 ||
                    hyphyMPIOptimizerMode == _hyphyLFMPIModeSiteTemplate && theDataFilters.lLength>1 && slaveNodes>=theDataFilters.lLength) {

                if (hyphyMPIOptimizerMode == _hyphyLFMPIModeSiteTemplate) {
                    if (templateKind != _hyphyLFComputationalTemplateBySite) {
                        ReportWarning ("[MPI] By site template optimizer cannot work with a not by site computational template");
                        return;
                    }
                } else if (templateKind != _hyphyLFComputationalTemplateNone) {
                    ReportWarning ("[MPI] Partition/Auto Parallel Optimizer does not presently work with computational templates");
                    return;
                }


                long             senderID = 0,
                                 fromPart = 0,
                                 toPart   = 0;

                transferrableVars = 0;
                _List             masterNodeMapL;
                _AVLListX         masterNodeMap (&masterNodeMapL);

                for (long k = 0; k < indexInd.lLength; k++) {
                    masterNodeMap.Insert(LocateVar(indexInd.lData[k])->GetName()->makeDynamic(),indexInd.lData[k],false);
                }
                for (long k = 0; k < indexDep.lLength; k++) {
                    masterNodeMap.Insert(LocateVar(indexDep.lData[k])->GetName()->makeDynamic(),indexDep.lData[k],false);
                }

                if ( hyphyMPIOptimizerMode == _hyphyLFMPIModeAuto ) {
                    _SimpleList    * optimalOrder = (_SimpleList*)(optimalOrders (0));
                    _DataSetFilter const * aDSF   = GetIthFilter(0L);

                    hyFloat    minPatternsPerNode = 0.;
                    checkParameter (minimumSitesForAutoParallelize, minPatternsPerNode, 50.);

                    // adjust slaveNodes as needed
                    slaveNodes = MIN(slaveNodes, ceil(optimalOrder->lLength/minPatternsPerNode));

                    long          sitesPerNode  =  optimalOrder->lLength / slaveNodes,
                                  overFlow         =  optimalOrder->lLength % slaveNodes,
                                  unitLength  =  aDSF->GetUnitLength();


                    _SimpleList const * dupMap    = &aDSF->duplicateMap,
                                * orOrder   = &aDSF->theOriginalOrder;

                    totalNodeCount     = slaveNodes + 1;


                    ReportWarning    (_String ("InitMPIOptimizer (autoLF) with:") & (long)optimalOrder->lLength & " site patterns on " & (long)slaveNodes
                                      & " MPI computational nodes. " & sitesPerNode & " site patterns per node (+1 for "
                                      & overFlow & " nodes)");

                    MPISwitchNodesToMPIMode (slaveNodes);
                    for (long i = 1; i<=slaveNodes; i++) {
                        toPart = sitesPerNode;
                        if (overFlow) {
                            // add an extra site when needed
                            overFlow--;
                            toPart++;
                        }

                        if (fromPart+toPart > optimalOrder->lLength || i == slaveNodes)
                            // ensure that all remaining sites are added to the final block
                        {
                            toPart = optimalOrder->lLength-fromPart;
                        }

                        _SimpleList     map     (optimalOrder->lLength, 0, 0),
                                        subset;

                        for (long i2 = 0; i2 < toPart; i2++)
                            // lay the sites out in optimal order
                        {
                            map.lData[optimalOrder->lData[fromPart+i2]] = 1;
                        }

                        //ReportWarning    (_String((_String*)map.toStr()));

                        for (long i2 = 0; i2 < dupMap->lLength; i2++)
                            if (map.lData[dupMap->lData[i2]])
                                for (long cp = 0; cp < unitLength; cp++ ) {
                                    subset << orOrder->lData[i2*unitLength+cp];
                                }

                        ReportWarning (_String((_String*)subset.toStr()));
                        fromPart += toPart;

                        _StringBuffer          sLF (8192UL);
                        SerializeLF      (sLF,_hyphyLFSerializeModeVanilla,nil,&subset);
                        sLF.TrimSpace()  ;

                        MPISendString (sLF,i);
                        parallelOptimizerTasks.AppendNewInstance (new _SimpleList);
                    }
                }
                // no autoParallelize
                else {
                    long     perNode        =  (hyphyMPIOptimizerMode == _hyphyLFMPIModeAuto)?1:theDataFilters.lLength / slaveNodes,
                             overFlow       =  (hyphyMPIOptimizerMode == _hyphyLFMPIModeAuto)?0:theDataFilters.lLength % slaveNodes;

                    if (perNode == 0L) { // more nodes than filters
                        slaveNodes         = theDataFilters.lLength;
                        totalNodeCount     = slaveNodes + 1;
                        perNode            = 1L;
                        overFlow           = 0L;
                    }

                    _Matrix partition_weights (theDataFilters.lLength, 2, false, true);
                    
                    if (theDataFilters.lLength > slaveNodes) {
                        
                        for (unsigned long i = 0UL; i < theDataFilters.lLength; i++) {
                            //fprintf (stderr, "\nComputing block %ld\n", i);
                            TimeDifference timer;
                            ComputeBlock(i);
                            hyFloat timeDiff   = timer.TimeSinceStart();
                            partition_weights.Store (i, 0, timeDiff);
                        }
                        
                    } else {
                        for (unsigned long i = 0UL; i < theDataFilters.lLength; i++) {
                            partition_weights.Store (i, 0, 1.);
                        }
                    }
                    
                    _Constant * sum = (_Constant*)partition_weights.Sum();
                    partition_weights *= (slaveNodes/sum->Value());
                    for (unsigned long i = 0UL; i < theDataFilters.lLength; i++) {
                        partition_weights.Store (i, 1, i);
                    }
                    
                    sum->SetValue(0.);
                    _Matrix * sorted_by_weight = (_Matrix*)partition_weights.SortMatrixOnColumn(sum);
                    DeleteObject (sum);
                    

                    MPISwitchNodesToMPIMode (slaveNodes);
                    
                    long current_index = 0L;
                    
                    ReportWarning    (_String ("InitMPIOptimizer with:") & (long)theDataFilters.lLength & " partitions on " & (long)slaveNodes
                                      & " MPI computational nodes. ");
                    
                    for (long i = 1L; i<totalNodeCount; i++) {
                        hyFloat sum = 0.;
                        _SimpleList my_part;
                        do {
                            sum += (*sorted_by_weight) (current_index, 0);
                            my_part << round ((*sorted_by_weight) (current_index, 1));
                            current_index++;
                        } while (sum < 1. && theDataFilters.lLength - current_index >= (totalNodeCount - i));
                        
                        ReportWarning    (_String ("InitMPIOptimizer sending partitions ") & _String ((_String*)my_part.toStr()) & " to node " & i);
                        
                        
                        //fprintf (stderr, "%s\n", _String ((_String*)my_part.toStr()).getStr());
                        
                        
                        _StringBuffer     sLF (8192L);
                        SerializeLF       (sLF,_hyphyLFSerializeModeVanilla,&my_part);
                        sLF.TrimSpace     ();
                        
                        MPISendString    (sLF,i);
                        parallelOptimizerTasks.AppendNewInstance (new _SimpleList);
                    }
                    
                    
                    DeleteObject (sorted_by_weight);
                    
                }


                for (long i = 1; i<=parallelOptimizerTasks.lLength; i++) {
                    _String         *mapString      =   MPIRecvString (-1,senderID);
                    _List           varNames       =   mapString->Tokenize (";");
                    _SimpleList*    indexedSlaveV   =   (_SimpleList*)parallelOptimizerTasks(senderID-1);

                    for (long i2 = 0; i2 < varNames.lLength; i2++) {
                        long vi2    = masterNodeMap.Find(varNames.GetItem(i2));
                        if (vi2 < 0)
                            HandleApplicationError (_String ("[MPI] InitMPIOptimizer: Failed to map independent variable ")
                                       & ((_String*)varNames.GetItem(i2))->Enquote()
                                       & " for MPI node " & senderID &". Had variable string:"
                                       & mapString->Enquote());


                        (*indexedSlaveV) << masterNodeMap.GetXtra(vi2);
                    }

                    if (indexedSlaveV->lLength > transferrableVars) {
                        transferrableVars = indexedSlaveV->lLength;
                    }

                    DeleteObject (mapString);
                }

                _Matrix::CreateMatrix (&varTransferMatrix, 1, transferrableVars, false, true, false);
                ReportWarning(_String("[MPI] InitMPIOptimizer:Finished with the setup. Maximum of ") & transferrableVars & " transferrable parameters.");
            }
        }
    }
#endif //__HYPHYMPI__
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::CleanupMPIOptimizer (void)
{
#ifdef __HYPHYMPI__
    if (hyphyMPIOptimizerMode!=_hyphyLFMPIModeNone && parallelOptimizerTasks.lLength) {
        varTransferMatrix.theData[0] = -1.e101;


        for (long i=0; i<parallelOptimizerTasks.lLength; i++) {
            ReportMPIError(MPI_Send(varTransferMatrix.theData, ((_SimpleList*)parallelOptimizerTasks(i))->lLength , MPI_DOUBLE, i+1, HYPHY_MPI_VARS_TAG, MPI_COMM_WORLD),true);
            MPISendString (kEmptyString, i+1);
        }

        if (hyphyMPIOptimizerMode == _hyphyLFMPIModeREL) {
            resTransferMatrix.Clear();
        }
        varTransferMatrix.Clear();
        parallelOptimizerTasks.Clear();
    }
    hyphyMPIOptimizerMode = _hyphyLFMPIModeNone;
#endif
}

//_______________________________________________________________________________________
void            _LikelihoodFunction::SetupLFCaches              (void) {
    // need to decide which data represenation to use,
    // large trees short alignments
    // an acceptable cache size etc
    categID = 0;
    conditionalInternalNodeLikelihoodCaches = new hyFloat*   [theTrees.lLength];
    branchCaches                            = new hyFloat*   [theTrees.lLength];
    siteScalingFactors                      = new hyFloat*   [theTrees.lLength];
    conditionalTerminalNodeStateFlag        = new long*         [theTrees.lLength];
    overallScalingFactors.Populate                        (theTrees.lLength, 0,0);
    overallScalingFactorsBackup.Populate                  (theTrees.lLength, 0,0);
    matricesToExponentiate.Clear();

#ifdef MDSOCL
	OCLEval = new _OCLEvaluator[theTrees.lLength]();
#endif

    evalsSinceLastSetup = 0L;

    for (unsigned long i=0UL; i<theTrees.lLength; i++) {
        _TheTree * cT = GetIthTree(i);
        _DataSetFilter const *theFilter = GetIthFilter(i);

        conditionalInternalNodeLikelihoodCaches[i] = nil;
        conditionalTerminalNodeStateFlag       [i] = nil;
        siteScalingFactors                     [i] = nil;
        branchCaches                           [i] = nil;

        if (!theFilter->IsNormalFilter()) {
            siteCorrections < new _SimpleList;
            siteCorrectionsBackup < new _SimpleList;
            conditionalTerminalNodeLikelihoodCaches< new _Vector;
            continue;
        }

        unsigned long patternCount   = theFilter->GetPatternCount(),
             stateSpaceDim    = theFilter->GetDimension (),
             leafCount      = cT->GetLeafCount(),
             iNodeCount        = cT->GetINodeCount(),
             atomSize      = theFilter->GetUnitLength();

        long ambig_resolution_count = 1L;

        if (leafCount > 1UL) {
            conditionalInternalNodeLikelihoodCaches[i] = new hyFloat [patternCount*stateSpaceDim*iNodeCount*cT->categoryCount];
            branchCaches[i]                            = new hyFloat [2*patternCount*stateSpaceDim*cT->categoryCount];
        }

        siteScalingFactors[i]                          = new hyFloat [patternCount*iNodeCount*cT->categoryCount];
        conditionalTerminalNodeStateFlag[i]            = new long       [patternCount*MAX(2,leafCount)];

        cachedBranches < new _SimpleList (cT->categoryCount,-1,0);
        if (cT->categoryCount == 1UL) {
            siteCorrections < new _SimpleList (patternCount,0,0);
            siteCorrectionsBackup < new _SimpleList (patternCount,0,0);
        } else {
            siteCorrections < new _SimpleList (cT->categoryCount*patternCount,0,0);
            siteCorrectionsBackup < new _SimpleList (cT->categoryCount*patternCount,0,0);
        }

        InitializeArray(siteScalingFactors[i] , patternCount*iNodeCount*cT->categoryCount, 1.);

        // now process filter characters by site / column

        _List        foundCharactersAux;
        _AVLListX    foundCharacters (&foundCharactersAux);
        _String      aState ((unsigned long)atomSize);

        char  const ** columnBlock      = new char const*[atomSize];
        hyFloat      * translationCache  = new hyFloat [stateSpaceDim];
        _Vector  * ambigs            = new _Vector();

        for (unsigned long siteID = 0UL; siteID < patternCount; siteID ++) {
            siteScalingFactors[i][siteID] = 1.;
            for (unsigned long k = 0UL; k < atomSize; k++) {
                columnBlock[k] = theFilter->GetColumn(siteID*atomSize+k);
            }

            long uptoL = Maximum (2UL,leafCount);

            for (long leafID = 0; leafID < uptoL; leafID ++) {
                long mappedLeaf  = theFilter->theNodeMap.lData[leafID],
                     translation;

                for (long k = 0; k < atomSize; k++) {
                    aState.set_char (k, columnBlock[k][mappedLeaf]);
                }

                translation = foundCharacters.Find (&aState);
                if (translation < 0L) {
                    translation = theFilter->Translate2Frequencies (aState, translationCache, true);
                    if (translation < 0L) {
                        for (unsigned long j = 0UL; j < stateSpaceDim; j++) {
                            ambigs->Store(translationCache[j]);
                        }
                        translation = -ambig_resolution_count++;
                    }
                    foundCharacters.Insert (new _String(aState), translation);
                } else {
                    translation = foundCharacters.GetXtra (translation);
                }
                conditionalTerminalNodeStateFlag [i][leafID*patternCount + siteID] = translation;
            }
        }
        conditionalTerminalNodeLikelihoodCaches < ambigs;
        delete [] columnBlock;
        delete [] translationCache;

#ifdef MDSOCL
		OCLEval[i].init(patternCount, theFilter->GetDimension(), conditionalInternalNodeLikelihoodCaches[i]);
#endif
    }
}

//extern long marginalLFEvals, marginalLFEvalsAmb;

//_______________________________________________________________________________________

void        _LikelihoodFunction::LoggerLogL (hyFloat logL) {
  if (optimizatonHistory) {

    #ifdef  _COMPARATIVE_LF_DEBUG_CHECK
          if (_comparative_lf_debug_matrix && fabs ((*_comparative_lf_debug_matrix)[_comparative_lf_index] - logL) > 0.001) {
              char buffer [512];
              snprintf (buffer, 511, "[Divergent results obtained at EVAL %d] [logged] %15.12g vs [computed] %15.12g\n", _comparative_lf_index, (*_comparative_lf_debug_matrix)[_comparative_lf_index], logL);
              StringToConsole (buffer);
              for (unsigned long var_id = 0UL; var_id < indexInd.lLength; var_id++) {
                  StringToConsole (*GetIthIndependentName(var_id));
                  StringToConsole (_String (" = ") & GetIthIndependent(var_id));
                  NLToConsole();
              }

              raise(SIGTRAP);
          }
          _comparative_lf_index++;
    #endif
    #ifdef  _COMPARATIVE_LF_DEBUG_DUMP
          if (_comparative_lf_debug_matrix) {
              (*_comparative_lf_debug_matrix) << logL;
          }
    #endif


    *((_Vector*) this->optimizatonHistory->GetByKey("LogL")) << logL
    << ((_AssociativeList*)this->optimizatonHistory->GetByKey("Phases"))->Length();
  }
}

//_______________________________________________________________________________________

void        _LikelihoodFunction::LoggerAddGradientPhase (hyFloat precision) {
  if (optimizatonHistory) {
    _AssociativeList* new_phase = new _AssociativeList;
    (*new_phase) < (_associative_list_key_value){"type", new _FString ("Gradient descent")}
                 < (_associative_list_key_value){"precision", new _Constant (precision)};


     *((_AssociativeList*) this->optimizatonHistory->GetByKey("Phases")) < (_associative_list_key_value){nil, new_phase};
  }
}

//_______________________________________________________________________________________

void        _LikelihoodFunction::LoggerAddCoordinatewisePhase (hyFloat shrinkage, char convergence_mode) {
  if (optimizatonHistory) {
    _String phase_kind;
    switch (convergence_mode) {
      case 0:
        phase_kind = "Normal";
        break;
      case 1:
        phase_kind = "Accelerated convergence";
        break;
      case 2:
        phase_kind = "Slow convergence";
        break;
      case 3:
        phase_kind = "Very slow convergence";
        break;
    }

    _AssociativeList* new_phase = new _AssociativeList;
    (*new_phase) < (_associative_list_key_value){"type", new _FString ("Directional pass")}
    < (_associative_list_key_value){"shrinkage", new _Constant (shrinkage)}
    < (_associative_list_key_value){"mode", new _FString(phase_kind)};


    *((_AssociativeList*) this->optimizatonHistory->GetByKey("Phases")) < (_associative_list_key_value){nil, new_phase};
  }
}

//_______________________________________________________________________________________

void        _LikelihoodFunction::LoggerAllVariables (void) {
  if (optimizatonHistory) {
    _AssociativeList* variables = ((_AssociativeList*)this->optimizatonHistory->GetByKey("Parameters"));
    for (unsigned long var_id = 0UL; var_id < indexInd.lLength; var_id++) {
      *((_Vector*) variables->GetByKey(*GetIthIndependentName(var_id))) << GetIthIndependent(var_id);
    }
  }
}

//_______________________________________________________________________________________

void        _LikelihoodFunction::LoggerSingleVariable        (unsigned long index, hyFloat logL, hyFloat bracket_precision, hyFloat brent_precision, hyFloat bracket_width, unsigned long bracket_evals, unsigned long brent_evals) {

  if (optimizatonHistory) {
    _AssociativeList* new_phase = new _AssociativeList;
    (*new_phase) < (_associative_list_key_value){"type", new _FString (*GetIthIndependentName(index))}
                < (_associative_list_key_value){"bracket precision", new _Constant (bracket_precision)}
                < (_associative_list_key_value){"brent precision", new _Constant (brent_precision)}
                < (_associative_list_key_value){"bracket width", new _Constant (bracket_width)}
                < (_associative_list_key_value){"bracket evals", new _Constant (bracket_evals)}
                < (_associative_list_key_value){"brent evals", new _Constant (brent_evals)}
                < (_associative_list_key_value){"brent evals", new _Constant (brent_evals)};

    *((_AssociativeList*) this->optimizatonHistory->GetByKey("Phases")) < (_associative_list_key_value){nil, new_phase};

    LoggerLogL (logL);
    *((_Vector*) (((_AssociativeList*)this->optimizatonHistory->GetByKey("Parameters")))->GetByKey(*GetIthIndependentName(index))) << GetIthIndependent(index);
  }
}

//_______________________________________________________________________________________________

void _LikelihoodFunction::SetupParameterMapping (void) {
    parameterTransformationFunction.Clear();
    parameterValuesAndRanges = new _Matrix (indexInd.lLength, 4, false, true);
    smoothingPenalty   = hy_env::EnvVariableGetNumber(kAddLFSmoothing, 0.0);
    smoothingReduction = hy_env::EnvVariableGetNumber(kReduceLFSmoothing, 0.8);
    if (smoothingPenalty < 0.0) {
        smoothingPenalty = 0.0;
    }
    if (smoothingReduction <= 0.0 || smoothingReduction >= 1.0) {
        smoothingReduction = 0.8;
    }
    
    for (unsigned long pIndex = 0; pIndex < indexInd.lLength; pIndex++) {
        _Variable* cv        = GetIthIndependentVar(pIndex);
        hyFloat thisLB    = cv->GetLowerBound(),
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

//_______________________________________________________________________________________

_Matrix*        _LikelihoodFunction::Optimize (_AssociativeList const * options) {
    // various optimization-only env variables
    
    
    static const _String
    
        kProduceOptimizationLog                        ("PRODUCE_OPTIMIZATION_LOG"),
        kMethodCoordinate                              ("coordinate-wise"),
        kMethodNedlerMead                              ("nedler-mead"),
        kMethodHybrid                                  ("hybrid"),
        kMethodGradientDescent                         ("gradient-descent");

        // optimization setting to produce a detailed log of optimization runs

    auto   get_optimization_setting = [options] (const _String& arg, const hyFloat def_value) -> hyFloat {
        if (options) {
            try {
                return options->GetNumberByKey(arg);
            } catch (const _String err) {
            }
        }
        return hy_env::EnvVariableGetNumber(arg, def_value);
    };

    auto   get_optimization_setting_string = [options] (const _String& arg) -> _FString* {
        if (options) {
            _FString * res = (_FString*)options->GetByKey(arg, STRING);
            if (res) {
                return res;
            }
        }
        return (_FString *)hy_env::EnvVariableGet(arg, STRING);
    };
    
    
    enum   optimization_mode_enum {
        kOptimizationCoordinateWise,
        kOptimizationNedlerMead,
        kOptimizationGradientDescent,
        kOptimizationHybrid
    } optimization_mode = kOptimizationHybrid;
    
    auto get_opt_method_string = [&] () -> const _String {
        if (optimization_mode == kOptimizationCoordinateWise)
            return kMethodCoordinate;
        if (optimization_mode == kOptimizationHybrid)
            return kMethodHybrid;
        if (optimization_mode == kOptimizationGradientDescent)
            return kMethodGradientDescent;
        if (optimization_mode == kOptimizationNedlerMead)
            return kMethodNedlerMead;
    };

    if (lockedLFID != -1) {
        HandleApplicationError ("Optimize() could not be executed, because another optimization is already in progress.");
        return new _Matrix (1,1,false,true);
    }

    /*if (options) {
        ObjectToConsole((_AssociativeList*)options);
    }*/
    
    char           buffer [1024];

    RescanAllVariables ();

    if (optimizatonHistory) {
      DeleteObject (optimizatonHistory);
      optimizatonHistory = nil;
    }

    bool         keepOptimizationLog = ! CheckEqual (get_optimization_setting (kProduceOptimizationLog, 0.0), 0);

    if (keepOptimizationLog) {
      optimizatonHistory = new _AssociativeList;
      (*optimizatonHistory) < (_associative_list_key_value){"LogL", new _Vector}
      /*
       2 values per entry:
       logL ; optimization stage (indexed from 0 to max)
       */
      < (_associative_list_key_value){"Phases", new _AssociativeList};
      /*
       0 - N-1 indices
       */
    }


    if (indexInd.empty()) {
        _Matrix * result = new _Matrix (2UL, Maximum(3UL,indexDep.lLength), false, true);
        PrepareToCompute();
        hyFloat logL = Compute();
        result->Store (1L,0L,logL);
        LoggerLogL    (logL);
        result->Store (1L,1L,0.);
        result->Store (1L,2L,0.);
        for (unsigned long i=0UL; i<indexDep.lLength; i++) {
            result->Store(0,i,this->GetIthDependent(i));
        }
        DoneComputing();
        return result;
    } else {
      if (keepOptimizationLog) {
        _AssociativeList * variable_traces = new _AssociativeList;
        for (unsigned long var_id = 0; var_id < indexInd.lLength; var_id++) {
          (*variable_traces) < (_associative_list_key_value){GetIthIndependentVar(var_id)->GetName()->get_str(), new _Vector};
        }
        (*optimizatonHistory) < (_associative_list_key_value){"Parameters", variable_traces};
      }
    }


    hyFloat     intermediateP,
                maxSoFar            = -INFINITY,
                bestVal,
                lastMax,
                currentPrecision    = 0.1,
                bP,
                percentDone         = 0.0;
    
    long        fnDim               = MaximumDimension(),
                evalsIn             = likeFuncEvalCallCount,
                exponentiationsIn   = matrix_exp_count;


    TimeDifference timer;

    hyFloat                 hardLimitOnOptimizationValue = hy_env::EnvVariableGetNumber(kOptimizationHardLimit, (hyFloat) INFINITY);

    hardLimitOnOptimizationValue = MAX (hardLimitOnOptimizationValue, 0.0);
    if (hardLimitOnOptimizationValue != INFINITY) {
        ReportWarning (_String("Set a hard time limit for optimization routines to ") & hardLimitOnOptimizationValue & " seconds\n");
    }

    isInOptimize    = true;
    lockedLFID      = likeFuncList._SimpleList::Find ((long)this);

    //RankVariables   ();
    VerbosityLevel  ();

#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0) {
#endif
        _FString * pfs = (_FString*)CheckReceptacle (&optimizationStatusFile, kEmptyString, false)->ComputeMatchingType(STRING);
        progressFileString = pfs ? & pfs->get_str() : nil;

        if (verbosity_level==1L) {
            UpdateOptimizationStatus (0,0,0,true,progressFileString);
        }
#ifdef __HYPHYMPI__
    }
#endif
#endif


    for (unsigned long tree_index = 0UL; tree_index <theTrees.lLength; tree_index ++) {
         GetIthTree (tree_index)->CountTreeCategories();
    }

    SetupLFCaches       ();
    SetupCategoryCaches ();
    computationalResults.Clear();


#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0) {
        InitMPIOptimizer ();
        if (parallelOptimizerTasks.lLength == 0) {
            hyphyMPIOptimizerMode = _hyphyLFMPIModeNone;
        }
    }

    if (hyphyMPIOptimizerMode == _hyphyLFMPIModeNone) {
#endif
        SetReferenceNodes();
        for (unsigned long tree_index = 0UL; tree_index <theTrees.lLength; tree_index ++) {
          _TheTree *ith_tree = GetIthTree (tree_index);
          ith_tree->SetUpMatrices(ith_tree->categoryCount);
        }
        hasBeenSetUp = 1;

#ifdef __HYPHYMPI__
    }
#endif

    computationalResults.Clear();

#if !defined __UNIX__ || defined __HEADLESS__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    SetStatusBarValue (0,maxSoFar,0);
#endif


#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0) {
#endif

    BenchmarkThreads (this);


#ifdef __HYPHYMPI__
    }
#endif


    bool            skipCG                  = ! CheckEqual (get_optimization_setting (kSkipConjugateGradient, 0.0), 0.0),
                    keepStartingPoint       = ! CheckEqual (get_optimization_setting (kUseLastResults, 0.0), 0.0),
                    go2Bound                = ! CheckEqual (get_optimization_setting (kAllowBoundary, 1.0), 0.0);
                    //use_distances           = ! CheckEqual (get_optimization_setting (kAllowBoundary, 1.0), 0.0);

    
    _FString * custom_convergence_callback_name = (_FString *)hy_env::EnvVariableGet(hy_env::lf_convergence_criterion, STRING);
    long custom_convergence_callback = custom_convergence_callback_name ?  FindBFFunctionName(*custom_convergence_callback_name->get_str()) : -1L;
    if (custom_convergence_callback >= 0) {
         if (GetBFFunctionArgumentCount (custom_convergence_callback) != 2L) {
             HandleApplicationError("Custom convergence criterion convergence function must have exactly two arguments: current log-L, and an dictionary with id -> value mapping");
             return new _Matrix;
         }
    }
    
    if (keepStartingPoint) {
        indexInd.Each ([this] (long v, unsigned long i) -> void {
            _Variable *iv = GetIthIndependentVar (i);
            if (iv->HasBeenInitialized()) {
                iv->MarkModified();
            }
        });
    }
    
    if (! CheckEqual (get_optimization_setting (kUseInitialDistanceGuess, 1.0), 0.0)) {
        GetInitialValues();
    }

    

    if (!keepStartingPoint) {
        hyFloat  global_starting_point = get_optimization_setting (kGlobalStartingPoint, 0.1);
        
        if (CheckEqual (get_optimization_setting (kRandomStartingPerturbations, 0.0), 0.0)) {
            for (unsigned long i=0UL; i<indexInd.lLength; i++) {
                _Variable *iv = GetIthIndependentVar (i);
                if (iv->IsGlobal()) {
                    if (iv->Value() < 1.e-10) {
                        iv->SetValue(global_starting_point);
                    }
                } else {
                    if (!iv->HasChanged()) {
                        iv->SetValue(global_starting_point);
                    }
                }
            }
        } else {
            for (unsigned long i=0UL; i<indexInd.lLength; i++) {
                _Variable *iv = GetIthIndependentVar (i);
                if (!iv->HasChanged()) {
                    iv->SetValue(global_starting_point*(1.0+genrand_real1()));
                }
            }
        }
    }

#if !defined __UNIX__ || defined __HEADLESS__
    SetStatusBarValue (5,maxSoFar,(likeFuncEvalCallCount-evalsIn)/timer.TimeSinceStart());
#endif

    CheckDependentBounds();

    currentPrecision    = get_optimization_setting (kStartingPrecision, 0.1);
    precision           = get_optimization_setting (kOptimizationPrecision, 0.001);
    maxItersPerVar      = get_optimization_setting (kMaximumIterationsPerVariable, 5000.);
    
    
    
    _FString * method_string = get_optimization_setting_string (kOptimizationMethod);
    
    //optimization_mode
    if (method_string) {
        _StringBuffer const * sm  = &method_string->get_str();
        if (*sm == kMethodCoordinate) {
            optimization_mode = kOptimizationCoordinateWise;
        } else if (*sm == kMethodNedlerMead) {
            optimization_mode = kOptimizationNedlerMead;
        } else if (*sm == kMethodHybrid) {
            optimization_mode = kOptimizationHybrid;
        } else if (*sm == kMethodGradientDescent) {
            optimization_mode = kOptimizationGradientDescent;
        }
    } else {
        switch ((long) get_optimization_setting (kOptimizationMethod, -1.)) {
            case 0:
                optimization_mode = kOptimizationCoordinateWise;
                break;
            case 3:
                optimization_mode = kOptimizationNedlerMead;
                break;
            case 4:
                optimization_mode = kOptimizationHybrid;
                break;
            case 7:
                optimization_mode = kOptimizationGradientDescent;
                break;
        }
        
    }
    
    _OptimiztionProgress progress_tracker;
    //checkParameter (optimizationMethod,optMethodP,4.0);

    
    ReportWarning  (_String("Optimization settings:\n\t") & kOptimizationMethod & " = " & get_opt_method_string () &
                    "\n\t" & kOptimizationPrecision & " = " & precision &
                    "\n\t" & kMaximumIterationsPerVariable & " = " & maxItersPerVar
                    &"\n\nInitial parameter values\n");
    
    
    indexInd.Each([this] (long index, unsigned long idx) -> void {
        ReportWarning (*this->GetIthIndependentName(idx) & " = " & this->GetIthIndependent(idx));
    });
    
    maxItersPerVar *= indexInd.countitems();

#if !defined __UNIX__ || defined __HEADLESS__
#ifdef __HYPHYMPI__
    if (hy_mpi_node_rank == 0)
#endif
        if (terminate_execution) {
            CleanUpOptimize();
            return new _Matrix (1,1,false,true);
        }
#endif


    SetupParameterMapping   ();
    _Matrix variableValues;
    GetAllIndependent (variableValues);


    if (optimization_mode == kOptimizationHybrid && indexInd.countitems() == 1) {
        optimization_mode = kOptimizationCoordinateWise;
    }

    
    if (optimization_mode != kOptimizationCoordinateWise) {
        bP = get_optimization_setting (kBracketingPersistence, 2.5);
    } else {
        bP = get_optimization_setting (kBracketingPersistence, 2.5);
    }

    if (optimization_mode == kOptimizationHybrid || optimization_mode == kOptimizationGradientDescent) { // gradient descent
        _Matrix bestSoFar;

        GetAllIndependent (bestSoFar);
        intermediateP = get_optimization_setting (kIntermediatePrecision, 0.1);
        
        if (verbosity_level>20) {
            snprintf (buffer, sizeof(buffer),"\nGradient Precision %g  Opt Precision = %g", intermediateP, precision);
            BufferToConsole (buffer);
        }

        if (optimization_mode != kOptimizationGradientDescent) {
            //ConjugateGradientDescent (0.5, bestSoFar, true, 10);
            if (gradientBlocks.nonempty()) {
                for (long b = 0; b < gradientBlocks.lLength; b++) {
                    maxSoFar = ConjugateGradientDescent (currentPrecision, bestSoFar,true,10,(_SimpleList*)(gradientBlocks(b)),maxSoFar);
                }
            } else {
                maxSoFar = ConjugateGradientDescent (currentPrecision, bestSoFar,true,10,nil,maxSoFar);
            }
        } else {
            hyFloat current_precision = MAX(1., precision);
            while (current_precision > precision) {
              ConjugateGradientDescent(current_precision, bestSoFar, true);
              current_precision *= 0.1;
            }
            ConjugateGradientDescent(precision, bestSoFar, true);
        }
#if !defined __UNIX__ || defined __HEADLESS__
#ifdef __HYPHYMPI__
        if (hy_mpi_node_rank == 0) {
#endif
            if (terminate_execution) {
                CleanUpOptimize();
                return new _Matrix (1,1,false,true);
            }
#ifdef __HYPHYMPI__
        }
#endif
#endif
        maxSoFar = Compute();
        if (optimization_mode != kOptimizationGradientDescent) {
            optimization_mode = kOptimizationCoordinateWise;
        }
        currentPrecision = kOptimizationGradientDescent==7?sqrt(precision):intermediateP;
        percentDone = 10.0;
    }


    if (optimization_mode == kOptimizationCoordinateWise) {
        bool      forward = false;

        hyFloat averageChange = 0,
                   divFactor     = -1.0 ,
                   oldAverage      = -1.0,
                   stdFactor   = 10,
                   loopCounter    = 0.0,
                   lastMaxValue,
                   averageChange2 = 0.0,
                   //currentPrecision2 = precision,
                   doShuffle,
                   useAdaptiveStep = 0.0;

        long       inCount       = 0,
                   termFactor,
                   lfCount         = likeFuncEvalCallCount;

        _SimpleList noChange,
                    glVars,
                    shuffledOrder;

        _List               *stepHistory = nil;
        _Vector      logLHistory;

        maxSoFar  = lastMaxValue = Compute();

        logLHistory.Store(maxSoFar);

        checkParameter (useAdaptiveVariableStep, useAdaptiveStep, 1.0);

        if (useAdaptiveStep>0.5) {
            stepHistory = new _List;
            for (unsigned long j=0UL; j<indexInd.lLength; j++) {
                _Vector      *varHistory = new _Vector;
                varHistory->Store(variableValues[j]);
                stepHistory->AppendNewInstance(varHistory);
            }
        }

        noChange << -1;

#if !defined __UNIX__ || defined __HEADLESS__
        SetStatusBarValue (percentDone,maxSoFar,(likeFuncEvalCallCount-evalsIn)*TimerDifferenceFunction(true));
#endif

        if (indexInd.lLength<8) {
            stdFactor = 8.0;
        } else if (indexInd.lLength<16) {
            stdFactor = 4.0;
        } else if (indexInd.lLength<32) {
            stdFactor = 2.5;
        } else {
            stdFactor = .5*(1.+sqrt(5.));
        }

        for (unsigned long j=0UL; j<indexInd.lLength; j++) {
            averageChange+=fabs (variableValues[j]-GetIthIndependent(j));
        }

        if (averageChange == 0.0) {
            averageChange = currentPrecision;
        }

        currentPrecision = precision>.1?precision:.1;


        termFactor = stdFactor+1;
        if (termFactor*2 > indexInd.lLength) {
            termFactor = indexInd.lLength >> 1;
        }

        if (termFactor<3) {
            termFactor = 3;
        }

        checkParameter (doShuffleOrder, doShuffle, 0.0);

        for (unsigned long j=0UL; j<indexInd.lLength; j++)
            if (LocateVar (indexInd.lData[j])->IsGlobal()) {
                glVars << j;
            }



        while (inCount<termFactor || smoothingTerm > 0.) {
            if (smoothingTerm > 0. && inCount == termFactor) {
              smoothingTerm = 0.;
              inCount = 0;
            }
            if (likeFuncEvalCallCount-lfCount>maxItersPerVar) {
                ReportWarning ("Optimization routines returning before requested precision goal met. The maximum iteration number specified by MAXIMUM_ITERATIONS_PER_VARIABLE has been reached");
                DeleteObject  (stepHistory);
                stepHistory = nil;
                break;
            }
            if (averageChange <1e-20) {
                averageChange = fabs(oldAverage)<1?fabs(oldAverage*oldAverage):.5;
                if (averageChange<1e-20) {
                    averageChange = 1e-8;
                }
            }

            hyFloat diffs [5] = {0.,0.,0.,0.,0.};
            char convergenceMode = 0;
            /* 0, normal
             * 1, accelerated (last cycle obtained a bigger LL drop than the one before)
             * 2, slow convergence (last three cycles within a factor of 2 LL drop of each other)
             * 3, really slow convergence (last five cycles within a factor of 2 LL drop of each other)
            */


            if (useAdaptiveStep > 0.5) {
                convergenceMode = 0;
                if (oldAverage<0.0) {
                    divFactor = stdFactor;
                } else {
                    divFactor           = MIN(16,MAX(stdFactor,oldAverage/averageChange));

                    long       steps    = logLHistory.get_used();
                    for (long k = 1; k <= MIN(5, steps-1); k++) {
                        diffs[k-1] = logLHistory.theData[steps-k] - logLHistory.theData[steps-k-1];
                        //printf ("%ld : %g\n", k, diffs[k-1]);
                    }
                    if (steps > 2 && diffs[0] >= diffs[1]) {
                        convergenceMode = 1;
                    }

                    if (diffs[0] < precision*0.001) {
                        convergenceMode = 3;
                    }

                    if (convergenceMode < 2) {
                        if (steps > 3) {
                            if (diffs[0] > 0. && diffs[1] > 0. && diffs[2] > 0.) {
                                if (diffs[0] / diffs[1] >= _HY_SLOW_CONVERGENCE_RATIO_INV && diffs[0] / diffs[1] <= _HY_SLOW_CONVERGENCE_RATIO &&
                                        diffs[1] / diffs[2] >= _HY_SLOW_CONVERGENCE_RATIO_INV && diffs[1] / diffs[2] <= _HY_SLOW_CONVERGENCE_RATIO) {
                                    convergenceMode = 2;
                                    if (steps > 4) {
                                        if (diffs [3] > 0.) {
                                            if (diffs[2] / diffs[3] >= _HY_SLOW_CONVERGENCE_RATIO_INV && diffs[2] / diffs[3] <= _HY_SLOW_CONVERGENCE_RATIO) {
                                                convergenceMode = 3;
                                            }
                                        } else {
                                            convergenceMode = 3;
                                        }
                                    }
                                }
                            } else {
                                convergenceMode = 2;
                            }
                        }

                    }

                    switch (convergenceMode) {
                    case 1:
                        divFactor = 1.;
                        break;
                    case 2:
                        divFactor = 4.;
                        break;
                    case 3:
                        divFactor = 10.;
                        break;
                        //default:
                        //  divFactor = 4.;
                    }

                }
            } else {
                if (oldAverage<0.0 || stdFactor*averageChange>oldAverage) {
                    divFactor = stdFactor;
                } else {
                    divFactor = oldAverage/averageChange;
                }
            }

            oldAverage = averageChange;

            averageChange  = 1e-25;


            if (verbosity_level>1) {
                snprintf (buffer, sizeof(buffer),"\n\nOptimization Pass %ld (%ld). LF evalutations : %ld\n", (long)loopCounter, inCount,likeFuncEvalCallCount-lfCount);
                BufferToConsole (buffer);
                if (useAdaptiveStep > 0.5 && logLHistory.get_used() > 2) {
                    snprintf (buffer, sizeof(buffer), "\nLast cycle logL change = %g\n", diffs[0]);
                    BufferToConsole (buffer);
                }
            }
#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
            else if (verbosity_level==1) {
                UpdateOptimizationStatus (maxSoFar,-1,1,true,progressFileString);
            }
#endif
            if (smoothingTerm > 0.) {
              smoothingTerm *= smoothingReduction;
            }

            _SimpleList nc2;

            long        ncp = 0,
                        jjj;

            averageChange2 = 0.0;

            if (doShuffle > 0.1) {
                nc2.Populate (indexInd.lLength, 0, 1);
                shuffledOrder.Subtract (nc2, glVars);
                shuffledOrder.Permute  (1);

                for (unsigned long j=0UL; j<glVars.lLength; j++) {
                    shuffledOrder << glVars.lData[j];
                }

                nc2.Clear();
                shuffledOrder.Flip();
            }

            hyFloat stepScale = 1.;

            if (useAdaptiveStep > 0.5) {
                stepScale = 1/divFactor;
                if (verbosity_level>5) {
                    snprintf (buffer, sizeof(buffer),"\n[BRACKET SHRINKAGE: %g]", divFactor);
                    BufferToConsole (buffer);
                    snprintf (buffer, sizeof(buffer),"\n[Convergence mode = %d]", convergenceMode);
                    BufferToConsole (buffer);
                    snprintf (buffer, sizeof(buffer),"\n[Unchanged variables = %ld]", noChange.lLength);
                    BufferToConsole (buffer);
                }

                if (hardLimitOnOptimizationValue < INFINITY && timer.TimeSinceStart() > hardLimitOnOptimizationValue) {
                    ReportWarning (_String("Optimization terminated before convergence because the hard time limit was exceeded."));
                    break;
                }

                if (convergenceMode > 2) {
                    if (hardLimitOnOptimizationValue < INFINITY && timer.TimeSinceStart() > hardLimitOnOptimizationValue) {
                        ReportWarning (_String("Optimization terminated before convergence because the hard time limit was exceeded."));
                        break;
                    }
                    _Matrix             bestMSoFar;
                    GetAllIndependent   (bestMSoFar);
                    hyFloat prec = Minimum (diffs[0], diffs[1]);

                    prec = Minimum (Maximum (prec, precision), 1.);

                    if (gradientBlocks.nonempty()) {
                        for (long b = 0; b < gradientBlocks.lLength; b++) {
                            maxSoFar = ConjugateGradientDescent (prec, bestMSoFar,true,10,(_SimpleList*)(gradientBlocks(b)),maxSoFar);
                        }
                    } else {
                        maxSoFar = ConjugateGradientDescent (prec, bestMSoFar,true,10,nil,maxSoFar);
                    }

                    GetAllIndependent   (bestMSoFar);
                    for (unsigned long k = 0UL; k < indexInd.lLength; k++) {
                        ((_Vector*)(*stepHistory)(k))->Store (bestMSoFar.theData[k]);
                    }

                    stepScale = 1.;
                    logLHistory.Store(maxSoFar);
                }
            }

            LoggerAddCoordinatewisePhase (divFactor, convergenceMode);

            for (jjj=forward?0:indexInd.lLength-1; forward?(jjj<indexInd.lLength):jjj>=0; forward?jjj++:jjj--) {
                if (hardLimitOnOptimizationValue < INFINITY && timer.TimeSinceStart() > hardLimitOnOptimizationValue) {
                    break;
                }

                unsigned long current_index = doShuffle > 0.1 ? shuffledOrder.lData[jjj] : jjj;

                bool amIGlobal = GetIthIndependentVar(current_index)->IsGlobal();

#ifdef __HYPHYMPI__
                if (hy_mpi_node_rank == 0)
#endif
                    if (terminate_execution) {
                        CleanUpOptimize();
                        return new _Matrix (1,1,false,true);
                    }

                bestVal = GetIthIndependent(current_index);
                lastMax = maxSoFar;

                if (current_index==noChange.lData[ncp]) {
                    if (ncp<noChange.lLength-1) {
                        ncp++;
                    }
                    averageChange2+=bestVal;
                    continue;
                }

                _Vector     *vH = nil;
                hyFloat         precisionStep = 0.,
                                   brackStep;


                if (useAdaptiveStep>0.5) {
                    vH  = (_Vector*)(*stepHistory)(current_index);
                    //hyFloat    suggestedPrecision  = currentPrecision*(1.+198./(1.+exp(sqrt(loopCounter))));

                    long stepsSoFar = vH->get_used();

                    if (stepsSoFar>1) {
                        hyFloat  lastParameterValue          = vH->theData[stepsSoFar-1],
                                    previousParameterValue       = vH->theData[stepsSoFar-2];

                        //stepScale   = 0.25;

                        brackStep     = fabs(lastParameterValue-previousParameterValue);
                        if (brackStep == 0.0) {
                            long k = Maximum(stepsSoFar-3L,0L);

                            for (; k && brackStep == 0.0; k--) {
                                previousParameterValue          = vH->theData[k],
                                lastParameterValue              = vH->theData[k+1];
                                brackStep                       = fabs(lastParameterValue-previousParameterValue);
                            }

                            if (k == 0) {
                                brackStep = MIN(0.001,precision*0.001);
                            }
                        }

                        precisionStep = brackStep*stepScale;

                        if (inCount) {
                            precisionStep = lastParameterValue*precision;
                        }

                        precisionStep = MAX(precisionStep,lastParameterValue*precision);


                        if (precisionStep == 0.0) {
                            precisionStep = precision;
                            brackStep = 2.*precisionStep;
                        } else {
                            brackStep *= 2.;
                        }
                        //if (IsIthParameterGlobal (j))
                        //  precisionStep *= 2.0;

                        if (verbosity_level>50) {
                            snprintf (buffer, sizeof(buffer),"\n[BRACKET STEP: current = %g: previous = %g (diff = %g). bracket = %g, prec = %g]",
                                     lastParameterValue,
                                     previousParameterValue,
                                     lastParameterValue - previousParameterValue,
                                     brackStep,
                                     precisionStep);

                            BufferToConsole (buffer);
                        }

                    } else {
                        brackStep     = vH->theData[0];
                        precisionStep = MAX(0.001,brackStep * 0.1);
                    }

                    /*if (brackStep < 1e-4) {
                        brackStep = 1e-4;
                    }
                    if (precisionStep < 1e-6) {
                        precisionStep = 1e-6;
                    }*/

                } else {
                    if (amIGlobal)
                        brackStep = pow(currentPrecision/**(bestVal>1.?pow(e,long(log(bestVal))):1.)*/,
                                        .5+loopCounter/(indexInd.lLength*(loopCounter/indexInd.lLength+1)));
                    else {
                        brackStep = currentPrecision;
                    }
                }

                long brackStepSave = bracketFCount,
                     oneDStepSave  = oneDFCount;

                hyFloat         lastLogL = maxSoFar;

                if (useAdaptiveStep>0.5) {
                    if (convergenceMode < 2) {
                        LocateTheBump (current_index,precisionStep, maxSoFar, bestVal, go2Bound, precisionStep);
                    } else {
                        LocateTheBump (current_index,precisionStep, maxSoFar, bestVal, go2Bound, convergenceMode == 2? precisionStep*0.25: precisionStep*0.0625);
                    }
                } else {
                    LocateTheBump (current_index,brackStep, maxSoFar, bestVal, go2Bound);
                }

                hyFloat  cj = GetIthIndependent(current_index),
                            ch = fabs(bestVal-cj);



                if (useAdaptiveStep>0.5) {
                    if (cj != 0.) {
                        averageChange += fabs (ch/cj);
                    }
                    if ((ch < precisionStep*0.01 /*|| lastLogL - maxSoFar < 1e-6*/) && inCount == 0) {
                        nc2 << current_index;
                    }
                } else {
                    averageChange  += ch;
                    averageChange2 += cj;
                    if (ch<currentPrecision/indexInd.lLength) {
                        nc2 << current_index;
                    }
                }


                if (vH) {
                    vH->Store(cj);
                }

                variableValues[current_index] = ch;

                if (verbosity_level>1) {
                    snprintf (buffer, sizeof(buffer),"\nindex = %ld\tlog(L) = %14.10g\t param value = %10.6g ( diff = %10.6g, bracket = %10.6g, precision %10.6g) EVALS: %ld (BRACKET), %ld (BRENT) ",
                             current_index,
                             maxSoFar,
                             cj,
                             ch,
                             brackStep,
                             precisionStep,
                             bracketFCount-brackStepSave,
                             oneDFCount - oneDStepSave
                             );
                    BufferToConsole (buffer);
                    StringToConsole (*GetIthIndependentVar(current_index)->GetName());
                    BufferToConsole (CheckEqual(GetIthIndependentBound (current_index, true), cj)? ("[Lower bound]") : (CheckEqual(GetIthIndependentBound (current_index, false),cj) ? "[Upper bound]" : ""));
                }
#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
                else if (verbosity_level==1) {
                    UpdateOptimizationStatus (maxSoFar,-1,1,true,progressFileString);
                }
#endif

            }

            if (!forward || doShuffle > 0.1) {
                nc2.Sort();
            }

            averageChange2/=indexInd.lLength;
            // mean of parameter values
            averageChange/=(hyFloat)(indexInd.lLength-nc2.lLength+1);
            // mean of change in parameter values during the last step

            if (glVars.lLength == indexInd.lLength) {
                noChange.Clear();
                noChange.Duplicate (&nc2);
            } else {
                noChange.Subtract  (nc2,glVars);
            }

            if (noChange.lLength==0) {
                noChange << -1;
            }



            forward = true;
            percentDone = progress_tracker.PushValue (maxSoFar);
            
            if (verbosity_level==1) {
                UpdateOptimizationStatus (maxSoFar,percentDone,1,true,progressFileString);
            }

            logLHistory.Store(maxSoFar);
            loopCounter += 1.;

            if (verbosity_level>5) {
                snprintf (buffer, sizeof(buffer),"\nAverage Variable Change: %g, percent done: %g, divFactor: %g, oldAverage/averageChange: %g", averageChange, percentDone,divFactor,oldAverage/averageChange);
                BufferToConsole (buffer);
                snprintf (buffer, sizeof(buffer),"\nDiff: %g, Precision: %16.12g, termFactor: %ld", maxSoFar-lastMaxValue, precision, termFactor);
                BufferToConsole (buffer);
                snprintf (buffer, sizeof(buffer),"\nSmoothing term: %g", smoothingTerm);
                BufferToConsole (buffer);

            }

            currentPrecision = averageChange/divFactor;

            if(currentPrecision<1e-6) {
                currentPrecision = 1e-6;
            } else if (currentPrecision>averageChange2) {
                currentPrecision = averageChange2;
            }
            
            if (custom_convergence_callback >= 0) {
                _List arguments;
                _AssociativeList * parameters = new _AssociativeList;
                
                GetIndependentVars().Each ([&] (long value, unsigned long i) -> void {
                    parameters-> MStore (*GetIthIndependentName(i), new _Constant (GetIthIndependent(i, false)));
                });
                arguments < new _Constant (maxSoFar) < parameters;

                HBLObjectRef convegence_check = custom_convergence_callback_name->Call (&arguments, nil);
                
                if (convegence_check->Value () <= precision/termFactor) {
                    inCount ++;
                }
                else {
                    inCount = 0;
                }
                DeleteObject (convegence_check);
            } else {
                if (maxSoFar-lastMaxValue<=precision/termFactor) {
                    inCount++;
                } else {
                    inCount = 0;
                }
            }

            lastMaxValue = maxSoFar;

            if (useAdaptiveStep < 0.5) {
                  if (!skipCG && loopCounter&& indexInd.lLength>1 && ( (((long)loopCounter)%indexInd.lLength)==0 )) {
                      _Matrix             bestMSoFar;
                      GetAllIndependent   (bestMSoFar);
                      maxSoFar = ConjugateGradientDescent (currentPrecision, bestMSoFar);
                      logLHistory.Store(maxSoFar);
                  }
            }

            if (hardLimitOnOptimizationValue < INFINITY && timer.TimeSinceStart() > hardLimitOnOptimizationValue) {
                ReportWarning (_String("Optimization terminated before convergence because the hard time limit was exceeded."));
                break;
            }

        }

        ReportWarning (_String("Optimization finished in ") & loopCounter & " loop passes.\n" & _String ((long)likeFuncEvalCallCount-evalsIn) & " likelihood evaluation calls and " & _String ((long)matrix_exp_count - exponentiationsIn) & " matrix exponentiations calls were made\n");

        if (optimization_mode == kOptimizationGradientDescent) {
            _Matrix bestMSoFar (indexInd.lLength,1,false,true);
            GetAllIndependent(bestMSoFar);
            maxSoFar = ConjugateGradientDescent (currentPrecision*.01, bestMSoFar);
        }
        DeleteObject (stepHistory);

    } else if (optimization_mode== kOptimizationNedlerMead) {
        SimplexMethod (precision, get_optimization_setting (kMaximumIterations, 10000000), maxItersPerVar);
    }
    
    if (keepOptimizationLog) {
        ((_Vector*)optimizatonHistory->GetByKey("LogL"))->Trim();
        _AssociativeList* variable_traces = ((_AssociativeList*)optimizatonHistory->GetByKey("Parameters"));
        
        for (unsigned long var_id = 0; var_id < indexInd.lLength; var_id++) {
            ((_Vector*)variable_traces->GetByKey(*GetIthIndependentName(var_id)))->Trim();
        }
        CheckReceptacleAndStore(AppendContainerName("trace", GetObjectNameByType(HY_BL_LIKELIHOOD_FUNCTION,lockedLFID, false)), "", false, optimizatonHistory, true);
        DeleteObject (optimizatonHistory);
        optimizatonHistory = nil;
    }

    _Matrix result (2,indexInd.lLength+indexDep.lLength<3?3:indexInd.lLength+indexDep.lLength, false, true);

    //forceRecomputation = true;
    result.Store (1,0,Compute());
    //forceRecomputation = false;
    result.Store (1,1,indexInd.lLength);
    result.Store (1,2,CountObjects(kLFCountGlobalVariables));



  #if !defined __UNIX__ || defined __HEADLESS__
    SetStatusBarValue (-1,maxSoFar,(likeFuncEvalCallCount-evalsIn)/TimerDifferenceFunction(true));
  #endif

    CleanUpOptimize();

    HBLObjectRef pm;
    for (unsigned long i=0UL; i<indexInd.lLength; i++) {
        result.Store(0,i,GetIthIndependent(i));
    }
    for (unsigned long i=0UL; i<indexDep.lLength; i++) {
        result.Store(0,i+indexInd.lLength,GetIthDependent(i));
    }

    return (_Matrix*)result.makeDynamic();
}
//_______________________________________________________________________________________
    
void    _LikelihoodFunction::_TerminateAndDump(const _String &error, bool sig_term) {
  
    FILE * out = doFileOpen ("/tmp/hyphy.dump", "w");
  
    _String err ("Internal error ");
  
    if (out) {
       this->DoneComputing();
      _StringBuffer          sLF (8192L);
      SerializeLF      (sLF,_hyphyLFSerializeModeVanilla);
      sLF.TrimSpace();
      fwrite ((void*)sLF.get_str(), 1, sLF.length(), out);
      fclose (out);
    }

    HandleApplicationError (err & '\n' & error, true);

}
//_______________________________________________________________________________________

void _LikelihoodFunction::CleanUpOptimize (void) {
    categID = 0;
    CleanupParameterMapping ();
    //printf ("Done OPT LF eval %d MEXP %d\n", likeFuncEvalCallCount, matrix_exp_count);
#ifdef __HYPHYMPI__
    if (hyphyMPIOptimizerMode==_hyphyLFMPIModeNone) {
#endif
        for (long i=0; i<theTrees.lLength; i++) {
            _TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
            cT->CleanUpMatrices();
        }

        DeleteCaches (false);

        if (mstCache) {
            hyFloat      umst = 0.0;
            checkParameter (useFullMST,umst,0.0);
            if (umst>.5)
                for (long kk=0; kk<mstCache->cacheSize.lLength; kk++) {
                    long cS = mstCache->cacheSize.lData[kk];
                    if ((cS>0)&&(mstCache->resultCache[kk])) {
                        hyFloat ** c1 = (hyFloat**)mstCache->resultCache[kk];
                        for (long k2 = 0; k2<cS; k2++) {
                            delete c1[k2];
                        }
                        delete c1;

                        long ** c2 = (long**)mstCache->statesCache[kk];
                        for (long k3 = 0; k3<cS; k3++) {
                            delete c2[k3];
                        }
                        delete c2;

                        char ** c3 = (char**)mstCache->statesNCache[kk];
                        for (long k4 = 0; k4<cS; k4++) {
                            delete c3[k4];
                        }
                        delete c3;

                        ((_SimpleList*)leafSkips(kk))->Clear();
                        ((_SimpleList*)leafSkips(kk))->Duplicate (mstCache->stashedLeafOrders(kk));
                        //printf ("\n%s\n", ((_String*)((_SimpleList*)leafSkips(kk))->toStr())->getStr());
                    }
                }
            mstCache->resultCache.Clear();
            mstCache->statesCache.Clear();
            mstCache->statesNCache.Clear();
            mstCache->stashedLeafOrders.Clear();
        }
#ifdef __HYPHYMPI__
    } else {
        CleanupMPIOptimizer();
    }

#endif


#if defined __UNIX__ && ! defined __HEADLESS__ && !defined __HYPHYQT__ && !defined __HYPHY_GTK__
    if (verbosity_level==1) {
        UpdateOptimizationStatus (0,0,2,true,progressFileString);
    }
#endif
    //printf ("\n%d\n",matrix_exp_count);
    /*printf ("\nOptimization spool:\n");
     for (i=0; i<indexInd.lLength; i++)
     {
     _Variable* vspool = LocateVar (indexInd(i));
     printf ("%d %s %g\n", i, vspool->GetName()->getStr(), vspool->Compute()->Value());
     }*/

    setParameter (likeFuncCountVar,likeFuncEvalCallCount);
    isInOptimize = false;
    DoneComputing();
    hasBeenOptimized = true;
    hasBeenSetUp     = 0;
    lockedLFID       = -1;
    DeleteObject     (nonConstantDep);
    nonConstantDep = nil;
}

//_______________________________________________________________________________________
bool CheckEqual (hyFloat a, hyFloat b, hyFloat tolerance) {
    if (a!=0.0) {
        a = (a>b)?(a-b)/a:(b-a)/a;
        return a>0.0 ? a<=tolerance : a>=-tolerance;
    }
    return (b<=tolerance)&&(b>=-tolerance);
}
    
//_______________________________________________________________________________________
bool CheckRange (hyFloat value, hyFloat lb, hyFloat ub, bool exceptions) {
    if (value >= lb && value <= ub) {
        return true;
    }
    if (exceptions) {
        throw (_String ("Value ") & value & " is outside the expected interval [" & lb &
               "-" & ub & "]");
    }
    return false;
}
    
//_______________________________________________________________________________________
bool CheckArgumentType (HBLObjectRef object, long type, bool exceptions) {
    
    if (object->ObjectClass() & type) {
        return true;
    }
    
    if (exceptions) {
        throw (_String ("Unsupported agrument type ") & FetchObjectNameFromType (object->ObjectClass()).Enquote());
    }
    return false;
}

//_______________________________________________________________________________________
hyFloat _LikelihoodFunction::SetParametersAndCompute (long index, hyFloat value, _Matrix* baseLine, _Matrix* direction, bool skip_compute)
{
    if (index >= 0) {
        SetIthIndependent (index,value);
    } else {
        /*if (value < 0.) {
            _TerminateAndDump (_String ("Negative value for parameter ") & *GetIthIndependentName(index));
            return -INFINITY;
        }*/
        _Matrix newValue (*baseLine);
        newValue.AplusBx (*direction, value);
        SetAllIndependent (&newValue);

    }
    if (skip_compute) {
        return -INFINITY;
    }
    hyFloat logL = Compute();
    //if (index >=0)
    //  printf ("[SetParametersAndCompute %g = %g]\n", value, logL);

    return logL;

    //return Compute();
}

//_______________________________________________________________________________________

void     _LikelihoodFunction::UnregisterListeners (void) {
    unsigned long partition_count = CountObjects(kLFCountPartitions);
    for (unsigned long i = 0UL; i < partition_count; i++) {
        UnregisterChangeListenerForDataFilter(theDataFilters.GetElement(i), this);
    }
}
    
//_______________________________________________________________________________________
void _LikelihoodFunction::GetAllIndependent (_Matrix & storage) const {
    storage.Clear();
    _Matrix::CreateMatrix (&storage, indexInd.lLength,1,false,true, false);
    for (unsigned long i=0UL; i<indexInd.lLength; i++) {
        storage.theData[i]=GetIthIndependent(i);
    }

}


//_______________________________________________________________________________________
long    _LikelihoodFunction::Bracket (long index, hyFloat& left, hyFloat& middle, hyFloat& right,  hyFloat& leftValue, hyFloat& middleValue, hyFloat& rightValue, hyFloat& initialStep, bool go2Bound, _Matrix* gradient)
{
    _Variable* curVar = index >= 0 ? GetIthIndependentVar (index) : nil;

    bool       movingLeft = false,
               first      = true,
               successful = false;

    hyFloat lowerBound = curVar?GetIthIndependentBound(index,true):0.,
               upperBound = curVar?GetIthIndependentBound(index,false):0.,
               practicalUB,
               magR = 2.,//1.61803,
               //r,q,u,d,
               leftStep  = initialStep*.5,
               rightStep = initialStep*.5,
               saveL     = index<0?middle:NAN,
               saveM     = index<0?NAN:middle,
               saveR     = NAN,
               saveLV    = index<0?middleValue:0.0,
               saveMV    = index<0?0.0:middleValue,
               saveRV    = 0.0,
               stash_middle;



    _Matrix            currentValues;
    if (index < 0) {
        GetAllIndependent                      (currentValues);
        GetGradientStepBound                   (*gradient, lowerBound, upperBound);
        if (upperBound < 1e-10) {
            long freezeCount = 0;
            GetGradientStepBound                   (*gradient, lowerBound, upperBound, &freezeCount);
            //printf ("[FREEZE %ld/%g]\n", freezeCount, upperBound);
            if (freezeCount == 0 || freezeCount == indexInd.lLength || upperBound < 1e-10) {
                return -2;
            }
        }
        lowerBound                             = 0.;
    }

    practicalUB = upperBound>DEFAULTPARAMETERUBOUND?DEFAULTPARAMETERUBOUND:upperBound;
    long               funcCounts = likeFuncEvalCallCount;


    if (index >= 0) {
        middle  =  GetIthIndependent (index);
    } else {
        middle  = initialStep;
    }

    if (lowerBound>middle || upperBound<middle) {
        middle = (lowerBound+practicalUB) * .5;
    }

    if (CheckEqual(middle,lowerBound)) {
        leftStep = initialStep * .2;
        middle   = lowerBound+leftStep;
        saveL    = lowerBound;
        saveLV = middleValue;
    }

    if (middle == upperBound) {
        rightStep = initialStep * .2;
        middle    = upperBound-rightStep;
        saveR     = upperBound;
        saveRV = middleValue;
    }

    if (index < 0) {
        leftStep = middle;
    }

    if (saveM != middle) {
        saveM = NAN;
    }



    /*if (index < 0)
    {
        printf                                 ("[Bracket bounds %g - %g (%g)/%g]\n", lowerBound, upperBound, practicalUB, middle);
        for (unsigned long i = 0; i < indexInd.lLength; i++) {
            printf ("%s = %.16g\n", GetIthIndependentVar(i)->GetName()->sData, gradient->theData[i]);
        }
    }*/

    if (verbosity_level > 100) {
        char buf [512];
        snprintf (buf, sizeof(buf), "\n\t[_LikelihoodFunction::Bracket (index %ld, eval %ld) INITIAL BRACKET %15.12g <= %15.12g (current %15.12g) <= %15.12g]", index, likeFuncEvalCallCount, middle-leftStep, middle, index>=0?GetIthIndependent (index):0.0, middle+rightStep);
        BufferToConsole (buf);
    }

    /*
    if (likeFuncEvalCallCount > 0) {
      printf ("\n\n\nCHECK INDEX 6\n\n\n");
      SetIthIndependent(6L, GetIthIndependent(6L));
    }
    */

    while (1) {

        while (middle-leftStep < lowerBound) {
            if (verbosity_level > 100) {
              char buf [512];
              snprintf (buf, sizeof(buf), "\n\t[_LikelihoodFunction::Bracket (index %ld) HANDLING LEFT BOUNDARY CASES] : LB = %g, current try = %g, current evaluated midpoint value = %g (%s)", index, lowerBound, middle-leftStep, middleValue, first ? "first" : "NOT first");
              BufferToConsole (buf);
            }

            leftStep*=.125;
            if ( leftStep<initialStep*.1 && index >0 || index < 0 && leftStep < STD_GRAD_STEP) {
                if (!first) {
                    if (go2Bound) {
                        middle = lowerBound==0.0 ? PERTURBATION_OF_ZERO : lowerBound;
                        if (verbosity_level > 100) {
                          char buf [512];
                          snprintf (buf, sizeof(buf), "\n\t[_LikelihoodFunction::Bracket LEFT BOUNDARY (index %ld) UPDATED middle from %15.12g to %15.12g, LogL = %15.12g]", index, middle, middleValue);
                          BufferToConsole (buf);
                        }
                        middleValue = stash_middle = SetParametersAndCompute (index, middle, &currentValues, gradient);
                        
                    }
                    return -2;
                } else {
                    middle=MIN(lowerBound+initialStep*.1,upperBound-rightStep);
                    first = false;
                }
            }
        }


        while (rightStep+middle > upperBound) {
            rightStep*=.125;
            if (rightStep<initialStep*.1 && index >0 || index < 0 && rightStep < STD_GRAD_STEP) {
                if (!first) {
                    if (go2Bound) {
                        middleValue = stash_middle = SetParametersAndCompute (index, middle=upperBound, &currentValues, gradient);
                    }
                    return -2;
                } else {
                    middle=MAX(upperBound-initialStep*.1,lowerBound+leftStep);
                    first = false;
                }
            }

        }



        if (CheckEqual(middle,saveL)) {
           stash_middle = middleValue = saveLV;
        } else if (CheckEqual(middle,saveR)) {
            stash_middle = middleValue = saveRV;
        } else if (CheckEqual(middle,saveM)) {
            stash_middle = middleValue = saveMV;
        } else {
            middleValue = stash_middle = SetParametersAndCompute (index, middle, &currentValues, gradient);
            if (verbosity_level > 100) {
               char buf [512];
               snprintf (buf, sizeof(buf), "\n\t[_LikelihoodFunction::Bracket (index %ld) UPDATED middle to %15.12g, LogL = %15.12g]", index, middle, middleValue);
               BufferToConsole (buf);
            }
        }

        left = middle-leftStep;

        if (CheckEqual(left,saveL)) {
            leftValue = saveLV;
        } else if (CheckEqual(left,saveR)) {
            leftValue = saveRV;
        } else if (CheckEqual(left,saveM)) {
            leftValue = saveMV;
        } else {
            leftValue = SetParametersAndCompute (index, left, &currentValues, gradient);
            if (verbosity_level > 100) {
                char buf [512];
                snprintf (buf, sizeof(buf), "\n\t[_LikelihoodFunction::Bracket (index %ld) UPDATED left to %15.12g, LogL = %15.12g]", index, left, leftValue);
                BufferToConsole (buf);
            }
        }


        right = middle+rightStep;
        if (CheckEqual(right,saveL)) {
            rightValue = saveLV;
        } else if (CheckEqual(right,saveR)) {
            rightValue = saveRV;
        } else if (CheckEqual(right,saveM)) {
            rightValue = saveMV;
        } else {
            rightValue = SetParametersAndCompute (index, right, &currentValues, gradient);
            if (verbosity_level > 100) {
                char buf [512];
                snprintf (buf, sizeof(buf), "\n\t[_LikelihoodFunction::Bracket (index %ld) UPDATED right to %15.12g, LogL = %15.12g]", index, right, rightValue);
                BufferToConsole (buf);
            }
        }

        if (verbosity_level > 100) {
            char buf [512];
            snprintf (buf, 512, "\n\t[_LikelihoodFunction::Bracket (index %ld): BRACKET %g (LogL : %15.12g diff: %15.12g) - %g (logL: %15.12g) - %g (LogL : %15.12g diff: %15.12g)]", index, left, leftValue, leftValue-middleValue, middle, middleValue, right, rightValue, rightValue-middleValue);
            BufferToConsole (buf);
        }

        saveL       = left;
        saveLV      = leftValue;

        saveM       = middle;
        saveMV      = middleValue;

        saveR       = right;
        saveRV      = rightValue;


        if (rightValue<=middleValue && leftValue<=middleValue) {
            successful = true;
            break;
        }

        // case 1: /
        if (rightValue>=middleValue && middleValue>=leftValue) {
            if (movingLeft) {
                rightStep/=magR;
            } else {
                leftStep = rightStep;
                rightStep*=magR;
            }
            middle     = right;
            movingLeft = false;
        } else // case 2
            if (rightValue<=middleValue && middleValue<=leftValue) {
                if (!movingLeft && !first) {
                    leftStep/=magR;
                } else {
                    rightStep =  leftStep;
                    leftStep  *= magR;
                }
                if (index < 0) {
                    if (CheckEqual (left, lowerBound)) {
                        middle    = (middle-lowerBound)*0.5;
                        leftStep  = middle;
                        rightStep = middle;
                    } else {
                        middle     = left;
                    }
                } else {
                    middle     = left;
                }
                movingLeft = true;
            } else {
                if (movingLeft) {
                    middle = left;
                } else {
                    middle = right;
                }
            }

        if (middle>=practicalUB) {
            stash_middle = middleValue         = SetParametersAndCompute (index, middle = practicalUB, &currentValues, gradient);
            //if (index < 0) printf ("\nmiddle>=practicalUB\n");
            break;
        }
        /*
        if (middle-lowerBound < STD_GRAD_STEP*0.5) {
            middleValue         = SetParametersAndCompute (index, middle = lowerBound+STD_GRAD_STEP*0.5, &currentValues, gradient);
            //if (index < 0) printf ("\nmiddle-lowerBound < STD_GRAD_STEP*0.5\n");
            break;
        }
        */
        first = false;

    }

    if (curVar) {
        if (CheckAndSetIthIndependent(index,middle)) {
           middleValue = stash_middle; // Compute ()

           if (verbosity_level > 100) {
              char buf [256];
              snprintf (buf, 256, "\n\t[_LikelihoodFunction::Bracket (index %ld) recomputed the value to midpoint: L(%20.16g) = %20.16g [@%20.16g -> %20.16g:@%20.16g -> %20.16g]]", index, middle, middleValue, left, leftValue,right, rightValue);
               BufferToConsole (buf);
            }
         }
    } else {
        middleValue         = stash_middle;
        SetParametersAndCompute (index, middle, &currentValues, gradient, false);
    }

    if (successful ) {
        /** SLKP 20180709 need to have a more permissive check, because sometimes if the change is too small
                  (or involves a paremeter that has very little effect on the LF), recomputation could be within numerical error
         
        **/
        if (rightValue - middleValue > 1e-9 || leftValue - middleValue > 1e-9) {
         char buf[256], buf2[512];
         snprintf (buf, 256, " \n\tERROR: [_LikelihoodFunction::Bracket (index %ld) recomputed the value to midpoint: L(%g) = %g [@%g -> %g:@%g -> %g]]", index, middle, middleValue, left, leftValue,right, rightValue);
         snprintf (buf2, 512, "\n\t[_LikelihoodFunction::Bracket (index %ld) BRACKET %s: %20.16g <= %20.16g >= %20.16g. steps, L=%g, R=%g, values %15.12g : %15.12g - %15.12g]", index, successful ? "SUCCESSFUL" : "FAILED", left,middle,right, leftStep, rightStep, leftValue - middleValue, middleValue, rightValue - middleValue);
         _TerminateAndDump (_String (buf) & "\n" & buf2 &  "\nParameter name " & (index >= 0 ? *GetIthIndependentName(index) : "line optimization"));
        }
        successful = rightValue<=middleValue && leftValue<=middleValue;
    }
    
    if (verbosity_level > 100) {
        char buf [256];
        snprintf (buf, 256, "\n\t[_LikelihoodFunction::Bracket (index %ld) BRACKET %s: %15.12g <= %15.12g <= %15.12g. Steps, L=%g, R=%g, values (delta-left) %20.16g : %20.16g : (delta-right) %20.16g]", index, successful ? "SUCCESSFUL" : "FAILED", left,middle,right, leftStep, rightStep, middleValue-leftValue, middleValue, middleValue-rightValue);
        BufferToConsole (buf);
    }


    bracketFCount+=likeFuncEvalCallCount-funcCounts;
    bracketCount++;
    return successful ? 0 : -1;
}
//_______________________________________________________________________________________

void    _LikelihoodFunction::CheckStep (hyFloat& tryStep, _Matrix vect, _Matrix* selection) {
    for (unsigned long index = 0; index<indexInd.lLength; index++) {

        hyFloat  Bound,
                    currentValue,
                    locValue = vect.theData[index];

        if (fabs(locValue)<1e-14) {
            Bound = GetIthIndependentBound(index,false);
            locValue = 0.0;
        } else if (locValue<0) {
            Bound = GetIthIndependentBound(index,true);
        } else {
            Bound = GetIthIndependentBound(index,false);
        }


        if (selection) {
            currentValue = selection->theData[index];
        } else {
            currentValue = GetIthIndependent (index);
        }

        if (Bound>1000.) {
            Bound = 1000.;
        }
        if (locValue>=0) {
            while (currentValue+locValue*tryStep>Bound-1e-8) {
                tryStep/=5;
                if (tryStep<1e-8) {
                    tryStep = 0.;
                    return;
                }
            }
        } else {
            while (currentValue+locValue*tryStep<Bound+1e-8) {
                tryStep/=5;
                if (tryStep<1e-8) {
                    tryStep = 0.;
                    return;
                }
            }
        }
    }
}

//_______________________________________________________________________________________

HBLObjectRef   _LikelihoodFunction::CovarianceMatrix (_SimpleList* parameterList) {
    if (indexInd.empty()) {
        return nil;
    }

    hyFloat  h = 1.e-5,//STD_GRAD_STEP,
                functionValue,
                t1,
                t2,
                cm; // small h and L(x_opt)

    checkParameter (covariancePrecision,cm,1.0);

    long        parameter_count = parameterList?parameterList->countitems():indexInd.countitems(),
                j;

    PrepareToCompute();

    functionValue = Compute();
    _Variable      *thisVar;

    bool        useIndirectIndexing = false;

    if (parameterList) {
        useIndirectIndexing = true;
    } else {
        parameterList = &indexInd;
    }

    if (cm<1.) {
        // use likelihood profile with the appropriate signifcance level
        _Matrix     * sigLevels = new _Matrix (parameter_count,7,false,true);

        // find the appropriate significance level for chi2

        _String xxc = _HYGenerateANameSpace () & (".x");
        thisVar = CheckReceptacle (&xxc, kEmptyString);
        thisVar->SetBounds (-1.e30,1.e30);

        try {
            if (cm>0.0) {
                _String     fString = _String ("CChi2(") & xxc & ",1)-" & cm;
                _Formula    CChi2Fla (fString,nil);
                t1 = CChi2Fla.Brent (thisVar,0,0);
                if (fabs(CChi2Fla.Compute()->Value())>1.e-6) {
                    throw _String ("Failed to compute chi-square significance level");
                }
            } else {
                if (cm<0.0) {
                    t1 = -2.*cm;
                } else {
                    throw (_String("Must have a non-zero ") & covariancePrecision);
                }
            }
        } catch (const _String& e) {
            DoneComputing();
            DeleteVariable(thisVar->get_index(), true, false);
            HandleApplicationError (e);
            return    sigLevels;
        }

        thisVar->SetNumericValue (t1);

        /* ____________________________________ NEW CODE BY AFYP _____________________________________ */
        long        mastodon    = likeFuncList._SimpleList::Find((long)this);
        _String *   myName;

        if (mastodon < 0) {
            mastodon    = scfgList._SimpleList::Find((long)this);
            myName      = (_String *) scfgNamesList (mastodon);

        } else {
            // generate HBL
            myName      = (_String*)likeFuncNamesList (mastodon);
        }

        _String xxv = "_xxv";
        
        _String     fString     = _String("function _profileFit(") & xxv & ",_variableIndex){SetParameter(" &*myName&",_variableIndex," &xxv&");LFCompute("
                                  //    &*myName&(",_xxres);  fprintf (stdout,\"\\n\",_xxv_,\" \",_xxres);  return _xxres;}");
                                  &*myName & _String (",_xxres);") //"fprintf (stdout,\"\\n\",") & xxv & ",\" \",_xxres);
                                  &"return _xxres;}";

        //printf ("%s\n", fString.get_str());

        /* ___________________________________ ! NEW CODE BY AFYP ____________________________________ */


        t1 = functionValue-.5*t1;

        _ExecutionList exL (fString);
        exL.Execute ();

        for (parameter_count=0; parameter_count<parameterList->lLength; parameter_count++) {
            j = useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count;
            t2 = GetIthIndependent (j);
            _Variable* thisVar2 = LocateVar (indexInd.lData[j]);
            thisVar->SetBounds (thisVar2->GetLowerBound(),thisVar2->GetUpperBound());
            sigLevels->Store (parameter_count,1,t2);

            char buffer[255];
            snprintf (buffer, sizeof(buffer),"%.14g",t1);

            fString = _String("_profileFit(") & xxc & "," & j & ")-(" & buffer& ')';
            _Formula    FitFla (fString,nil);
            
            if (CheckEqual(t2,thisVar2->GetLowerBound())) {
                sigLevels->Store (parameter_count,0,t2);
            } else {
                h = FitFla.Brent (thisVar,t2+1,t2,t2*0.0001+0.000001);
                //sigLevels.Store (i,0,MAX(h,thisVar2->GetLowerBound()));
                 if (h <= thisVar2->GetLowerBound()) {
                  _String lf_buffer ("_xxres");
                  snprintf (buffer, sizeof(buffer),"%.14g",FetchVar(LocateVarByName(lf_buffer))->Value() + 0.0001);
                  lf_buffer = _String("_profileFit(_xx_,") & j & ")-(" & buffer& ')';
                  _Formula try_again (lf_buffer,nil);
                  h = try_again.Brent (thisVar,t2+1,t2,t2*0.0001+0.000001);
                  sigLevels->Store (parameter_count,5,h);
                  sigLevels->Store (parameter_count,0,thisVar2->GetLowerBound());

                } else {
                  sigLevels->Store (parameter_count,0,h);
                  sigLevels->Store (parameter_count,5,h);
                }
           }

            snprintf (buffer, sizeof(buffer),"%.14g",(*sigLevels) (parameter_count,0));
            _String checkLFDIFF = _String("CChi2(2*(-_profileFit(") & buffer & "," & j & ")+(" & functionValue & ")),1)";
            HBLObjectRef lf_diff = (HBLObjectRef) _FString (checkLFDIFF, false).Evaluate(_hyDefaultExecutionContext);
            sigLevels->Store (parameter_count,3,lf_diff->Value());
            DeleteObject (lf_diff);


            if (CheckEqual(t2,thisVar2->GetUpperBound())) {
                sigLevels->Store (parameter_count,2,t2);
            } else {
                //_List store_evals;
                h = FitFla.Brent (thisVar,t2,t2,t2*0.0001+0.000001);//, &store_evals);
                if (h >= thisVar2->GetUpperBound()) {
                  _String lf_buffer ("_xxres");
                  snprintf (buffer, sizeof(buffer),"%.14g",FetchVar(LocateVarByName(lf_buffer))->Value() + 0.0001);
                  lf_buffer = _String("_profileFit(_xx_,") & j & ")-(" & buffer& ')';
                  _Formula try_again (lf_buffer,nil);
                  h = try_again.Brent (thisVar, t2,t2,t2*0.0001+0.000001);
                  sigLevels->Store (parameter_count,6,h);
                  sigLevels->Store (parameter_count,2,thisVar2->GetUpperBound());

                } else {
                  sigLevels->Store (parameter_count,2,h);
                  sigLevels->Store (parameter_count,6,h);
                }
            }

             snprintf (buffer, sizeof(buffer),"%.14g",(*sigLevels) (parameter_count,2));
             checkLFDIFF =_String("CChi2(2*(-_profileFit(") & buffer & "," & j & ")+(" & functionValue & ")),1)";
             lf_diff = (HBLObjectRef) _FString (checkLFDIFF, false).Evaluate(_hyDefaultExecutionContext);
             sigLevels->Store (parameter_count,4,lf_diff->Value());
             DeleteObject (lf_diff);

            SetIthIndependent (j,t2);
        }

        DoneComputing();
        DeleteVariable(thisVar->get_index(), true, false);
        return  sigLevels;
    }

    _Matrix     hessian    (parameter_count,parameter_count,false,true),
                funcValues (parameter_count,5,false,true),
                *iMap = nil;

    _AssociativeList * mapMethod = nil;

    hyFloat      uim = 0.0;
    checkParameter  (useIntervalMapping, uim, 0.0);

    if (uim > 0.5) {
        iMap = new _Matrix (parameter_count,3,false,true);
        mapMethod = new _AssociativeList;

    }
    // y,x',x''
    // first check for boundary values and move the parameter values a bit if needed
    for (parameter_count=0; parameter_count<parameterList->lLength; parameter_count++) {
        long     dIndex = useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count;
        thisVar = LocateVar (indexInd.lData[dIndex]);
        t1 = thisVar->Value();
        funcValues.Store (parameter_count,3,t1);

        hyFloat locH = 1./131072.; //*(t1>10.?exp(log(10.)*((long)log(t1)/log(10.))):1.);

        if (locH<1e-7) {
            locH = 1e-7;
        }

        /*hyFloat locH = 1./1024.,
                   tryH = locH * 0.25,
                   dv1,
                   dv2,
                   lastErr = 1e100;

        SetIthIndependent (dIndex, t1+locH);
        dv1 = Compute();
        SetIthIndependent (dIndex, t1-locH);
        dv1 += Compute()-2.*functionValue;
        dv1 /= 4.*locH*locH;

        while (tryH >= 1e-12)
        {
            SetIthIndependent (dIndex, t1+tryH);
            dv2 = Compute ();
            SetIthIndependent (dIndex, t1-tryH);
            dv2 += Compute ()-2.*functionValue;
            dv2 /= 4.*tryH*tryH;
            hyFloat err = fabs (dv2-dv1) / MAX (fabs(dv2), fabs(dv1));
            if (err > lastErr)
                break;
            else
                lastErr = err;
            dv1 = dv2;
            locH = tryH;
            tryH *= 0.25;
        }*/


        funcValues.Store (parameter_count,4,(locH+t1)-t1);

        if (t1+locH > thisVar->GetUpperBound()) {
            SetIthIndependent (dIndex,thisVar->GetUpperBound()-2.0*locH);
        } else if (t1-locH < thisVar->GetLowerBound()) {
            SetIthIndependent (dIndex,thisVar->GetLowerBound()+2.0*locH);
        }

        if (uim > 0.5) {
            hyFloat      lb  = thisVar->GetLowerBound(),
                            ub  = thisVar->GetUpperBound(),
                            y ,
                            dy2,
                            dy;

            _FString        varKeyName (*thisVar->GetName());
            _Constant       varMapMethod;

            if (CheckEqual (lb, DEFAULTPARAMETERLBOUND) && CheckEqual (ub, DEFAULTPARAMETERUBOUND))
                // use exp<->log map
            {
                y = log (t1-lb);
                dy = exp(y);
                dy2 = dy;
                varMapMethod.SetValue(1.);
            }
            // use logistic map
            else {
                y   = (t1-lb)/(ub-lb);
                t1  = y/(1-y);
                y   = log (t1);
                dy  = (ub-lb)*(t1/(t1+1)/(t1+1));
                dy2 = dy*((1-t1*t1)/(t1+1)/(t1+1));
                varMapMethod.SetValue(2.);
            }
            iMap->Store(parameter_count,0,y);
            iMap->Store(parameter_count,1,dy);
            iMap->Store(parameter_count,2,dy2);
            mapMethod->MStore (&varKeyName, &varMapMethod, true);
        }
    }



#if defined __MAC__ || defined __WINDOZE__ || defined __HYPHYQT__ || defined __HYPHY_GTK__
    hyFloat totalCount    = 2*parameterList->lLength;

    long       finishedCount = 0;

    if (cm>1.1) {
        totalCount += 2*parameterList->lLength*(parameterList->lLength-1);
    } else {
        totalCount += (parameterList->lLength*(parameterList->lLength-1)/2);
    }

    BenchmarkThreads(this);
#endif

    // fill in funcValues with L(...,x_i\pm h,...) and 1st derivatives and get 2nd derivatives
    for (parameter_count=0; parameter_count<parameterList->lLength; parameter_count++) {
        long              pIdx = useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count;

        hyFloat        pVal = GetIthIndependent (pIdx),
                          d1,
                          locH = funcValues (parameter_count,4);

        SetIthIndependent (pIdx,pVal-locH); // - step
        t1 = Compute();
        funcValues.Store (parameter_count,0,t1);
        SetIthIndependent (pIdx,pVal+locH); // + step
        t2 = Compute();
        funcValues.Store (parameter_count,1,t2);          // reset value
        SetIthIndependent (pIdx,pVal);
        d1 = (t2-t1)/(2.0*locH);
        // central 1st derivative
        funcValues.Store (parameter_count,2,d1);

        t1  = ((t1-functionValue)+(t2-functionValue))/(locH*locH);
        // Standard central second derivative

        if (uim < 0.5) {
            hessian.Store (parameter_count,parameter_count,-t1);
        } else {
            hessian.Store (parameter_count,parameter_count,-(t1*(*iMap)(parameter_count,1)*(*iMap)(parameter_count,1)+(*iMap)(parameter_count,2)*d1));
        }

#ifndef __UNIX__
        finishedCount += 2;
        if (TimerDifferenceFunction(true)>1.) {
            SetStatusBarValue (finishedCount/totalCount*100.,1,0);
            TimerDifferenceFunction (false);
        }
#endif
    }


    if (cm>1.1) {
        // fill in off-diagonal elements using the f-la
        // f_xy = 1/4h^2 (f(x+h,y+h)-f(x+h,y-h)+f(x-h,y-h)-f(x-h,y+h))

        for (parameter_count=0; parameter_count<parameterList->lLength-1; parameter_count++) {
            long        iidx = useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count;

            hyFloat  ival  = GetIthIndependent(iidx),
                        locHi = 1/8192.;//funcValues (i,4);

            for (j=parameter_count+1; j<parameterList->lLength; j++) {
                long        jidx = useIndirectIndexing?parameterList->lData[j]:j;

                hyFloat  jval  = GetIthIndependent(jidx),
                            locHj = locHi, //funcValues (j,4),
                            a, // f (x+h,y+h)
                            b, // f (x+h,y-h)
                            c, // f (x-h,y-h)
                            d; // f (x-h,y+h)

                SetIthIndependent (iidx,ival+locHi);
                SetIthIndependent (jidx,jval+locHj);
                a = Compute();
                SetIthIndependent (jidx,jval-locHj);
                b = Compute();
                SetIthIndependent (iidx,ival-locHi);
                c = Compute();
                SetIthIndependent (jidx,jval+locHj);
                d = Compute();

                t2 = (a-b-d+c)/(4*locHi*locHj);

                if (uim > 0.5) {
                    t2 *= (*iMap)(parameter_count,1)*(*iMap)(j,1);
                }

                hessian.Store (parameter_count,j,-t2);
                hessian.Store (j,parameter_count,-t2);
                SetIthIndependent (iidx,ival);
                SetIthIndependent (jidx,jval);
#ifndef __UNIX__
                finishedCount += 4;
                if (TimerDifferenceFunction(true)>1.) {
                    SetStatusBarValue (finishedCount/totalCount*100.,1,0);
                    TimerDifferenceFunction (false);
                }
#endif
            }
        }

    } else {
        // fill in off-diagonal elements using the f-la
        // f_xy = 1/h^2 (f(x+h,y+h)-f(x)-f_x h -f_y h -.5h^2(f_xx+f_yy))

        if (CheckEqual(cm,1.)) {
            for (parameter_count=0; parameter_count<parameterList->lLength-1; parameter_count++) {
                hyFloat t3 = GetIthIndependent(useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count),
                           t5 = hessian(parameter_count,parameter_count),
                           t6 = funcValues(parameter_count,2);

                SetIthIndependent (useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count,t3+h);
                for (j=parameter_count+1; j<parameterList->lLength; j++) {
                    hyFloat t4 = GetIthIndependent(useIndirectIndexing?parameterList->lData[j]:j);
                    SetIthIndependent (useIndirectIndexing?parameterList->lData[j]:j,t4+h);
                    t1 = Compute();
                    t2 = (t1-functionValue-(t6+funcValues(j,2)-.5*(t5+hessian(j,j))*h)*h)/(h*h);
                    hessian.Store (parameter_count,j,-t2);
                    hessian.Store (j,parameter_count,-t2);
                    SetIthIndependent (useIndirectIndexing?parameterList->lData[j]:j,t4);
#ifndef __UNIX__
                    finishedCount ++;
                    if (TimerDifferenceFunction(true)>1.) {
                        SetStatusBarValue (finishedCount/totalCount*100.,1,0);
                        TimerDifferenceFunction (false);
                    }
#endif
                }
                SetIthIndependent (useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count,t3);
            }
        }
    }
    // undo changes to var values if needed

    DoneComputing();

    if (iMap) {
        DeleteObject (iMap);
        setParameter (intervalMappingMethod, mapMethod);
        DeleteObject (mapMethod);
    }

    for (parameter_count=0; parameter_count<parameterList->lLength; parameter_count++) {
        t1 = funcValues(parameter_count,3);
        t2 = GetIthIndependent (useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count);
        if (!CheckEqual (t1,t2)) {
            SetIthIndependent (useIndirectIndexing?parameterList->lData[parameter_count]:parameter_count,t1);
        }
    }
    return hessian.Inverse();
    //return (_Matrix*)hessian.makeDynamic();
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::GetGradientStepBound (_Matrix& gradient,hyFloat& left, hyFloat& right, long * freezeCount)
{
    left = right = DEFAULTPARAMETERUBOUND;

    if (freezeCount) {
        *freezeCount = 0;
    }

    for (unsigned long i = 0; i < indexInd.lLength; i++) {
        hyFloat directionalStep = gradient.theData[i];
        if (directionalStep) {
            hyFloat currentValue = GetIthIndependent (i),
                       ub                       = GetIthIndependentBound (i,false)-currentValue,
                       lb                       = currentValue-GetIthIndependentBound (i,true);

            //if (ub < 1e-10 || lb < 1e-10)
            //  printf ("[HIT BOUNDARY AT %s; %g, %g]\n", cv->GetName()->sData, lb, ub);

            if (directionalStep > 0.0) {
                ub /= directionalStep;
                lb /= directionalStep;
            } else {
                currentValue = -lb/directionalStep;
                lb = -ub/directionalStep;
                ub = currentValue;
            }

            left  = MIN(left, lb);
            if (ub < 1e-6 && freezeCount) {
                (*freezeCount) ++;
                gradient.theData[i] = 0.;
            } else {
                right = MIN(right, ub);
            }
        }
    }

    if (left < 1.-8) {
        left = 0.;
    }
    if (right < 1.-8) {
        right = 0.;
    }

    left = -left;
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::ComputeGradient (_Matrix& gradient,  hyFloat& gradientStep, _Matrix& values,_SimpleList& freeze, long order, bool normalize)
{
    hyFloat funcValue;

    //CheckStep     (gradientStep,unit,&values);
    /*if (order>1)
    {
        _Matrix nG (unit);
        nG*=-1;
        CheckStep (gradientStep,nG,&values);
    }
    if (gradientStep==0)
        return;*/

    if (order==1) {
        funcValue = Compute();
        for (long index=0; index<indexInd.lLength; index++) {
            if (freeze.Find(index)!=-1) {
                gradient[index]=0.;
            } else {
                //_Variable  *cv            = GetIthIndependentVar (index);
                hyFloat currentValue = GetIthIndependent(index),
                           ub           = GetIthIndependentBound(index,false)-currentValue,
                           lb            = currentValue-GetIthIndependentBound(index,true),
                           testStep    = MAX(currentValue * gradientStep,gradientStep);

                if (testStep >= ub) {
                  if (testStep < lb) {
                    testStep = -testStep;
                  } else {
                    if (ub > lb) {
                      testStep = ub;
                    } else {
                      if (lb >= ub)
                        testStep = lb == 0. ? 0. : -lb;
                    }
                  }
                }


                if (testStep) {
                    SetIthIndependent(index,currentValue+testStep);
                    gradient[index]=(Compute()-funcValue)/testStep;
                    SetIthIndependent(index,currentValue);
                } else {
                    gradient[index]= 0.;
                }

                /*if (verbosity_level > 50)
                {
                    printf ("[GRADIENT @ %s, [%g-%g-%g], %g. der = %g]\n", cv->GetName()->sData, lb, currentValue, ub, testStep, gradient[index]);
                }*/
            }
        }
    } else {
        for (long index=0; index<indexInd.lLength; index++) {
            if (freeze.Find(index)!=-1) {
                gradient[index]=0.;
            } else {
                SetIthIndependent(index,GetIthIndependent(index)-gradientStep);
                hyFloat temp = Compute();
                SetIthIndependent(index,GetIthIndependent(index)+2*gradientStep);
                gradient[index]=(Compute()-temp)/gradientStep*0.5;
                SetIthIndependent(index,GetIthIndependent(index)-gradientStep);
            }
        }

    }


    // normalize the gradient
    if (normalize) {
        funcValue = 0.;
        for (long index=0; index<indexInd.lLength; index++) {
            funcValue+=gradient.theData[index]*gradient.theData[index];
        }

        //printf ("%10.10g\n",  funcValue);

        if (CheckEqual (funcValue,0.0)) {
            return;
        }

        funcValue = 1./sqrt(funcValue);
        for (long index=0; index<indexInd.lLength; index++) {
            gradient[index]*=funcValue;
        }
    }
}
//_______________________________________________________________________________________

bool    _LikelihoodFunction::SniffAround (_Matrix& values, hyFloat& bestSoFar, hyFloat& step)
{
    for (long index = 0; index<indexInd.lLength; index++) {

        hyFloat lowerBound       = GetIthIndependentBound(index, true),
                   tryStep          = step,
                   funcValue,
                   upperBound       = GetIthIndependentBound(index, false),
                   practicalUB      = upperBound>1000?1000:upperBound,
                   val              = GetIthIndependent     (index);

        // try moving backwards
        while (val-tryStep<lowerBound+1e-8) {
            tryStep/=2;
            if (tryStep<1e-8) {
                break;
            }
        }
        if (tryStep>=1e-8) {
            SetIthIndependent(index,val-tryStep);
            funcValue = Compute();
            if (funcValue>bestSoFar) {
                bestSoFar = funcValue;
                values[index]=val-tryStep;
                return true;
            }

        }
        tryStep = step  ;
        while (val+tryStep>practicalUB-1e-8) {
            tryStep/=2;
            if (tryStep<1e-8) {
                break;
            }
        }
        if (tryStep>=1e-8) {
            SetIthIndependent(index,val+tryStep);
            funcValue = Compute();
            if (funcValue>bestSoFar) {
                bestSoFar = funcValue;
                values[index]=val-tryStep;
                return true;
            }

        }
        SetIthIndependent(index,val);


    }
    return false;
}

//_______________________________________________________________________________________

hyFloat    _LikelihoodFunction::ConjugateGradientDescent (hyFloat precision, _Matrix& bestVal, bool localOnly, long iterationLimit, _SimpleList* only_these_parameters, hyFloat check_value) {

    hyFloat     gradientStep     = STD_GRAD_STEP,
                maxSoFar          = Compute(),
                initial_value     = maxSoFar,
                currentPrecision = localOnly?precision:.01;

    if (check_value != INFINITY) {
        if (!CheckEqual(check_value, maxSoFar)) {
            _String errorStr = _String("Internal error in _LikelihoodFunction::ConjugateGradientDescent. The function evaluated at current parameter values [") & maxSoFar & "] does not match the last recorded LF maximum [" & check_value & "]";
            if (check_value - 0.01 > maxSoFar) {
                if (optimizatonHistory) {
                    ReportWarning (_String ((_String*)optimizatonHistory->toStr()));
                }
                _TerminateAndDump (errorStr);
                return check_value;
            }
            //return;
        }
    }


    _SimpleList freeze;

    if (only_these_parameters) {
        only_these_parameters->Sort();
        _SimpleList all (indexInd.lLength,0,1);
        freeze.Intersect (all, *only_these_parameters);
    }


    char        buffer[1024];

    if (verbosity_level>1) {
        snprintf (buffer, sizeof(buffer),"\nConjugate Gradient Pass %d, precision %g, gradient step %g, max so far %15.12g\n",0,precision,gradientStep,maxSoFar);
        BufferToConsole (buffer);
    }


    const unsigned long dim = bestVal.GetHDim() * bestVal.GetVDim();
    _Matrix     gradient (bestVal),
                current_direction,
                previous_direction,
                previous_gradient;
    
    hyFloat gradL;
    ComputeGradient     (gradient, gradientStep, bestVal, freeze, 1, false);

    gradL = gradient.AbsValue ();


    if (gradL > 0.0) {
        current_direction   = gradient;
        
        for (long index = 0; index< MAX (dim, 10) && index < iterationLimit; index++, currentPrecision*=0.25) {
            hyFloat current_maximum = maxSoFar;

            if (currentPrecision < 0.00001) {
                currentPrecision = 0.00001;
            }

            hyFloat line_search_precision = localOnly?precision:currentPrecision;
            GradientLocateTheBump(line_search_precision, maxSoFar, bestVal, current_direction);
            
            LoggerAddGradientPhase (line_search_precision);
            LoggerAllVariables ();
            LoggerLogL (maxSoFar);

            if (verbosity_level>1) {
                snprintf (buffer, sizeof(buffer),"Conjugate Gradient Pass %ld, precision %g, gradient step %g, max so far %15.12g\n",index+1,precision,gradientStep,maxSoFar);
                BufferToConsole (buffer);
            }
            //if (localOnly) {
            if (fabs(maxSoFar-current_maximum)<=precision) {
                break;
            }
            
            previous_gradient = gradient;
            ComputeGradient (gradient, gradientStep, bestVal, freeze, 1, false);
            
            if (CheckEqual(gradL,0.0)) {
                break;
            }
 
            previous_direction = current_direction;
            hyFloat beta = 0., scalar_product = 0.;

            // use Polak–Ribière direction
            for (unsigned long i = 0UL; i < dim; i++) {
                scalar_product += previous_gradient.theData[i] * previous_gradient.theData[i];
                beta += gradient.theData[i] * ( gradient.theData[i] - previous_gradient.theData[i]);
            }
            
            // use Dai–Yuan
            /*
             for (unsigned long i = 0UL; i < dim; i++) {
               beta += gradient.theData[i] * gradient.theData[i];
               scalar_product += previous_direction.theData[i] * ( gradient.theData[i] - previous_gradient.theData[i]);
             }
             beta = -beta;
            */

            // Hestenes-Stiefel
            /*
             for (unsigned long i = 0UL; i < dim; i++) {
               beta += gradient.theData[i] * ( gradient.theData[i] - previous_gradient.theData[i]);
               scalar_product += previous_direction.theData[i] * ( gradient.theData[i] - previous_gradient.theData[i]);
             }
             beta = -beta;
            */

            beta /= scalar_product;
            beta = MAX (beta, 0.0);
            previous_direction = current_direction;
            previous_gradient = previous_direction;
            previous_gradient *= beta;
            current_direction  = gradient;
            current_direction += previous_gradient;
        }
    }

    SetAllIndependent (&bestVal);
    if (maxSoFar < initial_value && CheckEqual(maxSoFar, initial_value,  kMachineEpsilon * 100.) == false) {
        HandleApplicationError (_String("Internal optimization error in _LikelihoodFunction::ConjugateGradientDescent. Worsened likelihood score from ") & initial_value & " to " & maxSoFar);
    }

    if (verbosity_level>1) {
        BufferToConsole("\n");
    }

    return maxSoFar;
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::GradientDescent (hyFloat& gPrecision, _Matrix& bestVal)
{

    hyFloat      currentPrecision = 0.1,
                    gradientStep     = STD_GRAD_STEP,
                    temp,
                    tryStep,
                    maxSoFar = Compute(),
                    bestTry;

    _SimpleList     leastChange,
                    freeze,
                    countLC;

    _Matrix         gradient (bestVal);

    long            vl = verbosity_level,
                    index;

    while (currentPrecision>=gPrecision && freeze.lLength<indexInd.lLength) {

        gradientStep = STD_GRAD_STEP;
        if (terminate_execution) {
            return;
        }

        char    buffer[128];
        ComputeGradient (gradient, gradientStep, bestVal, freeze, 1);
        if (gradientStep==0) {
            break;
        }
        bool done = false;
        tryStep = currentPrecision;
        while(!done) {
            CheckStep (tryStep,gradient, &bestVal);
            if (tryStep == 0.0) {
                long wereFrozen = freeze.lLength;
                for (index=0; index<indexInd.lLength; index++) {
                    if (freeze.Find(index)>=0) {
                        continue;
                    }

                    /*if ((GetIthIndependent(index)-GetIthIndependentBound(index,true)<1.0e-20)
                            && (gradient(0,index)<0.0)) {
                        freeze<<index;
                        break;
                    }
                    if ((-GetIthIndependent(index)+GetIthIndependentBound(index,false)<1.0e-20)
                            && (gradient(0,index)>0.0)) {
                        freeze<<index;
                        break;
                    }*/
                }

                 if (freeze.lLength==wereFrozen) {
                    return;
                }
                break;
            }
            for (index=0; index<indexInd.lLength; index++) {
                SetIthIndependent(index,bestVal(0,index)+tryStep*gradient(index,0));
            }
            temp = Compute();

            if (temp>maxSoFar) {
                if (vl>=5) {
                    snprintf (buffer, sizeof(buffer),"\nMoving down along the gradient with step %g value %g", tryStep, temp);
                    BufferToConsole (buffer);
                }
                _Matrix delta;
                delta = gradient;
                delta*=tryStep;
                if ((temp-maxSoFar)/fabs(maxSoFar)<gPrecision) {
                    bestVal += delta;
                    return;
                }

                maxSoFar = temp;
                bestVal += delta;
                //see which variable changed the least
                temp = INFINITY;
                long  suspect,f;
                for (long i = 0; i<indexInd.lLength; i++) {
                    if (fabs(delta(i,0))<temp) {
                        if (freeze.Find(i)!=-1) {
                            continue;
                        }
                        temp = fabs(delta(i,0));
                        suspect = i;
                    }

                }

                f = leastChange.Find(suspect);
                if (f==-1) {
                    leastChange<<suspect;
                    countLC<<1;
                }

                if (tryStep<currentPrecision/10) {
                    currentPrecision/=10;
                } else {
                    tryStep*=2;
                }
                done = true;
            } else {
                tryStep/=10;
                if (tryStep<gPrecision) {
                    for (index=0; index<indexInd.lLength; index++) {
                        SetIthIndependent(index,bestVal(0,index));
                    }
                    if (leastChange.lLength) {
                        SortLists(&leastChange, &countLC);
                        if (vl>=5) {
                            snprintf (buffer, sizeof(buffer),"\nFreezing Variable %ld",leastChange(leastChange.lLength-1));
                            BufferToConsole (buffer);
                        }
                        freeze<<leastChange(leastChange.lLength-1);
                        leastChange.Delete (leastChange.lLength-1);
                        countLC.Delete (countLC.lLength-1);
                        currentPrecision = sqrt(gPrecision);
                        done = true;
                        break;
                    }
                    currentPrecision=0;
                    break;
                }
                if (vl>=5) {
                    snprintf (buffer, sizeof(buffer),"\nShrinking step to %g (%g %g)", tryStep, tryStep, gPrecision);
                    BufferToConsole (buffer);
                }
            }
        }
    }
}



//_______________________________________________________________________________________

  void    _LikelihoodFunction::GradientLocateTheBump (hyFloat gPrecision, hyFloat& maxSoFar, _Matrix& bestVal, _Matrix& gradient) {
    DetermineLocalUpdatePolicy           ();
    hyFloat  leftValue   = maxSoFar,
    middleValue = maxSoFar,
    rightValue  = maxSoFar,
    initialValue = maxSoFar,
    bp          = gPrecision*.1,
    left = 0., right = 0., middle = 0.;
      
      
    _Matrix                          left_vector        ;
    GetAllIndependent               (left_vector);
      
    _Matrix         right_vector    (left_vector),
                    middle_vector   (left_vector),
                    prior_parameter_values    (left_vector);
    // _GrowingVector brentHistory;

    middle_vector = bestVal;
      
    int  outcome = Bracket(-1, left,middle,right,leftValue, middleValue, rightValue,bp, true, &gradient);

    if (middleValue < initialValue) {
      SetAllIndependent (&bestVal);
      FlushLocalUpdatePolicy();
      return;
    }

    if (outcome >=0 && (leftValue > middleValue || rightValue > middleValue)) {
      HandleApplicationError (_String ("Internal error in  _LikelihoodFunction::GradientLocateTheBump: bracket reported successful (") & (long)outcome & "), but likelihood values are inconsistent with it. " & leftValue & " / " & middleValue & " / " & rightValue & " initial value = " & maxSoFar);
      return;
    }

    //printf ("[LogL = %.20g GRADIENT BRACKET %g/%.20g, %g/%.20g, %g/%.20g; %d]\n",maxSoFar,lV,leftValue,ms,middleValue,rV,rightValue, outcome);

    left_vector.AplusBx   (gradient, left);
    middle_vector.AplusBx (gradient, middle);
    right_vector.AplusBx  (gradient, right);
      
    bool reset = false;

    if (outcome!=-1) { // successfull bracket
                       // set up left, right, middle

      if (outcome == -2) {
        if (verbosity_level > 50) {
          char buf [256];
          snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump BRACKET == -2 clause]");
          BufferToConsole (buf);
        }
          
        if (middleValue>maxSoFar) {
            maxSoFar = middleValue;
            bestVal  = middle_vector;
            SetAllIndependent (&middle_vector);
        } else {
           SetAllIndependent (&bestVal);
        }

        FlushLocalUpdatePolicy();
        return;
      }



      if (outcome == indexInd.countitems()) {
        reset = true;
      } else {
          hyFloat U,V,W,X=middle,E=0.,FX,FW,FV,XM,R,Q,P,ETEMP,D=0.,FU;
          _Matrix current_best_vector (prior_parameter_values);
          current_best_vector.AplusBx(gradient, middle);
          W = .0;
          V = .0;
          FX = -middleValue;
          FV = FX;
          FW = FX;
          outcome = 0;
          while (outcome < 20) {
          // brentHistory.Store (-FX - initialValue);
              bool parabolic_step = false;
              XM = .5*(left+right);

              if (verbosity_level > 50) {
                  char buf [256];
                  snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump (current max = %20.16g) GOLDEN RATIO INTERVAL CHECK: %g <= %g (%g = %g) <= %g, span = %g]", middleValue, left, XM, X, fabs(X-XM), right, right-left);
                  BufferToConsole (buf);
              }
              
              hyFloat tol1 = fabs (X) * MIN (gPrecision, 1e-7) + kMachineEpsilon,
              tol2 = 2.*tol1;

              if (fabs(X-XM) <= tol2 && outcome > 0) {
                break;
              }

              if (fabs(E)>tol1) {
                  R = (X-W)*(FX-FV);
                  Q = (X-V)*(FX-FW);
                  P = (X-V)*Q-(X-W)*R;
                  Q = 2.0 * (Q-R);
                  if (Q>0.) {
                    P = -P;
                  }
                  Q = fabs(Q);
                  ETEMP = E;
                  E = D;
                  if (!(fabs (P) > fabs (.5*Q*ETEMP) || P <= Q * (left-X) || P >= Q *( right-X))) {
                      parabolic_step = true;
                      D = P/Q;
                      U = X+D;
                      if (U - left < tol2 || right - U < tol2) {
                        D = (XM - X >= 0.) ? tol1 : -tol1;
                      }
                  }
              }
          

              if (!parabolic_step) {
                E = (X >= XM ? left : right) - X;
                D = GOLDEN_RATIO_C * E;
              }
              U = fabs (D) >= tol1 ? X + D : X + (D > 0. ? tol1 : -tol1);
              FU = -SetParametersAndCompute (-1,U,&prior_parameter_values,&gradient);
              if (verbosity_level > 50) {
                  char buf [256];
                  snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump GOLDEN RATIO TRY: param %20.16g, log L %20.16g (%20.16g)]",  U, -FU, -FX);
                  BufferToConsole (buf);
              }

              if (FU<=FX) { // accept the move
                    if (verbosity_level > 50) {
                    char buf [256];
                    snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump (eval %ld) ACCEPT new try, confirm value %20.16g (delta = %20.16g)", likeFuncEvalCallCount,  U, FX-FU);
                    BufferToConsole (buf);
                    }
                    current_best_vector = prior_parameter_values;
                    current_best_vector.AplusBx(gradient, U);
                    if (U>=X) {
                      left = X;
                    } else {
                      right = X;
                    }
                    V = W;
                    FV = FW;
                    W = X;
                    FW = FX;
                    X = U;
                    FX = FU;
              } else {
                    if (verbosity_level > 50) {
                     char buf [256];
                     snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump (eval %ld) REJECT new try (%20.16g) (delta = %20.16g)", likeFuncEvalCallCount, U, FX-FU);
                     BufferToConsole (buf);
                    }
                    if (U<X) {
                      left = U;
                    } else {
                      right = U;
                    }
                    if (FU<=FW || W==X ) {
                      V = W;
                      FV = FW;
                      W = U;
                      FW = FU;
                    } else {
                      if (FU<=FV || V==X || V==W) {
                        V = U;
                        FV = FU;
                      }
                    }
              }
              outcome++;
        }

        if (verbosity_level > 50) {
         char buf [256];
         snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump GOLDEN RATIO SEARCH SUCCESSFUL: precision %g, parameter moved from %15.12g to %15.12g, Log L new/old = %15.12g/%15.12g ]\n\n", gPrecision, X, -FX, middleValue, maxSoFar);
          BufferToConsole (buf);
        }
        middleValue = -FX;
        //brentHistory.Store (0.);
        if (middleValue < maxSoFar) {
            if (verbosity_level > 50) {
                char buf [256];
                snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump RESETTING THE VALUE (worse log likelihood obtained; current value %20.16g, best value %20.16g) ]\n\n", middleValue, maxSoFar);
                BufferToConsole (buf);
            }
            SetAllIndependent (&bestVal);
            maxSoFar = middleValue;
        } else {
          SetAllIndependent (&current_best_vector);
          maxSoFar    = Compute();
          bestVal     = current_best_vector;
          if (verbosity_level > 50) {
            char buf [256];
            snprintf (buf, 256, "\n\t[_LikelihoodFunction::GradientLocateTheBump moving parameter value (should trigger LL update) %15.12g ||L2|| move ]\n\n", (current_best_vector-bestVal).AbsValue());
            BufferToConsole (buf);
          }
        }

        if (maxSoFar < initialValue && !CheckEqual (maxSoFar, initialValue, 100. * kMachineEpsilon)) {
          HandleApplicationError (_String ("Internal error in  _LikelihoodFunction::GradientLocateTheBump: in the Brent loop iteration ") & long(outcome) & ". " & _String (maxSoFar, "%18.16g") & " / " & _String (initialValue,"%18.16g") & ".\n");

          return;
        }

        //bestVal = middle;
        //maxSoFar = middleValue;
      }
      //middle = X;
    }

    else {
      reset = true;
      if (verbosity_level>1) {
        BufferToConsole ("Line optimization unsuccessful\n");
      }
      if (leftValue>middleValue) {
        middleValue = leftValue;
        middle_vector = left_vector;
      }
      if (rightValue>middleValue) {
        middleValue = rightValue;
        middle_vector = right_vector;
      }

      if (middleValue>maxSoFar) {
        SetAllIndependent (&middle_vector);
        maxSoFar = middleValue;
        reset = false;
      }
    }

    if (reset)
      SetAllIndependent (&bestVal);

    FlushLocalUpdatePolicy();
  }

//_______________________________________________________________________________________

void    _LikelihoodFunction::LocateTheBump (long index,hyFloat gPrecision, hyFloat& maxSoFar, hyFloat& bestVal, bool go2Bound, hyFloat bracketSetting) {
    hyFloat left,
               right,
               middle           = bestVal,
               leftValue,
               middleValue       = maxSoFar,
               rightValue,
               bp               = 2.*gPrecision,
               brentPrec        = bracketSetting>0.?bracketSetting:gPrecision;

    DetermineLocalUpdatePolicy           ();

    /*if (optimizatonHistory && ((_AssociativeList*)this->optimizatonHistory->GetByKey("Phases"))->Length() == 2171) {
        verbosity_level = 1000;
    } else {
        verbosity_level = 1;
    }*/

    unsigned long        inCount = likeFuncEvalCallCount;
    int outcome = Bracket (index,left,middle,right,leftValue, middleValue, rightValue,bp, go2Bound);
    unsigned long        bracketCount = likeFuncEvalCallCount - inCount;

    if (outcome != -1) { // successfull bracket
        if (right - left < 4*brentPrec) {
           brentPrec = (right-left) * 0.2;
           //printf ("\nResetting brentPrec to %g\n", brentPrec, "\n");
        }
        /*if (index >= 0) {
         // try Newton Raphson
         
         if (middle - left > STD_GRAD_STEP && right - middle > STD_GRAD_STEP) {
         
         
         _Parameter last_value, current_value = middle, current_fx = middleValue;
         
         auto store_max = [&] (_Parameter x, _Parameter fx) -> void {
         if (fx > maxSoFar) {
         maxSoFar = fx;
         bestVal  = x;
         }
         };
         
         do  {
         last_value = current_value;
         
         _Parameter x_plus_h  = last_value + STD_GRAD_STEP,
         x_minus_h = last_value - STD_GRAD_STEP,
         fx_plus_h  = SetParametersAndCompute(index, x_plus_h),
         fx_minus_h = SetParametersAndCompute(index, x_minus_h),
         dFdX = (fx_plus_h - fx_minus_h) / (2. * STD_GRAD_STEP),
         d2FdX2 = ((fx_plus_h - current_fx) + (fx_minus_h - current_fx)) / (STD_GRAD_STEP * STD_GRAD_STEP);
         
         
         store_max (x_plus_h, fx_plus_h);
         store_max (x_minus_h, fx_minus_h);
         
         if (CheckEqual(d2FdX2, 0.0)) {
         current_value = last_value;
         } else {
         current_value = last_value - dFdX / d2FdX2;
         }
         
         if (current_value < left || current_value > right) {
         break;
         }
         
         //printf ("\n\nf(%20.16g) = %20.16g; f(%20.16g) = %20.16g\n", x_plus_h, fx_plus_h, x_minus_h, fx_minus_h);
         //printf ("f(%g) = %g; dF = %g, dF2 = %g\n", last_value, current_fx, dFdX, d2FdX2);
         
         if (CheckAndSetIthIndependent(index, current_value)) {
         current_fx = SetParametersAndCompute (index, current_value);
         }
         //printf (" == move by %g with value %20.16g\n", current_value-last_value, current_fx);
         
         store_max (current_value, current_fx);
         
         
         } while (fabs (current_value - last_value) >= brentPrec && current_value - left > STD_GRAD_STEP && right - current_value > STD_GRAD_STEP);
         
         if (fabs (current_value - last_value) < brentPrec && current_value - left > STD_GRAD_STEP && right - current_value > STD_GRAD_STEP) {
         CheckAndSetIthIndependent(index, bestVal);
         //printf (" == SUCCESS\n");
         FlushLocalUpdatePolicy            ();
         return;
         }
         CheckAndSetIthIndependent(index, bestVal);
         middle = bestVal;
         middleValue = maxSoFar;
         //printf (" == FAILURE %20.16g -> %20.16g\n", middle, middleValue);
         
         }
         }*/
        
        hyFloat U,V,W,X=middle,E=0.,FX,FW,FV,XM,R,Q,P,ETEMP,D=0.,FU;
        W       = middle;
        V       = middle;
        FX      = -middleValue;
        FV      = FX;
        FW      = FX;
        outcome = 0;


        while (outcome < 20) {
            XM = .5*(left+right);

            bool parabolic_step = false;

            if (verbosity_level > 50) {
                char buf [256];
                snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (index %ld) (current max = %15.12g) GOLDEN RATIO INTERVAL CHECK: %g <= %g (%g = %g) <= %g, span = %g]", index, bestVal, left, XM, X, fabs(X-XM), right, right-left);
                BufferToConsole (buf);
            }

            if (fabs(X-XM) <= brentPrec && outcome > 0) {
              break;
            }

            hyFloat tol1 = fabs (X) * Minimum (brentPrec, 1e-7) + kMachineEpsilon,
                       tol2 = 2.*tol1;



            if (fabs(E)>tol1) {
                R = (X-W)*(FX-FV);
                Q = (X-V)*(FX-FW);
                P = (X-V)*Q-(X-W)*R;
                Q = 2.0 * (Q-R);
                if (Q>0.) {
                    P = -P;
                }
                Q = fabs(Q);
                ETEMP = E;
                E = D;
                if (!(fabs (P) > fabs (.5*Q*ETEMP) || P <= Q * (left-X) || P >= Q *( right-X))) {
                  parabolic_step = true;
                  D = P/Q;
                  U = X+D;
                  if (U - left< tol2 || right - U < tol2) {
                    D = (XM - X >= 0.) ? tol1 : -tol1;
                  }
                }
            }


            if (!parabolic_step) {
              E = (X >= XM ? left : right) - X;
              D = GOLDEN_RATIO_C * E;
            }

            U = fabs (D) >= tol1 ? X + D : X + (D > 0. ? tol1 : -tol1);

          //U = X + D;
            SetIthIndependent (index,U);
            FU = -Compute();

            if (verbosity_level > 50) {
                char buf [256];
                snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (index %ld) GOLDEN RATIO TRY: param %20.16g, log L %20.16g]", index, U, -FU);
                BufferToConsole (buf);
            }

            if (FU<=FX) { // value at U is the new minimum
                if (verbosity_level > 50) {
                    char buf [256];
                    snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (eval %ld) ACCEPT new try, confirm value %20.16g (delta = %20.16g)", likeFuncEvalCallCount,  GetIthIndependent(index), FX-FU);
                    BufferToConsole (buf);
                }

                if (U>=X) {
                    left = X;
                } else {
                    right = X;
                }
                V = W;
                FV = FW;
                W = X;
                FW = FX;
                X = U;
                FX = FU;
            } else { // value at X remains the minimum
                if (verbosity_level > 50) {
                    char buf [256];
                    snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (eval %ld) REJECT new try (%20.16g) (delta = %20.16g)", likeFuncEvalCallCount, U, FX-FU);
                    BufferToConsole (buf);
                }

                if (U < X) {
                    left = U;
                } else {
                    right = U;
                }
                if (FU<=FW || W==X) {
                    V = W;
                    FV = FW;
                    W = U;
                    FW = FU;
                } else {
                    if (FU<=FV || V==X || V==W) {
                        V = U;
                        FV = FU;
                    }
                }
            }
            outcome++;

        }

        if (verbosity_level > 50) {
            char buf [256];
            snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (index %ld) GOLDEN RATIO SEARCH SUCCESSFUL: precision %g, parameter moved from %15.12g to %15.12g, Log L new/old = %15.12g/%15.12g ]\n\n", index, brentPrec, bestVal, X, -FX, maxSoFar);
            BufferToConsole (buf);
        }
        middleValue = -FX;
        middle      = X;

        if (middleValue<maxSoFar) {
            if (verbosity_level > 50) {
              char buf [256];
              snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (index %ld) RESETTING THE VALUE (worse log likelihood obtained; current value %20.16g, best value %20.16g) ]\n\n", index, GetIthIndependent(index), bestVal);
              BufferToConsole (buf);
            }
            if (CheckEqual(GetIthIndependent(index), bestVal) && fabs (middleValue-maxSoFar) > 1e-9) {            
                char buf[256];
                snprintf (buf, 256, " \n\tERROR: [_LikelihoodFunction::LocateTheBump (index %ld) current value %20.16g (parameter = %20.16g), best value %20.16g (parameter = %20.16g)); delta = %20.16g ]\n\n", index, middleValue, GetIthIndependent(index), maxSoFar, bestVal, maxSoFar - middleValue);
                _TerminateAndDump (_String (buf) & "\n" &  "\nParameter name " & *GetIthIndependentName(index));
            }
            SetIthIndependent(index,bestVal);
        } else {
            if (!CheckEqual(GetIthIndependent(index),middle)) {
                if (verbosity_level > 50) {
                    char buf [256];
                    snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (index %ld) moving parameter value (should trigger LL update) %15.12g to %15.12g ]\n\n", index, GetIthIndependent(index), middle);
                    BufferToConsole (buf);
                }
                SetIthIndependent (index,middle);
            } else {
                if (verbosity_level > 50) {
                    char buf [256];
                    snprintf (buf, 256, "\n\t[_LikelihoodFunction::LocateTheBump (index %ld) KEEPS parameter value (no LL update) %15.12g == %15.12g ]\n\n", index, GetIthIndependent(index), middle);
                    BufferToConsole (buf);
                }
            }
            maxSoFar = middleValue;
        }
    }

    if (index >= 0) {
      LoggerSingleVariable (index, maxSoFar, bp, brentPrec, outcome != -1 ? right-left : -1., bracketCount, likeFuncEvalCallCount-inCount-bracketCount);
    }

    oneDFCount += likeFuncEvalCallCount-inCount-bracketCount;
    oneDCount ++;
    FlushLocalUpdatePolicy            ();
}

//_______________________________________________________________________________________
bool        _LikelihoodFunction::checkPermissibility (_Matrix&m, long row)
{
    for (unsigned long j=0; j<indexInd.lLength; j++) {
        hyFloat junk = m (row,j);

        _Variable *v = LocateVar (indexInd(j));

        if (junk<v->GetLowerBound()) {
            return FALSE;
        } else if (junk>v->GetUpperBound()) {
            return FALSE;
        }
    }

    return true;
}

//_______________________________________________________________________________________
hyFloat      _LikelihoodFunction::computeAtAPoint (_Matrix&m, long row)
{

    if (!checkPermissibility(m,row)) {
        return -1e300;
    }

    for (long k=0; k < indexInd.lLength; k++) {
        SetIthIndependent (k, m(row,k));
    }

    return Compute();
}

//_______________________________________________________________________________________
hyFloat      _LikelihoodFunction::replaceAPoint (_Matrix&m, long row, _Matrix&p, hyFloat& nV, _Matrix& fv)
{
    for (long k=0; k < indexInd.lLength; k++) {
        m.Store (row,k,p(0,k));
    }
    fv.Store (0,row,nV);
    return nV;
}


//_______________________________________________________________________________________
hyFloat      _LikelihoodFunction::SimplexMethod (hyFloat& gPrecision, unsigned long iterations, unsigned long max_evaluations) {
    
#define DEFAULT_STEP 0.05
#define DEFAULT_STEP_OFF_BOUND 0.00025
    
    /** "Implementing the Nelder-Mead simplex algorithm with adaptive parameters"
         DOI 10.1007/s10589-010-9329-3
     
        "CONVERGENCE PROPERTIES OF THE NELDER–MEAD SIMPLEX METHOD IN LOW DIMENSIONS"
         SIAM J. OPTIM. Vol. 9, No. 1, pp. 112–147
     **/
     
    // the dimension of the problem
    
    _OptimiztionProgress progress_tracker;
    
    long  N = indexInd.countitems(),
          lf_evaluations = 0L;
    
    hyFloat N_inv = 1./N,
            simplex_alpha = 1.,
            simplex_beta  = 1. + 2. * N_inv,
            simplex_gamma = .75 - 0.5*N_inv,
            simplex_delta = 1. - N_inv;
    
    _Constant zero (0.);
    
    // allocate the points of the simplex + scratch
    _Matrix * simplex = new _Matrix [N+1],
    centroid (N, 1, false, true),
    function_values = _Matrix (N+1, 2, false, true);

    auto replace_point = [&] (_Matrix const& new_point, hyFloat new_value, long index, long sorted_index) -> void {
        simplex[index] = new_point;
        function_values.Store (sorted_index, 0, new_value);
    };
    
    auto set_point_and_compute = [this,&lf_evaluations] (_Matrix& v) -> hyFloat {
        lf_evaluations++;
        SetAllIndependent (&v);
        return -Compute();
    };
    
    auto resort_values = [&zero] (_Matrix& v) -> void {
        _Matrix *sorted = (_Matrix *)v.SortMatrixOnColumn(&zero);
        v = *sorted;
        DeleteObject (sorted);
    };

    GetAllIndependent(simplex[0]);
    
    // current FEASIBLE point
    function_values.Store (0,0, set_point_and_compute (simplex[0]));
                            
    for (long i = 0L; i < N; i++) {
        simplex[i+1] = simplex[0];
        
        
        hyFloat ith_coordinate = GetIthIndependent(i),
                lb = GetIthIndependentBound(i, true),
                ub = GetIthIndependentBound(i, false);
        
#ifdef NEDLER_MEAD_DEBUG
        printf ("\nCoordinate %ld, value %g, range [%g, %g]\n", i, ith_coordinate, lb, ub);
#endif

        if (CheckEqual(ith_coordinate, lb)) {
            simplex[i+1][i] += MIN(DEFAULT_STEP_OFF_BOUND, (ub-lb)*0.5);
        }
        if (ub - ith_coordinate > DEFAULT_STEP) {
            simplex[i+1][i] += DEFAULT_STEP;
        } else {
            if (ub - ith_coordinate > ith_coordinate - lb) {
                simplex[i+1][i] += (ub - ith_coordinate)*0.5;
            } else {
                simplex[i+1][i] -= (ith_coordinate - lb)*0.5;
            }
        }
        
#ifdef NEDLER_MEAD_DEBUG
       ObjectToConsole(&simplex[i+1]);
#endif
        //SetAllIndependent(&simplex[i+1]);
        function_values.Store(i+1,0, set_point_and_compute (simplex[i+1]));
        function_values.Store(i+1,1, i+1);
    }
    
    resort_values (function_values);
#ifdef NEDLER_MEAD_DEBUG
    printf ("\n**********\nSimplex iteration in\n");
    ObjectToConsole(&function_values);
#endif
    for (long it_count = 0L; it_count <= iterations && lf_evaluations <= max_evaluations; it_count ++) {
        /** compute the centroid of all EXCEPT the WORST point **/
        
       long worst_index = function_values (N,1);
       for (long k = 0L; k < N; k++) {
            hyFloat kth_coordinate = 0.;
            for (long c = 0L; c <= N; c++) {
                if (c != worst_index) {
                    kth_coordinate += simplex[c](k,0);
                }
            }
            centroid[k] = kth_coordinate * N_inv;
        }
        
        /** reflection
         x[r] = ¯x +α(¯x −x[n+1]).
        */
        
        _Matrix new_point (centroid);
        _Matrix t (centroid);
        t -= simplex [worst_index];
        t *= simplex_alpha;
        new_point += t;
        hyFloat trial_value = set_point_and_compute (new_point);

#ifdef NEDLER_MEAD_DEBUG
        printf ("\tTrial point value %g\n", trial_value);
#endif
        bool  do_shrink = false;

        if (trial_value >= function_values(0,0) && trial_value < function_values (N-1,0)) {
            // accept the reflection
#ifdef NEDLER_MEAD_DEBUG
            printf ("### ACCEPT REFLECTION replace %d\n", (long)function_values(N,1));
#endif
            replace_point (new_point, trial_value, (long)function_values(N,1), N);
        } else {
            if (trial_value < function_values(0,0)) { // try expansion
                // x[e] = ¯x +β(x[r] − ¯x)
#ifdef NEDLER_MEAD_DEBUG
                printf ("--- Trying expansion\n");
#endif
                _Matrix expansion (centroid),
                        t (new_point);
                t -= centroid;
                t *= simplex_beta;
                expansion += t;
                hyFloat expansion_value = set_point_and_compute (expansion);
                if (expansion_value < trial_value) {
#ifdef NEDLER_MEAD_DEBUG
                    printf ("### ACCEPT EXPANSION\n");
#endif
                    replace_point (expansion, expansion_value, (long)function_values(N,1), N);
                } else {
#ifdef NEDLER_MEAD_DEBUG
                    printf ("### REJECT EXPANSION, ACCEPT REFLECTION\n");
#endif
                    replace_point (new_point, trial_value, (long)function_values(N,1), N);
                }
            } else {
                //if (trial_value >= function_values (N,0)) {
                    long worst_index = function_values (N,1);
#ifdef NEDLER_MEAD_DEBUG
                    printf ("--- Trying CONTRACTION\n");
#endif
                    // x[ic] = x ̄− γ (x ̄−x[n+1])
                    // x[oc] = ¯x ̄+γ alpha (x ̄−x[n+1])
                    _Matrix cv (centroid),
                            t (centroid);
                
                    t -= simplex [worst_index];
                    t *= simplex_gamma;
                    bool inside = false;
                    if (trial_value >= function_values (N,0)) { // INSIDE
#ifdef NEDLER_MEAD_DEBUG
                        printf ("--- INSIDE contraction\n");
#endif
                        inside = true;
                        cv -= t;
                    } else {
#ifdef NEDLER_MEAD_DEBUG
                        printf ("--- OUTSIDE contraction\n");
#endif
                         t *= simplex_alpha;
                        cv += t;
                    }
                
                    hyFloat contraction_value = set_point_and_compute (cv);
                
#ifdef NEDLER_MEAD_DEBUG
                    printf ("--Contraction value = %g (relected = %g, worst = %g)\n", contraction_value, trial_value, function_values (N,0));
#endif
                    if ((inside && contraction_value < function_values (N,0)) || (!inside && contraction_value < trial_value))  {
#ifdef NEDLER_MEAD_DEBUG
                        printf ("### ACCEPT contraction\n");
#endif
                        replace_point (cv, contraction_value, worst_index, N);
                    } else {
#ifdef NEDLER_MEAD_DEBUG
                        printf ("### REJECT contraction\n");
#endif
                        do_shrink = true;
                    }
            }
        }
        if (do_shrink) {
            long best_idx = function_values(0,1);
#ifdef NEDLER_MEAD_DEBUG
            printf ("### PEFRORM SHRINKAGE\n");
            BufferToConsole("\nBest point\n");
            ObjectToConsole(&simplex[best_idx]);
#endif
            //x[i] = x[1] + δ(x[i] − x[1]).
            for (long i = 1L; i <= N; i++) {
                long idx = function_values(i,1);
                _Matrix t (simplex[idx]);
                t -= simplex[best_idx];
                t *= simplex_delta;
#ifdef NEDLER_MEAD_DEBUG
                BufferToConsole("\nOld point\n");
                ObjectToConsole(&simplex[idx]);
#endif
                simplex[idx] = simplex[best_idx];
                simplex[idx] += t;
#ifdef NEDLER_MEAD_DEBUG
                BufferToConsole("\nMoved point\n");
                ObjectToConsole(&simplex[idx]);
#endif
                function_values.Store (i, 0, set_point_and_compute(simplex[idx]));
            }
        }
        resort_values (function_values);
        
        if (verbosity_level==1) {
            UpdateOptimizationStatus (-function_values (0,0),progress_tracker.PushValue (-function_values (0,0)),1,true,progressFileString);
        } else {
            if (verbosity_level>5) {
                char buffer [2048];
                snprintf (buffer, sizeof(buffer),"Simplex iteration %10ld; current max %12.6g, current value spread %12.6g [precision %g]\n", it_count, -function_values(0,0), fabs (function_values (N,0) - function_values (1,0)), gPrecision);
                BufferToConsole (buffer);
            }
        }
        
        if (fabs (function_values (N,0) - function_values (1,0)) < gPrecision) {
            break;
        }
        
#ifdef NEDLER_MEAD_DEBUG
        for (long i = 0; i < N; i++) {
            for (long j = i + 1L; j <= N; j++) {
                if (simplex[i].Equal(&simplex[j])) {
                    printf ("COLLAPSED SIMPLEX %ld == %ld\n\n", i, j);
                    exit (1);
                }
            }
        }
        printf ("\n\n**********\nSimplex iteration %ld \n", it_count);
        ObjectToConsole(&function_values);
        BufferToConsole("\nBest point\n");
        ObjectToConsole(&simplex[(long)function_values(0,1)]);
#endif
    }
    
    
    long best_idx = function_values(0,1);
    
    N_inv = -set_point_and_compute (simplex[best_idx]);
    
    delete [] simplex;
    return N_inv;
    
    /*
    
    _Matrix     points (indexInd.lLength+1, indexInd.lLength, false, true), //the matrix of points
                functionalEvaluations (1, indexInd.lLength+1,false, true),
                scratch1 (1, indexInd.lLength,false, true),
                scratch2 (1, indexInd.lLength,false, true),
                bestSoFar (1, indexInd.lLength+1,false, true),
                temp;

    hyFloat  lastBError      = 0.,
                lastError        = 0.0,
                simplexAlpha     = 1.0,
                simplexBeta    = 0.5,
                simplexGamma    = 2.0,
                testValue,
                minError      = 1.e300,
                bumpingFactor    = 1.;

    long        iterationsCount = 0,
                nBumps          = 0;


    // first populate the matrix of initial points
    long j,k;

    for (k=0; k<indexInd.lLength; k++) {
        _Variable *v = LocateVar (indexInd (k));

        hyFloat lowP     = v->GetLowerBound(),
                   highP     = v->GetUpperBound(),
                   span     = highP-lowP;

        if (span>1) {
            lowP += 0.1;
            highP = lowP+1;
        } else {
            lowP += span/10;
            highP += lowP+span/2;
        }

        for (j=0; j<k; j++) {
            points.Store (j,k,lowP);
        }

        for (j=k+1; j<=indexInd.lLength; j++) {
            points.Store (j,k,lowP);
        }

        points.Store (k,k,highP);
    }


    for (j=0; j<=indexInd.lLength; j++) {
        functionalEvaluations.Store (0,j,computeAtAPoint(points,j));
    }

    do {
        // compute the reflection
        long            indexMin = 0,
                        indexMax = 0,
                        index2Min;

        hyFloat      max  = -1e300,
                        min  = 1e300,
                        min2 = 1e300,
                        error = 0;

        for (j=0; j<=indexInd.lLength; j++) {
            if (functionalEvaluations(0,j)>max) {
                indexMax = j;
                max = functionalEvaluations(0,j);
            }
            if (functionalEvaluations(0,j)<min) {
                if (min<min2) {
                    min2 = min;
                    index2Min = indexMin;
                }
                indexMin = j;
                min = functionalEvaluations(0,j);
            } else {
                if (functionalEvaluations(0,j)<min2) {
                    index2Min = j;
                    min2 = functionalEvaluations(0,j);
                }
            }
        }




        error = functionalEvaluations(0,indexMax)-functionalEvaluations(0,indexMin);
        if (verbosity_level>1) {
            char    buffer[64];
            snprintf (buffer, sizeof(buffer),"\n Error = %15.15g", error);
            BufferToConsole (buffer);
        }
        if (error<minError) {
            minError = error;
            for (k=0; k < indexMax; k++) {
                bestSoFar.Store(0,k,points(indexMax,k));
            }
        }

        // find the centroid of all the points excluding the worst!
        for (j=0; j<indexInd.lLength; j++) {
            hyFloat junk = 0;
            for (k=0; k < indexMin; k++) {
                junk+= points(k,j);
            }

            for (k=indexMin+1; k <= indexInd.lLength; k++) {
                junk+= points(k,j);
            }
            _Variable *v = LocateVar (indexInd(j));

            junk/=indexInd.lLength;
            if (junk<v->GetLowerBound()) {
                scratch1.Store (0,j,v->GetLowerBound());
            } else if (junk>v->GetUpperBound()) {
                scratch1.Store (0,j,v->GetUpperBound());
            } else {
                scratch1.Store (0,j,junk);
            }

        }



        if (error<gPrecision) {
            break;
        }
        if (fabs(lastError-error)<gPrecision*error) {
            iterationsCount++;
            if (iterationsCount%10==0) {
                if ((error>lastBError)||(iterationsCount%(30)==0)) {
                    //perturb the coodinates of simplex vertices to escape a vicios loop
                    nBumps++;
                    for (k=0; k<indexInd.lLength; k++) {

                        _Variable *v = LocateVar (indexInd (k));
                        hyFloat lowP = v->GetLowerBound() , highP = v->GetUpperBound(), trial;

                        if (fabs(lastBError-error)>gPrecision) {
                            bumpingFactor*=10;
                        } else {
                            bumpingFactor = 1;
                        }
                        //trial = bestSoFar(0,k)+bumpingFactor*(genrand_int32()-RAND_MAX_32)/(hyFloat)RAND_MAX_32*error;
//                      else
//                          trial = scratch1(0,k)+100*(genrand_int32()-RAND_MAX_32)/(hyFloat)RAND_MAX_32*error;
                        trial = lowP+(highP-lowP)*genrand_real1()*error;

                        if ((trial<(lowP+gPrecision))||(trial>(highP-gPrecision))) {
                            trial = lowP+ genrand_real1 () * (highP-lowP);
                        }
                        scratch1.Store (0,k,trial);
                        lastBError = error;

                    }
                    if (verbosity_level>1) {
                        char    buffer[64];
                        snprintf (buffer, sizeof(buffer),"\nBUMPING... with factor %g", bumpingFactor);
                        BufferToConsole (buffer);
                    }
                }

            }
            if (nBumps>25) {
                char str[255];
                snprintf (str, sizeof(str),"Simplex Method Failed to Converge to Desired Precision. Precision attained: %15.15g", minError);
                ReportWarning (_String(str));
                break;
            }
        } else {
            iterationsCount = 0;
        }
        lastError = error;

        for (k=0; k < indexInd.lLength; k++) {
            scratch2.Store (0,k,points(indexMin,k));
        }

        // reflect the worst point across the centroid
        _Matrix dummy;
        temp = scratch2;
        dummy = scratch1;
        dummy -= temp;
        dummy *=(1.0+simplexAlpha);
        temp += dummy;

        testValue = computeAtAPoint (temp);


        if ((testValue>functionalEvaluations(0,indexMin))&&(testValue<functionalEvaluations(0,indexMax))) {
            replaceAPoint (points, indexMin, temp, testValue, functionalEvaluations);
        } else if (testValue>functionalEvaluations(0,indexMax)) {
            replaceAPoint (points, indexMin, temp, testValue, functionalEvaluations);
            scratch2 =temp;
            scratch2*=simplexGamma;
            scratch1*=(1-simplexGamma);
            scratch2+=scratch1;
            hyFloat testValue3 = computeAtAPoint (scratch2);

            if (testValue3>functionalEvaluations(0,index2Min)) {
                replaceAPoint (points, index2Min, scratch2,testValue3, functionalEvaluations);
            }
        } else {
            for (k=0; k < indexInd.lLength; k++) {
                scratch2.Store (0,k,points(indexMin,k));
            }
            temp=scratch2;
            temp*=simplexBeta;
            scratch1*=(1-simplexBeta);
            temp+=scratch1;

            hyFloat testValue2 = computeAtAPoint (temp);
            if (testValue2>=testValue) {
                replaceAPoint (points, indexMin, temp, testValue, functionalEvaluations);
            } else {
                for (j=0; j<indexInd.lLength; j++)
                    for (k=0; k<=indexInd.lLength; k++) {
                        points.Store (j,k,(points(j,k)+points(indexMax,k))/2);
                    }
            }
        }
    } while (1);

    for (k=0; k < indexInd.lLength; k++) {
        SetIthIndependent (k, bestSoFar(0,k));
    }

    return true;
    */
}


//_______________________________________________________________________________________

void    _LikelihoodFunction::Anneal (hyFloat&)
// simple simulated annealing approach
{
    /*hyFloat jumpMagnitude , decreaseRate, startingRate, currentValue = -1e300, trialValue, bestValue = -1e300;
    hyFloat iterationsPerDecrease, decisionOnConvergence, terminateAfter;
    _Matrix jumpVector (1, indexInd.lLength,false, true),
            currentPoint (1, indexInd.lLength,false, true),
            trialPoint (1, indexInd.lLength,false, true),
            bestPoint (1, indexInd.lLength,false, true);
    long i,count = 0, smCount = 0;

    checkParameter (_String("SA_JUMP_MAGNITUDE"),jumpMagnitude,0.1);
    checkParameter (_String("SA_DECREASE_RATE"),decreaseRate,1.25);
    checkParameter (_String("SA_STARTING_RATE"),startingRate,1);
    checkParameter (_String("SA_ITERATIONS_PER_DECREASE"),iterationsPerDecrease,500);
    checkParameter (_String("SA_CONVERGENCE_THERESHOLD"),decisionOnConvergence,100);
    checkParameter (_String("SA_TERMINATION_THERESHOLD"),terminateAfter,1000000);


    for (i=0; i<indexInd.lLength; i++)
    {
        trialPoint.Store(0,i,GetIthIndependent(0));
    }


    currentValue = computeAtAPoint (currentPoint);
    while (1)
    {
        count++;
        if (count>terminateAfter)
        {
            computeAtAPoint (bestPoint);
            ReportWarning (_String("Simulated annealing terminated before convergence was attained"));
            return;
        }
        if (count%(long)iterationsPerDecrease==0)
        {
            startingRate *= decreaseRate;
            iterationsPerDecrease/=decreaseRate;
        }
        trialValue = computeAtAPoint (trialPoint);
        // compare the trial value and the current value and decide what to do;
        if (fabs(trialValue-currentValue)<gPrecision)
        {
            smCount++;
        }
        else
            smCount=0;

        if (smCount > decisionOnConvergence)
        {
            computeAtAPoint (bestPoint);
            return;
        }

        if (trialValue > currentValue)
        {
            if (VerbosityLevel()>1)
            {
                printf ("\nJump Up to %g accepted", trialValue);
            }
            currentValue = trialValue;
            currentPoint = trialPoint;
        }
        else
        {
            hyFloat decisionProb = exp((trialValue-currentValue)*startingRate);
            if (((hyFloat)genrand_int32())/RAND_MAX_32<decisionProb)
            {
                currentValue = trialValue;
                currentPoint = trialPoint;
                if (VerbosityLevel()>1)
                {
                    printf ("\nJump down to %g accepted", currentValue);
                }
            }
            else
                if (VerbosityLevel()>900)
                {
                    printf ("\nJump to %g rejected", trialValue);
                }

        }

        // see if we have the max

        if (currentValue>bestValue)
        {
            bestValue = currentValue;
            bestPoint = currentPoint;
        }

        // now generate a random jump
        hyFloat distance = jumpMagnitude*genrand_int32()/(hyFloat)RAND_MAX_32, dist1 = 0;

        // generate an arbitrary point in the allowed domain


        for (i=0; i<indexInd.lLength; i++)
        {
            _Variable*v = LocateVar (indexInd(i));
            jumpVector.Store (0,i,v->GetLowerBound()+(v->GetUpperBound()-v->GetLowerBound())*genrand_int32()/(hyFloat)RAND_MAX_32-currentPoint(0,i));
            dist1+=jumpVector(0,i)*jumpVector(0,i);
        }

        jumpVector *= (distance/sqrt(dist1)); // normalize the jump;

        trialPoint = jumpVector;
        trialPoint += currentPoint;


    }       */
    HandleApplicationError("Simulated Annealing is yet to be implemented. Sorry about that.");
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::RescanAllVariables (void)
{
    indexCat.Clear ();
    indexDep.Clear ();
    indexInd.Clear ();
    computationalResults.Clear();
    indVarsByPartition.Clear  ();
    depVarsByPartition.Clear  ();
    ScanAllVariables();
}

//_______________________________________________________________________________________
long    _LikelihoodFunction::DependOnTree (_String const & treeName) const {
    long f = LocateVarByName (treeName);
    if (f>=0L) {
        return theTrees.Find(variableNames.GetXtra(f));
    }
    return -1L;
}

//_______________________________________________________________________________________
long    _LikelihoodFunction::DependOnDS (long ID) const {

    void * data_set_pointer = dataSetList.GetItem (ID);

    for (long k = 0L; k<theDataFilters.lLength; k++) {
      if (GetIthFilter (k) -> GetData() == data_set_pointer) {
        return k;
      }
    }
    return -1L;
}

//_______________________________________________________________________________________
long    _LikelihoodFunction::DependOnModel (_String const& modelTitle) const {
    if (modelTitle.nonempty()) {
        long modelIndex = FindModelName(modelTitle);
        if (modelIndex != HY_NO_MODEL) {
            for (long k=0; k<theTrees.lLength; k++) {
                _TreeIterator ti (GetIthTree(k), _HY_TREE_TRAVERSAL_POSTORDER);
                while (_CalcNode* iterator = ti.Next()) {
                    if (iterator->GetModelIndex() == modelIndex) {
                        return k;
                    }

                }
            }
        }
    }
    return -1;
}

//_______________________________________________________________________________________
_CategoryVariable*  _LikelihoodFunction::FindCategoryVar (long index)
{
    _CategoryVariable* res = nil;
    if (index>=0 && index<blockDependancies.lLength) {
        res = (_CategoryVariable*)LocateVar(indexCat(HighestBit(blockDependancies.lData[index])));
    }
    return res;
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::ScanAllVariables (void) {
    _SimpleList   allVariables,
                  covCat,
                  cpCat,
                  treeSizes,
                  rankVariablesSupp;


    _AVLListX    rankVariables (&rankVariablesSupp);
    _AVLList     avl (&allVariables);
    
    // scan parameters that may be present in frequency vectors and compute tree sizes
    // rank these global variables as influential
    
    theProbabilities.Each([this, &treeSizes,&rankVariables,&avl] (long var_index, unsigned long index) -> void {
        long iNodeCount, lNodeCount;
        GetIthTree (index)->EdgeCount(iNodeCount, lNodeCount);
        treeSizes << (iNodeCount + lNodeCount);
        LocateVar (var_index)->ScanForVariables (avl,true,&rankVariables,treeSizes.GetElement(-1) << 16);
    });
    

    if (computingTemplate) {
        computingTemplate->ScanFForVariables (avl,true,false,true, false, &rankVariables, treeSizes.Sum() << 16);
    }

    avl.ReorderList();

    if (templateKind<0) { // remove dependance of the function from the HMM model specifier
        allVariables.Delete (allVariables.Find(-templateKind-1));
        rankVariables.Delete ((BaseRef)(-templateKind-1));
    }

    _AVLList iia (&indexInd),
             iid (&indexDep);

    for (unsigned  long i=0UL; i<allVariables.countitems(); i++) {
        long variableIndex = allVariables.get(i);
        _Variable* theV = (_Variable*)LocateVar(variableIndex);
        if (theV->IsCategory()) {
            if (((_CategoryVariable*)theV)->IsUncorrelated()) {
                if (((_CategoryVariable*)theV)->IsConstantOnPartition()) {
                    indexCat << variableIndex;
                } else {
                    cpCat << variableIndex;
                }
            } else {
                covCat<< variableIndex;
            }
            continue;
        }
        if (theV->IsIndependent()) {
            iia.Insert ((BaseRef)variableIndex);
        } else {
            iid.Insert ((BaseRef)variableIndex);
        }
    }
    
    theTrees.Each ([this, &iia, &iid, &rankVariables, &treeSizes] (long tree_index, unsigned long index) -> void {
        _TheTree * tree_var = (_TheTree*)LocateVar(tree_index);
        tree_var->ScanAndAttachVariables ();
        tree_var->ScanForGVariables (iia, iid,&rankVariables, treeSizes.GetElement (index) << 16);
        tree_var->ScanContainerForVariables  (iia, iid, &rankVariables, 1L + treeSizes.GetElement (index));
        tree_var->ScanForDVariables (iid, iia);
        tree_var->SetUp();
    });
    

    /*for (unsigned long i=0; i<theTrees.lLength; i++) {
        ((_TheTree*)(LocateVar(theTrees(i))))->ScanAndAttachVariables ();
        ((_TheTree*)(LocateVar(theTrees(i))))->ScanForGVariables (iia, iid,&rankVariables, treeSizes.GetElement (i) << 16);
    }

    for (unsigned long i=0; i<theTrees.lLength; i++) {
        _TheTree * cT = ((_TheTree*)(LocateVar(theTrees(i))));
        cT->ScanContainerForVariables  (iia, iid, &rankVariables, 1 + treeSizes.GetElement (i));
        cT->ScanForDVariables (iid, iia);
        cT->SetUp();
    }*/

    iia.ReorderList ();
    iid.ReorderList ();
 


    /* mod 10/28/2005; Need to make two passes on trees; one to check for all category variables defined on the tree
       and second to lay out dependency bits (after all the category variables have been concatenated) */

    bool haveHMM                 = false,
         haveConstantOnPartition = false;

    for (unsigned long i=0; i<theTrees.lLength; i++) {
        _SimpleList localCategVars;
        {
            _AVLList ca (&localCategVars);
            ((_TheTree*)(LocateVar(theTrees(i))))->ScanForCVariables (ca);
            ca.ReorderList();
        }

        for (unsigned long i2 = 0; i2 < localCategVars.lLength; i2++) {
            _CategoryVariable * cvRef = (_CategoryVariable*)LocateVar (localCategVars.lData[i2]);
            haveHMM                     = haveHMM || cvRef->IsHiddenMarkov();
            haveConstantOnPartition     = haveConstantOnPartition || cvRef->IsConstantOnPartition();
            if (cvRef->IsUncorrelated()) {
                if (cvRef->IsConstantOnPartition()) {
                    indexCat >> cvRef->get_index();
                } else {
                    cpCat >> cvRef->get_index();
                }
            } else {
                covCat >> cvRef->get_index();
            }
        }
    }


    if ((haveHMM || haveConstantOnPartition) && templateKind != _hyphyLFComputationalTemplateNone && templateKind != _hyphyLFComputationalTemplateByPartition) {
        HandleApplicationError ("Non-partition based computational templates in likelihood functions cannot be combined with dependence on HMM or constant-on-partition random variables.");
        return;
    }

    indexCat << cpCat;
    indexCat << covCat;

    if (indexCat.lLength>sizeof(long)*8) {
        HandleApplicationError (_String ("The number of category variables exceeded ") &_String((long)sizeof(long)*8) );
        return;
    }

    for (unsigned long i=0; i<theTrees.lLength; i++) {
        _SimpleList categVars;
        {
            _AVLList ca (&categVars);
            ((_TheTree*)(LocateVar(theTrees(i))))->ScanForCVariables (ca);
            ca.ReorderList();
        }

        if (categVars.lLength) {
            long dep = 0;
            for (long k=categVars.lLength-1; k>=0; k--) {
                SetNthBit(dep,indexCat.Find (categVars(k)));
            }
            blockDependancies<<dep;
        } else {
            blockDependancies<<0;
        }
    }

    if (indexCat.lLength) // scan global variables upon which category variables depend
        // which may have been missed
        // indexInd will be sorted at this point
    {
        _SimpleList sortedCats (indexCat);
        sortedCats.Sort();
        if (sortedCats.CountCommonElements(indexInd,true)) {
            _SimpleList newList;
            newList.Subtract (indexInd,sortedCats);
            indexInd.Clear();
            indexInd.Duplicate(&newList);
        }

        {
            _SimpleList   auxL;
            _AVLList      auxa (&auxL);

            for (unsigned long i=0; i<indexCat.lLength; i++) {
                _CategoryVariable *theCV = (_CategoryVariable*)LocateVar(indexCat(i));
                theCV->ScanForGVariables (auxa);
            }

            auxa.ReorderList();
            if (auxL.lLength) {
                _SimpleList nl;
                nl.Union (indexInd,auxL);
                if (nl.lLength > indexInd.lLength) {
                    indexInd.Clear();
                    indexInd.Duplicate (&nl);
                }
            }
        }
    }

    hyFloat l = DEFAULTLOWERBOUND*(1.0-kMachineEpsilon),
               u = DEFAULTUPPERBOUND*(1.0-kMachineEpsilon);



    for (unsigned long i=0; i<indexInd.lLength; i++) {
        _Variable *_cv = GetIthIndependentVar(i);
        if (_cv->GetLowerBound()<=l) {
            _cv->SetBounds(DEFAULTPARAMETERLBOUND,_cv->GetUpperBound());
        }
        if (_cv->GetUpperBound()>=u) {
            _cv->SetBounds(_cv->GetLowerBound(),DEFAULTPARAMETERUBOUND);
        }
        //printf ("%s -> %d\n", _cv->theName->sData, rankVariables.GetXtra(rankVariables.Find ((BaseRef)indexInd.lData[i])));
    }

    for (unsigned  long i=0; i<indexDep.lLength; i++) {
        _Variable *_cv = GetIthDependentVar(i);
        if (_cv->GetLowerBound()<=l) {
            _cv->SetBounds(DEFAULTPARAMETERLBOUND,_cv->GetUpperBound());
        }
        if (_cv->GetUpperBound()>=u) {
            _cv->SetBounds(_cv->GetLowerBound(),DEFAULTPARAMETERUBOUND);
        }
    }

    _SimpleList pidx (1,0,0);
    for (unsigned long p = 0; p < theTrees.lLength; p++) {
        pidx.lData[0] = p;
        _SimpleList iv,dv,cv;
        ScanAllVariablesOnPartition (pidx, iv, dv, cv, true);
        indVarsByPartition && & iv;
        depVarsByPartition && & dv;
    }

    RankVariables(&rankVariables);

    //for (long iid = 0; iid < indexInd.lLength; iid++)
    //  printf ("%ld: %s\n", iid, LocateVar(indexInd.lData[iid])->GetName()->sData);

}

//_______________________________________________________________________________________
void    _LikelihoodFunction::ScanAllVariablesOnPartition (_SimpleList& pidx, _SimpleList& iind, _SimpleList& idep, _SimpleList &icat, bool treeOnly)
{
    _SimpleList   allVariables,
                  covCat,
                  cpCat;


    if (!treeOnly) {
        _AVLList avl (&allVariables);
        for (long i = 0; i < pidx.lLength; i++) {
            LocateVar(theProbabilities(pidx(i)))->ScanForVariables (avl,true);
        }

        if (computingTemplate) {
            computingTemplate->ScanFForVariables (avl,true,false,true);
        }

        avl.ReorderList();
    }

    if (treeOnly == false && templateKind<0) {
        allVariables.Delete (allVariables.Find(-templateKind-1));
    }

    {
        _AVLList iia (&iind),
                 iid (&idep);

        if (treeOnly == false) {
            for (long i=0; i<allVariables.lLength; i++) {
                _Variable* theV = ((_Variable*)LocateVar(allVariables(i)));
                if (theV->IsCategory()) {
                    if (((_CategoryVariable*)theV)->IsUncorrelated()) {
                        if (((_CategoryVariable*)theV)->IsConstantOnPartition()) {
                            icat << allVariables(i);
                        } else {
                            cpCat << allVariables(i);
                        }
                    } else {
                        covCat<<allVariables(i);
                    }
                    continue;
                }
                if (theV->IsIndependent()) {
                    iia.Insert ((BaseRef)allVariables(i));
                } else {
                    iid.Insert ((BaseRef)allVariables(i));
                }
            }

            indexCat << cpCat;
            indexCat << covCat;
        }

        for (long i2=0; i2<pidx.lLength; i2++) {
            _TheTree * cT = (_TheTree*)(LocateVar(theTrees.lData[pidx.lData[i2]]));
            cT->ScanForGVariables (iia, iid);
        }

        for (long i=0; i<pidx.lLength; i++) {
            _TheTree * cT = (_TheTree*)(LocateVar(theTrees.lData[pidx.lData[i]]));
            cT->ScanContainerForVariables (iia,iid);
            cT->ScanForDVariables (iid, iia);
        }

        iia.ReorderList ();
        iid.ReorderList ();
    }


    for (long i=0; i<pidx.lLength; i++) {
        _SimpleList categVars;
        _AVLList ca (&categVars);
        ((_TheTree*)(LocateVar(theTrees.lData[pidx.lData[i]])))->ScanForCVariables (ca);
        ca.ReorderList();

        if (categVars.lLength)
            for (long k=categVars.lLength-1; k>=0; k--) {
                long f = icat.Find (categVars(k));
                if (f==-1) {
                    icat<<categVars(k);
                }
            }
    }

    if (icat.lLength) {
        for (long i=0; i<iind.lLength; i++) {
            long t = icat.Find (iind.lData[i]);
            if (t>=0) {
                iind.Delete(i);
                i--;
            }
        }
        {
            _SimpleList   auxL;
            _AVLList      auxa (&auxL);

            for (long i=0; i<icat.lLength; i++) {
                _CategoryVariable *theCV = (_CategoryVariable*)LocateVar(icat(i));
                theCV->ScanForGVariables (auxa);
            }

            auxa.ReorderList();
            if (auxL.lLength) {
                _SimpleList nl;
                nl.Union (iind,auxL);
                if (nl.lLength > iind.lLength) {
                    iind.Clear();
                    iind.Duplicate (&nl);
                }
            }
        }
    }
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::UpdateIndependent (long index, bool purgeResults, _SimpleList* whichList, _SimpleList* secondList) {
    _SimpleList * theList = &indexInd;
    if (whichList) {
        theList = whichList;
    } else {
        secondList = &indexDep;
    }

    long f = theList->Find (index);
    if (f!=-1) {
        theList->Delete(f);
        (*secondList)<<index;

        _SimpleList   newVars;
        {
            _AVLList  al (&newVars);
            LocateVar(index)->ScanForVariables(al,true);
            al.ReorderList();
        }

        for (f=0; f<newVars.countitems(); f++) {
            _Variable* cv = LocateVar(newVars.lData[f]);
            if (cv->IsIndependent()&& theList->Find(newVars.lData[f])==-1) {
                (*theList) << newVars.lData[f];
            }
        }

        if (theList!=whichList) {
            for (f = 0; f < indVarsByPartition.lLength; f++) {
                UpdateIndependent (index, false, (_SimpleList*)indVarsByPartition(f), (_SimpleList*)depVarsByPartition(f));
            }
        }

        if (purgeResults) {
            computationalResults.Clear();
        }
    }

}

//_______________________________________________________________________________________
void    _LikelihoodFunction::UpdateDependent (long index) {
    long f = indexDep.Find (index);
    if (f >= 0L) {
        indexDep.Delete(f);
        indexInd<<index;
        for (unsigned long k = 0; k<depVarsByPartition.lLength; k++) {
            f = ((_SimpleList*)depVarsByPartition(k))->Find(index);
            if (f >= 0) {
                ((_SimpleList*)depVarsByPartition(k))->Delete(f);
                (*(_SimpleList*)indVarsByPartition(k)) << index;
            }
        }
    }

}


//_______________________________________________________________________________________

void    _LikelihoodFunction::Cleanup (void) {
    //DeleteCaches();

    Clear();
    DeleteObject (parameterValuesAndRanges);

#ifdef MDSOCL
	for (int i = 0; i < theTrees.lLength; i++)
	{
		OCLEval[i].~_OCLEvaluator();
	}
	delete [] OCLEval;
#endif
}


//_______________________________________________________________________________________

void    _LikelihoodFunction::DeleteCaches (bool all)
{
    if (all) {
        DeleteObject (siteResults);
        siteResults = nil;
        DeleteObject (bySiteResults);
        bySiteResults = nil;
    }

    conditionalTerminalNodeLikelihoodCaches.Clear();
    cachedBranches.Clear();
    siteCorrections.Clear();
    siteCorrectionsBackup.Clear();
    siteScalerBuffer.Clear();

//  computedLocalUpdatePolicy.Clear();
//  treeTraversalMasks.Clear();
//  matricesToExponentiate.Clear();
//  overallScalingFactors.Clear();
//  localUpdatePolicy.Clear();

    if (conditionalInternalNodeLikelihoodCaches) {
        for (long k = 0; k < theTrees.lLength; k++)
            if (conditionalInternalNodeLikelihoodCaches[k]) {
                delete [] (conditionalInternalNodeLikelihoodCaches[k]);
            }
        delete [] conditionalInternalNodeLikelihoodCaches ;
        conditionalInternalNodeLikelihoodCaches = nil;
    }
    if (branchCaches) {
        for (long k = 0; k < theTrees.lLength; k++)
            if (branchCaches[k]) {
                delete [] branchCaches[k];
            }
        delete [] branchCaches;
        branchCaches = nil;
    }
    if (conditionalTerminalNodeStateFlag) {
        for (long k = 0; k < theTrees.lLength; k++)
            if (conditionalTerminalNodeStateFlag[k]) {
                delete [] conditionalTerminalNodeStateFlag[k];
            }
        delete [] conditionalTerminalNodeStateFlag;
        conditionalTerminalNodeStateFlag = nil;
    }
    if (siteScalingFactors) {
        for (long k = 0; k < theTrees.lLength; k++)
            if (siteScalingFactors[k]) {
                delete [] siteScalingFactors[k];
            }
        delete [] siteScalingFactors;
        siteScalingFactors = nil;
    }
}


//_______________________________________________________________________________________

void    _LikelihoodFunction::Setup (bool check_reversibility)
{
    hyFloat kp       = 0.0;
    //RankVariables();
    checkParameter      (useFullMST,kp,0.0);

    if (kp>.5 && !mstCache) {
        mstCache = new MSTCache;
    }

    if (theTrees.lLength==optimalOrders.lLength) {
        //check to see if we need to recompute the
        // optimal summation order
        checkParameter (keepOptimalOrder,kp,0.0);
        if (kp) {
            for (unsigned long i=0; i<theTrees.lLength; i++) {
                _SimpleList*    s = (_SimpleList*)optimalOrders(i),
                                *   l = (_SimpleList*)leafSkips(i);

                _DataSetFilter const * df   = GetIthFilter(i);
                _Matrix       *glFreqs      = (_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue();
                _TheTree      *t            = ((_TheTree*)LocateVar(theTrees.lData[i]));

                t->InitializeTreeFrequencies (glFreqs, true);
                if (s->lLength!=df->GetPatternCount()) {
                    s->Clear();
                    l->Clear();
                    OptimalOrder (i,*s);
                    df->MatchStartNEnd (*s,*l);
                }
            }
            return;
        }
    }

    optimalOrders.Clear();
    leafSkips.Clear();

    treeTraversalMasks.Clear();

    if (!check_reversibility) {
        if (canUseReversibleSpeedups.countitems() != theTrees.lLength) {
            check_reversibility = true;
            canUseReversibleSpeedups.Clear();
        }
    } else {
        canUseReversibleSpeedups.Clear();
    }
    _SimpleList alreadyDoneModelsL;
    _AVLListX   alreadyDoneModels (&alreadyDoneModelsL);

    hyFloat assume_rev = hy_env::EnvVariableGetNumber(hy_env::assume_reversible, 0.);
 
    for (unsigned long i=0UL; i<theTrees.lLength; i++) {
        _Matrix         *glFreqs = GetIthFrequencies(i);
        _DataSetFilter const* df      = GetIthFilter(i);
        _TheTree        *t       = GetIthTree (i);
        t->InitializeTreeFrequencies (glFreqs, true);
        _SimpleList        *s = new _SimpleList,
        *l = new _SimpleList;

        treeTraversalMasks.AppendNewInstance(new _SimpleList (t->GetINodeCount() * df->GetPatternCount() / _HY_BITMASK_WIDTH_ + 1,0,0));
        OptimalOrder      (i,*s);
        df->MatchStartNEnd(*s,*l);
        optimalOrders.AppendNewInstance(s);
        leafSkips.AppendNewInstance(l);

        if (check_reversibility) {
          _SimpleList treeModels;
          t->CompileListOfModels (treeModels);
          bool isReversiblePartition = true;
          if (assume_rev > 0.5) {
              ReportWarning (_String ("Partition ") & long(i) & " is ASSUMED to have a reversible model");
          } else {
              if (assume_rev < -0.5) {
                isReversiblePartition = false;
                ReportWarning (_String ("Partition ") & long(i) & " is ASSUMED to have a non-reversible model");
              } else {
                for (unsigned long m = 0; m < treeModels.lLength && isReversiblePartition; m++) {
                    long alreadyDone = alreadyDoneModels.Find ((BaseRef)treeModels.lData[m]);
                    if (alreadyDone>=0) {
                        alreadyDone = alreadyDoneModels.GetXtra (alreadyDone);
                    } else {
                        alreadyDone = IsModelReversible (treeModels.lData[m]);
                        alreadyDoneModels.Insert ((BaseRef)treeModels.lData[m], alreadyDone);
                    }
                    isReversiblePartition = isReversiblePartition && alreadyDone;
                }
                ReportWarning (_String ("Partition ") & (long)i & " reversible model flag computed as " & (long)isReversiblePartition);
              }
          }
          canUseReversibleSpeedups << isReversiblePartition;
        }
    }

}

//_______________________________________________________________________________________

bool    _LikelihoodFunction::HasPartitionChanged (long index) {

    return ListAny (*(_SimpleList*)indVarsByPartition(index),
                    [] (const long value, const unsigned long index) -> bool {
                        return LocateVar(value)->HasChanged();
                       }
                    );

}

//#define _HY_GPU_EXAMPLE_CALCULATOR

//_______________________________________________________________________________________

hyFloat  _LikelihoodFunction::ComputeBlock (long index, hyFloat* siteRes, long currentRateClass, long branchIndex, _SimpleList * branchValues)
// compute likelihood over block index i
/*
    to optimize
        -no need to recurse the entire tree to decide if it had changed;
         should cache variable->block dependancies for rapid lookup
*/
{

    // set up global matrix frequencies

    _SimpleList               *sl = (_SimpleList*)optimalOrders.GetItem(index);

    _Matrix                   *glFreqs          = GetIthFrequencies(index);
    _DataSetFilter            const *df         = GetIthFilter(index);
    _TheTree                   *t               = GetIthTree(index);
    bool                       canClear         = true,
                               rootFreqsChange  = forceRecomputation?true:glFreqs->HasChanged();

    if (currentRateClass >=0 && t->HasForcedRecomputeList()) {
        canClear = TotalRateClassesForAPartition(index) == currentRateClass+1;
    }

    t->InitializeTreeFrequencies          ((_Matrix*)glFreqs->ComputeNumeric());

    usedCachedResults                   = false;

    if (computingTemplate&&templateKind) {
        if (!(forceRecomputation||!siteArrayPopulated||HasPartitionChanged(index)||rootFreqsChange)) {
            usedCachedResults = true;
            //printf ("\n[CACHED]\n");
            return            -1e300;
        }
    } else {
        if (!forceRecomputation && computationalResults.get_used()==optimalOrders.lLength && !siteRes && !HasPartitionChanged(index) && !rootFreqsChange) {
            usedCachedResults = true;
            //printf ("\n[CACHED]\n");
            return computationalResults.theData[index];
        }
    }

    if (conditionalInternalNodeLikelihoodCaches) {

        long        catID            = siteRes?currentRateClass:-1;

        if (conditionalInternalNodeLikelihoodCaches[index]) {
            // not a 2 sequence analysis
            long blockID    = df->GetPatternCount()*t->GetINodeCount(),
                 patternCnt = df->GetPatternCount();

            _SimpleList         *tcc  = (_SimpleList*)treeTraversalMasks(index);

            hyFloat          *inc  = (currentRateClass<1)?conditionalInternalNodeLikelihoodCaches[index]:
                                        conditionalInternalNodeLikelihoodCaches[index] + currentRateClass*df->GetDimension()*blockID,
                                *ssf  = (currentRateClass<1)?siteScalingFactors[index]: siteScalingFactors[index] + currentRateClass*blockID,
                                *bc   = (currentRateClass<1)?branchCaches[index]: (branchCaches[index] + currentRateClass*patternCnt*df->GetDimension()*2);

            long  *scc = nil,
                  *sccb = nil;

            if (siteRes) {
                sccb = ((_SimpleList*)siteCorrectionsBackup(index))->lData + ((currentRateClass<1)?0:patternCnt*currentRateClass);
                scc =  ((_SimpleList*)siteCorrections(index))->lData       + ((currentRateClass<1)?0:patternCnt*currentRateClass);
            }

            _SimpleList changedBranches, *branches;
            _List       changedModels,   *matrices;
            long        doCachedComp     = 0,     // whether or not to use a cached branch calculation when only one
                        // local tree parameter is being adjusted at a time

                        ciid          = MAX(0,currentRateClass),
                        *cbid            = &(((_SimpleList*)cachedBranches(index))->lData[ciid]);

            if (computedLocalUpdatePolicy.lLength && branchIndex < 0) {
                branches = (_SimpleList*)(*((_List*)localUpdatePolicy(index)))(ciid);
                matrices = (_List*)      (*((_List*)matricesToExponentiate(index)))(ciid) ;

                long nodeID = ((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid];
                if (nodeID == 0 || nodeID == 1) {
                    long snID = -1;

                    if (nodeID == 1) {
                        if (matrices->lLength == 2) {
                            branches->Clear();
                            matrices->Clear();

                            snID = t->DetermineNodesForUpdate          (*branches, matrices,catID,*cbid,canClear);
                        }
                    } else {
                        snID = t->DetermineNodesForUpdate          (*branches, matrices,catID,*cbid,canClear);
                    }

#ifdef _UBER_VERBOSE_LF_DEBUG
                    fprintf (stderr, "\nCached %s (nodeID = %lD)/New %s (touched matrices %ld) Eval id = %ld\n", *cbid >= 0 ? t->GetNodeFromFlatIndex (*cbid)->GetName()->getStr() : "None", nodeID, snID >= 0 ? t->GetNodeFromFlatIndex (snID)->GetName()->getStr() : "None", matrices->lLength, likeFuncEvalCallCount);
#endif
                    if (snID != *cbid) {
                        RestoreScalingFactors (index, *cbid, patternCnt, scc, sccb);
                        *cbid = -1;
                        if (snID >= 0 && canUseReversibleSpeedups.lData[index]) {
                            ((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid] = snID+3;
                            doCachedComp = -snID-1;
                        } else {
                            ((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid] = nodeID + 1;
                        }
                    } else {
                    // 20120718: SLKP added this branch to reuse the old cache if the branch that is being computed
                    // is the same as the one cached last time, e.g. sequentially iterating through all local parameters
                    // of a given branch.
                        if (snID >= 0 && canUseReversibleSpeedups.lData[index]) {
                            doCachedComp = ((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid] = snID+3;
                          } else {
                            ((_SimpleList*)computedLocalUpdatePolicy(index))->lData[ciid] = nodeID + 1;
                        }
                    }
                } else {
                    doCachedComp = nodeID;
                }

            } else {
                RestoreScalingFactors       (index, *cbid, patternCnt, scc, sccb);


                t->DetermineNodesForUpdate  (changedBranches,&changedModels,catID,(branchIndex >=0 )?
                                             (branchIndex<t->GetINodeCount()?branchIndex+t->GetLeafCount():branchIndex):*cbid,canClear);
                *cbid                       = -1;
                branches                    = &changedBranches;
                matrices                    = &changedModels;
            }

            if (evalsSinceLastSetup == 0) {
                branches->Populate (t->GetINodeCount()+t->GetLeafCount()-1,0,1);
            }

#ifdef _UBER_VERBOSE_LF_DEBUG
            fprintf (stderr, "%d matrices, %d branches marked for rate class %d\n", matrices->lLength, branches->lLength, catID);
            if (matrices->lLength == 0) {
                fprintf (stderr, "Hmm\n");
            }
#endif
            if (matrices->lLength) {
                t->ExponentiateMatrices(*matrices, GetThreadCount(),catID);
            }



            hyFloat sum  = 0.;

            if (doCachedComp >= 3) {
#ifdef _UBER_VERBOSE_LF_DEBUG
                fprintf (stderr, "CACHE compute branch %d\n",doCachedComp-3);
#endif
                sum = t->ComputeLLWithBranchCache (*sl,
                                                   doCachedComp-3,
                                                   bc,
                                                   df,
                                                   0,
                                                   df->GetPatternCount (),
                                                   catID,
                                                   siteRes)
                      - _logLFScaler * overallScalingFactors.lData[index];
                return sum;
            }

            long np = 1;
#ifdef _OPENMP
            np           = MIN(GetThreadCount(),omp_get_max_threads());
#endif
            long sitesPerP    = df->GetPatternCount() / np + 1L;

#ifdef _UBER_VERBOSE_LF_DEBUG
                fprintf (stderr, "NORMAL compute lf \n");
#endif

            hyFloat* thread_results = new hyFloat[np];
            #pragma omp  parallel for default(shared) schedule(static,1) private(blockID) num_threads (np) if (np>1)
              for (blockID = 0; blockID < np; blockID ++) {
                thread_results[blockID] = t->ComputeTreeBlockByBranch (*sl,
                                                    *branches,
                                                    tcc,
                                                    df,
                                                    inc,
                                                    conditionalTerminalNodeStateFlag[index],
                                                    ssf,
                                                    (_Vector*)conditionalTerminalNodeLikelihoodCaches(index),
                                                    overallScalingFactors.lData[index],
                                                    blockID * sitesPerP,
                                                    (1+blockID) * sitesPerP,
                                                    catID,
                                                    siteRes,
                                                    scc,
                                                    branchIndex,
                                                    branchIndex >= 0 ? branchValues->lData: nil);
              }


            if (np > 1) {
              hyFloat correction = 0.;
              for (blockID = 0; blockID < np; blockID ++)  {
                thread_results[blockID] -= correction;
                hyFloat temp_sum = sum +  thread_results[blockID];
                correction = (temp_sum - sum) - thread_results[blockID];
                sum = temp_sum;
              }


            } else {
              sum = thread_results[0];

            }

            delete [] thread_results;
            /*
            #pragma omp  parallel for default(shared) schedule(static,1) private(blockID) num_threads (np) reduction(+:sum) if (np>1)
            for (blockID = 0; blockID < np; blockID ++) {
                sum += t->ComputeTreeBlockByBranch (*sl,
                                                    *branches,
                                                    tcc,
                                                    df,
                                                    inc,
                                                    conditionalTerminalNodeStateFlag[index],
                                                    ssf,
                                                    (_GrowingVector*)conditionalTerminalNodeLikelihoodCaches(index),
                                                    overallScalingFactors.lData[index],
                                                    blockID * sitesPerP,
                                                    (1+blockID) * sitesPerP,
                                                    catID,
                                                    siteRes,
                                                    scc,
                                                    branchIndex,
                                                    branchIndex >= 0 ? branchValues->lData: nil);
            }*/

            sum -= _logLFScaler * overallScalingFactors.lData[index];
            

            if (doCachedComp < 0) {
                //printf ("Cache check in %d %d\n", doCachedComp, overallScalingFactors[index]);
                doCachedComp = -doCachedComp-1;
                //printf ("Set up %d\n", doCachedComp);
                *cbid = doCachedComp;


                overallScalingFactorsBackup.lData[index] = overallScalingFactors.lData[index];
                if (sccb)
                    for (long recoverIndex = 0; recoverIndex < patternCnt; recoverIndex++) {
                        sccb[recoverIndex] = scc[recoverIndex];
                    }

                /*for (unsigned long p_id = 0; p_id < indexInd.lLength; p_id++) {
                  printf ("%ld %s = %15.12g\n", p_id, GetIthIndependentVar(p_id)->GetName()->sData, (*parameterValuesAndRanges)(p_id,0));
                }*/

                 #pragma omp  parallel for default(shared) schedule(static,1) private(blockID) num_threads (np) if (np>1)
                for (blockID = 0; blockID < np; blockID ++) {
                    t->ComputeBranchCache (*sl,doCachedComp, bc, inc, df,
                                           conditionalTerminalNodeStateFlag[index],
                                           ssf,
                                           scc,
                                           (_Vector*)conditionalTerminalNodeLikelihoodCaches(index),
                                           overallScalingFactors.lData[index],
                                           blockID * sitesPerP,
                                           (1+blockID) * sitesPerP,
                                           catID,tcc,siteRes);
                }

                // check results

                if (sum > -INFINITY) {
                   hyFloat checksum = t->ComputeLLWithBranchCache (*sl,
                                                     doCachedComp,
                                                     bc,
                                                     df,
                                                     0,
                                                     df->GetPatternCount (),
                                                     catID,
                                                     siteRes)
                  - _logLFScaler * overallScalingFactors.lData[index];

                  if (fabs ((checksum-sum)/sum) > 1.e-10 * df->GetPatternCount ()) {
                    /*hyFloat check2 = t->ComputeTreeBlockByBranch (*sl,
                                                                     *branches,
                                                                     tcc,
                                                                     df,
                                                                     inc,
                                                                     conditionalTerminalNodeStateFlag[index],
                                                                     ssf,
                                                                     (_GrowingVector*)conditionalTerminalNodeLikelihoodCaches(index),
                                                                     overallScalingFactors.lData[index],
                                                                     0,
                                                                     df->GetPatternCount(),
                                                                     catID,
                                                                     siteRes,
                                                                     scc,
                                                                     branchIndex,
                                                                     branchIndex >= 0 ? branchValues->lData: nil);*/

                    _String* node_name =   GetIthTree (index)->GetNodeFromFlatIndex(doCachedComp)->GetName();

                    _TerminateAndDump (_String("Internal error in ComputeBranchCache (branch ") & *node_name &
                                         +                                       " ) reversible model cached likelihood = "& _String (checksum, "%20.16g") & ", directly computed likelihood = " & _String (sum, "%20.16g") &
                                         +                                       ". This is most likely because a non-reversible model was incorrectly auto-detected (or specified by the model file in environment variables).");

                     return -INFINITY;
                  }
                }

                // need to update siteRes when computing cache and changing scaling factors!
            }
            return sum;
        } else if (conditionalTerminalNodeStateFlag[index] || !df->IsNormalFilter()) {
            // two sequence analysis

            _SimpleList bc;
            _List       mc;
            t->DetermineNodesForUpdate         (bc, &mc);
            if (mc.lLength) {
                t->ExponentiateMatrices(mc, GetThreadCount(),catID);
            }

            //branchedGlobalCache.Clear(false);
            //matricesGlobalCache.Clear(false);
            if (df->IsNormalFilter())
                return t->ComputeTwoSequenceLikelihood (*sl,
                                                        df,
                                                        conditionalTerminalNodeStateFlag[index],
                                                        (_Vector*)conditionalTerminalNodeLikelihoodCaches(index),
                                                        0,
                                                        df->GetPatternCount (),
                                                        catID,
                                                        siteRes);
            else {
                return t->Process3TaxonNumericFilter ((_DataSetFilterNumeric*)df);
            }

        }

    } else {
        HandleApplicationError ("Dude -- lame! No cache. I can't compute like that with the new LF engine.");
    }

    return 0.0;
}

//_______________________________________________________________________________________
long        _LikelihoodFunction::CostOfPath  (_DataSetFilter const* df, _TheTree const* t, _SimpleList& sl, _SimpleList* tcc) const {
    long res = 0L;
    for (long i=1L; i<sl.countitems(); i++) {
        res+=t->ComputeReleafingCost (df,sl.get (i-1L),sl.get(i), tcc, i);
    }
    return res;
}




//_______________________________________________________________________________________

//#define __SORTING_DEBUG
#ifdef   __SORTING_DEBUG
#include "HYChartWindow.h"
#endif

void    countingTraverse (node<long>* startingNode, long& totalLength, long currentSize, long& maxSize, bool add2Size)
{
    if (startingNode->parent) {
        totalLength+=startingNode->in_object;
    }

    if (add2Size) {
        currentSize++;
    }

    if (currentSize>maxSize) {
        maxSize = currentSize;
    }

    for (long k=1; k<startingNode->nodes.length; k++) {
        countingTraverse (startingNode->go_down(k), totalLength, currentSize, maxSize, true);
    }

    if (startingNode->nodes.length) {
        countingTraverse (startingNode->go_down(startingNode->nodes.length), totalLength, currentSize, maxSize, false);
    }
}

//_______________________________________________________________________________________

void    countingTraverseArbRoot (node<long>* startingNode, node<long>* childNode, long& totalLength, long currentSize, long& maxSize)
{
    if (childNode)
        for (long k=1; k<=startingNode->nodes.length; k++) {
            node<long>* cNode = startingNode->go_down(k);
            if (cNode!=childNode) {
                countingTraverse (cNode, totalLength, currentSize, maxSize, true);
            }
        }
    else
        for (long k=1; k<=startingNode->nodes.length; k++) {
            countingTraverse (startingNode->go_down(k), totalLength, currentSize, maxSize, true);
        }

    if (startingNode->parent) {
        totalLength+=startingNode->in_object;
        countingTraverseArbRoot (startingNode->parent, startingNode, totalLength, currentSize, maxSize);
    }
}
//_______________________________________________________________________________________

long    findAvailableSlot (_SimpleList& cs, long& ff)
{
    for (long k = ff; k<cs.lLength; k++)
        if (cs.lData[k] == -1) {
            ff = k+1;
            return k;
        }
    {
        for (long k = 0; k<ff; k++)
            if (cs.lData[k] == -1) {
                ff = k+1;
                return k;
            }
    }
    cs << -1;
    ff = 0;
    return cs.lLength-1;
}

//_______________________________________________________________________________________

void    setComputingArrays (node<long>* startingNode, node<long>* childNode, _SimpleList& co, _SimpleList& so, _SimpleList &cs, _SimpleList& ro, _SimpleList& po, long& ff)
{
    bool isFirstNode = (co.lLength == 0);

    co << startingNode->in_object;

    if (isFirstNode || startingNode->nodes.length) {
        long id = findAvailableSlot(cs,ff);
        cs.lData[id] = startingNode->in_object;
        startingNode->in_object = id;
        so << id;
        if (isFirstNode) {
            ro << -1;
            po << -1;
        }
    } else {
        so << -1;
    }

    if (childNode) { // above root node
        ro << cs.lData[childNode->in_object];
        po << (long)childNode;

        if (startingNode->parent) {
            for (long k=1; k<=startingNode->nodes.length; k++) {
                node<long>* cNode = startingNode->go_down(k);
                if (cNode!=childNode) {
                    setComputingArrays (startingNode->go_down(k), nil, co,so,cs,ro,po,ff);
                }
            }
            cs.lData[childNode->in_object] = -1;
            setComputingArrays (startingNode->parent, startingNode, co,so,cs,ro,po,ff);
        } else { // actual tree root
            bool passedUpPath = false;

            for (long k=1; k<=startingNode->nodes.length; k++) {
                node<long>* cNode = startingNode->go_down(k);
                if (cNode!=childNode) {
                    if ((k==startingNode->nodes.length)||((k==startingNode->nodes.length-1)&&(!passedUpPath))) {
                        cs.lData[startingNode->in_object] = -1;
                    }
                    setComputingArrays (startingNode->go_down(k), nil, co,so,cs,ro,po,ff);
                } else {
                    passedUpPath = true;
                }
            }
        }
    } else { // "below" root node
        if ((!isFirstNode)&&startingNode->parent) {
            ro << startingNode->parent->in_object;
            po << (long)startingNode->parent;
        }

        for (long k=1; k<startingNode->nodes.length; k++) {
            setComputingArrays (startingNode->go_down(k), nil, co,so,cs,ro,po,ff);
        }

        if (!isFirstNode)
            if (startingNode->nodes.length) {
                cs.lData[startingNode->in_object] = -1;
            }

        if (startingNode->nodes.length) {
            setComputingArrays (startingNode->go_down(startingNode->nodes.length), nil, co,so,cs,ro,po,ff);
        }

        if (isFirstNode&&startingNode->parent) {
            cs.lData[startingNode->in_object] = -1;
            setComputingArrays (startingNode->parent,startingNode, co,so,cs,ro,po,ff);
        }

    }
}

//_______________________________________________________________________________________

void        _LikelihoodFunction::OptimalOrder    (long index, _SimpleList& sl) {

    _DataSetFilter const* df = GetIthFilter (index);

    long            partition           = -1,
                    totalSites           = 0,
                    completedSites         = 0,
                    startpt,
                    endpt,
                    j,
                    max                 = -1,
                    k,
                    intI,
                    totalLength;



    hyFloat      skipo = 1.0;
    _TheTree        *t = (_TheTree*)LocateVar(theTrees(index));
    checkParameter  (optimizeSummationOrder,skipo,1.0);

    if (!skipo || df->GetPatternCount()==1 || t->IsDegenerate() || !df->IsNormalFilter ()) { // do not optimize
        for (k = 0; k < df->GetPatternCount(); k++) {
            sl<<k;
        }
        return;
    }
    SetStatusLine ("Optimizing data ordering");
    checkParameter (optimizePartitionSize,skipo,0.0);
    totalSites = df->GetPatternCount();
    if (skipo) { //  partition the sequence into smaller subseqs. for optimization
        partition = (long)skipo;
        if ((partition<=0)||(partition>totalSites)) {
            partition = totalSites;
        }
    } else { // revert to default partitioning
        partition = totalSites>1500?1500:totalSites;
    }

    long    vLevel       = VerbosityLevel(),
            globalLength = 0;

    if (vLevel>5) {
        char buffer [128];
        snprintf (buffer, sizeof(buffer),"\nOptimizing Column Order for block %ld", index);
        BufferToConsole (buffer);
    }

    _SimpleList   partitionSites, distances, edges;


    while (completedSites<totalSites) {
        if (totalSites-completedSites<partition) {
            partition = totalSites-completedSites;
        }

        intI = 0; // internal index for partition

        // populate sites allowed

        if (df->GetUnitLength()==1) {
            for (k=completedSites+1; k<completedSites+partition; k++) {
                partitionSites<<k;
                distances<<t->ComputeReleafingCostChar (df, completedSites, k);
                edges<<completedSites;
            }
        } else {
            for (k=completedSites+1; k<completedSites+partition; k++) {
                partitionSites<<k;
                distances<<t->ComputeReleafingCost (df, completedSites, k);
                edges<<completedSites;
            }
        }

        node<long>*   spanningTreeRoot;
        _SimpleList   spanningTreePointers,
                      spanningTreeSites;

        if (mstCache) {
            spanningTreeRoot = new node<long>;
            spanningTreeSites << 0;
            for (k=completedSites+1; k<completedSites+partition; k++) {
                spanningTreePointers << (long)spanningTreeRoot;
                spanningTreeSites << 0;
            }
            spanningTreeSites.lData[0] = (long)spanningTreeRoot;
        }

        sl << completedSites;

        while (intI<partition-1) {
            // search for the shortest branch to add to the tree.
            max = 0x0fffffff;
            long * pl  = distances.lData;

            for (k=0; k<distances.lLength; k++,pl++)
                if (*pl<=max) {
                    max = *pl;
                    endpt = partitionSites.lData[k];
                    startpt = k;
                }

            // delete k from allowed sites
            partitionSites.Delete(startpt);
            distances.Delete(startpt);
            // insert the endpt into the tree before the pt that the edge came from

            node<long>*   spanningTreeNode;
            if (mstCache) {
                spanningTreeNode = new node<long>;
                spanningTreeNode->in_object = max;
                ((node<long>*)spanningTreePointers.lData[startpt])->add_node (*spanningTreeNode);
                spanningTreeSites.lData[endpt-completedSites] = (long)spanningTreeNode;
                spanningTreePointers.Delete(startpt);
            }

            j = sl.Find(edges.lData[startpt],completedSites);
            sl.InsertElement ((BaseRef)endpt,j+1,false,false);
            edges.Delete(startpt);

            // make one more pass and update the distances if needed
            if (df->GetUnitLength()==1) {
                for (k=0; k<distances.lLength; k++) {
                    j = t->ComputeReleafingCostChar (df,endpt, partitionSites.lData[k]);
                    if (j<distances.lData[k]) {
                        distances.lData[k]=j;
                        edges.lData[k] = endpt;

                        if (mstCache) {
                            spanningTreePointers.lData[k] = (long)spanningTreeNode;
                        }
                    }
                }
            } else {
                for (k=0; k<distances.lLength; k++) {
                    j = t->ComputeReleafingCost (df,endpt, partitionSites.lData[k]);
                    if (j<distances.lData[k]) {
                        distances.lData[k]=j;
                        edges.lData[k] = endpt;

                        if (mstCache) {
                            spanningTreePointers.lData[k] = (long)spanningTreeNode;
                        }
                    }
                }
            }

            intI++;
            if (intI%50==0) {
#if !defined __UNIX__ || defined __HEADLESS__
                SetStatusBarValue ((intI+completedSites)*100/totalSites,1.0, 0.0);
#endif
                //yieldCPUTime();
            }
        }



        partitionSites.Clear();
        distances.Clear();
        edges.Clear();

        if (mstCache) {
            // print out the info on the spanning tree
            // keep track of the "deepest path"
            long        maxLevel    = 0;

            totalLength = 0;
            countingTraverse (spanningTreeRoot,totalLength,0,maxLevel,true);

            globalLength += totalLength;


            long        level,
                        nc          =   0;


            node_iterator<long> ni (spanningTreeRoot, _HY_TREE_TRAVERSAL_POSTORDER);

            while (node<long>* iterator = ni.Next()) {
              if (iterator != spanningTreeRoot) {
                long maxLevel2 = 0L;
                     totalLength = 0L;
                countingTraverseArbRoot (spanningTreeRoot, nil, totalLength, 1, maxLevel2);
                maxLevel = MIN (maxLevel, maxLevel2);
              }
            }

            _SimpleList   computingOrder,
                          storageOrder,
                          cacheSlots,
                          referenceOrder,
                          parentOrder;

            for (level=0; level<maxLevel; level++) {
                cacheSlots << -1;
            }

            for (level=0; level<spanningTreeSites.lLength; level++) {
                ((node<long>*)spanningTreeSites.lData[level])->in_object = level+completedSites;
            }

            setComputingArrays (spanningTreeRoot, nil, computingOrder, storageOrder, cacheSlots, referenceOrder, parentOrder, nc);

            for (level=0; level<spanningTreeSites.lLength; level++) {
                ((node<long>*)spanningTreeSites.lData[level])->in_object = level+completedSites;
            }

            for (level=1; level<parentOrder.lLength; level++) {
                parentOrder.lData[level] = ((node<long>*)parentOrder.lData[level])->in_object;
            }

            if (completedSites) {
                level = mstCache->computingOrder.lLength-1;
                (*(_SimpleList*)mstCache->computingOrder(level)) << computingOrder;
                (*(_SimpleList*)mstCache->storageOrder(level)) << storageOrder;
                (*(_SimpleList*)mstCache->referenceOrder(level)) << referenceOrder;
                (*(_SimpleList*)mstCache->parentOrder(level)) << parentOrder;
                if (cacheSlots.lLength > mstCache->cacheSize.lData[level]) {
                    mstCache->cacheSize.lData[level] = cacheSlots.lLength;
                }
            } else {
                mstCache->computingOrder    && & computingOrder;
                mstCache->storageOrder      && & storageOrder;
                mstCache->referenceOrder    && & referenceOrder;
                mstCache->parentOrder       && & parentOrder;
                mstCache->cacheSize         << cacheSlots.lLength;
            }

            spanningTreeRoot->delete_tree();
            delete (spanningTreeRoot); // dmalloc fix 06162005

        }

        completedSites+=partition;
        if (vLevel>5) {
            char   buffer[64];
            snprintf (buffer, sizeof(buffer),"\n%ld %% done", (long)(completedSites*100/totalSites));
            BufferToConsole (buffer);
        }
    }

    _SimpleList straight (sl.lLength, 0, 1),
                * tcc = nil;

#ifdef _SLKP_LFENGINE_REWRITE_
    if (treeTraversalMasks.lLength > index) {
        tcc =  (_SimpleList*) treeTraversalMasks(index);
    }

#endif

    hyFloat strl = CostOfPath (df,t,straight),
               optl = CostOfPath (df,t,sl,tcc);

    if (vLevel>500) {
        _String* opPath = (_String*)sl.toStr();
        BufferToConsole("\nSite ordering:");
        StringToConsole(*opPath);
        DeleteObject(opPath);
    }

    char buffer [512];
    snprintf (buffer, sizeof(buffer),"\nPseudo-optimal path's cost %ld vs %ld for 1..k=> a %g x improvement",(long)optl,(long)strl,strl/(double)optl );
    ReportWarning (buffer);
    if (vLevel>5) {
        BufferToConsole (buffer);
    }

    if (mstCache) {
        long     memOverhead = mstCache->cacheSize.lData[mstCache->cacheSize.lLength-1];
        if (memOverhead) {
            memOverhead *= (t->GetINodeCount()*(sizeof(hyFloat)*t->GetCodeBase()+sizeof(long)+sizeof (char))+t->GetLeafCount()*(sizeof(hyFloat)+sizeof(long)))/1024;
            snprintf (buffer, sizeof(buffer),"\nIf using full MST heurisitcs: %ld vs %ld for 1..k=> a %g x (relative %g x) improvement with %ld KB memory overhead",globalLength,(long)strl,strl/(double)globalLength,optl/(double)globalLength,memOverhead);
            ReportWarning (buffer);
            if (vLevel>5) {
                BufferToConsole (buffer);
            }
        }
    }

}
//_______________________________________________________________________________________

void    _LikelihoodFunction::ComputePruningEfficiency (long& full, long& saved) {
    full = 0;
    saved = 0;
    for (long i=0; i<theTrees.lLength; i++) {
        _TheTree    *cT = ((_TheTree*)(LocateVar(theTrees(i))));
        _SimpleList *l  = (_SimpleList*)leafSkips (i);
        HBLObjectRef   lc  = cT->TipCount();

        long         leafCount = lc->Value(),
                     iCount;

        DeleteObject (lc);
        lc = cT->BranchCount ();
        iCount = lc->Value();
        DeleteObject (lc);

        saved += leafCount+iCount;
        full  += (leafCount+iCount) * (l->lLength+1);

        for (long k=0; k<l->lLength; k++) {
            unsigned long j = l->lData[k],
                          p1 = j&0xffff,
                          p2 = ((j>>16)&0xffff);

            saved += leafCount - p1 - (leafCount - 1 - p2);
            saved += iCount - cT->GetLeftINodes().lData[p1];
        }
    }
}

//_______________________________________________________________________________________

unsigned long    _LikelihoodFunction::CountObjects (_LikelihoodFunctionCountType kind) const {
    switch (kind) {
      case kLFCountGlobalVariables: {
          unsigned long res = 0UL;
          indexInd.Each ( [&res] (long var_id, unsigned long) -> void {
              if (LocateVar (var_id)->IsGlobal()) {
                  ++res;
              }
          });
          
          return res;
      }
      case kLFCountLocalCariables:
          return indexInd.countitems() - CountObjects (kLFCountGlobalVariables);
      case kLFCountDependentVariables:
          return indexDep.countitems();
      case kLFCountCategoryVariables:
          return indexCat.countitems();
    }

    return theTrees.countitems();
}

//_______________________________________________________________________________________

//______________________________________________________________________________
void _LikelihoodFunction::SerializeLF(_StringBuffer & rec, char opt,
                                      _SimpleList * partitionList,
                                      _SimpleList * exportPart) {

    // first - check how many data sets the LF depends on;
    // if only one - then spool it out into a NEXUS string;
    // if more than one - output a NEXUS file where
    // data sets after the first one are embedded as strings.
    RescanAllVariables();

    if (partitionList) {

        // validate the list
        _String errMsg;
        if (partitionList->lLength == 0 ||
            partitionList->lLength > theDataFilters.lLength) {
            errMsg = "The partition list for sub-export passed to SerializeLF "
                     "was either kEmptyString or too long";
        } else {
            // check for duplicates and index overrun
            partitionList->Sort();
            if (partitionList->lData[partitionList->lLength - 1] >=
                theDataFilters.lLength) {
                errMsg = "The partition list contained invalid partition "
                         "indices (too high) in call to SerializeLF";
            } else {
                for (unsigned long k2 = 1; k2 < partitionList->lLength; k2++) {
                    if (partitionList->lData[k2] ==
                        partitionList->lData[k2 - 1]) {
                        errMsg = "The partition list contained duplicate "
                                 "partition indices in call to SerializeLF";
                        break;
                    }
                }
            }
        }

        if (errMsg.nonempty()) {
            HandleApplicationError(errMsg);
            return;
        }
    }

    _String *lfName = (_String *)likeFuncNamesList(
        likeFuncList._SimpleList::Find((long) this));

    ReportWarning(_String("Serializing ") & *lfName & " with mode " &
                  _String((long) opt));

    if (opt == _hyphyLFSerializeModeShortMPI) {
        _String resVarName = *lfName & "_MLES";
        _Variable *resVar = FetchVar(LocateVarByName(resVarName));
        if (resVar) {
            rec.AppendAnAssignmentToBuffer(
                &resVarName, (_String *)resVar->Compute()->toStr());
        }

        resVarName = *lfName & "_MLE_VALUES";

        rec.AppendAnAssignmentToBuffer(&resVarName, new _String (kEmptyAssociativeList));
        rec.AppendVariableValueAVL(&resVarName, indexInd);
        rec.AppendVariableValueAVL(&resVarName, indexDep);

        return;
    }


    _SimpleList *redirector = nil;
    if (partitionList) {
        redirector = new _SimpleList;
        for (unsigned long pidx = 0; pidx < partitionList->lLength;
             pidx = pidx + 1) {
            (*redirector) << theDataFilters.lData[partitionList->lData[pidx]];
        }
    } else {
        redirector = &theDataFilters;
    }

    _SimpleList taggedDataSets, dataSetsByFilter;

    _AVLListX taggedDS(&taggedDataSets);

    for (unsigned long idx = 0; idx < redirector->lLength; idx++) {

        long tIdx = dataSetList._SimpleList::Find(
            (long)(GetDataFilter(redirector->get(idx))
                       ->GetData()));


        tIdx = taggedDS.Insert((BaseRef) tIdx, taggedDS.countitems());

        if (tIdx < 0L) {
            tIdx = -tIdx - 1L;
        }

        dataSetsByFilter << tIdx;

    }

    if (taggedDS.countitems() > 1 && exportPart) {
        HandleApplicationError(_String("Can't represent likelihood function ") &
                               ((_String *)likeFuncNamesList(
                                                             likeFuncList._SimpleList::Find((long) this)))->Enquote() &
                               " as a single file with a prespecified partition, because it "
                               "depends on multiple datasets. This is an internal error.");
        return;
    }

    _SimpleList indexedDataSets(taggedDS.countitems(), 0, 0);

    for (unsigned long idx2 = 0; idx2 < taggedDataSets.lLength; idx2++) {
        indexedDataSets.lData[taggedDS.xtraD.lData[idx2]] =
            taggedDataSets.lData[idx2];
    }

    _List involvedSites, avlSupport, dV;

    _SimpleList unique_sites;

    for (unsigned long idx3 = 0; idx3 < indexedDataSets.lLength; idx3++) {
        unique_sites << 0;
        _SimpleList *esl = new _SimpleList;
        avlSupport.AppendNewInstance(esl);
        involvedSites.AppendNewInstance(new _AVLListX(esl));
        dV.AppendNewInstance(new _SimpleList);
    }

    // create a dummy data filter and spool the data to a file
  
  
  
    if (partitionList || !exportPart) {
        for (unsigned long idx = 0; idx < redirector->lLength; idx++) {

            _SimpleList const *originalOrderFF =
                &GetDataFilter(redirector->get(idx))->theOriginalOrder;

            _AVLListX *involvedSitesL =
                (_AVLListX *)(involvedSites(dataSetsByFilter.lData[idx]));

            long *unique_sitesL =
                unique_sites.lData + dataSetsByFilter.lData[idx];

            for (unsigned long idx2 = 0; idx2 < originalOrderFF->lLength;
                 idx2++) {
                if (involvedSitesL->Insert(
                        (BaseRef)(originalOrderFF->lData[idx2]),
                        *unique_sitesL) >= 0) {
                    (*unique_sitesL)++;
                }
            }
        }

        for (unsigned long idx4 = 0; idx4 < indexedDataSets.lLength; idx4++) {

            _AVLListX *involvedSitesL = (_AVLListX *)involvedSites(idx4);
            _SimpleList tcache, *dVL = (_SimpleList *)dV(idx4);

            long iv, k = involvedSitesL->Traverser(tcache, iv,
                                                   involvedSitesL->GetRoot());

            for (; k >= 0; k = involvedSitesL->Traverser(tcache, iv)) {
                (*dVL) << (long) involvedSitesL->Retrieve(k);
                involvedSitesL->SetXtra(k, dVL->lLength - 1);
            }
        }
    } else if (exportPart) {
        ((_SimpleList *)dV(dataSetsByFilter(0)))->Duplicate(exportPart);
    }

    bool stashed_nfold = hy_env::EnvVariableTrue(hy_env::skip_omissions);
    hy_env::EnvVariableSet(hy_env::skip_omissions, new _Constant (HY_CONSTANT_FALSE), false);
    for (unsigned long idx5 = 0; idx5 < indexedDataSets.lLength; idx5++) {

        _DataSetFilter *entireSet = new _DataSetFilter();
        _SimpleList dH;
        entireSet->SetFilter(
            (_DataSet *)dataSetList(indexedDataSets.lData[idx5]), 1, dH,
            *(_SimpleList *)dV(idx5));
        if (idx5) {
            stashParameter(hy_env::data_file_print_format, 0.0, true);
        } else {
            stashParameter(hy_env::data_file_print_format, 4.0, true);
        }

        _String *cs = (_String *)entireSet->toStr();

        if (idx5) {
            if (idx5 == 1) {
                rec << "\n\nBEGIN HYPHY;\n\n";
            }
            rec << "_tdsstring_ = \"";
            rec.SanitizeAndAppend(*cs);
            rec << "\";\nDataSet ";
            rec << (_String *)dataSetNamesList(indexedDataSets.lData[idx5]);
            rec << " = ReadFromString (_tdsstring_);\n_tdsstring_=0;\n";
        } else {
            rec.AppendNewInstance(cs);
        }
        DeleteObject(entireSet);
        stashParameter(hy_env::data_file_print_format, 0.0, false);
    }
    hy_env::EnvVariableSet(hy_env::skip_omissions, stashed_nfold ? new  HY_CONSTANT_TRUE : new  HY_CONSTANT_FALSE , false);
  
    if (indexedDataSets.lLength == 1) {
        rec << "\n\nBEGIN HYPHY;\n\n";
    }

    if (opt == _hyphyLFSerializeModeOptimize) {
         if (! hy_env::EnvVariableTrue(kUseLastResults) && hy_env::EnvVariableTrue(kUseInitialDistanceGuess)) {
            for (unsigned long i = 0; i < indexInd.lLength; i++) {
                LocateVar(indexInd.lData[i])->MarkDone();
            }

            GetInitialValues();
        }

        _FString *mpiPrefix = (_FString *)FetchObjectFromVariableByType(
            &mpiPrefixCommand, STRING);

        if (mpiPrefix) {
            rec << mpiPrefix->get_str();
        }
    }

    // write out all globals
    _StringBuffer glVars(1024L), locVars(1024L);

    char str[4096];

    _SimpleList *indepVarList = nil, *depVarList = nil, *catVarList = nil;

    if (partitionList) {
        indepVarList = new _SimpleList;
        depVarList = new _SimpleList;
        catVarList = new _SimpleList;

        ScanAllVariablesOnPartition(*partitionList, *indepVarList, *depVarList,
                                    *catVarList);
    } else {
        indepVarList = &indexInd;
        depVarList = &indexDep;
        catVarList = &indexCat;
    }

    ExportIndVariables(glVars, locVars, indepVarList);
    ExportDepVariables(glVars, locVars, depVarList);

    glVars.TrimSpace();
    locVars.TrimSpace();

    rec << glVars;

    // write out all categs
    if (opt == _hyphyLFSerializeModeCategoryAsGlobal) {
        for (long idx = 0; idx < catVarList->lLength; idx++) {
            _CategoryVariable *theC =
                (_CategoryVariable *)LocateVar(catVarList->lData[idx]);
            snprintf(str, sizeof(str), "\nglobal %s;",
                     theC->GetName()->get_str());
            rec << str;
        }
    } else {
        ExportCatVariables(rec, catVarList);
    }

    //write out all the models
    _SimpleList *redirectorT = nil, dH, dV2, dU;

    if (partitionList) {
        redirectorT = new _SimpleList;
        for (long pidx = 0; pidx < partitionList->lLength; pidx = pidx + 1) {
            (*redirectorT) << theTrees.lData[partitionList->lData[pidx]];
        }
    } else {
        redirectorT = &theTrees;
    }

    for (long idx = 0; idx < redirectorT->lLength; idx++) {
        _SimpleList dT;
        ((_TheTree *)LocateVar(redirectorT->lData[idx]))
            ->CompileListOfModels(dT);

        if (dT.lLength == 1) {
            dV2 << dT.lData[0];
        } else {
            dV2 << -1;
        }

        dH.Union(dU, dT);
        dU.Duplicate(&dH);
    }

    {
        _SimpleList modelAux;
        _AVLList modelDupChecker(&modelAux);
        for (long idx = 0; idx < dH.lLength; idx++) {
            SerializeModel(rec, dH.lData[idx], &modelDupChecker);
        }

    }


    // write out all the trees, including model definitions if needed
    //hyFloat stashIM = hy_env::EnvVariableTrue(hy_env::accept_rooted_trees);
    //checkParameter(tryNumericSequenceMatch, stashIM, 0.0);
  

    rec.AppendAnAssignmentToBuffer(&hy_env::try_numeric_sequence_match,
                                   new _String((hyFloat)hy_env::EnvVariableTrue(hy_env::try_numeric_sequence_match)));
    rec.AppendAnAssignmentToBuffer(&hy_env::accept_rooted_trees,
                                 new _String((hyFloat)hy_env::EnvVariableTrue(hy_env::accept_rooted_trees)));

    bool include_model_info = hy_env::EnvVariableTrue(hy_env::include_model_spec);
 
    for (long idx = 0; idx < redirectorT->lLength; idx++) {
        if (dV2.lData[idx] >= 0) {
            rec << "\nUseModel (";
            rec << *((_String *)modelNames(dV2.lData[idx]));
            rec << ");\n";
            hy_env::EnvVariableSet(hy_env::include_model_spec, new HY_CONSTANT_FALSE, false);
        } else {
            hy_env::EnvVariableSet(hy_env::include_model_spec, new HY_CONSTANT_TRUE, false);
        }

        rec << "Tree ";
        rec << LocateVar(redirectorT->lData[idx])->GetName();
        rec << '=';
        rec.AppendNewInstance((_String *)(
            (_TheTree *)LocateVar(redirectorT->lData[idx]))->toStr());
        rec << '\n';
    }

    hy_env::EnvVariableSet(hy_env::include_model_spec, include_model_info ? new HY_CONSTANT_TRUE : new HY_CONSTANT_FALSE, false);
    

    rec << locVars;

    // spool out the specs for datafilter
    rec << "\nDataSet ";
    rec << *((_String *)dataSetNamesList(indexedDataSets.lData[0]));
    rec << " = ReadDataFile(";
    rec << useNexusFileData;
    rec << ");\n";

    dU.Clear();
    _AVLList writtenDF(&dU);

    for (long idx = 0; idx < redirector->lLength; idx++) {
        if (writtenDF.Insert((BaseRef) redirector->lData[idx]) >= 0) {
            _DataSetFilter const *theDF = GetDataFilter (redirector->get(idx));
            _String *horPart = nil;

            if (partitionList || !exportPart) {
                _SimpleList remappedOO;
                _AVLListX *involvedSitesL =
                    (_AVLListX *)involvedSites(dataSetsByFilter.lData[idx]);
                for (long idx2 = 0; idx2 < theDF->theOriginalOrder.lLength;
                     idx2++) {
                    remappedOO
                        << involvedSitesL->GetXtra(involvedSitesL->Find(
                               (BaseRef) theDF->theOriginalOrder.lData[idx2]));
                }
                horPart = (_String *)remappedOO.ListToPartitionString();
            } else if (exportPart) {
                horPart =
                    new _String(_String("0-") & (long) exportPart->lLength - 1);
            }
            // else
            // horPart = (_String*)theDF->theOriginalOrder.ListToPartitionString();

            rec     << "DataSetFilter "
                    << *GetFilterName(redirector->get(idx))
                    << " = CreateFilter("
                    << *((_String *)dataSetNamesList(
                       indexedDataSets.lData[dataSetsByFilter.lData[idx]]))
                    << ','
                    << _String((long) theDF->GetUnitLength())
                    << ',';

            if (horPart) {
              rec << horPart->Enquote('"');
              DeleteObject(horPart);
            }
            horPart = (_String *)theDF->theNodeMap.ListToPartitionString();

            rec << ','
                << horPart->Enquote('"');

            DeleteObject(horPart);

            horPart = theDF->GetExclusions();
            if (horPart->nonempty()) {
                rec << ',' << horPart->Enquote('"');
            }
            DeleteObject(horPart);

            rec << ");\n";
        }
    }

    // write out the global variable for enforcing reversible models
    rec.AppendAnAssignmentToBuffer(&hy_env::assume_reversible, new _String(hy_env::EnvVariableGetNumber(hy_env::assume_reversible, 0.)));

    rec << "LikelihoodFunction " << *lfName << " = (";

    long dsID = 0;
    {
        for (long idx = 0; idx < redirector->lLength; idx++) {
            if (dsID) {
                rec << ',';
            } else {
                dsID = 1;
            }

            rec << *GetFilterName(redirector->get(idx))
             << ','
             << *LocateVar(redirectorT->lData[idx])->GetName();
        }
    }
    if (computingTemplate && templateKind == 1) {
        rec << ",\"";
        rec << (_String *)computingTemplate->toStr(kFormulaStringConversionNormal);
        rec << '"';
    }

    if (opt == _hyphyLFSerializeModeOptimize) {
        rec << ");\n";
        rec.AppendAnAssignmentToBuffer(&kOptimizationPrecision, new _String(hy_env::EnvVariableGetNumber(kOptimizationPrecision, 0.01)));
        rec.AppendAnAssignmentToBuffer(&kOptimizationMethod, new _String(hy_env::EnvVariableGetNumber(kOptimizationMethod, 4.)));
        rec.AppendAnAssignmentToBuffer(&kUseInitialDistanceGuess, new _String(hy_env::EnvVariableGetNumber(kUseInitialDistanceGuess, 1.)));
        rec.AppendAnAssignmentToBuffer(&kUseLastResults, new _String(hy_env::EnvVariableGetNumber(kUseLastResults, 1.)));
        
        hyFloat hard_limit = hy_env::EnvVariableGetNumber(kOptimizationHardLimit, -1.);
        if (hard_limit > 0.) {
            rec.AppendAnAssignmentToBuffer(&kOptimizationHardLimit,
                                           new _String(hard_limit));
        }
        rec << "Optimize(";
        rec << lfName;
        rec << "_MLES,";
        rec << lfName;
        rec << ");\n";
        rec << mpiMLELFValue;
        rec << '=';
        rec << lfName;
        rec << "_MLES[1][0];\n";
        rec << lf2SendBack;
        rec << "=\"";
        rec << lfName;
        rec << "\";\n";
        rec.AppendAnAssignmentToBuffer(&shortMPIReturn, new _String(hy_env::EnvVariableGetNumber(shortMPIReturn, 0.)));
    } else if (opt == _hyphyLFSerializeModeLongMPI) {
        rec << ");\n";
        rec.AppendAnAssignmentToBuffer(
            &mpiMLELFValue, new _String(FetchVar(LocateVarByName(mpiMLELFValue))
                                            ->Compute()->Value()));

        _String resVarName = *lfName & "_MLES";
        _Variable *resVar = FetchVar(LocateVarByName(resVarName));
        if (resVar) {
            rec.AppendAnAssignmentToBuffer(
                &resVarName, (_String *)resVar->Compute()->toStr());
        }
    } else {
        rec << ");";
    }

    _FString *haveExtra =
        (_FString *)FetchObjectFromVariableByType(&lfExtraLFExportCode, STRING);
    if (haveExtra) {
        rec << haveExtra->get_str();
    }

    rec << "\n\nEND;";

    if (partitionList) {
        BatchDelete(redirector, indepVarList,depVarList,catVarList,redirectorT);
    }
}

//_______________________________________________________________________________________

BaseRef _LikelihoodFunction::toStr (unsigned long) {
    hyFloat longOrShort,
               value = 0.0;

    checkParameter(likefuncOutput,longOrShort,2.0);

    if (longOrShort < 4.0) {
        PrepareToCompute();
        value = Compute();
    }


    if (longOrShort>5.5) { // serialized self-contained LF
        _StringBuffer   *sLF = new _StringBuffer (8192L);
        SerializeLF     (*sLF);
        sLF->TrimSpace  ();
        return          sLF;
    }

    _StringBuffer * res = new _StringBuffer (1024UL);
    
    auto _spool_bounds = [] (_StringBuffer& res, _Variable* const the_var) -> void {
        if (!CheckEqual(the_var->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
            res << *the_var->GetName() << ":>" << _String (the_var->GetLowerBound(), "%.16g") << ';';
        }
        
        if (!CheckEqual(the_var->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
            res << *the_var->GetName() << ":<" << _String (the_var->GetUpperBound(), "%.16g") << ';';
        }
    };

    if (longOrShort>=4.0) { // just spool out the names of parameters
        _Variable * thisVar;

        for (unsigned long i=0UL; i<indexInd.countitems(); i++) {
            thisVar = GetIthIndependentVar(i);
            
            *res << '\n';
            if (thisVar->IsGlobal()) {
                *res << "global ";
            }
            *res << *thisVar->GetName() << '=' << _String (GetIthIndependent(i), "%.16g") << ';';
            _spool_bounds (*res, thisVar);
        }
        for (unsigned long i=0UL; i<indexDep.countitems(); i++) {
            thisVar = GetIthDependentVar(i);

            res->AppendAnAssignmentToBuffer(thisVar->GetName(), thisVar->GetFormulaString(kFormulaStringConversionNormal), kAppendAnAssignmentToBufferFree | (thisVar->IsGlobal() ? kAppendAnAssignmentToBufferGlobal : 0));
            _spool_bounds (*res, thisVar);
        }
        
        if (longOrShort>4.0) {
            longOrShort = 2.5;
        } else {
            DoneComputing();
            res->TrimSpace ();
            return res;
        }

    }
    
    if (longOrShort==3.0) { // just spool out the names of parameters
        *res << "Log Likelihood = " << _String (value, "%g") << ";\n\nIndependent Parameters List\n";
        
        for (unsigned long i=0UL; i<indexInd.countitems(); i++) {
            *res <<  "\n  Parameter " << _String ((long)(i+1)) << " : " << *GetIthIndependentName(i);
        }
        if (indexDep.countitems()) {
            *res << "\n\nConstrained Parameters List\n";
            for (unsigned long i=0UL; i<indexDep.countitems(); i++) {
                *res <<  "\n  Parameter " << _String ((long)(i+1)) << " : " << *GetIthDependentName(i) << " = " <<
                    GetIthDependentVar(i)->GetFormulaString(kFormulaStringConversionNormal);
            }
        }
    } else if (longOrShort>1.1) {
        if (longOrShort!=2.5) {
            * res << "Log Likelihood = " << _String (value, "%15.15g") << ";";
            
            bool globals = false;
            
            for (unsigned long i = 0UL; i<indexInd.countitems(); i++) {
                _Variable* v = GetIthIndependentVar(i);
                if (v->IsGlobal()) {
                    if (!globals) {
                        globals = true;
                        *res<< "\nGlobal Parameters:\n";
                    }
                    res->AppendAnAssignmentToBuffer(GetIthIndependentName(i), (_String*)v->toStr());
                }
            }

            for (unsigned long i = 0UL; i<indexDep.countitems(); i++) {
                _Variable* v = GetIthDependentVar(i);
                if (v->IsGlobal()) {
                    if (!globals) {
                        globals = true;
                        *res<< "\nGlobal Parameters:\n";
                    }
                    _String value ((_String*)v->GetFormulaString(kFormulaStringConversionNormal));
                    value = value & " = " & _String ((_String*)v->toStr());
                    res->AppendAnAssignmentToBuffer(GetIthDependentName(i), &value, kAppendAnAssignmentToBufferPlain);
                }
            }
    
        }
        // traverse the trees now
        for (unsigned long treeCounter = 0UL; treeCounter<theTrees.countitems(); treeCounter++) {
            _TheTree* currentTree = GetIthTree (treeCounter);
            
            long level = 0, lastLevel = 0,l1,l2,j;
            *res<<"\nTree "
                <<currentTree->GetName()
                <<'=';

            _TreeIterator ti (currentTree, _HY_TREE_TRAVERSAL_POSTORDER);

            _CalcNode* currentNode=ti.Next(),
                     * nextNode;

            level    = ti.Depth();
            nextNode = ti.Next();

            // decide if we can use expected substituion measure
            bool    useExpectedSubstitutions = longOrShort<3;
 
            while (nextNode) {
                _SimpleList iV, dV2;
                
                if (level>lastLevel) {
                    if (lastLevel) {
                        *res<<',';
                    }
                    res->AppendNCopies('(', level-lastLevel);

                } else if (level<lastLevel) {
                    res->AppendNCopies(')', lastLevel-level);
                } else {
                    *res<<',';
                }
                *res<< currentTree->GetNodeName (currentNode);
                if (!useExpectedSubstitutions) {
                    *res<<'{';
                    
                    _AVLList    iVA (&iV),
                                dVA (&dV2);

                    currentNode->ScanContainerForVariables(iVA,dVA);
                    currentNode->ScanForDVariables(dVA,iVA);

                    iVA.ReorderList ();
                    dVA.ReorderList ();
                    
                    iV.Each([&] (long var_index, unsigned long i) -> void {
                        _Variable * v = LocateVar (var_index);
                        if (i) {
                            *res << ',';
                        }
                        *res << v->ContextFreeName() << " = ";
                        res->AppendNewInstance((_String*)v->toStr());
                    });

                    dV2.Each([&] (long var_index, unsigned long i) -> void {
                        _Variable * v = LocateVar (var_index);
                        if (i) {
                            *res << ',';
                        }
                        *res << v->ContextFreeName() << " = ";
                        (res->AppendNewInstance(v->GetFormulaString(kFormulaStringConversionNormal)) << " = ")
                           .AppendNewInstance((_String*)v->toStr());
                    });
                    
                    *res<<'}';
                } else {
                    *res<<':';
                    if (!currentNode->ComputeModelMatrix()) {
                        *res <<  "NAN";
                    } else {
                        *res << _String ((hyFloat)fabs(currentNode->ComputeBranchLength()));
                    }

                }
                lastLevel = level;
                level = ti.Depth();
                currentNode = nextNode;
                nextNode=ti.Next() ;

            }
            res->AppendNCopies(')', lastLevel - level);

            *res<<';';
        }
        *res << '\n';
    } else if (longOrShort>.1) {
        * res << "Log Likelihood = " << _String (value, "%15.15g") << ";\n";

        for (unsigned long i = 0UL; i<indexInd.countitems(); i++) {
            res->AppendAnAssignmentToBuffer(GetIthIndependentName(i), (_String*)GetIthIndependentVar(i)->toStr());
        }
        
        for (unsigned long i = 0UL; i<indexDep.countitems(); i++) {
            _Variable* v = GetIthDependentVar(i);
            _String value ((_String*)v->GetFormulaString(kFormulaStringConversionNormal), kAppendAnAssignmentToBufferPlain);
            value = value & " = " & _String ((_String*)v->toStr());
            res->AppendAnAssignmentToBuffer(GetIthDependentName(i), &value);
        }
    } else {
        * res << _String (value, "%15.15g");
    }

    if (longOrShort < 4.0) {
        DoneComputing();
    }

    res->TrimSpace();
    return res;

}

//_______________________________________________________________________________________

void    _LikelihoodFunction::StateCounter (long functionCallback) const {
  HandleApplicationError ("This feature has not yet been implemented in the new LF engine framework");
}

    //_______________________________________________________________________________________

  void    _LikelihoodFunction::Simulate (_DataSet &target, _List& theExclusions, _Matrix* catValues, _Matrix* catNames, _Matrix* spawnValues, _String const* storeIntermediates) const {
      // will step thru multiple trees of the project and simulate  a dataset from the likelihood function

    enum {
      kLFSimulateCategoriesNone,
      kLFSimulateCategoriesDiscrete,
      kLFSimulateCategoriesContinuous
    } category_simulation_mode = kLFSimulateCategoriesNone;

      //char catSim = 0;
      // 0 - not at all
      // 1 - discrete,
      // 2 - continuous


    _SimpleList     continuous_category_variables,
    discrete_category_variables,
    HMM_category_variables,
    HMM_state;

    bool            simulate_column_wise = false;

    if (indexCat.lLength) {
      checkParameter (categorySimulationMethod,categorySimMethod,2.0);

      if (categorySimMethod > 1.5) {
        category_simulation_mode = kLFSimulateCategoriesContinuous;
      } else {
        if (categorySimMethod > 0.5) {
          category_simulation_mode = kLFSimulateCategoriesDiscrete;
        }
      }


      if (category_simulation_mode  == kLFSimulateCategoriesContinuous) {
        for (unsigned long cat_index = 0UL; cat_index < indexCat.lLength; cat_index ++) {
          _CategoryVariable* ith_category = GetIthCategoryVar(cat_index);
          if (ith_category->IsHiddenMarkov()) {
            HMM_category_variables << cat_index;
            HMM_state << 0L;
          } else if (ith_category->GetCumulative().IsEmpty()) {
            ReportWarning(_String (__PRETTY_FUNCTION__) & " will treat category variable " & ith_category->GetName()->Enquote() &
                          " as discrete since no cumulative distribution is available.");
            discrete_category_variables << cat_index;
            HMM_state << cat_index;
          } else {
            continuous_category_variables << cat_index;
          }
        }
      } else {
        for (unsigned long cat_index = 0UL; cat_index < indexCat.lLength; cat_index ++) {
          GetIthCategoryVar(cat_index)->Refresh (true);
          discrete_category_variables << cat_index;
        }
      }

      if (catNames && indexCat.lLength) {

        catNames->Clear();
        _Matrix::CreateMatrix (catNames,indexCat.lLength,1,false,true,false);
        catNames->Convert2Formulas();

        _SimpleList* all_arrays[3] = { & HMM_category_variables, & discrete_category_variables, & continuous_category_variables };

        unsigned long through_index = 0UL;

        for (_SimpleList * array : all_arrays) {
          for (unsigned long i = 0UL; i < array->lLength; i++, through_index++) {
            catNames->StoreFormula (through_index,0L,*new _Formula (new _FString (*LocateVar (array->GetElement(i))->GetName())), false, false);
          }
        }

      }
    }

    _DataSetFilter const *first_filter = GetIthFilter(0);

    unsigned long species_count = first_filter->NumberSpecies();

    for (unsigned long sequence_index = 0UL; sequence_index < species_count; sequence_index ++) {
      target.AddName(*first_filter->GetSequenceName(sequence_index));
    }

    unsigned long    internal_node_count = 0UL;

    if (storeIntermediates && storeIntermediates->nonempty() == 0) {
      GetIthTree (0L)->AddNodeNamesToDS (&target,false,true,0); // only add internal node names
      internal_node_count = target.GetNames().lLength - species_count;
    }
    target.SetTranslationTable (first_filter->GetData());

    unsigned long sequences_to_simulate = target.GetNames().lLength,
    site_offset_raw       = 0UL,
      // raw offset (e.g. in nucleotides for codon data
    site_offset          = 0UL,
    total_sites          = 0UL;

    for (unsigned long i = 0UL; i<theTrees.lLength; i++) {
      total_sites += GetIthFilter(i)->GetSiteCountInUnits();
    }

    if (catValues && indexCat.lLength) {
      catValues->Clear();
      _Matrix::CreateMatrix (catValues,indexCat.lLength,total_sites,false,true,false);
    }

    TimeDifference timer;


    bool       column_wise  = false;


    for (unsigned long partition_index = 0UL; partition_index<theTrees.lLength; partition_index++) { // loop thru the trees one at a time
                                                                                                     // will propagate a complete column thru the tree
                                                                                                     // first generate a root sequence
      _SimpleList user_exclusions_numeric;

      _DataSetFilter const * this_filter  = GetIthFilter(partition_index);
      if (this_filter->NumberSpecies()+internal_node_count != sequences_to_simulate) {
        ReportWarning (_String ("Ignoring partition ") & _String ((long) (partition_index + 1L)) & " of the likelihood function since it has a different number of sequences/tree than the first part.");
        continue;
      }
      hyFloat * this_freqs = ((_Matrix*)GetIthFrequencies(partition_index)->ComputeNumeric())->fastIndex();
      _TheTree *this_tree = GetIthTree (partition_index);


      unsigned long    this_site_count  = this_filter->GetSiteCountInUnits(),
      leaf_count       = 0L,
      good_sites       = 0L,
      filter_dimension = this_filter->GetDimension(true),
      sites_per_unit    = this_filter->GetUnitLength(),
      this_raw_site_count = this_filter->GetSiteCount();


      if (theExclusions.lLength>partition_index) {
        _List*      user_exclusions = (_List*)theExclusions(partition_index);


        if (user_exclusions->countitems() > 0) {
          const _TranslationTable* this_translation_table = this_filter->GetTranslationTable();

          for (unsigned long state = 0UL; state < user_exclusions->lLength; state++) {
            long resolved_state = this_translation_table->MultiTokenResolutions(*(_String*)user_exclusions->GetItem(state), NULL);
            if (resolved_state < 0L) {
              ReportWarning (_String ("Excluded character ") & ((_String*)user_exclusions->GetItem(state))->Enquote() & " does not represent a unique state and will therefore be ignored in data simulation");
            } else {
              user_exclusions_numeric << state;
            }
          }

          column_wise = true;
        }
      }

      if (category_simulation_mode != kLFSimulateCategoriesNone) {
        column_wise = true;
      }

      if (column_wise) {

        _TheTree * this_tree = GetIthTree (partition_index);
        this_tree->SetUpMatrices(1);
        while (good_sites < this_site_count) {

          if (category_simulation_mode != kLFSimulateCategoriesNone ) {

            for (unsigned long hmm_category_index =0UL;
                 hmm_category_index < HMM_category_variables.lLength;
                 hmm_category_index ++) { // deal with HMM

              _CategoryVariable* hmm_cat = GetIthCategoryVar(HMM_category_variables(hmm_category_index));

              _Matrix* category_weight_matrix = hmm_cat->GetWeights();
              hyFloat* category_weights;

              unsigned long category_count = hmm_cat->GetNumberOfIntervals();

              if (good_sites == 0L) {
                category_weights = category_weight_matrix->fastIndex();
              } else {
                _Matrix * hmm = hmm_cat->ComputeHiddenMarkov();
                category_weights = hmm->theData+hmm->GetVDim()* HMM_state(hmm_category_index);
              }

              unsigned long root_state = DrawFromDiscrete(category_weights, category_count);

              hmm_cat->SetIntervalValue(root_state);
              if (good_sites > 0L) {
                HMM_state [hmm_category_index] = root_state;
              }

              if (catValues) {
                catValues->Store (hmm_category_index,site_offset+good_sites,hmm_cat->Compute()->Value());
              }
            }

            for (unsigned long discrete_category_index = 0UL;
                 discrete_category_index < discrete_category_variables.lLength;
                 discrete_category_index++) {


              _CategoryVariable* discrete_cat = GetIthCategoryVar(discrete_category_variables(discrete_category_index));

              unsigned long category_value = DrawFromDiscrete(discrete_cat->GetWeights()->fastIndex(), discrete_cat->GetNumberOfIntervals());

              discrete_cat->SetIntervalValue(category_value);

              if (catValues) {
                catValues->Store (discrete_category_index+HMM_category_variables.lLength,site_offset+good_sites,discrete_cat->Compute()->Value());
              }
            }

            for (unsigned long continuous_category_index = 0UL;
                 continuous_category_index < continuous_category_variables.lLength;
                 continuous_category_index++) { // use discrete values here

              _CategoryVariable* continuous_cat_var = GetIthCategoryVar(continuous_category_variables(continuous_category_index));

              hyFloat  category_value = continuous_cat_var->GetCumulative().Newton(continuous_cat_var->GetDensity(),MAX (genrand_real2(), 1e-30),continuous_cat_var->GetMinX(),continuous_cat_var->GetMaxX(),hy_x_variable);

              continuous_cat_var->SetValue(new _Constant (category_value), false);
              if (catValues) {
                catValues->Store (continuous_category_index+discrete_category_variables.lLength+HMM_category_variables.lLength,site_offset+good_sites,category_value);
              }
            }
          } // end category initialization block

          unsigned long root_state;

          if (spawnValues) {
            root_state = spawnValues->theData[site_offset+good_sites];
          } else {
            root_state = DrawFromDiscrete(this_freqs, filter_dimension);
          }


          _SimpleList ancestral_values,
                      leaf_values;

          if (SingleBuildLeafProbs    (this_tree->GetRoot(), root_state, leaf_values, user_exclusions_numeric,this_tree, true, this_filter, storeIntermediates?&ancestral_values:nil)) {
            good_sites++;
              //add this site to the simulated dataset



            _String simulated_unit (this_filter->ConvertCodeToLetters(this_filter->CorrectCode(leaf_values(0)), sites_per_unit));

            for (unsigned long character_index = 0UL; character_index < sites_per_unit; character_index ++) {
              target.AddSite (simulated_unit (character_index));
              leaf_count ++;
            }

            for (unsigned long sequence_index = 1UL; sequence_index < species_count; sequence_index++) {
              simulated_unit = this_filter->ConvertCodeToLetters(this_filter->CorrectCode(leaf_values(sequence_index)), sites_per_unit);
              for (unsigned long character_index = 0UL; character_index < sites_per_unit; character_index ++) {
                target.Write2Site (site_offset_raw + leaf_count - sites_per_unit + character_index, simulated_unit (character_index));
              }
            }


            target.ResetIHelper();
            for (unsigned long character_index = 0UL; character_index < sites_per_unit; character_index ++) {
              target.Compact(site_offset_raw + leaf_count - sites_per_unit + character_index);
            }

            if (storeIntermediates && storeIntermediates->nonempty() == 0UL) {
              for (unsigned long internal_node_index = 0UL; internal_node_index < internal_node_count; internal_node_index++) {
                simulated_unit = this_filter->ConvertCodeToLetters(this_filter->CorrectCode(ancestral_values(internal_node_index)), sites_per_unit);
                for (unsigned long character_index = 0UL; character_index < sites_per_unit; character_index ++) {
                  target.Write2Site (site_offset_raw + leaf_count - sites_per_unit + character_index, simulated_unit (character_index));
                }
                target.ResetIHelper();
                for (unsigned long character_index = 0UL; character_index < sites_per_unit; character_index ++) {
                  target.Compact(site_offset_raw + leaf_count - sites_per_unit + character_index);
                }

              }
            }
          }

          hyFloat time_elapsed = timer.TimeSinceStart();


          if (time_elapsed > .25) {
#if !defined __UNIX__ || defined __HEADLESS__
            SetStatusBarValue (100.*(site_offset_raw+good_sites)/total_sites, 1, 0);
#endif
            timer.Start();
          }
        }
        this_tree->CleanUpMatrices();

      } else {// end simulate column by column


        unsigned long *  simulated_sequence = new unsigned long [this_site_count];

          // generate a random "spawning vector"

        if (spawnValues) { // use supplied starting values
          for (unsigned long site_index = 0UL;  site_index < this_site_count; site_index++) {
            simulated_sequence[site_index] = spawnValues->theData[site_index+site_offset];
          }
        } else {
          for (unsigned long site_index = 0UL;  site_index < this_site_count; site_index++) {
            simulated_sequence[site_index] = DrawFromDiscrete(this_freqs, filter_dimension);
          }
        }

          // now proceed down the tree branches to get the values of the probabilities at the leaves
          // this is done recursively

        _DataSet * ancestral_sequences = nil;

        if (storeIntermediates) {
          if (storeIntermediates->nonempty()) {
            FILE * file_for_ancestral_sequences = doFileOpen (storeIntermediates->get_str(),"w");
            if (!file_for_ancestral_sequences) {
              HandleApplicationError (_String ("Failed to open ") & storeIntermediates->Enquote() & " for writing.");
              target.Finalize();
              return;
            } else {
              ancestral_sequences = new _DataSet (file_for_ancestral_sequences);
              _TheTree *datree = (_TheTree*)LocateVar(theTrees(0));
              datree->AddNodeNamesToDS (ancestral_sequences,false,true,0);
            }
          } else {
            ancestral_sequences = new _DataSet (this_site_count);
          }

          ancestral_sequences->SetTranslationTable (this_filter->GetData());
        }

        BuildLeafProbs (this_tree->GetRoot(), simulated_sequence, this_site_count, target, this_tree, leaf_count, true, sites_per_unit, this_filter, site_offset_raw ,ancestral_sequences);

        if (ancestral_sequences) {
          ancestral_sequences->Finalize();
          if (storeIntermediates->nonempty() == 0) {
            for (unsigned long sequence_index = 0UL; sequence_index < internal_node_count; sequence_index ++ ) {
              for (unsigned long raw_site_index = 0UL; raw_site_index < this_raw_site_count; raw_site_index++) {
                target.Write2Site(site_offset_raw + raw_site_index, ancestral_sequences->GetSite (raw_site_index)->char_at(sequence_index));
              }
            }
          }
          DeleteObject (ancestral_sequences);
        }

        delete [] simulated_sequence;
      } // end over sequence-wise simulation

      site_offset_raw   += this_raw_site_count;
      site_offset += this_site_count;

    } // end loop over partitions

    target.Finalize();
    target.SetNoSpecies(target.GetNames().lLength);
  }

//_______________________________________________________________________________________

void    _LikelihoodFunction::BuildLeafProbs (node<long>& curNode, long unsigned * baseVector, unsigned long vecSize, _DataSet& target, _TheTree* curTree, unsigned long& leafCount, bool isRoot, long baseLength, _DataSetFilter const* dsf, long DSOffset, _DataSet* intNodes) const
/* SLKP TODO check that this works with the new category spec route */

{
  unsigned long * curVector = nil;

    _CalcNode* ccurNode = (_CalcNode*)LocateVar (curNode.get_data());

    if (!isRoot) {

        curVector = new unsigned long [vecSize];
        // first "mutate" the parent vector

        if (ccurNode->NeedNewCategoryExponential(-1)) {
            ccurNode->RecomputeMatrix(0,1);
        }

        hyFloat * baseI = ccurNode->GetCompExp()->fastIndex();

        unsigned long const m = ccurNode->GetCompExp()->GetVDim();

        for (unsigned long i = 0UL; i<vecSize; i++) {
          curVector[i] = DrawFromDiscrete( baseI + baseVector[i]*m, m);
        }
    } else {
        // handle the degenerate tree case
        if (curNode.nodes.length == 1) {
            for (unsigned long k = 0UL; k<vecSize; k++) {
                _String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(baseVector[k]), baseLength);
                for (unsigned long m = 0UL; m<letterValue.length(); m++) {
                    target.AddSite (letterValue.char_at(m));
                }
            }
            leafCount++;
            BuildLeafProbs (*curNode.go_down(1), baseVector, vecSize, target, curTree, leafCount, false, baseLength, dsf,DSOffset,intNodes);
            return;
        }
    }

    // now scan the "children" and pass on the parameters as needed

    if (curNode.nodes.length) {
        if (intNodes) {
            bool writeOrAdd = intNodes->lLength;
            for (unsigned long k = 0UL; k<vecSize; k++) {
                _String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(curVector?curVector[k]:baseVector[k]), baseLength);
                if (writeOrAdd)
                    for (unsigned long m = 0UL; m<letterValue.length(); m++) {
                        intNodes->Write2Site (letterValue.length()*k+m, letterValue.char_at(m));
                    }
                else
                  for (unsigned long m = 0UL; m<letterValue.length(); m++) {
                      intNodes->AddSite (letterValue.char_at(m));
                  }
            }
        }
        for (unsigned long k = 1UL; k<=curNode.get_num_nodes(); k++) {
            BuildLeafProbs (*curNode.go_down(k), curVector?curVector:baseVector, vecSize, target, curTree, leafCount, false, baseLength, dsf,DSOffset,intNodes);
        }
    } else
        // reached a leaf
    {
        // attach a row to the new data set
        long siteCount = DSOffset;
        if (!leafCount) {
            for (unsigned long k = 0UL; k<vecSize; k++) {
                _String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(curVector[k]), baseLength);
                for (unsigned long m = 0UL; m<letterValue.length(); m++) {
                    target.AddSite (letterValue.char_at(m));
                }
            }
            leafCount++;
        } else {
            for (unsigned long k = 0UL; k<vecSize; k++) {
                _String letterValue = dsf->ConvertCodeToLetters (dsf->CorrectCode(curVector[k]), baseLength);
                for (unsigned long m = 0UL; m<letterValue.length(); m++) {
                    target.Write2Site (siteCount++, letterValue.char_at(m));
                }
            }
        }
    }

    if (!isRoot) {
        ccurNode->FreeUpMemory (0);
    }

    if (curVector) {
        delete [] curVector;
    }
}

//_______________________________________________________________________________________

bool    _LikelihoodFunction::SingleBuildLeafProbs (node<long>& curNode, long parentState, _SimpleList& target, _SimpleList& theExc, _TheTree* curTree, bool isRoot, _DataSetFilter const* dsf, _SimpleList * iNodes) const {

    long myState = parentState;

    if (!isRoot) {

      _CalcNode* ccurNode = (_CalcNode*)LocateVar (curNode.get_data());

      if (ccurNode->NeedNewCategoryExponential(-1)) {
        ccurNode->RecomputeMatrix(0,1);
      }

      unsigned long matrix_dimension = ccurNode->GetCompExp()->GetVDim();


      hyFloat * fastI = ccurNode->GetCompExp()->fastIndex()+parentState*matrix_dimension;
      myState = DrawFromDiscrete(fastI, matrix_dimension);

      if (! curNode.is_leaf()) {
        if (iNodes) {
          if (theExc.Find(myState)!=kNotFound) {
            return false;
          }
          (*iNodes)<<myState;
          //return true;
        }
      } else {
        if (theExc.Find(myState)!=-kNotFound) {
          return false;
        }
        target<<myState;
        return true;
      }
    } else {
      if (curNode.nodes.length == 1) { // two taxon sumulation
        target << parentState;
      } else if (iNodes) {
        (*iNodes)<<parentState;
      }
    }

    // now scan the "children" and pass on the parameters as needed

    for (long k = 1; k<=curNode.get_num_nodes(); k++) {
      if(!SingleBuildLeafProbs (*curNode.go_down(k), myState, target, theExc, curTree, false, dsf, iNodes)) {
        return false;
      }
    }
    return true;
  }


//_______________________________________________________________________________________


void    _LikelihoodFunction::SetNthBit (long& reference, char n)
{
    long bitshifter = 1;
    bitshifter = bitshifter<<n;
    reference = reference|bitshifter;
}

//_______________________________________________________________________________________

bool    _LikelihoodFunction::CheckNthBit (long& reference, char n)
{
    unsigned long bitshifter = 1;
    bitshifter = bitshifter<<n;
    return reference&bitshifter;
}

//_______________________________________________________________________________________

char    _LikelihoodFunction::HighestBit (long reference)
{
    unsigned long bitshifter = 1,count = sizeof(long)*8-1;
    bitshifter = bitshifter<<(sizeof(long)*8-1);
    while(!(reference&bitshifter)) {
        bitshifter = bitshifter>>1;
        count--;
    }
    return count;
}

//_______________________________________________________________________________________

long    _LikelihoodFunction::HasHiddenMarkov (long reference, bool hmm) const {
    unsigned long bitshifter = 1, count = sizeof(long)*8-1;
    long     hMarkov = -1;
    bitshifter = bitshifter<<(sizeof(long)*8-1);
    while(bitshifter) {
        if (bitshifter&reference) {
            _CategoryVariable* thisC = (_CategoryVariable*)LocateVar(indexCat.lData[count]);
            if (hmm) {
                if (thisC->IsHiddenMarkov ()) {
                    hMarkov = count;
                }
            } else if (thisC->IsConstantOnPartition ()) {
                return count;
            }
        }
        bitshifter = bitshifter>>1;
        count--;
    }
    return hMarkov;
}

//_______________________________________________________________________________________

char    _LikelihoodFunction::LowestBit (long reference)
{
    unsigned long bitshifter = 1,count = 0;
    while(!(reference&bitshifter)) {
        bitshifter = bitshifter<<1;
        count++;
    }
    return count;
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::BuildIncrements (long reference, _SimpleList& incList)
{
    unsigned long currentIncr = 1;
    for (long i=0; i<indexCat.lLength; i++) {
        if (CheckNthBit(reference,i)) {
            incList<<currentIncr;
            currentIncr*=((_CategoryVariable*)LocateVar(indexCat(i)))->GetNumberOfIntervals();
        } else {
            incList<<0;
        }
    }
}

//_______________________________________________________________________________________

long    _LikelihoodFunction::MaximumDimension (void)
{
    long maxDim = 0;
    for (long i=0; i<theTrees.lLength; i++) {
        _Matrix *cM = (_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue();
        long myDim = cM->GetHDim()>cM->GetVDim()?cM->GetHDim():cM->GetVDim();
        if (myDim>maxDim) {
            maxDim = myDim;
        }
    }
    return maxDim;
}

//_______________________________________________________________________________________

long    _LikelihoodFunction::SequenceCount (long partID)
{
    if ((partID>=0)&&(partID<theTrees.lLength)) {
        _TheTree * cT = ((_TheTree*)(LocateVar(theTrees(partID))));
        HBLObjectRef  seqs = cT->TipCount();
        long res = seqs->Value();
        DeleteObject (seqs);
        return res;
    }
    return -1;
}

//_______________________________________________________________________________________

unsigned long    _LikelihoodFunction::SiteCount (void) const {

    unsigned long res = 0UL;

    for (unsigned long i=0UL; i<theDataFilters.lLength; i++) {
        res += GetIthFilter (i)->GetSiteCountInUnits();
    }

    return res;
}

//_______________________________________________________________________________________



//_______________________________________________________________________________________

void    _LikelihoodFunction::PrepareToCompute (bool disableClear) {
    if (hasBeenSetUp == 0L) {
        long categCount = 1L;


        for (unsigned long i=0UL; i<theTrees.lLength; i++) {
            _TheTree  * cT = GetIthTree(i);
            long locCC = cT->CountTreeCategories();
            StoreIfGreater( categCount, locCC );
            cT->SetUpMatrices (locCC);
        }

        for (unsigned long i=0UL; i<theProbabilities.lLength; i++) {
            ((_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue())->MakeMeSimple();
        }

        SetupCategoryCaches   ();
        SetupLFCaches         ();
        SetReferenceNodes     ();

        if (disableClear) {
            hasBeenSetUp = 0x6FFFFFFF;
        } else {
            hasBeenSetUp ++;
        }
        siteArrayPopulated = false;

    } else {
        hasBeenSetUp++;
    }
}

//_______________________________________________________________________________________

void    _LikelihoodFunction::DoneComputing (bool force)
{
    if (hasBeenSetUp == 1 || (hasBeenSetUp > 0 && force)) {
        for (unsigned long i=0UL; i<theTrees.lLength; i++) {
            GetIthTree(i)->CleanUpMatrices();
        }
        if (mstCache) {
            mstCache->resultCache.Clear();
            mstCache->statesCache.Clear();
        }
        for (unsigned long i=0UL; i<theProbabilities.lLength; i++) {
            ((_Matrix*)LocateVar(theProbabilities.lData[i])->GetValue())->MakeMeGeneral();
        }

        DeleteObject (siteResults);
        siteResults = 0;

        DeleteCaches        (false);
        categoryTraversalTemplate.Clear();
        hasBeenSetUp       = 0;
        siteArrayPopulated = false;
    } else if (hasBeenSetUp) {
        hasBeenSetUp --;
    }
}

//_______________________________________________________________________________________
void    _LikelihoodFunction::FillInConditionals(long partIndex) {
    if (partIndex >= 0) {

        long    catCounter = 0;
        _SimpleList             pcats;
        PartitionCatVars        (pcats,partIndex);
        catCounter            = pcats.lLength;

        _TheTree   *tree        = GetIthTree (partIndex);
        _DataSetFilter const *dsf = GetIthFilter (partIndex);;

        _SimpleList* tcc            = (_SimpleList*)treeTraversalMasks(partIndex);
        if (tcc) {
            long shifter = dsf->GetDimension()*dsf->GetPatternCount()*tree->GetINodeCount();
            for (long cc = 0; cc <= catCounter; cc++) {
                tree->FillInConditionals(dsf, conditionalInternalNodeLikelihoodCaches[partIndex] + cc*shifter, tcc);
            }
        }
    } else {
        for (long i = 0; i < theTrees.lLength; i++) {
            FillInConditionals (i);
        }
    }
}
//_______________________________________________________________________________________
void    _LikelihoodFunction::RankVariables(_AVLListX* tagger)
{
    _SimpleList varRank (indexInd.lLength,0,0),
                holder;

    gradientBlocks.Clear();

    if (tagger) {
        for (unsigned long k=0; k<indexInd.lLength; k++) {
            long idx = tagger->Find((BaseRef)indexInd.lData[k]);
            if (idx < 0) {
                ReportWarning (_String("Internal error in '_LikelihoodFunction::RankVariables': missing parameter name ") & *LocateVar(indexInd.lData[k])->theName);
            } else {
                varRank.lData[k] = -tagger->GetXtra(idx);
            }
        }
    }
    else {
        for (unsigned long k=0; k<indexInd.lLength; k++) {
            if (LocateVar(indexInd.lData[k])->IsGlobal()) {
                varRank<<10000;
            } else {
                varRank<<10050;
            }
            //printf ("%s -> %ld\n", LocateVar(indexInd.lData[k])->theName->sData, varRank.GetElement(-1));
        }

        for (unsigned long k=0; k<indexDep.lLength; k++) {
            holder.Clear();
            {
                _AVLList   al (&holder);
                LocateVar (indexDep.lData[k])->ScanForVariables(al,true);
                al.ReorderList ();
            }
            for (unsigned long j=0; j<holder.lLength; j++) {
                long f = indexInd.Find(holder.lData[j]);
                if (f>=0) {
                    varRank.lData[f]--;
                }
            }
        }
    }

    SortLists (&varRank,&indexInd);
    gradientBlocks.Clear();

    // enforce user provided rankings

    _AssociativeList * variableGrouping = (_AssociativeList*)FetchObjectFromVariableByType(&userSuppliedVariableGrouping, ASSOCIATIVE_LIST);
    if (variableGrouping) {

        _SimpleList  hist,
                     supportList;

        _AVLListX    existingRanking (&supportList);
        
        for (unsigned long vi = 0; vi < indexInd.lLength; vi ++ ) {
            existingRanking.Insert((BaseRef)indexInd.lData[vi], vi, true);
        }
        long  offset = 1;
        bool  re_sort = false;

        for (AVLListXLIteratorKeyValue key_value : variableGrouping->ListIterator()) {
            HBLObjectRef anEntry = (HBLObjectRef)key_value.get_object();
            if (anEntry->ObjectClass() == MATRIX) {
                _Matrix *variableGroup = (_Matrix*) anEntry;
                if (variableGroup -> IsAStringMatrix()) {
                    unsigned long dimension = variableGroup->GetHDim() * variableGroup->GetVDim ();

                    _SimpleList thisBlock;
                    for (unsigned long variable_id = 0; variable_id < dimension; variable_id ++) {
                        _String variableID ((_String*)variableGroup->GetFormula (variable_id,-1)->Compute()->toStr());
                        long variableIndex = LocateVarByName(variableID);
                        if (variableIndex >= 0) {
                            existingRanking.UpdateValue((BaseRef)variableIndex, -offset - dimension + variable_id, 1);
                            thisBlock << variableIndex;
                            //printf ("%s<%ld>\n",variableID.sData, variableIndex );
                            re_sort = true;
                        }
                    }
                    if (thisBlock.lLength) {
                        gradientBlocks && & thisBlock;
                    }
                    offset += dimension;
                }
            }
        }
        if (re_sort) {
            _SimpleList new_ranks;

            for (unsigned long vi = 0; vi < indexInd.lLength; vi ++ ) {
                new_ranks << existingRanking.GetXtra(existingRanking.Find ((BaseRef)indexInd.lData[vi]));
            }
            SortLists (&new_ranks,&indexInd);


            if (gradientBlocks.lLength) {
                _SimpleList  aux_list,
                             included (indexInd.lLength, 0,0),
                             not_listed;

                _AVLListX    indexIndToGlobalID (&aux_list);

                for (unsigned long vi = 0; vi < indexInd.lLength; vi ++ ) {
                    indexIndToGlobalID.Insert((BaseRef)indexInd.lData[vi], vi, true);
                    //printf ("[%ld]\n",indexInd.lData[vi] );
                }

                for (long b = 0; b < gradientBlocks.countitems(); b++) {
                    _SimpleList *a_block = (_SimpleList*)(gradientBlocks(b));
                    for (long i = 0; i < a_block->countitems(); i++){
                        long t = indexIndToGlobalID.Find ((BaseRef)a_block->lData[i]);
                        //printf ("%ld %ld\n",i,t);
                        if (t >= 0) {
                            t = indexIndToGlobalID.GetXtra(t);
                            //printf ("%ld %ld %ld/ %ld /%ld\n", b, i, a_block->lData[i], indexIndToGlobalID.Find ((BaseRef)a_block->lData[i]), t);
                            a_block->lData[i] = t;
                            included.lData  [t] = 1;
                        } else {
                            a_block->Delete (i--);
                        }
                   }
                    if (a_block->lLength == 0) {
                        gradientBlocks.Delete(b--);
                    }
                }

                if (gradientBlocks.lLength) {
                    for (long t = 0; t < included.lLength; t++) {
                        if (included.lData[t]==0) {
                            not_listed << t;
                        }
                    }
                    if (not_listed.lLength) {
                        gradientBlocks && & not_listed;
                    }
                }

                /*for (long b = 0; b < gradientBlocks.lLength; b++) {
                    _SimpleList *a_block = (_SimpleList*)(gradientBlocks(b));
                    for (long i = 0; i < a_block->lLength; i++){
                        printf ("Block %ld variable %s\n", b, LocateVar(indexInd.lData[a_block->lData[i]])->GetName()->sData);
                    }
                }*/

            }

        }
    }

}


//_______________________________________________________________________________________

hyFloat _CustomFunction::Compute (void) {
    likeFuncEvalCallCount++;
    _SimpleList const * iv = &GetIndependentVars ();
    for (unsigned long i=0UL; i<iv->lLength; i++) {
        hyFloat result = GetIthIndependent(i);

        if (result<GetIthIndependentBound (i,true) || result>GetIthIndependentBound (i,false)) {
            return -INFINITY;
        }
    }

    HBLObjectRef res = myBody.Compute();
    if (res) {
        return res->Value();
    }
    return 0.0;
}
