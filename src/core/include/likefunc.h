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

#ifndef __LIKELIHOODF__
#define __LIKELIHOODF__

//#pragma once
#include "category.h"
#include "calcnode.h"

#ifdef __HYALTIVEC__
#define   A_LARGE_NUMBER          1.e35
#else
#define   A_LARGE_NUMBER          1.e100
#endif

#define   DEFAULTPARAMETERLBOUND  0.0
#define   DEFAULTPARAMETERUBOUND  10000.0

/* various run modes for PopulateConditionalProbabilities */

#define   _hyphyLFConditionProbsRawMatrixMode       0
#define   _hyphyLFConditionProbsScaledMatrixMode    1
#define   _hyphyLFConditionProbsWeightedSum         2
#define   _hyphyLFConditionProbsMaxProbClass        3
#define   _hyphyLFConditionProbsClassWeights        4
#define   _hyphyLFConditionMPIIterate               5

/* computational template kinds for the likelihood function */

#define   _hyphyLFComputationalTemplateNone         0
#define   _hyphyLFComputationalTemplateBySite       1
#define   _hyphyLFComputationalTemplateByPartition  2

/* conditional likelihood matrix reconstruction modes */

#define   _hyphyLFConstructCategoryMatrixConditionals   0
#define   _hyphyLFConstructCategoryMatrixClasses        1
#define   _hyphyLFConstructCategoryMatrixWeights        2
#define   _hyphyLFConstructCategoryMatrixPosteriors     3
#define   _hyphyLFConstructCategoryMatrixSiteProbabilities      4

/* likelihood seialization model */

#define   _hyphyLFSerializeModeVanilla                  0
#define   _hyphyLFSerializeModeOptimize                 1
#define   _hyphyLFSerializeModeLongMPI                  2
#define   _hyphyLFSerializeModeCategoryAsGlobal         4
#define   _hyphyLFSerializeModeShortMPI                 5

/* likelihood function parallelize mode */

#define   _hyphyLFMPIModeNone                           0
#define   _hyphyLFMPIModePartitions                     1
#define   _hyphyLFMPIModeSiteTemplate                   2
#define   _hyphyLFMPIModeREL                            3
#define   _hyphyLFMPIModeAuto                           4

/* partition category variable types */

#define   _hyphyCategoryNormal                          0x1
#define   _hyphyCategoryHMM                             0x2
#define   _hyphyCategoryCOP                             0x4

/* interval mapping functions */

#define   _hyphyIntervalMapID                           0x0 // identify
#define   _hyphyIntervalMapExpit                        0x1 // expit -- maps [-infty,infty) to [0,1] 
// 1/(1+exp[-x])
#define   _hyphyIntervalMapSqueeze                      0x2 // maps [0,infty) to [0,1)
// x / (1+x)


//_______________________________________________________________________________________

struct  MSTCache {
    _List               computingOrder,
                        storageOrder,
                        referenceOrder,
                        parentOrder,
                        stashedLeafOrders;

    _SimpleList         statesNCache,
                        resultCache,
                        statesCache,
                        cacheSize;

};

//_______________________________________________________________________________________

class   _LikelihoodFunction: public BaseObj
{

public:

    // constructors

    _LikelihoodFunction (void); // default - doesn't do much
    _LikelihoodFunction (_String&, _VariableContainer*); // from triplets
    _LikelihoodFunction (_LikelihoodFunction&); // stack copy
    void        Init      (void);

    bool        Construct (_List&,_VariableContainer*);

    virtual ~_LikelihoodFunction (void) {
        Cleanup();   // destructor
    }

    virtual BaseRef     toStr (void);

    virtual BaseRef     makeDynamic (void);      // dynamic copy of this object

    virtual void        Duplicate (BaseRef);         // duplicate an object into this one

    _SimpleList&GetIndependentVars (void); // return a list of all indepenent variables
    _SimpleList&GetDependentVars   (void); // get all dependent vars of this object
    _SimpleList&GetCategoryVars    (void); // get all category variables
    void        GetGlobalVars      (_SimpleList&);

    long        CountObjects       (char);     // return object count
    // 0 - partitions
    // 1 - global variables
    // 2 - local independents
    // 3 - dependents
    // 4 - category variables


    _Parameter  GetIthIndependent           (long);     // get the value of i-th independent variable
    _Parameter  GetIthDependent             (long);     // get the value of i-th dependent variable
    void        GetAllIndependent           (_Matrix&); // store all indepenent values in a matrix
    _Variable*  GetIthIndependentVar        (long);     // get the variable object of i-th independent variable
    _Variable*  GetIthDependentVar          (long);     // get the variable object of i-th dependent variable
    _Parameter  GetIthIndependentBound      (long, bool isLower = true);
    // get the lower / upper bound for the i-th indepdendent variable

    void        SetIthIndependent (long, _Parameter);           // set the value of i-th independent variable
    bool        CheckAndSetIthIndependent (long, _Parameter);   // set the value of i-th independent variable
    void        SetIthDependent           (long, _Parameter);   // set the value of i-th dependent variable
    bool        IsIthParameterGlobal      (long);

    long        SetAllIndependent         (_Matrix*);


    void        UpdateIndependent (long,bool,_SimpleList* = nil,_SimpleList* = nil);
    void        UpdateDependent (long);
    _AssociativeList*
    SimulateCodonNeutral (_Matrix*, _Matrix*, long);

    bool        PreCompute      (void);
    void        PostCompute     (void);
    virtual
    _Parameter  Compute         (void);

    void        PrepareToCompute (bool = false);
    void        DoneComputing    (bool = false);
    virtual
    _Matrix*    Optimize ();
    _Matrix*    ConstructCategoryMatrix     (const _SimpleList&, char, bool = true, _String* = nil);

    _Parameter  SimplexMethod               (_Parameter& precision);
    void        Anneal                      (_Parameter& precision);

    void        Simulate                    (_DataSet &,_List&, _Matrix* = nil, _Matrix* = nil, _Matrix* = nil, _String* = nil);

    void        ReconstructAncestors        (_DataSet &, _SimpleList&, _String&, bool = false, bool = false, bool = false);
    // 20090224: added an argument to allow the marginal state reconstruction
    // 20091009: added an argument to allow the reconstruction of leaves

    long        MaximumDimension            (void);

    virtual _PMathObj   CovarianceMatrix            (_SimpleList* = nil);

    // compute  covariance matrix  based on the Hessian
    // optional list of parameters to estimate the conditional covariance for


    virtual     void        RescanAllVariables      (void);

    long        DependOnTree            (_String&);
    long        DependOnModel           (_String&);
    long        DependOnDS              (long);
    bool        DependOnDF              (long ID) {
        return theDataFilters.Find(ID)>=0;
    }
    bool        MapTreeTipsToData       (long, bool leafScan = false);
    bool        UpdateFilterSize        (long);
    void        VoidOldResults          (void) {
        computationalResults.ZeroUsed();
    }
    _CategoryVariable*
    FindCategoryVar         (long);
    // return the category variable for a given partition
    void        RankVariables           (_AVLListX* tagger = nil);
    _SimpleList&GetTheTrees             (void) {
        return theTrees;
    }
    _SimpleList&GetTheFilters           (void) {
        return theDataFilters;
    }
    _SimpleList&GetBaseFreqs            (void) {
        return theProbabilities;
    }
    
    void        FillInConditionals      (long = -1);
    
    void        Setup                   (void);
    bool&       HasBeenOptimized (void) {
        return hasBeenOptimized;
    }
    void        ComputePruningEfficiency(long&, long&);

    long        SequenceCount           (long);
    long        SiteCount               (void);
    void        Rebuild                 (void);
    void        SerializeLF             (_String&, char=0, _SimpleList* = nil, _SimpleList* = nil);
    _Formula*   HasComputingTemplate    (void) {
        return computingTemplate;
    }
    void        StateCounter            (long);
    void        MPI_LF_Compute          (long, bool = false);

#if defined _SLKP_LFENGINE_REWRITE_
#if defined _OPENMP
    void        SetThreadCount            (long tc) {
        if (tc != lfThreadCount) {
            lfThreadCount = tc;
            FillInConditionals ();
        }
    }
    long        GetThreadCount            (void) {
        return lfThreadCount;
    }
#else
    long        GetThreadCount            (void) {
        return 1;
    }
#endif
#endif

    bool            ProcessPartitionList        (_SimpleList&, _Matrix*, _String);
    // given a matrix argument (argument 2; can be nil to include all)
    // populate a sorted list (argument 1)
    // of partitions indexed in the matrix (e.g. {{1,3}} would include the 2-nd and 4-th partitions
    // the third argument is the ID of the calling function for error reporting
    // returns true if at least one of the requested partitions is valid
    // otherwise returns false

    long            TotalRateClassesForAPartition (long, char = 0);
    // given a partition index (assuming that category caches have been set up)
    // returns how many total rate classes there are for this partition
    // if the index is < 0, return the total number of categories for the LF as a whole

    // char flag: 0 - all categories
    //            1 - only HMM
    //            2 - only constant on partition


    bool            HasBeenSetup                (void) {
        return hasBeenSetUp > 0;
    }

    /**
       @return  an associative array with a variety of metrics (see GetString on HyPhy wiki)
                about this LF.
       @author  SLKP
       @version 20110608
    */

    _AssociativeList*CollectLFAttributes         (void);
protected:

    _Matrix*        PairwiseDistances       (long index);
    void            CheckDependentBounds    (void);
    void            AllocateSiteResults     (void);
    void            ZeroSiteResults         (void);
    // 20090211: A utility function to reset site results.
    void            GetInitialValues        (void);
    bool            checkPermissibility     (_Matrix&m, long row);
    _Parameter      computeAtAPoint         (_Matrix&m, long row = 0);
    _Parameter      replaceAPoint           (_Matrix&m, long row, _Matrix&p, _Parameter& nV, _Matrix& fv);

    void            ScanAllVariablesOnPartition
    (_SimpleList&, _SimpleList&, _SimpleList&, _SimpleList&, bool = false);
    // internal function to scan all the on a set of partitions variables in
    // 20091103: SLKP
    // the boolean flag (if true), requests that only tree-associated variables
    // (no templates or frequencies) be scanned

    virtual void            ScanAllVariables        (void);
    // internal function to scan all the variables in

    void            OptimalOrder            (long, _SimpleList&);
    // determine the optimal order of compuation for a block

    _Parameter      ComputeBlock            (long, _Parameter* siteResults = nil, long currentRateClass = -1, long = -1, _SimpleList* = nil);
    // 20090224: SLKP
    // added the option to pass an interior branch (referenced by the 3rd argument in the same order as flatTree)
    // and a set of values for each site pattern (indexed left to right) in the 4th argument

    void            SetReferenceNodes       (void);
    // compute likelihood over block index i

    // internal Setup, called from within constructors

//      void            DumpingOrder            (long index, _SimpleList& sl);
    // used for FreeUpMemory

    void            Cleanup                 (void);
    // internal cleanup, called from within destructors

    void            Clear                   (void);


    long            Bracket                 (long , _Parameter& , _Parameter& , _Parameter& ,
            _Parameter& , _Parameter& , _Parameter& , _Parameter&, _Matrix* = nil);
    //long          GradientBracketOneVar (_Matrix&, _Matrix& , _Matrix& , _Matrix&,  _Parameter& ,
    //                                      _Parameter&, _Parameter&, _Parameter&, bool retry = false);
    void            LocateTheBump         (long,_Parameter , _Parameter& , _Parameter&, _Parameter = -1.);
    void            GradientLocateTheBump (_Parameter, _Parameter&, _Matrix&, _Matrix&);
    void            GradientDescent       (_Parameter& , _Matrix& );
    void            ConjugateGradientDescent
    (_Parameter , _Matrix& , bool localOnly = false, long = 0x7fffffff,_SimpleList* only_these_parameters = nil, _Parameter check_lf = A_LARGE_NUMBER);

    _Parameter      SetParametersAndCompute
    (long, _Parameter, _Matrix* = nil, _Matrix* = nil);

    long            CostOfPath            (_DataSetFilter*, _TheTree* , _SimpleList&, _SimpleList* = nil);

    void            BuildLeafProbs        (node<long>& , long*, long&, _DataSet&, _TheTree*, long&, bool, long, _DataSetFilter*, long, _DataSet* = nil);
    bool            SingleBuildLeafProbs  (node<long>&, long, _SimpleList&, _SimpleList&, _TheTree*, bool,_DataSetFilter*, _SimpleList* = nil);
    void            CodonNeutralSimulate  (node<long>&, long, bool,_Matrix*,_Matrix*, _Parameter&, _Parameter&);

    bool            HasBlockChanged       (long);
    long            BlockLength           (long);
    void            PartitionCatVars      (_SimpleList&, long);
    // 20090210: extract variable indices for category variables in i-th partition
    // and append them to _SimpleList


    static  void            RandomizeList               (_SimpleList&, long);
    static  void            CheckFibonacci              (_Parameter);



    long            PartitionLengths            (char = 0,  _SimpleList const* = nil);
    /*
        SLKP: 20090317

        argument = 0
            return the length (sites) of the longest part in this
            likelihood function
        argument = 1
            return the cumulative length (sites) of all parts in this
            likelihood function

        SLKP: 20090608

            argument 2 provides an optional subcollection of partititons to operate on;
                the default is to operate on all
    */


private:


    void            ComputeParameterPenalty     (void);
    
    
    
    bool            SendOffToMPI                (long);
    void            InitMPIOptimizer            (void);
    void            CleanupMPIOptimizer         (void);
    void            ComputeBlockInt1            (long,_Parameter&,_TheTree*,_DataSetFilter*, char);
    void            CheckStep                   (_Parameter&, _Matrix, _Matrix* selection = nil);
    void            GetGradientStepBound        (_Matrix&, _Parameter &, _Parameter &, long* = nil);
    void            ComputeGradient             (_Matrix&, _Matrix&,  _Parameter&, _Matrix&, _SimpleList&,
            long, bool normalize = true);
    bool            SniffAround                 (_Matrix& , _Parameter& , _Parameter&);
    long            HasPrecisionBeenAchieved    (_Parameter funcValue = 2.*A_LARGE_NUMBER, bool = false);
    void            RecurseCategory             (long,long,long,long,_Parameter
#ifdef _SLKP_LFENGINE_REWRITE_
            ,_SimpleList* = nil, char = 0, _Parameter* = nil,
            long = -1, _SimpleList* = nil
#endif
                                                );
    void            RecurseConstantOnPartition  (long, long, long, long, _Parameter, _Matrix&);



    void            SetNthBit                   (long&,char);
    bool            CheckNthBit                 (long&,char);
    void            BuildIncrements             (long, _SimpleList&);
    char            HighestBit                  (long);
    char            LowestBit                   (long);
    long            HasHiddenMarkov             (long, bool hmm = true);
    _Matrix*        RemapMatrix                 (_Matrix*, const _SimpleList&);
    /* given a matrix where each column corresponds to a site pattern
     and a list of partitions that are covered
     this function returns a new matrix where each pattern is mapped
     to the corresponding sites in the alignment
    */
    void            CleanUpOptimize             (void);
    void            ComputeBlockForTemplate     (long, bool = false);
    void            ComputeBlockForTemplate2    (long, _Parameter*, _Parameter*, long);
    void            DeleteCaches                (bool = true);
    void            PopulateConditionalProbabilities
    (long index, char runMode, _Parameter* buffer, _SimpleList& scalers, long = -1, _SimpleList* = nil);
    void            ComputeSiteLikelihoodsForABlock
    (long, _Parameter*, _SimpleList&, long = -1, _SimpleList* = nil,  char = 0);

    // this function computes a list of site probabilities for the i-th block (1st parameter)
    // stores them in pattern (left to right) order (2nd argument)
    // and writes scaling factors for each site into the third argument
    // arguments 4 (branch index) and 5 (pattern assignments)
    // allows the calculation of the probability vector while setting a specific interior branch
    // to a given sequence

    _Parameter          SumUpHiddenMarkov (const _Parameter *, _Matrix&, _Matrix&, _SimpleList *, const _SimpleList*, long);
    /*
        SLKP 20090420

        given conditional probs                         (1)
        the transition matrix for the HMM               (2)
        the initial probabilities for the HMM           (3)
        the duplicate map (or nil if no need to remap)  (4)
        the vector of by-site scaling factors as        (5)
            a concatenated _SimpleList if (4) is not nil or
            as a _List of _SimpleLists if (4) is nil
        partition length                                (6)

        compute the log likelihood of the partition using the forward HMM algorithm with scaling
     */

    void                    RunViterbi (_Matrix & , const _Parameter * , _Matrix& , _Matrix& , _SimpleList * ,  const _SimpleList* , long );
    /* Viterbi decoding for HMM; parameter meanings as in SumUpHiddenMarkov,
       except the first, which will store the optimal path to be returned */

    _Parameter              SumUpSiteLikelihoods        (long, const _Parameter*, const _SimpleList&);
    /*
     SLKP 20090318

     given a partition index (argument 1),
     a vector of _pattern_ likelihoods (argument 2) and
     a list of pattern scaling factors (argument 3)

     compute the log likelihood of the partition

     */


    void            UpdateBlockResult           (long, _Parameter);
    /*
        SLKP 20090318

        update cached log-likelihood value for a given
        partition
    */


    _List*          RecoverAncestralSequencesMarginal
    (long, _Matrix&,_List&, bool = false);
    void            DetermineLocalUpdatePolicy  (void);
    void            FlushLocalUpdatePolicy      (void);
    void            RestoreScalingFactors       (long, long, long, long*, long *);
    void            SetupLFCaches               (void);
    void            SetupCategoryCaches         (void);
    bool            HasPartitionChanged         (long);
    void            SetupParameterMapping       (void);
    void            CleanupParameterMapping     (void);

    _SimpleList     theTrees,
                    theDataFilters,
                    theProbabilities,
                    indexInd,
                    indexDep,
                    indexCat,
                    *nonConstantDep,
                    blockDependancies,
                    parameterTransformationFunction;
    /* 20110718: SLKP this list holds the index of the parameter interval mapping function
        used during optimization */

    _GrowingVector  computationalResults;

    _List           optimalOrders,
                    leafSkips,
                    categoryTraversalTemplate,
                    /*SLKP: 20090225

                         This list contains as many entries (themselves of type _List) as there are partitions
                         The entry will be empty for a partition without category variables
                         For a partition with N category variables the entry will contain a

                         1). _List of references to category variables themselves in the order that
                         they appear in the blockDependancies

                         2). _SimpleList of category counts for each variable C_1, C_2, ... C_N + a final entry with
                         C_1 * C_2 * ... * C_N -- the total number of rate categories for this partition

                         3). _SimpleList of incremental loop offsets for each variable, e.g. if the
                         _SimpleList in 2). is 2,3,4, then this _SimpleList is
                         12,4,1

                      SLKP: 20100416

                         4). A _SimpleList including _relative_ indexing for the buffer array
                             for HMM and "constant on partition" variables similar to 3, but
                             ignoring all other variables. The last entry is the total number
                             of HMM (or similar) rate classes on a partition

                         5). A _SimpleList storing the type of each cat variable in a partition (flag)
                             the last element is the cumulative flag composed of bit toggles to indicate
                             what type of category variables there are for a given partition:
                                _hyphyCategoryNormal
                                _hyphyCategoryHMM
                                _hyphyCategoryCOP
                     */
                    indVarsByPartition,
                    depVarsByPartition;



    long            evalsSinceLastSetup,
                    // this is necessary to force internal caches to be filled in at least once
                    // for trees where branches have different numbers of rate categories
                    templateKind,

                    /*
                        _hyphyLFComputationalTemplateNone:
                            no computational template: simply return the summed up log-likelihood
                            for all partitions in the likelihood function

                        _hyphyLFComputationalTemplateBySite:
                            all partitons have the same number of sites
                            the template is in terms of SITE_LIKELIHOOD
                            meant to implement model/tree mixture constructs

                        _hyphyLFComputationalTemplateByPartition
                            the template is in terms of BLOCK_LIKELIHOOD
                            meant to implement random effects on the level of partitions
                            models

                        <0                                          : (-1-variable index) if a hidden markov process on site likelihoods
                        >_hyphyLFComputationalTemplateByPartition   : a _hyphyLFComputationalTemplateByPartition shifted index of a user
                                                                      function that will assemble the likelihood
                                                                      from site conditionals, written into matrices
                                                                      SITE_LIKELIHOODS_0, SITE_LIKELIHOODS_1, ..
                    */
                    hasBeenSetUp;

    _Matrix         *siteResults,
                    *bySiteResults,
                    *parameterValuesAndRanges;

    /* 20110718 SLKP
        this Nx4 (N is the number of adjustable parameters) matrix
        stores variable bounds and values for each independent variable;
        this is used for internal normalization transforms during optimization
     */

    bool            hasBeenOptimized,
                    siteArrayPopulated;

    _Formula*       computingTemplate;
    MSTCache*       mstCache;
    
    _Parameter      smoothingTerm,
                    smoothingReduction,
                    smoothingPenalty;

#ifdef  _SLKP_LFENGINE_REWRITE_
    /*
       these variables store conditional likelihoods for every paritition
       internal node (post-order traversal)
       and leaves;

       -conditionalInternalNodeLikelihoodCaches: 3D vector

            1st coordinate - the node index (in post-order traversal)
            2nd coordinate - the site (unique pattern) index (left-to-right)
            3rd coordiante - i-th marginal (for the i-th character; 0-filterCharDimension)

            stores the probability for the subtree below this node given that i-th character
            is present at the node

            for partitions with category variables this vector is copied enough times to
            store every combination of category values

        -conditionalTerminalNodeStateFlag: 2D vector
            1st coordinate - the leaf index (in post-order traversal)
            2nd coordinate - the site (unique pattern) index (left-to-right)

            stores a non-negative integer for a leaf with a simple (non-ambiguous) character
            stores a negative value to index (-value-1) into the vector of conditionalTerminalNodeLikelihoodCaches
            and read off filterCharDimension characters from there
    */

    _Parameter**        conditionalInternalNodeLikelihoodCaches,
               **     siteScalingFactors,
               **     branchCaches;

    _List               conditionalTerminalNodeLikelihoodCaches;
    long      **        conditionalTerminalNodeStateFlag;

    _SimpleList         overallScalingFactors,
                        overallScalingFactorsBackup,
                        // this (and siteCorrectionsBackup)
                        // is needed to store site/partition scaling factors when treeTraversalMasks
                        // and branch caching are used, otherwise scaling done by ComputeBranchCaches
                        // will fubar the scaling for ComputeTreeBlockByBranch depending on site
                        // patterns
                        canUseReversibleSpeedups,
                        // a partition will be tagged with '1' if its tree has only
                        // time-reversible models
                        siteScalerBuffer
                        // used for LF with category variables to
                        // store site-by-site scaling factors
                        ;

    _List               localUpdatePolicy,
                        matricesToExponentiate,
                        treeTraversalMasks,
                        computedLocalUpdatePolicy,
                        siteCorrections,
                        siteCorrectionsBackup,
                        cachedBranches,
                        // for models with categories, a list of site by site scaling operation counts
                        partScalingCache,
                        // used to store site by site scalers in computations that are performed
                        // on a site-by-site basis; includes scratch cache for remapping
                        gradientBlocks
                        ;

#ifdef  _OPENMP
    long                lfThreadCount;
#endif

#endif
};

//_______________________________________________________________________________________

class   _CustomFunction: public _LikelihoodFunction
{

public:

    _CustomFunction         (_String*);

    virtual     _Parameter  Compute                 (void);
    virtual     void        RescanAllVariables      (void) {}



    _Formula myBody;
};


//_______________________________________________________________________________________

extern  bool forceRecomputation,
        isInOptimize;

extern  long lockedLFID;

extern  _String // declare shared global keys/settings
globalStartingPoint             ,
randomStartingPerturbations    ,
optimizationPrecision          ,
startingPrecision              ,
optimizationMethod                 ,
useLastResults                     ,
allowBoundary                  ,
bracketingPersistence          ,
intermediatePrecision          ,
keepOptimalOrder               ,
skipOmissions                  ,
optimizeSummationOrder             ,
optimizePartitionSize          ,
maximumIterationsPerVariable   ,
optimizationPrecisionMethod    ,
relativePrecision              ,
likefuncOutput                 ,
dataFilePrintFormat                ,
dataFileDefaultWidth           ,
dataFileGapWidth               ,
categorySimulationMethod       ,
useInitialDistanceGuess            ,
randomSeed                     ,
assignedSeed                   ,
covariancePrecision                ,
cacheSubtrees                  ,
likeFuncCountVar               ,
doShuffleOrder                 ,
forceDistanceEstimates         ,
useDuplicateMatrixCaching      ,
siteWiseMatrix                 ,
blockWiseMatrix                    ,
useFullMST                     ,
stateCountMatrix               ,
wStateCountMatrix              ,
tryNumericSequenceMatch            ,
allowSequenceMismatch          ,
shortMPIReturn                 ,
mpiPrefixCommand               ,
skipConjugateGradient          ,
useIntervalMapping             ,
intervalMappingMethod          ,
useAdaptiveVariableStep            ,
storeRootSupportFlag           ,
supportMatrixVariable          ,
optimizationStatusFile         ,
autoParalellizeLF              ,
addLFSmoothing                 ,
reduceLFSmoothing              ;



bool    CheckEqual                  (_Parameter, _Parameter);
void    StateCounterResultHandler   (_Formula&, _SimpleList*,long&,long&,long,_Matrix&,_Matrix&);

_LikelihoodFunction*
FindLikeFuncByName           (_String&);

extern  bool                usedCachedResults;

extern _Parameter           _lfScalerUpwards,
       _lfScalingFactorThreshold,
       _logLFScaler;

extern  _GrowingVector      _scalerMultipliers,
        _scalerDividers;

_Parameter                  acquireScalerMultiplier (long);
_Parameter                  myLog                   (_Parameter);
_Parameter                  mapParameterToInverval  (_Parameter, char, bool);

#ifdef  __HYPHYMPI__
extern                  _Matrix     resTransferMatrix;
extern                  long        hyphyMPIOptimizerMode;
extern                  _String     mpiLoopSwitchToOptimize,
                        mpiLoopSwitchToBGM;

extern                  _SimpleList mpiNodesThatCantSwitch;

long                    RetrieveMPICount            (char);
void                    MPISwitchNodesToMPIMode     (long);


#endif

#endif









