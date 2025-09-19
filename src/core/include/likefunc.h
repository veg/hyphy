/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

// #pragma once
#include "associative_list.h"
#include "category.h"
#include "matrix.h"
#include "tree.h"
#include "vector.h"

// #define   NUMERICAL_INFINITY          INFINITY

#define DEFAULTPARAMETERLBOUND 0.0
#define DEFAULTPARAMETERUBOUND 10000.0

/** @name Run modes for PopulateConditionalProbabilities
    @{
*/
#define _hyphyLFConditionProbsRawMatrixMode 0
#define _hyphyLFConditionProbsScaledMatrixMode 1
#define _hyphyLFConditionProbsWeightedSum 2
#define _hyphyLFConditionProbsMaxProbClass 3
#define _hyphyLFConditionProbsClassWeights 4
#define _hyphyLFConditionMPIIterate 5
///@}

/** @name Computational template kinds for the likelihood function
    @{
*/
#define _hyphyLFComputationalTemplateNone 0
#define _hyphyLFComputationalTemplateBySite 1
#define _hyphyLFComputationalTemplateByPartition 2
///@}

/** @name Conditional likelihood matrix reconstruction modes
    @{
*/
#define _hyphyLFConstructCategoryMatrixConditionals 0
#define _hyphyLFConstructCategoryMatrixClasses 1
#define _hyphyLFConstructCategoryMatrixWeights 2
#define _hyphyLFConstructCategoryMatrixPosteriors 3
#define _hyphyLFConstructCategoryMatrixSiteProbabilities 4
#define _hyphyLFConstructCategoryMatrixPartitions 5
///@}

/** @name Likelihood serialization model
    @{
*/
#define _hyphyLFSerializeModeVanilla 0
#define _hyphyLFSerializeModeOptimize 1
#define _hyphyLFSerializeModeLongMPI 2
#define _hyphyLFSerializeModeCategoryAsGlobal 4
#define _hyphyLFSerializeModeShortMPI 5
///@}

/** @name Likelihood function parallelize mode
    @{
*/
#define _hyphyLFMPIModeNone 0
#define _hyphyLFMPIModePartitions 1
#define _hyphyLFMPIModeSiteTemplate 2
#define _hyphyLFMPIModeREL 3
#define _hyphyLFMPIModeAuto 4
///@}

/** @name Partition category variable types
    @{
*/
#define _hyphyCategoryNormal 0x1
#define _hyphyCategoryHMM 0x2
#define _hyphyCategoryCOP 0x4
///@}

/** @name Interval mapping functions
    @{
*/
#define _hyphyIntervalMapID 0x0    // identify
#define _hyphyIntervalMapExpit 0x1 // expit -- maps [-infty,infty) to [0,1]
// 1/(1+exp[-x])
#define _hyphyIntervalMapSqueeze 0x2 // maps [0,infty) to [0,1)
// x / (1+x)
///@}

//_______________________________________________________________________________________

/** @brief A cache for MST computation */
struct MSTCache {
  _List computingOrder, storageOrder, referenceOrder, parentOrder,
      stashedLeafOrders;

  _SimpleList statesNCache, resultCache, statesCache, cacheSize;
};

/** @brief Likelihood function count types */
enum _LikelihoodFunctionCountType {
  kLFCountPartitions,
  kLFCountGlobalVariables,
  kLFCountLocalCariables,
  kLFCountDependentVariables,
  kLFCountCategoryVariables
};

//_______________________________________________________________________________________

/**
 * @brief A class for likelihood functions
 *
 */
class _LikelihoodFunction : public BaseObj {

public:
  /**
   * @brief Construct a new _LikelihoodFunction object
   *
   */
  _LikelihoodFunction(void);
  /**
   * @brief Construct a new _LikelihoodFunction object
   *
   * @param s The string to construct from
   * @param vc The variable container
   */
  _LikelihoodFunction(_String &, _VariableContainer *); // from triplets
  /**
   * @brief Construct a new _LikelihoodFunction object
   *
   * @param lf The likelihood function to copy
   */
  _LikelihoodFunction(_LikelihoodFunction const &);     // stack copy
  /**
   * @brief The assignment operator
   *
   * @param lf The likelihood function to copy
   * @return const _LikelihoodFunction&
   */
  const _LikelihoodFunction &operator=(_LikelihoodFunction const &);
  /**
   * @brief Initialize the likelihood function
   *
   */
  void Init(void);

  /**
   * @brief Construct the likelihood function
   *
   * @param l The list of triplets
   * @param vc The variable container
   * @return true if the construction was successful, false otherwise
   */
  bool Construct(_List &, _VariableContainer *);

  /**
   * @brief Destroy the _LikelihoodFunction object
   *
   */
  virtual ~_LikelihoodFunction(void) {
    Cleanup(); // destructor
  }

  /**
   * @brief Convert the likelihood function to a string
   *
   * @param ul The format to use
   * @return BaseRef The string representation of the likelihood function
   */
  virtual BaseRef toStr(unsigned long = 0UL);

  /**
   * @brief Make a dynamic copy of the likelihood function
   *
   * @return BaseRef The dynamic copy of the likelihood function
   */
  virtual BaseRef makeDynamic(void) const; // dynamic copy of this object

  /**
   * @brief Duplicate the likelihood function
   *
   * @param brc The likelihood function to duplicate
   */
  virtual void Duplicate(BaseRefConst); // duplicate an object into this one

  /**
   * @brief Get the Independent Vars object
   *
   * @return _SimpleList const&
   */
  _SimpleList const &
  GetIndependentVars(void) const; // return a list of all indepenent variables
  /**
   * @brief Get the Dependent Vars object
   *
   * @return _SimpleList const&
   */
  _SimpleList const &
  GetDependentVars(void) const; // get all dependent vars of this object
  /**
   * @brief Get the Category Vars object
   *
   * @return _SimpleList const&
   */
  _SimpleList const &GetCategoryVars(void) const; // get all category variables
  /**
   * @brief Get the Global Vars object
   *
   * @param sl The simple list to store the global variables in
   */
  void GetGlobalVars(_SimpleList &) const;

  /**
   * @brief Count the number of objects of a given type
   *
   * @param type The type of object to count
   * @return unsigned long The number of objects of the given type
   */
  unsigned long
      CountObjects(_LikelihoodFunctionCountType) const; // return object count
  // 0 - partitions
  // 1 - global variables
  // 2 - local independents
  // 3 - dependents
  // 4 - category variables

  /**
   * @brief Get the Ith Independent object
   *
   * @param i The index of the independent variable
   * @param b Whether to check for bounds
   * @return hyFloat The value of the independent variable
   */
  hyFloat GetIthIndependent(
      long, bool = true) const; // get the value of i-th independent variable
  /**
   * @brief Get the Ith Independent Name object
   *
   * @param i The index of the independent variable
   * @return const _String* The name of the independent variable
   */
  const _String *GetIthIndependentName(
      long) const; // get the name of i-th independent variable
  /**
   * @brief Get the Ith Dependent Name object
   *
   * @param i The index of the dependent variable
   * @return const _String* The name of the dependent variable
   */
  const _String *
  GetIthDependentName(long) const; // get the name of i-th independent variable
  /**
   * @brief Get the Ith Dependent object
   *
   * @param i The index of the dependent variable
   * @return hyFloat The value of the dependent variable
   */
  hyFloat
  GetIthDependent(long) const; // get the value of i-th dependent variable
  /**
   * @brief Get the All Independent object
   *
   * @param m The matrix to store the independent variables in
   */
  void
  GetAllIndependent(_Matrix &) const; // store all indepenent values in a matrix
  /**
   * @brief Get the Ith Independent Var object
   *
   * @param i The index of the independent variable
   * @return _Variable* The independent variable
   */
  _Variable *GetIthIndependentVar(
      long) const; // get the variable object of i-th independent variable
  /**
   * @brief Get the Ith Dependent Var object
   *
   * @param i The index of the dependent variable
   * @return _Variable* The dependent variable
   */
  _Variable *GetIthDependentVar(
      long) const; // get the variable object of i-th dependent variable
  /**
   * @brief Get the Ith Category Var object
   *
   * @param i The index of the category variable
   * @return _CategoryVariable* The category variable
   */
  _CategoryVariable *GetIthCategoryVar(
      long) const; // get the variable object of i-th category variable
  /**
   * @brief Get the Ith Independent Bound object
   *
   * @param i The index of the independent variable
   * @param isLower Whether to get the lower or upper bound
   * @return hyFloat The bound of the independent variable
   */
  hyFloat GetIthIndependentBound(long, bool isLower = true) const;
  // get the lower / upper bound for the i-th indepdendent variable

  /**
   * @brief Obtain the g'(x) for the chain rule differentiation under transformed variables
   *
   * @param i The index of the independent variable
   * @param f The value of the independent variable
   * @return hyFloat The derivative correction
   */
  hyFloat DerivativeCorrection(
      long, hyFloat) const; // obtain the g'(x) for the chain rule
                            // differentiation under transformed variables
  /**
   * @brief Set the Ith Independent object
   *
   * @param i The index of the independent variable
   * @param f The value to set
   */
  void SetIthIndependent(long,
                         hyFloat); // set the value of i-th independent variable
  /**
   * @brief Check and set the Ith Independent object
   *
   * @param i The index of the independent variable
   * @param f The value to set
   * @return true if the value was set, false otherwise
   */
  bool CheckAndSetIthIndependent(
      long, hyFloat); // set the value of i-th independent variable
  /**
   * @brief Check if the Ith parameter is global
   *
   * @param i The index of the parameter
   * @return true if the parameter is global, false otherwise
   */
  bool IsIthParameterGlobal(long) const;

  /**
   * @brief Set all independent variables
   *
   * @param m The matrix of independent variables
   * @return long The number of independent variables set
   */
  long SetAllIndependent(_Matrix *);

  /**
   * @brief Update an independent variable
   *
   * @param i The index of the independent variable
   * @param b Whether to update the variable
   * @param sl1 The first simple list
   * @param sl2 The second simple list
   */
  void UpdateIndependent(long, bool, _SimpleList * = nil, _SimpleList * = nil);
  /**
   * @brief Update a dependent variable
   *
   * @param i The index of the dependent variable
   */
  void UpdateDependent(long);
  /**
   * @brief Update a dependent variable
   *
   * @param avl The AVL list of dependent variables to update
   */
  void UpdateDependent(_AVLList const &);

  /**
   * @brief Pre-compute the likelihood function
   *
   * @return true if the pre-computation was successful, false otherwise
   */
  bool PreCompute(void);
  /**
   * @brief Post-compute the likelihood function
   *
   */
  void PostCompute(void);
  /**
   * @brief Compute the likelihood function
   *
   * @return hyFloat The value of the likelihood function
   */
  virtual hyFloat Compute(void);

  /**
   * @brief Prepare the likelihood function for computation
   *
   * @param b Whether to prepare the function
   */
  void PrepareToCompute(bool = false);
  /**
   * @brief Finish the computation of the likelihood function
   *
   * @param b Whether to finish the computation
   */
  void DoneComputing(bool = false);
  /**
   * @brief Optimize the likelihood function
   *
   * @param options The optimization options
   * @return _Matrix* The result of the optimization
   */
  virtual _Matrix *Optimize(_AssociativeList const *options = nil);
  /**
   * @brief Construct a category matrix
   *
   * @param sl The simple list of categories
   * @param u The unsigned integer
   * @param b Whether to construct the matrix
   * @param s The string
   * @return _Matrix* The category matrix
   */
  _Matrix *ConstructCategoryMatrix(const _SimpleList &, unsigned, bool = true,
                                   _String * = nil);
  /**
   * @brief Get the name of the likelihood function
   *
   * @return const _String The name of the likelihood function
   */
  const _String GetMyName(void) const;

  /**
   * @brief The Simplex method for optimization
   *
   * @param precision The precision of the optimization
   * @param max_iterations The maximum number of iterations
   * @param max_evals The maximum number of evaluations
   * @return hyFloat The result of the optimization
   */
  hyFloat SimplexMethod(hyFloat &precision,
                        unsigned long max_iterations = 100000UL,
                        unsigned long max_evals = 0xFFFFFF);
  /**
   * @brief Anneal the likelihood function
   *
   * @param precision The precision of the annealing
   */
  void Anneal(hyFloat &precision);

  /**
   * @brief Simulate a dataset
   *
   * @param ds The dataset to simulate from
   * @param l The list of parameters
   * @param m1 The first matrix
   * @param m2 The second matrix
   * @param m3 The third matrix
   * @param s The string
   */
  void Simulate(_DataSet &, _List &, _Matrix * = nil, _Matrix * = nil,
                _Matrix * = nil, _String const * = nil) const;

  /**
   * @brief Reconstruct ancestral sequences
   *
   * @param ds The dataset
   * @param sl The simple list of nodes
   * @param s The string
   * @param b1 The first boolean
   * @param b2 The second boolean
   * @param b3 The third boolean
   */
  void ReconstructAncestors(_DataSet &, _SimpleList &, _String &, bool = false,
                            bool = false, bool = false);
  // 20090224: added an argument to allow the marginal state reconstruction
  // 20091009: added an argument to allow the reconstruction of leaves

  /**
   * @brief Get the maximum dimension of the likelihood function
   *
   * @return long The maximum dimension
   */
  long MaximumDimension(void);

  /**
   * @brief Get the covariance matrix
   *
   * @param sl The simple list of parameters
   * @return HBLObjectRef The covariance matrix
   */
  virtual HBLObjectRef CovarianceMatrix(_SimpleList * = nil);

  // compute  covariance matrix  based on the Hessian
  // optional list of parameters to estimate the conditional covariance for

  /**
   * @brief Rescan all variables in the likelihood function
   *
   * @param obtain_variable_mapping Whether to obtain the variable mapping
   */
  virtual void RescanAllVariables(bool obtain_variable_mapping = false);

  /**
   * @brief Check if the likelihood function depends on a given tree
   *
   * @param s The name of the tree
   * @return long The index of the tree if it is a dependency, -1 otherwise
   */
  long DependOnTree(_String const &) const;
  /**
   * @brief Check if the likelihood function depends on a given model
   *
   * @param s The name of the model
   * @return long The index of the model if it is a dependency, -1 otherwise
   */
  long DependOnModel(_String const &) const;
  /**
   * @brief Check if the likelihood function depends on a given dataset
   *
   * @param l The index of the dataset
   * @return long The index of the dataset if it is a dependency, -1 otherwise
   */
  long DependOnDS(long) const;
  /**
   * @brief Check if the likelihood function depends on a given data filter
   *
   * @param ID The ID of the data filter
   * @return true if the likelihood function depends on the data filter, false otherwise
   */
  bool DependOnDF(long ID) const { return theDataFilters.Find(ID) >= 0; }
  /**
   * @brief Map the tree tips to the data
   *
   * @param l The index of the tree
   * @param errorString The error string to store any errors in
   * @param leafScan Whether to scan the leaves
   * @return true if the mapping was successful, false otherwise
   */
  bool MapTreeTipsToData(long, _String *errorString, bool leafScan = false);
  /**
   * @brief Void the old results
   */
  void VoidOldResults(void) { computationalResults.ZeroUsed(); }
  /**
   * @brief Find a category variable
   *
   * @param l The index of the category variable
   * @return _CategoryVariable* The category variable
   */
  _CategoryVariable *FindCategoryVar(long);
  // return the category variable for a given partition
  /**
   * @brief Rank the variables
   *
   * @param tagger The tagger to use
   */
  void RankVariables(_AVLListX *tagger = nil);

  /**
   * @brief Get the Ith Tree object
   *
   * @param i The index of the tree
   * @return _TheTree* The tree
   */
  _TheTree *GetIthTree(long) const;
  /**
   * @brief Get the Ith Filter object
   *
   * @param i The index of the filter
   * @return _DataSetFilter const* The filter
   */
  _DataSetFilter const *GetIthFilter(long) const;
  /**
   * @brief Get the Ith Filter Mutable object
   *
   * @param i The index of the filter
   * @return _DataSetFilter* The filter
   */
  _DataSetFilter *GetIthFilterMutable(long) const;
  /**
   * @brief Get the Ith Filter Name object
   *
   * @param i The index of the filter
   * @return _String const* The name of the filter
   */
  _String const *GetIthFilterName(long) const;

  /**
   * @brief Get the Ith Frequencies object
   *
   * @param i The index of the frequencies
   * @return _Matrix* The frequencies
   */
  _Matrix *GetIthFrequencies(long) const;
  /**
   * @brief Get the Ith Frequencies Name object
   *
   * @param i The index of the frequencies
   * @return _String const* The name of the frequencies
   */
  _String const *GetIthFrequenciesName(long) const;

  /**
   * @brief Fill in the conditional probabilities
   *
   * @param l The index of the partition
   */
  void FillInConditionals(long = -1);

  /**
   * @brief Setup the likelihood function
   *
   * @param check_reversibility Whether to check for reversibility
   */
  void Setup(bool check_reversibility = true);

  /**
   * @brief Get the number of sequences
   *
   * @param l The index of the partition
   * @return long The number of sequences
   */
  long SequenceCount(long);
  /**
   * @brief Get the number of sites
   *
   * @return unsigned long The number of sites
   */
  unsigned long SiteCount(void) const;
  /**
   * @brief Rebuild the likelihood function
   *
   * @param rescan_parameters Whether to rescan the parameters
   */
  void Rebuild(bool rescan_parameters = false);
  /**
   * @brief Serialize the likelihood function
   *
   * @param sb The string buffer to serialize to
   * @param c The character to use
   * @param sl1 The first simple list
   * @param sl2 The second simple list
   */
  virtual void SerializeLF(_StringBuffer &, char = 0, _SimpleList * = nil,
                           _SimpleList * = nil);
  /**
   * @brief Check if the likelihood function has a computing template
   *
   * @return _Formula* The computing template
   */
  _Formula *HasComputingTemplate(void) const { return computingTemplate; }
  /**
   * @brief Compute the likelihood function using MPI
   *
   * @param l The index of the partition
   * @param b Whether to use MPI
   */
  void MPI_LF_Compute(long, bool = false);

#if defined _OPENMP
  /**
   * @brief Set the Thread Count object
   *
   * @param tc The thread count
   * @param fill_in Whether to fill in the conditionals
   */
  void SetThreadCount(long tc, bool fill_in = true) {
    if (tc != lfThreadCount) {
      lfThreadCount = tc;
      if (fill_in) {
        FillInConditionals();
      }
    }
  }
  /**
   * @brief Get the Thread Count object
   *
   * @return long The thread count
   */
  long GetThreadCount(void) const { return lfThreadCount; }
#else
  long GetThreadCount(void) const { return 1; }
  void SetThreadCount(long, bool = false) {};
#endif

  /**
   * @brief Compute the dependency lists
   *
   * @param receptacle The list to store the dependency lists in
   * @param max_dep The maximum number of dependencies to compute
   */
  void ComputeDependencyLists(_List &receptacle,
                              long max_dep = 0x7fffffff) const;
  /**
      This function computes a mapping from the index of each independent
     variable to the list of all dependent variables that it affects. index =>
     {n1,n2,...n3} special cases: index => null (no dependencies) index => {}
     (all variables, or >max_dep)

      @param max_dep the maximum number of dependent variables affected by a
     paritcular independent to generate an explicit list
  */

  /**
   * @brief Process the partition list
   *
   * @param sl The simple list of partitions
   * @param m The matrix of partitions
   * @return true if the processing was successful, false otherwise
   */
  bool ProcessPartitionList(_SimpleList &, _Matrix *) const;
  // given a matrix argument (argument 2; can be nil to include all)
  // populate a sorted list (argument 1)
  // of partitions indexed in the matrix (e.g. {{1,3}} would include the 2-nd
  // and 4-th partitions the third argument is the ID of the calling function
  // for error reporting returns true if at least one of the requested
  // partitions is valid otherwise returns false

  /**
   * @brief Get the total number of rate classes for a partition
   *
   * @param l The index of the partition
   * @param c The character flag
   * @return long The total number of rate classes
   */
  long TotalRateClassesForAPartition(long, char = 0);
  // given a partition index (assuming that category caches have been set up)
  // returns how many total rate classes there are for this partition
  // if the index is < 0, return the total number of categories for the LF as a
  // whole

  // char flag: 0 - all categories
  //            1 - only HMM
  //            2 - only constant on partition

  /**
   * @brief Check if the likelihood function has been set up
   *
   * @return true if the likelihood function has been set up, false otherwise
   */
  bool HasBeenSetup(void) { return hasBeenSetUp > 0; }

  /**
     @return  an associative array with a variety of metrics (see GetString on
     HyPhy wiki) about this LF.
     @author  SLKP
     @version 20110608
  */

  _AssociativeList *CollectLFAttributes(void) const;
  /**
   * @brief Unregister the listeners
   */
  void UnregisterListeners(void);
  /**
   * @brief Determine the local update policy
   */
  void DetermineLocalUpdatePolicy(void);
  /**
   * @brief Flush the local update policy
   */
  void FlushLocalUpdatePolicy(void);
  /**
   * @brief Terminate and dump the likelihood function
   *
   * @param error The error message
   * @param sig_term Whether the termination was due to a signal
   */
  void _TerminateAndDump(const _String &error, bool sig_term = false);

protected:
  void AllocateTemplateCaches(void);
  bool CheckIthPartition(unsigned long index, _String *errorString,
                         _String const * = nil, _String const * = nil,
                         _String const * = nil);

  _Matrix *PairwiseDistances(long index);
  void CheckDependentBounds(void);
  void AllocateSiteResults(void);
  void ZeroSiteResults(void);
  // 20090211: A utility function to reset site results.
  void GetInitialValues(void) const;
  bool checkPermissibility(_Matrix &m, long row);
  hyFloat computeAtAPoint(_Matrix &m, long row = 0);
  hyFloat replaceAPoint(_Matrix &m, long row, _Matrix &p, hyFloat &nV,
                        _Matrix &fv);

  void ScanAllVariablesOnPartition(_SimpleList &, _SimpleList &, _SimpleList &,
                                   _SimpleList &, bool = false);
  // internal function to scan all the on a set of partitions variables in
  // 20091103: SLKP
  // the boolean flag (if true), requests that only tree-associated variables
  // (no templates or frequencies) be scanned

  virtual void ScanAllVariables(void);
  // internal function to scan all the variables in

  void OptimalOrder(long, _SimpleList &, const _SimpleList *clone = nil);
  // determine the optimal order of compuation for a block

  hyFloat ComputeBlock(long, hyFloat *siteResults = nil,
                       long currentRateClass = -1, long = -1,
                       _SimpleList * = nil);
  // 20090224: SLKP
  // added the option to pass an interior branch (referenced by the 3rd argument
  // in the same order as flatTree) and a set of values for each site pattern
  // (indexed left to right) in the 4th argument

  // internal Setup, called from within constructors

  //      void            DumpingOrder            (long index, _SimpleList& sl);
  // used for FreeUpMemory

  void Cleanup(void);
  // internal cleanup, called from within destructors

  void Clear(void);

  long Bracket(long, hyFloat &, hyFloat &, hyFloat &, hyFloat &, hyFloat &,
               hyFloat &, hyFloat &, bool, _Matrix * = nil);
  // long          GradientBracketOneVar (_Matrix&, _Matrix& , _Matrix& ,
  // _Matrix&,  hyFloat& ,
  //                                       hyFloat&, hyFloat&, hyFloat&, bool
  //                                       retry = false);
  void LocateTheBump(long, hyFloat, hyFloat &, hyFloat &, bool, hyFloat = -1.);
  hyFloat GradientLocateTheBump(hyFloat, hyFloat &, _Matrix &, _Matrix &);
  void GradientDescent(hyFloat &, _Matrix &);
  hyFloat ConjugateGradientDescent(hyFloat, _Matrix &, bool localOnly = false,
                                   unsigned long = 0x7fffffff,
                                   _SimpleList *only_these_parameters = nil,
                                   hyFloat check_lf = -INFINITY,
                                   hyFloat min_improvement_to_continue = 0.);

  hyFloat SetParametersAndCompute(long, hyFloat, _Matrix * = nil,
                                  _Matrix * = nil, bool skip_compute = false);

  long CostOfPath(_DataSetFilter const *, _TheTree const *, _SimpleList &,
                  _SimpleList * = nil) const;

  void BuildLeafProbs(node<long> &, unsigned long *, unsigned long, _DataSet &,
                      _TheTree *, unsigned long &, bool, long,
                      _DataSetFilter const *, long, _DataSet * = nil) const;
  bool SingleBuildLeafProbs(node<long> &, long, _SimpleList &, _SimpleList &,
                            _TheTree *, bool, _DataSetFilter const *,
                            _SimpleList * = nil) const;

  bool HasBlockChanged(long) const;
  unsigned long BlockLength(long) const;
  void PartitionCatVars(_SimpleList &, long);
  // 20090210: extract variable indices for category variables in i-th partition
  // and append them to _SimpleList

  static void RandomizeList(_SimpleList &, long);
  static void CheckFibonacci(hyFloat);

  unsigned long PartitionLengths(char = 0, _SimpleList const * = nil) const;
  /*
      SLKP: 20090317

      argument = 0
          return the length (sites) of the longest part in this
          likelihood function
      argument = 1
          return the cumulative length (sites) of all parts in this
          likelihood function

      SLKP: 20090608

          argument 2 provides an optional subcollection of partititons to
     operate on; the default is to operate on all
  */

  void ComputeParameterPenalty(void);

  bool SendOffToMPI(long, void * = nil);
  void InitMPIOptimizer(void);
  void CleanupMPIOptimizer(void);
  void ComputeBlockInt1(long, hyFloat &, _TheTree *, _DataSetFilter *, char);
  void CheckStep(hyFloat &, const _Matrix &, _Matrix *selection = nil);
  void GetGradientStepBound(_Matrix &, hyFloat &, hyFloat &, long * = nil);
  void ComputeGradient(_Matrix &, hyFloat &, _Matrix &, _SimpleList &, long,
                       bool normalize = true);
  bool SniffAround(_Matrix &, hyFloat &, hyFloat &);
  void RecurseCategory(long, long, long, long, hyFloat, _SimpleList * = nil,
                       char = 0, hyFloat * = nil, long = -1,
                       _SimpleList * = nil);
  void RecurseConstantOnPartition(long, long, long, long, hyFloat, _Matrix &);

  static void SetNthBit(long &, unsigned char);
  static bool CheckNthBit(long &, unsigned char);
  void BuildIncrements(long, _SimpleList &);
  static unsigned char HighestBit(long);
  static unsigned char LowestBit(long);
  long HasHiddenMarkov(long, bool hmm = true) const;
  _Matrix *RemapMatrix(_Matrix *, const _SimpleList &) const;
  /* given a matrix where each column corresponds to a site pattern
   and a list of partitions that are covered
   this function returns a new matrix where each pattern is mapped
   to the corresponding sites in the alignment
  */
  void CleanUpOptimize(void);
  void ComputeBlockForTemplate(long, bool = false);
  void ComputeBlockForTemplate2(long, hyFloat *, hyFloat *, long);
  void DeleteCaches(bool = true);
  void PopulateConditionalProbabilities(long index, char runMode,
                                        hyFloat *buffer, _SimpleList &scalers,
                                        long = -1, _SimpleList * = nil);
  void ComputeSiteLikelihoodsForABlock(long, hyFloat *, _SimpleList &,
                                       long = -1, _SimpleList * = nil,
                                       char = 0);

  // this function computes a list of site probabilities for the i-th block (1st
  // parameter) stores them in pattern (left to right) order (2nd argument) and
  // writes scaling factors for each site into the third argument arguments 4
  // (branch index) and 5 (pattern assignments) allows the calculation of the
  // probability vector while setting a specific interior branch to a given
  // sequence

  hyFloat SumUpHiddenMarkov(const hyFloat *, _Matrix &, _Matrix &,
                            _SimpleList const *, const _SimpleList *, long);
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

      compute the log likelihood of the partition using the forward HMM
     algorithm with scaling
   */

  hyFloat SumUpConstantOnPartition(const hyFloat *, _Matrix const &,
                                   _SimpleList const &,
                                   _SimpleList const &) const;
  /*
      SLKP 20230113

      given conditional probs for each pattern        (1)
      the matrix of category weights                  (2)
      the vector of pattern frequencies               (3)
      the vector of by-site scaling factors as        (4)

      compute the log likelihood of the partition which is : log (sum [cat]
     (partition likelihood | cat) * Prob (cat))
   */

  void RunViterbi(_Matrix &, const hyFloat *, _Matrix &, _Matrix &,
                  _SimpleList const *, const _SimpleList *, long);
  /* Viterbi decoding for HMM; parameter meanings as in SumUpHiddenMarkov,
     except the first, which will store the optimal path to be returned */

  hyFloat SumUpSiteLikelihoods(long, const hyFloat *, const _SimpleList &);
  /*
   SLKP 20090318

   given a partition index (argument 1),
   a vector of _pattern_ likelihoods (argument 2) and
   a list of pattern scaling factors (argument 3)

   compute the log likelihood of the partition

   */

  /** optimization logger functions **/

  void LoggerLogL(hyFloat logL);
  void LoggerAddGradientPhase(hyFloat precision, hyFloat beta,
                              hyFloat scalar_product);
  void LoggerAddCoordinatewisePhase(hyFloat shrinkage, char convergence_mode);
  void LoggerAllVariables();
  void LoggerSingleVariable(unsigned long index, hyFloat logL,
                            hyFloat bracket_precision, hyFloat brent_precision,
                            hyFloat bracket_width, hyFloat movement,
                            unsigned long bracket_evals,
                            unsigned long brent_evals, unsigned long exp_count);

  void UpdateBlockResult(unsigned long, hyFloat);
  /*
      SLKP 20090318

      update cached log-likelihood value for a given
      partition
  */

  void CompileConstraints(void);
  void UncompileConstraints(void);

  /**
      20240305   SLKP
          These two functions populate / depopulate 'compiled_constraints'

   */

  _List *RecoverAncestralSequencesMarginal(long, _Matrix &, _List const &,
                                           bool = false);
  void RestoreScalingFactors(long, long, long, long *, long *);
  void SetupLFCaches(void);
  void SetupCategoryCaches(void);
  bool HasPartitionChanged(long);
  void SetupParameterMapping(void);
  void CleanupParameterMapping(void);

  void ConvertDependenciesToSimpleFormulas();
  /**
      Traverse the list of dependent variables and generate simple formulas if
     they can be created.
  */

  _SimpleList theTrees, theDataFilters, theProbabilities, indexInd, indexDep,
      indexCat, *nonConstantDep, blockDependancies,
      parameterTransformationFunction;
  /* 20110718: SLKP this list holds the index of the parameter interval mapping
     function used during optimization */

  _Vector computationalResults;

  _List optimalOrders, leafSkips, categoryTraversalTemplate,
      /*SLKP: 20090225

           This list contains as many entries (themselves of type _List) as
        there are partitions The entry will be empty for a partition without
        category variables For a partition with N category variables the entry
        will contain a

           1). _List of references to category variables themselves in the order
        that they appear in the blockDependancies

           2). _SimpleList of category counts for each variable C_1, C_2, ...
        C_N + a final entry with C_1 * C_2 * ... * C_N -- the total number of
        rate categories for this partition

           3). _SimpleList of incremental loop offsets for each variable, e.g.
        if the _SimpleList in 2). is 2,3,4, then this _SimpleList is 12,4,1

        SLKP: 20100416

           4). A _SimpleList including _relative_ indexing for the buffer array
               for HMM and "constant on partition" variables similar to 3, but
               ignoring all other variables. The last entry is the total number
               of HMM (or similar) rate classes on a partition

           5). A _SimpleList storing the type of each cat variable in a
        partition (flag) the last element is the cumulative flag composed of bit
        toggles to indicate what type of category variables there are for a
        given partition: _hyphyCategoryNormal _hyphyCategoryHMM
                  _hyphyCategoryCOP
       */
      indVarsByPartition, depVarsByPartition, *variable_to_node_map;

  long evalsSinceLastSetup,
      // this is necessary to force internal caches to be filled in at least
      // once for trees where branches have different numbers of rate categories
      templateKind,

      /*
          _hyphyLFComputationalTemplateNone:
              no computational template: simply return the summed up
         log-likelihood for all partitions in the likelihood function

          _hyphyLFComputationalTemplateBySite:
              all partitons have the same number of sites
              the template is in terms of SITE_LIKELIHOOD
              meant to implement model/tree mixture constructs

          _hyphyLFComputationalTemplateByPartition
              the template is in terms of BLOCK_LIKELIHOOD
              meant to implement random effects on the level of partitions
              models

          <0                                          : (-1-variable index) if a
         hidden markov process on site likelihoods
          >_hyphyLFComputationalTemplateByPartition   : a
         _hyphyLFComputationalTemplateByPartition shifted index of a user
                                                        function that will
         assemble the likelihood from site conditionals, written into matrices
                                                        SITE_LIKELIHOODS_0,
         SITE_LIKELIHOODS_1, ..
      */
      hasBeenSetUp;

  _Matrix *siteResults, *bySiteResults, *parameterValuesAndRanges;

  /* 20110718 SLKP
      this Nx4 (N is the number of adjustable parameters) matrix
      stores variable bounds and values for each independent variable;
      this is used for internal normalization transforms during optimization
   */

  bool siteArrayPopulated;

  _Formula *computingTemplate;
  MSTCache *mstCache;

  hyFloat smoothingTerm, smoothingReduction, smoothingPenalty, errorTolerance;

  /*
     these variables store conditional likelihoods for every paritition
     internal node (post-order traversal)
     and leaves;

     -conditionalInternalNodeLikelihoodCaches: 3D vector

          1st coordinate - the node index (in post-order traversal)
          2nd coordinate - the site (unique pattern) index (left-to-right)
          3rd coordiante - i-th marginal (for the i-th character;
     0-filterCharDimension)

          stores the probability for the subtree below this node given that i-th
     character is present at the node

          for partitions with category variables this vector is copied enough
     times to store every combination of category values

      -conditionalTerminalNodeStateFlag: 2D vector
          1st coordinate - the leaf index (in post-order traversal)
          2nd coordinate - the site (unique pattern) index (left-to-right)

          stores a non-negative integer for a leaf with a simple (non-ambiguous)
     character stores a negative value to index (-value-1) into the vector of
     conditionalTerminalNodeLikelihoodCaches and read off filterCharDimension
     characters from there
  */

  hyFloat **conditionalInternalNodeLikelihoodCaches, **siteScalingFactors,
      **branchCaches;

  _List conditionalTerminalNodeLikelihoodCaches;
  long **conditionalTerminalNodeStateFlag;

  _SimpleList overallScalingFactors, overallScalingFactorsBackup,
      // this (and siteCorrectionsBackup)
      // is needed to store site/partition scaling factors when
      // treeTraversalMasks and branch caching are used, otherwise scaling done
      // by ComputeBranchCaches will fubar the scaling for
      // ComputeTreeBlockByBranch depending on site patterns
      canUseReversibleSpeedups,
      // a partition will be tagged with '1' if its tree has only
      // time-reversible models
      siteScalerBuffer,
      // used for LF with category variables to
      // store site-by-site scaling factors
      *_variables_changed_during_last_compute;

  _AVLList *variables_changed_during_last_compute;

  _List localUpdatePolicy, matricesToExponentiate, treeTraversalMasks,
      computedLocalUpdatePolicy, siteCorrections, siteCorrectionsBackup,
      cachedBranches,
      // for models with categories, a list of site by site scaling operation
      // counts
      partScalingCache,
      // used to store site by site scalers in computations that are performed
      // on a site-by-site basis; includes scratch cache for remapping
      gradientBlocks;

  _AssociativeList *optimizatonHistory;

  _CompiledMatrixData *compiled_constraints;
  /**
      SLKP: 20240305
      Use this object for storing compiled
  */

#ifdef _OPENMP
  long lfThreadCount;
#endif
};

//_______________________________________________________________________________________

/**
 * @brief A class for custom functions
 *
 */
class _CustomFunction : public _LikelihoodFunction {

public:
  /**
   * @brief Construct a new _CustomFunction object
   *
   * @param s The string to construct from
   * @param vc The variable container
   */
  _CustomFunction(const _String &, _VariableContainer const *context = nil);

  /**
   * @brief Compute the custom function
   *
   * @return hyFloat The value of the custom function
   */
  virtual hyFloat Compute(void);
  /**
   * @brief Rescan all variables in the custom function
   *
   * @param b Whether to rescan the variables
   */
  virtual void RescanAllVariables(bool = false) {}
  /**
   * @brief Serialize the custom function
   *
   * @param res The string buffer to serialize to
   * @param c The character to use
   * @param sl1 The first simple list
   * @param sl2 The second simple list
   */
  virtual void SerializeLF(_StringBuffer &res, char = 0, _SimpleList * = nil,
                           _SimpleList * = nil) {
    res.AppendNewInstance(
        (_String *)myBody.toStr(kFormulaStringConversionNormal));
  }

private:
  _Formula myBody;
};

//_______________________________________________________________________________________

extern bool forceRecomputation, isInOptimize;

extern long lockedLFID;
/*
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
dataFileDefaultWidth           ,
dataFileGapWidth               ,
categorySimulationMethod       ,
useInitialDistanceGuess            ,
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
reduceLFSmoothing              ;*/

/**
 * @brief Handle the result of a state counter
 *
 * @param f The formula
 * @param sl The simple list
 * @param l1 The first long
 * @param l2 The second long
 * @param l3 The third long
 * @param m1 The first matrix
 * @param m2 The second matrix
 */
void StateCounterResultHandler(_Formula &, _SimpleList *, long &, long &, long,
                               _Matrix &, _Matrix &);

/**
 * @brief Find a likelihood function by name
 *
 * @param s The name of the likelihood function
 * @return _LikelihoodFunction* The likelihood function
 */
_LikelihoodFunction *FindLikeFuncByName(_String &);

/**
 * @brief Do something for each likelihood function
 *
 * @tparam ACTION The action to perform
 * @param cb The callback function
 */
template <typename ACTION> void DoForEachLikelihoodFunction(ACTION cb) {
  for (unsigned long i = 0UL; i < likeFuncNamesList.lLength; i++) {
    if (((_String *)likeFuncNamesList.GetItem(i))->nonempty()) {
      cb((_LikelihoodFunction *)likeFuncList.GetItem(i), i);
    }
  }
}

extern bool usedCachedResults;

extern hyFloat _lfScalerUpwards, _lfScalingFactorThreshold, _logLFScaler;

extern _Vector _scalerMultipliers, _scalerDividers;

/**
 * @brief Acquire a scaler multiplier
 *
 * @param l The index of the scaler multiplier
 * @return hyFloat The scaler multiplier
 */
hyFloat acquireScalerMultiplier(long);
/**
 * @brief A custom log function
 *
 * @param f The value to take the log of
 * @return hyFloat The log of the value
 */
hyFloat myLog(hyFloat);
/**
 * @brief Add a scaler
 *
 * @param f The value to add
 * @param l1 The first long
 * @param l2 The second long
 * @return long The result of the addition
 */
long addScaler(hyFloat, long, long);
/**
 * @brief Map a parameter to an interval
 *
 * @param f The parameter value
 * @param l The interval
 * @param b Whether to map the parameter
 * @return hyFloat The mapped parameter
 */
hyFloat mapParameterToInverval(hyFloat, long const, bool);
/**
 * @brief Obtain the derivative correction
 *
 * @param f The value
 * @param l The interval
 * @return hyFloat The derivative correction
 */
hyFloat obtainDerivativeCorrection(hyFloat, long const);

#ifdef __HYPHYMPI__
extern _Matrix resTransferMatrix;
extern long hyphyMPIOptimizerMode;
extern _String mpiLoopSwitchToOptimize, mpiLoopSwitchToBGM;

extern _SimpleList mpiNodesThatCantSwitch;

/**
 * @brief Retrieve the MPI count
 *
 * @param c The character
 * @return long The MPI count
 */
long RetrieveMPICount(char);
/**
 * @brief Switch the MPI nodes to MPI mode
 *
 * @param l The MPI mode
 */
void MPISwitchNodesToMPIMode(long);

#endif

#endif
