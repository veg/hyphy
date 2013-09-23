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

#ifndef __THETREE__
#define __THETREE__

#include "treetopology.h"
#include "growingvector.h"
#include "datasetfilternumeric.h"

#define _HY2TREE(X) (dynamic_cast<_TheTree*>(X))


class _TheTree : public _TreeTopology {

  // theModel matrix of _TheTree contains the column matrix of probabilities,
  // which is computed
  // based on the DataSetFilter passed on to the tree at the initialization
  // stage

public:

  _TheTree();            // default constructor - doesn't do much
  _TheTree(_String name, _String &parms,
           bool = true); // builds a tree from a string
  _TheTree(_String name, _TreeTopology *); // builds a tree from a tree topology
  _TheTree(_String name, _TheTree *);      // builds a tree from another tree

  virtual ~_TheTree(void);
  virtual bool HasChanged(void);
  virtual void MarkDone(void);
  bool HasChanged2(void);

  _CalcNode *DepthWiseTraversal(bool = false);
  //performs a post-order traversal
  _CalcNode *DepthWiseTraversalRight(bool = false);
  //performs a post-order tree traversal going right first
  _CalcNode *DepthWiseTraversalLevel(long &, bool = false);
  //performs a post-order tree traversal
  //storing current node depth
  _CalcNode *StepWiseTraversal(bool = false);
  //performs a pre-order  tree traversal
  _CalcNode *StepWiseTraversalLevel(long &, bool = false);
  //performs a pre-order wise tree traversal
  //storing current node depth

  _CalcNode *LeafWiseTraversal(bool = false);
  //iterate through the leaves (left-to-right)

  virtual bool FinalizeNode(node<long> *, long, _String &, _String &, _String &,
                            _String * = NULL);
  virtual BaseRef makeDynamic(void);

  virtual BaseRef makeDynamicCopy(_String *);
  node<long> *DuplicateTreeStructure(node<long> *, _String *, bool);
  virtual BaseRef toStr(void);
  virtual unsigned long ObjectClass(void) { return TREE; }

  virtual _PMathObj
  Execute(long, _PMathObj = nil, _PMathObj = nil,
          _hyExecutionContext *context = _hyDefaultExecutionContext,
          _PMathObj = nil);
  virtual _PMathObj TEXTreeString(_PMathObj);
  virtual _PMathObj PlainTreeString(_PMathObj, _PMathObj);

  virtual void GetNodeName(node<long> *, _String &, bool = false);
  virtual void GetBranchLength(node<long> *, _String &, bool = false);
  virtual void GetBranchLength(node<long> *, _Parameter &);
  virtual void GetBranchValue(node<long> *, _String &);
  virtual _String *GetBranchSpec(node<long> *);
  virtual void GetBranchVarValue(node<long> *, _String &, long);

  void InitializeTreeFrequencies(_Matrix *, bool = false);

  _Parameter ReleafTreeAndCheck(_DataSetFilter *, long, bool,
                                long categID = -1);
  _Parameter ReleafTreeAndCheckChar4(_DataSetFilter *, long, bool,
                                     long categID = -1);

  _Parameter ReleafTree(_DataSetFilter *, long, long, long, long);
  _Parameter ReleafTreeDegenerate(_DataSetFilter *, long);

  _Parameter ReleafTreeCache(_DataSetFilter *, long, long, long, long, long);
#if USE_SCALING_TO_FIX_UNDERFLOW
  _Parameter ThreadReleafTreeCache(_DataSetFilter *, long, long, long, long,
                                   long, long offset = 0, long fixAttempt = 0,
                                   _Parameter = 690.);
  _Parameter doScaling(_DataSetFilter *, long, long, long, long, _Parameter,
                       bool, bool);
#else
  _Parameter ThreadReleafTreeCache(_DataSetFilter *, long, long, long, long,
                                   long, long offset = 0);
#endif

  _Parameter Process3TaxonNumericFilter(_DataSetFilterNumeric *, long = 0);

  void ThreadMatrixUpdate(long, bool);
  void SerialMatrixUpdate(long, bool);
  void MatrixCacheUpdate(void);

  _Parameter ThreadReleafTreeCharCache(_DataSetFilter *, long, long, long, long,
                                       long, long offset = 0);
  _Parameter ReleafTreeCharDegenerate(_DataSetFilter *, long);
  _Parameter ReleafTreeChar4(_DataSetFilter *, long, long, long, long, long);
  _Parameter ReleafTreeChar4Degenerate(_DataSetFilter *, long);

  _Parameter ThreadReleafTreeChar4(_DataSetFilter *, long, long, long, long,
                                   long, long offset = 0);
  _Parameter ReleafTreeChar4(_DataSetFilter *, long, long, long, long);

  _Parameter Probij(long, long, _CalcNode *);
  _Parameter PruneTree(long categID = -1);
  _Parameter PruneTreeChar(long categID = -1);
  _Parameter PruneTreeCharCache(long categID = -1);
  _Parameter PruneTreeChar4(long categID = -1);
  _Parameter PruneTreeChar4Cache(long categID = -1);

  _List *RecoverAncestralSequences(_DataSetFilter *, _SimpleList &, _List &,
                                   _Parameter *, _Parameter *, long, long *,
                                   _GrowingVector *, bool = false);
  void RecoverNodeSupportStates(_DataSetFilter *, long, long, _Matrix &);
  void RecoverNodeSupportStates2(node<long> *, _Parameter *, _Parameter *,
                                 long);
  _List *SampleAncestors(_DataSetFilter *, node<long> *);
  void PurgeTree(void);

  long ComputeReleafingCost(_DataSetFilter *, long, long, _SimpleList * = nil,
                            long = 0);
  long ComputeReleafingCostChar(_DataSetFilter *, long, long);
  void DumpingOrder(_DataSetFilter *, _SimpleList &);
  void SetTreeCodeBase(long);
  long IsLinkedToALF(long &);

  bool HasCache(void) { return topLevelNodes.lLength > 0; }

  long GetLeafCount(void) { return flatLeaves.lLength; }

  long GetINodeCount(void) { return flatNodes.lLength; }

  void ScanAndAttachVariables(void);
  void ScanForVariables(_AVLList &l, _AVLList &l2, _AVLListX *tagger = nil,
                        long weight = 0);
  void ScanForDVariables(_AVLList &l, _AVLList &l2);
  void ScanForGVariables(_AVLList &, _AVLList &, _AVLListX *tagger = nil,
                         long weight = 0);
  void ScanForCVariables(_AVLList &);
  void MolecularClock(_String &, _List &);

  void SetUp(void);
  void SetUpMatrices(long);
  void CleanUpMatrices(void);
  void BuildTopLevelCache(void);
  void KillTopLevelCache(void);
  void SetCompMatrices(long);

  virtual void ClearConstraints(void);

  bool FindScalingVariables(_SimpleList &);
  bool HaveStringBranchLengths(void);
  void AssignLabelsToBranches(node<nodeCoord> *, _String *, bool);

  node<nodeCoord> *AlignedTipsMapping(bool first = false,
                                      bool respectRoot = true);

  void AlignNodes(node<nodeCoord> *);

  node<nodeCoord> *ScaledBranchMapping(node<nodeCoord> *, _String *, long,
                                       long &, char);

  node<nodeCoord> *RadialBranchMapping(node<long> *, node<nodeCoord> *,
                                       _String *, _Parameter, long &,
                                       _Parameter &, char);

  void ScaledBranchReMapping(node<nodeCoord> *, _Parameter);
  char &RootedFlag(void) { return rooted; }

  nodeCoord TreeTEXRecurse(node<nodeCoord> *, _String &, _Parameter, _Parameter,
                           long, long);
  void TreePSRecurse(node<nodeCoord> *, _String &, _Parameter, _Parameter, long,
                     long, long, long, _AssociativeList * = nil, char = 0,
                     _Parameter * = nil);

  bool AllBranchesHaveModels(long);
  void ScanSubtreeVars(_List &, char, _CalcNode *);
  void BuildINodeDependancies(void);
  void AllocateResultsCache(long);
  long CountTreeCategories(void);
  void CompileListOfModels(_SimpleList &);

  void MarkMatches(_DataSetFilter *, long, long);
  long GetLowerBoundOnCost(_DataSetFilter *);
  long GetLowerBoundOnCostWithOrder(_DataSetFilter *, _SimpleList *);
  _SimpleList &GetLeftINodes(void) { return leftiNodes; }
  bool MatchLeavesToDF(_SimpleList &, _DataSetFilter *, bool);
  virtual void RemoveModel(void);
  _String *TreeUserParams(void);

  _String CompareSubTrees(_TheTree *, node<long> *);
  _String FindMaxCommonSubTree(_TheTree *, long &, _List *);
  void WeightedCharacterDifferences(_Parameter, _Matrix *, _Matrix *,
                                    long = -1);
  void AddNodeNamesToDS(_DataSet *, bool, bool, char);
  // if the
  _Parameter PSStringWidth(_String &);

  _Parameter DetermineBranchLengthGivenScalingParameter(long, _String &, char);

  _AVLListX *ConstructNodeToIndexMap(bool);
  // 20090206: SLKP
  // makes an AVL of with keys storing memory addresses of node<long> tree nodes
  // and values showing the order in either flatLeaves (bool = false) or
  // flatNodes (bool = true)

  void MapPostOrderToInOderTraversal(_SimpleList &, bool = true);
  // 20090306: SLKP
  // 20100511: SLKP
  // construct a post-order -> in-order traveral map for internal nodes
  // bool = true (internal nodes), bool = false (leaf nodes)

  void AddBranchToForcedRecomputeList(long idx) {
    forceRecalculationOnTheseBranches << idx;
  }
  void ClearForcedRecomputeList(void) {
    forceRecalculationOnTheseBranches.Clear();
  }
  bool HasForcedRecomputeList(void) {
    return forceRecalculationOnTheseBranches.lLength;
  }
  // 20090306: SLKP
  // the previous two functions are used to manipulate the list of
  // branches that are marked as 'dirty' for LF computation purposes
  // because of external manipulations of the cache (e.g. computing the LF with
  // one of the interior
  // nodes set to

  void SetupCategoryMapsForNodes(_List &, _SimpleList &, _SimpleList &);
  /* 20090325: SLKP
     a wrapper function to set up category variable maps for nodes;
     see comments for remapMyCategories in _CalcNode
   */

  _Parameter VerySimpleLikelihoodEvaluator(_SimpleList &updateNodes,
                                           _DataSetFilter *theFilter,
                                           _Parameter *iNodeCache,
                                           long *lNodeFlags,
                                           _GrowingVector *lNodeResolutions);

#ifdef MDSOCL
  _Parameter OCLLikelihoodEvaluator(_SimpleList &updateNodes,
                                    _DataSetFilter *theFilter,
                                    _Parameter *iNodeCache, long *lNodeFlags,
                                    _GrowingVector *lNodeResolutions,
                                    _OCLEvaluator &OCLEval);
#endif

#ifdef _SLKP_LFENGINE_REWRITE_
  void SampleAncestorsBySequence(_DataSetFilter *, _SimpleList &, node<long> *,
                                 _AVLListX *, _Parameter *, _List &,
                                 _SimpleList *, _List &, _Parameter *, long);

  _Parameter ComputeTreeBlockByBranch(_SimpleList &, _SimpleList &,
                                      _SimpleList *, _DataSetFilter *,
                                      _Parameter *, long *, _Parameter *,
                                      _GrowingVector *, long &, long, long,
                                      long = -1, _Parameter * = nil,
                                      long * = nil, long = -1, long * = nil);
  long DetermineNodesForUpdate(_SimpleList &, _List * = nil, long = -1,
                               long = -1, bool = true);
  void ExponentiateMatrices(_List &, long, long = -1);
  void FillInConditionals(_DataSetFilter *, _Parameter *, _SimpleList *);

  void ComputeBranchCache(_SimpleList &, long nodeID, _Parameter *cache,
                          _Parameter *iNodeCache, _DataSetFilter *theFilter,
                          long *lNodeFlags, _Parameter *scalingAdjustments,
                          long *siteCorrectionCounts,
                          _GrowingVector *lNodeResolutions, long &overallScaler,
                          long siteFrom, long siteTo, long catID,
                          _SimpleList * = nil, _Parameter * = nil);

  _Parameter ComputeLLWithBranchCache(_SimpleList &siteOrdering, long brID,
                                      _Parameter *cache,
                                      _DataSetFilter *theFilter, long siteFrom,
                                      long siteTo, long catID,
                                      _Parameter *storageVec = nil);

  _Parameter ComputeTwoSequenceLikelihood(
      _SimpleList &siteOrdering, _DataSetFilter *theFilter, long *lNodeFlags,
      _GrowingVector *lNodeResolutions, long siteFrom, long siteTo, long catID,
      _Parameter *storageVec = nil);
#endif

  // --------------------------

  long *nodeStates;
  char *nodeMarkers;

  _Parameter *rootIChildrenCache, *marginalLikelihoodCache;

  _AVLListXL *aCache;

  long categoryCount;

protected:

  bool IntPopulateLeaves(_DataSetFilter *, long, long);

  _Parameter ConditionalBranchLikelihood(node<long> *, node<long> *,
                                         _Parameter *, _Parameter *, long,
                                         long);

  _Parameter ConditionalNodeLikelihood(node<long> *, node<long> *, _Parameter *,
                                       _Parameter *, long, long);

  _List *MapCBaseToCharacters(_DataSetFilter *, bool = true);

  virtual void PreTreeConstructor(bool);
  virtual void PostTreeConstructor(bool);

  // all of the following members exist to speed-up the pruning algorithm
  // the are created by calling the function set up
  _SimpleList flatLeaves, flatNodes, leftiNodes, flatParents, topLevelNodes,
      topLevelLeftL, topLevelRightL, forceRecalculationOnTheseBranches,
      nodesToUpdate;

};

#endif
