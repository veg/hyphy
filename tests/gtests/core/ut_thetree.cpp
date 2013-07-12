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

#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//Generate necessary includes from the respective implementation file
#include "hy_globals.h"
#include "thetree.h"
#include "calcnode_globals.h"
#include "growingvector.h"
#include "math.h"
#include "float.h"
#include "ctype.h"
#include "string.h"
#include "calcnode.h"
#include "scfg.h"
#include "parser.h"
#include "category.h"
#include "batchlan.h"
#include "likefunc.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__TheTreeTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__TheTreeTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _Parametertest = new _Parameter();
    _CalcNodetest = new _CalcNode();
    chartest = new char();
    _Listtest = new _List();
    _DataSetFilterNumerictest = new _DataSetFilterNumeric();
    //node<long>test = new node<long>();
    _TheTreetest = new _TheTree();
    _DataSettest = new _DataSet(FILEtest);
    _hyExecutionContexttest = new _hyExecutionContext(_VariableContainertest, _Stringtest);
    _Matrixtest = new _Matrix();
    _AVLListXtest = new _AVLListX(_SimpleListtest);
    _SimpleListtest = new _SimpleList();
    _DataSetFiltertest = new _DataSetFilter(_DataSettest, *chartest, *_Stringtest);
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
    //node<nodeCoord>test = new node<nodeCoord>();
    _PMathObjtest = new _PMathObj();
    _VariableContainertest = new _VariableContainer();
  }

  virtual ~DISABLED__TheTreeTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
    delete longtest;
    delete _Stringtest;
    delete _Parametertest;
    delete _CalcNodetest;
    delete chartest;
    delete _Listtest;
    delete _DataSetFilterNumerictest;
    //delete node<long>test;
    delete _TheTreetest;
    delete _DataSettest;
    delete _hyExecutionContexttest;
    delete _Matrixtest;
    delete _AVLListXtest;
    delete _SimpleListtest;
    delete _DataSetFiltertest;
    delete booltest;
    delete _AVLListtest;
    //delete node<nodeCoord>test;
    delete _PMathObjtest;
    delete _VariableContainertest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  _Parameter* _Parametertest;
  _CalcNode* _CalcNodetest;
  char* chartest;
  _List* _Listtest;
  _DataSetFilterNumeric* _DataSetFilterNumerictest;
  //node<long>* node<long>test;
  _TheTree* _TheTreetest;
  _DataSet* _DataSettest;
  _hyExecutionContext* _hyExecutionContexttest;
  _Matrix* _Matrixtest;
  _AVLListX* _AVLListXtest;
  _SimpleList* _SimpleListtest;
  _DataSetFilter* _DataSetFiltertest;
  bool* booltest;
  _AVLList* _AVLListtest;
  //node<nodeCoord>* node<nodeCoord>test;
  _PMathObj* _PMathObjtest;
  _VariableContainer* _VariableContainertest;
};


TEST_F(DISABLED__TheTreeTest, AddNodeNamesToDSTest) {

  //_TheTreetest->AddNodeNamesToDS(_DataSettest, *booltest, *booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, AlignNodesTest) {

  //_TheTreetest->AlignNodes(node<nodeCoord>test);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, AlignedTipsMappingTest) {

  //node<nodeCoord>* resultnode<nodeCoord> = _TheTreetest->AlignedTipsMapping(*booltest, *booltest);
  //EXPECT_EQ (resultnode<nodeCoord>*, 0);

}


TEST_F(DISABLED__TheTreeTest, AllBranchesHaveModelsTest) {

  bool resultbool = _TheTreetest->AllBranchesHaveModels(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, AllocateResultsCacheTest) {

  _TheTreetest->AllocateResultsCache(*longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, AllocateUnderflowScalersTest) {

  //_TheTreetest->AllocateUnderflowScalers(*longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, AssignLabelsToBranchesTest) {

  //_TheTreetest->AssignLabelsToBranches(node<nodeCoord>test);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, BuildINodeDependanciesTest) {

  _TheTreetest->BuildINodeDependancies();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, BuildTopLevelCacheTest) {

  _TheTreetest->BuildTopLevelCache();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, CleanUpMatricesTest) {

  _TheTreetest->CleanUpMatrices();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ClearConstraintsTest) {

  _TheTreetest->ClearConstraints();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, CompareSubTreesTest) {

  //_String result_String = _TheTreetest->CompareSubTrees(_TheTreetest, node<long>test);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(DISABLED__TheTreeTest, CompileListOfModelsTest) {

  _TheTreetest->CompileListOfModels(*_SimpleListtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ComputeBranchCacheTest) {

  //_TheTreetest->ComputeBranchCache();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ComputeLLWithBranchCacheTest) {

  //_Parameter result_Parameter = _TheTreetest->ComputeLLWithBranchCache(*_SimpleListtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ComputeReleafingCostTest) {

  //long resultlong = _TheTreetest->ComputeReleafingCost(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, ComputeReleafingCostCharTest) {

  //long resultlong = _TheTreetest->ComputeReleafingCostChar(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, ComputeTreeBlockByBranchTest) {

  //_Parameter result_Parameter = _TheTreetest->ComputeTreeBlockByBranch();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ComputeTwoSequenceLikelihoodTest) {

  //_Parameter result_Parameter = _TheTreetest->ComputeTwoSequenceLikelihood();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ConditionalBranchLikelihoodTest) {

  //_Parameter result_Parameter = _TheTreetest->ConditionalBranchLikelihood(node<long>test);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ConditionalNodeLikelihoodTest) {

  //_Parameter result_Parameter = _TheTreetest->ConditionalNodeLikelihood(node<long>test);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ConstructNodeToIndexMapTest) {

  _AVLListX* result_AVLListX = _TheTreetest->ConstructNodeToIndexMap(*booltest);
  //EXPECT_EQ (result_AVLListX*, 0);

}


TEST_F(DISABLED__TheTreeTest, CountTreeCategoriesTest) {

  long resultlong = _TheTreetest->CountTreeCategories();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, DeallocateUnderflowScalersTest) {

  //_TheTreetest->DeallocateUnderflowScalers();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, DepthWiseTraversalTest) {

  _CalcNode* result_CalcNode = _TheTreetest->DepthWiseTraversal(*booltest);
  //EXPECT_EQ (result_CalcNode*, 0);

}


TEST_F(DISABLED__TheTreeTest, DepthWiseTraversalLevelTest) {

  _CalcNode* result_CalcNode = _TheTreetest->DepthWiseTraversalLevel(*longtest, *booltest);
  //EXPECT_EQ (result_CalcNode*, 0);

}


TEST_F(DISABLED__TheTreeTest, DepthWiseTraversalRightTest) {

  _CalcNode* result_CalcNode = _TheTreetest->DepthWiseTraversalRight(*booltest);
  //EXPECT_EQ (result_CalcNode*, 0);

}


TEST_F(DISABLED__TheTreeTest, DetermineBranchLengthGivenScalingParameterTest) {

  //_Parameter result_Parameter = _TheTreetest->DetermineBranchLengthGivenScalingParameter();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, DetermineNodesForUpdateTest) {

  long resultlong = _TheTreetest->DetermineNodesForUpdate(*_SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, DumpingOrderTest) {

  _TheTreetest->DumpingOrder(_DataSetFiltertest, *_SimpleListtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, DuplicateTreeStructureTest) {

  //node<long>* resultnode<long> = _TheTreetest->DuplicateTreeStructure(node<long>test);
  //EXPECT_EQ (resultnode<long>*, 0);

}


TEST_F(DISABLED__TheTreeTest, ExecuteTest) {

  _PMathObj result_PMathObj = _TheTreetest->Execute(*longtest, *_PMathObjtest, *_PMathObjtest, _hyExecutionContexttest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__TheTreeTest, ExponentiateMatricesTest) {

  _TheTreetest->ExponentiateMatrices(*_Listtest, *longtest, *longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, FillInConditionalsTest) {

  //_TheTreetest->FillInConditionals(_DataSetFiltertest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, FinalizeNodeTest) {

  //bool resultbool = _TheTreetest->FinalizeNode(node<long>test, *longtest, *_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, FindMaxCommonSubTreeTest) {

  //_String result_String = _TheTreetest->FindMaxCommonSubTree(_TheTreetest, *longtest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(DISABLED__TheTreeTest, FindScalingVariablesTest) {

  bool resultbool = _TheTreetest->FindScalingVariables(*_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, GetBranchLengthTest) {

  //_TheTreetest->GetBranchLength(node<long>test, *_Parametertest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, GetBranchLength1Test) {

  //_TheTreetest->GetBranchLength(node<long>test, *_Stringtest, *booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, GetBranchSpecTest) {

  //_String* result_String = _TheTreetest->GetBranchSpec(node<long>test);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(DISABLED__TheTreeTest, GetBranchValueTest) {

  //_TheTreetest->GetBranchValue(node<long>test, *_Stringtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, GetBranchVarValueTest) {

  //_TheTreetest->GetBranchVarValue(node<long>test, *_Stringtest, *longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, GetLowerBoundOnCostTest) {

  long resultlong = _TheTreetest->GetLowerBoundOnCost(_DataSetFiltertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, GetLowerBoundOnCostWithOrderTest) {

  //long resultlong = _TheTreetest->GetLowerBoundOnCostWithOrder(_DataSetFiltertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, GetNodeNameTest) {

  //_TheTreetest->GetNodeName(node<long>test, *_Stringtest, *booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, HasChangedTest) {

  bool resultbool = _TheTreetest->HasChanged();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, HasChanged2Test) {

  bool resultbool = _TheTreetest->HasChanged2();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, HaveStringBranchLengthsTest) {

  bool resultbool = _TheTreetest->HaveStringBranchLengths();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, InitializeTreeFrequenciesTest) {

  _TheTreetest->InitializeTreeFrequencies(_Matrixtest, *booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, IntPopulateLeavesTest) {

  //bool resultbool = _TheTreetest->IntPopulateLeaves(_DataSetFiltertest, *longtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, IsLinkedToALFTest) {

  long resultlong = _TheTreetest->IsLinkedToALF(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__TheTreeTest, KillTopLevelCacheTest) {

  _TheTreetest->KillTopLevelCache();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, LeafWiseTraversalTest) {

  _CalcNode* result_CalcNode = _TheTreetest->LeafWiseTraversal(*booltest);
  //EXPECT_EQ (result_CalcNode*, 0);

}


TEST_F(DISABLED__TheTreeTest, MapCBaseToCharactersTest) {

  //_List* result_List = _TheTreetest->MapCBaseToCharacters(_DataSetFiltertest, *booltest);
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(DISABLED__TheTreeTest, MapPostOrderToInOderTraversalTest) {

  _TheTreetest->MapPostOrderToInOderTraversal(*_SimpleListtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, MarkDoneTest) {

  _TheTreetest->MarkDone();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, MarkMatchesTest) {

  //_TheTreetest->MarkMatches(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, MatchLeavesToDFTest) {

  //bool resultbool = _TheTreetest->MatchLeavesToDF(*_SimpleListtest, _DataSetFiltertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__TheTreeTest, MatrixCacheUpdateTest) {

  _TheTreetest->MatrixCacheUpdate();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, MolecularClockTest) {

  _TheTreetest->MolecularClock(*_Stringtest, *_Listtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, OCLLikelihoodEvaluatorTest) {

  //_Parameter result_Parameter = _TheTreetest->OCLLikelihoodEvaluator(*_SimpleListtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PSStringWidthTest) {

  _Parameter result_Parameter = _TheTreetest->PSStringWidth(*_Stringtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PlainTreeStringTest) {

  _PMathObj result_PMathObj = _TheTreetest->PlainTreeString(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__TheTreeTest, PostTreeConstructorTest) {

  //_TheTreetest->PostTreeConstructor(*booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, PreTreeConstructorTest) {

  //_TheTreetest->PreTreeConstructor(*booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ProbijTest) {

  _Parameter result_Parameter = _TheTreetest->Probij(*longtest, *longtest, _CalcNodetest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, Process3TaxonNumericFilterTest) {

  _Parameter result_Parameter = _TheTreetest->Process3TaxonNumericFilter(_DataSetFilterNumerictest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PruneTreeTest) {

  _Parameter result_Parameter = _TheTreetest->PruneTree(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PruneTreeCharTest) {

  _Parameter result_Parameter = _TheTreetest->PruneTreeChar(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PruneTreeChar4Test) {

  _Parameter result_Parameter = _TheTreetest->PruneTreeChar4(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PruneTreeChar4CacheTest) {

  _Parameter result_Parameter = _TheTreetest->PruneTreeChar4Cache(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, PurgeTreeTest) {

  _TheTreetest->PurgeTree();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, RadialBranchMappingTest) {

  //node<nodeCoord>* resultnode<nodeCoord> = _TheTreetest->RadialBranchMapping();
  //EXPECT_EQ (resultnode<nodeCoord>*, 0);

}


TEST_F(DISABLED__TheTreeTest, RecoverAncestralSequencesTest) {

  //_List* result_List = _TheTreetest->RecoverAncestralSequences();
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(DISABLED__TheTreeTest, RecoverNodeSupportStatesTest) {

  //_TheTreetest->RecoverNodeSupportStates(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, RecoverNodeSupportStates2Test) {

  //_TheTreetest->RecoverNodeSupportStates2(node<long>test);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeTest) {

  //_Parameter result_Parameter = _TheTreetest->ReleafTree(_DataSetFiltertest, *longtest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeAndCheckTest) {

  //_Parameter result_Parameter = _TheTreetest->ReleafTreeAndCheck(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeAndCheckChar4Test) {

  //_Parameter result_Parameter = _TheTreetest->ReleafTreeAndCheckChar4(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeCacheTest) {

  //_Parameter result_Parameter = _TheTreetest->ReleafTreeCache(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeChar4Test) {

  //_Parameter result_Parameter = _TheTreetest->ReleafTreeChar4(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeChar4DegenerateTest) {

  _Parameter result_Parameter = _TheTreetest->ReleafTreeChar4Degenerate(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeCharDegenerateTest) {

  _Parameter result_Parameter = _TheTreetest->ReleafTreeCharDegenerate(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ReleafTreeDegenerateTest) {

  _Parameter result_Parameter = _TheTreetest->ReleafTreeDegenerate(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, RemoveModelTest) {

  _TheTreetest->RemoveModel();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SampleAncestorsBySequenceTest) {

  //_TheTreetest->SampleAncestorsBySequence();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScaledBranchMappingTest) {

  //node<nodeCoord>* resultnode<nodeCoord> = _TheTreetest->ScaledBranchMapping(node<nodeCoord>test);
  //EXPECT_EQ (resultnode<nodeCoord>*, 0);

}


TEST_F(DISABLED__TheTreeTest, ScaledBranchReMappingTest) {

  //_TheTreetest->ScaledBranchReMapping(node<nodeCoord>test, *_Parametertest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScanAndAttachVariablesTest) {

  _TheTreetest->ScanAndAttachVariables();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScanForCVariablesTest) {

  _TheTreetest->ScanForCVariables(*_AVLListtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScanForDVariablesTest) {

  _TheTreetest->ScanForDVariables(*_AVLListtest, *_AVLListtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScanForGVariablesTest) {

  _TheTreetest->ScanForGVariables(*_AVLListtest, *_AVLListtest, _AVLListXtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScanForVariablesTest) {

  _TheTreetest->ScanForVariables(*_AVLListtest, *_AVLListtest, _AVLListXtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ScanSubtreeVarsTest) {

  _TheTreetest->ScanSubtreeVars(*_Listtest, *chartest, _CalcNodetest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SerialMatrixUpdateTest) {

  _TheTreetest->SerialMatrixUpdate(*longtest, *booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SetCompMatricesTest) {

  _TheTreetest->SetCompMatrices(*longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SetTreeCodeBaseTest) {

  _TheTreetest->SetTreeCodeBase(*longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SetUpTest) {

  _TheTreetest->SetUp();
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SetUpMatricesTest) {

  _TheTreetest->SetUpMatrices(*longtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, SetupCategoryMapsForNodesTest) {

  //_TheTreetest->SetupCategoryMapsForNodes(*_Listtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, StepWiseTraversalTest) {

  _CalcNode* result_CalcNode = _TheTreetest->StepWiseTraversal(*booltest);
  //EXPECT_EQ (result_CalcNode*, 0);

}


TEST_F(DISABLED__TheTreeTest, StepWiseTraversalLevelTest) {

  _CalcNode* result_CalcNode = _TheTreetest->StepWiseTraversalLevel(*longtest, *booltest);
  //EXPECT_EQ (result_CalcNode*, 0);

}


TEST_F(DISABLED__TheTreeTest, TEXTreeStringTest) {

  _PMathObj result_PMathObj = _TheTreetest->TEXTreeString(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__TheTreeTest, ThreadMatrixUpdateTest) {

  _TheTreetest->ThreadMatrixUpdate(*longtest, *booltest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, ThreadReleafTreeCacheTest) {

  //_Parameter result_Parameter = _TheTreetest->ThreadReleafTreeCache(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ThreadReleafTreeChar4Test) {

  //_Parameter result_Parameter = _TheTreetest->ThreadReleafTreeChar4(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, ThreadReleafTreeCharCacheTest) {

  //_Parameter result_Parameter = _TheTreetest->ThreadReleafTreeCharCache(_DataSetFiltertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, TreePSRecurseTest) {

  //_TheTreetest->TreePSRecurse(node<nodeCoord>test, *_Stringtest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, TreeTEXRecurseTest) {

  //nodeCoord resultnodeCoord = _TheTreetest->TreeTEXRecurse(node<nodeCoord>test, *_Stringtest);
  //EXPECT_EQ (resultnodeCoord, 0);

}


TEST_F(DISABLED__TheTreeTest, TreeUserParamsTest) {

  _String* result_String = _TheTreetest->TreeUserParams();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(DISABLED__TheTreeTest, VerySimpleLikelihoodEvaluatorTest) {

  //_Parameter result_Parameter = _TheTreetest->VerySimpleLikelihoodEvaluator();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__TheTreeTest, WeightedCharacterDifferencesTest) {

  //_TheTreetest->WeightedCharacterDifferences(*_Parametertest);
  //EXPECT_EQ (_TheTreetest, 0);

}


TEST_F(DISABLED__TheTreeTest, makeDynamicTest) {

  BaseRef resultBaseRef = _TheTreetest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__TheTreeTest, makeDynamicCopyTest) {

  BaseRef resultBaseRef = _TheTreetest->makeDynamicCopy(_Stringtest);
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__TheTreeTest, toStrTest) {

  BaseRef resultBaseRef = _TheTreetest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
