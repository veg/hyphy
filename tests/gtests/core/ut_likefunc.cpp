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
#include "hy_globals.h"
#include "helperfunctions.h"
#include "thetree.h"
#include "datasetfilter.h"
#include "executionlist.h"
#include "scfg.h"
    #include <pthread.h>

namespace {

// The fixture for testing class Foo.
class DISABLED__LikelihoodFunctionTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__LikelihoodFunctionTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    chartest = new char();
    _Listtest = new _List();
    _TheTreetest = new _TheTree();
    _DataSettest = new _DataSet(FILEtest);
    BaseReftest = new BaseRef();
    _Matrixtest = new _Matrix();
    _AVLListXtest = new _AVLListX(_SimpleListtest);
    _SimpleListtest = new _SimpleList();
    _DataSetFiltertest = new _DataSetFilter(_DataSettest, *chartest, *_Stringtest);
    _VariableContainertest = new _VariableContainer();
    booltest = new bool();
    _Parametertest = new _Parameter();
    _LikelihoodFunctiontest = new _LikelihoodFunction(*_LikelihoodFunctiontest);
  }

  virtual ~DISABLED__LikelihoodFunctionTest() {
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
    delete chartest;
    delete _Listtest;
    delete _TheTreetest;
    delete _DataSettest;
    delete BaseReftest;
    delete _Matrixtest;
    delete _AVLListXtest;
    delete _SimpleListtest;
    delete _DataSetFiltertest;
    delete _VariableContainertest;
    delete booltest;
    delete _Parametertest;
    delete _LikelihoodFunctiontest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  char* chartest;
  _List* _Listtest;
  _TheTree* _TheTreetest;
  _DataSet* _DataSettest;
  BaseRef* BaseReftest;
  _Matrix* _Matrixtest;
  _AVLListX* _AVLListXtest;
  _SimpleList* _SimpleListtest;
  _DataSetFilter* _DataSetFiltertest;
  _VariableContainer* _VariableContainertest;
  bool* booltest;
  _Parameter* _Parametertest;
  _LikelihoodFunction* _LikelihoodFunctiontest;
};


TEST_F(DISABLED__LikelihoodFunctionTest, AllocateSiteResultsTest) {

  //_LikelihoodFunctiontest->AllocateSiteResults();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, AnnealTest) {

  _LikelihoodFunctiontest->Anneal(*_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, BlockLengthTest) {

  //long resultlong = _LikelihoodFunctiontest->BlockLength(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, BracketTest) {

  //long resultlong = _LikelihoodFunctiontest->Bracket();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, BuildIncrementsTest) {

  //_LikelihoodFunctiontest->BuildIncrements(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, BuildLeafProbsTest) {

  //_LikelihoodFunctiontest->BuildLeafProbs();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CheckAndSetIthIndependentTest) {

  bool resultbool = _LikelihoodFunctiontest->CheckAndSetIthIndependent(*longtest, *_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CheckDependentBoundsTest) {

  //_LikelihoodFunctiontest->CheckDependentBounds();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CheckFibonacciTest) {

  //_LikelihoodFunctiontest->CheckFibonacci(*_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CheckNthBitTest) {

  //bool resultbool = _LikelihoodFunctiontest->CheckNthBit(*longtest, *chartest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CheckStepTest) {

  //_LikelihoodFunctiontest->CheckStep(*_Parametertest, *_Matrixtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CleanUpOptimizeTest) {

  //_LikelihoodFunctiontest->CleanUpOptimize();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CleanupTest) {

  //_LikelihoodFunctiontest->Cleanup();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CleanupMPIOptimizerTest) {

  //_LikelihoodFunctiontest->CleanupMPIOptimizer();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CleanupParameterMappingTest) {

  //_LikelihoodFunctiontest->CleanupParameterMapping();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ClearTest) {

  //_LikelihoodFunctiontest->Clear();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CodonNeutralSimulateTest) {

  //_LikelihoodFunctiontest->CodonNeutralSimulate();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CollectLFAttributesTest) {

  _AssociativeList* result_AssociativeList = _LikelihoodFunctiontest->CollectLFAttributes();
  //EXPECT_EQ (result_AssociativeList*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputeTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->Compute();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputeBlockTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->ComputeBlock();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputeBlockForTemplateTest) {

  //_LikelihoodFunctiontest->ComputeBlockForTemplate(*longtest, *booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputeBlockForTemplate2Test) {

  //_LikelihoodFunctiontest->ComputeBlockForTemplate2();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputeGradientTest) {

  //_LikelihoodFunctiontest->ComputeGradient();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputePruningEfficiencyTest) {

  _LikelihoodFunctiontest->ComputePruningEfficiency(*longtest, *longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ComputeSiteLikelihoodsForABlockTest) {

  //_LikelihoodFunctiontest->ComputeSiteLikelihoodsForABlock();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ConjugateGradientDescentTest) {

  //_LikelihoodFunctiontest->ConjugateGradientDescent();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ConstructTest) {

  bool resultbool = _LikelihoodFunctiontest->Construct(*_Listtest, _VariableContainertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CostOfPathTest) {

  //long resultlong = _LikelihoodFunctiontest->CostOfPath(_DataSetFiltertest, _TheTreetest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CountObjectsTest) {

  long resultlong = _LikelihoodFunctiontest->CountObjects(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, CovarianceMatrixTest) {

  _PMathObj result_PMathObj = _LikelihoodFunctiontest->CovarianceMatrix(_SimpleListtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DeleteCachesTest) {

  //_LikelihoodFunctiontest->DeleteCaches(*booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DependOnDSTest) {

  long resultlong = _LikelihoodFunctiontest->DependOnDS(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DependOnModelTest) {

  long resultlong = _LikelihoodFunctiontest->DependOnModel(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DependOnTreeTest) {

  long resultlong = _LikelihoodFunctiontest->DependOnTree(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DetermineLocalUpdatePolicyTest) {

  //_LikelihoodFunctiontest->DetermineLocalUpdatePolicy();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DoneComputingTest) {

  _LikelihoodFunctiontest->DoneComputing(*booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, DuplicateTest) {

  _LikelihoodFunctiontest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, FillInConditionalsTest) {

  _LikelihoodFunctiontest->FillInConditionals(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, FindCategoryVarTest) {

  _CategoryVariable* result_CategoryVariable = _LikelihoodFunctiontest->FindCategoryVar(*longtest);
  //EXPECT_EQ (result_CategoryVariable*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, FlushLocalUpdatePolicyTest) {

  //_LikelihoodFunctiontest->FlushLocalUpdatePolicy();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetAllIndependentTest) {

  _LikelihoodFunctiontest->GetAllIndependent(*_Matrixtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetCategoryVarsTest) {

  _SimpleList result_SimpleList = _LikelihoodFunctiontest->GetCategoryVars();
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetDependentVarsTest) {

  _SimpleList result_SimpleList = _LikelihoodFunctiontest->GetDependentVars();
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetGlobalVarsTest) {

  _LikelihoodFunctiontest->GetGlobalVars(*_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetGradientStepBoundTest) {

  //_LikelihoodFunctiontest->GetGradientStepBound();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetIndependentVarsTest) {

  _SimpleList result_SimpleList = _LikelihoodFunctiontest->GetIndependentVars();
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetInitialValuesTest) {

  //_LikelihoodFunctiontest->GetInitialValues();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetIthDependentTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->GetIthDependent(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetIthDependentVarTest) {

  _Variable* result_Variable = _LikelihoodFunctiontest->GetIthDependentVar(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetIthIndependentTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->GetIthIndependent(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetIthIndependentBoundTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->GetIthIndependentBound(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GetIthIndependentVarTest) {

  _Variable* result_Variable = _LikelihoodFunctiontest->GetIthIndependentVar(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GradientDescentTest) {

  //_LikelihoodFunctiontest->GradientDescent(*_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, GradientLocateTheBumpTest) {

  //_LikelihoodFunctiontest->GradientLocateTheBump();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, HasBlockChangedTest) {

  //bool resultbool = _LikelihoodFunctiontest->HasBlockChanged(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, HasHiddenMarkovTest) {

  //long resultlong = _LikelihoodFunctiontest->HasHiddenMarkov(*longtest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, HasPartitionChangedTest) {

  //bool resultbool = _LikelihoodFunctiontest->HasPartitionChanged(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, HasPrecisionBeenAchievedTest) {

  //long resultlong = _LikelihoodFunctiontest->HasPrecisionBeenAchieved(*_Parametertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, HighestBitTest) {

  //char resultchar = _LikelihoodFunctiontest->HighestBit(*longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, InitTest) {

  _LikelihoodFunctiontest->Init();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, InitMPIOptimizerTest) {

  //_LikelihoodFunctiontest->InitMPIOptimizer();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, IsIthParameterGlobalTest) {

  bool resultbool = _LikelihoodFunctiontest->IsIthParameterGlobal(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, LocateTheBumpTest) {

  //_LikelihoodFunctiontest->LocateTheBump();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, LowestBitTest) {

  //char resultchar = _LikelihoodFunctiontest->LowestBit(*longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, MPIDISABLED__LF_ComputeTest) {

  _LikelihoodFunctiontest->MPI_LF_Compute(*longtest, *booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, MapTreeTipsToDataTest) {

  bool resultbool = _LikelihoodFunctiontest->MapTreeTipsToData(*longtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, MaximumDimensionTest) {

  long resultlong = _LikelihoodFunctiontest->MaximumDimension();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, OptimalOrderTest) {

  //_LikelihoodFunctiontest->OptimalOrder(*longtest, *_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, OptimizeTest) {

  _Matrix* result_Matrix = _LikelihoodFunctiontest->Optimize();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, PartitionCatVarsTest) {

  //_LikelihoodFunctiontest->PartitionCatVars(*_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, PartitionLengthsTest) {

  //long resultlong = _LikelihoodFunctiontest->PartitionLengths(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, PopulateConditionalProbabilitiesTest) {

  //_LikelihoodFunctiontest->PopulateConditionalProbabilities();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, PostComputeTest) {

  _LikelihoodFunctiontest->PostCompute();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, PreComputeTest) {

  bool resultbool = _LikelihoodFunctiontest->PreCompute();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, PrepareToComputeTest) {

  _LikelihoodFunctiontest->PrepareToCompute(*booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ProcessPartitionListTest) {

  //bool resultbool = _LikelihoodFunctiontest->ProcessPartitionList(*_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RandomizeListTest) {

  //_LikelihoodFunctiontest->RandomizeList(*_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RankVariablesTest) {

  _LikelihoodFunctiontest->RankVariables(_AVLListXtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RebuildTest) {

  _LikelihoodFunctiontest->Rebuild();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ReconstructAncestorsTest) {

  //_LikelihoodFunctiontest->ReconstructAncestors(*_DataSettest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RecoverAncestralSequencesMarginalTest) {

  //_List* result_List = _LikelihoodFunctiontest->RecoverAncestralSequencesMarginal();
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RecurseCategoryTest) {

  //_LikelihoodFunctiontest->RecurseCategory();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RecurseConstantOnPartitionTest) {

  //_LikelihoodFunctiontest->RecurseConstantOnPartition();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RemapMatrixTest) {

  //_Matrix* result_Matrix = _LikelihoodFunctiontest->RemapMatrix(_Matrixtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RescanAllVariablesTest) {

  _LikelihoodFunctiontest->RescanAllVariables();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RestoreScalingFactorsTest) {

  //_LikelihoodFunctiontest->RestoreScalingFactors(*longtest, *longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, RunViterbiTest) {

  //_LikelihoodFunctiontest->RunViterbi(*_Matrixtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ScanAllVariablesOnPartitionTest) {

  //_LikelihoodFunctiontest->ScanAllVariablesOnPartition();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SendOffToMPITest) {

  //bool resultbool = _LikelihoodFunctiontest->SendOffToMPI(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SequenceCountTest) {

  long resultlong = _LikelihoodFunctiontest->SequenceCount(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SerializeLFTest) {

  _LikelihoodFunctiontest->SerializeLF(*_Stringtest, *chartest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetAllIndependentTest) {

  long resultlong = _LikelihoodFunctiontest->SetAllIndependent(_Matrixtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetIthDependentTest) {

  _LikelihoodFunctiontest->SetIthDependent(*longtest, *_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetIthIndependentTest) {

  _LikelihoodFunctiontest->SetIthIndependent(*longtest, *_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetNthBitTest) {

  //_LikelihoodFunctiontest->SetNthBit(*longtest, *chartest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetParametersAndComputeTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->SetParametersAndCompute();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetupTest) {

  _LikelihoodFunctiontest->Setup();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetupLFCachesTest) {

  //_LikelihoodFunctiontest->SetupLFCaches();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SetupParameterMappingTest) {

  //_LikelihoodFunctiontest->SetupParameterMapping();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SimplexMethodTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->SimplexMethod(*_Parametertest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SimulateTest) {

  //_LikelihoodFunctiontest->Simulate();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SimulateCodonNeutralTest) {

  //_AssociativeList* result_AssociativeList = _LikelihoodFunctiontest->SimulateCodonNeutral();
  //EXPECT_EQ (result_AssociativeList*, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SingleBuildLeafProbsTest) {

  //bool resultbool = _LikelihoodFunctiontest->SingleBuildLeafProbs();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SiteCountTest) {

  long resultlong = _LikelihoodFunctiontest->SiteCount();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SniffAroundTest) {

  //bool resultbool = _LikelihoodFunctiontest->SniffAround(*_Matrixtest, *_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, StateCounterTest) {

  _LikelihoodFunctiontest->StateCounter(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, SumUpHiddenMarkovTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->SumUpHiddenMarkov();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, TotalRateClassesForAPartitionTest) {

  long resultlong = _LikelihoodFunctiontest->TotalRateClassesForAPartition(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, UpdateBlockResultTest) {

  //_LikelihoodFunctiontest->UpdateBlockResult(*longtest, *_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, UpdateDependentTest) {

  _LikelihoodFunctiontest->UpdateDependent(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, UpdateFilterSizeTest) {

  bool resultbool = _LikelihoodFunctiontest->UpdateFilterSize(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, UpdateIndependentTest) {

  _LikelihoodFunctiontest->UpdateIndependent(*longtest, *booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, ZeroSiteResultsTest) {

  //_LikelihoodFunctiontest->ZeroSiteResults();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, checkPermissibilityTest) {

  //bool resultbool = _LikelihoodFunctiontest->checkPermissibility(*_Matrixtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, computeAtAPointTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->computeAtAPoint(*_Matrixtest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, makeDynamicTest) {

  BaseRef resultBaseRef = _LikelihoodFunctiontest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, replaceAPointTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->replaceAPoint();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__LikelihoodFunctionTest, toStrTest) {

  BaseRef resultBaseRef = _LikelihoodFunctiontest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
