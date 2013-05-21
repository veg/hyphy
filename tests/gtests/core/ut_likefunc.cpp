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
class _LikelihoodFunctionTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _LikelihoodFunctionTest() {
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

  virtual ~_LikelihoodFunctionTest() {
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


TEST_F(_LikelihoodFunctionTest, AllocateSiteResultsTest) {

  //_LikelihoodFunctiontest->AllocateSiteResults();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, AnnealTest) {

  _LikelihoodFunctiontest->Anneal(*_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, BlockLengthTest) {

  //long resultlong = _LikelihoodFunctiontest->BlockLength(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, BracketTest) {

  //long resultlong = _LikelihoodFunctiontest->Bracket();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, BuildIncrementsTest) {

  //_LikelihoodFunctiontest->BuildIncrements(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, BuildLeafProbsTest) {

  //_LikelihoodFunctiontest->BuildLeafProbs();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CheckAndSetIthIndependentTest) {

  bool resultbool = _LikelihoodFunctiontest->CheckAndSetIthIndependent(*longtest, *_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, CheckDependentBoundsTest) {

  //_LikelihoodFunctiontest->CheckDependentBounds();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CheckFibonacciTest) {

  //_LikelihoodFunctiontest->CheckFibonacci(*_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CheckNthBitTest) {

  //bool resultbool = _LikelihoodFunctiontest->CheckNthBit(*longtest, *chartest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, CheckStepTest) {

  //_LikelihoodFunctiontest->CheckStep(*_Parametertest, *_Matrixtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CleanUpOptimizeTest) {

  //_LikelihoodFunctiontest->CleanUpOptimize();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CleanupTest) {

  //_LikelihoodFunctiontest->Cleanup();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CleanupMPIOptimizerTest) {

  //_LikelihoodFunctiontest->CleanupMPIOptimizer();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CleanupParameterMappingTest) {

  //_LikelihoodFunctiontest->CleanupParameterMapping();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ClearTest) {

  //_LikelihoodFunctiontest->Clear();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CodonNeutralSimulateTest) {

  //_LikelihoodFunctiontest->CodonNeutralSimulate();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, CollectLFAttributesTest) {

  _AssociativeList* result_AssociativeList = _LikelihoodFunctiontest->CollectLFAttributes();
  //EXPECT_EQ (result_AssociativeList*, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputeTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->Compute();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputeBlockTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->ComputeBlock();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputeBlockForTemplateTest) {

  //_LikelihoodFunctiontest->ComputeBlockForTemplate(*longtest, *booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputeBlockForTemplate2Test) {

  //_LikelihoodFunctiontest->ComputeBlockForTemplate2();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputeGradientTest) {

  //_LikelihoodFunctiontest->ComputeGradient();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputePruningEfficiencyTest) {

  _LikelihoodFunctiontest->ComputePruningEfficiency(*longtest, *longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ComputeSiteLikelihoodsForABlockTest) {

  //_LikelihoodFunctiontest->ComputeSiteLikelihoodsForABlock();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ConjugateGradientDescentTest) {

  //_LikelihoodFunctiontest->ConjugateGradientDescent();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ConstructTest) {

  bool resultbool = _LikelihoodFunctiontest->Construct(*_Listtest, _VariableContainertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, CostOfPathTest) {

  //long resultlong = _LikelihoodFunctiontest->CostOfPath(_DataSetFiltertest, _TheTreetest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, CountObjectsTest) {

  long resultlong = _LikelihoodFunctiontest->CountObjects(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, CovarianceMatrixTest) {

  _PMathObj result_PMathObj = _LikelihoodFunctiontest->CovarianceMatrix(_SimpleListtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_LikelihoodFunctionTest, DeleteCachesTest) {

  //_LikelihoodFunctiontest->DeleteCaches(*booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, DependOnDSTest) {

  long resultlong = _LikelihoodFunctiontest->DependOnDS(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, DependOnModelTest) {

  long resultlong = _LikelihoodFunctiontest->DependOnModel(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, DependOnTreeTest) {

  long resultlong = _LikelihoodFunctiontest->DependOnTree(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, DetermineLocalUpdatePolicyTest) {

  //_LikelihoodFunctiontest->DetermineLocalUpdatePolicy();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, DoneComputingTest) {

  _LikelihoodFunctiontest->DoneComputing(*booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, DuplicateTest) {

  _LikelihoodFunctiontest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, FillInConditionalsTest) {

  _LikelihoodFunctiontest->FillInConditionals(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, FindCategoryVarTest) {

  _CategoryVariable* result_CategoryVariable = _LikelihoodFunctiontest->FindCategoryVar(*longtest);
  //EXPECT_EQ (result_CategoryVariable*, 0);

}


TEST_F(_LikelihoodFunctionTest, FlushLocalUpdatePolicyTest) {

  //_LikelihoodFunctiontest->FlushLocalUpdatePolicy();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, GetAllIndependentTest) {

  _LikelihoodFunctiontest->GetAllIndependent(*_Matrixtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, GetCategoryVarsTest) {

  _SimpleList result_SimpleList = _LikelihoodFunctiontest->GetCategoryVars();
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(_LikelihoodFunctionTest, GetDependentVarsTest) {

  _SimpleList result_SimpleList = _LikelihoodFunctiontest->GetDependentVars();
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(_LikelihoodFunctionTest, GetGlobalVarsTest) {

  _LikelihoodFunctiontest->GetGlobalVars(*_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, GetGradientStepBoundTest) {

  //_LikelihoodFunctiontest->GetGradientStepBound();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, GetIndependentVarsTest) {

  _SimpleList result_SimpleList = _LikelihoodFunctiontest->GetIndependentVars();
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(_LikelihoodFunctionTest, GetInitialValuesTest) {

  //_LikelihoodFunctiontest->GetInitialValues();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, GetIthDependentTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->GetIthDependent(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, GetIthDependentVarTest) {

  _Variable* result_Variable = _LikelihoodFunctiontest->GetIthDependentVar(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(_LikelihoodFunctionTest, GetIthIndependentTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->GetIthIndependent(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, GetIthIndependentBoundTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->GetIthIndependentBound(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, GetIthIndependentVarTest) {

  _Variable* result_Variable = _LikelihoodFunctiontest->GetIthIndependentVar(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(_LikelihoodFunctionTest, GradientDescentTest) {

  //_LikelihoodFunctiontest->GradientDescent(*_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, GradientLocateTheBumpTest) {

  //_LikelihoodFunctiontest->GradientLocateTheBump();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, HasBlockChangedTest) {

  //bool resultbool = _LikelihoodFunctiontest->HasBlockChanged(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, HasHiddenMarkovTest) {

  //long resultlong = _LikelihoodFunctiontest->HasHiddenMarkov(*longtest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, HasPartitionChangedTest) {

  //bool resultbool = _LikelihoodFunctiontest->HasPartitionChanged(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, HasPrecisionBeenAchievedTest) {

  //long resultlong = _LikelihoodFunctiontest->HasPrecisionBeenAchieved(*_Parametertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, HighestBitTest) {

  //char resultchar = _LikelihoodFunctiontest->HighestBit(*longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_LikelihoodFunctionTest, InitTest) {

  _LikelihoodFunctiontest->Init();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, InitMPIOptimizerTest) {

  //_LikelihoodFunctiontest->InitMPIOptimizer();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, IsIthParameterGlobalTest) {

  bool resultbool = _LikelihoodFunctiontest->IsIthParameterGlobal(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, LocateTheBumpTest) {

  //_LikelihoodFunctiontest->LocateTheBump();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, LowestBitTest) {

  //char resultchar = _LikelihoodFunctiontest->LowestBit(*longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_LikelihoodFunctionTest, MPI_LF_ComputeTest) {

  _LikelihoodFunctiontest->MPI_LF_Compute(*longtest, *booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, MapTreeTipsToDataTest) {

  bool resultbool = _LikelihoodFunctiontest->MapTreeTipsToData(*longtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, MaximumDimensionTest) {

  long resultlong = _LikelihoodFunctiontest->MaximumDimension();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, OptimalOrderTest) {

  //_LikelihoodFunctiontest->OptimalOrder(*longtest, *_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, OptimizeTest) {

  _Matrix* result_Matrix = _LikelihoodFunctiontest->Optimize();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(_LikelihoodFunctionTest, PartitionCatVarsTest) {

  //_LikelihoodFunctiontest->PartitionCatVars(*_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, PartitionLengthsTest) {

  //long resultlong = _LikelihoodFunctiontest->PartitionLengths(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, PopulateConditionalProbabilitiesTest) {

  //_LikelihoodFunctiontest->PopulateConditionalProbabilities();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, PostComputeTest) {

  _LikelihoodFunctiontest->PostCompute();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, PreComputeTest) {

  bool resultbool = _LikelihoodFunctiontest->PreCompute();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, PrepareToComputeTest) {

  _LikelihoodFunctiontest->PrepareToCompute(*booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ProcessPartitionListTest) {

  //bool resultbool = _LikelihoodFunctiontest->ProcessPartitionList(*_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, RandomizeListTest) {

  //_LikelihoodFunctiontest->RandomizeList(*_SimpleListtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RankVariablesTest) {

  _LikelihoodFunctiontest->RankVariables(_AVLListXtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RebuildTest) {

  _LikelihoodFunctiontest->Rebuild();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ReconstructAncestorsTest) {

  //_LikelihoodFunctiontest->ReconstructAncestors(*_DataSettest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RecoverAncestralSequencesMarginalTest) {

  //_List* result_List = _LikelihoodFunctiontest->RecoverAncestralSequencesMarginal();
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(_LikelihoodFunctionTest, RecurseCategoryTest) {

  //_LikelihoodFunctiontest->RecurseCategory();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RecurseConstantOnPartitionTest) {

  //_LikelihoodFunctiontest->RecurseConstantOnPartition();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RemapMatrixTest) {

  //_Matrix* result_Matrix = _LikelihoodFunctiontest->RemapMatrix(_Matrixtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(_LikelihoodFunctionTest, RescanAllVariablesTest) {

  _LikelihoodFunctiontest->RescanAllVariables();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RestoreScalingFactorsTest) {

  //_LikelihoodFunctiontest->RestoreScalingFactors(*longtest, *longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, RunViterbiTest) {

  //_LikelihoodFunctiontest->RunViterbi(*_Matrixtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ScanAllVariablesOnPartitionTest) {

  //_LikelihoodFunctiontest->ScanAllVariablesOnPartition();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SendOffToMPITest) {

  //bool resultbool = _LikelihoodFunctiontest->SendOffToMPI(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, SequenceCountTest) {

  long resultlong = _LikelihoodFunctiontest->SequenceCount(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, SerializeLFTest) {

  _LikelihoodFunctiontest->SerializeLF(*_Stringtest, *chartest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SetAllIndependentTest) {

  long resultlong = _LikelihoodFunctiontest->SetAllIndependent(_Matrixtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, SetIthDependentTest) {

  _LikelihoodFunctiontest->SetIthDependent(*longtest, *_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SetIthIndependentTest) {

  _LikelihoodFunctiontest->SetIthIndependent(*longtest, *_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SetNthBitTest) {

  //_LikelihoodFunctiontest->SetNthBit(*longtest, *chartest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SetParametersAndComputeTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->SetParametersAndCompute();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, SetupTest) {

  _LikelihoodFunctiontest->Setup();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SetupLFCachesTest) {

  //_LikelihoodFunctiontest->SetupLFCaches();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SetupParameterMappingTest) {

  //_LikelihoodFunctiontest->SetupParameterMapping();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SimplexMethodTest) {

  _Parameter result_Parameter = _LikelihoodFunctiontest->SimplexMethod(*_Parametertest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, SimulateTest) {

  //_LikelihoodFunctiontest->Simulate();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SimulateCodonNeutralTest) {

  //_AssociativeList* result_AssociativeList = _LikelihoodFunctiontest->SimulateCodonNeutral();
  //EXPECT_EQ (result_AssociativeList*, 0);

}


TEST_F(_LikelihoodFunctionTest, SingleBuildLeafProbsTest) {

  //bool resultbool = _LikelihoodFunctiontest->SingleBuildLeafProbs();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, SiteCountTest) {

  long resultlong = _LikelihoodFunctiontest->SiteCount();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, SniffAroundTest) {

  //bool resultbool = _LikelihoodFunctiontest->SniffAround(*_Matrixtest, *_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, StateCounterTest) {

  _LikelihoodFunctiontest->StateCounter(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, SumUpHiddenMarkovTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->SumUpHiddenMarkov();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, TotalRateClassesForAPartitionTest) {

  long resultlong = _LikelihoodFunctiontest->TotalRateClassesForAPartition(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_LikelihoodFunctionTest, UpdateBlockResultTest) {

  //_LikelihoodFunctiontest->UpdateBlockResult(*longtest, *_Parametertest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, UpdateDependentTest) {

  _LikelihoodFunctiontest->UpdateDependent(*longtest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, UpdateFilterSizeTest) {

  bool resultbool = _LikelihoodFunctiontest->UpdateFilterSize(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, UpdateIndependentTest) {

  _LikelihoodFunctiontest->UpdateIndependent(*longtest, *booltest);
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, ZeroSiteResultsTest) {

  //_LikelihoodFunctiontest->ZeroSiteResults();
  //EXPECT_EQ (_LikelihoodFunctiontest, 0);

}


TEST_F(_LikelihoodFunctionTest, checkPermissibilityTest) {

  //bool resultbool = _LikelihoodFunctiontest->checkPermissibility(*_Matrixtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_LikelihoodFunctionTest, computeAtAPointTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->computeAtAPoint(*_Matrixtest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, makeDynamicTest) {

  BaseRef resultBaseRef = _LikelihoodFunctiontest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_LikelihoodFunctionTest, replaceAPointTest) {

  //_Parameter result_Parameter = _LikelihoodFunctiontest->replaceAPoint();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_LikelihoodFunctionTest, toStrTest) {

  BaseRef resultBaseRef = _LikelihoodFunctiontest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
