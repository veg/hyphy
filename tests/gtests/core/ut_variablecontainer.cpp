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
#include "defines.h"
#include "variablecontainer.h"
#include "operation.h"
#include "likefunc.h"
#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"

namespace {

// The fixture for testing class Foo.
class _VariableContainerTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _VariableContainerTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _Listtest = new _List();
    BaseReftest = new BaseRef();
    _AVLListXLtest = new _AVLListXL(_SimpleListtest);
    _SimpleListtest = new _SimpleList();
    _VariableContainertest = new _VariableContainer();
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
  }

  virtual ~_VariableContainerTest() {
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
    delete _Listtest;
    delete BaseReftest;
    delete _AVLListXLtest;
    delete _SimpleListtest;
    delete _VariableContainertest;
    delete booltest;
    delete _AVLListtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  _List* _Listtest;
  BaseRef* BaseReftest;
  _AVLListXL* _AVLListXLtest;
  _SimpleList* _SimpleListtest;
  _VariableContainer* _VariableContainertest;
  bool* booltest;
  _AVLList* _AVLListtest;
};


TEST_F(_VariableContainerTest, CheckAndAddUserExpressionTest) {

  long resultlong = _VariableContainertest->CheckAndAddUserExpression(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_VariableContainerTest, ClearTest) {

  _VariableContainertest->Clear();
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, ClearConstraintsTest) {

  _VariableContainertest->ClearConstraints();
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, CompileListOfDependentsTest) {

  _VariableContainertest->CompileListOfDependents(*_SimpleListtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, CopyMatrixParametersTest) {

  _VariableContainertest->CopyMatrixParameters(_VariableContainertest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, CountAllTest) {

  long resultlong = _VariableContainertest->CountAll();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_VariableContainerTest, CountIndependentsTest) {

  long resultlong = _VariableContainertest->CountIndependents();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_VariableContainerTest, DuplicateTest) {

  _VariableContainertest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, GetExplicitFormModelTest) {

  _Formula* result_Formula = _VariableContainertest->GetExplicitFormModel();
  //EXPECT_EQ (result_Formula*, 0);

}


TEST_F(_VariableContainerTest, GetFreqMatrixTest) {

  _Matrix* result_Matrix = _VariableContainertest->GetFreqMatrix();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(_VariableContainerTest, GetIthDependentTest) {

  _Variable* result_Variable = _VariableContainertest->GetIthDependent(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(_VariableContainerTest, GetIthIndependentTest) {

  _Variable* result_Variable = _VariableContainertest->GetIthIndependent(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(_VariableContainerTest, GetIthParameterTest) {

  _Variable* result_Variable = _VariableContainertest->GetIthParameter(*longtest);
  //EXPECT_EQ (result_Variable*, 0);

}


TEST_F(_VariableContainerTest, GetListOfModelParametersTest) {

  _VariableContainertest->GetListOfModelParameters(*_Listtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, GetModelDimensionTest) {

  long resultlong = _VariableContainertest->GetModelDimension();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_VariableContainerTest, GetModelMatrixTest) {

  _Matrix* result_Matrix = _VariableContainertest->GetModelMatrix(_Listtest, _SimpleListtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(_VariableContainerTest, GetModelNameTest) {

  _String result_String = _VariableContainertest->GetModelName();
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_VariableContainerTest, GetSaveableListOfUserParametersTest) {

  _String* result_String = _VariableContainertest->GetSaveableListOfUserParameters();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_VariableContainerTest, HasChangedTest) {

  bool resultbool = _VariableContainertest->HasChanged();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, HasExplicitFormModelTest) {

  bool resultbool = _VariableContainertest->HasExplicitFormModel();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, HasLocalsTest) {

  bool resultbool = _VariableContainertest->HasLocals();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, InitializeVarContTest) {

  //_VariableContainertest->InitializeVarCont(*_Stringtest, *_Stringtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, IsConstantTest) {

  bool resultbool = _VariableContainertest->IsConstant();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, IsModelVarTest) {

  bool resultbool = _VariableContainertest->IsModelVar(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, KillUserExpressionTest) {

  _VariableContainertest->KillUserExpression(*longtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, MarkDoneTest) {

  _VariableContainertest->MarkDone();
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, MatchParametersToListTest) {

  _VariableContainertest->MatchParametersToList(*_Listtest, *booltest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, NeedToExponentiateTest) {

  bool resultbool = _VariableContainertest->NeedToExponentiate(*booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, RemoveDependanceTest) {

  bool resultbool = _VariableContainertest->RemoveDependance(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, ScanAndAttachVariablesTest) {

  _VariableContainertest->ScanAndAttachVariables();
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, ScanForDVariablesTest) {

  _VariableContainertest->ScanForDVariables(*_AVLListtest, *_AVLListtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, ScanForGVariablesTest) {

  _VariableContainertest->ScanForGVariables(*_AVLListtest, *_AVLListtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, ScanForVariablesTest) {

  _VariableContainertest->ScanForVariables(*_AVLListtest, *_AVLListtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, ScanModelBasedVariablesTest) {

  //_VariableContainertest->ScanModelBasedVariables(*_Stringtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, SetDependanceTest) {

  long resultlong = _VariableContainertest->SetDependance(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_VariableContainerTest, SetMDependanceTest) {

  bool resultbool = _VariableContainertest->SetMDependance(*_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_VariableContainerTest, SetModelTest) {

  _VariableContainertest->SetModel(*longtest, _AVLListXLtest);
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, SortVarsTest) {

  //_VariableContainertest->SortVars();
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, TrimMemoryTest) {

  _VariableContainertest->TrimMemory();
  //EXPECT_EQ (_VariableContainertest, 0);

}


TEST_F(_VariableContainerTest, makeDynamicTest) {

  BaseRef resultBaseRef = _VariableContainertest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_VariableContainerTest, toStrTest) {

  BaseRef resultBaseRef = _VariableContainertest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
