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
#include "variable.h"
#include "operation.h"
#include "likefunc.h"
#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"
#include "thetree.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__VariableTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__VariableTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _Variabletest = new _Variable(*_Stringtest, *_Stringtest, *booltest);
    _Formulatest = new _Formula(*_PMathObjtest, *booltest);
    BaseReftest = new BaseRef();
    _AVLListXtest = new _AVLListX(_SimpleListtest);
    _SimpleListtest = new _SimpleList();
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
    _Parametertest = new _Parameter();
    _PMathObjtest = new _PMathObj();
  }

  virtual ~DISABLED__VariableTest() {
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
    delete _Variabletest;
    delete _Formulatest;
    delete BaseReftest;
    delete _AVLListXtest;
    delete _SimpleListtest;
    delete booltest;
    delete _AVLListtest;
    delete _Parametertest;
    delete _PMathObjtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  _Variable* _Variabletest;
  _Formula* _Formulatest;
  BaseRef* BaseReftest;
  _AVLListX* _AVLListXtest;
  _SimpleList* _SimpleListtest;
  bool* booltest;
  _AVLList* _AVLListtest;
  _Parameter* _Parametertest;
  _PMathObj* _PMathObjtest;
};


TEST_F(DISABLED__VariableTest, CheckAndSetTest) {

  _Variabletest->CheckAndSet(*_Parametertest, *booltest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, CheckFForDependenceTest) {

  bool resultbool = _Variabletest->CheckFForDependence(*longtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__VariableTest, ClearConstraintsTest) {

  _Variabletest->ClearConstraints();
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, CompileListOfDependentsTest) {

  _Variabletest->CompileListOfDependents(*_SimpleListtest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, ComputeTest) {

  _PMathObj result_PMathObj = _Variabletest->Compute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__VariableTest, ComputeReferenceTest) {

  _PMathObj result_PMathObj = _Variabletest->ComputeReference(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__VariableTest, ContextFreeNameTest) {

  _String result_String = _Variabletest->ContextFreeName();
  //EXPECT_EQ (result_String, 0);

}


TEST_F(DISABLED__VariableTest, DuplicateTest) {

  _Variabletest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, EnsureTheValueIsInBoundsTest) {

  _Variabletest->EnsureTheValueIsInBounds();
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, HasChangedTest) {

  bool resultbool = _Variabletest->HasChanged(*booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__VariableTest, InitializeTest) {

  _Variabletest->Initialize();
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, IsConstantTest) {

  bool resultbool = _Variabletest->IsConstant();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__VariableTest, IsVariableTest) {

  bool resultbool = _Variabletest->IsVariable();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__VariableTest, MarkDoneTest) {

  _Variabletest->MarkDone();
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, ParentObjectNameTest) {

  _String result_String = _Variabletest->ParentObjectName();
  //EXPECT_EQ (result_String, 0);

}


TEST_F(DISABLED__VariableTest, PostMarkChangedTest) {

  _Variabletest->PostMarkChanged();
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, PreMarkChangedTest) {

  _Variabletest->PreMarkChanged();
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, ScanForVariablesTest) {

  _Variabletest->ScanForVariables(*_AVLListtest, *booltest, _AVLListXtest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, SetBoundsTest) {

  _Variabletest->SetBounds(*_Parametertest, *_Parametertest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, SetFormulaTest) {

  _Variabletest->SetFormula(*_Formulatest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, SetNumericValueTest) {

  _Variabletest->SetNumericValue(*_Parametertest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, SetValueTest) {

  _Variabletest->SetValue(*_PMathObjtest, *booltest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, makeDynamicTest) {

  BaseRef resultBaseRef = _Variabletest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__VariableTest, toFileStrTest) {

  _Variabletest->toFileStr(FILEtest);
  //EXPECT_EQ (_Variabletest, 0);

}


TEST_F(DISABLED__VariableTest, toStrTest) {

  BaseRef resultBaseRef = _Variabletest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
