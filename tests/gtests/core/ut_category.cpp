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
#include  "category.h"
#include  "math.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__CategoryVariableTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__CategoryVariableTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleListtest = new _SimpleList();
    longtest = new long();
    _CategoryVariabletest = new _CategoryVariable(*_Stringtest, _Listtest, _VariableContainertest);
    _Stringtest = new _String(FILEtest);
    _Listtest = new _List();
    BaseReftest = new BaseRef();
    _Matrixtest = new _Matrix();
    _VariableContainertest = new _VariableContainer();
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
  }

  virtual ~DISABLED__CategoryVariableTest() {
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
    delete _SimpleListtest;
    delete longtest;
    delete _CategoryVariabletest;
    delete _Stringtest;
    delete _Listtest;
    delete BaseReftest;
    delete _Matrixtest;
    delete _VariableContainertest;
    delete booltest;
    delete _AVLListtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _SimpleList* _SimpleListtest;
  _CategoryVariable* _CategoryVariabletest;
  _String* _Stringtest;
  _List* _Listtest;
  BaseRef* BaseReftest;
  _Matrix* _Matrixtest;
  _VariableContainer* _VariableContainertest;
  bool* booltest;
  _AVLList* _AVLListtest;
};


TEST_F(DISABLED__CategoryVariableTest, ChangeNumberOfIntervalsTest) {

  _CategoryVariabletest->ChangeNumberOfIntervals(*longtest);
  //EXPECT_EQ (_CategoryVariabletest, 0);

}

TEST_F(DISABLED__CategoryVariableTest, ComputeHiddenMarkovTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->ComputeHiddenMarkov();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__CategoryVariableTest, ComputeHiddenMarkovFreqsTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->ComputeHiddenMarkovFreqs();
  //EXPECT_EQ (result_Matrix*, 0);

}

TEST_F(DISABLED__CategoryVariableTest, DuplicateTest) {

  _CategoryVariabletest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_CategoryVariabletest, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetCurrentStateTest) {

  long resultlong = _CategoryVariabletest->GetCurrentState();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetHiddenMarkovTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->GetHiddenMarkov();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetHiddenMarkovFreqsTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->GetHiddenMarkovFreqs();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetIntervalEndsTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->GetIntervalEnds();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetIntervalValueTest) {

  _Parameter result_Parameter = _CategoryVariabletest->GetIntervalValue(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetIntervalWeightTest) {

  _Parameter result_Parameter = _CategoryVariabletest->GetIntervalWeight(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetValuesTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->GetValues();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__CategoryVariableTest, GetWeightsTest) {

  _Matrix* result_Matrix = _CategoryVariabletest->GetWeights(*booltest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__CategoryVariableTest, HaveParametersChangedTest) {

  bool resultbool = _CategoryVariabletest->HaveParametersChanged(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__CategoryVariableTest, IsConstantTest) {

  bool resultbool = _CategoryVariabletest->IsConstant();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__CategoryVariableTest, IsGlobalTest) {

  bool resultbool = _CategoryVariabletest->IsGlobal();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__CategoryVariableTest, IsLayeredTest) {

  bool resultbool = _CategoryVariabletest->IsLayered();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__CategoryVariableTest, IsUncorrelatedTest) {

  bool resultbool = _CategoryVariabletest->IsUncorrelated();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__CategoryVariableTest, MeanTest) {

  _Parameter result_Parameter = _CategoryVariabletest->Mean();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__CategoryVariableTest, ScanForGVariablesTest) {

  _CategoryVariabletest->ScanForGVariables(*_AVLListtest);
  //EXPECT_EQ (_CategoryVariabletest, 0);

}


TEST_F(DISABLED__CategoryVariableTest, ScanForVariablesTest) {

  _CategoryVariabletest->ScanForVariables(*_AVLListtest, *booltest);
  //EXPECT_EQ (_CategoryVariabletest, 0);

}


TEST_F(DISABLED__CategoryVariableTest, SerializeCategoryTest) {

  _CategoryVariabletest->SerializeCategory(*_Stringtest);
  //EXPECT_EQ (_CategoryVariabletest, 0);

}


TEST_F(DISABLED__CategoryVariableTest, SetIntervalValueTest) {

  _Parameter result_Parameter = _CategoryVariabletest->SetIntervalValue(*longtest, *booltest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__CategoryVariableTest, makeDynamicTest) {

  BaseRef resultBaseRef = _CategoryVariabletest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__CategoryVariableTest, toStrTest) {

  BaseRef resultBaseRef = _CategoryVariabletest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
