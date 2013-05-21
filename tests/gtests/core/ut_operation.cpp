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
#include "parser.h"
#include "polynoml.h"
#include "batchlan.h"
#include "hy_globals.h"
#include "executionlist.h"

namespace {

// The fixture for testing class Foo.
class _OperationTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _OperationTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _PMathObjtest = new _MathObject();
    _Stringtest = new _String(FILEtest);
    _Stacktest = new _Stack();
    _Operationtest = new _Operation(_PMathObjtest);
    BaseReftest = new BaseRef();
    _VariableContainertest = new _VariableContainer();
    booltest = new bool();
  }

  virtual ~_OperationTest() {
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
    delete _Stacktest;
    delete _Operationtest;
    delete BaseReftest;
    delete _PMathObjtest;
    delete _VariableContainertest;
    delete booltest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _PMathObj _PMathObjtest;
  _String* _Stringtest;
  _Stack* _Stacktest;
  _Operation* _Operationtest;
  BaseRef* BaseReftest;
  _VariableContainer* _VariableContainertest;
  bool* booltest;
};


TEST_F(_OperationTest, CanResultsBeCachedTest) {

  bool resultbool = _Operationtest->CanResultsBeCached(_Operationtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, DuplicateTest) {

  _Operationtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Operationtest, 0);

}


TEST_F(_OperationTest, EqualOpTest) {

  bool resultbool = _Operationtest->EqualOp(_Operationtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, ExecuteTest) {

  bool resultbool = _Operationtest->Execute(*_Stacktest, _VariableContainertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, ExecutePolynomialTest) {

  bool resultbool = _Operationtest->ExecutePolynomial(*_Stacktest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, HasChangedTest) {

  bool resultbool = _Operationtest->HasChanged();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, InitializeTest) {

  _Operationtest->Initialize();
  //EXPECT_EQ (_Operationtest, 0);

}


TEST_F(_OperationTest, IsAFunctionCallTest) {

  bool resultbool = _Operationtest->IsAFunctionCall();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, IsAVariableTest) {

  bool resultbool = _Operationtest->IsAVariable(*booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, IsConstantTest) {

  bool resultbool = _Operationtest->IsConstant();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, ReportOperationExecutionErrorTest) {

  //bool resultbool = _Operationtest->ReportOperationExecutionError(*_Stringtest, _Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_OperationTest, StackDepthTest) {

  _Operationtest->StackDepth(*longtest);
  //EXPECT_EQ (_Operationtest, 0);

}


TEST_F(_OperationTest, makeDynamicTest) {

  BaseRef resultBaseRef = _Operationtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_OperationTest, toStrTest) {

  BaseRef resultBaseRef = _Operationtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
