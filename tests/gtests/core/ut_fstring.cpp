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
#include "hy_globals.h"
#include "fstring.h"
#include "constant.h"
#include "matrix.h"
#include "calcnode.h"
#include "thetree.h"
#include "growingvector.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__FStringTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__FStringTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _VariableContainertest = new _VariableContainer();
    _hyExecutionContexttest = new _hyExecutionContext(_VariableContainertest, _Stringtest);
    BaseReftest = new BaseRef();
    booltest = new bool();
    _PMathObjtest = new _PMathObj();
    _FStringtest = new _FString(*_Stringtest, *booltest);
  }

  virtual ~DISABLED__FStringTest() {
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
    //delete test;
    delete _Stringtest;
    delete _hyExecutionContexttest;
    delete _VariableContainertest;
    delete BaseReftest;
    delete booltest;
    delete _PMathObjtest;
    delete _FStringtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  _VariableContainer* _VariableContainertest;
  _hyExecutionContext* _hyExecutionContexttest;
  BaseRef* BaseReftest;
  bool* booltest;
  _PMathObj* _PMathObjtest;
  _FString* _FStringtest;
};


TEST_F(DISABLED__FStringTest, AddTest) {

  _PMathObj result_PMathObj = _FStringtest->Add(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, AddOnTest) {

  long resultlong = _FStringtest->AddOn(*_PMathObjtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__FStringTest, AreEqualTest) {

  _PMathObj result_PMathObj = _FStringtest->AreEqual(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, AreEqualCISTest) {

  _PMathObj result_PMathObj = _FStringtest->AreEqualCIS(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, CharAccessTest) {

  _PMathObj result_PMathObj = _FStringtest->CharAccess(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, CountGlobalObjectsTest) {

  _PMathObj result_PMathObj = _FStringtest->CountGlobalObjects();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, DereferenceTest) {

  _PMathObj result_PMathObj = _FStringtest->Dereference(*booltest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, DifferentiateTest) {

  _PMathObj result_PMathObj = _FStringtest->Differentiate(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, DuplicateTest) {

  _FStringtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_FStringtest, 0);

}


TEST_F(DISABLED__FStringTest, EqualAmbTest) {

  _PMathObj result_PMathObj = _FStringtest->EqualAmb(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, EqualRegExpTest) {

  _PMathObj result_PMathObj = _FStringtest->EqualRegExp(*_PMathObjtest, *booltest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, EvaluateTest) {

  _PMathObj result_PMathObj = _FStringtest->Evaluate(_hyExecutionContexttest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, ExecuteTest) {

  //_PMathObj result_PMathObj = _FStringtest->Execute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, FileExistsTest) {

  _PMathObj result_PMathObj = _FStringtest->FileExists();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, GreaterTest) {

  _PMathObj result_PMathObj = _FStringtest->Greater(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, GreaterEqTest) {

  _PMathObj result_PMathObj = _FStringtest->GreaterEq(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, JoinTest) {

  _PMathObj result_PMathObj = _FStringtest->Join(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, LessTest) {

  _PMathObj result_PMathObj = _FStringtest->Less(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, LessEqTest) {

  _PMathObj result_PMathObj = _FStringtest->LessEq(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, MapStringToVectorTest) {

  _PMathObj result_PMathObj = _FStringtest->MapStringToVector(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, NotEqualTest) {

  _PMathObj result_PMathObj = _FStringtest->NotEqual(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, ReplaceReqExpTest) {

  _PMathObj result_PMathObj = _FStringtest->ReplaceReqExp(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, RerootTreeTest) {

  _PMathObj result_PMathObj = _FStringtest->RerootTree();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__FStringTest, makeDynamicTest) {

  BaseRef resultBaseRef = _FStringtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__FStringTest, toStrTest) {

  BaseRef resultBaseRef = _FStringtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
