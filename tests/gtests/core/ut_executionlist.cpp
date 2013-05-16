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
#include "executionlist.h"
#include "likefunc.h"
#include "batchlan.h"
#include "string.h"
#include "ctype.h"
#include "polynoml.h"
#include "time.h"
#include "scfg.h"
#include "HYNetInterface.h"
#include "hy_globals.h"
#include "batchlan_globals.h"
#include "bayesgraph.h"
#include "hyphy_qt_helpers.h"

namespace {

// The fixture for testing class Foo.
class _ExecutionListTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _ExecutionListTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _ExecutionListtest = new _ExecutionList();
    BaseReftest = new BaseRef();
    _SimpleListtest = new _SimpleList();
    booltest = new bool();
  }

  virtual ~_ExecutionListTest() {
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
    delete _ExecutionListtest;
    delete BaseReftest;
    delete _SimpleListtest;
    delete booltest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  _ExecutionList* _ExecutionListtest;
  BaseRef* BaseReftest;
  _SimpleList* _SimpleListtest;
  bool* booltest;
};


TEST_F(_ExecutionListTest, AddNameSpaceToIDTest) {

  _String result_String = _ExecutionListtest->AddNameSpaceToID(*_Stringtest, _Stringtest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_ExecutionListTest, BuildListTest) {

  bool resultbool = _ExecutionListtest->BuildList(*_Stringtest, _SimpleListtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ExecutionListTest, DuplicateTest) {

  _ExecutionListtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_ExecutionListtest, 0);

}


TEST_F(_ExecutionListTest, ExecuteTest) {

  _PMathObj result_PMathObj = _ExecutionListtest->Execute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_ExecutionListTest, ExecuteAndCleanTest) {

  long resultlong = _ExecutionListtest->ExecuteAndClean(*longtest, _Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ExecutionListTest, ExecuteSimpleTest) {

  _ExecutionListtest->ExecuteSimple();
  //EXPECT_EQ (_ExecutionListtest, 0);

}


TEST_F(_ExecutionListTest, FetchFromStdinRedirectTest) {

  _String* result_String = _ExecutionListtest->FetchFromStdinRedirect();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_ExecutionListTest, GetFileNameTest) {

  _String result_String = _ExecutionListtest->GetFileName();
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_ExecutionListTest, GetNameSpaceTest) {

  _String* result_String = _ExecutionListtest->GetNameSpace();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_ExecutionListTest, ReportAnExecutionErrorTest) {

  _ExecutionListtest->ReportAnExecutionError(*_Stringtest);
  //EXPECT_EQ (_ExecutionListtest, 0);

}


TEST_F(_ExecutionListTest, ResetFormulaeTest) {

  _ExecutionListtest->ResetFormulae();
  //EXPECT_EQ (_ExecutionListtest, 0);

}


TEST_F(_ExecutionListTest, ResetNameSpaceTest) {

  _ExecutionListtest->ResetNameSpace();
  //EXPECT_EQ (_ExecutionListtest, 0);

}


TEST_F(_ExecutionListTest, SetNameSpaceTest) {

  _ExecutionListtest->SetNameSpace(*_Stringtest);
  //EXPECT_EQ (_ExecutionListtest, 0);

}


TEST_F(_ExecutionListTest, TrimNameSpaceFromIDTest) {

  _String result_String = _ExecutionListtest->TrimNameSpaceFromID(*_Stringtest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_ExecutionListTest, TryToMakeSimpleTest) {

  bool resultbool = _ExecutionListtest->TryToMakeSimple();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ExecutionListTest, makeDynamicTest) {

  BaseRef resultBaseRef = _ExecutionListtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ExecutionListTest, toStrTest) {

  BaseRef resultBaseRef = _ExecutionListtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
