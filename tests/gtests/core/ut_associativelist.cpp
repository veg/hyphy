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
#include "associativelist.h"
#include "batchlan.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__AssociativeListTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__AssociativeListTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _AssociativeListtest = new _AssociativeList();
    _Listtest = new _List();
    BaseReftest = new BaseRef();
    booltest = new bool();
    _PMathObjtest = new _PMathObj();
  }

  virtual ~DISABLED__AssociativeListTest() {
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
    delete _AssociativeListtest;
    delete _Listtest;
    delete BaseReftest;
    delete booltest;
    delete _PMathObjtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  _AssociativeList* _AssociativeListtest;
  _List* _Listtest;
  BaseRef* BaseReftest;
  bool* booltest;
  _PMathObj* _PMathObjtest;
};


TEST_F(DISABLED__AssociativeListTest, ComputeTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->Compute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, DeleteByKeyTest) {

  _AssociativeListtest->DeleteByKey(*_PMathObjtest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, DuplicateTest) {

  _AssociativeListtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, ExecuteTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->Execute(*longtest, *_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, FillInListTest) {

  _AssociativeListtest->FillInList(*_Listtest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, GetByKeyTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->GetByKey(*_Stringtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, GetByKey1Test) {

  _PMathObj result_PMathObj = _AssociativeListtest->GetByKey(*_Stringtest, *longtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, GetByKey2Test) {

  _PMathObj result_PMathObj = _AssociativeListtest->GetByKey(*longtest, *longtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, GetKeysTest) {

  _List* result_List = _AssociativeListtest->GetKeys();
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(DISABLED__AssociativeListTest, MAccessTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->MAccess(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, MCoordTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->MCoord(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, MIteratorTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->MIterator(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, MStoreTest) {

  _AssociativeListtest->MStore(*_PMathObjtest, *_PMathObjtest, *booltest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, MStore1Test) {

  _AssociativeListtest->MStore(*_Stringtest, *_PMathObjtest, *booltest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, MStore2Test) {

  _AssociativeListtest->MStore(*_Stringtest, *_Stringtest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, MergeTest) {

  _AssociativeListtest->Merge(*_PMathObjtest);
  //EXPECT_EQ (_AssociativeListtest, 0);

}


TEST_F(DISABLED__AssociativeListTest, ParseStringRepresentationTest) {

  bool resultbool = _AssociativeListtest->ParseStringRepresentation(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__AssociativeListTest, SerializeTest) {

  _String* result_String = _AssociativeListtest->Serialize(*_Stringtest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(DISABLED__AssociativeListTest, SumTest) {

  _PMathObj result_PMathObj = _AssociativeListtest->Sum();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__AssociativeListTest, makeDynamicTest) {

  BaseRef resultBaseRef = _AssociativeListtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__AssociativeListTest, toStrTest) {

  BaseRef resultBaseRef = _AssociativeListtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
