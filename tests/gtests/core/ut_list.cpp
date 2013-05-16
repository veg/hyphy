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
#include "list.h"
#include "hy_strings.h"
#include "errorfns.h"
#include "parser.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

namespace {

// The fixture for testing class Foo.
class _ListTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _ListTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleListtest = new _SimpleList();
    longtest = new long();
    booltest = new bool();
    BaseReftest = new BaseRef();
    _Listtest = new _List();
  }

  virtual ~_ListTest() {
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
    delete booltest;
    delete BaseReftest;
    delete _Listtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _SimpleList* _SimpleListtest;
  long* longtest;
  bool* booltest;
  BaseRef* BaseReftest;
  _List* _Listtest;
};


TEST_F(_ListTest, AppendNewInstanceTest) {

  _Listtest->AppendNewInstance(*BaseReftest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, BinaryFindTest) {

  long resultlong = _Listtest->BinaryFind(*BaseReftest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, BinaryInsertTest) {

  long resultlong = _Listtest->BinaryInsert(*BaseReftest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, ClearTest) {

  _Listtest->Clear(*booltest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, CompareTest) {

  long resultlong = _Listtest->Compare(*BaseReftest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, Compare1Test) {

  long resultlong = _Listtest->Compare(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, CountTest) {

  long resultlong = _Listtest->Count();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, DeleteTest) {

  _Listtest->Delete(*longtest, *booltest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, DeleteListTest) {

  _Listtest->DeleteList(*_SimpleListtest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, DuplicateTest) {

  _Listtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, EqualTest) {

  bool resultbool = _Listtest->Equal(*_Listtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ListTest, FindTest) {

  long resultlong = _Listtest->Find(*BaseReftest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, FindStringTest) {

  long resultlong = _Listtest->FindString(*BaseReftest, *longtest, *booltest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, FreeUpMemoryTest) {

  long resultlong = _Listtest->FreeUpMemory(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_ListTest, GetItemTest) {

  BaseRef resultBaseRef = _Listtest->GetItem(*longtest);
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ListTest, InsertElementTest) {

  _Listtest->InsertElement(*BaseReftest, *longtest, *booltest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, IntersectTest) {

  _Listtest->Intersect(*_Listtest, *_Listtest, _SimpleListtest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, JoinTest) {

  BaseRef resultBaseRef = _Listtest->Join(*BaseReftest, *longtest, *longtest);
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ListTest, MapTest) {

  _Listtest->Map(*_Listtest, *_SimpleListtest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, PlaceTest) {

  _Listtest->Place(*BaseReftest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, ReplaceTest) {

  _Listtest->Replace(*longtest, *BaseReftest, *booltest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, bumpNInstTest) {

  _Listtest->bumpNInst();
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, makeDynamicTest) {

  BaseRef resultBaseRef = _Listtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ListTest, operatorAmpTest) {

  _List result_List = _Listtest->operatorAmp(*BaseReftest);
  //EXPECT_EQ (result_List, 0);

}


TEST_F(_ListTest, operatorAmp1Test) {

  _List result_List = _Listtest->operatorAmp(*_Listtest);
  //EXPECT_EQ (result_List, 0);

}


TEST_F(_ListTest, operatorDoubleAmpTest) {

  _Listtest->operatorDoubleAmp(chartest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, operatorParenthsTest) {

  BaseRef resultBaseRef = _Listtest->operatorParenths();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ListTest, operatorDoubleLessTest) {

  _Listtest->operatorDoubleLess(*BaseReftest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, operatorDoubleLess1Test) {

  _Listtest->operatorDoubleLess(*_Listtest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, operatorEqualTest) {

  _List result_List = _Listtest->operatorEqual(*_Listtest);
  //EXPECT_EQ (result_List, 0);

}


TEST_F(_ListTest, operatorDoubleEqualTest) {

  bool resultbool = _Listtest->operatorDoubleEqual(*_Listtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_ListTest, operatorBracketsTest) {

  BaseRef resultBaseRef = _Listtest->operatorBrackets(*longtest);
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_ListTest, toFileStrTest) {

  _Listtest->toFileStr(FILEtest);
  //EXPECT_EQ (_Listtest, 0);

}


TEST_F(_ListTest, toStrTest) {

  BaseRef resultBaseRef = _Listtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
