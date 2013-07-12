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
#include "avllist.h"
#include "hy_strings.h"
#include "errorfns.h"
#include "parser.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

namespace {

// The fixture for testing class Foo.
class DISABLED__AVLListTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__AVLListTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleListtest = new _SimpleList();
    longtest = new long();
    booltest = new bool();
    BaseReftest = new BaseRef();
    _AVLListtest = new _AVLList(_SimpleListtest);
  }

  virtual ~DISABLED__AVLListTest() {
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
    delete _AVLListtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _SimpleList* _SimpleListtest;
  long* longtest;
  bool* booltest;
  BaseRef* BaseReftest;
  _AVLList* _AVLListtest;
};


TEST_F(DISABLED__AVLListTest, ClearTest) {

  _AVLListtest->Clear(*booltest);
  //EXPECT_EQ (_AVLListtest, 0);

}


TEST_F(DISABLED__AVLListTest, ConsistencyCheckTest) {

  _AVLListtest->ConsistencyCheck();
  //EXPECT_EQ (_AVLListtest, 0);

}


TEST_F(DISABLED__AVLListTest, DeleteTest) {

  _AVLListtest->Delete(*BaseReftest, *booltest);
  //EXPECT_EQ (_AVLListtest, 0);

}


TEST_F(DISABLED__AVLListTest, FindTest) {

  long resultlong = _AVLListtest->Find(*BaseReftest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, Find1Test) {

  long resultlong = _AVLListtest->Find(*BaseReftest, *_SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, FindBestTest) {

  char resultchar = _AVLListtest->FindBest(*BaseReftest, *longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(DISABLED__AVLListTest, FindLongTest) {

  long resultlong = _AVLListtest->FindLong(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, FirstTest) {

  long resultlong = _AVLListtest->First();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, GetByIndexTest) {

  long resultlong = _AVLListtest->GetByIndex(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, HasDataTest) {

  bool resultbool = _AVLListtest->HasData(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__AVLListTest, InsertTest) {

  long resultlong = _AVLListtest->Insert(*BaseReftest, *longtest, *booltest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, InsertDataTest) {

  long resultlong = _AVLListtest->InsertData(*BaseReftest, *longtest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, LastTest) {

  long resultlong = _AVLListtest->Last();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, NextTest) {

  long resultlong = _AVLListtest->Next(*longtest, *_SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, PrevTest) {

  long resultlong = _AVLListtest->Prev(*longtest, *_SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, ReorderListTest) {

  _AVLListtest->ReorderList(_SimpleListtest);
  //EXPECT_EQ (_AVLListtest, 0);

}


TEST_F(DISABLED__AVLListTest, RetrieveTest) {

  BaseRef resultBaseRef = _AVLListtest->Retrieve(*longtest);
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__AVLListTest, TraverserTest) {

  long resultlong = _AVLListtest->Traverser(*_SimpleListtest, *longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, countitemsTest) {

  long resultlong = _AVLListtest->countitems();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__AVLListTest, toStrTest) {

  BaseRef resultBaseRef = _AVLListtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
