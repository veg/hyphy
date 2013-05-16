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
#include "hy_strings.h"
#include "errorfns.h"
#include "list.h"
#include "simplelist.h"
#include "parser.h"
#include "defines.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

namespace {

// The fixture for testing class Foo.
class _SimpleListTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _SimpleListTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleListtest = new _SimpleList();
    longtest = new long();
    unsignedtest = new unsigned();
    booltest = new bool();
    BaseReftest = new BaseRef();
  }

  virtual ~_SimpleListTest() {
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
    delete unsignedtest;
    delete booltest;
    delete BaseReftest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _SimpleList* _SimpleListtest;
  long* longtest;
  unsigned* unsignedtest;
  bool* booltest;
  BaseRef* BaseReftest;
};


TEST_F(_SimpleListTest, BinaryFindTest) {

  long resultlong = _SimpleListtest->BinaryFind(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, BinaryInsertTest) {

  long resultlong = _SimpleListtest->BinaryInsert(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, BubbleSortTest) {

  _SimpleListtest->BubbleSort();
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, ClearTest) {

  _SimpleListtest->Clear(*booltest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, ClearFormulasInListTest) {

  _SimpleListtest->ClearFormulasInList();
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, CompareTest) {

  long resultlong = _SimpleListtest->Compare(*BaseReftest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, Compare1Test) {

  long resultlong = _SimpleListtest->Compare(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, CountCommonElementsTest) {

  long resultlong = _SimpleListtest->CountCommonElements(*_SimpleListtest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, CountingSortTest) {

  _SimpleList* result_SimpleList = _SimpleListtest->CountingSort(*longtest, _SimpleListtest);
  //EXPECT_EQ (result_SimpleList*, 0);

}


TEST_F(_SimpleListTest, DebugVarListTest) {

  _SimpleListtest->DebugVarList();
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, DeleteTest) {

  _SimpleListtest->Delete(*longtest, *booltest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, DeleteDuplicatesTest) {

  _SimpleListtest->DeleteDuplicates();
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, DeleteListTest) {

  _SimpleListtest->DeleteList(*_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, DisplaceTest) {

  _SimpleListtest->Displace(*longtest, *longtest, *longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, DuplicateTest) {

  _SimpleListtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, ElementTest) {

  long resultlong = _SimpleListtest->Element(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, EqualTest) {

  bool resultbool = _SimpleListtest->Equal(*_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_SimpleListTest, FilterRangeTest) {

  _SimpleListtest->FilterRange(*longtest, *longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, FindTest) {

  long resultlong = _SimpleListtest->Find(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, FindSteppingTest) {

  long resultlong = _SimpleListtest->FindStepping(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, FlipTest) {

  _SimpleListtest->Flip();
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, GetElementTest) {

  long resultlong = _SimpleListtest->GetElement(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, InitializeTest) {

  _SimpleListtest->Initialize(*booltest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, InsertElementTest) {

  _SimpleListtest->InsertElement(*BaseReftest, *longtest, *booltest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, IntersectTest) {

  _SimpleListtest->Intersect(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, ListToPartitionStringTest) {

  BaseRef resultBaseRef = _SimpleListtest->ListToPartitionString();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_SimpleListTest, MaxTest) {

  long resultlong = _SimpleListtest->Max();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, MergeTest) {

  _SimpleListtest->Merge(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, MinTest) {

  long resultlong = _SimpleListtest->Min();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, NChooseKTest) {

  bool resultbool = _SimpleListtest->NChooseK(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_SimpleListTest, NChooseKInitTest) {

  bool resultbool = _SimpleListtest->NChooseKInit(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_SimpleListTest, NormalizeCoordinatesTest) {

  _SimpleListtest->NormalizeCoordinates(*longtest, *longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, OffsetTest) {

  _SimpleListtest->Offset(*longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, PermuteTest) {

  _SimpleListtest->Permute(*longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, PermuteWithReplacementTest) {

  _SimpleListtest->PermuteWithReplacement(*longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, PopTest) {

  long resultlong = _SimpleListtest->Pop();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, PopulateTest) {

  _SimpleListtest->Populate(*longtest, *longtest, *longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, QuickSortTest) {

  _SimpleListtest->QuickSort(*longtest, *longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, RecursiveIndexSortTest) {

  _SimpleListtest->RecursiveIndexSort(*longtest, *longtest, _SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, RequestSpaceTest) {

  _SimpleListtest->RequestSpace(*longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, SortTest) {

  _SimpleListtest->Sort(*booltest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, SubsetTest) {

  _SimpleList* result_SimpleList = _SimpleListtest->Subset(*longtest, *booltest);
  //EXPECT_EQ (result_SimpleList*, 0);

}


TEST_F(_SimpleListTest, SubtractTest) {

  _SimpleListtest->Subtract(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, SumTest) {

  long resultlong = _SimpleListtest->Sum();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, SwapTest) {

  _SimpleListtest->Swap(*longtest, *longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, TrimMemoryTest) {

  _SimpleListtest->TrimMemory();
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, UnionTest) {

  _SimpleListtest->Union(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, XORTest) {

  _SimpleListtest->XOR(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, countitemsTest) {

  long resultlong = _SimpleListtest->countitems();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, makeDynamicTest) {

  BaseRef resultBaseRef = _SimpleListtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_SimpleListTest, operatorAmpTest) {

  _SimpleList result_SimpleList = _SimpleListtest->operatorAmp(*_SimpleListtest);
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(_SimpleListTest, operatorParenthsTest) {

  long resultlong = _SimpleListtest->operatorParenths();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, operatorDoubleLessTest) {

  _SimpleListtest->operatorDoubleLess(*_SimpleListtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, operatorDoubleLess1Test) {

  _SimpleListtest->operatorDoubleLess(*longtest);
  //EXPECT_EQ (_SimpleListtest, 0);

}


TEST_F(_SimpleListTest, operatorEqualTest) {

  _SimpleList result_SimpleList = _SimpleListtest->operatorEqual(*_SimpleListtest);
  //EXPECT_EQ (result_SimpleList, 0);

}


TEST_F(_SimpleListTest, operatorDoubleGreaterTest) {

  bool resultbool = _SimpleListtest->operatorDoubleGreater(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_SimpleListTest, operatorBracketsTest) {

  long resultlong = _SimpleListtest->operatorBrackets(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_SimpleListTest, toStrTest) {

  BaseRef resultBaseRef = _SimpleListtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
