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
#include "math.h"
#include "ctype.h"
#include "string.h"
#include "calcnode.h"
#include "calcnode_globals.h"
#include "thetree.h"
#include "scfg.h"
#include "parser.h"
#include "category.h"
#include "batchlan.h"
#include "likefunc.h"
#include "float.h"
#include "hy_globals.h"

namespace {

// The fixture for testing class Foo.
class _CalcNodeTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _CalcNodeTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleListtest = new _SimpleList();
    longtest = new long();
    unsignedtest = new unsigned();
    _Stringtest = new _String(FILEtest);
    _CalcNodetest = new _CalcNode();
    _Listtest = new _List();
    inttest = new int();
    BaseReftest = new BaseRef();
    _Matrixtest = new _Matrix();
    _AVLListXLtest = new _AVLListXL(_SimpleListtest);
    _VariableContainertest = new _VariableContainer();
    booltest = new bool();
//    node<long>test = new node<long>();
  }

  virtual ~_CalcNodeTest() {
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
    delete unsignedtest;
    delete _Stringtest;
    delete _CalcNodetest;
    delete _Listtest;
    delete inttest;
    delete BaseReftest;
    delete _Matrixtest;
    delete _AVLListXLtest;
    delete _VariableContainertest;
    delete booltest;
    delete _SimpleListtest;
//    delete node<long>test;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  unsigned* unsignedtest;
  _SimpleList* _SimpleListtest;
  _String* _Stringtest;
  _CalcNode* _CalcNodetest;
  _List* _Listtest;
  int* inttest;
  BaseRef* BaseReftest;
  _Matrix* _Matrixtest;
  _AVLListXL* _AVLListXLtest;
  _VariableContainer* _VariableContainertest;
  bool* booltest;
////  node<long>* node<long>test;
};


TEST_F(_CalcNodeTest, BranchLengthTest) {

  _Parameter result_Parameter = _CalcNodetest->BranchLength();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_CalcNodeTest, CheckForReferenceNodeTest) {

  long resultlong = _CalcNodetest->CheckForReferenceNode();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_CalcNodeTest, ComputeTest) {

  _PMathObj result_PMathObj = _CalcNodetest->Compute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_CalcNodeTest, ComputeModelMatrixTest) {

  _Matrix* result_Matrix = _CalcNodetest->ComputeModelMatrix(*booltest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(_CalcNodeTest, ConvertFromSimpleMatrixTest) {

  _CalcNodetest->ConvertFromSimpleMatrix();
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, ConvertToSimpleMatrixTest) {

  long resultlong = _CalcNodetest->ConvertToSimpleMatrix();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_CalcNodeTest, DuplicateTest) {

  _CalcNodetest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, FreeUpMemoryTest) {

  long resultlong = _CalcNodetest->FreeUpMemory(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_CalcNodeTest, GetCompExpTest) {

  _Matrix* result_Matrix = _CalcNodetest->GetCompExp(*longtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(_CalcNodeTest, HasChangedTest) {

  bool resultbool = _CalcNodetest->HasChanged();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_CalcNodeTest, InitializeCNTest) {

  _CalcNodetest->InitializeCN(*_Stringtest, *inttest, _VariableContainertest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, LocateMeInTreeTest) {

//  node<long>* resultnode<long> = _CalcNodetest->LocateMeInTree();
//  //EXPECT_EQ (resultnode<long>*, 0);

}


TEST_F(_CalcNodeTest, MatchSubtreeTest) {

  bool resultbool = _CalcNodetest->MatchSubtree(_CalcNodetest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_CalcNodeTest, NeedToExponentiateTest) {

  bool resultbool = _CalcNodetest->NeedToExponentiate(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_CalcNodeTest, ParentTreeTest) {

  _VariableContainer* result_VariableContainer = _CalcNodetest->ParentTree();
  //EXPECT_EQ (result_VariableContainer*, 0);

}


TEST_F(_CalcNodeTest, RecomputeMatrixTest) {

  bool resultbool = _CalcNodetest->RecomputeMatrix(*longtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_CalcNodeTest, RecurseMCTest) {

  //_Formula* result_Formula = _CalcNodetest->RecurseMC(*longtest, node<long>test);
  //EXPECT_EQ (result_Formula*, 0);

}


TEST_F(_CalcNodeTest, RemoveModelTest) {

  _CalcNodetest->RemoveModel();
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, ReplaceModelTest) {

  //TODO
  //_CalcNodetest->ReplaceModel(*_Stringtest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, SetCodeBaseTest) {

  _CalcNodetest->SetCodeBase(*inttest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, SetCompExpTest) {

  _CalcNodetest->SetCompExp(_Matrixtest, *longtest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, SetCompMatrixTest) {

  _CalcNodetest->SetCompMatrix(*longtest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, SetDependanceTest) {

  long resultlong = _CalcNodetest->SetDependance(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_CalcNodeTest, SetModelTest) {

  _CalcNodetest->SetModel(*longtest, _AVLListXLtest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, SetupCategoryMapTest) {
  //TODO
  //_CalcNodetest->SetupCategoryMap(*_Listtest);
  //EXPECT_EQ (_CalcNodetest, 0);

}


TEST_F(_CalcNodeTest, makeDynamicTest) {

  BaseRef resultBaseRef = _CalcNodetest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_CalcNodeTest, operatorBracketsTest) {

  //TODO
  //_Parameter result_Parameter = _CalcNodetest->operatorBrackets(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_CalcNodeTest, toStrTest) {

  BaseRef resultBaseRef = _CalcNodetest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
