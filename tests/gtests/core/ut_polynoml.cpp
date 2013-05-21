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
#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "polynoml.h"

namespace {

// The fixture for testing class Foo.
class _PolynomialTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _PolynomialTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _MathObjecttest = new _MathObject();
    BaseReftest = new BaseRef();
    _AVLListXtest = new _AVLListX(_SimpleListtest);
    _SimpleListtest = new _SimpleList();
    _Polynomialtest = new _Polynomial(*_Parametertest);
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
    _Parametertest = new _Parameter();
  }

  virtual ~_PolynomialTest() {
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
    delete _MathObjecttest;
    delete BaseReftest;
    delete _AVLListXtest;
    delete _SimpleListtest;
    delete _Polynomialtest;
    delete booltest;
    delete _AVLListtest;
    delete _Parametertest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _MathObject* _MathObjecttest;
  BaseRef* BaseReftest;
  _AVLListX* _AVLListXtest;
  _SimpleList* _SimpleListtest;
  _Polynomial* _Polynomialtest;
  bool* booltest;
  _AVLList* _AVLListtest;
  _Parameter* _Parametertest;
};


TEST_F(_PolynomialTest, AddTest) {

  _MathObject* result_MathObject = _Polynomialtest->Add(_MathObjecttest);
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, CheckTermTest) {

  _Polynomialtest->CheckTerm();
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, ComputeTest) {

  _MathObject* result_MathObject = _Polynomialtest->Compute();
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, ComputePTest) {

  //_Parameter result_Parameter = _Polynomialtest->ComputeP(_Parametertest, _Parametertest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_PolynomialTest, ComputePolynomialTest) {

  _Parameter result_Parameter = _Polynomialtest->ComputePolynomial();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_PolynomialTest, Convert2ComputationFormTest) {

  _Polynomialtest->Convert2ComputationForm(_SimpleListtest, _SimpleListtest);
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, Convert2OperationFormTest) {

  //_Polynomialtest->Convert2OperationForm();
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, DuplicateTest) {

  _Polynomialtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, EqualTest) {

  bool resultbool = _Polynomialtest->Equal(_MathObjecttest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialTest, ExecuteTest) {

  //_PMathObj result_PMathObj = _Polynomialtest->Execute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(_PolynomialTest, HasChangedTest) {

  bool resultbool = _Polynomialtest->HasChanged();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialTest, IsANumberTest) {

  _MathObject* result_MathObject = _Polynomialtest->IsANumber(*booltest);
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, IsMaxElementTest) {

  bool resultbool = _Polynomialtest->IsMaxElement(*_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialTest, IsObjectEmptyTest) {

  bool resultbool = _Polynomialtest->IsObjectEmpty();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialTest, MinusTest) {

  _MathObject* result_MathObject = _Polynomialtest->Minus();
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, MultTest) {

  _MathObject* result_MathObject = _Polynomialtest->Mult(_MathObjecttest);
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, PlusTest) {

  _MathObject* result_MathObject = _Polynomialtest->Plus(_MathObjecttest, *booltest);
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, RaiseTest) {

  _MathObject* result_MathObject = _Polynomialtest->Raise(_MathObjecttest);
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, RankTermsTest) {

  _Polynomialtest->RankTerms(_SimpleListtest);
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, ScanForVariablesTest) {

  _Polynomialtest->ScanForVariables(*_AVLListtest, *booltest, _AVLListXtest);
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, SubTest) {

  _MathObject* result_MathObject = _Polynomialtest->Sub(_MathObjecttest);
  //EXPECT_EQ (result_MathObject*, 0);

}


TEST_F(_PolynomialTest, makeDynamicTest) {

  BaseObj* resultBaseObj = _Polynomialtest->makeDynamic();
  //EXPECT_EQ (resultBaseObj*, 0);

}


TEST_F(_PolynomialTest, toFileStrTest) {

  _Polynomialtest->toFileStr(FILEtest);
  //EXPECT_EQ (_Polynomialtest, 0);

}


TEST_F(_PolynomialTest, toStrTest) {

  BaseObj* resultBaseObj = _Polynomialtest->toStr();
  //EXPECT_EQ (resultBaseObj*, 0);

}


}
