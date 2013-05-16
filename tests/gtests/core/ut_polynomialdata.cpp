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
#include "defines.h"
#include "math.h"
#include "polynoml.h"
#include "polynomialdata.h"

namespace {

// The fixture for testing class Foo.
class _PolynomialDataTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _PolynomialDataTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    BaseReftest = new BaseRef();
    _Parametertest = new _Parameter();
    _PolynomialDatatest = new _PolynomialData(*_PolynomialDatatest);
  }

  virtual ~_PolynomialDataTest() {
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
    delete BaseReftest;
    delete _Parametertest;
    delete _PolynomialDatatest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  BaseRef* BaseReftest;
  _Parameter* _Parametertest;
  _PolynomialData* _PolynomialDatatest;
};


TEST_F(_PolynomialDataTest, AddTermTest) {

  _PolynomialDatatest->AddTerm(*_Parametertest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, AddTerm1Test) {

  _PolynomialDatatest->AddTerm(longtest, *_Parametertest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, AddTerm2Test) {

  _PolynomialDatatest->AddTerm(longtest, *_Parametertest, longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, BinaryRaiseTest) {

  _Parameter result_Parameter = _PolynomialDatatest->BinaryRaise(*_Parametertest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_PolynomialDataTest, ChopTermsTest) {

  _PolynomialDatatest->ChopTerms();
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, CompareTermsTest) {

  char resultchar = _PolynomialDatatest->CompareTerms(longtest, longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_PolynomialDataTest, CompareTerms1Test) {

  char resultchar = _PolynomialDatatest->CompareTerms(longtest, longtest, longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_PolynomialDataTest, CompareTerms2Test) {

  char resultchar = _PolynomialDatatest->CompareTerms(longtest, longtest, longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_PolynomialDataTest, DeleteTermTest) {

  _PolynomialDatatest->DeleteTerm(*longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, DuplicateTest) {

  _PolynomialDatatest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, FindTermTest) {

  long resultlong = _PolynomialDatatest->FindTerm(longtest, longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_PolynomialDataTest, GetTermTest) {

  long* resultlong = _PolynomialDatatest->GetTerm(*longtest);
  //EXPECT_EQ (resultlong*, 0);

}


TEST_F(_PolynomialDataTest, IsFirstANumberTest) {

  bool resultbool = _PolynomialDatatest->IsFirstANumber();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialDataTest, MultiplyTermsTest) {

  _PolynomialDatatest->MultiplyTerms(longtest, longtest, longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, RaiseTermTest) {

  _PolynomialDatatest->RaiseTerm(longtest, *longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, RearrangeTermTest) {

  _PolynomialDatatest->RearrangeTerm(longtest, longtest, longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, ResortTermsTest) {

  _PolynomialDatatest->ResortTerms(longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, SumOfPowersTest) {

  long resultlong = _PolynomialDatatest->SumOfPowers(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_PolynomialDataTest, WeightedSumOfPowersTest) {

  long resultlong = _PolynomialDatatest->WeightedSumOfPowers(*longtest, _Parametertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_PolynomialDataTest, WriteTermTest) {

  _PolynomialDatatest->WriteTerm(longtest, *longtest);
  //EXPECT_EQ (_PolynomialDatatest, 0);

}


TEST_F(_PolynomialDataTest, checkMeTest) {

  bool resultbool = _PolynomialDatatest->checkMe();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialDataTest, checkTermTest) {

  bool resultbool = _PolynomialDatatest->checkTerm(*_Parametertest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_PolynomialDataTest, makeDynamicTest) {

  BaseRef resultBaseRef = _PolynomialDatatest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
