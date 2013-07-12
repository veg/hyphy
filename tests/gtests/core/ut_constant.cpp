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
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include "string.h"
#include "stdlib.h"
#include "time.h"
#include "constant.h"
#include "parser.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__ConstantTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__ConstantTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");
    _Parametertest = new _Parameter(-4.5);
    BaseReftest = new BaseRef();
    _Constanttest = new _Constant(*_Parametertest);
    _PMathObjtest = new _PMathObj();
  }

  virtual ~DISABLED__ConstantTest() {
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
    delete BaseReftest;
    delete _Parametertest;
    delete _Constanttest;
    delete _PMathObjtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  BaseRef* BaseReftest;
  _Parameter* _Parametertest;
  _Constant* _Constanttest;
  _PMathObj* _PMathObjtest;
};


TEST_F(DISABLED__ConstantTest, AbsTest) {

  _PMathObj result_PMathObj = _Constanttest->Abs();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, AddTest) {

  _PMathObj result_PMathObj = _Constanttest->Add(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, ArctanTest) {

  _PMathObj result_PMathObj = _Constanttest->Arctan();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, AreEqualTest) {

  _PMathObj result_PMathObj = _Constanttest->AreEqual(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, AssignTest) {

  _Constanttest->Assign(*_PMathObjtest);
  //EXPECT_EQ (_Constanttest, 0);

}


TEST_F(DISABLED__ConstantTest, BetaTest) {

  _PMathObj result_PMathObj = _Constanttest->Beta(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, CChi2Test) {

  _PMathObj result_PMathObj = _Constanttest->CChi2(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, CGammaDistTest) {

  _PMathObj result_PMathObj = _Constanttest->CGammaDist(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, CosTest) {

  _PMathObj result_PMathObj = _Constanttest->Cos();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, DivTest) {

  _PMathObj result_PMathObj = _Constanttest->Div(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, DuplicateTest) {

  _Constanttest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Constanttest, 0);

}


TEST_F(DISABLED__ConstantTest, EqualTest) {

  bool resultbool = _Constanttest->Equal(*_PMathObjtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__ConstantTest, ErfTest) {

  _PMathObj result_PMathObj = _Constanttest->Erf();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, ExpTest) {

  _PMathObj result_PMathObj = _Constanttest->Exp();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, FormatNumberStringTest) {

  _PMathObj result_PMathObj = _Constanttest->FormatNumberString(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, GammaTest) {

  _PMathObj result_PMathObj = _Constanttest->Gamma();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, GammaDistTest) {

  _PMathObj result_PMathObj = _Constanttest->GammaDist(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, GreaterTest) {

  _PMathObj result_PMathObj = _Constanttest->Greater(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, GreaterEqTest) {

  _PMathObj result_PMathObj = _Constanttest->GreaterEq(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, IBetaTest) {

  _PMathObj result_PMathObj = _Constanttest->IBeta(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, IGammaTest) {

  _PMathObj result_PMathObj = _Constanttest->IGamma(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, InitializeTest) {

  _Constanttest->Initialize();
  //EXPECT_EQ (_Constanttest, 0);

}


TEST_F(DISABLED__ConstantTest, InvChi2Test) {

  _PMathObj result_PMathObj = _Constanttest->InvChi2(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LAndTest) {

  _PMathObj result_PMathObj = _Constanttest->LAnd(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LNotTest) {

  _PMathObj result_PMathObj = _Constanttest->LNot();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LOrTest) {

  _PMathObj result_PMathObj = _Constanttest->LOr(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LessTest) {

  _PMathObj result_PMathObj = _Constanttest->Less(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LessEqTest) {

  _PMathObj result_PMathObj = _Constanttest->LessEq(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LnGammaTest) {

  _PMathObj result_PMathObj = _Constanttest->LnGamma();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, LogTest) {

  _PMathObj result_PMathObj = _Constanttest->Log();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, MaxTest) {

  _PMathObj result_PMathObj = _Constanttest->Max(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, MinTest) {

  _PMathObj result_PMathObj = _Constanttest->Min(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, MinusTest) {

  _PMathObj result_PMathObj = _Constanttest->Minus();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, MultTest) {

  _PMathObj result_PMathObj = _Constanttest->Mult(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, NotEqualTest) {

  _PMathObj result_PMathObj = _Constanttest->NotEqual(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, RaiseTest) {

  _PMathObj result_PMathObj = _Constanttest->Raise(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, RandomTest) {

  _PMathObj result_PMathObj = _Constanttest->Random(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, SinTest) {

  _PMathObj result_PMathObj = _Constanttest->Sin();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, SqrtTest) {

  _PMathObj result_PMathObj = _Constanttest->Sqrt();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, SubTest) {

  _PMathObj result_PMathObj = _Constanttest->Sub(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, SumTest) {

  _PMathObj result_PMathObj = _Constanttest->Sum();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, TanTest) {

  _PMathObj result_PMathObj = _Constanttest->Tan();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, TimeTest) {

  _PMathObj result_PMathObj = _Constanttest->Time();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, ValueTest) {

  _Parameter result_Parameter = _Constanttest->Value();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__ConstantTest, ZCDFTest) {

  _PMathObj result_PMathObj = _Constanttest->ZCDF();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, lDivTest) {

  _PMathObj result_PMathObj = _Constanttest->lDiv(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, longDivTest) {

  _PMathObj result_PMathObj = _Constanttest->longDiv(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__ConstantTest, makeDynamicTest) {

  BaseRef resultBaseRef = _Constanttest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__ConstantTest, toStrTest) {

  BaseRef resultBaseRef = _Constanttest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
