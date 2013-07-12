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
#include "defines.h"
#include "formula.h"
#include "formulaparsingcontext.h"
#include "parser.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__FormulaTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__FormulaTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleFormulaDatumtest = new _SimpleFormulaDatum();
    longtest = new long();
    _Stringtest = new _String(FILEtest);
    _Listtest = new _List();
    inttest = new int();
    _Parametertest = new _Parameter();
    _Formulatest = new _Formula(*_PMathObjtest, *booltest);
    BaseReftest = new BaseRef();
    _Variabletest = new _Variable(*_Stringtest, *_Stringtest, *booltest);
    _SimpleListtest = new _SimpleList();
    _VariableContainertest = new _VariableContainer();
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
    //node<long>test = new node<long>();
    _PMathObjtest = new _PMathObj();
  }

  virtual ~DISABLED__FormulaTest() {
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
    delete _SimpleFormulaDatumtest;
    delete longtest;
    delete _Stringtest;
    delete _Listtest;
    delete inttest;
    delete _Parametertest;
    delete _Formulatest;
    delete BaseReftest;
    delete _Variabletest;
    delete _SimpleListtest;
    delete _VariableContainertest;
    delete booltest;
    delete _AVLListtest;
    //delete node<long>test;
    delete _PMathObjtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _SimpleFormulaDatum* _SimpleFormulaDatumtest;
  long* longtest;
  _String* _Stringtest;
  _List* _Listtest;
  int* inttest;
  _Parameter* _Parametertest;
  _Formula* _Formulatest;
  BaseRef* BaseReftest;
  _Variable* _Variabletest;
  _SimpleList* _SimpleListtest;
  _VariableContainer* _VariableContainertest;
  bool* booltest;
  _AVLList* _AVLListtest;
  //node<long>* node<long>test;
  _PMathObj* _PMathObjtest;
};


//TEST_F(DISABLED__FormulaTest, AmISimpleTest) {

//  bool resultbool = _Formulatest->AmISimple(*longtest, *_SimpleListtest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, BrentTest) {

//  _Parameter result_Parameter = _Formulatest->Brent(_Variabletest, *_Parametertest, *_Parametertest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, CheckFForDependenceTest) {

//  bool resultbool = _Formulatest->CheckFForDependence(*longtest, *booltest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, CheckSimpleTermTest) {

//  //bool resultbool = _Formulatest->CheckSimpleTerm(*_PMathObjtest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, ClearTest) {

//  _Formulatest->Clear();
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, ComputeTest) {

//  _PMathObj result_PMathObj = _Formulatest->Compute(*longtest, _VariableContainertest);
//  //EXPECT_EQ (result_PMathObj, 0);

//}


//TEST_F(DISABLED__FormulaTest, ComputeSimpleTest) {

//  //_Parameter result_Parameter = _Formulatest->ComputeSimple(_SimpleFormulaDatumtest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, ConstructPolynomialTest) {

//  _PMathObj result_PMathObj = _Formulatest->ConstructPolynomial();
//  //EXPECT_EQ (result_PMathObj, 0);

//}


//TEST_F(DISABLED__FormulaTest, ConvertFromSimpleTest) {

//  _Formulatest->ConvertFromSimple(*_SimpleListtest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, ConvertFromTreeTest) {

//  //_Formulatest->ConvertFromTree();
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, ConvertMatrixArgumentsToSimpleOrComplexFormTest) {

//  _Formulatest->ConvertMatrixArgumentsToSimpleOrComplexForm(*booltest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, ConvertToSimpleTest) {

//  bool resultbool = _Formulatest->ConvertToSimple(*_SimpleListtest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, ConvertToTreeTest) {

//  //_Formulatest->ConvertToTree(*booltest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, DependsOnVariableTest) {

//  bool resultbool = _Formulatest->DependsOnVariable(*longtest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, DereferenceTest) {

//  _Variable* result_Variable = _Formulatest->Dereference(*booltest);
//  //EXPECT_EQ (result_Variable*, 0);

//}


//TEST_F(DISABLED__FormulaTest, DifferentiateTest) {

//  _Formula* result_Formula = _Formulatest->Differentiate(*_Stringtest, *booltest);
//  //EXPECT_EQ (result_Formula*, 0);

//}


//TEST_F(DISABLED__FormulaTest, DuplicateTest) {

//  _Formulatest->Duplicate(*BaseReftest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, DuplicateFormulaTest) {

//  //node<long>* resultnode<long> = _Formulatest->DuplicateFormula(node<long>test, *_Formulatest);
//  //EXPECT_EQ (resultnode<long>*, 0);

//}


//TEST_F(DISABLED__FormulaTest, DuplicateReferenceTest) {

//  _Formulatest->DuplicateReference(_Formulatest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, EqualFormulaTest) {

//  bool resultbool = _Formulatest->EqualFormula(_Formulatest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, ExtractMatrixExpArgumentsTest) {

//  long resultlong = _Formulatest->ExtractMatrixExpArguments(_Listtest);
//  //EXPECT_EQ (resultlong, 0);

//}


//TEST_F(DISABLED__FormulaTest, GetIthTermTest) {

//  _Operation* result_Operation = _Formulatest->GetIthTerm(*longtest);
//  //EXPECT_EQ (result_Operation*, 0);

//}


//TEST_F(DISABLED__FormulaTest, GetTheMatrixTest) {

//  _PMathObj result_PMathObj = _Formulatest->GetTheMatrix();
//  //EXPECT_EQ (result_PMathObj, 0);

//}


//TEST_F(DISABLED__FormulaTest, HasChangedTest) {

//  bool resultbool = _Formulatest->HasChanged(*booltest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, HasChangedSimpleTest) {

//  bool resultbool = _Formulatest->HasChangedSimple(*_SimpleListtest);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, InitializeTest) {

//  _Formulatest->Initialize();
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, IntegralTest) {

//  _Parameter result_Parameter = _Formulatest->Integral(_Variabletest, *_Parametertest, *_Parametertest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, InternalDifferentiateTest) {

//  //node<long>* resultnode<long> = _Formulatest->InternalDifferentiate(node<long>test);
//  //EXPECT_EQ (resultnode<long>*, 0);

//}


//TEST_F(DISABLED__FormulaTest, InternalSimplifyTest) {

//  //bool resultbool = _Formulatest->InternalSimplify(node<long>test);
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, IsAConstantTest) {

//  bool resultbool = _Formulatest->IsAConstant();
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, IsArrayAccessTest) {

//  bool resultbool = _Formulatest->IsArrayAccess();
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, IsConstantTest) {

//  bool resultbool = _Formulatest->IsConstant();
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, IsEmptyTest) {

//  bool resultbool = _Formulatest->IsEmpty();
//  //EXPECT_EQ (resultbool, 0);

//}


//TEST_F(DISABLED__FormulaTest, LocalizeFormulaTest) {

//  //_Formulatest->LocalizeFormula(*_Formulatest, *_Stringtest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, MeanIntegralTest) {

//  //_Parameter result_Parameter = _Formulatest->MeanIntegral(_Variabletest, *_Parametertest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, NewtonTest) {

//  //_Parameter result_Parameter = _Formulatest->Newton(*_Formulatest, *_Parametertest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, Newton1Test) {

//  //_Parameter result_Parameter = _Formulatest->Newton(*_Formulatest, _Variabletest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, Newton2Test) {

//  //_Parameter result_Parameter = _Formulatest->Newton(_Variabletest, *_Parametertest);
//  //EXPECT_EQ (result_Parameter, 0);

//}


//TEST_F(DISABLED__FormulaTest, NumberOperationsTest) {

//  long resultlong = _Formulatest->NumberOperations();
//  //EXPECT_EQ (resultlong, 0);

//}


//TEST_F(DISABLED__FormulaTest, ObjectClassTest) {

//  long resultlong = _Formulatest->ObjectClass();
//  //EXPECT_EQ (resultlong, 0);

//}


//TEST_F(DISABLED__FormulaTest, PatchFormulasTogetherTest) {

//  //_Formula result_Formula = _Formulatest->PatchFormulasTogether(*_Formulatest);
//  //EXPECT_EQ (result_Formula, 0);

//}


//TEST_F(DISABLED__FormulaTest, ScanFForTypeTest) {

//  _Formulatest->ScanFForType(*_SimpleListtest, *inttest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, ScanFForVariablesTest) {

//  _Formulatest->ScanFForVariables(*_AVLListtest, *booltest);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, SimplifyConstantsTest) {

//  _Formulatest->SimplifyConstants();
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, internalToStrTest) {

//  //_Formulatest->internalToStr(*_Stringtest, node<long>test);
//  //EXPECT_EQ (_Formulatest, 0);

//}


//TEST_F(DISABLED__FormulaTest, makeDynamicTest) {

//  BaseRef resultBaseRef = _Formulatest->makeDynamic();
//  //EXPECT_EQ (resultBaseRef, 0);

//}


//TEST_F(DISABLED__FormulaTest, operatorPointerTest) {

//  //_Formula result_Formula = _Formulatest->operatorPointer(*_Formulatest);
//  //EXPECT_EQ (result_Formula, 0);

//}


//TEST_F(DISABLED__FormulaTest, operatorUnaryPlusTest) {

//  //_Formula result_Formula = _Formulatest->operatorUnaryPlus(*_Formulatest);
//  //EXPECT_EQ (result_Formula, 0);

//}


//TEST_F(DISABLED__FormulaTest, operatorUnaryNegTest) {

//  //_Formula result_Formula = _Formulatest->operatorUnaryNeg(*_Formulatest);
//  //EXPECT_EQ (result_Formula, 0);

//}


//TEST_F(DISABLED__FormulaTest, operatorDivTest) {

//  //_Formula result_Formula = _Formulatest->operatorDiv(*_Formulatest);
//  //EXPECT_EQ (result_Formula, 0);

//}


//TEST_F(DISABLED__FormulaTest, operatorCarotTest) {

//  //_Formula result_Formula = _Formulatest->operatorCarot(*_Formulatest);
//  //EXPECT_EQ (result_Formula, 0);

//}


//TEST_F(DISABLED__FormulaTest, toStrTest) {

//  BaseRef resultBaseRef = _Formulatest->toStr(_Listtest, *booltest);
//  //EXPECT_EQ (resultBaseRef, 0);

//}


}
