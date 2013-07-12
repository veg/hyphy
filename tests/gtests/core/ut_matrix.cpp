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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include "batchlan.h"
#include "polynoml.h"
#include "likefunc.h"
#include "avllist.h"
#include "avllistx.h"
#include "avllistxl.h"
#include "hy_globals.h"
#include "associativelist.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__MatrixTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__MatrixTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    chartest = new char();
    _Listtest = new _List();
    _Formulatest = new _Formula(*_PMathObjtest, *booltest);
    _MathObjecttest = new _MathObject();
    BaseReftest = new BaseRef();
    _Constanttest = new _Constant(*_Parametertest);
    _SimpleListtest = new _SimpleList();
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
    _Parametertest = new _Parameter();
    _PMathObjtest = new _PMathObj();
    _Matrixtest = new _Matrix();
  }

  virtual ~DISABLED__MatrixTest() {
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
    delete chartest;
    delete _Listtest;
    delete _Formulatest;
    delete _MathObjecttest;
    delete BaseReftest;
    delete _Constanttest;
    delete _SimpleListtest;
    delete booltest;
    delete _AVLListtest;
    delete _Parametertest;
    delete _PMathObjtest;
    delete _Matrixtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  char* chartest;
  _List* _Listtest;
  _Formula* _Formulatest;
  _MathObject* _MathObjecttest;
  BaseRef* BaseReftest;
  _Constant* _Constanttest;
  _SimpleList* _SimpleListtest;
  bool* booltest;
  _AVLList* _AVLListtest;
  _Parameter* _Parametertest;
  _PMathObj* _PMathObjtest;
  _Matrix* _Matrixtest;
};


TEST_F(DISABLED__MatrixTest, AbsTest) {

  _PMathObj result_PMathObj = _Matrixtest->Abs();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, AbsValueTest) {

  _Parameter result_Parameter = _Matrixtest->AbsValue();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, AddTest) {

  //_Matrixtest->Add(*_Matrixtest, *_Matrixtest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, AddObjTest) {

  _PMathObj result_PMathObj = _Matrixtest->AddObj(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, AddWithThresholdTest) {

  //bool resultbool = _Matrixtest->AddWithThreshold(*_Matrixtest, *_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, AgreeObjectsTest) {

  //_Matrixtest->AgreeObjects(*_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, AmISparseTest) {

  bool resultbool = _Matrixtest->AmISparse();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, AmISparseFastTest) {

  //bool resultbool = _Matrixtest->AmISparseFast(*_Matrixtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, AplusBxTest) {

  _Matrixtest->AplusBx(*_Matrixtest, *_Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, BalanceTest) {

  //_Matrixtest->Balance();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, BranchLengthExpressionTest) {

  _String* result_String = _Matrixtest->BranchLengthExpression(_Matrixtest, *booltest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(DISABLED__MatrixTest, CheckCoordinatesTest) {

  bool resultbool = _Matrixtest->CheckCoordinates(*longtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, CheckDimensionsTest) {

  //bool resultbool = _Matrixtest->CheckDimensions(*_Matrixtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, CheckIfSparseEnoughTest) {

  _Matrixtest->CheckIfSparseEnough(*booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, CholeskyDecomposeTest) {

  _PMathObj result_PMathObj = _Matrixtest->CholeskyDecompose();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, ClearTest) {

  _Matrixtest->Clear();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ClearFormulaeTest) {

  //_Matrixtest->ClearFormulae();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ClearObjectsTest) {

  //_Matrixtest->ClearObjects();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, CompareRowsTest) {

  long resultlong = _Matrixtest->CompareRows(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__MatrixTest, CompressSparseMatrixTest) {

  _Matrixtest->CompressSparseMatrix(*booltest, _Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ComputeTest) {

  _PMathObj result_PMathObj = _Matrixtest->Compute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, ComputeNumericTest) {

  _PMathObj result_PMathObj = _Matrixtest->ComputeNumeric(*booltest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, ComputeRowAndColSumsTest) {

  _List* result_List = _Matrixtest->ComputeRowAndColSums();
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(DISABLED__MatrixTest, Convert2FormulasTest) {

  _Matrixtest->Convert2Formulas();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ConvertFormulas2PolyTest) {

  //_Matrixtest->ConvertFormulas2Poly(*booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ConvertNumbers2PolyTest) {

  //_Matrixtest->ConvertNumbers2Poly();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ConvertToSimpleListTest) {

  _Matrixtest->ConvertToSimpleList(*_SimpleListtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, CopyABlockTest) {

  _Matrixtest->CopyABlock(_Matrixtest, *longtest, *longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, DirichletDeviateTest) {

  _PMathObj result_PMathObj = _Matrixtest->DirichletDeviate();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, DuplicateTest) {

  _Matrixtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, EigenDecompTest) {

  //_Matrixtest->EigenDecomp(*_Matrixtest, *_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, EigensystemTest) {

  _PMathObj result_PMathObj = _Matrixtest->Eigensystem();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, EqualTest) {

  bool resultbool = _Matrixtest->Equal(*_PMathObjtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, EvaluateTest) {

  _PMathObj result_PMathObj = _Matrixtest->Evaluate(*booltest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, EvaluateSimpleTest) {

  _PMathObj result_PMathObj = _Matrixtest->EvaluateSimple();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, ExecuteTest) {

  //_PMathObj result_PMathObj = _Matrixtest->Execute();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, ExpNumberOfSubsTest) {

  _Parameter result_Parameter = _Matrixtest->ExpNumberOfSubs(_Matrixtest, *booltest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, ExponentiateTest) {

  _Matrix* result_Matrix = _Matrixtest->Exponentiate();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__MatrixTest, ExportMatrixExpTest) {

  _Matrixtest->ExportMatrixExp(_Matrixtest, FILEtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ExtractElementsByEnumerationTest) {

  _Matrix* result_Matrix = _Matrixtest->ExtractElementsByEnumeration(_SimpleListtest, _SimpleListtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__MatrixTest, FillInListTest) {

  _Matrixtest->FillInList(*_Listtest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, FisherExactTest) {

  _Parameter result_Parameter = _Matrixtest->FisherExact(*_Parametertest, *_Parametertest, *_Parametertest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, GaussianDeviateTest) {

  _PMathObj result_PMathObj = _Matrixtest->GaussianDeviate(*_Matrixtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, GetFormulaTest) {

  _Formula* result_Formula = _Matrixtest->GetFormula(*longtest, *longtest);
  //EXPECT_EQ (result_Formula*, 0);

}


TEST_F(DISABLED__MatrixTest, HasChangedTest) {

  bool resultbool = _Matrixtest->HasChanged();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, HashTest) {

  long resultlong = _Matrixtest->Hash(*longtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__MatrixTest, HashBackTest) {

  //long resultlong = _Matrixtest->HashBack(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__MatrixTest, ImportMatrixExpTest) {

  bool resultbool = _Matrixtest->ImportMatrixExp(FILEtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, IncreaseStorageTest) {

  //bool resultbool = _Matrixtest->IncreaseStorage();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, InitMxVarTest) {

  //_Matrixtest->InitMxVar(*_SimpleListtest, *_Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, InitializeTest) {

  _Matrixtest->Initialize();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, InverseTest) {

  _PMathObj result_PMathObj = _Matrixtest->Inverse();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, InverseWishartDeviateTest) {

  _PMathObj result_PMathObj = _Matrixtest->InverseWishartDeviate(*_Matrixtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, IsAStringMatrixTest) {

  bool resultbool = _Matrixtest->IsAStringMatrix();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, IsAVectorTest) {

  bool resultbool = _Matrixtest->IsAVector(*chartest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, IsConstantTest) {

  bool resultbool = _Matrixtest->IsConstant();
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, IsMaxElementTest) {

  bool resultbool = _Matrixtest->IsMaxElement(*_Parametertest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, IsNonEmptyTest) {

  //inlinebool resultinlinebool = _Matrixtest->IsNonEmpty(*longtest);
  //EXPECT_EQ (resultinlinebool, 0);

}


TEST_F(DISABLED__MatrixTest, IsReversibleTest) {

  bool resultbool = _Matrixtest->IsReversible(_Matrixtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, KDISABLED__MeansTest) {

  _PMathObj result_PMathObj = _Matrixtest->K_Means(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, LUDecomposeTest) {

  _PMathObj result_PMathObj = _Matrixtest->LUDecompose();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, LUSolveTest) {

  _PMathObj result_PMathObj = _Matrixtest->LUSolve(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, LogTest) {

  _PMathObj result_PMathObj = _Matrixtest->Log();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MAccessTest) {

  _PMathObj result_PMathObj = _Matrixtest->MAccess(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MCoordTest) {

  _PMathObj result_PMathObj = _Matrixtest->MCoord(*_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MResolveTest) {

  bool resultbool = _Matrixtest->MResolve(*_PMathObjtest, *_PMathObjtest, *longtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, MStoreTest) {

  _Matrixtest->MStore(*_PMathObjtest, *_PMathObjtest, *_Formulatest, *longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MStore1Test) {

  _Matrixtest->MStore(*_PMathObjtest, *_PMathObjtest, *_PMathObjtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MStore2Test) {

  _Matrixtest->MStore(*longtest, *longtest, *_Formulatest, *longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MStore3Test) {

  _Matrixtest->MStore(*longtest, *longtest, *_PMathObjtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MakeMeGeneralTest) {

  _Matrixtest->MakeMeGeneral();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MakeMeSimpleTest) {

  _Matrixtest->MakeMeSimple();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MakeTreeFromParentTest) {

  _Matrix* result_Matrix = _Matrixtest->MakeTreeFromParent(*longtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__MatrixTest, MaxElementTest) {

  _Parameter result_Parameter = _Matrixtest->MaxElement(*chartest, longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, MaxRelErrorTest) {

  _Parameter result_Parameter = _Matrixtest->MaxRelError(*_Matrixtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, MinElementTest) {

  _Parameter result_Parameter = _Matrixtest->MinElement(*chartest, longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, MultByFreqsTest) {

  _PMathObj result_PMathObj = _Matrixtest->MultByFreqs(*longtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MultElementsTest) {

  _PMathObj result_PMathObj = _Matrixtest->MultElements(*_PMathObjtest, *booltest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MultObjTest) {

  _PMathObj result_PMathObj = _Matrixtest->MultObj(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MultbySTest) {

  //_Matrixtest->MultbyS(*_Matrixtest, *booltest, _Matrixtest, _Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, MultinomialSampleTest) {

  _PMathObj result_PMathObj = _Matrixtest->MultinomialSample(_Constanttest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, MultiplyTest) {

  //_Matrixtest->Multiply(*_Matrixtest, *_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, Multiply1Test) {

  //_Matrixtest->Multiply(*_Matrixtest, *_Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, NeighborJoinTest) {

  _Matrix* result_Matrix = _Matrixtest->NeighborJoin(*booltest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__MatrixTest, NonZeroEntriesTest) {

  _Matrixtest->NonZeroEntries(*_SimpleListtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, PathLogLikelihoodTest) {

  //_PMathObj result_PMathObj = _Matrixtest->PathLogLikelihood(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, PoissonLLTest) {

  _PMathObj result_PMathObj = _Matrixtest->PoissonLL(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, PopulateConstantMatrixTest) {

  _Matrixtest->PopulateConstantMatrix(*_Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ProcessFormulasTest) {

  //bool resultbool = _Matrixtest->ProcessFormulas(*longtest, *_SimpleListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__MatrixTest, ProfileMeanFitTest) {

  //_PMathObj result_PMathObj = _Matrixtest->ProfileMeanFit(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, RandomTest) {

  _PMathObj result_PMathObj = _Matrixtest->Random(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, RecursiveIndexSortTest) {

  //_Matrixtest->RecursiveIndexSort(*longtest, *longtest, _SimpleListtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ResizeTest) {

  _Matrixtest->Resize(*longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, RetrieveNumericTest) {

  _PMathObj result_PMathObj = _Matrixtest->RetrieveNumeric();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, RowAndColumnMaxTest) {

  //_Matrixtest->RowAndColumnMax(*_Parametertest, *_Parametertest, _Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ScanForVariablesTest) {

  _Matrixtest->ScanForVariables(*_AVLListtest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, ScanForVariables2Test) {

  _Matrixtest->ScanForVariables2(*_AVLListtest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SchurTest) {

  //_Matrixtest->Schur();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SerializeTest) {

  _Matrixtest->Serialize(*_Stringtest, *_Stringtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SimplexHelper1Test) {

  //_Matrixtest->SimplexHelper1(*longtest, *_SimpleListtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SimplexHelper2Test) {

  //_Matrixtest->SimplexHelper2(*longtest, *longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SimplexHelper3Test) {

  //_Matrixtest->SimplexHelper3(*longtest, *longtest, *longtest, *longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SimplexSolveTest) {

  _Matrix* result_Matrix = _Matrixtest->SimplexSolve(*_Parametertest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__MatrixTest, SortMatrixOnColumnTest) {

  _PMathObj result_PMathObj = _Matrixtest->SortMatrixOnColumn(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, SqrTest) {

  _Matrixtest->Sqr(_Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, StoreTest) {

  _Matrixtest->Store(*longtest, *longtest, *_Parametertest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, StoreFormulaTest) {

  _Matrixtest->StoreFormula(*longtest, *longtest, *_Formulatest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, StoreObjectTest) {

  _Matrixtest->StoreObject(*longtest, *longtest, _MathObjecttest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, StoreObject1Test) {

  _Matrixtest->StoreObject(*longtest, _MathObjecttest, *booltest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SubObjTest) {

  _PMathObj result_PMathObj = _Matrixtest->SubObj(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, SubtractTest) {

  //_Matrixtest->Subtract(*_Matrixtest, *_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SumTest) {

  _PMathObj result_PMathObj = _Matrixtest->Sum();
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, SwapTest) {

  _Matrixtest->Swap(*_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, SwapRowsTest) {

  _Matrixtest->SwapRows(*longtest, *longtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, TransposeTest) {

  _Matrixtest->Transpose();
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, UpdateDiagTest) {

  _Matrixtest->UpdateDiag(*longtest, *longtest, _MathObjecttest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, WishartDeviateTest) {

  _PMathObj result_PMathObj = _Matrixtest->WishartDeviate(*_Matrixtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, WishartDeviate1Test) {

  _PMathObj result_PMathObj = _Matrixtest->WishartDeviate(*_Matrixtest, *_Matrixtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, branchLengthStencilTest) {

  //_Matrix* result_Matrix = _Matrixtest->branchLengthStencil();
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__MatrixTest, computePFDRTest) {

  //_Parameter result_Parameter = _Matrixtest->computePFDR(*_Parametertest, *_Parametertest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, makeDynamicTest) {

  BaseRef resultBaseRef = _Matrixtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__MatrixTest, operatorParenthsTest) {

  //_Parameter result_Parameter = _Matrixtest->operatorParenths();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, operatorPointerTest) {

  //_Matrix result_Matrix = _Matrixtest->operatorPointer(*_Matrixtest);
  //EXPECT_EQ (result_Matrix, 0);

}


TEST_F(DISABLED__MatrixTest, operatorPointer1Test) {

  //_Matrix result_Matrix = _Matrixtest->operatorPointer(*_Parametertest);
  //EXPECT_EQ (result_Matrix, 0);

}


TEST_F(DISABLED__MatrixTest, operatorMultEqualTest) {

  //_Matrixtest->operatorMultEqual(*_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, operatorUnaryPlusTest) {

  //_Matrix result_Matrix = _Matrixtest->operatorUnaryPlus(*_Matrixtest);
  //EXPECT_EQ (result_Matrix, 0);

}


TEST_F(DISABLED__MatrixTest, operatorPlusEqualTest) {

  //_Matrixtest->operatorPlusEqual(*_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, operatorUnaryNegTest) {

  //_Matrix result_Matrix = _Matrixtest->operatorUnaryNeg(*_Matrixtest);
  //EXPECT_EQ (result_Matrix, 0);

}


TEST_F(DISABLED__MatrixTest, operatorSubEqualTest) {

  //_Matrixtest->operatorSubEqual(*_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, operatorEqualTest) {

  //_Matrixtest->operatorEqual(*_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, operatorEqual1Test) {

  //_Matrixtest->operatorEqual(_Matrixtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, operatorBracketsTest) {

  //_Parameter result_Parameter = _Matrixtest->operatorBrackets(*longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(DISABLED__MatrixTest, pFDRTest) {

  _PMathObj result_PMathObj = _Matrixtest->pFDR(*_PMathObjtest);
  //EXPECT_EQ (result_PMathObj, 0);

}


TEST_F(DISABLED__MatrixTest, toFileStrTest) {

  _Matrixtest->toFileStr(FILEtest);
  //EXPECT_EQ (_Matrixtest, 0);

}


TEST_F(DISABLED__MatrixTest, toStrTest) {

  BaseRef resultBaseRef = _Matrixtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
