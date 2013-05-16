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
#include "site.h"
#include "dataset.h"
#include "batchlan.h"
#include "ctype.h"

namespace {

// The fixture for testing class Foo.
class _DataSetTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _DataSetTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    _Stringtest = new _String(FILEtest);
    chartest = new char();
    _PMathObjtest = new _PMathObj();
    _DataSettest = new _DataSet(FILEtest);
    _Formulatest = new _Formula(*_PMathObjtest, *booltest);
    _TranslationTabletest = new _TranslationTable(*_TranslationTabletest);
    _SimpleListtest = new _SimpleList();
    booltest = new bool();
    _Parametertest = new _Parameter();
  }

  virtual ~_DataSetTest() {
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
    delete _DataSettest;
    delete _PMathObjtest;
    delete _Formulatest;
    delete _TranslationTabletest;
    delete _SimpleListtest;
    delete booltest;
    delete _Parametertest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  _String* _Stringtest;
  char* chartest;
  _DataSet* _DataSettest;
  _PMathObj* _PMathObjtest;
  _Formula* _Formulatest;
  _TranslationTable* _TranslationTabletest;
  _SimpleList* _SimpleListtest;
  bool* booltest;
  _Parameter* _Parametertest;
};


TEST_F(_DataSetTest, AddNameTest) {

  _DataSettest->AddName(*_Stringtest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, AddSiteTest) {

  _DataSettest->AddSite(*chartest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, CheckAlphabetConsistencyTest) {

  _Parameter result_Parameter = _DataSettest->CheckAlphabetConsistency();
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_DataSetTest, CheckCompatibilityTest) {

  //TODO:
  //_TranslationTable* result_TranslationTable = _DataSettest->CheckCompatibility(*_SimpleListtest);
  //EXPECT_EQ (result_TranslationTable*, 0);

}


TEST_F(_DataSetTest, CheckMappingTest) {

  _DataSettest->CheckMapping(*longtest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, ClearTest) {

  _DataSettest->Clear();
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, CombineTest) {

  _DataSet* result_DataSet = _DataSettest->Combine(*_SimpleListtest);
  //EXPECT_EQ (result_DataSet*, 0);

}


TEST_F(_DataSetTest, CompactTest) {

  _DataSettest->Compact(*longtest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, ComputeSizeTest) {

  long resultlong = _DataSettest->ComputeSize();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_DataSetTest, ConcatenateTest) {

  _DataSet* result_DataSet = _DataSettest->Concatenate(*_SimpleListtest);
  //EXPECT_EQ (result_DataSet*, 0);

}


TEST_F(_DataSetTest, ConvertRepresentationsTest) {

  _DataSettest->ConvertRepresentations();
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, FinalizeTest) {

  _DataSettest->Finalize();
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, FindAllSitesLikeThisOneTest) {

  _DataSettest->FindAllSitesLikeThisOne(*longtest, *_SimpleListtest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, GetCharDimensionTest) {

  long resultlong = _DataSettest->GetCharDimension();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_DataSetTest, GetFreqTypeTest) {

  long resultlong = _DataSettest->GetFreqType(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_DataSetTest, GetNoTypesTest) {

  long resultlong = _DataSettest->GetNoTypes();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_DataSetTest, HarvestFrequenciesTest) {
  //TODO
  //_Matrix* result_Matrix = _DataSettest->HarvestFrequencies(*chartest, *chartest, *booltest);
  //EXPECT_EQ (result_Matrix*, 0);
}


TEST_F(_DataSetTest, MatchIndicesTest) {

  //TODO
  //_DataSettest->MatchIndices(*_Formulatest, *_SimpleListtest, *booltest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, ProcessPartitionTest) {

  //TODO
  //_DataSettest->ProcessPartition(*_Stringtest, *_SimpleListtest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, ResetIHelperTest) {

  _DataSettest->ResetIHelper();
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, SetTranslationTableTest) {

  _DataSettest->SetTranslationTable(_DataSettest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, SetTranslationTable1Test) {

  _DataSettest->SetTranslationTable(_TranslationTabletest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, Write2SiteTest) {

  _DataSettest->Write2Site(*longtest, *chartest);
  //EXPECT_EQ (_DataSettest, 0);

}

TEST_F(_DataSetTest, makeDynamicTest) {

  BaseRef resultBaseRef = _DataSettest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_DataSetTest, operatorParenthsTest) {

  //TODO
  //char resultchar = _DataSettest->operatorParenths();
  //EXPECT_EQ (resultinlinechar, 0);

}


TEST_F(_DataSetTest, toFileStrTest) {

  _DataSettest->toFileStr(FILEtest);
  //EXPECT_EQ (_DataSettest, 0);

}


TEST_F(_DataSetTest, toStrTest) {

  BaseRef resultBaseRef = _DataSettest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
