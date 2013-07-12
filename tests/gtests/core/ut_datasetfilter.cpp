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
#include "datasetfilter.h"
#include "hy_globals.h"
#include "likefunc.h"
#include "site.h"

namespace {

// The fixture for testing class Foo.
class DISABLED__DataSetFilterTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  DISABLED__DataSetFilterTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    unsignedtest = new unsigned();
    _Stringtest = new _String(FILEtest);
    chartest = new char();
    _Listtest = new _List();
    _DataSettest = new _DataSet(FILEtest);
    _Matrixtest = new _Matrix();
    _SimpleListtest = new _SimpleList();
    _DataSetFiltertest = new _DataSetFilter(_DataSettest, *chartest, *_Stringtest);
    booltest = new bool();
    _AVLListtest = new _AVLList(_SimpleListtest);
    _Parametertest = new _Parameter();
  }

  virtual ~DISABLED__DataSetFilterTest() {
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
    delete chartest;
    delete _Listtest;
    delete _DataSettest;
    delete _Matrixtest;
    delete _SimpleListtest;
    delete _DataSetFiltertest;
    delete booltest;
    delete _AVLListtest;
    delete _Parametertest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  unsigned* unsignedtest;
  _String* _Stringtest;
  char* chartest;
  _List* _Listtest;
  _DataSet* _DataSettest;
  _Matrix* _Matrixtest;
  _SimpleList* _SimpleListtest;
  _DataSetFilter* _DataSetFiltertest;
  bool* booltest;
  _AVLList* _AVLListtest;
  _Parameter* _Parametertest;
};


TEST_F(DISABLED__DataSetFilterTest, CompareTwoSitesTest) {

  //TODO
  //bool resultbool = _DataSetFiltertest->CompareTwoSites(*longtest, *longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__DataSetFilterTest, CompareTwoSitesCharTest) {

  //TODO
  //bool resultbool = _DataSetFiltertest->CompareTwoSitesChar(*longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__DataSetFilterTest, ComputePairwiseDifferencesTest) {

  _Matrix* result_Matrix = _DataSetFiltertest->ComputePairwiseDifferences(*longtest, *longtest, *chartest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, ComputePairwiseDifferences1Test) {

  _DataSetFiltertest->ComputePairwiseDifferences(*_Matrixtest, *longtest, *longtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, ComputePatternToSiteMapTest) {

  _List* result_List = _DataSetFiltertest->ComputePatternToSiteMap();
  //EXPECT_EQ (result_List*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, ConvertCodeToLettersBufferedTest) {
  //TODO
  //_DataSetFiltertest->ConvertCodeToLettersBuffered(*longtest, *chartest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, CopyFilterTest) {

  _DataSetFiltertest->CopyFilter(_DataSetFiltertest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, CorrectCodeTest) {

  long resultlong = _DataSetFiltertest->CorrectCode(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, CountAndResolveTest) {

  _SimpleList* result_SimpleList = _DataSetFiltertest->CountAndResolve(*longtest, _Parametertest);
  //EXPECT_EQ (result_SimpleList*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, FilterDeletionsTest) {

  _DataSetFiltertest->FilterDeletions(_SimpleListtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, FindAllSitesLikeThisOneTest) {

  //TODO
  //_DataSetFiltertest->FindAllSitesLikeThisOne(*longtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, FindSpeciesNameTest) {

  long resultlong = _DataSetFiltertest->FindSpeciesName(*_Listtest, *_SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, FindUniqueSequencesTest) {

  //TODO
  //long resultlong = _DataSetFiltertest->FindUniqueSequences(*_SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, FreeUpMemoryTest) {

  long resultlong = _DataSetFiltertest->FreeUpMemory(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, FreezeTest) {

  _DataSetFiltertest->Freeze(*longtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GenerateConsensusStringTest) {

  _String result_String = _DataSetFiltertest->GenerateConsensusString(_SimpleListtest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GetCharTest) {

  char resultchar = _DataSetFiltertest->GetChar(*longtest, *longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GetDimensionTest) {

  long resultlong = _DataSetFiltertest->GetDimension(*booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GetExclusionsTest) {

  _String* result_String = _DataSetFiltertest->GetExclusions();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GetFilterCharactersTest) {

  _Matrix* result_Matrix = _DataSetFiltertest->GetFilterCharacters(*booltest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GetOriginalToShortMapTest) {

  long resultlong = _DataSetFiltertest->GetOriginalToShortMap(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GetSequenceCharactersTest) {

  _String* result_String = _DataSetFiltertest->GetSequenceCharacters(*longtest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GrabSiteTest) {

  _DataSetFiltertest->GrabSite(*longtest, *longtest, chartest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, GrabSite1Test) {

  //TODO
  //_DataSetFiltertest->GrabSite(*longtest, *longtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, HarvestFrequenciesTest) {

  _Matrix* result_Matrix = _DataSetFiltertest->HarvestFrequencies(*chartest, *chartest, *booltest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, HasDeletionsTest) {

  bool resultbool = _DataSetFiltertest->HasDeletions(*longtest, _AVLListtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__DataSetFilterTest, HasExclusionsTest) {

  //TODO
  //long resultlong = _DataSetFiltertest->HasExclusions(*longtest, _SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, IsConstantTest) {

  bool resultbool = _DataSetFiltertest->IsConstant(*longtest, *booltest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(DISABLED__DataSetFilterTest, LookupConversionTest) {

  long resultlong = _DataSetFiltertest->LookupConversion(*chartest, _Parametertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, MapStringToCharIndexTest) {

  long resultlong = _DataSetFiltertest->MapStringToCharIndex(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, MatchStartNEndTest) {

  _DataSetFiltertest->MatchStartNEnd(*_SimpleListtest, *_SimpleListtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, PairFilterTest) {

  //TODO
  //_DataSetFilter* result_DataSetFilter = _DataSetFiltertest->PairFilter(*longtest, *longtest);
  //EXPECT_EQ (result_DataSetFilter*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, PairwiseCompareTest) {

  _Matrix* result_Matrix = _DataSetFiltertest->PairwiseCompare(_SimpleListtest, _SimpleListtest);
  //EXPECT_EQ (result_Matrix*, 0);

}


TEST_F(DISABLED__DataSetFilterTest, PatternToSiteMapperTest) {

  //TODO
  //_DataSetFiltertest->PatternToSiteMapper(**sourcetest, **targettest, *chartest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, RetrieveStateTest) {

  //TODO
  //_DataSetFiltertest->RetrieveState(*longtest, *longtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, SetDimensionsTest) {

  _DataSetFiltertest->SetDimensions();
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, SetExclusionsTest) {

  _DataSetFiltertest->SetExclusions(_Stringtest, *booltest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, SetFilterTest) {

  //TODO
  //_DataSetFiltertest->SetFilter(_DataSettest, *chartest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, SetMapTest) {

  _DataSetFiltertest->SetMap(*_Stringtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, SetupConversionTest) {

  _DataSetFiltertest->SetupConversion();
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, SiteFrequencyTest) {

  long resultlong = _DataSetFiltertest->SiteFrequency(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, Translate2FrequenciesTest) {

  //TODO
  //long resultlong = _DataSetFiltertest->Translate2Frequencies(*_Stringtest, _Parametertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, Translate2Frequencies1Test) {

  //TODO
  //long resultlong = _DataSetFiltertest->Translate2Frequencies(*chartest, _Parametertest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(DISABLED__DataSetFilterTest, UnFreezeTest) {

  _DataSetFiltertest->UnFreeze(*longtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, XferwCorrectionTest) {

  //TODO
  //_DataSetFiltertest->XferwCorrection(*_Matrixtest, _Parametertest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, XferwCorrection1Test) {

  //TODO
  //_DataSetFiltertest->XferwCorrection(_Parametertest, _Parametertest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, XferwCorrection2Test) {

  //TODO
  //_DataSetFiltertest->XferwCorrection(longtest, _Parametertest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, internalToStrTest) {

  //TODO
  //_DataSetFiltertest->internalToStr(FILEtest, *_Stringtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, makeDynamicTest) {

  BaseRef resultBaseRef = _DataSetFiltertest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(DISABLED__DataSetFilterTest, operatorParenthsTest) {

  //TODO
  //_String result_String = _DataSetFiltertest->operatorParenths();
  //EXPECT_EQ (result_String, 0);

}


TEST_F(DISABLED__DataSetFilterTest, toFileStrTest) {

  _DataSetFiltertest->toFileStr(FILEtest);
  //EXPECT_EQ (_DataSetFiltertest, 0);

}


TEST_F(DISABLED__DataSetFilterTest, toStrTest) {

  BaseRef resultBaseRef = _DataSetFiltertest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
