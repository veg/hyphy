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
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "hy_strings.h"
#include "list.h"
#include "translationtable.h"
#include "site.h"

namespace {

// The fixture for testing class Foo.
class _TranslationTableTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _TranslationTableTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    unsignedtest = new unsigned();
    _Stringtest = new _String(FILEtest);
    BaseReftest = new BaseRef();
    _TranslationTabletest = new _TranslationTable(*_TranslationTabletest);
  }

  virtual ~_TranslationTableTest() {
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
    delete BaseReftest;
    delete _TranslationTabletest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  unsigned* unsignedtest;
  _String* _Stringtest;
  BaseRef* BaseReftest;
  _TranslationTable* _TranslationTabletest;
};


TEST_F(_TranslationTableTest, AddBaseSetTest) {

  _TranslationTabletest->AddBaseSet(*_Stringtest);
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, AddTokenCodeTest) {

  _TranslationTabletest->AddTokenCode(*chartest, *_Stringtest);
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, CheckTypeTest) {

  bool resultbool = _TranslationTabletest->CheckType(*chartest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TranslationTableTest, CheckValidAlphabetTest) {

  bool resultbool = _TranslationTabletest->CheckValidAlphabet(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TranslationTableTest, ClearTest) {

  _TranslationTabletest->Clear();
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, CodeToLetterTest) {

  char resultchar = _TranslationTabletest->CodeToLetter(longtest);
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_TranslationTableTest, ConvertCodeToLettersTest) {

  _String result_String = _TranslationTabletest->ConvertCodeToLetters(*longtest, *chartest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_TranslationTableTest, DetectTypeTest) {

  char resultchar = _TranslationTabletest->DetectType();
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_TranslationTableTest, DuplicateTest) {

  _TranslationTabletest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, GetDefaultAlphabetTest) {

  _String* result_String = _TranslationTabletest->GetDefaultAlphabet(*longtest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_TranslationTableTest, GetGapCharTest) {

  char resultchar = _TranslationTabletest->GetGapChar();
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_TranslationTableTest, GetSkipCharTest) {

  char resultchar = _TranslationTabletest->GetSkipChar();
  //EXPECT_EQ (resultchar, 0);

}


TEST_F(_TranslationTableTest, IsCharLegalTest) {

  bool resultbool = _TranslationTabletest->IsCharLegal(*chartest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TranslationTableTest, LengthOfAlphabetTest) {

  long resultlong = _TranslationTabletest->LengthOfAlphabet();
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TranslationTableTest, MergeTablesTest) {

  _TranslationTable* result_TranslationTable = _TranslationTabletest->MergeTables(_TranslationTabletest);
  //EXPECT_EQ (result_TranslationTable*, 0);

}


TEST_F(_TranslationTableTest, PrepareForChecksTest) {

  _TranslationTabletest->PrepareForChecks();
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, RetrieveCharactersTest) {

  _String* result_String = _TranslationTabletest->RetrieveCharacters();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_TranslationTableTest, SetStandardTypeTest) {

  _TranslationTabletest->SetStandardType(*chartest);
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, SplitTokenCodeTest) {

  _TranslationTabletest->SplitTokenCode(*longtest);
  //EXPECT_EQ (_TranslationTabletest, 0);

}


TEST_F(_TranslationTableTest, TokenCodeTest) {

  bool resultbool = _TranslationTabletest->TokenCode(*chartest, longtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TranslationTableTest, TokenCode1Test) {

  long resultlong = _TranslationTabletest->TokenCode(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TranslationTableTest, makeDynamicTest) {

  BaseRef resultBaseRef = _TranslationTabletest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
