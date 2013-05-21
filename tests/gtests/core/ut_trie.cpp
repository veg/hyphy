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
#include "trie.h"
#include "helperfunctions.h"

namespace {

// The fixture for testing class Foo.
class _TrieTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _TrieTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    _SimpleListtest = new _SimpleList();
    chartest = new char('h');
    longtest = new long(1);
    booltest = new bool();
    BaseReftest = new BaseRef();
    _Stringtest = new _String(FILEtest);
    _Trietest = new _Trie(_Stringtest);
  }

  virtual ~_TrieTest() {
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
    delete booltest;
    delete BaseReftest;
    delete _Stringtest;
    delete _Trietest;
    delete chartest;
    delete longtest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  _SimpleList* _SimpleListtest;
  bool* booltest;
  char* chartest;
  long* longtest;
  BaseRef* BaseReftest;
  _String* _Stringtest;
  _Trie* _Trietest;
};


TEST_F(_TrieTest, AlphabetTest) {

  _String result_String = _Trietest->Alphabet();
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_TrieTest, ClearTest) {

  _Trietest->Clear(*booltest);
  //EXPECT_EQ (_Trietest, 0);

}


TEST_F(_TrieTest, DeleteTest) {

  bool resultbool = _Trietest->Delete(*_Stringtest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TrieTest, Delete1Test) {

  bool resultbool = _Trietest->Delete(chartest);
  //EXPECT_EQ (resultbool, 0);

}


TEST_F(_TrieTest, Delete2Test) {

  //long resultlong = _Trietest->Delete(*_Listtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, DumpRawTest) {

  //_Trietest->DumpRaw();
  //EXPECT_EQ (_Trietest, 0);

}


TEST_F(_TrieTest, DuplicateTest) {

  _Trietest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_Trietest, 0);

}


TEST_F(_TrieTest, FindTest) {

  long resultlong = _Trietest->Find(*_Stringtest, _SimpleListtest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, Find1Test) {

  long resultlong = _Trietest->Find(*chartest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, FindNextLetterTest) {

  //long resultlong = _Trietest->FindNextLetter(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, FindNextUnusedIndexTest) {

  //long resultlong = _Trietest->FindNextUnusedIndex(*booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, GetValueTest) {

  long resultlong = _Trietest->GetValue(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, GetValueFromStringTest) {

  long resultlong = _Trietest->GetValueFromString(*_Stringtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, InsertTest) {

  long resultlong = _Trietest->Insert(*_Stringtest, *longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, Insert1Test) {

  long resultlong = _Trietest->Insert(chartest, *longtest, *booltest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, Insert2Test) {

  //long resultlong = _Trietest->Insert(*_Listtest, _SimpleListtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, InsertNextLetterTest) {

  //long resultlong = _Trietest->InsertNextLetter(*chartest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_TrieTest, RetrieveKeyByPayloadTest) {

  _String result_String = _Trietest->RetrieveKeyByPayload(*longtest);
  //EXPECT_EQ (result_String, 0);

}


TEST_F(_TrieTest, RetrieveStringFromPathTest) {

  //_String* result_String = _Trietest->RetrieveStringFromPath(*_SimpleListtest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_TrieTest, SetAlphabetTest) {

  //_Trietest->SetAlphabet(_Stringtest, *booltest);
  //EXPECT_EQ (_Trietest, 0);

}


TEST_F(_TrieTest, UpdateValueTest) {

  _Trietest->UpdateValue(*longtest, *longtest);
  //EXPECT_EQ (_Trietest, 0);

}


TEST_F(_TrieTest, makeDynamicTest) {

  BaseRef resultBaseRef = _Trietest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_TrieTest, toStrTest) {

  BaseRef resultBaseRef = _Trietest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
