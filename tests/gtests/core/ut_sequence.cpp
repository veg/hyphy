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
#include "sequence.h"
#include "errorfns.h"
#include "stdio.h"
#include "string.h"
#include "helperfunctions.h"
#include "stdio.h"

namespace {

// The fixture for testing class Foo.
class _CStringTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  _CStringTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

    FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    longtest = new long();
    unsignedtest = new unsigned();
    _Stringtest = new _String(FILEtest);
    chartest = new char();
    _CStringtest = new _CString(*_Stringtest);
    BaseReftest = new BaseRef();
    booltest = new bool();
  }

  virtual ~_CStringTest() {
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
    delete _CStringtest;
    delete BaseReftest;
    delete booltest;
    fclose (FILEtest);
  }

  FILE* FILEtest;
  long* longtest;
  unsigned* unsignedtest;
  _String* _Stringtest;
  char* chartest;
  _CString* _CStringtest;
  BaseRef* BaseReftest;
  bool* booltest;
};


TEST_F(_CStringTest, BestCompressTest) {

  _Parameter result_Parameter = _CStringtest->BestCompress(*chartest, *longtest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_CStringTest, DecompressTest) {

  _String* result_String = _CStringtest->Decompress();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_CStringTest, DecompressFrequencyTest) {

  _String* result_String = _CStringtest->DecompressFrequency();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_CStringTest, DecompressLZWTest) {

  _String* result_String = _CStringtest->DecompressLZW();
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_CStringTest, DuplicateTest) {

  _CStringtest->Duplicate(*BaseReftest);
  //EXPECT_EQ (_CStringtest, 0);

}


TEST_F(_CStringTest, FinalizeTest) {

  _CStringtest->Finalize();
  //EXPECT_EQ (_CStringtest, 0);

}


TEST_F(_CStringTest, FreeUpMemoryTest) {

  long resultlong = _CStringtest->FreeUpMemory(*longtest);
  //EXPECT_EQ (resultlong, 0);

}


TEST_F(_CStringTest, FrequencyCompressTest) {

  _Parameter result_Parameter = _CStringtest->FrequencyCompress(*chartest, *booltest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_CStringTest, LZWCompressTest) {

  _Parameter result_Parameter = _CStringtest->LZWCompress(*chartest);
  //EXPECT_EQ (result_Parameter, 0);

}


TEST_F(_CStringTest, SelectAlphaTest) {

  _String* result_String = _CStringtest->SelectAlpha(*chartest);
  //EXPECT_EQ (result_String*, 0);

}


TEST_F(_CStringTest, makeDynamicTest) {

  BaseRef resultBaseRef = _CStringtest->makeDynamic();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(_CStringTest, operatorDoubleLessTest) {

  _CStringtest->operatorDoubleLess(_Stringtest);
  //EXPECT_EQ (_CStringtest, 0);

}


TEST_F(_CStringTest, operatorDoubleLess1Test) {

  _CStringtest->operatorDoubleLess(*chartest);
  //EXPECT_EQ (_CStringtest, 0);

}


}
