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

#include <iostream>
#include "gtest/gtest.h"
#include <stdio.h>
#include <stdlib.h>

//Generate necessary includes from the respective implementation file
#include "baseobj.h"

namespace {

// The fixture for testing class Foo.
class BaseObjTest : public ::testing::Test {

protected:
  // You can remove any or all of the following functions if its body
  // is empty.

    BaseObjTest() {
    // You can do set-up work for each test here.
    // Create objects of every type needed. Performance doesn't matter.

        //FILEtest = fopen ("./tests/gtests/res/HIV_gp120.nex" , "r");

    BaseObjtest = new BaseObj();
  }

  virtual ~BaseObjTest() {
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
    delete BaseObjtest;
      //fclose (FILEtest);
  }

  BaseObj* BaseObjtest;
};


TEST_F(BaseObjTest, makeDynamicTest) {

  BaseObj* resultBaseObj = BaseObjtest->makeDynamic();
  //EXPECT_EQ (resultBaseObj*, 0);

}


TEST_F(BaseObjTest, toErrStrTest) {

  BaseRef resultBaseRef = BaseObjtest->toErrStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


TEST_F(BaseObjTest, toFileStrTest) {

    //BaseObjtest->toFileStr(FILEtest);
  //EXPECT_EQ (BaseObjtest, 0);

}


TEST_F(BaseObjTest, toStrTest) {

  BaseRef resultBaseRef = BaseObjtest->toStr();
  //EXPECT_EQ (resultBaseRef, 0);

}


}
