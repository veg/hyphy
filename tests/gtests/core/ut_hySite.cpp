/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon@cfenet.ubc.ca)
  Steven Weaver (sweaver@ucsd.edu)

Module Developers:
	Lance Hepler (nlhepler@gmail.com)
	Martin Smith (msmith@ucsd.edu)

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

#ifndef _HY_GENSITE_
#define _HY_GENSITE_

#include <iostream>
#include "gtest/gtest.h"
#include "site.h"
//#include <string>

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Range;
using ::testing::Combine;

namespace {

// The fixture for testing class Foo.
class _SiteTest : public ::testing::Test {

  protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _SiteTest() {
      // You can do set-up work for each test here.
    }

    virtual ~_SiteTest() {
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
    }

    // Objects declared here can be used by all tests in the test case for Foo.
};


TEST_F(_SiteTest, ConstructorsTest)
{
  // Normal
  _String test_s("hyphy");
  _Site test(test_s);
  EXPECT_EQ(5, test.s_length);
  EXPECT_EQ(-1, test.getRefNo());
  EXPECT_FALSE(test.isComplete());
  test.setRefNo(162);
  EXPECT_FALSE(test.isComplete());
  EXPECT_EQ(162, test.getRefNo());
  EXPECT_STREQ("hyphy", test);
  test.complete();
  EXPECT_TRUE(test.isComplete());

  _Site test2;
  test2 << "hyphy";
  EXPECT_EQ(5, test2.s_length);
  EXPECT_EQ(-1, test2.getRefNo());
  EXPECT_FALSE(test2.isComplete());
  test2.setRefNo(162);
  EXPECT_FALSE(test2.isComplete());
  EXPECT_EQ(162, test2.getRefNo());
  EXPECT_STREQ("hyphy", test2);
  test2.complete();
  EXPECT_TRUE(test2.isComplete());

  char test_c = 'h';
  _Site test3(test_c);
  EXPECT_EQ(1, test3.s_length);
  EXPECT_EQ(-1, test3.getRefNo());
  EXPECT_FALSE(test3.isComplete());
  test3.setRefNo(162);
  EXPECT_FALSE(test3.isComplete());
  EXPECT_EQ(162, test3.getRefNo());
  EXPECT_STREQ("h", test3);
  test3.complete();
  EXPECT_TRUE(test3.isComplete());

  long test_i = 162;
  _Site test4(test_i);
  EXPECT_EQ(0, test4.s_length);
  EXPECT_EQ(162, test4.getRefNo());
  EXPECT_FALSE(test4.isComplete());
  test4.setRefNo(163);
  EXPECT_FALSE(test4.isComplete());
  EXPECT_EQ(163, test4.getRefNo());
  test4.complete();
  EXPECT_TRUE(test4.isComplete());
}

TEST_F(_SiteTest, CompleteTest)
{
  // Normal
  _Site test;
  EXPECT_FALSE(test.isComplete());
  test.complete();
  EXPECT_TRUE(test.isComplete());

  _Site test1;
  EXPECT_FALSE(test1.isComplete());
  test1.complete();
  EXPECT_TRUE(test1.isComplete());
  test1.setRefNo(163);
  EXPECT_FALSE(test1.isComplete());
  test1.complete();
  EXPECT_TRUE(test1.isComplete());
}

TEST_F(_SiteTest, duplicateTest)
{
  // mod original
  _Site test;
  test << "hyphy";
  EXPECT_EQ(5, test.s_length);
  EXPECT_STREQ("hyphy", test);
  EXPECT_EQ(-1, test.getRefNo());
  _Site test2;
  test2.duplicate(&test);
  EXPECT_EQ(5, test2.s_length);
  EXPECT_STREQ("hyphy", test2);
  EXPECT_EQ(-1, test2.getRefNo());
  test.setRefNo(162);
  EXPECT_EQ(162, test.getRefNo());
  EXPECT_EQ(-1, test2.getRefNo());
  EXPECT_FALSE(test.isComplete());
  EXPECT_FALSE(test2.isComplete());
  test.complete();
  EXPECT_TRUE(test.isComplete());
  EXPECT_FALSE(test2.isComplete());
  test2.complete();
  EXPECT_TRUE(test2.isComplete());

  // mod duplicate
  _Site test3;
  test3 << "hyphy";
  EXPECT_EQ(5, test3.s_length);
  EXPECT_STREQ("hyphy", test3);
  EXPECT_EQ(-1, test3.getRefNo());
  _Site test4;
  test4.duplicate(&test3);
  EXPECT_EQ(5, test4.s_length);
  EXPECT_STREQ("hyphy", test4);
  EXPECT_EQ(-1, test4.getRefNo());
  test4.setRefNo(162);
  EXPECT_EQ(162, test4.getRefNo());
  EXPECT_EQ(-1, test3.getRefNo());
  EXPECT_FALSE(test3.isComplete());
  EXPECT_FALSE(test4.isComplete());
  test4.complete();
  EXPECT_TRUE(test4.isComplete());
  EXPECT_FALSE(test3.isComplete());
  test3.complete();
  EXPECT_TRUE(test4.isComplete());
}

TEST_F(_SiteTest, clearTest)
{
  // Normal
  _Site test;
  EXPECT_FALSE(test.isComplete());
  test.complete();
  EXPECT_TRUE(test.isComplete());
  test << "hyphy";
  EXPECT_STREQ("hyphy", test);
  test.clear();
  EXPECT_FALSE(test.isComplete());
  EXPECT_EQ(-1, test.getRefNo());
  EXPECT_STREQ(NULL, test);
}

TEST_F(_SiteTest, getRefNoTest)
{
  _Site test;
  EXPECT_EQ(-1, test.getRefNo());
  test << "hyphy";
  EXPECT_EQ(-1, test.getRefNo());
  test.complete();
  EXPECT_EQ(-1, test.getRefNo());
  EXPECT_TRUE(test.isComplete());
  EXPECT_EQ(-1, test.getRefNo());
  test.setRefNo(162);
  EXPECT_FALSE(test.isComplete());
  EXPECT_EQ(162, test.getRefNo());
  test.complete();
  EXPECT_EQ(162, test.getRefNo());
  EXPECT_TRUE(test.isComplete());
}

TEST_F(_SiteTest, isCompleteTest)
{
  _Site test;
  EXPECT_EQ(-1, test.getRefNo());
  test << "hyphy";
  EXPECT_EQ(-1, test.getRefNo());
  EXPECT_FALSE(test.isComplete());
  test.setRefNo(162);
  EXPECT_FALSE(test.isComplete());
  test.complete();
  EXPECT_TRUE(test.isComplete());
  test.setRefNo(162);
  EXPECT_FALSE(test.isComplete());
  test.complete();
  EXPECT_TRUE(test.isComplete());
}

TEST_F(_SiteTest, setRefNoTest)
{
  _Site test;
  test << "hyphy";
  // Normal use
  EXPECT_EQ(-1, test.getRefNo());
  EXPECT_EQ(-1, test.getRefNo());
  EXPECT_FALSE(test.isComplete());
  test.setRefNo(162);
  EXPECT_EQ(162, test.getRefNo());
  EXPECT_FALSE(test.isComplete());
  test.complete();
  EXPECT_EQ(162, test.getRefNo());
  EXPECT_TRUE(test.isComplete());
  test.setRefNo(162);
  EXPECT_FALSE(test.isComplete());
  EXPECT_EQ(162, test.getRefNo());
  test.complete();
  EXPECT_EQ(162, test.getRefNo());
  EXPECT_TRUE(test.isComplete());

  int tests[5] = {162, 2, 1, 0, -1,};
  for (int i = 0; i < 5; i++) {
    test.setRefNo(tests[i]);
    EXPECT_EQ(tests[i], test.getRefNo());
    EXPECT_FALSE(test.isComplete());
    test.complete();
    EXPECT_EQ(tests[i], test.getRefNo());
    EXPECT_TRUE(test.isComplete());
    test.setRefNo(tests[i]);
    EXPECT_FALSE(test.isComplete());
    EXPECT_EQ(tests[i], test.getRefNo());
    test.complete();
    EXPECT_EQ(tests[i], test.getRefNo());
    EXPECT_TRUE(test.isComplete());
  }
}

} // namespace


#endif
