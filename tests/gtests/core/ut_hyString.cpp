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
#include <math.h>


//Generate necessary includes from the respective implementation file
#include "hy_strings.h"
#include "helperfunctions.h"

namespace {
  // The fixture for testing class Foo.
  class _hyStringTest : public ::testing::Test {

  protected:

    _hyStringTest() {
      
    }

    virtual ~_hyStringTest() {
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
    
  public:
    // Per-test-case set-up.
    // Called before the first test in this test case.
    // Can be omitted if not needed.
    static void SetUpTestCase() {
        //shared_resource_ = new ...;
    }

    // Per-test-case tear-down.
    // Called after the last test in this test case.
    // Can be omitted if not needed.
    static void TearDownTestCase() {
        // delete shared_resource_;
        // shared_resource_ = NULL;
    }

    
  };

}


TEST_F (_hyStringTest, StringConstuctorTests) {
  _String blank,
          literal_string ("The quick brown fox jumps over the lazy dog"),
          literal_w_string (L"The quick brown fox jumps over the lazy dog");
  
  EXPECT_EQ (0UL, blank.Length ()) << "Default string constructor should make an empty string";
  EXPECT_EQ (NULL, blank.getStr()) << "Default string constructor should set string data to NULL";
  
  _String * test_string = new _String (literal_string);
  
  EXPECT_TRUE (literal_string.Equal(test_string)) << "Full stack copy failed";
  delete test_string;
  test_string = new _String (literal_string, 5, 10);
  EXPECT_STREQ ("uick b", test_string->getStr()) << "Partial stack copy failed";
  delete test_string;
  test_string = new _String (literal_string, -1, 256);
  EXPECT_TRUE (literal_string.Equal(test_string)) << "Stack copy with ranges implying a full copy failed";
  delete test_string;
  test_string = new _String (blank, 2, 5);
  EXPECT_TRUE (test_string->Length () == 0UL && test_string->getStr() != NULL) << "Stack copy constructor from an empty string must return an empty string regardless of the range";
  
  _String single_char ('W');
  EXPECT_TRUE (single_char.Length () == 1UL && single_char.getChar (0) == 'W') << "Single character constructor failed";
  
  FILE * temp_file = tmpfile();
  fprintf (temp_file, "%s", literal_string.getStr());
  fflush(temp_file);
  test_string = new _String (temp_file);
  EXPECT_TRUE (literal_string.Equal(test_string)) << "String constructor from file did not work correctly";
  delete test_string;
  fclose(temp_file);
  
  temp_file = NULL;
  test_string = new _String (temp_file);
  EXPECT_TRUE (blank.Equal(test_string)) << "String constructor from NULL file is expected to return an empty string";
  delete test_string;
  
  EXPECT_TRUE (literal_string.Equal (&literal_w_string)) << "Strings constructed from wchar_t* and char* did not match";
  
  char buffer[256];
  for (unsigned long i = 0; i < 50; i++) {
    _Parameter dn = genrand_real2() * pow (10., 200 - genrand_int32() % 400);
    _String double_string (dn);
    snprintf(buffer, 255, PRINTF_FORMAT_STRING, dn);
    EXPECT_STREQ(buffer, double_string.getStr()) << "Incorrect _String (Parameter) constructor for " << dn;
  }
  
  literal_w_string = blank;
  blank = literal_string;
  EXPECT_TRUE (literal_string.Equal (&blank)) << "operator = failed";
  
}

TEST_F (_hyStringTest, StringAccessorTests) {
  _String s ("I am the great Cornholio");
  const _String sc (s);
  
  EXPECT_EQ('t', s[5]) << "String [] access (valid range) failed";
  EXPECT_EQ(0, s[-1]) << "String [] access (-1) failed";
  EXPECT_EQ(0, s[s.Length()]) << "String [] access (length + 1) failed";
  EXPECT_EQ('t', sc(5)) << "String () access (valid range) failed";
  EXPECT_EQ(0, sc(-1)) << "String () access (-1) failed";
  EXPECT_EQ(0, sc(sc.Length())) << "String () access (length + 1) failed";
  EXPECT_EQ('t', sc.getChar(5)) << "String getChar access (valid range) failed";
  EXPECT_EQ(0, sc.getChar(-1)) << "String getChar access (-1) failed";
  EXPECT_EQ(0, sc.getChar(sc.Length())) << "String getChar access (length + 1) failed";
}

TEST_F (_hyStringTest, StringLogicalOps) {
  const _String s1 ("beavis"),
                s2 ("butthead"),
                s3 ("beaviss");
  
  EXPECT_TRUE (s1 < s2) << "beavis > butthead";
  EXPECT_TRUE (s1 <= s2) << "beavis > butthead";
  EXPECT_TRUE (s3 > s1) << "beavis > beaviss";
  EXPECT_TRUE (s3 >= s1) << "beavis > beaviss";
  EXPECT_TRUE (s2 > empty) << "butthead << \"\"";
  EXPECT_TRUE (s1 == s1) << "beavis != beavis";
  EXPECT_TRUE (s1 != s2) << "beavis == butthead";
  
}


