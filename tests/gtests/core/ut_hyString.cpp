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


TEST_F (_hyStringTest, StringConstuctorAndConverstionTests) {
  _String blank,
          literal_string ("The quick brown fox jumps over the lazy dog"),
          literal_w_string (L"The quick brown fox jumps over the lazy dog");
  
  EXPECT_EQ (0UL, blank.Length ()) << "Default string constructor should make an empty string";
  EXPECT_EQ (NULL, blank.getStr()) << "Default string constructor should set string data to NULL";
  
  _String * test_string = new _String (literal_string);
  
  EXPECT_TRUE (literal_string.Equal(test_string)) << "Full stack copy failed";
  delete test_string;
  test_string = new _String (literal_string, 5, 10);
  EXPECT_STREQ ("uick b", (const char*)*test_string) << "Partial stack copy (or char*) failed";
  delete test_string;
  test_string = new _String (literal_string, -1, 256);
  EXPECT_TRUE (literal_string.Equal(test_string)) << "Stack copy with ranges implying a full copy failed";
  delete test_string;
  test_string = new _String (blank, 2, 5);
  EXPECT_TRUE (test_string->Length () == 0UL && test_string->getStr() != NULL) << "Stack copy constructor from an empty string must return an empty string regardless of the range";
  delete test_string;
  
  _String single_char ('W');
  EXPECT_TRUE (single_char.Length () == 1UL && single_char.getChar (0) == 'W') << "Single character constructor failed";
  
  FILE * temp_file = tmpfile();
  fprintf (temp_file, "%s", literal_string.getStr());
  fflush(temp_file);
  test_string = new _String (temp_file);
  EXPECT_TRUE (literal_string.Equal(test_string)) << "String constructor from file did not work correctly";
  test_string->AddAReference();
  _String test_copy (test_string);
  EXPECT_TRUE (test_string->SingleReference()) << "Incorrect reference count handling for CopyDynamicString";
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
    EXPECT_STREQ(buffer, (const char*)double_string) << "Incorrect _String (Parameter) constructor (or char* conversion) for " << dn;
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

TEST_F (_hyStringTest, StringManipulationOps) {
  
  for (long k = 0; k < 1024L; k++) {
    _String s (_String::Random (genrand_int32() % 32)),
            f (s);
    
    f.Flip();
    f.Flip();
    
    EXPECT_TRUE (s == f) << "s.Flip().Flip() roundtrip failed for '" << s << "'";
    
    unsigned long string_length = s.Length();
    long     breakpoint   = genrand_int32 () % MAX (string_length + 10L, 1L) - 5L;
                            
                            
    _String reassembled = breakpoint >= 0L ? (s.Cut (0, breakpoint) & s.Cut (breakpoint+1, -1L)) :
                                             s.Cut (0, breakpoint);
                                             
    EXPECT_TRUE (s == reassembled) << "Failed string cut test for '" << (const char*) s << "' at " << breakpoint;
    if (breakpoint >= 0L) {
      reassembled = s.Chop (breakpoint+1, -1) & s.Chop (0, breakpoint);    
      EXPECT_TRUE (s == reassembled) << "Failed string chop test for '" << (const char*) s << "' at " << breakpoint 
            << " had '" << (const char*) reassembled << "'" ;
      _String p1 (s), p2 (s);
      p1.Delete (0, breakpoint);
      p2.Delete (breakpoint+1, -1);
      EXPECT_TRUE (s == (p2 & p1)) << "Failed to reassemble a string from partial delections for '" << (const char*) s << "' at " << breakpoint << " had '" << (const char*) reassembled << "'" ;
      
      reassembled.Insert ('X', breakpoint);
      breakpoint = MIN (breakpoint, string_length);
      EXPECT_STREQ ("X", reassembled.Cut (breakpoint, breakpoint).getStr()) << "Insert character failed for '" << (const char*) s << "' at " << breakpoint;
      
    }
    
  } 
  
  EXPECT_EQ (0UL, (empty & empty).Length()) << "empty & empty must be empty";
  EXPECT_EQ (0UL, empty.Chop (0,0).Length()) << "empty.Cut (anything) must be empty";
}

TEST_F (_hyStringTest, StringCaseChange) {
  for (long k = 0; k < 1024L; k++) {
    _String s (_String::Random (genrand_int32() % 32)),
            l (s.LoCase()),
            u (s.UpCase());
            
     
    if (l != u) {
      EXPECT_EQ (l, u.LoCase()) << "Case conversion failed";
    }
    
  }
}

TEST_F (_hyStringTest, StringSearching) {

  _String numbers ("0123456789"), 
          letters ("abcdefghijklmnopqrstuvwxyz");
                    

  for (long k = 0; k < 1024L; k++) {
    _String number_string (_String::Random (genrand_int32() % 32, &numbers)),
            letter_string (_String::Random (genrand_int32() % 32, &letters)),
            upcase_letter_string (letter_string.UpCase()),
            combined = letter_string & number_string & upcase_letter_string & number_string;
            
            
    if (number_string.Length () == 0) {
      EXPECT_EQ (HY_NOT_FOUND, combined.Find (number_string)) << "Finding an empty string must fail";
    } else {
      long point1 = letter_string.Length(),
           point2 = letter_string.Length() * 2 + number_string.Length();
           
      EXPECT_EQ (point1, combined.Find (number_string)) << "Find string (whole string) failed";
      EXPECT_EQ (point2 , combined.Find (number_string, letter_string.Length() + number_string.Length())) << "Find string (partial string) failed";
      
      EXPECT_EQ (point2, combined.FindBackwards (number_string)) << "Find string backwards (whole string) failed";
      EXPECT_EQ (HY_NOT_FOUND, combined.FindBackwards (upcase_letter_string, 0, point1 + number_string.Length() - 1 )) << "Incorrect result for '" << combined << ".FindBackwards('" <<
                    upcase_letter_string << "')";
                    
      EXPECT_EQ (HY_NOT_FOUND, combined.FindBackwards (number_string, combined.Length(), 0)) << "Invalid range returned a valid Find result";
      
      
    }
    if (letter_string.Length()) {
        EXPECT_TRUE (combined.ContainsSubstring (upcase_letter_string)) << "Incorrect result for '" << combined << ".ContainsSubstring('" <<
                    upcase_letter_string << "')";
         EXPECT_EQ (HY_NOT_FOUND, combined.Find (letter_string, letter_string.Length())) << "Incorrectly found uppercase string";
         EXPECT_EQ (letter_string.Length() + number_string.Length(), combined.FindAnyCase (letter_string, letter_string.Length())) << "Failed with a case insensitive search for " << letter_string << " in " << combined;
 
    }
      
  }
}

