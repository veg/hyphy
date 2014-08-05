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
#include "hy_list_reference.h"
#include "hy_list_numeric.h"
#include "hy_globals.h"

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

_String numbers ("0123456789"),
        letters ("abcdefghijklmnopqrstuvwxyz"),
        alnum  = numbers & letters & letters.UpCase(),
        alpha_ = letters & letters.UpCase() & '_',
        alnum_ = alnum & '_',
        punctuation ("`!@#$%^&*()+=/,';][\\"),
        whitespaces (" \t\n\r\v\f");

const char test_case [] = "The quick brown fox jumps over the lazy dog";


TEST_F (_hyStringTest, StringConstuctorAndConverstionTests) {
  _String blank,
          literal_string (test_case),
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
  test_string->DuplicateErasing(&literal_string);
  EXPECT_TRUE (test_string->Equal(&literal_string)) << "DuplicateErasing did not correctly copy the source string";
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
  
  test_string = (_String*) literal_string.toStr();
  EXPECT_TRUE (literal_string.Equal (test_string)) << "toStr () did not return the same string back";
  delete test_string;
  
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
  
    // check setChar
  
  s.setChar (0, '!');
  EXPECT_EQ ('!', s.getChar (0)) << "setChar failed with a valid index";
  s.setChar (0, 'I');
  s.setChar (s.Length() + 1, 'X');
  s.setChar (-1, 'Y');
  EXPECT_EQ (sc, s) << "setChar with an invalid index range modified the string";
  
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
  
  _String a_literal = test_case,
          another_literal;
  a_literal.Trim (5,2);
  EXPECT_EQ (empty, a_literal) << "Trim with an invalid range must produce an empty string";
  a_literal = test_case;
  a_literal.Trim(2,5);
  EXPECT_EQ (a_literal, _String (test_case, 2,5)) << "Trim with a valid range failed compared to the substring constructor";
  a_literal = test_case;
  a_literal.Trim(0,HY_NOT_FOUND);
  EXPECT_EQ (a_literal, _String(test_case)) << "Trim which should leave the string intact failed";
  
  another_literal = a_literal;
  a_literal.StripQuotes();
  EXPECT_EQ (a_literal, another_literal) << "Strip quotes on a string not enclosed in \" marks changed the string";
  another_literal = _String ('[') & a_literal & ']';
  another_literal.StripQuotes('[',']');
  EXPECT_EQ(a_literal, another_literal) << "StripQuotes '[string]' did not return 'string'";
  
  another_literal = "aa";
  another_literal.StripQuotes('a','a');
  EXPECT_EQ (another_literal, empty) << "StripQuotes on '' did not return an empty string";
  
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
      
      EXPECT_EQ (combined.Find (number_string(0)), point1) << "Find (char) failed";
      EXPECT_EQ (combined.Find ('0'), combined.Find (_String ("0"))) <<
        "Find 'char' and Find String(char) returned unequal results";
      
      EXPECT_EQ (HY_NOT_FOUND, number_string.Find ('a')) << "Found a letter in a numeric string";
      EXPECT_EQ (combined.contains ('z'), combined.Find ('z') != HY_NOT_FOUND) << "contains (char) did not match Find results";
      EXPECT_EQ (combined.contains ("X"), combined.Find ("X") != HY_NOT_FOUND) << "contains (_String) did not match Find results";
      EXPECT_EQ (number_string.Equal('1'), number_string == _String("1")) << "equal (char) did not match '==' results";
      
      if (letter_string.sLength) {
        EXPECT_TRUE   (combined.startswith(letter_string)) << "AB starts with A failed (match case), A == '" << letter_string << "'";
        EXPECT_FALSE  (letter_string.startswith(combined)) << "NOT AB startswith A failed";
        EXPECT_TRUE   (combined.startswith(upcase_letter_string, false)) << "aB startswith A failed (ignore case)";
        EXPECT_FALSE  (combined.startswith(upcase_letter_string, true)) << "NOT aB startswith A failed (match case)";
      }
      EXPECT_TRUE (combined.endswith (number_string)) << "AB endswith B failed (match case)";
      EXPECT_FALSE(number_string.endswith (combined)) << "NOT B endswith AB failed (match case)" << number_string << "/" << combined;
      
      
    }
    if (letter_string.Length()) {
        EXPECT_TRUE (combined.ContainsSubstring (upcase_letter_string)) << "Incorrect result for '" << combined << ".ContainsSubstring('" <<
                    upcase_letter_string << "')";
         EXPECT_EQ (HY_NOT_FOUND, combined.Find (letter_string, letter_string.Length())) << "Incorrectly found uppercase string";
         EXPECT_EQ (letter_string.Length() + number_string.Length(), combined.FindAnyCase (letter_string, letter_string.Length())) << "Failed with a case insensitive search for " << letter_string << " in " << combined;
 
    }
      
  }
}



TEST_F (_hyStringTest, EqualWithWildChar) {
  _String match_me (test_case),
          random_string (_String::Random(10L, &numbers));
  
  EXPECT_TRUE (match_me.EqualWithWildChar("*qui*ck****dog*", '*'));
  EXPECT_TRUE  (match_me.EqualWithWildChar(test_case, '*'));
  EXPECT_TRUE  (match_me.EqualWithWildChar(test_case, 'a'));
  EXPECT_TRUE  (match_me.EqualWithWildChar("*quick*", '*'));
  EXPECT_TRUE  (match_me.EqualWithWildChar("*quick*dog", '*'));
  EXPECT_TRUE  (match_me.EqualWithWildChar("?", '?'));
  EXPECT_TRUE  (empty.EqualWithWildChar(empty, '?'));
  EXPECT_FALSE (match_me.EqualWithWildChar("*quick*ogd", '*'));
  EXPECT_FALSE (match_me.EqualWithWildChar("*quick", '*'));
  EXPECT_TRUE  (match_me.EqualWithWildChar("*quick***", '*'));
  EXPECT_FALSE (match_me.EqualWithWildChar(empty, '?'));
  EXPECT_FALSE (match_me.EqualWithWildChar(random_string,'*'));
  
}

TEST_F (_hyStringTest, StringReplace) {
  _String case1 ("ABABABAB"),
          case2 ("XXXXXXXXXXY"),
          r;
  
  EXPECT_EQ (empty.Replace ("A","B",false).Length(), 0L) << "Replace called on an empty string must return an empty string";
  EXPECT_EQ (case1.Replace (empty, "C", true), case1) << "Replace an empty string must do nothing";
  
  r = case1.Replace ("B",empty,false);
  EXPECT_EQ (r, _String("AABABAB")) << "Failed replace the first 'B' with an empty string in '" << case1 << "' / " << r;
  
  r = case1.Replace ("B",empty,true);
  EXPECT_EQ (r, _String("AAAA")) << "Failed replace all 'B' with an empty string in '" << case1 << "' / " << r;
  
  r = case1.Replace ("B","CDE",true);
  EXPECT_EQ (r, _String("ACDEACDEACDEACDE")) << "Failed replace all 'B' with an 'CDE' string in '" << case1 << "' / " << r;

  r = case1.Replace ("AB",empty,true);
  EXPECT_EQ (r, empty) << "Failed replace all 'AB' with an '' in '" << case1 << "' / " << r;

  r = case1.Replace ("ABA","123",true);
  EXPECT_EQ (r, _String("123B123B")) << "Failed replace all 'ABA' with '123' in '" << case1 << "' / " << r;

  r = case2.Replace ("XXXXY","Beavis",true);
  EXPECT_EQ (r, _String("XXXXXXBeavis")) << "Failed replace all 'XXXXY' with 'Beavis' in '" << case2 << "' / " << r;

  r = case2.Replace ("XXXXY","Beavis",false);
  EXPECT_EQ (r, _String("XXXXXXBeavis")) << "Failed replace the only 'XXXXY' with 'Beavis' in '" << case2 << "' / " << r;
  
}

TEST_F (_hyStringTest, StringFormatters) {
  EXPECT_STREQ(_String::FormatTimeString(0L), "00:00:00");
  EXPECT_STREQ(_String::FormatTimeString(71L), "00:01:11");
  EXPECT_STREQ(_String::FormatTimeString(3599L), "00:59:59");
  EXPECT_STREQ(_String::FormatTimeString(3600L), "01:00:00");
  EXPECT_STREQ(_String::FormatTimeString(3600L*101L - 1L), "100:59:59");
}

TEST_F (_hyStringTest, Adler32) {
  const _String referenceString (_String::Random (16L + genrand_int32() % 128L));
  _String mutableReference (referenceString);
  
  long checksum = referenceString.Adler32(),
       hits     = 0L;
  
  for (long loop = 0L; loop < 4096L; loop ++) {
    long idx = genrand_int32() % referenceString.Length();
    mutableReference.setChar(idx, referenceString(idx) + 1);
    if (checksum == mutableReference.Adler32()) {
      hits ++;
    }
    mutableReference.setChar (idx, referenceString(idx));
  }
  
  EXPECT_LE(hits, 2L) << "Too many Adler32 matches for one-character-away strings";
}

TEST_F (_hyStringTest, Tokenize) {
  for (long k = 0L; k < 4096L; k++) {
    
    _String number_string (_String::Random (genrand_int32() % 32, &numbers)),
            letter_string (_String::Random (genrand_int32() % 32, &letters)),
            number_string2 (_String::Random (genrand_int32() % 32, &numbers)),
            joined = number_string & '|' & letter_string & '|' & number_string2;
    
    _List   single_token (number_string.Tokenize ("|"));
    ASSERT_TRUE (single_token.Length () == 1 && * single_token(0) == number_string) << "A.Tokenize (char not in A) != A";

    _List   multiple_tokens (joined.Tokenize ("|"));
    
    ASSERT_EQ (multiple_tokens.Length(), 3UL) << "A|B|C.Tokenize ('|') failed '" <<  joined.getStr() << "'";
    
    ASSERT_TRUE (* multiple_tokens(0) == number_string && * multiple_tokens(2) == number_string2
                 && * multiple_tokens(1) == letter_string) << "Incorrect tokens returned A|B|C.Tokenize ('|') ";

    multiple_tokens = joined.Tokenize (letter_string);
    
    if (letter_string.Length()) {
      ASSERT_EQ (multiple_tokens.Length (), 2UL) << "Expected the ABC.Tokenize (B) to return [A,B] for string'" <<  joined.getStr() << "'";
       ASSERT_TRUE (* multiple_tokens(0) == (number_string & '|') && * multiple_tokens(1) ==  (_String('|') & number_string2))
                    << "Incorrect tokens returned by the ABC.Tokenize (B)";
    } else {
    
      ASSERT_TRUE (multiple_tokens.Length () == 1 && * multiple_tokens(0) == joined) << "A.Tokenize (empty) != A";
    }
  

  }
}

TEST_F (_hyStringTest, toNum) {
    ASSERT_EQ (0.0, empty.toNum()) << "Empty string must convert to 0.";
    ASSERT_EQ (-1.23e10, _String ("-1.23e+10").toNum()) << " '-1.23e+10' conversion to double failed";
    ASSERT_EQ (1.e-8, _String ("0.00000001").toNum()) << " '0.0000001' conversion to double failed";
}


TEST_F (_hyStringTest, Sort) {
    
    for (unsigned long k = 0UL; k < 1024UL; k++) {
        
        _SimpleList index_list;
        
        const   _String random_string (_String::Random (genrand_int32() % 384)),
                sorted (random_string.Sort()),
                sorted_index (random_string.Sort(&index_list));
        
        for (unsigned long i = 0UL; i < random_string.Length(); i++) {
            if (i) {
                ASSERT_GE (sorted(i), sorted (i-1UL)) << "Characters are not sorted following the applicaiton of 'Sort' (no index) to '" << sorted.getStr() << "'";
                ASSERT_GE (sorted_index(i), sorted_index (i-1UL)) << "Characters are not sorted following the applicaiton of 'Sort' (index) to '" << sorted_index.getStr() << "'";
            }
            ASSERT_EQ (sorted_index (i), random_string (index_list(i))) << "The character index at " << i << " is incorrect following the applicaiton of 'Sort' (index) to '"  << sorted_index.getStr() << "'";
        }
        
    }
}

TEST_F (_hyStringTest, Identifiers) {
    // some standard cases
    _String id1 ("hyphy3"),
            id2 ("hyphy3.new"),
            id3 ("0hyphy"),
            id4 ("hy-phy"),
            id5 ("hy_phy"),
            id6 ("hyphy.test.dots"),
            id7 ("hyphy.fail..dots");
    
    ASSERT_TRUE (id1.IsValidIdentifier()) << id1.getStr() << " declared an invalid identifier";
    ASSERT_TRUE (id5.IsValidIdentifier()) << id5.getStr() << " declared an invalid identifier";
    ASSERT_TRUE (id6.IsValidIdentifier()) << id6.getStr() << " declared an invalid identifier";

    ASSERT_FALSE (id3.IsValidIdentifier()) << id3.getStr() << " declared a valid identifier";
    ASSERT_FALSE (id4.IsValidIdentifier()) << id4.getStr() << " declared a valid identifier";
    ASSERT_FALSE (id7.IsValidIdentifier()) << id7.getStr() << " declared a valid identifier";
    ASSERT_FALSE (id2.IsValidIdentifier(false)) << id2.getStr() << " declared a valid identifier (no compounds allowed)";
    ASSERT_FALSE (empty.IsValidIdentifier()) << "An empty string declared a valid identifier";

    
    for (unsigned long k = 0UL; k < 1024UL; k++) {
        const _String randomIdent  = _String::Random (1, &alpha_) & _String::Random (genrand_int32() % 32, &alnum_),
                      randomIdent2 = _String::Random (1, &alpha_) & _String::Random (genrand_int32() % 32, &alnum_),
                      randomIdent3 = _String::Random (1, &alpha_) & _String::Random (genrand_int32() % 32, &alnum_),
                      compound2 = randomIdent & '.' & randomIdent2,
                      compound3 = randomIdent & '.' & randomIdent2  & '.' & randomIdent3,
                      partial1 = randomIdent & _String::Random (1 + genrand_int32() % 32, &punctuation),
                      partial2 = compound3 & _String::Random (1 + genrand_int32() % 32, &punctuation);
        
        
        
        
        ASSERT_TRUE (randomIdent.IsValidIdentifier()) << "IsValidIdentifier failed on '" << randomIdent.getStr() << "'";
        ASSERT_TRUE (compound3.IsValidIdentifier()) << "IsValidIdentifier failed on '" << compound3.getStr() << "'" ;
        ASSERT_FALSE(compound3.IsValidIdentifier(false)) << "IsValidIdentifier (no compounds) failed on '" << compound3.getStr() << "'";
        ASSERT_FALSE(partial1.IsValidIdentifier(false)) << "IsValidIdentifier failed on '" << partial1.getStr() << "'";
        
        ASSERT_EQ (compound2.ShortenVarID(randomIdent), randomIdent2) << "'A.B'.ShortenVarID ('A') failed on '" << compound2.getStr() << "'";
        ASSERT_EQ (compound2.ShortenVarID(randomIdent2), compound2) << "'A.B'.ShortenVarID ('B') failed on '" << compound2.getStr() << "'";
        ASSERT_EQ ((randomIdent&randomIdent2).ShortenVarID(randomIdent), randomIdent&randomIdent2) << "'AB'.ShortenVarID ('A') failed on '" << (randomIdent&randomIdent2).getStr() << "'";
    }
}

TEST_F (_hyStringTest, SpaceFunctions) {
  for (unsigned long k = 0UL; k < 1024UL; k++) {
    unsigned long segments = 1UL + genrand_int32() % 32;
    
    _StringBuffer with_spaces,
                  single_spaces,
                  no_spaces;
    
    bool last_space = false;
    for (unsigned long p = 0UL; p < segments; p++) {
      bool is_space = genrand_real2() > 0.5;
      _String random_s = _String::Random (1L + genrand_int32() % 5, is_space ? &whitespaces : &alnum_);
      if (is_space) {
        with_spaces << random_s;
        if (!last_space) {
          single_spaces << ' ';
        }
      } else {
        no_spaces << random_s;
        single_spaces << random_s;
        with_spaces << random_s;
      }
      last_space = is_space;
    }
    
    ASSERT_EQ (with_spaces.CompressSpaces(),single_spaces) << "Compress spaces failed on '" << with_spaces.getStr() << "' return value '" << with_spaces.CompressSpaces() << "'. Expected '" << single_spaces << "'";
    ASSERT_EQ (single_spaces.CompressSpaces(),single_spaces) << "Compress (should do nothing) spaces failed on '" << single_spaces.getStr() << '\'';
    ASSERT_EQ (with_spaces.KillSpaces (),no_spaces) << "Kill spaces failed on '" << with_spaces.getStr() << '\'';
    ASSERT_EQ (HY_NOT_FOUND, no_spaces.FirstSpaceIndex()) << "Found a space character in a string that should have none";
    ASSERT_EQ (0, no_spaces.FirstNonSpaceIndex());
    ASSERT_EQ (no_spaces.Length() - 1UL, no_spaces.FirstNonSpaceIndex(0, HY_NOT_FOUND, HY_STRING_DIRECTION_BACKWARD));
    
    
  }
}
