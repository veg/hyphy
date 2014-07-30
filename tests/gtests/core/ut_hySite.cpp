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

/*
TEST_F(_StringBufferTest, ConstructorConstCharTest)
{
  // Normal
  const char* test_s = "hyphy";
  _StringBuffer test(test_s);
  EXPECT_EQ(5, test.s_length);
  EXPECT_STREQ("hyphy", test);

  // One Space
  const char* test_s2= " ";
  _StringBuffer test2(test_s2);
  EXPECT_EQ(1, test2.s_length);
  EXPECT_STREQ(" ", test2);

  // Empty
  const char* test_s3 = "";
  _StringBuffer test3(test_s3);
  EXPECT_EQ(0, test3.s_length);
  EXPECT_STREQ(NULL, test3);

  const char* test_s4 = "";
  _StringBuffer test4(test_s4);
  test4 << "hyphy";
  EXPECT_EQ(5, test4.s_length);
  EXPECT_STREQ("hyphy", test4);
}

TEST_F(_StringBufferTest, ConstructorLengthTest)
{
  // Normal
  _StringBuffer test(5);
  EXPECT_EQ(0, test.s_length);
  EXPECT_STREQ("", test);

  _StringBuffer test2(0UL);
  EXPECT_EQ(0, test2.s_length);
  EXPECT_STREQ("", test2);

  _StringBuffer test3(0UL);
  test3 << "hyphy";
  EXPECT_EQ(5, test3.s_length);
  EXPECT_STREQ("hyphy", test3);
}


TEST_F(_StringBufferTest, InitializeTest)
{
  _StringBuffer test(new _String("hyphy"));
  EXPECT_STREQ("hyphy", test);
  EXPECT_EQ(5, test.s_length);
  test.Initialize();
  EXPECT_STREQ(NULL, test);
  EXPECT_EQ(0, test.s_length);
}

TEST_F(_StringBufferTest, ConstructorEmptyTest)
{
  _StringBuffer test;
  EXPECT_EQ(0, test.s_length);
  EXPECT_STREQ("", test);

  _StringBuffer test2;
  test2 << "hyphy";
  EXPECT_EQ(5, test2.s_length);
  EXPECT_STREQ("hyphy", test2);

  _StringBuffer test3;
  test3 << "";
  EXPECT_EQ(0, test3.s_length);
  EXPECT_STREQ("", test3);
}

TEST_F(_StringBufferTest, ConstructorStringBufferTest)
{
  // Empty buffer from empty buffer
  _StringBuffer test;
  _StringBuffer test2(test);
  EXPECT_EQ(0, test2.s_length);
  EXPECT_STREQ("", test2);

  // Empty buffer from empty "" buffer
  _StringBuffer test3;
  test3 << "";
  _StringBuffer test4(test3);
  EXPECT_EQ(0, test4.s_length);
  EXPECT_STREQ("", test4);

  // buffer from normal buffer
  _StringBuffer test5;
  test5 << "hyphy";
  _StringBuffer test6(test5);
  EXPECT_EQ(5, test6.s_length);
  EXPECT_STREQ("hyphy", test6);

  // modify the original, test duplicate
  _StringBuffer test7;
  test7 << "hyphy";
  _StringBuffer test8(test7);
  test7 << "hyphy";
  EXPECT_EQ(5, test8.s_length);
  EXPECT_STREQ("hyphy", test8);
  EXPECT_STREQ("hyphyhyphy", test7);

  // modify the duplicate, test original
  _StringBuffer test9;
  test9 << "hyphy";
  _StringBuffer test10(test9);
  test9 << "xxxxx";
  test10 << "yyyyy";
  EXPECT_EQ(10, test9.s_length);
  EXPECT_EQ(10, test10.s_length);
  EXPECT_STREQ("hyphyxxxxx", test9);
  EXPECT_STREQ("hyphyyyyyy", test10);
}

TEST_F(_StringBufferTest, DuplicateTest)
{
  // Modify the original, test duplicate
  _StringBuffer test;
  _StringBuffer test2("hyphy");
  test.Duplicate(&test2);
  test2 << "hyphy";
  EXPECT_EQ(5, test.s_length);
  EXPECT_STREQ("hyphy", test);
  EXPECT_STREQ("hyphyhyphy", test2);

  // modify the duplicate, test the original
  _StringBuffer test3;
  _StringBuffer test4("hyphy");
  test3.Duplicate(&test4);
  test3 << "hyphy";
  EXPECT_EQ(5, test4.s_length);
  EXPECT_EQ(10, test3.s_length);
  EXPECT_STREQ("hyphy", test4);
  EXPECT_STREQ("hyphyhyphy", test3);

  // overwrite content with duplicate
  _StringBuffer test5("hyphy");
  _StringBuffer test6("data");
  test6.Duplicate(&test5);
  EXPECT_EQ(5, test5.s_length);
  EXPECT_EQ(5, test6.s_length);
  EXPECT_STREQ("hyphy", test5);
  EXPECT_STREQ("hyphy", test6);

  // test with a string instead of a string buffer
  _String test7("hyphy");
  _StringBuffer test8;
  test8.Duplicate(&test7);
  EXPECT_EQ(5, test7.s_length);
  EXPECT_EQ(5, test8.s_length);
  EXPECT_STREQ("hyphy", test7);
  EXPECT_STREQ("hyphy", test8);
}

TEST_F(_StringBufferTest, makeDynamicTest)
{
  _StringBuffer test("hyphy");
  _StringBuffer* test2 = (_StringBuffer*)(test.makeDynamic());
  *test2<<"hyphy";
  _StringBuffer test3(*test2);
  test3<<"hyphy";
  EXPECT_STREQ("hyphy", test);
  EXPECT_STREQ("hyphyhyphy", *test2);
  EXPECT_STREQ("hyphyhyphyhyphy", test3);
  delete test2;
}

TEST_F(_StringBufferTest, AppendStringPointerTest)
{
  _StringBuffer test;
  _String* test_s = new _String("hyphy");
  test << test_s;
  EXPECT_STREQ("hyphy", test);
}

TEST_F(_StringBufferTest, AppendStringReferenceTest)
{
  _StringBuffer test;
  _String test_s("hyphy");
  test << test_s;
  EXPECT_STREQ("hyphy", test);

  _StringBuffer test2;
  _String test_s2("");
  test2 << test_s2;
  EXPECT_STREQ("", test2);

  _StringBuffer test3;
  _String test_s3;
  test3 << test_s3;
  EXPECT_STREQ("", test3);
}

TEST_F(_StringBufferTest, AppendCharArrayTest)
{
  _StringBuffer test;
  const char* test_ca = "hyphy";
  test << test_ca;
  EXPECT_STREQ("hyphy", test);

  _StringBuffer test2;
  const char* test_ca2 = "";
  test2 << test_ca2;
  EXPECT_STREQ("", test2);
}

TEST_F(_StringBufferTest, AppendCharTest)
{
  _StringBuffer test;
  const char test_c = 'h';
  test << test_c;
  EXPECT_STREQ("h", test);

  _StringBuffer test2;
  const char test_c2 = '\0';
  test2 << test_c2;
  EXPECT_STREQ("", test2);

  _StringBuffer test3;
  for (int i = 0; i < 1000; i++) {
    test3 << 'h';
  }
  std::string test3_e = std::string(1000, 'h');
  EXPECT_STREQ(test3_e.data(), test3);

}

// There is still a delete error if you pass it a non-dynamic string which is
// later deleted.
// Also, why is test_s not a NULL ptr?
TEST_F(_StringBufferTest, AppendNewInstanceTest)
{
  _StringBuffer test;
  _String* test_s = new _String("hyphy");
  test.appendNewInstance(test_s);
  EXPECT_STREQ("hyphy", test);
}

TEST_F(_StringBufferTest, sanitizeForSQLTest)
{
  _StringBuffer test("AB");
  const char test_c = '\'';
  test.sanitizeForSQLAndAppend(test_c);
  EXPECT_STREQ("AB''", test);
  test.sanitizeForSQLAndAppend('h');
  EXPECT_STREQ("AB''h", test);

  _StringBuffer test2("AB");
  test2.sanitizeForSQLAndAppend(new _String("hyphy\'"));
  EXPECT_STREQ("ABhyphy''", test2);
}

TEST_F(_StringBufferTest, sanitizeForHTMLTest)
{
  _StringBuffer test("AB");
  test.sanitizeForHTMLAndAppend('"');
  EXPECT_STREQ("AB&quot;", test);
  test.sanitizeForHTMLAndAppend('\'');
  EXPECT_STREQ("AB&quot;&apos;", test);
  test.sanitizeForHTMLAndAppend('<');
  EXPECT_STREQ("AB&quot;&apos;&lt;", test);
  test.sanitizeForHTMLAndAppend('>');
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;", test);
  test.sanitizeForHTMLAndAppend('&');
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;", test);
  test.sanitizeForHTMLAndAppend('^');
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;^", test);

  _StringBuffer test2("AB");
  test2.sanitizeForHTMLAndAppend(new _String("hyphy\"\'<>&^hyphy"));
  EXPECT_STREQ("ABhyphy&quot;&apos;&lt;&gt;&amp;^hyphy", test2);
}

TEST_F(_StringBufferTest, sanitizeTest)
{
  _StringBuffer test("AB");
  test.sanitizeAndAppend('<');
  EXPECT_STREQ("AB<", test);

  _StringBuffer test2("AB");
  test2.sanitizeAndAppend('\n');
  EXPECT_STREQ("AB\\n", test2);
  test2.sanitizeAndAppend('\t');
  EXPECT_STREQ("AB\\n\\t", test2);
  test2.sanitizeAndAppend("\"");
  EXPECT_STREQ("AB\\n\\t\\\"", test2);
  test2.sanitizeAndAppend("\\");
  EXPECT_STREQ("AB\\n\\t\\\"\\\\", test2);
  test2.sanitizeAndAppend("!");
  EXPECT_STREQ("AB\\n\\t\\\"\\\\!", test2);

  _StringBuffer test3("AB");
  test3.sanitizeAndAppend(new _String("hyphy\n\t\"\\!hyphy"));
  EXPECT_STREQ("ABhyphy\\n\\t\\\"\\\\!hyphy", test3);
}

TEST_F(_StringBufferTest, sanitizePSTest)
{
  _StringBuffer test("AB");
  test.sanitizeForPostScriptAndAppend('(');
  EXPECT_STREQ("AB\\(", test);
  test.sanitizeForPostScriptAndAppend(')');
  EXPECT_STREQ("AB\\(\\)", test);
  test.sanitizeForPostScriptAndAppend('%');
  EXPECT_STREQ("AB\\(\\)\\%", test);
  test.sanitizeForPostScriptAndAppend('@');
  EXPECT_STREQ("AB\\(\\)\\%@", test);

  _StringBuffer test2("AB");
  test2.sanitizeForPostScriptAndAppend(new _String("hyphy()%@hyphy"));
  EXPECT_STREQ("ABhyphy\\(\\)\\%@hyphy", test2);
}

TEST_F(_StringBufferTest, sanitizeRETest)
{
  _StringBuffer test("AB");
  test.sanitizeForRegExAndAppend('[');
  EXPECT_STREQ("AB\\[", test);
  test.sanitizeForRegExAndAppend('^');
  EXPECT_STREQ("AB\\[\\^", test);
  test.sanitizeForRegExAndAppend('$');
  EXPECT_STREQ("AB\\[\\^\\$", test);
  test.sanitizeForRegExAndAppend('.');
  EXPECT_STREQ("AB\\[\\^\\$\\.", test);
  test.sanitizeForRegExAndAppend('|');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|", test);
  test.sanitizeForRegExAndAppend('?');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?", test);
  test.sanitizeForRegExAndAppend('*');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*", test);
  test.sanitizeForRegExAndAppend('+');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+", test);
  test.sanitizeForRegExAndAppend('(');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(", test);
  test.sanitizeForRegExAndAppend(')');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)", test);
  test.sanitizeForRegExAndAppend('\\');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\", test);
  test.sanitizeForRegExAndAppend('@');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\@", test);

  _StringBuffer test2("AB");
  test2.sanitizeForRegExAndAppend(new _String("hyphy[^$.|?*+()\\@hyphy"));
  EXPECT_STREQ("ABhyphy\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\@hyphy", test2);
}

TEST_F(_StringBufferTest, AppendAnAssignmentToBufferTest)
{
  _StringBuffer test;
  test << "hyphy";
  test.appendAnAssignmentToBuffer(new _String("12"), new _String("12"), false, false, false);
  EXPECT_STREQ("hyphy12=12;\n", test);

  _StringBuffer test2;
  test2 << "hyphy";
  _String* test2_s = new _String("12");
  test2.appendAnAssignmentToBuffer(new _String("12"), test2_s, true, false, false);
  EXPECT_STREQ("hyphy12=12;\n", test2);

  _StringBuffer test3;
  test3 << "hyphy";
  _String* test3_s = new _String("12");
  test3.appendAnAssignmentToBuffer(new _String("12"), test3_s, true, true, false);
  EXPECT_STREQ("hyphy12=\"12\";\n", test3);

  _StringBuffer test4;
  test4 << "hyphy";
  _String* test4_s = new _String("12");
  test4.appendAnAssignmentToBuffer(new _String("12"), test4_s, true, true, true);
  EXPECT_STREQ("hyphy12:=\"12\";\n", test4);
}

TEST_F(_StringBufferTest, AppendSubstringTest)
{
  _StringBuffer test(new _String("hyphy"));
  _StringBuffer test2;
  test2.appendSubstring(test, 1, 2);
  EXPECT_STREQ("yp", test2);

  _StringBuffer test3;
  test3.appendSubstring(test, 1, 1);
  EXPECT_STREQ("y", test3);

  _StringBuffer test4;
  test4.appendSubstring(test, 1, 0);
  EXPECT_STREQ("", test4);

  _StringBuffer test5;
  test5.appendSubstring(test, -1, 2);
  EXPECT_STREQ("hyp", test5);

  _StringBuffer test6;
  test6.appendSubstring(test, 1, 12);
  EXPECT_STREQ("yphy", test6);
}
*/

} // namespace


#endif
