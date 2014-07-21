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

#ifndef _HY_STRING_BUFFER
#define _HY_STRING_BUFFER

#include <iostream>
#include "gtest/gtest.h"
#include "hy_string_buffer.h"

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Range;
using ::testing::Combine;

namespace {

// The fixture for testing class Foo.
class _StringBufferTest : public ::testing::Test {

  protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _StringBufferTest() {
      // You can do set-up work for each test here.
    }

    virtual ~_StringBufferTest() {
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


/******************************************/
TEST_F(_StringBufferTest, ConstructorStringTest)
{
  // Normal
  _StringBuffer test(new _String("hyphy"));
  EXPECT_EQ(5, test.sLength);
  EXPECT_STREQ("hyphy", test);

  // One Space
  _StringBuffer test2(new _String(" "));
  EXPECT_EQ(1, test2.sLength);
  EXPECT_STREQ(" ", test2);

  // Empty
  // I'm not sure I agree with this expectation...
  _StringBuffer test3(new _String(""));
  EXPECT_EQ(0, test3.sLength);
  EXPECT_STREQ(NULL, test3);

  _StringBuffer test4(new _String(""));
  test4 << "hyphy";
  EXPECT_EQ(5, test4.sLength);
  EXPECT_STREQ("hyphy", test4);
}

/******************************************/
TEST_F(_StringBufferTest, ConstructorConstCharTest)
{
  // Normal
  const char* test_s = "hyphy";
  _StringBuffer test(test_s);
  EXPECT_EQ(5, test.sLength);
  EXPECT_STREQ("hyphy", test);

  // One Space
  const char* test_s2= " ";
  _StringBuffer test2(test_s2);
  EXPECT_EQ(1, test2.sLength);
  EXPECT_STREQ(" ", test2);

  // Empty
  const char* test_s3 = "";
  _StringBuffer test3(test_s3);
  EXPECT_EQ(0, test3.sLength);
  EXPECT_STREQ(NULL, test3);

  const char* test_s4 = "";
  _StringBuffer test4(test_s4);
  test4 << "hyphy";
  EXPECT_EQ(5, test4.sLength);
  EXPECT_STREQ("hyphy", test4);
}

/******************************************/
TEST_F(_StringBufferTest, ConstructorLengthTest)
{
  // Normal
  _StringBuffer test(5);
  EXPECT_EQ(0, test.sLength);
  EXPECT_STREQ("", test);

  _StringBuffer test2(0UL);
  EXPECT_EQ(0, test2.sLength);
  EXPECT_STREQ("", test2);

  _StringBuffer test3(0UL);
  test3 << "hyphy";
  EXPECT_EQ(5, test3.sLength);
  EXPECT_STREQ("hyphy", test3);
}


/******************************************/
TEST_F(_StringBufferTest, InitializeTest)
{
  _StringBuffer test(new _String("hyphy"));
  EXPECT_STREQ("hyphy", test);
  EXPECT_EQ(5, test.sLength);
  test.Initialize();
  EXPECT_STREQ(NULL, test);
  EXPECT_EQ(0, test.sLength);
}

/******************************************/
TEST_F(_StringBufferTest, ConstructorEmptyTest)
{
  _StringBuffer test;
  EXPECT_EQ(0, test.sLength);
  EXPECT_STREQ("", test);

  _StringBuffer test2;
  test2 << "hyphy";
  EXPECT_EQ(5, test2.sLength);
  EXPECT_STREQ("hyphy", test2);

  _StringBuffer test3;
  test3 << "";
  EXPECT_EQ(0, test3.sLength);
  EXPECT_STREQ("", test3);
}

/******************************************/
TEST_F(_StringBufferTest, ConstructorStringBufferTest)
{
  // Empty buffer from empty buffer
  _StringBuffer test;
  _StringBuffer test2(test);
  EXPECT_EQ(0, test2.sLength);
  EXPECT_STREQ("", test2);

  // Empty buffer from empty "" buffer
  _StringBuffer test3;
  test3 << "";
  _StringBuffer test4(test3);
  EXPECT_EQ(0, test4.sLength);
  EXPECT_STREQ("", test4);

  // buffer from normal buffer
  _StringBuffer test5;
  test5 << "hyphy";
  _StringBuffer test6(test5);
  EXPECT_EQ(5, test6.sLength);
  EXPECT_STREQ("hyphy", test6);

  // modify the original, test duplicate
  _StringBuffer test7;
  test7 << "hyphy";
  _StringBuffer test8(test7);
  test7 << "hyphy";
  EXPECT_EQ(5, test8.sLength);
  EXPECT_STREQ("hyphy", test8);
  EXPECT_STREQ("hyphyhyphy", test7);

  // modify the duplicate, test original
  _StringBuffer test9;
  test9 << "hyphy";
  _StringBuffer test10(test9);
  test9 << "xxxxx";
  test10 << "yyyyy";
  EXPECT_EQ(10, test9.sLength);
  EXPECT_EQ(10, test10.sLength);
  EXPECT_STREQ("hyphyxxxxx", test9);
  EXPECT_STREQ("hyphyyyyyy", test10);
}

/******************************************/
TEST_F(_StringBufferTest, DuplicateTest)
{
  // Modify the original, test duplicate
  _StringBuffer test;
  _StringBuffer test2("hyphy");
  test.Duplicate(&test2);
  test2 << "hyphy";
  EXPECT_EQ(5, test.sLength);
  EXPECT_STREQ("hyphy", test);
  EXPECT_STREQ("hyphyhyphy", test2);

  // modify the duplicate, test the original
  _StringBuffer test3;
  _StringBuffer test4("hyphy");
  test3.Duplicate(&test4);
  test3 << "hyphy";
  EXPECT_EQ(5, test4.sLength);
  EXPECT_EQ(10, test3.sLength);
  EXPECT_STREQ("hyphy", test4);
  EXPECT_STREQ("hyphyhyphy", test3);

  // overwrite content with duplicate
  _StringBuffer test5("hyphy");
  _StringBuffer test6("data");
  test6.Duplicate(&test5);
  EXPECT_EQ(5, test5.sLength);
  EXPECT_EQ(5, test6.sLength);
  EXPECT_STREQ("hyphy", test5);
  EXPECT_STREQ("hyphy", test6);

  // test with a string instead of a string buffer
  _String test7("hyphy");
  _StringBuffer test8;
  test8.Duplicate(&test7);
  EXPECT_EQ(5, test7.sLength);
  EXPECT_EQ(5, test8.sLength);
  EXPECT_STREQ("hyphy", test7);
  EXPECT_STREQ("hyphy", test8);
}

/******************************************/
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

/******************************************/
TEST_F(_StringBufferTest, AppendStringPointerTest)
{
  _StringBuffer test;
  _String* test_s = new _String("hyphy");
  test << test_s;
  EXPECT_STREQ("hyphy", test);
}

/******************************************/
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

/******************************************/
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

/******************************************/
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
}

/******************************************/
// There is still a delete error if you pass it a non-dynamic string which is
// later deleted.
// Also, why is test_s not a NULL ptr?
TEST_F(_StringBufferTest, AppendNewInstanceTest)
{
  _StringBuffer test;
  _String* test_s = new _String("hyphy");
  test.AppendNewInstance(test_s);
  EXPECT_STREQ("hyphy", test);
  EXPECT_STREQ("", *test_s);
}

/******************************************/
TEST_F(_StringBufferTest, sanitizeForSQLTest)
{
  _StringBuffer test3("AB");
  const char test_c = '\'';
  test3.sanitizeForSQLAndAppend(test_c);
  EXPECT_STREQ("AB''", test3);
}

/******************************************/
TEST_F(_StringBufferTest, sanitizeForHTMLTest)
{
  _StringBuffer test5("AB");
  test5.sanitizeForHTMLAndAppend('"');
  EXPECT_STREQ("AB&quot;", test5);
  test5.sanitizeForHTMLAndAppend('\'');
  EXPECT_STREQ("AB&quot;&apos;", test5);
  test5.sanitizeForHTMLAndAppend('<');
  EXPECT_STREQ("AB&quot;&apos;&lt;", test5);
  test5.sanitizeForHTMLAndAppend('>');
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;", test5);
  test5.sanitizeForHTMLAndAppend('&');
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;", test5);
  test5.sanitizeForHTMLAndAppend('^');
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;^", test5);
}

/******************************************/
TEST_F(_StringBufferTest, sanitizeTest)
{
  _StringBuffer test("AB");
  test.sanitizeAndAppend('<');
  EXPECT_STREQ("AB<", test);

  _StringBuffer test4("AB");
  test4.sanitizeAndAppend('\n');
  EXPECT_STREQ("AB\\n", test4);
  test4.sanitizeAndAppend('\t');
  EXPECT_STREQ("AB\\n\\t", test4);
  test4.sanitizeAndAppend("\"");
  EXPECT_STREQ("AB\\n\\t\\\"", test4);
  test4.sanitizeAndAppend("\\");
  EXPECT_STREQ("AB\\n\\t\\\"\\\\", test4);
  test4.sanitizeAndAppend("!");
  EXPECT_STREQ("AB\\n\\t\\\"\\\\!", test4);
}

/******************************************/
TEST_F(_StringBufferTest, sanitizePSTest)
{
  _StringBuffer test2("AB");
  test2.sanitizeForPostScriptAndAppend('(');
  EXPECT_STREQ("AB\\(", test2);
  test2.sanitizeForPostScriptAndAppend(')');
  EXPECT_STREQ("AB\\(\\)", test2);
  test2.sanitizeForPostScriptAndAppend('%');
  EXPECT_STREQ("AB\\(\\)\\%", test2);
}

/******************************************/
TEST_F(_StringBufferTest, sanitizeRETest)
{
  _StringBuffer test6("AB");
  test6.sanitizeForRegExAndAppend('[');
  EXPECT_STREQ("AB\\[", test6);
  test6.sanitizeForRegExAndAppend('^');
  EXPECT_STREQ("AB\\[\\^", test6);
  test6.sanitizeForRegExAndAppend('$');
  EXPECT_STREQ("AB\\[\\^\\$", test6);
  test6.sanitizeForRegExAndAppend('.');
  EXPECT_STREQ("AB\\[\\^\\$\\.", test6);
  test6.sanitizeForRegExAndAppend('|');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|", test6);
  test6.sanitizeForRegExAndAppend('?');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?", test6);
  test6.sanitizeForRegExAndAppend('*');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*", test6);
  test6.sanitizeForRegExAndAppend('+');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+", test6);
  test6.sanitizeForRegExAndAppend('(');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(", test6);
  test6.sanitizeForRegExAndAppend(')');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)", test6);
  test6.sanitizeForRegExAndAppend('\\');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\", test6);
  test6.sanitizeForRegExAndAppend('@');
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\@", test6);
}

/******************************************/
TEST_F(_StringBufferTest, EscapeAndAppendCharTest)
{
  _StringBuffer test("AB");
  test.EscapeAndAppend('<',HY_ESCAPE_NORMAL);
  EXPECT_STREQ("AB<", test);

  _StringBuffer test2("AB");
  test2.EscapeAndAppend('(',HY_ESCAPE_POSTSCRIPT);
  EXPECT_STREQ("AB\\(", test2);
  test2.EscapeAndAppend(')',HY_ESCAPE_POSTSCRIPT);
  EXPECT_STREQ("AB\\(\\)", test2);
  test2.EscapeAndAppend('%',HY_ESCAPE_POSTSCRIPT);
  EXPECT_STREQ("AB\\(\\)\\%", test2);

  _StringBuffer test3("AB");
  const char test_c = '\'';
  test3.EscapeAndAppend(test_c,HY_ESCAPE_SQLITE);
  EXPECT_STREQ("AB''", test3);

  _StringBuffer test4("AB");
  test4.EscapeAndAppend('\n',HY_ESCAPE_UNUSED);
  EXPECT_STREQ("AB\\n", test4);
  test4.EscapeAndAppend('\t',HY_ESCAPE_UNUSED);
  EXPECT_STREQ("AB\\n\\t", test4);
  test4.EscapeAndAppend("\"",HY_ESCAPE_UNUSED);
  EXPECT_STREQ("AB\\n\\t\\\"", test4);
  test4.EscapeAndAppend("\\",HY_ESCAPE_UNUSED);
  EXPECT_STREQ("AB\\n\\t\\\"\\\\", test4);
  test4.EscapeAndAppend("!",HY_ESCAPE_UNUSED);
  EXPECT_STREQ("AB\\n\\t\\\"\\\\!", test4);

  _StringBuffer test5("AB");
  test5.EscapeAndAppend('"',HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;", test5);
  test5.EscapeAndAppend('\'',HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;&apos;", test5);
  test5.EscapeAndAppend('<',HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;&apos;&lt;", test5);
  test5.EscapeAndAppend('>',HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;", test5);
  test5.EscapeAndAppend('&',HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;", test5);
  test5.EscapeAndAppend('^',HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;^", test5);

  _StringBuffer test6("AB");
  test6.EscapeAndAppend('[',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[", test6);
  test6.EscapeAndAppend('^',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^", test6);
  test6.EscapeAndAppend('$',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$", test6);
  test6.EscapeAndAppend('.',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.", test6);
  test6.EscapeAndAppend('|',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|", test6);
  test6.EscapeAndAppend('?',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?", test6);
  test6.EscapeAndAppend('*',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*", test6);
  test6.EscapeAndAppend('+',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+", test6);
  test6.EscapeAndAppend('(',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(", test6);
  test6.EscapeAndAppend(')',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)", test6);
  test6.EscapeAndAppend('\\',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\", test6);
  test6.EscapeAndAppend('@',HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\@", test6);
}

/******************************************/
TEST_F(_StringBufferTest, EscapeAndAppendStringTest)
{
  _StringBuffer test("AB");
  test.EscapeAndAppend(new _String("<>()!"),HY_ESCAPE_NORMAL);
  EXPECT_STREQ("AB<>()!", test);

  _StringBuffer test2("AB");
  test2.EscapeAndAppend(new _String("()%@!"),HY_ESCAPE_POSTSCRIPT);
  EXPECT_STREQ("AB\\(\\)\\%@!", test2);

  _StringBuffer test3("AB");
  test3.EscapeAndAppend(new _String("\'\'@#"),HY_ESCAPE_SQLITE);
  EXPECT_STREQ("AB''''@#", test3);

  _StringBuffer test4("AB");
  test4.EscapeAndAppend(new _String("\n\t\"\\@!"),HY_ESCAPE_UNUSED);
  EXPECT_STREQ("AB\\n\\t\\\"\\\\@!", test4);

  _StringBuffer test5("AB");
  test5.EscapeAndAppend(new _String("\"\'<>&^#!"),HY_ESCAPE_HTML);
  EXPECT_STREQ("AB&quot;&apos;&lt;&gt;&amp;^#!", test5);

  _StringBuffer test6("AB");
  test6.EscapeAndAppend(new _String("[]^$.|?*+()\\@!#"),HY_ESCAPE_REGEXP);
  EXPECT_STREQ("AB\\[]\\^\\$\\.\\|\\?\\*\\+\\(\\)\\\\@!#", test6);
}

/******************************************/
TEST_F(_StringBufferTest, AppendAnAssignmentToBufferTest)
{
  _StringBuffer test;
  test << "hyphy";
  test.AppendAnAssignmentToBuffer(new _String("12"), new _String("12"), false, false, false);
  EXPECT_STREQ("hyphy12=12;\n", test);

  _StringBuffer test2;
  test2 << "hyphy";
  _String* test2_s = new _String("12");
  test2.AppendAnAssignmentToBuffer(new _String("12"), test2_s, true, false, false);
  EXPECT_STREQ("hyphy12=12;\n", test2);
  EXPECT_STREQ("", *test2_s);

  _StringBuffer test3;
  test3 << "hyphy";
  _String* test3_s = new _String("12");
  test3.AppendAnAssignmentToBuffer(new _String("12"), test3_s, true, true, false);
  EXPECT_STREQ("hyphy12=\"12\";\n", test3);
  EXPECT_STREQ("", *test3_s);

  _StringBuffer test4;
  test4 << "hyphy";
  _String* test4_s = new _String("12");
  test4.AppendAnAssignmentToBuffer(new _String("12"), test4_s, true, true, true);
  EXPECT_STREQ("hyphy12:=\"12\";\n", test4);
  EXPECT_STREQ("", *test4_s);
}

/******************************************/
/*
TEST_F(_StringBufferTest, AppendVariableValueTest)
{
  _StringBuffer test;
  _String test_s("hyphy");
  //_SimpleList test_l;
  //test_l << "one";
  //test.AppendVariableValueAVL(test_s, test_l);
}
*/

} // namespace


#endif
