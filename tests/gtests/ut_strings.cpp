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

//#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include "ut_strings.h"


#include "hy_strings.h"
#include "list.h"
#include "simplelist.h"

#include "baseobj.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Range;
using ::testing::Combine;

namespace
{

// The fixture for testing class Foo.
class StringTest : public ::testing::Test
{
protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    StringTest() {
        // You can do set-up work for each test here.
    }

    virtual ~StringTest() {
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


/* 20110825: SLKP made this a global string.*/
_String globalTest1 ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");

/******************************************/

TEST_F(StringTest,InitializeTest)
{
    _String test("hyphy"),
            empty;

    test.Initialize();
    empty.Initialize();

    EXPECT_EQ (test.sData, (char*)NULL);
    EXPECT_EQ (test.sLength, 0);

    EXPECT_EQ (empty.sData, (char*)NULL);
    EXPECT_EQ (empty.sLength, 0);
}

/******************************************/

TEST_F(StringTest,DuplicateTest)
{
    /* 20110825: SLKP code coverage complete

        - added an extreme case (empty string)
        - changed ASSERT to EXPECT
    */
    _String test    = _String ("hyphy"),
            empty;
    _String dupe1,
            dupe2;

    dupe1.Duplicate(&test);
    EXPECT_STREQ(test.getStr(), dupe1.getStr());

    dupe2.Duplicate(&empty);
    EXPECT_EQ (dupe2.sData, (char*)NULL);
    EXPECT_EQ (dupe2.sLength, 0);
}

/******************************************/

TEST_F(StringTest,DuplicateErasingTest)
{
    /* 20110825: SLKP code coverage complete

             - changed ASSERT to EXPECT
             - added a test to ensure that the pointer
               to the string got reallocated after
               DuplicateErasing is called
    */
    _String test = _String ("hyphy"),
            dupe = _String ("old_hyphy");

    char*   oldSData = dupe.sData;
    dupe.DuplicateErasing(&test);

    test = empty;

    EXPECT_STREQ ("hyphy", dupe.getStr());
}

/******************************************/


TEST_F(StringTest,makeDynamicTest)
{
    //What is the difference between this and dupicate?

    /* 20110825: SLKP code coverage complete
          The idea of make dynamic is to convert a stack object into a heap object, see code below
     */

    _String stackString (globalTest1, 5, -1), // this helps test one of the constructors
            *heapString = (_String*)stackString.makeDynamic();

    EXPECT_STREQ (stackString.getStr(), heapString->getStr());
    stackString = empty; // overwrite the stack object
    EXPECT_EQ (heapString->sLength, globalTest1.sLength-5);
    DeleteObject (heapString);

}

/******************************************/

TEST_F(StringTest,getCharTest)
{

    /* 20110825: SLKP code coverage complete */

    _String test (globalTest1);
    EXPECT_EQ('e', test.getChar(5));

    //Default return is 0
    _String test2 = empty;
    EXPECT_EQ(0, test2.getChar(5));

}

/******************************************/


TEST_F(StringTest,setCharTest)
{
    _String test (globalTest1);
    test.setChar(5,'d');
    EXPECT_EQ('d', test.getChar(5));

    //Should be 0
    _String test2 = empty;
    test2.setChar(5,'d');
    EXPECT_EQ('\0', test2.getChar(5));
}

/******************************************/

TEST_F(StringTest,CopyDynamicStringTest)
{

    /* 20110825: SLKP code coverage complete

          - added a test case when the source string has object counter > 1, i.e. it can't be
            simply moved to the destination but has to be copied

    */


    _String* test = new _String ("hyphy");
    _String dupe = _String("old_hyphy");
    dupe.CopyDynamicString(test);
    EXPECT_STREQ("hyphy", dupe.getStr());

    _String* test2 = new _String ("beavis");
    test2->AddAReference(); // add 1 to object counter
    _String dupe2, dupe3, blank;
    dupe2.CopyDynamicString(test2);
    blank.AddAReference();
    dupe3.CopyDynamicString(&blank);
    EXPECT_STREQ("beavis", dupe2.getStr());
    EXPECT_STREQ("beavis", test2->getStr());
    EXPECT_STREQ("", dupe3.getStr());
    DeleteObject (test2);
}

/******************************************/


TEST_F(StringTest, LengthTest)
{
    _String test = new _String ("You're asking me to run MCMC without reporting any tests.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ(99, test.Length());

    _String test2 = new _String ("");
    EXPECT_EQ(0, test2.Length());
}

TEST_F(StringTest,InsertTest)
{
    _String test = "AAGGCCTTA";
    _String expected_test = "CAAGGCCTTA";
    test.Insert('C',0);
    EXPECT_STREQ(expected_test.getStr(), test.getStr());

    _String test2 = "";
    _String expected_test2 = "";
    test2.Insert('C',5);
    EXPECT_STREQ(expected_test2.getStr(), test2.getStr());

    test = "AAGGCCTTA";
    expected_test = "AAGGCCTTAC";
    test.Insert('C',-11);
    EXPECT_STREQ(expected_test.getStr(), test.getStr());

}

TEST_F(StringTest,DeleteTest)
{
    _String test = "AAGGCCTTA";
    test.Delete(3,4);
    EXPECT_STREQ("AAGCTTA", test.getStr());

    /*  //ERROR: Program crashes if you go beyond end of string.
     *
     *  _String test2 = "AAGG";
     *  test2.Delete(1,14);
     *  EXPECT_STREQ("A", test2.getStr());
     *
     */

    _String test2 = "AAGG";
    test2.Delete(1,-5);
    EXPECT_STREQ("A", test2.getStr());

    _String test3 = "AAGG";
    test3.Delete(-1,-5);
    EXPECT_STREQ("", test3.getStr());

}

TEST_F(StringTest,AppendNewInstanceTest)
{
    //TODO: Debug for coverage
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String test = orig&to_append;
    EXPECT_STREQ(expected.getStr(), test.getStr());

    orig = _String ("");
    to_append = _String("");
    expected = _String("");
    test = orig&to_append;
    EXPECT_STREQ(expected.getStr(), test.getStr());

}


TEST_F(StringTest,EscapeAndAppendTest)
{
    _String result = _String("AAGG");
    _String expected = _String("AAGG\\\\(&lt;\\\[''");


    result.EscapeAndAppend('\\',2);
    result.EscapeAndAppend('(',1);
    result.EscapeAndAppend('<',4);
    result.EscapeAndAppend('[',5);
    result.EscapeAndAppend('\'',2);

    EXPECT_EQ(expected.getStr()[expected.sLength - 1], 
              result.getStr()[expected.sLength - 1]);

}

TEST_F(StringTest,EscapeAndAppend2Test)
{
    _String str = _String("AAGG");
    _String append = _String("\"AAGG\"");
    _String expected = _String("AAGG&quot;AAGG&quot;");

    str.EscapeAndAppend(append,4);

    EXPECT_EQ(expected.getStr()[expected.sLength - 1], 
              str.getStr()[expected.sLength - 1]);
}

TEST_F(StringTest,FinalizeTest)
{
    _String orig = _String ("hyphy");
    orig.Finalize();
    EXPECT_EQ(0, orig[orig.Length()]);
}

TEST_F(StringTest,getStrTest)
{
    _String test = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_STREQ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n", test.getStr());

    _String test2 = _String ("");
    EXPECT_STREQ("", test2.getStr());

}

TEST_F(StringTest,ChopTest)
{

    _String test = globalTest1;
    _String expected = _String("'re asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r1 = test.Chop(0,2);
    EXPECT_STREQ(expected.getStr(), r1.getStr());

    _String test2 = globalTest1;
    _String expected2 = _String("'re asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r2 = test2.Chop(-1,2);
    EXPECT_STREQ(expected2.getStr(), r2.getStr());

    _String test3 = _String ("");
    _String r3 = test3.Chop(0,2);
    EXPECT_STREQ("", r3.getStr());

    _String test4 = globalTest1;
    _String r4 = test4.Chop(5,4);
    EXPECT_STREQ("", r4.getStr());

    _String test5 = globalTest1;
    _String expected5 = _String("Y're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r5 = test5.Chop(1,2);
    EXPECT_STREQ(expected5.getStr(), r5.getStr());

    _String test6 = globalTest1;
    _String r6 = test6.Chop(0,-1);
    EXPECT_STREQ("", r6.getStr());

    //Error: Should return empty, but it returns "AA"
    _String test7 = _String ("ABBA");
    _String substr7 = test7.Chop(11,2);
    EXPECT_STREQ("", substr7.getStr());
 

    //Error: Memory allocation error
    //_String test3 = _String ("ABBA");
    //_String substr3 = test3.Chop(2,20);
    //EXPECT_STREQ("BA", substr3.getStr());
   
}

TEST_F(StringTest,CutTest)
{
    _String test = globalTest1;
    _String substr = test.Cut(0,2);
    EXPECT_STREQ("You", substr.getStr());

    _String test2 = _String ("AABBCC");
    _String substr2 = test2.Cut(4,12);
    EXPECT_STREQ("CC", substr2.getStr());

    _String test3 = _String ("");
    _String substr3 = test3.Cut(4,13);
    EXPECT_STREQ("", substr3.getStr());

    _String test4 = _String ("AABBCC");
    _String substr4 = test4.Cut(-1,-1);
    EXPECT_STREQ("AABBCC", substr4.getStr());

    _String test5 = _String ("AABBCC");
    _String substr5 = test5.Cut(5,4);
    EXPECT_STREQ("", substr5.getStr());

    _String test6 = _String ("");
    _String substr6 = test6.Cut(4,13);
    EXPECT_STREQ("", substr6.getStr());

}

TEST_F(StringTest,FlipTest)
{
    _String result = _String ("ABC");
    result.Flip();
    EXPECT_STREQ("CBA",result.getStr());
}

TEST_F(StringTest,Adler32Test)
{
    //TODO: Coverage
    _String result = new _String ("Wikipedia");
    EXPECT_EQ(300286872, result.Adler32());
}

TEST_F(StringTest,TrimTest)
{
    _String test = globalTest1; 
    test.Trim(7,12);
    EXPECT_STREQ("asking", test.getStr());

    _String test2("");
    test2.Trim(7,12);
    EXPECT_STREQ("", test2.getStr());

    //Is there a reason from should act in such a manner?
    _String test3 = _String ("AABBCC");
    test3.Trim(13,-1);
    EXPECT_STREQ("C", test3.getStr());

    _String test4 = _String ("AABBCC");
    test4.Trim(0,-1);
    EXPECT_STREQ("AABBCC", test4.getStr());

    _String test5 = _String ("AABBCC");
    test5.Trim(-1,30, true);
    EXPECT_STREQ("AABBCC", test5.getStr());

}

TEST_F(StringTest,FirstNonSpaceIndexTest)
{
    //Error: curious that it would give me an 8
    _String test = _String ("    hyphy");
    EXPECT_EQ(4, test.FirstNonSpaceIndex());

    test = _String ("    hyphy");
    EXPECT_EQ(8, test.FirstNonSpaceIndex(-1));

    _String test2 = _String ("    hyphy ");
    EXPECT_EQ(8, test2.FirstNonSpaceIndex(0,-1,-1));

    _String test3 = _String ("");
    EXPECT_EQ(-1, test3.FirstNonSpaceIndex());

}

TEST_F(StringTest,KillSpacesTest)
{
    _String result = _String ("  h   y   p    ");
    _String r2;
    result.KillSpaces(r2);
    EXPECT_STREQ("hyp", r2.getStr());

    _String test2 = _String ("hyp");
    _String result_string2;
    test2.KillSpaces(result_string2);
    EXPECT_STREQ("hyp", result_string2.getStr());

    _String test3 = _String ("");
    _String result_string3;
    test3.KillSpaces(result_string3);
    EXPECT_STREQ("", result_string3.getStr());
}

TEST_F(StringTest,CompressSpacesTest)
{
    _String test = _String ("ABB   and    CCD");
    test.CompressSpaces();
    EXPECT_STREQ("ABB and CCD",test.getStr());

    _String test2 = _String ("AB C");
    test2.CompressSpaces();
    EXPECT_STREQ("AB C",test2.getStr());

    _String test3 = _String ("");
    test3.CompressSpaces();
    EXPECT_STREQ("",test3.getStr());

}

TEST_F(StringTest,FirstSpaceIndexTest)
{
    _String test = _String ("AA BB");
    EXPECT_EQ(2, test.FirstSpaceIndex());

    _String test2 = _String ("");
    EXPECT_EQ(-1, test2.FirstSpaceIndex());

    _String test3 = _String ("A BBB");
    EXPECT_EQ(1, test3.FirstSpaceIndex(-1,-1,-1));

    _String test4 = _String (" A BBB");
    EXPECT_EQ(0, test4.FirstSpaceIndex());
}

TEST_F(StringTest,FirstNonSpaceTest)
{
    _String test = _String ("  AA BB");
    EXPECT_EQ('A', test.FirstNonSpace());

    _String test2 = _String ("AABB ");
    EXPECT_EQ('B', test2.FirstNonSpace(0,-1,-1));
}

TEST_F(StringTest,FirstNonSpace2Test)
{
    _String test = _String ("  AA BB");
    EXPECT_EQ('A', test.FirstNonSpace());

    _String test2 = _String ("AABB ");
    EXPECT_EQ('B', test2.FirstNonSpace(0,-1,-1));
}


TEST_F(StringTest,FindTest)
{
    _String test = _String ("AABBBCCAADDDAABBBEEEFFFF");
    EXPECT_EQ(19, test.Find(_String("EF"), -1));
    EXPECT_EQ(-1, test.Find(_String("EFF"),8,9));
}

TEST_F(StringTest,FindCharTest)
{
    _String test = _String ("AABBBCCAADDDAABBBEEEFFFF");
    EXPECT_EQ(17, test.Find('E', -1));
    EXPECT_EQ(-1, test.Find('E',9,8));

    test = _String ("");
    EXPECT_EQ(-1, test.Find('E'));
}

TEST_F(StringTest,FindKMPTest)
{
    //TODO
}


TEST_F(StringTest,FindEndOfIdentTest)
{
    _String test = _String ("iden12&iden34");
    EXPECT_EQ(5, test.FindEndOfIdent(0,-1,'.'));

    _String test2 = _String ("iden12");
    EXPECT_EQ(5, test2.FindEndOfIdent(0,-1,'.'));

    //Error, returns -2
    _String test3 = _String ("");
    EXPECT_EQ(-1, test3.FindEndOfIdent(-1,-1,'.'));

    _String test4 = _String ("iden12&iden34");
    EXPECT_EQ(5, test4.FindEndOfIdent(0,-1,'.'));

    _String test5 = _String ("iden12__&iden34");
    EXPECT_EQ(5, test5.FindEndOfIdent(0,-1,'.'));
}

TEST_F(StringTest,FindAnyCaseTest)
{
    _String result = _String ("AABBCCDD");
    EXPECT_EQ(2, result.FindAnyCase("BBcCDD",-1));

    _String test = _String ("AABBCCDD");
    EXPECT_EQ(-1, test.FindAnyCase("cBcCDD"));

    _String test2 = _String ("");
    EXPECT_EQ(-1, test2.FindAnyCase("cBcCDD"));

    test2 = _String("AABBCCDD"); 
    EXPECT_EQ(-1, test2.FindAnyCase("AABCC", 6,5));

    test2 = _String("AABBCCDD"); 
    EXPECT_EQ(-1, test2.FindAnyCase("AABCC", 5,6));
}

TEST_F(StringTest,ContainsSubstringTest)
{

    _String test = _String ("AABBCCDD");
    _String substr = _String ("CC");
    EXPECT_EQ(true, test.ContainsSubstring(substr));

    _String test2 = _String ("AABBCCDD");
    _String substr2 = _String ("cC");
    EXPECT_FALSE(test2.ContainsSubstring(substr2));

    _String test3 = _String ("AABBCCDD");
    _String substr3 = _String ("");
    EXPECT_EQ(true, test3.ContainsSubstring(substr3));

    _String test4 = _String ("");
    _String substr4 = _String ("AABCC");
    EXPECT_EQ(false, test4.ContainsSubstring(substr4));

    //ERROR: returns true
    _String test5 = _String ("AA");
    _String substr5 = _String ("AABCC");
    EXPECT_EQ(false, test5.ContainsSubstring(substr5));

}

TEST_F(StringTest,FindBackwardsTest)
{
    _String test = _String ("AABBCCDD");
    EXPECT_EQ(-1, test.FindBackwards("DC",0,3));

    _String test2 = _String ("AABBCCDD");
    EXPECT_EQ(5, test2.FindBackwards("CD",-1,-1));

    _String test3 = _String ("");
    EXPECT_EQ(-1, test3.FindBackwards("CD",-1,-1));

    _String test4 = _String ("AABBCCDD");
    EXPECT_EQ(-1, test4.FindBackwards("CD",5,4));

    _String test5 = _String ("AABBCCDD");
    EXPECT_EQ(-1, test5.FindBackwards("CDAACC",8,9));

}

TEST_F(StringTest,FindBinaryTest)
{
    _String haystack = _String ("AABBCDDD");
    char needle = 'C';

    long loc = haystack.FindBinary(needle);

    EXPECT_EQ(4,loc);


    _String haystack2 = _String ("AABBCDDD");
    char needle2 = 'F';

    long loc2 = haystack2.FindBinary(needle2);

    EXPECT_EQ(-1,loc2);

}

TEST_F(StringTest,EqualsTest)
{
    _String* test = new _String ("AABBCCDD");
    _String* r2 = new _String ("AABBCCDD");
    EXPECT_EQ(true, test->Equal(r2));

    _String test2 = _String("AADCC");
    EXPECT_EQ(false, test2.Equal(r2));

    _String test3("AABBCCDD");
    _String r3("AABBCCDD");
    EXPECT_TRUE(test3==r3);

    delete test;
    delete r2;
}

TEST_F(StringTest,CompareTest)
{
    //house precedes household
    //Household precedes house
    //composer precedes computer
    //H2O precedes HOTEL

    //ERROR: This returns true
    _String result = _String ("household");
    _String* substr = new _String ("house");

    _String result2 = _String ("household");
    _String* substr2 = new _String ("household");
    EXPECT_EQ(0, result2.Compare(substr2));
}

TEST_F(StringTest,EqualWithWildCharTest)
{

    _String test = _String ("AABBCCDD");
    _String* t = new _String ("EEBBCCDD");

    EXPECT_EQ(true, test.EqualWithWildChar(t, 'E'));

    _String test2 = _String ("AAFBCCDD");
    EXPECT_EQ(false, test2.EqualWithWildChar(t, 'E'));

    _String test3 = _String ("EEBBBCDD");
    EXPECT_EQ(false, test3.EqualWithWildChar(t, 'B'));

    _String test4 = _String ("EEBBBCDB");
    EXPECT_EQ(false, test4.EqualWithWildChar(t, 'B'));

    delete t;
}

TEST_F(StringTest,GreaterTest)
{
    //TODO: Debug Coverage 
    //This is a lexicographic comparison
    _String test = _String ("house");
    _String t = _String ("household");

    EXPECT_EQ(false, test>t);

    //This is a lexicographic comparison
    _String test2 = _String ("");
    _String t2 = _String ("household");

    EXPECT_EQ(false, test2>t2);
}

TEST_F(StringTest,LessTest)
{

    //TODO: Debug Coverage 

    _String test = _String ("house");
    _String t = _String ("household");
    EXPECT_EQ(true, test<t);

    _String test2 = _String ("");
    _String t2 = _String ("");
    EXPECT_EQ(false, test2<t2);

}

TEST_F(StringTest,containsTest)
{

    _String t = _String ("household");
    _String test = _String ("house");
    EXPECT_EQ(true, t.contains(test));

    _String t2 = _String ("");
    _String test2 = _String ("");
    EXPECT_EQ(false, t2.contains(test2));

    _String t3 = _String ("household");
    char test3 = 'o'; 
    EXPECT_EQ(true, t3.contains(test3));

}

TEST_F(StringTest,beginswithTest)
{
    //Why not have an overloaded function instead of beginsWith and startswith?
    _String test = _String ("household");
    _String t = _String ("house");
    EXPECT_EQ(true, test.beginswith(t, true));

    _String test2 = _String ("household");
    _String t2 = _String ("hold");
    EXPECT_EQ(false, test2.beginswith(t2, false));

    _String test3 = _String ("household");
    _String t3 = _String ("House");
    EXPECT_EQ(false, test3.beginswith(t3, true));

    _String test4 = _String ("household");
    _String t4 = _String ("House");
    EXPECT_EQ(true, test4.beginswith(t4, false));

    _String test5 = _String ("");
    _String t5 = _String ("");
    EXPECT_EQ(true, test5.beginswith(t5, false));

    _String test6 = _String ("AA");
    _String t6 = _String ("AABB");
    EXPECT_EQ(false, test6.beginswith(t6, false));

}

TEST_F(StringTest,startswithTest)
{
    //Why not have an overloaded function instead of beginsWith and startswith?
    _String result = _String ("household");
    _String substr = _String ("house");
    EXPECT_EQ(true, result.startswith(substr));

    result = _String ("household");
    substr = _String ("louse");
    EXPECT_EQ(false, result.startswith(substr));

    result = _String ("house");
    substr = _String ("household");
    EXPECT_EQ(false, result.startswith(substr));
}

TEST_F(StringTest,endswithTest)
{
    _String result = _String ("household");
    _String r2 = _String ("hold");
    EXPECT_EQ(true, result.endswith(r2));

    result = _String ("household");
    r2 = _String ("HOLD");
    EXPECT_EQ(true, result.endswith(r2, false));

    result = _String ("household");
    r2 = _String ("HlLD");
    EXPECT_EQ(false, result.endswith(r2, false));

    result = _String ("household");
    r2 = _String ("HOLD");
    EXPECT_EQ(false, result.endswith(r2));

    result = _String ("hold");
    r2 = _String ("household");
    EXPECT_EQ(false, result.endswith(r2));
}

TEST_F(StringTest,FormatTimeStringTest)
{
    //Takes seconds
    long time_diff = 459132;
    _String result = new _String("127:32:12");
    _String r2 = new _String("hyphy");
    r2.FormatTimeString(time_diff);
    EXPECT_STREQ(result.getStr(), r2.getStr());

    //Takes seconds
    time_diff = 0;
    result = new _String("00:00:00");
    r2 = new _String("hyphy");
    r2.FormatTimeString(time_diff);
    EXPECT_STREQ(result.getStr(), r2.getStr());
}

TEST_F(StringTest,ReplaceTest)
{

    _String orig_string = _String("household");
    _String to_replace = _String("hold");
    _String replacer = _String("house");
    _String result = _String("househouse");

    _String real_result = orig_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ(result.getStr(), real_result.getStr());

    _String empty_string = _String("");
    real_result = empty_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ("", real_result.getStr());

    _String short_string = _String("hi");
    real_result = short_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ("hi", real_result.getStr());

    short_string = _String("hi");
    to_replace = _String("");
    real_result = short_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ("hi", real_result.getStr());


    orig_string = _String("household");
    to_replace = _String("ho");
    replacer = _String("mo");
    result = _String("mousemold");

    real_result = orig_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ(result.getStr(), real_result.getStr());

    orig_string = _String("household");
    to_replace = _String("ffff");
    replacer = _String("mo");
    result = _String("household");
    real_result = orig_string.Replace(to_replace, replacer, true);

    EXPECT_STREQ(result.getStr(), real_result.getStr());

    orig_string = _String("household");
    to_replace = _String("ho");
    replacer = _String("mo");
    result = _String("mousehold");
    real_result = orig_string.Replace(to_replace, replacer, false);

    EXPECT_STREQ(result.getStr(), real_result.getStr());

    orig_string = _String("household");
    to_replace = _String("ffff");
    replacer = _String("mo");
    result = _String("household");
    real_result = orig_string.Replace(to_replace, replacer, false);

    EXPECT_STREQ(result.getStr(), real_result.getStr());
}


TEST_F(StringTest,TokenizeTest)
{
    _String test_string = _String("house,condo,hyphy");
    _String* sub_string = new _String(",");

    _List* result_list = test_string.Tokenize(sub_string);
    _String* result = (_String*)result_list->lData[0];

    EXPECT_STREQ("house", result->getStr());


//    _String test_string2 = _String("house,condo,hyphy");
//    _String* sub_string2 = new _String("");

//    _List* result_list2 = test_string2.Tokenize(sub_string2);
//    _String* result2 = (_String*)result_list2->lData[0];
//    EXPECT_STREQ("house,condo,hyphy", result2->getStr());


}

TEST_F(StringTest,toNumTest)
{
    _String test_string = _String("3.14");
    double result = test_string.toNum();
    double expected = 3.14;
    EXPECT_EQ(expected, result);

    test_string = _String("");
    result = test_string.toNum();
    expected = 0.;
    EXPECT_EQ(expected, result);
}


TEST_F(StringTest,UpCaseTest)
{
    _String result = _String("HOUSE");
    _String insert = _String("house");
    insert.UpCase();

    EXPECT_STREQ(result.getStr(), insert.getStr());
}

TEST_F(StringTest,LoCaseTest)
{
    _String insert = _String("HOUSE");
    _String result = _String("house");
    insert.LoCase();

    EXPECT_STREQ(result.getStr(), insert.getStr());
}

TEST_F(StringTest,ProcessTreeBranchLengthTest)
{
    //All this does is find the toNum, if it begins with a ':', just skip that.
    _String test_string = _String(":3.14");
    double result = test_string.ProcessTreeBranchLength();
    double expected = 3.14;
    EXPECT_EQ(expected, result);

    test_string = _String("3.14");
    result = test_string.ProcessTreeBranchLength();
    expected = 3.14;
    EXPECT_EQ(expected, result);

    test_string = _String("1e-12");
    result = test_string.ProcessTreeBranchLength();
    expected = 1e-10;
    EXPECT_EQ(expected, result);
}

TEST_F(StringTest,ExtractEnclosedExpressionTest)
{
    //TODO: Coverage and Debug
    //returns position
    _String test_string = _String("[hyp[house]hy]");
    long i = 0;
    long j = 0;
    j = test_string.ExtractEnclosedExpression(i,'[',']',true, true);
    EXPECT_EQ(13, j);
}

TEST_F(StringTest,IsALiteralArgumentTest)
{
    _String test_string = _String("\"house\"");
    bool result = test_string.IsALiteralArgument(true);
    EXPECT_EQ(true, result);

    test_string = _String("hi");
    result = test_string.IsALiteralArgument(true);
    EXPECT_EQ(false, result);
}

TEST_F(StringTest,StripQuotesTest)
{
    //Only strips the outer quotes
    _String insert = _String("\"So this\"");
    _String result = _String("So this");
    insert.StripQuotes();
    EXPECT_STREQ(result.getStr(), insert.getStr());
}

TEST_F(StringTest,IsValidIdentifierTest)
{
    // Valid Identifier must be "greater than 0, doesn't start with non-alpha character and doesn't start wtih "_" if strict, if not strict, accept number"
    // Also cannot have keyword
    _String test_string = _String("house");
    EXPECT_EQ(true, test_string.IsValidIdentifier(true));

    test_string = _String("$house");
    EXPECT_EQ(false, test_string.IsValidIdentifier(true));

    test_string = _String("");
    EXPECT_EQ(false, test_string.IsValidIdentifier(true));

    test_string = _String("$house");
    EXPECT_EQ(false, test_string.IsValidIdentifier(false));
}

TEST_F(StringTest,IsValidRefIdentifierTest)
{
    //Same as IsValidIdentifier, but ends with a &
    _String test_string = _String("house&");
    EXPECT_EQ(true, test_string.IsValidRefIdentifier());
}

/*TEST_F(StringTest,ProcessParameterTest) { */
////This will be properly tested with batchlan
//_String initial_path = _String("/home/sergei/hyphy");
//initial_path.ProcessParameter();

//EXPECT_STREQ("/home/sergei/hyphy", initial_path);
/*}*/

//TEST_F(StringTest,ProcessFileNameTest) {
////This will be properly tested with batchlan
//_String* test_string = new _String("/Users/stevenweaver/Documents/sergei/hyphy/trunk/UnitTests/mtDNA.fas");
//test_string->ProcessFileName();
//EXPECT_EQ(2, 2);
//}

TEST_F(StringTest,PathCompositionTest)
{
    _String initial_path = _String("/home/sergei/hyphy");
    _String change_path = _String("../trunk");
    _String actual_path = initial_path.PathComposition(change_path);
    _String result_path = _String("/home/sergei/trunk");
    EXPECT_STREQ(result_path.getStr(), actual_path.getStr());

    initial_path = _String("/home/sergei/hyphy");
    change_path = _String("/trunk");
    actual_path = initial_path.PathComposition(change_path);
    result_path = _String("/trunk");
    EXPECT_STREQ(result_path.getStr(), actual_path.getStr());

    initial_path = _String("/home/sergei/hyphy");
    change_path = _String("/trunk/");
    actual_path = initial_path.PathComposition(change_path);
    result_path = _String("/trunk/");
    EXPECT_STREQ(result_path.getStr(), actual_path.getStr());

}

TEST_F(StringTest,PathSubtractionTest)
{
    _String initial_path = _String("/home/sergei/");
    _String sub_path = _String("/home/sergei/hyphy");
    _String result_path = initial_path.PathSubtraction(sub_path,'A');
    EXPECT_STREQ("hyphy", result_path.getStr());

    initial_path = _String("/home/sergei/");
    sub_path = _String("/home/sergei/documents/hyphy");
    result_path = initial_path.PathSubtraction(sub_path,'A');
    EXPECT_STREQ("documents/hyphy", result_path.getStr());

    initial_path = _String("/home/steven");
    sub_path = _String("/home/sergei/documents/hyphy");
    result_path = initial_path.PathSubtraction(sub_path,'A');
    EXPECT_STREQ("sergei/documents/hyphy", result_path.getStr());
}

TEST_F(StringTest,ConvertToAnIdentTest)
{
    //Takes a String and converts it to a valid hyphy ident test
    _String initial = _String("$house");
    initial.ConvertToAnIdent();
    EXPECT_STREQ("_house", initial.getStr());

    initial.ConvertToAnIdent(false);
    EXPECT_STREQ("_house", initial.getStr());

    initial = _String("_house");
    initial.ConvertToAnIdent();
    EXPECT_STREQ("_house", initial.getStr());

    initial.ConvertToAnIdent(false);
    EXPECT_STREQ("_house", initial.getStr());

    _String initial2 = _String("_ho_use");
    initial2.ConvertToAnIdent();
    EXPECT_STREQ("_ho_use", initial2.getStr());
}

TEST_F(StringTest,ShortenVarIDTest)
{
    //TODO Debug Coverage
    _String initial = _String("house.room");
    _String container = _String("house");
    _String result = initial.ShortenVarID(container);
    EXPECT_STREQ("room", result.getStr());

    initial = _String("house.room");
    container = _String("friend");
    result = initial.ShortenVarID(container);
    EXPECT_STREQ("house.room", result.getStr());
}

TEST_F(StringTest,RegExpMatchOnceTest)
{
    _String initial = new _String("hyphy");
    _String* pattern = new _String("hyph");
    _SimpleList matched_pairs;
    initial.RegExpMatchOnce(pattern, matched_pairs, false, false);
    EXPECT_EQ(0,matched_pairs.lData[0]);
}

TEST_F(StringTest,RegExpMatchTest)
{

    _String initial = new _String("hyphy");
    _String* pattern = new _String("phy");
    _SimpleList matched_pairs;
    initial.RegExpMatchOnce(pattern, matched_pairs, false, false);

    EXPECT_EQ(2,matched_pairs.lData[0]);
}

TEST_F(StringTest,RegExpMatchAllTest)
{

    _String initial = new _String("hyphy");
    _String* pattern = new _String("phy");
    _SimpleList matched_pairs;
    int errNo = 0;

    Ptr regex = PrepRegExp (pattern, errNo, false);
    initial.RegExpMatchAll(regex, matched_pairs);

    EXPECT_EQ(2,matched_pairs.lData[0]);
}

TEST_F(StringTest,LempelZivProductionHistoryTest)
{
    //{1,0,01,1110,1100, 0010}
    _String initial = _String("1001111011000010");
    _SimpleList rec;
    long result = initial.LempelZivProductionHistory(&rec);

    EXPECT_EQ(6, result);
}

TEST_F(StringTest,SortTest)
{
    _String initial = _String("hgfedcba");
    _SimpleList* index = new _SimpleList();
    _String* result = initial.Sort(index);
    EXPECT_STREQ("abcdefgh", result->getStr());


    initial = _String();
    _SimpleList index2;
    _String* result2 = initial.Sort(&index2);
    EXPECT_STREQ(nil, result2->getStr());

    delete index;
    delete result;
    delete result2;
}

TEST_F(StringTest,FindTerminatorTest)
{
    _String initial = _String("dither");
    _String terminator = _String("h");
    long result = initial.FindTerminator(0,terminator);
    EXPECT_EQ(3, result);

    initial = _String("\"[dither]{dither}(dither)\"dither");
    terminator = _String("h");
    result = initial.FindTerminator(0,terminator);
    EXPECT_EQ(29, result);
}

TEST_F(StringTest,AppendAnAssignmentToBufferTest)
{

    _String initial = _String("dither");
    _String append = _String("test");
    _String append2 = _String("34");
    initial.AppendAnAssignmentToBuffer(&append, &append2, false, false, false);

    _String expected = _String("dithertest=34;\n");
    EXPECT_EQ(expected.getStr()[expected.sLength - 1], 
              initial.getStr()[expected.sLength -1]);

    initial = _String("dither");
    append = _String("12");
    _String* pAppend = new _String("34");
    initial.AppendAnAssignmentToBuffer(&append, pAppend, true, true, true);

    _String expected2 = _String("dither12:=\"34\";\n");
    EXPECT_EQ(expected2.getStr()[expected2.sLength - 1], 
              initial.getStr()[expected2.sLength - 1]);

}

/*
 *TEST_F(StringTest,AppendVariableValueAVLTest) {
 *
 *    _String initial = _String("life");
 *    _String* append = new _String("answer");
 *
 *    unsigned long num = 5;
 *    _SimpleList matched_pairs(num);
 *
 *    long int1 = 4;
 *    long int2 = 2;
 *
 *    matched_pairs.InsertElement((BaseRef)int1,-1,false,false);
 *    matched_pairs.InsertElement((BaseRef)int1,-1,false,false);
 *
 *    initial.AppendVariableValueAVL(append, matched_pairs);
 *
 *    EXPECT_STREQ("life[answer]=4,2", initial.getStr());
 *}
 */


//Operator Tests
TEST_F(StringTest,BracketTest)
{
    //[]
    _String result = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ('e', result[5L]);
}

TEST_F(StringTest,ParanthTest)
{
    //()
    _String result = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ('e', result(5));

    EXPECT_EQ(0, result(500));
}

TEST_F(StringTest,EqualTest)
{
    //=
    _String result = _String ("hyphy");
    _String dupe = result;
    EXPECT_STREQ(result.getStr(), dupe.getStr());
}

TEST_F(StringTest, AmpersandTest)
{

    // &
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String result = orig&to_append;

    EXPECT_STREQ(expected.getStr(), result.getStr());
}

TEST_F(StringTest,DoubleLessTest)
{
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String result = orig&to_append;

    EXPECT_STREQ(expected.getStr(), result.getStr());
}

TEST_F(StringTest,DoubleEqualTest)
{
    _String* result = new _String ("AABBCCDD");
    _String* r2 = new _String ("AABBCCDD");
    EXPECT_EQ(true, result->Equal(r2));

    delete result;
    delete r2;
}

TEST_F(StringTest,GreaterOpTest)
{
    // >
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(false, result>r2);
}

TEST_F(StringTest,LesserOpTest)
{
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(true, result<r2);
}

TEST_F(StringTest,LesserEqualOpTest)
{
    //<=
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(true, result<=r2);
}

TEST_F(StringTest,GreaterEqualOpTest)
{
    // >=
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(false, result>=r2);
}

TEST_F(StringTest,NotEqualOpTest)
{
    // !=
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(true, result!=r2);
}

//Diagnostic Testing

//I should probably move this out into its own separate test
/*
 *class _KMPDiagnosticTest : public TestWithParam< ::std::tr1::tuple<int, int> > {
 *    protected:
 *    // You can remove any or all of the following functions if its body
 *    // is empty.
 *
 *    _KMPDiagnosticTest() {
 *        // You can do set-up work for each test here.
 *    }
 *
 *    virtual ~_KMPDiagnosticTest() {
 *        // You can do clean-up work that doesn't throw exceptions here.
 *    }
 *
 *    // If the constructor and destructor are not enough for setting up
 *    // and cleaning up each test, you can define the following methods:
 *
 *    virtual void SetUp() {
 *        // Code here will be called immediately after the constructor (right
 *        // before each test).
 *        FILE * test;
 *        test = fopen ("/Users/stevenweaver/Documents/sergei/hyphy/trunk/UnitTests/mtDNA.fas","r");
 *        buffer = new _String(test);
 *        fclose(test);
 *
 *        substr = *buffer;
 *
 *        rand_start = pow(2,::std::tr1::get<0>(GetParam()));
 *        rand_length = pow(10,::std::tr1::get<1>(GetParam()));
 *
 *        substr.Trim(rand_start, rand_start+rand_length);
 *        buffer->buildKmpTable(substr);
 *    }
 *
 *    virtual void TearDown() {
 *        // Code here will be called immediately after each test (right
 *        // before the destructor).
 *        delete buffer;
 *    }
 *
 *    // Objects declared here can be used by all tests in the test case for Foo.
 *    _String* buffer;
 *    _String substr;
 *    long rand_start;
 *    long rand_length;
 *};
 *
 *INSTANTIATE_TEST_CASE_P(SpeedTest,
 *                        _KMPDiagnosticTest,
 *                        Combine(Range(1,10), Range(1,10)));
 *
 *TEST_P(_KMPDiagnosticTest, KMPDiagnostic) {
 *    // Inside a test, access the test parameter with the GetParam() method
 *    // of the TestWithParam<T> class:
 *    //substr
 *    int i=0;
 *    while(i<25){
 *        buffer->FindKMP(substr);
 *        ++i;
 *    }
 *
 *    EXPECT_EQ(rand_start, buffer->FindKMP(substr));
 *}
 *
 *TEST_P(_KMPDiagnosticTest, FindDiagnostic) {
 *    // Inside a test, access the test parameter with the GetParam() method
 *    // of the TestWithParam<T> class:
 *    int i=0;
 *    while(i<25){
 *        buffer->Find(substr);
 *        ++i;
 *    }
 *
 *    //substr
 *    EXPECT_EQ(rand_start, buffer->Find(substr));
 *}
 */

}  // namespace
