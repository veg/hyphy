/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
 Steven Weaver (sweaver@temple.edu)
 
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

    EXPECT_EQ (test.length(), 0);
    EXPECT_EQ (empty.length(), 0);
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
    EXPECT_STREQ(test.get_str(), dupe1.get_str());

    dupe2.Duplicate(&empty);
    EXPECT_EQ (dupe2.length(), 0);
}

/******************************************/


TEST_F(StringTest,makeDynamicTest)
{

    /* 20110825: SLKP code coverage complete
          The idea of make dynamic is to convert a stack object into a heap object, see code below
     */

    _String stackString (globalTest1, 5, -1), // this helps test one of the constructors
            *heapString = (_String*)stackString.makeDynamic(),
            empty;

    empty.Initialize();

    EXPECT_STREQ (stackString.get_str(), heapString->get_str());
    stackString = empty; // overwrite the stack object
    EXPECT_EQ (heapString->length(), globalTest1.length()-5);
    DeleteObject (heapString);

}

/******************************************/

TEST_F(StringTest,getCharTest)
{

    /* 20110825: SLKP code coverage complete */

    _String test (globalTest1),
            empty;

    empty.Initialize();

    EXPECT_EQ('e', test.get_char(5));

    //Default return is 0
    _String test2 = empty;
    EXPECT_EQ(0, test2.get_char(5));

}

/******************************************/


TEST_F(StringTest,setCharTest)
{
    _String test (globalTest1),
            empty;

    empty.Initialize();

    test.set_char(5,'d');
    EXPECT_EQ('d', test.get_char(5));

    //Should be 0
    _String test2 = empty;
    test2.set_char(5,'d');
    EXPECT_EQ('\0', test2.get_char(5));
}

/******************************************/

TEST_F(StringTest, lengthTest)
{
    _String test = new _String ("You're asking me to run MCMC without reporting any tests.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ(99, test.length());

    _String test2 = new _String ("");
    EXPECT_EQ(0, test2.length());
}

TEST_F(StringTest,InsertTest)
{
    _String test = "AAGGCCTTA";
    _String expected_test = "CAAGGCCTTA";
    test.Insert('C',0);
    EXPECT_STREQ(expected_test.get_str(), test.get_str());

    //_String test2 = "";
    //_String expected_test2 = "";
    //test2.Insert('C',5);
    //EXPECT_STREQ(expected_test2.get_str(), test2.get_str());

    test = "AAGGCCTTA";
    expected_test = "AAGGCCTTAC";
    test.Insert('C',-11);
    EXPECT_STREQ(expected_test.get_str(), test.get_str());

}

TEST_F(StringTest,DeleteTest)
{
    _String test = "AAGGCCTTA";
    test.Delete(3,4);
    EXPECT_STREQ("AAGCTTA", test.get_str());

    /*  //ERROR: Program crashes if you go beyond end of string.
     *
     *  _String test2 = "AAGG";
     *  test2.Delete(1,14);
     *  EXPECT_STREQ("A", test2.get_str());
     *
     */

    _String test2 = "AAGG";
    test2.Delete(1,-5);
    EXPECT_STREQ("A", test2.get_str());

    _String test3 = "AAGG";
    test3.Delete(-1,-5);
    EXPECT_STREQ("", test3.get_str());

}

TEST_F(StringTest,AppendNewInstanceTest)
{
    //TODO: Debug for coverage
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String test = orig&to_append;
    EXPECT_STREQ(expected.get_str(), test.get_str());

    orig = _String ("");
    to_append = _String("");
    expected = _String("");
    test = orig&to_append;
    EXPECT_STREQ("", expected.get_str());

}


//TEST_F(StringTest,EscapeAndAppendTest)
//{
//    _String result = _String("AAGG");
//    _String expected = _String("AAGG\\\\(&lt;\\\[''");


//    result.EscapeAndAppend('\\',2);
//    result.EscapeAndAppend('(',1);
//    result.EscapeAndAppend('<',4);
//    result.EscapeAndAppend('[',5);
//    result.EscapeAndAppend('\'',2);

//    EXPECT_EQ(expected.get_str()[expected.length() - 1], 
//              result.get_str()[expected.length() - 1]);

//}

//TEST_F(StringTest,EscapeAndAppend2Test)
//{
//    _String str = _String("AAGG");
//    _String append = _String("\"AAGG\"");
//    _String expected = _String("AAGG&quot;AAGG&quot;");

//    str.EscapeAndAppend(append,4);

//    EXPECT_EQ(expected.get_str()[expected.length() - 1], 
//              str.get_str()[expected.length() - 1]);
//}

TEST_F(StringTest,get_strTest)
{
    _String test = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_STREQ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n", test.get_str());

    _String test2 = _String ("");
    EXPECT_STREQ("", test2.get_str());

}

TEST_F(StringTest,ChopTest)
{

    _String test = globalTest1;
    _String expected = _String("'re asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r1 = test.Chop(0,2);
    EXPECT_STREQ(expected.get_str(), r1.get_str());

    _String test2 = globalTest1;
    _String expected2 = _String("'re asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r2 = test2.Chop(-1,2);
    EXPECT_STREQ(expected2.get_str(), r2.get_str());

    _String test3 = _String ("");
    _String r3 = test3.Chop(0,2);
    EXPECT_STREQ("", r3.get_str());

    // TODO: No longer works
    //_String test4 = globalTest1;
    //_String r4 = test4.Chop(5,4);
    //EXPECT_STREQ("", r4.get_str());

    _String test5 = globalTest1;
    _String expected5 = _String("Y're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r5 = test5.Chop(1,2);
    EXPECT_STREQ(expected5.get_str(), r5.get_str());

    _String test6 = globalTest1;
    _String r6 = test6.Chop(0,-1);
    EXPECT_STREQ("", r6.get_str());

    //Error: Should return empty, but it returns "AA"
    //_String test7 = _String ("ABBA");
    //_String substr7 = test7.Chop(11,2);
    //EXPECT_STREQ("", substr7.get_str());
 

    //Error: Memory allocation error
    //_String test3 = _String ("ABBA");
    //_String substr3 = test3.Chop(2,20);
    //EXPECT_STREQ("BA", substr3.get_str());
   
}

TEST_F(StringTest,CutTest)
{
    _String test = globalTest1;
    _String substr = test.Cut(0,2);
    EXPECT_STREQ("You", substr.get_str());

    _String test2 = _String ("AABBCC");
    _String substr2 = test2.Cut(4,12);
    EXPECT_STREQ("CC", substr2.get_str());

    _String test3 = _String ("");
    _String substr3 = test3.Cut(4,13);
    EXPECT_STREQ("", substr3.get_str());

    _String test4 = _String ("AABBCC");
    _String substr4 = test4.Cut(-1,-1);
    EXPECT_STREQ("AABBCC", substr4.get_str());

    _String test5 = _String ("AABBCC");
    _String substr5 = test5.Cut(5,4);
    EXPECT_STREQ("", substr5.get_str());

    _String test6 = _String ("");
    _String substr6 = test6.Cut(4,13);
    EXPECT_STREQ("", substr6.get_str());

}

TEST_F(StringTest,FlipTest)
{
    _String result = _String ("ABC");
    result.Flip();
    EXPECT_STREQ("CBA",result.get_str());
}

TEST_F(StringTest,Adler32Test)
{
    //TODO: Coverage
    _String result = new _String ("Wikipedia");
    EXPECT_EQ(300286872, result.Adler32());
}

TEST_F(StringTest,FirstNonSpaceIndexTest)
{
    //Error: curious that it would give me an 8
    _String test = _String ("    hyphy");
    EXPECT_EQ(4, test.FirstNonSpaceIndex());

    test = _String ("    hyphy");
    EXPECT_EQ(4, test.FirstNonSpaceIndex(-1));

    _String test2 = _String ("    hyphy ");
    EXPECT_EQ(4, test2.FirstNonSpaceIndex(0,-1,kStringDirectionForward));

    _String test3 = _String ("");
    EXPECT_EQ(-1, test3.FirstNonSpaceIndex());

}

TEST_F(StringTest,KillSpacesTest)
{
    _String result = _String ("  h   y   p    ");
    EXPECT_STREQ("hyp", result.KillSpaces());

    _String test2 = _String ("hyp");
    EXPECT_STREQ("hyp", test2.KillSpaces());

    _String test3 = _String ("");
    EXPECT_STREQ("", test3.KillSpaces());
}

TEST_F(StringTest,CompressSpacesTest)
{
    _String test = _String ("ABB   and    CCD");
    EXPECT_STREQ("ABB and CCD",test.CompressSpaces());

    _String test2 = _String ("AB C");
    EXPECT_STREQ("AB C",test2.CompressSpaces());

    _String test3 = _String ("");
    EXPECT_STREQ("",test3.CompressSpaces());

}

TEST_F(StringTest,FirstSpaceIndexTest)
{
    _String test = _String ("AA BB");
    EXPECT_EQ(2, test.FirstSpaceIndex());

    _String test2 = _String ("");
    EXPECT_EQ(-1, test2.FirstSpaceIndex());

    _String test3 = _String ("A BBB");
    EXPECT_EQ(1, test3.FirstSpaceIndex(-1,-1,kStringDirectionForward));

    _String test4 = _String (" A BBB");
    EXPECT_EQ(0, test4.FirstSpaceIndex());
}

TEST_F(StringTest,FirstNonSpaceTest)
{
    _String test = _String ("  AA BB");
    EXPECT_EQ('A', test.FirstNonSpace());

    _String test2 = _String ("AABB ");
    EXPECT_EQ('B', test2.FirstNonSpace(0,-1,kStringDirectionBackward));
}

TEST_F(StringTest,FirstNonSpace2Test)
{
    _String test = _String ("  AA BB");
    EXPECT_EQ('A', test.FirstNonSpace());

    _String test2 = _String ("AABB ");
    EXPECT_EQ('B', test2.FirstNonSpace(0,-1,kStringDirectionBackward));
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

TEST_F(StringTest,EqualsTest)
{
    _String* test = new _String ("AABBCCDD");
    _String* r2 = new _String ("AABBCCDD");
    EXPECT_EQ(true, test->Equal(r2));

    //_String test2 = _String("AADCC");
    //EXPECT_EQ(false, test2.Equal(r2));

    //_String test3("AABBCCDD");
    //_String r3("AABBCCDD");
    //EXPECT_TRUE(test3==r3);

}

TEST_F(StringTest,CompareTest)
{
    //house precedes household
    //Household precedes house
    //composer precedes computer
    //H2O precedes HOTEL

    _String result = _String ("household");
    _String* substr = new _String ("house");
    EXPECT_EQ(1, result.Compare(substr));

    _String result2 = _String ("household");
    _String* substr2 = new _String ("household");
    EXPECT_EQ(0, result2.Compare(substr2));
}

TEST_F(StringTest,EqualWithWildCharTest)
{

    _String test = _String ("AABBCCDD");
    _String* t = new _String ("EEBBCCDD");

    EXPECT_EQ(true, test.EqualWithWildChar(t, 'E'));

    //_String test2 = _String ("AAFBCCDD");
    //EXPECT_EQ(false, test2.EqualWithWildChar(t, 'E'));

    //_String test3 = _String ("EEBBBCDD");
    //EXPECT_EQ(false, test3.EqualWithWildChar(t, 'B'));

    //_String test4 = _String ("EEBBBCDB");
    //EXPECT_EQ(false, test4.EqualWithWildChar(t, 'B'));

    //delete t;
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

TEST_F(StringTest,beginswithTest)
{
    //Why not have an overloaded function instead of beginsWith and BeginsWith?
    _String test = _String ("household");
    _String t = _String ("house");
    EXPECT_EQ(true, test.BeginsWith(t, true));

    _String test2 = _String ("household");
    _String t2 = _String ("hold");
    EXPECT_EQ(false, test2.BeginsWith(t2, false));

    _String test3 = _String ("household");
    _String t3 = _String ("House");
    EXPECT_EQ(false, test3.BeginsWith(t3, true));

    _String test4 = _String ("household");
    _String t4 = _String ("House");
    EXPECT_EQ(true, test4.BeginsWith(t4, false));

    _String test5 = _String ("");
    _String t5 = _String ("");
    EXPECT_EQ(true, test5.BeginsWith(t5, false));

    _String test6 = _String ("AA");
    _String t6 = _String ("AABB");
    EXPECT_EQ(false, test6.BeginsWith(t6, false));

}

TEST_F(StringTest,endswithTest)
{
    _String result = _String ("household");
    _String r2 = _String ("hold");
    EXPECT_EQ(true, result.EndsWith(r2));

    result = _String ("household");
    r2 = _String ("HOLD");
    EXPECT_EQ(true, result.EndsWith(r2, false));

    result = _String ("household");
    r2 = _String ("HlLD");
    EXPECT_EQ(false, result.EndsWith(r2, false));

    result = _String ("household");
    r2 = _String ("HOLD");
    EXPECT_EQ(false, result.EndsWith(r2));

    result = _String ("hold");
    r2 = _String ("household");
    EXPECT_EQ(false, result.EndsWith(r2));
}

TEST_F(StringTest,FormatTimeStringTest)
{
    //Takes seconds
    long time_diff = 459132;
    _String result = new _String("127:32:12");
    _String r2 = new _String("hyphy");
    EXPECT_STREQ(result.get_str(), r2.FormatTimeString(time_diff));

    //Takes seconds
    time_diff = 0;
    result = new _String("00:00:00");
    r2 = new _String("hyphy");
    EXPECT_STREQ(result.get_str(), r2.FormatTimeString(time_diff));
}

TEST_F(StringTest,ReplaceTest)
{

    _String orig_string = _String("household");
    _String to_replace = _String("hold");
    _String replacer = _String("house");
    _String result = _String("househouse");

    _String real_result = orig_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ(result.get_str(), real_result.get_str());

    _String empty_string = _String("");
    real_result = empty_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ("", real_result.get_str());

    _String short_string = _String("hi");
    real_result = short_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ("hi", real_result.get_str());

    short_string = _String("hi");
    to_replace = _String("");
    real_result = short_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ("hi", real_result.get_str());


    orig_string = _String("household");
    to_replace = _String("ho");
    replacer = _String("mo");
    result = _String("mousemold");

    real_result = orig_string.Replace(to_replace, replacer, true);
    EXPECT_STREQ(result.get_str(), real_result.get_str());

    orig_string = _String("household");
    to_replace = _String("ffff");
    replacer = _String("mo");
    result = _String("household");
    real_result = orig_string.Replace(to_replace, replacer, true);

    EXPECT_STREQ(result.get_str(), real_result.get_str());

    orig_string = _String("household");
    to_replace = _String("ho");
    replacer = _String("mo");
    result = _String("mousehold");
    real_result = orig_string.Replace(to_replace, replacer, false);

    EXPECT_STREQ(result.get_str(), real_result.get_str());

    orig_string = _String("household");
    to_replace = _String("ffff");
    replacer = _String("mo");
    result = _String("household");
    real_result = orig_string.Replace(to_replace, replacer, false);

    EXPECT_STREQ(result.get_str(), real_result.get_str());
}


TEST_F(StringTest,TokenizeTest)
{
    _String test_string = _String("house,condo,hyphy");
    _String* sub_string = new _String(",");

    _List result_list = test_string.Tokenize(sub_string);
    _String* result = (_String*)result_list.list_data[0];

    EXPECT_STREQ("house", result->get_str());


//    _String test_string2 = _String("house,condo,hyphy");
//    _String* sub_string2 = new _String("");

//    _List* result_list2 = test_string2.Tokenize(sub_string2);
//    _String* result2 = (_String*)result_list2->list_data[0];
//    EXPECT_STREQ("house,condo,hyphy", result2->get_str());


}

TEST_F(StringTest,toNumTest)
{
    _String test_string = _String("3.14");
    double result = test_string.to_float();
    double expected = 3.14;
    EXPECT_EQ(expected, result);

    test_string = _String("");
    result = test_string.to_float();
    expected = 0.;
    EXPECT_EQ(expected, result);
}


TEST_F(StringTest,UpCaseTest)
{
    _String result = _String("HOUSE");
    _String insert = _String("house");
    EXPECT_STREQ(result.get_str(), insert.ChangeCase(kStringUpperCase));
}

TEST_F(StringTest,LoCaseTest)
{
    _String insert = _String("HOUSE");
    _String result = _String("house");
    EXPECT_STREQ(result.get_str(), insert.ChangeCase(kStringLowerCase));
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
    EXPECT_STREQ(result.get_str(), insert.get_str());
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


TEST_F(StringTest,ConvertToAnIdentTest)
{
    //Takes a String and converts it to a valid hyphy ident test
    _String initial = _String("$house");

    EXPECT_STREQ("_house", initial.ConvertToAnIdent());
    EXPECT_STREQ("_house", initial.ConvertToAnIdent(false));

    initial = _String("_house");
    EXPECT_STREQ("_house", initial.ConvertToAnIdent());
    EXPECT_STREQ("_house", initial.ConvertToAnIdent(false));

    _String initial2 = _String("_ho_use");
    EXPECT_STREQ("_ho_use", initial2.ConvertToAnIdent());
}

TEST_F(StringTest,LempelZivProductionHistoryTest)
{
    //{1,0,01,1110,1100, 0010}
    _String initial = _String("1001111011000010");
    _SimpleList rec;
    long result = initial.LempelZivProductionHistory(&rec);

    EXPECT_EQ(6, result);
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
    EXPECT_STREQ(result.get_str(), dupe.get_str());
}

TEST_F(StringTest, AmpersandTest)
{

    // &
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String result = orig&to_append;

    EXPECT_STREQ(expected.get_str(), result.get_str());
}

TEST_F(StringTest,DoubleLessTest)
{
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String result = orig&to_append;

    EXPECT_STREQ(expected.get_str(), result.get_str());
}

TEST_F(StringTest,DoubleEqualTest)
{
    _String* result = new _String ("AABBCCDD");
    _String* r2 = new _String ("AABBCCDD");
    EXPECT_EQ(true, result->Equal(r2));

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

}  // namespace
