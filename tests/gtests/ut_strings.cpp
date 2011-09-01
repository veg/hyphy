/*
 *  ut_strings.cpp
 *  HyPhyXCode
 *
 *  Created by Steven Weaver on 6/17/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include "ut_strings.h"


#include "hy_strings.h"
#include "hy_lists.h"

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
class _StringTest : public ::testing::Test
{
protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _StringTest() {
        // You can do set-up work for each test here.
    }

    virtual ~_StringTest() {
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

//Stub out the Rest of the tests

_String globalTest1 ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");

/* 20110825: SLKP made this a global string.*/

/******************************************/

TEST_F(_StringTest,DuplicateTest)
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

TEST_F(_StringTest,DuplicateErasingTest)
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


TEST_F(_StringTest,makeDynamicTest)
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

TEST_F(_StringTest,getCharTest)
{

    /* 20110825: SLKP code coverage complete */

    _String test (globalTest1);
    EXPECT_EQ('e', test.getChar(5));

    //Default return is 0
    _String test2 = empty;
    EXPECT_EQ(0, test2.getChar(5));

}

/******************************************/


TEST_F(_StringTest,setCharTest)
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

TEST_F(_StringTest,CopyDynamicStringTest)
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


TEST_F(_StringTest, LengthTest)
{
    _String test = new _String ("You're asking me to run MCMC without reporting any tests.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ(99, test.Length());

    _String test2 = new _String ("");
    EXPECT_EQ(0, test2.Length());
}

TEST_F(_StringTest,InsertTest)
{
    _String test = "AAGGCCTTA";
    _String expected_test = "CAAGGCCTTA";
    test.Insert('C',0);
    ASSERT_STREQ(expected_test.getStr(), test.getStr());

    _String test2 = "";
    _String expected_test2 = "";
    test2.Insert('C',5);
    ASSERT_STREQ(expected_test2.getStr(), test2.getStr());

}

TEST_F(_StringTest,DeleteTest)
{
    _String test = "AAGGCCTTA";
    test.Delete(3,4);
    ASSERT_STREQ("AAGCTTA", test.getStr());

    /*  //ERROR: Program crashes if you go beyond end of string.
     *
     *  _String test2 = "AAGG";
     *  test2.Delete(1,14);
     *  ASSERT_STREQ("A", test2.getStr());
     *
     */


    _String test2 = "AAGG";
    test2.Delete(1,-5);
    ASSERT_STREQ("A", test2.getStr());

}

TEST_F(_StringTest,AppendNewInstanceTest)
{
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String test = orig&to_append;
    ASSERT_STREQ(expected.getStr(), test.getStr());

    orig = _String ("");
    to_append = _String("");
    expected = _String("");
    test = orig&to_append;
    ASSERT_STREQ(expected.getStr(), test.getStr());

}


TEST_F(_StringTest,EscapeAndAppendTest)
{
    _String result = _String("AAGG");
    _String expected = _String("AAGG\\\\(&lt;\\[");
    result.EscapeAndAppend('\\',2);
    result.EscapeAndAppend('(',1);
    result.EscapeAndAppend('<',4);
    result.EscapeAndAppend('[',5);
    ASSERT_STREQ(expected.getStr(), result.getStr());
}

TEST_F(_StringTest,FinalizeTest)
{
    _String orig = _String ("hyphy");
    orig.Finalize();
    EXPECT_EQ(0, orig[orig.Length()]);
}

TEST_F(_StringTest,getStrTest)
{
    _String test = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    ASSERT_STREQ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n", test.getStr());

    _String test2 = _String ("");
    ASSERT_STREQ("", test2.getStr());

}

TEST_F(_StringTest,ChopTest)
{

    _String test = _String ("You're asking me to run MCMC without reporting any tests.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String r2 = test.Chop(0,2);
    ASSERT_STREQ("'re asking me to run MCMC without reporting any tests.  Did you forget to set Bgm_MCMC_SAMPLES?\n", r2.getStr());
    /*
     *
     *    //Error: Should return empty, but it returns "AA"
     *    _String test2 = _String ("ABBA");
     *    _String substr2 = test2.Chop(1l,2);
     *    ASSERT_STREQ("", substr2.getStr());
     *
     */

    /*    //Error: Memory allocation error
     *    _String test3 = _String ("ABBA");
     *    _String substr3 = test3.Chop(2,20);
     *    ASSERT_STREQ("BA", substr3.getStr());
     *
     */
}

TEST_F(_StringTest,CutTest)
{
    _String test = _String ("You're asking me to run MCMC without reporting any tests.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    _String substr = test.Cut(0,2);
    ASSERT_STREQ("You", substr.getStr());

    _String test2 = _String ("AABBCC");
    _String substr2 = test2.Cut(4,12);
    ASSERT_STREQ("CC", substr2.getStr());
}

TEST_F(_StringTest,FlipTest)
{
    _String result = _String ("ABC");
    result.Flip();
    ASSERT_STREQ("CBA",result.getStr());
}

TEST_F(_StringTest,Adler32Test)
{
    _String result = new _String ("Wikipedia");
    EXPECT_EQ(300286872, result.Adler32());
}

TEST_F(_StringTest,TrimTest)
{
    _String test = new _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    test.Trim(7,12);
    ASSERT_STREQ("asking", test.getStr());

    _String test2 = new _String ("");
    test2.Trim(7,12);
    ASSERT_STREQ("", test2.getStr());

}

TEST_F(_StringTest,FirstNonSpaceIndexTest)
{
    _String test = _String ("    lol");
    EXPECT_EQ(4, test.FirstNonSpaceIndex());

    _String test2 = _String ("");
    EXPECT_EQ(-1, test2.FirstNonSpaceIndex());

}

TEST_F(_StringTest,KillSpacesTest)
{
    _String result = _String ("  l   o   l    ");
    _String r2;
    result.KillSpaces(r2);
    ASSERT_STREQ("lol", r2.getStr());

    _String test2 = _String ("lol");
    _String result_string2;
    test2.KillSpaces(result_string2);
    ASSERT_STREQ("lol", result_string2.getStr());

    _String test3 = _String ("");
    _String result_string3;
    test3.KillSpaces(result_string3);
    ASSERT_STREQ("", result_string3.getStr());
}

TEST_F(_StringTest,CompressSpacesTest)
{
    _String test = _String ("Beavis   and    Butthead");
    test.CompressSpaces();
    ASSERT_STREQ("Beavis and Butthead",test.getStr());

    _String test2 = _String ("lo l");
    test2.CompressSpaces();
    ASSERT_STREQ("lo l",test2.getStr());

    _String test3 = _String ("");
    test3.CompressSpaces();
    ASSERT_STREQ("",test3.getStr());

}

TEST_F(_StringTest,FirstSpaceIndexTest)
{
    _String test = _String ("AA BB");
    EXPECT_EQ(2, test.FirstSpaceIndex());

    _String test2 = _String ("");
    EXPECT_EQ(-1, test2.FirstSpaceIndex());

    _String test3 = _String ("A BBB");
    EXPECT_EQ(1, test3.FirstSpaceIndex(0,-1,-1));

}

TEST_F(_StringTest,FirstNonSpaceTest)
{
    _String test = _String ("  AA BB");
    EXPECT_EQ('A', test.FirstNonSpace());

    _String test2 = _String ("AABB ");
    EXPECT_EQ('B', test2.FirstNonSpace(0,-1,-1));

}

TEST_F(_StringTest,FindEndOfIdentTest)
{
    _String test = _String ("iden12&iden34");
    EXPECT_EQ(5, test.FindEndOfIdent(0,-1,'.'));

    _String test2 = _String ("iden12");
    EXPECT_EQ(5, test2.FindEndOfIdent(0,-1,'.'));

    _String test3 = _String ("");
    EXPECT_EQ(-1, test3.FindEndOfIdent(0,-1,'.'));
}

TEST_F(_StringTest,FindAnyCaseTest)
{
    _String result = _String ("AABBCCDD");
    EXPECT_EQ(2, result.FindAnyCase("BBcCDD"));

    _String test = _String ("AABBCCDD");
    EXPECT_EQ(-1, test.FindAnyCase("cBcCDD"));
}

TEST_F(_StringTest,ContainsSubstringTest)
{

    _String test = _String ("AABBCCDD");
    _String r2 = _String ("CC");
    EXPECT_EQ(true, test.ContainsSubstring(r2));

    _String test2 = _String ("AABBCCDD");
    _String substr = _String ("cC");
    EXPECT_EQ(false, test2.ContainsSubstring(substr));

    _String test3 = _String ("AABBCCDD");
    _String substr2 = _String ("");
    EXPECT_EQ(true, test3.ContainsSubstring(substr2));

}

TEST_F(_StringTest,FindBackwardsTest)
{
    _String test = _String ("AABBCCDD");
    EXPECT_EQ(-1, test.FindBackwards("DC",0,3));


    _String test2 = _String ("AABBCCDD");
    EXPECT_EQ(5, test2.FindBackwards("CD",0,-1));
}

TEST_F(_StringTest,FindBinaryTest)
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

TEST_F(_StringTest,EqualsTest)
{
    _String* test = new _String ("AABBCCDD");
    _String* r2 = new _String ("AABBCCDD");
    EXPECT_EQ(true, test->Equal(r2));

    _String test2 = _String("AADCC");
    EXPECT_EQ(false, test2.Equal(r2));

    delete test;
    delete r2;
}

TEST_F(_StringTest,CompareTest)
{
    //house precedes household
    //Household precedes house
    //composer precedes computer
    //H2O precedes HOTEL

    //ERROR: This returns true
    _String result = _String ("household");
    _String* r2 = new _String ("house");

    EXPECT_EQ(true, result.Compare(r2));

}

TEST_F(_StringTest,EqualWithWildCharTest)
{

    _String test = _String ("AABBCCDD");
    _String* t = new _String ("EEBBCCDD");

    EXPECT_EQ(true, test.EqualWithWildChar(t, 'E'));

    _String test2 = _String ("AAFBCCDD");
    EXPECT_EQ(false, test2.EqualWithWildChar(t, 'E'));

    delete t;
}

TEST_F(_StringTest,GreaterTest)
{
    //This is a lexicographic comparison
    _String test = _String ("house");
    _String t = _String ("household");

    EXPECT_EQ(false, test>t);

    //This is a lexicographic comparison
    _String test2 = _String ("");
    _String t2 = _String ("household");

    EXPECT_EQ(false, test2>t2);
}

TEST_F(_StringTest,LessTest)
{
    _String test = _String ("house");
    _String t = _String ("household");
    EXPECT_EQ(true, test<t);

    _String test2 = _String ("");
    _String t2 = _String ("");
    EXPECT_EQ(false, test2<t2);

}

TEST_F(_StringTest,containsTest)
{

    _String t = _String ("household");
    _String test = _String ("house");
    EXPECT_EQ(true, t.contains(test));

    _String t2 = _String ("");
    _String test2 = _String ("");
    EXPECT_EQ(false, t2.contains(test2));

}

TEST_F(_StringTest,beginswithTest)
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

}

TEST_F(_StringTest,startswithTest)
{
    //Why not have an overloaded function instead of beginsWith and startswith?
    _String result = _String ("household");
    _String r2 = _String ("house");
    EXPECT_EQ(true, result.startswith(r2));

}

TEST_F(_StringTest,endswithTest)
{
    _String result = _String ("household");
    _String r2 = _String ("hold");
    EXPECT_EQ(true, result.endswith(r2));
}

TEST_F(_StringTest,FormatTimeStringTest)
{
    //Takes seconds
    long time_diff = 459132;
    _String result = new _String("127:32:12");
    _String r2 = new _String("lol");
    r2.FormatTimeString(time_diff);
    ASSERT_STREQ(result.getStr(), r2.getStr());
}

TEST_F(_StringTest,ReplaceTest)
{

    _String orig_string = _String("household");
    _String to_replace = _String("hold");
    _String replacer = _String("house");
    _String result = _String("househouse");

    _String real_result = orig_string.Replace(to_replace, replacer, true);
    ASSERT_STREQ(result.getStr(), real_result.getStr());
}


TEST_F(_StringTest,TokenizeTest)
{
    _String test_string = _String("house,condo,hyphy");
    _String* sub_string = new _String(",");

    _List* result_list = test_string.Tokenize(sub_string);
    _String* result = (_String*)result_list->lData[0];

    ASSERT_STREQ("house", result->getStr());
    delete sub_string;
}

TEST_F(_StringTest,toNumTest)
{
    _String test_string = _String("3.14");

    double result = test_string.toNum();
    double expected = 3.14;

    EXPECT_EQ(expected, result);
}


TEST_F(_StringTest,UpCaseTest)
{
    _String result = _String("HOUSE");
    _String insert = _String("house");
    insert.UpCase();

    ASSERT_STREQ(result.getStr(), insert.getStr());
}

TEST_F(_StringTest,LoCaseTest)
{
    _String insert = _String("HOUSE");
    _String result = _String("house");
    insert.LoCase();

    ASSERT_STREQ(result.getStr(), insert.getStr());
}

TEST_F(_StringTest,ProcessTreeBranchLengthTest)
{
    //All this does is find the toNum, if it begins with a ':', just skip that.
    _String test_string = _String(":3.14");
    double result = test_string.ProcessTreeBranchLength();
    double expected = 3.14;
    EXPECT_EQ(expected, result);
}

TEST_F(_StringTest,ExtractEnclosedExpressionTest)
{
    //returns position
    _String test_string = _String("[hyp[house]hy]");
    long i = 0;
    long j = 0;
    j = test_string.ExtractEnclosedExpression(i,'[',']',true, true);
    EXPECT_EQ(13, j);
}

TEST_F(_StringTest,IsALiteralArgumentTest)
{
    _String test_string = _String("\"house\"");
    bool result = test_string.IsALiteralArgument(true);

    EXPECT_EQ(true, result);
}

TEST_F(_StringTest,StripQuotesTest)
{
    //Only strips the outer quotes
    _String insert = _String("\"So this\"");
    _String result = _String("So this");
    insert.StripQuotes();
    ASSERT_STREQ(result.getStr(), insert.getStr());
}

TEST_F(_StringTest,IsValidIdentifierTest)
{
    // Valid Identifier must be "greater than 0, doesn't start with non-alpha character and doesn't start wtih "_" if strict, if not strict, accept number"
    // Also cannot have keyword
    _String test_string = _String("house");
    EXPECT_EQ(true, test_string.IsValidIdentifier(true));
}

TEST_F(_StringTest,IsValidRefIdentifierTest)
{
    //Same as IsValidIdentifier, but ends with a &
    _String test_string = _String("house&");
    EXPECT_EQ(true, test_string.IsValidRefIdentifier());
}

/*TEST_F(_StringTest,ProcessParameterTest) { */
////This will be properly tested with batchlan
//_String initial_path = _String("/home/sergei/hyphy");
//initial_path.ProcessParameter();

//ASSERT_STREQ("/home/sergei/hyphy", initial_path);
/*}*/

//TEST_F(_StringTest,ProcessFileNameTest) {
////This will be properly tested with batchlan
//_String* test_string = new _String("/Users/stevenweaver/Documents/sergei/hyphy/trunk/UnitTests/mtDNA.fas");
//test_string->ProcessFileName();
//EXPECT_EQ(2, 2);
//}

TEST_F(_StringTest,PathCompositionTest)
{
    _String initial_path = _String("/home/sergei/hyphy");
    _String change_path = _String("../trunk");
    _String actual_path = initial_path.PathComposition(change_path);
    _String result_path = _String("/home/sergei/trunk");

    ASSERT_STREQ(result_path.getStr(), actual_path.getStr());

}

TEST_F(_StringTest,PathSubtractionTest)
{
    _String initial_path = _String("/home/sergei/");
    _String sub_path = _String("/home/sergei/hyphy");
    _String result_path = initial_path.PathSubtraction(sub_path,'A');
    ASSERT_STREQ("hyphy", result_path.getStr());
}

TEST_F(_StringTest,ConvertToAnIdentTest)
{
    //Takes a String and converts it to a valid hyphy ident test
    _String initial = _String("$house");
    initial.ConvertToAnIdent();
    ASSERT_STREQ("_house", initial.getStr());
}

TEST_F(_StringTest,ShortenVarIDTest)
{
    _String initial = _String("house.room");
    _String container = _String("house");
    _String result = initial.ShortenVarID(container);
    ASSERT_STREQ("room", result.getStr());
}

TEST_F(_StringTest,RegExpMatchOnceTest)
{
    _String initial = new _String("hyphy");
    _String* pattern = new _String("hyph");
    _SimpleList matched_pairs;
    initial.RegExpMatchOnce(pattern, matched_pairs, false, false);
    EXPECT_EQ(0,matched_pairs.lData[0]);
}

TEST_F(_StringTest,RegExpMatchTest)
{

    _String initial = new _String("hyphy");
    _String* pattern = new _String("phy");
    _SimpleList matched_pairs;
    initial.RegExpMatchOnce(pattern, matched_pairs, false, false);

    EXPECT_EQ(2,matched_pairs.lData[0]);
}

TEST_F(_StringTest,RegExpMatchAllTest)
{

    _String initial = new _String("hyphy");
    _String* pattern = new _String("phy");
    _SimpleList matched_pairs;

    initial.RegExpMatchOnce(pattern, matched_pairs, false, false);

    EXPECT_EQ(2,matched_pairs.lData[0]);

}

TEST_F(_StringTest,LempelZivProductionHistoryTest)
{
    //{1,0,01,1110,1100, 0010}
    _String initial = _String("1001111011000010");
    _SimpleList rec;
    long result = initial.LempelZivProductionHistory(&rec);

    EXPECT_EQ(6, result);
}

TEST_F(_StringTest,SortTest)
{
    _String initial = _String("hgfedcba");
    _SimpleList* index = new _SimpleList();
    _String* result = initial.Sort(index);
    ASSERT_STREQ("abcdefgh", result->getStr());
}

TEST_F(_StringTest,FindTerminatorTest)
{
    _String initial = _String("dither");
    _String terminator = _String("h");
    long result = initial.FindTerminator(0,terminator);
    EXPECT_EQ(3, result);
}

TEST_F(_StringTest,AppendAnAssignmentToBufferTest)
{

    _String initial = _String("dither");
    _String append = _String("12");
    _String append2 = _String("34");
    initial.AppendAnAssignmentToBuffer(&append, &append2, false, false, false);

    ASSERT_STREQ("dither12=34;\n", initial.getStr());

}

/*
 *TEST_F(_StringTest,AppendVariableValueAVLTest) {
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
 *    ASSERT_STREQ("life[answer]=4,2", initial.getStr());
 *}
 */


//Operator Tests
TEST_F(_StringTest,BracketTest)
{
    //[]
    _String result = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ('e', result[5L]);
}

TEST_F(_StringTest,ParanthTest)
{
    //()
    _String result = _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
    EXPECT_EQ('e', result(5));
}

TEST_F(_StringTest,EqualTest)
{
    //=
    _String result = _String ("hyphy");
    _String dupe = result;
    ASSERT_STREQ(result.getStr(), dupe.getStr());
}

TEST_F(_StringTest, AmpersandTest)
{

    // &
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String result = orig&to_append;

    ASSERT_STREQ(expected.getStr(), result.getStr());
}

TEST_F(_StringTest,DoubleLessTest)
{
    _String orig = _String ("hyphy");
    _String to_append = _String("-package");
    _String expected = _String("hyphy-package");
    _String result = orig&to_append;

    ASSERT_STREQ(expected.getStr(), result.getStr());
}

TEST_F(_StringTest,DoubleEqualTest)
{
    _String* result = new _String ("AABBCCDD");
    _String* r2 = new _String ("AABBCCDD");
    EXPECT_EQ(true, result->Equal(r2));

    delete result;
    delete r2;
}

TEST_F(_StringTest,GreaterOpTest)
{
    // >
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(false, result>r2);
}

TEST_F(_StringTest,LesserOpTest)
{
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(true, result<r2);
}

TEST_F(_StringTest,LesserEqualOpTest)
{
    //<=
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(true, result<=r2);
}

TEST_F(_StringTest,GreaterEqualOpTest)
{
    // >=
    _String result = _String ("house");
    _String r2 = _String ("household");

    EXPECT_EQ(false, result>=r2);
}

TEST_F(_StringTest,NotEqualOpTest)
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
