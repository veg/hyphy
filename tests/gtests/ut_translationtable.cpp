/*
 *  ut_translationtable.cpp
 *  HyPhyXCode
 *
 *
 */

#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include "ut_strings.h"

#include "hy_strings.h"
#include "translationtable.h"

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
class TranslationTableTest : public ::testing::Test
{
protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    TranslationTableTest() {
        // You can do set-up work for each test here.
    }

    virtual ~TranslationTableTest() {
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

TEST_F(TranslationTableTest,DefaultConstructorTest) {
    _TranslationTable test;
           
    EXPECT_EQ (HY_TRANSLATION_TABLE_STANDARD_NUCLEOTIDE, test.DetectType());
}
    
/******************************************/

TEST_F(TranslationTableTest,StringConstructorTest) {
    _String    test1           ("ABC"),
               test2           (*_TranslationTable::GetDefaultAlphabet(20));
    
    _TranslationTable test_tt  (test1),
                      test_tt2 (test2);
    
    EXPECT_EQ (HY_TRANSLATION_TABLE_NONSTANDARD, test_tt.DetectType());
    EXPECT_EQ (3, test_tt.LengthOfAlphabet());
    EXPECT_STREQ(test_tt.RetrieveCharacters()->sData, test1.sData);

    EXPECT_EQ (HY_TRANSLATION_TABLE_STANDARD_PROTEIN, test_tt2.DetectType());
    EXPECT_EQ (20, test_tt2.LengthOfAlphabet());
    EXPECT_STREQ(test_tt2.RetrieveCharacters()->sData, _TranslationTable::GetDefaultAlphabet(20)->sData);
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
