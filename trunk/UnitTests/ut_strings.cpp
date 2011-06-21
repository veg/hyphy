/*
 *  ut_strings.cpp
 *  HyPhyXCode
 *
 *  Created by Steven Weaver on 6/17/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "gtest/gtest.h"
#include "ut_strings.h"

#include <iostream>
#include "hy_strings.h"

namespace {

// The fixture for testing class Foo.
class _StringTest : public ::testing::Test {
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

//// Tests that Foo does Xyz.
//TEST_F(_StringTest, DoesXyz) {
  //// Exercises the Xyz feature of Foo.
//}

//Stub out the Rest of the tests
TEST_F(_StringTest,makeDynamicTest) { 
  //What is the difference between this and dupicate?
  //_String* result = new _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
  //EXPECT_EQ(101, result->Length());

}

TEST_F(_StringTest,getCharTest) { 
  _String* result = new _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
  EXPECT_EQ('e', result->getChar(5));
}


TEST_F(_StringTest,setCharTest) { 
  _String* result = new _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
  result->getChar(5,'d')
  EXPECT_EQ('d', result->getChar(5));
}

TEST_F(_StringTest,DuplicateErasingTest) { 
//Not sure what this does

}

TEST_F(_StringTest, LengthTest) {
  _String* result = new _String ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
  EXPECT_EQ(101, result->Length());
}

TEST_F(_StringTest,InsertTest) { 

  _String result = "AAGGCCTTA";
  result.Insert('C',0)
  EXPECT_EQ("CAAGGCCTTA", result);
}

TEST_F(_StringTest,DeleteTest) { 
  _String result = "AAGGCCTTA";
  result.Delete(3,4)
  EXPECT_EQ("AAGCTTA", result);
}

TEST_F(_StringTest,AppendNewInstanceTest) { 
//This checks the << operator
//TODO: I don't know what this does

}

TEST_F(_StringTest,EscapeAndAppendTest) { 
//It escapes the character and then appends it to the end of the string.

  _String result = "AAGGCCTTA";
  result.Delete(3,4)
  EXPECT_EQ("AAGCTTA", result);

}

TEST_F(_StringTest,EscapeAndAppendTest) { 

  _String result = "AAGGCCTTA";
  result->EscapeAndAppend;
  EXPECT_EQ("AAGCTTA", result);

}

TEST_F(_StringTest,FinalizeTest) { 

}

TEST_F(_StringTest,getStrTest) { 

}

TEST_F(_StringTest,toStrTest) { }
TEST_F(_StringTest,ChopTest) { }
TEST_F(_StringTest,CutTest) { }
TEST_F(_StringTest,FlipTest) { }
TEST_F(_StringTest,Adler32Test) { }
TEST_F(_StringTest,TrimTest) { }
TEST_F(_StringTest,FirstNonSpaceIndexTest) { }
TEST_F(_StringTest,KillSpacesTest) { }
TEST_F(_StringTest,CompressSpacesTest) { }
TEST_F(_StringTest,FirstSpaceIndexTest) { }
TEST_F(_StringTest,FirstNonSpaceTest) { }
TEST_F(_StringTest,FindEndOfIdentTest) { }
TEST_F(_StringTest,FindTest) { }
TEST_F(_StringTest,FindAnyCaseTest) { }
TEST_F(_StringTest,ContainsSubstringTest) { }
TEST_F(_StringTest,FindTest) { }
TEST_F(_StringTest,FindBackwardsTest) { }
TEST_F(_StringTest,FindBinaryTest) { }
TEST_F(_StringTest,EqualTest) { }
TEST_F(_StringTest,CompareTest) { }
TEST_F(_StringTest,EqualWithWildCharTest) { }
TEST_F(_StringTest,GreaterTest) { }
TEST_F(_StringTest,LessTest) { }
TEST_F(_StringTest,containsTest) { }
TEST_F(_StringTest,beginswithTest) { }
TEST_F(_StringTest,startswithTest) { }
TEST_F(_StringTest,endswithTest) { }
TEST_F(_StringTest,FormatTimeStringTest) { }
TEST_F(_StringTest,ReplaceTest) { }
TEST_F(_StringTest,TokenizeTest) { }
TEST_F(_StringTest,toNumTest) { }
TEST_F(_StringTest,UpCaseTest) { }
TEST_F(_StringTest,LoCaseTest) { }
TEST_F(_StringTest,ProcessTreeBranchLengthTest) { }
TEST_F(_StringTest,IsALiteralArgumentTest) { }
TEST_F(_StringTest,StripQuotesTest) { }
TEST_F(_StringTest,IsValidIdentifierTest) { }
TEST_F(_StringTest,IsValidRefIdentifierTest) { }
TEST_F(_StringTest,ProcessParameterTest) { }
TEST_F(_StringTest,ProcessFileNameTest) { }
TEST_F(_StringTest,PathCompositionTest) { }
TEST_F(_StringTest,GetPlatformDirectoryCharTest) { }
TEST_F(_StringTest,PathSubtractionTest) { }
TEST_F(_StringTest,ConvertToAnIdentTest) { }
TEST_F(_StringTest,ShortenVarIDTest) { }
TEST_F(_StringTest,GetVersionStringTest) { }
TEST_F(_StringTest,GetTimeStampTest) { }
TEST_F(_StringTest,RegExpMatchOnceTest) { }
TEST_F(_StringTest,RegExpMatchTest) { }
TEST_F(_StringTest,RegExpMatchAllTest) { }
TEST_F(_StringTest,PrepRegExpTest) { }
TEST_F(_StringTest,GetRegExpErrorTest) { }
TEST_F(_StringTest,FlushRegExpTest) { }
TEST_F(_StringTest,LempelZivProductionHistoryTest) { }
TEST_F(_StringTest,SortTest) { }
TEST_F(_StringTest,ExtractEnclosedExpressionTest) { }
TEST_F(_StringTest,FindTerminatorTest) { }
TEST_F(_StringTest,AppendAnAssignmentToBufferTest) { }
TEST_F(_StringTest,AppendVariableValueAVLTest) { }
}  // namespace
