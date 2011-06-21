#include <iostream>
#include "hy_strings.h"
#include "gtest/gtest.h"

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

// Tests that the Foo::Bar() method does Abc.
TEST_F(_StringTest, LengthTest) {
  _String* result;
  _String oops ("You're asking me to run MCMC without reporting any results.  Did you forget to set Bgm_MCMC_SAMPLES?\n");
  EXPECT_EQ(128, 128);
}


//// Tests that Foo does Xyz.
//TEST_F(_StringTest, DoesXyz) {
  //// Exercises the Xyz feature of Foo.
//}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}