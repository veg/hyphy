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
#include <time.h>

//Generate necessary includes from the respective implementation file
#include "hy_stack.h"
#include "hy_list_reference.h"
#include "hy_strings.h"

namespace {

// The fixture for testing class Foo.
template <typename DATA>
class _hyStackTest : public ::testing::Test {

protected:

  _hyStackTest() {
    
  }

  virtual ~_hyStackTest() {
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

typedef ::testing::Types<char, long, double> _hyStackTypes;
TYPED_TEST_CASE(_hyStackTest, _hyStackTypes);
TYPED_TEST (_hyStackTest, ConstuctorTests) {
  _hyStack <TypeParam> null_list;
  ASSERT_EQ (0UL, null_list.Length()) << "Non-empty null list";
}

TYPED_TEST(_hyStackTest, MethodTests) {
  _hyStack <TypeParam> multiple_element_stack,
                       null_stack;

  TypeParam array[4] = {(TypeParam)1L, (TypeParam)4L, (TypeParam)9L, (TypeParam)16L};
  int array_size = sizeof(array)/sizeof(TypeParam); 

  //push
  for (int i = 0; i < array_size; i++) {
    multiple_element_stack.push(array[i]);
  }

  ASSERT_EQ (array_size, multiple_element_stack.Length()) << "All elements weren't pushed to the stack";

  //pop
  TypeParam result =  multiple_element_stack.pop();
  ASSERT_EQ (array[array_size - 1 ], result) << "Incorrect element popped from the stack";
  ASSERT_EQ (array_size - 1 , multiple_element_stack.Length()) << "Element was not popped off the stack";

  //reset
  multiple_element_stack.reset();
  ASSERT_EQ (0, multiple_element_stack.Length()) << "Stack should be empty after reset";

}

// The fixture for testing class Foo.
template <typename DATA>
class _hyStackReferenceTest : public ::testing::Test {

protected:

  _hyStackReferenceTest() {
    
  }

  virtual ~_hyStackReferenceTest() {
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

typedef ::testing::Types<_String> _hyStackReferenceTypes;
TYPED_TEST_CASE(_hyStackReferenceTest, _hyStackReferenceTypes);
TYPED_TEST(_hyStackReferenceTest, ReferenceConstructorTests) {

  _hyStack <TypeParam, _hyListReference<TypeParam> > string_stack;

 for (long k = 0; k < 50; k ++) {
    string_stack.push(new TypeParam (k));
  }


  ASSERT_EQ (50, string_stack.Length()) << "All elements weren't pushed to the stack";

  //pop
  TypeParam result =  string_stack.pop();
  EXPECT_TRUE(result.Equal(TypeParam (49L)));
  ASSERT_EQ (49, string_stack.Length()) << "Incorrect element popped from the stack";

  //reset
  string_stack.reset();
  ASSERT_EQ (0, string_stack.Length()) << "Stack should be empty after reset";

}

