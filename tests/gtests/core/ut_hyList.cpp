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

//Generate necessary includes from the respective implementation file
#include "hy_list.h"

namespace {

class _testPayload {
  public:
  
  _testPayload (unsigned long p) { data = p; unused = 0UL; }
  _testPayload (const _testPayload& o) { data = o.data; }
  
  bool operator == (const _testPayload & o) const { return data == o.data;} 
  
  unsigned long data, unused;
};

// The fixture for testing class Foo.
template <typename DATA>
class _hyListTest : public ::testing::Test {

protected:

  _hyListTest() {
    
  }

  virtual ~_hyListTest() {
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

_hyList <DATA> test_list;

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

TYPED_TEST_CASE_P(_hyListTest);


TYPED_TEST_P (_hyListTest, ConstuctorTests) {
  
   TypeParam array [5] = {(TypeParam)1, (TypeParam)4, (TypeParam)9, (TypeParam)16, (TypeParam)25};
  
  _hyList <TypeParam> null_list,
          single_element_list ((TypeParam)16),
          multiple_element_list (5,array),
          full_stack_copy (multiple_element_list),
          partial_stack_copy (multiple_element_list,2,HY_LIST_INSERT_AT_END);
          
  ASSERT_EQ (0UL, null_list.countitems()) << "Non-empty null list";
  ASSERT_EQ (1UL, single_element_list.countitems()) << "Single element list has wrong length";
  ASSERT_EQ (5UL, multiple_element_list.countitems()) << "Array of elements list has wrong length";
  ASSERT_EQ (5UL, full_stack_copy.countitems()) << "Stack copy list has wrong length";
  ASSERT_EQ (3UL, partial_stack_copy.countitems()) << "Partial stack copy list has wrong length";
  
  EXPECT_EQ (single_element_list (0), full_stack_copy(3));   
  for (unsigned long i = 0UL; i < multiple_element_list.countitems(); i++) {
    EXPECT_EQ (full_stack_copy (i), full_stack_copy[i]);
  }
  
  EXPECT_TRUE(full_stack_copy == multiple_element_list) << "Failed list == list";
  EXPECT_FALSE(full_stack_copy == partial_stack_copy) << "Failed list != list";
  
  EXPECT_EQ (partial_stack_copy (0), multiple_element_list (2));
  EXPECT_EQ (partial_stack_copy (2), full_stack_copy (4));
  
  partial_stack_copy.Clear (false);
  EXPECT_EQ (0UL, partial_stack_copy.countitems());
  EXPECT_LT (0UL, partial_stack_copy.allocated());
  partial_stack_copy.Clear (true);
  EXPECT_EQ (0UL, partial_stack_copy.allocated());
  
  _hyList <TypeParam> * dynamicList1 = new _hyList<TypeParam> (multiple_element_list),
  * dynamicList2 = dynamic_cast <_hyList<TypeParam>*> (multiple_element_list.makeDynamic());
  
  EXPECT_TRUE(*dynamicList1 == multiple_element_list) << "Failed new (list) == list";
  EXPECT_TRUE(*dynamicList2 == multiple_element_list) << "Failed makeDynamic (list) == list";
  
  partial_stack_copy.Duplicate (&multiple_element_list);
  EXPECT_TRUE(partial_stack_copy == multiple_element_list) << "Failed Duplicate(list) == list";
  
  delete dynamicList1;
  DeleteObject (dynamicList2);
  
  
  
}

TYPED_TEST_P (_hyListTest, AccessAndManipulationTests) {
  _hyList <TypeParam>  sequence,
                       sequence2,
                       sequence3,
                       sequence4;
  
  for (long i = 1L; i <= 10L; i++) {
    sequence.append ((TypeParam) i);
    sequence2.append ((TypeParam) i);
  }
  
  EXPECT_EQ (10UL, sequence.countitems()) << "A sequence of appends failed to generate the right list length";
  EXPECT_GE (sequence.allocated(), sequence.countitems()) << "The allocated space is less than the used space";
  
  for (long i = 1L; i <= 10L; i++) {
    sequence << (TypeParam) (i*i);
  }
  
  for (unsigned long idx = 0UL; idx < sequence.countitems(); idx++) {
      EXPECT_EQ (sequence.AtIndex (idx), sequence(idx)) << "AtIndex != operator ()";
  }
  
  EXPECT_EQ (20UL, sequence.countitems()) << "A sequence of << calls failed to generate the right list length";
  
  sequence >> sequence (0);
  EXPECT_EQ (20UL, sequence.countitems()) << ">> added an already existing element";
  
  sequence >> (TypeParam)42;
  EXPECT_EQ (21UL, sequence.countitems()) << ">> failed to add a new element";
  
  sequence << sequence2;
  EXPECT_EQ (31UL, sequence.countitems()) << "A << call with a list argument failed to generate the right list length";
  
  EXPECT_EQ (sequence[100], sequence2.Element (-1)) << "Accessing past the end of the list or with negative coordinates failed";
  
  sequence3 = sequence3 & sequence4;
  EXPECT_EQ (0UL, sequence3.countitems()) << "The concatenation of two empty lists is not empty";
  sequence3 = sequence & sequence2;
  sequence4 = sequence;
  sequence4 << sequence2;
  
  EXPECT_TRUE (sequence3 == sequence4) << "Failed test (C = A & B) == (C = A) << B";
  
  sequence4.SetItem (3, (TypeParam)123);
  sequence3.SetItem (8, (TypeParam)123);
  
  EXPECT_EQ (sequence4.AtIndex (3), sequence3.AtIndex (8)) << "SetItem failed to assign equal values to list items";
  
  
}
  
TYPED_TEST_P (_hyListTest, FindAddDelete) {
  
  _hyList <TypeParam>  list,
                       list2;
  
  TypeParam key1 = (TypeParam)4,
            key2 = (TypeParam)12;
  
  for (long i = 1L; i <= 10L; i++) {
    list << (TypeParam) i;
  }
  
  list2.Clone (&list);
  EXPECT_EQ(3L, list.Find (key1)) << "Failed to find the correct element";
  EXPECT_EQ(HY_NOT_FOUND, list.Find (key2)) << "Incorrectly found an element";
  EXPECT_EQ(HY_NOT_FOUND, list.FindStepping (key1, 2)) << "Incorrectly found an element during a stepping search";
  EXPECT_EQ(3L, list.FindStepping (key1, 2, 1)) << "Found an element during a stepping search";
  
    
}


REGISTER_TYPED_TEST_CASE_P (_hyListTest, ConstuctorTests, AccessAndManipulationTests, FindAddDelete);

}

typedef ::testing::Types<char, long, double, _testPayload> _hyListTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(_typedList, _hyListTest, _hyListTestTypes);
