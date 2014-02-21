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
#include "hy_list_numeric.h"

namespace {

class _testNumericPayload {
  public:
  _testNumericPayload (void) { data [0] = 0; data [1] = 0;}
  _testNumericPayload (unsigned long p) { data[0] = p; data[1] = p; }
  _testNumericPayload (unsigned long p, unsigned long p2) { data[0] = p; data[1] = p2; }
  _testNumericPayload (const _testNumericPayload& o) { data[0] = o.data[1]; data[1] = o.data[1]; }
  
  bool operator == (const _testNumericPayload & o) const { return data[0] == o.data[0] && data[1] == o.data[1];} 
  bool operator != (const _testNumericPayload & o) const { return ! ((*this) == o);} 
  bool operator <  (const _testNumericPayload & o) const { if (data[0] < o.data[0]) return true; if (data [0] == o.data[0]) {return data[1] < o.data[1];} return false;} 
  bool operator <= (const _testNumericPayload & o) const { return (*this) < o || (*this == o); } 
  bool operator >  (const _testNumericPayload & o) const { return  ! ((*this) <= o); } 
  bool operator >= (const _testNumericPayload & o) const { return  ! ((*this) < o); } 
  operator unsigned long (void) {return data[0];}
  
  unsigned long data[2];
};

// The fixture for testing class Foo.
template <typename DATA>
class _hyListNumericTest : public ::testing::Test {

protected:

  _hyListNumericTest() {
    
  }

  virtual ~_hyListNumericTest() {
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
  
  _hyListNumeric <DATA> make_random_list (const unsigned long size, const unsigned long range, const unsigned long offset = 0UL) {
    _hyListNumeric <DATA> random_list;
    for (unsigned long item = 0UL; item < size; item ++) {
      random_list.append ((DATA) (offset + genrand_int32() % range));
    }
    return random_list;
  }
  
  _StringBuffer dump_to_stream_as_longs (const _hyListNumeric<DATA>& data) {
    _StringBuffer result;
    char buffer [256];
    result << '[';
    for (unsigned long item = 0UL; item < data.countitems(); item++) {
      if (item) {
        result << ',';
      }
      snprintf(buffer,255,"%lu", (unsigned long)data.AtIndex (item));
      result << buffer;
    }
    result << ']';
    return result;
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

TYPED_TEST_CASE_P(_hyListNumericTest);

TYPED_TEST_P (_hyListNumericTest, ConstuctorTests) {
  
  _hyListNumeric <TypeParam> null_list,
                               single_element_list((TypeParam)16),
                               arithmetic_series_list(10UL, (TypeParam)4UL, (TypeParam)4UL),
                               full_stack_copy(arithmetic_series_list),
                               partial_stack_copy(arithmetic_series_list,2,HY_LIST_INSERT_AT_END);
          
  ASSERT_EQ (0UL, null_list.countitems()) << "Non-empty null list";
  ASSERT_EQ (1UL, single_element_list.countitems()) << "Single element list has wrong length";
  ASSERT_EQ (10UL, arithmetic_series_list.countitems()) << "Arithmetic Series of elements list has wrong length";
  ASSERT_EQ (40UL, arithmetic_series_list[9]) << "Incorrect arithmetic calculation";
  ASSERT_EQ (10UL, full_stack_copy.countitems()) << "Stack copy list has wrong length";
  ASSERT_EQ (8UL, partial_stack_copy.countitems()) << "Partial stack copy list has wrong length";

  

  EXPECT_EQ (single_element_list (0), arithmetic_series_list[3]);   

  for (unsigned long i = 0UL; i < arithmetic_series_list.countitems(); i++) {
    EXPECT_EQ (full_stack_copy(i), arithmetic_series_list[i]);
  }
  
  EXPECT_TRUE(full_stack_copy == arithmetic_series_list) << "Failed list == list";
  EXPECT_FALSE(full_stack_copy == partial_stack_copy) << "Failed list != list";
  
  EXPECT_EQ (partial_stack_copy(0), arithmetic_series_list(2));
  EXPECT_EQ (partial_stack_copy(2), full_stack_copy(4));
  
  partial_stack_copy.Clear(false);
  EXPECT_EQ (0UL, partial_stack_copy.countitems());
  EXPECT_LT (0UL, partial_stack_copy.allocated());
  partial_stack_copy.Clear(true);
  EXPECT_EQ (0UL, partial_stack_copy.allocated());
}

TYPED_TEST_P (_hyListNumericTest, MethodTests) {
  
  _hyListNumeric <TypeParam> null_list,
                               single_element_list((TypeParam)16),
                               arithmetic_series_list(10UL, (TypeParam)4UL, (TypeParam)4UL),
                               full_stack_copy(arithmetic_series_list),
                               partial_stack_copy(arithmetic_series_list,2,HY_LIST_INSERT_AT_END);
             
    long direct_index [6] = {4L,8L,12L,16L,20L,24L};

                                
    // test some basic element tests
    ASSERT_EQ ( direct_index[4], arithmetic_series_list.Element(4)) << "Failed retrieving the correct element";
    ASSERT_EQ ( direct_index[1], arithmetic_series_list.Element(-9)) << "Failed retrieving the correct element";
    ASSERT_EQ ( _HY_LIST_NUMERIC_INVALID_VALUE_, arithmetic_series_list.Element(12)) << "Should be null";
    ASSERT_EQ ( _HY_LIST_NUMERIC_INVALID_VALUE_, null_list.Element(0)) << "Failed retrieving the correct element";
    
    // summation
    ASSERT_EQ ( 220, arithmetic_series_list.Sum()) << "Failed summation of a list";
    ASSERT_EQ ( 16, single_element_list.Sum()) << "Failed summation of a list with a single element";
    ASSERT_EQ ( 0, null_list.Sum()) << "Failed summation of an empty list";

    // toStr tests
    //TODO: Fix string doubles and longs to strings
    _StringBuffer* arithmetic_string = (_StringBuffer*)arithmetic_series_list.toStr();
    _StringBuffer* single_element_string = (_StringBuffer*)single_element_list.toStr();
    _StringBuffer* null_string = (_StringBuffer*)null_list.toStr();
    EXPECT_STREQ("", arithmetic_string->getStr()) << "multiple numeric list to string failed";
    EXPECT_STREQ("", single_element_string->getStr()) << "single numeric list to string failed";
    EXPECT_STREQ("{}", null_string->getStr()) << "empty numeric list to string failed";

    // list to partition string test
    //TODO: Fix string doubles and longs to strings
    _StringBuffer* arithmetic_partition_string = (_StringBuffer*)arithmetic_series_list.ListToPartitionString();
    _StringBuffer* single_element_partition_string = (_StringBuffer*)single_element_list.ListToPartitionString();
    _StringBuffer* null_partition_string = (_StringBuffer*)null_list.ListToPartitionString();
    EXPECT_STREQ("", arithmetic_partition_string->getStr()) << "multiple numeric list to partition string failed";
    EXPECT_STREQ("", single_element_partition_string->getStr()) << "single numeric list to partition string failed";
    EXPECT_STREQ("", null_partition_string->getStr()) << "empty numeric list to partition string failed";

    // offset test
    arithmetic_series_list.Offset(24);
    single_element_list.Offset(24);
    null_list.Offset(24);
    ASSERT_EQ ( 36, arithmetic_series_list.Element(2)) << "Failed offset of a list";
    ASSERT_EQ ( 40, single_element_list.Element(0)) << "Failed offset of a list with a single element";
    ASSERT_EQ ( _HY_LIST_NUMERIC_INVALID_VALUE_, null_list.Element(0)) << "Failed offset of an empty list";

    // counting sort
    //for (int i; i < 12; i++) {
    //  if (i % 6)
    //    single_element_list.BinaryInsert((TypeParam)2UL);
    //  if (i % 4) 
    //    single_element_list.BinaryInsert((TypeParam)3UL);
    //  if (i % 3)
    //    single_element_list.BinaryInsert((TypeParam)4UL);
    //  if (i % 2)
    //    single_element_list.BinaryInsert((TypeParam)6UL);
    //  single_element_list.BinaryInsert((TypeParam)12UL);
    //}

    //single_element_list.CountingSort(10);
    //ASSERT_EQ ( 12UL, single_element_list[0]) << "Failed summation of a list";
    //ASSERT_EQ ( 6UL, single_element_list[1]) << "Failed summation of a list";
    //ASSERT_EQ ( 4UL, single_element_list[2]) << "Failed summation of a list with a single element";
}

REGISTER_TYPED_TEST_CASE_P (_hyListNumericTest, ConstuctorTests, MethodTests);

typedef ::testing::Types<long, double> _hyListNumericTestTypes;
//typedef ::testing::Types<long> _hyListNumericTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(_typedList, _hyListNumericTest, _hyListNumericTestTypes);
}
