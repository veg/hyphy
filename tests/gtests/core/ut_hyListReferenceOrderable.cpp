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
#include "hy_list_reference_orderable.h"

namespace {

class _testPayload: public BaseObj {

  public:
  
  _testPayload (void) { data = 0UL;}
  _testPayload (const char* p) { data = atoi (p); }
  _testPayload (const long p) { data = p; }
  _testPayload (const _testPayload& o) { data = o.data; }
  virtual BaseObj *makeDynamic(void) const { return new _testPayload (*this); }
  virtual void Duplicate(BaseObj const * ref) { data = ((_testPayload*)ref)->data; }
  virtual long Compare (const _testPayload * o) {return data > o->data ? 1 : (data < o->data ? -1 : 0);}
  virtual bool Equal (_testPayload const & o) {return data == o.data;}

  long data;
};

// The fixture for testing class Foo.
template <typename DATA>
class _hyListReferenceOrderableTest : public ::testing::Test {

protected:

  _hyListReferenceOrderableTest() {
    
  }

  virtual ~_hyListReferenceOrderableTest() {
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


TYPED_TEST_CASE_P(_hyListReferenceOrderableTest);


TYPED_TEST_P (_hyListReferenceOrderableTest, ConstructorTests) {
  
   TypeParam array [5] = {(TypeParam)"1", (TypeParam)"10", (TypeParam)"2", (TypeParam)"3", (TypeParam)"7"},
             *array2 [5] = {array, array + 1L, array + 2L, array + 3L, array + 4L},
             testValue = (TypeParam)"3",
             * dynamic_object = new TypeParam ("3.14");
             
  {
  
    TypeParam         *another_dynamic_object = new TypeParam ("2.7172");
    
    _hyListReferenceOrderable <TypeParam> null_list,
                                 single_element_list (testValue),
                                 multiple_element_list (5,array2),
                                 full_stack_copy (multiple_element_list),
                                 partial_stack_copy (multiple_element_list,2,HY_LIST_INSERT_AT_END),
                                 * dynamic_list = new _hyListReferenceOrderable <TypeParam> (dynamic_object),
                                 * another_dynamic_list = new _hyListReferenceOrderable <TypeParam> (another_dynamic_object);
            
    ASSERT_EQ (0UL, null_list.countitems()) << "Non-empty null list";
    ASSERT_EQ (1UL, single_element_list.countitems()) << "Single element list has wrong length";
    ASSERT_EQ (5UL, multiple_element_list.countitems()) << "Array of elements list has wrong length";
    ASSERT_EQ (5UL, full_stack_copy.countitems()) << "Stack copy list has wrong length";
    ASSERT_EQ (3UL, partial_stack_copy.countitems()) << "Partial stack copy list has wrong length";
    
    EXPECT_TRUE (single_element_list (0)->Equal (*full_stack_copy(3)));
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
    EXPECT_EQ (dynamic_object, dynamic_list->Element (0));
    EXPECT_TRUE (dynamic_object->SingleReference()) << "Incorrect reference count following list operations";
    
    EXPECT_TRUE (another_dynamic_list->SingleReference()) << "Incorrect reference count following list operations";
    
    full_stack_copy.Clear();
    delete dynamic_list;
    delete another_dynamic_list;
  }
  
  for (unsigned long k = 0; k < 5L; k++) {
    EXPECT_TRUE (array[k].SingleReference()) << "Incorrect reference count following list operations";
  }
  EXPECT_TRUE (testValue.SingleReference()) << "Incorrect reference count following list operations";
}

TYPED_TEST_P (_hyListReferenceOrderableTest, SortingTests) {
  

   TypeParam array [6] = {(TypeParam)5L, (TypeParam)2L, (TypeParam)3L, (TypeParam)2L, (TypeParam)9L, (TypeParam)1L},
             *array_ref [6] = {array, array + 1L, array + 2L, array + 3L, array + 4L, array + 5L},
             array_sorted [6] = {(TypeParam)1L, (TypeParam)2L, (TypeParam)2L, (TypeParam)3L, (TypeParam)5L, (TypeParam)9L},
             *array_sorted_ref [6] = {array_sorted, array_sorted + 1L, array_sorted + 2L, array_sorted + 3L, array_sorted + 4L, array_sorted + 5L};
             
             
   long      direct_index [6] = {0L,1L,2L,3L,4L,5L};
             
   _hyListReferenceOrderable  <TypeParam> unsorted      (6UL, array_ref),
                                          sorted        (6UL, array_sorted_ref),
                                          sort_me       (unsorted),
                                          reverse_large,
                                          direct_large, 
                                          empty;
                                
   _hyListOrderable <long>      bubble_sort_indexer (6UL, direct_index),
                                quick_sort_indexer (6UL, direct_index);
                          
    // test some basic comparison operations 
    
    ASSERT_EQ (HY_COMPARE_LESS, sorted.Compare (0,4)) << "Failed a<b element comparison";
    ASSERT_EQ (HY_COMPARE_EQUAL, sorted.Compare (1,2)) << "Failed a==b element comparison";
    ASSERT_EQ (HY_COMPARE_GREATER, sorted.Compare (4,0)) << "Failed a>b element comparison";
    
    ASSERT_TRUE (sorted.IsSorted()) << "IsSorted returned false for a sorted list";
    ASSERT_FALSE (unsorted.IsSorted()) << "IsSorted returned true for an unsorted list";
    
    empty.Sort();
    ASSERT_TRUE (empty.IsSorted()) << "Empty list claimed to be unsorted";
    
    
    sort_me.BubbleSort (&bubble_sort_indexer);
    for (unsigned long i = 0UL; i < sort_me.countitems(); i++) {
      EXPECT_TRUE (sorted (i) -> Equal (*sort_me[i])) << "Bubble sort failed at element " << (i+1);
      EXPECT_TRUE (unsorted (bubble_sort_indexer.AtIndex(i)) ->Equal (*sort_me.AtIndex (i))) << "Sorted indices for bubble sort are incorrect at position " << i << " mapping to index " << bubble_sort_indexer.AtIndex(i);
    }
    sort_me.Clone (&unsorted);
    sort_me.QuickSort (0, sort_me.countitems(), &quick_sort_indexer);
    for (unsigned long i = 0UL; i < sort_me.countitems(); i++) {
      EXPECT_TRUE (sorted (i) ->Equal (*sort_me[i])) << "Quick sort failed at element " << (i+1);
      EXPECT_TRUE (unsorted (quick_sort_indexer.AtIndex(i)) -> Equal (*sort_me.AtIndex (i))) << "Sorted indices for quick sort are incorrect at position " << i << " mapping to index " << quick_sort_indexer.AtIndex(i);
    }
    
    _hyListOrderable <long> index_list;
    
    for (unsigned long i = 0; i < 10; i++) {
      direct_large.AppendNewInstance (new TypeParam(i));
      reverse_large.AppendNewInstance (new TypeParam(9UL-i));
      index_list.append (i);
    }
    
    _hyListReferenceOrderable <TypeParam> copy_reverse (reverse_large),
                                 copy_direct  (direct_large);
                                 
    copy_direct.Sort();
    ASSERT_EQ(copy_direct, direct_large) << "Sorting an already sorted list failed";
    reverse_large.QuickSort(0, reverse_large.countitems()-1, &index_list);
    for (unsigned long i = 0UL; i < copy_reverse.countitems(); i++) {
      EXPECT_TRUE (direct_large.AtIndex(i)->Equal (*reverse_large.AtIndex(i))) << "Sort reverse list in ascending order at element " << (i+1);
    }
    reverse_large.Sort(false, &index_list);
    ASSERT_EQ (reverse_large, copy_reverse) << "Sort reverse list in descending order failed";
    
    // random tests

    
}


REGISTER_TYPED_TEST_CASE_P (_hyListReferenceOrderableTest, ConstructorTests, SortingTests);
  
}

typedef ::testing::Types<_testPayload, _String> _hyListReferenceOrderableTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(_typedList, _hyListReferenceOrderableTest, _hyListReferenceOrderableTestTypes);
