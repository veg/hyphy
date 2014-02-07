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
#include "hy_list.h"

namespace {

class _testPayload {
  public:
  
  _testPayload (void) { data = 0UL; unused = 0UL;}
  _testPayload (unsigned long p) { data = p; unused = 0UL; }
  _testPayload (const _testPayload& o) { data = o.data; unused = o.unused; }
  
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
  
  TypeParam e1 = list2.Pop(),
            e2 = list2.Pop();
            
  EXPECT_EQ(list.Element (9UL), e1) << "Popping last list element failed to return the right value";
  EXPECT_EQ(list.Element (8UL), e2) << "Popping last list element failed to return the right value";
  EXPECT_EQ (list.countitems()-2, list2.countitems()) << "Popping two elements did not shoren the list by two";
  
  list2.Clear();
  for (long i = 0L; i < 3*HY_LIST_ALLOCATION_CHUNK; i++) {
    list2 << (TypeParam) i;
  }
  
  list2.Delete (0, false);
  for (long i = HY_LIST_ALLOCATION_CHUNK; i < 2*HY_LIST_ALLOCATION_CHUNK; i++) {
    list2.Delete (HY_LIST_ALLOCATION_CHUNK, true);
  }
      
  EXPECT_EQ(2*HY_LIST_ALLOCATION_CHUNK-1, list2.countitems() ) << "List deletetion operation did not return the right length list";
  for (long i = 0L; i < HY_LIST_ALLOCATION_CHUNK-1; i++) {
    EXPECT_EQ((TypeParam)(i+1), list2.Element (i)) << "List deletion operation failed";
  }
  for (long i = HY_LIST_ALLOCATION_CHUNK; i < 2*HY_LIST_ALLOCATION_CHUNK-1; i++) {
    EXPECT_EQ((TypeParam)(HY_LIST_ALLOCATION_CHUNK+i+1), list2.Element (i)) << "List deletion operation failed";
  }
  
  list2.TrimMemory ();
  EXPECT_EQ (list2.countitems(), list2.allocated()) << "Trim memory did not correctly trim the allocated buffer";
  
  _hyList <long> indicesToDeleteList;
  for (long i = 0; i < HY_LIST_ALLOCATION_CHUNK - 1; i+=2) { 
    indicesToDeleteList.append (i);
  }
  indicesToDeleteList << list2.countitems()-1;
  indicesToDeleteList << list2.countitems();
  
  long run1 = HY_LIST_ALLOCATION_CHUNK/2;
  
  list2.DeleteList (&indicesToDeleteList); 
  indicesToDeleteList.Clear();

  for (long i = 0; i < list2.countitems(); i++) {
    EXPECT_EQ ((TypeParam)(i < run1 ? 2*(i+1) : 2*HY_LIST_ALLOCATION_CHUNK + (i-run1) + 1) , list2.Element (i)) << "DeleteList did not correctly delete the elements";
    indicesToDeleteList << i;
  }
  
  list2.RequestSpace (HY_LIST_ALLOCATION_CHUNK * 2);
  list2.DeleteList (&indicesToDeleteList); 
  EXPECT_EQ(0UL, list2.countitems()) << "Should have seen an empty list after all the elements have been deleted";
  EXPECT_EQ(0UL, list2.allocated()) << "Should have seen an empty list allocation after all the elements have been deleted";
  
  
}

TYPED_TEST_P (_hyListTest, InPlaceOperations) {
  _hyList <TypeParam> list,
                      list2;
                      
  for (long i = 0L; i < 100L; i++) {
    list.append ((TypeParam) i);
  }
  
  list2 = list;
  list.Flip();
  list.Flip();
  EXPECT_TRUE(list == list2) << "Flip is not reversible for an even sized list";
  list.InsertElement ((TypeParam)256, 1L);
  list2 = list;
  list.Flip();
  list.Flip();
  EXPECT_TRUE(list == list2) << "Flip is not reversible for an odd sized list";

  list.InsertElement ((TypeParam)256, 99L);
  EXPECT_EQ(list.Element(1L), list.Element (99L)) << "Insertion of the same element in two different places failed";
  
  list2 = list;
  list.Displace (10,20,2);
  EXPECT_FALSE (list == list2) << "Displace produced a list equal to the original";
  list.Displace (12,22,-2);
  EXPECT_TRUE (list == list2) << "Displace roundtrip failed";
  list.Displace (90,95,10);
  list.Displace (96,HY_LIST_INSERT_AT_END,-6);
  EXPECT_TRUE (list == list2) << "Displace roundtrip failed when hitting the right edge";
  list.Displace (5,10,-20);
  list.Displace (HY_LIST_INSERT_AT_END,5,5);
  EXPECT_TRUE (list == list2) << "Displace roundtrip failed when hitting the left edge";
  
  list2.Swap (5,10);
  EXPECT_EQ (list(10), list2(5)) << "Element swap failed";
  EXPECT_EQ (list(5), list2(10)) << "Element swap failed";
  
  
}

TYPED_TEST_P (_hyListTest, SubsetsAndSelectors) {
  time_t timer;
  init_genrand (time(&timer));
  
  TypeParam array [4] = {(TypeParam)1L, (TypeParam)4L, (TypeParam)9L, (TypeParam)16L};

  _hyList <TypeParam> list (4, array),
                      list2 = list;
  
    
  unsigned long iteration_limiter = 1000000UL,
                current_iteration = 0UL,
                base_array_size = list.countitems();
                
  char     *has_replacement,
           with_replacement [] = " with replacement",
           without_replacement [] = " without replacement";
                
  
  for (long i = 0; i < 2; i++) {
  
    long array_size;
  
    if (i == 0) {
      has_replacement = with_replacement;
      array_size = base_array_size * base_array_size;
    } else {
      has_replacement = without_replacement;
      array_size = base_array_size * (base_array_size-1);
    }
    
    _hyList <unsigned long> counter;
  
    counter.append_multiple (0UL, array_size);
         
    while (counter.Find (0UL) != HY_NOT_FOUND && current_iteration < iteration_limiter) {
      _hyList <TypeParam> * subset = list.Subset (2, i == 0);
      
      ASSERT_EQ (2UL, subset->countitems()) << " Incorrect subset list length";
      unsigned long i0 = list.Find (subset->Element (0)),
                    i1 = list.Find (subset->Element (1));
                  
      ASSERT_NE (HY_NOT_FOUND, i0) << "Subset list element not found in the original list";
      ASSERT_NE (HY_NOT_FOUND, i1) << "Subset list element not found in the original list";
      
      if (i == 1) {
        ASSERT_NE (i0, i1) << "Sampling without replacement generated items with replacement";
      }
      
      i0 = base_array_size * i0 + i1;
      if (i == 1) {
        i0 -= 1 + (i0 > 5) + (i0 > 10);
      }
      //printf ("%ld :: %ld\n", i0, counter.Find (0UL));
      counter[i0] ++;
      current_iteration++;
      delete subset;
    }
    
    ASSERT_LT (current_iteration, iteration_limiter) << "Not all possible subsets have been generated " << has_replacement << ". Has not generated " << counter.Find (0UL);
  }
  
  
  
  for (unsigned long repl = 0UL; repl < 2UL; repl ++) {
    if (repl) {
      has_replacement = with_replacement;
    } else {
      has_replacement = without_replacement;
    }
    for (unsigned long block = 1UL; block <= 2UL; block++) {
      list2 = list;
      long order [4] = {0L, 0L, 0L, 0L};
      for (unsigned long i = 0L; i < 100; i++) {
        _hyList <long> counter;
        counter.append_multiple (HY_NOT_FOUND, base_array_size);
        
        if (repl) {
          list2.PermuteWithReplacement (block);
        } else {
          list2.Permute (block);
        }
        for (unsigned long k = 0UL; k < list2.countitems(); k++) {
          long idx = list.Find (list2.Element (k)); 
          ASSERT_NE (HY_NOT_FOUND, idx) << "A permuted list contains elements different from the original list, block size " << block << has_replacement << " index " << k; 
          counter[idx] = k;
        }
        
        if (repl == 0) {
          ASSERT_EQ (HY_NOT_FOUND, counter.Find (HY_NOT_FOUND)) << "A permuted list is missing some of the original elements, block size " << block << has_replacement;
        }
        
        if (block == 2UL) {
          for (unsigned long k = 0UL; k < counter.countitems(); k+=2) {
            if (repl && counter (k) == HY_NOT_FOUND) {
              continue;
            }
            ASSERT_EQ (counter(k) + 1, counter (k+1)) << "Failed to preserve within block-ordering in Permute (2) in block " << k << has_replacement;
          }
          if (repl) {
            order [(counter(0UL) > 0) * 2 + (counter (2UL) > 0)] ++;
          } else {
            order[counter(0UL) > 0] ++;
          }
        }
      
      }
      
      if (block > 1) {
        ASSERT_LE (0, order[0]) << "Missing 0 permutation order" << has_replacement;
        ASSERT_LE (0, order[1]) << "Missing 1 permutation order" << has_replacement;
        if (repl) {
          ASSERT_LE (0, order[2]) << "Missing 2 permutation order" << has_replacement;
          ASSERT_LE (0, order[3]) << "Missing 3 permutation order" << has_replacement;
        }
      }
    }
  }
  
  // N choose K testing
  
  _hyList <TypeParam> chosen_items;
  _hyList <long> NKtest;
  
  ASSERT_TRUE (list.NChooseKInit (&NKtest, chosen_items, 3)) << "Failed to initialize N choose K iteration";
  
  long expected_orders [] [3] = { {0,1,2} , {0,1,3}, { 0, 2, 3 }, {1,2,3}};
  
  for (long i = 0; i < 4; i++) {  
    if (i < 3) {
      ASSERT_TRUE (list.NChooseK (&NKtest, chosen_items)) << "Failed to get the set " << (i+1) << " in N choose K"; 
    } else {
      ASSERT_FALSE (list.NChooseK (&NKtest, chosen_items)) << "Should stop after set " << (i+1) << " in N choose K"; 
    }
    for (long j = 0; j < 3; j++) {
      EXPECT_EQ (expected_orders[i][j], list.Find (chosen_items (j))) << "Incorrect element " << (j) << " in set " << (i+1) << " in 4 choose 3";
    }
   }
  
}

REGISTER_TYPED_TEST_CASE_P (_hyListTest, ConstuctorTests, AccessAndManipulationTests, FindAddDelete, InPlaceOperations, SubsetsAndSelectors);

}

typedef ::testing::Types<char, long, double, _testPayload> _hyListTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(_typedList, _hyListTest, _hyListTestTypes);
