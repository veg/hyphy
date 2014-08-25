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
#include "hy_avllist.h"
#include "hy_list_reference.h"

namespace {
  
  class _testPayload {
  public:
    
    _testPayload (void) { data = 0UL; unused = 0UL;}
    _testPayload (unsigned long p) { data = p; unused = 0UL; }
    _testPayload (const _testPayload& o) { data = o.data; unused = o.unused; }
    
    bool operator == (const _testPayload & o) const { return data == o.data;}
    
    operator unsigned long (void) {return data;}
    
    unsigned long data, unused;
  };
  
  // The fixture for testing class Foo.
  template <typename DATA>
  class _hyAVLListTest : public ::testing::Test {
    
  protected:
    
    _hyAVLListTest() {
      
    }
    
    virtual ~_hyAVLListTest() {
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
    
    _StringBuffer dump_to_stream_as_longs (const _hyList <DATA>& data) {
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
  
  TYPED_TEST_CASE_P(_hyAVLListTest);
  
  
  TYPED_TEST_P (_hyAVLListTest, ConstuctorTests) {
    TypeParam array [5] = {(TypeParam)1, (TypeParam)4, (TypeParam)9, (TypeParam)16, (TypeParam)25},
    array2 [5] = {(TypeParam)1, (TypeParam)0, (TypeParam)9, (TypeParam)16, (TypeParam)25};
    
    _hyList <TypeParam> multiple_element_defs (5,array),
    multiple_element_defs2 (5, array2);
    
    
    _AVLList <TypeParam>  test_list,
    single_item_list ((TypeParam)42),
    multiple_elements_list (multiple_element_defs),
    multiple_elements_list2 (multiple_element_defs2),
    
    stack_copy (multiple_elements_list);
    
    ASSERT_EQ (0UL, test_list.Length()) << "Default constructor must create an empty _AVLList";
    ASSERT_EQ (1UL, single_item_list.Length()) << "Default single item constructor must create a single item _AVLList";
    ASSERT_EQ (multiple_element_defs.Length(), multiple_elements_list.Length()) << "Default list of items constructor must create an _AVLList of the same length as the argument list";
    
    ASSERT_EQ (multiple_elements_list,stack_copy) << " Failed _AVLList A(B) => A==B";
    ASSERT_NE (multiple_elements_list,multiple_elements_list2) << " Failed A != B";
    
    
  }
  
  
  TYPED_TEST_P (_hyAVLListTest, TraversalTests) {
    long up_to = 1L << 16;
    
    _AVLList <TypeParam> arithmetic_series;
    
    for (long index = -up_to; index < up_to; index += 4) {
      ASSERT_LE(0L,arithmetic_series.Insert ((TypeParam)index)) << "Incorrect return value for insert (missing element reported as found";
      ASSERT_LE(0L,arithmetic_series.Insert ((TypeParam)(0x7FFFFF - index))) << "Incorrect return value for insert (missing element reported as found";
    }
    
    
    for (long i = 0; i < 2; i++) {
      _AVLListTraversalHistory history;
      long k = arithmetic_series.Traverser (history, i);
      TypeParam last = * arithmetic_series.AtIndex (k),
      current;
      while (1) {
        k = arithmetic_series.Traverser (history, i);
        if (k != HY_NOT_FOUND) {
          current = * arithmetic_series.AtIndex (k);
          if (i) {
            ASSERT_LT (current, last) << "Keys returned by AVL traversal in reverse are not sorted";
          } else {
            ASSERT_GT (current, last) << "Keys returned by AVL traversal are not sorted";
          }
          last = current;
        } else {
          break;
        }
      }
    }
    
    
  }
  
  TYPED_TEST_P (_hyAVLListTest, InsertTests) {
    
    for (long replicates = 0L; replicates < 128L; replicates++) {
      _AVLList <TypeParam> test_list;
      long up_to = genrand_int32() % 10000;
      for (long item = 0; item < up_to; item ++) {
        test_list.Insert ((TypeParam) genrand_int32());
        
      }
      const char * err = test_list.ConsistencyCheck();
      ASSERT_EQ (NULL, err) << "Failed an internal avl consistency check: " << err;
    }
    
    
    {
      _AVLList <TypeParam> long_direct,
      long_reverse;
      
      long upper_bound = 1L << 20;
      for (long k = 0; k < upper_bound; k++) {
        long_direct.Insert ((TypeParam)  (k));
        long_reverse.Insert ((TypeParam)  (upper_bound-1-k));
      }
      //long_direct.EchoList();
      //long_reverse.EchoList();
      
      ASSERT_EQ (long_direct,long_reverse) << " Failed _AVLList (A) == _AVLList (reverse[A])";
      
    }
  }
  
  TYPED_TEST_P (_hyAVLListTest, DeleteTests) {
    
    
    _AVLList <TypeParam> sequence;
    
    long upto = 1L << 3;
    
    for (long k = 0L; k < upto; k++) {
      sequence.Insert ((TypeParam)k);
    }
    
    _AVLList <TypeParam> sequence_copy (sequence);
    
    //sequence.EchoList();
    for (long k = 0L; k < upto; k++) {
      ASSERT_LE (0L,sequence.Delete ((TypeParam)k)) << "Missing key " << k << " during delete";
      const char * err = sequence.ConsistencyCheck();
      ASSERT_EQ (NULL, err) << "Failed an internal avl consistency check: '" << err << "' after deleting key " << k;
      //printf ("Deleted %ld\n", k);
    }
    ASSERT_EQ (0UL, sequence.Length()) << "AVL Length must be zero after all elements have been deleted.";
    
    sequence = sequence_copy;
    //sequence_copy.EchoList();
    for (long k = upto - 1; k >= 0; k--) {
      ASSERT_LE (0L,sequence.Delete ((TypeParam)k)) << "Missing key " << k << " during reverse order delete";
      const char * err = sequence.ConsistencyCheck();
      ASSERT_EQ (NULL, err) << "Failed an internal avl consistency check: '" << err << "' after deleting (reverse order) key " << k;
    }
    
 
    for (long replicates = 0L; replicates < 128L; replicates++) {
      _AVLList <TypeParam> test_list;
      _hyList <TypeParam>  contents;
      long up_to = genrand_int32() % 10000;
      
      _SimpleList random_vector (up_to, 0, 1);
      random_vector.Permute (1);

      for (long item = 0; item < up_to; item ++) {
        TypeParam value = (TypeParam) genrand_int32();
        test_list.Insert (value);
        contents << value;
      }
        
      for (long item = 0; item < up_to/2; item ++) {
        ASSERT_LE (0L,test_list.Delete (contents(random_vector(item)))) << "Missing key " << item << " during randomized delete";
      }
      const char * err = test_list.ConsistencyCheck();
      ASSERT_EQ (NULL, err) << "Failed an internal avl consistency check: '" << err << "' after deleting (random order) half of the keys ";
    }
    

  }
  
  TYPED_TEST_P (_hyAVLListTest, CombinedTests) {
    
    
    _AVLList <TypeParam> sequence;
    for (long replicates = 0L; replicates < 128L; replicates++) {
      _AVLList <TypeParam> small_values,
                           large_values,
                           all_values;
      
      long up_to = genrand_int32() % 512;
      for (long item = 0; item < up_to; item ++) {
        TypeParam low  = (TypeParam) genrand_int32() % 0xFFFF,
                  high = 0xFFFFF + (TypeParam) genrand_int32() % 0xFFFFF;
        
        small_values.Insert (low);
        large_values.Insert (high);
        all_values.Insert (low);
        all_values.Insert (high);
      }
      
      _AVLListTraversalHistory hist;
      _AVLList <TypeParam> all_copy (all_values);
      
      long k = small_values.Traverser (hist);
      while (k != HY_AVL_LEAF) {
        ASSERT_LE (0L, all_values.Delete (*small_values.AtIndex (k))) << "Failed to find an existing element of the AVLList";
        ASSERT_EQ (HY_NOT_FOUND, large_values.Find (*small_values.AtIndex (k))) << "Found a missing element of the AVLList";
        ASSERT_EQ (HY_NOT_FOUND, all_values.Find (*small_values.AtIndex (k))) << "Found an element that should have been deleted";
        k = small_values.Traverser (hist);
      }
      
      
      k = small_values.Traverser (hist);
      while (k != HY_AVL_LEAF) {
        ASSERT_LE (0L, all_values.Insert (*small_values.AtIndex (k))) << "Found an element that should have been deleted";
        //all_values.Insert (*small_values.AtIndex (k));
        k = small_values.Traverser (hist);
     }
     
      //all_values.EchoList();
      //all_copy.EchoList();
      
     ASSERT_EQ (all_values, all_copy) << "AB, delete A, add A, != AB";

      
      
    }
  }
  
  
  
  REGISTER_TYPED_TEST_CASE_P (_hyAVLListTest, CombinedTests, ConstuctorTests, TraversalTests, InsertTests, DeleteTests);
  
}

typedef ::testing::Types<long> _hyAVLListTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(_typedList, _hyAVLListTest, _hyAVLListTestTypes);
