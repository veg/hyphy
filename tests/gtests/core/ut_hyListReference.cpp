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
#include "hy_list_reference.h"

namespace {

class _testPayload: public BaseObj {

  public:
  
  _testPayload (void) { data = 0UL;}
  _testPayload (const char* p) { data = atoi (p); }
  _testPayload (const long p) { data = p; }
  _testPayload (const _testPayload& o) { data = o.data; }
  virtual BaseObj *makeDynamic(void) const { return new _testPayload (*this); }
  virtual void Duplicate(BaseObj const * ref) { data = ((_testPayload*)ref)->data; }
  virtual bool Equal (const _testPayload * o) {return data == o->data;}

  long data;
};

// The fixture for testing class Foo.
template <typename DATA>
class _hyListReferenceTest : public ::testing::Test {

protected:

  _hyListReferenceTest() {
    
  }

  virtual ~_hyListReferenceTest() {
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


TYPED_TEST_CASE_P(_hyListReferenceTest);


TYPED_TEST_P (_hyListReferenceTest, ConstuctorTests) {
  
   TypeParam array [5] = {(TypeParam)"1", (TypeParam)"10", (TypeParam)"2", (TypeParam)"3", (TypeParam)"7"},
             *array2 [5] = {array, array + 1L, array + 2L, array + 3L, array + 4L},
             testValue = (TypeParam)"3",
             * dynamic_object = new TypeParam ("3.14");
             
  {
  
    TypeParam         *another_dynamic_object = new TypeParam ("2.7172");
    
    _hyListReference <TypeParam> null_list,
                                 single_element_list (testValue),
                                 multiple_element_list (5,array2),
                                 full_stack_copy (multiple_element_list),
                                 partial_stack_copy (multiple_element_list,2,HY_LIST_INSERT_AT_END),
                                 * dynamic_list = new _hyListReference <TypeParam> (dynamic_object),
                                 * another_dynamic_list = new _hyListReference <TypeParam> (another_dynamic_object);
            
    ASSERT_EQ (0UL, null_list.countitems()) << "Non-empty null list";
    ASSERT_EQ (1UL, single_element_list.countitems()) << "Single element list has wrong length";
    ASSERT_EQ (5UL, multiple_element_list.countitems()) << "Array of elements list has wrong length";
    ASSERT_EQ (5UL, full_stack_copy.countitems()) << "Stack copy list has wrong length";
    ASSERT_EQ (3UL, partial_stack_copy.countitems()) << "Partial stack copy list has wrong length";
    
    EXPECT_TRUE (single_element_list (0)->Equal (full_stack_copy(3)));   
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
    EXPECT_EQ (dynamic_object, dynamic_list->Element(0));
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
  //EXPECT_TRUE (dynamic_object->SingleReference()) << "Incorrect reference count following list operations _list (object *)";
  
  //DeleteObject(dynamic_object);
}

TYPED_TEST_P (_hyListReferenceTest, ListManipulationTests) {
  _hyListReference <TypeParam> list1, 
                               list2,
                               list3,
                               list4;
                               
   for (long k = 0; k < 50; k ++) {
      TypeParam * test_object = new TypeParam (k),
                * test_object2;
                
      list1.AppendNewInstance (test_object);
      EXPECT_TRUE (test_object->SingleReference()) << "Append New Instance does not maintain reference counts";
      list2 << test_object; 
      EXPECT_FALSE (test_object->SingleReference()) << "operator << does not increment reference counts";
      list3 && test_object; 
      test_object2 = list3.Element (-1);
      EXPECT_NE (list2.Element (-1), test_object2) << "operator && does not make a dynamic object copy";
      EXPECT_TRUE (test_object2->SingleReference()) << "The object added by && does not have a single reference count";
   }
   
   list4 = list1;
   list4 << list2;
   
   list3.Clear();
   EXPECT_EQ (list4.countitems(), list1.countitems() + list2.countitems()) << "Incorrect list length obtained by <<";
   
   _hyListReference <TypeParam> list5 = list1 & list4;
   EXPECT_EQ (list5.countitems(), list1.countitems() + list4.countitems()) << "Incorrect list length obtained by list & list";
   
   TypeParam an_object (64L);
   EXPECT_EQ (list5.Find (&an_object), HY_NOT_FOUND);
   list5 = list5 & &an_object;
   
   EXPECT_EQ (list5.countitems(), list1.countitems() + list4.countitems() + 1) << "Incorrect list length obtained by list & item";
   EXPECT_EQ (list5.Find (&an_object), list5.countitems() - 1);
   
   
   
}

REGISTER_TYPED_TEST_CASE_P (_hyListReferenceTest, ConstuctorTests, ListManipulationTests);
  
}

typedef ::testing::Types<_testPayload, _String> _hyListReferenceTestTypes;
INSTANTIATE_TYPED_TEST_CASE_P(_typedList, _hyListReferenceTest, _hyListReferenceTestTypes);
