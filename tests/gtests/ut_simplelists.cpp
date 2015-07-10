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

//#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include "ut_simplelists.h"

#include "hy_strings.h"
#include "simplelist.h"

#include "baseobj.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using ::testing::TestWithParam;
using ::testing::Values;
using ::testing::Range;
//using ::testing::Combine;


namespace
{

// The fixture for testing class Foo.
class SimpleListTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    SimpleListTest() {
        // You can do set-up work for each test here.
    }

    virtual ~SimpleListTest() {
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

};

TEST_F(SimpleListTest,_StackCopyConstructorListTest){
    _SimpleList sl;
    sl.Populate(4,1,2);

    _SimpleList sl2(sl,(long)-1,(long)-1); 
    EXPECT_EQ(0,sl2.lLength);

    _SimpleList sl3(sl,(long)1,(long)3); 
    EXPECT_EQ(3,sl3[0]);
}

TEST_F(SimpleListTest,_LengthConstructorTest){
    _SimpleList sl((long)7);
    EXPECT_EQ(8,sizeof(sl.lData));
}

TEST_F(SimpleListTest, _PopulateTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    EXPECT_EQ(4,sl.lLength);
    EXPECT_EQ(7,sl[3]);
}

//TEST_F(SimpleListTest, _NormalizeCoordinatesTest){
//}

TEST_F(SimpleListTest, _BracketOpTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    EXPECT_EQ(1,sl[0]);
    EXPECT_EQ(7,sl[3]);

    sl.lLength = 0;
    EXPECT_EQ(1,sl[7]);

}

TEST_F(SimpleListTest, _ParenthOpTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(7,sl(3));

    //TODO: We need a way to test for warnErrors
    //EXPECT_EQ(0,sl(10));
}

TEST_F(SimpleListTest, _EqualOpTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList sl2((long)0); 
    sl2 = sl;
    EXPECT_EQ(3,sl[1]);
}

TEST_F(SimpleListTest, _DoubleLessOpTest){

    _SimpleList sl;
    sl.Populate(4,1,2);

    _SimpleList sl2;
    sl2.Populate(4,1,2);

    sl << sl2;
    EXPECT_EQ(3,sl[5]);

}

TEST_F(SimpleListTest, _OffsetTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl.Offset(3);
    EXPECT_EQ(10,sl[3]);
}

TEST_F(SimpleListTest, _ElementTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(7,sl.Element(3));
    EXPECT_EQ(3,sl.Element(-3));
    EXPECT_EQ(0,sl.Element(-5));
}

TEST_F(SimpleListTest, _PopTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(7,sl.Pop());
    EXPECT_EQ(3,sl.lLength);

    //Returns a 0
    _SimpleList sl2;
    EXPECT_EQ(0,sl2.Pop());
}

TEST_F(SimpleListTest, _countitemsTest){
    //TODO: Not camelcased
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(4,sl.countitems());
}

TEST_F(SimpleListTest, _EqualTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList sl2; 
    sl2.Populate(4,1,2);

    _SimpleList sl3;
    sl3.Populate(4,1,3);

    _SimpleList sl4;
    sl4.Populate(9,1,3);

    _SimpleList sl5;
    sl5.Populate(4,15,3);


    EXPECT_EQ(true,sl.Equal(sl2));
    EXPECT_EQ(false,sl.Equal(sl3));
    EXPECT_EQ(false,sl.Equal(sl4));
    EXPECT_EQ(false,sl.Equal(sl5));
}

TEST_F(SimpleListTest, _MergeTest){
    //TODO: Coverage Testing

    //Takes 4 parameters
    _SimpleList sl; 

    _SimpleList l1;
    l1.Populate(4,1,1);

    _SimpleList l2;
    l2.Populate(4,5,1);

    _SimpleList l3;
    l3.Populate(12,1,1);

    //List that are automatically going to be filled
    _SimpleList m1; 
    _SimpleList m2; 

    sl.Merge(l1,l2,&m1,&m2);
    EXPECT_EQ(8,sl[8]);

    sl.Merge(l2,l1,&m1,&m2);
    EXPECT_EQ(8,sl[8]);

    sl.Merge(l2,l3,&m1,&m2);
    EXPECT_EQ(9,sl[8]);

    sl.Merge(l1,l1,&m1,&m2);
    EXPECT_EQ(4,sl[8]);

}


//TEST_F(SimpleListTest,AmpersandOpTest){
//    //TODO: Assignment doesn't work approprately
//    _SimpleList sl((long)4,(long)1,(long)2); 
//    _SimpleList sl2((long)4,(long)11,(long)2); 
//    _SimpleList sl3(sl & sl2);
//
//    EXPECT_EQ(11,sl[5]);
//
//}

TEST_F(SimpleListTest,DoubleGreaterOpTest){
    //Does the same as lesser than, but no dupes and returns bool
    _SimpleList sl; 
    bool r1, r2;

    sl.Populate(4,1,2);

    r1 = sl >> 1; 
    r2 = sl >> 12; 
    

    EXPECT_FALSE(r1);
    EXPECT_TRUE(r2);
    EXPECT_EQ(12,sl[5]);
}


TEST_F(SimpleListTest, _ListToPartitionStringTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _String* returned_string = (_String*)sl.ListToPartitionString();
    EXPECT_STREQ("1,3,5,7", returned_string->getStr());

    _SimpleList sl2;
    sl2.Populate(4,1,1);

    _String* returned_string2 = (_String*)sl2.ListToPartitionString();
    EXPECT_STREQ("1-4", returned_string2->getStr());

    _SimpleList sl3;
    sl3.Populate(4,1,1);
    sl3 << (long)7;
    sl3 << (long)9;
    sl3 << (long)10;
    sl3 << (long)11;

    _String* returned_string3 = (_String*)sl3.ListToPartitionString();
    EXPECT_STREQ("1-4,7,9-11", returned_string3->getStr());
}

TEST_F(SimpleListTest, _RequestSpaceTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl.RequestSpace(10);
    EXPECT_EQ(1,sl[0]);
}

TEST_F(SimpleListTest, _toStrTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    _String* returned_string = (_String*)sl.toStr();
    EXPECT_STREQ("{1,3,5,7}", returned_string->getStr());

    _SimpleList sl2;
    sl2.Populate(100,1,1);
    _String* returned_string2 = (_String*)sl2.toStr();
    EXPECT_STREQ("{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100}", returned_string2->getStr());

    _SimpleList sl3;
    _String* returned_string3 = (_String*)sl3.toStr();
    EXPECT_STREQ("{}", returned_string3->getStr());
}

TEST_F(SimpleListTest, _makeDynamicTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList* dynamic_list = (_SimpleList*)sl.makeDynamic();
    EXPECT_EQ(4,dynamic_list->lLength);

}

TEST_F(SimpleListTest, _MinTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    long min = sl.Min();
    EXPECT_EQ(1, min);
}

TEST_F(SimpleListTest, _MaxTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    long min = sl.Min();
    EXPECT_EQ(1, min);
}

TEST_F(SimpleListTest, _ClearFormulasInListTest){
    //TODO: Debug Coverage
    _SimpleList sl;
    sl.Populate(4,1,2);

    long max = sl.Max();
    EXPECT_EQ(7, max);
}

TEST_F(SimpleListTest, _DebugVarListTest){
    //We don't have to do a UnitTest for this
}

TEST_F(SimpleListTest, _CountingSortTest){
    _SimpleList sl; 
    _SimpleList ordering_list;
    sl.Populate(4,1,2);

    _SimpleList* returned_list;

    returned_list = sl.CountingSort(20, &ordering_list);
    long ret = returned_list->Element(1);
    EXPECT_EQ(1,1);
}

TEST_F(SimpleListTest, _BinaryInsertTest){

    _SimpleList sl,sl2;
    long pos;

    sl.Populate(4,1,2);
    pos = sl.BinaryInsert(4);
    EXPECT_EQ(2, pos);

    pos = sl2.BinaryInsert(4);
    EXPECT_EQ(0, pos);

    pos = sl.BinaryInsert(1);
    EXPECT_EQ(4, pos);
}

TEST_F(SimpleListTest, _FindTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    long pos = sl.Find(1,3);
    EXPECT_EQ( -1, pos);

    pos = sl.Find(3,0);
    EXPECT_EQ( 1, pos);

}

TEST_F(SimpleListTest, _FindSteppingTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    long pos = sl.FindStepping(1,1,3);
    EXPECT_EQ( -1, pos);

    pos = sl.FindStepping(3,2,0);
    EXPECT_EQ( -1, pos);

    pos = sl.FindStepping(5,2,0);
    EXPECT_EQ( 2, pos);

}

TEST_F(SimpleListTest, _FilterRangeTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    sl.FilterRange(2,4);
    EXPECT_EQ(3, sl[0]);

    sl.FilterRange(4,2);
    EXPECT_EQ(0, sl.lLength);
}

TEST_F(SimpleListTest, _BinaryFindTest){

    _SimpleList sl; 
    long pos;

    sl.Populate(4,1,3);
    pos = sl.BinaryFind(3,0);
    EXPECT_EQ( -3, pos);

}

TEST_F(SimpleListTest, _SortTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl.Sort(false);
    EXPECT_EQ( 7, sl[0]);

    sl.Sort(true);
    EXPECT_EQ( 1, sl[0]);

    _SimpleList sl2; 
    sl2.Populate(20,1,2);
    sl2.Sort(false);
    EXPECT_EQ( 39, sl2[0]);

    //Create a worst case scenario for Quick Sort (Already sorted)
    sl2.Sort(false);
    EXPECT_EQ( 39, sl2[0]);
}


TEST_F(SimpleListTest, _CompareTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    long resp = sl.Compare((long)0,(long)1);
    EXPECT_EQ( -1, resp);

    resp = sl.Compare((long)1,(long)0);
    EXPECT_EQ( 1, resp);

    resp = sl.Compare((long)1,(long)1);
    EXPECT_EQ( 0, resp);
}

TEST_F(SimpleListTest, _ClearTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl.Clear(true);
    EXPECT_EQ(0, sl.lLength);
}

TEST_F(SimpleListTest, _DeleteTest){
    //TODO: Debug Coverage
    _SimpleList sl;
    sl.Populate(4,1,2);
    sl.Delete(0,true);
    EXPECT_EQ(3, sl[0]);
}

TEST_F(SimpleListTest, _TrimMemoryTest){

    _SimpleList sl;
    sl.Populate(4,1,2);
    sl.Delete(0,true);
    EXPECT_EQ(3, sl[0]);

}

TEST_F(SimpleListTest, _DeleteDuplicatesTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl << 5;
    sl.Sort(true); 
    sl.DeleteDuplicates();
    

    EXPECT_EQ(5, sl[2]);
    EXPECT_EQ(7, sl[3]);
}

TEST_F(SimpleListTest, _DeleteListTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList to_delete; 
    to_delete.Populate(3,0,1);

    sl.DeleteList(to_delete);
    EXPECT_EQ(7, sl[0]);
}

TEST_F(SimpleListTest, _DisplaceTest){
    _SimpleList sl, sl2, sl3;
    sl.Populate(10,1,1);
    sl2.Populate(10,1,1);
    sl3.Populate(10,1,1);

    //TODO: Doesn't work as expected
    sl.Displace(1,5,1); 
    EXPECT_EQ(7, sl[1]);
    EXPECT_EQ(5, sl[5]);

    //Don't do anything
    sl2.Displace(11,10,1); 
    EXPECT_EQ(1, sl2[0]);

    //This shouldn't do anything
    sl3.Displace(-1,25,0);
    EXPECT_EQ(1, sl3[0]);

}

TEST_F(SimpleListTest, _PermuteWithReplacementTest){

    //TODO: Go through with debugger
    _SimpleList sl, sl2;
    sl.Populate(10,1,1);
    sl2.Populate(10,1,1);

    sl.PermuteWithReplacement(2);

    //Cannot be zero
    //sl2.PermuteWithReplacement(0);
    sl2.PermuteWithReplacement(-1);

    EXPECT_EQ(1, 1);

}

TEST_F(SimpleListTest, _PermuteTest){

    //TODO: Go through with debugger
    _SimpleList sl, sl2; 
    sl.Populate(10,1,1);
    sl2.Populate(10,1,1);

    sl.Permute(2);

    //Cannot be zero
    //sl2.Permute(0);

    //Cannot be negative
    //sl2.Permute(-1);

    EXPECT_EQ(1, 1);

}

TEST_F(SimpleListTest, _NChooseKTest){

    _SimpleList sl, sl2;
    _SimpleList state, state2;
    _SimpleList store, store2;

    long stride = 3;
    long large_stride = 42;

    bool algorithm = true;

    sl.Populate(10,1,1);

    bool init_test = sl.NChooseKInit(state,store,large_stride,algorithm);
    EXPECT_EQ(false, init_test);
    init_test = sl.NChooseKInit(state,store,stride,algorithm);
    bool choose_test = sl.NChooseK(state,store);


    EXPECT_EQ(true, init_test);
    EXPECT_EQ(true, choose_test);

    //TODO:This doesn't seem very tuple-ly
    EXPECT_EQ(1,store[0]);
    EXPECT_EQ(1,store[0]);

}

TEST_F(SimpleListTest, _SwapTest){

    _SimpleList sl; 
    sl.Populate(10,1,1);
    sl.Swap(0,1);

    EXPECT_EQ(2,sl[0]);
    EXPECT_EQ(1,sl[1]);
}

TEST_F(SimpleListTest, _FlipTest){

    _SimpleList sl; 
    sl.Populate(10,1,1);
    sl.Flip();

    EXPECT_EQ(10,sl[0]);
    EXPECT_EQ(1,sl[9]);

}

//TEST_F(SimpleListTest, _RecursiveIndexSortTest){
////TODO

//}

TEST_F(SimpleListTest, _UnionTest){

    _SimpleList sl, sl2; 
    _SimpleList l1, l2; 

    sl.Populate(10,1,1);
    l1.Populate(10,1,1);
    l2.Populate(20,1,1);

    sl.Union(l1, l2);
    EXPECT_EQ(20,sl[19]);

    //For code coverage
    sl2.Union(l2, l1);
    EXPECT_EQ(20,sl[19]);

}

TEST_F(SimpleListTest, _IntersectTest){

    _SimpleList sl, sl2; 
    _SimpleList l1, l2; 

    sl.Populate(10,1,1);
    l1.Populate(10,1,1);
    l2.Populate(10,1,2);

    sl.Intersect(l1, l2);
    EXPECT_EQ(9,sl[4]);

    //For code coverage
    sl2.Intersect(l2, l1);
    EXPECT_EQ(9,sl2[4]);
}

TEST_F(SimpleListTest, _CountCommonElementsTest){

    _SimpleList sl; 

    _SimpleList l1; 
    l1.Populate(10,1,1);
    _SimpleList l2; 
    l2.Populate(10,1,2);

    sl.Union(l1, l2);
    EXPECT_EQ(5,l1.CountCommonElements(l2,false));

}

TEST_F(SimpleListTest, _XORTest){

    _SimpleList sl, sl2; 

    _SimpleList l1, l2; 

    sl.Populate(10,1,1);
    l1.Populate(10,1,1);
    l2.Populate(10,1,2);

    sl.XOR(l1, l2);
    EXPECT_EQ(10,sl[4]);

    //For coverage
    sl2.XOR(l2, l1);
    EXPECT_EQ(10,sl2[4]);
}

TEST_F(SimpleListTest, _SubtractTest){

    _SimpleList sl, sl2; 

    _SimpleList l1, l2; 


    sl.Populate(10,1,1);
    l1.Populate(10,1,1);
    l2.Populate(10,1,2);

    sl.Subtract(l1, l2);
    EXPECT_EQ(10,sl[9]);

    sl2.Subtract(l1, l2);
    EXPECT_EQ(10,sl2[9]);


}
}
