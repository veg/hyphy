/*
 *  ut_lists.cpp
 *  HyPhyXCode
 *
 *  Created by Steven Weaver on 6/17/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <tr1/tuple>
#include <iostream>
#include "gtest/gtest.h"
#include "ut_simplelists.h"

#include "hy_strings.h"
#include "hy_lists.h"

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
class _SimpleListTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _SimpleListTest() {
        // You can do set-up work for each test here.
    }

    virtual ~_SimpleListTest() {
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


TEST_F(_SimpleListTest, _PopulateTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    EXPECT_EQ(4,sl.lLength);
    EXPECT_EQ(7,sl[3]);
}

TEST_F(_SimpleListTest, _NormalizeCoordinatesTest){
// TODO: This function is questionable. 
// It doesn't seem to really normalize anything and it's not used
}

TEST_F(_SimpleListTest, _OffsetTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(7,sl[3]);
}

TEST_F(_SimpleListTest, _ElementTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(7,sl.Element(3));
}

TEST_F(_SimpleListTest, _PopTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(7,sl.Pop());
    EXPECT_EQ(3,sl.lLength);
}

TEST_F(_SimpleListTest, _countitemsTest){
    //TODO: Not camelcased
    _SimpleList sl; 
    sl.Populate(4,1,2);
    EXPECT_EQ(4,sl.countitems());
}

TEST_F(_SimpleListTest, _EqualTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList sl2; 
    sl.Populate(4,1,2);

    _SimpleList sl3; 
    sl.Populate(4,1,3);

    EXPECT_EQ(true,sl.Equal(sl2));
    EXPECT_EQ(false,sl.Equal(sl3));
}

TEST_F(_SimpleListTest, _MergeTest){
//TODO: This seems like it could be optimized
}


TEST_F(_SimpleListTest,AmpersandOpTest){}
TEST_F(_SimpleListTest,DoubleLessOpTest){}
TEST_F(_SimpleListTest,DoubleLess2OpTest){}
TEST_F(_SimpleListTest,DoubleGreaterOpTest){}


TEST_F(_SimpleListTest, _ListToPartitionStringTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _String* returned_string = (_String*)sl.ListToPartitionString();
    EXPECT_STREQ("1,3,5,7", returned_string->getStr());
}

TEST_F(_SimpleListTest, _RequestSpaceTest){
//TODO
}

TEST_F(_SimpleListTest, _toStrTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _String* returned_string = (_String*)sl.toStr();
    EXPECT_STREQ("{1,3,5,7}", returned_string->getStr());
}

TEST_F(_SimpleListTest, _makeDynamicTest){
//TODO
}

TEST_F(_SimpleListTest, _MinTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    long min = sl.Min();
    EXPECT_EQ(1, min);
}

TEST_F(_SimpleListTest, _MaxTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    long min = sl.Min();
    EXPECT_EQ(1, min);
}

TEST_F(_SimpleListTest, _ClearFormulasInListTest){
    _SimpleList sl;
    sl.Populate(4,1,2);

    long max = sl.Max();
    EXPECT_EQ(7, max);
}

TEST_F(_SimpleListTest, _DebugVarListTest){
    //TODO, We don't have to do a UnitTest for this
}

TEST_F(_SimpleListTest, _CountingSortTest){
    //TODO
    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList* returned_list;
    _SimpleList* ordering_list;

    //returned_list = sl.CountingSort(20, ordering_list);
    //long ret = returned_list->Element(1);
    EXPECT_EQ(1,1);
}

TEST_F(_SimpleListTest, _BinaryInsertTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    long pos = sl.BinaryInsert(4);
    
    EXPECT_EQ(2, pos);

}

TEST_F(_SimpleListTest, _FindTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    long pos = sl.Find(1,3);
    EXPECT_EQ( -1, pos);

    pos = sl.Find(3,0);
    EXPECT_EQ( 1, pos);

}

TEST_F(_SimpleListTest, _FindSteppingTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    long pos = sl.FindStepping(1,1,3);
    EXPECT_EQ( -1, pos);

    pos = sl.FindStepping(3,2,0);
    EXPECT_EQ( -1, pos);

    pos = sl.FindStepping(5,2,0);
    EXPECT_EQ( 2, pos);

}

TEST_F(_SimpleListTest, _FilterRangeTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    sl.FilterRange(2,4);
    EXPECT_EQ(3, sl[0]);
}

TEST_F(_SimpleListTest, _BinaryFindTest){

    _SimpleList sl; 
    sl.Populate(4,1,3);
    long pos = sl.BinaryFind(3,0);
    EXPECT_EQ( -3, pos);

}

TEST_F(_SimpleListTest, _SortTest){
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
}


TEST_F(_SimpleListTest, _CompareTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);

    long resp = sl.Compare((long)0,(long)1);
    EXPECT_EQ( -1, resp);

    resp = sl.Compare((long)1,(long)0);
    EXPECT_EQ( 1, resp);

    resp = sl.Compare((long)1,(long)1);
    EXPECT_EQ( 0, resp);
}

TEST_F(_SimpleListTest, _ClearTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl.Clear(true);
    EXPECT_EQ(0, sl.lLength);
}

TEST_F(_SimpleListTest, _DeleteTest){
    _SimpleList sl;
    sl.Populate(4,1,2);
    sl.Delete(0,true);
    EXPECT_EQ(3, sl[0]);
}

//TEST_F(_SimpleListTest, _TrimMemoryTest){
////TODO

//}

TEST_F(_SimpleListTest, _DeleteDuplicatesTest){
    _SimpleList sl; 
    sl.Populate(4,1,2);
    sl << 5;
    sl.Sort(true); 
    sl.DeleteDuplicates();
    

    EXPECT_EQ(5, sl[2]);
    EXPECT_EQ(7, sl[3]);
}

TEST_F(_SimpleListTest, _DeleteListTest){

    _SimpleList sl; 
    sl.Populate(4,1,2);

    _SimpleList to_delete; 
    to_delete.Populate(3,0,1);

    sl.DeleteList(to_delete);
    EXPECT_EQ(7, sl[0]);
}

TEST_F(_SimpleListTest, _DisplaceTest){

    _SimpleList sl; 
    sl.Populate(10,1,1);
    sl.Displace(0,5,1); 

    EXPECT_EQ(7, sl[0]);
    EXPECT_EQ(5, sl[5]);

}

TEST_F(_SimpleListTest, _PermuteWithReplacementTest){

    //TODO
    _SimpleList sl; 
    sl.Populate(10,1,1);
    sl.PermuteWithReplacement(2);
    EXPECT_EQ(1, 1);

}

TEST_F(_SimpleListTest, _PermuteTest){

    //TODO
    _SimpleList sl; 
    sl.Populate(10,1,1);

    sl.Permute(2);
    EXPECT_EQ(1, 1);

}

TEST_F(_SimpleListTest, _NChooseKTest){

    _SimpleList sl; 
    sl.Populate(10,1,1);

    _SimpleList state; 
    _SimpleList store; 
    long stride = 3;
    bool algorithm = true;

    bool init_test = sl.NChooseKInit(state,store,stride,algorithm);
    bool choose_test = sl.NChooseK(state,store);

    EXPECT_EQ(true, init_test);
    EXPECT_EQ(true, choose_test);

    //This doesn't seem very tuple-ly
    EXPECT_EQ(1,store[0]);
    EXPECT_EQ(1,store[0]);

}

TEST_F(_SimpleListTest, _SwapTest){

    _SimpleList sl; 
    sl.Populate(10,1,1);
    sl.Swap(0,1);

    EXPECT_EQ(2,sl[0]);
    EXPECT_EQ(1,sl[1]);
}

TEST_F(_SimpleListTest, _FlipTest){

    _SimpleList sl; 
    sl.Populate(10,1,1);
    sl.Flip();

    EXPECT_EQ(10,sl[0]);
    EXPECT_EQ(1,sl[9]);

}

//TEST_F(_SimpleListTest, _RecursiveIndexSortTest){
////TODO

//}

TEST_F(_SimpleListTest, _UnionTest){

    _SimpleList sl; 

    _SimpleList l1; 
    l1.Populate(10,1,1);
    _SimpleList l2; 
    l2.Populate(20,1,1);

    sl.Union(l1, l2);
    EXPECT_EQ(20,sl[19]);

}

TEST_F(_SimpleListTest, _IntersectTest){

    _SimpleList sl; 

    _SimpleList l1; 
    l1.Populate(10,1,1);
    _SimpleList l2; 
    l2.Populate(10,1,2);

    sl.Intersect(l1, l2);
    EXPECT_EQ(9,sl[4]);

}

TEST_F(_SimpleListTest, _CountCommonElementsTest){

    _SimpleList sl; 

    _SimpleList l1; 
    l1.Populate(10,1,1);
    _SimpleList l2; 
    l2.Populate(10,1,2);

    sl.Union(l1, l2);
    EXPECT_EQ(5,l1.CountCommonElements(l2,false));

}

TEST_F(_SimpleListTest, _XORTest){

    _SimpleList sl; 

    _SimpleList l1; 
    l1.Populate(10,1,1);
    _SimpleList l2; 
    l2.Populate(10,1,2);

    sl.XOR(l1, l2);
    EXPECT_EQ(10,sl[4]);

}

TEST_F(_SimpleListTest, _SubtractTest){

    _SimpleList sl; 

    _SimpleList l1; 
    l1.Populate(10,1,1);
    _SimpleList l2; 
    l2.Populate(10,1,2);

    sl.Subtract(l1, l2);
    EXPECT_EQ(10,sl[9]);

}
}
