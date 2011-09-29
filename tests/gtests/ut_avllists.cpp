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
#include "ut_avllists.h"
#include "hy_strings.h"
#include "baseobj.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


_List _AVcreateStrList() {

    _List str_list;

    str_list.AppendNewInstance(new _String("zero"));
    str_list.AppendNewInstance(new _String("one"));
    str_list.AppendNewInstance(new _String("two"));
    str_list.AppendNewInstance(new _String("three"));
    str_list.AppendNewInstance(new _String("four"));
    str_list.AppendNewInstance(new _String("five"));
    str_list.AppendNewInstance(new _String("six"));

    return str_list;
}

namespace
{

// The fixture for testing class Foo.
class _AVLListTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _AVLListTest() {
        // You can do set-up work for each test here.
    }

    virtual ~_AVLListTest() {
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

class _AVLListXTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _AVLListXTest() {
        // You can do set-up work for each test here.
    }

    virtual ~_AVLListXTest() {
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

class _AVLListXLTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _AVLListXLTest() {
        // You can do set-up work for each test here.
    }

    virtual ~_AVLListXLTest() {
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

TEST_F(_AVLListTest,FindTest){
    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 

    //sl.Populate(4,1,2);
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.Find((BaseRef)3);
    EXPECT_EQ(3,info);
 }

TEST_F(_AVLListTest,FindLongTest){
    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 

    //sl.Populate(4,1,2);
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.FindLong(3);
    EXPECT_EQ(3,info);
}

/*
 *TEST_F(_AVLListTest,FindBestTest){
 *    //TODO. It seems like it would only keep going 
 *    //until we have a leftChild or 0
 *
 *}
 *
 *TEST_F(_AVLListTest,FindTest){
 *
 *    //AVL List always takes a pointer to a SimpleList
 *    long index; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *    index = al.Find(5);
 *    
 *    EXPECT_EQ(2,index);
 *
 *}
 *
 *TEST_F(_AVLListTest,NextTest){
 *    //TODO
 *
 *
 *}
 *
 *TEST_F(_AVLListTest,FirstTest){
 *    //TODO
 *    //AVL List always takes a pointer to a SimpleList
 *    //Keeps going down through left children until lData is 0
 *
 *    long index; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    index = al.First();
 *    EXPECT_EQ(1,index);
 *}
 *
 *TEST_F(_AVLListTest,LastTest){
 *    //TODO
 *    //AVL List always takes a pointer to a SimpleList
 *    //Keeps going down through right children until lData is 0
 *
 *    long index; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    index = al.Last();
 *    EXPECT_EQ(4,index);
 *}
 *
 *TEST_F(_AVLListTest,GetByIndexTest){
 *    //AVL List always takes a pointer to a SimpleList
 *    //Keeps going down through right children until lData is 0
 *    long info; 
 *    _SimpleList sl; 
 *
 *    //sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    for(int i=0; i<=10; i++) {
 *        al.Insert((BaseRef)i, 0, true, false);    
 *    }
 *    
 *    info = al.GetByIndex(3);
 *    EXPECT_EQ(7,info);
 * 
 *}
 *
 *TEST_F(_AVLListTest,PrevTest){
 *    //TODO
 *}
 *
 *TEST_F(_AVLListTest,ReorderListTest){
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *    EXPECT_EQ(1,1);
 *}
 *
 *TEST_F(_AVLListTest,ConsistencyCheckTest){
 *    //Checks to see if it is a valid AVL tree
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *    al.ConsistencyCheckTest();
 *
 *    EXPECT_EQ(1,1);
 *
 *}
 *
 *TEST_F(_AVLListTest,TraverserTest){
 *    //TODO
 *}
 *
 *TEST_F(_AVLListTest,toStrTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    EXPECT_STREQ("hi",(String*)al.toStr());
 *
 *}
 *
 *TEST_F(_AVLListTest,RetrieveTest) {
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    info = al.GetByIndexTest(3);
 *    EXPECT_EQ(7,info);
 *
 *}
 *
 *TEST_F(_AVLListTest,ClearTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    al.Clear(true);
 *
 *    EXPECT_EQ(-1,al.root);
 *
 *}
 *
 *TEST_F(_AVLListTest,InsertTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *    al.Insert((long)1,4,false);
 *    EXPECT_EQ(5,al.lLength);
 *
 *}
 *
 *TEST_F(_AVLListTest,InsertDataTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    long a = al.InsertData((long)1);
 *    EXPECT_EQ(1,a);
 *}
 *
 *TEST_F(_AVLListTest,HasDataTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    EXPECT_EQ(true, al.HasData());
 *}
 *
 *TEST_F(_AVLListTest,DeleteTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *
 *    al.Delete(0,true);
 *    EXPECT_EQ(3, al[0]);
 *
 *}
 *
 *TEST_F(_AVLListTest,countitemsTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ReorderList();
 *
 *    al.countitemsTest();
 *    EXPECT_EQ(3, al[0]);
 *
 *}
 *
 *
 *TEST_F(_AVLListXTest,ClearTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    //A call to Reorderlist may be
 *    alx.ReorderList();
 *
 *    alx.Clear(true);
 *
 *    EXPECT_EQ(-1,alx.root);
 *
 *}
 *
 *TEST_F(_AVLListXTest,toStrTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    EXPECT_STREQ("hi",(String*)alx.toStr());
 *
 *}
 *
 *TEST_F(_AVLListXLTest,toStrTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    EXPECT_STREQ("hi",(String*)alx.toStr());
 *
 *}
 *
 *TEST_F(_AVLListXLTest,ClearTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    //A call to Reorderlist may be
 *    alx.ReorderList();
 *
 *    alx.Clear(true);
 *
 *    EXPECT_EQ(-1,alx.root);
 *
 *}
 *
 *TEST_F(_AVLListXTest,InsertDataTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    //A call to Reorderlist may be
 *    alx.ReorderList();
 *
 *    long a = alx.InsertData((long)1);
 *    EXPECT_EQ(1,a);
 *
 *}
 *
 *TEST_F(_AVLListXLTest,InsertDataTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListXL alxl(&sl);
 *
 *    //A call to Reorderlist may be
 *    alxl.ReorderList();
 *
 *    long a = alxl.InsertData((long)1);
 *    EXPECT_EQ(1,a);
 *
 *}
 *
 *TEST_F(_AVLListXLTest,DeleteXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListXL alxl(&sl);
 *
 *    //A call to Reorderlist may be
 *    alxl.ReorderList();
 *
 *    alxl.DeleteXtra((long)0);
 *    EXPECT_EQ(3,a);
 *
 *}
 *
 *TEST_F(_AVLListXTest,DeleteXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    //A call to Reorderlist may be
 *    alx.ReorderList();
 *
 *    alx.DeleteXtra((long)0);
 *    EXPECT_EQ(3,a);
 *
 *}
 *
 *TEST_F(_AVLListXTest,PopulateFromListTest){
 *
 *    _List str_list = createStrList();
 *    _AVLListX alx;
 *
 *    alx.PopulateFromList(str_list)
 *
 *    _String* return_str = (_String*)list[4];
 *    EXPECT_STREQ("four", return_str->getStr());
 *}
 *
 *
 *TEST_F(_AVLListXTest,GetXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    //A call to Reorderlist may be
 *    alx.ReorderList();
 *
 *    long x = alx.GetXtra((long)0);
 *    EXPECT_EQ(1,x);
 *
 *}
 *
 *TEST_F(_AVLListXLTest,GetXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListXL alxl(&sl);
 *
 *    //A call to Reorderlist may be
 *    alxl.ReorderList();
 *
 *    long x = alxl.GetXtra((long)0);
 *    EXPECT_EQ(1,x);
 *
 *}
 *
 *TEST_F(_AVLListXTest,SetXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListX alx(&sl);
 *
 *    //A call to Reorderlist may be
 *    alx.ReorderList();
 *
 *    alx.SetXtra(0,13);
 *
 *    long x = alx.GetXtra((long)0);
 *    EXPECT_EQ(13,x);
 *
 *}
 *
 *TEST_F(_AVLListXLTest,SetXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLListXL alxl(&sl);
 *
 *    //A call to Reorderlist may be
 *    alxl.ReorderList();
 *
 *    alxl.SetXtra(0,13);
 *    long x = alxl.GetXtra((long)0);
 *    EXPECT_EQ(13,x);
 *
 *}
 */
}
