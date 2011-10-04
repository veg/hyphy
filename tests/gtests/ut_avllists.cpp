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

_List createAVLStrList() {

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
class AVLListTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    AVLListTest() {
        // You can do set-up work for each test here.
    }

    virtual ~AVLListTest() {
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

class AVLListXTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    AVLListXTest() {
        // You can do set-up work for each test here.
    }

    virtual ~AVLListXTest() {
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

class AVLListXLTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    AVLListXLTest() {
        // You can do set-up work for each test here.
    }

    virtual ~AVLListXLTest() {
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

TEST_F(AVLListTest,FindTest){
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

TEST_F(AVLListTest,FindLongTest){
    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.FindLong(3);
    EXPECT_EQ(3,info);
}

TEST_F(AVLListTest,FindBestTest){
    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.FindLong(3);
    EXPECT_EQ(3,info);
}

TEST_F(AVLListTest,NextTest){

    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 
    _SimpleList hist;

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.Next(3, hist);
    EXPECT_EQ(4,info);

}

TEST_F(AVLListTest,PrevTest){

    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 
    _SimpleList hist;

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.Prev(3, hist);
    EXPECT_EQ(2,info);

}

TEST_F(AVLListTest,FirstTest){
    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through left children until lData is 0

    long index; 
    _SimpleList sl; 
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    index = al.First();
    EXPECT_EQ(0,index);

}

TEST_F(AVLListTest,LastTest){

    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through left children until lData is 0
    long index; 
    _SimpleList sl; 
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    index = al.Last();
    EXPECT_EQ(10,index);
}

TEST_F(AVLListTest,GetByIndexTest){
    //AVL List always takes a pointer to a SimpleList
    //Keeps going down through right children until lData is 0
    long info; 
    _SimpleList sl; 

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    info = al.GetByIndex(3);
    EXPECT_EQ(3,info);
}

TEST_F(AVLListTest,ReorderListTest){

    long info; 
    _SimpleList sl,sl2;
    _AVLList al(&sl);
    
    al.Insert((BaseRef)5, 0, true, false);
    al.Insert((BaseRef)3, 0, true, false);
    al.Insert((BaseRef)7, 0, true, false);
    al.Insert((BaseRef)1, 0, true, false);
    al.Insert((BaseRef)9, 0, true, false);
    
    //A call to Reorderlist may be
    al.ReorderList(&sl2);
    EXPECT_EQ(3,sl2[0]);
}

/*
 *TEST_F(AVLListTest,ConsistencyCheckTest){
 *    //Checks to see if it is a valid AVL tree
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    sl.Populate(4,1,2);
 *    _AVLList al(&sl);
 *
 *    //A call to Reorderlist may be
 *    al.ConsistencyCheckTest();
 *
 *    EXPECT_EQ(1,1);
 *
 *}
 */


TEST_F(AVLListTest,TraverserTest){

    long info; 
    long t = 5;
    long r = 0;

    _SimpleList sl,sl2;
    _AVLList al(&sl);

    al.Insert((BaseRef)5, 0, true, false);
    al.Insert((BaseRef)3, 0, true, false);
    al.Insert((BaseRef)7, 0, true, false);
    al.Insert((BaseRef)1, 0, true, false);
    al.Insert((BaseRef)9, 0, true, false);

    //A call to Reorderlist may be
    al.Traverser(sl2, t, r);
    EXPECT_EQ(0,sl2[0]);
}

TEST_F(AVLListTest,toStrTest){

    long info;
    _SimpleList sl;

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);
    }

    _String* return_str = (_String*)al.toStr();
    EXPECT_STREQ("0\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n", return_str->getStr());
}

TEST_F(AVLListTest,RetrieveTest) {
    
    long info; 

    _SimpleList sl; 
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    info = (long)al.Retrieve(3);
    EXPECT_EQ(3,info);
}

TEST_F(AVLListTest,ClearTest){

    long info; 
    _SimpleList sl; 

    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    al.Clear(false);

    EXPECT_EQ(-1,al.root);

}

TEST_F(AVLListTest,InsertTest){
    long info; 
    _SimpleList sl; 
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    al.Insert((BaseRef)13,4,false);
    EXPECT_EQ(13,al.dataList->lData[11]);
    
    al.Insert((BaseRef)13,4,false);
    EXPECT_FALSE(al.dataList->lData[12]==13);
}

TEST_F(AVLListTest,InsertDataTest){
    long info;
    _SimpleList sl;
    _AVLList al(&sl);

    long a = al.InsertData((BaseRef)1, 1, false);
    EXPECT_EQ(0,a);
    
    a = al.InsertData((BaseRef)1, 1, false);
    EXPECT_EQ(1,a);
}

TEST_F(AVLListTest,HasDataTest){
    long info; 
    _SimpleList sl; 
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    EXPECT_EQ(true, al.HasData(2));
    //Do not understand the fail condition
    //EXPECT_EQ(false, al.HasData(20));
}


TEST_F(AVLListTest,DeleteTest){
    long info; 
    _SimpleList sl; 
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }

    al.Delete((BaseRef)1,false);
    EXPECT_EQ(0, al.dataList->lData[1]);
}

TEST_F(AVLListTest,countitemsTest){

    long info; 
    _SimpleList sl;
    _AVLList al(&sl);

    for(int i=0; i<=10; i++) {
        al.Insert((BaseRef)i, 0, true, false);    
    }
    
    EXPECT_EQ(11, al.countitems());
}

/*
 *TEST_F(AVLListXTest,toStrTest){
 *    //_AVLListX is supposed to be strings
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    _AVLListX al(&sl);
 *
 *    al.InsertData((BaseRef)"zero", 0, true, false);    
 *    al.InsertData((BaseRef)"one", 0, true, false);    
 *    al.InsertData((BaseRef)"two", 0, true, false);    
 *    al.InsertData((BaseRef)"three", 0, true, false);    
 *    al.Insert((BaseRef)"three", 0, true, false);    
 *
 *    _String* return_str = (_String*)al.toStr();
 *    EXPECT_STREQ("0\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n", return_str->getStr());
 *
 *}
 *
 *TEST_F(AVLListXLTest,toStrTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    _AVLListXL al(&sl);
 *
 *    for(int i=0; i<=10; i++) {
 *        al.Insert((BaseRef)i, 0, true, false);    
 *    }
 *
 *    _String* return_str = (_String*)al.toStr();
 *    EXPECT_STREQ("0\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n", return_str->getStr());
 *
 *}
 */

TEST_F(AVLListXLTest,ClearTest){

    long info; 
    _SimpleList sl; 
    _AVLList alxl(&sl);

    for(int i=0; i<=10; i++) {
        alxl.Insert((BaseRef)i, 0, true, false);    
    }

    alxl.Clear(false);
    EXPECT_EQ(-1,alxl.root);

}

TEST_F(AVLListXTest,InsertDataTest){

    long info;
    _SimpleList sl;
    _AVLListX al(&sl);

    long a = al.InsertData((BaseRef)1, 1, false);
    EXPECT_EQ(0,a);

}

TEST_F(AVLListXLTest,InsertDataTest){

    long info;
    _SimpleList sl;
    _AVLListX al(&sl);

    long a = al.InsertData((BaseRef)1, 1, false);
    EXPECT_EQ(0,a);

}

/*
 *TEST_F(AVLListXLTest,DeleteXtraTest){
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
 *    alxl.DeleteXtra(0);
 *    EXPECT_EQ(3,a);
 *
 *}
 *
 *TEST_F(AVLListXTest,DeleteXtraTest){
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
 */

/*
 *TEST_F(AVLListXTest,PopulateFromListTest){
 *
 *    //It seems as though if you don't use a pointer to 
 *    //a list, you are going to segfault
 *    //segfaults on command line but not in xcode 
 *
 *    _SimpleList sl;
 *    _AVLListX alx(&sl);
 *    
 *    _String test_string = _String("house,condo,hyphy");
 *    _String* sub_string = new _String(",");
 *    
 *    _List* result_list = test_string.Tokenize(sub_string);
 *    
 *    alx.PopulateFromList(*result_list);
 *    _String* return_str = (_String*)alx.dataList->lData[2];
 *    EXPECT_STREQ("hyphy", return_str->getStr());
 *}
 */


/*
 *TEST_F(AVLListXTest,SetAndGetXtraTest){
 *    //It seems as though if you don't use a pointer to 
 *    //a list, you are going to segfault
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    _AVLListX alx(&sl);
 *    _String test_string = _String("house,condo,hyphy");
 *    _String* sub_string = new _String(",");
 *    _List* result_list = test_string.Tokenize(sub_string);
 *    alx.PopulateFromList(*result_list);
 *
 *    alx.SetXtra(0,13);
 *    long x = alx.GetXtra((long)0);
 *    EXPECT_EQ(13,x);
 *
 *}
 */

/*
 *TEST_F(AVLListXLTest,SetAndGetXtraTest){
 *
 *    long info; 
 *    _SimpleList sl; 
 *
 *    _AVLListXL alxl(&sl);
 *
 *    for(int i=0; i<=10; i++) {
 *        alxl.Insert((BaseRef)i, 0, true, false);    
 *    }
 *
 *    alxl.SetXtra(0,13);
 *    long x = alxl.GetXtra((long)0);
 *    EXPECT_EQ(13,x);
 *}
 */

}
