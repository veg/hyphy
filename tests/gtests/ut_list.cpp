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
#include "ut_lists.h"

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

//_List createStrList();
//_SimpleList createIntList(int n, _SimpleList& int_list);

_List createStrList() {

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

void createIntList(int n, _SimpleList& int_list) {

    for(int i=0;i<=n;i++) {
        int_list << i;
    }
}

namespace
{

// The fixture for testing class Foo.
class _ListTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _ListTest() {
        // You can do set-up work for each test here.
    }

    virtual ~_ListTest() {
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


//does nothing
TEST_F(_ListTest,_NothingListTest){
    _List list = _List();
    EXPECT_EQ(list.lLength,0);
}

TEST_F(_ListTest,_LengthConstructorListTest){
    _List str_list = createStrList();
    unsigned long l = 7;
    EXPECT_EQ(str_list.lLength,l);
}

TEST_F(_ListTest,_StackCopyConstructorListTest){
    _List str_list = createStrList();
    _List list = _List(str_list, 0, 4); 

    _String* return_str = (_String*)list[4];
    EXPECT_STREQ("four", return_str->getStr());
}

TEST_F(_ListTest,_SubStrConstructorListTest){

    _String* str = new _String("one,two,three");
    _List list = _List((BaseRef)str, ','); 
    _String* return_string = (_String*)list[0];

    //TODO: Contemplate whether we should be getting information from this class in this manner
    EXPECT_STREQ("one", return_string->getStr());
}

TEST_F(_ListTest,_DataConstructorListTest){
    //one member list
    _String* string = new _String("one");
    _List list = _List((BaseRef)string);
    _String* return_string = (_String*)list[0];

    EXPECT_STREQ("one", return_string->getStr());

}


TEST_F(_ListTest,paranthesisTest){

    _List str_list = createStrList();
    _String* return_string = (_String*)str_list(4);

    EXPECT_STREQ("four", return_string->getStr());
}

TEST_F(_ListTest,EqualOpTest){

    _List str_list = createStrList();
    _List list = str_list;

    _String* return_string = (_String*)list[4];

    EXPECT_STREQ("four", return_string->getStr());
}

TEST_F(_ListTest,EqualTest){

    _List str_list = createStrList();
    _List list = str_list; 

    EXPECT_EQ(true,list.Equal(str_list));
    EXPECT_EQ(false,str_list.Equal(str_list));
}

TEST_F(_ListTest,AmpersandOpTest){
    //Append Operator
    _List str_list = createStrList();

    _String string ("one,two,three");
    _List* list = string.Tokenize(','); 

    _String append_string("four,five,six");
    _List* append_list = append_string.Tokenize(','); 

    _String expected_string("one,two,three,four,five,six");
    _List* expected_list = expected_string.Tokenize(','); 

    _List result_list = *list & *append_list;

    _String* result_string = (_String*)result_list[0];

    EXPECT_STREQ("one", result_string->getStr());
} 

TEST_F(_ListTest,DoubleAmpersandOpTest){

    _List str_list = createStrList();

    _String* str = new _String("one");
    str_list && (_String*)str;

    _String* result_string = (_String*)str_list[str_list.lLength-1];
    
    EXPECT_STREQ("one",result_string->getStr());

}


TEST_F(_ListTest,DoubleLessOpTest){

    //Insert Element
    _List str_list = createStrList();
    _String* str = new _String("one");
    str_list << str;

    _String* result_string = (_String*)str_list[str_list.lLength-1];
    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(_ListTest,DoubleLess2OpTest){
    //Append Operator

    _List str_list = createStrList();
    _String string ("one,two,three");

    _List* list = string.Tokenize(','); 

    _String append_string("four,five,six");
    _List* append_list = append_string.Tokenize(','); 

    _String* expected_string = new _String("one,two,three,four,five,six");

    *list << *append_list;

    _List result_list = *list;

    _String* result_string = (_String*)result_list[5];
    EXPECT_STREQ(expected_string->getStr(), result_string->getStr());

}

TEST_F(_ListTest,AppendNewInstanceTest){

    //Insert Element
    _List str_list = createStrList();
    str_list.AppendNewInstance(new _String("one"));
    _String* result_string = (_String*)str_list[5];

    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(_ListTest,PlaceTest){
    //Place Test
    _List str_list = createStrList();
    str_list.Place(new _String("one"));

    _String* result_string = (_String*)str_list[5];
    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(_ListTest,InsertElementTest){

    _List str_list = createStrList();

    //Place Test
    str_list.InsertElement(new _String("one"),3,true);
    _String* result_string = (_String*)str_list[3];
    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(_ListTest,getStrTest){
    _List str_list = createStrList();

    _String* str = new _String("{1,2,3,4,5,6,7,8,9,10,11,12,13,14}");
    _String* result_string = (_String*)str_list.toStr();

    EXPECT_STREQ(str->getStr(), result_string->getStr());
}

//TEST_F(_ListTest,toFileStrTest){
    ////TODO

//}


TEST_F(_ListTest,bumpNInstTest){

    _List str_list = createStrList();
    str_list.bumpNInst();

    _String* result_string = (_String*)str_list[3];
    EXPECT_STREQ("three", result_string->getStr());
}

TEST_F(_ListTest,FindTest){
    //This cast the object as a string and then checks if it's equal. 
    _List str_list = createStrList();

    int index = str_list.Find((BaseRef)"two");
    EXPECT_EQ(2, index);

    index = str_list.Find((BaseRef)"seventeen");
    EXPECT_EQ(-1, index);
}

TEST_F(_ListTest,FindStringTest){
    //Find the position of a search string in the list of strings (ONLY)
    long upTo = -1;
    _List str_list = createStrList();

    int index = str_list.FindString((BaseRef)"two", 0, false, upTo); 
    EXPECT_EQ(2, index);

    index = str_list.FindString((BaseRef)"two", 0, true, upTo); 
    EXPECT_EQ(2, index);
} 

TEST_F(_ListTest,JoinTest){
    //Append Operator
    _String* str = new _String("one,two,three");
    _List* list = new _List((BaseRef)str, ','); 
    //_String* return_string = (_String*)list->Join((BaseRef)"1");
    EXPECT_STREQ("1", str->getStr());
}

TEST_F(_ListTest,BinaryFindTest){
    //Find the position of a search string in the list of strings (ONLY)
    //Same as FindString
    int upTo = -1;
    _List str_list = createStrList();

    int index = str_list.FindString((BaseRef)"two", 0, false, upTo); 
    EXPECT_EQ(2, index);

    index = str_list.FindString((BaseRef)"two", 0, true, upTo); 
    EXPECT_EQ(2, index);
} 

TEST_F(_ListTest,BinaryInsertTest){
    //Binary Insert 
    //TODO
}

TEST_F(_ListTest,CompareTest){
    _List str_list = createStrList();

    _String* test = new _String("hyphy");
    EXPECT_EQ(false, str_list.Compare((BaseRef)test,1));
}

TEST_F(_ListTest,FreeUpMemoryTest){
    _List str_list = createStrList();
    EXPECT_EQ(16,str_list.FreeUpMemory(16));
} 

TEST_F(_ListTest,ClearTest){

    _List str_list = createStrList();
    str_list.Clear(true);

    EXPECT_EQ(0,str_list.lLength);
} 


TEST_F(_ListTest,DeleteTest){

    _List str_list = createStrList();
    str_list.Delete(0);
    _String* return_string = (_String*)str_list[1];
    EXPECT_STREQ("two", return_string->getStr());
} 

TEST_F(_ListTest,DeleteListTest){
    //See also Clear
    //It looks like this deletes a subset of a list based on a list of indices
    _List str_list = createStrList();
    _SimpleList to_delete; 
    createIntList(3, to_delete);

    str_list.DeleteList(to_delete);
    EXPECT_EQ(3,str_list.lLength);
} 

TEST_F(_ListTest,ReplaceTest){
    _List str_list = createStrList();
    str_list.Replace(1, new _String("help"), false);
    _String* return_string = (_String*)str_list[1];
    EXPECT_STREQ("help", return_string->getStr());
} 

TEST_F(_ListTest,IntersectTest){
// compute the union of two sorted lists
// each repeat appears exactly once

    _List str_list = createStrList();

    _List l1 = str_list;
    _List l2 = str_list;
    _SimpleList* idx; 
    _SimpleList* idx2;

    _List l3;
    //l3.Intersect(l1, l2, idx, idx2);

    _String* return_string = (_String*)str_list[1];
    EXPECT_STREQ("one",return_string->getStr());
}

}
