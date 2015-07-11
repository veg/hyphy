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
#include "ut_lists.h"

#include "hy_strings.h"
#include "list.h"

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
class ListTest : public ::testing::Test
{

protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    ListTest() {
        // You can do set-up work for each test here.
    }

    virtual ~ListTest() {
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


TEST_F(ListTest,_NothingListTest){
    _List list = _List();
    EXPECT_EQ(list.lLength,0);
}

TEST_F(ListTest,_LengthConstructorListTest){
    _List str_list = createStrList();
    unsigned long l = 7;
    EXPECT_EQ(str_list.lLength,l);
}

TEST_F(ListTest,_StackCopyConstructorListTest){

    _List str_list = createStrList();
    _List list = _List(str_list, 0, 4); 

    _String* return_str = (_String*)list[4];
    EXPECT_STREQ("four", return_str->getStr());

}

TEST_F(ListTest,_SubStrConstructorListTest){

    _String* str = new _String(",one,two,three");
    _List list = _List((BaseRef)str, ','); 
    _String* return_string = (_String*)list[1];
    EXPECT_STREQ("one", return_string->getStr());
}

TEST_F(ListTest,_DataConstructorListTest){
    //one member list
    _String* string = new _String("one");
    _List list = _List((BaseRef)string);
    _String* return_string = (_String*)list[0];

    EXPECT_STREQ("one", return_string->getStr());

}

TEST_F(ListTest,paranthesisTest){

    _List str_list = createStrList();
    _String* return_string = (_String*)str_list(4);

    EXPECT_STREQ("four", return_string->getStr());
}

TEST_F(ListTest,EqualOpTest){

    _List str_list = createStrList();
    _List list = _List();
    list = str_list;

    _String* return_string = (_String*)list[4];

    EXPECT_STREQ("four", return_string->getStr());
}

TEST_F(ListTest,EqualTest){

    _List str_list = createStrList();
    _List list = str_list; 
    _List l2;

    EXPECT_EQ(true,list.Equal(str_list));
    EXPECT_EQ(true,list==str_list);

    EXPECT_EQ(false,str_list.Equal(l2));
    list.AppendNewInstance(new _String("zero"));
    EXPECT_EQ(false,str_list.Equal(list));
}

TEST_F(ListTest,AmpersandOpTest){

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


    _List* l2  = new _List(); 
    _List* al2 = new _List();
    _List rl2 = *l2 & *al2;
    EXPECT_EQ(0, rl2.lLength);
}

TEST_F(ListTest,Ampersand2OpTest){

    _List str_list = createStrList();
    _String* str = new _String("one");

    _List result_list = str_list & (_String*)str;

    _String* result_string = (_String*)result_list[result_list.lLength-1];

    EXPECT_STREQ("one",result_string->getStr());
}

TEST_F(ListTest,DoubleAmpersandOpTest){

    _List str_list = createStrList();

    _String* str = new _String("one");
    str_list && (_String*)str;

    _String* result_string = (_String*)str_list[str_list.lLength-1];
    
    EXPECT_STREQ("one",result_string->getStr());

}

TEST_F(ListTest,DoubleAmpersand2OpTest){
    _List str_list = createStrList();
    str_list && "one"; 
    _String* result_string = (_String*)str_list[str_list.lLength-1];
    EXPECT_STREQ("one",result_string->getStr());
}

TEST_F(ListTest,DoubleLessOpTest){

    //Insert Element
    _List str_list = createStrList();
    _String* str = new _String("one");
    str_list << str;

    _String* result_string = (_String*)str_list[str_list.lLength-1];
    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(ListTest,DoubleLess2OpTest){
    //Append Operator

    _List str_list = createStrList();
    _String string ("one,two,three");

    _List* list = string.Tokenize(','); 

    _String append_string("four,five,six");
    _List* append_list = append_string.Tokenize(','); 

    _String* expected_string = new _String("six");

    *list << *append_list;

    _List result_list = *list;

    _String* result_string = (_String*)result_list[5];
    EXPECT_STREQ(expected_string->getStr(), result_string->getStr());

}

TEST_F(ListTest,AppendNewInstanceTest){

    //Insert Element
    _List str_list = createStrList();
    str_list.AppendNewInstance(new _String("one"));
    _String* result_string = (_String*)str_list[7];

    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(ListTest,PlaceTest){
    //Place Test

    //TODO: Figure out when laLength would be larger than lLength
    _List str_list = createStrList();
    str_list.Place(new _String("one"));

    _String* result_string = (_String*)str_list[7];
    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(ListTest,InsertElementTest){

    _List str_list = createStrList();
    str_list.InsertElement(new _String("one"),3,true);
    _String* result_string = (_String*)str_list[3];
    EXPECT_STREQ("one", result_string->getStr());
}

TEST_F(ListTest,getStrTest){
    _List str_list = createStrList();

    _String* str = new _String("{zero,one,two,three,four,five,six}");
    _String* result_string = (_String*)str_list.toStr();

    EXPECT_STREQ(str->getStr(), result_string->getStr());
}

TEST_F(ListTest, makeDynamicTest){
    _List str_list = createStrList();
    _List* dynamic_list = (_List*)str_list.makeDynamic();
    EXPECT_EQ(7,dynamic_list->lLength);
}

TEST_F(ListTest,toFileStrTest){
    //Doesn't need unit tested
}


TEST_F(ListTest,bumpNInstTest){
    //TODO
    _List str_list = createStrList();
    str_list.bumpNInst();

    _String* result_string = (_String*)str_list[3];
    EXPECT_STREQ("three", result_string->getStr());
}

TEST_F(ListTest,FindTest){
    //This cast the object as a string and then checks if it's equal. 
    //Why do we send in a BaseRef when we could just pass in a _String if it's strings only.
    _List str_list = createStrList();
    
    _String* needle = new _String("two");

    int index = str_list.Find((BaseRef)needle);
    EXPECT_EQ(2, index);

    index = str_list.Find((BaseRef)"seventeen");
    EXPECT_EQ(-1, index);
}

TEST_F(ListTest,FindStringTest){
    //Find the position of a search string in the list of strings (ONLY)
    long upTo = -1;
    _List str_list = createStrList();
    
    _String* needle = new _String("two");

    int index = str_list.FindString((BaseRef)needle, 0, false, upTo); 
    EXPECT_EQ(2, index);

    index = str_list.FindString((BaseRef)needle, 0, true, upTo); 
    EXPECT_EQ(2, index);

    _String* needle2 = new _String("france");
    index = str_list.FindString((BaseRef)needle2, 0, true, upTo); 
    EXPECT_EQ(-1, index);

} 

TEST_F(ListTest,JoinTest){
    //Append Operator
    _String* str = new _String("one,two,three");
    _String* spacer = new _String(";");
    _List* list = new _List((BaseRef)str, ','); 
    _String* return_string = (_String*)list->Join((BaseRef)spacer);
    EXPECT_STREQ("one;two;three", return_string->getStr());
}

TEST_F(ListTest,BinaryFindTest){

    //Find the position of a search string in the list of strings (ONLY)
    int upTo = -1;
    _List str_list = createStrList();
    _String* needle = new _String("six");
    int index = str_list.BinaryFind((BaseRef)needle); 
    EXPECT_EQ(-4, index);

}

TEST_F(ListTest,BinaryInsertTest){

    //TODO: Debug Coverage

    _List str_list = createStrList();

    //Place Test
    str_list.BinaryInsert((BaseRef)new _String("one"));

    _String* result_string = (_String*)str_list[3];
    EXPECT_STREQ("three", result_string->getStr());
}

TEST_F(ListTest,CompareTest){
    _List str_list = createStrList();

    _String* test = new _String("hyphy");
    EXPECT_EQ(-1, str_list.Compare((BaseRef)test,1));
}

TEST_F(ListTest,Compare2Test){
    _List str_list = createStrList();
    EXPECT_EQ(1, str_list.Compare((long)0,(long)1));
}

TEST_F(ListTest,FreeUpMemoryTest){
    _List str_list = createStrList();
    EXPECT_EQ(0,str_list.FreeUpMemory(16));
} 

TEST_F(ListTest,ClearTest){

    _List str_list = createStrList();
    str_list.InsertElement(new _String("one"),3,true);
    str_list.Clear(true);

    EXPECT_EQ(0,str_list.lLength);
} 


TEST_F(ListTest,DeleteTest){

    _List str_list = createStrList();
    str_list.Delete(0);
    _String* return_string = (_String*)str_list[1];
    EXPECT_STREQ("two", return_string->getStr());
} 

TEST_F(ListTest,DeleteListTest){
    //See also Clear
    //It looks like this deletes a subset of a list based on a list of indices
    _List str_list = createStrList();
    _SimpleList to_delete; 
    createIntList(3, to_delete);

    str_list.DeleteList(to_delete);
    EXPECT_EQ(3,str_list.lLength);
} 


TEST_F(ListTest,ReplaceTest){

    _List str_list = createStrList();

    str_list.Replace(1, new _String("help"), false);
    _String* return_string = (_String*)str_list[1];
    EXPECT_STREQ("help", return_string->getStr());

    str_list.Replace(1, new _String("two"), false);
    _String* return_string2 = (_String*)str_list[1];
    EXPECT_STREQ("two", return_string2->getStr());
} 


TEST_F(ListTest,IntersectTest){
    // compute the union of two sorted lists
    // each repeat appears exactly once

    _List str_list = createStrList();

    _List l1 = str_list;
    _List l2 = str_list;
    _SimpleList idx; 
    _SimpleList idx2;

    _List l3;
    l3.Intersect(l1, l2, &idx, &idx2);

    _String* return_string = (_String*)str_list[1];
    EXPECT_STREQ("one",return_string->getStr());

}

}
