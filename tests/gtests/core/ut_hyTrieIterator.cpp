/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon@cfenet.ubc.ca)
  Steven Weaver (sweaver@ucsd.edu)

Module Developers:
	Lance Hepler (nlhepler@gmail.com)
	Martin Smith (msmith@ucsd.edu)

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


#include "hy_list_reference.h"
#include "trieiterator.h"
#include "gtest/gtest.h"

namespace {

// The fixture for testing class Foo.
class _TrieIteratorTest : public ::testing::Test {

  protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _TrieIteratorTest() {
      // You can do set-up work for each test here.
    }

    virtual ~_TrieIteratorTest() {
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

    // Objects declared here can be used by all tests in the test case for Foo.
};


TEST_F(_TrieIteratorTest, MethodTests) {

  _Trie the_trie;

  the_trie.Insert("to", 1L);
  the_trie.Insert("inn", 2L);
  the_trie.Insert("tea", 3L);
  the_trie.Insert("ted", 4L);
  the_trie.Insert("ten", 5L);

  _TrieIterator the_trieiterator(the_trie);

  _String* result = the_trieiterator.First();
  EXPECT_STREQ("to", *result);
  result = the_trieiterator.Next();
  EXPECT_STREQ("tea", *result);
  result = the_trieiterator.Next();
  EXPECT_STREQ("ted", *result);
  result = the_trieiterator.Next();
  EXPECT_STREQ("ten", *result);
  result = the_trieiterator.Next();
  EXPECT_STREQ("inn", *result);

  result = the_trieiterator.Previous();
  EXPECT_STREQ("ten", *result);
  result = the_trieiterator.Previous();
  EXPECT_STREQ("ted", *result);
  result = the_trieiterator.Previous();
  EXPECT_STREQ("tea", *result);
  result = the_trieiterator.Previous();
  EXPECT_STREQ("to", *result);

  result = the_trieiterator.Previous();
  EXPECT_EQ(NULL, result);

  result = the_trieiterator.Next();
  EXPECT_EQ(NULL, result);


  result = the_trieiterator.Last();
  EXPECT_STREQ("inn", *result) << "the_trieiterator.First() should be first element of random_str_list";

  _Trie null_trie;
  _TrieIterator null_trieiterator(null_trie);

  _String* null = null_trieiterator.First();
  EXPECT_EQ(NULL, null);
  null = null_trieiterator.Last();
  EXPECT_EQ(NULL, null);


}


}
