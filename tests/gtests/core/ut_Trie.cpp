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


#include "trie.h"
#include "hy_list_reference.h"
#include "gtest/gtest.h"

namespace {

// The fixture for testing class Foo.
class _TrieTest : public ::testing::Test {

  protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    _TrieTest() {
      // You can do set-up work for each test here.
    }

    virtual ~_TrieTest() {
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


TEST_F(_TrieTest, ConstructorsTest)
{
  // Normal
  _StringBuffer test_s ("hyphy");

  _Trie the_test(&test_s);
  EXPECT_EQ ('\0', the_test.alphabet()[0]) << "Should be \\0";
  EXPECT_EQ ('h', the_test.alphabet()[1]) << "Should be h";
  EXPECT_EQ ('p', the_test.alphabet()[2]) << "Should be p";
  EXPECT_EQ ('y', the_test.alphabet()[3]) << "Should be y";

  //Insert an item that has an invalid character
  EXPECT_EQ(HY_TRIE_INVALID_LETTER, the_test.Insert("@handle#rocks", 1L)) << "Insertion of a string with invalid characters was allowed to proceed";
  EXPECT_NE(HY_TRIE_INVALID_LETTER, the_test.Insert ("hppy", 0L)) << "Insertion of a string with all valid characters failed";

}

TEST_F(_TrieTest, MethodTests)
{

  _Trie the_trie;
  the_trie.Insert("Dirichlet", 1L);
  the_trie.Insert("Gaussian", 2L);
  the_trie.Insert("Wishart", 3L);
  the_trie.Insert("InverseWishart", 4L);
  the_trie.Insert("Multinomial", 5L);
  EXPECT_STREQ("{\"Dirichlet\":1\n\"Gaussian\":2,\n\"Wishart\":3,\n\"InverseWishart\":4,\n\"Multinomial\":5,\n}", (char *)the_trie.toStr().getStr());

  //Clear the trie
  the_trie.clear();
  EXPECT_STREQ("{}", (char *)the_trie.toStr().getStr()) << "String should be empty after clearing";

  //Insert again
  the_trie.Insert("Dirichlet", 1L);
  the_trie.Insert("Gaussian", 2L);
  the_trie.Insert("Wishart", 3L);
  the_trie.Insert("InverseWishart", 4L);
  the_trie.Insert("Multinomial", 5L);
  EXPECT_STREQ("{\"Dirichlet\":1\n\"Gaussian\":2,\n\"Wishart\":3,\n\"InverseWishart\":4,\n\"Multinomial\":5,\n}", (char *)the_trie.toStr().getStr());

  the_trie.clear();
  the_trie.Insert("apple", 1L);
  the_trie.Insert("apropos", 2L);
  the_trie.Insert("banana", 3L);
  the_trie.Insert("orange", 4L);
  EXPECT_STREQ("banana", the_trie.RetrieveKeyByPayload(3L));
  EXPECT_EQ(3L, the_trie.GetValueFromString("banana"));

  EXPECT_EQ(26, the_trie.countitems());

  bool deleted = the_trie.Delete("banana");
  EXPECT_EQ(true, deleted);
  EXPECT_GT(0, the_trie.Find("banana"));
  EXPECT_STREQ("0", the_trie.RetrieveKeyByPayload(3L));



  // Dreate a random strings and ensure the count is correct
  the_trie.clear();

  long random_list_length = 16L + genrand_int32() % 128L;

  _hyListReference <_String> random_str_list;
  _hyListNumeric <long> random_long_list;
  _String* alph = new _String("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

  for (long i = 0; i <= random_list_length; i++) {
    long random_length = 16L + genrand_int32() % 128L;
    _String* random_string = new _String(_String::Random (random_length, alph));
    random_str_list.AppendNewInstance(random_string);
    random_long_list.append(random_length);
    the_trie.Insert(*random_string, random_length);
  }

  long random_index = genrand_int32() % random_long_list.Length();
  EXPECT_LT(random_long_list.Sum(), the_trie.countitems() - 1);
  _String str(*random_str_list.Element(random_index));
  EXPECT_LT(0, the_trie.Find(str));
  EXPECT_GT(random_long_list.Sum(), the_trie.Find(str));



  // Delete every item from a list
  deleted = the_trie.Delete(random_str_list);
  EXPECT_EQ(true, deleted);

   

  //Insert the entire list at once
  the_trie.clear();

  _Trie second_trie;  
  long how_many = second_trie.Insert(random_str_list, random_long_list);
  EXPECT_EQ(random_list_length, how_many - 1);
  EXPECT_LT(random_long_list.Sum(), second_trie.countitems() - 1);
  EXPECT_LT(0, second_trie.Find(str));
  EXPECT_GT(random_long_list.Sum(), second_trie.Find(str));

  delete alph;
  random_str_list.Clear();

}
  
  TEST_F(_TrieTest, LargeTrieTest) {
    _String alph ("ACGT");
    _Trie nucleotide_trie (&alph);
    long i = 0L;
    
    for (; i <= (1<<20L); i++) {
      nucleotide_trie.Insert(_String::Random (20, &alph), i);
    }
    
    //ASSERT_EQ (i, nucleotide_trie.Length());
    
  }



}
