/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef _TrieIteratorHeader_
#define _TrieIteratorHeader_

#include "trie.h"
#include "hy_string_buffer.h"
#include "function_templates.h"

class _Trie; // forward declaration

//_____________________________________________________________________________

class TrieIteratorKeyValue {
private:
  _String             key;
  long                value;
    
public:
  TrieIteratorKeyValue (_String * k, long val) : key (k), value (val) {}

  _String const& get_key    (void) const {return key;}
  long           get_value (void) const {return value;}
  
};

class TrieIterator {
  /** 
   * A C++ range compiant iterator for _Trie
   */
  
private:
  _Trie        const * the_trie;
  _SimpleList  traversal_history,
               *root_list;
    
  
  bool         Next ();
  
public:
  TrieIterator (_Trie const * trie);
  /**
    * Build an iterator on the underlying list
    * 
    * @param the_list the list to iterate over
   */
  

  TrieIterator (TrieIterator const & the_list);
  /**
   * Copy constructor
   *
   * @param the_list the iterator to copy from
   */

   TrieIterator& begin (void);
  
  /**
   * Return an iterator pointing to the first element of the list
   * @return the iterator poiting the to the first element fo the list
  */
  
  TrieIterator& end   (void);

  /**
   * Return an iterator pointing past the end of the list
   * @return the iterator poiting past the end of the list
   */

  TrieIterator& operator ++ (void);

  /**
   * Return an iterator pointing past the end of the list
   * @return the iterator poiting past the end of the list
  */

  
  bool      operator == (TrieIterator const & compare);
  
  /**
   * Compare two iterators
   * @return true if the iterators are the same
   */

  TrieIterator&       operator = (TrieIterator const & assign_from);
  
  /**
   * Assignment operator
   * @param assign_from assign from this object
   * @return *this
   */

   bool      operator != (TrieIterator const & compare);
  
  /**
   * Compare two iterators
   * @return true if the iterators are not the same
   */

   TrieIteratorKeyValue const    operator * (void);
  
  /**
   * Return an object pointed to by the iterator
   * @return the object pointed to by the iterator
   */

};

//_____________________________________________________________________________

#endif
