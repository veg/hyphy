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

#ifndef _AVLListXLIteratorHeader_
#define _AVLListXLIteratorHeader_

#include "avllistxl.h"

#define  AVL_LISTXL_ITERATOR_ENDINDEX (-1L)

//_____________________________________________________________________________

class AVLListXLIteratorKeyValue {
private:
  long      index;
  BaseRef   object;
  const     _String * key;
public:
  AVLListXLIteratorKeyValue (long idx, const _String* k, BaseRef obj) : index (idx), object (obj), key (k) {}
  long get_index (void) const {return index;}
  BaseRef get_object (void) const {return object;}
  _String const* get_key (void) {return key;}
  
};

class AVLListXLIterator {
  /** 
   * A C++ range compiant iterator for _AVLListXL
   */
  
private:
  _AVLListXL const * the_list;
  long         current_location,
               traversal_lookup;
  _SimpleList  traversal_history;
  
  
public:
  AVLListXLIterator (_AVLListXL const * the_list);
  /**
    * Build an iterator on the underlying list
    * 
    * @param the_list the list to iterate over
   */
  

  AVLListXLIterator (AVLListXLIterator const & the_list);
  /**
   * Copy constructor
   *
   * @param the_list the iterator to copy from
   */

   AVLListXLIterator& begin (void);
  
  /**
   * Return an iterator pointing to the first element of the list
   * @return the iterator poiting the to the first element fo the list
  */
  
  AVLListXLIterator& end   (void);

  /**
   * Return an iterator pointing past the end of the list
   * @return the iterator poiting past the end of the list
   */

  AVLListXLIterator& operator ++ (void);

  /**
   * Return an iterator pointing past the end of the list
   * @return the iterator poiting past the end of the list
  */

  
  bool      operator == (AVLListXLIterator const & compare);
  
  /**
   * Compare two iterators
   * @return true if the iterators are the same
   */

  AVLListXLIterator&       operator = (AVLListXLIterator const & assign_from);
  
  /**
   * Assignment operator
   * @param assign_from assign from this object
   * @return *this
   */

   bool      operator != (AVLListXLIterator const & compare);
  
  /**
   * Compare two iterators
   * @return true if the iterators are not the same
   */

   AVLListXLIteratorKeyValue const    operator * (void);
  
  /**
   * Return an object pointed to by the iterator
   * @return the object pointed to by the iterator
   */
    
    bool is_done (void) const {
        return current_location == AVL_LISTXL_ITERATOR_ENDINDEX;
    }

};

//_____________________________________________________________________________

#endif
