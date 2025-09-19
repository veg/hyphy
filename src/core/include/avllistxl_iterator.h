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

/**
 * @brief A key-value pair for the AVLListXLIterator
 *
 */
class AVLListXLIteratorKeyValue {
private:
  long      index;
  BaseRef   object;
  const     _String * key;
public:
  /**
   * @brief Construct a new AVLListXLIteratorKeyValue object
   *
   * @param idx
   * @param k
   * @param obj
   */
  AVLListXLIteratorKeyValue (long idx, const _String* k, BaseRef obj) : index (idx), object (obj), key (k) {}
  /**
   * @brief Get the index
   *
   * @return long
   */
  long get_index (void) const {return index;}
  /**
   * @brief Get the object
   *
   * @return BaseRef
   */
  BaseRef get_object (void) const {return object;}
  /**
   * @brief Get the key
   *
   * @return _String const*
   */
  _String const* get_key (void) {return key;}
  
};

/**
 * @brief A C++ range compliant iterator for _AVLListXL
 *
 */
class AVLListXLIterator {
  
private:
  _AVLListXL const * the_list;
  long         current_location,
               traversal_lookup;
  _SimpleList  traversal_history;
  
  
public:
  /**
   * @brief Construct a new AVLListXLIterator object
   *
   * @param the_list
   */
  AVLListXLIterator (_AVLListXL const * the_list);
  
  /**
   * @brief Construct a new AVLListXLIterator object
   *
   * @param the_list
   */
  AVLListXLIterator (AVLListXLIterator const & the_list);

  /**
   * @brief Return an iterator pointing to the first element of the list
   *
   * @return AVLListXLIterator&
   */
   AVLListXLIterator& begin (void);
  
  /**
   * @brief Return an iterator pointing past the end of the list
   *
   * @return AVLListXLIterator&
   */
  AVLListXLIterator& end   (void);

  /**
   * @brief Increment the iterator
   *
   * @return AVLListXLIterator&
   */
  AVLListXLIterator& operator ++ (void);
  
  /**
   * @brief Compare two iterators
   *
   * @param compare
   * @return true
   * @return false
   */
  bool      operator == (AVLListXLIterator const & compare);
  
  /**
   * @brief Assignment operator
   *
   * @param assign_from
   * @return AVLListXLIterator&
   */
  AVLListXLIterator&       operator = (AVLListXLIterator const & assign_from);
  
  /**
   * @brief Compare two iterators
   *
   * @param compare
   * @return true
   * @return false
   */
   bool      operator != (AVLListXLIterator const & compare);
  
  /**
   * @brief Dereference the iterator
   *
   * @return AVLListXLIteratorKeyValue
   */
   AVLListXLIteratorKeyValue const    operator * (void);
  
    /**
     * @brief Check if the iterator is at the end
     *
     * @return true
     * @return false
     */
    bool is_done (void) const {
        return current_location == AVL_LISTXL_ITERATOR_ENDINDEX;
    }

};

//_____________________________________________________________________________

#endif
