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

#ifndef _AVLListXIteratorHeader_
#define _AVLListXIteratorHeader_

#include "avllistx.h"

class _AVLListX; // forward declaration

//_____________________________________________________________________________

/**
 * @brief A key-value pair for the AVLListXIterator
 *
 */
class AVLListXIteratorKeyValue {
private:
  long      index;
  long      value;
    
public:
  /**
   * @brief Construct a new AVLListXIteratorKeyValue object
   *
   * @param idx
   * @param val
   */
  AVLListXIteratorKeyValue (long idx, long val) : index (idx), value (val) {}
  /**
   * @brief Get the index
   *
   * @return long
   */
  long get_index (void) const {return index;}
  /**
   * @brief Get the value
   *
   * @return long
   */
  long get_value (void) const {return value;}
  
};

/**
 * @brief A C++ range compliant iterator for _AVLListX
 *
 */
class AVLListXIterator {
  
private:
  _AVLListX const * the_list;
  long         current_location,
               traversal_lookup;
  _SimpleList  traversal_history;
  
  
public:
  /**
   * @brief Construct a new AVLListXIterator object
   *
   * @param the_list
   */
  AVLListXIterator (_AVLListX const * the_list);
  
  /**
   * @brief Construct a new AVLListXIterator object
   *
   * @param the_list
   */
  AVLListXIterator (AVLListXIterator const & the_list);

  /**
   * @brief Return an iterator pointing to the first element of the list
   *
   * @return AVLListXIterator&
   */
   AVLListXIterator& begin (void);
  
  /**
   * @brief Return an iterator pointing past the end of the list
   *
   * @return AVLListXIterator&
   */
  AVLListXIterator& end   (void);

  /**
   * @brief Increment the iterator
   *
   * @return AVLListXIterator&
   */
  AVLListXIterator& operator ++ (void);
  
  /**
   * @brief Compare two iterators
   *
   * @param compare
   * @return true
   * @return false
   */
  bool      operator == (AVLListXIterator const & compare);
  
  /**
   * @brief Assignment operator
   *
   * @param assign_from
   * @return AVLListXIterator&
   */
  AVLListXIterator&       operator = (AVLListXIterator const & assign_from);
  
  /**
   * @brief Compare two iterators
   *
   * @param compare
   * @return true
   * @return false
   */
   bool      operator != (AVLListXIterator const & compare);
  
  /**
   * @brief Dereference the iterator
   *
   * @return AVLListXIteratorKeyValue
   */
   AVLListXIteratorKeyValue const    operator * (void);
  
};

//_____________________________________________________________________________

#endif
