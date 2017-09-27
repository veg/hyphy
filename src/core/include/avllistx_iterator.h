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

class AVLListXIteratorKeyValue {
private:
  long      index;
  long      value;
    
public:
  AVLListXIteratorKeyValue (long idx, long val) : index (idx), value (val) {}
  long get_index (void) const {return index;}
  long get_value (void) const {return value;}
  
};

class AVLListXIterator {
  /** 
   * A C++ range compiant iterator for _AVLListX
   */
  
private:
  _AVLListX const * the_list;
  long         current_location,
               traversal_lookup;
  _SimpleList  traversal_history;
  
  
public:
  AVLListXIterator (_AVLListX const * the_list);
  /**
    * Build an iterator on the underlying list
    * 
    * @param the_list the list to iterate over
   */
  

  AVLListXIterator (AVLListXIterator const & the_list);
  /**
   * Copy constructor
   *
   * @param the_list the iterator to copy from
   */

   AVLListXIterator& begin (void);
  
  /**
   * Return an iterator pointing to the first element of the list
   * @return the iterator poiting the to the first element fo the list
  */
  
  AVLListXIterator& end   (void);

  /**
   * Return an iterator pointing past the end of the list
   * @return the iterator poiting past the end of the list
   */

  AVLListXIterator& operator ++ (void);

  /**
   * Return an iterator pointing past the end of the list
   * @return the iterator poiting past the end of the list
  */

  
  bool      operator == (AVLListXIterator const & compare);
  
  /**
   * Compare two iterators
   * @return true if the iterators are the same
   */

  AVLListXIterator&       operator = (AVLListXIterator const & assign_from);
  
  /**
   * Assignment operator
   * @param assign_from assign from this object
   * @return *this
   */

   bool      operator != (AVLListXIterator const & compare);
  
  /**
   * Compare two iterators
   * @return true if the iterators are not the same
   */

   AVLListXIteratorKeyValue const    operator * (void);
  
  /**
   * Return an object pointed to by the iterator
   * @return the object pointed to by the iterator
   */

};

//_____________________________________________________________________________

#endif
