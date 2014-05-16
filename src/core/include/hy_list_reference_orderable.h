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

#ifndef _HY_LIST_REFERENCE_
#define _HY_LIST_REFERENCE_

#include "hy_list_orderable.h"
#include "hy_list_reference.h"

/*
 
 A resizable list which stores PAYLOAD* types and tracks
 references
 
 PAYLOAD objects are assumed to derive from BaseObject and 
 implement the Compare and Equal operations (only if sorting 
 is needed)
 
 */

template <typename PAYLOAD>
class _hyListReferenceOrderable : public virtual _hyListOrderable <PAYLOAD*>, public virtual _hyListReference <PAYLOAD*> {

public:

  /**
  * A constructor.
  * A simple constructor that does nothing
  */
  _hyListReferenceOrderable();


  /**
   * A _hyListReference constructor for X items
   * @param items how many items should the storage be allocated for
   */
  _hyListReferenceOrderable(const unsigned long items);

  /**
   * A _hyListReference constructor from a single item
   * The reference counter for the item will be incremented
   * @param item the element to add to the list
   */
  _hyListReferenceOrderable(const PAYLOAD&);


  /**
   * A _hyListReference constructor from a single item
   * The reference counter for the item will NOT be incremented
   * @param item the element to add to the list
   */
  _hyListReferenceOrderable(const PAYLOAD*);
  
  /**
  * Stack copy contructor
  * @param l List to be copied
  * @param from Beginning index to copy from
  * @param to Last index to copy to
  */
  _hyListReferenceOrderable(const _hyListReference<PAYLOAD> &, const long = 0, const long = HY_LIST_INSERT_AT_END);


  /**
   * Data constructor from a static list of PAYLOAD objects
   * @param const unsigned long the number of PAYLOAD arguments supplied
   * to the constructor in the second argument
   * @param const PAYLOAD: the static list of items to pass to the constructor
   */
  
  _hyListReferenceOrderable (const unsigned long, const PAYLOAD * []);

  /**
   * Construct a list of substrings from the original string separated by char
   * ONLY SPECIALIZED FOR PAYLOAD == _String
   * \n\n \b Example: \code  list = _hyListReference<_String>(
   * _String("one,two,three"), ','); returns ['one','two','three']; \endcode
   * @param the_string The substring to be parsed, remember to cast it as a BaseRef
   * @param separator The separator for the string
   */
  
  //_hyListReference(const _String& the_string, const char separator);

  virtual _hyListReferenceOrderable (void);
  
  /**
   * Add a new reference to the list WITHOUT incrementing the reference counter
   * \n\n \b Example: \code  list.AppendNewInstance (new _String ("something"))
   * @param the_ref The reference to be added.
  */
  
  /**
   * Compares two elements of the list
   * @param i The index to compare
   * @param j The second index to compare
   * @return -1 if i<j, 0 if i==j, or 1 if i>j
   */
   virtual long Compare(const long, const long) const;
  
    //   virtual BaseRef makeDynamic(void) const;

  
};

#include "hy_list_reference_orderable.cpp"

#endif

