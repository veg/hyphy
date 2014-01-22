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

#ifndef _HY_LIST_
#define _HY_LIST_
//#pragma once
#include "baseobj.h"
#include "hy_string_buffer.h"

/*
  
  A resizable list which stores PAYLOAD types
  
  assumes that PAYLOAD understands the equal (==), assign (=), 
  and stack copy operations  

*/ 

#define HY_LIST_ALLOCATION_CHUNK 8UL
#define HY_LIST_INSERT_AT_END    (-1L)

template <typename PAYLOAD>
class _hyList : public virtual BaseObj {

protected:
  //memory allocated enough for this many slots
  unsigned long laLength;
  PAYLOAD       *lData;
  unsigned long lLength; //actual length
  
  void          CompactList (void);
  void          ResizeList  (void);

  //Methods

  /*
  ==============================================================
  Constructors
  ==============================================================
  */

public:

  //does nothing
  _hyList();

  //_hyList constructor
  _hyList(const unsigned long);

  // stack copy contructor
  _hyList(const _hyList <PAYLOAD>&, const long = 0UL, const long = HY_LIST_INSERT_AT_END);

  // data constructor (1 member list)
  _hyList(const PAYLOAD);

  //destructor
  virtual ~_hyList(void);

  /**
   * Data constructor list of longs supplied as a variable
   * @param long the first string to add to the list
   * @param const unsigned long the number of additional long arguments supplied
   * to the constructor
   * @param 2-N: long to be added to the list
   */
  _hyList(const PAYLOAD, const unsigned long, ...);

  /*
      ==============================================================
      Operator Overloads
      ==============================================================
  */

  // element location functions - read/write
  PAYLOAD &operator[](const long);

  // element location functions - read only
  PAYLOAD operator()(const unsigned long);

  // assignment operator
  const _hyList<PAYLOAD> operator=(const _hyList<PAYLOAD>&);

  // append operator
  const _hyList<PAYLOAD> operator&(const _hyList<PAYLOAD>);

  // append an instance to this
  
  virtual void append    (const PAYLOAD);
  
  virtual void operator<<(const PAYLOAD);

  // append number to this if it's not in the list (search first). List assumed
  // unsorted.
  virtual bool operator>>(const PAYLOAD);

  virtual void operator<<(const _hyList <PAYLOAD>&);
  
  

  /*
  ==============================================================
  Methods
  ==============================================================
  */

  /**
   * Clear the current list and make a copy from the argment
   * Argument is a pointer to make it possible to overload Clone
   * @param clone_from the object to clone from
   * @return None.
   */

  virtual void Clone (const _hyList<PAYLOAD>* clone_from);
  
  /**
  * Retrieve the element in position index if index if positive or
  * length + index if index is negative
  * Example: _hyList(1,3,5,7).GetElement(1) = 3,
  * _hyList(1,3,5,7).GetElement(-1) = 7
  * @param index The index of the elemnt to retrieve
  * @return the value of the element at the specified index.
  */
  inline PAYLOAD GetElement(const long index) const;

  /**
  * Clear the list (set lLength to 0)
  * @param compact Free allocated memory if true
  */

  virtual void          Clear(bool = true);
  
  /**
   * Element storage functions - assuming the space is already allocated
   * used to avoid (*list)[3] = x which are hard to read
   * @param index where to store the item (not range-checked)
   * @param item the stuff to store
   */
  virtual void SetItem(const unsigned long, PAYLOAD);


  /**
  * Return the number of elements in the list
  * @return list length
  */

  inline unsigned long countitems (void) const;


  /**
  * delete the item at a given poisiton
  * @param index the index (0-based) of the element to delete
  * @param compact_list Free allocated memory if true
  */  
  void         Delete(const long, bool compact_list = true);

  virtual void    Duplicate  (BaseRefConst);
  virtual BaseRef makeDynamic(void);

  /**
  * Delete a sorted list of indices in a sorted list
  * Example: _hyList(1,3,5,7).DeleteList([0,1,2]) = [7]
  * @param toDelete _hyList of indices to delete
  * @return Nothing. Acts on the List object it was called from.
  */
  virtual void DeleteList(const _hyList<long> *);

  /**
  * Shift a range of elements in the array
  * Example: _hyList(1,3,5,7).Displace([1,2,1]) = [1,7,3,5];  _hyList(1,3,5,7).Displace([1,2,-1]) = [3,5,1,7]
  * @param from start shifting at this index
  * @param end end shifting at this index (inclusive)
  * @param delta offset the elements by this much
  * @return Nothing. Acts on the List object it was called from.
  */
  void Displace(long, long, long);

  /**
  * Much like [] and () except negative indices return offsets from the end.
  * Note: Invalid indices trigger an error message
  * Example: _hyList(1,3,5,7).Element(-1) = [7]
  * @param index Which item you want.
  * @return A long
  */
  
  virtual PAYLOAD Element(const long);

  /**
  * Checks if list is identical to other list
  * Example: _hyList([4, 1, 2]).Equal(_hyList([4, 1, 2]) = true
  * @return true if equal.
  */
  bool Equal(const _hyList <PAYLOAD>&) const;


  /**
   * Checks if list if list element at 'index' is equal to a fixed value
   * @param index the list index of the element to test
   * @param value the value to compare the result to
   * @return true if equal.
   */
  virtual bool ItemEqualToValue (unsigned long index, const _hyList <PAYLOAD>& value) const;

  /**
  * Find the position of an item in an unsorted list using linear search
  * Example: _hyList(1,3,5,7).Find(3,3) = -1
  * @param s The integer to find
  * @param startAt Index to start at
  * @return HY_NOT_FOUND if not found, index if found
  */
  virtual long Find(const PAYLOAD, const long startAt = 0L) const;

  /**
  * Same as find, but steps over indices
  * Example: _hyList(1,3,5,7).Find(3,3) = -1
  * @param s The integer to find
  * @param step The number to skip between searches
  * @param startAt Index to start at
  * @return -1 if not found, index if found
  */
  virtual long FindStepping(const PAYLOAD, const long, const long = 0L) const;

  /**
  * Flips the list
  * Example: _hyList sl(1,2,3).Flip() = [3,2,1]
  * @return Nothing. Acts on the List object it was called from.
  */
  void Flip(void); //flip the order of list elements

  /**
  * Initializes the list
  * @param doMemAlloc If true, perform default memory allocation
  * @return Nothing. Acts on the List object it was called from. */
  virtual void Initialize(bool = true);

  /**
  * Insert an element at a specific point
  * Example: _hyList sl.Populate(4, 1, 2).InsertElement(1,1,?,false)
  * @param item The variable to insert
  * @param insert_at The position to insert into
  * @return Nothing. Acts on the List object it was called from.
  */
  virtual void InsertElement(const PAYLOAD, const long insertAt = -1);



  /**
  * Permute elements in blocks of given size
  * Example: _hyList(1,3,5,7).Permute(2) = [5, 7, 1, 3] // random output
  * @param block_size Treat groups with this many elements as units to base the permuations on
  * @return Nothing. Acts on the List object it was called from.
  */
  void Permute(const unsigned long);

  /**
  * Permute elements in blocks of given size with replacement
  * Example: _hyList(1,3,5,7).PermuteWithReplacement(2) = [5, 7, 5, 7] // random output
  * @param block_size Treat groups with this many elements as units to base the permuations on
  * @return Nothing. Acts on the List object it was called from.
  */
  void PermuteWithReplacement(const unsigned long);

  /**
   * Select a number of list elements at random (either with or w/o replacement)
   * Example: _hyList(1,3,5,7).Subset(2) = (1,7)
   * @param size How many elements to select (values >= lLength are rest to
   * lLength)
   * @param replacement Sample with our without replacement
   * @return Return the list of sampled elements
   */
  _hyList <PAYLOAD> *Subset(unsigned long size, const bool replacement = false);

  /**
  * Retrive the last value and shorten the list by 1
  * Example: _hyList(1,3,5,7).Pop() = 7
  * @return Return last value from the list
  */
  virtual PAYLOAD Pop();

  void Subtract(_hyList <PAYLOAD> &, _hyList <PAYLOAD> &);

  /**
  * Swaps two positions
  * Example: _hyList sl[1,3,5].Swap(0, 1) = [3,1]
  * @param i First index to swap
  * @param j Second index to swap with
  * @return Nothing. Acts on the List object it was called from.
  */
  inline void    Swap    (const unsigned long, const unsigned long); //swap two elements
  virtual BaseRef toStr(void);

  /**
  * Request space for this many elements
  * @param slots The number of elements to allocate the space for
  */
  
  void    RequestSpace(const unsigned long);
  /**
  * Trim excess memory from the list allocation 
  */
  void    TrimMemory  (void);

  /**
  * Initialize the function to select all k-element subsets of a given simple
  * list
  * @param state  a state-storing simple list; will be approximately the same
  * length as (*this) _SimpleList;
  * DO NOT MANIPULATE this list outside NChooseKInit; it must persist between
  * calls to NChooseK
  * @param store the receptacle list that will store k-tuples
  * @param stride how many elements to choose; must be <= lLength
  * @param algorithm which algorithm to use for k-tuple generation; false -
  * lexicographic (in the sense of the original list order)
  * : true - 'revolving door' method - TBA
  * @return true if successfully initialized
  */
  bool NChooseKInit(_hyList <long> *, _hyList <PAYLOAD>&, const unsigned long, bool = false);

  /**
  * Select the next k-tuple
  * Example: SimpleList(1,3,5,7).DeleteList([0,1,2]) = [7]
  * @param state the state-storing list previously populated by NChooseKInit
  * @param store the receptacle that will store k-tuples
  * @return [bool] true is more k-tuples are available; [false] if the last one
  * has just been stored
  */
  bool NChooseK(_hyList <long> *, _hyList <PAYLOAD>&);

};

#include "hy_list.cpp"

#endif
