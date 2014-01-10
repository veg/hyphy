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

#ifndef _SIMPLELIST_
#define _SIMPLELIST_
//#pragma once
#include "string.h"
#include "baseobj.h"


/*
  
  A resizable list which stores _SIMPLE_ types
  which require no custom constructor or descructor action, 
  supports comparison operators and '+'
*/ 

class _String;

#define HY_LIST_ALLOCATION_CHUNK 8UL
#define HY_LIST_INSERT_AT_END    (-1L)

template <typename PAYLOAD>
class _SimpleList : public virtual BaseObj {

protected:
  //memory allocated enough for this many slots
  unsigned long laLength;
  PAYLOAD       *lData;
  unsigned long lLength; //actual length
  
  void      CompactList (void);
  void      ResizeList  (void);

  //Methods

  /*
  ==============================================================
  Constructors
  ==============================================================
  */

public:

  //does nothing
  _SimpleList();

  //length constructor
  _SimpleList(const unsigned long);

  // stack copy contructor
  _SimpleList(const _SimpleList <PAYLOAD>&, const long = 0, const long = HY_LIST_INSERT_AT_END);

  // data constructor (1 member list)
  _SimpleList(const PAYLOAD);

  // arithmetic series populator: size, first item, step
  _SimpleList(const unsigned long, const PAYLOAD, const PAYLOAD);

  //destructor
  virtual ~_SimpleList(void);

  /**
   * Data constructor list of longs supplied as a variable
   * @param long the first string to add to the list
   * @param const unsigned long the number of additional long arguments supplied
   * to the constructor
   * @param 2-N: long to be added to the list
   */
  _SimpleList(const PAYLOAD, const unsigned long, ...);

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
  virtual const _SimpleList<PAYLOAD> operator=(const _SimpleList<PAYLOAD>);

  // append operator
  virtual const _SimpleList<PAYLOAD> operator&(const _SimpleList<PAYLOAD>);

  // append an instance to this
  
  virtual void append    (const PAYLOAD);
  
  virtual void operator<<(const PAYLOAD);

  // append number to this if it's not in the list (search first). List assumed
  // unsorted.
  virtual bool operator>>(const PAYLOAD);

  virtual void operator<<(const _SimpleList <PAYLOAD>&);

  /*
  ==============================================================
  Methods
  ==============================================================
  */

  /**
  * Retrieve the element in position index if index if positive or
  * length + index if index is negative
  * Example: SimpleList(1,3,5,7).GetElement(1) = 3,
  * SimpleList(1,3,5,7).GetElement(-1) = 7
  * @param index The index of the elemnt to retrieve
  * @return the value of the element at the specified index.
  */
  inline PAYLOAD GetElement(const long index) const;

  /**
  * Find the position of a search string in the list of strings (ONLY)
  * Example: SimpleList(1,3,5,7).Find(3,3) = -1
  * @param s The integer to find
  * @param startAt Index to start at
  * @return -1 if not found, index if found
  */
  virtual long BinaryFind(const PAYLOAD, const long startAt = 0) const;

  // insert an element into the sorted list preserving the sortedness
  long BinaryInsert(const PAYLOAD);

  void Clear(bool = true);


  /**
  * Compares two elements of the list
  * Example: SimpleList sl(1,3,5).Compare(0,1) = -1
  * @param i The index to compare
  * @param j The second index to compare
  * @return -1 if i<j, 0 if i==j, or 1 if i>j
  */
  virtual long Compare(const long, const long);

  long CountCommonElements(const _SimpleList<PAYLOAD> &, bool at_least_one = false) const;

  /**
  //Lists length
  * Example: SimpleList SimpleList([4, 1, 2]).countitems() = 4
  * @return Unsigned long of item length
  */
  unsigned long countitems (void) const;


  // delete the item at a given poisiton
  void Delete(const long, bool compact_list = true);

  virtual void Duplicate  (BaseRef);

  /**
  * Delete all duplicates in a sorted list
  * Example: SimpleList(1,3,3,5,7).DeleteDuplicates() = [1,3,5,7]
  * @return Nothing. Acts on the List object it was called from.
  */
  void DeleteDuplicates(void);

  /**
  * Delete list of indices in a sorted list
  * Example: SimpleList(1,3,5,7).DeleteList([0,1,2]) = [7]
  * @param toDelete SimpleList of indices to
  * @return Nothing. Acts on the List object it was called from.
  */
  virtual void DeleteList(const _SimpleList<long> &);

  // shift a range of elements in the array
  void Displace(long, long, long);

  /**
  * Much like [] and () except negative indices return offsets from the end.
  * Invalid indices return 0;
  * Example: SimpleList(1,3,5,7).Element(1) = [3]
  * @param index Which item you want.
  * @return A long
  */
  PAYLOAD Element(const long);

  /**
  * Checks if list is identical to other list
  * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4
  * @return true if equal.
  */
  bool Equal(const _SimpleList <PAYLOAD>&);

  /**
  * Retain all those elements that are between (strictly) the 1st and the 2nd
  * argument
  * Example: SimpleList(1,3,5,7).FilterRange(2,4) = [3,5,7]
  * @param lb Start of new list
  * @param ub End of new list
  * @return Nothing. Operates on class that called it.
  */
  virtual void FilterRange(const PAYLOAD, const PAYLOAD);

  /**
  * Find the position of a search string in the list of strings (ONLY)
  * Example: SimpleList(1,3,5,7).Find(3,3) = -1
  * @param s The integer to find
  * @param startAt Index to start at
  * @return -1 if not found, index if found
  */
  virtual long Find(const PAYLOAD, const long startAt = 0L) const;

  /**
  * Same as find, but steps over indices
  * Example: SimpleList(1,3,5,7).Find(3,3) = -1
  * @param s The integer to find
  * @param step The number to skip between searches
  * @param startAt Index to start at
  * @return -1 if not found, index if found
  */
  virtual long FindStepping(const PAYLOAD, const long, const long = 0L) const;

  /**
  * Flips the list
  * Example: SimpleList sl(1,2,3).Flip() = [3,2,1]
  * @return Nothing. Acts on the List object it was called from.
  */
  void Flip(void); //flip the order of list elements

  virtual void Initialize(bool = true);

  /**
  * Insert an element at a specific point
  * Example: SimpleList sl.Populate(4, 1, 2).InsertElement(1,1,?,false)
  * @param br The variable to insert
  * @param insertAt The position to insert into
  * @param store
  * @param pointer
  * @return Nothing. Acts on the List object it was called from.
  */
  virtual void InsertElement(const PAYLOAD, const long insertAt = -1);

  void Intersect(_SimpleList <PAYLOAD>&, _SimpleList<PAYLOAD> &);

    // only specialized for long
  _String* ListToPartitionString(void);

  virtual BaseRef makeDynamic(void);

  /**
  * SLKP: 20090508
  * Return the maximum value in the list
  * Example: _SimpleList([4, 1, 2]).Min() = 1
  * @return maximum value in the list
  */
  PAYLOAD Max(void) const;

  /**
  * SLKP: 20090508
  * Return the sum of all values in the list
  * Example: _SimpleList([4, 1, 2]).Sum() = 7
  * @return the sum of all values in the list
  */
  PAYLOAD Sum(void) const;

  /**
  * Populate a Simple List with integers incrementally.
  * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
  * @param s The substring to find
  * @param startat The index to start searching from
  * @param increment by Pass true for a case sensitive search
  * @return Nothing. Acts on the List object it was called from.
  */
  void Merge(_SimpleList <PAYLOAD> &l1, _SimpleList <PAYLOAD> &l2, _SimpleList <long> *mergeResults = nil,
             _SimpleList <long> *mergeResults2 = nil);

  /**
  * SLKP: 20090508
  * Return the minimum value in the list
  * Example: _SimpleList([4, 1, 2]).Min() = 1
  * @return minimum value in the list
  */
  PAYLOAD Min(void) const;

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
  bool NChooseKInit(_SimpleList <PAYLOAD> &, _SimpleList <PAYLOAD>&, unsigned long, bool = false);

  /**
  * Select the next k-tuple
  * Example: SimpleList(1,3,5,7).DeleteList([0,1,2]) = [7]
  * @param state the state-storing list previously populated by NChooseKInit
  * @param store the receptacle that will store k-tuples
  * @return [bool] true is more k-tuples are available; [false] if the last one
  * has just been stored
  */
  bool NChooseK(_SimpleList <PAYLOAD> &, _SimpleList <PAYLOAD> &);

  /**
  * SLKP: 20090316
  * Given a range [from,to] and a given list,
  * make the range conform to the list(e.g. resolve negative to and/or from
  * coordinates)
  * clip the range to fit the list etc
  * Example: SimpleList sl.NormalizeCoordinates(4, 1, 2) = [1, 3, 5, 7]
  * @param from The substring to find
  * @param to The index to start searching from
  * @param refLength The third argument is the length of the list to normalize
  * with respect to.
  * @return Nothing. Acts on the List object it was called from.
  */
  static void NormalizeCoordinates(long &, long &, const unsigned long);

  /**
  * Add a number to each entry in the array
  * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
  * @param shift Number to add
  * @return Nothing. Acts on the List object it was called from.
  */
  void Offset(PAYLOAD);

  /**
  * TODO:Permute elements in blocks of given size
  * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
  * @param shift Number to add
  * @return Nothing. Acts on the List object it was called from.
  */
  void Permute(const long);

  /**
  * TODO:Permute elements in blocks of given size with possible replacement
  * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
  * @param shift Number to add
  * @return Nothing. Acts on the List object it was called from.
  */
  void PermuteWithReplacement(const long);

  /**
   * Select a number of list elements at random (either with or w/o replacement)
   * Example: SimpleList(1,3,5,7).Subset(2) = (1,7)
   * @param size How many elements to select (values >= lLength are rest to
   * lLength)
   * @param select Sample with our without replacement
   * @return Return the list of sampled elements
   */
  _SimpleList <PAYLOAD> *Subset(unsigned long size, const bool replacement = false);

  /**
  * Retrive the last value and shorten the list by 1
  * Example: SimpleList(1,3,5,7).Pop() = 7
  * @return Return last value from the list
  */
  long Pop();

  /**
  * Populate a Simple List with integers incrementally.
  * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
  * @param s The substring to find
  * @param startat The index to start searching from
  * @param increment by Pass true for a case sensitive search
  * @return Nothing. Acts on the List object it was called from.
  */
  void Populate(const unsigned long, const PAYLOAD, const PAYLOAD);

  /**
  * TODO
  * Example: SimpleList sl(1,2,3).Flip() = [3,2,1]
  * @return Nothing. Acts on the List object it was called from.
  */
  void RecursiveIndexSort(long from, long to, _SimpleList <long> *index);

  /**
  * Request space for a given # of elements
  * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4
  * @return true if equal.
  */
  void RequestSpace(const long);

  void Subtract(_SimpleList <PAYLOAD> &, _SimpleList <PAYLOAD> &);

  /**
  * Swaps two positions
  * Example: SimpleList sl[1,3,5].Swap(0, 1) = [3,1]
  * @param i First index to swap
  * @param j Second index to swap with
  * @return Nothing. Acts on the List object it was called from.
  */
  void Swap(const long, const long); //swap two elements

  virtual BaseRef toStr(void);

  /**
  *
  */
  void TrimMemory(void);

  /**
  * Performs union of two SimpleLists
  * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
  * @param shift Number to add
  * @return Nothing. Acts on the List object it was called from.
  */
  void Union(_SimpleList <PAYLOAD> &, _SimpleList <PAYLOAD> &);

  void XOR(_SimpleList <PAYLOAD> &, _SimpleList <PAYLOAD> &);

  /**
  * Sorts List
  * Example: SimpleList sl.Sort([5,4,3,2,1]) = [1, 2, 3, 4, 5]
  * @param ascending true if ascending, false for descending sort
  * @return Nothing. Acts on the List object it was called from.
  */
  void Sort(bool ascending = true);

  /**
  * SLKP: 20090508
  * Implements a counting sort procedure, ASSUMING that all
  * list values are in [0, upperBound-1]; if the 1st argument is <0, it is
  * automatically
  * determined
  * @return a pointer to the sorted list
  * if the second argument is not nil, then
  * the new_order->old_order mapping is returned in the array pointed to
  *
  */
  
    // this only specialized for 'long'
  
  _SimpleList <PAYLOAD> *CountingSort(PAYLOAD, _SimpleList <long> * = nil);

  void BubbleSort(void);
  void QuickSort(long, long);

  long *quickArrayAccess(void) { return (long *)lData; }
};

//TODO:Why is this a global function? If it needs to be, should be in
//helpers.cpp
template<typename PAYLOAD>
void SortLists(_SimpleList<PAYLOAD> *, _SimpleList <PAYLOAD> *);

#include "simplelist.cpp"

#endif
