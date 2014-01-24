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

#ifndef _HY_LIST_ORDERABLE_
#define _HY_LIST_ORDERABLE_
//#pragma once
#include "hy_list.h"

#define HY_COMPARE_EQUAL 0L
#define HY_COMPARE_LESS  (-1L)
#define HY_COMPARE_GREATER  1L



/*
  
  A resizable list which stores PAYLOAD types
  
  assumes that PAYLOAD understands the equal (==), assign (=), 
  and stack copy operations, as well as all ordering operations 
  (<,>,<=,>=)  

*/ 


template <typename PAYLOAD>
class _hyListOrderable : public _hyList<PAYLOAD> {
  public:
  
    //does nothing
  _hyListOrderable();
  
    //_hyList constructor
  _hyListOrderable(const unsigned long);
  
    // stack copy contructor
  _hyListOrderable(const _hyListOrderable <PAYLOAD>&, const long = 0UL, const long = HY_LIST_INSERT_AT_END);
  
    // data constructor (1 member list)
  _hyListOrderable(const PAYLOAD);
  
  
  /**
   * Data constructor list of longs supplied as a variable
   * @param long the first string to add to the list
   * @param const unsigned long the number of additional long arguments supplied
   * to the constructor
   * @param 2-N: long to be added to the list
   */
  _hyListOrderable(const PAYLOAD, const unsigned long, ...);

  /**
    * Find the position of an item
    * Example: _hyListOrderable(1,3,5,7).Find(5,1) = 2
    * @param s The value to find
    * @param startAt Index to start at
    * @return HY_NOT_FOUND if not found, index if found
    */
    virtual long BinaryFind(const PAYLOAD, const long startAt = 0) const;

    /**
    * Insert the item into a sorted list, preserving sortedness
    * Example: _hyListOrderable(1,3,5,7).Insert(4) = [1,3,4,5,7] (returns 2)
    * @param item The value to find
    * @return the index (0-based) where the item was inserted
    */
    long BinaryInsert(const PAYLOAD);


    /**
    * Counts the number of elements in common between [SORTED] *this and another [SORTED] argument
    * Example: _hyListOrderable(1,3,5,7).CountCommonElements([2,5,7,8]) = 2
    * @param l1 The second sorted list
    * @param at_least_one if true, returns 0/1 if the interesction of the lists is non-empty
    * @return the number of elements in common or 0/1 is at_least_one is true
    */
    long CountCommonElements(const _hyListOrderable<PAYLOAD> &, bool at_least_one = false) const;


    /**
    * Delete list of duplicate items in a sorted list
    * Example: _hyListOrderable(1,3,3,7).DeleteDuplicates() = [1,3,7]
    * @param toDelete SimpleList of indices to
    * @return Nothing. Acts on the  object it was called from.
    */
    void DeleteDuplicates(void);


    /**
    * Compares two elements of the list
    * Example: _hyListOrderable sl(1,3,5).Compare(0,1) = -1
    * @param i The index to compare
    * @param j The second index to compare
    * @return -1 if i<j, 0 if i==j, or 1 if i>j
    */
    virtual long Compare(const long, const long) const;


    /**
     * Compares the element at a given index to a fixed value
     * Example: _hyListOrderable sl(1,3,5).Compare(0,1) = -1
     * @param item The index of the item to compare
     * @param value The value to compare the index to
     * @return -1 if i<j, 0 if i==j, or 1 if i>j
     */
    virtual long CompareToValue (const long item, const PAYLOAD& value) const;

  /**
    * Retain all those elements that are between (strictly) the 1st and the 2nd
    * argument (lb < value < ub)
    * Example: _hyListOrderable(1,3,5,7).FilterRange(2,5) = [3]
    * @param lb Lower Bound
    * @param ub Upper Bound
    * @return Nothing. Operates on class that called it.
    */
    virtual void FilterRange(const PAYLOAD, const PAYLOAD);

    /**
    * Return the minimum value in the list
    * Example: _hyListOrderable([4, 1, 2]).Min() = 1
    * @return minimum value in the list; undefined on empty sets
    */
    PAYLOAD Min(void) const;

    /**
    * Return the maximum value in the list
    * Example: _hyListOrderable([4, 1, 2]).Max() = 4
    * @return maximum value in the list; undefined on empty sets
    */
    PAYLOAD Max(void) const;
    
    /**
    * Computes the union of two SORTED orderable lists, store the result in *this
    * Example: _hyListOrderable.Union ([1,4,6], [1,2,4,7]) = [1,2,4,6,7]
    * @param l1 the first list
    * @param l2 the second list
    * @return Nothing. Acts on the List object it was called from.
    */
    void Union(const _hyListOrderable <PAYLOAD> &, const _hyListOrderable <PAYLOAD> &);

    /**
    * Computes the symmetric difference of two SORTED orderable lists, store the result in *this
    * Example: _hyListOrderable.XOR ([1,4,6], [1,2,4,7]) = [2,6,7]
    * @param l1 the first list
    * @param l2 the second list
    * @return Nothing. Acts on the List object it was called from.
    */
    void XOR(const _hyListOrderable <PAYLOAD> &, const _hyListOrderable <PAYLOAD> &);

    /**
    * Computes the intersection of two SORTED orderable lists, store the result in *this
    * Example: _hyListOrderable.Intersect ([1,4,6], [1,2,4,7]) = [1,4]
    * @param l1 the first list
    * @param l2 the second list
    * @return Nothing. Acts on the List object it was called from.
    */
    void Intersect (const _hyListOrderable <PAYLOAD> &, const _hyListOrderable <PAYLOAD> &);


    /**
    * Computes the difference of two SORTED orderable lists, store the result in *this
    * Example: _hyListOrderable.Subtract ([1,4,6], [1,2,4,7]) = [6]
    * @param l1 the first list
    * @param l2 the second list
    * @return Nothing. Acts on the List object it was called from.
    */
    void Subtract (const _hyListOrderable <PAYLOAD> &, const _hyListOrderable <PAYLOAD> &);

    /**
    * Runs the (quadratic) bubble-sort algoritm
    * If provided to the function, the list passed as the argument will have its elements ordered the same way
    * as the ordered list being sorted (this is useful to keep track of the sorted order of element indices)
    * @param index_list
    */
    void BubbleSort (_hyListOrderable<long>* = nil);

    /**
    * Runs the (N log(N)) quick-sort algoritm on a range of list elements
    * If provided to the function, the list passed as the third argument will have its elements ordered the same way
    * as the ordered list being sorted (this is useful to keep track of the sorted order of element indices)
    * Recursively calls self
    * @param from start the sort at this index
    * @param to end the sort at this index
    * @param index_list
    */
    void QuickSort (const unsigned long, const unsigned long, _hyListOrderable<long>* = nil);
    
    /**
    * Sorts List
    * Example:  Sort([5,4,3,2,1]) = [1, 2, 3, 4, 5]
    * If provided to the function, the list passed as the third argument will have its elements ordered the same way
    * as the ordered list being sorted (this is useful to keep track of the sorted order of element indices)
    * @param ascending true if ascending, false for descending sort
    * @param index_list
    * @return Nothing. Acts on the List object it was called from.
    */
    void Sort(bool ascending = true, _hyListOrderable<long>* = nil);
    
    /**
    * Merge two sorted lists into *this
    * Example:  Sort([5,4,3,2,1]) = [1, 2, 3, 4, 5]
    * If provided to the function, the list passed as the third argument will have its elements ordered the same way
    * as the ordered list being sorted (this is useful to keep track of the sorted order of element indices)
    * @param l1 first list to merge
    * @param l2 second list to merge
    * @param mergeResults1 if not NULL then store the indices of l1 entries in the merged array
    * @param mergeResults2 if not NULL then store the indices of l2 entries in the merged array
    */
    void Merge(const _hyListOrderable <PAYLOAD> &l1, const _hyListOrderable <PAYLOAD> &l2, _hyListOrderable <long> *mergeResults1 = nil,
             _hyListOrderable <long> *mergeResults2 = nil);

};

#include "hy_list_orderable.cpp"

#endif
