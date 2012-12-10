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

#define  MEMORYSTEP 8

class _SimpleList:public BaseObj
{
    friend class _AVLList;

    protected:
        //memory allocated enough for this many slots
        unsigned long laLength;

    public:
        //Data Members
        long* lData;
        unsigned long lLength;//actual length

        //Methods

        /*
        ==============================================================
        Constructors
        ==============================================================
        */

        //does nothing
        _SimpleList();

        //length constructor
        _SimpleList(unsigned long);

        // stack copy contructor
        _SimpleList(_SimpleList&,long=0,long=-1);

        // data constructor (1 member list)
        _SimpleList(long);

        // arithmetic series populator: size, first item, step
        _SimpleList(long,long,long);

        //destructor
        virtual ~_SimpleList(void);


        /**
         * Data constructor list of longs supplied as a variable
         * @param long the first string to add to the list
         * @param const unsigned long the number of additional long arguments supplied to the constructor
         * @param 2-N: long to be added to the list
         */
        _SimpleList(const long, const unsigned long, ...);

    /*
        ==============================================================
        Operator Overloads
        ==============================================================
        */

        // element location functions - read/write
        long& operator [] (const long);

        // element location functions - read only
        long operator () (const unsigned long);

        // assignment operator
        virtual _SimpleList operator = (_SimpleList);

        // append operator
        virtual _SimpleList operator & (_SimpleList);

        // append number to this
        virtual void operator << (long);

        // append number to this if it's not in the list (search first). List assumed unsorted.
        virtual bool operator >> (long);

        virtual void operator << (_SimpleList&);


        /*
        ==============================================================
        Methods
        ==============================================================
        */

        /**
        * Retrieve the element in position index if index if positive or 
        * length + index if index is negative
        * Example: SimpleList(1,3,5,7).GetElement(1) = 3, SimpleList(1,3,5,7).GetElement(-1) = 7 
        * @param index The index of the elemnt to retrieve 
        * @return the value of the element at the specified index.
        */
        long GetElement (const long index);

        /**
        * Find the position of a search string in the list of strings (ONLY)
        * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
        * @param s The integer to find
        * @param startAt Index to start at 
        * @return -1 if not found, index if found
        */
        virtual long BinaryFind(long, long startAt=0);

        // insert an element into the sorted list preserving the sortedness
        long BinaryInsert(long);

        void Clear(bool = true);

        /* SLKP: 20110209
         * An UGLY hack to automate clearing lists that have pointers to formulas in them 
        */
        void ClearFormulasInList(void);

        /**
        * Compares two parts of the list 
        * Example: SimpleList sl(1,3,5).Compare(0,1) = -1 
        * @param i The index to compare 
        * @param j The second index to compare
        * @return -1 if i<j, 0 if i==j, or 1 if i>j 
        */
        virtual long Compare(long,long);
        virtual long Compare(BaseRef,long);

        long CountCommonElements(_SimpleList&, bool=false);

        /**
        //Lists length
        * Example: SimpleList SimpleList([4, 1, 2]).countitems() = 4 
        * @return Unsigned long of item length    
        */
        unsigned long countitems(void);

        /**
        * SLKP: 20090611    
        * Print the names of variables whose indices are
        * contained in the list
        * @return Nothing. Prints out to screen 
        */
        void DebugVarList(void);

        // delete the item at a given poisiton
        void Delete(long, bool=true);

        virtual void Duplicate(BaseRef);

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
        virtual void DeleteList(const _SimpleList&);

        // shift a range of elements in the array
        void Displace(long,long,long);

        /**
        * Much like [] and () except negative indices return offsets from the end. Invalid indices return 0;
        * Example: SimpleList(1,3,5,7).Element(1) = [3]
        * @param index Which item you want.
        * @return A long    
        */
        long Element(long);

        /**
        * Checks if list is identical to other list 
        * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4 
        * @return true if equal. 
        */
        bool Equal(_SimpleList&);

        /**
        * Retain all those elements that are between (strictly) the 1st and the 2nd argument
        * Example: SimpleList(1,3,5,7).FilterRange(2,4) = [3,5,7] 
        * @param lb Start of new list 
        * @param ub End of new list 
        * @return Nothing. Operates on class that called it. 
        */
        virtual void FilterRange(long, long);

        /**
        * Find the position of a search string in the list of strings (ONLY)
        * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
        * @param s The integer to find
        * @param startAt Index to start at 
        * @return -1 if not found, index if found
        */
        virtual long Find(long, long startAt = 0);

        /**
        * Same as find, but steps over indices 
        * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
        * @param s The integer to find
        * @param step The number to skip between searches 
        * @param startAt Index to start at 
        * @return -1 if not found, index if found
        */
        virtual long FindStepping(long, long, long=0);

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
        virtual void InsertElement(BaseRef br,long insertAt=-1, bool store=true, bool pointer=true);

        void Intersect(_SimpleList&, _SimpleList&);

        BaseRef ListToPartitionString(void);

        virtual BaseRef makeDynamic(void);

        /**
        * SLKP: 20090508
        * Return the maximum value in the list 
        * Example: _SimpleList([4, 1, 2]).Min() = 1 
        * @return maximum value in the list 
        */
        long Max(void);


        /**
        * SLKP: 20090508
        * Return the sum of all values in the list 
        * Example: _SimpleList([4, 1, 2]).Sum() = 7 
        * @return the sum of all values in the list 
        */
        long Sum (void);

        /**
        * Populate a Simple List with integers incrementally.
        * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
        * @param s The substring to find
        * @param startat The index to start searching from
        * @param increment by Pass true for a case sensitive search 
        * @return Nothing. Acts on the List object it was called from. 
        */
        void Merge(_SimpleList& l1, _SimpleList& l2, _SimpleList* mergeResults = nil, _SimpleList* mergeResults2 = nil);

        /**
        * SLKP: 20090508
        * Return the minimum value in the list 
        * Example: _SimpleList([4, 1, 2]).Min() = 1 
        * @return minimum value in the list 
        */
        long Min(void);


        /**
        * Initialize the function to select all k-element subsets of a given simple list
        * @param state  a state-storing simple list; will be approximately the same length as (*this) _SimpleList;
        * DO NOT MANIPULATE this list outside NChooseKInit; it must persist between calls to NChooseK
        * @param store the receptacle list that will store k-tuples
        * @param stride how many elements to choose; must be <= lLength
        * @param algorithm which algorithm to use for k-tuple generation; false - lexicographic (in the sense of the original list order)
        * : true - 'revolving door' method - TBA
        * @return true if successfully initialized
        */
        bool NChooseKInit(_SimpleList&, _SimpleList&, unsigned long, bool = false);

        /**
        * Select the next k-tuple
        * Example: SimpleList(1,3,5,7).DeleteList([0,1,2]) = [7] 
        * @param state the state-storing list previously populated by NChooseKInit
        * @param store the receptacle that will store k-tuples 
        * @return [bool] true is more k-tuples are available; [false] if the last one has just been stored
        */
        bool NChooseK(_SimpleList&, _SimpleList&);

        /**
        * SLKP: 20090316
        * Given a range [from,to] and a given list,
        * make the range conform to the list(e.g. resolve negative to and/or from coordinates)
        * clip the range to fit the list etc
        * Example: SimpleList sl.NormalizeCoordinates(4, 1, 2) = [1, 3, 5, 7]
        * @param from The substring to find
        * @param to The index to start searching from
        * @param refLength The third argument is the length of the list to normalize with respect to.
        * @return Nothing. Acts on the List object it was called from. 
        */
        static void NormalizeCoordinates(long&, long&, const unsigned long);


        /**
        * Add a number to each entry in the array
        * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
        * @param shift Number to add 
        * @return Nothing. Acts on the List object it was called from. 
        */
        void Offset(long);

        /**
        * TODO:Permute elements in blocks of given size
        * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
        * @param shift Number to add 
        * @return Nothing. Acts on the List object it was called from. 
        */
        void Permute(long);
    

        /**
        * TODO:Permute elements in blocks of given size with possible replacement
        * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
        * @param shift Number to add 
        * @return Nothing. Acts on the List object it was called from. 
        */
        void PermuteWithReplacement(long);

    
        /**
         * Select a number of list elements at random (either with or w/o replacement)
         * Example: SimpleList(1,3,5,7).Subset(2) = (1,7)
         * @param size How many elements to select (values >= lLength are rest to lLength)
         * @param select Sample with our without replacement 
         * @return Return the list of sampled elements
         */
        _SimpleList* Subset (unsigned long size, const bool replacement = false);
    
    
        /**
        * Retrive the last value and shorted the list by 1
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
        void Populate(long, long, long); 

        /**
        * TODO
        * Example: SimpleList sl(1,2,3).Flip() = [3,2,1]
        * @return Nothing. Acts on the List object it was called from. 
        */
        void RecursiveIndexSort(long from, long to, _SimpleList* index);


        /**
        * Request space for a given # of elements 
        * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4 
        * @return true if equal. 
        */
        void RequestSpace(long);

        void Subtract(_SimpleList&, _SimpleList&);

        /**
        * Swaps two positions  
        * Example: SimpleList sl[1,3,5].Swap(0, 1) = [3,1]
        * @param i First index to swap 
        * @param j Second index to swap with
        * @return Nothing. Acts on the List object it was called from. 
        */
        void Swap(long, long); //swap two elements

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
        void Union(_SimpleList&, _SimpleList&);

        void XOR(_SimpleList&, _SimpleList&);


        /**
        * Sorts List 
        * Example: SimpleList sl.Sort([5,4,3,2,1]) = [1, 2, 3, 4, 5]
        * @param ascending true if ascending, false for descending sort 
        * @return Nothing. Acts on the List object it was called from. 
        */
        void Sort(bool ascending=true);


        /**
        * SLKP: 20090508
        * Implements a counting sort procedure, ASSUMING that all
        * list values are in [0, upperBound-1]; if the 1st argument is <0, it is automatically
        * determined
        * @return a pointer to the sorted list
        * if the second argument is not nil, then
        * the new_order->old_order mapping is returned in the array pointed to
        *
        */
        _SimpleList* CountingSort(long, _SimpleList* = nil);


        void BubbleSort(void);
        void QuickSort(long, long);

        long* quickArrayAccess(void) {
            return (long*)lData;
        }
};

//TODO:Why is this a global function? If it needs to be, should be in helpers.cpp
void SortLists(_SimpleList*, _SimpleList*);

#endif
