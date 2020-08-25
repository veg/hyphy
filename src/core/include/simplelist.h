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

#ifndef _SIMPLELIST_
#define _SIMPLELIST_
//#pragma once
#include "string.h"
#include "baseobj.h"

#define  MEMORYSTEP 8

class _SimpleList:public BaseObj {
    friend class _AVLList;
    
    public:
    //Data Members
        long* list_data;
        unsigned long lLength;//actual length

    protected:
        //memory allocated enough for this many slots
        unsigned long laLength;
        // default static buffer to hold smaller arrays and avoid dynamic mallocs
        long     static_data [MEMORYSTEP];
    
        void     _UpdateStorageType (void);
        void     _EnsureCorrectStorageType (void);
        void     _CopyStatic (void);

    public:
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

        //length constructor with storage already preallocated somewhere else
        //e.g. using `alloca`
        _SimpleList(unsigned long, long* preallocated_storage);
    
        // stack copy contructor
        _SimpleList(_SimpleList const&,long=0,long=-1);
    
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
        long operator () (const unsigned long) const;

        // assignment operator
        virtual _SimpleList const & operator = (_SimpleList const&);

        // append operator
        virtual _SimpleList operator & (_SimpleList const&);

        // append number to this
        _SimpleList& operator << (long);

        // append number to this if it's not in the list (search first). List assumed unsorted.
        virtual bool operator >> (long);

        virtual void operator << (_SimpleList const &);
  


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
        long GetElement (const long index) const;

        /**
        * Find the position of a search string in the list of strings (ONLY)
        * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
        * @param s The integer to find
        * @param startAt Index to start at 
        * @return -1 if not found, index if found
        */
        virtual long BinaryFind(long, long startAt=0) const;

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
        virtual hyComparisonType Compare(long,long) const;
        virtual hyComparisonType Compare(BaseObj const*,long) const;

        long CountCommonElements(_SimpleList const&, bool=false) const;

        /**
        //Lists length
        * Example: SimpleList SimpleList([4, 1, 2]).countitems() = 4 
        * @return Unsigned long of item length    
        */
        inline unsigned long countitems(void) const {return lLength;}
  
  
        /**
         * Append a range to the current list
         * @param how_many : the number of items in the range
         * @param start: the first element
         * @param step : step from i to i+1 element
         
         */
  
        void   AppendRange (unsigned long how_many, long start, long step);

        /**
         //Is the list empty
         * Example: SimpleList SimpleList([4, 1, 2]).empty() = false
         * @return True if the list is empty
         */
        inline bool empty (void) const {return lLength == 0UL;}

        /**
           Does this list have dynamic memory allocation, or does it have a static array?
         * @return True if the list has dynamically managed memory
         */
        inline bool is_dynamic (void) const {return list_data != static_data;}

        /**
         //Is the list non-empty
         * Example: SimpleList SimpleList([4, 1, 2]).nonempty() = true
         * @return True if the list is non-empty
         */
        inline bool nonempty (void) const {return lLength > 0UL;}

        /**
        * SLKP: 20090611
        * Print the names of variables whose indices are
        * contained in the list
        * @return Nothing. Prints out to screen 
        */
        void DebugVarList(void);

        // delete the item at a given poisiton
        void Delete(long, bool=true);

        virtual void Duplicate(BaseRefConst);

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

        /** shift a range of elements in the array
         * Move a block of elements in the list 
         * Example: SimpleList (1,3,5,7).Displace (1,2,1) = (1,7,3,5)
         * @param from start the block here (0-based, inclusive)
         * @param to end the block here (0-based, inclusive)
         * @param delta number of slots to move (positive or negative, range checked)
         */
        void Displace(long,long,long);
        
        /** Adjust the argument for skipped elements in the [0-max] range
         * so that the argument is increased by however many elements in the list
         * are less than it. The list must be sorted
         * Example: SimpleList (2,4,5).SkipCorrect (2) = 3 
         * Example: SimpleList (2,4,5).SkipCorrect (4) = 7
         
         * @param index the index to correct
         * @return the corrected index
         */
        long SkipCorrect (long) const;

        /** Adjust the argument for skipped elements in the [0-max] range
         * so that the argument is remapped to the range with elements in this
         * list excluded. (*this) list must be sorted
         * Example: SimpleList (2,4,5).SkipCorrect (3) = 2
         
         * @param index the index to correct
         * @param excluded the value to return if index is in the list
         * @return the corrected index
         */
        long CorrectForExclusions (long index, long excluded = -1) const;

        /** Adjust the [sorted] list of indcies argument for skipped elements in the [0-max] range
         * so that the arguments is remapped to the range with elements in this
         * list excluded. (*this) list must be sorted
         * Example: SimpleList (2,4,5).CorrectForExclusions ([2,3],2) = 1 (and the list is now [2])
         
         * @param index the list of indices to correct; corrected indices are written here
         * @param count the length of indices
         * @return the number of entires in the indices list that are not excluded
         */
  
        long CorrectForExclusions (long * indices, long count) const;
  
        /**
        * Much like [] and () except negative indices return offsets from the end. Invalid indices return 0;
        * Example: SimpleList(1,3,5,7).Element(1) = [3]
        * @param index Which item you want.
        * @return A long    
        */
        long Element(long) const;

        /**
         * no range checking element access
         * Example: SimpleList(1,3,5,7).Get (1) = [3]
         * @param index Which item you want.
         * @return A long
         */
        inline long get (long index) const {return list_data [index];}

        /**
        * Checks if list is identical to other list
        * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4 
        * @return true if equal. 
        */
        bool Equal(_SimpleList const&) const;

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
        virtual long Find(long, long startAt = 0) const;

        template <typename FILTER> long FindOnCondition (FILTER condition, long startAt = 0) const {
          for (unsigned long i = startAt; i<lLength; i++) {
            if ( condition (((long*)(list_data))[i], i )) {
              return i;
            }
          }
          return kNotFound;
        }

    
        template <typename MAPPER> void Each (MAPPER&& mapper, long startAt = 0) const {
          for (unsigned long i = startAt; i<lLength; i++) {
            mapper ( ((long*)(list_data))[i], i );
          }
        }

        template <typename MAPPER> bool Any (MAPPER&& mapper, long startAt = 0) const {
            for (unsigned long i = startAt; i<lLength; i++) {
                if (mapper ( ((long*)(list_data))[i], i )) return true;
            }
            return false;
        }

        template <typename CONDITION> bool Every (CONDITION&& predicate, long startAt = 0) const {
            for (unsigned long i = startAt; i<lLength; i++) {
                if (!predicate ( ((long*)(list_data))[i], i )) return false;
            }
            return true;
        }

        template <typename FILTER> _SimpleList const FilterIndex (FILTER condition) const {
          _SimpleList filtered;
          for (unsigned long i = 0UL; i<lLength; i++) {
            if (condition (((long*)(list_data))[i] , i)) {
              filtered << i;
            }
          }
          return filtered;
        }

        template <typename FILTER> _SimpleList const Filter (FILTER condition, unsigned long start_index = 0UL) const {
            _SimpleList filtered;
            for (unsigned long i = start_index; i<lLength; i++) {
                long value = ((long*)(list_data))[i];
                if (condition (value , i)) {
                    filtered << value;
                }
            }
            return filtered;
        }

        template <typename FUNCTOR> _SimpleList const MapList (FUNCTOR transform, unsigned long start_index = 0UL) const {
            _SimpleList mapped;
            mapped.RequestSpace(countitems() - (long)start_index);
            for (unsigned long i = start_index; i<lLength; i++) {
                mapped << transform (((long*)(list_data))[i] , i);
            }
            return mapped;
        }
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

        /**
            Initialize the list object
            @param bool if set to true, make a default memory allocation for the elements
         
            @return nothing; acts on this
         */
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

        BaseRef ListToPartitionString(void) const;

        virtual BaseRef makeDynamic(void) const;

        /**
        * SLKP: 20090508
        * Return the maximum value in the list 
        * Example: _SimpleList([4, 1, 2]).Min() = 1 
        * @return maximum value in the list 
        */
        long Max(void) const;


        /**
        * SLKP: 20090508
        * Return the sum of all values in the list 
        * Example: _SimpleList([4, 1, 2]).Sum() = 7 
        * @return the sum of all values in the list 
        */
        long Sum (void) const;

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
        long Min(void) const;


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
         * TODO: Sample (without replacement)
         * @param size how many elements to sample
         * @return resampled list
         */
        _SimpleList const Sample (unsigned long size) const;

        /**
         * Draw a random element from the list
         * @return the index of the sampled element or -1 if the list is empty
         */
        long Choice () const;

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
        * @param discard; if provided, discard this many items off the end of the list before popping
        * @return Return last value from the list
        */
        long Pop(unsigned long discard = 0UL);

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
         * Populate a SimpleList from another list using a transform
         * @param source start with this list
         * @param action (target array [long* __restrict__], source array [long const* __restrict__], index [unsigned long]) -> ()
         * @return *this. Acts on the _List object it was called from.
         */
    
        template <typename functor> _SimpleList& Populate(_SimpleList const& source, functor action) {
            this->RequestSpace (source.lLength);
            for (unsigned long idx = 0UL; idx < source.lLength; idx++) {
                action (this->list_data, source.list_data, idx);
            }
            this->lLength = source.lLength;
            return *this;
        }
    
        /**
         * Some functors for the populator above
         * Assuming that the list on N elements stores a 1-1 map from [0,N-1] to [0, N-1], invert it
         * so that if source (x) = y becomes target (y) = x
         */
        static void action_invert (long * __restrict__ target, long * __restrict__ source, unsigned long idx) {
            target [source[idx]] = idx;
        }
        /**
         * Some functors for the populator above
         * so that if target [x] -> source (target (x))
         */
        static void action_compose (long * __restrict__ target, long * __restrict__ source, unsigned long idx) {
            target [idx] = source [target[idx]];
        }

        /**
        * TODO
        * Example: SimpleList sl(1,2,3).Flip() = [3,2,1]
        * @return Nothing. Acts on the List object it was called from. 
        */
        void RecursiveIndexSort(long from, long to, _SimpleList* index);


        /**
        * Request space for a given # of elements 
        */
        void RequestSpace(long);
    
        virtual bool is_numeric_list (void) const { return true;}

        void Subtract(_SimpleList const &, _SimpleList const&);

        /**
        * Swaps two positions  
        * Example: SimpleList sl[1,3,5].Swap(0, 1) = [3,1]
        * @param i First index to swap 
        * @param j Second index to swap with
        * @return Nothing. Acts on the List object it was called from. 
        */
        void Swap(long, long); //swap two elements

        virtual BaseRef toStr(unsigned long = 0UL);

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
         * Maps the integer argument to the value of the corresponding 
         * element in the list. Performs range checks, and returns default value 
         * if the mapping failed. Equivalent to 
         * return index >= 0 && index < lLength ? lData[index] : map_failed
         
         * @param index the index to map
         * @param map_failed the value to return if the index is out of range
         * @return the mapped value
         */
  
         long Map (long index, long map_failed = -1L) const;
 


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
    
        _SimpleList* CountingSort(long, _SimpleList* = nil, bool wantResult = true);

        void BubbleSort(void);
        void QuickSort(long, long);

        long* quickArrayAccess(void) {
            return (long*)list_data;
        }
};

//TODO:Why is this a global function? If it needs to be, should be in helpers.cpp
void SortLists(_SimpleList*, _SimpleList*);

#endif
