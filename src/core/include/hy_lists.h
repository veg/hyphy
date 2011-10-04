/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2011
Primary Development:
  Sergei L Kosakovsky Pond (sergeilkp@mac.com)
Significant contributions from:
  Spencer V Muse (muse@stat.ncsu.edu)
  Simon DW Frost (sdfrost@ucsd.edu)
  Art FY Poon    (apoon@biomail.ucsd.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#ifndef _HLIST_
#define _HLIST_
//#pragma once
#include "string.h"
#include "baseobj.h"

#define  MEMORYSTEP 8

// for storing longs

//_____________________________________________________________________________

class _SimpleList:public BaseObj
{

    // contructor/destructor methods
public:

    _SimpleList                     ();
    //does nothing
    _SimpleList                     (unsigned long);
    //length constructor
    _SimpleList                     (_SimpleList&, long = 0, long = -1);
    // stack copy contructor
    _SimpleList                     (long);
    // data constructor (1 member list)

    _SimpleList                     (long, long, long);
    // arithmetic series populator: size, first item, step

    virtual     ~_SimpleList                    (void);
    //destructor

    long& operator []               (long);
    // element location functions - read/write

    long operator ()                (unsigned long);
    // element location functions - read only

    /**
    * Checks if list is identical to other list 
    * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4 
    * @return true if equal. 
    */
    bool        Equal                           (_SimpleList&);

    void        TrimMemory                      (void);

    // 
    /**
    * Request space for a given # of elements 
    * Example: _SimpleList([4, 1, 2]).Equal(_SimpleList([4, 1, 2]) = 4 
    * @return true if equal. 
    */
    void        RequestSpace                    (long);

    virtual     _SimpleList operator =          (_SimpleList);
    // assignment operator


    /**
    //Lists length
    * Example: SimpleList SimpleList([4, 1, 2]).countitems() = 4 
    * @return Unsigned long of item length    
    */
    unsigned long   countitems      (void);

    virtual     _SimpleList operator &          (_SimpleList);
    // append operator

    virtual     void operator <<                (long);
    // append number to this

    virtual     bool operator >>                (long);
    // append number to this if it's not in the list (search first). List assumed unsorted.

    virtual     void operator <<                (_SimpleList&);

    /**
    * Insert an element at a specific point
    * Example: SimpleList sl.Populate(4, 1, 2).InsertElement(1,1,?,false)
    * @param br The variable to insert 
    * @param insertAt The position to insert into
    * @param store 
    * @param pointer 
    * @return Nothing. Acts on the List object it was called from. 
    */
    virtual     void InsertElement              (BaseRef br, long insertAt = -1, bool store = true, bool pointer = true);

    void Clear                      (bool = true);

    /**
    * Much like [] and () except negative indices return offsets from the end. Invalid indices return 0;
    * Example: SimpleList(1,3,5,7).Element(1) = [3]
    * @param index Which item you want.
    * @return A long    
    */
    long Element                    (long);

    /**
    * Retrive the last value and shorted the list by 1
    * Example: SimpleList(1,3,5,7).Pop() = 7 
    * @return Return last value from the list
    */
    long Pop                        ();

    /**
    * Find the position of a search string in the list of strings (ONLY)
    * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
    * @param s The integer to find
    * @param startAt Index to start at 
    * @return -1 if not found, index if found
    */
    virtual     long  Find                      (long, long startAt = 0);

    // 

    /**
    * Retain all those elements that are between (strictly) the 1st and the 2nd argument
    * Example: SimpleList(1,3,5,7).FilterRange(2,4) = [3,5,7] 
    * @param lb Start of new list 
    * @param ub End of new list 
    * @return Nothing. Operates on class that called it. 
    */
    virtual     void  FilterRange               (long, long);

    /**
    * Same as find, but steps over indices 
    * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
    * @param s The integer to find
    * @param step The number to skip between searches 
    * @param startAt Index to start at 
    * @return -1 if not found, index if found
    */
    virtual     long  FindStepping              (long, long, long = 0);


    /**
    * Find the position of a search string in the list of strings (ONLY)
    * Example: SimpleList(1,3,5,7).Find(3,3) = -1 
    * @param s The integer to find
    * @param startAt Index to start at 
    * @return -1 if not found, index if found
    */
    virtual     long  BinaryFind                (long, long startAt = 0);

    // insert an element into the sorted list preserving the sortedness
    long  BinaryInsert              (long);

    // delete the item at a given poisiton
    void  Delete                    (long, bool = true);

    /**
    * Delete all duplicates in a sorted list
    * Example: SimpleList(1,3,3,5,7).DeleteDuplicates() = [1,3,5,7] 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  DeleteDuplicates          (void);

    /**
    * Delete list of indices in a sorted list
    * Example: SimpleList(1,3,5,7).DeleteList([0,1,2]) = [7] 
    * @param toDelete SimpleList of indices to  
    * @return Nothing. Acts on the List object it was called from. 
    */
    virtual     void  DeleteList                (const _SimpleList&);

    void  Displace                  (long,long,long);
    // shift a range of elements in the array


    /**
    * Add a number to each entry in the array
    * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
    * @param shift Number to add 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Offset                    (long);

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
    bool  NChooseKInit              (_SimpleList&, _SimpleList&, unsigned long, bool = false);

    /**
    * Select the next k-tuple
    * Example: SimpleList(1,3,5,7).DeleteList([0,1,2]) = [7] 
    * @param state the state-storing list previously populated by NChooseKInit
    * @param store the receptacle that will store k-tuples 
    * @return [bool] true is more k-tuples are available; [false] if the last one has just been stored
    */
    bool NChooseK                   (_SimpleList&, _SimpleList&);

    /**
    * TODO:Permute elements in blocks of given size
    * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
    * @param shift Number to add 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Permute                   (long);

    /**
    * TODO:Permute elements in blocks of given size with possible replacement
    * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
    * @param shift Number to add 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  PermuteWithReplacement    (long);

    /**
    * Performs union of two SimpleLists
    * Example: SimpleList(1,3,5,7).Offset(2) = [3, 5, 7, 9]
    * @param shift Number to add 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Union                     (_SimpleList&, _SimpleList&);


    void  Intersect                 (_SimpleList&, _SimpleList&);

    void  XOR                       (_SimpleList&, _SimpleList&);

    void  Subtract                  (_SimpleList&, _SimpleList&);

    long  CountCommonElements       (_SimpleList&, bool = false);


    /**
    * Populate a Simple List with integers incrementally.
    * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
    * @param s The substring to find
    * @param startat The index to start searching from
    * @param increment by Pass true for a case sensitive search 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Merge                     (_SimpleList& l1, _SimpleList& l2, _SimpleList* mergeResults = nil, _SimpleList* mergeResults2 = nil);

    /**
    * Swaps two positions  
    * Example: SimpleList sl[1,3,5].Swap(0, 1) = [3,1]
    * @param i First index to swap 
    * @param j Second index to swap with
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Swap                      (long, long); //swap two elements

    /**
    * Populate a Simple List with integers incrementally.
    * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
    * @param s The substring to find
    * @param startat The index to start searching from
    * @param increment by Pass true for a case sensitive search 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Populate                  (long, long, long); 

    /**
    * Flips the list 
    * Example: SimpleList sl(1,2,3).Flip() = [3,2,1]
    * @return Nothing. Acts on the List object it was called from. 
    */
    void  Flip                      (void); //flip the order of list elements

    virtual     BaseRef toStr                   (void);

    /**
    * TODO
    * Example: SimpleList sl(1,2,3).Flip() = [3,2,1]
    * @return Nothing. Acts on the List object it was called from. 
    */
    void        RecursiveIndexSort              (long from, long to, _SimpleList* index);

    BaseRef     ListToPartitionString           (void);

    virtual     BaseRef makeDynamic             (void);

    virtual     void    Initialize              (bool = true);

    virtual     void    Duplicate               (BaseRef);

    /**
    * Sorts List 
    * Example: SimpleList sl.Sort([5,4,3,2,1]) = [1, 2, 3, 4, 5]
    * @param ascending true if ascending, false for descending sort 
    * @return Nothing. Acts on the List object it was called from. 
    */
    void                Sort                    (bool ascending = true);

    /**
    * SLKP: 20090611    
    * Print the names of variables whose indices are
    * contained in the list
    * @return Nothing. Prints out to screen 
    */
    void                DebugVarList            (void);

    void                ClearFormulasInList     (void);
    /* SLKP: 20110209
        An UGLY hack to automate clearing lists that have pointers to formulas in them */

    /**
    * SLKP: 20090508
    * Return the minimum value in the list 
    * Example: _SimpleList([4, 1, 2]).Min() = 1 
    * @return minimum value in the list 
    */
    long                Min                     (void);

    /**
    * SLKP: 20090508
    * Return the maximum value in the list 
    * Example: _SimpleList([4, 1, 2]).Min() = 1 
    * @return maximum value in the list 
    */
    long                Max                     (void);

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

    _SimpleList*        CountingSort            (long, _SimpleList* = nil);

    void                BubbleSort              (void);
    void                QuickSort               (long, long);

    long*               quickArrayAccess        (void) {
        return (long*)lData;
    }

    /**
    * Compares two parts of the list 
    * Example: SimpleList sl(1,3,5).Compare(0,1) = -1 
    * @param i The index to compare 
    * @param j The second index to compare
    * @return -1 if i<j, 0 if i==j, or 1 if i>j 
    */

    virtual     long    Compare                 (long,long);
    virtual     long    Compare                 (BaseRef,long);

    /**
    * SLKP: 20090316
    * Given a range [from,to] and a given list,
    * make the range conform to the list (e.g. resolve negative to and/or from coordinates)
    * clip the range to fit the list etc
    * Example: SimpleList sl.NormalizeCoordinates(4, 1, 2) = [1, 3, 5, 7]
    * @param from The substring to find
    * @param to The index to start searching from
    * @param refLength The third argument is the length of the list to normalize with respect to.
    * @return Nothing. Acts on the List object it was called from. 
    */
    static      void    NormalizeCoordinates    (long&, long&, const unsigned long);

    friend class _AVLList;

protected:

    // data fields
    unsigned long laLength;
    //memory allocated enough for this many slots

public:

    long*         lData;
    unsigned long lLength;//actual length

};

//_____________________________________________________________________________

class _List:public _SimpleList
{

    // contructor/destructor methods
public:

    /**
    * A constructor.
    * A simple constructor that does nothing
    */
    _List ();


    /**
    * Length constructor.
    * @param sL Length of the string
    */
    _List (unsigned long);


    /**
    * Stack copy contructor
    * @param l List to be copied
    * @param from Beginning index to copy from 
    * @param to Last index to copy to  
    */
    _List (const _List&, long = 0, long = -1);

    /**
    * Construct a list of substrings from the original string separated by char
    * \n\n \b Example: \code _List list = _List((BaseRef)new _String("one,two,three"), ','); \endcode
    * @param ss The substring to be parsed, remember to cast it as a BaseRef
    * @param sep The separator for the string  
    */
    _List (BaseRef, char);

    /**
    * Data constructor (1 member list)
    * \n\n \b Example: \code _List list = _List((BaseRef)new _String("one")); \endcode
    * @param br The object to be changed
    */
    _List (BaseRef);

    /**
    * The deconstructor
    */
    virtual     ~_List(void);

    /**
    * Element location function - read/write
    */
    BaseRef& operator [] (long);

    /**
    * Element location functions - read only
    */
    BaseRef operator () (unsigned long);

    /**
    * Element location functions - read only
    */
    virtual     _List operator = (_List&);

    /**
    * Append operator
    * \n\n \b Example: \code _List result_list = list & append_list; \endcode 
    * @return New concatenated list
    */
    _List operator & (_List&);

    /**
    * Append operator
    * \n\n \b Example: \code _List result_list = list && new _String("one"); \endcode 
    * @return Nothing. Acts on list that is being operated on
    */
    void operator && (BaseRef);

    /**
    * Append operator
    * TODO: need to check this works
    * \n\n \b Example: \code _List result_list = list && "one"; \endcode 
    */
    void operator && (const char*);

    /**
    * Append reference to *this
    * \n\n \b Example: \code _List result_list << new _String("one"); \endcode 
    * @return Nothing. Operates on the _List.
    */
    void operator << (BaseRef);

    /**
    * Appends existing list to *this
    * \n\n \b Example: \code _List result_list << existing_list \endcode 
    * @param l2 The list to be appended
    * @return Nothing. Operates on the _List.
    * @sa AppendNewInstance()
    */
    void operator << (_List&);

    /**
    * Append reference to *this
    * \n\n \b Example: \code _List result_list << existing_list \endcode 
    * @param br The object to be appended
    * @return Nothing. Operates on the _List.
    * @sa AppendNewInstance()
    */
    void AppendNewInstance (BaseRef);

    /**
    * Checks if Lists are identical to each other. Must be _String castable 
    * \n\n \b Example: \code list1.Equal(list2) \endcode 
    * @return bool, true if identical.
    * @sa AppendNewInstance()
    */
    bool Equal       (_List&);

    /**
    * Identical to << operator. Places new value at the end of the list.
    * \n\n \b Example: \code list1.Place(new _String("one")) \endcode 
    * @return Nothing, manipulates *this.
    * @sa InsertElement()
    */
    void Place (BaseRef);


    /**
    * Replace an item
    * \n\n \b Example: \code list.Replace(1, new _String("one"), false); \endcode 
    * @param index The location in the list to be replaced
    * @param newObj The object to be inserted
    * @param dup Allows a duplication
    * @return Nothing, manipulates *this.
    */
    void  Replace (long, BaseRef,bool dup = true);


    virtual     long     FreeUpMemory (long);

    /**
    * Find the position of a search string in the list of strings (ONLY)
    * SLKP: 20100811
    * \n Equivalent to Python's join using the argument as the spacer
    * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
    * @param spacer What you want to be the spacer 
    * @return A pointer to the new list 
    * @sa Find()
    */
    BaseRef Join       (BaseRef);


    /**
    */
    void bumpNInst (void);

    virtual     void Clear     (bool = true);

    // append operator
    _List operator & (BaseRef);

    /**
    * Find the position of a search string in the list of strings (ONLY)
    * \n\n \b Example: \code 
    * _String* needle = new _String("two") 
    * _List("zero","one,"two").Find((BaseRef)needle)
    * \endcode
    * @param s The integer to find
    * @return -1 if not found, index if found
    */
    virtual     long  Find (BaseRef, long startat = 0);

    virtual     long  FindPointer (BaseRef b, long startat = 0) {
        return _SimpleList::Find ((long)b, startat);
    }

    /**
    * Find the position of a search string in the list of strings (ONLY)
    * \n Faster than the Find(), since it assumes string entries
    * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
    * @param s The substring to find
    * @param startat The index to start searching from
    * @param caseSensitive Pass true for a case sensitive search 
    * @param upTo Upper limit for search index. 
    * @return -1 if not found, the index if it is found.
    * @sa Find()
    * @sa BinaryFind()
    */
    virtual     long  FindString          (BaseRef, long startat = 0, bool caseSensitive = true, long upTo = -1);

    /**
    * Find the position of a search string in the list of strings (ONLY)
    * \n Faster than the Find(), since it assumes string entries
    * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
    * @param s The substring to find
    * @param startat The index to start searching from
    * @param caseSensitive Pass true for a case sensitive search 
    * @param upTo Upper limit for search index. 
    * @return -1 if not found, the index if it is found.
    * @sa Find()
    * @sa BinaryFind()
    */
    virtual     long  BinaryFind          (BaseRef );

    // insert an element into the sorted list preserving the sortedness
    long  BinaryInsert        (BaseRef);

    // delete the item at a given poisiton
    void  Delete (long );

    // delete the item at a given poisiton
    virtual     void  DeleteList (const _SimpleList&);

    /**
    * Populate a Simple List with integers incrementally.
    * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
    * @param s The substring to find
    * @param startat The index to start searching from
    * @param increment by Pass true for a case sensitive search 
    * @return Nothing. Acts on the List object it was called from. 
    */
    virtual     void InsertElement (BaseRef br, long insertAt = -1, bool store = true);


    virtual     BaseRef toStr           (void);

    virtual     void    toFileStr       (FILE*);

    virtual     BaseRef makeDynamic     (void);

    virtual     void    Duplicate       (const BaseRef);

    virtual     long    Compare         (long,long);
    virtual     long    Compare         (BaseRef,long);

    void    Intersect       (_List&, _List&, _SimpleList* = nil, _SimpleList* = nil);


};

//_____________________________________________________________________________

class _AVLList: public BaseObj
{

public:

    _AVLList    (_SimpleList*);
    virtual     ~_AVLList   (void) {}


    long                Find                (BaseRef);
    long                FindLong            (long);
    // 20100623: a shortcut function to look for integers only
    // avoids calling ::Compare
    char                FindBest            (BaseRef, long&);
    long                Find                (BaseRef, _SimpleList&);
    long                Next                (long,    _SimpleList&);
    long                Prev                (long,    _SimpleList&);
    long                First               (void);
    long                Last                (void);
    long                GetByIndex          (const long);
    long                Insert              (BaseRef, long = 0, bool = true, bool = false);
    // the 1st bool flag is to say whether to dup the object being inserted
    // the 2nd bool flag (if the first flag is false) if set to true,
    // will cause failed inserts (key already exists) to delete the key
    void                Delete              (BaseRef, bool = false);
    virtual  void       ReorderList         (_SimpleList* = nil);
    BaseRef             Retrieve            (long);
    virtual  long       InsertData          (BaseRef, long, bool);
    unsigned long       countitems          (void);
    virtual  void       Clear               (bool = false);
    virtual  BaseRef    toStr               (void);
    virtual  bool       HasData             (long);
    virtual  long       Traverser           (_SimpleList&, long &, long = -1);
    virtual  long       GetRoot             (void) {
        return root;
    }
    virtual  void       DeleteXtra          (long) {};
    virtual  void       DeleteAll           (bool cL) {
        Clear(cL);
        DeleteObject (dataList);
    }

    // data members

    void                ConsistencyCheck    (void);

    _SimpleList *dataList,
                leftChild,
                rightChild,
                balanceFactor,
                emptySlots;

    long        root;

};

//_____________________________________________________________________________

class _AVLListX: public _AVLList
{

public:

    _AVLListX   (_SimpleList*);
    virtual     ~_AVLListX  (void) {}



    virtual  BaseRef    toStr               (void);
    virtual  long       InsertData          (BaseRef, long, bool);
    long        GetXtra             (long);
    void        SetXtra             (long,long);
    virtual  void       Clear               (bool = false);
    virtual  void       DeleteXtra          (long);
    virtual  void       PopulateFromList    (_List&);
    /* SLKP: 20090817
        add key: index values from the list of strings
     */

    _SimpleList         xtraD;

};

//_____________________________________________________________________________

class _AVLListXL: public _AVLList
{

public:

    _AVLListXL  (_SimpleList*);
    virtual     ~_AVLListXL (void) {}


    virtual  BaseRef    toStr               (void);
    virtual  long       InsertData          (BaseRef, long,bool);
    BaseRef GetXtra             (long);
    void        SetXtra             (long,BaseRef,bool);
    virtual  void       Clear               (bool = false);
    virtual  void       DeleteXtra          (long);

    _List               xtraD;

};

//_____________________________________________________________________________

void        SortLists (_SimpleList*, _SimpleList*);

extern      _List     pathNames;

#endif
