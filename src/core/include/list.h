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
#include "simplelist.h"

#define  MEMORYSTEP 8

//_____________________________________________________________________________

//_____________________________________________________________________________

class _List:public _SimpleList
{

    public:

        /**
        * A constructor.
        * A simple constructor that does nothing
        */
        _List();


        /**
        * Length constructor.
        * @param sL Length of the string
        */
        _List(unsigned long);


        /**
        * Stack copy contructor
        * @param l List to be copied
        * @param from Beginning index to copy from 
        * @param to Last index to copy to  
        */
        _List(const _List&,long=0,long=-1);

        /**
        * Construct a list of substrings from the original string separated by char
        * \n\n \b Example: \code _List list = _List((BaseRef)new _String("one,two,three"), ','); \endcode
        * @param ss The substring to be parsed, remember to cast it as a BaseRef
        * @param sep The separator for the string  
        */
        _List(BaseRef,char);

        /**
        * Data constructor (1 member list)
        * \n\n \b Example: \code _List list = _List((BaseRef)new _String("one")); \endcode
        * @param br The object to be changed
        */
        _List(BaseRef);

        /**
        * The deconstructor
        */
        virtual ~_List(void);

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
        virtual _List operator = (_List&);

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

        // append operator
        _List operator & (BaseRef);


        //Return number of elements 
        unsigned long Count();


        /**
        * Append reference to *this
        * \n\n \b Example: \code _List result_list << existing_list \endcode 
        * @param br The object to be appended
        * @return Nothing. Operates on the _List.
        * @sa AppendNewInstance()
        */
        void AppendNewInstance(BaseRef);

        /**
        * Checks if Lists are identical to each other. Must be _String castable 
        * \n\n \b Example: \code list1.Equal(list2) \endcode 
        * @return bool, true if identical.
        * @sa AppendNewInstance()
        */
        bool Equal(_List&);

        /**
        * Identical to << operator. Places new value at the end of the list.
        * \n\n \b Example: \code list1.Place(new _String("one")) \endcode 
        * @return Nothing, manipulates *this.
        * @sa InsertElement()
        */
        void Place(BaseRef);


        /**
        * Replace an item
        * \n\n \b Example: \code list.Replace(1, new _String("one"), false); \endcode 
        * @param index The location in the list to be replaced
        * @param newObj The object to be inserted
        * @param dup Allows a duplication
        * @return Nothing, manipulates *this.
        */
        void Replace(long,BaseRef,bool dup=true);


        virtual long FreeUpMemory(long);

        /**
        * Find the position of a search string in the list of strings (ONLY)
        * SLKP: 20100811
        * \n Equivalent to Python's join using the argument as the spacer
        * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
        * @param spacer What you want to be the spacer 
        * @return A pointer to the new list 
        * @sa Find()
        */
        BaseRef Join(BaseRef);

        /**
        */
        void bumpNInst(void);

        virtual void Clear(bool=true);


        /**
        * Find the position of a search string in the list of strings (ONLY)
        * \n\n \b Example: \code 
        * _String* needle = new _String("two") 
        * _List("zero","one,"two").Find((BaseRef)needle)
        * \endcode
        * @param s The integer to find
        * @return -1 if not found, index if found
        */
        virtual long Find(BaseRef, long startat = 0);

        virtual long FindPointer(BaseRef b, long startat = 0) {

            return _SimpleList::Find((long)b, startat);

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
        virtual long FindString(BaseRef,long startat=0,bool caseSensitive=true,long upTo=-1);

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
        virtual long BinaryFind(BaseRef);

        // insert an element into the sorted list preserving the sortedness
        long BinaryInsert(BaseRef);

        // delete the item at a given poisiton
        void Delete(long);

        // delete the item at a given poisiton
        virtual void DeleteList(const _SimpleList&);

        /**
        * Populate a Simple List with integers incrementally.
        * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
        * @param s The substring to find
        * @param startat The index to start searching from
        * @param increment by Pass true for a case sensitive search 
        * @return Nothing. Acts on the List object it was called from. 
        */
        virtual void InsertElement(BaseRef br,long insertAt=-1, bool store=true);

        virtual BaseRef toStr(void);
        virtual void toFileStr(FILE*);
        virtual BaseRef makeDynamic(void);
        virtual void Duplicate(const BaseRef);
        virtual long Compare(long,long);
        virtual long Compare(BaseRef,long);
        void Intersect(_List&, _List&, _SimpleList* = nil, _SimpleList* = nil);
};

//TODO:Can we avoid using this extern?
extern _List pathNames;

#endif
