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

#ifndef _HLIST_
#define _HLIST_
//#pragma once
#include "simplelist.h"

#define  MEMORYSTEP 8

//_____________________________________________________________________________

class _String;
// forward declaration

//_____________________________________________________________________________

class _List:public _SimpleList {

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
          variadic initializer
         
         */
        template <typename... T>
        _List (const BaseRef& first, const T&... args) {
          AppendObjectsToList(first);
          AppendObjectsToList (args...);
        };


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
        BaseRef operator () (const unsigned long);

  
        /**
         * Element location functions - read only
         * used to avoid (*list)(3) which are hard to read
         * checks for valid index range and returns NULL if outside the range
         */
        virtual BaseRef GetItemRangeCheck     (const unsigned long) const;

        /**
        * Element location functions - read only
        */
        const _List& operator = (const _List&);

        template <typename MAPPER> void ForEach (MAPPER mapper, long startAt = 0) const {
          for (unsigned long i = startAt; i<lLength; i++) {
            mapper ( ((BaseRef*)(list_data))[i], i );
          }
        }
    
        template <typename MAPPER> const _List Map (MAPPER mapper, long startAt = 0, long endAt = -1L) const {
            _List result;
            if (endAt < 0) {
                endAt = lLength;
            } else if (endAt > lLength) {
                endAt = lLength;
            }
            for (unsigned long i = startAt; i<endAt; i++) {
                result < mapper ( ((BaseRef*)(list_data))[i], i );
            }
            return result;
        }

        /**
        * Append operator
        * \n\n \b Example: \code _List result_list = list & append_list; \endcode 
        * @return New concatenated list
        */
        const _List operator & (_List const&) const;

        /**
        * Append operator
        * \n\n \b Example: \code _List result_list = list && new _String("one"); \endcode 
        * @return Nothing. Acts on list that is being operated on
        */
        _List& operator && (BaseRef);

        /**
        * Append operator
        * \n\n \b Example: \code _List result_list = list && "one"; \endcode 
        */
        _List& operator && (const char*);

        /**
        * Append reference to *this (<< also increments the reference counter)
        * \n\n \b Example: \code _List result_list << new _String("one"); \endcode 
        * @return Nothing. Operates on the _List.
        */
        _List& operator << (BaseRef);
        _List& operator <  (BaseRef);

        /**
        * Appends existing list to *this (<< also increments the reference counters)
        * \n\n \b Example: \code _List result_list << existing_list \endcode 
        * @param l2 The list to be appended
        * @return Nothing. Operates on the _List.
        * @sa AppendNewInstance()
        */
        _List& operator << (_List const&);
        _List& operator < (_List const&);
        _List& operator < (const char *);

        /**
        * Append operator
        */
        const _List operator & (BaseRef) const;

        /**
        * @sa Equal()
        */
        bool operator == (_List const&) const;
  
  
        template <typename... Args>
        void AppendObjectsToList (BaseRef const& first, const Args&... args) {
          AppendNewInstance (first);
          AppendObjectsToList (args...);
        }
  
        void AppendObjectsToList (BaseRef const& first) {
          AppendNewInstance ((BaseRef)first);
        }


        /**
        * Append reference to *this
        * \n\n \b Example: \code _List result_list << existing_list \endcode 
        * @param br The object to be appended
        * @return The length of the list following element addition
        * @sa AppendNewInstance()
        */
        long AppendNewInstance(BaseRef);
  
        /**
         * Append a variable number of arguments without increasing ref counts
         
         * @sa _List constructors()
         */
        void AppendNewInstance (BaseObj* ref, const unsigned long number, ...);

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
        virtual long BinaryFindObject (BaseObj const *, long startAt = 0) const;

        template <typename FILTER> long FindOnCondition (FILTER condition, long startAt = 0) const {
            for (unsigned long i = startAt; i<lLength; i++) {
                if ( condition (((BaseRefConst*)(list_data))[i], i) ) {
                    return i;
                }
            }
            return kNotFound;
        }

        /**
        * Insert an element into the sorted list preserving the sortedness
        */
        long BinaryInsert(BaseRef);

        /**
        */
        void bumpNInst(void);

        /**
        */
        virtual hyComparisonType Compare(long,long) const;

        /**
        */
        virtual hyComparisonType Compare(BaseObj const *,long) const;

        /**
        * Return number of elements 
        */
        unsigned long Count();

        /**
        */
        virtual void Clear(bool=true);

        /**
        * Delete the item at a given poisiton
        * if the bool flag is false, then only remove the reference to the object 
        */
        void Delete(long, bool = true);

        /**
        * Delete all items from the starting position to the end of the list
        * if the bool flag is false, then only remove the reference to the object
        */
        void DeleteTail (long, bool = true);

        /**
        */
        virtual void Duplicate(BaseRefConst);

        /**
        * Delete the item at a given poisiton
        */
        virtual void DeleteList(const _SimpleList&);

        /**
        * Checks if Lists are identical to each other. Must be _String castable 
        * \n\n \b Example: \code list1.Equal(list2) \endcode 
        * @return bool, true if identical.
        * @sa AppendNewInstance()
        */
        bool Equal(_List const&) const;

        /**
        * Find the position of a search string in the list of strings (ONLY)
        * \n\n \b Example: \code 
        * _String* needle = new _String("two") 
        * _List("zero","one,"two").Find((BaseRef)needle)
        * \endcode
        * @param s The integer to find
        * @return -1 if not found, index if found
        */
        virtual long FindObject (BaseRefConst, long startat = 0) const;

        /**
        */
        virtual long FindPointer(BaseRef b, long startat = 0) {
            return _SimpleList::Find((long)b, startat);
       }

        /**
         * Element location functions - read only
         * used to avoid (*list)(3) which are hard to read
         */
        template <typename INDEX> BaseRef GetItem     (INDEX index) const {
          return ((BaseRef*)list_data)[index];
        }

        template<typename INDEX, typename ...INDICES> BaseRef GetItem (INDEX i0, INDICES... rest) const {
           return ((_List*)GetItem (i0))->GetItem (rest...);
         }
 
 

        /**
        * Populate a Simple List with integers incrementally.
        * Example: SimpleList sl.Populate(4, 1, 2) = [1, 3, 5, 7]
        * @param s The substring to find
        * @param startat The index to start searching from
        * @param increment by Pass true for a case sensitive search 
        * @return Nothing. Acts on the List object it was called from. 
        */
        virtual void InsertElement(BaseRef br,long insertAt=-1, bool store=true, bool pointer=true);

        /**
        */
        void Intersect(_List&, _List&, _SimpleList* = nil, _SimpleList* = nil);

        /**
        * Find the position of a search string in the list of strings (ONLY)
        * SLKP: 20100811
        * \n Equivalent to Python's join using the argument as the spacer
        * \n\n \b Example: \code _String ("AABBCC").Find("B")\endcode
        * @param spacer What you want to be the spacer 
        * @param startAt start at this list element 
        * @param endAt end at this list element 
        * @return A pointer to the new string 
        * @sa Find()
        */
        _String* Join(_String const& spacer, long startAt = 0, long endAt = -1) const;

        /**
        * Identical to << operator. Places new value at the end of the list.
        * \n\n \b Example: \code list1.Place(new _String("one")) \endcode 
        * @return Nothing, manipulates *this.
        * @sa InsertElement()
        */
        void Place(BaseRef);
        
        
        /**
        * Map the values of the first list to the corresponding indices in the second list (treated as STRINGS).
        * \n\n \b Example: \code ("a","b","d").Map (["b","c","a"], mapping); mapping = (2,0,-1)\endcode 
        * @param target The target list of the mapping
        * @param mapping The list that will store the mapping
        * @return Nothing, manipulates mapping.
        */

        void Map (_List& target, _SimpleList& mapping);

        /**
        */
        virtual BaseRef makeDynamic(void) const;

        /**
        * Replace an item
        * \n\n \b Example: \code list.Replace(1, new _String("one"), false); \endcode 
        * @param index The location in the list to be replaced
        * @param newObj The object to be inserted
        * @param dup Allows a duplication
        * @return Nothing, manipulates *this.
        */
        void Replace(long,BaseRef,bool dup=true);

        /**
        */
        virtual BaseRef toStr(unsigned long = 0UL);

        /**
        */
        virtual void toFileStr(FILE*, unsigned long = 0UL);
  
        /**
         Generate a string that is not present in the list.
         The string will look like 'base_[autogenerated number]'
         
         @param base the prefix to use for the name
         @sorted is the list sorted?
         
         @return a unique string
         
         */
        const _String GenerateUniqueNameForList (_String const& base, bool sorted) const;
    
        virtual bool is_numeric_list (void) const { return false;}


};

//TODO:Can we avoid using this extern?
extern _List pathNames;

#endif
