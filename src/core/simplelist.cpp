/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)
 
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

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

#include "global_things.h"
#include "hy_strings.h"
#include "hy_string_buffer.h"
#include "list.h"
#include "simplelist.h"
#include "parser.h"
#include "mersenne_twister.h"
#include "function_templates.h"

using namespace hy_global;


/*
==============================================================
Constructors
==============================================================
*/

// Does nothing
_SimpleList::_SimpleList () : static_data {0L} {
    Initialize(false);
}

//Data constructor (1 member list)
_SimpleList::_SimpleList (long br) : static_data {0L} {
    lLength     = 1;
    laLength    = MEMORYSTEP;
    list_data   = (long*)&static_data;
    static_data[0] = br;
}

//Length constructor and populator
_SimpleList::_SimpleList (long l, long start, long step) : static_data {0L} {
    Initialize (false);
    Populate   (l,start,step);
}

//Length constructor
_SimpleList::_SimpleList (unsigned long l) : static_data {0L} {
    lLength = 0UL;
    if (l <= MEMORYSTEP) {
        laLength = MEMORYSTEP;
        list_data   = (long*)&static_data;
    } else {
        laLength = (l/MEMORYSTEP + 1)*MEMORYSTEP;
        list_data = (long*)MemAllocate (laLength * sizeof(hyPointer));
        memset (list_data,0,laLength * sizeof(hyPointer));
    }
}

//Length constructor
_SimpleList::_SimpleList (unsigned long l, long * preallocated_storage) : static_data {0L} {
    lLength = 0UL;
    if (preallocated_storage) {
        laLength = l;
        list_data = preallocated_storage;
        AddAReference();
    } else {
        if (l <= MEMORYSTEP) {
            laLength = MEMORYSTEP;
            list_data   = (long*)&static_data;
        } else {
            laLength = (l/MEMORYSTEP + 1)*MEMORYSTEP;
            list_data = (long*)MemAllocate (laLength * sizeof(hyPointer));
        }
        memset (list_data,0,laLength * sizeof(hyPointer));
    }
}

//Stack copy contructor
_SimpleList::_SimpleList (_SimpleList const & l, long from, long to) : static_data {0L} {
  Initialize           (false);
  if (from == 0L && to == -1L) {
    to = l.lLength;
  } else {
    NormalizeCoordinates (from, to, l.lLength);
  }
  if (to > from) {
      long upto = to-from ;
      RequestSpace(upto);
      for (lLength = 0UL; lLength < upto; lLength++) {
          list_data [lLength] = l.list_data[from+lLength];
      }
  }
}

// Data constructor (variable number of long constants)
_SimpleList::_SimpleList (const long value1, const unsigned long number, ...) : static_data {0L} {
    Initialize (true);
    va_list vl;
    
    (*this) << value1;
    
    va_start(vl,number);
    for (unsigned long arg_id =0;arg_id<number;arg_id++) {
        const long this_arg =va_arg(vl,long);
        (*this) << this_arg;
    }
    va_end(vl);
}

//Destructor
_SimpleList::~_SimpleList(void) {
    if (CanFreeMe()) {
        if (list_data && list_data != (long*)static_data) {
            free (list_data);
        }
    } else {
        RemoveAReference();
    }
}

/*
==============================================================
Operator Overloads
==============================================================
*/

//Element location functions (0,llength - 1)
long& _SimpleList::operator [] (const long i) {
    return list_data[i];
}

//Element location functions (0,llength - 1)
long _SimpleList::operator () (const unsigned long i) const {
    //if (lLength == 0) return 0;
    //Is there a reason why this is commented out?
    //if (i>=lLength) i = lLength-1;
    if(i < lLength) {
        return list_data [i];
    }
    HandleApplicationError ("List index out of range");
    return -1;
}

//Assignment operator
const _SimpleList& _SimpleList::operator = (_SimpleList const &l) {
    if (l.laLength && list_data && laLength >= l.lLength && laLength <= l.laLength && CanFreeMe()) {
        // reuse existing memory allocation
      lLength  = l.lLength;
      if (lLength) {
        memcpy (list_data,l.list_data,lLength*sizeof (hyPointer));
      }
    } else { // allocate new memory
      Clear();
      lLength  = l.lLength;
      laLength = l.laLength;
      if (laLength >= MEMORYSTEP) {
          list_data = (long*)MemAllocate (laLength*sizeof (hyPointer));
      } else {
          list_data = static_data;
      }
      if (lLength) {
        memcpy (list_data,l.list_data,lLength*sizeof (hyPointer));
      }

    }

    return *this;
}

//Append operator
_SimpleList _SimpleList::operator & (_SimpleList const& l) {
    _SimpleList res (l.lLength + lLength);
    
    if (!res.laLength) {
        return res;
    }

    if (list_data && lLength) {
        memcpy(res.list_data,list_data,lLength*sizeof(void*));
    }

    if (l.list_data&&l.lLength) {
        memcpy((char*)res.list_data+lLength*sizeof(void*),l.list_data, l.lLength*sizeof(void*));
    }

    res.lLength = l.lLength + lLength;

    return res;
}

_SimpleList& _SimpleList::operator << (long br) {
  _SimpleList::InsertElement ((BaseRef)br, -1, false, false);
  return *this;
}

bool _SimpleList::operator >> (long br) {
    if (Find(br) == -1) {
        InsertElement ((BaseRef)br, -1, false, false);
        return true;
    }
    return false;
}

void _SimpleList::operator << (_SimpleList const & source) {
    for (unsigned long k=0UL; k<source.lLength; k++) {
        (*this) << source.list_data [k];
    }
}

/*
==============================================================
Methods
==============================================================
*/


//Element location functions (0,llength - 1), negative values return 
// elements from the end of the list

long _SimpleList::GetElement (const long index) const {
    if (index >= 0L) {
        if ((const unsigned long) index < lLength) {
            return list_data [index];
        }
    } 
    if ((const unsigned long) (-index) <= lLength) {
        return list_data [lLength + index];
    }
    HandleApplicationError (_String("List index '") & (long)((const unsigned long) (-index)) & "' out of range in _SimpleList::GetElement on list of length " & long (lLength));
    return 0;
}


long  _SimpleList::BinaryFind (long s, long startAt) const {
    
    if (lLength == 0L) {
        return -2;
    }
    
    long top    =   lLength-1,
         bottom  =   startAt,
         middle;


    while (top>bottom) {
        middle = (top+bottom) >> 1;
        if ( s<list_data [middle]) {
            top = middle==top? top-1 : middle;
        } else if (s > list_data [middle]) {
            bottom = middle == bottom ? bottom+1 : middle;
        } else {
            return middle;
        }
    }

    middle     = top;
    long comp  = list_data[middle]-s;
    if (!comp) {
        return middle;
    }

    return comp<0?-middle-3:-middle-2;
}

long  _SimpleList::BinaryInsert (long n) {
    if (lLength == 0L) {
        (*this) << n;
        return 0;
    }

    long pos = -BinaryFind (n)-2;

    if (pos<0L) {
        return -pos+2;
    }

    if (pos < lLength && list_data [pos]<n) {
        pos++;
    }

    InsertElement ((BaseRef)n,pos,false,false);

    return pos>=lLength?lLength-1:pos;
}

void _SimpleList::ClearFormulasInList(void) {
    for (unsigned long k = 0UL; k < lLength; k++)
        if (list_data[k]) {
            delete (_Formula*) list_data [k];
        }
}

long _SimpleList::Sum (void) const{
    long sum = 0L;
    for (unsigned long k = 0UL; k < lLength; k++) {
       sum += list_data [k];
    }
    return sum;
}

hyComparisonType  _SimpleList::Compare (long i, long j) const {
    long    v1 = list_data[i],
            v2 = list_data[j];

    if (v1<v2) {
        return kCompareLess;
    } else if (v1==v2) {
        return kCompareEqual;
    } else {
        return kCompareGreater;
    }


    //return ((long*)list_data)[i]-((long*)list_data)[j];
}

hyComparisonType  _SimpleList::Compare (BaseObj const *i, long j) const {
    long    v1 = (long)i,
            v2 = list_data[j];


    if (v1<v2) {
        return kCompareLess;
    } else if (v1==v2) {
        return kCompareEqual;
    } else {
        return kCompareGreater;
    }

    //return (long)i-((long*)list_data)[j];
}

// Compute the number of shared of two sorted lists
long    _SimpleList::CountCommonElements (_SimpleList const& l1, bool yesNo) const {
    long  c1    = 0,
          c2    = 0,
          res   = 0;


    while (c1<l1.lLength && c2<lLength) {
        while (l1.list_data[c1]<list_data[c2]) {
            c1++;
            if (c1==l1.lLength) {
                break;
            }
        }
        if (c1==l1.lLength) {
            break;
        }

        while (l1.list_data[c1]==list_data[c2]) {
            c2++;
            if (yesNo) {
                return 1;
            } else {
                res++;
            }
            if (c1==l1.lLength || c2==lLength) {
                break;
            }
        }
        if (c1==l1.lLength || c2==lLength) {
            break;
        }
        while (list_data[c2]<l1.list_data[c1]) {
            c2++;
            if (c2==lLength) {
                break;
            }
        }
    }

    return res;
}



_SimpleList*  _SimpleList::CountingSort (long upperBound, _SimpleList* ordering, bool wantResult) {
    static const unsigned long kAllocaLimit = 4096UL;
    
    if (ordering) {
        ordering->Clear();
    }

    if (lLength) {
        if (upperBound <= 0) {
            upperBound = Max()+1;
        }

        long * storage = upperBound < kAllocaLimit ? (long*)alloca (sizeof (long) * upperBound) : nil;
        
        if (storage) {
            memset (storage, 0, sizeof (long) * upperBound);
        }
        
        _SimpleList buffer     (upperBound, storage);
        _SimpleList *result    =  (wantResult || !ordering) ? new _SimpleList (lLength) : nil;
        
        for (long pass1 = 0; pass1 < lLength; pass1 ++) {
            buffer.list_data[list_data[pass1]] ++;
        }
        for (long pass2 = 1; pass2 < upperBound; pass2 ++) {
            buffer.list_data[pass2] += buffer.list_data[pass2-1];
        }
        
        if (!wantResult && ordering) {
            ordering->Populate (lLength, 0, 0);
            for (long pass3 = lLength-1; pass3 >=0; pass3--) {
                ordering->list_data[--buffer.list_data[list_data[pass3]]] = pass3;
            }
        } else {
            if (ordering) {
                ordering->Populate (lLength, 0, 0);
                for (long pass3 = lLength-1; pass3 >=0; pass3--) {
                    result->list_data[--buffer.list_data[list_data[pass3]]] = list_data[pass3];
                    ordering->list_data[buffer.list_data[list_data[pass3]]] = pass3;
                }
            } else
                for (long pass3 = lLength-1; pass3 >=0; pass3--) {
                    result->list_data[--buffer.list_data[list_data[pass3]]] = list_data[pass3];
                }
            
            result->lLength = lLength;
        }
        
        

        return result;
    }
    return wantResult ? new _SimpleList : nil;
}

void  _SimpleList::Clear (bool completeClear) {
    if (CanFreeMe()) {
        lLength = 0UL;
        if (completeClear) {
            laLength = MEMORYSTEP;
            if (is_dynamic()) {
                free (list_data);
                list_data = static_data;
            }
        }
    } else {
        RemoveAReference();
    }
}

void _SimpleList::DebugVarList(void) {
    printf ("\nVariable list dump:\n");
    for  (unsigned long e = 0UL; e < lLength; e++) {
        if (list_data[e] >= 0L) {
            _Variable * theV = LocateVar (list_data[e]);
            if (theV) {
                printf ("[%s]\n", theV->GetName()->get_str());
                continue;
            }
        }
        printf ("[Empty]\n");
    }
}

void  _SimpleList::_EnsureCorrectStorageType (void) {
    if (laLength >= MEMORYSTEP) {
        if (list_data == static_data) {
            list_data = (long*)MemAllocate(laLength*sizeof(hyPointer));
            _CopyStatic ();

        } else {
            list_data = (long*)MemReallocate ((long*)list_data, laLength*sizeof(hyPointer));
        }
    } else {
        if (list_data != static_data) {
            for (long k = 0L; k < lLength; k++) {
                static_data [k] = list_data[k];
            }
            free (list_data);
            list_data = static_data;
        }
    }
}

void  _SimpleList::_UpdateStorageType (void) {
    if (laLength-lLength>MEMORYSTEP) {
        laLength -= ((laLength-lLength)/MEMORYSTEP)*MEMORYSTEP;
        _EnsureCorrectStorageType ();
    }
}

void  _SimpleList::_CopyStatic(void) {
    for (long k = 0L; k < MEMORYSTEP; k++) {
        list_data[k] = static_data[k];
    }
}

//Delete item at index (>=0)
void  _SimpleList::Delete (long index, bool compact) {
    if (index>=0 && index<lLength) {
        lLength--;
        if (lLength > index) {
            for (long k = index; k < lLength; k++) {
                list_data[k] = list_data[k+1];
            }
            //memmove ((hyPointer)list_data+sizeof(BaseRef)*(index),(hyPointer)list_data+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
        }
    }
    if (compact) {
        _UpdateStorageType ();
    }

}

//Delete duplicates from a sorted list
void  _SimpleList::DeleteDuplicates (void) {
    if (lLength>1L) {
        /*_SimpleList noDups;

        long    lastValue = list_data [0]+1;
        for (unsigned long k=0; k<lLength; k++) {
            long thisValue = list_data[k];
            if (thisValue!=lastValue) {
                noDups << thisValue;
                lastValue = thisValue;
            }
        }

        if (noDups.lLength != lLength) {
            Duplicate (&noDups);
        }*/
        long last_retained   = 0L,
             current_store   = 1L;
        
        for (long k = 1L; k < lLength; k++) {
            if (list_data[k] != list_data[last_retained]) {
                last_retained = k;
                list_data[current_store++] = list_data[k];
            }
        }
        
        if (lLength != current_store) {
            lLength = current_store;
            TrimMemory();
        }
    }

}



//Delete items from a sorted list
void  _SimpleList::DeleteList (const _SimpleList& toDelete) {
    if (toDelete.lLength) {
        unsigned long k = 0;
        for (unsigned long i = 0UL; i<lLength; i++) {
            if (k<toDelete.lLength && i==toDelete.list_data[k])
                //if (k<toDelete.lLength)
            {
                k++;
            } else {
                list_data[i-k] = list_data[i];
            }
        }
        lLength -= toDelete.lLength;
    }
    _UpdateStorageType ();
}


long _SimpleList::SkipCorrect (long index) const {
  for (unsigned long k=0UL; k < lLength; k++)
    if (index >= list_data[k]) {
      index++;
    }
  return index;
}

long _SimpleList::CorrectForExclusions(long index, long excluded_value) const {
  long correction = 0L;
  for (unsigned long k=0UL; k < lLength && index >= list_data[k]; k++) {
    if (index == list_data[k]) {
      return excluded_value;
    }
    correction++;
  }
  return index - correction;
}

long _SimpleList::CorrectForExclusions(long * index,  long length) const {
  long exclusion_index = 0L,
       mapped = 0UL;
  
  for (long k=0UL; k < length; k++) {
    if (exclusion_index < lLength && index[k] >= list_data[exclusion_index]) {
      if (index[k] > list_data[exclusion_index]) {
        k--;
      }
      exclusion_index++;
      continue;
    }
    index[mapped++] = index[k] - exclusion_index;
  }
  return mapped;
}


//Shift the range from start to end
void  _SimpleList::Displace (long start, long end, long delta) {
    if (start<0) {
        start = 0;
    } else if (start>=lLength) {
        start = lLength-1;
    }

    if (end<0) {
        end = lLength-1;
    } else if (end>=lLength) {
        end = lLength-1;
    }

    if ((end-start>=0) && delta && (end-start<lLength-1))
    {
        if (delta>0 && lLength-end <= delta) { // shift up
            delta = lLength-end-1;
        } else if (start-delta<0){
            delta = start;
        }

        if (delta) {
            long i,j,delta2 = end-start+1;
            _SimpleList swapList ((unsigned long)(end-start+1));

            for (i=start; i<=end; i++) {
                swapList << list_data[i];
            }

            if (delta>0) {

                for (i=end+1; i<=end+delta; i++) {
                    list_data[i-delta2] = list_data[i];
                }

            } else {
                for (i=start-1; i>=start+delta; i--) {
                    list_data[i+delta2] = list_data[i];
                }
            }

            for (i=start+delta, j=0; i<=end+delta; i++,j++) {
                list_data[i] = swapList.list_data[j];
            }
        }

    }
}

void _SimpleList::Duplicate(BaseRefConst theRef) {
    _SimpleList const* l  = (_SimpleList const*)theRef;
    lLength         = l->lLength;
    laLength        = l->laLength;
    //list_data       = l->list_data;
    if (l->is_dynamic()) {
        list_data = (long*)MemAllocate (laLength*sizeof (hyPointer));
        memcpy (list_data, l->list_data, lLength*sizeof (hyPointer));
    } else {
        list_data = static_data;
        for (long k = 0L; k < lLength; k++) {
            list_data [k] = l->list_data[k];
        }
    }
}

//Element location functions (0,llength - 1)
//Negative indices return offsets from the end of the list
long _SimpleList::Element(long index) const {
    if (index >= 0 && index < lLength) {
        return list_data[index];
    }

    else if(index < 0L && -index <= lLength) {
        return list_data[(long)lLength+index];
    }

    return 0L;
}

bool _SimpleList::Equal(_SimpleList const& l2) const
{
    if (lLength == l2.lLength) {      
      for (unsigned long i=0UL; i<lLength; i++)
          if (list_data[i] != l2.list_data[i]) {
              return false;
          }
      
      return true;
    }
    return false;
}

long  _SimpleList::Find (long s, long startAt) const {
    for (unsigned long i = startAt; i<lLength; i++) {
        if ( list_data[i] == s ) {
            return i;
        }
    }
    return kNotFound;
}

long  _SimpleList::FindStepping (long s, long step, long startAt) {
    for (unsigned long i = startAt; i<lLength; i+=step) {
        if (list_data [i] == s) {
            return i;
        }
    }

    return kNotFound;
}

void  _SimpleList::FilterRange (long lb, long ub) {
    if (ub <= lb) {
        Clear();
    } else {
        /*_SimpleList toDelete;
        for (long k = 0; k < lLength; k++)
            if (list_data[k] <= lb || list_data[k] >= ub) {
                toDelete << k;
            }
        DeleteList (toDelete);*/
        unsigned long current_store = 0UL;
        for (long k = 0; k < lLength; k++) {
            if (list_data [k] >= lb || list_data[k] <= ub) {
                list_data[current_store++] = list_data[k];
            }
        }
        
        if (lLength != current_store) {
            lLength = current_store;
            TrimMemory();
        }
        
        
    }
}

void _SimpleList::Flip () {
    for (long k=0L, l=lLength-1; k<l; k++,l--) {
        Exchange(list_data[k], list_data[l]);
    }
}

void _SimpleList::Initialize(bool doMemAlloc) {
    BaseObj::Initialize();
    lLength = 0UL;
    /*if (doMemAlloc) {
        laLength = MEMORYSTEP;
        list_data = (long*)MemAllocate (laLength * sizeof(hyPointer));
    } else {
        laLength = 0;
        list_data    = nil;
    }*/
    laLength = MEMORYSTEP;
    list_data = static_data;
}

//Append & store operator
void _SimpleList::InsertElement (BaseRef br, long insertAt, bool store, bool pointer) {
    lLength++;
    if (lLength>laLength) {
        unsigned long incBy = ((MEMORYSTEP << 2) > lLength)? MEMORYSTEP: (lLength >> 2);

        laLength+=incBy;

        if (is_dynamic()) {
            list_data = (long*)MemReallocate((char*)list_data, laLength*sizeof(long));
        } else {
            list_data = (long*) MemAllocate(laLength*sizeof(long));
            _CopyStatic();
        }

    }
  
    if (insertAt==-1) {
        insertAt = lLength-1;
    } else {
        //insertAt = insertAt>=lLength?lLength:insertAt;
        insertAt = insertAt>=lLength?lLength-1:insertAt;
        long     moveThisMany = (laLength-insertAt-1);
        if (moveThisMany < 32L)
            for (long k=insertAt+moveThisMany; k> insertAt ; k--) {
                list_data[k] = list_data[k-1];
            }
        else {
            memmove (((char**)list_data)+(insertAt+1), ((char**)list_data)+insertAt, moveThisMany*sizeof(void*));
        }

      
    }
    if (store) {
      ((BaseRef*)list_data)[insertAt]=br->makeDynamic();
    } else {
      ((BaseRef*)list_data)[insertAt]=br;
      if (pointer) {
        br->AddAReference();
      }
    }


}

//Convert a list into a partition style string
BaseRef _SimpleList::ListToPartitionString () const {
    _StringBuffer *result = new _StringBuffer ((unsigned long)64);
  
    for (unsigned long k=0UL; k<lLength; k++) {
        unsigned long m;
        for (m=k+1UL; m<lLength; m++)
            if (list_data[m]-list_data[m-1]!=1) {
                break;
            }
        if (m>k+2UL) {
            (*result) << _String (list_data[k])
                      << '-'
                      << _String (list_data[m-1UL]);

            if (m<lLength) {
                (*result) << ',';
            }
            k = m-1UL;
        } else {
            (*result) << _String (list_data[k]);
            if (k<lLength-1) {
                (*result) << ',';
            }
        }
    }
    return result;
}

BaseRef _SimpleList::makeDynamic(void) const {
    _SimpleList * Res = new _SimpleList;
    Res->Duplicate (this);
    return Res;
}

long _SimpleList::Max(void) const{
    long res = LONG_MIN;
    for  (long e = 0L; e < lLength; e++)
        if (list_data[e] > res) {
            res = list_data[e];
        }
    return res;
}

long _SimpleList::Min(void) const {
    long res = LONG_MAX;
    for  (long e = 0L; e < lLength; e++)
        if (list_data[e] < res) {
            res = list_data[e];
        }
    return res;
}

//Merge 2 lists (sorted)

void _SimpleList::Merge(_SimpleList& l1, _SimpleList& l2, _SimpleList* mergeResults1, _SimpleList* mergeResults2) {
  
  this->Clear();
  if (mergeResults1) {
    mergeResults1->Clear();
  }
  if (mergeResults2) {
    mergeResults2->Clear();
  }

  enum    machine_states {INIT, ADVANCE1, ADVANCE2, FLUSH1, FLUSH2} advancing = INIT;

  unsigned long   pos1 = 0,
                  pos2 = 0,
                  nt1 = l1.countitems(),
                  nt2 = l2.countitems();


  bool  keep_going = true;

  while (keep_going) {
    switch (advancing) {
        
      case ADVANCE1: { // advancing in the 1st list
        pos1++;
        if (pos1==nt1) {
          advancing = FLUSH2;
          continue;
        }
        
        long cmp = l1.list_data[pos1] - l2.list_data [pos2];
        //if (l1lData[pos1] <= l2->list_data[pos2]) {
        if ( cmp <= 0L) {
          if (mergeResults1) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<this->lLength;
          }
          
          (*this) << (l1.list_data[pos1]);
          
          //        if (mergeResults2 && l1->list_data[pos1] < l2->list_data[pos2]) {
          if (cmp < 0L) {
            if (mergeResults2 && mergeResults2->countitems () == 0UL && pos1 >= pos2) {
              for (unsigned long i=0UL; i<pos2; i++) {
                (*mergeResults2) << i;
              }
            }
            continue;
          }
        }
        
        if (cmp > 0L) {
          //if (l1->list_data[pos1] > l2->list_data[pos2]) {
          advancing = ADVANCE2;
          if (mergeResults1 && mergeResults1->countitems () == 0UL ) {
            for (unsigned long i=0UL; i<pos1; i++) {
              //printf("MERGE line number %d in file %s\n", __xLINE__, __FILE__);
              (*mergeResults1) << i;
            }
          }
          
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength;
          }
          
          //printf ("this->append (l2->AtIndex(pos2))\n");
          (*this)  << l2.list_data[pos2];
          continue;
        }
        break;
      }
        
      case ADVANCE2: { // advancing in the 2nd list
        pos2++;
        if (pos2==nt2) {
          advancing = FLUSH1;
          continue;
        }
        
        long cmp = l2.list_data[pos2] - l1.list_data[pos1];
        
        //if (l2->list_data[pos2] <= l1->list_data[pos1]) {
        if (cmp <= 0L) {
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength;
          }
          
          (*this) << l2.list_data[pos2];
          
          //if (l2->list_data[pos2] < l1->list_data[pos1] && mergeResults1 && !doMerge1 && pos2>=pos1) {
          if (cmp < 0L ) {
            if(mergeResults1 && mergeResults1->countitems() == 0 && pos2>=pos1) {
              for (unsigned long i=0UL; i<pos1; i++) {
                //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
                (*mergeResults1)<<i;
              }
            }
            continue;
          }
        }
        
        //if (l2->list_data[pos2] > l1->list_data[pos1]) {
        if (cmp > 0L) {
          advancing = ADVANCE1;
          if (mergeResults2 && mergeResults2->countitems() == 0) {
            for (unsigned long i=0UL; i<pos2; i++) {
              (*mergeResults2)<<i;
            }
          }
          if (mergeResults1) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<this->lLength;
          }
          (*this) << l1.list_data[pos1];
          continue;
        }
        break;
      }
      case FLUSH2: { // flush out the 2nd list
        if (mergeResults1 && pos2<nt2 ) {
          for (unsigned long i=pos1; i<nt1; i++) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<i;
          }
        }
        if (mergeResults2)
          while (pos2<nt2) {
            (*mergeResults2)<<this->lLength;
            (*this) << (l2.list_data[pos2++]);
          }
        else
          while (pos2<nt2) {
            (*this) << (l2.list_data[pos2++]);
          }
        keep_going = false;
        break;
      }
      case FLUSH1: { // flush out the 1st list
        if (mergeResults2 && pos1<nt1) {
          for (unsigned long i=pos2; i<nt2; i++) {
            (*mergeResults2)<<i;
          }
        }
        if (mergeResults1)
          while (pos1<nt1) {
            //printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
            (*mergeResults1)<<this->lLength;
            (*this) << (l1.list_data[pos1++]);
          }
        else
          while (pos1<nt1) {
            (*this) << (l1.list_data[pos1++]);
          }
        keep_going = false;
        break;
      }
      case INIT: { // just starting
        if (!nt1) { // first list is empty!
          advancing = FLUSH2;
          continue;
        }
        if (!nt2) { // second list is empty!
          advancing = FLUSH1;
          continue;
        }
        
        //if (l1->list_data[pos1] <= l2->list_data[pos2]) { // begin with the first list
        
        long cmp = l1.list_data[pos1] - l2.list_data[pos2];
        
        if (cmp <= 0L) { // begin with the first list
          if (mergeResults1) {
            (*mergeResults1)<<this->lLength;
          }
          (*this) << (l1.list_data[0]);
          advancing = ADVANCE1;
          if (cmp != 0L) {
            continue;
          }
        } else {
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength;
          }
          (*this) << (l2.list_data[0]);
          advancing = ADVANCE2;
          continue;
        }
        break;
      }
    } // end SWITCH
    
    if (keep_going) {
      if (advancing == ADVANCE1) { // moving up in the second term
        pos1++;
        if (pos1==nt1) {
          pos2++;
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength-1UL;
          }
          advancing = FLUSH2;
          continue;
        } else {
          advancing = ADVANCE2;
          if (mergeResults2) {
            (*mergeResults2)<<this->lLength-1UL;
          }
        }
      } else {
        pos2++;
        if (pos2==nt2) {
          pos1++;
          if (mergeResults1) {
            (*mergeResults1)<<this->lLength-1UL;
          }
          advancing = FLUSH1;
          continue;
        } else {
          advancing = ADVANCE1;
          if (mergeResults1) {
            (*mergeResults1)<<this->lLength-1UL;
          }
        }
      }
    }
  }
}

// Together with the next function
// Implements algorithm NEXKSB from p.27 of http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
bool  _SimpleList::NChooseKInit (_SimpleList& state, _SimpleList& store, unsigned long stride, bool algorithm) {
    if (stride <= lLength && lLength) {
        state.Clear();
        state.RequestSpace (stride+3);
        state << stride;
        store.Clear();
        store.RequestSpace (stride);
        return true;
    }
    return false;
}

// Implements algorithm NEXKSB from p.27 of http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
bool  _SimpleList::NChooseK (_SimpleList& state, _SimpleList& store) {
    if (state.lLength == 1) { // first pass
        state << 0;              // m
        state << state.list_data[0]; // h
        state.lLength = state.list_data[0]+3;
        store.lLength = state.list_data[0];
        if (store.lLength == 0) {
            return false;
        }
    } else {
        if (state.list_data[1] < lLength - state.list_data[2]) {
            state.list_data[2] = 0;
        }
        state.list_data[2] ++;
        state.list_data[1] = state.list_data[3 + state.list_data[0] - state.list_data[2]] + 1;
    }
    for (long j=1; j <= state.list_data[2]; j++) {
        long  anIndex   = j+state.list_data[0]-state.list_data[2],
              anIndex2    = state.list_data[1]+j;
        state.list_data[anIndex+2]   = anIndex2-1;
        store.list_data[anIndex-1]   = list_data[anIndex2-1];
    }
    return state.list_data[3] < lLength-state.list_data[0];
}

//Coordinate normalizer
void _SimpleList::NormalizeCoordinates (long& from, long& to, const unsigned long refLength) {
    if (to < 0) {
        to = refLength+to;
    } else {
        to = to < refLength-1 ? to : refLength - 1;
    }
    if (from < 0) {
        from = refLength+from;
    }
}

void _SimpleList::Offset (long shift) {
    for (unsigned long k=0UL; k<lLength; k++) {
        list_data[k] += shift;
    }
}

_SimpleList* _SimpleList::Subset (unsigned long size, bool replacement) {
    _SimpleList* result = new _SimpleList;
    if (size > 0) {
        size = MIN(size, lLength);
        if (replacement) {
            for (long k = 0; k < size; k++) {
                (*result) << list_data[genrand_int32()%lLength];
            }
        } else {
            (*result) << (*this);
            for (long k = 0; k < size; k++) {
                long idx = list_data[genrand_int32()%(lLength-k)];
                long t = result->list_data[k];
                result->list_data[k] = result->list_data[idx];
                result->list_data[idx] = t;
            }
            result->lLength = size;
            result->TrimMemory();
        }
    }
    return result;
}

// Create a permutation of the list's elements
void  _SimpleList::Permute (long blockLength) {
    
    unsigned long blockCount = lLength/blockLength;

    if (blockCount > 1) {
        if (blockLength>1) {
            for (unsigned long k=0; k<blockCount-1; k=k+1) {
                unsigned long k2 = genrand_real2()*(blockCount-k);
                if (k2) {
                    k2 += k;
                    k2 *= blockLength;

                    for (long j = 0; j<blockLength; j++) {
                        Exchange(list_data[k2+j], list_data[k*blockLength+j]);
                    }
                }
            }

        } else {
            for (unsigned long k=0L; k<blockCount-1; k=k+1) {
                unsigned long k2 = genrand_real2()*(blockCount-k);
                if (k2) {
                    k2+=k;
                    Exchange (list_data[k2], list_data[k]);
                }
            }
        }
    }
}

  // Create a permutation of the list's elements
long  _SimpleList::Choice () const {
  if (lLength) {
    return genrand_int32() % lLength;
  }
  return kNotFound;
}

// Create a permutation of the list's elements
_SimpleList const  _SimpleList::Sample (unsigned long size) const {
    // TODO SLKP 20171026: this is new, need to check correctness
    if (size >= lLength) {
        return *this;
    }
    
    _SimpleList result (size, 0, 0),
                tracker (size, 0, 0);
    
    for (unsigned long k=0; k < size; k++) {
        unsigned long k2 = k + genrand_real2()*(lLength-k);
        result.list_data[k] = list_data[k2];
        tracker[k] = k2;
        if (k2 != k) {
            Exchange (list_data[k2], list_data[k]);
        }
    }
    
    // reshuffle moved elements LIFO to ensure correct ordering in *this
    
    for (long k = size - 1; k >= 0; k--) {
        Exchange (list_data[k], list_data[tracker.get(k)]);
    }
    
    return result;
}

// Create a permutation of the list's elements with possible repetitions
void  _SimpleList::PermuteWithReplacement (long blockLength)
{
    unsigned long blockCount = lLength/blockLength;
    _SimpleList   result ((unsigned long)(blockCount*blockLength));
    if (blockLength>1)
        for (long i = 0; i<blockCount; i++) {
            unsigned long sample = (unsigned long)(genrand_real2()*blockCount);
            sample *= blockLength;
            for (long j = 0; j<blockLength; j++,sample++) {
                result<<list_data[sample];
            }
        } else {
            for (long i = 0; i<blockCount; i++) {
                unsigned long sample = genrand_real2()*blockCount;
                result<<list_data[sample];
            }
        }

    Clear();
    Duplicate(&result);

}

long _SimpleList::Pop (unsigned long discard) {
    if (lLength > discard) {
        return list_data[lLength -= (discard+1UL)];
    }
    return 0L;
}

//Length constructor and populator
void _SimpleList::Populate (long l, long start, long step) {
    RequestSpace (l);
    for (long k = 0L; k < l; k++, start+=step) {
        list_data[k] = start;
    }
    lLength = l;
}


 void _SimpleList::AppendRange(unsigned long how_many, long start, long step)  {
  RequestSpace (how_many + lLength);
  for (unsigned long k = 0UL; k < how_many; k++, start+=step) {
    list_data [lLength + k] = start;
  }
  lLength = how_many + lLength;
}

void _SimpleList::RecursiveIndexSort (long from, long to, _SimpleList* index) {
    long middle = (from+to) >> 1, middleV = list_data[middle],
         bottommove = 1, topmove = 1, temp,i, imiddleV = (*index)(middle);

    if (middle)
        while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle) != kCompareLess)) {
            bottommove++;
        }

    if (from<to)
        while ((middle+topmove<=to)&&(Compare(middle+topmove,middle) != kCompareGreater)) {
            topmove++;
        }

    // now shuffle
    for (i=from; i<middle-bottommove; i++) {
        if (Compare(i,middle) != kCompareLess) {
            Exchange (list_data[middle-bottommove], list_data[i]);
            Exchange (index->list_data[middle-bottommove], index->list_data[i]);
            bottommove++;
            while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle)!= kCompareLess)) {
                bottommove++;
            }
        }
    }

    for (i=middle+topmove+1; i<=to; i++) {
        if (Compare(i,middle)!= kCompareGreater) {
            Exchange (list_data[middle+topmove], list_data[i]);
            Exchange (index->list_data[middle+topmove], index->list_data[i]);
            topmove++;
            while ((middle+topmove<=to)&&(Compare(middle+topmove,middle)!= kCompareGreater)) {
                topmove++;
            }
        }
    }

    if (topmove==bottommove) {
        for (i=1; i<bottommove; i++) {
            Exchange (list_data[middle+i], list_data[middle-i]);
            Exchange (index->list_data[middle+i], index->list_data[middle-i]);
        }
    } else if (topmove>bottommove) {
        long shift = topmove-bottommove;
        for (i=1; i<bottommove; i++) {
            Exchange (list_data[middle+i+shift], list_data[middle-i]);
            Exchange (index->list_data[middle+i+shift], index->list_data[middle-i]);
        }
        
        for (i=0; i<shift; i++) {
            list_data[middle+i]=list_data[middle+i+1];
            index->list_data [middle+i]= index->list_data [middle+i+1];
        }
        middle+=shift;
        list_data[middle]=middleV;
        index->list_data[middle]=imiddleV;
    } else {
        long shift = bottommove-topmove;
        for (i=1; i<topmove; i++) {
            Exchange (list_data[middle-i-shift], list_data[middle+i]);
            Exchange (index->list_data[middle-i-shift], index->list_data[middle+i]);
        }
        for (i=0; i<shift; i++) {
            list_data[middle-i]=list_data[middle-i-1];
            index->list_data[middle-i]=index->list_data[middle-i-1];
        }
        middle-=shift;
        list_data[middle]=middleV;
        index->list_data[middle]=imiddleV;
    }
    
    if (to>middle+1) {
        RecursiveIndexSort (middle+1,to, index);
    }
    if (from<middle-1) {
        RecursiveIndexSort (from,middle-1, index);
    }
}

//Append & store operator
void _SimpleList::RequestSpace (long slots) {
    if (slots>laLength) {
        laLength=(slots/MEMORYSTEP+1)*MEMORYSTEP;
        _EnsureCorrectStorageType ();
     }
}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void    _SimpleList::Subtract (_SimpleList const& l1, _SimpleList const& l2) {
    if (lLength) {
        Clear();
    }

    long  c1 = 0L,
          c2 = 0L;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (c1<l1.lLength && l1.list_data[c1]<l2.list_data[c2]) {
            (*this) << l1.list_data[c1++];
        }
        if (c1==l1.lLength) {
            break;
        }
        while ( c1<l1.lLength && c2<l2.lLength  && l1.list_data[c1]==l2.list_data[c2] ) {
            c1++;
            c2++;
        }
        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }
        while (c2<l2.lLength && l2.list_data[c2]<l1.list_data[c1]) {
            c2++;
        }
    }

    while (c1<l1.lLength) {
        (*this) << l1.list_data[c1++];
    }
}

void _SimpleList::Swap (long i, long j) {
    if ( i>=lLength || j>=lLength ) {
        return;
    }
    Exchange(list_data[i], list_data[j]);
}

//Char* conversion
BaseRef _SimpleList::toStr(unsigned long) {
    if (lLength) {
        unsigned long ma  = lLength*(1+log10((double)lLength));

        _StringBuffer * s = new _StringBuffer (MAX(32UL,ma));

        (*s) << "{";

        char c[32];
        for (unsigned long i = 0UL; i<lLength-1L; i++) {
            snprintf (c, sizeof(c),"%ld",list_data[i]);
            (*s) << c << ',';
        }
        snprintf (c, sizeof(c),"%ld",list_data[lLength-1L]);
        (*s) << c << '}';

        return s;
    } else {
        return new _String ("{}");
    }
}


//Delete item at index (>=0)
void  _SimpleList::TrimMemory (void) {
    if (laLength>lLength) {
        laLength = Maximum(lLength, (unsigned long)MEMORYSTEP);
        _EnsureCorrectStorageType();
    }
}

/*
==============================================================
Sort Methods
==============================================================
*/

void  _SimpleList::Sort (bool ascending) {
    if (lLength > 0UL) {
      if (lLength<10UL) { // use bubble sort
          BubbleSort();
      } else {
          QuickSort(0UL,lLength-1UL);
      }

      if (!ascending) {
          for (unsigned long i=0UL, j = lLength - 1UL; i<j; i++,j--) {
              Exchange(list_data[i], list_data[j]);
          }
      }
    }
}

void  _SimpleList::BubbleSort (void) {
    bool done = lLength == 0UL;
    while (!done) {
        done = true;
        unsigned long upper_bound = lLength - 1L;
        for (unsigned long i=0UL; i < upper_bound; i++) {
          if (Compare(i,i+1) == kCompareGreater) {
            done = false;
            Exchange (list_data[i], list_data[i+1]);
            upper_bound = i;
          }
        }
    }
}

void  _SimpleList::QuickSort (long from, long to)
{
    long middle = (from+to) >> 1,
         middleV = ((long*)list_data)[middle],
         top = to,
         bottommove = 1,
         topmove = 1,
          i;

    if (middle)
        //while ((middle-bottommove>=from)&&(((long*)list_data)[middle-bottommove]>middleV))
        while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle) == kCompareGreater)) {
            bottommove++;
        }

    if (from<to)
        //while ((middle+topmove<=to)&&(((long*)list_data)[middle+topmove]<middleV))
        while ((middle+topmove<=to)&&(Compare(middle+topmove,middle) == kCompareLess)) {
            topmove++;
        }
    // now shuffle
    for (i=from; i<middle-bottommove; i++) {
        if (Compare(i,middle) == kCompareGreater) {
            Exchange (list_data[middle-bottommove],list_data[i]);
            bottommove++;

            //while ((middle-bottommove>=from)&&(((long*)list_data)[middle-bottommove]>middleV))
            while ((middle-bottommove>=from)&&(Compare(middle-bottommove,middle) == kCompareGreater)) {
                bottommove++;
            }
        }
    }

    for (i=middle+topmove+1; i<=top; i++) {
        if (Compare(i,middle) == kCompareLess) {
            Exchange (list_data[middle+topmove],list_data[i]);
            topmove++;

            //while ((middle+topmove<=to)&&(((long*)list_data)[middle+topmove]<middleV))
            while ((middle+topmove<=to)&&(Compare(middle+topmove,middle) == kCompareLess)) {
                topmove++;
            }
        }
    }

    if (topmove==bottommove) {
        for (i=1; i<bottommove; i++) {
            Exchange (list_data[middle+i],list_data[middle-i]);
        }
    } else if (topmove>bottommove) {
        long shift = topmove-bottommove;
        for (i=1; i<bottommove; i++) {
          Exchange (list_data[middle+i+shift],list_data[middle-i]);
        }
        for (i=0; i<shift; i++) {
             list_data[middle+i] = list_data[middle+i+1];
        }
        middle+=shift;
        list_data[middle]=middleV;
    } else {
        long shift = bottommove-topmove;
        for (i=1; i<topmove; i++) {
            Exchange (list_data[middle-i-shift],list_data[middle+i]);
        }
        for (i=0; i<shift; i++) {
            list_data[middle-i] = list_data[middle-i-1];
        }
        middle-=shift;
        list_data[middle]=middleV;
    }
    if (to>middle+1) {
        QuickSort (middle+1,top);
    }
    if (from<middle-1) {
        QuickSort (from,middle-1);
    }
}

//TODO: This is a global. Should it be here?
void SortLists (_SimpleList* ref, _SimpleList* index) {
  if (ref->lLength > 0UL) {
      if (ref->lLength!=index->lLength) {
          return;
      }
      if ((*ref).lLength<=10UL) {
          bool done = false;

 
          while (!done) {
              done = true;
              unsigned long upper_bound = ref->lLength - 1L;
              for (unsigned long i=0UL; i < upper_bound; i++) {
                if (ref->Compare(i,i+1) == kCompareGreater) {
                  done = false;
                  Exchange(ref->list_data[i], ref->list_data[i+1]);
                  Exchange(index->list_data[i], index->list_data[i+1]);
                  upper_bound = i;
                }
              }
          }
      } else {
          (*ref).RecursiveIndexSort (0, (*ref).lLength-1UL,index);
      }
  }
}

/*
==============================================================
Set Methods
==============================================================
*/

// Compute the intersection of two sorted lists
void    _SimpleList::Intersect (_SimpleList& l1, _SimpleList& l2) {
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (l1.list_data[c1]<l2.list_data[c2]) {
            c1++;
            if (c1==l1.lLength) {
                break;
            }
        }
        if (c1==l1.lLength) {
            break;
        }

        while (l1.list_data[c1]==l2.list_data[c2]) {
            (*this) << l1.list_data[c1++];
            c2++;
            if (c1==l1.lLength || c2==l2.lLength) {
                break;
            }
        }
        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }
        while (l2.list_data[c2]<l1.list_data[c1]) {
            c2++;
            if (c2==l2.lLength) {
                break;
            }
        }
    }
}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void    _SimpleList::Union (_SimpleList& l1, _SimpleList& l2) {
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (l1.list_data[c1]<l2.list_data[c2]) {
            (*this) << l1.list_data[c1++];
            if (c1==l1.lLength) {
                break;
            }
        }

        if (c1==l1.lLength) {
            break;
        }

        while (l1.list_data[c1]==l2.list_data[c2]) {
            (*this) << l1.list_data[c1++];
            c2++;
            if (c1==l1.lLength || c2==l2.lLength) {
                break;
            }
        }

        if (c1==l1.lLength || c2==l2.lLength) {
            break;
        }

        while (l2.list_data[c2]<l1.list_data[c1]) {
            (*this) << l2.list_data[c2++];
            if (c2==l2.lLength) {
                break;
            }
        }
    }

    while (c1<l1.lLength) {
        (*this) << l1.list_data[c1++];
    }
    while (c2<l2.lLength) {
        (*this) << l2.list_data[c2++];
    }

}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void    _SimpleList::XOR (_SimpleList& l1, _SimpleList& l2)
{
    if (lLength) {
        Clear();
    }

    long  c1 = 0,
          c2 = 0;

    while (c1<l1.lLength && c2<l2.lLength) {
        while (c1<l1.lLength && l1.list_data[c1]<l2.list_data[c2]) {
            (*this) << l1.list_data[c1++];
        }
        if (c1==l1.lLength) {
            break;
        }
        while (c1<l1.lLength && c2<l2.lLength && l1.list_data[c1]==l2.list_data[c2] ) {
            c1++;
            c2++;
        }
        if (c1==l1.lLength||c2==l2.lLength) {
            break;
        }

        while (c2<l2.lLength && l2.list_data[c2]<l1.list_data[c1]) {
            (*this) << l2.list_data[c2++];
        }
    }

    while (c1<l1.lLength) {
        (*this) << l1.list_data[c1++];
    }
    while (c2<l2.lLength) {
        (*this) << l2.list_data[c2++];
    }
}

long _SimpleList::Map (long index, long map_failed) const {
  return index >= 0L && index < lLength ? list_data[index]: map_failed;
}
