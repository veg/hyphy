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

#include "hy_strings.h"
#include "errorfns.h"
#include "list.h"
#include "simplelist.h"
#include "legacy_parser.h"
#include "defines.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

/*
==============================================================
Constructors
==============================================================
*/

// Does nothing
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::_SimpleList() {
    Initialize(false);
}

//Data constructor (1 member list)
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::_SimpleList(const PAYLOAD item) {
  lLength     = 1UL;
  laLength    = HY_LIST_ALLOCATION_CHUNK;
  lData       = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
  lData       = item;
}

//Length constructor and populator
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::_SimpleList(const unsigned long l, const PAYLOAD start, const PAYLOAD step ){
  Initialize  (false);
  Populate    (l, start, step);
}

//Length constructor
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::_SimpleList(unsigned long l) {
  lLength  = 0UL;
  laLength = (l / HY_LIST_ALLOCATION_CHUNK + 1) * HY_LIST_ALLOCATION_CHUNK;
  lData    = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
  memset     (lData, 0, laLength * sizeof (PAYLOAD));
}

//Stack copy contructor
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::_SimpleList(const _SimpleList <PAYLOAD> &l, const long from, const long to) {
  if (from == 0 && to == HY_LIST_INSERT_AT_END) { // copy the whole thing
    Duplicate(&l);
  } else {
    Initialize();
    NormalizeCoordinates(from, to, l.lLength);
    RequestSpace(to - from);
    long upto = to - from;
    for (long k = 0; k < upto; k++) {
      lData[k] = l.lData[from + k];
    }
  }
}

// Data constructor (variable number of long constants)
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::_SimpleList(const PAYLOAD value1, const unsigned long number, ...) {
  Initialize(true);
  va_list vl;

  append(value1);

  va_start(vl, number);
  for (unsigned long arg_id = 0; arg_id < number; arg_id++) {
    const PAYLOAD this_arg = (PAYLOAD)va_arg(vl, PAYLOAD);
    append(this_arg);
  }
  va_end(vl);
}

//Destructor
template<typename PAYLOAD>
_SimpleList<PAYLOAD>::~_SimpleList(void) {
  if (CanFreeMe()) {
    if (lData) {
      free(lData);
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
template<typename PAYLOAD>
PAYLOAD &_SimpleList<PAYLOAD>::operator[](const long i) {
  if (lLength == 0) {
    return lData[0];
  }

  const unsigned long in = (const unsigned long) i;
  if (in > lLength - 1) {
    return lData[lLength - 1];
  }

  return lData[in];
}
//Element location functions (0,llength - 1)
template<typename PAYLOAD>
PAYLOAD _SimpleList<PAYLOAD>::operator()(const unsigned long i) {
  if (i < lLength) {
    return lData[i];
  }
  WarnError("List index out of range");
  return lData[0];
}

//Assignment operator
template<typename PAYLOAD>
_SimpleList<PAYLOAD> const _SimpleList<PAYLOAD>::operator=(const _SimpleList<PAYLOAD> l) {
  Clear();
  lLength  = l.lLength;
  laLength = l.laLength;
  if (laLength) {
    lData = (long *)MemAllocate(laLength * sizeof(PAYLOAD));
    if (lLength) {
      memcpy(lData, l.lData, lLength * sizeof(PAYLOAD));
    }
  }

  return *this;
}

//Append operator
template<typename PAYLOAD>
const _SimpleList<PAYLOAD> _SimpleList<PAYLOAD>::operator&(const _SimpleList<PAYLOAD> l) {
  _SimpleList<PAYLOAD> res(l.lLength + lLength);
  if (res.laLength == 0UL) {
    return res;
  }

  if (lData && lLength) {
    memcpy(res.lData, lData, lLength * sizeof(PAYLOAD));
  }

  if (l.lData && l.lLength) {
    memcpy((char *)res.lData + lLength * sizeof(PAYLOAD), l.lData,
           l.lLength * sizeof(PAYLOAD));
  }

  res.lLength = l.lLength + lLength;

  return res;
}


template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::append(const PAYLOAD item) {
  InsertElement(item, HY_LIST_INSERT_AT_END);
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::operator<<(const PAYLOAD item) {
  append (item);
}

template<typename PAYLOAD>
bool _SimpleList<PAYLOAD>::operator>>(const PAYLOAD item) {
  if (Find(item) == HY_NOT_FOUND) {
    InsertElement(item, HY_LIST_INSERT_AT_END);
    return true;
  }
  return false;
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::operator<<(const _SimpleList<PAYLOAD> &source) {
  for (unsigned long k = 0; k < source.lLength; k++) {
    append(source.GetElement(k));
  }
}

/*
==============================================================
Methods
==============================================================
*/

//Element location functions (0,llength - 1), negative values return
// elements from the end of the list

template<typename PAYLOAD>
PAYLOAD _SimpleList<PAYLOAD>::GetElement(const long index) const{
  if (index >= 0L) {
    if ((const unsigned long) index < lLength) {
      return lData[index];
    }
  }
  if ((const unsigned long)(-index) <= lLength) {
    return lData[lLength + index];
  }
  WarnError(_String("List index '") & (long)((const unsigned long)(-index)) &
            "' out of range in _SimpleList::GetElement on list of length " &
            long(lLength));
  
  return lData[0UL];
}

template<typename PAYLOAD>
long _SimpleList<PAYLOAD>::BinaryFind(const PAYLOAD item, const long startAt) const {
  long top = lLength - 1L,
       bottom = startAt,
       middle;

  if (top == -1L) {
    return -2L;
  }

  while (top > bottom) {
    middle = (top + bottom) >> 1;
    if (item < GetElement (middle)) {
      top = middle == top ? top - 1L : middle;
    } else if (item > GetElement (middle)) {
      bottom = middle == bottom ? bottom + 1 : middle;
    } else {
      return middle;
    }
  }

  middle = top;
  if (GetElement (middle) == item) {
    return middle;
  } else {
    if (GetElement (middle) < item) {
      return -middle - 3;
    }
  }

  return -middle - 2;
}

template<typename PAYLOAD>
long _SimpleList<PAYLOAD>::BinaryInsert(const PAYLOAD item) {
  if (lLength == 0UL) {
    append (item);
    return 0L;
  }

  long pos = -BinaryFind(item) - 2;

  if (pos < 0L) {
    return -pos + 2;
  }

  if (lData[pos] < item) {
    pos++;
  }

  InsertElement(item, pos);

  return pos >= lLength ? lLength - 1 : pos;
}


template<typename PAYLOAD>
PAYLOAD _SimpleList<PAYLOAD>::Sum(void) const {
  PAYLOAD sum = 0;
  for (unsigned long k = 0; k < lLength; k++) {
    sum += lData[k];
  }
  return sum;
}

template<typename PAYLOAD>
long _SimpleList<PAYLOAD>::Compare(const long i, const long j) {
  PAYLOAD v1 = GetElement (i), v2 = GetElement (j);

  if (v1 < v2) {
    return -1;
  } else if (v1 == v2) {
    return 0;
  }
  return 1;
}


// Compute the number of shared of two sorted lists
template<typename PAYLOAD>
long _SimpleList<PAYLOAD>::CountCommonElements(const _SimpleList &l1, bool at_least_one) const {
  long c1 = 0L, c2 = 0L, res = 0L;

  while (c1 < l1.lLength && c2 < lLength) {
    while (l1.lData[c1] < lData[c2]) {
      c1++;
      if (c1 == l1.lLength) {
        break;
      }
    }
    if (c1 == l1.lLength) {
      break;
    }

    while (l1.lData[c1] == lData[c2]) {
      c2++;
      if (at_least_one) {
        return 1;
      } else {
        res++;
      }
      if (c1 == l1.lLength || c2 == lLength) {
        break;
      }
    }
    if (c1 == l1.lLength || c2 == lLength) {
      break;
    }
    while (lData[c2] < l1.lData[c1]) {
      c2++;
      if (c2 == lLength) {
        break;
      }
    }
  }
  return res;
}

//List length
template<typename PAYLOAD>
unsigned long _SimpleList<PAYLOAD>::countitems(void) const { return lLength; }


template<typename PAYLOAD>
_SimpleList<PAYLOAD> *_SimpleList<PAYLOAD>::CountingSort(PAYLOAD, _SimpleList<long>*) {
  return new _SimpleList <PAYLOAD>;
}

template<>
_SimpleList<long> *_SimpleList<long>::CountingSort(long upperBound, _SimpleList<long> *ordering) {
  if (ordering) {
    ordering->Clear();
  }

  if (lLength) {
    if (upperBound < 0) {
      upperBound = Max() + 1;
    }

    _SimpleList<long> buffer,
               *result = new _SimpleList<long>(upperBound);
    
    buffer.Populate(upperBound, 0L, 0L);
    
    for (long pass1 = 0; pass1 < lLength; pass1++) {
      buffer.lData[lData[pass1]]++;
    }
    for (long pass2 = 1; pass2 < upperBound; pass2++) {
      buffer.lData[pass2] += buffer.lData[pass2 - 1];
    }
    if (ordering) {
      ordering->Populate(lLength, 0, 0);
      for (long pass3 = lLength - 1; pass3 >= 0; pass3--) {
        result->lData[--buffer.lData[lData[pass3]]] = lData[pass3];
        ordering->lData[buffer.lData[lData[pass3]]] = pass3;
      }
    } else
      for (long pass3 = lLength - 1; pass3 >= 0; pass3--) {
        result->lData[--buffer.lData[lData[pass3]]] = lData[pass3];
      }
    result->lLength = lLength;

    return result;
  }
  return new _SimpleList<long>;
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Clear(bool completeClear) {
  if (CanFreeMe()) {
    lLength = 0;
    if (completeClear) {
      laLength = 0;
      if (lData) {
        free(lData);
      }
      lData = nil;
    }
  } else {
    RemoveAReference();
  }
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::CompactList(void) {
  if (laLength - lLength > HY_LIST_ALLOCATION_CHUNK) {
    laLength -= ((laLength - lLength) / HY_LIST_ALLOCATION_CHUNK) * HY_LIST_ALLOCATION_CHUNK;
    if (laLength) {
      lData = (PAYLOAD *)MemReallocate(lData, laLength * sizeof(PAYLOAD));
    } else {
      free(lData);
      lData = nil;
    }
  }
 
}


//Delete item at index (>=0)
template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Delete(const long index, bool compact_list) {
  if (index >= 0L && (const unsigned long)index < lLength) {
    lLength--;
    if (lLength - index) {
      memmove(lData + sizeof(PAYLOAD) * (index),
              lData + sizeof(PAYLOAD) * (index + 1),
              sizeof(PAYLOAD) * (lLength - index));
    }
  }
  if (compact_list) {
    CompactList ();
  }

}

//Delete duplicates from a sorted list
template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::DeleteDuplicates(void) {
  if (lLength > 1) {
    _SimpleList <PAYLOAD> noDups;

    PAYLOAD lastValue = lData[0];
    noDups << lastValue;
    
    for (unsigned long k = 1; k < lLength; k++) {
      PAYLOAD thisValue = lData[k];
      if (thisValue != lastValue) {
        noDups << thisValue;
        lastValue = thisValue;
      }
    }

    if (noDups.lLength != lLength) {
      Duplicate(&noDups);
    }
  }
}

//Delete items from a sorted list
template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::DeleteList(const _SimpleList<long> &indices_to_delete) {
  if (indices_to_delete.lLength) {
    unsigned long k = 0UL;
    for (unsigned long i = 0UL; i < lLength; i++) {
      if (k < indices_to_delete.lLength && i == indices_to_delete.lData[k]) {
        k++;
      } else {
        lData[i - k] = lData[i];
      }
    }
    lLength -= indices_to_delete.lLength;
  }
  CompactList();
}

//Shift the range from start to end
template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Displace(long start, long end, long delta) {
  if (start <= 0L) {
    start = 0L;
  } else if (start >= lLength) {
    start = lLength - 1L;
  }

  if (end < 0L) {
    end = lLength - 1L;
  } else if (end >= lLength) {
    end = lLength - 1L;
  }

  if (end - start >= 0L && delta && end - start < lLength - 1) {
    if (delta > 0 && lLength - end <= delta) { // shift up
      delta = lLength - end - 1;
    } else if (start - delta < 0) {
      delta = start;
    }

    if (delta) {
      long i, j, delta2 = end - start + 1;
      _SimpleList <PAYLOAD> swapList((unsigned long)(end - start + 1));

      for (i = start; i <= end; i++) {
        swapList << lData[i];
      }

      if (delta > 0) {

        for (i = end + 1; i <= end + delta; i++) {
          lData[i - delta2] = lData[i];
        }

      } else {
        for (i = start - 1; i >= start + delta; i--) {
          lData[i + delta2] = lData[i];
        }
      }

      for (i = start + delta, j = 0; i <= end + delta; i++, j++) {
        lData[i] = swapList.lData[j];
      }
    }
  }
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Duplicate(BaseRef theRef) {
  _SimpleList<PAYLOAD> *l = dynamic_cast<_SimpleList<PAYLOAD> >(theRef);
  lLength   = l->lLength;
  laLength  = l->laLength;
  lData     = l->lData;
  if (lData) {
    lData = (long *)MemAllocate(laLength * sizeof(PAYLOAD));
    memcpy(lData, l->lData, lLength * sizeof(PAYLOAD));
  }
}

//Element location functions (0,llength - 1)
//Negative indices return offsets from the end of the list
template<typename PAYLOAD>
PAYLOAD _SimpleList<PAYLOAD>::Element(const long index) {
  if (index >= 0L && index < lLength) {
    return lData[index];
  } else if (-index <= lLength) {
    return lData[lLength - (-index)];
  }
  return 0L;
}

template<typename PAYLOAD>
bool _SimpleList<PAYLOAD>::Equal(const _SimpleList<PAYLOAD> &l2) {
  if (lLength != l2.lLength) {
    return false;
  }

  for (unsigned long i = 0UL; i < lLength; i++)
    if (lData[i] != l2.lData[i]) {
      return false;
    }

  return true;
}

template<typename PAYLOAD>
long _SimpleList<PAYLOAD>::Find(const PAYLOAD item, const long startAt) const {
  for (unsigned long i = startAt; i < lLength; i++) {
    if (lData[i] == item) {
      return i;
    }
  }
  return HY_NOT_FOUND;
}

template<typename PAYLOAD>
long _SimpleList<PAYLOAD>::FindStepping(const PAYLOAD item, const long step,
                               const long startAt) const {
  for (unsigned long i = startAt; i < lLength; i += step)
    if (lData[i] == item) {
      return i;
    }

  return HY_NOT_FOUND;
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::FilterRange(const PAYLOAD lb, const PAYLOAD ub) {
  if (ub <= lb) {
    Clear();
  } else {
    _SimpleList <long> toDelete;
    for (unsigned long k = 0UL; k < lLength; k++)
      if (lData[k] <= lb || lData[k] >= ub) {
        toDelete << k;
      }
    DeleteList(toDelete);
  }
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Flip() {
  if (lLength > 0UL) {
    for (unsigned long k = 0UL, l = lLength - 1UL; k < l; k++, l--) {
      PAYLOAD pt = lData[k];
      lData[k] = lData[l];
      lData[l] = pt;
    }
  }
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Initialize(bool doMemAlloc) {
  BaseObj::Initialize();
  lLength = 0UL;
  if (doMemAlloc) {
    laLength = HY_LIST_ALLOCATION_CHUNK;
    lData    = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
  } else {
    laLength = 0;
    lData    = nil;
  }
}

template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::ResizeList(void) {
  if (lLength > laLength) {
    unsigned long incBy = (HY_LIST_ALLOCATION_CHUNK * 5UL > lLength) ? HY_LIST_ALLOCATION_CHUNK : lLength / 5UL;
    
    laLength += incBy;
    
    if (lData) {
      lData = (PAYLOAD *)MemReallocate(lData, laLength * sizeof(PAYLOAD));
    } else {
      lData = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
    }
    
    if (!lData) {
      checkPointer(lData);
    }
  }
}
//Append & store operator
template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::InsertElement(const PAYLOAD item, long insert_at) {
  
  lLength++;
  ResizeList();
  
  if (insert_at == HY_LIST_INSERT_AT_END) {
    lData[lLength-1UL] = item;
  } else {
    insert_at = insert_at >= lLength ? lLength - 1UL : insert_at;
    long moveThisMany = (laLength - insert_at - 1L);
    if (moveThisMany < 32L)
      for (long k = insert_at + moveThisMany; k > insert_at; k--) {
        lData[k] = lData[k - 1];
      }
    else {
      memmove(lData + (insert_at + 1), lData + insert_at,
              moveThisMany * sizeof(PAYLOAD));
    }
  }
}

//Convert a list into a partition style string
template<typename PAYLOAD>
_String* _SimpleList<PAYLOAD>::ListToPartitionString() {
  return new _String;
}
  
template<>
_String* _SimpleList<long>::ListToPartitionString() {
    _String *result = new _String((unsigned long) 64, true), conv;
    
    for (unsigned long k = 0; k < lLength; k++) {
      long m;
      for (m = k + 1; m < lLength; m++)
        if (lData[m] - lData[m - 1] != 1) {
          break;
        }
      if (m > k + 2) {
        conv = lData[k];
        (*result) << &conv;
        (*result) << '-';
        conv = lData[m - 1];
        (*result) << &conv;
        if (m < lLength) {
          (*result) << ',';
        }
        k = m - 1;
      } else {
        conv = lData[k];
        (*result) << &conv;
        if (k < lLength - 1) {
          (*result) << ',';
        }
      }
    }
    result->Finalize();
    return result;
}

template<typename PAYLOAD>
BaseRef _SimpleList<PAYLOAD>::makeDynamic(void) {
  _SimpleList <PAYLOAD> *res = new _SimpleList <PAYLOAD>;
  res->nInstances = 1;
  res->lData = nil;
  res->Duplicate(this);
  return res;
}

template<typename PAYLOAD>
PAYLOAD _SimpleList<PAYLOAD>::Max(void) const{
  PAYLOAD res = lData[0];
  for (unsigned long e = 1UL; e < lLength; e++)
    if (lData[e] > res) {
      res = lData[e];
    }
  return res;
}

template<typename PAYLOAD>
PAYLOAD _SimpleList<PAYLOAD>::Min(void) const {
  long res = lData[0];
  for (unsigned long e = 0; e < lLength; e++)
    if (lData[e] < res) {
      res = lData[e];
    }
  return res;
}

//Merge 2 lists (sorted)
template<typename PAYLOAD>
void _SimpleList<PAYLOAD>::Merge(_SimpleList<PAYLOAD> &l1, _SimpleList<PAYLOAD> &l2,
                        _SimpleList<long> *mergeResults1,
                        _SimpleList<long> *mergeResults2) {
  Clear();
  if (mergeResults1) {
    mergeResults1->Clear();
  }
  if (mergeResults2) {
    mergeResults2->Clear();
  }
  char advancing = -1;

  PAYLOAD * list1 = l1.lData,
          * list2 = l2.lData;
  
  long pos1 = 0,
       pos2 = 0,
       nt1 = l1.lLength,
       nt2 = l2.lLength;

  if (mergeResults1 && mergeResults2) {
    bool doMerge1 = false, doMerge2 = false;
    while (1) {             // stuff left to do
      if (advancing == 0) { // advancing in the 1st list
        pos1++;
        if (pos1 == nt1) {
          advancing = 2;
          continue;
        }
        list1++;
        if (*list1 <= *list2) {
          if (doMerge1) {
            (*mergeResults1) << lLength;
          }

          append (*list1);
          if (*list1 < *list2) {
            if (!doMerge2) {
              if (pos1 >= pos2) {
                doMerge2 = true;
                for (long i = 0L; i < pos2; i++) {
                  (*mergeResults2) << i;
                }
              }
            }
            continue;
          }
        }
        
        if (*list1 > *list2) {
          advancing = 1;
          if (!doMerge1) {
            for (long i = 0; i < pos1; i++) {
              (*mergeResults1) << i;
            }
            doMerge1 = true;
          }
          if (doMerge2) {
            (*mergeResults2) << lLength;
          }
          append(*list2);
          continue;
        }

      } else if (advancing == 1) { // advancing in the 2nd list
        pos2++;
        if (pos2 == nt2) {
          advancing = 3;
          continue;
        }
        list2++;
          //c = *list2 - *list1;
        if (*list2 <= *list1) {
          if (doMerge2) {
            (*mergeResults2) << lLength;
          }
          append (*list2);
          if (*list2 < *list1) {
            if (!doMerge1) {
              if (pos2 >= pos1) {
                doMerge1 = true;
                for (long i = 0; i < pos1; i++) {
                  (*mergeResults1) << i;
                }
              }
            }
            continue;
          }
        }
        if (c > 0) {
          advancing = 0;
          if (!doMerge2) {
            for (i = 0; i < pos2; i++) {
              (*mergeResults2) << i;
            }
            doMerge2 = true;
          }
          if (doMerge1) {
            (*mergeResults1) << lLength;
          }
          (*this) << *list1;
          continue;
        }
      } else if (advancing == 2) { // flush out the 2nd list
        if (!doMerge1 && (pos2 < nt2)) {
          for (i = 0; i < nt1; i++) {
            (*mergeResults1) << i;
          }
        }
        if (doMerge2)
          while (pos2 < nt2) {
            (*mergeResults2) << lLength;
            (*this) << *list2;
            list2++;
            pos2++;
          }
        else
          while (pos2 < nt2) {
            (*this) << *list2;
            list2++;
            pos2++;
          }
        break;
      } else if (advancing == 3) { // flush out the 1st list
        if (!doMerge2 && (pos1 < nt1)) {
          for (i = 0; i < nt2; i++) {
            (*mergeResults2) << i;
          }
        }
        if (doMerge1)
          while (pos1 < nt1) {
            (*mergeResults1) << lLength;
            (*this) << *list1;
            list1++;
            pos1++;
          }
        else
          while (pos1 < nt1) {
            (*this) << *list1;
            list1++;
            pos1++;
          }
        break;
      } else if (advancing == -1) { // just starting
        if (!nt1) {                 // first list is empty!
          advancing = 2;
          continue;
        }
        if (!nt2) { // second list is empty!
          advancing = 3;
          continue;
        }
        c = *list1 - *list2;
        if (c <= 0) { // begin with the first list
          (*this) << *list1;
          advancing = 0;
          if (c) {
            doMerge2 = true;
            continue;
          }
        } else {
          (*this) << *list2;
          advancing = 1;
          doMerge1 = true;
          continue;
        }

      }

      if (advancing == 0) { // moving up in the second term
        pos1++;
        if (pos1 == nt1) {
          list2++;
          pos2++;
          if (doMerge2) {
            (*mergeResults2) << lLength - 1;
          }
          advancing = 2;
          continue;
        } else {
          advancing = 1;
          if (doMerge2) {
            (*mergeResults2) << lLength - 1;
          }
          list1++;
        }
      } else {
        pos2++;
        if (pos2 == nt2) {
          list1++;
          pos1++;
          if (doMerge1) {
            (*mergeResults1) << lLength - 1;
          }
          advancing = 3;
          continue;
        } else {
          list2++;
          if (doMerge1) {
            (*mergeResults1) << lLength - 1;
          }
          advancing = 0;
        }
      }
    }
  } else {
    while (1) {             // stuff left to do
      if (advancing == 0) { // advancing in the 1st list
        pos1++;
        if (pos1 == nt1) {
          advancing = 2;
          continue;
        }
        list1++;
        c = *list1 - *list2;
        if (c <= 0) {
          (*this) << *list1;
          if (c < 0) {
            continue;
          }
        }
        if (c > 0) {
          advancing = 1;
          (*this) << *list2;
          continue;
        }

      } else if (advancing == 1) { // advancing in the 2nd list
        pos2++;
        if (pos2 == nt2) {
          advancing = 3;
          continue;
        }
        list2++;
        c = *list2 - *list1;
        if (c <= 0) {
          (*this) << *list2;
          if (c < 0) {
            continue;
          }
        }
        if (c > 0) {
          advancing = 0;
          (*this) << *list1;
          continue;
        }
      } else if (advancing == 2) { // flush out the 2nd list
        while (pos2 < nt2) {
          (*this) << *list2;
          list2++;
          pos2++;
        }
        break;
      } else if (advancing == 3) { // flush out the 2nd list
        while (pos1 < nt1) {
          (*this) << *list1;
          list1++;
          pos1++;
        }
        break;
      } else if (advancing == -1) { // just starting
        if (!nt1) {                 // first list is empty!
          advancing = 2;
          continue;
        }
        if (!nt2) { // second list is empty!
          advancing = 3;
          continue;
        }
        c = *list1 - *list2;
        if (c <= 0) { // begin with the first list
          (*this) << *list1;
          advancing = 0;
          if (c) {
            continue;
          }
        } else {
          (*this) << *list2;
          advancing = 1;
          continue;
        }

      }

      if (advancing == 0) { // moving up in the second term
        pos1++;
        if (pos1 == nt1) {
          list2++;
          pos2++;
          advancing = 2;
          continue;
        } else {
          advancing = 1;
          list1++;
        }
      } else {
        pos2++;
        if (pos2 == nt2) {
          list1++;
          pos1++;
          advancing = 3;
          continue;
        } else {
          list2++;
          advancing = 0;
        }
      }

    }
  }
}

// Together with the next function
// Implements algorithm NEXKSB from p.27 of
// http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
bool _SimpleList::NChooseKInit(_SimpleList &state, _SimpleList &store,
                               unsigned long stride, bool algorithm) {
  if (stride <= lLength && lLength) {
    state.Clear();
    state.RequestSpace(stride + 3);
    state << stride;
    store.Clear();
    store.RequestSpace(stride);
    return true;
  }
  return false;
}

// Implements algorithm NEXKSB from p.27 of
// http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
bool _SimpleList::NChooseK(_SimpleList &state, _SimpleList &store) {
  if (state.lLength == 1) {  // first pass
    state << 0;              // m
    state << state.lData[0]; // h
    state.lLength = state.lData[0] + 3;
    store.lLength = state.lData[0];
    if (store.lLength == 0) {
      return false;
    }
  } else {
    if (state.lData[1] < lLength - state.lData[2]) {
      state.lData[2] = 0;
    }
    state.lData[2]++;
    state.lData[1] = state.lData[3 + state.lData[0] - state.lData[2]] + 1;
  }
  for (long j = 1; j <= state.lData[2]; j++) {
    long anIndex = j + state.lData[0] - state.lData[2],
         anIndex2 = state.lData[1] + j;
    state.lData[anIndex + 2] = anIndex2 - 1;
    store.lData[anIndex - 1] = lData[anIndex2 - 1];
  }
  return state.lData[3] < lLength - state.lData[0];
}

//Coordinate normalizer
void _SimpleList::NormalizeCoordinates(long &from, long &to,
                                       const unsigned long refLength) {
  if (to < 0) {
    to = refLength + to;
  } else {
    to = to < refLength - 1 ? to : refLength - 1;
  }
  if (from < 0) {
    from = refLength + from;
  }
}

void _SimpleList::Offset(long shift) {
  if (lData) {
    for (long k = 0; k < lLength; k++) {
      lData[k] += shift;
    }
  }
}

_SimpleList *_SimpleList::Subset(unsigned long size, bool replacement) {
  _SimpleList *result = new _SimpleList;
  if (size > 0) {
    size = MIN(size, lLength);
    if (replacement) {
      for (long k = 0; k < size; k++) {
        (*result) << lData[genrand_int32() % lLength];
      }
    } else {
      (*result) << (*this);
      for (long k = 0; k < size; k++) {
        long idx = lData[genrand_int32() % (lLength - k)];
        long t = result->lData[k];
        result->lData[k] = result->lData[idx];
        result->lData[idx] = t;
      }
      result->lLength = size;
      result->TrimMemory();
    }
  }
  return result;
}

// Create a permutation of the list's elements
void _SimpleList::Permute(long blockLength) {
  unsigned long blockCount = lLength / blockLength;

  if (blockLength > 1) {
    /*_SimpleList     result ((unsigned long)(blockCount*blockLength));
    while (blockCount)
    {
        unsigned long sample = (unsigned long)(genrand_real2()*blockCount);
        sample *= blockLength;
        for (long j = 0; j<blockLength; j++)
        {
            result<<lData[sample];
            Delete(sample);
        }
        blockCount --;
    }
    Duplicate(&result); */

    for (unsigned long k = 0; k < blockCount - 1; k = k + 1) {
      unsigned long k2 = genrand_real2() * (blockCount - k);
      if (k2) {
        k2 += k;
        k2 *= blockLength;

        for (long j = 0; j < blockLength; j++) {
          long t = lData[k2 + j];
          lData[k2 + j] = lData[k * blockLength + j];
          lData[k * blockLength + j] = t;
        }
      }
    }

  } else {
    for (unsigned long k = 0; k < blockCount - 1; k = k + 1) {
      unsigned long k2 = genrand_real2() * (blockCount - k);
      if (k2) {
        k2 += k;
        long t = lData[k2];
        lData[k2] = lData[k];
        lData[k] = t;
      }
    }
    /*{
        while (blockCount)
        {
            unsigned long sample = genrand_real2()*blockCount;
            result<<lData[sample];
            Delete(sample);
            blockCount --;
        }
    }*/
  }
}

// Create a permutation of the list's elements with possible repetitions
void _SimpleList::PermuteWithReplacement(long blockLength) {
  unsigned long blockCount = lLength / blockLength;
  _SimpleList result((unsigned long)(blockCount * blockLength));
  if (blockLength > 1)
    for (long i = 0; i < blockCount; i++) {
      unsigned long sample = (unsigned long)(genrand_real2() * blockCount);
      sample *= blockLength;
      for (long j = 0; j < blockLength; j++, sample++) {
        result << lData[sample];
      }
    }
  else {
    for (long i = 0; i < blockCount; i++) {
      unsigned long sample = genrand_real2() * blockCount;
      result << lData[sample];
    }
  }

  Clear();
  Duplicate(&result);

}

long _SimpleList::Pop(void) {
  if (lLength > 0) {
    lLength--;
    return lData[lLength];
  }

  return 0;
}


//Length constructor and populator
void _SimpleList::Populate(long l, long start, long step) {
  RequestSpace(l);
  for (long k = 0; k < l; k++, start += step) {
    lData[k] = start;
  }

  lLength = l;
}

void _SimpleList::RecursiveIndexSort(long from, long to, _SimpleList *index) {
  long middle = (from + to) / 2, middleV = lData[middle], bottommove = 1,
       topmove = 1, temp, i, imiddleV = (*index)(middle);
  long *idata = (*index).quickArrayAccess();
  if (middle)
    while ((middle - bottommove >= from) &&
           (Compare(middle - bottommove, middle) >= 0)) {
      bottommove++;
    }

  if (from < to)
    while ((middle + topmove <= to) &&
           (Compare(middle + topmove, middle) <= 0)) {
      topmove++;
    }

  // now shuffle
  for (i = from; i < middle - bottommove; i++) {
    if (Compare(i, middle) >= 0) {
      temp = lData[middle - bottommove];
      lData[middle - bottommove] = lData[i];
      lData[i] = temp;
      temp = idata[middle - bottommove];
      idata[middle - bottommove] = idata[i];
      idata[i] = temp;
      bottommove++;
      while ((middle - bottommove >= from) &&
             (Compare(middle - bottommove, middle) >= 0)) {
        bottommove++;
      }
    }
  }

  for (i = middle + topmove + 1; i <= to; i++) {
    if (Compare(i, middle) <= 0) {
      temp = lData[middle + topmove];
      lData[middle + topmove] = lData[i];
      lData[i] = temp;
      temp = idata[middle + topmove];
      idata[middle + topmove] = idata[i];
      idata[i] = temp;
      topmove++;
      while ((middle + topmove <= to) &&
             (Compare(middle + topmove, middle) <= 0)) {
        topmove++;
      }
    }
  }

  if (topmove == bottommove) {
    for (i = 1; i < bottommove; i++) {
      temp = lData[middle + i];
      lData[middle + i] = lData[middle - i];
      lData[middle - i] = temp;
      temp = idata[middle + i];
      idata[middle + i] = idata[middle - i];
      idata[middle - i] = temp;
    }
  } else if (topmove > bottommove) {
    long shift = topmove - bottommove;
    for (i = 1; i < bottommove; i++) {
      temp = lData[middle + i + shift];
      lData[middle + i + shift] = lData[middle - i];
      lData[middle - i] = temp;
      temp = idata[middle + i + shift];
      idata[middle + i + shift] = idata[middle - i];
      idata[middle - i] = temp;
    }
    for (i = 0; i < shift; i++) {
      lData[middle + i] = lData[middle + i + 1];
      idata[middle + i] = idata[middle + i + 1];
    }
    middle += shift;
    lData[middle] = middleV;
    idata[middle] = imiddleV;
  } else {
    long shift = bottommove - topmove;
    for (i = 1; i < topmove; i++) {
      temp = lData[middle - i - shift];
      lData[middle - i - shift] = lData[middle + i];
      lData[middle + i] = temp;
      temp = idata[middle - i - shift];
      idata[middle - i - shift] = idata[middle + i];
      idata[middle + i] = temp;
    }
    for (i = 0; i < shift; i++) {
      lData[middle - i] = lData[middle - i - 1];
      idata[middle - i] = idata[middle - i - 1];
    }
    middle -= shift;
    lData[middle] = middleV;
    idata[middle] = imiddleV;
  }
  if (to > middle + 1) {
    RecursiveIndexSort(middle + 1, to, index);
  }
  if (from < middle - 1) {
    RecursiveIndexSort(from, middle - 1, index);
  }
}

//Append & store operator
void _SimpleList::RequestSpace(long slots) {
  if (slots > laLength) {
    laLength = (slots / MEMORYSTEP + 1) * MEMORYSTEP;
    if (lData) {
      checkPointer(lData = (long *)MemReallocate((char *)lData,
                                                 laLength * sizeof(void *)));
    } else {
      checkPointer(lData = (long *)MemAllocate(laLength * sizeof(void *)));
    }
  }
}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void _SimpleList::Subtract(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (c1 < l1.lLength && l1.lData[c1] < l2.lData[c2]) {
      (*this) << l1.lData[c1++];
    }
    if (c1 == l1.lLength) {
      break;
    }
    while (c1 < l1.lLength && c2 < l2.lLength && l1.lData[c1] == l2.lData[c2]) {
      c1++;
      c2++;
    }
    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }
    while (c2 < l2.lLength && l2.lData[c2] < l1.lData[c1]) {
      c2++;
    }
  }

  while (c1 < l1.lLength) {
    (*this) << l1.lData[c1++];
  }
}

void _SimpleList::Swap(long i, long j) {
  if (i >= lLength || j >= lLength) {
    return;
  }

  void *pt = ((void **)lData)[j];
  ((void **)lData)[j] = ((void **)lData)[i];
  ((void **)lData)[i] = pt;
}

//Char* conversion
BaseRef _SimpleList::toStr(void) {
  if (lLength) {
    unsigned long ssi = _String::storageIncrement,
                  ma = lLength * (1 + log10((double) lLength));

    if (ma > ssi) {
      _String::storageIncrement = ma;
    }

    _String *s = new _String(10L, true);

    checkPointer(s);

    (*s) << "{";

    for (unsigned long i = 0; i < lLength; i++) {
      char c[32];
      snprintf(c, sizeof(c), "%ld", ((long *)lData)[i]), (*s) << c;
      if (i < lLength - 1) {
        (*s) << ',';
      }
    }

    (*s) << '}';

    s->Finalize();
    _String::storageIncrement = ssi;
    return s;
  } else {
    return new _String("{}");
  }
}

//Delete item at index (>=0)
void _SimpleList::TrimMemory(void) {
  if (laLength > lLength) {
    laLength = lLength;
    if (laLength) {
      if (lData) {
        lData = (long *)MemReallocate((char *)lData, laLength * sizeof(Ptr));
      } else {
        lData = (long *)MemAllocate(laLength * sizeof(Ptr));
      }
      if (!lData) {
        checkPointer(lData);
      }
    } else {
      if (lData) {
        free(lData);
        lData = nil;
      }
    }
  }
}

/*
==============================================================
Sort Methods
==============================================================
*/

void _SimpleList::Sort(bool ascending) {
  if (lLength < 10) { // use bubble sort
    BubbleSort();
  } else {
    QuickSort(0, lLength - 1);
  }

  if (!ascending) {
    long swap, i, j;
    for (i = 0, j = lLength - 1; i < j; i++, j--) {
      swap = ((long *)lData)[i];
      ((long *)lData)[i] = ((long *)lData)[j];
      ((long *)lData)[j] = swap;
    }
  }
}

void _SimpleList::BubbleSort(void) {
  bool done = false;
  long swap, i, j;
  while (!done) {
    done = true;
    for (i = lLength - 1, j = i - 1; i > 0; i--, j--) {
      if (Compare(i, j) < 0) {
        done = false;
        swap = ((long *)lData)[i];
        ((long *)lData)[i] = ((long *)lData)[j];
        ((long *)lData)[j] = swap;
      }
    }
  }
}

void _SimpleList::QuickSort(long from, long to) {
  long middle = (from + to) / 2, middleV = ((long *)lData)[middle], top = to,
       bottommove = 1, topmove = 1, temp, i;

  if (middle)
    //while
    //((middle-bottommove>=from)&&(((long*)lData)[middle-bottommove]>middleV))
    while ((middle - bottommove >= from) &&
           (Compare(middle - bottommove, middle) > 0)) {
      bottommove++;
    }

  if (from < to)
    //while ((middle+topmove<=to)&&(((long*)lData)[middle+topmove]<middleV))
    while ((middle + topmove <= to) &&
           (Compare(middle + topmove, middle) < 0)) {
      topmove++;
    }
  // now shuffle
  for (i = from; i < middle - bottommove; i++) {
    if (Compare(i, middle) > 0) {
      temp = ((long *)lData)[middle - bottommove];
      ((long *)lData)[middle - bottommove] = ((long *)lData)[i];
      ((long *)lData)[i] = temp;
      bottommove++;

      //while
      //((middle-bottommove>=from)&&(((long*)lData)[middle-bottommove]>middleV))
      while ((middle - bottommove >= from) &&
             (Compare(middle - bottommove, middle) > 0)) {
        bottommove++;
      }
    }
  }

  for (i = middle + topmove + 1; i <= top; i++) {
    if (Compare(i, middle) < 0) {
      temp = ((long *)lData)[middle + topmove];
      ((long *)lData)[middle + topmove] = ((long *)lData)[i];
      ((long *)lData)[i] = temp;
      topmove++;

      //while ((middle+topmove<=to)&&(((long*)lData)[middle+topmove]<middleV))
      while ((middle + topmove <= to) &&
             (Compare(middle + topmove, middle) < 0)) {
        topmove++;
      }
    }
  }

  if (topmove == bottommove) {
    for (i = 1; i < bottommove; i++) {
      temp = ((long *)lData)[middle + i];
      ((long *)lData)[middle + i] = ((long *)lData)[middle - i];
      ((long *)lData)[middle - i] = temp;
    }
  } else if (topmove > bottommove) {
    long shift = topmove - bottommove;
    for (i = 1; i < bottommove; i++) {
      temp = ((long *)lData)[middle + i + shift];
      ((long *)lData)[middle + i + shift] = ((long *)lData)[middle - i];
      ((long *)lData)[middle - i] = temp;
    }
    for (i = 0; i < shift; i++) {
      ((long *)lData)[middle + i] = ((long *)lData)[middle + i + 1];
    }
    middle += shift;
    ((long *)lData)[middle] = middleV;
  } else {
    long shift = bottommove - topmove;
    for (i = 1; i < topmove; i++) {
      temp = ((long *)lData)[middle - i - shift];
      ((long *)lData)[middle - i - shift] = ((long *)lData)[middle + i];
      ((long *)lData)[middle + i] = temp;
    }
    for (i = 0; i < shift; i++) {
      ((long *)lData)[middle - i] = ((long *)lData)[middle - i - 1];
    }
    middle -= shift;
    ((long *)lData)[middle] = middleV;
  }
  if (to > middle + 1) {
    QuickSort(middle + 1, top);
  }
  if (from < middle - 1) {
    QuickSort(from, middle - 1);
  }
}

//TODO: This is a global. Should it be here?
void SortLists(_SimpleList *ref, _SimpleList *index) {
  if ((*ref).lLength != index->lLength) {
    return;
  }
  if ((*ref).lLength <= 10) {
    bool done = false;

    while (!done) {
      done = true;
      for (long i = 1; i < (*ref).lLength; i++) {
        if (ref->Compare(i - 1, i) > 0) {
          long swap;
          swap = ((long *)ref->lData)[i];
          ((long *)ref->lData)[i] = ((long *)ref->lData)[i - 1];
          ((long *)ref->lData)[i - 1] = swap;
          swap = ((long *)index->lData)[i];
          ((long *)index->lData)[i] = ((long *)index->lData)[i - 1];
          ((long *)index->lData)[i - 1] = swap;
          done = false;
        }
      }
    }
  } else {
    (*ref).RecursiveIndexSort(0, (*ref).lLength - 1, index);
  }
}



/*
==============================================================
Set Methods
==============================================================
*/

// Compute the intersection of two sorted lists
void _SimpleList::Intersect(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (l1.lData[c1] < l2.lData[c2]) {
      c1++;
      if (c1 == l1.lLength) {
        break;
      }
    }
    if (c1 == l1.lLength) {
      break;
    }

    while (l1.lData[c1] == l2.lData[c2]) {
      (*this) << l1.lData[c1++];
      c2++;
      if (c1 == l1.lLength || c2 == l2.lLength) {
        break;
      }
    }
    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }
    while (l2.lData[c2] < l1.lData[c1]) {
      c2++;
      if (c2 == l2.lLength) {
        break;
      }
    }
  }
}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void _SimpleList::Union(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }

  long c1 = 0, c2 = 0;

  while (c1 < l1.lLength && c2 < l2.lLength) {
    while (l1.lData[c1] < l2.lData[c2]) {
      (*this) << l1.lData[c1++];
      if (c1 == l1.lLength) {
        break;
      }
    }

    if (c1 == l1.lLength) {
      break;
    }

    while (l1.lData[c1] == l2.lData[c2]) {
      (*this) << l1.lData[c1++];
      c2++;
      if (c1 == l1.lLength || c2 == l2.lLength) {
        break;
      }
    }

    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }

    while (l2.lData[c2] < l1.lData[c1]) {
      (*this) << l2.lData[c2++];
      if (c2 == l2.lLength) {
        break;
      }
    }
  }

  while (c1 < l1.lLength) {
    (*this) << l1.lData[c1++];
  }
  while (c2 < l2.lLength) {
    (*this) << l2.lData[c2++];
  }

}

// Compute the union of two sorted lists
// Each repeat appears exactly once
void _SimpleList::XOR(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }

  long c1 = 0, c2 = 0;

  while ((c1 < l1.lLength) && (c2 < l2.lLength)) {
    while (c1 < l1.lLength && l1.lData[c1] < l2.lData[c2]) {
      (*this) << l1.lData[c1++];
    }
    if (c1 == l1.lLength) {
      break;
    }
    while (c1 < l1.lLength && c2 < l2.lLength && l1.lData[c1] == l2.lData[c2]) {
      c1++;
      c2++;
    }
    if (c1 == l1.lLength || c2 == l2.lLength) {
      break;
    }

    while (c2 < l2.lLength && l2.lData[c2] < l1.lData[c1]) {
      (*this) << l2.lData[c2++];
    }
  }

  while (c1 < l1.lLength) {
    (*this) << l1.lData[c1++];
  }
  while (c2 < l2.lLength) {
    (*this) << l2.lData[c2++];
  }
}
