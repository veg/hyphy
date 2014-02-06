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

#include <stdarg.h>
#include <string.h>
#include "helperfunctions.h"
#include "errorfns.h"
#include "hy_string_buffer.h"



/*
==============================================================
Constructors
==============================================================
*/

// Does nothing
template<typename PAYLOAD>
_hyList<PAYLOAD>::_hyList()
{
  Initialize(false);
}

//Data constructor (1 member list)
template<typename PAYLOAD>
_hyList<PAYLOAD>::_hyList(const PAYLOAD item)
{
  lLength     = 1UL;
  laLength    = HY_LIST_ALLOCATION_CHUNK;
  lData       = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
  lData[0]    = item;
}


//Stack copy contructor
template<typename PAYLOAD>
_hyList<PAYLOAD>::_hyList(const _hyList <PAYLOAD> &l, const long from, const long to)
{
  Initialize ();
  Clone (&l, from, to);
}

// Data constructor (variable number of long constants)
template<typename PAYLOAD>
_hyList<PAYLOAD>::_hyList(const unsigned long number, const PAYLOAD items[])
{
  Initialize(true);
  for (unsigned long arg_id = 0UL; arg_id < number; arg_id++) {
    append(items[arg_id]);
  }
}

//Destructor
template<typename PAYLOAD>
_hyList<PAYLOAD>::~_hyList(void)
{
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
PAYLOAD &_hyList<PAYLOAD>::operator[](const long i)
{
  if (lLength == 0) {
    warnError ("_hyList[] called on an empty list");
  }

  const unsigned long in = (const unsigned long) i;
  if (in > lLength - 1) {
    return lData[lLength - 1];
  }

  return lData[in];
}
//Element location functions (0,llength - 1)
template<typename PAYLOAD>
PAYLOAD _hyList<PAYLOAD>::operator()(const unsigned long i) const
{
  if (i < lLength) {
    return lData[i];
  }
  warnError("List index out of range");
  return lData[0];
}

//Assignment operator
template<typename PAYLOAD>
_hyList<PAYLOAD> const _hyList<PAYLOAD>::operator=(const _hyList<PAYLOAD>& l)
{
  Clear();
  Clone (&l);
  return *this;
}

//Append operator
template<typename PAYLOAD>
const _hyList<PAYLOAD> _hyList<PAYLOAD>::operator&(const _hyList<PAYLOAD> l)
{
  _hyList<PAYLOAD> res(l.lLength + lLength);
  if (res.laLength == 0UL) {
    return res;
  }

  if (lData && lLength) {
    memcpy((Ptr)res.lData, (Ptr)lData, lLength * sizeof(PAYLOAD));
  }

  if (l.lData && l.lLength) {
    memcpy((Ptr)&(res.lData [lLength]), (Ptr)l.lData,
           l.lLength * sizeof(PAYLOAD));
  }

  res.lLength = l.lLength + lLength;

  return res;
}



template<typename PAYLOAD>
void _hyList<PAYLOAD>::operator<<(const PAYLOAD item)
{
  append (item);
}

template<typename PAYLOAD>
bool _hyList<PAYLOAD>::operator>>(const PAYLOAD item)
{
  if (Find(item) == HY_NOT_FOUND) {
    InsertElement(item, HY_LIST_INSERT_AT_END);
    return true;
  }
  return false;
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::operator<<(const _hyList<PAYLOAD> &source)
{
  for (unsigned long k = 0; k < source.lLength; k++) {
    append(source.AtIndex(k));
  }
}


template<typename PAYLOAD>
bool _hyList<PAYLOAD>::operator == (const _hyList<PAYLOAD> &source) const
{
  return this->Equal (source);
}

/*
==============================================================
Methods
==============================================================
*/


template<typename PAYLOAD>
void _hyList<PAYLOAD>::Clone(const _hyList<PAYLOAD>* clone_from, const long from, const long to) {
  if (from == 0UL && to == HY_LIST_INSERT_AT_END) {
    lLength  = clone_from->lLength;
    RequestSpace (clone_from->laLength);
    if (lLength) {
      memcpy((Ptr)lData, (Ptr)clone_from->lData, lLength * sizeof(PAYLOAD));
    }
  } else {
    long f = from, t = to;
    NormalizeCoordinates(f, t, clone_from->lLength);
    long upto = t - f + 1;
    RequestSpace(upto);
    for (long k = 0L; k < upto; k++) {
      lData[k] = clone_from->lData[from + k];
    }
    lLength = upto;
  }
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::append(const PAYLOAD item)
{
  InsertElement(item, HY_LIST_INSERT_AT_END);
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::RequestSpace(const unsigned long slots)
{
  if (slots > laLength) {
    laLength = (slots / HY_LIST_ALLOCATION_CHUNK + 1) * HY_LIST_ALLOCATION_CHUNK;
    if (lData) {
      checkPointer(lData = (PAYLOAD *)MemReallocate((Ptr)lData,
                           laLength * sizeof(PAYLOAD)));
    } else {
      checkPointer(lData = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD)));
    }
  }
}

template<typename PAYLOAD>
PAYLOAD _hyList<PAYLOAD>::AtIndex(const unsigned long index) const
{
  return lData[index];
}

template<typename PAYLOAD>
unsigned long _hyList<PAYLOAD>::countitems(void) const
{
  return lLength;
}

template<typename PAYLOAD>
unsigned long _hyList<PAYLOAD>::allocated(void) const
{
  return laLength;
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::Clear(bool completeClear)
{
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
void _hyList<PAYLOAD>::CompactList(void)
{
  if (laLength - lLength > HY_LIST_ALLOCATION_CHUNK) {
    laLength -= ((laLength - lLength) / HY_LIST_ALLOCATION_CHUNK) * HY_LIST_ALLOCATION_CHUNK;
    if (laLength) {
      lData = (PAYLOAD *)MemReallocate((Ptr)lData, laLength * sizeof(PAYLOAD));
    } else {
      free(lData);
      lData = nil;
    }
  }

}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::TrimMemory(void)
{
  if (laLength > lLength) {
    laLength = lLength;
    if (laLength) {
      if (lData) {
        lData = (PAYLOAD *)MemReallocate(lData, laLength * sizeof(PAYLOAD));
      } else {
        lData = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
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

template<typename PAYLOAD>
void _hyList<PAYLOAD>::Initialize(bool doMemAlloc)
{
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
BaseRef _hyList<PAYLOAD>::makeDynamic(void) const
{
  _hyList <PAYLOAD> *res = new _hyList <PAYLOAD>;
  res->Duplicate(this);
  return res;
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::ResizeList(void)
{
  if (lLength > laLength) {
    unsigned long incBy = (HY_LIST_ALLOCATION_CHUNK * 5UL > lLength) ? HY_LIST_ALLOCATION_CHUNK : lLength / 5UL;

    laLength += incBy;

    if (lData) {
      lData = (PAYLOAD *)MemReallocate((Ptr)lData, laLength * sizeof(PAYLOAD));
    } else {
      lData = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
    }

    if (!lData) {
      checkPointer(lData);
    }
  }
}

//Delete item at index (>=0)
template<typename PAYLOAD>
void _hyList<PAYLOAD>::Delete(const long index, bool compact_list)
{
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


template<typename PAYLOAD>
void _hyList<PAYLOAD>::DeleteList(const _hyList<long> *indices_to_delete)
{
  
  unsigned long del_list_length = indices_to_delete->countitems();
  
  if (del_list_length) {
    unsigned long k = 0UL;
    
    for (unsigned long i = 0UL; i < lLength; i++) {
      if (k < del_list_length && i == indices_to_delete->AtIndex(k)) {
        k++;
      } else {
        lData[i - k] = lData[i];
      }
    }
    lLength -= del_list_length;
  }
  CompactList();
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::Displace(long start, long end, long delta)
{
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
      _hyList <PAYLOAD> swapList((unsigned long)(end - start + 1));

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
void _hyList<PAYLOAD>::Duplicate(BaseRefConst theRef)
{
  const _hyList<PAYLOAD> *l = dynamic_cast<const _hyList<PAYLOAD>* >(theRef);
  lLength   = l->lLength;
  laLength  = l->laLength;
  lData     = l->lData;
  if (lData) {
    lData = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
    memcpy((Ptr)lData, (Ptr)l->lData, lLength * sizeof(PAYLOAD));
  }
}

//Element location functions (0,llength - 1)
//Negative indices return offsets from the end of the list
template<typename PAYLOAD>
PAYLOAD _hyList<PAYLOAD>::Element(const long index)
{
  if (index >= 0L && index < lLength) {
    return lData[index];
  } else if (-index <= lLength) {
    return lData[lLength - (-index)];
  }
  return lData[0];
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::SetItem(const unsigned long index, PAYLOAD item)
{
  lData[index] = item;
}

template<typename PAYLOAD>
bool _hyList<PAYLOAD>::Equal(const _hyList<PAYLOAD> &l2) const
{
  if (lLength != l2.lLength) {
    return false;
  }

  for (unsigned long i = 0UL; i < lLength; i++)
    if (ItemEqualToValue (i, l2 (i)) {
      return false;
    }

  return true;
}

template<typename PAYLOAD>
bool _hyList<PAYLOAD>::ItemEqualToValue(unsigned long index, const PAYLOAD& value) const
{
  return lData[index] == value;
}


template<typename PAYLOAD>
long _hyList<PAYLOAD>::Find(const PAYLOAD item, const long startAt) const
{
  for (unsigned long i = startAt; i < lLength; i++) {
    if (ItemEqualToValue (i, item)) {
      return i;
    }
  }
  return HY_NOT_FOUND;
}

template<typename PAYLOAD>
long _hyList<PAYLOAD>::FindStepping(const PAYLOAD item, const long step,
                                    const long startAt) const
{
  for (unsigned long i = startAt; i < lLength; i += step)
    if (ItemEqualToValue (i, item)) {
      return i;
    }

  return HY_NOT_FOUND;
}


template<typename PAYLOAD>
void _hyList<PAYLOAD>::Flip()
{
  if (lLength > 0UL) {
    for (unsigned long k = 0UL, l = lLength - 1UL; k < l; k++, l--) {
      PAYLOAD pt = lData[k];
      lData[k] = lData[l];
      lData[l] = pt;
    }
  }
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::InsertElement(const PAYLOAD item, const long insert_at)
{

  lLength++;
  ResizeList();

  if (insert_at == HY_LIST_INSERT_AT_END) {
    lData[lLength-1UL] = item;
  } else {
    long insert_here = insert_at >= lLength ? lLength - 1UL : insert_at;
    long moveThisMany = (laLength - insert_here - 1L);
    if (moveThisMany < 32L)
      for (long k = insert_here + moveThisMany; k > insert_here; k--) {
        lData[k] = lData[k - 1];
      }
    else {
      memmove(lData + (insert_here + 1), lData + insert_here,
              moveThisMany * sizeof(PAYLOAD));
    }
  }
}


template<typename PAYLOAD>
_hyList<PAYLOAD> *_hyList<PAYLOAD>::Subset(unsigned long size, bool replacement)
{
  _hyList<PAYLOAD> *result = new _hyList<PAYLOAD>;
  if (size > 0UL) {
    size = MIN(size, lLength);
    if (replacement) {
      for (long k = 0; k < size; k++) {
        append (lData[genrand_int32() % lLength]);
      }
    } else {
      (*result) << (*this);
      for (unsigned long k = 0UL; k < size; k++) {
        long idx = lData[genrand_int32() % (lLength - k)];
        PAYLOAD t;
        SWAP (result->lData[k], result->lData[idx], t);
      }
      result->lLength = size;
      result->TrimMemory();
    }
  }
  return result;
}

// Create a permutation of the list's elements
template<typename PAYLOAD>
void _hyList<PAYLOAD>::Permute(const unsigned long blockLength)
{
  unsigned long blockCount = lLength / blockLength;

  if (blockLength > 1UL) {
    for (unsigned long k = 0; k < blockCount - 1; k ++) {
      unsigned long k2 = genrand_real2() * (blockCount - k);
      if (k2) {
        k2 += k;
        k2 *= blockLength;

        for (long j = 0; j < blockLength; j++) {
          PAYLOAD t;
          SWAP (lData[k2 + j],lData[k * blockLength + j],t);
        }
      }
    }
  } else {
    for (unsigned long k = 0; k < blockCount - 1; k = k + 1) {
      unsigned long k2 = genrand_real2() * (blockCount - k);
      if (k2) {
        k2 += k;
        PAYLOAD t;
        SWAP (lData[k], lData[k2], t);
      }
    }
  }
}

// Create a permutation of the list's elements with possible repetitions
template<typename PAYLOAD>
void _hyList<PAYLOAD>::PermuteWithReplacement(const unsigned long blockLength)
{
  unsigned long blockCount = lLength / blockLength;
  _hyList<PAYLOAD> result((unsigned long)(blockCount * blockLength));
  if (blockLength > 1UL)
    for (unsigned long i = 0; i < blockCount; i++) {
      unsigned long sample = (unsigned long)(genrand_real2() * blockCount);
      sample *= blockLength;
      for (unsigned long j = 0; j < blockLength; j++, sample++) {
        result << lData[sample];
      }
    }
  else {
    for (unsigned long i = 0; i < blockCount; i++) {
      unsigned long sample = genrand_real2() * blockCount;
      result << lData[sample];
    }
  }

  Clear();
  Duplicate(&result);
}

template<typename PAYLOAD>
PAYLOAD _hyList<PAYLOAD>::Pop(void)
{
  if (lLength > 0UL) {
    lLength--;
    return lData[lLength];
  }

  warnError ("_hyList::Pop called on an empty list");
  return lData[0UL];
}



template<typename PAYLOAD>
void _hyList<PAYLOAD>::Swap(const unsigned long i, const unsigned long j)
{
  if (i < lLength && j < lLength) {
      PAYLOAD t;
      SWAP (lData[i], lData[j], t);
  }
}

template<typename PAYLOAD>
BaseRef _hyList<PAYLOAD>::toStr(void)
{
  _StringBuffer * stringified = new _StringBuffer ();
  (*stringified) << "_hyList with ";
  (*stringified) << (long)lLength;
  (*stringified) << " elements.";
  
  return stringified;
}

// Together with the next function
// Implements algorithm NEXKSB from p.27 of
// http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
template<typename PAYLOAD>
bool _hyList<PAYLOAD>::NChooseKInit(_hyList <long> *state, _hyList <PAYLOAD>& store,
                               const unsigned long stride, bool algorithm) {
  if (stride <= lLength && lLength) {
    state->Clear();
    state->RequestSpace(stride + 3);
    state->append(stride);
    store.Clear();
    store.RequestSpace(stride);
    return true;
  }
  return false;
}

// Implements algorithm NEXKSB from p.27 of
// http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
template<typename PAYLOAD>
bool _hyList<PAYLOAD>::NChooseK(_hyList <long> *state, _hyList <PAYLOAD>& store) {
  if (state->lLength == 1) {  // first pass
    state->append (0L);              // m
    state->append (state->lData[0]); // h
    state->lLength = state->lData[0] + 3;
    store.lLength = state->lData[0];
    
    if (store.lLength == 0) {
      return false;
    }
  } else {
    if (state->lData[1] < lLength - state->lData[2]) {
      state->lData[2] = 0L;
    }
    state->lData[2]++;
    state->lData[1] = state->lData[3 + state->lData[0] - state->lData[2]] + 1L;
  }
  for (long j = 1L; j <= state->lData[2]; j++) {
    long anIndex = j + state->lData[0] - state->lData[2],
         anIndex2 = state->lData[1] + j;
    state->lData[anIndex + 2] = anIndex2 - 1;
    store.lData[anIndex - 1] = lData[anIndex2 - 1];
  }
  return state->lData[3] < lLength - state->lData[0];
}

