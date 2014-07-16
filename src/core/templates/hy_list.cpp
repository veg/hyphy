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
#include <stdlib.h>
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
_hyList<PAYLOAD>::_hyList(const PAYLOAD item) {
  lLength     = 1UL;
  laLength    = HY_LIST_ALLOCATION_CHUNK;
  lData       = (PAYLOAD *)MemAllocate(laLength * sizeof(PAYLOAD));
  lData[0]    = item;
}


//Stack copy contructor
template<typename PAYLOAD>
_hyList<PAYLOAD>::_hyList(const _hyList <PAYLOAD> &l, const long from, const long to) {
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
  _hyList<PAYLOAD> res;
  unsigned long combined_length = l.lLength + lLength;
  if (combined_length == 0UL) {
    return res;
  }
  res.RequestSpace (combined_length);

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
void _hyList<PAYLOAD>::append_multiple(const PAYLOAD item, const unsigned long copies)
{
  RequestSpace (laLength + copies);
  for (unsigned long i = 0UL; i < copies; i++) {
    append (item);
  }
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::RequestSpace(const unsigned long slots, bool set_length)
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
  if (set_length) {
    lLength = slots;
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
        lData = (PAYLOAD *)MemReallocate((Ptr)lData, laLength * sizeof(PAYLOAD));
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
      for (unsigned long el = index; el < lLength; el ++)  {
        lData[el] = lData[el+1];
      }
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
      if (k < del_list_length && i == (unsigned long)indices_to_delete->AtIndex(k)) {
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

  if (lLength > 0UL) {

    unsigned long ustart, 
                  uend;

    if (start < 0L) {
      ustart = 0L;
    } else if ((unsigned long)start >= lLength) {
      ustart = lLength - 1L;
    } else {
      ustart = start;
    }

    if (end < 0L || (unsigned long)end >= lLength) {
      uend = lLength - 1L;
    } else {
      uend = end;
    }
    
    
    if (uend >= ustart && delta != 0L && uend - ustart + 1UL < lLength) {
      if (delta > 0L && lLength <= delta + uend) { // shift up
        delta = lLength - uend - 1L;
      } else if (start + delta < 0L) {
        delta = -ustart;
      }
      
      if (delta) {
        unsigned long span = uend - ustart + 1L;
        
        
        if (delta > 0L) {
          if (span > delta) {
            Displace (uend+1UL, uend + delta, -span);
            return;
          }
        } else {
          if (span > -delta) {
             Displace (ustart + delta, ustart - 1UL, span);
             return;
          }
        }
        
        
        _hyList <PAYLOAD> swapList;
        swapList.RequestSpace (span);

        for (unsigned long i = ustart; i <= uend; i++) {
          swapList.append (lData[i]);
        }

        if (delta > 0L) {
          for (unsigned long i = ustart; i <= ustart + delta; i++) {
            lData[i] = lData[i+span];
          }

        } else {
          for (long i = (long)ustart - 1L; i >= (long)ustart + delta; i--) {
            lData[i+span] = lData[i];
          }
        }
  
        for (unsigned long i = ustart + delta, j = 0UL; i <= uend + delta; i++, j++) {
          lData[i] = swapList.lData[j];
        }
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
PAYLOAD _hyList<PAYLOAD>::Element(const long index) const
{
  if (index >= 0L && (unsigned long)index < lLength) {
    return lData[index];
  } else if ((unsigned long)(-index) <= lLength) {
    return lData[(long)lLength + index];
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
    if (!ItemEqualToValue (i, l2 (i))) {
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
      PAYLOAD pt;
      SWAP (lData[k], lData[l], pt);
    }
  }
}

template<typename PAYLOAD>
void _hyList<PAYLOAD>::InsertElement(const PAYLOAD item, const long insert_at)
{

  lLength++;
  ResizeList();
  
 
  if (insert_at == HY_LIST_INSERT_AT_END || insert_at < 0L) {
    lData[lLength-1UL] = item;
  } else {
    long insert_here = insert_at >= lLength ? lLength - 1UL : insert_at;
    long moveThisMany = (laLength - insert_here - 1L);
    if (moveThisMany < 32L)
      for (long k = insert_here + moveThisMany; k > insert_here; k--) {
        lData[k] = lData[k - 1];
      }
    else {
      memmove(&lData[insert_here + 1], &lData[insert_here],
              moveThisMany * sizeof(PAYLOAD));
    }
    lData[insert_here] = item;
  }
}


template<typename PAYLOAD>
_hyList<PAYLOAD> *_hyList<PAYLOAD>::Subset(const unsigned long size, const bool replacement) const
{
  _hyList<PAYLOAD> *result = new _hyList<PAYLOAD>;
  if (size > 0UL) {
    unsigned long sample_size = MIN(size, lLength);
    if (replacement) {
      for (unsigned long k = 0UL; k < sample_size; k++) {
        result->append (lData[genrand_int32() % lLength]);
      }
    } else {
      (*result) << (*this);
      for (unsigned long k = 0UL; k < sample_size; k++) {
        long idx = genrand_int32() % (lLength - k);
        PAYLOAD t;
        SWAP (result->lData[k], result->lData[k+idx], t);
      }
      result->lLength = sample_size;
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
    for (unsigned long k = 0UL; k < blockCount - 1; k ++) {
      unsigned long k2 = genrand_int32() % (blockCount - k);
      if (k2) {
        k2 += k;
        k2 *= blockLength;

        for (unsigned long j = 0UL; j < blockLength; j++) {
          PAYLOAD t;
          SWAP (lData[k2 + j],lData[k * blockLength + j],t);
        }
      }
    }
  } else {
    for (unsigned long k = 0; k < blockCount - 1; k ++) {
      unsigned long k2 = genrand_int32() % (blockCount - k);
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
  _hyList<PAYLOAD> result;
  result.RequestSpace (blockCount * blockLength);
  if (blockLength > 1UL)
    for (unsigned long i = 0UL; i < blockCount; i++) {
      unsigned long sample = genrand_int32() % blockCount;
      sample *= blockLength;
      for (unsigned long j = 0UL; j < blockLength; j++) {
        result << lData[sample + j];
      }
    }
  else {
    for (unsigned long i = 0UL; i < blockCount; i++) {
      unsigned long sample = genrand_int32() % blockCount;
      result << lData[sample];
    }
  }
  Clone (&result);
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
BaseRef _hyList<PAYLOAD>::toStr(void) const
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
    state->append_multiple (HY_NOT_FOUND, stride+3);
    state->SetItem (0L, stride);
    store.Clear();
    store.append_multiple (lData[0], stride);
    return true;
  }
  return false;
}

// Implements algorithm NEXKSB from p.27 of
// http://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf
template<typename PAYLOAD>
bool _hyList<PAYLOAD>::NChooseK(_hyList <long> *state, _hyList <PAYLOAD>& store) {
  if (state->AtIndex (1) == HY_NOT_FOUND) {  // first pass
    state->SetItem (1L ,0L);                 // m
    state->SetItem (2L, state->AtIndex (0)); // h
    
    if (store.countitems() == 0UL) {
      return false;
    }
  } else {
    if (state->AtIndex(1) < lLength - state->AtIndex(2)) {
      state->SetItem(2UL,0L);
    }
    (*state)[2]++;
    (*state)[1] = state->AtIndex(3 + state->AtIndex(0) - state->AtIndex(2)) + 1L;
  }
  
  
  for (long j = 1L; j <= state->AtIndex(2); j++) {
    long anIndex = j + state->AtIndex (0) - state->AtIndex(2),
         anIndex2 = state->AtIndex(1) + j -1 ;
    (*state)[anIndex + 2] = anIndex2;
    store[anIndex - 1] = lData[anIndex2];
  }
  
  //printf ("[%ld %ld %ld %ld]\n", state->AtIndex(0),  state->AtIndex(1),  state->AtIndex(2),  state->AtIndex(3));
  
  return state->AtIndex (3) < lLength - state->AtIndex (0);
}

