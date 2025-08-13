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

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "function_templates.h"
#include "global_things.h"
#include "hy_string_buffer.h"
#include "hy_strings.h"
#include "list.h"
#include "mersenne_twister.h"
#include "parser.h"
#include "simplelist.h"

using namespace hy_global;

/*
==============================================================
Constructors
==============================================================
*/

/**
 * @brief Default constructor. Initializes an empty list.
 */
_SimpleList::_SimpleList() : static_data{0L} { Initialize(false); }

/**
 * @brief Constructor for a list with a single element.
 * @param br The initial element to add to the list.
 */
_SimpleList::_SimpleList(long br) : static_data{0L} {
  lLength = 1;
  laLength = MEMORYSTEP;
  list_data = (long *)&static_data;
  static_data[0] = br;
}

/**
 * @brief Constructor that creates and populates a list with a sequence of
 * numbers.
 * @param l The number of elements in the list.
 * @param start The starting value of the sequence.
 * @param step The increment between consecutive elements.
 * @code
 *  _SimpleList list(5, 0, 2); // Creates a list with elements {0, 2, 4, 6, 8}
 * @endcode
 */
_SimpleList::_SimpleList(long l, long start, long step) : static_data{0L} {
  Initialize(false);
  Populate(l, start, step);
}

/**
 * @brief Constructor that pre-allocates memory for a given number of elements.
 * @param l The number of elements to allocate space for.
 */
_SimpleList::_SimpleList(unsigned long l) : static_data{0L} {
  lLength = 0UL;
  if (l <= MEMORYSTEP) {
    laLength = MEMORYSTEP;
    list_data = (long *)&static_data;
  } else {
    laLength = (l / MEMORYSTEP + 1) * MEMORYSTEP;
    list_data = static_cast<long *>(MemAllocate(laLength * sizeof(hyPointer)));
    memset(list_data, 0, laLength * sizeof(hyPointer));
  }
}

/**
 * @brief Constructor that uses a pre-allocated buffer for storage.
 * @param l The number of elements the pre-allocated storage can hold.
 * @param preallocated_storage A pointer to the pre-allocated buffer. If NULL,
 * new memory is allocated.
 */
_SimpleList::_SimpleList(unsigned long l, long *preallocated_storage)
    : static_data{0L} {
  lLength = 0UL;
  if (preallocated_storage) {
    laLength = l;
    list_data = preallocated_storage;
    AddAReference();
  } else {
    if (l <= MEMORYSTEP) {
      laLength = MEMORYSTEP;
      list_data = (long *)&static_data;
    } else {
      laLength = (l / MEMORYSTEP + 1) * MEMORYSTEP;
      list_data =
          static_cast<long *>(MemAllocate(laLength * sizeof(hyPointer)));
    }
    memset(list_data, 0, laLength * sizeof(hyPointer));
  }
}

/**
 * @brief Copy constructor that copies a sub-list of another list.
 * @param l The source list to copy from.
 * @param from The starting index of the sub-list to copy (inclusive).
 * @param to The ending index of the sub-list to copy (exclusive).
 */
_SimpleList::_SimpleList(_SimpleList const &l, long from, long to)
    : static_data{0L} {
  Initialize(false);
  if (from == 0L && to == -1L) {
    to = l.lLength;
  } else {
    NormalizeCoordinates(from, to, l.lLength);
  }
  if (to > from) {
    unsigned long upto = to - from;
    RequestSpace(upto);
    memcpy(list_data, l.list_data + from, upto * sizeof(long));
    lLength = upto;
  }
}

/**
 * @brief Constructor that initializes the list with a variable number of
 * arguments.
 * @param value1 The first value to add.
 * @param number The number of additional arguments.
 * @param ... The additional long integer values.
 * @code
 * _SimpleList list(10, 2, 20, 30); // Creates a list with elements {10, 20, 30}
 * @endcode
 */
_SimpleList::_SimpleList(const long value1, const unsigned long number, ...)
    : static_data{0L} {
  Initialize(true);
  va_list vl;

  (*this) << value1;

  va_start(vl, number);
  for (unsigned long arg_id = 0; arg_id < number; arg_id++) {
    const long this_arg = va_arg(vl, long);
    (*this) << this_arg;
  }
  va_end(vl);
}

/**
 * @brief Destructor. Frees the memory allocated by the list.
 */
_SimpleList::~_SimpleList(void) {
  if (CanFreeMe()) {
    if (list_data && list_data != static_cast<long *>(static_data)) {
      free(list_data);
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

/**
 * @brief Access element by index (l-value).
 * @param i The index of the element to access.
 * @return A reference to the element at the specified index.
 * @note No bounds checking is performed.
 */
long &_SimpleList::operator[](const long i) { return list_data[i]; }

/**
 * @brief Access element by index (r-value).
 * @param i The index of the element to access.
 * @return The value of the element at the specified index.
 * @note Performs bounds checking and calls HandleApplicationError on failure.
 */
long _SimpleList::operator()(const unsigned long i) const {

  if (i < lLength) {
    return list_data[i];
  }
  HandleApplicationError("List index out of range");
  return -1;
}

/**
 * @brief Assignment operator.
 * @param l The list to assign from.
 * @return A reference to this list.
 */
const _SimpleList &_SimpleList::operator=(_SimpleList const &l) {
  if (l.laLength && list_data && laLength >= l.lLength &&
      laLength <= l.laLength && CanFreeMe()) {
    // reuse existing memory allocation
    lLength = l.lLength;
    if (lLength) {
      memcpy(list_data, l.list_data, lLength * sizeof(hyPointer));
    }
  } else { // allocate new memory
    Clear();
    lLength = l.lLength;
    laLength = l.laLength;
    if (laLength >= MEMORYSTEP) {
      list_data =
          static_cast<long *>(MemAllocate(laLength * sizeof(hyPointer)));
    } else {
      list_data = static_data;
    }
    if (lLength) {
      memcpy(list_data, l.list_data, lLength * sizeof(hyPointer));
    }
  }

  return *this;
}

/**
 * @brief Concatenation operator.
 * @param l The list to concatenate with.
 * @return A new list containing the elements of both lists.
 */
_SimpleList const _SimpleList::operator&(_SimpleList const &l) {
  _SimpleList res(l.lLength + lLength);

  if (!res.laLength) {
    return res;
  }

  if (list_data && lLength) {
    memcpy(res.list_data, list_data, lLength * sizeof(void *));
  }

  if (l.list_data && l.lLength) {
    memcpy((char *)res.list_data + lLength * sizeof(void *), l.list_data,
           l.lLength * sizeof(void *));
  }

  res.lLength = l.lLength + lLength;

  return res;
}

/**
 * @brief Append element operator.
 * @param br The element to append.
 * @return A reference to this list.
 */
_SimpleList &_SimpleList::operator<<(long br) {
  _SimpleList::InsertElement((BaseRef)br, -1, false, false);
  return *this;
}

/**
 * @brief Append element if not already present.
 * @param br The element to append.
 * @return `true` if the element was appended, `false` otherwise.
 */
bool _SimpleList::operator>>(long br) {
  if (Find(br) == -1) {
    InsertElement((BaseRef)br, -1, false, false);
    return true;
  }
  return false;
}

/**
 * @brief Append another list.
 * @param source The list to append.
 */
void _SimpleList::operator<<(_SimpleList const &source) {
  for (unsigned long k = 0UL; k < source.lLength; k++) {
    (*this) << source.list_data[k];
  }
}

/*
==============================================================
Methods
==============================================================
*/

/**
 * @brief Get an element at a specific index.
 * @param index The index of the element to retrieve. Negative indices count
 * from the end of the list.
 * @return The value of the element at the specified index.
 * @note Performs bounds checking and calls HandleApplicationError on failure.
 */
long _SimpleList::GetElement(const long index) const {
  if (index >= 0L) {
    if (static_cast<const unsigned long>(index) < lLength) {
      return list_data[index];
    }
  }
  if (static_cast<const unsigned long>(-index) <= lLength) {
    return list_data[lLength + index];
  }
  HandleApplicationError(
      _String("List index '") &
      static_cast<long>(static_cast<const unsigned long>(-index)) &
      "' out of range in _SimpleList::GetElement on list of length " &
      static_cast<long>(lLength));
  return 0;
}

/**
 * @brief Performs a binary search for a value in a sorted list.
 * @param s The value to search for.
 * @param startAt The index to start the search from.
 * @return If the value is found, its index is returned. Otherwise, a negative
 * value is returned, where `-index - 2` is the insertion point.
 */
long _SimpleList::BinaryFind(long s, long startAt) const {

  if (lLength == 0L) {
    return -2;
  }

  long top = lLength - 1, bottom = startAt, middle;

  while (top >= bottom) {
    middle = bottom + ((top - bottom) >> 1);
    long middle_value = list_data[middle];
    if (s < middle_value) {
      top = middle - 1;
    } else if (s > middle_value) {
      bottom = middle + 1;
    } else {
      return middle;
    }
  }

  return -bottom - 2;

  // exiting here; because no match was found and middle = top = bottom
  // if middle_value < s, then bottom = middle + 1, so we return -middle - 3
  // if middle_value < s then bottom = middle - 1, otherwise bottom = middle

  // comp < 0, means middle_value < s
  // return comp < 0 ? -middle - 3 : -middle - 2;
}

/**
 * @brief Inserts a value into a sorted list, maintaining the sort order.
 * @param n The value to insert.
 * @return The index where the value was inserted.
 */
long _SimpleList::BinaryInsert(long n) {
  if (lLength == 0L) {
    (*this) << n;
    return 0L;
  }

  long pos = -BinaryFind(n) - 2;

  if (pos < 0L) {
    return -pos + 2;
  }

  if (pos < (long)lLength && list_data[pos] < n) {
    pos++;
  }

  InsertElement((BaseRef)n, pos, false, false);

  return pos >= (long)lLength ? lLength - 1 : pos;
}

/**
 * @brief Deletes all _Formula objects stored in the list.
 * @note This function assumes that the list contains pointers to _Formula
 * objects.
 */
void _SimpleList::ClearFormulasInList(void) {
  for (unsigned long k = 0UL; k < lLength; k++)
    if (list_data[k]) {
      delete (_Formula *)list_data[k];
    }
}

/**
 * @brief Calculates the sum of all elements in the list.
 * @return The sum of all elements.
 */
long _SimpleList::Sum(void) const {
  long sum = 0L;
  for (unsigned long k = 0UL; k < lLength; k++) {
    sum += list_data[k];
  }
  return sum;
}

/**
 * @brief Compares two elements in the list.
 * @param i The index of the first element.
 * @param j The index of the second element.
 * @return kCompareLess if the first element is smaller, kCompareEqual if they
 * are equal, and kCompareGreater if the first element is larger.
 */
hyComparisonType _SimpleList::Compare(long i, long j) const {
  long v1 = list_data[i], v2 = list_data[j];

  if (v1 < v2) {
    return kCompareLess;
  } else if (v1 == v2) {
    return kCompareEqual;
  } else {
    return kCompareGreater;
  }

  // return ((long*)list_data)[i]-((long*)list_data)[j];
}

/**
 * @brief Compares a BaseObj pointer with an element in the list.
 * @param i A pointer to a BaseObj.
 * @param j The index of the element in the list to compare with.
 * @return kCompareLess if the pointer value is smaller, kCompareEqual if they
 * are equal, and kCompareGreater if the pointer value is larger.
 */
hyComparisonType _SimpleList::Compare(BaseObj const *i, long j) const {
  long v1 = reinterpret_cast<long>(i), v2 = list_data[j];

  if (v1 < v2) {
    return kCompareLess;
  } else if (v1 == v2) {
    return kCompareEqual;
  } else {
    return kCompareGreater;
  }

  // return (long)i-((long*)list_data)[j];
}

/**
 * @brief Counts the number of common elements between two sorted lists.
 * @param l1 The other sorted list to compare with.
 * @param yesNo If true, the function returns 1 as soon as a common element is
 * found. Otherwise, it counts all common elements.
 * @return The number of common elements.
 */
long _SimpleList::CountCommonElements(_SimpleList const &l1, bool yesNo) const {
  unsigned long res = 0;

  unsigned long i = 0, j = 0;
  while (i < l1.lLength && j < lLength) {
    if (l1.list_data[i] < list_data[j]) {
      i++;
    } else if (list_data[j] < l1.list_data[i]) {
      j++;
    } else {
      if (yesNo)
        return 1;
      res++;
      i++;
      j++;
    }
  }

  return res;
}

/**
 * @brief Performs a counting sort on the list.
 * @param upperBound The upper bound of the values in the list. If 0 or less, it
 * is computed from the list.
 * @param ordering A list to store the original indices of the sorted elements.
 * @param wantResult If true, a new sorted list is returned.
 * @return If `wantResult` is true, a new sorted list. Otherwise, `nil`.
 */
_SimpleList *_SimpleList::CountingSort(long upperBound, _SimpleList *ordering,
                                       bool wantResult) {
  static const long kAllocaLimit = 4096L;

  if (ordering) {
    ordering->Clear();
  }

  if (lLength) {
    if (upperBound <= 0) {
      upperBound = Max() + 1;
    }

    long *storage = upperBound < kAllocaLimit
                        ? static_cast<long *>(alloca(sizeof(long) * upperBound))
                        : nil;

    if (storage) {
      memset(storage, 0, sizeof(long) * upperBound);
    }

    _SimpleList buffer(upperBound, storage);
    _SimpleList *result =
        (wantResult || !ordering) ? new _SimpleList(lLength) : nil;

    for (unsigned long pass1 = 0; pass1 < lLength; pass1++) {
      buffer.list_data[list_data[pass1]]++;
    }
    for (long pass2 = 1; pass2 < upperBound; pass2++) {
      buffer.list_data[pass2] += buffer.list_data[pass2 - 1];
    }

    if (!wantResult && ordering) {
      ordering->Populate(lLength, 0, 0);
      for (long pass3 = lLength - 1; pass3 >= 0; pass3--) {
        ordering->list_data[--buffer.list_data[list_data[pass3]]] = pass3;
      }
    } else {
      if (ordering) {
        ordering->Populate(lLength, 0, 0);
        for (long pass3 = lLength - 1; pass3 >= 0; pass3--) {
          result->list_data[--buffer.list_data[list_data[pass3]]] =
              list_data[pass3];
          ordering->list_data[buffer.list_data[list_data[pass3]]] = pass3;
        }
      } else
        for (long pass3 = lLength - 1; pass3 >= 0; pass3--) {
          result->list_data[--buffer.list_data[list_data[pass3]]] =
              list_data[pass3];
        }

      result->lLength = lLength;
    }
    if (storage) {
      // this is to remove a static analysis warning
      buffer.list_data = nil;
    }
    return result;
  }
  return wantResult ? new _SimpleList : nil;
}

/**
 * @brief Clears the list.
 * @param completeClear If true, the allocated memory is freed and the list is
 * reset to its initial state. Otherwise, only the length is reset to 0.
 */
void _SimpleList::Clear(bool completeClear) {
  if (CanFreeMe()) {
    lLength = 0UL;
    if (completeClear) {
      laLength = MEMORYSTEP;
      if (is_dynamic()) {
        free(list_data);
        list_data = static_data;
      }
    }
  } else {
    RemoveAReference();
  }
}

/**
 * @brief Dumps the contents of a variable list to the console for debugging.
 * @note This function assumes that the list contains indices of variables.
 */
void _SimpleList::DebugVarList(void) {
  printf("\nVariable list dump:\n");
  for (unsigned long e = 0UL; e < lLength; e++) {
    if (list_data[e] >= 0L) {
      _Variable *theV = LocateVar(list_data[e]);
      if (theV) {
        printf("[%s]\n", theV->GetName()->get_str());
        continue;
      }
    }
    printf("[Empty]\n");
  }
}

/**
 * @brief Ensures that the list uses the correct storage type (static or
 * dynamic) based on its allocated length.
 */
void _SimpleList::_EnsureCorrectStorageType(void) {
  if (laLength >= MEMORYSTEP) {
    if (list_data == static_data) {
      list_data =
          static_cast<long *>(MemAllocate(laLength * sizeof(hyPointer)));
      _CopyStatic();

    } else {
      list_data = (long *)MemReallocate((long *)list_data,
                                        laLength * sizeof(hyPointer));
    }
  } else {
    if (list_data != static_data) {
      memcpy(static_data, list_data, lLength * sizeof(long));
      free(list_data);
      list_data = static_data;
    }
  }
}

/**
 * @brief Updates the storage type of the list if there is a significant amount
 * of unused allocated space.
 */
void _SimpleList::_UpdateStorageType(void) {
  if (laLength - lLength > MEMORYSTEP) {
    laLength -= ((laLength - lLength) / MEMORYSTEP) * MEMORYSTEP;
    _EnsureCorrectStorageType();
  }
}

/**
 * @brief Copies the contents of the static buffer to the dynamic buffer.
 */
void _SimpleList::_CopyStatic(void) {
  for (long k = 0L; k < MEMORYSTEP; k++) {
    list_data[k] = static_data[k];
  }
}

/**
 * @brief Deletes an element at a specific index.
 * @param index The index of the element to delete.
 * @param compact If true, the storage is compacted after deletion.
 */
void _SimpleList::Delete(long index, bool compact) {
  if (index_in_range(index)) {
    lLength--;
    if ((long)lLength > index) {
      for (unsigned long k = index; k < lLength; k++) {
        list_data[k] = list_data[k + 1];
      }
      // memmove
      // ((hyPointer)list_data+sizeof(BaseRef)*(index),(hyPointer)list_data+sizeof(BaseRef)*(index+1),sizeof(BaseRef)*(lLength-index));
    }
  }
  if (compact) {
    _UpdateStorageType();
  }
}

/**
 * @brief Deletes duplicate values from a sorted list.
 */
void _SimpleList::DeleteDuplicates(void) {
  if (lLength > 1L) {
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
    long last_retained = 0L, current_store = 1L;

    for (unsigned long k = 1; k < lLength; k++) {
      if (list_data[k] != list_data[last_retained]) {
        last_retained = k;
        list_data[current_store++] = list_data[k];
      }
    }

    if ((long)lLength != current_store) {
      lLength = current_store;
      TrimMemory();
    }
  }
}

/**
 * @brief Deletes a list of indices from this list.
 * @param toDelete A sorted list of indices to delete.
 */
void _SimpleList::DeleteList(const _SimpleList &toDelete) {
  if (toDelete.lLength) {
    unsigned long k = 0;
    for (unsigned long i = 0UL; i < lLength; i++) {
      if (k < toDelete.lLength && (long)i == toDelete.list_data[k])
      // if (k<toDelete.lLength)
      {
        k++;
      } else {
        list_data[i - k] = list_data[i];
      }
    }
    lLength -= toDelete.lLength;
  }
  _UpdateStorageType();
}

/**
 * @brief Corrects an index by skipping over the elements in this list.
 * @param index The index to correct.
 * @return The corrected index.
 */
long _SimpleList::SkipCorrect(long index) const {
  for (unsigned long k = 0UL; k < lLength; k++)
    if (index >= list_data[k]) {
      index++;
    }
  return index;
}

/**
 * @brief Corrects an index based on a list of excluded indices.
 * @param index The index to correct.
 * @param excluded_value The value to return if the index is in the exclusion
 * list.
 * @return The corrected index, or `excluded_value` if the index is excluded.
 */
long _SimpleList::CorrectForExclusions(long index, long excluded_value) const {
  long correction = 0L;
  for (unsigned long k = 0UL; k < lLength && index >= list_data[k]; k++) {
    if (index == list_data[k]) {
      return excluded_value;
    }
    correction++;
  }
  return index - correction;
}

/**
 * @brief Corrects an array of indices based on a list of excluded indices.
 * @param index A pointer to the array of indices to correct.
 * @param length The length of the array.
 * @return The new length of the corrected array.
 */
long _SimpleList::CorrectForExclusions(long *index, long length) const {
  long exclusion_index = 0L, mapped = 0UL;

  for (long k = 0UL; k < length; k++) {
    if (exclusion_index < (long)lLength &&
        index[k] >= list_data[exclusion_index]) {
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

/**
 * @brief Displaces a range of elements within the list.
 * @param start The starting index of the range.
 * @param end The ending index of the range.
 * @param delta The amount to displace the range by.
 */
void _SimpleList::Displace(long start, long end, long delta) {
  if (start < 0) {
    start = 0;
  } else if (start >= (long)lLength) {
    start = lLength - 1;
  }

  if (end < 0) {
    end = lLength - 1;
  } else if (end >= (long)lLength) {
    end = lLength - 1;
  }

  if ((end - start >= 0) && delta && (end - start < (long)lLength - 1)) {
    if (delta > 0 && (long)lLength - end <= delta) { // shift up
      delta = (long)lLength - end - 1;
    } else if (start - delta < 0) {
      delta = start;
    }

    if (delta) {
      long i, j, delta2 = end - start + 1;
      _SimpleList swapList((unsigned long)(end - start + 1));

      for (i = start; i <= end; i++) {
        swapList << list_data[i];
      }

      if (delta > 0) {

        for (i = end + 1; i <= end + delta; i++) {
          list_data[i - delta2] = list_data[i];
        }

      } else {
        for (i = start - 1; i >= start + delta; i--) {
          list_data[i + delta2] = list_data[i];
        }
      }

      for (i = start + delta, j = 0; i <= end + delta; i++, j++) {
        list_data[i] = swapList.list_data[j];
      }
    }
  }
}

/**
 * @brief Duplicates the contents of another SimpleList into this list.
 * @param theRef A constant reference to the SimpleList to duplicate.
 */
void _SimpleList::Duplicate(BaseRefConst theRef) {
  _SimpleList const *l = static_cast<const _SimpleList *>(theRef);
  lLength = l->lLength;
  laLength = l->laLength;
  // list_data       = l->list_data;
  if (l->is_dynamic()) {
    list_data = static_cast<long *>(MemAllocate(laLength * sizeof(hyPointer)));
    memcpy(list_data, l->list_data, lLength * sizeof(hyPointer));
  } else {
    list_data = static_data;
    for (unsigned long k = 0L; k < lLength; k++) {
      list_data[k] = l->list_data[k];
    }
  }
}

/**
 * @brief Retrieves an element from the list.
 * @param index The index of the element to retrieve. Negative indices count
 * from the end of the list.
 * @return The element at the specified index, or 0L if the index is out of
 * bounds.
 */
long _SimpleList::Element(long index) const {
  if (index_in_range(index)) {
    return list_data[index];
  }

  else if (index < 0L && -index <= (long)lLength) {
    return list_data[static_cast<long>(lLength) + index];
  }

  return 0L;
}

/**
 * @brief Compares this list with another list for equality.
 * @param l2 The other list to compare with.
 * @return True if the lists are equal (same length and all elements are the
 * same), false otherwise.
 */
bool _SimpleList::Equal(_SimpleList const &l2) const {
  if (lLength == l2.lLength) {
    for (unsigned long i = 0UL; i < lLength; i++)
      if (list_data[i] != l2.list_data[i]) {
        return false;
      }

    return true;
  }
  return false;
}

/**
 * @brief Finds the first occurrence of a value in the list.
 * @param s The value to search for.
 * @param startAt The index to start the search from.
 * @return The index of the first occurrence of the value, or kNotFound if not
 * found.
 */
long _SimpleList::Find(long s, long startAt) const {
  for (unsigned long i = startAt; i < lLength; i++) {
    if (list_data[i] == s) {
      return i;
    }
  }
  return kNotFound;
}

/**
 * @brief Finds a value in the list by stepping through elements.
 * @param s The value to search for.
 * @param step The step size to use when iterating through the list.
 * @param startAt The index to start the search from.
 * @return The index of the first occurrence of the value, or kNotFound if not
 * found.
 */
long _SimpleList::FindStepping(long s, long step, long startAt) {
  for (unsigned long i = startAt; i < lLength; i += step) {
    if (list_data[i] == s) {
      return i;
    }
  }

  return kNotFound;
}

/**
 * @brief Filters the list, keeping only elements within a specified range.
 * @param lb The lower bound of the range (inclusive).
 * @param ub The upper bound of the range (exclusive).
 */
void _SimpleList::FilterRange(long lb, long ub) {
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
    for (unsigned long k = 0; k < lLength; k++) {
      if (list_data[k] >= lb || list_data[k] <= ub) {
        list_data[current_store++] = list_data[k];
      }
    }

    if (lLength != current_store) {
      lLength = current_store;
      TrimMemory();
    }
  }
}

/**
 * @brief Reverses the order of elements in the list.
 */
void _SimpleList::Flip() {
  for (long k = 0L, l = lLength - 1; k < l; k++, l--) {
    Exchange(list_data[k], list_data[l]);
  }
}

/**
 * @brief Initializes the list.
 * @param doMemAlloc If true, memory is allocated for the list.
 */
void _SimpleList::Initialize(bool) {
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

/**
 * @brief Inserts an element into the list at a specified position.
 * @param br The element to insert.
 * @param insertAt The index at which to insert the element. If -1, the element
 * is appended.
 * @param store If true, the element is duplicated and stored. Otherwise, the
 * element is stored directly.
 * @param pointer If true, a reference to the element is added.
 */
void _SimpleList::InsertElement(BaseRef br, long insertAt, bool store,
                                bool pointer) {
  lLength++;
  if (lLength > laLength) {
    unsigned long incBy =
        ((MEMORYSTEP << 2) > lLength) ? MEMORYSTEP : (lLength >> 2);

    laLength += incBy;

    if (is_dynamic()) {
      list_data =
          (long *)MemReallocate((char *)list_data, laLength * sizeof(long));
    } else {
      list_data = static_cast<long *>(MemAllocate(laLength * sizeof(long)));
      _CopyStatic();
    }
  }

  if (insertAt == -1) {
    insertAt = lLength - 1;
  } else {
    // insertAt = insertAt>=lLength?lLength:insertAt;
    insertAt = insertAt >= (long)lLength ? lLength - 1 : insertAt;
    long moveThisMany = (laLength - insertAt - 1);
    if (moveThisMany < 32L)
      for (long k = insertAt + moveThisMany; k > insertAt; k--) {
        list_data[k] = list_data[k - 1];
      }
    else {
      memmove(((char **)list_data) + (insertAt + 1),
              ((char **)list_data) + insertAt, moveThisMany * sizeof(void *));
    }
  }
  if (store) {
    ((BaseRef *)list_data)[insertAt] = br->makeDynamic();
  } else {
    ((BaseRef *)list_data)[insertAt] = br;
    if (pointer) {
      br->AddAReference();
    }
  }
}

/**
 * @brief Converts the list into a partition-style string representation.
 * @return A BaseRef to a _StringBuffer containing the partition string.
 * @code
 *  _SimpleList list;
 *  list << 1 << 2 << 3 << 5 << 6 << 8;
 *  // list.ListToPartitionString() would return "{1-3,5-6,8}"
 * @endcode
 */
BaseRef _SimpleList::ListToPartitionString() const {
  _StringBuffer *result = new _StringBuffer((unsigned long)64);

  for (unsigned long k = 0UL; k < lLength; k++) {
    unsigned long m;
    for (m = k + 1UL; m < lLength; m++)
      if (list_data[m] - list_data[m - 1] != 1) {
        break;
      }
    if (m > k + 2UL) {
      (*result) << _String(list_data[k]) << '-' << _String(list_data[m - 1UL]);

      if (m < lLength) {
        (*result) << ',';
      }
      k = m - 1UL;
    } else {
      (*result) << _String(list_data[k]);
      if (k < lLength - 1) {
        (*result) << ',';
      }
    }
  }
  return result;
}

/**
 * @brief Creates a dynamic copy of this list.
 * @return A BaseRef to a new dynamically allocated _SimpleList object.
 */
BaseRef _SimpleList::makeDynamic(void) const {
  _SimpleList *Res = new _SimpleList;
  Res->Duplicate(this);
  return Res;
}

/**
 * @brief Finds the maximum value in the list.
 * @return The maximum value in the list.
 */
long _SimpleList::Max(void) const {
  long res = LONG_MIN;
  for (unsigned long e = 0L; e < lLength; e++)
    if (list_data[e] > res) {
      res = list_data[e];
    }
  return res;
}

/**
 * @brief Finds the minimum value in the list.
 * @return The minimum value in the list.
 */
long _SimpleList::Min(void) const {
  long res = LONG_MAX;
  for (unsigned long e = 0L; e < lLength; e++)
    if (list_data[e] < res) {
      res = list_data[e];
    }
  return res;
}

/**
 * @brief Merges two sorted lists into this list.
 * @param l1 The first sorted list.
 * @param l2 The second sorted list.
 * @param mergeResults1 An optional list to store the indices of elements from
 * l1 in the merged list.
 * @param mergeResults2 An optional list to store the indices of elements from
 * l2 in the merged list.
 */
void _SimpleList::Merge(_SimpleList &l1, _SimpleList &l2,
                        _SimpleList *mergeResults1,
                        _SimpleList *mergeResults2) {

  this->Clear();
  if (mergeResults1) {
    mergeResults1->Clear();
  }
  if (mergeResults2) {
    mergeResults2->Clear();
  }

  enum machine_states {
    INIT,
    ADVANCE1,
    ADVANCE2,
    FLUSH1,
    FLUSH2
  } advancing = INIT;

  unsigned long pos1 = 0, pos2 = 0, nt1 = l1.countitems(),
                nt2 = l2.countitems();

  bool keep_going = true;

  while (keep_going) {
    switch (advancing) {

    case ADVANCE1: { // advancing in the 1st list
      pos1++;
      if (pos1 == nt1) {
        advancing = FLUSH2;
        continue;
      }

      long cmp = l1.list_data[pos1] - l2.list_data[pos2];
      // if (l1lData[pos1] <= l2->list_data[pos2]) {
      if (cmp <= 0L) {
        if (mergeResults1) {
          // printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
          (*mergeResults1) << this->lLength;
        }

        (*this) << (l1.list_data[pos1]);

        //        if (mergeResults2 && l1->list_data[pos1] <
        //        l2->list_data[pos2]) {
        if (cmp < 0L) {
          if (mergeResults2 && mergeResults2->countitems() == 0UL &&
              pos1 >= pos2) {
            for (unsigned long i = 0UL; i < pos2; i++) {
              (*mergeResults2) << i;
            }
          }
          continue;
        }
      }

      if (cmp > 0L) {
        // if (l1->list_data[pos1] > l2->list_data[pos2]) {
        advancing = ADVANCE2;
        if (mergeResults1 && mergeResults1->countitems() == 0UL) {
          for (unsigned long i = 0UL; i < pos1; i++) {
            // printf("MERGE line number %d in file %s\n", __xLINE__, __FILE__);
            (*mergeResults1) << i;
          }
        }

        if (mergeResults2) {
          (*mergeResults2) << this->lLength;
        }

        // printf ("this->append (l2->AtIndex(pos2))\n");
        (*this) << l2.list_data[pos2];
        continue;
      }
      break;
    }

    case ADVANCE2: { // advancing in the 2nd list
      pos2++;
      if (pos2 == nt2) {
        advancing = FLUSH1;
        continue;
      }

      long cmp = l2.list_data[pos2] - l1.list_data[pos1];

      // if (l2->list_data[pos2] <= l1->list_data[pos1]) {
      if (cmp <= 0L) {
        if (mergeResults2) {
          (*mergeResults2) << this->lLength;
        }

        (*this) << l2.list_data[pos2];

        // if (l2->list_data[pos2] < l1->list_data[pos1] && mergeResults1 &&
        // !doMerge1 && pos2>=pos1) {
        if (cmp < 0L) {
          if (mergeResults1 && mergeResults1->countitems() == 0 &&
              pos2 >= pos1) {
            for (unsigned long i = 0UL; i < pos1; i++) {
              // printf("MERGE line number %d in file %s\n", __LINE__,
              // __FILE__);
              (*mergeResults1) << i;
            }
          }
          continue;
        }
      }

      // if (l2->list_data[pos2] > l1->list_data[pos1]) {
      if (cmp > 0L) {
        advancing = ADVANCE1;
        if (mergeResults2 && mergeResults2->countitems() == 0) {
          for (unsigned long i = 0UL; i < pos2; i++) {
            (*mergeResults2) << i;
          }
        }
        if (mergeResults1) {
          // printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
          (*mergeResults1) << this->lLength;
        }
        (*this) << l1.list_data[pos1];
        continue;
      }
      break;
    }
    case FLUSH2: { // flush out the 2nd list
      if (mergeResults1 && pos2 < nt2) {
        for (unsigned long i = pos1; i < nt1; i++) {
          // printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
          (*mergeResults1) << i;
        }
      }
      if (mergeResults2)
        while (pos2 < nt2) {
          (*mergeResults2) << this->lLength;
          (*this) << (l2.list_data[pos2++]);
        }
      else
        while (pos2 < nt2) {
          (*this) << (l2.list_data[pos2++]);
        }
      keep_going = false;
      break;
    }
    case FLUSH1: { // flush out the 1st list
      if (mergeResults2 && pos1 < nt1) {
        for (unsigned long i = pos2; i < nt2; i++) {
          (*mergeResults2) << i;
        }
      }
      if (mergeResults1)
        while (pos1 < nt1) {
          // printf("MERGE line number %d in file %s\n", __LINE__, __FILE__);
          (*mergeResults1) << this->lLength;
          (*this) << (l1.list_data[pos1++]);
        }
      else
        while (pos1 < nt1) {
          (*this) << (l1.list_data[pos1++]);
        }
      keep_going = false;
      break;
    }
    case INIT: {  // just starting
      if (!nt1) { // first list is empty!
        advancing = FLUSH2;
        continue;
      }
      if (!nt2) { // second list is empty!
        advancing = FLUSH1;
        continue;
      }

      // if (l1->list_data[pos1] <= l2->list_data[pos2]) { // begin with the
      // first list

      long cmp = l1.list_data[pos1] - l2.list_data[pos2];

      if (cmp <= 0L) { // begin with the first list
        if (mergeResults1) {
          (*mergeResults1) << this->lLength;
        }
        (*this) << (l1.list_data[0]);
        advancing = ADVANCE1;
        if (cmp != 0L) {
          continue;
        }
      } else {
        if (mergeResults2) {
          (*mergeResults2) << this->lLength;
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
        if (pos1 == nt1) {
          pos2++;
          if (mergeResults2) {
            (*mergeResults2) << this->lLength - 1UL;
          }
          advancing = FLUSH2;
          continue;
        } else {
          advancing = ADVANCE2;
          if (mergeResults2) {
            (*mergeResults2) << this->lLength - 1UL;
          }
        }
      } else {
        pos2++;
        if (pos2 == nt2) {
          pos1++;
          if (mergeResults1) {
            (*mergeResults1) << this->lLength - 1UL;
          }
          advancing = FLUSH1;
          continue;
        } else {
          advancing = ADVANCE1;
          if (mergeResults1) {
            (*mergeResults1) << this->lLength - 1UL;
          }
        }
      }
    }
  }
}

/**
 * @brief Initializes the state for the NChooseK algorithm.
 * @param state A SimpleList to store the state of the algorithm.
 * @param store A SimpleList to store the current combination.
 * @param stride The number of elements to choose (k).
 * @param algorithm A boolean flag (unused in this implementation).
 * @return True if initialization is successful, false otherwise.
 * @note This function is part of the implementation of algorithm NEXKSB from
 * "Combinatorial Algorithms" by Albert Nijenhuis and Herbert S. Wilf.
 */
bool _SimpleList::NChooseKInit(_SimpleList &state, _SimpleList &store,
                               unsigned long stride, bool) {
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

/**
 * @brief Generates the next combination for NChooseK algorithm.
 * @param state A SimpleList storing the state of the algorithm (modified by
 * this function).
 * @param store A SimpleList to store the current combination (modified by this
 * function).
 * @return True if a new combination is generated, false if all combinations
 * have been generated.
 * @note This function is part of the implementation of algorithm NEXKSB from
 * "Combinatorial Algorithms" by Albert Nijenhuis and Herbert S. Wilf.
 */
bool _SimpleList::NChooseK(_SimpleList &state, _SimpleList &store) {
  if (state.lLength == 1) {      // first pass
    state << 0;                  // m
    state << state.list_data[0]; // h
    state.lLength = state.list_data[0] + 3;
    store.lLength = state.list_data[0];
    if (store.lLength == 0) {
      return false;
    }
  } else {
    if (state.list_data[1] < (long)lLength - state.list_data[2]) {
      state.list_data[2] = 0;
    }
    state.list_data[2]++;
    state.list_data[1] =
        state.list_data[3 + state.list_data[0] - state.list_data[2]] + 1;
  }
  for (long j = 1; j <= state.list_data[2]; j++) {
    long anIndex = j + state.list_data[0] - state.list_data[2],
         anIndex2 = state.list_data[1] + j;
    state.list_data[anIndex + 2] = anIndex2 - 1;
    store.list_data[anIndex - 1] = list_data[anIndex2 - 1];
  }
  return state.list_data[3] < (long)lLength - state.list_data[0];
}

/**
 * @brief Normalizes 'from' and 'to' coordinates based on a reference length.
 * @param from The starting coordinate (modified by this function).
 * @param to The ending coordinate (modified by this function).
 * @param refLength The reference length to normalize against.
 */
void _SimpleList::NormalizeCoordinates(long &from, long &to,
                                       const unsigned long refLength) {
  if (to < 0) {
    to = refLength + to;
  } else {
    to = to < (long)refLength - 1 ? to : (long)refLength - 1;
  }
  if (from < 0) {
    from = refLength + from;
  }
}

/**
 * @brief Offsets all elements in the list by a given shift value.
 * @param shift The value to add to each element.
 */
void _SimpleList::Offset(long shift) {
  for (unsigned long k = 0UL; k < lLength; k++) {
    list_data[k] += shift;
  }
}

/**
 * @brief Creates a subset of the list.
 * @param size The desired size of the subset.
 * @param replacement If true, sampling is done with replacement. Otherwise,
 * without replacement.
 * @return A new SimpleList containing the subset.
 */
_SimpleList *_SimpleList::Subset(unsigned long size, bool replacement) {
  _SimpleList *result = new _SimpleList;
  if (size > 0UL && lLength > 0UL) {
    size = MIN(size, lLength);
    if (replacement) {
      for (unsigned long k = 0; k < size; k++) {
        (*result) << list_data[genrand_int32() % lLength];
      }
    } else {
      (*result) << (*this);
      for (unsigned long k = 0; k < size; k++) {
        long idx = list_data[genrand_int32() % (lLength - k)];
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

/**
 * @brief Creates a random permutation of the list's elements.
 * @param blockLength The size of blocks to permute. If 1, individual elements
 * are permuted.
 */
void _SimpleList::Permute(long blockLength) {

  unsigned long blockCount = lLength / blockLength;

  if (blockCount > 1) {
    if (blockLength > 1) {
      for (unsigned long k = 0; k < blockCount - 1; k = k + 1) {
        unsigned long k2 = (unsigned long)(genrand_real2() * (blockCount - k));
        if (k2) {
          k2 += k;
          k2 *= blockLength;

          for (long j = 0; j < blockLength; j++) {
            Exchange(list_data[k2 + j], list_data[k * blockLength + j]);
          }
        }
      }

    } else {
      for (unsigned long k = 0L; k < blockCount - 1; k = k + 1) {
        unsigned long k2 = (unsigned long)(genrand_real2() * (blockCount - k));
        if (k2) {
          k2 += k;
          Exchange(list_data[k2], list_data[k]);
        }
      }
    }
  }
}

/**
 * @brief Returns a random index from the list.
 * @return A random index within the list's bounds, or kNotFound if the list is
 * empty.
 */
long _SimpleList::Choice() const {
  if (lLength) {
    return genrand_int32() % lLength;
  }
  return kNotFound;
}

/**
 * @brief Samples elements from the list with replacement.
 * @param size The number of elements to sample.
 * @return A new SimpleList containing the sampled elements.
 */
_SimpleList const _SimpleList::Sample(unsigned long size) const {
  // TODO SLKP 20171026: this is new, need to check correctness
  if (size >= lLength) {
    return *this;
  }

  _SimpleList result(size, 0, 0), tracker(size, 0, 0);

  for (unsigned long k = 0; k < size; k++) {
    unsigned long k2 = k + (unsigned long)(genrand_real2() * (lLength - k));
    result.list_data[k] = list_data[k2];
    tracker[k] = k2;
    if (k2 != k) {
      Exchange(list_data[k2], list_data[k]);
    }
  }

  // reshuffle moved elements LIFO to ensure correct ordering in *this

  for (long k = size - 1; k >= 0; k--) {
    Exchange(list_data[k], list_data[tracker.get(k)]);
  }

  return result;
}

/**
 * @brief Creates a permutation of the list's elements with possible
 * repetitions.
 * @param blockLength The size of blocks to permute.
 */
void _SimpleList::PermuteWithReplacement(long blockLength) {
  unsigned long blockCount = lLength / blockLength;
  _SimpleList result(static_cast<unsigned long>(blockCount * blockLength));
  if (blockLength > 1)
    for (unsigned long i = 0; i < blockCount; i++) {
      unsigned long sample =
          static_cast<unsigned long>(genrand_real2() * blockCount);
      sample *= blockLength;
      for (long j = 0; j < blockLength; j++, sample++) {
        result << list_data[sample];
      }
    }
  else {
    for (unsigned long i = 0; i < blockCount; i++) {
      unsigned long sample = (unsigned long)(genrand_real2() * blockCount);
      result << list_data[sample];
    }
  }

  Clear();
  Duplicate(&result);
}

/**
 * @brief Removes and returns an element from the end of the list.
 * @param discard The number of elements to discard from the end before popping.
 * @return The popped element, or 0L if the list is empty or `discard` is too
 * large.
 */
long _SimpleList::Pop(unsigned long discard) {
  if (lLength > discard) {
    return list_data[lLength -= (discard + 1UL)];
  }
  return 0L;
}

/**
 * @brief Populates the list with a sequence of numbers.
 * @param l The number of elements to populate.
 * @param start The starting value of the sequence.
 * @param step The increment between consecutive elements.
 */
void _SimpleList::Populate(long l, long start, long step) {
  RequestSpace(l);

  if (start == 0 && step == 0 && l >= 128) {
    memset(list_data, 0, sizeof(long) * l);
  } else {
    for (long k = 0L; k < l; k++, start += step) {
      list_data[k] = start;
    }
  }
  lLength = l;
}

/**
 * @brief Appends a range of numbers to the list.
 * @param how_many The number of elements to append.
 * @param start The starting value of the sequence.
 * @param step The increment between consecutive elements.
 */
void _SimpleList::AppendRange(unsigned long how_many, long start, long step) {
  RequestSpace(how_many + lLength);
  for (unsigned long k = 0UL; k < how_many; k++, start += step) {
    list_data[lLength + k] = start;
  }
  lLength = how_many + lLength;
}

/**
 * @brief Recursively sorts a portion of the list and updates an index list
 * accordingly.
 * @param from The starting index of the portion to sort.
 * @param to The ending index of the portion to sort.
 * @param index A pointer to a SimpleList that stores the original indices of
 * the elements.
 */
void _SimpleList::RecursiveIndexSort(long from, long to, _SimpleList *index) {
  long middle = (from + to) >> 1, middleV = list_data[middle], bottommove = 1,
       topmove = 1, i, imiddleV = (*index)(middle);

  if (middle)
    while ((middle - bottommove >= from) &&
           (Compare(middle - bottommove, middle) == kCompareGreater)) {
      bottommove++;
    }

  if (from < to)
    while ((middle + topmove <= to) &&
           (Compare(middle + topmove, middle) == kCompareLess)) {
      topmove++;
    }

  // now shuffle
  for (i = from; i < middle - bottommove; i++) {
    if (Compare(i, middle) != kCompareLess) {
      Exchange(list_data[middle - bottommove], list_data[i]);
      Exchange(index->list_data[middle - bottommove], index->list_data[i]);
      bottommove++;
      while ((middle - bottommove >= from) &&
             (Compare(middle - bottommove, middle) == kCompareGreater)) {
        bottommove++;
      }
    }
  }

  for (i = middle + topmove + 1; i <= to; i++) {
    if (Compare(i, middle) != kCompareGreater) {
      Exchange(list_data[middle + topmove], list_data[i]);
      Exchange(index->list_data[middle + topmove], index->list_data[i]);
      topmove++;
      while ((middle + topmove <= to) &&
             (Compare(middle + topmove, middle) == kCompareLess)) {
        topmove++;
      }
    }
  }

  if (topmove == bottommove) {
    for (i = 1; i < bottommove; i++) {
      Exchange(list_data[middle + i], list_data[middle - i]);
      Exchange(index->list_data[middle + i], index->list_data[middle - i]);
    }
  } else if (topmove > bottommove) {
    long shift = topmove - bottommove;
    for (i = 1; i < bottommove; i++) {
      Exchange(list_data[middle + i + shift], list_data[middle - i]);
      Exchange(index->list_data[middle + i + shift],
               index->list_data[middle - i]);
    }

    for (i = 0; i < shift; i++) {
      list_data[middle + i] = list_data[middle + i + 1];
      index->list_data[middle + i] = index->list_data[middle + i + 1];
    }
    middle += shift;
    list_data[middle] = middleV;
    index->list_data[middle] = imiddleV;
  } else {
    long shift = bottommove - topmove;
    for (i = 1; i < topmove; i++) {
      Exchange(list_data[middle - i - shift], list_data[middle + i]);
      Exchange(index->list_data[middle - i - shift],
               index->list_data[middle + i]);
    }
    for (i = 0; i < shift; i++) {
      list_data[middle - i] = list_data[middle - i - 1];
      index->list_data[middle - i] = index->list_data[middle - i - 1];
    }
    middle -= shift;
    list_data[middle] = middleV;
    index->list_data[middle] = imiddleV;
  }

  if (to > middle + 1) {
    RecursiveIndexSort(middle + 1, to, index);
  }
  if (from < middle - 1) {
    RecursiveIndexSort(from, middle - 1, index);
  }
}

/**
 * @brief Requests additional space for the list.
 * @param slots The number of slots to request.
 */
void _SimpleList::RequestSpace(unsigned long slots) {
  if (slots > laLength) {
    laLength = (slots / MEMORYSTEP + 1) * MEMORYSTEP;
    _EnsureCorrectStorageType();
  }
}

/**
 * @brief Computes the set difference of two sorted lists (this = l1 - l2).
 * @param l1 The first list.
 * @param l2 The second list (elements to subtract).
 */
void _SimpleList::Subtract(_SimpleList const &l1, _SimpleList const &l2) {
  if (lLength) {
    Clear();
  }

  auto add_to_result = [this](long element) -> void { (*this) << element; };
  auto noop = [](long) -> void {};

  IterateOverTwoSortedLists(l1, l2, noop, add_to_result, noop);
}

/**
 * @brief Swaps two elements in the list.
 * @param i The index of the first element.
 * @param j The index of the second element.
 */
void _SimpleList::Swap(long i, long j) {
  if (i >= (long)lLength || j >= (long)lLength) {
    return;
  }
  Exchange(list_data[i], list_data[j]);
}

/**
 * @brief Converts the list to a string representation.
 * @param unused An unused parameter.
 * @return A BaseRef to a _StringBuffer containing the string representation of
 * the list.
 * @code
 *  _SimpleList list;
 *  list << 1 << 2 << 3;
 *  // list.toStr() would return "{1,2,3}"
 * @endcode
 */
BaseRef _SimpleList::toStr(unsigned long) {
  if (lLength) {
    unsigned long ma = (unsigned long)(lLength * (1. + log10((double)lLength)));

    _StringBuffer *s = new _StringBuffer(MAX(32UL, ma));

    (*s) << "{";

    char c[32];
    for (unsigned long i = 0UL; i < lLength - 1L; i++) {
      snprintf(c, sizeof(c), "%ld", list_data[i]);
      (*s) << c << ',';
    }
    snprintf(c, sizeof(c), "%ld", list_data[lLength - 1L]);
    (*s) << c << '}';

    return s;
  } else {
    return new _String("{}");
  }
}

/**
 * @brief Trims the allocated memory to match the current length of the list, or
 * to MEMORYSTEP if the list is smaller.
 */
void _SimpleList::TrimMemory(void) {
  if (laLength > lLength) {
    laLength = Maximum(lLength, (unsigned long)MEMORYSTEP);
    _EnsureCorrectStorageType();
  }
}

/*
==============================================================
Sort Methods
==============================================================
*/

/**
 * @brief Sorts the elements in the list.
 * @param ascending If true, sorts in ascending order. Otherwise, sorts in
 * descending order.
 */
void _SimpleList::Sort(bool ascending) {
  if (lLength > 0UL) {
    if (lLength < 10UL) { // use bubble sort
      BubbleSort();
    } else {
      QuickSort(0UL, lLength - 1UL);
    }

    if (!ascending) {
      for (unsigned long i = 0UL, j = lLength - 1UL; i < j; i++, j--) {
        Exchange(list_data[i], list_data[j]);
      }
    }
  }
}

/**
 * @brief Sorts the elements in the list using the Bubble Sort algorithm.
 */
void _SimpleList::BubbleSort(void) {
  bool done = lLength == 0UL;
  while (!done) {
    done = true;
    unsigned long upper_bound = lLength - 1L;
    for (unsigned long i = 0UL; i < upper_bound; i++) {
      if (Compare(i, i + 1) == kCompareGreater) {
        done = false;
        Exchange(list_data[i], list_data[i + 1]);
        upper_bound = i;
      }
    }
  }
}

/**
 * @brief Sorts a portion of the list using the Quick Sort algorithm.
 * @param from The starting index of the portion to sort.
 * @param to The ending index of the portion to sort.
 */
/**
 * @brief Partitions a sub-array using the Hoare partition scheme.
 * @param from The starting index of the sub-array.
 * @param to The ending index of the sub-array.
 * @return An index that splits the array into two partitions. Note: This index
 * is NOT the final sorted position of the pivot. All elements to the
 * left of the returned index are <= the pivot, and all elements to the
 * right are >= the pivot.
 */
long _SimpleList::Partition(long from, long to) {
  // We'll use the middle element as the pivot to be consistent with the
  // original and to avoid worst-case behavior with already-sorted data.
  long pivot_idx = from + (to - from) / 2;
  long i = from - 1;
  long j = to + 1;

  while (true) {
    // Find an element on the left side that should be on the right side
    // (i.e., an element >= pivot)
    do {
      i++;
    } while (Compare(i, pivot_idx) == kCompareLess);

    // Find an element on the right side that should be on the left side
    // (i.e., an element <= pivot)
    do {
      j--;
    } while (Compare(j, pivot_idx) == kCompareGreater);

    // If the pointers have crossed, the partition is done.
    if (i >= j) {
      return j; // Return the split point
    }

    // Swap the two elements that are in the wrong partition.
    // Important: If the pivot itself is swapped, we must update its index.
    if (i == pivot_idx) {
      pivot_idx = j;
    } else if (j == pivot_idx) {
      pivot_idx = i;
    }
    Exchange(list_data[i], list_data[j]);
  }
}

/**
 * @brief Sorts the list using the QuickSort algorithm with Hoare partitioning.
 * @param from The starting index for the sort.
 * @param to The ending index for the sort.
 */
void _SimpleList::QuickSort(long from, long to) {
  // Base case: if the segment has 0 or 1 elements, it's already sorted.
  if (from < to) {
    // Get the split point from the partition function.
    long pi = Partition(from, to);

    // Recursively sort the two sub-arrays.
    // Because Hoare's scheme doesn't place the pivot at its final spot,
    // the partition index 'pi' is included in the first sub-array.
    QuickSort(from, pi);
    QuickSort(pi + 1, to);
  }
}

/**
 * @brief Sorts a list and updates an index list accordingly.
 * @param ref A pointer to the SimpleList to sort.
 * @param index A pointer to a SimpleList that stores the original indices of
 * the elements.
 */
void SortLists(_SimpleList *ref, _SimpleList *index) {
  if (ref->lLength > 0UL) {
    if (ref->lLength != index->lLength) {
      return;
    }
    if ((*ref).lLength <= 10UL) {
      bool done = false;

      while (!done) {
        done = true;
        unsigned long upper_bound = ref->lLength - 1L;
        for (unsigned long i = 0UL; i < upper_bound; i++) {
          if (ref->Compare(i, i + 1) == kCompareGreater) {
            done = false;
            Exchange(ref->list_data[i], ref->list_data[i + 1]);
            Exchange(index->list_data[i], index->list_data[i + 1]);
            upper_bound = i;
          }
        }
      }
    } else {
      (*ref).RecursiveIndexSort(0, (*ref).lLength - 1UL, index);
    }
  }
}

/*
==============================================================
Set Methods
==============================================================
*/

/**
 * @brief Computes the intersection of two sorted lists.
 * @param l1 The first sorted list.
 * @param l2 The second sorted list.
 */
void _SimpleList::Intersect(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }

  auto add_to_result = [this](long element) -> void { (*this) << element; };
  auto noop = [](long) -> void {};

  IterateOverTwoSortedLists(l1, l2, add_to_result, noop, noop);
}

/**
 * @brief Computes the union of two sorted lists. Each element appears exactly
 * once.
 * @param l1 The first sorted list.
 * @param l2 The second sorted list.
 */
void _SimpleList::Union(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }
  auto add_to_result = [this](long element) -> void { (*this) << element; };
  IterateOverTwoSortedLists(l1, l2, add_to_result, add_to_result,
                            add_to_result);
}

/**
 * @brief Computes the symmetric difference (XOR) of two sorted lists.
 * @param l1 The first sorted list.
 * @param l2 The second sorted list.
 */
void _SimpleList::XOR(_SimpleList &l1, _SimpleList &l2) {
  if (lLength) {
    Clear();
  }
  auto add_to_result = [this](long element) -> void { (*this) << element; };
  auto noop = [](long) -> void {};

  IterateOverTwoSortedLists(l1, l2, noop, add_to_result, add_to_result);
}

/**
 * @brief Maps an index to a value in the list.
 * @param index The index to map.
 * @param map_failed The value to return if the index is out of bounds.
 * @return The value at the given index, or `map_failed` if the index is out of
 * bounds.
 */
long _SimpleList::Map(long index, long map_failed) const {
  return index_in_range(index) ? list_data[index] : map_failed;
}
