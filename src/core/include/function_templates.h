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

#pragma once

/**
 * @file function_templates.h
 * @brief A collection of useful function templates.
 */

#include "mersenne_twister.h"
#include "parser.h"
#include "variablecontainer.h"

/**
 * @brief Check for a parameter in a variable container and assign a default value if not found.
 *
 * @tparam ARG_TYPE The type of the parameter.
 * @param name The name of the parameter.
 * @param dest The destination variable to store the parameter value.
 * @param def The default value if the parameter is not found.
 * @param theP The variable container to search in.
 */
template <typename ARG_TYPE>
void checkParameter(_String const &name, ARG_TYPE &dest, const ARG_TYPE def,
                    const _VariableContainer *theP = nil) {
  _Variable *v = FetchVar(
      LocateVarByName(WrapInNamespace(name, theP ? theP->GetName() : nil)));
  dest = v ? (ARG_TYPE)v->Value() : def;
}

/**
 * @brief Store a value if it is greater than the current maximum.
 *
 * @tparam ARG_TYPE The type of the values.
 * @param current_max The current maximum value.
 * @param value_to_check The value to check.
 * @return true if the value was stored, false otherwise.
 */
template <typename ARG_TYPE>
bool StoreIfGreater(ARG_TYPE &current_max, ARG_TYPE const &value_to_check) {
  if (value_to_check > current_max) {
    current_max = value_to_check;
    return true;
  }
  return false;
}

/**
 * @brief Returns the maximum of two values.
 *
 * @tparam ARG_TYPE The type of the values.
 * @param a The first value.
 * @param b The second value.
 * @return The maximum of the two values.
 */
template <typename ARG_TYPE>
ARG_TYPE Maximum(ARG_TYPE const a, ARG_TYPE const b) {
  if (a > b) {
    return a;
  }
  return b;
}

/**
 * @brief Returns the minimum of two values.
 *
 * @tparam ARG_TYPE The type of the values.
 * @param a The first value.
 * @param b The second value.
 * @return The minimum of the two values.
 */
template <typename ARG_TYPE>
ARG_TYPE Minimum(ARG_TYPE const a, ARG_TYPE const b) {
  if (a > b) {
    return b;
  }
  return a;
}

/**
 * @brief Exchanges the values of two variables.
 *
 * @tparam ARG_TYPE The type of the variables.
 * @param a The first variable.
 * @param b The second variable.
 */
template <typename ARG_TYPE> void Exchange(ARG_TYPE &a, ARG_TYPE &b) {
  ARG_TYPE t = a;
  a = b;
  b = t;
}

/**
 * @brief Store a value if it is less than the current minimum.
 *
 * @tparam ARG_TYPE The type of the values.
 * @param current_min The current minimum value.
 * @param value_to_check The value to check.
 * @return true if the value was stored, false otherwise.
 */
template <typename ARG_TYPE>
bool StoreIfLess(ARG_TYPE &current_min, ARG_TYPE const &value_to_check) {
  if (value_to_check < current_min) {
    current_min = value_to_check;
    return true;
  }
  return false;
}

/**
 * @brief Computes the power of a number.
 *
 * @tparam ARG_TYPE The type of the number.
 * @param base The base.
 * @param exponent The exponent.
 * @return The result of the power operation.
 */
template <typename ARG_TYPE>
ARG_TYPE ComputePower(ARG_TYPE base, unsigned long exponent) {
  ARG_TYPE result = 1.;
  unsigned long mask = 1L << (sizeof(unsigned long) * 8 - 2);
  // left shift to left-most position of binary sequence for long integer
  // e.g. 100...0 (30 zeroes for signed long)

  while (exponent && ((exponent & mask) == 0)) {
    mask >>= 1; // bitwise AND, right-shift mask until overlaps with first '1'
  }

  while (mask) {
    result *= result;
    if (exponent & mask) {
      result *= base;
    }
    mask >>= 1;
  }
  return result;
}

/**
 * @brief Checks if any element in an array satisfies a condition.
 *
 * @tparam ARG_TYPE The type of the array elements.
 * @tparam LAMBDA The type of the condition.
 * @param array The array to check.
 * @param dimension The dimension of the array.
 * @param condition The condition to check.
 * @return true if any element satisfies the condition, false otherwise.
 */
template <typename ARG_TYPE, typename LAMBDA>
bool ArrayAny(ARG_TYPE const *array, unsigned long dimension,
              LAMBDA &&condition) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    if (condition(array[i], i)) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Checks if all elements in an array satisfy a condition.
 *
 * @tparam ARG_TYPE The type of the array elements.
 * @tparam LAMBDA The type of the condition.
 * @param array The array to check.
 * @param dimension The dimension of the array.
 * @param condition The condition to check.
 * @return true if all elements satisfy the condition, false otherwise.
 */
template <typename ARG_TYPE, typename LAMBDA>
bool ArrayAll(ARG_TYPE const *array, unsigned long dimension,
              LAMBDA &&condition) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    if (!condition(array[i], i)) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Checks if any element in a list satisfies a condition.
 *
 * @tparam LAMBDA The type of the condition.
 * @param list The list to check.
 * @param condition The condition to check.
 * @return true if any element satisfies the condition, false otherwise.
 */
template <typename LAMBDA> bool ListAny(_SimpleList &list, LAMBDA &&condition) {
  unsigned long const list_length = list.countitems();
  for (unsigned long i = 0UL; i < list_length; i++) {
    if (condition(list.get(i), i)) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Applies a transformation to each element of an array.
 *
 * @tparam ARG_TYPE The type of the array elements.
 * @tparam LAMBDA The type of the transformation.
 * @param array The array to transform.
 * @param dimension The dimension of the array.
 * @param transform The transformation to apply.
 */
template <typename ARG_TYPE, typename LAMBDA>
void ArrayForEach(ARG_TYPE *array, unsigned long dimension,
                  LAMBDA &&transform) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    array[i] = transform(array[i], i);
  }
}

/**
 * @brief Initializes an array with a given value.
 *
 * @tparam ARG_TYPE The type of the array elements.
 * @param array The array to initialize.
 * @param dimension The dimension of the array.
 * @param value The value to initialize the array with.
 */
template <typename ARG_TYPE>
void InitializeArray(ARG_TYPE *array, unsigned long dimension,
                     ARG_TYPE &&value) {
// #pragma clang loop unroll_count(8)
#pragma GCC unroll 4
  for (unsigned long i = 0UL; i < dimension; i++) {
    array[i] = value;
  }
}

/**
 * @brief Copies an array.
 *
 * @tparam ARG_TYPE The type of the array elements.
 * @param to The destination array.
 * @param from The source array.
 * @param dimension The dimension of the arrays.
 */
template <typename ARG_TYPE>
void CopyArray(ARG_TYPE *to, ARG_TYPE const *from, unsigned long dimension) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    to[i] = from[i];
  }
}

/**
 * @brief Copies an array with an offset.
 *
 * @tparam ARG_TYPE The type of the array elements.
 * @param to The destination array.
 * @param from The source array.
 * @param dimension The dimension of the arrays.
 * @param offset The offset to apply to the source array.
 */
template <typename ARG_TYPE>
void CopyArrayWithOffset(ARG_TYPE *to, ARG_TYPE const *from,
                         unsigned long dimension, long offset) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    to[i] = from[i + offset];
  }
}

/**
 * @brief Deconstruct a number into 'places' digits according to the supplied radix
 *
 * @tparam ARG_TYPE The type of the number.
 * @param composition The number to deconstruct.
 * @param places The number of places.
 * @param radix The radix to use.
 * @return _SimpleList The list of digits.
 *
 * @example SplitIntoDigits (5,3,2) will return (higher to lower significe digits)
 * 1,0,1 (e.g. 101 in binary)
 */
template <typename ARG_TYPE>
const _SimpleList SplitIntoDigits(ARG_TYPE composition, unsigned long places,
                                  unsigned long radix) {
  _SimpleList result(places, 0, 0);

  ARG_TYPE remainder = composition;
  unsigned long index = 0;

  while (remainder > 0 && index < places) {
    result.list_data[places - index - 1] = (remainder % radix);
    remainder /= radix;
  }

  return result;
}

/**
 * @brief Reconstruct a number from digits according to the supplied radix.
 *
 * @tparam ARG_TYPE The type of the number.
 * @param digits The digits to use.
 * @param places The number of places.
 * @param radix The radix to use.
 * @return The reconstructed number.
 *
 * @example CombineDigits ([5,3,2], 3, 4) will return
 * 2 + 3*4 + 5*16 = 94
 */
template <typename ARG_TYPE>
const ARG_TYPE CombineDigits(ARG_TYPE const *digits, unsigned long places,
                             unsigned long radix) {

  ARG_TYPE number = 0, multiplier = 1;

  for (long digit = places - 1; digit >= 0L; digit--) {
    number += multiplier * digits[digit];
    multiplier *= radix;
  }

  return number;
}

/**
 * @brief Delete (or decrease ref count) of an object pointer. If the object was deleted, set the pointer to NULL.
 *
 * @tparam ARG_TYPE The type of the object.
 * @param object The object to delete.
 */
template <typename ARG_TYPE> void DeleteAndZeroObject(ARG_TYPE &object) {
  DeleteObject(object);
  object = NULL;
}

/**
 * @brief Draw a random index from a discrete distribution.
 *
 * @tparam ARG_TYPE The type of the probabilities.
 * @param cdf The cumulative distribution function.
 * @param dimension The dimension of the distribution.
 * @return The drawn index.
 */
template <typename ARG_TYPE>
unsigned long DrawFromDiscrete(ARG_TYPE const *cdf, unsigned long dimension) {
  unsigned long index = 0UL;
  ARG_TYPE sum_so_far = cdf[0], random_draw = genrand_real2();

  while (sum_so_far < random_draw && index < dimension) {
    sum_so_far += cdf[++index];
  }

  return index;
}

/**
 * @brief Draw a random index from a discrete distribution using a generator function.
 *
 * @tparam FUNCTOR The type of the generator function.
 * @param generator The generator function.
 * @param dimension The dimension of the distribution.
 * @return The drawn index.
 */
template <typename FUNCTOR>
unsigned long DrawFromDiscreteGenerator(FUNCTOR &&generator,
                                        unsigned long dimension) {
  unsigned long index = 0UL;

  auto sum_so_far = generator(0);
  double random_draw = genrand_real2();

  while (sum_so_far < random_draw && index < dimension) {
    sum_so_far += generator(++index);
  }

  return index;
}

/**
 * @brief Deletes a single object.
 *
 * @tparam ARG_TYPE The type of the object.
 * @param first The object to delete.
 */
template <typename ARG_TYPE> void BatchDelete(ARG_TYPE first) { delete first; }

/**
 * @brief Deletes a single object using DeleteObject.
 *
 * @tparam ARG_TYPE The type of the object.
 * @param first The object to delete.
 */
template <typename ARG_TYPE> void BatchDeleteObject(ARG_TYPE first) {
  DeleteObject(first);
}

/**
 * @brief Deletes multiple objects.
 *
 * @tparam ARG_TYPE The type of the first object.
 * @tparam Args The types of the other objects.
 * @param first The first object to delete.
 * @param args The other objects to delete.
 */
template <typename ARG_TYPE, typename... Args>
void BatchDelete(ARG_TYPE first, const Args &...args) {
  delete first;
  BatchDelete(args...);
}

/**
 * @brief Deletes multiple objects using DeleteObject.
 *
 * @tparam ARG_TYPE The type of the first object.
 * @tparam Args The types of the other objects.
 * @param first The first object to delete.
 * @param args The other objects to delete.
 */
template <typename ARG_TYPE, typename... Args>
void BatchDeleteObject(ARG_TYPE first, const Args &...args) {
  DeleteObject(first);
  BatchDeleteObject(args...);
}

/**
 * @brief Deletes a single array.
 *
 * @tparam ARG_TYPE The type of the array.
 * @param first The array to delete.
 */
template <typename ARG_TYPE> void BatchDeleteArray(ARG_TYPE first) {
  delete[] first;
}

/**
 * @brief Deletes multiple arrays.
 *
 * @tparam ARG_TYPE The type of the first array.
 * @tparam Args The types of the other arrays.
 * @param first The first array to delete.
 * @param args The other arrays to delete.
 */
template <typename ARG_TYPE, typename... Args>
void BatchDeleteArray(ARG_TYPE first, const Args &...args) {
  delete[] first;
  BatchDeleteArray(args...);
}
