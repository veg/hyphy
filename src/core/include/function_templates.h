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


template <typename ARG_TYPE>
void        checkParameter  (_String const& name, ARG_TYPE& dest, const ARG_TYPE def, const _VariableContainer* theP = nil) {
  _Variable *v = FetchVar(LocateVarByName(WrapInNamespace(name, theP ? theP->GetName() : nil)));
  dest = v ? v->Value() : def;
}

template <typename ARG_TYPE>
bool        StoreIfGreater (ARG_TYPE& current_max, ARG_TYPE const & value_to_check) {
  if (value_to_check > current_max) {
    current_max = value_to_check;
    return true;
  }
  return false;
}

template <typename ARG_TYPE>
ARG_TYPE        Maximum (ARG_TYPE const a, ARG_TYPE const b) {
  if (a > b) {
    return a;
  }
  return b;
}

template <typename ARG_TYPE>
ARG_TYPE        Minimum (ARG_TYPE const a, ARG_TYPE const b) {
  if (a > b) {
    return b;
  }
  return a;
}

template <typename ARG_TYPE>
bool        StoreIfLess (ARG_TYPE& current_min, ARG_TYPE const & value_to_check) {
  if (value_to_check < current_min) {
    current_min = value_to_check;
    return true;
  }
  return false;
}

template <typename ARG_TYPE>
ARG_TYPE        ComputePower  (ARG_TYPE base, unsigned long exponent) {
    ARG_TYPE    result = 1;
    unsigned long mask   = 1L<<(sizeof(unsigned long)*8-2);
            // left shift to left-most position of binary sequence for long integer
            // e.g. 100...0 (30 zeroes for signed long)
    
    while ((exponent & mask) == 0) {
      mask >>= 1;    // bitwise AND, right-shift mask until overlaps with first '1'
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

template <typename ARG_TYPE, typename LAMBDA>
bool      ArrayAny (ARG_TYPE const* array, unsigned long dimension, LAMBDA&& condition) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    if (condition (array[i],i)) {
      return true;
    }
  }
  return false;
}

template <typename ARG_TYPE, typename LAMBDA>
void      ArrayForEach (ARG_TYPE* array, unsigned long dimension, LAMBDA&& transform) {
  for (unsigned long i = 0UL; i < dimension; i++) {
    transform (array[i], i);
  }
}


template <typename ARG_TYPE>
void InitializeArray (ARG_TYPE* array, unsigned long dimension, ARG_TYPE&& value) {
  for (unsigned long i = 0UL; i < dimension; i++) {
     array[i] = value;
  }
}

template <typename ARG_TYPE>
const _SimpleList SplitIntoDigits (ARG_TYPE composition, unsigned long places, unsigned long radix) {
  /**
   Deconstruct a number into 'places' digits according to the supplied radix
   
   SplitIntoDigits (5,3,2) will return (higher to lower significe digits)
   1,0,1 (e.g. 101 in binary)
   
   */
  
  _SimpleList result (places, 0, 0);
  
  ARG_TYPE remainder = composition;
  unsigned long index = 0;
  
  while (remainder > 0 && index < places) {
    result.lData[places-index-1] = (remainder % radix);
    remainder /= radix;
  }
  
  return result;
}

template <typename ARG_TYPE>
const ARG_TYPE CombineDigits (ARG_TYPE const* digits, unsigned long places, unsigned long radix) {
  /**
   Reconstruct a number from digits according to the supplied radix.
   
   CombineDigits ([5,3,2], 3, 4) will return 
    2 + 3*4 + 5*16 = 94
   */
  
  
  ARG_TYPE number = 0,
           multiplier = 1;
  
  for (long digit = places-1; digit >= 0L; digit --  ) {
    number += multiplier * digits[digit];
    multiplier *= radix;
  }
  
  return number;
}

template <typename ARG_TYPE>
unsigned long DrawFromDiscrete (ARG_TYPE const * cdf, unsigned long dimension) {
  /** 
    assuming that cdf is an array of probabilities summing to 1,
    draw a random index from the distribution
   
  */
  
  unsigned long index  = 0UL;
  ARG_TYPE sum_so_far  = cdf[0],
           random_draw = genrand_real2 ();
  
  while (sum_so_far < random_draw && index < dimension) {
    sum_so_far += cdf[++index];
  }
  
  return index;
}
