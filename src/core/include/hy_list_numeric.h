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

#ifndef _HY_LIST_NUMERIC_
#define _HY_LIST_NUMERIC_
//#pragma once
#include "hy_list_orderable.h"

/*
  
  A resizable list which stores PAYLOAD types
  
  assumes that PAYLOAD understands the equal (==), assign (=), 
  and stack copy operations, as well as all ordering operations 
  (<,>,<=,>=), the addition/subtraction operation, and can be 
  cast to from _HY_LIST_NUMERIC_INVALID_VALUE_ and 0L

*/ 

#define _HY_LIST_NUMERIC_INVALID_VALUE_ 0L


template <typename PAYLOAD>
class _hyListNumeric : public _hyListOrderable<PAYLOAD> {
  public:
  
        //does nothing
      _hyListNumeric();
      
        // stack copy contructor
      _hyListNumeric(const _hyListNumeric <PAYLOAD>&, const long = 0UL, const long = HY_LIST_INSERT_AT_END);
      
        // data constructor (1 member list)
      _hyListNumeric(const PAYLOAD);

        // data constructor from _hyListOrderable
  
      _hyListNumeric(const _hyListOrderable<PAYLOAD>&);

    // arithmetic series populator: size, first item, step
      _hyListNumeric(const unsigned long, const PAYLOAD, const PAYLOAD);


      /**
       * Data constructor from a static list of PAYLOAD objects
       * @param const unsigned long the number of PAYLOAD arguments supplied
       * to the constructor in the second argument
       * @param const PAYLOAD: the static list of items to pass to the constructor
       */
      
      _hyListNumeric (const unsigned long, const PAYLOAD []);

      /**
      * Populate a _hyListNumeric with integers incrementally.
      * Example: sl.Populate(4, 1, 2) = [1, 3, 5, 7]
      * @param s The substring to find
      * @param startat The index to start searching from
      * @param increment by Pass true for a case sensitive search
      * @return Nothing. Acts on the List object it was called from.
      */
      void Populate(const unsigned long, const PAYLOAD, const PAYLOAD);


      /* instead of throwing an error, returns the missing value */
      virtual PAYLOAD Element(const long) const;

      /**
      * SLKP: 20090508
      * Return the sum of all values in the list
      * Example: _SimpleList([4, 1, 2]).Sum() = 7
      * @return the sum of all values in the list
      */
      PAYLOAD Sum(void) const;


      /**
      * Add a value to each entry in the array
      * Example: _hyListNumeric(1,3,5,7).Offset(2) = [3, 5, 7, 9]
      * @param shift Value to add
      * @return Nothing. Acts on the List object it was called from.
      */
      void Offset(const PAYLOAD);

  
  
  
      /**
      * Convert (a sorted) list into a partition string for consumption by datafilters etc
      * Example: _hyListNumeric(1,2,3,7,8,9,10).ListToPartitionString() = "1-3,7-10"
      * @return The partition string
      */
      // define this only for "long" PAYLOADS
      _StringBuffer * ListToPartitionString(void) const;


      /**
      * SLKP: 20090508
      * Implements a counting sort procedure, ASSUMING that all
      * list values are in [0, upperBound-1]; if the 1st argument is 
      * _HY_LIST_NUMERIC_INVALID_VALUE_, it is
      * automatically determined. This only really makes sense for integer types
      * @return a pointer to the sorted list
      * if the second argument is not nil, then
      * the new_order->old_order mapping is returned in the array pointed to
      *
      */
      
      _hyListNumeric <PAYLOAD> CountingSort(PAYLOAD, _hyListNumeric <long> * = nil);
      virtual BaseRef toStr() const;
  
private:
  
      /**
       * Appends a number to the passed string buffer
       * @return void
       */
      void AppendNumtoStr(_StringBuffer*, const PAYLOAD) const;

};

  //#include "hy_list_numeric.cpp"

  // FUNCTION DEFINITIONS -- moved here so that gcov can view them

// Does nothing
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric() {
}

  //Data constructor (1 member list)
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const PAYLOAD item) : _hyListOrderable<PAYLOAD> (item) {
}


  //Stack copy contructor
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const _hyListNumeric <PAYLOAD> &l, const long from, const long to)
{
  this->Clone(&l, from, to); }

  // Data constructor from array
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const unsigned long number, const PAYLOAD items[]) :
_hyListOrderable<PAYLOAD> (number, items) {
}

  // Does nothing
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const _hyListOrderable<PAYLOAD>& source):_hyListOrderable<PAYLOAD> (source) {
}

template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric (const unsigned long l, const PAYLOAD start, const PAYLOAD step) {
  this->Initialize(false);
  Populate (l, start, step);
}


template<typename PAYLOAD>
PAYLOAD _hyListNumeric<PAYLOAD>::Element(const long index) const {
  if (this->lLength) {
    if (index >= 0L && (unsigned long)index < this->lLength) {
      return this->lData[index];
    } else if ((unsigned long)(-index) <= this->lLength) {
      return this->lData[this->lLength - (-index)];
    }
  }
  return _HY_LIST_NUMERIC_INVALID_VALUE_;
}

template<typename PAYLOAD>
PAYLOAD _hyListNumeric<PAYLOAD>::Sum(void) const {
  PAYLOAD sum = 0L;
  for (unsigned long k = 0UL; k < this->lLength; k++) {
    sum += this->lData[k];
  }
  return sum;
}

template<typename PAYLOAD>
void _hyListNumeric<PAYLOAD>::Offset(const PAYLOAD shift) {
  for (unsigned long k = 0UL; k < this->lLength; k++) {
    this->lData[k] += shift;
  }
}

template<typename PAYLOAD>
void _hyListNumeric<PAYLOAD>::Populate (const unsigned long l, const PAYLOAD start, const PAYLOAD step) {
  this->RequestSpace (l);
  PAYLOAD current_value = start;
  for (unsigned long k = 0UL; k < l; k++, current_value+=step) {
    this->lData[k] = current_value;
  }
  
  this->lLength = l;
}


template<typename PAYLOAD>
void _hyListNumeric<PAYLOAD>::AppendNumtoStr(_StringBuffer* s, PAYLOAD num) const{
  
  char buffer[128];
  snprintf(buffer, 127, "%ld", long (num));
  (*s) << buffer;
}


  //Char* conversion
template<typename PAYLOAD>
BaseRef _hyListNumeric<PAYLOAD>::toStr(void) const {
  
  if (this->countitems()) {
    _StringBuffer* s = new _StringBuffer();
    
    (*s) << '{';
    
    for (unsigned long i = 0UL; i<this->lLength; i++) {
      this->AppendNumtoStr(s, this->lData[i]);
      if (i<this->lLength-1) {
        (*s) << ',';
      }
    }
    (*s) << '}';
    
    return s;
  } else {
    return new _String ("{}");
  }
}


template<typename PAYLOAD>
_StringBuffer* _hyListNumeric<PAYLOAD>::ListToPartitionString (void) const {
  
  _StringBuffer *result = new _StringBuffer (64UL);
  
  for (unsigned long k=0UL; k<this->lLength; k++) {
    unsigned long m;
    for (m=k+1UL; m<this->lLength; m++)
      if (this->lData[m]-this->lData[m-1]!=1L) {
        break;
      }
    if (m>k+2) {
      this->AppendNumtoStr(result, this->lData[k]);
      (*result) << '-';
      this->AppendNumtoStr(result, this->lData[m-1]);
      if (m<this->lLength) {
        (*result) << ',';
      }
      k = m-1;
    } else {
      this->AppendNumtoStr(result, this->lData[k]);
      if (k<this->lLength-1) {
        (*result) << ',';
      }
    }
  }
  return result;
}

template<typename PAYLOAD>
_hyListNumeric <PAYLOAD>  _hyListNumeric<PAYLOAD>::CountingSort (PAYLOAD upperBound, _hyListNumeric <long> * ordering)
{
  if (ordering) {
    ordering->Clear();
  }
  
  _hyListNumeric<PAYLOAD> result;

  if (this->lLength) {
    if (upperBound == _HY_LIST_NUMERIC_INVALID_VALUE_) {
      upperBound = this->Max()+1UL;
    }
    
    _hyListNumeric<PAYLOAD> count;
    
    count.RequestSpace (upperBound+1UL);
    result.RequestSpace (this->lLength, true);
    
    
    count.Populate(upperBound, 0UL, 0UL);
    
    for (unsigned long pass1 = 0UL; pass1 < this->lLength; pass1++) {
      count[this->lData[pass1]]++;
    }
    
    for (unsigned long pass2 = 1UL; pass2 < upperBound; pass2++) {
      count[pass2] += count[pass2-1];
    }
    
    if (ordering) {
      ordering->Populate (this->lLength, 0UL, 0UL);
      for (long pass3 = this->lLength-1; pass3 >=0L; pass3--) {
        result[--count[this->Element(pass3)]] = this->Element(pass3);
        (*ordering)[count.Element(this->Element(pass3))] = pass3;
      }
    } else {
      for (long pass3 = this->lLength-1; pass3 >= 0L; pass3--) {
        result[--count[this->Element(pass3)]] = this->Element(pass3);
      }
    }
  }
  return result;
}

#endif
