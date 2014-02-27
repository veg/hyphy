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

/*
 ==============================================================
 Constructors
 ==============================================================
 */ // Does nothing
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric() {
}

  //Data constructor (1 member list)
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const PAYLOAD item) : _hyListOrderable<PAYLOAD> (item) {
}


  //Length constructor
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(unsigned long l) : _hyListOrderable<PAYLOAD> (l) {
  this->lData = (PAYLOAD *)MemAllocate(l * sizeof(PAYLOAD));
}

  //Stack copy contructor
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const _hyListNumeric <PAYLOAD> &l, const long from, const long to)
{
  this->Clone(&l, from, to); }

  // Data constructor (variable number of long constants)
template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric(const PAYLOAD value1, const unsigned long number, ...)
{
  this->Initialize(true);
  va_list vl;
  
  this->append(value1);
  
  va_start(vl, number);
  for (unsigned long arg_id = 0; arg_id < number; arg_id++) {
    const PAYLOAD this_arg = (PAYLOAD)va_arg(vl, PAYLOAD);
    this->append(this_arg);
  }
  va_end(vl);
}


template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric (const unsigned long l, const PAYLOAD start, const PAYLOAD step) {
  this->Initialize(false);
  Populate (l, start, step);
}


template<typename PAYLOAD>
PAYLOAD _hyListNumeric<PAYLOAD>::Element(const long index) {
  if (this->lLength) {
    if (index >= 0L && index < this->lLength) {
      return this->lData[index];
    } else if (-index <= this->lLength) {
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

    _StringBuffer* char_list = new _StringBuffer();

    if(num != 0UL) {
      for(PAYLOAD val = num; val; val /= 10) {
          (*char_list) << "0123456789"[val % 10];
      }
    } else {
      (*char_list) << "0";
    }

    for(int j=0; j < char_list->sLength; j++) {
        (*s) << char_list->sData[char_list->sLength - j - 1];
    }

    return;
}

//Char* conversion
template<typename PAYLOAD>
BaseRef _hyListNumeric<PAYLOAD>::toStr(void) {

  if (this->lLength) {
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
_String* _hyListNumeric<PAYLOAD>::ListToPartitionString (void) const {

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
_hyListNumeric <PAYLOAD>*  _hyListNumeric<PAYLOAD>::CountingSort (PAYLOAD upperBound, _hyListNumeric <long> * ordering)
{
    if (ordering) {
        ordering->Clear();
    }

    if (this->lLength) {
        if (upperBound == _HY_LIST_NUMERIC_INVALID_VALUE_) {
            upperBound = this->Max()+1UL;
        } else {
          upperBound = upperBound;
        }

        _hyListNumeric<PAYLOAD> *count  =  new _hyListNumeric <PAYLOAD> ((long)upperBound);
        _hyListNumeric<PAYLOAD> *result =  new _hyListNumeric <PAYLOAD> (this->lLength + 1);
                  
        count->Populate(upperBound, 0UL, 0UL);
        
        for (unsigned long pass1 = 0UL; pass1 < this->lLength; pass1++) {
            count->lData[this->lData[pass1]]++;
        }

        for (unsigned long pass2 = 1UL; pass2 < upperBound; pass2++) {
            count->lData[pass2] += count->lData[pass2-1];
        }

        if (ordering) {
            ordering->Populate (this->lLength, 0UL, 0UL);
            for (long pass3 = this->lLength-1; pass3 >=0L; pass3--) {
                result->lData[--count->lData[this->lData[pass3]]] = this->lData[pass3];
                ordering->lData[count->lData[this->lData[pass3]]] = pass3;
            }
        } else {
            for (long pass3 = this->lLength-1; pass3 >= 0L; pass3--) {
                result->lData[--count->lData[this->lData[pass3]]] = this->lData[pass3];
            }
        }
        result->lLength = this->lLength;
        return result;
    }
    return new _hyListNumeric <PAYLOAD>;
}

