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

template<typename PAYLOAD>
_hyListNumeric<PAYLOAD>::_hyListNumeric (const unsigned long l, const PAYLOAD start, const PAYLOAD step) {
  this->Initialize(false);
  Populate (l, start, step);
}


template<typename PAYLOAD>
PAYLOAD _hyListNumeric<PAYLOAD>::Element(const long index)
{
  if (index >= 0L && index < this->lLength) {
    return this->lData[index];
  } else if (-index <= this->lLength) {
    return this->lData[this->lLength - (-index)];
  }
  return _HY_LIST_NUMERIC_INVALID_VALUE_;
}

template<typename PAYLOAD>
PAYLOAD _hyListNumeric<PAYLOAD>::Sum(void) const
{
    PAYLOAD sum = 0L;
    for (unsigned long k = 0UL; k < this->lLength; k++) {
       sum += this->lData[k];
    }
    return this->sum;
}

template<typename PAYLOAD>
void _hyListNumeric<PAYLOAD>::Offset(const PAYLOAD shift)
{
    for (unsigned long k = 0UL; k < this->lLength; k++) {
       this->lData[k] += shift;
    }
}

template<typename PAYLOAD>
void _hyListNumeric<PAYLOAD>::Populate (const unsigned long l, const PAYLOAD start, const PAYLOAD step)
{
    this->RequestSpace (l);
    PAYLOAD current_value = start;
    for (unsigned long k = 0UL; k < l; k++, current_value+=step) {
        this->lData[k] = current_value;
    }

    this->lLength = l;
}


//Char* conversion
template<typename PAYLOAD>
BaseRef _hyListNumeric<PAYLOAD>::toStr(void)
{
  if (this->lLength) {
      _String * s = new _String (10L, true);
      (*s) << '{';

      for (unsigned long i = 0UL; i<this->lLength; i++) {
          (*s) << this->lData[i];
          if (i<this->lLength-1) {
              (*s) << ',';
          }
      }
      (*s) << '}';

      s->Finalize();
      return s;
  } else {
      return new _String ("{}");
  }
}


template<typename PAYLOAD>
_String* _hyListNumeric<PAYLOAD>::ListToPartitionString (void) const
{
    _String *result = new _String ((unsigned long)64,true),
    conv;

    for (unsigned long k=0UL; k<this->lLength; k++) {
        unsigned long m;
        for (m=k+1UL; m<this->lLength; m++)
            if (this->lData[m]-this->lData[m-1]!=1L) {
                break;
            }
        if (m>k+2) {
            conv = this->lData[k];
            (*result) << & conv;
            (*result) << '-';
            conv = this->lData[m-1];
            (*result) << & conv;
            if (m<this->lLength) {
                (*result) << ',';
            }
            k = m-1;
        } else {
            conv = this->lData[k];
            (*result) << &conv;
            if (k<this->lLength-1) {
                (*result) << ',';
            }
        }
    }
    (*result).Finalize();
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
        }

        _hyListNumeric<PAYLOAD> buffer,
                    * result =  new _hyListNumeric <PAYLOAD> (this->lLength);
                  
        buffer.Populate (upperBound, 0L, 0L);
        
        for (unsigned long pass1 = 0UL; pass1 < this->lLength; pass1 ++) {
            buffer.lData[this->lData[pass1]] ++;
        }
        for (unsigned long pass2 = 1UL; pass2 < upperBound; pass2 ++) {
            buffer.lData[pass2] += buffer.lData[pass2-1];
        }
        
        if (ordering) {
            ordering->Populate (this->lLength, 0L, 0L);
            for (long pass3 = this->lLength-1; pass3 >=0L; pass3--) {
                result->lData[--buffer.lData[this->lData[pass3]]] = this->lData[pass3];
                ordering->lData[buffer.lData[this->lData[pass3]]] = pass3;
            }
        } else
            for (long pass3 = this->lLength-1; pass3 >=0L; pass3--) {
                result->lData[--buffer.lData[this->lData[pass3]]] = this->lData[pass3];
            }
        result->lLength = this->lLength;

        return result;
    }
    return new _hyListNumeric <PAYLOAD>;
}