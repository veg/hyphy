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

#ifndef     __POLYNOMIALDATA__
#define     __POLYNOMIALDATA__


#define  POLY_DATA_INCREMENT 10

#include "parser.h"

class _PolynomialData : public BaseObj
{

public:

    _PolynomialData (void);
    _PolynomialData (long);
    _PolynomialData (_PolynomialData&);
    _PolynomialData (long,long, _Parameter*);

    virtual ~_PolynomialData ();

    virtual BaseObj*    makeDynamic(void);
    virtual void        Duplicate  (BaseRef);

    inline  _Parameter*         GetCoeff (void) {
        return theCoeff;
    }
    inline  _Parameter&         GetCoeff (long index) {
        return theCoeff[index];
    }

    long    *           GetTerm (long);
    long                GetNoTerms (void) {
        return actTerms;
    }
    void                AddTerm (long*, _Parameter);
    void                AddTerm (long*, _Parameter, long*, long);
    void                AddTerm (_Parameter);
    void                WriteTerm (long*,long);
    void                DeleteTerm (long);
    bool                IsFirstANumber (void);
    inline long         NumberOfTerms (void) {
        return actTerms;
    }
    long                SumOfPowers (long);
    long                WeightedSumOfPowers (long,_Parameter*);

    // temp!

    bool                checkMe (void);

    friend class _Polynomial;

    void                MultiplyTerms (long*, long*, long*);
    void                RaiseTerm     (long*, long);
    static  _Parameter  BinaryRaise   (_Parameter, long);
    static  void        RearrangeTerm (long*, long*, long*,long);
    char                CompareTerms  (long*, long*);
    char                CompareTerms  (long*, long*, long*, long);
    char                CompareTerms  (long*, long*, long*, long*, long, long);
    long                FindTerm      (long*, long*, long start = 0);
    void                ResortTerms   (long*);
    void                ChopTerms     (void);
    bool                checkTerm     (_Parameter, long);


protected:

    _Parameter*     theCoeff;
    long*           thePowers;
    long            numberVars, actTerms, allocTerms;

};

#endif
