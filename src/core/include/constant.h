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

#ifndef     __CONSTANT__
#define     __CONSTANT__

#include "mathobj.h"

class _Constant : public _MathObject   // a numerical constant
{

public:

    _Constant (_Parameter);
    _Constant (_String&);
    _Constant (void);
    ~_Constant (void) {}

    virtual _PMathObj Add           (_PMathObj);
    virtual _PMathObj Sub           (_PMathObj);
    virtual _PMathObj Minus         (void) ;
    virtual _PMathObj Sum           (void) ;
    virtual _PMathObj Mult          (_PMathObj);
    virtual _PMathObj Div           (_PMathObj);
    virtual _PMathObj lDiv          (_PMathObj);
    virtual _PMathObj longDiv       (_PMathObj);
    virtual _PMathObj Raise         (_PMathObj);
    virtual void      Assign        (_PMathObj);
    virtual bool      Equal         (_PMathObj);
    virtual _PMathObj Abs           (void);
    virtual _PMathObj Sin           (void);
    virtual _PMathObj Cos           (void);
    virtual _PMathObj Tan           (void);
    virtual _PMathObj Exp           (void);
    virtual _PMathObj Log           (void);
    virtual _PMathObj Sqrt          (void);
    virtual _PMathObj Time          (void);
    virtual _PMathObj Arctan        (void);
    virtual _PMathObj Gamma         (void);
    virtual _PMathObj LnGamma       (void);         /* <- added by afyp, February 8, 2007 */
    virtual _PMathObj Beta          (_PMathObj);
    virtual _PMathObj Min           (_PMathObj);
    virtual _PMathObj Max           (_PMathObj);
    virtual _PMathObj GammaDist     (_PMathObj,_PMathObj);
    virtual _PMathObj CGammaDist    (_PMathObj,_PMathObj);
    virtual _PMathObj IBeta         (_PMathObj,_PMathObj);
    virtual _PMathObj IGamma        (_PMathObj);
    virtual _PMathObj CChi2         (_PMathObj);
    virtual _PMathObj InvChi2       (_PMathObj);
    virtual _PMathObj Erf           (void);
    virtual _PMathObj ZCDF          (void);
    virtual _PMathObj Less          (_PMathObj);
    virtual _PMathObj Greater       (_PMathObj);
    virtual _PMathObj LessEq        (_PMathObj);
    virtual _PMathObj GreaterEq     (_PMathObj);
    virtual _PMathObj AreEqual      (_PMathObj);
    virtual _PMathObj NotEqual      (_PMathObj);
    virtual _PMathObj LAnd          (_PMathObj);
    virtual _PMathObj LOr           (_PMathObj);
    virtual _PMathObj LNot          ();
    virtual _PMathObj Random        (_PMathObj);
    virtual _Parameter
    Value       (void);
    virtual _PMathObj FormatNumberString
    (_PMathObj,_PMathObj);
    virtual _PMathObj Compute       (void) {
        return this;
    };

    virtual   void    Initialize            (void);
    virtual   void    Duplicate             (BaseRef);
    virtual   BaseRef makeDynamic           (void);
    virtual   BaseRef toStr                 (void);
    virtual   unsigned long    ObjectClass           (void) {
        return NUMBER;
    }
    virtual   void    SetValue              (_Parameter pl) {
        theValue = pl;
    }

public:
    _Parameter theValue;
};

#endif
