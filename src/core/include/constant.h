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

#ifndef     __CONSTANT__
#define     __CONSTANT__

#include "mathobj.h"
#include "global_things.h"

#define  _HY_CONSTANT_PREALLOCATE_SLOTS 16384


class _Constant : public _MathObject {

private:
    template <class T> HBLObjectRef _check_type_and_compute (HBLObjectRef operand, T functor) {
        if (operand) {
            if (operand->ObjectClass() == NUMBER) {
                return new _Constant (functor (Value (), ((_Constant*)operand)->Value()));
            } else {
                hy_global::HandleApplicationError (_String("<'constant' operation 'X'>, where 'X' is not a number. \nconstant = ") & (_String((_String*)toStr())) & "\n'X' = " & (_String((_String*)operand->toStr())));
            }
        } else {
            hy_global::HandleApplicationError (_String("<'constant' operation 'null'>, where constant = ") & (_String((_String*)toStr())));
        }
        return new _MathObject;
    }
    
    template <class T> HBLObjectRef _check_type_and_compute_3 (HBLObjectRef operand, HBLObjectRef operand2, T functor) {
        if (operand && operand2 && operand->ObjectClass() == NUMBER && operand2->ObjectClass() == NUMBER) {
            return new _Constant (functor (Value (), ((_Constant*)operand)->Value(), ((_Constant*)operand2)->Value()));
        }
        hy_global::HandleApplicationError ("Not a numeric 'X' type in a <'constant' operation 'X'> call");
        return new _MathObject;
    }

public:

    _Constant (hyFloat);
    _Constant (_String&);
    _Constant (void);
    ~_Constant (void) {}

    virtual HBLObjectRef Add           (HBLObjectRef);
    virtual HBLObjectRef Sub           (HBLObjectRef);
    virtual HBLObjectRef Minus         (void) ;
    virtual HBLObjectRef Sum           (void) ;
    virtual HBLObjectRef Mult          (HBLObjectRef);
    virtual HBLObjectRef Div           (HBLObjectRef);
    virtual HBLObjectRef lDiv          (HBLObjectRef);
    virtual HBLObjectRef longDiv       (HBLObjectRef);
    virtual HBLObjectRef Raise         (HBLObjectRef);
    virtual bool         Equal         (HBLObjectRef);
    virtual HBLObjectRef Abs           (void);
    virtual HBLObjectRef Sin           (void);
    virtual HBLObjectRef Cos           (void);
    virtual HBLObjectRef Tan           (void);
    virtual HBLObjectRef Exp           (void);
    virtual HBLObjectRef Log           (void);
    virtual HBLObjectRef Sqrt          (void);
    virtual HBLObjectRef Time          (void);
    virtual HBLObjectRef Arctan        (void);
    virtual HBLObjectRef Gamma         (void);
    virtual HBLObjectRef LnGamma       (void);         /* <- added by afyp, February 8, 2007 */
    virtual HBLObjectRef Beta          (HBLObjectRef);
    virtual HBLObjectRef Min           (HBLObjectRef);
    virtual HBLObjectRef Max           (HBLObjectRef);
    virtual HBLObjectRef GammaDist     (HBLObjectRef,HBLObjectRef);
    virtual HBLObjectRef CGammaDist    (HBLObjectRef,HBLObjectRef);
    virtual HBLObjectRef IBeta         (HBLObjectRef,HBLObjectRef);
    virtual HBLObjectRef IGamma        (HBLObjectRef);
    virtual HBLObjectRef CChi2         (HBLObjectRef);
    virtual HBLObjectRef InvChi2       (HBLObjectRef);
    virtual HBLObjectRef Erf           (void);
    virtual HBLObjectRef ZCDF          (void);
    virtual HBLObjectRef Less          (HBLObjectRef);
    virtual HBLObjectRef Greater       (HBLObjectRef);
    virtual HBLObjectRef LessEq        (HBLObjectRef);
    virtual HBLObjectRef GreaterEq     (HBLObjectRef);
    virtual HBLObjectRef AreEqual      (HBLObjectRef);
    virtual HBLObjectRef NotEqual      (HBLObjectRef);
    virtual HBLObjectRef LAnd          (HBLObjectRef);
    virtual HBLObjectRef LOr           (HBLObjectRef);
    virtual HBLObjectRef LNot          ();
    virtual HBLObjectRef Random        (HBLObjectRef);
    virtual hyFloat
    Value       (void);
    virtual HBLObjectRef FormatNumberString
    (HBLObjectRef,HBLObjectRef);
    virtual HBLObjectRef Compute       (void) {
        return this;
    };

    virtual   void    Initialize            (bool = true);
    virtual   void    Duplicate             (BaseRefConst);
    virtual   BaseRef makeDynamic           (void) const;
    virtual   BaseRef toStr                 (unsigned long = 0UL);
    
    virtual   unsigned long    ObjectClass  (void) const {
        return NUMBER;
    }
    
    virtual   void    SetValue              (hyFloat pl) {
        theValue = pl;
    }
    
    void * operator new       (size_t size);
    void   operator delete    (void * p);
  
    static  _SimpleList                    free_slots;
    static unsigned char                       preallocated_buffer[];

public:
    hyFloat theValue;
    
    
};

hyFloat _ln_gamma (hyFloat alpha);
hyFloat  gaussDeviate (void);
hyFloat  exponDeviate (void);
hyFloat  gammaDeviate (hyFloat a, hyFloat scale = 1.);
hyFloat  chisqDeviate (double df);

#endif
