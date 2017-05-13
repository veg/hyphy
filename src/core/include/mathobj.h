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

#ifndef     __MATHOBJ__
#define     __MATHOBJ__

#include "baseobj.h"
#include "defines.h"
#include "list.h"
#include "avllistx.h"
#include "hy_strings.h"
#include "_hyExecutionContext.h"



class   _MathObject : public BaseObj  { //abstract math operations class
  
protected:
  _MathObject* _extract_argument (_List * arguments, unsigned long index, bool fill_in) const;


public:

    virtual _MathObject* Add        (_MathObject*)     ;
    virtual _MathObject* Sub        (_MathObject*)     ;
    virtual _MathObject* Minus      (void)             ;
    virtual _MathObject* Sum        (void)             ;
    virtual _MathObject* Mult       (_MathObject*)     ;
    virtual _MathObject* Div        (_MathObject*)     ;
    virtual _MathObject* lDiv       (_MathObject*)     ;
    virtual _MathObject* longDiv    (_MathObject*)     ;
    virtual _MathObject* Raise      (_MathObject*)     ;
    virtual void         Assign     (_MathObject*)     ;
    virtual bool         Equal      (_MathObject*)     ;
    virtual _MathObject* Abs        (void)             ;
    virtual _MathObject* Sin        (void)             ;
    virtual _MathObject* Cos        (void)             ;
    virtual _MathObject* Tan        (void)             ;
    virtual _MathObject* Exp        (void)             ;
    virtual _MathObject* Log        (void)             ;
    virtual _MathObject* Sqrt       (void)             ;
    virtual _MathObject* Gamma      (void)             ;
    virtual _MathObject* Erf        (void)             ;
    virtual _MathObject* LnGamma    (void)             ;
    virtual _MathObject* Beta       (_MathObject*)     ;
    virtual _MathObject* IGamma     (_MathObject*)     ;
    virtual _MathObject* CChi2      (_MathObject*)     ;
    virtual _MathObject* IBeta      (_MathObject*,_MathObject*) ;
    virtual _MathObject* Simplex    (void)             ;
    
    virtual _MathObject* Simplify    (void)             ;
    
    virtual _MathObject* Min        (_MathObject*)     ;
    virtual _MathObject* Max        (_MathObject*)     ;
    virtual _MathObject* InvChi2    (_MathObject*)     ;
    virtual _MathObject* ZCDF       (void)             ;
    virtual _MathObject* Time       (void)             ;
    virtual _MathObject* Arctan     (void)             ;
    virtual _MathObject* Less       (_MathObject*)     ;
    virtual _MathObject* Random     (_MathObject*)     ;
    virtual _MathObject* Greater    (_MathObject*)     ;
    virtual _MathObject* LessEq     (_MathObject*)     ;
    virtual _MathObject* GreaterEq  (_MathObject*)     ;
    virtual _MathObject* AreEqual   (_MathObject*)     ;
    virtual _MathObject* NotEqual   (_MathObject*)     ;
    virtual _MathObject* LAnd       (_MathObject*)     ;
    virtual _MathObject* LOr        (_MathObject*)     ;
    virtual _MathObject* GammaDist  (_MathObject*,_MathObject*) ;
    virtual _MathObject* CGammaDist (_MathObject*,_MathObject*) ;
    virtual _MathObject* LNot       (void)             ;
    virtual _MathObject* TipCount   (void)             ;
    virtual _MathObject* BranchCount (void)            ;
    virtual _MathObject* TipName     (_MathObject*)    ;
    virtual _MathObject* BranchName  (_MathObject*)    ;
    virtual _MathObject* BranchLength(_MathObject*)    ;
    virtual _MathObject* RerootTree  (_MathObject*)    ;
    virtual _MathObject* TEXTreeString(_MathObject*) const ;
    virtual _MathObject* Type                          (void);
    virtual _MathObject* PlainTreeString(_MathObject*,_MathObject*) ;
    virtual _MathObject* FormatNumberString (_MathObject*,_MathObject*) ;
    virtual hy_float   Value (void)              ;
    
    virtual _MathObject* Compute (void)            {
        return this;
    }
    virtual void         ScanForVariables (_AVLList&,bool = false, _AVLListX* = nil, long = 0)
    {}

    virtual      BaseRef makeDynamic               (void) const;
    virtual      void    Duplicate                 (BaseRefConst);
    
    virtual bool         IsVariable (void)         {
        return false;
    }
    virtual bool         IsObjectEmpty (void)      {
        return true;
    }
    virtual bool         IsPrintable (void)        {
        return false;
    }

    virtual bool         IsIndependent (void)       {
        return true;
    }
    virtual unsigned long  ObjectClass (void)       {
        return HY_UNDEFINED;
    }
    // returns a unique ID for this object
    // 0 - undefined
    // 1 - number
    // 4 - matrix

  
    virtual _MathObject* ExecuteSingleOp (long opCode, _List* arguments = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    // execute this operation with the list of Args

    virtual bool         HasChanged (bool = false) {
        return false;
    }

    virtual   bool       IsConstant (void) {
        return true;
    }
    
    
};

// pointer to a math object
typedef _MathObject* _PMathObj ;


#endif
