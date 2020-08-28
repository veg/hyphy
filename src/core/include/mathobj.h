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
#include "avllistx.h"
#include "hy_strings.h"
#include "_hyExecutionContext.h"



class   _MathObject : public BaseObj  { //abstract math operations class
  
protected:
  _MathObject* _extract_argument (_List * arguments, unsigned long index, bool fill_in) const;


public:

    virtual _MathObject* Add        (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* Sub        (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* Minus      (_MathObject* cache = nil)             ;
    virtual _MathObject* Sum        (_MathObject* cache = nil)             ;
    virtual _MathObject* Mult       (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* Div        (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* lDiv       (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* longDiv    (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* Raise      (_MathObject*, _MathObject* cache = nil)     ;
    virtual bool         Equal      (_MathObject*)     ;
    virtual _MathObject* Abs        (_MathObject* cache = nil)             ;
    virtual _MathObject* Sin        (_MathObject* cache = nil)             ;
    virtual _MathObject* Cos        (_MathObject* cache = nil)             ;
    virtual _MathObject* Tan        (_MathObject* cache = nil)             ;
    virtual _MathObject* Exp        (_MathObject* cache = nil)             ;
    virtual _MathObject* Log        (_MathObject* cache = nil)             ;
    virtual _MathObject* Sqrt       (_MathObject* cache = nil)             ;
    virtual _MathObject* Gamma      (_MathObject* cache = nil)             ;
    virtual _MathObject* Erf        (_MathObject* cache = nil)             ;
    virtual _MathObject* LnGamma    (_MathObject* cache = nil)             ;
    virtual _MathObject* Beta       (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* IGamma     (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* CChi2      (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* IBeta      (_MathObject*,_MathObject*, _MathObject* cache = nil) ;
    virtual _MathObject* Simplex    (_MathObject* cache = nil)             ;
    
    virtual _MathObject* Simplify    (_MathObject* cache = nil)             ;
    
    virtual _MathObject* Min        (_MathObject*,_MathObject* cache = nil)     ;
    virtual _MathObject* Max        (_MathObject*,_MathObject* cache = nil)     ;
    virtual _MathObject* InvChi2    (_MathObject*,_MathObject* cache = nil)     ;
    virtual _MathObject* ZCDF       (_MathObject* cache = nil)             ;
    virtual _MathObject* Time       (_MathObject* cache = nil)             ;
    virtual _MathObject* Arctan     (_MathObject* cache = nil)             ;
    virtual _MathObject* Less       (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* Random     (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* Greater    (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* LessEq     (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* GreaterEq  (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* AreEqual   (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* NotEqual   (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* LAnd       (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* LOr        (_MathObject*, _MathObject* cache = nil)     ;
    virtual _MathObject* GammaDist  (_MathObject*,_MathObject*, _MathObject* cache = nil) ;
    virtual _MathObject* CGammaDist (_MathObject*,_MathObject*, _MathObject* cache = nil) ;
    virtual _MathObject* LNot       (_MathObject* cache = nil)             ;
    virtual _MathObject* TipCount   (_MathObject* cache = nil)             ;
    virtual _MathObject* BranchCount (_MathObject* cache = nil)            ;
    virtual _MathObject* TipName     (_MathObject*, _MathObject* cache = nil)    ;
    virtual _MathObject* BranchName  (_MathObject*, _MathObject* cache = nil)    ;
    virtual _MathObject* BranchLength(_MathObject*, _MathObject* cache = nil)    ;
    virtual _MathObject* RerootTree  (_MathObject*, _MathObject* cache = nil)    ;
    virtual _MathObject* TEXTreeString(_MathObject*, _MathObject* cache = nil) const ;
    virtual _MathObject* Type                          (_MathObject* cache = nil);
    virtual _MathObject* PlainTreeString(_MathObject*,_MathObject*, _MathObject* cache = nil) ;
    virtual _MathObject* FormatNumberString (_MathObject*,_MathObject*, _MathObject* cache = nil) ;
    virtual hyFloat   Value (void)              ;
    
    virtual _MathObject* Compute (void)            {
        return this;
    }
    virtual void         ScanForVariables (_AVLList&,bool = false, _AVLListX* = nil, long = 0) const
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
    virtual unsigned long  ObjectClass (void) const       {
        return HY_UNDEFINED;
    }
    // returns a unique ID for this object
    // 0 - undefined
    // 1 - number
    // 4 - matrix

  
    virtual _MathObject* ExecuteSingleOp (long opCode, _List* arguments = nil, _hyExecutionContext* context = _hyDefaultExecutionContext, _MathObject* cache = nil);
    // execute this operation with the list of Args

    virtual bool         HasChanged (bool = false) {
        return false;
    }

    virtual   bool       IsConstant (void) {
        return true;
    }
    
    private :
        _MathObject* _null_handler ();
    
    
    
};

// pointer to a math object

typedef _MathObject* HBLObjectRef ;
typedef _MathObject const * HBLObjectRefConst ;

HBLObjectRef   _returnConstantOrUseCache (hyFloat value, HBLObjectRef cache);

#endif
