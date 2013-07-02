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

#ifndef     __FSTRING__
#define     __FSTRING__

#include "mathobj.h"
#include "_hyExecutionContext.h"

//__________________________________________________________________________________

class _FString : public _MathObject   // strings encountered in formulas
{

public:

    _FString (_String&, bool = true);
    _FString (long);
    _FString (_String*);
    _FString (void);
    virtual  ~_FString ();
//  ~_Constant (void);

    virtual BaseRef   makeDynamic       (void);
    virtual void      Duplicate         (BaseRef);
    virtual _PMathObj Add               (_PMathObj);
    virtual long      AddOn             (_PMathObj);
    virtual _PMathObj AreEqual          (_PMathObj);
    virtual _PMathObj AreEqualCIS       (_PMathObj);
    virtual _PMathObj Less              (_PMathObj);
    virtual _PMathObj LessEq            (_PMathObj);
    virtual _PMathObj Greater           (_PMathObj);
    virtual _PMathObj GreaterEq         (_PMathObj);
    virtual _PMathObj NotEqual          (_PMathObj);
    virtual _PMathObj RerootTree        (void);
    virtual _PMathObj EqualAmb          (_PMathObj);
    virtual _PMathObj EqualRegExp       (_PMathObj,bool = false);
    virtual _PMathObj ReplaceReqExp     (_PMathObj);
    virtual _PMathObj CountGlobalObjects(void);
    virtual _PMathObj FileExists        (void);
    virtual _PMathObj Evaluate          (_hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual _PMathObj Join              (_PMathObj);
    virtual _PMathObj Differentiate     (_PMathObj);
    virtual unsigned long      ObjectClass       (void) {
        return STRING;
    }
    virtual _PMathObj Compute           (void) {
        return this;
    }
    _PMathObj         Dereference       (bool ignore_context, _hyExecutionContext* context = _hyDefaultExecutionContext, bool return_variable_ref = false);

    virtual _PMathObj MapStringToVector (_PMathObj);
    virtual _PMathObj CharAccess        (_PMathObj,_PMathObj);
    virtual _PMathObj Execute           (long opCode, _MathObject* p = nil , _MathObject* p2 = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual BaseRef   toStr             (void);

    virtual bool      IsVariable        (void) {
        return true;
    }

    virtual bool      HasChanged        (void) {
        return true;
    }

    virtual bool      IsEmpty           (void) {
        return !theString || theString->sLength == 0;
    }
    // SLKP 20100907: a simple utility function to check if the object is an empty string

    _String*          theString;

};

#endif
