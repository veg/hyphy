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

#ifndef     __FSTRING__
#define     __FSTRING__

#include "mathobj.h"
#include "hy_string_buffer.h"
#include "_hyExecutionContext.h"

//__________________________________________________________________________________

class _FString : public _MathObject   // strings encountered in formulas
{

public:

    _FString (_String const&, bool = true);
    _FString (long);
    _FString (_String*);
    _FString (const char*);
    _FString (void);
    _FString (const _FString&);
    _FString (_FString&&);
    virtual  ~_FString ();
//  ~_Constant (void);

    virtual BaseRef   makeDynamic       (void) const;
    virtual void      Duplicate         (BaseRefConst);
    virtual HBLObjectRef Add               (HBLObjectRef, HBLObjectRef cache);
    virtual long      AddOn             (HBLObjectRef);
    virtual HBLObjectRef AreEqual          (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef AreEqualCIS       (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef Less              (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef LessEq            (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef Greater           (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef GreaterEq         (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef NotEqual          (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef RerootTree        (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef EqualAmb          (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef EqualRegExp       (HBLObjectRef,bool, HBLObjectRef cache);
    virtual HBLObjectRef ReplaceReqExp     (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef CountGlobalObjects(HBLObjectRef cache);
    virtual HBLObjectRef FileExists        (HBLObjectRef cache);
    virtual HBLObjectRef Call              (_List*,_hyExecutionContext*,HBLObjectRef cache);
    virtual HBLObjectRef Sum               (HBLObjectRef cache);
    virtual HBLObjectRef Evaluate          (_hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual HBLObjectRef SubstituteAndSimplify
                                        (HBLObjectRef arguments,HBLObjectRef cache);
    virtual HBLObjectRef Join              (HBLObjectRef,HBLObjectRef cache);
    virtual HBLObjectRef Differentiate     (HBLObjectRef,HBLObjectRef cache);
    virtual unsigned long      ObjectClass       (void) const {
        return STRING;
    }
    
    
    /**
        Bind string contents directly to the string argument
        *this takes over reference management for the argument (which is assumed to be a "free-and-clear" object
     
     */
    
    void  SetStringContent (_StringBuffer * );
    template <class T> void  SetData          (T src) {
        if (the_string) {
            the_string->Clear();
            (*the_string) << src;
        } else {
            SetStringContent (new _StringBuffer (src));
        }
    }
    
    
    virtual HBLObjectRef Compute           (void) {
        return this;
    }
    HBLObjectRef         Dereference       (bool ignore_context, _hyExecutionContext* context = _hyDefaultExecutionContext, bool return_variable_ref = false);

    virtual HBLObjectRef MapStringToVector (HBLObjectRef, HBLObjectRef cache);
    virtual HBLObjectRef CharAccess        (HBLObjectRef,HBLObjectRef,HBLObjectRef cache);
    virtual HBLObjectRef ExecuteSingleOp   (long opCode, _List* args = nil, _hyExecutionContext* context = _hyDefaultExecutionContext, HBLObjectRef cache = nil);
    virtual BaseRef   toStr             (unsigned long = 0UL);

    virtual bool      IsVariable        (void) {
        return true;
    }

    virtual bool      HasChanged        (bool = false) {
        return true;
    }
  
    bool              has_data          (void) const {return the_string && the_string->nonempty();}
  
    hyComparisonType  Compare           (HBLObjectRef, bool convert_non_strings = true);
  
    inline _StringBuffer const&    get_str           (void) const {return *the_string;}
    _StringBuffer *    get_str_ref        (void)  {the_string->AddAReference(); return the_string;}

    virtual bool      empty           (void) const {
        return !the_string || the_string->empty();
    }
    // SLKP 20100907: a simple utility function to check if the object is an empty string

protected:
  
    _StringBuffer*          the_string;
  
};

template <class T> HBLObjectRef _returnStringOrUseCache (T source, HBLObjectRef cache) {
    if (cache && cache->ObjectClass() == STRING) {
        ((_FString*)cache)->SetData(source);
    }
    return new _FString (source);
}

#endif
