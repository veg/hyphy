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

#ifndef     __VARIABLE__
#define     __VARIABLE__

#include "constant.h"
#include "hy_strings.h"
#include "list.h"
#include "avllist.h"
#include "operation.h"
#include "formula.h"


class _Variable : public _Constant {

    friend class _Operation;

public:

    _Variable (void);
    _Variable (_String const&, bool isG = false); // name
    _Variable (_String const&, _String const&, bool isG = false); // name and formula

    virtual ~_Variable (void);

    virtual   void          Initialize (bool = true);
    virtual   void          Duplicate (BaseRefConst);
    virtual   BaseRef       makeDynamic(void) const;
    virtual   BaseRef       toStr (unsigned long = 0UL);
    virtual    void         toFileStr (FILE*, unsigned long = 0UL);

    virtual   void          MarkDone (void);

    virtual     HBLObjectRef   Compute (void);       // compute or return the value
                HBLObjectRef   ComputeMatchingType (long);
                // return a value if the type is matched, otherwise nil
    virtual     bool        IsVariable (void); //
    virtual     bool        IsIndependent (void) {
        return (varFormula&&varFormula->theFormula.lLength)?
               false:
               (varValue?varValue->IsIndependent():true);
    }
    virtual     bool        IsConstant (void);
    //void        SetValue (HBLObjectRef, bool = true, bool = true); // set the value of the variable
    void        SetValue (HBLObjectRef, bool duplicate_value, bool do_checks, _AVLList* keep_track_of_changes);
    void        SetValue (hyFloat); // set the value of the variable
    void        SetNumericValue (hyFloat);
    void        CheckAndSet (hyFloat, bool, _AVLList* keep_track_of_changed);
    // set the value of the variable
    // bool flag is used to indicate that out of bounds values should be rejected

    HBLObjectRef   GetValue (void) {
        return varValue;   // get the value of the variable
    }
    void        SetFormula (_Formula&); // set the variable to a new formula
  
    const     _Formula * get_constraint (void) const {
      return varFormula;
    }

    
    void *      operator new       (size_t size);
    void        operator delete    (void * p);

    virtual     bool        HasChanged      (bool = false);
    virtual     void        PreMarkChanged  ();
    virtual     void        PostMarkChanged ();
    virtual     bool        IsGlobal (void) {
        return varFlags & HY_VARIABLE_GLOBAL;
    }
    virtual     bool        IsCategory (void) {
        return false;
    }

    virtual unsigned long        ObjectClass (void) const;
    
    void        SetIndex (long i) {
        theIndex = i;
    }
    long        get_index (void) const {
        return theIndex;
    }
    virtual void        ScanForVariables (_AVLList& l, bool globals = false, _AVLListX* tagger = nil, long weight = 0) const;
    virtual     bool        IsContainer (void) {
        return false;
    }

    void        SetBounds (hyFloat lb, hyFloat ub);
    void        EnsureTheValueIsInBounds (void);
    bool        IsValueInBounds (hyFloat v)
                           { return v >= lowerBound && v <= upperBound; }

    hyFloat  GetLowerBound (void) {
        return lowerBound;
    }
    hyFloat  GetUpperBound (void) {
        return upperBound;
    }

    virtual     void        ClearConstraints    (void);
    virtual     bool        CheckFForDependence (long, bool = false);
    virtual     bool        CheckFForDependence (_AVLList const&, bool = false);
    virtual     bool        HasBeenInitialized (void) const {return !(varFlags & HY_VARIABLE_NOTSET);}
    virtual     void        MarkModified  (void) {varFlags = varFlags | HY_VARIABLE_CHANGED;}

    _String const     ContextFreeName                 (void) const;
    _StringBuffer&    ContextFreeName                 (_StringBuffer&) const;
    _String const    ParentObjectName                 (void) const;
 
    _String*    GetName                         (void) const{
        return theName;
    }
    _String*    GetFormulaString        (_hyFormulaStringConversionMode mode) {
        return varFormula?(_String*)varFormula->toStr(mode):new _String;
    }

    virtual     void        CompileListOfDependents (_SimpleList&);
    HBLObjectRef   ComputeReference        (_MathObject const *) const;


    friend      void        ResetVariables          (void);
    //friend      _Variable*  LocateVar               (long);
    friend      void        InsertVar               (_Variable*);
    bool        has_been_set                        (void) const {return !(HY_VARIABLE_NOTSET & varFlags);}

public:

    _String*   theName;

    HBLObjectRef  varValue;

    long       theIndex; // index of this variable in the global variable pool

    // the class of this variable - i.e global, local, category or random
    int       varFlags;

    hyFloat lowerBound,
               upperBound;
    // dynamic lower and upper bounds here

    _Formula*  varFormula;

};

long    DereferenceVariable (long index, _MathObject const *  context, char reference_type);
long    DereferenceString   (HBLObjectRef, _MathObject const * context, char reference_type);
_String const WrapInNamespace (_String const&, _String const*);


#endif
