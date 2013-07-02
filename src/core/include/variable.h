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

#ifndef     __VARIABLE__
#define     __VARIABLE__

#include "constant.h"
#include "hy_strings.h"
#include "list.h"
#include "avllist.h"
#include "operation.h"
#include "formula.h"


class _Variable : public _Constant
{

    friend class _Operation;

public:

    _Variable (void);
    _Variable (_String&, bool isG = false); // name
    _Variable (_String&, _String&, bool isG = false); // name and formula

    virtual ~_Variable (void);

    virtual   void          Initialize (void);
    virtual   void          Duplicate (BaseRef);
    virtual   BaseRef       makeDynamic(void);
    virtual   BaseRef       toStr (void);
    virtual    void         toFileStr (FILE*);

    virtual   void          MarkDone (void);

    virtual     _PMathObj   Compute (void);       // compute or return the value
    virtual     bool        IsVariable (void); //
    virtual     bool        IsIndependent (void) {
        return (varFormula&&varFormula->theFormula.lLength)?
               false:
               (varValue?varValue->IsIndependent():true);
    }
    virtual     bool        IsConstant (void);
    void        SetValue (_PMathObj, bool = true); // set the value of the variable
    void        SetNumericValue (_Parameter);
    void        CheckAndSet (_Parameter, bool = false);
    // set the value of the variable
    // bool flag is used to indicate that out of bounds values should be rejected

    _PMathObj   GetValue (void) {
        return varValue;   // get the value of the variable
    }
    void        SetFormula (_Formula&); // set the variable to a new formula

    virtual     bool        HasChanged      (bool = false);
    virtual     void        PreMarkChanged  ();
    virtual     void        PostMarkChanged ();
    virtual     bool        IsGlobal (void) {
        return varFlags & HY_VARIABLE_GLOBAL;
    }
    virtual     bool        IsCategory (void) {
        return false;
    }
    virtual     long        GetAVariable (void) {
        return theIndex;
    }
    virtual unsigned long        ObjectClass (void) {
        return varValue?varValue->ObjectClass():((varFormula&&varFormula->theFormula.lLength)?varFormula->ObjectClass():1);
    }
    void        SetIndex (long i) {
        theIndex = i;
    }
    long        GetIndex (void) {
        return theIndex;
    }
    virtual     void        ScanForVariables (_AVLList& l, bool globals = false, _AVLListX* tagger = nil, long weight = 0);
    virtual     bool        IsContainer (void) {
        return false;
    }

    void        SetBounds (_Parameter lb, _Parameter ub);
    void        EnsureTheValueIsInBounds (void);
    bool        IsValueInBounds (_Parameter v)
                           { return v >= lowerBound && v <= upperBound; }

    _Parameter  GetLowerBound (void) {
        return lowerBound;
    }
    _Parameter  GetUpperBound (void) {
        return upperBound;
    }

    virtual     void        ClearConstraints    (void);
    virtual     bool        CheckFForDependence (long, bool = false);

    _String     ContextFreeName                 (void);
    _String     ParentObjectName                 (void);
 
    _String*    GetName                         (void) {
        return theName;
    }
    _String*    GetFormulaString        (void) {
        return varFormula?(_String*)varFormula->toStr():(_String*)empty.makeDynamic();
    }

    virtual     void        CompileListOfDependents (_SimpleList&);
    virtual     _PMathObj   ComputeReference        (_PMathObj);


    friend      void        ResetVariables          (void);
    friend      _Variable*  LocateVar               (long);
    friend      void        InsertVar               (_Variable*);

public:

    _String*   theName;

    _PMathObj  varValue;

    long       theIndex; // index of this variable in the global variable pool

    // the class of this variable - i.e global, local, category or random
    char       varFlags;

    _Parameter lowerBound,
               upperBound;
    // dynamic lower and upper bounds here

    _Formula*  varFormula;

};

long    DereferenceVariable (long index, _PMathObj context, char reference_type);
long    DereferenceString   (_PMathObj, _PMathObj context, char reference_type);


#endif
