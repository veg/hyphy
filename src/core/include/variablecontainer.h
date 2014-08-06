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

#ifndef     __VARIABLECONTAINER__
#define     __VARIABLECONTAINER__

#include "variable.h"
#include "baseobj.h"
#include "list.h"
#include "avllistx.h"
#include "avllistxl.h"
#include "operation.h"

//__________________________________________________________________________________

// this class defines a computational (or storage) class which, as a variable, may contain
// other variables locally.

class _Matrix;
class   _VariableContainer: public _Variable
{

    friend class _Operation;
    friend class _Variable;

public:

    _VariableContainer (void);
    _VariableContainer (_String theName, _String theTmplt = "", _VariableContainer* theP = nil);
    // name, matrix constructor, the parent (if there is one)
    virtual ~_VariableContainer(void);

    void                    InitializeVarCont       (_String&, _String&, _VariableContainer*, _AVLListXL* = nil);
    void                    ScanModelBasedVariables (_String&, _AVLListXL*);
    virtual     void        SetModel                (long,_AVLListXL*);

    virtual     void        MarkDone (void);

    // variable access/operation functions

    virtual     bool        IsContainer                 (void) {
        return true;
    }

    virtual     bool        HasChanged                  (void);
    virtual     bool        NeedToExponentiate          (bool = false);

    void        ScanAndAttachVariables      (void);

    virtual     void        ScanForVariables            (_AVLList&,_AVLList&, _AVLListX* tagger = nil, long weight = 0);
    virtual     void        ScanForDVariables           (_AVLList&,_AVLList&);
    virtual     void        ScanForGVariables           (_AVLList&,_AVLList&, _AVLListX* tagger = nil, long weight = 0);

    virtual     bool        IsModelVar                  (long);
    virtual     bool        IsConstant                  (void);
    virtual     BaseRef     makeDynamic                 (void);
    virtual     void        Duplicate                   (BaseRef);

    virtual     BaseRef     toStr                       (void);

    bool        HasLocals                   (void);

    virtual     bool        RemoveDependance            (long);
    virtual     long        SetDependance               (long);
    bool        SetMDependance              (_SimpleList&);

    void        Clear                       (void);
    virtual     void        ClearConstraints            (void);

    long        CountIndependents           (void);
    long        CountAll                    (void);

    virtual     _Variable*  GetIthIndependent           (long);
    virtual     _Variable*  GetIthDependent             (long);
    virtual     _Variable*  GetIthParameter             (long);

    long        CheckAndAddUserExpression   (_String&, long startWith = 0);
    void        KillUserExpression          (long);
    virtual     void        CompileListOfDependents     (_SimpleList&);

    void        MatchParametersToList       (_List&, bool doAll = false, bool indOnly = false);
    _Matrix*    GetModelMatrix              (_List* = nil, _SimpleList* = nil);
    _Matrix*    GetFreqMatrix               (void);
    bool        HasExplicitFormModel        (void);
    _Formula*   GetExplicitFormModel        (void);

    long        GetModelIndex               (void) {
        return theModel;
    }
    
    _String*    GetModelName                (void);
    
    long        GetModelDimension           (void);
    /* 20100316 SLKP
        return the dimension of the model; needed to handle the case
        of explicit model exponentials
     */

    void        CopyMatrixParameters                (_VariableContainer*, bool match_by_name = false);
    void        GetListOfModelParameters            (_List&);
    _String*    GetSaveableListOfUserParameters     (void);
    void        TrimMemory                          (void);
    _VariableContainer* GetTheParent                (void) {
        return theParent;
    }

protected: // data members

    _SimpleList         *iVariables,
                        *dVariables,
                        *gVariables;

    void               SortVars (void);

    long                theModel;   // model template for the container
    _VariableContainer  *theParent; // a higher level container, if there is one.

};

#endif
