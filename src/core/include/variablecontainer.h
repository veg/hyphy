#ifndef     __VARIABLECONTAINER__
#define     __VARIABLECONTAINER__

#include "baseobj.h"
#include "avllist.h"
#include "operation.h"
#include "variable.h"

//__________________________________________________________________________________

// this class defines a computational (or storage) class which, as a variable, may contain
// other variables locally.

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

    virtual     void        MarkDone (void);

    // variable access/operation functions

    virtual     bool        IsContainer                 (void) {
        return true;
    }

    virtual     bool        HasChanged                  (void);
    virtual     bool        NeedToExponentiate          (bool = false);

    void        ScanAndAttachVariables      (void);

    virtual     void        ScanForVariables            (_AVLList&,_AVLList&);
    virtual     void        ScanForDVariables           (_AVLList&,_AVLList&);
    virtual     void        ScanForGVariables           (_AVLList&,_AVLList&);

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
    _Matrix*    GetModelMatrix              (void);
    _Matrix*    GetFreqMatrix               (void);
    bool        HasExplicitFormModel        (void);
    _Formula*   GetExplicitFormModel        (void);

    long        GetModelIndex               (void) {
        return theModel;
    }
    long        GetModelDimension           (void);
    /* 20100316 SLKP
        return the dimension of the model; needed to handle the case
        of explicit model exponentials
     */

    void        CopyMatrixParameters        (_VariableContainer*);
    void        GetListOfModelParameters    (_List&);
    _String*    GetSaveableListOfUserParameters
    (void);
    void        TrimMemory                  (void);
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


