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
    virtual     long        ObjectClass (void) {
        return varValue?varValue->ObjectClass():((varFormula&&varFormula->theFormula.lLength)?varFormula->ObjectClass():1);
    }
    void        SetIndex (long i) {
        theIndex = i;
    }
    long        GetIndex (void) {
        return theIndex;
    }
    virtual     void        ScanForVariables (_AVLList& l, bool globals = false) {
        if (varValue) {
            varValue->ScanForVariables (l, globals);
        }
        if (varFormula && varFormula->theFormula.lLength) {
            varFormula->ScanFForVariables(l,globals);
        }
    }

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

    _String*    GetName                 (void) {
        return theName;
    }
    _String*    GetFormulaString        (void) {
        return varFormula?(_String*)varFormula->toStr():(_String*)empty.makeDynamic();
    }

    virtual     void        CompileListOfDependents (_SimpleList&);


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

#endif
