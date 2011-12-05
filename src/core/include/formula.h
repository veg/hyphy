#ifndef     __FORMULAS__
#define     __FORMULAS__

#include "baseobj.h"
#include "classes.h"
#include "defines.h"
#include "avllist.h"
#include "stack.h"
#include "operation.h"

union       _SimpleFormulaDatum {
    _Parameter value;
    Ptr        reference;
};

class _Variable;
class _VariableContainer;
class   _Formula   // a computational formula
{

    friend class _Variable;
    friend class _VariableContainer;

public:
    _Formula (void);
    _Formula (_String&,_VariableContainer* theParent=nil,bool errors=true);
    _Formula (_PMathObj, bool isAVar = false);
    ~_Formula (void);
    _PMathObj   Compute             (long = 0, _VariableContainer* = nil);
    // compute the value of the formula
    // 1st argument : execute from this instruction onwards
    // see the commend for ExecuteFormula for the second argument

    bool        IsEmpty             (void); // is there anything in the formula
    long        NumberOperations    (void); // how many ops in the formula?

    friend  long        Parse               (_Formula*, _String&, long&, _VariableContainer* = nil, _Formula* = nil, bool flagErrors = true, bool* isVolatile = nil); // the parser
    friend  long        ExecuteFormula      (_Formula*, _Formula*, long, long, _VariableContainer* = nil);
    // the execution block for "compiled formulae
    /*
     SLKP 20100119: added an execution name space to allow correct scoping of "pass-by-reference"
                    arguments when calling ExecuteAFile within a namespace.

                    e.g. in

                    function foo (var&)
                    {
                        ...
                    }

                    foo ("varID");

                    varID may need to be prefixed by a namespace ID.
     */

    bool        CheckFormula        (void); // check to see if this formula is valid and compute the obj class

    _MathObject*ConstructPolynomial (void);

    virtual void        Initialize          (void);
    virtual void        Duplicate           (BaseRef);
    void        DuplicateReference  (_Formula*);
    virtual BaseRef     makeDynamic         (void);
    virtual BaseRef     toStr               (_List* matchNames = nil, bool = false);

    virtual long        ObjectClass         (void);


    virtual void        ScanFForVariables   (_AVLList&l, bool includeGlobals = false, bool includeAll = false, bool includeCateg = true, bool = false);
    virtual void        ScanFForType        (_SimpleList&,  int);
    /* SLKP 20100716:
            A simple utility function to retrieve all variables of a given type
     */

    virtual bool        CheckFForDependence (long, bool checkAll = false);
    _List&      GetList             (void) {
        return theFormula;
    }

    bool        HasChanged          (bool = false); // does  the formula need recomputing
    bool        HasChangedSimple    (_SimpleList&);
    bool        EqualFormula        (_Formula*);
    bool        IsAConstant         (void); //  does this formula include variables, or is it just a constant?
    bool        IsConstant          (void); //  does this formula depend on something other that constants and fixed parameters?
    bool        DependsOnVariable   (long);
    /*
        SLKP 20090315: added a missing utility function
        given a variable index as an argument, returns true if
        the formula depends on a it; false otherwise
    */
    _Operation* GetIthTerm          (long);
    /*
        SLKP 20090315: added a missing utility function
        given an index (i) as the argument, the function retrieves
        the i-th term of the formula
    */
    void        Clear               (void);
    _PMathObj   GetTheMatrix        (void);

    bool        AmISimple           (long& stackDepth, _SimpleList& variableIndex);
    void        ConvertToSimple     (_SimpleList& variableIndex);
    void        ConvertFromSimple   (_SimpleList& variableIndex);
    void        SimplifyConstants   (void);

    _Parameter  ComputeSimple       (_SimpleFormulaDatum* stack, _SimpleFormulaDatum* varValues) ;

    _Parameter  Newton              (_Formula&, _Variable*,  _Parameter, _Parameter, _Parameter);
    _Parameter  Newton              (_Formula&, _Parameter, _Parameter, _Parameter, _Variable*);
    _Parameter  Newton              (_Variable*,  _Parameter, _Parameter, _Parameter, _Parameter);
    _Parameter  Newton              (_Variable*,_Parameter, _Parameter, _Parameter);

    _Parameter  Brent               (_Variable*, _Parameter, _Parameter, _Parameter = 1.e-7, _List* = nil, _Parameter = 0.);

    _Parameter  Integral            (_Variable*,_Parameter, _Parameter, bool inifinite = false);
    _Parameter  MeanIntegral        (_Variable*,_Parameter, _Parameter, bool inifinite = false);
    _Formula*   Differentiate       (_String, bool = true);
    node<long>* InternalDifferentiate
    (node<long>*, long,_SimpleList&, _SimpleList&, _Formula&);

    bool        InternalSimplify    (node<long>*);

    void        LocalizeFormula     (_Formula&, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv);
    void        ConvertMatrixArgumentsToSimpleOrComplexForm (bool);

protected:

    void        internalToStr       (_String& result,node<long>*, char opLevel, _List* matchNames, _Operation* = nil);
    void        ConvertToTree       (void);
    void        ConvertFromTree     (void);
    bool        CheckSimpleTerm     (_PMathObj);
    node<long>* DuplicateFormula    (node<long>*,_Formula&);

    _List       theFormula,
                *resultCache;

    _Stack      theStack;
    node<long>* theTree; // this formula converted to a tree for operation purposes
    // such as simplification, differentiation and printing.
    // trees store numbers referencing operations inside
    // "theFormula"

};

#endif
