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

#ifndef     __FORMULAS__
#define     __FORMULAS__

#include "baseobj.h"
#include "classes.h"
#include "defines.h"
#include "avllistx.h"
#include "stack.h"
#include "operation.h"

class _Variable;
class _VariableContainer;


union       _SimpleFormulaDatum {
    _Parameter value;
    Ptr        reference;
};


class _FormulaParsingContext {

    long                 assignment_ref_id;
    char                 assignment_ref_type;
    bool                 is_volatile;
    _String            * err_msg;
    _VariableContainer * formula_scope;
    
    public:
        _FormulaParsingContext (_String* = nil, _VariableContainer* = nil);
        bool&       isVolatile (void)                   { return is_volatile; }
        long&       assignmentRefID (void)              { return assignment_ref_id; }
        char&       assignmentRefType (void)            { return assignment_ref_type;} 
        _String*    errMsg (void)                       { return err_msg; }
        _VariableContainer* formulaScope (void)         { return formula_scope; }     
        _String     contextualizeRef (_String&);
};


class   _Formula   // a computational formula
{

    friend class _Variable;
    friend class _VariableContainer;

public:
    _Formula (void);
    _Formula (_String&,_VariableContainer* theParent=nil,_String* errorString = nil);
    _Formula (_PMathObj, bool isAVar = false);
    virtual ~_Formula (void);
    _PMathObj   Compute             (long = 0, _VariableContainer* = nil, _List* additionalCacheArguments = nil, _String *errMsg = nil);
    // compute the value of the formula
    // 1st argument : execute from this instruction onwards
    // see the commend for ExecuteFormula for the second argument

    bool        IsEmpty             (void); // is there anything in the formula
    long        NumberOperations    (void); // how many ops in the formula?

    friend  long        Parse               (_Formula*, _String&, _FormulaParsingContext&, _Formula*); 
    // the main expression parser
    
    friend  long        ExecuteFormula      (_Formula*, _Formula*, long, long, _VariableContainer*, char);
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
    void        DuplicateReference  (const _Formula*);
    virtual BaseRef     makeDynamic         (void);
    virtual BaseRef     toStr               (_List* matchNames = nil, bool = false);

    virtual long        ObjectClass         (void);


    virtual void        ScanFForVariables   (_AVLList&l, bool includeGlobals = false, bool includeAll = false, bool includeCateg = true, bool skipMatrixAssignments = false, _AVLListX* tagger = nil, long weight = 0);
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
    bool        IsArrayAccess       (void); // check to see if this formula performs a matrix access
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
    bool        ConvertToSimple     (_SimpleList& variableIndex);
    void        ConvertFromSimple   (_SimpleList& variableIndex);
    void        SimplifyConstants   (void);
    _Variable * Dereference         (bool, _hyExecutionContext* = _hyDefaultExecutionContext);

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

    void        LocalizeFormula           (_Formula&, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv);
    void        ConvertMatrixArgumentsToSimpleOrComplexForm (bool);
    long        ExtractMatrixExpArguments        (_List*);
    
    virtual     _Formula operator + (const _Formula&);
    virtual     _Formula operator - (const _Formula&);
    virtual     _Formula operator * (const _Formula&);
    virtual     _Formula operator / (const _Formula&);
    virtual     _Formula operator ^ (const _Formula&);
    
    _Formula&        PatchFormulasTogether (_Formula&, const _Formula&, const char op_code);

protected:

    void        internalToStr       (_String& result,node<long>*, char opLevel, _List* matchNames, _Operation* = nil);
    void        ConvertToTree       (bool err_msg = true);
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
