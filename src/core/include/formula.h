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

#ifndef     __FORMULAS__
#define     __FORMULAS__

#include "baseobj.h"
#include "classes.h"
#include "defines.h"
#include "avllistxl.h"
#include "stack.h"
#include "operation.h"


#include "hy_string_buffer.h"
#include "formula_parsing_context.h"

class _Variable;
class _VariableContainer;
class _Polynomial;


union       _SimpleFormulaDatum {
    hyFloat value;
    hyPointer        reference;
};


enum _hyFormulaStringConversionMode  {
  kFormulaStringConversionNormal = 0L,
  kFormulaStringConversionSubstiteValues = 2L,
  kFormulaStringConversionReportRanges = 3L
};

class   _Formula {

    friend class _Variable;
    friend class _VariableContainer;
    
protected:

    unsigned    long    call_count;
    HBLObjectRef        recursion_calls;
    _List*              resultCache;
    _Stack              theStack;
    _List               theFormula;

    node<long>* theTree; // this formula converted to a tree for operation purposes
    // such as simplification, differentiation and printing.
    // trees store numbers referencing operations inside
    // "theFormula"


public:
    _Formula (void);
    _Formula (_String const&,_VariableContainer const* theParent=nil,_String* errorString = nil);
    _Formula (const _Polynomial *);
    
    long     ParseFormula (_String const&,_VariableContainer const* theParent=nil,_String* errorString = nil);
    
    _Formula (HBLObjectRef, bool isAVar = false);
    _Formula (_Formula const & rhs);
    const _Formula & operator = (_Formula const & rhs);
    virtual ~_Formula (void);
    HBLObjectRef   Compute             (long = 0, _VariableContainer const* = nil, _List* additionalCacheArguments = nil, _String *errMsg = nil, long object_type = HY_ANY_OBJECT);
    // compute the value of the formula
    // 1st argument : execute from this instruction onwards
    // see the commend for ExecuteFormula for the second argument

    bool        IsEmpty             (void) const; // is there anything in the formula
    long        NumberOperations    (void) const; // how many ops in the formula?

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


    _MathObject*ConstructPolynomial (void);

    virtual void        Initialize          (void);
    virtual void        Duplicate           (_Formula const *);
    void        DuplicateReference          (const _Formula*);
    virtual BaseRef     makeDynamic         (void) const;
    virtual BaseRef     toStr               (_hyFormulaStringConversionMode mode, _List* matchNames = nil, bool = false);
    _StringBuffer const     toRPN               (_hyFormulaStringConversionMode mode, _List* matchNames = nil);

    virtual long        ObjectClass         (void);


    virtual void        ScanFForVariables   (_AVLList&l, bool includeGlobals = false, bool includeAll = false, bool includeCateg = true, bool skipMatrixAssignments = false, _AVLListX* tagger = nil, long weight = 0) const;
    virtual void        ScanFForType        (_SimpleList&,  int);
    /* SLKP 20100716:
            A simple utility function to retrieve all variables of a given type
     */

    virtual bool        CheckFForDependence (long, bool checkAll = false);
    virtual bool        CheckFForDependence (_AVLList const & indices, bool checkAll = false);
    
    _List&      GetList             (void) {
        return theFormula;
    }

    bool        HasChanged          (bool = false); // does  the formula need recomputing
    bool        HasChangedSimple    (_SimpleList&);
    bool        EqualFormula        (_Formula*);
    bool        IsAConstant         (void); //  does this formula include variables, or is it just a constant?
    bool        IsConstant          (bool strict = false); //  does this formula depend on something other that constants and fixed parameters?
    bool        DependsOnVariable   (long);
    bool        IsArrayAccess       (void); // check to see if this formula performs a matrix access
    /*
        SLKP 20090315: added a missing utility function
        given a variable index as an argument, returns true if
        the formula depends on a it; false otherwise
    */
    _Operation* GetIthTerm          (long) const;
    /*
        SLKP 20090315: added a missing utility function
        given an index (i) as the argument, the function retrieves
        the i-th term of the formula
    */

    _Operation* ItemAt          (long) const;
    /*
     Same as GetIthTerm, but no range checking
     */

    unsigned long Length            (void) const {return theFormula.lLength;}

    void        Clear               (void);
    HBLObjectRef   GetTheMatrix        (void);

    void        PushTerm            (BaseRef);

    /* 20151008: if the argument is a _List, then treat as a list of _Operations and push them onto this formula (increment reference counters as well)
                 otherwise assume it's a MathObject and push it to this forumla (+1 reference counter)
                 dynamic_cast is used to determine what type of object this is

    */

    bool        AmISimple           (long& stack_depth, _AVLList& variable_index);
    long        StackDepth          (long start_at = 0L, long end_at = -1L) const;
      /**
        starting at operation 'start_at', counting up to 'end_at' (-1 == the end),
        evaluate how many values would be on the stack after the execution of these commands
       */
    bool        ConvertToSimple     (_AVLList& variableIndex);
    void        ConvertFromSimple   (_AVLList const& variableIndex);
    void        ConvertFromSimpleList   (_SimpleList const& variableIndex);
    void        SimplifyConstants   (void);
    _Variable * Dereference         (bool, _hyExecutionContext* = _hyDefaultExecutionContext);

    hyFloat  ComputeSimple       (_SimpleFormulaDatum* stack, _SimpleFormulaDatum* varValues) ;

    hyFloat  Newton              (_Formula&, _Variable*,  hyFloat, hyFloat, hyFloat);
    hyFloat  Newton              (_Formula&, hyFloat, hyFloat, hyFloat, _Variable*);
    hyFloat  Newton              (_Variable*,  hyFloat, hyFloat, hyFloat, hyFloat);
    hyFloat  Newton              (_Variable*,hyFloat, hyFloat, hyFloat);

    hyFloat  Brent               (_Variable*, hyFloat, hyFloat, hyFloat = 1.e-7, _List* = nil, hyFloat = 0.);

    hyFloat  Integral            (_Variable*,hyFloat, hyFloat, bool inifinite = false);
    hyFloat  MeanIntegral        (_Variable*,hyFloat, hyFloat, bool inifinite = false);
    _Formula*   Differentiate       (_String const&, bool = true, bool convert_from_tree = true);
    node<long>* InternalDifferentiate
    (node<long>*, long,_SimpleList const &, _Formula  * const *, _Formula&);

    bool        InternalSimplify    (node<long>*);

    void        LocalizeFormula           (_Formula&, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv);
    void        ConvertMatrixArgumentsToSimpleOrComplexForm (bool);
    long        ExtractMatrixExpArguments        (_List*);

    virtual     _Formula const operator + (const _Formula&);
    virtual     _Formula const operator - (const _Formula&);
    virtual     _Formula const operator * (const _Formula&);
    virtual     _Formula const operator / (const _Formula&);
    virtual     _Formula const operator ^ (const _Formula&);

    static      _Formula*        PatchFormulasTogether (const _Formula& op1, const _Formula& op2, const char op_code);
    static      _Formula*        PatchFormulasTogether (const _Formula& op1, HBLObjectRef op2, const char op_code);
    
    void        ScanFormulaForHBLFunctions (_AVLListX& collection , bool recursive);
  
  
    /** A compute and forget utility function.
        Parse an expression, optionally check to see that it's of the right type and return the value
        
      @param expression : the string to parse
      @param use_exceptions : if true, throw const _String exceptions, otherwise handle errors directly 
      @param requested_type: return nil if the computed value is not of this type
      @param formula parsing contrext
     
      @return expression value or nil; the value needs to be managed by the caller
      Revision history
          20170921 : SLKP initial implementation 
     
    */
     
    static      HBLObjectRef          ParseAndCompute (_String const& expression, bool use_exceptions = false, long requested_type = HY_ANY_OBJECT, _hyExecutionContext * context = nil);


protected:

    void        SubtreeToString     (_StringBuffer & result, node<long>* top_node, int op_level, _List* match_names, _Operation* this_node_op, _hyFormulaStringConversionMode mode = kFormulaStringConversionNormal);
    void        ConvertToTree       (bool err_msg = true);
    void        ConvertFromTree     (void);
    bool        CheckSimpleTerm     (HBLObjectRef);
    node<long>* DuplicateFormula    (node<long>*,_Formula&) const;


};

#endif
