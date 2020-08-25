/*

 HyPhy - Hypothesis Testing Using Phylogenies.

 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon42@uwo.ca)
 Steven Weaver (sweaver@temple.edu)

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

#include <math.h>
#include <float.h>
#include "defines.h"
#include "formula.h"


#include "parser.h"
#include "batchlan.h"
#include "function_templates.h"
#include "global_things.h"
#include "global_object_lists.h"
#include "polynoml.h"

using namespace hy_global;
using namespace hyphy_global_objects;

//Constants
hyFloat const sqrtPi = 1.77245385090551603,
              twoOverSqrtPi = 2./sqrtPi;

//__________________________________________________________________________________


_Formula::_Formula (void) {
    Initialize();
}

//__________________________________________________________________________________

_Formula::_Formula (HBLObjectRef p, bool is_a_var) {
    theTree     = nil;
    resultCache = nil;
    recursion_calls = nil;
    call_count = 0UL;

    if (!is_a_var) {
        theFormula.AppendNewInstance (new _Operation (p));
    } else {
        _Variable* v = (_Variable*)p;
        theFormula.AppendNewInstance (new _Operation (true,*v->GetName(),v->IsGlobal(), nil));
    }
}

//__________________________________________________________________________________
_Formula::_Formula (_String const &s, _VariableContainer const* theParent, _String* reportErrors) {
    ParseFormula (s, theParent, reportErrors);
}

//__________________________________________________________________________________

_Formula::_Formula(const _Polynomial* source) {
    Initialize();
    _PolynomialData * poly_terms = source->GetTheTerms();
    
    _List poly_vars;
    for (long i = 0L; i < source->GetNoVariables(); i++) {
        poly_vars << source->GetIthVariable(i);
    }
    
    for (long i = 0L; i < poly_terms->NumberOfTerms(); i++) {
        hyFloat c = poly_terms->GetCoeff(i);
        /*theFormula.AppendNewInstance (new _Operation (new _Constant (poly_terms->GetCoeff(i))));*/
        bool started_term = false;
        if (!CheckEqual(c, 1.0)) {
            theFormula.AppendNewInstance (new _Operation (new _Constant (poly_terms->GetCoeff(i))));
            started_term = true;
        }
        long *exponents = poly_terms->GetTerm(i);
        poly_vars.ForEach ([this, exponents, &started_term] (BaseRefConst v, unsigned long idx) -> void {
            if (exponents[idx] != 0L) {
                this->theFormula.AppendNewInstance (new _Operation (*((_Variable const*)v)));
                if (exponents[idx] != 1L) {
                    this->theFormula.AppendNewInstance (new _Operation (new _Constant (exponents[idx])));
                    this->theFormula.AppendNewInstance (new _Operation (HY_OP_CODE_POWER, 2));
                }
                if (started_term) {
                    this->theFormula.AppendNewInstance (new _Operation (HY_OP_CODE_MUL, 2));
                } else {
                    started_term = true;
                }
            }
        });
        if (!started_term) {
            this->theFormula.AppendNewInstance (new _Operation (new _Constant (1.0)));
        }
        if (i) {
            this->theFormula.AppendNewInstance (new _Operation (HY_OP_CODE_ADD, 2));
        }
    }
}

//__________________________________________________________________________________
void _Formula::Initialize (void) {
    theTree     = nil;
    resultCache = nil;
    call_count = 0UL;
    recursion_calls = nil;
}

//__________________________________________________________________________________
_Formula::_Formula (_Formula const& rhs) {
    *this = rhs;
}

//__________________________________________________________________________________
const _Formula& _Formula::operator =  (_Formula const& rhs) {
    if (this != &rhs) {
        Initialize();
        Duplicate (&rhs);
    }
    return *this;
}

//__________________________________________________________________________________
void _Formula::Duplicate  (_Formula const * f_cast) {

    theFormula.Duplicate       (& f_cast->theFormula);
    theStack.theStack.Duplicate(& f_cast->theStack.theStack);
    call_count = f_cast->call_count;
    if (f_cast->recursion_calls) {
      recursion_calls = (HBLObjectRef)f_cast->recursion_calls->makeDynamic();
    } else {
      recursion_calls = nil;
    }
    if (f_cast->theTree) {
        theTree = f_cast->theTree->duplicate_tree();
    } else {
        theTree = nil;
    }

    if (f_cast->resultCache) {
        resultCache = (_List*)f_cast->resultCache->makeDynamic();
    } else {
        resultCache = nil;
    }
}

//__________________________________________________________________________________
void _Formula::DuplicateReference  (const _Formula* f) {
    for (unsigned long i=0; i<f->theFormula.countitems(); i++) {
        _Operation *ith_term = f->GetIthTerm(i);
        // -2 is used to code for value substitutions (x__)
        if (ith_term->IsValueSubstitution()) {
            theFormula.AppendNewInstance(new _Operation ((HBLObjectRef)LocateVar (ith_term->GetAVariable())->Compute()->makeDynamic()));
        } else {
            theFormula && ith_term;
        }
    }
}

  //__________________________________________________________________________________

HBLObjectRef _Formula::ParseAndCompute (_String const& expression, bool use_exceptions, long type, _hyExecutionContext * context) {
  
  _String error_message;
  _Formula f;
  
  long parse_result = f.ParseFormula (expression, context ? context -> GetContext() : nil , &error_message);
  
  try {
    if (error_message.nonempty()) {
      throw (_String ("Failed to parse ") & expression.Enquote () & " with the following error: " & error_message);
    }
    if (parse_result != HY_FORMULA_EXPRESSION) {
      throw (expression.Enquote () & " did not parse to a simple expression");
    }
    if (f.IsEmpty ()) {
      throw (expression.Enquote () & " parsed to an empty expression");
    }
    if (!(type == HY_ANY_OBJECT || f.ObjectClass() == type)) {
        // TODO SLKP 20170704: ObjectClass will compute the expression with current values which may fail
      throw (expression.Enquote () & " did not evaluate to a " & FetchObjectNameFromType (type));
    }
  } catch (_String const & e) {
    if (use_exceptions) {
      throw (e);
    } else {
      HandleApplicationError (e);
      return nil;
    }
  }
  
  HBLObjectRef result = f.Compute();
  result->AddAReference();
  return result;
}


//__________________________________________________________________________________
BaseRef _Formula::makeDynamic (void) const {
    _Formula * res = new _Formula;
    res->Duplicate(this);
    return (BaseRef)res;
}
//__________________________________________________________________________________

_Formula::~_Formula (void) {
    Clear();
}

//__________________________________________________________________________________
void    _Formula::Clear (void) {
    if (theTree) {
        theTree->delete_tree();
        delete theTree;
    }
    theTree = nil;
    if (resultCache) {
        DeleteAndZeroObject(resultCache);
    }

    theFormula.Clear();
    if (recursion_calls) {
      delete (recursion_calls);
    }

//  theStack.Clear();
}
//__________________________________________________________________________________
_StringBuffer const _Formula::toRPN (_hyFormulaStringConversionMode mode, _List* matched_names) {
    _StringBuffer r;
    if (theFormula.countitems()) {
        SubtreeToString (r,nil,0,nil,GetIthTerm(0UL), mode);
        for (unsigned long k=1UL; k<theFormula.countitems(); k++) {
            r <<'|';
            SubtreeToString (r,nil,0,nil,GetIthTerm(k), mode);
        }
    }
    return r;
}
//__________________________________________________________________________________
BaseRef _Formula::toStr (_hyFormulaStringConversionMode mode, _List* matched_names, bool drop_tree) {
    ConvertToTree(false);

    _StringBuffer * result = new _StringBuffer (64UL);

    long          savepd = print_digit_specification;
    print_digit_specification          = 0L;

    if (theTree) { // there is something to do
        SubtreeToString (*result, theTree, -1, matched_names, nil, mode);
    } else {
        if (theFormula.countitems()) {
            (*result) << "RPN:" << toRPN (mode, matched_names);
        }
    }

    print_digit_specification = savepd;
    result->TrimSpace();
    if (theTree && drop_tree) {
        theTree->delete_tree();
        delete theTree;
        theTree = nil;
    }
  
    return result;
}
//__________________________________________________________________________________
node<long>* _Formula::DuplicateFormula (node<long>* src, _Formula& tgt) const {
    
    node<long>* copied = new node<long>;
    
    for (long k=1L; k<=src->get_num_nodes(); k++) {
        copied->add_node (*DuplicateFormula (src->go_down (k), tgt));
    }

    copied->in_object = tgt.theFormula.lLength;
    tgt.theFormula && theFormula.GetItem(src->in_object);

    return     copied;
}

//__________________________________________________________________________________
_Formula* _Formula::Differentiate (_String const & var_name, bool bail, bool convert_from_tree) {
   
    _Variable * dx = FetchVar(LocateVarByName(var_name));
    
    if (!dx) { // create the dX variable if it's not here ye
        return new _Formula (new _Constant (0.0));
    }

    long dx_id = dx->get_index();

     ConvertToTree    ();
    
     //printf ("\n **** Diff %s on %s\n\n", _String ((_String*)toStr(kFormulaStringConversionNormal)).get_str(), var_name.get_str());

    _SimpleList  var_refs = PopulateAndSort ([&] (_AVLList & parameter_list) -> void {
                            this->ScanFForVariables (parameter_list, true, true, true);
    });
    
    if (var_refs.empty()) { // handle the case of constant expressions here
        return new _Formula (new _Constant (0.0));
    }

    _Formula*     res = new _Formula ();
    _Formula    ** dydx = new _Formula* [var_refs.countitems()] {0};// stores precomputed derivatives for all the
    
    auto dydx_cleanup = [&] () -> void {
        for (unsigned long k=0UL; k < var_refs.countitems(); k++) {
            if (dydx[k]) {
                delete dydx[k];
            }
        }
        delete [] dydx;
    };
   
    node<long>*           dTree = nil;

    try {
        for (unsigned long k=0UL; k < var_refs.countitems(); k++) {
            _Variable* thisVar = LocateVar (var_refs.GetElement(k));
            _Formula * dYdX;

            if (thisVar->IsIndependent()) {
                dYdX = new _Formula ((*thisVar->GetName() == var_name)?new _Constant (1.0):new _Constant (0.0));
            } else {
                dYdX = thisVar->varFormula->Differentiate (var_name, bail, false);
                if (dYdX->IsEmpty()) {
                    delete dYdX;
                    dydx_cleanup ();
                    return res;
                }
            }
            dYdX->ConvertToTree();
            dydx [k] = dYdX;
          }

            // SortLists             (&varRefs, &dydx);
            // this is already sorted coming from PopulateAndSort
            

        if (!(dTree = InternalDifferentiate (theTree, dx_id, var_refs, dydx, *res))) {
            throw (_String ("Differentiation of ") & _String((_String*)toStr(kFormulaStringConversionNormal)) & " failed.");
        }
    } catch (_String const &e) {
        dydx_cleanup ();
        
        if (bail) {
            HandleApplicationError (e);
            res->Clear();
            return       res;
        } else {
            delete res;
            return nil;
        }
    }

    dydx_cleanup ();

    //res->theFormula.AppendNewInstance (new _Operation(new _Constant (0.0))) ; // why is this here?
    res->theTree         = dTree;
    res->InternalSimplify (dTree);
    //printf ("%s\n", _String ((_String*)res->theFormula.toStr(kFormulaStringConversionNormal)).get_str());
    
    if (convert_from_tree) {
      res->ConvertFromTree  ();
    }
    //printf ("%s\n", _String ((_String*)res->theFormula.toStr(kFormulaStringConversionNormal)).get_str());

    // try polynomial simplification
    _Polynomial*    is_poly = (_Polynomial*)res->ConstructPolynomial();
    if (is_poly) {
        _Formula * simplified_polynomial = new _Formula (is_poly);
        //printf ("%s\n", _String ((_String*)simplified_polynomial->toStr(kFormulaStringConversionNormal)).get_str());
        delete res;
        
        //printf ("\n RESULT : %s \n", _String ((_String*)simplified_polynomial->toStr(kFormulaStringConversionNormal)).get_str());

        return simplified_polynomial;
    } else {
        //printf ("%s\n", _String ((_String*)res->toStr(kFormulaStringConversionNormal)).get_str());
    }
    
    return res;

}

//__________________________________________________________________________________
node<long>* _Formula::InternalDifferentiate (node<long>* currentSubExpression, long varID, _SimpleList const & varRefs, _Formula  * const * dydx, _Formula& tgt) {
    _Operation * op = GetIthTerm(currentSubExpression->in_object);
    
    if (op->theData!=-1) {
        long k     = varRefs.BinaryFind (op->GetAVariable());
        if (k<0L) {
            return nil;
        }
        
        _Formula const* dYdX = dydx[k];
        return dYdX->DuplicateFormula (dYdX->theTree, tgt);
    }
    
    if (op->theNumber) {
        _Formula src (new _Constant (0.0));
        src.ConvertToTree ();
        return   src.DuplicateFormula (src.theTree, tgt);
    }
    
    node<long>* newNode = new node<long>;
    node <long> *  created_nodes [7]{nil};
    created_nodes [0] = newNode;
    
    try {
        switch (op->opCode) {
            case HY_OP_CODE_MUL: {
                /**
                 d (X*Y) = X*dY + Y*dX
                 */
                
                created_nodes [1]  = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes [1]) {
                    throw (2);
                }
                
                node<long>*       YtimesDX = new node<long>;
                created_nodes[3] = YtimesDX;
                YtimesDX->add_node (*created_nodes [1],*DuplicateFormula (currentSubExpression->go_down(2),tgt));
                YtimesDX->in_object = tgt.theFormula.AppendNewInstance (new _Operation (HY_OP_CODE_MUL,2)) - 1;
                

                created_nodes [2]  = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
                if (!created_nodes [2]) {
                    throw (3);
                }
                node<long>*       XtimesDY = new node<long>;
                XtimesDY->add_node (*created_nodes [2],*DuplicateFormula (currentSubExpression->go_down(1),tgt));
                XtimesDY->in_object = tgt.theFormula.AppendNewInstance (new _Operation (HY_OP_CODE_MUL,2)) - 1;
                newNode->add_node  (*YtimesDX,*XtimesDY);
                newNode->in_object  = tgt.theFormula. AppendNewInstance (new _Operation (HY_OP_CODE_ADD,2)) - 1;
                
                return          newNode;
            }
            break;
                
            case HY_OP_CODE_ADD: // +
            case HY_OP_CODE_SUB: { // -
                
                // TODO SLKP 20170426: clean up the rest of the function
                
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes[1]) {
                    throw (1);
                }
                
                long      isUnary = (currentSubExpression->get_num_nodes()==1);
                
                if (!isUnary) {
                    created_nodes[2] = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
                    if (!created_nodes[2]) {
                        throw (2);
                    }
                }
                
                
                newNode->add_node (*created_nodes[1]);
                if (!isUnary) {
                    newNode->add_node (*created_nodes[2]);
                }
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (*op)) - 1;
                return          newNode;
            }
            break;
                
            case HY_OP_CODE_DIV: { // /
                
                /*
                    d (X/Y) = (y * dX - x * dY ) * Y^{-2}
                */
                created_nodes [1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes [1]) {
                    throw (2);
                }
                
                
                node<long>*       yDx = new node<long>;
                created_nodes[3] = yDx;
                yDx->add_node (*created_nodes [1],*DuplicateFormula (currentSubExpression->go_down(2),tgt)); // dX * Y
                yDx->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL,2))-1;
                
                created_nodes [2] = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
                if (!created_nodes [2]) {
                    throw (4);
                }

                node<long>*       y_raise_m2 = new node<long>;
                node<long>*       yDx_minus_xDy = new node<long>;
                node<long>*       xDy = new node<long>;
                node<long>*       m2  = new node<long>;

                xDy->add_node (*created_nodes [2],*DuplicateFormula (currentSubExpression->go_down(1),tgt)); // dY * X
                xDy->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL,2))-1;
                
                yDx_minus_xDy->add_node  (*yDx,*xDy);
                yDx_minus_xDy->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_SUB,2))-1; // numerator
                
                y_raise_m2->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt), *m2);
                m2->in_object =  tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (-2.)))-1;
                y_raise_m2->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_POWER,2))-1; // y^{-2}

                newNode->add_node(*yDx_minus_xDy, *y_raise_m2);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL,2)) - 1;
                
                return          newNode;
            }
                break;
                
            case HY_OP_CODE_ARCTAN: { // Arctan
                
                // Arctax (y)' = y' / (1+y^2)
                
                created_nodes [1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes[1]) {
                    throw (1);
                }
                
                node<long>*       one   = new node<long>,
                          *       two   = new node<long>,
                          *       denom = new node<long>,
                          *       y_squared = new node<long>;
                
                y_squared->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt), *two);
                two->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (2.))) - 1;
                y_squared->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_POWER,2)) - 1;
                denom->add_node(*y_squared, *one);
                one->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (1.))) - 1;
                denom->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_ADD,2)) - 1;
                newNode->add_node(*created_nodes[1], *denom);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_DIV,2)) - 1;
                
                return          newNode;
            }
                break;
                
            case HY_OP_CODE_COS: {
                // d Cos (f[x]) = - f'[x] * Sin (x)
                
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (! created_nodes[1] ) {
                    throw (1);
                }
                
                node <long> * func_fx = new node <long>;
                func_fx->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
                func_fx->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_SIN,1)) - 1;
                node <long> * product_term = new node <long>;
                product_term->add_node (*created_nodes[1], *func_fx);
                product_term->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL,2)) - 1;
                newNode->add_node(*product_term);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_SUB,1)) - 1;

                return          newNode;
            }
                break;
                
            case HY_OP_CODE_ERF: { // Erf
                
                // d ERF (f [X]) = Exp(-f[X]^2) / 2Pi
                
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes[1]) {
                    throw (1);
                }
                
                node <long> * f_squared = new node<long>,
                            * two = new node<long>;
                
                f_squared->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt), *two);
                two->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (2.))) - 1L;
                f_squared->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_POWER, 2)) - 1L;
 
                node <long> * minus_f_squared = new node<long>;
                minus_f_squared->add_node(*f_squared);
                minus_f_squared->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_SUB, 1)) - 1L;

                node <long> * exp_minus_f_squared = new node<long>;
                exp_minus_f_squared->add_node(*minus_f_squared);
                exp_minus_f_squared->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_EXP, 1)) - 1L;

                
                node <long> * two_pi        = new node<long>,
                            * div_by_two_pi = new node<long>;

                div_by_two_pi->add_node(*exp_minus_f_squared, *two_pi);
                two_pi->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (twoOverSqrtPi))) - 1L;
                div_by_two_pi ->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_DIV, 2)) - 1L;
                
                newNode->add_node(*created_nodes[1], *div_by_two_pi);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL, 2)) - 1L;
                
                return          newNode;
            }
                break;
                
            case HY_OP_CODE_EXP: // HY_OP_CODE_EXP
            case HY_OP_CODE_SIN: { // HY_OP_CODE_SIN
                
                // d Exp (f[X]) = f'[X] Exp (f[X])
                // d Sin (f[X]) = f'[X] Cos (f[X])
                
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (! created_nodes[1] ) {
                    throw (1);
                }
                
                node <long> * func_fx = new node <long>;
                func_fx->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
                func_fx->in_object = tgt.theFormula.AppendNewInstance(new _Operation (op->opCode==HY_OP_CODE_SIN ? HY_OP_CODE_COS : HY_OP_CODE_EXP,1)) - 1;
                newNode->add_node (*created_nodes[1], *func_fx);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL,2)) - 1;
                
                return           newNode;
            }
                break;
                
            case HY_OP_CODE_LOG: {
                // Log (f(x))' = f'(x) * f(x)^{-1}
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes[1]) {
                    throw (2);
                }
                
                node <long> * n1 = new node<long>, * n2 = new node<long>;
                n2->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt), *n1); // f(x) , {-1}
                n1->in_object = tgt.theFormula.AppendNewInstance (new _Operation (new _Constant (-1))) - 1; // push -1
                n2->in_object = tgt.theFormula.AppendNewInstance (new _Operation (HY_OP_CODE_POWER,2)) - 1; // push '^'
                newNode->add_node (*created_nodes[1], *n2); // f'(x) , f(x) (^) {-1}
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL,2))-1; // push *
                
                return           newNode;
            }
            break;
                
            case HY_OP_CODE_SQRT: {
                
                // d (Sqrt (f[x])) = 0.5 * f'[x] / Sqrt (f[x]))
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes[1]) {
                    throw (1);
                }
                
                node <long>* half           = new node<long>;
                node <long>* half_times_dfx = new node<long>;
                
                half_times_dfx->add_node (*created_nodes[1],*half);
                half->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (0.5))) - 1;
                half_times_dfx->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL, 2)) - 1;
                
                node <long> * sqrt_fx = new node<long>;
                sqrt_fx->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
                sqrt_fx->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_SQRT, 1)) - 1;

                newNode->add_node  (*half_times_dfx,*sqrt_fx);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_DIV, 2)) - 1;                
                return          newNode;
            }
                break;
                
            case HY_OP_CODE_TAN: {
                // Tan (f[X]) = f'[X] / Cos ^ 2 (f[X])
                
                created_nodes[1] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                
                if (!created_nodes[1]) {
                    throw (1);
                }
 
                node <long>* two                    = new node<long>;
                node <long>* cos_f                  = new node<long>;
                node <long>* cos_f_squared         = new node<long>;

                cos_f->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
                cos_f->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_COS, 1)) - 1;
                cos_f_squared->add_node(*cos_f, *two);
                two->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (2.))) - 1;
                cos_f_squared->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_POWER, 2)) - 1;
                
                newNode->add_node(*created_nodes[1], *cos_f_squared);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_DIV, 2)) - 1;
                
                return          newNode;
            }
                break;
                
            case HY_OP_CODE_POWER: // ^
                // f[x]^g[x] (g'[x] Log[f[x]] + f'[x]g[x] * (f[x]^{-1}))
            {
                node <long> * f_raise_g = new node <long>;
                f_raise_g->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt),*DuplicateFormula (currentSubExpression->go_down(2),tgt));
                f_raise_g->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_POWER, 2)) - 1;
                created_nodes[1] = f_raise_g;
                created_nodes[2] = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
                if (!created_nodes[2]) {
                    throw (3);
                }
                node <long> * log_f = new node <long>;
                created_nodes[3] = log_f;
                log_f->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
                log_f->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_LOG, 1)) - 1;
                node <long> * dg_times_log_fx = new node <long>;
                created_nodes[4] = dg_times_log_fx;
                dg_times_log_fx->add_node (*created_nodes[2], *log_f);
                dg_times_log_fx->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL, 2)) - 1;
                created_nodes[5] = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);
                if (!created_nodes[5]) {
                    throw (6);
                }
                node <long> * df_times_g = new node <long>;
                df_times_g->add_node (*created_nodes[5],*DuplicateFormula (currentSubExpression->go_down(2),tgt));
                df_times_g->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL, 2)) - 1;
                node <long> * inverse_f = new node <long>;
                node <long> * minus1 = new node<long>;
                inverse_f->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt),*minus1);
                minus1->in_object = tgt.theFormula.AppendNewInstance(new _Operation (new _Constant (-1.))) - 1;
                inverse_f->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_POWER, 2)) - 1;
                node <long> * product_term = new node<long>; // f'[x]g[x] [*] (f[x]^{-1}
                product_term->add_node (*df_times_g, *inverse_f);
                product_term->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL, 2)) - 1;
                node <long> * sum_term = new node<long>; // g'[x] Log[f[x]] <+> f'[x]g[x] * (f[x]^{-1})
                sum_term->add_node (*dg_times_log_fx, *product_term);
                sum_term->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_ADD, 2)) - 1;
                newNode->add_node (*f_raise_g, *sum_term);
                newNode->in_object = tgt.theFormula.AppendNewInstance(new _Operation (HY_OP_CODE_MUL, 2)) - 1;

                
                return          newNode;
            }
        }
    } catch (int cleanup_length) {
        for (long i = 0L; i < cleanup_length; i++) {
            if (created_nodes [i]) {
                created_nodes[i]->delete_tree(true);
            }
        }
        return nil;
    }
    delete (newNode);
    return nil;
}



//__________________________________________________________________________________
bool _Formula::InternalSimplify (node<long>* top_node) {
// returns true if the subexpression at
// and below startnode is constant
    _Operation* op = GetIthTerm(top_node->get_data());
    long        n_children = top_node->get_num_nodes();

    if  (n_children == 0L) {
      return !op->IsAVariable();
    }


    bool        all_constant  = true,
                left_constant  = true,
                right_constant = (n_children>1L);

    long        prune_this_child = -1;

    hyFloat     evaluated_to      = 0.0;
    HBLObjectRef   replace_with      = nil;


    //printf ("InternalSimplify %x\n", startNode);

    for  (unsigned long k=1UL; k<=n_children; k++) {
        if (k==1UL) {
            left_constant = InternalSimplify (top_node->go_down(k));
        } else if (k==2UL) {
            right_constant = InternalSimplify (top_node->go_down(k));
        } else {
          if (!InternalSimplify (top_node->go_down(k))) {
            all_constant = false;
          }
        }
    }

    all_constant = all_constant && left_constant && (n_children==1 || right_constant);

    if (op->opCode > HY_OP_CODE_NONE) {
        if (all_constant) { // this executes the subxpression starting at the current node
            _Stack scrap;
            for  (unsigned long k=1UL; k<=n_children; k++) {
                ((_Operation*)theFormula (top_node->go_down(k)->get_data()))->Execute (scrap);
            }
            op->Execute (scrap);
            replace_with = (HBLObjectRef)scrap.Pop();//->makeDynamic();
        } else {
            if (left_constant||right_constant) {

                HBLObjectRef constant_value = ((_Operation*)theFormula (top_node->go_down(left_constant?1:2)->get_data()))->GetANumber();


                if (constant_value->ObjectClass() == NUMBER) {
                  evaluated_to  = constant_value->Value();

                  switch (op->opCode) {
                      case HY_OP_CODE_MUL: { // *
                          if (CheckEqual (evaluated_to,0.0)) { // *0 => 0
                                                         //printf ("*0\n");
                              replace_with = new _Constant (0.0);
                              break;
                          }
                          if (CheckEqual (evaluated_to,1.0)) { // x*1 => x
                                                         //printf ("*1\n");
                              prune_this_child = left_constant?1:2;
                              break;
                          }
                      }
                      break;

                      case HY_OP_CODE_ADD: { // +
                          if (CheckEqual (evaluated_to,0.0)) { // x+0 => x
                                                         //printf ("+0\n");
                              prune_this_child = left_constant?1:2;
                          }
                          break;
                      }

                      case HY_OP_CODE_SUB: { // x-0 => x
                                             // 0-x => -x

                          if (CheckEqual (evaluated_to,0.0)) {
                              //printf ("-0\n");
                             prune_this_child = left_constant? -2 : 2;
                          }
                          break;
                      }

                      case HY_OP_CODE_DIV: { // /
                          if (left_constant&&CheckEqual (evaluated_to,0.0)) { // 0/x => 0
                              replace_with = new _Constant (0.0);
                              //printf ("0/\n");

                              break;
                          }
                          if (right_constant&&CheckEqual (evaluated_to,1.0)) { // x/1 => x
                                                                      //printf ("/1\n");
                              prune_this_child = 2;
                              break;
                          }
                      }
                      break;

                      case HY_OP_CODE_POWER: { // ^
                          if (left_constant&&CheckEqual (evaluated_to,1.0)) { // 1^x => 1
                                                                     //printf ("1^\n");
                              replace_with = new _Constant (1.0);
                              break;
                          }
                          if (right_constant&&CheckEqual (evaluated_to,1.0)) { // x^1 => ?
                                                                      //printf ("^1\n");
                              prune_this_child = 1;
                              break;
                          }
                      }
                      break;
                  }
                }
            }
        }
    }

    if (replace_with) {
        all_constant = true;
        for  (int k=1; k <= n_children; k++) {
            top_node->go_down(k)->delete_tree(true);
        }
        top_node->kill_all_nodes();
        top_node->in_object = theFormula.lLength;
        theFormula < (new _Operation(replace_with));
    } else {

      if (prune_this_child !=- 1L) {
          if (prune_this_child > 0L) {
              top_node->go_down(prune_this_child)->delete_tree(true);
              top_node->kill_node (prune_this_child);
              node <long>*    pruned_tree = top_node->go_down(1);

              top_node->kill_all_nodes();

              for (unsigned long k=1; k<=pruned_tree->get_num_nodes(); k++) {
                  top_node->add_node (*pruned_tree->go_down(k));
              }
              top_node->in_object = pruned_tree->in_object;
              delete (pruned_tree);

          } else { // 0-? => -?
              top_node->go_down(1)->delete_tree(true);
              top_node->kill_node(1);
              op->SetTerms(1);
             //startNode->kill_node(1);
          }
      }
    }
    
    // check to see if this op may be of the kind F (F^-1 (X)), like Log (Exp (..)))
    
    if (!all_constant && top_node->get_num_nodes() == 1) {
        long check_op_code = GetIthTerm(top_node->get_data())->opCode;
        node <long>* child_node = top_node->go_down(1);
        if (child_node->get_num_nodes() == 1) {
            long child_op_code = GetIthTerm(child_node->get_data())->opCode;
            bool cancel_pair = _Operation::AreOpsInverse (child_op_code, check_op_code);
            if (cancel_pair) {
                node<long> * remaining_trunk = child_node->go_down(1);
                delete (top_node->go_down(1));
                top_node->in_object = remaining_trunk->in_object;
                top_node->kill_all_nodes();
                for (long k = 1; k <= remaining_trunk->get_num_nodes(); k++) {
                    top_node->add_node (*remaining_trunk->go_down (k));
                }
                delete remaining_trunk;
                
            }
        }
    }
    
    return all_constant;
}


//__________________________________________________________________________________
void _Formula::SubtreeToString (_StringBuffer & result, node<long>* top_node, int op_level, _List* match_names, _Operation* this_node_op, _hyFormulaStringConversionMode mode) {
    
    if (!this_node_op) {
        if (!top_node) {
          return;
        }
        this_node_op = GetIthTerm (top_node->get_data());
    }

    // decide what to do about this operation

    if (this_node_op->IsAVariable(false) || this_node_op->IsValueSubstitution()) {
        // this operation is just a variable - add ident to string and return
        long node_variable_index = this_node_op->GetAVariable();
        if (node_variable_index < 0L) {
          return;
        }
      
        if (mode != kFormulaStringConversionNormal) {
          
            _Variable *node_variable = LocateVar(node_variable_index);
          
            if (mode == kFormulaStringConversionSubstiteValues) {
                 if  (hy_x_variable && (node_variable->get_index()==hy_x_variable->get_index())) {
                    result << hy_x_variable->GetName();
                    return;
                }
            }

            HBLObjectRef replace_with = node_variable->Compute();

            if (replace_with->ObjectClass () == NUMBER) {
                if (mode == kFormulaStringConversionReportRanges) {
                    result << node_variable->GetName()
                           << '[';
                    result.AppendNewInstance(new _String (replace_with->Value()));
                    result << ':';
                    result.AppendNewInstance(new _String (node_variable->GetLowerBound()));
                    result << '-';
                    result.AppendNewInstance(new _String (node_variable->GetUpperBound()));
                    result << ']';

                } else {
                    result.AppendNewInstance(new _String (replace_with->Value()));
                }
            } else if (replace_with->ObjectClass () == STRING) {
                result.AppendNewInstance((_String*)replace_with->toStr());
                if (this_node_op->IsValueSubstitution()) {
                    result << "__";
                }
            } else {
                result << node_variable->GetName();
                if (this_node_op->IsValueSubstitution()) {
                    result << "__";
                }
            }
        } else {
            _String * variable_name = LocateVar(node_variable_index)->GetName();
            if (match_names) {

                long  f = ((_List*)(*match_names)(0))->FindObject (variable_name);

                if (f == kNotFound) {
                    result<<variable_name;
                } else {
                    result << (_String*)((_List*)(*match_names)(1))->GetItem(f);
                }
            } else {
                result<<variable_name;
            }
            if (this_node_op->IsValueSubstitution()) {
                result << "__";
            }
        }
        return;
    }

    long node_op_count = this_node_op->GetNoTerms();
    if (node_op_count > 0) {
        // a built-in operation or a function call
        // check if it's a built-in binary operation
      
        long f = _Operation::BinOpCode(this_node_op->GetCode());
      
        if (f != kNotFound) {
            // indeed - a binary operation is what we have. check if need to wrap the return in parentheses
          
            if (!top_node || top_node->get_num_nodes() == 2 ) {
                unsigned char op_precedence   = opPrecedence(f),
                              op_precedence_right = op_precedence;

                if (associativeOps.Find(f) == kNotFound) {
                    op_precedence_right ++;
                }
                if (op_level >= 0) { // need to worry about op's precedence
                    bool need_parens = op_precedence < op_level;

                    if (need_parens && top_node ) { // put parentheses around the return of this expression
                        result<<'(';
                        SubtreeToString (result, top_node->go_down(1),op_precedence,match_names, nil, mode);
                        result <<&this_node_op->GetCode();
                        SubtreeToString (result, top_node->go_down(2),op_precedence_right,match_names, nil, mode);
                        result<<')';
                        return;
                    }
                }
              
                if (top_node) {
                    SubtreeToString (result, top_node->go_down(1),op_precedence,match_names,nil, mode);
                }
                result << &this_node_op->GetCode();
                if (top_node) {
                    SubtreeToString (result, top_node->go_down(2),op_precedence_right,match_names,nil, mode);
                }
                return;
            } else { // mixed binary-unary operation
                result<<&this_node_op->GetCode();
                if (top_node) {
                    result<<'(';
                    SubtreeToString (result, top_node->go_down(1),opPrecedence(f),match_names,nil, mode);
                    result<<')';
                }
                return;
            }
        } else {
             if (this_node_op->TheCode() != HY_OP_CODE_MACCESS) {
                result<<&this_node_op->GetCode();
                if (top_node) {
                    result<<'(';
                    SubtreeToString (result, top_node->go_down(1),-1,match_names,nil, mode);
                    for (long k=2UL; k <= node_op_count; k++) {
                        result << ',';
                        SubtreeToString (result, top_node->go_down(k),-1,match_names,nil, mode);
                    }
                    result<<')';
                }
            } else { // matrix element access - treat specially
                if (top_node) {
                    SubtreeToString (result, top_node->go_down(1),-1,match_names,nil,mode);
                    for (long k=2; k<=node_op_count; k++) {
                        result<<'[';
                        SubtreeToString (result, top_node->go_down(k),-1,match_names,nil, mode);
                        result<<']';
                    }
                }
            }
        }
        return;
    }
    if (node_op_count < 0L) {
        // a user-defined function
        long func_id = this_node_op->UserFunctionID();
        result<< & GetBFFunctionNameByIndex(func_id);
        if (top_node) {
            result<<'(';
            long argument_count = GetBFFunctionArgumentCount(func_id);
            SubtreeToString (result, top_node->go_down(1),-1,match_names,nil,mode);
            for (long k=2L; k<=argument_count; k++) {
                result<<',';
                SubtreeToString (result, top_node->go_down(k),-1,match_names,nil,mode);
            }
            result<<')';
        }
        return;
    }
  
    HBLObjectRef op_data = this_node_op->GetANumber();
    if (op_data) {
        _String* conv = (_String*)op_data->toStr();
        if (op_data->ObjectClass()==STRING) {
            (result <<'"').AppendNewInstance (conv) << '"';
        } else {
            if (op_data->ObjectClass() == NUMBER && op_data->Value() < 0.0) {
              (result <<'(').AppendNewInstance (conv) << ')';
            } else {
                result.AppendNewInstance (conv);
            }
        }
    } else {
        result << "<null>";
    }
}

//__________________________________________________________________________________
bool     _Formula::IsEmpty(void) const {
    return theFormula.empty();
}

//__________________________________________________________________________________
hyFloat   _Formula::Newton(_Formula& derivative, _Variable* unknown, hyFloat target_value, hyFloat left, hyFloat right) {
    // find a root of the formulaic expression, using Newton's method, given the derivative and a bracketed root.
    // will alternate between bisections and Newton iterations based on what is fatser
    // check that there is indeed a sign change on the interval
  
  auto set_and_compute = [&] (hyFloat x) -> hyFloat {
    unknown->SetValue(x);
    return Compute()->Value()-target_value;
  };
  
  hyFloat    func_left, func_right, // keeps track of function values in the current interval, [left,right]
  root_guess, func_root_guess,
  lastCorrection = 100.,
  newCorrection;

  func_left = set_and_compute (left);
  if (func_left==0.0) {
    return left;
  }
  unknown->SetValue(right);
  func_right = set_and_compute (right);
  if (func_right==0.0) {
    return right;
  }

  if (func_left*func_right>0.0) { // bracket fail
    ReportWarning (_String((_String*)toStr(kFormulaStringConversionReportRanges)) & " = "&_String(target_value)&" has no (or multiple) roots in ["&_String(left)&",Inf); " & func_left & "-" & func_right);
    return    left;
  }
  // else all is good we can start the machine
  bool useNewton = false;

  root_guess = (right+left) * .5;

  for (unsigned long iterCount  = 0L; fabs(right-left)/MAX(left,right) > kMachineEpsilon*10. && iterCount < 200UL; iterCount++) {
    func_root_guess = set_and_compute (root_guess);
    if (func_root_guess == 0.) {
      return root_guess;
    }
    // get the correction term from the derivative
    hyFloat df_dx = derivative.Compute()->Value(),
    adjusted_root_guess;

    useNewton = true;
    if (df_dx==0.0) {
      useNewton = false;
    } else {
      newCorrection = -func_root_guess/df_dx;

      if (fabs(newCorrection/func_root_guess)<kMachineEpsilon*2. || fabs(newCorrection)<kMachineEpsilon*2.) { // correction too small - the root has been found
          return root_guess;
      }

      if (fabs(newCorrection/lastCorrection)>4.) { // Newton correction too large - revert to bisection
        useNewton = false;
      }

      adjusted_root_guess = root_guess +newCorrection;
      if (adjusted_root_guess<=left || adjusted_root_guess >=right) {
        useNewton = false;
      } else {
        lastCorrection = newCorrection;
      }
    }

    if (useNewton) {
      root_guess = adjusted_root_guess;
    } else {
      if (func_root_guess==0.0) {
        return root_guess;
      }
      if (func_root_guess*func_left > 0.0) { // move to the right half
        func_left   = func_root_guess;
        left  = root_guess;
      } else {
        right = root_guess;
      }
      root_guess = (right+left) * .5;
    }
  }
  return root_guess;
}


//__________________________________________________________________________________
hyFloat   _Formula::Brent(_Variable* unknown, hyFloat a, hyFloat b, hyFloat tol, _List* store_evals, hyFloat rhs) {
// find a root of the formulaic expression, using Brent's method and a bracketed root.
    // check that there is indeed a sign change on the interval

    auto set_and_compute = [&] (hyFloat x) -> hyFloat {
      unknown->SetValue(x);
      hyFloat fx = Compute()->Value()-rhs;
      
      if (store_evals) {
        store_evals->AppendNewInstance(new _Constant (x));
        store_evals->AppendNewInstance(new _Constant (fx));
      }
      
      return fx;
    };

    hyFloat  fa = 0.0,fb = 0.0,fc,d = b-a,e = b-a ,min1,min2,xm,p,q,r,s,tol1,
                c = b;

    min1 = unknown->GetLowerBound();
    min2 = unknown->GetUpperBound();

    long        it = 0;

    if (a>b) { // autobracket to the left
        fb = set_and_compute (b);
      
        if (b<0.00001 && b>-0.00001) {
            a = b-0.0001;
        } else {
            a = b-fabs(b)*0.1;
        }

        if (a<min1) {
            a = min1+0.5*(b-min1);
        }
      
        fa = set_and_compute (a);
      
        for (long k=0L; k<50L && fb*fa >= 0.0; k++) {
 
            d  = (b-a)*GOLDEN_RATIO;
            b  = a;
            a -= d;

            if (a<min1) {
                if (b>min1) {
                    a = min1;
                } else {
                    break;
                }
            }

            fb = fa;
            fa = set_and_compute (a);
        }
    } else if (CheckEqual (a,b)) { // autobracket to the right
        fa = set_and_compute (a);
      
        a = b;

        if ((b<0.00001)&&(b>-0.00001)) {
            b = b+0.0001;
        } else {
            b = b+fabs(b)*0.1;
        }

        if (b>min2) {
            b = a+0.5*(min2-a);
        }

        fb = set_and_compute (b);
      
        for (long k=0L; k<50L && fb*fa >= 0.0; k++) {

            d  = (b-a)*GOLDEN_RATIO;
            a  = b;
            b += d;

            if (b>min2) {
                if (a<min2) {
                    b = min2;
                } else {
                    break;
                }
            }

            fa = fb;
            fb = set_and_compute (b);
        }
    }


    if (fa == 0.0) {
        fa = set_and_compute (a);
        if (fa == 0.0) {
            return a;
        }
    }

    if (fb == 0.0) {
        fb = set_and_compute (b);
        if (fb == 0.0) {
            return b;
        }
    }

    if (fa*fb<0.0) {
        fc = fb;
        c = b;

        for (it = 0; it < MAX_BRENT_ITERATES; it++) {
            if (fb*fc>0.0) {
                fc = fa;
                c  = a;
                e = d = b-a;
            }

            if (fabs (fc) < fabs (fb)) {
                a     = b;
                b     = c;
                c     = a;
                fa    = fb;
                fb    = fc;
                fc    = fa;
            }

            tol1 = 2.*fabs(b)*kMachineEpsilon+.5*tol;

            xm = .5*(c-b);

            if (fabs(xm)<=tol1 || fb == 0.0) {
                return b;
            }

            if (fabs(e)>=tol1 && fabs (fa) > fabs (fb)) {
                s = fb/fa;
                if (a==c) {
                    p = 2.*xm*s;
                    q = 1.-s;
                } else {
                    q = fa/fc;
                    r = fb/fc;
                    p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.0));
                    q = (q-1.)*(r-1.)*(s-1.);
                }
                if (p>0.0) {
                    q = -q;
                }

                if (p<0.0) {
                    p = -p;
                }

                min1 = 3.*xm*q-fabs(tol1*q);
                min2 = fabs (e*q);
                if (2.*p < (min1<min2?min1:min2)) {
                    e = d;
                    d = p/q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = b;
            fa = fb;
            if (fabs(d)>tol1) {
                b+=d;
            } else {
                if (xm > 0.) {
                    b += fabs (tol1);
                } else {
                    b -= fabs (tol1);
                }
            }
            fb = set_and_compute (b);
        }
    }


    /*for (long i = 0; i < theFormula.lLength; i++) {
      _Operation *op_i = GetIthTerm(i);
      printf ("%ld: %s\n", i+1, _String((_String*)op_i->toStr()).sData);
    }*/

   _String msg ((_String*)toStr(kFormulaStringConversionSubstiteValues));
    msg = msg & "=" & rhs;
    if (it < MAX_BRENT_ITERATES) {
        msg =   msg & " has no (or multiple) roots in ["&_String(a)&","&_String(b)&"]";
    } else {
        msg =   msg & " failed to converge to sufficient precision in " & MAX_BRENT_ITERATES &" iterates.";
    }
    ReportWarning (msg);
    return    b;
}
//__________________________________________________________________________________
hyFloat   _Formula::Newton(_Formula& derivative, hyFloat target_value, hyFloat left, hyFloat max_right, _Variable* unknown) {
// given a monotone function and a left bracket bound, found the right bracket bound and solve
  
    // check that there is indeed a sign change on the interval
    unknown->SetValue(left);
    hyFloat  t1 = Compute()->Value(), right = left, t2, step = 1.0;

    if (max_right-left < step * 100.) {
        step = (max_right-left) * 0.01;
    }
    if  (step==0.0) {
        return left;
    }
    do {
        right += step;
        if (right>max_right) { // function doesn't seem to have a root
            ReportWarning (_String((_String*)toStr(kFormulaStringConversionSubstiteValues))&"="&_String(target_value)&" has no (or multiple) roots in ["&_String(left)&","&_String(right)&")");
            return    left;
        }
        unknown->SetValue(right);
        t2 = Compute()->Value();
        step*=2;
        if (right+step>max_right)
            if (right<max_right) {
                step = max_right-right;
            }
    } while ((target_value-t1)*(target_value-t2)>0);
    return Newton (derivative, unknown, target_value, left, right);

}

//__________________________________________________________________________________
hyFloat   _Formula::Newton(_Variable* unknown, hyFloat target_value, hyFloat x_min, hyFloat left, hyFloat right) {
    // check that there is indeed a sign change on the interval
    hyFloat  t1,t2,t3,t4,t5,lastCorrection = 100, newCorrection;
    _String     msg;
    t1 =Integral(unknown, x_min, left)-target_value;
    if (t1==0.0) {
        return left;
    }
    t2 = t1+Integral(unknown, left, right);
    if (t2==0.0) {
        return right;
    }
    if (t1*t2>0.0) {
      ReportWarning (_String((_String*)toStr(kFormulaStringConversionSubstiteValues))&"="&_String(target_value)&" has no (or multiple) roots in ["&_String(left)&","&_String(right)&")");
      return    left;
    }
    // else all is good we can start the machine
    bool useNewton = false;

    t3 = (right + left) * 0.5;

    while (right-left>1e-6) { // stuff to do
        if (!useNewton) {
            t3 = (right+left) * 0.5;
        }
        unknown->SetValue(t3);
        t4 = Integral(unknown, x_min, t3)-target_value;
        // get the correction term from the derivative
        unknown->SetValue(t3);
         t5 = Compute()->Value();
        useNewton = true;
        if (t5==0.0) {
            useNewton = false;
        } else {
            newCorrection = -t4/t5;
            if (fabs(newCorrection)<1e-5) { // correction too small - the root has been found
                return t3;
            }
            if (fabs(newCorrection/lastCorrection)>4) { // Newton correction too large - revert to bisection
                useNewton = false;
            }
            t5 = t3+newCorrection;
            if ((t5<=left)||(t5>=right)) {
                useNewton = false;
            } else {
                lastCorrection = newCorrection;
            }
        }
        if (useNewton) {
            t3 = t5;
        } else {
            t4 = Integral(unknown, x_min, t3)-target_value;
            if (t4==0.0) {
                return t3;
            }
            if (t4*t1 >0) {
                t1 = t4;
                left = t3;
            } else {
                right = t3;
            }
        }

    }
    return t3;
}

//__________________________________________________________________________________
hyFloat   _Formula::Newton( _Variable* unknown, hyFloat target_value,hyFloat x_min, hyFloat left) {
// given a monotone function and a left bracket bound, found the right bracket bound and solve
    // check that there is indeed a sign change on the interval
    hyFloat  t1 = Integral(unknown, x_min, left), right = left, t2, step = 1.0;
    do {
        right += step;
        t2 = Integral(unknown, right-step, right);
        step*=2;
        if (right>=1e10) { // function doesn't seem to have a root
            ReportWarning (_String((_String*)toStr(kFormulaStringConversionSubstiteValues))&"="&_String(target_value)&" has no (or multiple) roots in ["&_String(left)&",Inf)");
            return    0.0;
        }
    } while ((target_value-t1)*(target_value-t2-t1)>=0);
    return Newton ( unknown, target_value,x_min, left, right);

}

//__________________________________________________________________________________

hyFloat   _Formula::Integral(_Variable* dx, hyFloat left, hyFloat right, bool infinite) {
// uses Romberg's intergation method
    if (infinite) { // tweak "right" here
        hyFloat value = 1.0, step = 1.0, right1= -1.;
        right = left;
        while (value>1e-8) {
            right+=step;
            dx->SetValue(right);
            value = fabs(Compute()->Value());
            if ( value < 1e-4 && right1 < 0.) { // SLKP 20181205 why is this here?
                right1 = right;
            }
            step *= 2;
            if (step > 100000) { // integrand decreasing too slowly
                HandleApplicationError (_String(*(_String*)toStr(kFormulaStringConversionNormal)).Enquote() & " decreases too slowly to be integrated over an infinite interval");
                return 0.0;
            }
        }
        
        
        if (right1<right-kMachineEpsilon) {
            return Integral(dx,left,right1,false)+Integral(dx,right1,right,false);
        } else {
            return Integral(dx,left,right1,false);
        }
    }

    hyFloat          precision_factor =  hy_env::EnvVariableGetNumber(hy_env::integration_precision_factor);
    long             max_iterations  =  hy_env::EnvVariableGetNumber(hy_env::integration_maximum_iterations);

    hyFloat ss,
               dss,
               *s = new hyFloat [max_iterations],
               *h = new hyFloat [max_iterations+1L];

 
    h[0]=1.0;

    long         interpolate_steps = 5L,
                 stack_depth = 0L;

    _SimpleList  fvidx_aux,
                 changing_vars,
                 idx_map;
  
  
    _AVLList    fvidx (&fvidx_aux);

    hyFloat   * ic = new hyFloat[interpolate_steps],
              * id = new hyFloat[interpolate_steps];

    _SimpleFormulaDatum
      * stack = nil,
      * vvals = nil;


    if (AmISimple (stack_depth,fvidx)) {
        stack = new _SimpleFormulaDatum [stack_depth];
        
        ConvertToSimple (fvidx);
        if (fvidx.countitems()) { // could be a constant expression with no variable dependancies
            vvals = new _SimpleFormulaDatum [fvidx.countitems()];
            //fvidx.ReorderList();
            long dx_index = dx->get_index();
            for (long vi = 0; vi < fvidx_aux.countitems(); vi++) {
                _Variable* variable_in_expression = LocateVar (fvidx_aux.Element(vi));
                if (variable_in_expression->CheckFForDependence (dx_index,true)) {
                    changing_vars << fvidx_aux.Element(vi);
                    idx_map << vi;
                }
                vvals[vi].value = variable_in_expression->Compute()->Value();
            }
            
            long depends_on_dx = fvidx.FindLong(dx_index);
            
            if (depends_on_dx >= 0) {
                changing_vars.InsertElement ((BaseRef)dx_index,0,false,false);
                idx_map.InsertElement ((BaseRef)fvidx.FindLong(dx_index),0,false,false);
            }
        } else {
            vvals = new _SimpleFormulaDatum [1L];
        }
    } else {
        stack_depth = -1L;
    }

    for (long step = 0L; step < max_iterations; step ++) {
      //printf ("%d\n", step);
      if (stack_depth >=0) { // compiled
          s[step] = TrapezoidLevelKSimple(*this, dx, left, right, step+1L, stack, vvals,changing_vars,idx_map);
        } else {
            s[step] = TrapezoidLevelK(*this, dx, left, right, step+1L);
        }
        if (step >= 4L) {
            ss = InterpolateValue(&h[step-4L],&s[step-4L],interpolate_steps,ic,id,0.0, dss);
            if (fabs(dss)<= precision_factor * fabs(ss)) {
              
                BatchDeleteArray (s,h,ic,id);
                if (stack_depth >=0L ) {
                    ConvertFromSimple(fvidx);
                    BatchDeleteArray (stack, vvals);
                }
                return ss;
            }
        }
        h[step+1] =  h[step]/9.0;
    }

    if (stack_depth >=0) {
        ConvertFromSimple(fvidx);
        BatchDeleteArray (stack, vvals);
    }

    ReportWarning (_String("Integral of ")&*_String((_String*)toStr(kFormulaStringConversionNormal)) & " over ["&_String(left)&","&_String(right)&"] converges slowly, loss of precision may occur. Change either INTEGRATION_PRECISION_FACTOR or INTEGRATION_MAX_ITERATES");
    BatchDeleteArray (s,h,ic,id);
    return ss;
}
  //__________________________________________________________________________________
hyFloat  InterpolateValue (hyFloat* theX, hyFloat* theY, long n, hyFloat *c , hyFloat *d, hyFloat x, hyFloat& err) {
    // Neville's algoruthm for polynomial interpolation (Numerical Recipes' rawinterp)
  hyFloat y,
  den,
  dif = 1e10,
  dift,
  ho,
  hp,
  w;
  
  long   ns = 0L;
  
  for (long i=0L; i<n; i++) {
    dift = fabs(x-theX[i]);
    if (dift<dif) {
      ns = i;
      dif = dift;
    }
    c[i] = d[i] = theY[i];
  }
  
  y = theY[ns--];
  
  for (long m=1L; m<n; m++) {
    for (long i=0L; i<=n-m-1; i++) {
      ho = theX[i]-x;
      hp = theX[i+m]-x;
      w = c[i+1]-d[i];
      den = w/(ho-hp);
      d[i] = hp*den;
      c[i] = ho*den;
    }
    err = 2*ns< (n-m)? c[ns+1]: d[ns--];
    y += err;
  }
  
  return y;
}



  //__________________________________________________________________________________

hyFloat  TrapezoidLevelKSimple (_Formula&f, _Variable* xvar, hyFloat left, hyFloat right, long k, _SimpleFormulaDatum * stack, _SimpleFormulaDatum* values, _SimpleList& changingVars, _SimpleList& varToStack) {
    
    hyFloat x,
    tnm,
    sum,
    del,
    ddel;
    
    static hyFloat s;
    
    //_Constant dummy;
    
    auto set_and_compute = [&] (hyFloat x) -> hyFloat {
        if (changingVars.countitems() == 1L) {
            values[varToStack.get(0)].value = x;
        } else {
            xvar->SetNumericValue  (x);
            changingVars.Each ([&] (long v, unsigned long vi) -> void {
                values[varToStack.get(vi)].value = LocateVar(v)->Compute()->Value();
            });
        }
        return f.ComputeSimple(stack, values);
    };
    
    
    if (k==1) {
        s = set_and_compute ((left+right)*0.5);
        return s;
    }
    
    long        it =1;
    
    for (long j=1L; j<k-1; j++) {
        it *= 3L;
    }
    
    tnm = it;
    del = (right-left)/(3.0*tnm);
    ddel = del+del;
    x   = left+del*.5;
    sum = 0.;
    for (long j=1L; j<=it; j++, x+=del) {
        sum += set_and_compute (x);
        x+=ddel;
        sum += set_and_compute (x);
    }
    s = (s+(right-left)*sum/tnm)/3.0;
    return s;
}

  //__________________________________________________________________________________
hyFloat  TrapezoidLevelK (_Formula&f, _Variable* xvar, hyFloat left, hyFloat right, long k) {
    hyFloat x,
    tnm,
    sum,
    del,
    ddel;
    
    static hyFloat s;
    
    long        it = 1L;
    
    if (k==1L) {
        xvar -> SetNumericValue((left + right) * 0.5);
        return s = f.Compute()->Value();
    }
    
    for (long j = 1L; j < k-1; j++) {
        it *= 3L;
    }
    
    
    tnm = it;
    del = (right-left)/(3.0*tnm);
    ddel = del+del;
    x   = left+del*.5;
    
    sum = 0.0;
    
    for (long j=1L; j<=it; j++, x+=del) {
        xvar->SetNumericValue(x);
        sum += f.Compute()->Value();
        x+=ddel;
        xvar->SetNumericValue(x);
        sum += f.Compute()->Value();
    }
    s = (s+(right-left)*sum/tnm)/3.0;
    return s;
}

//__________________________________________________________________________________
hyFloat   _Formula::MeanIntegral(_Variable* dx, hyFloat left, hyFloat right, bool infinite) {
    _Formula newF;
    newF.Duplicate(this);
    newF.theFormula < new _Operation (true, *(dx->GetName())) < new _Operation (HY_OP_CODE_MUL, 2);
    return newF.Integral (dx,left,right, infinite);
}

//__________________________________________________________________________________
long     _Formula::NumberOperations(void) const {
// number of operations in the formula
    return theFormula.countitems();
}

//__________________________________________________________________________________

long      _Formula::ExtractMatrixExpArguments (_List* storage) {
  long count = 0L;

    if (resultCache && resultCache->empty() == false) {
        long cacheID      = 0L;
        bool cacheUpdated = false;
            // whether or not a cached result was used


        for (unsigned long i=0UL; i + 1UL <theFormula.countitems(); i++) {
            _Operation* this_op = GetIthTerm(i),
                      * next_op = GetIthTerm(i + 1UL);

              if (! cacheUpdated && next_op->CanResultsBeCached(this_op)) {
                  /*
                   StringToConsole("\n----\n");
                  ObjectToConsole(FetchVar(LocateVarByName("k"))->Compute());
                  NLToConsole();*/

                  _Stack temp;
                  this_op->Execute (temp);
                  

                  _Matrix *currentArg  = (_Matrix*)temp.Pop(true),
                          *cachedArg   = (_Matrix*)((HBLObjectRef)(*resultCache)(cacheID)),
                          *diff        = nil;

                  if (cachedArg->ObjectClass() == MATRIX) {
                      diff =  (_Matrix*)cachedArg->SubObj(currentArg, nil);
                  }

                  if (diff && diff->MaxElement() <= 1e-12) {
                      cacheID += 2;
                      i ++;
                  } else {
                      cacheUpdated = true;
                      cacheID++;
                      if (next_op->CanResultsBeCached(this_op, true)) {
                          storage->AppendNewInstance(currentArg);
                          count ++;
                      }
                  }
                  DeleteObject (diff);
                  continue;
              }
              if (cacheUpdated) {
                  cacheID++;
                  cacheUpdated = false;
              }
        }
    }

    return count;
}

//__________________________________________________________________________________
_Variable * _Formula::Dereference (bool ignore_context, _hyExecutionContext* theContext) {
    _Variable * result = nil;
    HBLObjectRef computedValue = Compute (0, theContext->GetContext(), nil, theContext->GetErrorBuffer());
    if (computedValue && computedValue->ObjectClass() == STRING) {
        result =  (_Variable*)((_FString*)computedValue)->Dereference(ignore_context, theContext, true);
    }

    if (!result) {
        theContext->ReportError( (_String ("Failed to dereference '") & _String ((_String*)toStr(kFormulaStringConversionNormal)) & "' in the " & (ignore_context ? "global" : "local") & " context"));
    }

    return result;
}

  //unsigned long ticker = 0UL;

//__________________________________________________________________________________
HBLObjectRef _Formula::Compute (long startAt, _VariableContainer const * nameSpace, _List* additionalCacheArguments, _String* errMsg, long valid_type)
// compute the value of the formula
// TODO SLKP 20170925 Needs code review
{
    _Stack * scrap_here;
    if (theFormula.empty()) {
        theStack.theStack.Clear();
        theStack.Push (new _MathObject, false);
        scrap_here = &theStack;
    } else {
        bool wellDone = true;


        if (call_count++) {
          scrap_here = new _Stack;
        } else {
          scrap_here = &theStack;
          if (startAt == 0) {
              theStack.Reset();
          }
        }

        /*ticker++;
        if (ticker >= 1462440) {
          printf ("\n_Formula::Compute (%x, %d)  %ld terms, stack depth %ld\n", this, ticker, theFormula.lLength, theStack.theStack.lLength);
        }*/
      
        const unsigned long term_count = NumberOperations();

        if (startAt == 0L && resultCache && !resultCache->empty()) {
            long cacheID     = 0L;
                // where in the cache are we currently looking
            bool cacheUpdated = false;
                // whether or not a cached result was used

            for (unsigned long i=0; i< term_count ; i++) {
                _Operation* thisOp = ItemAt (i);
                if ( i + 1UL < term_count ) {
                    _Operation* nextOp  = ItemAt (i+1UL);

                    if (! cacheUpdated && nextOp->CanResultsBeCached(thisOp)) {
                        if (!thisOp->Execute(*scrap_here,nameSpace, errMsg)) {
                            wellDone = false;
                            break;
                        }

                        _Matrix *currentArg  = (_Matrix*)scrap_here->Pop(false),
                                *cachedArg   = (_Matrix*)((HBLObjectRef)(*resultCache)(cacheID)),
                                *diff        = nil;

                        if (cachedArg->ObjectClass() == MATRIX) {
                            diff =  (_Matrix*)cachedArg->SubObj(currentArg, nil);
                        }

                        bool    no_difference = diff && diff->MaxElement() <= 1e-12;

                        if (no_difference || (additionalCacheArguments && !additionalCacheArguments->empty() && nextOp->CanResultsBeCached(thisOp,true))) {
                            DeleteObject  (scrap_here->Pop  ());
                            if (no_difference) {
                                scrap_here->Push ((HBLObjectRef)(*resultCache)(cacheID+1));
                            } else {

                                scrap_here->Push ((HBLObjectRef)additionalCacheArguments->GetItem (0));
                                resultCache->Replace(cacheID,scrap_here->Pop(false),true);
                                resultCache->Replace(cacheID+1,(HBLObjectRef)additionalCacheArguments->GetItem (0),false);
                                additionalCacheArguments->Delete (0, false);
                                //printf ("_Formula::Compute additional arguments %ld\n", additionalCacheArguments->lLength);
                           }
                            cacheID += 2;
                            i ++;
                            //printf ("Used formula cache %s\n", _String((_String*)nextOp->toStr()).sData);
                        } else {
                            cacheUpdated = true;
                            resultCache->Replace(cacheID++,scrap_here->Pop(false),true);
                            //printf ("Updated formula cache %s\n", _String((_String*)nextOp->toStr()).sData);
                       }
                       DeleteObject (diff);
                       continue;
                    }
                }
                if (!thisOp->Execute(*scrap_here,nameSpace, errMsg)) { // does this always get executed?
                    wellDone = false;
                    break;
                }
                if (cacheUpdated) {
                    resultCache->Replace(cacheID++,scrap_here->Pop(false),true);
                    cacheUpdated = false;
                }
            }
        } else {

            for (unsigned long i=startAt; i< term_count; i++) {
                  if (!ItemAt (i)->Execute(*scrap_here, nameSpace, errMsg, call_count == 1)) {
                      wellDone = false;
                      break;
                  }

            }
        }
        if (scrap_here->StackDepth() != 1L || !wellDone) {
            _String errorText = _String ("'") & _String((_String*)toStr(kFormulaStringConversionNormal)) & _String("' evaluated with errors ");
            if (errMsg && errMsg->nonempty()) {
                errorText = errorText & " " & errMsg->Enquote('(',')');
            }

            if (scrap_here->StackDepth() > 1 && wellDone) {
                errorText = errorText & " Unconsumed values on the stack";
                for (long stack_id = scrap_here->StackDepth()-1; stack_id >= 0; stack_id --) {
                  errorText = errorText & "\n[" & (stack_id+1) & "]------------------\n" & (_String*) scrap_here->Pop(false)->toStr();
                }
                errorText & "\n------------------\n";
            }

            if (errMsg) {
                *errMsg = *errMsg & errorText;
            }
            else {
                HandleApplicationError (errorText);
            }
            scrap_here->Reset();
            scrap_here->Push (new _Constant (0.0), false);
         }
    }


    HBLObjectRef return_value = scrap_here->Pop(false);

    if (theFormula.lLength) {
         DeleteObject (recursion_calls);
         if (--call_count) {
          recursion_calls = return_value;
          return_value->AddAReference();
          delete scrap_here;

        } else {
          recursion_calls = nil;
        }
    }
    return valid_type == HY_ANY_OBJECT ? return_value : ((return_value->ObjectClass() & valid_type) ? return_value : nil);
}

//__________________________________________________________________________________
bool _Formula::CheckSimpleTerm (HBLObjectRef thisObj) {
    if (thisObj) {
        long oc = thisObj->ObjectClass();
        if (oc != NUMBER) {
            if (oc ==MATRIX) {
                _Matrix * mv = (_Matrix*)thisObj->Compute();
                if (!mv->SparseDataStructure() && mv->IsIndependent ()) {
                    return true;
                }
            }
        } else {
            return true;
        }
    }
    return false;
}

//__________________________________________________________________________________
_Formula const _Formula::operator+ (const _Formula& operand2) {
    return *PatchFormulasTogether (*this,operand2,HY_OP_CODE_ADD);
}

//__________________________________________________________________________________
_Formula const _Formula::operator- (const _Formula& operand2) {
    return *PatchFormulasTogether (*this,operand2,HY_OP_CODE_SUB);
}

//__________________________________________________________________________________
_Formula const _Formula::operator* (const _Formula& operand2) {
    return *PatchFormulasTogether (*this,operand2,HY_OP_CODE_MUL);
}

//__________________________________________________________________________________
_Formula const _Formula::operator/ (const _Formula& operand2) {
    return *PatchFormulasTogether (*this,operand2,HY_OP_CODE_DIV);
}

//__________________________________________________________________________________
_Formula const _Formula::operator^ (const _Formula& operand2) {
    return *PatchFormulasTogether (*this,operand2,HY_OP_CODE_POWER);
}

//__________________________________________________________________________________
_Formula* _Formula::PatchFormulasTogether (const _Formula& op1, const _Formula& op2, const char op_code) {
    _Formula * result = new _Formula;
    result->DuplicateReference(&op1);
    result->DuplicateReference(&op2);
    result->theFormula.AppendNewInstance(new _Operation (op_code, 2));
    return result;
}

//__________________________________________________________________________________
_Formula* _Formula::PatchFormulasTogether (const _Formula& op1, HBLObjectRef op2, const char op_code) {
    _Formula * result = new _Formula;
    result->DuplicateReference(&op1);
    result->theFormula.AppendNewInstance(new _Operation (op2));
    result->theFormula.AppendNewInstance(new _Operation (op_code, 2));
    return result;
}


//__________________________________________________________________________________
void _Formula::ConvertMatrixArgumentsToSimpleOrComplexForm (bool make_complex) {
  unsigned long term_count = NumberOperations();
  
  if (make_complex) {
    if (resultCache) {
      DeleteObject (resultCache);
      resultCache = nil;
    }
  } else {
    if (!resultCache) {
      resultCache = new _List();
      for (unsigned long i = 1UL ; i< term_count ; i++) {
        if ( ItemAt(i)->CanResultsBeCached(ItemAt (i-1))) {
          resultCache->AppendNewInstance(new _MathObject());
          resultCache->AppendNewInstance(new _MathObject());
        }
      }
    }
  }

  for (unsigned long i = 0UL ; i< term_count ; i++) {
    _Operation* this_op = ItemAt (i);
    _Matrix   * op_matrix = nil;

    if (this_op->theNumber) {
      if (this_op->theNumber->ObjectClass() == MATRIX) {
        op_matrix = (_Matrix*)this_op->theNumber;
      }
    } else {
      if (this_op->theData >= 0L) {
        _Variable* thisVar = LocateVar (this_op->theData);
        if (thisVar->ObjectClass() == MATRIX) {
          op_matrix = (_Matrix*)thisVar->GetValue();
        }
      }
    }

    if (op_matrix) {
      if (make_complex) {
        op_matrix->MakeMeGeneral();
      } else {
        op_matrix->MakeMeSimple();
      }
    }
  }
}

//__________________________________________________________________________________
long _Formula::StackDepth (long from, long to) const {
  _SimpleList::NormalizeCoordinates(from, to, NumberOperations());
  long result = 0L;

  for ( long i = from; i <= to; i++) {
     result += ItemAt(i)->StackDepth ();
  }

  return result;

}


//__________________________________________________________________________________
bool _Formula::AmISimple (long& stack_depth, _AVLList& variable_index) {
    if (theFormula.empty()) {
        return true;
    }

    long loc_depth = 0L;

    for (unsigned long i=0UL; i<theFormula.countitems(); i++) {
        _Operation* this_op =ItemAt (i);
        loc_depth++;
        if ( this_op->theData<-2 || this_op->numberOfTerms<0) {
            if (this_op->theData < -2 && i == 0UL) {
              variable_index.InsertNumber (this_op->GetAVariable());
              continue;
            }
            return false;
        }

        if (this_op->theNumber) {
            if (this_op->theNumber->ObjectClass() != NUMBER) {
                return false;
            }
        } else {
            if (this_op->theData >= 0) {
                _Variable* this_var = LocateVar (this_op->theData);
                if (this_var->ObjectClass()!=NUMBER) {
                    HBLObjectRef cv = this_var->GetValue();
                    if (!CheckSimpleTerm (cv)) {
                        return false;
                    }
                }
                variable_index.InsertNumber (this_op->GetAVariable());
            } else {
                long op_code = this_op->TheCode();
              
                if (simpleOperationCodes.Find(op_code)==kNotFound) {
                    return false;
                } else if ((op_code == HY_OP_CODE_MACCESS || op_code == HY_OP_CODE_MCOORD || op_code == HY_OP_CODE_MUL) && this_op->GetNoTerms() != 2) {
                    return false;
                }

                loc_depth -= this_op->GetNoTerms();
            }
        }
        if (loc_depth>stack_depth) {
            stack_depth = loc_depth;
        } else if (loc_depth==0L) {
             HandleApplicationError (_String("Invalid formula (no return value) passed to ") & __PRETTY_FUNCTION__ & " :" & _String ((_String*)toStr(kFormulaStringConversionNormal)).Enquote());
            return false;
        }
    }
    return true;
}

//__________________________________________________________________________________
bool _Formula::IsArrayAccess (void){
    if (!theFormula.empty()) {
        return (ItemAt (theFormula.countitems()-1)->TheCode() == HY_OP_CODE_MACCESS);
    }
    return false;
}

//__________________________________________________________________________________
bool _Formula::ConvertToSimple (_AVLList& variable_index) {
    bool has_volatiles = false;
    if (!theFormula.empty()) {
        for (unsigned long i=0UL; i<theFormula.countitems(); i++) {
            _Operation* this_op = ItemAt (i);
            if (this_op->theNumber) {
                continue;
            } else if (this_op->theData >= 0) {
                this_op->theData = variable_index.FindLong (this_op->theData);
            } else if (this_op->opCode == HY_OP_CODE_SUB && this_op->numberOfTerms == 1) {
                this_op->opCode = (long)MinusNumber;
            } else {
                if (this_op->opCode == HY_OP_CODE_MACCESS) {
                    this_op->numberOfTerms = -2;
                } else {
                  if (this_op->opCode == HY_OP_CODE_MCOORD) {
                    this_op->numberOfTerms = -3;
                  }
                }
                if (this_op->opCode == HY_OP_CODE_RANDOM || this_op->opCode == HY_OP_CODE_TIME)
                    has_volatiles = true;
                this_op->opCode = simpleOperationFunctions(simpleOperationCodes.Find(this_op->opCode));
            }

        }
    }
    return has_volatiles;
}

//__________________________________________________________________________________
void _Formula::ConvertFromSimple (_AVLList const& variableIndex) {
  ConvertFromSimpleList (*variableIndex.dataList);
}

//__________________________________________________________________________________
void _Formula::ConvertFromSimpleList (_SimpleList const& variableIndex) {
  if (theFormula.empty()) {
    return;
  }
  
  for (unsigned long i=0UL; i<theFormula.countitems(); i++) {
    _Operation* this_op = ItemAt (i);
    if (this_op->theNumber) {
      continue;
    } else {
      if (this_op->theData>-1) {
        this_op->theData = variableIndex.get (this_op->theData);
      } else if (this_op->opCode == (long)MinusNumber) {
        this_op->opCode = HY_OP_CODE_SUB;
      } else {
        if (this_op->opCode == (long)FastMxAccess) {
          this_op->numberOfTerms = 2;
        } else {
          if (this_op->opCode == (long)FastMxWrite) {
            this_op->numberOfTerms = 2;
          }
        }
        this_op->opCode = simpleOperationCodes(simpleOperationFunctions.Find(this_op->opCode));
      }
    }
  }
}

//__________________________________________________________________________________
hyFloat _Formula::ComputeSimple (_SimpleFormulaDatum* stack, _SimpleFormulaDatum* varValues)
{
    if (theFormula.nonempty()) {
        long stackTop = 0;
        unsigned long upper_bound = NumberOperations();

        for (unsigned long i=0UL; i<upper_bound; i++) {
            _Operation const* thisOp = ItemAt (i);
            if (thisOp->theNumber) {
                stack[stackTop++].value = thisOp->theNumber->Value();
                continue;
            } else {
                if (thisOp->theData>-1) {
                    stack[stackTop++] = varValues[thisOp->theData];
                } else {
                    stackTop--;
                    if (thisOp->numberOfTerms==2) {
                        hyFloat  (*theFunc) (hyFloat, hyFloat);
                        theFunc = (hyFloat(*)(hyFloat,hyFloat))thisOp->opCode;
                        if (stackTop<1L) {
                            HandleApplicationError ("Internal error in _Formula::ComputeSimple - stack underflow.)", true);
                            return 0.0;
                        }
                        stack[stackTop-1].value = (*theFunc)(stack[stackTop-1].value,stack[stackTop].value);
                    } else {
                      switch (thisOp->numberOfTerms) {
                        case -2 : {
                            hyFloat  (*theFunc) (hyPointer,hyFloat);
                            theFunc = (hyFloat(*)(hyPointer,hyFloat))thisOp->opCode;
                            if (stackTop<1L) {
                                HandleApplicationError ("Internal error in _Formula::ComputeSimple - stack underflow.)", true);
                                return 0.0;
                            }
                            stack[stackTop-1].value = (*theFunc)(stack[stackTop-1].reference,stack[stackTop].value);
                            break;
                          }
                        case -3 : {
                          void  (*theFunc) (hyPointer,hyFloat,hyFloat);
                          theFunc = (void(*)(hyPointer,hyFloat,hyFloat))thisOp->opCode;
                          if (stackTop != 2 || i != theFormula.lLength - 1) {
                            HandleApplicationError ("Internal error in _Formula::ComputeSimple - stack underflow or MCoord command is not the last one.)", true);

                            return 0.0;
                          }
                          //stackTop = 0;
                          // value, reference, index
                          (*theFunc)(stack[1].reference,stack[2].value, stack[0].value);
                          break;
                        }
                        default: {
                            hyFloat  (*theFunc) (hyFloat);
                            theFunc = (hyFloat(*)(hyFloat))thisOp->opCode;
                            stack[stackTop].value = (*theFunc)(stack[stackTop].value);
                            ++stackTop;
                        }
                      }

                    }
                }
            }
        }
        return stack[0].value;
    }
    return 0.0;
}

//__________________________________________________________________________________
bool _Formula::EqualFormula (_Formula* f) {
    if (theFormula.countitems() == f->theFormula.countitems()) {
        unsigned long const upper_bound = NumberOperations();
      
        for (unsigned long i=0UL; i<upper_bound; i++) {
            if (ItemAt (i) -> EqualOp (f->ItemAt(i)) == false) {
                return false;
            }
        }
        return true;
    }
    return false;
}

//__________________________________________________________________________________
HBLObjectRef _Formula::ConstructPolynomial (void) {
    theStack.Reset();
    bool wellDone = true;
    _String errMsg;

    for (unsigned long i=0UL; wellDone && i<theFormula.countitems(); i++) {
      wellDone = ItemAt (i)->ExecutePolynomial(theStack, nil, &errMsg);
        /*if (wellDone) {
            printf ("%ld (%ld)\n", i, theStack.StackDepth());
            for (long k = 0; k < theStack.StackDepth(); k++) {
                printf ("\t%s\n", _String((_String*)theStack.Peek (k)->toStr()).get_str());
            }
        }*/
    }
    
    if (theStack.StackDepth() == 1 && wellDone) {
        return theStack.Pop(false);
    }

    return nil;
}

//__________________________________________________________________________________
bool _Formula::HasChanged (bool ingoreCats) {
    unsigned long const upper_bound = NumberOperations();
  
    for (unsigned long i=0UL; i<upper_bound; i++) {
        long data_id;
        _Operation * this_op = ItemAt (i);
      
        if (this_op->IsAVariable()) {
            data_id = this_op->GetAVariable();
            if (data_id>=0) {
                if (((_Variable*)(((BaseRef*)(variablePtrs.list_data))[data_id]))->HasChanged(ingoreCats)) {
                    return true;
                }
            } else if (this_op->theNumber->HasChanged()) {
                return true;
            }
        } else if (this_op->opCode == HY_OP_CODE_BRANCHLENGTH||this_op->opCode == HY_OP_CODE_RANDOM||this_op->opCode == HY_OP_CODE_TIME)
            // time, random or branch length
        {
            return true;
        } else if (this_op->numberOfTerms<0L) {
            data_id = -this_op->numberOfTerms-2;
            if (IsBFFunctionIndexValid (data_id)) {
                if (GetBFFunctionType (data_id) == kBLFunctionSkipUpdate) {
                    continue;
                }
            }
            return true;
        }
    }
    return false;
}


//__________________________________________________________________________________

void _Formula::ScanFormulaForHBLFunctions (_AVLListX& collection , bool recursive) {


  auto handle_function_id = [&collection, recursive] (const long hbl_id) -> void {
    if (IsBFFunctionIndexValid(hbl_id)) {
      _String function_name = GetBFFunctionNameByIndex(hbl_id);

      if (collection.Find(&function_name) < 0) {
        collection.Insert (new _String (function_name), HY_BL_HBL_FUNCTION, false, false);
        if (recursive) {
          GetBFFunctionBody(hbl_id).BuildListOfDependancies(collection, true);
        }
      }
    }
  };

  ConvertToTree(false);

  if (theTree) {

    InternalSimplify(theTree);
    node_iterator<long> ni (theTree, _HY_TREE_TRAVERSAL_PREORDER);

    while (node<long>* iterator = ni.Next()) {
      _Operation *this_op = GetIthTerm(iterator->get_data());

      long hbl_id = -1L;

      if (this_op -> IsHBLFunctionCall()) {
        hbl_id = this_op -> GetHBLFunctionID();
      } else {
        if (this_op->opCode == HY_OP_CODE_CALL) {
          node <long>* function_to_call = iterator->go_down (1);
          _Operation * function_to_call_value = GetIthTerm (function_to_call->get_data());

          if (function_to_call->get_num_nodes() == 0 && function_to_call_value -> IsConstantOfType(STRING)) {
              hbl_id = FindBFFunctionName (((_FString*)function_to_call_value->theNumber->Compute())->get_str());
          } else {
              ReportWarning ("Cannot export Call function arguments which are run-time dependent");
          }
        } else {
          if (this_op->opCode == HY_OP_CODE_MACCESS) { // handle AVL iterators
            if (this_op->GetNoTerms() == 3) { // [][]
              if (iterator->go_down(2)->get_num_nodes() == 0 && iterator->go_down(3)->get_num_nodes() == 0) {
                _Operation* bracket_1 = GetIthTerm (iterator->go_down(2)->get_data());
                _Operation* bracket_2 = GetIthTerm (iterator->go_down(3)->get_data());
                if (bracket_1->IsConstantOfType(STRING) && bracket_2->IsConstantOfType(STRING)) {
                  handle_function_id (FindBFFunctionName (((_FString*)bracket_1->theNumber->Compute())->get_str()));
                  handle_function_id (FindBFFunctionName (((_FString*)bracket_2->theNumber->Compute())->get_str()));
                  continue;
                }
              }
              ReportWarning ("Potentially missed dependence on a function in [][]; arguments are run-time dependent");
            }
          }
        }
      }

      handle_function_id (hbl_id);



    }
  }
}

//__________________________________________________________________________________
bool _Formula::HasChangedSimple (_SimpleList& variableIndex) {
    unsigned long const upper_bound = NumberOperations();
    
    for (unsigned long i=0UL; i<upper_bound; i++) {
      _Operation * this_op = ItemAt (i);
        if (this_op->theNumber) {
            continue;
        } else if (this_op->theData >= 0) {
            if (((_Variable*)(((BaseRef*)(variablePtrs.list_data))[variableIndex.list_data[this_op->theData]]))->HasChanged(false)) {
                return true;
            }
        } else {
            if (this_op->opCode == (long) RandomNumber) {
                return true;
            }
        }
    }
    return false;
}

//__________________________________________________________________________________
void _Formula::ScanFForVariables (_AVLList&l, bool includeGlobals, bool includeAll, bool includeCategs, bool skipMatrixAssignments, _AVLListX* tagger, long weight) const {
    unsigned long const upper_bound = NumberOperations();
  
    for (unsigned long i=0UL; i<upper_bound; i++) {
        _Operation * this_op = ItemAt (i);
        if (this_op->IsAVariable()) {
            if (!includeGlobals)
                // This change was part of a commit that introduced an optimizer bug (suspected location:
                // src/core/batchlan2.cpp:2220). This change is suspicious as well (removed and undocumented condition).
                //if ((((_Variable*)LocateVar(theObj->GetAVariable()))->IsGlobal())||
                 //       (((_Variable*)LocateVar(theObj->GetIndex()))->ObjectClass()!=NUMBER)) {
                if (((_Variable*)LocateVar(this_op->GetAVariable()))->IsGlobal()) {
                    continue;
                }

            long f = this_op->GetAVariable();

            if (f>=0) {
                _Variable * v = LocateVar(f);

                if (v->IsCategory()&&includeCategs) {
                    v->ScanForVariables (l,includeGlobals,tagger, weight);
                }

                if(includeAll || v->ObjectClass()==NUMBER) {
                    l.Insert ((BaseRef)f);
                    if (tagger) {
                        tagger -> UpdateValue((BaseRef)f, weight, 0);
                    }
                }

                if (skipMatrixAssignments) {
                    if (v->ObjectClass()!=MATRIX || !this_op->AssignmentVariable()) {
                        v->ScanForVariables(l,includeGlobals,tagger, weight);
                    }
                } else if (!v->IsIndependent()) {
                    v->ScanForVariables(l,includeGlobals,tagger);
                }
            } else if (this_op->theNumber)
                if (this_op->theNumber->ObjectClass()==MATRIX) {
                    ((_Matrix*)this_op->theNumber)->ScanForVariables(l,includeGlobals,tagger, weight);
                }
        }
    }
}

//__________________________________________________________________________________
void _Formula::ScanFForType (_SimpleList &l, int type) {
  unsigned long const upper_bound = NumberOperations();
  
  for (unsigned long i=0UL; i<upper_bound; i++) {
    _Operation * this_op = ItemAt (i);
    if (this_op->IsAVariable()) {
      long f = this_op->GetAVariable();
      
      if (f>=0) {
        _Variable * v = LocateVar(f);
        
        if(v->ObjectClass()==type) {
          l << f;
        }
        
      }
    }
  }
}

//__________________________________________________________________________________
bool _Formula::CheckFForDependence (long varID, bool checkAll) {
  unsigned long const upper_bound = NumberOperations();
  
  for (unsigned long i=0UL; i<upper_bound; i++) {
    _Operation * this_op = ItemAt (i);
    if (this_op->IsAVariable()) {
      long f = this_op->GetAVariable();
      if (f>=0) {
        if (f == varID) {
          return true;
        }
        if (checkAll) {
            if (LocateVar(f)->CheckFForDependence(varID)) {
              return true;
            }
        }
      }
    }
  }
  return false;
}

//__________________________________________________________________________________
bool _Formula::CheckFForDependence (_AVLList const& indices, bool checkAll) {
  unsigned long const upper_bound = NumberOperations();
  
  for (unsigned long i=0UL; i<upper_bound; i++) {
    _Operation * this_op = ItemAt (i);
    if (this_op->IsAVariable()) {
      long f = this_op->GetAVariable();
      if (f>=0) {
        if (indices.FindLong (f) >= 0) {
          return true;
        }
        if (checkAll) {
             if (LocateVar(f)->CheckFForDependence(indices)) {
              return true;
            }
        }
      }
    }
  }
  return false;
}


  //__________________________________________________________________________________
void  _Formula::LocalizeFormula (_Formula& ref, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv) {
  unsigned long const upper_bound = ref.NumberOperations();
  
  for (unsigned long i=0UL; i<upper_bound; i++) {
    _Operation * this_op = ref.ItemAt (i);
    if (this_op->IsAVariable()) {
      long     vIndex = this_op->GetAVariable();
      _Variable* theV = LocateVar (vIndex);
      if (theV->IsGlobal()) {
        theFormula&& ref.theFormula(i);
        continue;
      }
      if (theV->IsContainer()) {
        continue;
      }
      _String  localized_name = parentName&"."&*(theV->GetName());
      
      _Variable * localized_var = CheckReceptacle (&localized_name, kEmptyString, false, false);
      
      if (theV->IsIndependent()) {
        iv<<localized_var->get_index();
        iiv<<vIndex;
      } else {
        dv<<localized_var->get_index();
        idv<<vIndex;
        
      }
      theFormula.AppendNewInstance( new _Operation (true, localized_name));
    } else {
      theFormula&& ref.theFormula(i);
    }
  }
}

  //__________________________________________________________________________________
bool _Formula::DependsOnVariable (long idx) {
  unsigned long const upper_bound = NumberOperations();
  
  for (unsigned long i=0UL; i<upper_bound; i++) {
    _Operation * this_op = ItemAt (i);
    if (this_op->IsAVariable() && this_op->GetAVariable() == idx) {
      return true;
    }
  }
  
  return false;
}

//__________________________________________________________________________________
_Operation* _Formula::GetIthTerm (long idx) const {
    if (idx >= 0 && idx < theFormula.lLength) {
        return (_Operation*)theFormula.GetItem(idx);
    }

    return nil;
}

//__________________________________________________________________________________
_Operation* _Formula::ItemAt (long idx) const {
  return (_Operation*)theFormula.GetItem(idx);
}

  //__________________________________________________________________________________
bool _Formula::IsAConstant (void) {
  
  unsigned long const upper_bound = NumberOperations();
  
  for (unsigned long i=0UL; i<upper_bound; i++) {
    if ( ItemAt (i)->IsAVariable()) {
      return false;
    }    
  }
  return true;

}

//__________________________________________________________________________________
bool _Formula::IsConstant (bool strict) {
  unsigned long const upper_bound = NumberOperations();

  for (unsigned long i=0UL; i<upper_bound; i++) {
    if (ItemAt (i)->IsConstant(strict) == false) {
      return false;
    }
  }
    
  return true;
}

//__________________________________________________________________________________
void _Formula::PushTerm (BaseRef object) {
  _List * test_list;
  if ((test_list = dynamic_cast<_List*>(object))) {
    theFormula << *test_list;
  } else {
    theFormula << object;
  }
}

//__________________________________________________________________________________
void _Formula::SimplifyConstants (void){
  ConvertToTree    ();
  InternalSimplify (theTree);
  ConvertFromTree ();
}

//__________________________________________________________________________________
HBLObjectRef _Formula::GetTheMatrix (void) {
    if (theFormula.countitems()==1) {
        _Operation* firstOp = ItemAt (0);
        HBLObjectRef   ret = firstOp->GetANumber();
        if (ret&&(ret->ObjectClass()==MATRIX)) {
            return ret;
        } else {
            if (firstOp->theData!=-1) {
                _Variable* firstVar = LocateVar (firstOp->GetAVariable());
                ret = firstVar->GetValue();
                if (ret&&(ret->ObjectClass()==MATRIX)) {
                    return ret;
                }
            }
        }
    }
    return nil;
}

//__________________________________________________________________________________
long _Formula::ObjectClass (void) {
    if (theStack.StackDepth()) {
        return ((HBLObjectRef)theStack.theStack.list_data[0])->ObjectClass();
    }

    HBLObjectRef res =   Compute();

    if (res) {
        return res->ObjectClass();
    }

    return HY_UNDEFINED;
}



//__________________________________________________________________________________
long _Formula::ParseFormula (_String const &s, _VariableContainer const* theParent, _String* reportErrors) {
    theTree     = nil;
    resultCache = nil;
    recursion_calls = nil;
    call_count = 0UL;

    _FormulaParsingContext fpc (reportErrors, theParent);

    _String formula_copy (s);

    long    return_value = Parse (this, formula_copy, fpc, nil);

    if (return_value != HY_FORMULA_EXPRESSION) {
        Clear();
    }

    return return_value;

}

//__________________________________________________________________________________
void    _Formula::ConvertToTree (bool err_msg) {
    if (!theTree && theFormula.countitems()) { // work to do
        _SimpleList nodeStack;

        unsigned long const upper_bound = NumberOperations();
        theTree = nil;
        try {
        
            for (unsigned long i=0UL; i<upper_bound; i++) {

                _Operation* currentOp = ItemAt(i);

                if (currentOp->theNumber || currentOp->theData >= 0L || currentOp->theData <= -2L) { // a data bit
                    node<long>* leafNode = new node<long>;
                    leafNode->init(i);
                    nodeStack<<(long)leafNode;
                } else { // an operation
                    long nTerms = currentOp->GetNoTerms();
                    if (nTerms<0L) {
                        nTerms = GetBFFunctionArgumentCount(currentOp->opCode);
                    }

                    if (nTerms>nodeStack.lLength) {
                        throw (_String ("Insufficient number of arguments for a call to ") & _String ((_String*)currentOp->toStr()) & " while converting " & toRPN(kFormulaStringConversionNormal).Enquote() & " to a parse tree");
                    }

                    node<long>* operationNode = new node<long>;
                    operationNode->init(i);
                    for (long j=0; j<nTerms; j++) {
                        operationNode->prepend_node(*((node<long>*)nodeStack.Pop()));
                    }
                    nodeStack<<(long)operationNode;
                }
            }
            if (nodeStack.lLength!=1) {
                throw ((_String)"The expression '" & toRPN (kFormulaStringConversionNormal) & "' has " & (long)nodeStack.lLength & " terms left on the stack after evaluation");
            } else {
                theTree = (node<long>*)nodeStack(0);
            }
        } catch (_String const& e) {
            if (err_msg) {
                HandleApplicationError(e);
            }
            nodeStack.Each ([this] (long v, unsigned long) -> void {
                node<long>* operationNode = (node<long>*)v;
                operationNode->delete_tree();
                delete (operationNode);
            });
            theTree = nil;
        }
    }
}
//__________________________________________________________________________________
void    _Formula::ConvertFromTree (void) {
    if (theTree) { // work to do
        _SimpleList termOrder;
        node_iterator<long> ni (theTree, _HY_TREE_TRAVERSAL_POSTORDER);

      while (node<long>* iterator = ni.Next()) {
            termOrder<<iterator->get_data();
        }

        if (termOrder.lLength!=theFormula.lLength) { // something has changed
            _List newFormula;
            for (long i=0; i<termOrder.lLength; i++) {
                newFormula<<theFormula(termOrder(i));
            }
            theFormula.Clear();
            theFormula.Duplicate(&newFormula);
            //ConvertToTree();
        }
        theTree->delete_tree();
        delete (theTree);
        theTree = nil;
    }
}

