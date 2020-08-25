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

#include "fstring.h"
#include "batchlan.h"
#include "global_object_lists.h"
#include "global_things.h"
#include "calcnode.h"
#include "function_templates.h"
#include "tree.h"
#include "tree_iterator.h"

/*extern long lastMatrixDeclared;
extern _AVLListX _HY_GetStringGlobalTypes;

extern _List likeFuncList,
             batchLanguageFunctionNames,
             dataSetList,
             scfgList;

extern _SimpleList modelMatrixIndices;
extern _String last_model_parameter_list;*/

using namespace hyphy_global_objects;
using namespace hy_global;

//__________________________________________________________________________________
_FString::_FString (void) {
    the_string = new _StringBuffer;
}

//__________________________________________________________________________________
_FString::~_FString () {
    if (CanFreeMe()) {
        DeleteObject (the_string);
    } else {
        RemoveAReference();
    }
}

//__________________________________________________________________________________
_FString::_FString (_String* data) {
    the_string = new _StringBuffer (data);
}

//__________________________________________________________________________________

_FString::_FString (const char* data) {
    the_string = new _StringBuffer (data);
}

//__________________________________________________________________________________
_FString::_FString (long in_data){
   the_string = new _StringBuffer (new _String (in_data));
}

//__________________________________________________________________________________
_FString::_FString (const _FString& source){
  the_string = new _StringBuffer (*source.the_string);
}

//__________________________________________________________________________________
_FString::_FString (_FString&& source){
  if (source.CanFreeMe()) {
    the_string = source.the_string;
    source.the_string = nil;
  } else {
    the_string = new _StringBuffer (*source.the_string);
  }
}

//__________________________________________________________________________________
_FString::_FString (_String const& data, bool meta) {
    if (meta) {
 
        the_string = new _StringBuffer(data.length());
        for (unsigned long k=0; k < data.length(); k++) {
            char c = data.char_at(k);
            if (c=='\\') {
                  c = data.char_at(k+1UL);
                    switch (c) {
                      case 'n':
                          *the_string<<'\n';
                          break;
                      case 't':
                          *the_string<<'\t';
                          break;
                      case '\\':
                          *the_string<<'\\';
                          break;
                      default:
                          *the_string<<'\\' <<c;
                    }
                  k++;
                  continue;
             }
             *the_string<<c;
        }
    } else {
        the_string = new _StringBuffer (data);
    }
}

//__________________________________________________________________________________
BaseRef _FString::makeDynamic (void) const {
  if (the_string) {
    return new _FString(*the_string, false);
  }
  return new _FString;
}

//__________________________________________________________________________________
void _FString::Duplicate (BaseRefConst o) {
    DeleteObject (the_string);
    the_string = (_StringBuffer*)((_FString const*)o)->get_str().makeDynamic();
}

//__________________________________________________________________________________
void _FString::SetStringContent (_StringBuffer * arg){
    DeleteObject (the_string);
    the_string = arg;
 }


//__________________________________________________________________________________
HBLObjectRef _FString::Add (HBLObjectRef p, HBLObjectRef cache) {
    _StringBuffer  * res;
  
    if (p->ObjectClass()==STRING) {
      _FString* add_this = (_FString*)p;
      res    = new _StringBuffer (get_str().length() + add_this->get_str().length());
      (*res) << get_str()  << add_this->get_str();
    } else {
      res    = new _StringBuffer (get_str().length());
      ((*res) << get_str()).AppendNewInstance((_String*)p->toStr());
    }
    return new _FString (res, false);
}

//__________________________________________________________________________________

HBLObjectRef _FString::Sum (HBLObjectRef cache) {
  return _returnConstantOrUseCache (get_str().to_float(), cache);
}

//__________________________________________________________________________________
long _FString::AddOn (HBLObjectRef p) {
    if (p->ObjectClass()==STRING) {
        *the_string << ((_FString*)p)->get_str();
        return ((_FString*)p)->get_str().length();
    } else if (p->ObjectClass()==NUMBER) {
        long s = p->Value();
        if (s) {
            DeleteObject (the_string);
            the_string = new _StringBuffer (s);
            return s;
        }
    } else {
        HandleApplicationError ("Invalid 2nd argument in call to string*number");
    }
    return 0;
}

//__________________________________________________________________________________
HBLObjectRef _FString::AreEqual (HBLObjectRef p, HBLObjectRef cache) {
    if (p->ObjectClass()==STRING) {
        return _returnConstantOrUseCache (get_str() == ((_FString*)p)->get_str(), cache);
    } else {
        return _returnConstantOrUseCache (0., cache);
    }
}

//__________________________________________________________________________________
HBLObjectRef _FString::AreEqualCIS (HBLObjectRef p, HBLObjectRef cache) {
  if (p->ObjectClass()==STRING) {
      return _returnConstantOrUseCache (get_str().CompareIgnoringCase(((_FString*)p)->get_str()) == kCompareEqual, cache);
  } else {
    return _returnConstantOrUseCache (0., cache);
  }
}

//__________________________________________________________________________________
HBLObjectRef _FString::Join (HBLObjectRef p, HBLObjectRef cache) {
    _List theStrings;

    if (p->ObjectClass()==MATRIX) {
        ((_Matrix*)(p->Compute()))->FillInList (theStrings,true);
    } else if (p->ObjectClass()==ASSOCIATIVE_LIST) {
        ((_AssociativeList*)(p->Compute()))->FillInList (theStrings);
    }

    return _returnStringOrUseCache ((_String*)theStrings.Join(get_str()), cache);
}

//__________________________________________________________________________________
HBLObjectRef _FString::EqualAmb (HBLObjectRef p, HBLObjectRef cache) {
    bool result;
    if (p->ObjectClass()==STRING) {
        result = the_string->EqualWithWildChar (((_FString*)p)->get_str());
    } else {
        _String  convStr      ((_String*)p->toStr());
        result = the_string->EqualWithWildChar(convStr);
    }

    return _returnConstantOrUseCache (result, cache);

}

//__________________________________________________________________________________

HBLObjectRef _FString::EqualRegExp (HBLObjectRef p, bool match_all, HBLObjectRef cache)
{
    if (p->ObjectClass()==STRING) {
       _SimpleList matches  = match_all ? the_string->RegExpAllMatches(((_FString*)p)->get_str(), true, true) :
                                          the_string->RegExpMatch     (((_FString*)p)->get_str(), true, true);
      
        if (matches.empty()) {
            matches << -1 << -1;
        }
        _Matrix * res = new _Matrix (matches);
        res->Transpose();
        return res;
    } else {
        HandleApplicationError ("Invalid 2nd argument in call to string$reg.exp.");
        return new _Matrix (2,1,false,true);
    }
}

//__________________________________________________________________________________
HBLObjectRef _FString::ReplaceReqExp (HBLObjectRef p, HBLObjectRef cache) {
    if (has_data ()) {
      
        if (p->ObjectClass()==MATRIX) {
            _Matrix * m = (_Matrix*)p;

            if (m->IsAStringMatrix() && m->GetHDim() * m->GetVDim() >= 2) {
              
                _FString*  pattern      = (_FString*)m->GetFormula(0, 0)->Compute(),
                         * replacement  = (_FString*)m->GetFormula(1,-1)->Compute();

  
                _SimpleList matches = get_str().RegExpAllMatches(pattern->get_str(), true, true);
              
                _FString * res;
              
                if (matches.empty() == false) {
                    _StringBuffer * replaced = new _StringBuffer (get_str().length() + 1L);
                    long      idx  = matches.Element (0L),
                              midx = 0;

                    for (long k=0; k<the_string->length();) {
                        if (k==idx) {
                            *replaced << replacement->get_str();
                            k = matches.Element (midx+1)+1;
                            midx += 2;
                            if (midx == matches.lLength) {
                                idx = -1;
                            } else {
                                idx = matches.Element (midx);
                            }
                        } else {
                            *replaced << the_string->char_at(k++);
                        }
                    }
                    return _returnStringOrUseCache (replaced, cache);
                    //res = new _FString (replaced);
                } else {
                    return _returnStringOrUseCache (get_str(), cache);
                    //res = new _FString (get_str(), false);
                }

                return res;
            }
        }

        HandleApplicationError ("Invalid 2nd argument in call to string^{{pattern,replacement}}");
    }
    return _returnStringOrUseCache(kEmptyString,cache);
}

//__________________________________________________________________________________

hyComparisonType _FString::Compare  (HBLObjectRef p, bool convert_non_strings) {
  hyComparisonType result = kCompareUndefined;
  if (p->ObjectClass()==STRING) {
    result = the_string->Compare(((_FString*)p)->get_str());
  } else {
    if (convert_non_strings) {
      _String  convStr      ((_String*)p->toStr());
      result = the_string->Compare(convStr);
    }
  }
  return result;

}

//__________________________________________________________________________________

HBLObjectRef _FString::NotEqual (HBLObjectRef p, HBLObjectRef cache) {
    return _returnConstantOrUseCache(Compare (p, true) != kCompareEqual, cache);
}

//__________________________________________________________________________________
HBLObjectRef _FString::Less (HBLObjectRef p, HBLObjectRef cache) {
  return _returnConstantOrUseCache(Compare (p, true) == kCompareLess, cache);
}

//__________________________________________________________________________________
HBLObjectRef _FString::LessEq (HBLObjectRef p, HBLObjectRef cache) {
  hyComparisonType r = Compare (p, true);
  return _returnConstantOrUseCache (r == kCompareLess || r == kCompareEqual, cache);
}

  //__________________________________________________________________________________
HBLObjectRef _FString::GreaterEq (HBLObjectRef p, HBLObjectRef cache) {
  hyComparisonType r = Compare (p, true);
  return _returnConstantOrUseCache (r == kCompareGreater || r == kCompareEqual, cache);
}

  //__________________________________________________________________________________
HBLObjectRef _FString::Greater (HBLObjectRef p, HBLObjectRef cache) {
  return _returnConstantOrUseCache (Compare (p, true) == kCompareGreater, cache);
}


//__________________________________________________________________________________
HBLObjectRef _FString::Differentiate (HBLObjectRef p, HBLObjectRef cache) {
    _Formula F;

    _String  *X,
             *DFDX = nil;

    bool     deleteX = false;

    if (p->ObjectClass()==STRING) {
        X = ((_FString*)p)->the_string;
    } else {
        deleteX = true;
        X = (_String*)p->toStr();
    }

    _String copyMe (*the_string);
  
    _FormulaParsingContext fpc;
    if (Parse (&F,copyMe, fpc, nil) == HY_FORMULA_EXPRESSION) {
        _Formula *DF = F.Differentiate (*X,true);
        if (DF) {
            DFDX = (_String*)DF->toStr(kFormulaStringConversionNormal);
        }
    }

    if (deleteX) {
        DeleteObject (X);
    }

    if (DFDX)
        return _returnStringOrUseCache(DFDX, cache);
    
    return _returnStringOrUseCache (kEmptyString, cache);
    
    //return new _FString (DFDX?DFDX:new _String());

}


//__________________________________________________________________________________
BaseRef  _FString::toStr(unsigned long) {
    the_string->AddAReference();
    return the_string;
}

//__________________________________________________________________________________
HBLObjectRef _FString::RerootTree (HBLObjectRef root, HBLObjectRef cache) {
    // TODO SKLP 20170921 This needs algorithmic review
  
    static const _String _internal_reroot_tree ("_INTERNAL_REROOT_TREE_");

  
    long     stashed_model_id = lastMatrixDeclared;
    lastMatrixDeclared      = HY_NO_MODEL;
    /* unset current model; do not want the internal tree to have an attached model */

    _TheTree    rTree (_internal_reroot_tree,get_str());
 
    if (rTree.IsDegenerate()) { // no need to reroot two-sequence trees
        lastMatrixDeclared = stashed_model_id;
        DeleteVariable  (_internal_reroot_tree);
        return _returnStringOrUseCache(get_str(), cache);
        //return new _FString (get_str(), false);
    }

    if (terminate_execution) {
        lastMatrixDeclared = stashed_model_id;
        DeleteVariable  (_internal_reroot_tree);
        return _returnStringOrUseCache(kEmptyString, cache);
        //return new _FString;
    }
  

    _CalcNode   *rerootAt = nil;
  
    node<long>* counted_descendants = rTree.theRoot->duplicate_tree(node_count_descendants);

    long        maxMin         = 0L,
                totalNodeCount = counted_descendants->in_object + 1L;
  
    hyFloat  bRatio  = 0.0;
  
    node_iterator<long> ni (counted_descendants, _HY_TREE_TRAVERSAL_POSTORDER);
    _TreeIterator ti (&rTree, _HY_TREE_TRAVERSAL_POSTORDER);

    while       (_CalcNode * iterator = ti.Next()) {
        node<long>* counter_tree = ni.Next();
        long      nodeMin    = totalNodeCount-counter_tree->in_object-1L;
        hyFloat thisRatio = nodeMin/(1L+counter_tree->in_object);

        if (thisRatio>1.0) {
            thisRatio = 1.0/thisRatio;
        }

        if (counter_tree->is_leaf()) {
          nodeMin = 1L;
        } else {
            for (int k = counter_tree->get_num_nodes(); k; k--) {
              long tt = counter_tree->go_down(k)->in_object;
              if (tt<nodeMin) {
                nodeMin = tt;
              }
            }
        }

        if (nodeMin>maxMin || (nodeMin==maxMin && thisRatio>bRatio)) {
            bRatio = thisRatio;
            maxMin = nodeMin;
            rerootAt = iterator;
            if (counter_tree->is_root()) {
                rerootAt = nil;
            }
        }
    }
  
    counted_descendants->delete_tree(true);

    HBLObjectRef res;
    if (rerootAt) {
        _FString    rAt  (rerootAt->ContextFreeName());
        res = (_FString*)rTree.RerootTree (&rAt, cache);
    } else {
        res = _returnStringOrUseCache (get_str(), cache);
    }

    DeleteVariable  (_internal_reroot_tree);
    lastMatrixDeclared = stashed_model_id;

    return      res;
}

//__________________________________________________________________________________

HBLObjectRef _FString::Evaluate (_hyExecutionContext* context) {
    if (has_data ()) {
        _String     s (get_str());
        _Formula    evaluator (s, (_VariableContainer*)context->GetContext());
        HBLObjectRef   evalTo = evaluator.Compute(0,(_VariableContainer*)context->GetContext());

        if (evalTo && !terminate_execution) {
            evalTo->AddAReference();
            return evalTo;
        }
    }
    return new _MathObject;
}

  //__________________________________________________________________________________

HBLObjectRef _FString::SubstituteAndSimplify(HBLObjectRef arguments, HBLObjectRef cache) {
  /**
   "arguments" is expected to be a dictionary of with key : value pairs like
    "x" : 3, 
    "y" : {{1,2}}
   
    etc
   */
  if (has_data()) {
    _String     s (get_str());
    _Formula    evaluator (s);
      
    _Polynomial*    is_poly = (_Polynomial*)evaluator.ConstructPolynomial();
    if (is_poly) {
        _Formula pf (is_poly);
        evaluator.Duplicate(&pf);
        //printf ("%s\n", _String ((_String*)simplified_polynomial->toStr(kFormulaStringConversionNormal)).get_str());
        
        //printf ("\n RESULT : %s \n", _String ((_String*)simplified_polynomial->toStr(kFormulaStringConversionNormal)).get_str());
    }
    
    
    if (!terminate_execution) {
      _AssociativeList* argument_substitution_map = (_AssociativeList*) (arguments->ObjectClass() == ASSOCIATIVE_LIST ? arguments : nil);
      if (argument_substitution_map) { // do direct argument substitution
        for (unsigned long expression_term = 0UL; expression_term < evaluator.Length(); expression_term++) {
          _Operation* current_term       = evaluator.GetIthTerm(expression_term);
          _Variable * variable_reference = current_term->RetrieveVar();
          if (variable_reference) {
            HBLObjectRef replacement = argument_substitution_map->GetByKey (*variable_reference->GetName());
            if (replacement) {
              current_term->SetAVariable(-1);
              current_term->SetNumber ((HBLObjectRef)replacement->makeDynamic());
            }
          }
        }
      }
      
      evaluator.SimplifyConstants();
      return _returnStringOrUseCache(((_String*)evaluator.toStr(kFormulaStringConversionNormal)), cache);
    }
  }
  return new _MathObject;
}

//__________________________________________________________________________________

HBLObjectRef _FString::Dereference(bool ignore_context, _hyExecutionContext* context, bool return_variable_ref) {
   HBLObjectRef result = nil;
   _String referencedVariable;
   if (has_data()) {
      referencedVariable = get_str();
      if (!ignore_context && context) {
          referencedVariable = AppendContainerName(referencedVariable, (_VariableContainer *)context->GetContext());
      }
      if (return_variable_ref) {
          return FetchVar (LocateVarByName(referencedVariable));
      }
      result = FetchObjectFromVariableByType(&referencedVariable, HY_ANY_OBJECT);
   }
    //printf ("\n\nDereferencing %s in this context %x\n\n", referencedVariable.sData, context);
    if (!result) {
        _String errM = _String("Failed to dereference '") & referencedVariable & "'";
        if (context) {
            context->ReportError(errM);
        } else {
            HandleApplicationError (errM);
        }
        result = new _FString;
    } else {
        result->AddAReference();
    }
    return result;
}

//__________________________________________________________________________________


HBLObjectRef _FString::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context, HBLObjectRef cache)   {
  
  switch (opCode) { // first check operations without arguments
    case HY_OP_CODE_NOT: // !
      return FileExists(cache);
    case HY_OP_CODE_ABS: // Abs
      return _returnConstantOrUseCache (get_str().length(), cache);
    case HY_OP_CODE_EVAL: // Eval
        return Evaluate(context);
    case HY_OP_CODE_EXP: // Exp
      return _returnConstantOrUseCache (get_str().LempelZivProductionHistory(nil), cache);
    case HY_OP_CODE_LOG: // Log - check sum
      return _returnConstantOrUseCache (get_str().Adler32(), cache);
    case HY_OP_CODE_INVERSE:  // Inverse
      return _returnStringOrUseCache(get_str().Reverse(), cache);
    case HY_OP_CODE_MCOORD: // MCoord
      return _returnStringOrUseCache(get_str(), cache);
    case HY_OP_CODE_TYPE: // Type
      return Type(cache);
    case HY_OP_CODE_REROOTTREE: // RerootTree
      return RerootTree (nil, cache);
    case HY_OP_CODE_ROWS: // Count Objects of given type
      return CountGlobalObjects(cache);
  }
  
  _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
  switch (opCode) { // next check operations without arguments or with at least one argument
    case HY_OP_CODE_MUL: // *
      if (arg0) {
        // NOT a dereference
        if (arg0->ObjectClass() == MATRIX) {
          return      MapStringToVector (arg0, cache);
        } else {
          return _returnConstantOrUseCache (AddOn(arg0), cache);
        }
      } else {
          return Dereference(false, context, cache);
      }
    case HY_OP_CODE_ADD: // +
      if (arg0) {
        return Add(arg0, cache);
      } else {
        return Sum(cache);
      }
    case HY_OP_CODE_POWER: {
      // Replace (^)
      if (arg0)
        return ReplaceReqExp (arg0, cache);
      return Dereference(true, context, cache);
    }
      
    case HY_OP_CODE_CALL: // call the function
      return Call (arguments, context, cache);
  }
  
  if (arg0) {
    switch (opCode) {
      case HY_OP_CODE_NEQ: // !=
        return NotEqual(arg0, cache);
      case HY_OP_CODE_IDIV: // $ match regexp
        return EqualRegExp(arg0, false, cache);
      case HY_OP_CODE_MOD: // % equal case insenstive
        return AreEqualCIS(arg0, cache);
            
      case HY_OP_CODE_MIN: {// resolve file path
          hyFloat local_resolve = 0.0;
          if (arg0->ObjectClass () == NUMBER) {
              local_resolve = arg0->Value();
          }
          _String value = get_str();
          ProcessFileName (value, false, false, (hyPointer) (context ? context->GetContext() : nil), false, nil, !CheckEqual(local_resolve,0.0));
          return _returnStringOrUseCache(value, cache);
      }
            
      case HY_OP_CODE_AND: { // && upcase or lowercase
        hyFloat pVal = 0.0;
        if (arg0->ObjectClass () == NUMBER) {
          pVal = arg0->Value();
        }
        
        if (pVal < 0.0) {
          return (HBLObjectRef)makeDynamic();
        } else {
          _String * t = nil;
          
          long      conversion_type = pVal;
          
          switch (conversion_type) {
            case 1L:
              t = new _String( get_str().ChangeCase(kStringUpperCase));
              break;
            case 2L: {
              t = &(new _StringBuffer)->SanitizeAndAppend(get_str());
              break;
            }
            case 3L: {
              t = &(new _StringBuffer)->SanitizeForPostScriptAndAppend(get_str());
              break;
            }
            case 4L: {
              t = &(new _StringBuffer)->SanitizeForSQLAndAppend(get_str());
              break;
            }
            case 5L: {
              t = &(new _StringBuffer)->SanitizeForHTMLAndAppend(get_str());
              break;
            }
             
            case 6L: {
              t = &(new _StringBuffer)->SanitizeForRegExAndAppend(get_str());
              break;
            }
                  
            case 7L: {
                t = &(*(new _StringBuffer) <<  get_str().ConvertToAnIdent(fIDAllowCompound));
                break;
            }
              
            default:
              t = new _String( get_str().ChangeCase(kStringLowerCase));
              break;
          }
          
          return _returnStringOrUseCache (t, cache);
        }
      }
        
      case HY_OP_CODE_DIV: // /
        return EqualAmb(arg0, cache);
      case HY_OP_CODE_LESS: // <
        return Less(arg0, cache);
      case HY_OP_CODE_LEQ: // <=
        return LessEq(arg0, cache);
      case HY_OP_CODE_EQ: // ==
        return AreEqual(arg0, cache);
      case HY_OP_CODE_GREATER: // >
        return Greater(arg0, cache);
      case HY_OP_CODE_GEQ: // >=
        return GreaterEq(arg0, cache);
      case HY_OP_CODE_DIFF: // Differentiate
        return Differentiate(arg0, cache);
      case HY_OP_CODE_JOIN: // JOIN
        return Join (arg0, cache);
      case HY_OP_CODE_SIMPLIFY: // Simplify an expression
        return SubstituteAndSimplify (arg0, cache);
      case HY_OP_CODE_OR: // Match all instances of the reg.ex (||)
        return EqualRegExp (arg0, true, cache);
    }
    
    _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
    
    switch (opCode) {
      case HY_OP_CODE_MACCESS: // MAccess
        return CharAccess(arg0,arg1,cache);
    }
    
    if (arg1) {
      switch (opCode) {
        case HY_OP_CODE_FORMAT: { // Format
          HBLObjectRef fv = nil;
          try {
           fv = _Formula::ParseAndCompute (get_str(), true, NUMBER, context);
          } catch (_String const &e) {
            
          }
          if (fv) {
            return ((_Constant*)fv)->FormatNumberString (arg0,arg1, cache);
          } else {
            ReportWarning (_String("Failed to evaluate ")& get_str() & " to a number in call to Format (string...)");
            return new _FString();
          }
        }
        break;
      }
    }
  }
  
  switch (opCode) {
    case HY_OP_CODE_NEQ: // !=
    case HY_OP_CODE_IDIV: // $ match regexp
    case HY_OP_CODE_MOD: // % equal case insenstive
    case HY_OP_CODE_AND:// && upcase or lowercase
    case HY_OP_CODE_DIV: // /
    case HY_OP_CODE_LESS: // <
    case HY_OP_CODE_LEQ: // <=
    case HY_OP_CODE_EQ: // ==
    case HY_OP_CODE_GREATER: // >
    case HY_OP_CODE_GEQ: // >=
    case HY_OP_CODE_DIFF: // Differentiate
    case HY_OP_CODE_JOIN: // Inverse
    case HY_OP_CODE_SIMPLIFY: // Simplify an expression
    case HY_OP_CODE_OR: // Match all instances of the reg.ex (||)
    case HY_OP_CODE_MACCESS: // MAccess
    case HY_OP_CODE_FORMAT:  // Format
      WarnWrongNumberOfArguments (this, opCode,context, arguments);
      break;
    default:
      WarnNotDefined (this, opCode,context);

  }
  
  return new _MathObject;
  
}

//__________________________________________________________________________________
HBLObjectRef   _FString::MapStringToVector (HBLObjectRef p, HBLObjectRef cache) {
    if (has_data() && p->ObjectClass () == MATRIX) {
        _Matrix         * factoringMatrix = (_Matrix *)p;

        if ((factoringMatrix->is_row () || factoringMatrix->is_column ()) && factoringMatrix->IsAStringMatrix()) {
            long            mapper [255],
                            keys    = factoringMatrix->GetHDim() * factoringMatrix->GetVDim(),
                            byRows   = factoringMatrix->is_column ();

            InitializeArray(mapper, 255 , -1L);
 
            for (long r = 0; r < keys; r++) {
                _FString* aKey = (_FString*)factoringMatrix->GetFormula(byRows?r:0,byRows?0:r)->Compute();
                if (aKey->get_str().length() == 1UL) {
                    unsigned char thisChar = aKey->get_str().get_uchar (0);
                    if (mapper[thisChar] < 0L) {
                        mapper[thisChar] = r;
                    }
                }
            }

            _SimpleList mapped;
          
            for (long s = 0; s < get_str().length(); s++) {
                mapped << mapper[get_str().get_uchar(s)];
            }

            return new _Matrix (mapped);
        }
    }

    return new _Matrix;
}

//__________________________________________________________________________________
HBLObjectRef   _FString::CharAccess (HBLObjectRef p,HBLObjectRef p2,HBLObjectRef cache)
{
    long index = p->Value();

    if (p2) {
        long index2 = p2->Value();
        return _returnStringOrUseCache(get_str().Cut (index,index2), cache);
        
        //return new _FString (new _String (get_str().Cut (index,index2)));
    } else if (index<get_str().length()) {
        //return new _FString (new _String (get_str().char_at(index)));
        char buffer[2];
        buffer[0] = get_str()(index);
        buffer[1] = 0;
        return _returnStringOrUseCache(buffer, cache);
     }

    return _returnStringOrUseCache(kEmptyString, cache);
    //return new _FString (kEmptyString, false);
}
//__________________________________________________________________________________
HBLObjectRef   _FString::FileExists (HBLObjectRef cache) {
    hyFloat exists = 0.0;
    if (has_data()) {
        _String cpy (get_str());
        if (ProcessFileName(cpy)) {
            // TODO use fstat instead
          FILE * test = doFileOpen (cpy, "rb");
          if (test) {
              exists = 1.0;
              fclose (test);
          }
        }
    }
    return _returnConstantOrUseCache(exists, cache);
}

//__________________________________________________________________________________
HBLObjectRef   _FString::Call (_List* arguments, _hyExecutionContext* context, HBLObjectRef cache) {
  long function_id = FindBFFunctionName (get_str(), NULL);
  if (function_id >= 0) {
       _Formula the_call;
    
      if (arguments) {
        for (long k = 0; k < arguments->countitems() ; k ++) {
          HBLObjectRef payload = (HBLObjectRef)arguments->GetItem (k);
          _Operation *arg_k = new _Operation (payload);
          payload->AddAReference();
          the_call.PushTerm(arg_k);
          arg_k->RemoveAReference();
        }
      }
      
      _Operation * function_call_term = new _Operation (function_id, -1L-(arguments?arguments->countitems():0L));
      the_call.PushTerm(function_call_term);
      DeleteObject (function_call_term);
      
      HBLObjectRef result = the_call.Compute();
      result->AddAReference();
      return result;
      
  } else {
    HandleApplicationError (_String ("The first argument ('") & get_str() & "') to 'Call' was not an HBL function name");
  }
  
  return new _MathObject;
    
}

//__________________________________________________________________________________
HBLObjectRef   _FString::CountGlobalObjects (HBLObjectRef cache)
{
    hyFloat res = 0.0;

    long      standardType = _HY_GetStringGlobalTypes.Find(&get_str());
    if (standardType >=0 ) {
        standardType = _HY_GetStringGlobalTypes.GetXtra (standardType);
    }

    switch (standardType) {
    case HY_BL_LIKELIHOOD_FUNCTION:
        return _returnConstantOrUseCache(likeFuncList.lLength, cache);
    case HY_BL_DATASET:
        return _returnConstantOrUseCache(dataSetList.lLength, cache);
    case HY_BL_DATASET_FILTER:
        return _returnConstantOrUseCache(CountObjectsByType (HY_BL_DATASET_FILTER), cache);
    case HY_BL_HBL_FUNCTION:
        return _returnConstantOrUseCache(GetBFFunctionCount(), cache);
    case HY_BL_TREE: {
        _SimpleList tc;
        long        si,
                    vi = variableNames.Traverser (tc,si,variableNames.GetRoot());

        for (; vi >= 0; vi = variableNames.Traverser (tc,si))
            if (((_Variable*)FetchVar(vi))->ObjectClass () == TREE) {
                res += 1.;
            }

        break;
    }

    case HY_BL_SCFG:
        return _returnConstantOrUseCache(scfgList.lLength, cache);
    case HY_BL_VARIABLE:
        return _returnConstantOrUseCache(variableNames.countitems(), cache);
    }

    if (standardType < 0) {
        if (get_str() == hy_env::last_model_parameter_list) {
            if (lastMatrixDeclared>=0) {
              res =  PopulateAndSort ([&] (_AVLList & parameter_list) -> void  {
                LocateVar (modelMatrixIndices.list_data[lastMatrixDeclared])->ScanForVariables(parameter_list, false);
              }).countitems();
            }
        }
    }
    return _returnConstantOrUseCache (res, cache);
}


