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

#include "defines.h"
#include "fstring.h"
#include "constant.h"
#include "matrix.h"
#include "calcnode.h"
#include "batchlan.h"
#include "global_object_lists.h"
#include "global_things.h"

extern long lastMatrixDeclared;
extern _AVLListX _HY_GetStringGlobalTypes;

extern _List likeFuncList,
             batchLanguageFunctionNames,
             dataSetList,
             scfgList;

extern _SimpleList modelMatrixIndices;
extern _String lastModelParameterList;

using namespace hyphy_global_objects;
using namespace hy_global;

_String internalRerootTreeID("_INTERNAL_REROOT_TREE_");
//__________________________________________________________________________________
_FString::_FString (void)
{
    theString = new _String;
}

//__________________________________________________________________________________
_FString::~_FString ()
{
    if (CanFreeMe()) {
        DeleteObject (theString);
    } else {
        RemoveAReference();
    }

}

//__________________________________________________________________________________
_FString::_FString (_String* data) {
    theString = data;
}

//__________________________________________________________________________________
_FString::_FString (long inData){
   theString = new _String (inData);
}

//__________________________________________________________________________________
_FString::_FString (_String const& data, bool meta) {
    if (meta) {
        unsigned long ssi = _String::storageIncrement;

        if (data.sLength>ssi) {
            _String::storageIncrement = data.sLength;
        }

        theString = new _String (data.sLength,true);
        for (unsigned long k=0; k<data.sLength; k++) {
            char c = data.sData[k];
            if (c=='\\') {
                if (k<data.sLength-1) {
                    c = data.sData[k+1];
                    switch (c) {
                    case 'n':
                        (*theString)<<'\n';
                        break;
                    case 't':
                        (*theString)<<'\t';
                        break;
                    case '\\':
                        (*theString)<<'\\';
                        break;
                    default:
                        (*theString)<<'\\';
                        (*theString)<<c;
                    }
                    k++;
                    continue;
                }
            }
            (*theString)<<c;
        }
        _String::storageIncrement = ssi;
        theString->Finalize();
    } else {
        theString = new _String (data);
    }
}

//__________________________________________________________________________________
BaseRef _FString::makeDynamic (void) const {
    return new _FString(*theString, false);
}

//__________________________________________________________________________________
void _FString::Duplicate (BaseRefConst o) {
    DeleteObject (theString);
    _FString const * m = (_FString const*)o;
    theString = (_String*)m->theString->makeDynamic();
}

//__________________________________________________________________________________
_PMathObj _FString::Add (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _String * res    = new _String (theString->sLength+theStr->theString->sLength,true);
        (*res) << theString << theStr->theString;
        res->Finalize();
        return new _FString (res);
     }
    _String* convStr = (_String*)p->toStr();
    _String res (*theString& *((_String*)convStr));
    DeleteObject (convStr);
    return new _FString (res, false);
}

//__________________________________________________________________________________

_PMathObj _FString::Sum (void) {
  return new _Constant (theString->toNum());
}

//__________________________________________________________________________________
long _FString::AddOn (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        (*theString) << theStr->theString;
        return theStr->theString->sLength;
    } else if (p->ObjectClass()==NUMBER) {
        long s = p->Value();
        if (s) {
            delete theString;
            theString = new _String (s, true);
            return s;
        } else {
            theString->Finalize();
        }
    } else {
        HandleApplicationError ("Invalid 2nd argument in call to string*number");
    }
    return 0;
}

//__________________________________________________________________________________
_PMathObj _FString::AreEqual (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Equal(theStr->theString);
        return new _Constant ((hy_float)equal);
    } else {
        /*_String* convStr = (_String*)p->toStr();
         bool     equal = theString->Equal(convStr);
         DeleteObject (convStr);
         return new _Constant ((hy_float)equal);*/
         return new HY_CONSTANT_FALSE;
    }
}

//__________________________________________________________________________________
_PMathObj _FString::AreEqualCIS (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _String   t1 (*theString),
                  t2 (*theStr->theString);
        t1.UpCase();
        t2.UpCase();
        bool     equal = t1.Equal(&t2);
        return new _Constant ((hy_float)equal);
    } else {
        return AreEqual (p);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Join (_PMathObj p)
{
    _List theStrings;

    if (p->ObjectClass()==MATRIX) {
        ((_Matrix*)(p->Compute()))->FillInList (theStrings,true);
    } else if (p->ObjectClass()==ASSOCIATIVE_LIST) {
        ((_AssociativeList*)(p->Compute()))->FillInList (theStrings);
    }

    return new _FString((_String*)theStrings.Join(theString));
}

//__________________________________________________________________________________
_PMathObj _FString::EqualAmb (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->EqualWithWildChar(theStr->theString);
        return new _Constant ((hy_float)equal);
    } else {
        _String  convStr      ((_String*)p->toStr());
        return   new _Constant(theString->EqualWithWildChar(&convStr));
    }
}

//__________________________________________________________________________________

_PMathObj _FString::EqualRegExp (_PMathObj p, bool matchAll)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _SimpleList matches;

        if (matchAll) {
            int errNo = 0;
            hy_pointer regex = PrepRegExp (theStr->theString, errNo, true);
            if (regex) {
                theString->RegExpMatchAll(regex, matches);
                FlushRegExp (regex);
            } else {
                HandleApplicationError (GetRegExpError (errNo));
            }
        } else {
            theString->RegExpMatchOnce(theStr->theString, matches, true, true);
        }

        if (matches.lLength == 0) {
            matches << -1;
            matches << -1;
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
_PMathObj _FString::ReplaceReqExp (_PMathObj p)
{
    if (theString && theString->sLength) {
        if (p->ObjectClass()==MATRIX) {
            _Matrix * m = (_Matrix*)p;

            if (m->IsAStringMatrix() && m->GetHDim() * m->GetVDim() >= 2) {
                _FString* theStr  = (_FString*)m->GetFormula(0,0)->Compute(),
                          * repWith = (_FString*)m->GetFormula(1,-1)->Compute();

                _SimpleList matches;

                int errNo = 0;
                hy_pointer regex = PrepRegExp (theStr->theString, errNo, true);

                if (!regex) {
                    HandleApplicationError (GetRegExpError (errNo));
                    return new _FString (kEmptyString);
                }

                theString->RegExpMatchAll(regex, matches);
                _FString * res;
                if (matches.lLength) {
                    _String * newString = new _String (theString->sLength+1,true);
                    long      idx  = matches.lData[0],
                              midx = 0;

                    for (long k=0; k<theString->sLength;) {
                        if (k==idx) {
                            (*newString) << repWith->theString;
                            k = matches.lData[midx+1]+1;
                            midx += 2;
                            if (midx == matches.lLength) {
                                idx = -1;
                            } else {
                                idx = matches.lData[midx];
                            }
                        } else {
                            (*newString) << theString->sData[k++];
                        }
                    }
                    newString->Finalize();
                    res = new _FString (newString);
                } else {
                    res = new _FString (*theString,false);
                }

                FlushRegExp (regex);
                return res;
            }
        }

        HandleApplicationError ("Invalid 2nd argument in call to string^{{pattern,replacement}}");
    }
    return new _FString (kEmptyString,false);
}

//__________________________________________________________________________________
_PMathObj _FString::NotEqual (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Equal(theStr->theString);
        return new _Constant ((hy_float)!equal);
    } else {
        //_String* convStr = (_String*)p->toStr();
        //bool     equal = theString->Equal(convStr);
        //DeleteObject (convStr);
        //return new _Constant ((hy_float)!equal);
        return new HY_CONSTANT_TRUE;
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Less (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Less(theStr->theString);
        return new _Constant ((hy_float)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Less(convStr);
        DeleteObject (convStr);
        return new _Constant ((hy_float)equal);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::LessEq (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Less(theStr->theString)||theString->Equal(theStr->theString);
        return new _Constant ((hy_float)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Less(convStr)||theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((hy_float)equal);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Differentiate (_PMathObj p)
{
    _Formula F;

    _String  *X,
             *DFDX = nil;

    bool     deleteX = false;

    if (p->ObjectClass()==STRING) {
        X = ((_FString*)p)->theString;
    } else {
        deleteX = true;
        X = (_String*)p->toStr();
    }

    _String copyMe (*theString);
    _FormulaParsingContext fpc;
    if (Parse (&F,copyMe, fpc, nil) == HY_FORMULA_EXPRESSION) {
        _Formula *DF = F.Differentiate (*X,true);
        if (DF) {
            DFDX = (_String*)DF->toStr();
        }

    }

    if (deleteX) {
        DeleteObject (X);
    }

    return new _FString (DFDX?DFDX:new _String());

}

//__________________________________________________________________________________
_PMathObj _FString::GreaterEq (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Greater(theStr->theString)||theString->Equal(theStr->theString);
        return new _Constant ((hy_float)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Greater(convStr)||theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((hy_float)equal);
    }
}

//__________________________________________________________________________________
_PMathObj _FString::Greater (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Greater(theStr->theString);
        return new _Constant ((hy_float)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Greater(convStr);
        DeleteObject (convStr);
        return new _Constant ((hy_float)equal);
    }
}

//__________________________________________________________________________________
BaseRef  _FString::toStr(unsigned long)
{
    theString->AddAReference();
    return theString;
}

//__________________________________________________________________________________
_PMathObj _FString::RerootTree (_PMathObj) {
  
    long     stashedModelID = lastMatrixDeclared;
    lastMatrixDeclared      = HY_NO_MODEL;
    /* unset current model; do not want the internal tree to have an attached model */

    _TheTree    rTree (internalRerootTreeID,*theString);

    if (rTree.IsDegenerate()) { // no need to reroot two-sequence trees
        lastMatrixDeclared = stashedModelID;
        DeleteVariable  (internalRerootTreeID);
        return new _FString (*theString, false);
    }

    if (terminate_execution) {
        lastMatrixDeclared = stashedModelID;
        DeleteVariable  (internalRerootTreeID);
        return new _FString;
    }
  

    _CalcNode   *rerootAt = nil;
  
    node<long>* counted_descendants = rTree.theRoot->duplicate_tree(node_count_descendants);

    long        maxMin         = 0L,
                totalNodeCount = counted_descendants->in_object + 1L;
  
    hy_float  bRatio  = 0.0;
  
    node_iterator<long> ni (counted_descendants, _HY_TREE_TRAVERSAL_POSTORDER);
    _TreeIterator ti (&rTree, _HY_TREE_TRAVERSAL_POSTORDER);

    while       (_CalcNode * iterator = ti.Next()) {
        node<long>* counter_tree = ni.Next();
        long      nodeMin    = totalNodeCount-counter_tree->in_object-1L;
        hy_float thisRatio = nodeMin/(1L+counter_tree->in_object);

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

    _FString* res;
    if (rerootAt) {
        _FString    rAt  (rerootAt->ContextFreeName());
        res = (_FString*)rTree.RerootTree (&rAt);
    } else {
        res = new _FString (*theString, false);
    }

    DeleteVariable  (internalRerootTreeID);
    lastMatrixDeclared = stashedModelID;

    return      res;
}

//__________________________________________________________________________________

_PMathObj _FString::Evaluate (_hyExecutionContext* context)
{
    if (theString && theString->sLength) {
        _String     s (*theString);
        _Formula    evaluator (s, (_VariableContainer*)context->GetContext());
        _PMathObj   evalTo = evaluator.Compute(0,(_VariableContainer*)context->GetContext());

        if (evalTo && !terminate_execution) {
            evalTo->AddAReference();
            return evalTo;
        }
    }
    return new _MathObject;
}

  //__________________________________________________________________________________

_PMathObj _FString::SubstituteAndSimplify(_PMathObj arguments) {
  /**
   "arguments" is expected to be a dictionary of with key : value pairs like
    "x" : 3, 
    "y" : {{1,2}}
   
    etc
   */
  if (theString && theString->sLength) {
    _String     s (*theString);
    _Formula    evaluator (s);
    
    
    if (!terminate_execution) {
      _AssociativeList* argument_substitution_map = (_AssociativeList*) (arguments->ObjectClass() == ASSOCIATIVE_LIST ? arguments : nil);
      if (argument_substitution_map) { // do direct argument substitution
        for (unsigned long expression_term = 0UL; expression_term < evaluator.Length(); expression_term++) {
          _Operation* current_term       = evaluator.GetIthTerm(expression_term);
          _Variable * variable_reference = current_term->RetrieveVar();
          if (variable_reference) {
            _PMathObj replacement = argument_substitution_map->GetByKey (*variable_reference->GetName());
            if (replacement) {
              current_term->SetAVariable(-1);
              current_term->SetNumber ((_PMathObj)replacement->makeDynamic());
            }
          }
        }
      }
      
      evaluator.SimplifyConstants();
      return new _FString ((_String*)evaluator.toStr());
    }
  }
  return new _MathObject;
}

//__________________________________________________________________________________

_PMathObj _FString::Dereference(bool ignore_context, _hyExecutionContext* context, bool return_variable_ref) {
    _String referencedVariable = *theString;
    if (!ignore_context && context) {
        referencedVariable = AppendContainerName(referencedVariable, (_VariableContainer *)context->GetContext());
    }
    if (return_variable_ref) {
        return FetchVar (LocateVarByName(referencedVariable));
    }
    _PMathObj result = FetchObjectFromVariableByType(&referencedVariable, HY_ANY_OBJECT); 
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


_PMathObj _FString::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context)   {
  
  switch (opCode) { // first check operations without arguments
    case HY_OP_CODE_NOT: // !
      return FileExists();
    case HY_OP_CODE_ABS: // Abs
      return new _Constant (theString->sLength);
    case HY_OP_CODE_EVAL: // Eval
        return Evaluate(context);
    case HY_OP_CODE_EXP: // Exp
      return new _Constant (theString->LempelZivProductionHistory(nil));
    case HY_OP_CODE_LOG: // Log - check sum
      return new _Constant (theString->Adler32());
    case HY_OP_CODE_INVERSE:  // Inverse
      return new _FString (new _String (theString->Reverse()));
    case HY_OP_CODE_MCOORD: // MCoord
      return new _FString (*theString, true);
    case HY_OP_CODE_TYPE: // Type
      return Type();
    case HY_OP_CODE_REROOTTREE: // RerootTree
      return RerootTree (nil);
    case HY_OP_CODE_ROWS: // Count Objects of given type
      return CountGlobalObjects();
  }
  
  _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
  switch (opCode) { // next check operations without arguments or with one argument
    case HY_OP_CODE_MUL: // *
      if (arg0) {
        // NOT a dereference
        if (arg0->ObjectClass() == MATRIX) {
          return      MapStringToVector (arg0);
        } else {
          return new _Constant(AddOn(arg0));
        }
      } else {
        return Dereference(false, context);
      }
    case HY_OP_CODE_ADD: // +
      if (arg0) {
        return Add(arg0);
      } else {
        return Sum();
      }
    case HY_OP_CODE_POWER: {
      // Replace (^)
      if (arg0)
        return ReplaceReqExp (arg0);
      return Dereference(true, context);
    }
      
    case HY_OP_CODE_CALL: // call the function
      return Call (arguments, context);
  }
  
  if (arg0) {
    switch (opCode) {
      case HY_OP_CODE_NEQ: // !=
        return NotEqual(arg0);
      case HY_OP_CODE_IDIV: // $ match regexp
        return EqualRegExp(arg0);
      case HY_OP_CODE_MOD: // % equal case insenstive
        return AreEqualCIS(arg0);
      case HY_OP_CODE_AND: { // && upcase or lowercase
        hy_float pVal = 0.0;
        if (arg0->ObjectClass () == NUMBER) {
          pVal = arg0->Value();
        }
        
        if (pVal < 0.0) {
          return (_PMathObj)makeDynamic();
        } else {
          _String * t = nil;
          
          if (CheckEqual(pVal,2.0) || CheckEqual(pVal,3.0) || CheckEqual(pVal,4.0) || CheckEqual(pVal,5.0) || CheckEqual(pVal,6.0)) {
            t = new _String (theString->sLength+1,true);
            t->EscapeAndAppend (*theString, CheckEqual(pVal,3.0) + 2*CheckEqual(pVal,4.0) + 4*CheckEqual(pVal,5.0) + 5*CheckEqual(pVal,6.0));
            t->Finalize();
          } else {
            t = new _String (*theString);
            if (CheckEqual(pVal,1.0)) {
              t->UpCase();
            } else {
              t->LoCase();
            }
          }
          
          return new _FString (t);
        }
      }
      case HY_OP_CODE_DIV: // /
        return EqualAmb(arg0);
      case HY_OP_CODE_LESS: // <
        return Less(arg0);
      case HY_OP_CODE_LEQ: // <=
        return LessEq(arg0);
      case HY_OP_CODE_EQ: // ==
        return AreEqual(arg0);
      case HY_OP_CODE_GREATER: // >
        return Greater(arg0);
      case HY_OP_CODE_GEQ: // >=
        return GreaterEq(arg0);
      case HY_OP_CODE_DIFF: // Differentiate
        return Differentiate(arg0);
      case HY_OP_CODE_JOIN: // JOIN
        return Join (arg0);
      case HY_OP_CODE_SIMPLIFY: // Simplify an expression
        return SubstituteAndSimplify (arg0);
      case HY_OP_CODE_OR: // Match all instances of the reg.ex (||)
        return EqualRegExp (arg0, true);
    }
    
    _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
    
    switch (opCode) {
      case HY_OP_CODE_MACCESS: // MAccess
        return CharAccess(arg0,arg1);
    }
    
    if (arg1) {
      switch (opCode) {
        case HY_OP_CODE_FORMAT: { // Format
          _String   cpyString (*theString);
          _Formula f (cpyString);
          _PMathObj fv = f.Compute();
          if (fv && fv->ObjectClass () == NUMBER) {
            return ((_Constant*)fv)->FormatNumberString (arg0,arg1);
          } else {
            ReportWarning (_String("Failed to evaluate ")& *theString & " to a number in call to Format (string...)");
            return new _FString();
          }
        }
          
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
_PMathObj   _FString::MapStringToVector (_PMathObj p)
{
    if (theString->sLength && p->ObjectClass () == MATRIX) {
        _Matrix         * factoringMatrix = (_Matrix *)p;

        if (factoringMatrix->IsAVector () && factoringMatrix->IsAStringMatrix()) {
            long            mapper [255],
                            keys    = factoringMatrix->GetHDim() * factoringMatrix->GetVDim(),
                            byRows   = factoringMatrix->IsAVector (HY_MATRIX_COLUMN_VECTOR);

            for (long c = 0; c < 255; c++) {
                mapper[c] = -1;
            }

            for (long r = 0; r < keys; r++) {
                _FString* aKey = (_FString*)factoringMatrix->GetFormula(byRows?r:0,byRows?0:r)->Compute();
                if (aKey->theString->sLength == 1) {
                    unsigned char thisChar = aKey->theString->sData[0];
                    if (mapper[thisChar] < 0) {
                        mapper[thisChar] = r;
                    }
                }
            }

            _SimpleList mapped;
            for (long s = 0; s < theString->sLength; s++) {
                mapped << mapper[theString->getUChar(s)];
            }

            return new _Matrix (mapped);
        }
    }

    return new _Matrix;
}

//__________________________________________________________________________________
_PMathObj   _FString::CharAccess (_PMathObj p,_PMathObj p2)
{
    unsigned long index = p->Value();
    _String res;

    if (p2) {
        unsigned long index2 = p2->Value();
        res = theString->Cut (index,index2);
    } else if (index<theString->sLength) {
        res = theString->sData[index];
    }

    return new _FString (res);
}
//__________________________________________________________________________________
_PMathObj   _FString::FileExists (void) {
    _Constant  * retValue = new _Constant (0.0);
    if (theString) {
        _String cpy (*theString);
        cpy.ProcessFileName();
        FILE * test = doFileOpen (cpy, "rb");
        if (test) {
            retValue->SetValue (1.0);
            fclose (test);
        }
    }
    return retValue;
}

//__________________________________________________________________________________
_PMathObj   _FString::Call (_List* arguments, _hyExecutionContext* context) {
  long function_id = FindBFFunctionName (*theString, NULL);
  if (function_id >= 0) {
       _Formula the_call;
    
      if (arguments) {
        for (long k = 0; k < arguments->countitems() ; k ++) {
          _PMathObj payload = (_PMathObj)arguments->GetItem (k);
          _Operation *arg_k = new _Operation (payload);
          payload->AddAReference();
          the_call.PushTerm(arg_k);
          arg_k->RemoveAReference();
        }
      }
      
      _Operation * function_call_term = new _Operation (function_id, -1L-(arguments?arguments->countitems():0L));
      the_call.PushTerm(function_call_term);
      DeleteObject (function_call_term);
      
      _PMathObj result = the_call.Compute();
      result->AddAReference();
      return result;
      
  } else {
    HandleApplicationError (_String ("The first argument ('") & *theString & "') to 'Call' was not an HBL function name");
  }
  
  return new _MathObject;
    
}

//__________________________________________________________________________________
_PMathObj   _FString::CountGlobalObjects (void)
{
    hy_float res = 0.0;

    long      standardType = _HY_GetStringGlobalTypes.Find(theString);
    if (standardType >=0 ) {
        standardType = _HY_GetStringGlobalTypes.GetXtra (standardType);
    }

    switch (standardType) {
    case HY_BL_LIKELIHOOD_FUNCTION:
        return new _Constant (likeFuncList.lLength);
    case HY_BL_DATASET:
        return new _Constant (dataSetList.lLength);
    case HY_BL_DATASET_FILTER:
        return new _Constant (CountObjectsByType (HY_BL_DATASET_FILTER));
    case HY_BL_HBL_FUNCTION:
        return new _Constant (batchLanguageFunctionNames.lLength);
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
        return new _Constant (scfgList.lLength);
    case HY_BL_VARIABLE:
        return new _Constant (variableNames.countitems());

    }

    if (standardType < 0) {
        if ((*theString)==lastModelParameterList) {
            if (lastMatrixDeclared>=0) {
                _SimpleList p;
                _Variable *theM = LocateVar (modelMatrixIndices.lData[lastMatrixDeclared]);
                {
                    _AVLList pA (&p);
                    theM->ScanForVariables (pA,false);
                    pA.ReorderList();
                }
                res = p.lLength;
            }
        }
    }
    return new _Constant (res);
}
