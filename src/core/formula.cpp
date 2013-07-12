#include <math.h>
#include <float.h>
#include "defines.h"
#include "formula.h"
#include "formulaparsingcontext.h"
#include "legacy_parser.h"
#include "batchlan.h"

//Constants
extern _Parameter twoOverSqrtPi;

extern _List batchLanguageFunctionNames, batchLanguageFunctionParameterLists;

extern _SimpleList BinOps, opPrecedence, FunctionArgumentCount,
    batchLanguageFunctionParameters, batchLanguageFunctionClassification,
    associativeOps;

extern _SimpleList simpleOperationCodes, simpleOperationFunctions;

extern _Variable *_x_, *_n_;

extern _Parameter machineEps;
extern _Parameter tolerance;

extern _String intPrecFact;
extern _String intMaxIter;

_Parameter maxRombergSteps = 8., integrationPrecisionFactor = 1.e-5;

_Formula::_Formula(void) {
  theTree = nil;
  resultCache = nil;
}

//______________________________________________________________________________
_Formula::_Formula(_PMathObj p, bool isAVar) {
  theTree = nil;
  resultCache = nil;
  if (!isAVar) {
    theFormula.AppendNewInstance(new _Operation(p));
  } else {
    _Variable *v = (_Variable *)p;
    theFormula.AppendNewInstance(
        new _Operation(  *v->GetName(), true, v->IsGlobal(), nil));
  }
}

//______________________________________________________________________________
void _Formula::Initialize(void) {}

//______________________________________________________________________________
void _Formula::Duplicate(BaseRef f) {
  _Formula *f_cast = (_Formula *)f;

  theFormula.Duplicate(&f_cast->theFormula);
  theStack.theStack.Duplicate(&f_cast->theStack.theStack);

  if (f_cast->theTree) {
    theTree = f_cast->theTree->duplicate_tree();
  } else {
    theTree = nil;
  }

  if (f_cast->resultCache) {
    resultCache = (_List *)f_cast->resultCache->makeDynamic();
  } else {
    resultCache = nil;
  }
}

//______________________________________________________________________________
void _Formula::DuplicateReference(const _Formula *f) {
  for (unsigned long i = 0; i < f->theFormula.lLength; i++) {
    _Operation *theO = ((_Operation **)f->theFormula.lData)[i];
    if (theO->IsDeferred ()) {
      _Operation * resolvedCopy = new _Operation (*theO);
      resolvedCopy->ResolveDeferredAction(_hyDefaultExecutionContext);
      theFormula.AppendNewInstance(resolvedCopy);
    } else {
      theFormula &&theO;
    }
  }
}

//______________________________________________________________________________
BaseRef _Formula::makeDynamic(void) {
  _Formula *res = new _Formula;
  res->Duplicate((BaseRef) this);
  return (BaseRef) res;
}

//______________________________________________________________________________
_Formula::~_Formula(void) { 
  Clear(); 
}

//______________________________________________________________________________
void _Formula::Clear(void) {
  DumpTree();
  if (resultCache) {
    DeleteObject(resultCache);
  }

  theFormula.Clear();
  //theStack.Clear();
}

//______________________________________________________________________________
BaseRef _Formula::toStr(_List *matchedNames, bool dropTree) {

  ConvertToTree(false);

  _String *result =
      (_String *)new _String((unsigned long) 16, true);

  long savepd = printDigits;
  printDigits = 0;

  if (theTree) { // there is something to do
    internalToStr(*result, theTree, -1, matchedNames);
  } else {
    if (theFormula.lLength) {
      (*result) << "RPN:";
      internalToStr(*result, nil, 0, nil, (_Operation *)(theFormula(0)));
      for (unsigned long k = 1; k < theFormula.lLength; k++) {
        (*result) << ',';
        internalToStr(*result, nil, 0, nil, (_Operation *)(theFormula(k)));
      }
    }
  }

  printDigits = savepd;
  result->Finalize();
  if (theTree && dropTree) {
    theTree->delete_tree();
    delete theTree;
    theTree = nil;
  }
  return result;

}

//______________________________________________________________________________
node<long> *_Formula::DuplicateFormula(node<long> *src, _Formula &tgt) {
  node<long> *resNode = new node<long>;
  checkPointer(resNode);

  tgt.theFormula &&(_Operation *)theFormula(src->in_object);

  resNode->in_object = tgt.theFormula.lLength - 1;

  for (long k = 1; k <= src->get_num_nodes(); k++) {
    resNode->add_node(*DuplicateFormula(src->go_down(k), tgt));
  }

  return resNode;
}

//______________________________________________________________________________
_Formula *_Formula::Differentiate(_String varName, bool bail) {
  long varID = LocateVarByName(varName);

  if (varID < 0) {
    return new _Formula(new _Constant(0.0));
  }

  varID = variableNames.GetXtra(varID);

  _Formula *res = new _Formula();
  ConvertToTree();

  _SimpleList varRefs, dydx;

  {
    _AVLList al(&varRefs);
    ScanFForVariables(al, true, true, true);
    al.ReorderList();
  }

  for (unsigned long k = 0; k < varRefs.lLength; k++) {
    _Variable *thisVar = LocateVar(varRefs.lData[k]);
    _Formula *dYdX;
    if (thisVar->IsIndependent()) {
      dYdX = new _Formula((thisVar->GetName()->Equal(&varName))
                              ? new _Constant(1.0)
                              : new _Constant(0.0));
      dYdX->ConvertToTree();
      dydx << (long) dYdX;
    } else {
      dYdX = thisVar->varFormula->Differentiate(varName);
      if (dYdX->theFormula.lLength == 0) {
        delete (dYdX);
        return res;
      }
      dydx << (long) dYdX;
    }
  }

  SortLists(&varRefs, &dydx);
  node<long> *dTree = nil;

  if (!(dTree = InternalDifferentiate(theTree, varID, varRefs, dydx, *res))) {
    for (unsigned long k = 0; k < dydx.lLength; k++) {
      delete ((_Formula *)dydx.lData[k]);
    }

    if (bail) {
      WarnError(_String("Differentiation of ") & _String((_String *)toStr()) &
                " failed.");
      res->Clear();
      return res;
    } else {
      delete res;
      return nil;
    }
  }

  for (unsigned long k = 0; k < dydx.lLength; k++) {
    delete ((_Formula *)dydx.lData[k]);
  }

  res->theFormula.AppendNewInstance(new _Operation(new _Constant(0.0)));
  res->theTree = dTree;
  res->InternalSimplify(dTree);
  res->ConvertFromTree();
  return res;

}

//______________________________________________________________________________
// returns true if the subexpression at
// and below startnode is constant

bool _Formula::InternalSimplify(node<long> *startNode) {

  long numChildren = startNode->get_num_nodes(), k, collapse2 = -1;
  bool isConstant = true, firstConst = true, secondConst = (numChildren > 1);

  _Parameter theVal = 0.0;

  _PMathObj newVal = nil;

  _Operation *op = (_Operation *)theFormula(startNode->get_data());

  if (numChildren == 0) {
    return !op->IsAVariable();
  }

  for (k = 1; k <= numChildren; k++) {
    InternalSimplify(startNode->go_down(k));
    if (k == 1) {
      firstConst = InternalSimplify(startNode->go_down(k));
    } else if (k == 2) {
      secondConst = InternalSimplify(startNode->go_down(k));
    } else {
      isConstant = isConstant && InternalSimplify(startNode->go_down(k));
    }
  }

  isConstant = isConstant && firstConst && (numChildren == 1 || secondConst);

  if (op->IsComputed()) {
    if (isConstant) {
      // this executes the subxpression starting at the current
      // node
      if (op->IsVolatileOp()) {
        isConstant = false;
      } else {
        _Stack scrap;
        for (k = 1; k <= numChildren; k++) {
          ((_Operation *)theFormula(startNode->go_down(k)->get_data()))
              ->Execute(scrap);
        }
        op->Execute(scrap);
        newVal = (_PMathObj) scrap.Pop();
      }
    } else {
      if (op->GetOpKind () == _HY_OPERATION_BUILTIN && (firstConst || secondConst)) {
        theVal =
            ((_Operation *)theFormula(startNode->go_down(firstConst ? 1 : 2)
                                          ->get_data()))->GetPayload()->Value();

        switch (op->GetReference()) {
        case HY_OP_CODE_MUL: {           // *
          if (CheckEqual(theVal, 0.0)) { // ?*0 => 0
            newVal = new _Constant(0.0);
            break;
          }
          if (CheckEqual(theVal, 1.0)) { // ?*1 => ?
            collapse2 = firstConst ? 2 : 1;
            break;
          }
        } break;

        case HY_OP_CODE_ADD: {           // +
          if (CheckEqual(theVal, 0.0)) { // ?+0 => ?
            collapse2 = firstConst ? 2 : 1;
            break;
          }
        }

        case HY_OP_CODE_SUB: { // -
          if (CheckEqual(theVal, 0.0)) {
            collapse2 = firstConst ? (-2) : 1;
            break;
          }
        }

        case HY_OP_CODE_DIV: {                         // /
          if (firstConst && CheckEqual(theVal, 0.0)) { // 0/? => 0
            newVal = new _Constant(0.0);
            break;
          }
          if (secondConst && CheckEqual(theVal, 1.0)) { // ?/1 => ?
            collapse2 = 1;
            break;
          }
        } break;

        case HY_OP_CODE_POWER: {                       // ^
          if (firstConst && CheckEqual(theVal, 1.0)) { // 1^? => 1
            newVal = new _Constant(1.0);
            break;
          }
          if (secondConst && CheckEqual(theVal, 1.0)) { // ?^1 => ?
            collapse2 = 1;
            break;
          }
        } break;
        }
      }
    }
  }

  if (newVal) {
    for (k = numChildren; k; k--) {
      startNode->go_down(k)->delete_tree(true);
      startNode->kill_node(k);
    }
    startNode->in_object = theFormula.lLength;
    theFormula.AppendNewInstance(new _Operation(newVal));
  }

  if (collapse2 != -1) {
    if (collapse2 > 0) {
      k = 3 - collapse2;

      startNode->go_down(k)->delete_tree(true);
      node<long> *replaceWith = startNode->go_down(collapse2);

      startNode->kill_node(1);
      startNode->kill_node(2);

      for (k = 1; k <= replaceWith->get_num_nodes(); k++) {
        startNode->add_node(*replaceWith->go_down(k));
      }

      startNode->in_object = replaceWith->in_object;

      delete (replaceWith);

    } else { 
      // 0-? => -?
      delete (startNode->go_down(1));
      startNode->kill_node(1);
    }
  }
  return isConstant;
}

  //______________________________________________________________________________
void _Formula::internalToStr(_String &result, node<long> *currentNode,
                             long opLevel, _List *matchNames,
                             _Operation *thisNodeOperation) {
  
  if (!thisNodeOperation) {
    thisNodeOperation = (_Operation *)theFormula(currentNode->get_data());
  }
  
    // decide what to do about this operation
  if (thisNodeOperation->IsAVariable(false)) {
      // this operation is just a variable - add ident to string and return
    if (subNumericValues) {
      if (subNumericValues == 2) {
        _Variable *theV = LocateVar(thisNodeOperation->GetAVariable());
        if (_x_ && (theV->GetAVariable() == _x_->GetAVariable())) {
          result << _x_->GetName();
          return;
        }
      }
      
      _Variable *thisVariable = LocateVar(thisNodeOperation->GetAVariable());
      _PMathObj subThisValue = thisVariable->Compute();
      
      if (subThisValue->ObjectClass() == NUMBER) {
        if (subNumericValues == 3) {
          result << LocateVar(thisNodeOperation->GetAVariable())->GetName();
          result << '[';
          result.AppendNewInstance(new _String(subThisValue->Value()));
          result << ':';
          result.AppendNewInstance(new _String(thisVariable->GetLowerBound()));
          result << '-';
          result.AppendNewInstance(new _String(thisVariable->GetUpperBound()));
          result << ']';
        } else {
          result.AppendNewInstance(new _String(subThisValue->Value()));
        }
      } else if (subThisValue->ObjectClass() == STRING) {
        result.AppendNewInstance((_String *)subThisValue->toStr());
      } else {
        result << LocateVar(thisNodeOperation->GetAVariable())->GetName();
      }
    } else {
      long variableIDX = thisNodeOperation->GetAVariable();
      if (variableIDX >= 0) {
        _String *vName = LocateVar(variableIDX)->GetName();
        if (matchNames) {
          _List *p1 = (_List *)(*matchNames)(0),
          *p2 = (_List *)(*matchNames)(1);
          
          long f = p1->Find(vName);
          
          if (f < 0) {
            result << vName;
          } else {
            result << (_String *)(*p2)(f);
          }
        } else {
          result << vName;
        }
      }
    }
    return;
  }
  
  switch (thisNodeOperation->GetOpKind()) {

    case _HY_OPERATION_BUILTIN: {
        // check if it's a built-in binary operation
      long precedence_level = thisNodeOperation->GetOperationPrecedence();
      
      if (precedence_level > _HY_OPERATION_MIN_PRECEDENCE && precedence_level < _HY_OPERATION_MAX_PRECEDENCE) {
          // indeed - a binary operation is what we have. check if need to wrap
          // the return in parentheses
        if (!currentNode || currentNode->get_num_nodes() == 2) {
          
          char  tOpLevel  = precedence_level, 
          tOpLevel2 = tOpLevel;
          
          if (!thisNodeOperation->IsAssociativeOp()) {
            tOpLevel2++;
          }
          
          if (opLevel > _HY_OPERATION_MIN_PRECEDENCE) { 
            
              // need to worry about op's precedence
            bool parens = precedence_level < opLevel;
            
            if (parens && currentNode) { 
                // put parentheses around the return of this expression
              result << '(';
              internalToStr(result, currentNode->go_down(1), tOpLevel, matchNames);
              result << thisNodeOperation->GetCode();
              internalToStr(result, currentNode->go_down(2), tOpLevel2, matchNames);
              result << ')';
              return;
            }
          }
          if (currentNode) {
            internalToStr(result, currentNode->go_down(1), tOpLevel, matchNames);
          }
          result << thisNodeOperation->GetCode();
          if (currentNode) {
            internalToStr(result, currentNode->go_down(2), tOpLevel2, matchNames);
          }
          return;
        } else { 
            // mixed binary-unary operation
          result << thisNodeOperation->GetCode();
          if (currentNode) {
            result << '(';
            internalToStr(result, currentNode->go_down(1), _HY_OPERATION_MIN_PRECEDENCE,
                          matchNames);
            result << ')';
          }
          return;
        }
      } else {
        long nOps = thisNodeOperation->GetAttribute();
        if (thisNodeOperation->GetReference() != HY_OP_CODE_MACCESS ) {
          result << thisNodeOperation->GetCode();
          if (currentNode) {
            result << '(';
            for (long k = 1; k <= nOps; k++) {
              if (k > 1) {
                result << ',';
              }
              internalToStr(result, currentNode->go_down(k), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
            }
            result << ')';
          }
        } else { 
            // matrix element access - treat specially
          if (currentNode) {
            internalToStr(result, currentNode->go_down(1), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
            for (long k = 2; k <= nOps; k++) {
              result << '[';
              internalToStr(result, currentNode->go_down(k), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
              result << ']';
            }
          }
        }
      }
      return;
    }
    case _HY_OPERATION_FUNCTION_CALL: {
      result << *_HBLObjectNameByType(HY_BL_HBL_FUNCTION, thisNodeOperation->GetReference());
      if (currentNode) {
        result << '(';
        for (long k = 1; k <= thisNodeOperation->GetAttribute(); k++) {
          if (k > 1) {
            result << ',';
          }
          internalToStr(result, currentNode->go_down(k), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
        }
        result << ')';
      }
      return;
    }
      
    case _HY_OPERATION_ASSIGNMENT_VALUE:
    case _HY_OPERATION_ASSIGNMENT_EXPRESSION:
    case _HY_OPERATION_ASSIGNMENT_UPPER_BOUND:
    case _HY_OPERATION_ASSIGNMENT_LOWER_BOUND: {

        long op_count = currentNode->get_num_nodes();
        printf("%d", op_count);
        if (currentNode) {
          long up_to;
          if(thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_VALUE
             || thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_EXPRESSION) {
            up_to = op_count - (thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_VALUE);
          } else {
            up_to = op_count - (thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_LOWER_BOUND);
          }

          if (up_to == 1) {
            internalToStr(result, currentNode->go_down(1), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
          } else {
            char delims [4] = {'(',')','[',']'};
            for (long k = 1; k <= up_to; k++) {
              result << delims[2*(k>1)];
              internalToStr(result, currentNode->go_down(k), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
              result << delims[2*(k>1)+1];
            }
          }
          if (thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_VALUE) {
            result << '=';
            internalToStr(result, currentNode->go_down(op_count), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
          } else if (thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_LOWER_BOUND 
                     || thisNodeOperation->GetOpKind() == _HY_OPERATION_ASSIGNMENT_UPPER_BOUND) {
            result << ":<";
            internalToStr(result, currentNode->go_down(op_count), _HY_OPERATION_MIN_PRECEDENCE, matchNames);
          } else {
            result << ":=";
            result << (_String*)((_Formula*)thisNodeOperation->GetPayload())->toStr();
          }
        
        }
      }
      return;
      
    default: {
      _PMathObj opValue = thisNodeOperation->GetPayload();
      if (opValue) {
        _String *conv = (_String *)opValue->toStr();
        if (opValue->ObjectClass() == STRING) {
          result << '"';
          result.AppendNewInstance(conv);
          result << '"';
        } else {
          if (opValue->ObjectClass() == NUMBER && opValue->Value() < 0.0) {
            result << '(';
            result.AppendNewInstance(conv);
            result << ')';
          } else {
            result.AppendNewInstance(conv);
          }
        }
      }  
    }
  }
}

//______________________________________________________________________________
bool _Formula::IsEmpty(void) {
  // is there anything in the formula
  return bool(!theFormula.lLength);
}

//______________________________________________________________________________
_Parameter _Formula::Newton(_Formula &derivative, _Variable *unknown,
                            _Parameter targetValue, _Parameter left,
                            _Parameter right) {
  // find a root of the formulaic expression, using Newton's method, given the
  // derivative and a bracketed root.
  // will alternate between bisections and Newton iterations based on what is
  // fatser
  // check that there is indeed a sign change on the interval

  _Constant dummy;
  _Parameter t1, t2, t3, t4, t5, lastCorrection = 100., newCorrection;
  _String msg;
  long iterCount = 0;
  dummy.SetValue(left);
  unknown->SetValue(&dummy);
  t1 = Compute()->Value() - targetValue;
  if (t1 == 0.0) {
    return left;
  }
  dummy.SetValue(right);
  unknown->SetValue(&dummy);
  t2 = Compute()->Value() - targetValue;
  if (t2 == 0) {
    return right;
  }
  if (t1 * t2 > 0.0) {
    subNumericValues = 2;
    _String *s = (_String *)toStr();
    subNumericValues = 0;
    _String msg = *s & "=" & _String(targetValue) &
                  " has no (or multiple) roots in [" & _String(left) & ",Inf)";
    ReportWarning(msg);
    DeleteObject(s);
    return left;
  }

  // else all is good we can start the machine
  bool useNewton = false;
  while ((fabs((right - left) / MAX(left, right)) > machineEps * 10.) && (iterCount < 200)) {
    // stuff to do
    iterCount++;
    if (!useNewton) {
      t3 = (right + left) / 2;
    }
    dummy.SetValue(t3);
    unknown->SetValue(&dummy);
    t4 = Compute()->Value() - targetValue;

    // get the correction term from the derivative
    dummy.SetValue(t3);
    unknown->SetValue(&dummy);
    t5 = derivative.Compute()->Value();
    useNewton = true;
    if (t5 == 0.0) {
      useNewton = false;
    } else {
      newCorrection = -t4 / t5;
      if (t3) {
        if (fabs(newCorrection / t3) < machineEps * 2.) { 
          // correction too small - the root has been found
          return t3;
        }
      } else if (fabs(newCorrection) < machineEps * 2.) { 
        // correction too small - the root has been found
        return t3;
      }
      if (fabs(newCorrection / lastCorrection) > 4) { 
        // Newton correction too large - revert to bisection
        useNewton = false;
      }
      t5 = t3 + newCorrection;
      if ((t5 <= left) || (t5 >= right)) {
        useNewton = false;
      } else {
        lastCorrection = newCorrection;
      }
    }
    if (useNewton) {
      t3 = t5;
    } else {
      dummy.SetValue(t3);
      unknown->SetValue(&dummy);
      t4 = Compute()->Value() - targetValue;
      if (t4 == 0.0) {
        return t3;
      }
      if (t4 * t1 > 0) {
        t1 = t4;
        left = t3;
      } else {
        t2 = t4;
        right = t3;
      }
    }

  }
  return t3;
}

//______________________________________________________________________________
// find a root of the formulaic expression, using Brent's method and a
// bracketed root.
_Parameter _Formula::Brent(_Variable *unknown, _Parameter a, _Parameter b,
                           _Parameter tol, _List *storeEvals, _Parameter rhs) {

  // check that there is indeed a sign change on the interval
  _Constant dummy;

  _Parameter fa = 0.0, fb = 0.0, fc, d, e, min1, min2, xm, p, q, r, s, tol1,
             c = b;

  min1 = unknown->GetLowerBound();
  min2 = unknown->GetUpperBound();

  long it = 0;

  if (a > b) { 
    // autobracket to the left
    dummy.SetValue(b);
    unknown->SetValue(&dummy);
    fb = Compute()->Value();
    if (storeEvals) {
      (*storeEvals) && &dummy;
      dummy.SetValue(fb);
      (*storeEvals) && &dummy;
    }

    if (b < 0.00001 && b > -0.00001) {
      a = b - 0.0001;
    } else {
      a = b - fabs(b) * 0.1;
    }

    if (a < min1) {
      a = min1 + 0.5 * (b - min1);
    }

    dummy.SetValue(a);
    unknown->SetValue(&dummy);
    fa = Compute()->Value() - rhs;
    if (storeEvals) {
      (*storeEvals) && &dummy;
      dummy.SetValue(fa);
      (*storeEvals) && &dummy;
    }

    for (long k = 0; k < 50; k++) {
      if (fb * fa < 0.0) {
        break;
      }

      d = (b - a) * GOLDEN_RATIO;
      b = a;
      a -= d;

      if (a < min1) {
        if (b > min1) {
          a = min1;
        } else {
          break;
        }
      }

      fb = fa;

      dummy.SetValue(a);
      unknown->SetValue(&dummy);
      fa = Compute()->Value() - rhs;
      if (storeEvals) {
        (*storeEvals) && &dummy;
        dummy.SetValue(fa);
        (*storeEvals) && &dummy;
      }
    }
  } else if (CheckEqual(a, b)) { 
    // autobracket to the right
    dummy.SetValue(a);
    unknown->SetValue(&dummy);
    fa = Compute()->Value() - rhs;

    if (storeEvals) {
      (*storeEvals) && &dummy;
      dummy.SetValue(fa);
      (*storeEvals) && &dummy;
    }
    a = b;

    if ((b < 0.00001) && (b > -0.00001)) {
      b = b + 0.0001;
    } else {
      b = b + fabs(b) * 0.1;
    }

    if (b > min2) {
      b = a + 0.5 * (min2 - a);
    }

    dummy.SetValue(b);
    unknown->SetValue(&dummy);
    fb = Compute()->Value() - rhs;

    if (storeEvals) {
      (*storeEvals) && &dummy;
      dummy.SetValue(fb);
      (*storeEvals) && &dummy;
    }

    for (long k = 0; k < 50; k++) {
      if (fb * fa < 0.0) {
        break;
      }

      d = (b - a) * GOLDEN_RATIO;
      a = b;
      b += d;

      if (b > min2) {
        if (a < min2) {
          b = min2;
        } else {
          break;
        }
      }

      fa = fb;

      dummy.SetValue(b);
      unknown->SetValue(&dummy);
      fb = Compute()->Value() - rhs;
      if (storeEvals) {
        (*storeEvals) && &dummy;
        dummy.SetValue(fb);
        (*storeEvals) && &dummy;
      }
    }
  }

  if (fa == 0.0) {
    dummy.SetValue(a);
    unknown->SetValue(&dummy);
    fa = Compute()->Value() - rhs;
    if (storeEvals) {
      (*storeEvals) && &dummy;
      dummy.SetValue(fa);
      (*storeEvals) && &dummy;
    }
    if (fa == 0.0) {
      return a;
    }
  }

  if (fb == 0.0) {
    dummy.SetValue(b);
    unknown->SetValue(&dummy);
    fb = Compute()->Value() - rhs;
    if (storeEvals) {
      (*storeEvals) && &dummy;
      dummy.SetValue(fb);
      (*storeEvals) && &dummy;
    }
    if (fb == 0) {
      return b;
    }
  }

  if (fa * fb < 0.0) {
    fc = fb;

    for (it = 0; it < MAX_BRENT_ITERATES; it++) {
      if (fb * fc > 0.0) {
        fc = fa;
        c = a;
        e = d = b - a;
      }

      if (fabs(fc) < fabs(fb)) {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
      }

      tol1 = 2. * fabs(b) * machineEps + .5 * tol;

      xm = .5 * (c - b);

      if (fabs(xm) <= tol1 || fb == 0.0) {
        return b;
      }

      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
        s = fb / fa;
        if (a == c) {
          p = 2. * xm * s;
          q = 1. - s;
        } else {
          q = fa / fc;
          r = fb / fc;
          p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.0));
          q = (q - 1.) * (r - 1.) * (s - 1.);
        }
        if (p > 0.0) {
          q = -q;
        }

        if (p < 0.0) {
          p = -p;
        }

        min1 = 3. * xm * q - fabs(tol1 * q);
        min2 = fabs(e * q);
        if (2. * p < (min1 < min2 ? min1 : min2)) {
          e = d;
          d = p / q;
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
      if (fabs(d) > tol1) {
        b += d;
      } else {
        b += fabs(tol1) * (1. * 2. * (xm >= .0));
      }

      dummy.SetValue(b);
      unknown->SetValue(&dummy);
      fb = Compute()->Value() - rhs;
      if (storeEvals) {
        (*storeEvals) && &dummy;
        dummy.SetValue(fb);
        (*storeEvals) && &dummy;
      }
    }
  }

  subNumericValues = 2;
  _String msg((_String *)toStr());
  subNumericValues = 0;
  msg = msg & "=" & rhs;
  if (it < MAX_BRENT_ITERATES) {
    msg = msg & " has no (or multiple) roots in [" & _String(a) & "," &
          _String(b) & "]";
  } else {
    msg = msg & " failed to converge to sufficient precision in " &
          MAX_BRENT_ITERATES & " iterates.";
  }

  ReportWarning(msg);
  return b;
}
//______________________________________________________________________________
_Parameter _Formula::Newton(_Formula &derivative, _Parameter targetValue,
                            _Parameter left, _Parameter max_right,
                            _Variable *unknown)
    // given a monotone function and a left bracket bound, found the right
    // bracket bound and solve
    {
  // check that there is indeed a sign change on the interval
  _Constant dummy;
  dummy.SetValue(left);
  unknown->SetValue(&dummy);
  _Parameter t1 = Compute()->Value(), right = left, t2, step = 1.0;

  if (max_right - left < step * 100) {
    step = (max_right - left) / 100;
  }
  if (step == 0.0) {
    return left;
  }
  do {
    right += step;
    if (right > max_right) { // function doesn't seem to have a root
      subNumericValues = 2;
      _String *s = (_String *)toStr();
      subNumericValues = 0;
      _String msg =
          *s & "=" & _String(targetValue) & " has no (or multiple) roots in [" &
          _String(left) & "," & _String(right) & ")";
      ReportWarning(msg);
      DeleteObject(s);
      return left;
    }
    dummy.SetValue(right);
    unknown->SetValue(&dummy);
    t2 = Compute()->Value();
    step *= 2;
    if (right + step > max_right)
      if (right < max_right) {
        step = max_right - right;
      }
  } while ((targetValue - t1) * (targetValue - t2) > 0);
  return Newton(derivative, unknown, targetValue, left, right);

}

//______________________________________________________________________________
_Parameter _Formula::Newton(_Variable *unknown, _Parameter targetValue,
                            _Parameter x_min, _Parameter left,
                            _Parameter right) {

  // check that there is indeed a sign change on the interval
  _Constant dummy;
  _Parameter t1, t2, t3, t4, t5, lastCorrection = 100, newCorrection;
  _String msg;
  t1 = Integral(unknown, x_min, left) - targetValue;
  if (t1 == 0.0) {
    return left;
  }
  t2 = t1 + Integral(unknown, left, right);
  if (t2 == 0) {
    return right;
  }
  if (t1 * t2 > 0.0) {
    subNumericValues = 2;
    _String *s = (_String *)toStr();
    subNumericValues = 0;
    _String msg =
        *s & "=" & _String(targetValue) & " has no (or multiple) roots in [" &
        _String(left) & "," & _String(right) & "]";
    ReportWarning(msg);
    DeleteObject(s);
    return left;
  }

  // else all is good we can start the machine
  bool useNewton = false;
  while (right - left > 1e-6) { 
    // stuff to do
    if (!useNewton) {
      t3 = (right + left) / 2;
    }
    dummy.SetValue(t3);
    unknown->SetValue(&dummy);
    t4 = Integral(unknown, x_min, t3) - targetValue;
    // get the correction term from the derivative
    dummy.SetValue(t3);
    unknown->SetValue(&dummy);
    t5 = Compute()->Value();
    useNewton = true;
    if (t5 == 0.0) {
      useNewton = false;
    } else {
      newCorrection = -t4 / t5;
      if (fabs(newCorrection) < 1e-5) { 
        // correction too small - the root has been found
        return t3;
      }
      if (fabs(newCorrection / lastCorrection) > 4) { 
        // Newton correction too large - revert to bisection
        useNewton = false;
      }
      t5 = t3 + newCorrection;
      if ((t5 <= left) || (t5 >= right)) {
        useNewton = false;
      } else {
        lastCorrection = newCorrection;
      }
    }
    if (useNewton) {
      t3 = t5;
    } else {
      t4 = Integral(unknown, x_min, t3) - targetValue;
      if (t4 == 0.0) {
        return t3;
      }
      if (t4 * t1 > 0) {
        t1 = t4;
        left = t3;
      } else {
        t2 = t4;
        right = t3;
      }
    }

  }
  return t3;
}

//______________________________________________________________________________
// given a monotone function and a left bracket bound, found the right
// bracket bound and solve
_Parameter _Formula::Newton(_Variable *unknown, _Parameter targetValue,
                            _Parameter x_min, _Parameter left) {

  // check that there is indeed a sign change on the interval
  _Parameter t1 = Integral(unknown, x_min, left), right = left, t2, step = 1.0;
  do {
    right += step;
    t2 = Integral(unknown, right - step, right);
    step *= 2;
    if (right >= 1e10) { 
      // function doesn't seem to have a root
      subNumericValues = 2;
      _String *s = (_String *)toStr();
      subNumericValues = 0;
      _String msg =
          *s & "=" & _String(targetValue) & " has no (or multiple) roots in [" &
          _String(left) & ",Inf)";
      WarnError(msg);
      DeleteObject(s);
      return 0.0;
    }
  } while ((targetValue - t1) * (targetValue - t2 - t1) >= 0);
  return Newton(unknown, targetValue, x_min, left, right);

}

// uses Romberg's intergation method
_Parameter _Formula::Integral(_Variable *dx, _Parameter left, _Parameter right,
                              bool infinite) {
  if (infinite) { 
    // tweak "right" here
    _Parameter value = 1.0, step = 1.0, right1 = -1;
    right = left;
    while (value > machineEps) {
      right += step;
      _Constant dummy(right);
      dx->SetValue(&dummy);
      value = fabs(Compute()->Value());
      if ((value < 0.001) && (right1 < 0)) {
        right1 = right;
      }
      step *= 2;
      if (step > 100000) { 
        // integrand decreasing too slowly
        _String msg, *s = (_String *)toStr();
        msg = *s & " decreases too slowly to be integrated over an infinite interval";
        DeleteObject(s);
        WarnError(msg);
        return 0.0;
      }
    }
    if (right1 < right - machineEps) {
      return Integral(dx, left, right1, false) +
             Integral(dx, right1, right, false);
    } else {
      return Integral(dx, left, right1, false);
    }
  }

  checkParameter(intPrecFact, integrationPrecisionFactor, 1e-4);
  checkParameter(intMaxIter, maxRombergSteps, 8);

  _Parameter ss, dss, *s, *h;

  s = new _Parameter[(long) maxRombergSteps];
  h = new _Parameter[(long)(maxRombergSteps + 1)];
  checkPointer(s);
  checkPointer(h);

  h[0] = 1.0;

  long interpolateSteps = 5, stackDepth = 0;

  _SimpleList fvidx, changingVars, idxMap;

  _Parameter *ic = new _Parameter[interpolateSteps],
             *id = new _Parameter[interpolateSteps];

  _SimpleFormulaDatum *stack = nil, *vvals = nil;

  checkPointer(ic);
  checkPointer(id);

  if (AmISimple(stackDepth, fvidx)) {
    stack = new _SimpleFormulaDatum[stackDepth];
    checkPointer(stack);
    vvals = new _SimpleFormulaDatum[fvidx.lLength];
    checkPointer(vvals);
    ConvertToSimple(fvidx);
    stackDepth = dx->GetAVariable();
    for (unsigned long vi = 0; vi < fvidx.lLength; vi++) {
      _Variable *checkvar = LocateVar(fvidx.lData[vi]);
      if (checkvar->CheckFForDependence(stackDepth, true)) {
        changingVars << fvidx.lData[vi];
        idxMap << vi;
      }
      vvals[vi].value = checkvar->Compute()->Value();
    }
    changingVars.InsertElement((BaseRef) stackDepth, 0, false, false);
    idxMap.InsertElement((BaseRef) fvidx.Find(stackDepth), 0, false, false);
  } else {
    stackDepth = -1;
  }

  for (long j = 0; j < (long) maxRombergSteps; j++) {
    if (stackDepth >= 0) {
      s[j] = TrapezoidLevelKSimple(*this, dx, left, right, j + 1, stack, vvals,
                                   changingVars, idxMap);
    } else {
      s[j] = TrapezoidLevelK(*this, dx, left, right, j + 1);
    }
    if (j >= 4) {
      ss = InterpolateValue(&h[j - 4], &s[j - 4], interpolateSteps, ic, id, 0.0,
                            dss);
      if (fabs(dss) <= integrationPrecisionFactor * fabs(ss)) {
        delete s;
        delete h;
        delete ic;
        delete id;
        if (stackDepth >= 0) {
          ConvertFromSimple();
          delete (stack);
          delete (vvals);
        }
        return ss;
      }
    }
    h[j + 1] = h[j] / 9.0;
  }

  if (stackDepth >= 0) {
    ConvertFromSimple();
    delete (stack);
    delete (vvals);
  }

  _String *str = (_String *)toStr(),
          msg =
              _String("Integral of ") & *str & " over [" & _String(left) & "," &
              _String(right) &
              "] converges slowly, loss of precision may occur. Change either "
              "INTEGRATION_PRECISION_FACTOR or INTEGRATION_MAX_ITERATES";

  DeleteObject(str);
  ReportWarning(msg);

  delete s;
  delete h;
  delete ic;
  delete id;
  return ss;
}

//______________________________________________________________________________
_Parameter _Formula::MeanIntegral(_Variable *dx, _Parameter left,
                                  _Parameter right, bool infinite) {
  _Formula newF;
  newF.Duplicate((BaseRef) this);
  newF.theFormula.AppendNewInstance (new _Operation( *(dx->GetName()), true, 
    dx->IsGlobal()));
  newF.theFormula.AppendNewInstance (new _Operation( HY_OP_CODE_MUL, 2));
  return newF.Integral(dx, left, right, infinite);
}

//______________________________________________________________________________
// number of operations in the formula
long _Formula::NumberOperations(void) {
  return theFormula.lLength;
}

//______________________________________________________________________________
long _Formula::ExtractMatrixExpArguments(_List *storage) {
  long count = 0;

  if (resultCache && resultCache->lLength) {
    long cacheID = 0;
    bool cacheUpdated = false;

    // whether or not a cached result was used
    for (unsigned long i = 0; i < theFormula.lLength; i++) {
      _Operation *thisOp((_Operation *)(((BaseRef **)theFormula.lData)[i]));
      if (i < theFormula.lLength - 1) {
        _Operation *nextOp(
            (_Operation *)(((BaseRef **)theFormula.lData)[i + 1]));

        if (!cacheUpdated && nextOp->CanResultsBeCached(thisOp)) {
          _Stack temp;
          thisOp->Execute(temp);

          _Matrix *currentArg = (_Matrix *)temp.Pop(true),
                  *cachedArg = (_Matrix *)((_PMathObj)(*resultCache)(cacheID)),
                  *diff = nil;

          if (cachedArg->ObjectClass() == MATRIX) {
            diff = (_Matrix *)cachedArg->SubObj(currentArg);
          }

          if (diff && diff->MaxElement() <= 1e-12) {
            cacheID += 2;
            i++;
          } else {
            cacheUpdated = true;
            cacheID++;
            if (nextOp->CanResultsBeCached(thisOp, true)) {
              storage->AppendNewInstance(currentArg);
              count++;
            }
          }
          DeleteObject(diff);
          continue;
        }
        if (cacheUpdated) {
          cacheID++;
          cacheUpdated = false;
        }
      }
    }
  }

  return count;
}

//______________________________________________________________________________
_Variable *_Formula::Dereference(bool ignore_context,
                                 _hyExecutionContext *theContext) {
  _Variable *result = nil;
  _PMathObj computedValue =
      Compute(0, theContext);
  if (computedValue && computedValue->ObjectClass() == STRING) {
    result = (_Variable *)((_FString *)computedValue)
        ->Dereference(ignore_context, theContext, true);
  }

  if (!result) {
    theContext->ReportError(
        (_String("Failed to dereference '") & _String((_String *)toStr()) &
         "' in the " & (ignore_context ? "global" : "local") & " context"));
  }

  return result;
}

//______________________________________________________________________________
// compute the value of the formula
_PMathObj _Formula::Compute(long startAt, _hyExecutionContext* execContext,
                            _List *additionalCacheArguments, 
                            unsigned long result_type) {
  if (theFormula.lLength == 0) {
    theStack.theStack.Clear();
    theStack.Push(new _MathObject, false);
  } else {
    bool wellDone = true;
    if (startAt == 0) {
      theStack.theStack.Clear();
    }

    if (startAt == 0 && resultCache && resultCache->lLength) {

      long cacheID = 0;

      // where in the cache are we currently looking
      bool cacheUpdated = false;

      // whether or not a cached result was used
      for (unsigned long i = 0; i < theFormula.lLength; i++) {
        _Operation *thisOp((_Operation *)(((BaseRef **)theFormula.lData)[i]));
        if (i < theFormula.lLength - 1) {
          _Operation *nextOp(
              (_Operation *)(((BaseRef **)theFormula.lData)[i + 1]));

          if (!cacheUpdated && nextOp->CanResultsBeCached(thisOp)) {
            if (!thisOp->Execute(theStack, execContext)) {
              wellDone = false;
              break;
            }

            _Matrix *currentArg = (_Matrix *)theStack.Pop(false),
                    *cachedArg =
                        (_Matrix *)((_PMathObj)(*resultCache)(cacheID)),
                    *diff = nil;

            if (cachedArg->ObjectClass() == MATRIX) {
              diff = (_Matrix *)cachedArg->SubObj(currentArg);
            }

            bool no_difference = diff && diff->MaxElement() <= 1e-12;

            if (no_difference || (additionalCacheArguments &&
                                  additionalCacheArguments->lLength &&
                                  nextOp->CanResultsBeCached(thisOp, true))) {
              DeleteObject(theStack.Pop());
              if (no_difference) {
                theStack.Push((_PMathObj)(*resultCache)(cacheID + 1));
              } else {

                theStack.Push((_PMathObj)(*additionalCacheArguments)(0));
                resultCache->Replace(cacheID, theStack.Pop(false), true);
                resultCache->Replace(cacheID + 1,
                                     (*additionalCacheArguments)(0), false);
                additionalCacheArguments->Delete(0, false);
                //printf ("_Formula::Compute additional arguments %ld\n",
                //additionalCacheArguments->lLength);
              }
              cacheID += 2;
              i++;
              //printf ("Used formula cache %s\n",
              //_String((_String*)nextOp->toStr()).sData);
            } else {
              cacheUpdated = true;
              resultCache->Replace(cacheID++, theStack.Pop(false), true);
              //printf ("Updated formula cache %s\n",
              //_String((_String*)nextOp->toStr()).sData);
            }
            DeleteObject(diff);
            continue;
          }
        }

        if (!thisOp->Execute(theStack, execContext)) { 
          // does this always get executed?
          wellDone = false;
          break;
        }

        if (cacheUpdated) {
          resultCache->Replace(cacheID++, theStack.Pop(false), true);
          cacheUpdated = false;
        }
      }

    } else {
      for (unsigned long i = startAt; i < theFormula.lLength; i++)
        if (!((_Operation *)(((BaseRef **)theFormula.lData)[i]))
                ->Execute(theStack, execContext)) {
          wellDone = false;
          break;
        }
    }
    if (theStack.theStack.lLength != 1 || !wellDone) {
      _String errorText =
          _String((_String *)toStr()) & _String(" contains errors.");
      if (execContext->GetErrorBuffer()) {
        *execContext->GetErrorBuffer() = *execContext->GetErrorBuffer() & errorText;
      } else {
        WarnError(errorText);
      }
      theStack.theStack.Clear();
      theStack.Push(new _Constant(0.0), false);
    }
  }

  _PMathObj r = theStack.Pop(false);
  if (result_type == HY_ANY_OBJECT || (r->ObjectClass() & result_type) != 0) {
    return r;
  }
  return nil;
}

//______________________________________________________________________________
bool _Formula::CheckSimpleTerm(_PMathObj thisObj) {
  if (thisObj) {
    long oc = thisObj->ObjectClass();
    if (oc != NUMBER) {
      if (oc == MATRIX) {
        _Matrix *mv = (_Matrix *)thisObj->Compute();
        if (mv->IsIndependent() && !mv->SparseDataStructure()) {
          return true;
        }
      }
    } else {
      return true;
    }
  }
  return false;
}

//______________________________________________________________________________
_Formula _Formula::operator+(const _Formula &operand2) {
  _Formula joint;
  return joint.PatchFormulasTogether(*this, operand2, HY_OP_CODE_ADD);
}

//______________________________________________________________________________
_Formula _Formula::operator-(const _Formula &operand2) {
  _Formula joint;
  return joint.PatchFormulasTogether(*this, operand2, HY_OP_CODE_SUB);
}

//______________________________________________________________________________
_Formula _Formula::operator*(const _Formula &operand2) {
  _Formula joint;
  return joint.PatchFormulasTogether(*this, operand2, HY_OP_CODE_MUL);
}

//______________________________________________________________________________
_Formula _Formula::operator/(const _Formula &operand2) {
  _Formula joint;
  return joint.PatchFormulasTogether(*this, operand2, HY_OP_CODE_DIV);
}

//______________________________________________________________________________
_Formula _Formula::operator^(const _Formula &operand2) {
  _Formula joint;
  return joint.PatchFormulasTogether(*this, operand2, HY_OP_CODE_POWER);
}

//______________________________________________________________________________
_Formula &_Formula::PatchFormulasTogether(_Formula &target,
                                          const _Formula &operand2,
                                          const char op_code) {
  target.Clear();
  target.DuplicateReference(this);
  target.DuplicateReference(&operand2);
  target.theFormula.AppendNewInstance(new _Operation( op_code, 2));
  return target;
}

//______________________________________________________________________________
void _Formula::ConvertMatrixArgumentsToSimpleOrComplexForm(bool makeComplex) {
  if (makeComplex) {
    if (resultCache) {
      DeleteObject(resultCache);
      resultCache = nil;
    }
  } else {
    if (!resultCache) {
      resultCache = new _List();
      for (unsigned long i = 1; i < theFormula.lLength; i++) {
        _Operation *thisOp =
            ((_Operation *)(((BaseRef **)theFormula.lData)[i]));
        if (thisOp->CanResultsBeCached(
                ((_Operation *)(((BaseRef **)theFormula.lData)[i - 1])))) {
          resultCache->AppendNewInstance(new _MathObject());
          resultCache->AppendNewInstance(new _MathObject());
        }
      }
    }
  }

  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    _PMathObj thisMatrix = ((_Operation *)(((BaseRef **)theFormula.lData)[i]))->FetchReferencedObject();

    if (thisMatrix && thisMatrix->ObjectClass() == MATRIX)
      if (makeComplex) {
        ((_Matrix*)thisMatrix)->MakeMeGeneral();
      } else {
        ((_Matrix*)thisMatrix)->MakeMeSimple();
      }
  }
}

//______________________________________________________________________________

_Operation* _Formula::getIthOp (unsigned long idx) const {
  return ((_Operation *)(((BaseRef **)theFormula.lData)[idx]));
}

//______________________________________________________________________________

long _Formula::PrepareLHS (long idx) {
  if (idx != HY_NOT_FOUND) {
    long reference = ((_Operation *)(((BaseRef **)theFormula.lData)[idx]))->PrepareLHS();
    if (reference != _HY_OPERATION_INVALID_REFERENCE) {
      // delete the MACCESS operation from the formula
      theFormula.Delete (idx, true);
      return reference;
    }
  }
  return _HY_OPERATION_INVALID_REFERENCE;
}


//______________________________________________________________________________
long _Formula::LValueIndex(const unsigned long start_at, const bool handle_matrix_op) const {
  
  _SimpleList    lvalues, indices;
  
  for (unsigned long op = start_at; op < theFormula.lLength; op++) {
    getIthOp (op) -> UpdateLValue (lvalues, indices, op);
  }
  
  for (long code = lvalues.lLength - 1; code >= 0L; code --) {
    if (lvalues.GetElement(code) != HY_NOT_FOUND) {
      // check if we need to toggle the matrix op code
      if (handle_matrix_op) {
        if (code > 0) {
          if (lvalues.GetElement (code-1) == lvalues.GetElement (code)) {
            ((_Operation *)(((BaseRef **)theFormula.lData)[indices.GetElement(code-1)]))->
            ToggleVarRef(true);
          }
        }
      }
      return indices.GetElement(code);
    }
  }
  
  return HY_NOT_FOUND;
}
//______________________________________________________________________________
bool _Formula::AmISimple(long &stackDepth, _SimpleList &variableIndex) {
  if (!theFormula.lLength) {
    return true;
  }
  
  long locDepth = 0;
  
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    _Operation *thisOp = getIthOp (i);
    locDepth++;
    
    switch (thisOp->GetOpKind()) {
      case _HY_OPERATION_VALUE:
        if (thisOp->GetPayload()->ObjectClass() != NUMBER) {
          return false;
        }
        break;
        
      case _HY_OPERATION_VAR:
      case _HY_OPERATION_VAR_OBJ: {
        long var_ref = thisOp->GetReference();
        _Variable *thisVar = LocateVar(var_ref);
        if (thisVar->ObjectClass() != NUMBER) {
          _PMathObj cv = thisVar->GetValue();
          if (!CheckSimpleTerm(cv)) {
            return false;
          }
        }
        if (variableIndex.Find(var_ref) == HY_NOT_FOUND) {
          variableIndex << var_ref;
        }
      }
        break;
        
      case _HY_OPERATION_BUILTIN: {
        long op_code = thisOp->GetReference();
        if (simpleOperationCodes.Find(op_code) == HY_NOT_FOUND) {
          return false;
        } 
        else {
          if ((op_code == HY_OP_CODE_MACCESS || op_code ==
               HY_OP_CODE_MUL) && thisOp->GetAttribute() != 2) {
            return false;
          }
          locDepth -= thisOp->GetAttribute();
        }
      }
        break;
        
      default:
        return false;
        
    }
    
    if (locDepth > stackDepth) {
      stackDepth = locDepth;
    } else if (locDepth == 0) {
      _String errStr =
      _String("Invalid formula passed to _Formula::AmISimple") &
      _String((_String *)toStr());
      WarnError(errStr);
      return false;
    }
  }
  return true;
}

//______________________________________________________________________________
bool _Formula::IsArrayAccess(void) {
  if (theFormula.lLength) {
    const _Operation * theOp = getIthOp(theFormula.lLength - 1);
    return theOp->GetOpKind() == _HY_OPERATION_BUILTIN && theOp->GetReference() == HY_OP_CODE_MACCESS;
  }
  return false;
}

//______________________________________________________________________________
bool _Formula::ConvertToSimple(_SimpleList &variableIndex) {
  bool has_volatiles = false;
  if (theFormula.lLength)
    for (unsigned long i = 0; i < theFormula.lLength; i++) {
      has_volatiles = has_volatiles || ((_Operation *)(((BaseRef **)theFormula.lData)[i]))->ToggleFastExec(true, &variableIndex);
    }
  return has_volatiles;
}

//______________________________________________________________________________
void _Formula::ConvertFromSimple(void) {

  if (!theFormula.lLength) {
    return;
  }

  for (unsigned long i = 0; i < theFormula.lLength; i++) {
      ((_Operation *)(((BaseRef **)theFormula.lData)[i]))->ToggleFastExec(false, NULL);
  }
}

//______________________________________________________________________________
_Parameter _Formula::ComputeSimple(_SimpleFormulaDatum *stack,
                                   _SimpleFormulaDatum *varValues) {
  if (!theFormula.lLength) {
    return 0.0;
  }

  long stackTop = 0;

  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    if (!((_Operation *)(((BaseRef **)theFormula.lData)[i]))->ExecuteFast(stack,varValues,stackTop, NULL)) {
      return 0.0;
    }
  }

  return stack->value;
}

//______________________________________________________________________________
bool _Formula::EqualFormula(_Formula *f) {
  if (theFormula.lLength == f->theFormula.lLength) {
    for (unsigned long i = 0; i < theFormula.lLength; i++) {
      if (!((_Operation *)(((BaseRef **)theFormula.lData)[i]))->EqualOp(
              ((_Operation *)(((BaseRef **)f->theFormula.lData)[i])))) {
        return false;
      }
    }
    return true;
  }
  return false;
}

//______________________________________________________________________________
// compute the value of the formula
_PMathObj _Formula::ConstructPolynomial(void) {
  theStack.Reset();
  bool wellDone = true;
  _String errMsg;

  for (unsigned long i = 0; i < theFormula.lLength; i++)
    if (!((_Operation *)((BaseRef **)theFormula.lData)[i])
            ->ExecutePolynomial(theStack, nil, &errMsg)) {
      wellDone = false;
      break;
    }

  if (theStack.theStack.lLength != 1 || !wellDone) {
    return nil;
  }

  return theStack.Pop(false);
}

//______________________________________________________________________________
bool _Formula::HasChanged(bool ignoreCats) {
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    if (((_Operation *)((BaseRef **)theFormula.lData)[i])->HasChanged(ignoreCats)) {
      return true;
    }
  }
  return false;
}

//______________________________________________________________________________
bool _Formula::HasChangedSimple(_SimpleList &variableIndex) {
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    if (((_Operation *)((BaseRef **)theFormula.lData)[i])->HasChanged(false, &variableIndex)){
      return true;
    } 
  }
  return false;
}

//______________________________________________________________________________
void _Formula::ScanFForVariables(_AVLList &l, bool includeGlobals,
                                 bool includeAll, bool includeCategs,
                                 bool skipMatrixAssignments, _AVLListX *tagger,
                                 long weight) {
  
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    _Operation *theObj = ((_Operation **)theFormula.lData)[i];
    if (theObj->IsAVariable(true)) {
      long f = theObj->GetAVariable();
      
      if (f != _HY_OPERATION_INVALID_REFERENCE) {
        
        if (!includeGlobals) {
            // This change was part of a commit that introduced an optimizer bug
            // (suspected location:
            // src/core/batchlan2.cpp:2220). This change is suspicious as well
            // (removed and undocumented condition).
            //if ((((_Variable*)LocateVar(theObj->GetAVariable()))->IsGlobal())||
            //       (((_Variable*)LocateVar(theObj->GetAVariable()))->ObjectClass()!=NUMBER))
            // {
          if ((_Variable *)LocateVar(f)->IsGlobal()) {
            continue;
          }
        }
        
        
        _Variable *v = LocateVar(f);
        
        if (v->IsCategory() && includeCategs) {
          v->ScanForVariables(l, includeGlobals, tagger, weight);
        }
        
        if (includeAll || v->ObjectClass() == NUMBER) {
          l.Insert((BaseRef) f);
          if (tagger) {
            tagger->UpdateValue((BaseRef) f, weight, 0);
          }
        }
        
        if (skipMatrixAssignments) {
          if (v->ObjectClass() != MATRIX || !theObj->AssignmentVariable()) {
            v->ScanForVariables(l, includeGlobals, tagger, weight);
          }
          
        } else if (!v->IsIndependent()) {
          v->ScanForVariables(l, includeGlobals, tagger);
        }
        
      } else {
        _PMathObj op_load = theObj->GetPayload();
        if (op_load && op_load->ObjectClass() == MATRIX) {
          ((_Matrix *)op_load)->ScanForVariables(l, includeGlobals, tagger, weight);
        }
      }
    }
  }
}

//______________________________________________________________________________
void _Formula::ScanFForType(_SimpleList &l, int type) {
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    _Operation *theObj = ((_Operation **)theFormula.lData)[i];
    if (theObj->IsAVariable()) {
      long f = theObj->GetAVariable();
      if (f >= 0) {
        _Variable *v = LocateVar(f);

        if (v->ObjectClass() == type) {
          l << f;
        }
      }
    }
  }
}

//______________________________________________________________________________
bool _Formula::CheckFForDependence(long varID, bool checkAll) {
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    _Operation *theObj = (_Operation *)theFormula(i);
    long f;
    if (theObj->IsAVariable()) {
      f = theObj->GetAVariable();
      if (f >= 0) {
        if (f == varID) {
          return true;
        }
        if (checkAll) {
          _Variable *v = LocateVar(f);
          if (!v->IsIndependent()) {
            if (v->CheckFForDependence(varID)) {
              return true;
            }
          }
        }
      }
    }
  }
  return false;
}

//______________________________________________________________________________
void _Formula::LocalizeFormula(_Formula &ref, _String &parentName,
                               _SimpleList &iv, _SimpleList &iiv,
                               _SimpleList &dv, _SimpleList &idv) {

  for (unsigned long i = 0; i < ref.theFormula.lLength; i++) {
    if (((_Operation *)ref.theFormula(i))->IsAVariable()) {
      long vIndex = ((_Operation *)ref.theFormula(i))->GetAVariable();
      _Variable *theV = LocateVar(vIndex);

      if (theV->IsGlobal()) {
        theFormula &&ref.theFormula(i);
        continue;
      }

      if (theV->IsContainer()) {
        continue;
      }

      _String fullName = parentName & "." & *(theV->GetName());
      long lvIndex = LocateVarByName(fullName);

      // local variable doesn't yet exist - create
      if (lvIndex == -1) {
        _Variable dummy(fullName);
        lvIndex = LocateVarByName(fullName);
        if (theV->IsIndependent()) {
          iv << lvIndex;
          iiv << vIndex;
        } else {
          dv << lvIndex;
          idv << vIndex;
        }
      }
      theFormula.AppendNewInstance (new _Operation ( fullName, true));
    } else {
      theFormula &&ref.theFormula(i);
    }
  }
}

//______________________________________________________________________________
bool _Formula::DependsOnVariable(long idx) {
  for (unsigned long f = 0; f < theFormula.lLength; f++) {
    _Operation *anOp = ((_Operation **)theFormula.lData)[f];
    if (anOp->IsAVariable() && anOp->GetAVariable() == idx) {
      return true;
    }
  }
  return false;
}

//______________________________________________________________________________
_Operation *_Formula::GetIthTerm(long idx) {
  if (idx >= 0 && idx < theFormula.lLength) {
    return ((_Operation **)theFormula.lData)[idx];
  }
  return nil;
}

//______________________________________________________________________________
bool _Formula::IsAConstant(bool deep) {
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    if (((_Operation *)((BaseRef **)theFormula.lData)[i])->IsAVariable(deep)) {
      return false;
    }
  }
  return true;
}

//______________________________________________________________________________
bool _Formula::IsConstant(void) {
  for (unsigned long i = 0; i < theFormula.lLength; i++) {
    if (((_Operation *)((BaseRef **)theFormula.lData)[i])->IsConstant() ==
        false) {
      return false;
    }
  }
  return true;
}

//______________________________________________________________________________
void _Formula::SimplifyConstants(void) {
  ConvertToTree();
  if (theTree) {
    InternalSimplify (theTree);
    ConvertFromTree (false);
  }
}

//______________________________________________________________________________
_PMathObj _Formula::GetTheMatrix(void) {
  if (theFormula.lLength == 1) {
    _Operation *firstOp = (_Operation *)theFormula(0);
    switch (firstOp->GetOpKind()) {
    
       case _HY_OPERATION_VALUE: {
          _PMathObj ret = firstOp->GetPayload();
          if (ret && (ret->ObjectClass() == MATRIX)) {
            return ret;
          } 
          break;
      }
      
      case _HY_OPERATION_VAR:
      case _HY_OPERATION_VAR_OBJ: {
        _Variable *firstVar = LocateVar(firstOp->GetReference());
        _PMathObj ret = firstVar->GetValue();
        if (ret && (ret->ObjectClass() == MATRIX)) {
          return ret;
        }
        break;
      }
    }
  }
  
  return nil;
}

//______________________________________________________________________________
long _Formula::FormulaType(void) const {
  
  if (theFormula.lLength) {
    const _Operation* last_op = getIthOp(theFormula.lLength-1);
    if (last_op -> GetOpKind() == _HY_OPERATION_ASSIGNMENT_VALUE) {
      return last_op->GetReference() > _HY_OPERATION_INVALID_REFERENCE ?
             HY_FORMULA_FORMULA_VALUE_ASSIGNMENT :
             HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT;
    }
  }
  return HY_FORMULA_EXPRESSION;
}

//______________________________________________________________________________
long _Formula::ObjectClass(void) {
  if (theStack.theStack.lLength) {
    return ((_PMathObj) theStack.theStack.lData[0])->ObjectClass();
  }

  _PMathObj res = Compute();

  if (res) {
    return res->ObjectClass();
  }

  return HY_UNDEFINED;
}

//______________________________________________________________________________
// the parser itself
_Formula::_Formula(_String &s, _VariableContainer *theParent,
                   _String *reportErrors) {

  theTree = nil;
  resultCache = nil;

  _FormulaParsingContext fpc(reportErrors, theParent);

  if (Parse(this, s, fpc) != HY_FORMULA_EXPRESSION) {
    Clear();
  }

}

//______________________________________________________________________________
void _Formula::ConvertToTree(bool err_msg) {
  
  if (!theTree && theFormula.lLength) { 
    
      // work to do
    _SimpleList nodeStack;
    bool can_convert = true;
    
    _Operation *currentOp;
    for (unsigned long i = 0; i < theFormula.lLength; i++) {
      currentOp = (_Operation *)theFormula(i);
      long opKind = currentOp->GetOpKind();
      
      switch (opKind) {
        case _HY_OPERATION_VALUE:
        case _HY_OPERATION_VAR:
        case _HY_OPERATION_VAR_OBJ:
        case _HY_OPERATION_REF:
        case _HY_OPERATION_NOOP:
        case _HY_OPERATION_DEFERRED_INLINE:
        case _HY_OPERATION_DICTIONARY:
        case _HY_OPERATION_SPARSE_MATRIX: {
          node<long> *leafNode = new node<long>;
          leafNode->init(i);
          nodeStack << (long) leafNode;
          break;
        }
          
        case  _HY_OPERATION_FUNCTION_CALL:
        case  _HY_OPERATION_DEFERRED_FUNCTION_CALL:
        case  _HY_OPERATION_BUILTIN:
        case  _HY_OPERATION_ASSIGNMENT_VALUE:
        case  _HY_OPERATION_ASSIGNMENT_EXPRESSION:
        case  _HY_OPERATION_ASSIGNMENT_UPPER_BOUND:
        case  _HY_OPERATION_ASSIGNMENT_LOWER_BOUND: {
          long nTerms = currentOp->OperandCount();
          
          //if (opKind == _HY_OPERATION_ASSIGNMENT_EXPRESSION
          //    || opKind == _HY_OPERATION_ASSIGNMENT_UPPER_BOUND
          //    || opKind == _HY_OPERATION_ASSIGNMENT_LOWER_BOUND ) {
          if (opKind == _HY_OPERATION_ASSIGNMENT_EXPRESSION) {
              // need to push the expression encoded in payload _Formula onto the tree
              nTerms --;
          }
          
          if (nTerms > nodeStack.lLength) {
            if (err_msg) {
              WarnError(
                        _String("Insufficient number of arguments for a call to ") &
                        _String((_String *)currentOp->toStr()));
            }
            theTree = nil;
            return;
          }

          node<long> *operationNode = new node<long>;
          operationNode->init(i);
          for (unsigned long j = 0; j < nTerms; j++) {
            operationNode->prepend_node(
                                        *((node<long> *)nodeStack(nodeStack.lLength - 1)));
            nodeStack.Pop();
          }
          nodeStack << (long) operationNode;
          break;
        }
          
        default:
          can_convert = false;
          break;
      }
    }
    
    if (can_convert == false || nodeStack.lLength != 1) {
      if (err_msg) {
        if (can_convert) {
          WarnError((_String) "The expression '" & _String((_String *)toStr()) &
                    "' has " & (long) nodeStack.lLength &
                    " terms left on the stack after evaluation");
        } else {
          WarnError((_String) "Failed to convert the formula to a parse tree: unsupported operation type");
        }
      }
      theTree = nil;
    } else {
      theTree = (node<long> *)nodeStack(0);
    }
  }
}

//______________________________________________________________________________
void _Formula::DumpTree (void) {
  if (theTree) {
      theTree->delete_tree();
      delete (theTree);
      theTree = nil;        
  }
}


//______________________________________________________________________________
void _Formula::ConvertFromTree(bool rebuild_tree) {
  if (theTree) { 
    // work to do
    _SimpleList termOrder;
    node<long> *currentNode = DepthWiseStepTraverser(theTree);
    while (currentNode) {
      termOrder << currentNode->get_data();
      currentNode = DepthWiseStepTraverser((node<long> *)nil);
    }
    if (termOrder.lLength != theFormula.lLength) { 
      _List newFormula;
      for (unsigned long i = 0; i < termOrder.lLength; i++) {
        newFormula << theFormula(termOrder(i));
      }
      theFormula.Clear();
      theFormula.Duplicate(&newFormula);
      DumpTree();
      if (rebuild_tree) {
        ConvertToTree();
      }
    } else {
      DumpTree();
    }
  }
}

//______________________________________________________________________________
node<long> *_Formula::InternalDifferentiate(node<long> *currentSubExpression,
                                            long varID, _SimpleList &varRefs,
                                            _SimpleList &dydx, _Formula &tgt) {

  _Operation *op = (_Operation *)theFormula(currentSubExpression->in_object);
  
  switch (op->GetOpKind()) {
    case _HY_OPERATION_VALUE: {
      _Formula src(new _Constant(0.0));
      src.ConvertToTree();
      return src.DuplicateFormula(src.theTree, tgt);
    }
      
    case _HY_OPERATION_VAR:
    case _HY_OPERATION_VAR_OBJ: {
      long k = varRefs.BinaryFind(op->GetReference());
      if (k < 0) {
        return nil;
      }
      _Formula *dYdX = (_Formula *)dydx(k);
      return dYdX->DuplicateFormula(dYdX->theTree, tgt);
    }
      
      
      
    case _HY_OPERATION_BUILTIN: {
      node<long> *newNode = (node<long> *)checkPointer(new node<long>);
      
      switch (op->GetReference()) {
        case HY_OP_CODE_MUL: {
          
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
          *b2 = InternalDifferentiate(currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
          
          if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
              b1->delete_tree(true);
            }
            if (b2) {
              b2->delete_tree(true);
            }
            return nil;
          }
          
          
          _Operation *newOp = new _Operation( HY_OP_CODE_ADD, 2),
          *newOp2 = new _Operation( HY_OP_CODE_MUL, 2),
          *newOp3 = new _Operation( HY_OP_CODE_MUL, 2);
          
          
          node<long> *newNode2 = (node<long> *)checkPointer(new node<long>);
          node<long> *newNode3 = (node<long> *)checkPointer(new node<long>);
          
          newNode2->add_node(*b1);
          newNode2->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(2), tgt));
          
          newNode3->add_node(*b2);
          newNode3->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode->add_node(*newNode2);
          newNode->add_node(*newNode3);
          
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp3);
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          return newNode;
          
        } break;
          
        case HY_OP_CODE_ADD:   // +
        case HY_OP_CODE_SUB: { // -
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt),
          *b2 = nil;
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          long isUnary = (currentSubExpression->get_num_nodes() == 1);
          
          if (!isUnary) {
            b2 = InternalDifferentiate(currentSubExpression->go_down(2), varID,
                                       varRefs, dydx, tgt);
            if (!b2) {
              b1->delete_tree(true);
              newNode->delete_tree(true);
              return nil;
            }
          }
          
          _Operation *newOp = new _Operation();
          checkPointer(newOp);
          newOp->Duplicate(op);
          newNode->add_node(*b1);
          if (!isUnary) {
            newNode->add_node(*b2);
          }
          newNode->in_object = tgt.theFormula.lLength;
          
          tgt.theFormula.AppendNewInstance(newOp);
          return newNode;
        } break;
          
        case HY_OP_CODE_DIV: { // /
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt),
          *b2 = InternalDifferentiate(currentSubExpression->go_down(2),
                                      varID, varRefs, dydx, tgt);
          
          if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
              b1->delete_tree(true);
            }
            if (b2) {
              b2->delete_tree(true);
            }
            return nil;
          }
          
          
          _Operation *newOp  = new _Operation( HY_OP_CODE_DIV, 2),
          *newOp2 = new _Operation( HY_OP_CODE_POWER, 2),
          *newOp3 = new _Operation( HY_OP_CODE_SUB, 2),
          *newOp4 = new _Operation( HY_OP_CODE_MUL, 2),
          *newOp5 = new _Operation( HY_OP_CODE_MUL, 2),
          *newOp6 = new _Operation(new _Constant(2.0));
          
          
          node<long> *newNode2 = new node<long>;
          node<long> *newNode3 = new node<long>;
          node<long> *newNode4 = new node<long>;
          node<long> *newNode5 = new node<long>;
          node<long> *newNode6 = new node<long>;
          
          
          newNode6->add_node(*b1);
          newNode6->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(2), tgt));
          
          newNode5->add_node(*b2);
          newNode5->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode4->add_node(*newNode6);
          newNode4->add_node(*newNode5);
          
          newNode2->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(2), tgt));
          newNode2->add_node(*newNode3);
          
          newNode->add_node(*newNode4);
          newNode->add_node(*newNode2);
          
          newNode6->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp5);
          newNode5->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp4);
          newNode4->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp3);
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp6);
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_ARCTAN: { 
          
            // Arctan
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          
          _Operation *newOp  = new _Operation( HY_OP_CODE_DIV, 2),
          *newOp2 = new _Operation( HY_OP_CODE_ADD, 2),
          *newOp3 = new _Operation(new _Constant(1.0)),
          *newOp4 = new _Operation( HY_OP_CODE_POWER, 2),
          *newOp5 = new _Operation(new _Constant(2.0));
          
          node<long> *newNode2 = new node<long>;
          node<long> *newNode3 = new node<long>;
          node<long> *newNode4 = new node<long>;
          node<long> *newNode5 = new node<long>;
          
          newNode4->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          newNode4->add_node(*newNode5);
          
          newNode2->add_node(*newNode3);
          newNode2->add_node(*newNode4);
          
          newNode->add_node(*b1);
          newNode->add_node(*newNode2);
          
          newNode5->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp5);
          newNode4->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp4);
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp3);
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_COS: { 
            // Cos
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          
          _Operation *newOp  = new _Operation( HY_OP_CODE_MUL, 2),
          *newOp2 = new _Operation( HY_OP_CODE_SUB, 1),
          *newOp3 = new _Operation( HY_OP_CODE_SIN, 1);
          
          node<long> *newNode2 = new node<long>;
          node<long> *newNode3 = new node<long>;
          
          newNode3->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode2->add_node(*newNode3);
          
          newNode->add_node(*newNode2);
          newNode->add_node(*b1);
          
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula << newOp3;
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula << newOp2;
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula << newOp;
          
          DeleteObject(newOp);
          DeleteObject(newOp2);
          DeleteObject(newOp3);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_ERF: { 
            // Erf
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          
          _Operation *newOp  = new _Operation( HY_OP_CODE_MUL, 2),
          *newOp2 = new _Operation( HY_OP_CODE_DIV, 2),
          *newOp3 = new _Operation(new _Constant(twoOverSqrtPi)),
          *newOp4 = new _Operation( HY_OP_CODE_EXP, 1),
          *newOp5 = new _Operation( HY_OP_CODE_SUB, 1),
          *newOp6 = new _Operation( HY_OP_CODE_POWER, 2),
          *newOp7 = new _Operation(new _Constant(2.0));
          
          
          node<long> *newNode2 = new node<long>;
          node<long> *newNode3 = new node<long>;
          node<long> *newNode4 = new node<long>;
          node<long> *newNode5 = new node<long>;
          node<long> *newNode6 = new node<long>;
          node<long> *newNode7 = new node<long>;
          
          newNode6->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          newNode6->add_node(*newNode7);
          
          newNode5->add_node(*newNode6);
          newNode4->add_node(*newNode5);
          newNode2->add_node(*newNode4);
          newNode2->add_node(*newNode3);
          
          newNode->add_node(*b1);
          newNode->add_node(*newNode2);
          
          newNode7->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp7);
          newNode6->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp6);
          newNode5->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp5);
          newNode4->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp4);
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp3);
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_EXP:   // HY_OP_CODE_EXP
        case HY_OP_CODE_SIN: { // HY_OP_CODE_SIN
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          long opC2;
          
          if (op->GetAVariable() == HY_OP_CODE_SIN) {
            opC2 = HY_OP_CODE_COS;
          } else {
            opC2 = HY_OP_CODE_EXP;
          }
          
          _Operation *newOp = new _Operation( HY_OP_CODE_MUL, 2L),
          *newOp2 = new _Operation( opC2, 1L);
          
          node<long> *newNode2 = new node<long>;
          
          newNode2->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode->add_node(*newNode2);
          newNode->add_node(*b1);
          
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_LOG: { 
            // Log
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          _Operation *newOp = new _Operation( HY_OP_CODE_DIV, 2L);
          
          
          newNode->add_node(*b1);
          newNode->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_SQRT: { 
            // Sqrt
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
  
          _Operation *newOp = new _Operation( HY_OP_CODE_DIV, 2L),
          *newOp2 = new _Operation( HY_OP_CODE_MUL, 2L),
          *newOp3 = new _Operation( HY_OP_CODE_SQRT, 1L),
          *newOp4 = new _Operation(new _Constant(2.0));
          
          
          node<long> *newNode2 = new node<long>;
          node<long> *newNode3 = new node<long>;
          node<long> *newNode4 = new node<long>;
          
          newNode3->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode2->add_node(*newNode4);
          newNode2->add_node(*newNode3);
          
          newNode->add_node(*b1);
          newNode->add_node(*newNode2);
          
          newNode4->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp4);
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp3);
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_TAN: { 
            // Tan
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          
          _Operation *newOp = new _Operation( HY_OP_CODE_DIV, 2L),
          *newOp2 = new _Operation( HY_OP_CODE_POWER , 2L),
          *newOp3 = new _Operation(new _Constant(2.0)),
          *newOp4 = new _Operation( HY_OP_CODE_COS, 1L);
          
          node<long> *newNode2 = new node<long>;
          node<long> *newNode3 = new node<long>;
          node<long> *newNode4 = new node<long>;
          
          
          newNode4->add_node(
                             *DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNode2->add_node(*newNode4);
          newNode2->add_node(*newNode3);
          
          newNode->add_node(*b1);
          newNode->add_node(*newNode2);
          
          newNode4->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp4);
          newNode3->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp3);
          newNode2->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp2);
          newNode->in_object = tgt.theFormula.lLength;
          tgt.theFormula.AppendNewInstance(newOp);
          
          return newNode;
        } break;
          
        case HY_OP_CODE_POWER: {
            // ^
            // f[x]^g[x] (g'[x] Log[f[x]] + f'[x]g[x]/f[x])
          node<long> *b1 = InternalDifferentiate(currentSubExpression->go_down(1),
                                                 varID, varRefs, dydx, tgt);
          if (!b1) {
            newNode->delete_tree(true);
            return nil;
          }
          node<long> *b2 = InternalDifferentiate(currentSubExpression->go_down(2),
                                                 varID, varRefs, dydx, tgt);
          if (!b2) {
            newNode->delete_tree(true);
            return nil;
          }
          
          long opCodes[7] = { HY_OP_CODE_MUL, HY_OP_CODE_POWER, HY_OP_CODE_ADD,
            HY_OP_CODE_DIV, HY_OP_CODE_MUL, HY_OP_CODE_MUL,
            HY_OP_CODE_LOG };
          long opArgs[7] = { 2, 2, 2, 2, 2, 2, 1 };
          
          node<long> *newNodes[7];
          newNodes[0] = newNode;
          for (long k = 1; k < 7; k++) {
            newNodes[k] = new node<long>;
          }
          
          newNodes[6]
          ->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNodes[5]->add_node(*b2);
          newNodes[5]->add_node(*newNodes[6]);
          
          newNodes[4]->add_node(*b1);
          newNodes[4]
          ->add_node(*DuplicateFormula(currentSubExpression->go_down(2), tgt));
          
          newNodes[3]->add_node(*newNodes[4]);
          newNodes[3]
          ->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
          
          newNodes[2]->add_node(*newNodes[5]);
          newNodes[2]->add_node(*newNodes[3]);
          
          newNodes[1]
          ->add_node(*DuplicateFormula(currentSubExpression->go_down(1), tgt));
          newNodes[1]
          ->add_node(*DuplicateFormula(currentSubExpression->go_down(2), tgt));
          
          newNode->add_node(*newNodes[1]);
          newNode->add_node(*newNodes[2]);
          
          for (long k = 6; k >= 0; k--) {
            newNodes[k]->in_object = tgt.theFormula.lLength;
            tgt.theFormula.AppendNewInstance(new _Operation( opCodes[k], opArgs[k]));
          }
          return newNode;
        }
      }
      
      delete (newNode);
    }
  }
  return nil;
}
