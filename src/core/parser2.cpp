/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

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

#include "legacy_parser.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "ctype.h"
#include "calcnode.h"
#include "time.h"
#include "likefunc.h"
#include "polynoml.h"
#include "float.h"
#include "batchlan.h"
#include "category.h"
#include "hy_globals.h"
#include "parser2013.h"
#include "Parser.h"

extern _SimpleList BinOps, opPrecedence, FunctionArgumentCount, associativeOps;

extern long expressionsParsed;

extern _AVLListX _HY_GetStringGlobalTypes;

extern _Parameter tolerance, sqrtPi, maxRombergSteps,
    integrationPrecisionFactor, machineEps;

extern _String intPrecFact, intMaxIter;

_Parameter verbosityLevel = 0.0;

extern _Parameter twoOverSqrtPi;

_SimpleList simpleOperationCodes, simpleOperationFunctions;

long subNumericValues = 0;

//______________________________________________________________________________
_Parameter AddNumbers(_Parameter x, _Parameter y) { return x + y; }
_Parameter SubNumbers(_Parameter x, _Parameter y) { return x - y; }
_Parameter MinusNumber(_Parameter x) { return -x; }
_Parameter MultNumbers(_Parameter x, _Parameter y) { return x * y; }
_Parameter DivNumbers(_Parameter x, _Parameter y) { return x / y; }
_Parameter LessThan(_Parameter x, _Parameter y) { return x < y; }
_Parameter GreaterThan(_Parameter x, _Parameter y) { return x > y; }
_Parameter LessThanE(_Parameter x, _Parameter y) { return x <= y; }
_Parameter GreaterThanE(_Parameter x, _Parameter y) { return x >= y; }
_Parameter Power(_Parameter x, _Parameter y) { 
    if (x==0.0) {
      if (y > 0.0) {
        return 0.0;
      } else {
        return 1.0;
      }
    }
    return pow(x, y); 
}
_Parameter MaxNumbers(_Parameter x, _Parameter y) { return x < y ? y : x; }
_Parameter MinNumbers(_Parameter x, _Parameter y) { return x < y ? x : y; }
_Parameter ExpNumbers(_Parameter x) { return exp(x); }
_Parameter LogNumbers(_Parameter x) { return log(x); }
_Parameter FastMxAccess(Ptr m, _Parameter index) {
  return ((_Parameter *)m)[(long) index];
}
_Parameter AndNumbers(_Parameter x, _Parameter y) {
  return x != 0.0 && y != 0.0;
}
_Parameter OrNumbers(_Parameter x, _Parameter y) {
  return x != 0.0 || y != 0.0;
}
_Parameter AbsNumber(_Parameter x) { return fabs(x); }

//______________________________________________________________________________
_Parameter RandomNumber(_Parameter l, _Parameter u) {
  _Parameter r = l;
  if (u > l) {
    r = genrand_int32();
    r /= RAND_MAX_32;
    r = l + (u - l) * r;
  }
  return r;
}

//______________________________________________________________________________
_Parameter EqualNumbers(_Parameter a, _Parameter b) {
  if (a != 0.0) {
    a = (a > b) ? (a - b) / a : (b - a) / a;
    return ((a > 0.) ? (a <= machineEps) : (a >= -machineEps));
  }
  return (b <= machineEps) && (b >= -machineEps);
}

//______________________________________________________________________________
void PopulateArraysForASimpleFormula(_SimpleList &vars,
                                     _SimpleFormulaDatum *values) {
  for (unsigned long k2 = 0; k2 < vars.lLength; k2++) {
    _PMathObj varValue = LocateVar(vars.lData[k2])->Compute();
    if (varValue->ObjectClass() == NUMBER) {
      values[k2].value = varValue->Value();
    } else {
      values[k2].reference = (Ptr)(_HY2MATRIX(varValue))->theData;
    }
  }
}

//______________________________________________________________________________
void WarnNotDefined(_PMathObj p, long opCode, _hyExecutionContext *context) {
  _FString *t = (_FString *)p->Type();
  context->ReportError(_String("Operation '") &
                       BuiltInFunctions.RetrieveKeyByPayload (opCode) &
                       "' is not implemented/defined for a " & *t->theString);
  DeleteObject(t);
}

//______________________________________________________________________________
_Parameter InterpolateValue(_Parameter *theX, _Parameter *theY, long n,
                            _Parameter *c, _Parameter *d, _Parameter x,
                            _Parameter &err) {
  _Parameter y, den, dif = 1e10, dift, ho, hp, w;

  long ns;

  for (long i = 0; i < n; i++) {
    dift = fabs(x - theX[i]);
    if (dift < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = d[i] = theY[i];
  }

  y = theY[ns];
  ns--;

  for (long m = 1; m < n; m++) {
    for (long i = 0; i <= n - m - 1; i++) {
      ho = theX[i] - x;
      hp = theX[i + m] - x;
      w = c[i + 1] - d[i];
      den = w / (ho - hp);
      d[i] = hp * den;
      c[i] = ho * den;
    }
    err = 2 * ns < (n - m) ? c[ns + 1] : d[ns--];
    y += err;
  }

  return y;
}

//______________________________________________________________________________
_Parameter TrapezoidLevelKSimple(_Formula &f, _Variable *xvar, _Parameter left,
                                 _Parameter right, long k,
                                 _SimpleFormulaDatum *stack,
                                 _SimpleFormulaDatum *values,
                                 _SimpleList &changingVars,
                                 _SimpleList &varToStack) {
  _Parameter x, tnm, sum, del, ddel;

  static _Parameter s;

  //_Constant dummy;

  long it, j;

  if (k == 1) {
    if (changingVars.lLength == 1) {
      values[varToStack.lData[0]].value = (left + right) * 0.5;
    } else {
      xvar->SetValue(new _Constant((left + right) * 0.5), false);
      for (long vi = 0; vi < changingVars.lLength; vi++) {
        values[varToStack.lData[vi]].value =
            LocateVar(changingVars.lData[vi])->Compute()->Value();
      }
    }
    s = f.ComputeSimple(stack, values);
    return s;
  }

  for (it = 1, j = 1; j < k - 1; j++) {
    it *= 3;
  }

  tnm = it;
  del = (right - left) / (3.0 * tnm);
  ddel = del + del;
  x = left + del * .5;
  for (sum = 0.0, j = 1; j <= it; j++, x += del) {
    if (changingVars.lLength == 1) {
      values[varToStack.lData[0]].value = x;
    } else {
      xvar->SetValue(new _Constant(x), false);
      for (long vi = 0; vi < changingVars.lLength; vi++) {
        values[varToStack.lData[vi]].value =
            LocateVar(changingVars.lData[vi])->Compute()->Value();
      }
    }
    sum += f.ComputeSimple(stack, values);

    x += ddel;

    if (changingVars.lLength == 1) {
      values[varToStack.lData[0]].value = x;
    } else {
      xvar->SetValue(new _Constant(x), false);
      for (long vi = 0; vi < changingVars.lLength; vi++) {
        values[varToStack.lData[vi]].value =
            LocateVar(changingVars.lData[vi])->Compute()->Value();
      }
    }
    sum += f.ComputeSimple(stack, values);
  }
  s = (s + (right - left) * sum / tnm) / 3.0;
  return s;
}

//______________________________________________________________________________
_Parameter TrapezoidLevelK(_Formula &f, _Variable *xvar, _Parameter left,
                           _Parameter right, long k) {
  _Parameter x, tnm, sum, del, ddel;

  static _Parameter s;

  _Constant dummy;

  long it, j;

  if (k == 1) {
    dummy.SetValue((left + right) / 2);
    xvar->SetValue(&dummy);
    s = f.Compute()->Value();
    return s;
  }

  for (it = 1, j = 1; j < k - 1; j++) {
    it *= 3;
  }

  tnm = it;
  del = (right - left) / (3.0 * tnm);
  ddel = del + del;
  x = left + del * .5;
  for (sum = 0.0, j = 1; j <= it; j++, x += del) {
    dummy.SetValue(x);
    xvar->SetValue(&dummy);
    sum += f.Compute()->Value();
    x += ddel;
    dummy.SetValue(x);
    xvar->SetValue(&dummy);
    sum += f.Compute()->Value();
  }
  s = (s + (right - left) * sum / tnm) / 3.0;
  return s;
}

//______________________________________________________________________________
long ExecuteFormula(_Formula *f, _Formula *f2, long code, long reference,
                    _VariableContainer *nameSpace, char assignment_type) {
  if (assignment_type != HY_STRING_DIRECT_REFERENCE && reference >= 0) {
    long dereferenced =
        DereferenceVariable(reference, nameSpace, assignment_type);
    if (dereferenced < 0) {
      WarnError(_String("Failed to dereference '") &
                *FetchVar(reference)->GetName() & "' in the " &
                ((assignment_type == HY_STRING_GLOBAL_DEREFERENCE) ? "global"
                                                                   : "local") &
                " context");
      return 0;
    }
    reference = dereferenced;
  }

  _hyExecutionContext localContext (nameSpace);

  if (code == HY_FORMULA_EXPRESSION ||
      code == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT ||
      code == HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT) {
    _PMathObj formulaValue = (code == HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT)
                                 ? f2->Compute(0, &localContext)
                                 : f->Compute(0, &localContext);
    if (!formulaValue) {
      return 0;
    }

    if (code == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
      LocateVar(reference)->SetValue(formulaValue);
      return 1;
    }

    if (code == HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT) {
      _hyExecutionContext localContext(nameSpace);
      _Variable *theV = f->Dereference(
          assignment_type == HY_STRING_GLOBAL_DEREFERENCE, &localContext);
      if (theV) {
        theV->SetValue(formulaValue);
      } else {
        return 0;
      }
      return 1;
    }
    return 0;
  }

  if (code == HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT) {
    _Formula fFixed;
    fFixed.DuplicateReference(f);
    LocateVar(reference)->SetFormula(fFixed);
    return 1;
  }

  if (code == HY_FORMULA_REFERENCE_FORMULA_ASSIGNMENT) {
    _hyExecutionContext localContext(nameSpace);
    _Variable *theV = f->Dereference(
        assignment_type == HY_STRING_GLOBAL_DEREFERENCE, &localContext);
    if (theV) {
      _Formula fFixed;
      fFixed.DuplicateReference(f2);
      theV->SetFormula(fFixed);
      return 1;
    } else {
      return 0;
    }

  }

  if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT ||
      code == HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT ||
      code == HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT ||
      code == HY_FORMULA_REFERENCE_LOWER_BOUND_ASSIGNMENT) {
    if (f2->IsEmpty()) {
      WarnError("Empty RHS in a constraint assignment.");
      return 0;
    }

    _PMathObj varObj = f2->Compute(0, &localContext);
    if (varObj->ObjectClass() != NUMBER) {
      WarnError("Not a numeric RHS in a constraint assignment.");
      return 0;
    }

    _Variable *theV;

    if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT ||
        code == HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT) {
      theV = LocateVar(reference);
    } else {
      _hyExecutionContext localContext(nameSpace);
      theV = f->Dereference(assignment_type == HY_STRING_GLOBAL_DEREFERENCE,
                            &localContext);
      if (!theV) {
        return 0;
      }
    }

    if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT ||
        code == HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT) {
      theV->SetBounds(varObj->Value(), theV->GetUpperBound());
    } else {
      theV->SetBounds(theV->GetLowerBound(), varObj->Value());
    }

    /*
        SLKP 20110301 if the new constraint makes the current variable value
        invalid, then the value will be modified to stay in bounds*/

    theV->EnsureTheValueIsInBounds();

  }

  if (code == HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT ||
      code == HY_FORMULA_FORMULA_VALUE_ASSIGNMENT ||
      code == HY_FORMULA_FORMULA_VALUE_INCREMENT) {
    _Formula newF;

    if (f2->IsEmpty()) {
      WarnError("Empty RHS in an assignment.");
      return 0;
    }

    if (code == HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT) {
      newF.DuplicateReference(f2);
    } else {
      newF.theFormula.AppendNewInstance(
          new _Operation(_HY2MATHOBJ (f2->Compute(0, &localContext)->makeDynamic())));
    }

    long stackD = -1, last0 = 0;


    for (long opID = 0; opID < f->theFormula.lLength - 1; opID++) {
      f->GetIthTerm (opID)->StackDepth(stackD);
      if (stackD == 0) {
        last0 = opID;
      }
    }

    _Matrix *mmx = nil;
    _AssociativeList *mma = nil;

    if (last0 > 0) {
      stackD = f->theFormula.lLength;
      f->theFormula.lLength = last0 + 1;
      _PMathObj lvalue = f->Compute(0, &localContext);
      f->theFormula.lLength = stackD;
      if (lvalue->ObjectClass() == MATRIX) {
        mmx = _HY2MATRIX (lvalue);
      }
      if (lvalue->ObjectClass() == ASSOCIATIVE_LIST) {
        mma = _HY2ASSOCIATIVE_LIST (lvalue);
      }
      last0++;
    } else {
      _Variable *mmo =
          LocateVar(f->GetIthTerm (0)->GetAVariable());

      if (mmo) {
        if (mmo->ObjectClass() == MATRIX) {
          mmx = _HY2MATRIX (mmo->GetValue());
          f->GetIthTerm (0)->ToggleVarRef(true);
        } else {
          if (mmo->ObjectClass() == ASSOCIATIVE_LIST) {
            mma = _HY2ASSOCIATIVE_LIST (mmo->GetValue());
            f->GetIthTerm (0)->ToggleVarRef(true);
          }
        }
      }
    }

    _PMathObj coordMx = nil;
    if (mma || mmx) {
      long expectedType = mmx ? MATRIX : STRING;
      coordMx = f->Compute(last0, &localContext);
      if (!coordMx || coordMx->ObjectClass() != expectedType) {
        if (mmx) {
          WarnError(_String("Matrix expected but not supplied."));
        } else {
          WarnError(_String("String key expected but not supplied."));
        }

        return 0;
      }
    } else {
      WarnError("Matrix/List LHS expected but not supplied.");
      return 0;
    }

    if (mmx) { // matrix LHS
      _Matrix *mcoord = _HY2MATRIX(coordMx);

      long hC = mcoord->theData[0], vC = mcoord->theData[1];

      if (mmx->CheckCoordinates(hC, vC)) {
        if (!ANALYTIC_COMPUTATION_FLAG) {
          mmx->MStore(hC, vC, newF, (code == HY_FORMULA_FORMULA_VALUE_INCREMENT)
                                        ? HY_OP_CODE_ADD
                                        : HY_OP_CODE_NONE);
        } else {
          _PMathObj newP = newF.ConstructPolynomial();
          if (!newP) {
            warnError(_String(
                "Can't assign non-polynomial entries to polynomial matrices."));
          } else {
            mmx->MStore(hC, vC, newP);
          }
        }
        mmx->CheckIfSparseEnough();
      }
    } else if (mma) { // Associative array LHS
      mma->MStore(
          coordMx, newF.Compute(), true,
          (code == HY_FORMULA_FORMULA_VALUE_INCREMENT) ? HY_OP_CODE_ADD
                                                       : HY_OP_CODE_NONE);
    }

    return 1;
  }
  return 0;
}

//______________________________________________________________________________

struct      characterChecker {
    characterChecker (_String s) {
        for (long r = 0; r<256; r++) {
            isAllowed [r] = false;
        }
        for (long r2 = 0; r2<s.sLength; r2++) {
            isAllowed [(unsigned char)s.sData[r2]] = true;
        }
    }
    bool     isAllowed [256];
}
alpha       ("ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz"),
            numeric     (".0123456789eE");

_String     globalToken ("global"),
            noneToken   ("None");

/*
_String     alpha   ("ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz"),
            numeric (".0123456789eE");
*/

#define CONTEXT_TRUNCATION 24

//______________________________________________________________________________
long HandleFormulaParsingError(_String errMsg, _String *saveError, _String &s,
                               long index) {
  if (index >= 0) {
    errMsg =
        errMsg & " in the following context: '" &
        s.Cut(MAX(0, index - CONTEXT_TRUNCATION), index) & "<ERROR HERE>" &
        s.Cut(index + 1, MIN(index + CONTEXT_TRUNCATION, s.sLength - 1)) & "'";
  }
  if (saveError) {
    *saveError = errMsg;
  } else {
    WarnError(errMsg);
  }
  return HY_FORMULA_FAILED;
}


//______________________________________________________________________________
long _parserHelperHandleInlineBoundCases(_String &s,
                                         _FormulaParsingContext &parsingContext,
                                         long i, _Variable *lhs_variable,
                                         _Formula *f, char deref,
                                         _Formula &newF) {
  _PMathObj varObj = newF.Compute();
  if (varObj->ObjectClass() != NUMBER) {
    return HandleFormulaParsingError(
        "Variable bound must evaluate to a number ", parsingContext.errMsg(), s,
        i);
  }

  long varID;
  
  _hyExecutionContext localContext (parsingContext.formulaScope(), parsingContext.errMsg());

  if (lhs_variable) {
    varID = DereferenceVariable(lhs_variable->GetAVariable(),
                                parsingContext.formulaScope(), deref);
  } else {
    varID = DereferenceString(f->Compute(0, &localContext),
                              parsingContext.formulaScope(), deref);
  }
  if (varID < 0) {
    return HandleFormulaParsingError("Failed to dereference ",
                                     parsingContext.errMsg(), s, i);
  }

  _Variable *theV = (_Variable *)LocateVar(varID);

  if (s.getChar(i) == '>') {
    theV->SetBounds(varObj->Value(), theV->GetUpperBound());
  } else {
    theV->SetBounds(theV->GetLowerBound(), varObj->Value());
  }
  return HY_FORMULA_EXPRESSION;
}

//______________________________________________________________________________
long _parserHelperHandleInlineAssignmentCases(
    _String &s, _FormulaParsingContext &parsingContext, long i,
    _Variable *lhs_variable, _Formula *f, char deref, _Formula &newF,
    bool twoToken) {

  long varID;

  _hyExecutionContext localContext (parsingContext.formulaScope(), parsingContext.errMsg());

  if (lhs_variable) {
    varID = DereferenceVariable(lhs_variable->GetAVariable(),
                                parsingContext.formulaScope(), deref);
  } else {
    varID = DereferenceString(f->Compute(0, &localContext),
                              parsingContext.formulaScope(), deref);
  }
  if (varID < 0) {
    return HandleFormulaParsingError("Failed to dereference ",
                                     parsingContext.errMsg(), s, i);
  }
  _Variable *theV = (_Variable *)LocateVar(varID);

  if (s.getChar(i - 1) != ':') {
    _PMathObj varObj = newF.Compute();
    if (!varObj) {
      return HandleFormulaParsingError("Invalid RHS in an assignment ",
                                       parsingContext.errMsg(), s, i);
    }
    if (twoToken && s.getChar(i - 1) == '+') {
      theV->SetValue(theV->Compute()->Execute(HY_OP_CODE_ADD, varObj));
    } else {
      theV->SetValue(varObj);
    }
  } else {
    theV->SetFormula(newF);
  }
  return HY_FORMULA_EXPRESSION;
}

//______________________________________________________________________________
// SLKP 20110908: added the concept of a 'volatile' formula, i.e. something that should be reparsed every time in ExecuteCase0
//                : currently those include
//                :    inline constructors (matrices, dictionaries)
//                :    `` substitutions in strings
// SLKP 20100817: decoupled return code from variable reference return
// returns:
// case                       | return value                                  | parsingContext.assignment_ref_id value

// parse failed               | HY_FORMULA_FAILED                             | undefined
// expresion (no LHS)         | HY_FORMULA_EXPRESSION                         | undefined
// z = x/y                    | HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT          | index of the LHS
// *|^(expr) = expr           | HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT         | undefined
// z := x/y                   | HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT        | index of the LHS
// object[a] := x/y           | HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT         | undefined
// object[a] = x/y            | HY_FORMULA_FORMULA_VALUE_ASSIGNMENT           | undefined
// z :< expr                  | HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT    | index of the LHS
// z :> expr                  | HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT    | index of the LHS
// Further, for (HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT,  HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT, 
// HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT, HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT, 
// HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT):
// case            |   parsingContext.assignment_ref_type value
// z op x/y        |   HY_STRING_DIRECT_REFERENCE
// *z op x/y       |   HY_STRING_LOCAL_DEREFERENCE
// ^z op x/y       |   HY_STRING_GLOBAL_DEREFERENCE
long Parse(_Formula *f, _String &s, _FormulaParsingContext &parsingContext,
           _Formula *f2) {

  expressionsParsed++;

  Scanner *scanner = new Scanner((unsigned char*)s.getStr(), s.sLength);
  
  Parser  *parser = new Parser(scanner, f, &parsingContext);
  parser->Parse();
  
  bool has_errors = parser->errors->count > 0;
  
  delete parser;
  delete scanner;
  

  return has_errors ? HY_FORMULA_FAILED : f->FormulaType();
}

//______________________________________________________________________________
long VerbosityLevel(void) {
  checkParameter(VerbosityLevelString, verbosityLevel, -1.0);
  return verbosityLevel;
}

//______________________________________________________________________________
void checkParameter(_String &name, _Parameter &dest, _Parameter def,
                    _VariableContainer *theP) {
  long f;
  if (theP) {
    _String ppn = *theP->GetName() & '.' & name;
    f = LocateVarByName(ppn);
  } else {
    f = LocateVarByName(name);
  }
  if (f < 0) {
    dest = def;
  } else {
    dest = FetchVar(f)->Value();
  }
}

//______________________________________________________________________________
void stashParameter(_String &name, _Parameter v, bool set) {

  static _Parameter stash = 0.0;

  long f = LocateVarByName(name);
  if (f >= 0) {
    _Variable *thisV = FetchVar(f);
    if (set) {
      stash = thisV->Value();
      _Constant dummy(v);
      thisV->SetValue(&dummy);
    } else {
      _Constant dummy(stash);
      thisV->SetValue(&dummy);
    }
  } else if (set) {
    stash = v;
    setParameter(name, v);
  }
}

//______________________________________________________________________________
void setParameter(_String &name, _Parameter def, _String *namespc) {
  if (namespc) {
    _String namespcd = AppendContainerName(name, namespc);
    setParameter(namespcd, def);
  } else {
    long f = LocateVarByName(name);
    if (f < 0) {
      _Variable cornholio(name);
      setParameter(name, def);
    } else {
      FetchVar(f)->SetValue(new _Constant(def), false);
    }
  }
}

//______________________________________________________________________________
void setParameter(_String &name, _PMathObj def, bool dup, _String *namespc) {
  if (namespc) {
    _String namespcd = AppendContainerName(name, namespc);
    setParameter(namespcd, def, dup);
  } else {
    long f = LocateVarByName(name);
    if (f < 0) {
      _Variable cornholio(name);
      setParameter(name, def, dup);
    } else {
      FetchVar(f)->SetValue(def, dup);
    }
  }
}

//______________________________________________________________________________
void ExportIndVariables(_String &glVars, _String &locVars,
                        _SimpleList *indepVarList) {
  _String *stIn, str;

  for (unsigned long i = 0; i < indepVarList->lLength; i++) {
    _Variable *thisVar = LocateVar(indepVarList->lData[i]);
    if (thisVar->IsGlobal()) {
      str = _String("\nglobal ") & *thisVar->GetName() & '=' &
            _String((_String *)parameterToString(thisVar->Compute()->Value())) &
            ';';
      stIn = &glVars;
    } else {
      str = _String("\n") & *thisVar->GetName() & '=' &
            _String((_String *)parameterToString(thisVar->Compute()->Value())) &
            ';';
      stIn = &locVars;
    }
    *stIn << str;
    if (!CheckEqual(thisVar->GetLowerBound(), DEFAULTPARAMETERLBOUND)) {
      str =
          _String("\n") & *thisVar->GetName() & ":>" &
          _String((_String *)parameterToString(thisVar->GetLowerBound())) & ';';
      *stIn << str;
    }
    if (!CheckEqual(thisVar->GetUpperBound(), DEFAULTPARAMETERUBOUND)) {
      str =
          _String("\n") & *thisVar->GetName() & ":<" &
          _String((_String *)parameterToString(thisVar->GetUpperBound())) & ';';
      *stIn << str;
    }
  }
}

//______________________________________________________________________________
void ExportDepVariables(_String &glVars, _String &locVars,
                        _SimpleList *depVarList) {
  if (depVarList->lLength) {
    _String *stIn, str;

    //first we have to reorder global variables, so that dependent global
    //  variables which depend
    //  on other dependent global variables are written afterwards (lest they be
    //  implicitly declared
    //  as local).
    //  The algorithm is very ugly, but since there are only a few global
    //  dependent variables (in general...) 

    _SimpleList _globalVariablesList, lfDepGlobs, tl1;

    _List dependancyLists;
    {
      for (unsigned long i = 0; i < depVarList->lLength; i++)
        if (LocateVar(depVarList->lData[i])->IsGlobal()) {
          lfDepGlobs << depVarList->lData[i];
        }
    }

    lfDepGlobs.Sort();

    for (unsigned long i = 0; i < depVarList->lLength; i++) {
      _Variable *thisVar = LocateVar(depVarList->lData[i]);
      if (thisVar->IsGlobal()) {
        _SimpleList globDependancyList, prunedList;

        _AVLList globDependancyListAVL(&globDependancyList);

        thisVar->ScanForVariables(globDependancyListAVL, true);

        globDependancyListAVL.ReorderList();

        prunedList.Intersect(globDependancyList, lfDepGlobs);

        if (prunedList.lLength) {
          _globalVariablesList << i;
          dependancyLists &&&prunedList;
          continue;
        }
        str = _String("\nglobal ") & *thisVar->GetName();
        stIn = &glVars;
      } else {
        str = _String("\n") & *thisVar->GetName();
        stIn = &locVars;
      }
      (*stIn) << str;
      (*stIn) << ":=";
      stIn->AppendNewInstance(thisVar->GetFormulaString());
      (*stIn) << ';';
      if (!CheckEqual(thisVar->GetLowerBound(), DEFAULTPARAMETERLBOUND)) {
        str = _String("\n") & *thisVar->GetName() & ":>" &
              _String((_String *)parameterToString(thisVar->GetLowerBound())) &
              ';';
        (*stIn) << str;
      }
      if (!CheckEqual(thisVar->GetUpperBound(), DEFAULTPARAMETERUBOUND)) {
        str = _String("\n") & *thisVar->GetName() & ":<" &
              _String((_String *)parameterToString(thisVar->GetUpperBound())) &
              ';';
        (*stIn) << str;
      }
    }

    if (_globalVariablesList.lLength)
        // check internal dependancies
        {
      _SimpleList writeOrder(_globalVariablesList.lLength, 0, 1),
          indexList(_globalVariablesList.lLength, 0, 1);

      for (unsigned long i2 = 0; i2 < _globalVariablesList.lLength; i2++) {
        long updatedIndex = writeOrder.lData[i2];
        _SimpleList *depList = _HY2SIMPLELIST (dependancyLists(i2));
        for (unsigned long i3 = 0; i3 < depList->lLength; i3++) {
          long i4 = _globalVariablesList.Find(depList->lData[i3]);
          if (i4 >= 0 && updatedIndex < writeOrder.lData[i4]) {
            updatedIndex = writeOrder.lData[i4] + 1;
          }
        }
        writeOrder.lData[i2] = updatedIndex;
      }

      SortLists(&writeOrder, &indexList);

      for (unsigned long i = 0; i < _globalVariablesList.lLength; i++) {
        _Variable *thisVar = LocateVar(
            depVarList->lData[_globalVariablesList.lData[indexList.lData[i]]]);
        str = _String("\nglobal ") & *thisVar->GetName();
        glVars << str;
        glVars << ":=";
        glVars << thisVar->GetFormulaString();
        glVars << ';';
        if (!CheckEqual(thisVar->GetLowerBound(), DEFAULTPARAMETERLBOUND)) {
          str =
              _String("\n") & *thisVar->GetName() & ":>" &
              _String((_String *)parameterToString(thisVar->GetLowerBound())) &
              ';';
          glVars << str;
        }
        if (!CheckEqual(thisVar->GetUpperBound(), DEFAULTPARAMETERUBOUND)) {
          str =
              _String("\n") & *thisVar->GetName() & ":<" &
              _String((_String *)parameterToString(thisVar->GetUpperBound())) &
              ';';
          glVars << str;
        }
      }
    }
  }
}

//______________________________________________________________________________
void ExportCatVariables(_String &rec, _SimpleList *catVarList) {
  _SimpleList nonInd;

  for (long idx = 0; idx < catVarList->lLength; idx++)
    if (((_CategoryVariable *)LocateVar(catVarList->lData[idx]))
            ->IsUncorrelated()) {
      ((_CategoryVariable *)LocateVar(catVarList->lData[idx]))
          ->SerializeCategory(rec);
    } else {
      nonInd << idx;
    }
  {
    for (long idx = 0; idx < nonInd.lLength; idx++) {
      ((_CategoryVariable *)LocateVar(catVarList->lData[nonInd.lData[idx]]))
          ->SerializeCategory(rec);
    }
  }
}

//______________________________________________________________________________
void SplitVariablesIntoClasses(_SimpleList &all, _SimpleList &i, _SimpleList &d,
                               _SimpleList &c) {
  for (long idx = 0; idx < all.lLength; idx++) {
    _Variable *thisVar = LocateVar(all.lData[idx]);
    if (thisVar->IsCategory()) {
      c << all.lData[idx];
    } else if (thisVar->IsIndependent()) {
      i << all.lData[idx];
    } else {
      d << all.lData[idx];
    }
  }
}
