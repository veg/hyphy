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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <time.h>


#include "parser.h"
#include "matrix.h"
#include "calcnode.h"
#include "likefunc.h"
#include "polynoml.h"
#include "batchlan.h"
#include "category.h"
#include "function_templates.h"
#include "global_things.h"

using namespace hy_global;


#define  GOLDEN_RATIO 1.618034
#define MAX_BRENT_ITERATES 100L

extern   _SimpleList BinOps,
         opPrecedence,
         FunctionArgumentCount,
         associativeOps;

extern   long           expressionsParsed;

extern    _List         FunctionNameList;

extern    _AVLListX     _HY_GetStringGlobalTypes;

extern    hyFloat    tolerance,
          maxRombergSteps,
          integrationPrecisionFactor,
          machineEps;


long      verbosityLevel = 0L;


_SimpleList simpleOperationCodes,
            simpleOperationFunctions;


long        subNumericValues = 0;

//__________________________________________________________________________________
hyFloat  AddNumbers  (hyFloat x, hyFloat y)
{
    return x+y;
}
hyFloat  SubNumbers  (hyFloat x, hyFloat y)
{
    return x-y;
}
hyFloat  MinusNumber (hyFloat x)
{
    return -x;
}
hyFloat  MultNumbers (hyFloat x, hyFloat y)
{
    return x*y;
}
hyFloat  DivNumbers  (hyFloat x, hyFloat y)
{
    return x/y;
}
hyFloat  LessThan    (hyFloat x, hyFloat y)
{
    return x<y;
}
hyFloat  GreaterThan (hyFloat x, hyFloat y)
{
    return x>y;
}
hyFloat  LessThanE   (hyFloat x, hyFloat y)
{
    return x<=y;
}
hyFloat  GreaterThanE(hyFloat x, hyFloat y)
{
    return x>=y;
}
hyFloat  Power       (hyFloat x, hyFloat y)
{ 
    if (x==0.0) {
      if (y > 0.0) {
        return 0.0;
      } else {
        return 1.0;
      }
    }
    return pow(x,y);
}
hyFloat  MaxNumbers  (hyFloat x, hyFloat y)
{
    return x<y?y:x;
}
hyFloat  MinNumbers  (hyFloat x, hyFloat y)
{
    return x<y?x:y;
}
hyFloat  ExpNumbers  (hyFloat x)
{
    return exp(x);
}
hyFloat  LogNumbers  (hyFloat x)
{
    return log(x);
}
hyFloat  FastMxAccess(hyPointer m, hyFloat index)
{
    return ((hyFloat*)m)[(unsigned long)index];
}

void  FastMxWrite(hyPointer m, hyFloat index, hyFloat value) {
  ((hyFloat*)m)[(unsigned long)index] = value;
}

hyFloat  AndNumbers  (hyFloat x, hyFloat y)
{
    return x != 0.0 && y != 0.0;
}
hyFloat  AbsNumber  (hyFloat x)
{
    return fabs (x);
}

//__________________________________________________________________________________

hyFloat  RandomNumber(hyFloat l, hyFloat u)
{
    hyFloat r = l;
    if (u>l) {
        r =l+(u-l)*genrand_real1();
    }
    return r;
}



//_______________________________________________________________________________________
hyFloat  EqualNumbers(hyFloat a, hyFloat b)
{
    if (a!=0.0) {
        a = (a>b)?(a-b)/a:(b-a)/a;
        return ((a>0.)?(a<=machineEps):(a>=-machineEps));
    }
    return (b<=machineEps)&&(b>=-machineEps);
}

//_______________________________________________________________________________________
void        PopulateArraysForASimpleFormula (_SimpleList& vars, _SimpleFormulaDatum* values)
{
    for (unsigned long k2 = 0; k2 < vars.lLength; k2++) {
        _PMathObj varValue = LocateVar (vars.lData[k2])->Compute();
        if (varValue->ObjectClass() == NUMBER) {
            values[k2].value = varValue->Value();
        } else {
            if (varValue->ObjectClass() == MATRIX) {
              values[k2].reference = (hyPointer)((_Matrix*)varValue)->theData;
            } else {
              HandleApplicationError ("Internal error in PopulateArraysForASimpleFormula", true);
            }
        }
    }
}


//__________________________________________________________________________________

void        WarnNotDefined (_PMathObj p, long opCode, _hyExecutionContext* context) {
    _FString * t = (_FString*)p->Type();
    context->ReportError  (_String("Operation '")&*(_String*)BuiltInFunctions(opCode)&"' is not implemented/defined for a " & *t->theString);
    DeleteObject (t);
}

//__________________________________________________________________________________

void        WarnWrongNumberOfArguments (_PMathObj p, long opCode, _hyExecutionContext* context, _List * args) {
  _FString * t = (_FString*)p->Type();
  context->ReportError  (_String("Operation '")&*(_String*)BuiltInFunctions(opCode)&"' was called with an incorrect number of arguments (" & (long) (args ? args->lLength : 0L) & ") for " & *t->theString);
  DeleteObject (t);
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

    for (unsigned long i=0; i<n; i++) {
        dift = fabs(x-theX[i]);
        if (dift<dif) {
            ns = i;
            dif = dift;
        }
        c[i] = d[i] = theY[i];
    }

    y = theY[ns--];

    for (unsigned long m=1; m<n; m++) {
        for (long i=0; i<=n-m-1; i++) {
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
hyFloat  TrapezoidLevelKSimple (_Formula&f, _Variable* xvar, hyFloat left, hyFloat right, long k, _SimpleFormulaDatum * stack, _SimpleFormulaDatum* values, _SimpleList& changingVars, _SimpleList& varToStack)
{
    hyFloat x,
               tnm,
               sum,
               del,
               ddel;

    static hyFloat s;

    //_Constant dummy;

    long        it,
                j;

    if (k==1) {
        if (changingVars.lLength == 1) {
            values[varToStack.lData[0]].value = (left+right)*0.5;
        } else {
            xvar->SetValue  (new _Constant ((left+right)*0.5), false);
            for (long vi = 0; vi < changingVars.lLength; vi++) {
                values[varToStack.lData[vi]].value = LocateVar(changingVars.lData[vi])->Compute()->Value();
            }
        }
        s = f.ComputeSimple(stack, values);
        return s;
    }

    for (it=1, j=1; j<k-1; j++) {
        it*=3;
    }

    tnm = it;
    del = (right-left)/(3.0*tnm);
    ddel = del+del;
    x   = left+del*.5;
    for (sum=0.0, j=1; j<=it; j++, x+=del) {
        if (changingVars.lLength == 1) {
            values[varToStack.lData[0]].value = x;
        } else {
            xvar->SetValue(new _Constant (x), false);
            for (long vi = 0; vi < changingVars.lLength; vi++) {
                values[varToStack.lData[vi]].value = LocateVar(changingVars.lData[vi])->Compute()->Value();
            }
        }
        sum += f.ComputeSimple(stack, values);

        x+=ddel;

        if (changingVars.lLength == 1) {
            values[varToStack.lData[0]].value = x;
        } else {
            xvar->SetValue(new _Constant (x), false);
            for (long vi = 0; vi < changingVars.lLength; vi++) {
                values[varToStack.lData[vi]].value = LocateVar(changingVars.lData[vi])->Compute()->Value();
            }
        }
        sum += f.ComputeSimple(stack, values);
    }
    s = (s+(right-left)*sum/tnm)/3.0;
    return s;
}

//__________________________________________________________________________________
hyFloat  TrapezoidLevelK (_Formula&f, _Variable* xvar, hyFloat left, hyFloat right, long k)
{
    hyFloat x,
               tnm,
               sum,
               del,
               ddel;

    static hyFloat s;

    _Constant dummy;

    long        it,
                j;

    if (k==1) {
        dummy.SetValue((left+right)/2);
        xvar->SetValue (&dummy);
        s = f.Compute()->Value();
        return s;
    }

    for (it=1, j=1; j<k-1; j++) {
        it*=3;
    }

    tnm = it;
    del = (right-left)/(3.0*tnm);
    ddel = del+del;
    x   = left+del*.5;
    for (sum=0.0, j=1; j<=it; j++, x+=del) {
        dummy.SetValue(x);
        xvar->SetValue(&dummy);
        sum += f.Compute()->Value();
        x+=ddel;
        dummy.SetValue(x);
        xvar->SetValue(&dummy);
        sum += f.Compute()->Value();
    }
    s = (s+(right-left)*sum/tnm)/3.0;
    return s;
}

//__________________________________________________________________________________

long       ExecuteFormula (_Formula*f , _Formula* f2, long code, long reference, _VariableContainer* nameSpace, char assignment_type)
{
    if (assignment_type != kStringDirectReference && reference >= 0) {
        long dereferenced = DereferenceVariable(reference, nameSpace, assignment_type);
        if (dereferenced < 0) {
            HandleApplicationError (_String ("Failed to dereference '") & *FetchVar(reference)->GetName() & "' in the " & ((assignment_type == kStringGlobalDeference) ? "global" : "local") & " context");
            return 0;
        }
        reference = dereferenced;
    }


    if (code == HY_FORMULA_EXPRESSION || code == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT || code == HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT) {
        _PMathObj  formulaValue = (code == HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT)?f2->Compute(0, nameSpace):f->Compute(0,nameSpace);
        if (!formulaValue) {
            return 0;
        }

        if (code == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
            LocateVar (reference)->SetValue (formulaValue);
            return 1;
        }

        if (code == HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT) {
            _hyExecutionContext localContext (nameSpace);
            _Variable * theV = f->Dereference(assignment_type == kStringGlobalDeference, &localContext);
            if (theV) {
                theV->SetValue (formulaValue);
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
        LocateVar (reference)->SetFormula (fFixed);
        return 1;
    }

    if (code == HY_FORMULA_REFERENCE_FORMULA_ASSIGNMENT) {
        _hyExecutionContext localContext (nameSpace);
        _Variable * theV = f->Dereference(assignment_type == kStringGlobalDeference, &localContext);
        if (theV) {
            _Formula fFixed;
            fFixed.DuplicateReference(f2);
            theV->SetFormula (fFixed);
            return 1;
        } else {
            return 0;
        }

    }

    if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT || code == HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT ||
        code == HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT || code == HY_FORMULA_REFERENCE_LOWER_BOUND_ASSIGNMENT ) {
        if (f2->IsEmpty()) {
            HandleApplicationError ("Empty RHS in a constraint assignment.");
            return 0;
        }

        _PMathObj varObj = f2->Compute(0, nameSpace);
        if (varObj->ObjectClass()!=NUMBER) {
            HandleApplicationError ("Not a numeric RHS in a constraint assignment.");
            return 0;
        }

        _Variable * theV;

        if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT || code == HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT) {
            theV= LocateVar (reference);
        } else {
            _hyExecutionContext localContext (nameSpace);
            theV = f->Dereference(assignment_type == kStringGlobalDeference, &localContext);
            if (!theV) {
                return 0;
            }
        }

        if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT || code == HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT) {
            theV->SetBounds(varObj->Value(),theV->GetUpperBound());
        } else {
            theV->SetBounds(theV->GetLowerBound(),varObj->Value());
        }


        /*
            SLKP 20110301 if the new constraint makes the current variable value
            invalid, then the value will be modified to stay in bounds*/

        theV->EnsureTheValueIsInBounds();

    }


    if ( code== HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT || code== HY_FORMULA_FORMULA_VALUE_ASSIGNMENT || code == HY_FORMULA_FORMULA_VALUE_INCREMENT) {
        _Formula newF;

        if (f2->IsEmpty()) {
            HandleApplicationError ("Empty RHS in an assignment.");
            return 0;
        }

        if (code == HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT) {
            newF.DuplicateReference(f2);
        } else {
            newF.theFormula.AppendNewInstance(new _Operation((_PMathObj)f2->Compute(0, nameSpace)->makeDynamic()));
        }

        long stackD = -1L,
             last0  = 0L;

        for (long opID = 0; opID < f->theFormula.lLength - 1; opID ++) {
            ((_Operation*)f->theFormula(opID)) -> StackDepth (stackD);
            if (stackD == 0L) {
                last0 = opID;
            }
        }

        _Matrix          * mmx = nil;
        _AssociativeList * mma = nil;

        if (last0 > 0) {
            stackD = f->theFormula.lLength;
            f->theFormula.lLength       = last0+1;
            _PMathObj   lvalue  = f->Compute(0, nameSpace);
            f->theFormula.lLength = stackD;
            if (lvalue->ObjectClass () == MATRIX) {
                mmx = (_Matrix*)lvalue;
            }
            if (lvalue->ObjectClass () == ASSOCIATIVE_LIST) {
                mma = (_AssociativeList*)lvalue;
            }
            last0++;
        } else {

            _Operation * lValue = f->GetIthTerm(0L);

            if (lValue && lValue->IsAVariable(false)) {
             _Variable* mmo = LocateVar(lValue->GetAVariable());
             if (mmo->ObjectClass () == MATRIX) {
                mmx = (_Matrix*)(mmo->GetValue());
                lValue->SetAVariable(-lValue->GetAVariable()-3);
              } else {

                if (mmo->ObjectClass () == ASSOCIATIVE_LIST) {
                  mma = (_AssociativeList*)(mmo->GetValue());
                  lValue->SetAVariable(-lValue->GetAVariable()-3);
                }
              }
            }
        }

        _PMathObj coordMx = nil;
        if (mma || mmx) {
            long expectedType = mmx?MATRIX:STRING;
            coordMx = f->Compute(last0);
            if (!coordMx || coordMx->ObjectClass() != expectedType) {
                if (mmx) {
                    HandleApplicationError ("Matrix expected but not supplied.");
                } else {
                    HandleApplicationError ("String key expected but not supplied.");
                }

                return 0;
            }
        } else {
            HandleApplicationError ("Matrix/List LHS expected but not supplied.");
            return 0;
        }


        if (mmx) { // matrix LHS
            _Matrix * mcoord = (_Matrix*)coordMx;

            long hC = mcoord->theData[0],
                 vC = mcoord->theData[1];

            if (mmx->CheckCoordinates (hC,vC)) {
                if (!ANALYTIC_COMPUTATION_FLAG) {
                    mmx->MStore (hC, vC, newF, (code==HY_FORMULA_FORMULA_VALUE_INCREMENT)?HY_OP_CODE_ADD:HY_OP_CODE_NONE);
                } else {
                    _PMathObj newP = newF.ConstructPolynomial();
                    if (!newP) {
                        HandleApplicationError ("Can't assign non-polynomial entries to polynomial matrices.");
                    } else {
                        mmx->MStore (hC,vC, newP);
                    }
                }
                mmx->CheckIfSparseEnough();
            }
        } else if (mma) { // Associative array LHS
            mma->MStore (coordMx, newF.Compute(), true, (code==HY_FORMULA_FORMULA_VALUE_INCREMENT)?HY_OP_CODE_ADD:HY_OP_CODE_NONE);
        }

        return 1;
    }
    return 0;
}

//__________________________________________________________________________________

struct      characterChecker {
    characterChecker (_String s) {
        for (long r = 0L; r<256L; r++) {
            isAllowed [r] = false;
        }
        for (long r2 = 0; r2<s.sLength; r2++) {
            isAllowed [s.get_uchar (r2)] = true;
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

//__________________________________________________________________________________
long        HandleFormulaParsingError (_String errMsg, _String* saveError, _String& s, long index) {
    if (index >= 0) {
        errMsg = errMsg & " in the following context: '"&s.Cut(MAX(0,index-CONTEXT_TRUNCATION),index)&"<ERROR HERE>"&s.Cut(index+1,MIN (index+CONTEXT_TRUNCATION, s.sLength-1)) & "'";
    }
    if (saveError) {
        *saveError = errMsg;
    } else {
        HandleApplicationError (errMsg);
    }
    return HY_FORMULA_FAILED;
}

//__________________________________________________________________________________
bool        checkLHS (_List* levelOps, _List* levelData, _String& errMsg, char & deref, _Formula * f, _Variable*& lhs, _FormulaParsingContext& context) {
    bool check = true;

    lhs = nil;

    if (f->IsEmpty()) {
    /* nothing has been added to the formula previously, so this should be a simple assignment
       to a variable */
        if (levelOps->lLength == 0) {
            if (levelData->lLength == 0) {
                errMsg = "An empty left-hand side";
                return false;
            }
        }
    } else {
        if (levelData->lLength > 0) {
             errMsg = "Only simple variable references [e.g. var = value or *ref = value or *(string expression) = value] can appear on the LHS of assignments";
             return false;
        }
    }
                  
    deref = kStringDirectReference;
    if (levelOps->lLength > 0) { // this is where 'f is non-empty' cases will go
        check = false;
        if (levelOps->lLength == 1) {
            char buffered_op = ((_Operation*)((*levelOps)(0)))->TheCode();
            if (buffered_op == HY_OP_CODE_MUL) {
                check = true; deref = kStringLocalDeference;
            } else {
                if (buffered_op == HY_OP_CODE_POWER) {
                    check = true; deref = kStringGlobalDeference;
                } else {
                    errMsg = "* and ^ are the two supported de-referencing operations";
                }
            }
        } else {
            errMsg = "Expressions (other than matrix/dict access) cannot appear on the left-hand side of assignments";
        }
    } else {
        if (levelData->lLength != 1) {
            errMsg = "The left hand side expression does not contain an object reference";
            check = false;
        }
    }
    if (check && levelData->lLength == 1) {
        _Operation * theOp = dynamic_cast<_Operation*>(levelData->GetItem(0));

        check = false;
        if (theOp) {
          if (theOp->IsAVariable(false)) {
            lhs = LocateVar(theOp->GetAVariable());
            check = true;
          } else {
            if (theOp->GetANumber() && theOp->GetANumber()->ObjectClass() == STRING) {
              // handle things like ^"id" = value
              if (deref != kStringDirectReference) {
                _hyExecutionContext cntx (context.formulaScope(), context.errMsg());
                lhs = (_Variable*) ((_FString*)theOp->GetANumber())->Dereference(deref == kStringGlobalDeference, &cntx, true);
                if (!lhs) {
                  errMsg = "The left-hand side of an assignment like ^\"id\" must reference an existing variable";
                  return false;
                } else {
                  deref = kStringDirectReference;
                  check = true;
                }
              }
            }
          }
        }

        if (!check) {
          errMsg = "The left-hand side of an assignment must be a variable (not a constant)";
        }
    }
    return check;
}

//__________________________________________________________________________________
long _parserHelperHandleInlineBoundCases (_String& s, _FormulaParsingContext& parsingContext, long i, _Variable* lhs_variable, _Formula * f, char deref, _Formula &newF) {
    _PMathObj varObj = newF.Compute();
    if (varObj->ObjectClass()!=NUMBER) {
        return HandleFormulaParsingError ("Variable bound must evaluate to a number ", parsingContext.errMsg(), s, i);
    }

    long varID;

    if (lhs_variable) {
        varID = DereferenceVariable(lhs_variable->GetAVariable(), parsingContext.formulaScope(), deref);
    } else {
        varID = DereferenceString(f->Compute(0, parsingContext.formulaScope(), nil, parsingContext.errMsg()), parsingContext.formulaScope(), deref);
    }
    if (varID < 0) {
        return HandleFormulaParsingError ("Failed to dereference ", parsingContext.errMsg(), s, i);
    }

    _Variable * theV = (_Variable*)LocateVar(varID);
        
    if (s.get_char(i)=='>') {
        theV->SetBounds(varObj->Value(),theV->GetUpperBound());
    } else {
        theV->SetBounds(theV->GetLowerBound(),varObj->Value());
    }
    return HY_FORMULA_EXPRESSION;
}

//__________________________________________________________________________________
long _parserHelperHandleInlineAssignmentCases (_String& s, _FormulaParsingContext& parsingContext, long i, _Variable* lhs_variable, _Formula * f, char deref, _Formula &newF, bool twoToken) {


    long varID;

    if (lhs_variable) {
        varID = DereferenceVariable(lhs_variable->GetAVariable(), parsingContext.formulaScope(), deref);
    } else {
        varID = DereferenceString(f->Compute(0, parsingContext.formulaScope(), nil, parsingContext.errMsg()), parsingContext.formulaScope(), deref);
    }
    if (varID < 0) {
        return HandleFormulaParsingError ("Failed to dereference ", parsingContext.errMsg(), s, i);
    }
    _Variable * theV = (_Variable*)LocateVar(varID);
    
    if (s.get_char(i-1) != ':') {
        _PMathObj varObj = newF.Compute();
        if (!varObj) {
            return HandleFormulaParsingError ("Invalid RHS in an assignment ", parsingContext.errMsg(), s, i);
        }
        if (twoToken && s.get_char(i-1) == '+') {
            _List arg;
            arg <<  varObj;
            theV->SetValue(theV->Compute()->ExecuteSingleOp(HY_OP_CODE_ADD,&arg));
        } else {
            theV->SetValue(varObj);
        }
    } else {
        theV->SetFormula (newF);
    }
    return HY_FORMULA_EXPRESSION;
}

//__________________________________________________________________________________

void        _parse_new_level (long & level, _List & operations, _List& operands, _List*& levelOps, _List*& levelData, _String& curOp, _SimpleList& functionCallTags, long function_tag = -1L) {

  level ++;
  operations.AppendNewInstance (new _List);
  operands.AppendNewInstance (new _List);

  levelOps  = (_List*)(operations(level));
  levelData = (_List*)(operands(level));

  functionCallTags << function_tag;
  
  curOp = kEmptyString;
}



//__________________________________________________________________________________
long        Parse (_Formula* f, _String& s, _FormulaParsingContext& parsingContext, _Formula* f2)
/* SLKP 20110908: added the concept of a 'volatile' formula, i.e. something that should be reparsed every time in ExecuteCase0
                : currently those include
                :    inline constructors (matrices, dictionaries)
                :    `` substitutions in strings


   SLKP 20100817: decoupled return code from variable reference return
*/

// returns:

/*

 case                       | return value                                  | parsingContext.assignment_ref_id value

 parse failed               | HY_FORMULA_FAILED                             | undefined
 expresion (no LHS)         | HY_FORMULA_EXPRESSION                         | undefined
 z = x/y                    | HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT          | index of the LHS
 *|^(expr) = expr           | HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT         | undefined
 z := x/y                   | HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT        | index of the LHS
 object[a] := x/y           | HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT         | undefined
 object[a] = x/y            | HY_FORMULA_FORMULA_VALUE_ASSIGNMENT           | undefined
 z :< expr                  | HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT    | index of the LHS
 z :> expr                  | HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT    | index of the LHS

 Further, for (HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT,  HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT,
 HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT, HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT,
 HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT):

 case            |   parsingContext.assignment_ref_type value

 z op x/y        |   kStringDirectReference
 *z op x/y       |   kStringLocalDeference
 ^z op x/y       |   kStringGlobalDeference


*/

{
  //static bool inAssignment = false;

    expressionsParsed++;

    _List           operations,
                    operands,
                    *levelOps,
                    *levelData;

    /*04252006*/

    _SimpleList     squareBrackets,
                    mergeMAccess,
                    mergeMAccessLevel,
                    functionCallTags
                    /*
                      for each context level, stores -1, if this level is not a function call,
                      or the first operation associated with the argument list for this function
                      call
                     */
                    ;



    long            level                 = -1;
    /* 04252006 mlevel = -1, */
    /* mcount = 0 ; */

    _String curOp;

    bool            impliedMult = false,
                    globalKey   = false,
                    twoToken    = false;

    char            storage     = 0;


    _parse_new_level (level, operations, operands, levelOps, levelData, curOp, functionCallTags);
    for (long i = 0; i<=s.sLength; i++) {
        storage = 0; // no implied ops by default

        if (isspace(s.get_char(i))) { // skip spaces and tabs
            continue;
        }
      
        char     lookAtMe = s.get_char(i);
        //printf ("at '%c', the formula looks like this %s\n", lookAtMe, (const char*)_String ((_String*)f->GetList().toStr()));

        if (i==s.sLength || lookAtMe == ')' || lookAtMe == ']' || lookAtMe == ',') {
            // closing ) or ]
            // or a parameter list
                /* 04252006 if (level == mlevel && s.get_char(i)!=']')*/
            if (squareBrackets.lLength && squareBrackets.lData[squareBrackets.lLength-1] == level && lookAtMe != ']') {
                return HandleFormulaParsingError ("Missing or unbalanced '[]' ", parsingContext.errMsg(), s, i);
            }

            if (lookAtMe != ',') {
              if (i != s.sLength) {
                level--;
              } else {
                if (level != 0L) {
                  level = -1L;
                }
              }
            }

            if (level<0) {
                return HandleFormulaParsingError ("Unbalanced '()' parentheses ", parsingContext.errMsg(), s, i);
            }

            if (lookAtMe ==',' && (level<1 || (squareBrackets.lLength && squareBrackets.lData[squareBrackets.lLength-1] == level))) {
                return HandleFormulaParsingError ("Parameter list is out of context ", parsingContext.errMsg(), s, i);
            }

            if (levelOps->lLength) { // there are some buffered operations left
                if (levelOps->lLength > 3 || levelData->lLength > 2) {
                    return HandleFormulaParsingError ("Syntax error ", parsingContext.errMsg(), s, i);
                }

                for (unsigned long i = 0UL; i<levelData->countitems(); i++) {
                    f->PushTerm (levelData->GetItem (i));
                }

                for (long k = levelOps->countitems()-1L; k>=0; k--) {
                    f->PushTerm (levelOps->GetItem(k));
                }

                levelOps->Clear();
            } else {
                if (levelData->lLength>1) {
                    return HandleFormulaParsingError ("Syntax error ", parsingContext.errMsg(), s, i);
                } else if (levelData->lLength) {
                  //f->theFormula << (*levelData)(0);    // mod 07072006 to not duplicate
                  f->PushTerm (levelData->GetItem(0));
                }

            }

            levelData->Clear();

            if (i<s.sLength && lookAtMe !=',' ) {
                operations.Delete (level+1);
                operands.Delete   (level+1);
                functionCallTags.Pop();

                levelOps          = (_List*)operations(level);
                levelData         = (_List*)operands(level);

                long function_call_pop = functionCallTags.Element (-1);

                if (function_call_pop >= 0) {
                  if (levelData->lLength == 0L && levelOps->lLength == 1L) {
                    long argument_count = f->StackDepth(function_call_pop);

                    _Operation* function_call = (_Operation*)levelOps->GetItem(0);
                    if (function_call->IsHBLFunctionCall()) {
                      function_call->SetTerms(-argument_count-1L);
                    } else {
                      function_call->SetTerms(argument_count);
                    }

                    f->PushTerm (function_call);

                    operations.Delete (level);
                    operands.Delete   (level);
                    functionCallTags.Pop();

                    levelOps          = (_List*)operations(--level);
                    levelData         = (_List*)operands(level);
                 } else {
                    return HandleFormulaParsingError ("Syntax error ", parsingContext.errMsg(), s, i);
                  }
                }

                if (lookAtMe !=']')
                    if ( BinOps.Find(s.get_char(i+1))==-1 && i<s.sLength-1 && s.get_char(i+1)!=')' && s.get_char(i+1)!=']' && s.get_char(i+1)!='[' && HalfOps.Find(s.get_char(i+1))==-1 && s.get_char(i+1)!=',') {
                        storage = s.get_char(i);
                        s.set_char(i,'*');
                    }
            }

            if (lookAtMe ==']') {
                if (!squareBrackets.lLength || squareBrackets.lData [squareBrackets.lLength-1] != level + 1) {
                    return HandleFormulaParsingError ("Unexpected ']' ", parsingContext.errMsg(), s, i);
                }
                squareBrackets.Delete(squareBrackets.lLength-1);
                curOp = *(_String*)BuiltInFunctions(HY_OP_CODE_MACCESS);
                if (mergeMAccess.lLength && mergeMAccess.lData[mergeMAccess.lLength-1] >= 0 && mergeMAccessLevel.lData[mergeMAccessLevel.lLength-1] == level) {
                    long mergeIndex              = mergeMAccess.lData[mergeMAccess.lLength-1];
                    _Operation * previousMaccess = (_Operation*) f->theFormula (mergeIndex);
                    if (previousMaccess->GetCode () != curOp) {
                        return HandleFormulaParsingError ("Internal error in Parse. Incorrect matrix access token code ", parsingContext.errMsg(), s, i);
                    }

                    if (previousMaccess->GetNoTerms() > 2) {
                        mergeMAccess.Delete (mergeMAccess.lLength-1,false);
                        mergeMAccessLevel.Delete (mergeMAccessLevel.lLength-1,false);
                        f->theFormula.AppendNewInstance(new _Operation (curOp ,2));
                    } else {
                        previousMaccess->SetTerms(3);
                        mergeMAccess.Delete (mergeMAccess.lLength-1,false);
                        mergeMAccessLevel.Delete (mergeMAccessLevel.lLength-1,false);
                        previousMaccess->AddAReference();
                        f->theFormula.Delete (mergeIndex);
                        f->theFormula.AppendNewInstance(previousMaccess);
                    }
                } else {
                    f->theFormula.AppendNewInstance(new _Operation (curOp ,2));
                }
            }

            if (!storage) {
                continue;
            }
        }


        if (s.get_char(i) == '=' && s.get_char(i+1) != '=' && (!twoToken || s.get_char(i-1)==':' || s.get_char (i-1) == '+')) { // assignment operator
            _String  errMsg;

            bool check               = !parsingContext.inAssignment(),
                 is_array_assignment = f->IsArrayAccess() && levelOps->lLength == 0UL;


            char deref = 0;

            _Variable *lhs_variable = nil;


            if (check) {
                if (is_array_assignment) {
                    (((_Operation*)((f->theFormula)(f->theFormula.lLength-1)))->TheCode()) = HY_OP_CODE_MCOORD;
                } else {
                    check = checkLHS (levelOps, levelData, errMsg, deref, f, lhs_variable, parsingContext);
                }
            } else {
                errMsg = "Can't assign within another assignment";
            }


            if (!check) {
                return HandleFormulaParsingError (errMsg, parsingContext.errMsg(), s, i);
            }

            parsingContext.inAssignment() = true;
            _String ss (s,i+1,-1); // this is the RHS
            _Formula  newF;

            if (Parse(&newF,ss,parsingContext, f2) != HY_FORMULA_EXPRESSION) {
                parsingContext.inAssignment() = false;
                return HY_FORMULA_FAILED;
            }
            parsingContext.inAssignment() = false;
            if (!is_array_assignment && lhs_variable)
                // normal variable assignment
            {
                if (!f2) { // immediate execution
                    if (_parserHelperHandleInlineAssignmentCases (s,parsingContext, i, lhs_variable, f,  deref, newF, twoToken) == HY_FORMULA_FAILED) {
                        return HY_FORMULA_FAILED;
                    }
                } else { // this gets called from ExecuteCase0...
                    if (twoToken && s.get_char(i-1) == '+') { // += gets handled here
                    
                        _Operation* self = new _Operation ();
                        self->SetAVariable(lhs_variable->GetAVariable());
                        newF.theFormula.InsertElement (self,0,false);
                        DeleteObject (self);
                        if (deref != kStringDirectReference) {
                             _Operation* ref = new _Operation (*(_String*)BuiltInFunctions(deref == kStringGlobalDeference ? HY_OP_CODE_POWER : HY_OP_CODE_MUL),1);
                             newF.theFormula.InsertElement (ref,1,false);
                             DeleteObject (ref);
                      }
                      newF.theFormula.AppendNewInstance (new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                    }
                    f->Duplicate((_Formula const*)&newF);
                }

                parsingContext.assignmentRefID()   = lhs_variable->GetAVariable();
                parsingContext.assignmentRefType() = deref;

                return (s.get_char(i-1)==':')?HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT:HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT;
            } else
                // matrix/associative array element assignment
            {
                if (is_array_assignment) {
                    long stackD = -1,
                         last0  = 0;

                    for (unsigned long opID = 0; opID < f->theFormula.lLength - 1; opID ++) {
                        ((_Operation*)f->theFormula(opID)) -> StackDepth (stackD);
                        if (stackD == 0) {
                            last0 = opID;
                        }
                    }

                    if (!f2) { // immediate execution
                        bool       anError = false;

                        if (newF.IsAConstant() || s.get_char(i-1) !=':') {
                            _PMathObj       currentValue = (_PMathObj)newF.Compute();
                            currentValue->AddAReference();
                            newF.theFormula.Clear();
                            newF.theFormula.AppendNewInstance (new _Operation(currentValue));
                        }


                        _Matrix          * mmx = nil;
                        _AssociativeList * mma = nil;

                        if (last0 > 0) {
                            stackD = f->theFormula.lLength;
                            f->theFormula.lLength   = last0+1;
                            _PMathObj   lvalue      = f->Compute();
                            f->theFormula.lLength   = stackD;

                            if (lvalue->ObjectClass () == MATRIX) {
                                mmx = (_Matrix*)lvalue;
                            } else if (lvalue->ObjectClass () == ASSOCIATIVE_LIST) {
                                mma = (_AssociativeList*)lvalue;
                            }

                            last0++;
                        } else {
                            _Variable* mmo = ((_Operation*)f->theFormula(0))->IsAVariable()?LocateVar(((_Operation*)f->theFormula(0))->GetAVariable()):nil;

                            if (mmo) {
                              if (mmo->ObjectClass () == MATRIX) {
                                mmx = (_Matrix*)(mmo->GetValue());
                                ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                              } else {

                                if (mmo->ObjectClass () == ASSOCIATIVE_LIST) {
                                  mma = (_AssociativeList*)(mmo->GetValue());
                                  ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                                }
                              }
                            }
                        }

                        if (mmx) {
                            _Matrix  *mcoord;
                            _PMathObj coordMx = f->Compute(last0);

                            if (!coordMx|| coordMx->ObjectClass()!=MATRIX) {
                                anError = true;
                            } else {
                                mcoord = (_Matrix*)coordMx;
                                _Constant hC ((*mcoord)[0]),
                                          vC ((*mcoord)[1]);

                                mmx->MStore (&hC, &vC, newF, (twoToken && s.get_char(i-1) =='+')?HY_OP_CODE_ADD:HY_OP_CODE_NONE);
                            }
                        } else if (mma) {
                            _PMathObj coordIdx = f->Compute(last0);

                            if (!coordIdx|| coordIdx->ObjectClass() != STRING ) {
                                anError = true;
                            } else {
                                mma->MStore (coordIdx, newF.Compute(),true, (twoToken && s.get_char(i-1) =='+')?HY_OP_CODE_ADD:HY_OP_CODE_NONE);
                            }
                        } else {
                            anError = true;
                        }


                        if (anError) {
                            return HandleFormulaParsingError ("Invalid matrix/associative list ident supplied ", parsingContext.errMsg(), s, i);
                        }

                        return HY_FORMULA_EXPRESSION;
                    } else {
                        bool isSimple = (s.get_char(i-1) != ':');
                        f2->Duplicate   (&newF);
                        if (last0 == 0) {
                            ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                        }
                        return isSimple?((s.get_char(i-1) == '+')?HY_FORMULA_FORMULA_VALUE_INCREMENT:HY_FORMULA_FORMULA_VALUE_ASSIGNMENT):HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT;
                    }
                } else {
                // *(expression) reference
                    if (f2) {
                        bool isSimple = (s.get_char(i-1) != ':');

                        if (twoToken && s.get_char(i-1) == '+') { // += gets handled here
                            newF.theFormula.InsertElement (new _Operation (*(_String*)BuiltInFunctions(deref == kStringGlobalDeference ? HY_OP_CODE_POWER : HY_OP_CODE_MUL),1), 0, false);
                            for (long lhs_ops = f->theFormula.lLength-1; lhs_ops >= 0; lhs_ops --) {  
                                newF.theFormula.InsertElement (f->theFormula(lhs_ops), 0, true);
                            }
                            newF.theFormula.AppendNewInstance (new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                        }
                        f2->Duplicate   (&newF);
                        parsingContext.assignmentRefType() = deref;
                        return isSimple?HY_FORMULA_REFERENCE_VALUE_ASSIGNMENT:HY_FORMULA_REFERENCE_FORMULA_ASSIGNMENT;

                    } else {
                        if (_parserHelperHandleInlineAssignmentCases (s,parsingContext, i, lhs_variable, f,  deref, newF, twoToken) == HY_FORMULA_FAILED) {
                            return HY_FORMULA_FAILED;
                        }
                    }
                }
            }
        }

        if ( s.get_char(i-1)==':' && (s.get_char(i)=='<' || s.get_char(i)=='>')) { // variable bounds
            _Variable * lhs = nil;
            _String errMsg;
            char    deref;

            if (parsingContext.inAssignment()||f->IsArrayAccess()||! checkLHS (levelOps, levelData, errMsg, deref, f, lhs, parsingContext)) {
               return HandleFormulaParsingError ("Can't set bounds like this ", parsingContext.errMsg(), s, i);
            }

            parsingContext.inAssignment() = true;

            _String ss (s,i+1,-1);
            _Formula newF;

            if (Parse(&newF,ss,parsingContext,f2) != HY_FORMULA_EXPRESSION) {
                parsingContext.inAssignment() = false;
                return HY_FORMULA_FAILED;
            }


            parsingContext.inAssignment() = false;



            if (lhs) {
                if (!f2) {
                    if (_parserHelperHandleInlineBoundCases (s,parsingContext,i,lhs,f, deref, newF) == HY_FORMULA_FAILED) {
                        return HY_FORMULA_FAILED;
                    }
                } else { // BOUND ASSIGNMENTS
                    f2->Duplicate   (&newF);

                    parsingContext.assignmentRefID()   = lhs->GetAVariable();
                    parsingContext.assignmentRefType() = deref;

                    return (s.get_char(i)=='>')?HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT:HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT;
                }
            } else {
                if (!f2) {
                    if (_parserHelperHandleInlineBoundCases (s,parsingContext,i,lhs,f, deref, newF) == HY_FORMULA_FAILED) {
                        return HY_FORMULA_FAILED;
                    }
                } else { // BOUND ASSIGNMENTS
                    f2->Duplicate   (&newF);
                    parsingContext.assignmentRefType() = deref;
                    return (s.get_char(i)=='>')?HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT:HY_FORMULA_REFERENCE_LOWER_BOUND_ASSIGNMENT;
                }

            }

            return HY_FORMULA_EXPRESSION;
        }


        if (s.get_char(i) == '{') // a matrix
            /* 20090803 SLKP:
                 fixed the code to deal with
            */
        {

            parsingContext.isVolatile() = true;
            
            int     j       = s.ExtractEnclosedExpression (i,'{','}',fExtractRespectQuote | fExtractRespectEscape);

            if (j<0) {
                return HandleFormulaParsingError ("Poorly formed matrix/associative array construct ", parsingContext.errMsg(), s, i);
            }

            _String matrixDef   (s,i,j);
            bool has_values = false;

            if (matrixDef.sLength == 2 || (has_values = matrixDef.FindTerminator(1, ":") >= 0L) || matrixDef.FirstNonSpaceIndex(1,-1) == matrixDef.sLength - 1) {
                _AssociativeList *theList = new _AssociativeList ();
                if (has_values) {
                    matrixDef.Trim (1,matrixDef.sLength-2);
                    if (!theList->ParseStringRepresentation (matrixDef,parsingContext)) {
                        return HandleFormulaParsingError ("Poorly formed associative array construct ", parsingContext.errMsg(), s, i);
                    }
                }

                levelData->AppendNewInstance (new _Operation (theList));
            } else {
                levelData->AppendNewInstance (new _Operation ( new _Matrix (matrixDef,false,parsingContext.formulaScope())));
            }

            i = j;
            continue;

        }

        if (s.get_char(i) == '[') { // opening [
            long  lastCode = -1;

            if (!f->IsEmpty()) {
                lastCode = ((_Operation*)((f->theFormula)(f->theFormula.lLength-1)))->TheCode();
            }

            if (lastCode == HY_OP_CODE_MACCESS && s.get_char(i-1) == ']') {
                mergeMAccess << f->theFormula.lLength-1;
                mergeMAccessLevel << level;
            } else {
                if (levelData->lLength == 0 && f->IsEmpty()) {
                   return HandleFormulaParsingError ("[..] must be preceded by an object to index ", parsingContext.errMsg(), s, i);
                }

                if (levelData->lLength) {
                    //f->theFormula.AppendNewInstance((*levelData)[levelData->lLength-1]);
                    f->PushTerm(levelData->GetItem(levelData->lLength-1));
                    levelData->Delete(levelData->lLength-1);
                }
            }


            _parse_new_level (level, operations, operands, levelOps, levelData, curOp, functionCallTags);
            squareBrackets << level;
            continue;
        }


        if (s.get_char(i) == '(') { // opening (
          // check to see if this is a function call

          _parse_new_level (level, operations, operands, levelOps, levelData, curOp, functionCallTags);
          continue;
        }

        if (s.get_char(i)=='"' || s.get_char (i) == '\'') { // a string literal
            long j             = 1,
                 inPlaceID     = -1;
          
          char terminator = s.get_char (i);

            _String * literal = new _String (16,true);
            _List * formula_list = nil;

            while (i+j<s.sLength) {
                char aChar = s.sData[i+j];
                if (aChar =='\\') {
                    if (i+j+1<s.sLength) {
                        if (s.sData[i+j+1]=='"' ||s.sData[i+j+1]=='`' ||  s.sData[i+j+1]=='\'') {
                            j++;
                            (*literal)<<s.sData[i+j++];
                        } else {
                            (*literal)<<s.sData[i+j++];
                            (*literal)<<s.sData[i+j++];
                        }
                    }
                    continue;
                }

                if (aChar == terminator && inPlaceID < 0) {
                    break;
                }


                if (aChar == '`') {
                    if (inPlaceID < 0) {
                        inPlaceID = ++j;
                    } else if (j == inPlaceID) {
                      //return HandleFormulaParsingError ("Attempted to string substitute an empty quotation ", parsingContext.errMsg(), s, i);
                      (*literal) << '`';
                      inPlaceID = -1L;
                      j++;

                    } else {
                        _String     inPlaceVID (s,i+inPlaceID,i+j-1);
                        _Formula    expressionProcessor;

                        long parse_result = Parse(&expressionProcessor, inPlaceVID, parsingContext, nil);

                        if (parse_result != HY_FORMULA_EXPRESSION) {
                          return HandleFormulaParsingError ("Not a valid/simple expression inside `` ", parsingContext.errMsg(), s, i);
                        }
                      
                      
                        if (expressionProcessor.IsConstant(true)) {
                          _PMathObj constant_literal = expressionProcessor.Compute (0, nil, nil, nil, STRING);
                          if (constant_literal) {
                            (*literal) << *((_FString*)constant_literal)->theString;
                          }
                          else {
                            return HandleFormulaParsingError ("Constant expression inside `` did not evaluate to a string ", parsingContext.errMsg(), s, i);
                          }
                        } else {
                          // push this expression down
                          if (!formula_list) {
                            formula_list = new _List;
                          }

                          literal->Finalize();
                          if (literal->sLength) {
                            formula_list->AppendNewInstance(new _Operation (new _FString(literal)));
                          } else {
                            DeleteObject (literal);
                          }
                          if (formula_list->lLength > 1L) {
                            formula_list->AppendNewInstance(new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                          }
                          literal = new _String (16,true);
                          (*formula_list) << expressionProcessor.theFormula;
                          if (formula_list->lLength > expressionProcessor.theFormula.lLength) {
                            formula_list->AppendNewInstance(new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                          }
                        }
                        inPlaceID = -1L;
                        //parsingContext.isVolatile() = true;
                        j++;


                    }

                } else {
                    if (inPlaceID < 0) {
                        (*literal)<<s.sData[i+j];
                    }
                    j++;
                }
            }


            literal->Finalize();

            if (formula_list) {
              if (literal->sLength) {
                formula_list->AppendNewInstance(new _Operation (new _FString(literal)));
                formula_list->AppendNewInstance(new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
              } else {
                DeleteObject (literal);
              }
              formula_list->AppendNewInstance(new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_MCOORD),1));

              levelData->AppendNewInstance(new _List(*formula_list));
              DeleteObject (formula_list);

            } else {
              levelData->AppendNewInstance (new _Operation (new _FString(*literal)));
              DeleteObject(literal);
            }


            if (inPlaceID >= 0) {
              return HandleFormulaParsingError ("Unterminated string substitution (``) inside a literal ", parsingContext.errMsg(), s, i);
            }

            i += j;
            continue;

        }

        if (alpha.isAllowed [(unsigned char)s.get_char(i)]) { // an identifier
            bool takeVarReference = false;

            if (twoToken) {
                char opChar = s.get_char(i-1);
                if (((_String*)BuiltInFunctions(HY_OP_CODE_REF))->Equal(opChar)) {
                    takeVarReference = true;
                    twoToken = false;
                } else {
                    _String thisOp (opChar);
                    levelOps->AppendNewInstance (new _Operation (thisOp,1L));
                }
            }
            
            impliedMult = (i && numeric.isAllowed [(unsigned char)s.get_char(i-1)]);

            long j = 1;
            while ( i+j<s.sLength && (alpha.isAllowed [(unsigned char)s.get_char(i+j)]|| numeric.isAllowed [(unsigned char)s.get_char(i+j)]) ) {
                j++;
            }


            curOp =  s.Cut(i,i+j-1L);
            i+=j-1L;

            if (curOp.Equal(&globalToken)) {
                if (takeVarReference) {
                    return HandleFormulaParsingError (_String("Cannot make a reference from a reserved word ") & globalToken, parsingContext.errMsg(), s, i);
                }
                globalKey = true;
                continue;
            }

            bool noneObject = false;
            if (curOp.Equal(&noneToken)) {
                 if (takeVarReference) {
                    return HandleFormulaParsingError (_String("Cannot make a reference from a reserved word ") & noneToken, parsingContext.errMsg(), s, i);
                }
               noneObject = true;
                globalKey  = true;
            }

            if (UnOps.FindKey(curOp)>=0) { // a standard function
                if (takeVarReference) {
                    return HandleFormulaParsingError ("Cannot make a reference from a built-in function", parsingContext.errMsg(), s, i);
                }

                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else { // a variable
                // check if this is a function defined  in the list of "standard functions"
                long bLang = noneObject?-1:FunctionNameList.BinaryFindObject (&curOp);
                if (bLang>=0) {
                    if (takeVarReference) {
                        return HandleFormulaParsingError ("Cannot make a reference from a built-in function", parsingContext.errMsg(), s, i);
                    }
                    // built-in function
                    _Operation * built_in_call = new _Operation (curOp,FunctionArgumentCount(bLang));
                    _parse_new_level (level, operations, operands, levelOps, levelData, curOp, functionCallTags, f->NumberOperations());
                    levelOps->AppendNewInstance (built_in_call);
                    continue;
                }

                // check if this is a function defined in the batch language

                if ((bLang =  noneObject?-1:FindBFFunctionName (curOp, parsingContext.formulaScope()))>=0) {
                    if (takeVarReference) {
                        return HandleFormulaParsingError ("Cannot make a reference from user-defined function", parsingContext.errMsg(), s, i);
                    }
                    // HBL function
                    _Operation * hbl_call = new _Operation (curOp,-bLang-1);
                    _parse_new_level (level, operations, operands, levelOps, levelData, curOp, functionCallTags, f->NumberOperations());
                    levelOps->AppendNewInstance (hbl_call);
                    continue;
                }

                long curOpl = curOp.sLength;
                if (curOpl>2 && curOp[curOpl-1]=='_' && curOp[curOpl-2]=='_') { // instant variable refrence
                    _String realVarName (curOp,0,curOpl-3);

                    realVarName = parsingContext.contextualizeRef (realVarName);


                    long realVarLoc = LocateVarByName (realVarName);
                    if (realVarLoc<0) { // bad instant variable reference
                        return HandleFormulaParsingError ("Attempted to evaluate an undeclared variable ", parsingContext.errMsg(), s, i);
                     }
                    if (!f2) { // 03/25/2004 ? Confused why the else
                        levelData->AppendNewInstance(new _Operation((_MathObject*)FetchVar (realVarLoc)->Compute()->makeDynamic()));
                    } else {
                        _Operation * variable_op = new _Operation (true, realVarName, globalKey, parsingContext.formulaScope());
                        variable_op -> SetTerms(-variableNames.GetXtra (realVarLoc)-1);
                        variable_op ->SetAVariable(-2);
                        (*levelData) < variable_op;
                    }
                } else {
                    if (noneObject)
                        levelData->AppendNewInstance (new _Operation (false, curOp));
                    else
                        if (parsingContext.formulaScope() && _hy_application_globals.Find(&curOp) >= 0) {
                            levelData->AppendNewInstance (new _Operation(true, curOp, globalKey, nil, takeVarReference));
                        } else {
                            levelData->AppendNewInstance (new _Operation(true, curOp, globalKey, parsingContext.formulaScope(), takeVarReference));
                        }
                }
                globalKey = false;
                if (impliedMult) {
                    storage = s.get_char(i);
                    s.set_char(i,((_String*)BuiltInFunctions(HY_OP_CODE_MUL))->get_char(0));
                } else if (s.get_char(i+1)=='(') {
                    if (!storage) {
                        storage = s.get_char(i);
                        s.set_char(i,((_String*)BuiltInFunctions(HY_OP_CODE_MUL))->get_char(0));
                    } else {
                        curOp = *(_String*)BuiltInFunctions(HY_OP_CODE_MUL);
                        levelOps->AppendNewInstance(new _Operation (curOp,2));
                    }
                }
                if (!storage) {
                    continue;
                }
            }
        }

        if (numeric.isAllowed [(unsigned char)s.get_char(i)]) {
            if (twoToken) {
                _String thisOp (s.get_char(i-1));
                levelOps->AppendNewInstance (new _Operation (thisOp,1L));
            }
            long j = 1;

            while ( i+j<s.sLength && (numeric.isAllowed [(unsigned char)s.get_char(i+j)] || ((s.get_char(i+j)=='-' || s.get_char(i+j)=='+' )&& tolower(s.get_char(i+j-1))=='e')) ) {
                j++;
            }

            curOp =  (s.Cut(i,i+j-1));
            i+=j-1;
            levelData->AppendNewInstance (new _Operation (false, curOp));
            if (i<s.sLength-1 && s.get_char(i+1)=='(') {
                storage = s.get_char(i);
                s.set_char(i,((_String*)BuiltInFunctions(HY_OP_CODE_MUL))->get_char(0));
            } else {
                continue;
            }
        }
   
        if ( BinOps.Find (s.get_char(i)) != -1L || (twoToken&&  _Operation::BinOpCode (s, i) != -1L)) {
          
            bool look_ahead = _Operation::BinOpCode (s, i+1) != -1L;

            if (!twoToken && look_ahead) {
                twoToken = true;
                continue;
            }

            if (twoToken|| look_ahead) {
                if (!twoToken) {
                    i++;
                }
                curOp = _String(s.get_char(i-1)) & s.get_char(i);
            } else {
                curOp = s.get_char(i);
            }

            long twoOrOne = 2;

            if (storage) {
                s.set_char(i,storage);
            }

            if (levelData->countitems()==0) {
                if (s[i-curOp.sLength]!=')' && storage!=')' && s[i-curOp.sLength] !=']') {
                    if (!twoToken && UnOps.FindKey (s.get_char(i)) >= 0) {
                        twoOrOne = 1;
                    } else {
                        return HandleFormulaParsingError ("Bad binary operator placement ", parsingContext.errMsg(), s, i);
                    }
                }
            }

            twoToken = false;

            if (levelData->countitems()) {
                if (storage) {
                    BaseRef newS = (*levelData)(levelData->countitems()-1)->makeDynamic();
                    for (unsigned long k = 0; k<levelData->countitems()-1; k++) {
                        f->PushTerm (levelData->GetItem (k));
                        //f->theFormula << ((*levelData)(k));
                    }

                    levelData->Clear();
                    levelData->AppendNewInstance (newS);
                } else {
                    for (unsigned long k = 0; k<levelData->countitems(); k++) {
                        f->PushTerm (levelData->GetItem (k));
                        //f->theFormula << ((*levelData)(k));
                    }
                    levelData->Clear();
                }
            }

            if (!levelOps->countitems()) {
                levelOps->AppendNewInstance (new _Operation (curOp,twoOrOne));
                if (terminate_execution) {
                    return HY_FORMULA_FAILED;
                }
                continue;
            }

            // check operation precedence

            long h,g;
            _String prevOp = *((((_Operation*)((*levelOps)(levelOps->countitems()-1)))->GetCode()));

            h = _Operation::BinOpCode (prevOp);
            g = _Operation::BinOpCode (curOp);

            if (h >= 0L) {
                h = opPrecedence (h);
            }
            g = opPrecedence (g);


           if (g>h && h!=-1) { // store the op, don't do it yet!
                levelOps->AppendNewInstance (new _Operation (curOp,twoOrOne));
                if (terminate_execution) {
                    return HY_FORMULA_FAILED;
                }
                continue;
            }

            // do the stored operations now

            for (int j = levelOps->countitems()-1; j>=0; j--) {
                _String  sss = (((_Operation*)((*levelOps)(levelOps->countitems()-1)))->GetCode());
                h = _Operation::BinOpCode (sss);
                if (h < 0L) {
                    h = 0xFFFF;
                } else {
                    h = opPrecedence (h);
                }

                if (h<g) {
                    break;
                }
                f->theFormula << levelOps->GetItem(j);
                levelOps->Delete(levelOps->lLength-1);
            }
            levelOps->AppendNewInstance (new _Operation (curOp,twoOrOne));
            if (terminate_execution) {
                return HY_FORMULA_FAILED;
            }
            continue;
        } else if (UnOps.FindKey (s.get_char(i)) >= 0) {
            if ((s.get_char(i)=='-' || s.get_char(i)=='+') && (!i|| s.get_char(i-1)=='(')) { // unary minus or plus
                curOp   = s.get_char(i);
                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else {
                if (HalfOps.contains(s.get_char(i))) {
                    twoToken = true;
                    continue;
                }
                return HandleFormulaParsingError ("Bad binary operator placement ", parsingContext.errMsg(), s, i);
            }
        } else {
            if (!HalfOps.contains(s.get_char(i))) {
                return HandleFormulaParsingError ("Unexpected symbol ", parsingContext.errMsg(), s, i);
            } else {
                twoToken = true;
            }

        }
    }
    return HY_FORMULA_EXPRESSION;
}
//__________________________________________________________________________________

long     VerbosityLevel (void)
{
    checkParameter (VerbosityLevelString, verbosityLevel, -1L);
    return verbosityLevel;
}


//__________________________________________________________________________________
void  stashParameter (_String const& name, hyFloat v, bool set)
{
    static  hyFloat stash = 0.0;

    long f = LocateVarByName (name);
    if (f>=0) {
        _Variable *thisV = FetchVar(f);
        if (set) {
            stash = thisV->Value();
            _Constant dummy (v);
            thisV->SetValue (&dummy);
        } else {
            _Constant dummy (stash);
            thisV->SetValue (&dummy);
        }
    } else if (set) {
        stash = v;
        setParameter (name,v);
    }
}


//__________________________________________________________________________________
void  setParameter (_String const & name, hyFloat def, _String* namespc)
{
    if (namespc) {
        _String namespcd = AppendContainerName(name,namespc);
        setParameter (namespcd,def);
    } else {
        long f = LocateVarByName (name);
        if (f<0) {
            _Variable cornholio(name);
            setParameter (name,def);
        } else {
            FetchVar(f)->SetValue(new _Constant (def), false);
        }
    }
}

//__________________________________________________________________________________

void  setParameter (_String const& name, _PMathObj def, _String* namespc, bool dup)
{
    if (namespc) {
        _String namespcd = AppendContainerName(name,namespc);
        setParameter (namespcd,def,nil, dup);
    } else {
        long f = LocateVarByName (name);
        if (f<0) {
            _Variable cornholio(name);
            setParameter (name,def,nil, dup);
        } else {
            FetchVar(f)->SetValue(def,dup);
        }
    }
}

//__________________________________________________________________________________

void ExportIndVariables (_StringBuffer& glVars, _StringBuffer& locVars, _SimpleList* indepVarList)
{
    _StringBuffer * stIn;
    _String str;

    for (unsigned long   i=0; i<indepVarList->lLength; i++) {
        _Variable *thisVar = LocateVar(indepVarList->lData[i]);
        if (thisVar->IsGlobal()) {
            str = _String ("\nglobal ") & *thisVar->GetName() & '=' & _String((_String*)parameterToString(thisVar->Compute()->Value())) & ';';
            stIn = &glVars;
        } else {
            str = _String ("\n") & *thisVar->GetName() & '=' & _String((_String*)parameterToString(thisVar->Compute()->Value())) & ';';
            stIn = &locVars;
        }
        *stIn << str;
        if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
            str = _String ("\n") & *thisVar->GetName() & ":>" & _String((_String*)parameterToString(thisVar->GetLowerBound())) & ';';
            *stIn << str;
        }
        if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
            str = _String ("\n") & *thisVar->GetName() & ":<" & _String((_String*)parameterToString(thisVar->GetUpperBound())) & ';';
            *stIn << str;
        }
    }
}

//__________________________________________________________________________________

void ExportDepVariables (_StringBuffer& glVars, _StringBuffer& locVars, _SimpleList* depVarList)
{
    if (depVarList->lLength) {
        _StringBuffer * stIn;
        _String str;

      /* first we have to reorder global variables, so that dependent global variables which depend
           on other dependent global variables are written afterwards (lest they be implicitly declared
           as local).
           The algorithm is very ugly, but since there are only a few global dependent variables (in general...) */

        _SimpleList     _globalVariablesList,
                        lfDepGlobs,
                        tl1;

        _List           dependancyLists;
        {
            for (unsigned long i=0; i<depVarList->lLength; i++)
                if (LocateVar(depVarList->lData[i])->IsGlobal()) {
                    lfDepGlobs << depVarList->lData[i];
                }
        }
        lfDepGlobs.Sort();

        for (unsigned long i=0; i<depVarList->lLength; i++) {
            _Variable * thisVar = LocateVar(depVarList->lData[i]);
            if (thisVar->IsGlobal()) {
                _SimpleList                 globDependancyList,
                                            prunedList;

                _AVLList                    globDependancyListAVL (&globDependancyList);

                thisVar->ScanForVariables (globDependancyListAVL,true);

                globDependancyListAVL.ReorderList ();

                prunedList.Intersect (globDependancyList,lfDepGlobs);

                if (prunedList.lLength) {
                    _globalVariablesList << i;
                    dependancyLists && & prunedList;
                    continue;
                }
                str = _String("\nglobal ") & *thisVar->GetName();
                stIn = &glVars;
            } else {
                str = _String("\n") & *thisVar->GetName();
                 stIn = &locVars;
            }
            (*stIn)<<str;
            (*stIn)<<":=";
            stIn->AppendNewInstance(thisVar->GetFormulaString());
            (*stIn)<<';';
            if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
                str = _String ("\n") & *thisVar->GetName() & ":>" & _String((_String*)parameterToString(thisVar->GetLowerBound())) & ';';
                (*stIn)<<str;
            }
            if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
                str = _String ("\n") & *thisVar->GetName() & ":<" & _String((_String*)parameterToString(thisVar->GetUpperBound())) & ';';
                (*stIn)<<str;
            }
        }

        if (_globalVariablesList.lLength)
            // check internal dependancies
        {
            _SimpleList writeOrder (_globalVariablesList.lLength,0,1),
                        indexList  (_globalVariablesList.lLength,0,1);


            for (unsigned long i2 = 0; i2 < _globalVariablesList.lLength; i2++) {
                long updatedIndex = writeOrder.lData[i2];
                _SimpleList * depList = (_SimpleList*)dependancyLists(i2);
                for (unsigned long i3 = 0; i3 < depList->lLength; i3 ++) {
                    long i4 = _globalVariablesList.Find (depList->lData[i3]);
                    if (i4 >= 0 && updatedIndex < writeOrder.lData[i4]) {
                        updatedIndex = writeOrder.lData[i4] + 1;
                    }
                }
                writeOrder.lData[i2] = updatedIndex;
            }

            SortLists (&writeOrder, &indexList);

            for (unsigned long i=0; i<_globalVariablesList.lLength; i++) {
                _Variable * thisVar = LocateVar(depVarList->lData[_globalVariablesList.lData[indexList.lData[i]]]);
                str = _String("\nglobal ") & *thisVar->GetName();
                glVars<<str;
                glVars<<":=";
                glVars<< thisVar->GetFormulaString();
                glVars<<';';
                if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
                    str = _String ("\n") & *thisVar->GetName() & ":>" & _String((_String*)parameterToString(thisVar->GetLowerBound())) & ';';
                    glVars<<str;
                }
                if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
                    str = _String ("\n") & *thisVar->GetName() & ":<" & _String((_String*)parameterToString(thisVar->GetUpperBound())) & ';';
                    glVars<<str;
                }
            }
        }
    }
}

//__________________________________________________________________________________

void ExportCatVariables (_StringBuffer & rec, _SimpleList* catVarList) {
    _SimpleList     nonInd;

    for (long idx = 0; idx < catVarList->lLength; idx++)
        if (((_CategoryVariable*)LocateVar(catVarList->lData[idx]))->IsUncorrelated()) {
            ((_CategoryVariable*)LocateVar(catVarList->lData[idx]))->SerializeCategory (rec);
        } else {
            nonInd << idx;
        }
    {
        for (long idx = 0; idx < nonInd.lLength; idx++) {
            ((_CategoryVariable*)LocateVar(catVarList->lData[nonInd.lData[idx]]))->SerializeCategory (rec);
        }
    }
}

//__________________________________________________________________________________

void SplitVariablesIntoClasses (_SimpleList& all, _SimpleList& i, _SimpleList& d, _SimpleList& c)
{
    for (long idx = 0; idx < all.lLength; idx++) {
        _Variable* thisVar = LocateVar (all.lData[idx]);
        if (thisVar->IsCategory()) {
            c << all.lData[idx];
        } else if (thisVar->IsIndependent()) {
            i << all.lData[idx];
        } else {
            d << all.lData[idx];
        }
    }
}
