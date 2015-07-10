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

#include "parser.h"
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

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


#define  GOLDEN_RATIO 1.618034
#define MAX_BRENT_ITERATES 100L

extern   _SimpleList BinOps,
         opPrecedence,
         FunctionArgumentCount,
         associativeOps;

extern   long           expressionsParsed;

extern    _List         FunctionNameList;

extern    _AVLListX     _HY_GetStringGlobalTypes;

extern    _Parameter    tolerance,
          sqrtPi,
          maxRombergSteps,
          integrationPrecisionFactor,
          machineEps;

extern    _String       intPrecFact ,
          intMaxIter;

_Parameter  verbosityLevel = 0.0;

extern _Parameter twoOverSqrtPi;

_SimpleList simpleOperationCodes,
            simpleOperationFunctions;


long        subNumericValues = 0;

//__________________________________________________________________________________
_Parameter  AddNumbers  (_Parameter x, _Parameter y)
{
    return x+y;
}
_Parameter  SubNumbers  (_Parameter x, _Parameter y)
{
    return x-y;
}
_Parameter  MinusNumber (_Parameter x)
{
    return -x;
}
_Parameter  MultNumbers (_Parameter x, _Parameter y)
{
    return x*y;
}
_Parameter  DivNumbers  (_Parameter x, _Parameter y)
{
    return x/y;
}
_Parameter  LessThan    (_Parameter x, _Parameter y)
{
    return x<y;
}
_Parameter  GreaterThan (_Parameter x, _Parameter y)
{
    return x>y;
}
_Parameter  LessThanE   (_Parameter x, _Parameter y)
{
    return x<=y;
}
_Parameter  GreaterThanE(_Parameter x, _Parameter y)
{
    return x>=y;
}
_Parameter  Power       (_Parameter x, _Parameter y)
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
_Parameter  MaxNumbers  (_Parameter x, _Parameter y)
{
    return x<y?y:x;
}
_Parameter  MinNumbers  (_Parameter x, _Parameter y)
{
    return x<y?x:y;
}
_Parameter  ExpNumbers  (_Parameter x)
{
    return exp(x);
}
_Parameter  LogNumbers  (_Parameter x)
{
    return log(x);
}
_Parameter  FastMxAccess(Ptr m, _Parameter index)
{
    return ((_Parameter*)m)[(long)index];
}
_Parameter  AndNumbers  (_Parameter x, _Parameter y)
{
    return x != 0.0 && y != 0.0;
}
_Parameter  AbsNumber  (_Parameter x)
{
    return fabs (x);
}

//__________________________________________________________________________________

_Parameter  RandomNumber(_Parameter l, _Parameter u)
{
    _Parameter r = l;
    if (u>l) {
        r=genrand_int32();
        r/=RAND_MAX_32;
        r =l+(u-l)*r;
    }
    return r;
}



//_______________________________________________________________________________________
_Parameter  EqualNumbers(_Parameter a, _Parameter b)
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
            values[k2].reference = (Ptr)((_Matrix*)varValue)->theData;
        }
    }
}


//__________________________________________________________________________________

void        WarnNotDefined (_PMathObj p, long opCode, _hyExecutionContext* context)
{
    _FString * t = (_FString*)p->Type();
    context->ReportError  (_String("Operation '")&*(_String*)BuiltInFunctions(opCode)&"' is not implemented/defined for a " & *t->theString);        
    DeleteObject (t);
}


//__________________________________________________________________________________
_Parameter  InterpolateValue (_Parameter* theX, _Parameter* theY, long n, _Parameter *c , _Parameter *d, _Parameter x, _Parameter& err)
{
    _Parameter y,
               den,
               dif = 1e10,
               dift,
               ho,
               hp,
               w;

    long   ns;

    for (long i=0; i<n; i++) {
        dift = fabs(x-theX[i]);
        if (dift<dif) {
            ns = i;
            dif = dift;
        }
        c[i] = d[i] = theY[i];
    }

    y = theY[ns];
    ns --;

    for (long m=1; m<n; m++) {
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
_Parameter  TrapezoidLevelKSimple (_Formula&f, _Variable* xvar, _Parameter left, _Parameter right, long k, _SimpleFormulaDatum * stack, _SimpleFormulaDatum* values, _SimpleList& changingVars, _SimpleList& varToStack)
{
    _Parameter x,
               tnm,
               sum,
               del,
               ddel;

    static _Parameter s;

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
_Parameter  TrapezoidLevelK (_Formula&f, _Variable* xvar, _Parameter left, _Parameter right, long k)
{
    _Parameter x,
               tnm,
               sum,
               del,
               ddel;

    static _Parameter s;

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
    if (assignment_type != HY_STRING_DIRECT_REFERENCE && reference >= 0) {
        long dereferenced = DereferenceVariable(reference, nameSpace, assignment_type);
        if (dereferenced < 0) {
            WarnError (_String ("Failed to dereference '") & *FetchVar(reference)->GetName() & "' in the " & ((assignment_type == HY_STRING_GLOBAL_DEREFERENCE) ? "global" : "local") & " context");
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
            _Variable * theV = f->Dereference(assignment_type == HY_STRING_GLOBAL_DEREFERENCE, &localContext);
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
        _Variable * theV = f->Dereference(assignment_type == HY_STRING_GLOBAL_DEREFERENCE, &localContext);
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
            WarnError ("Empty RHS in a constraint assignment.");
            return 0;
        }

        _PMathObj varObj = f2->Compute(0, nameSpace);
        if (varObj->ObjectClass()!=NUMBER) {
            WarnError ("Not a numeric RHS in a constraint assignment.");
            return 0;
        }

        _Variable * theV;
        
        if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT || code == HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT) {
            theV= LocateVar (reference);
        } else {
            _hyExecutionContext localContext (nameSpace);
            theV = f->Dereference(assignment_type == HY_STRING_GLOBAL_DEREFERENCE, &localContext);
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
            WarnError ("Empty RHS in an assignment.");
            return 0;
        }

        if (code == HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT) {
            newF.DuplicateReference(f2);
        } else {
            newF.theFormula.AppendNewInstance(new _Operation((_PMathObj)f2->Compute(0, nameSpace)->makeDynamic()));
        }

        long stackD = -1,
             last0  = 0;

        for (long opID = 0; opID < f->theFormula.lLength - 1; opID ++) {
            ((_Operation*)f->theFormula(opID)) -> StackDepth (stackD);
            if (stackD == 0) {
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
            _Variable* mmo = LocateVar(((_Operation*)f->theFormula(0))->GetAVariable());

            if (mmo)
                if (mmo->ObjectClass () == MATRIX) {
                    mmx = (_Matrix*)(mmo->GetValue());
                    ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                } else if (mmo->ObjectClass () == ASSOCIATIVE_LIST) {
                    mma = (_AssociativeList*)(mmo->GetValue());
                    ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                }
        }

        _PMathObj coordMx = nil;
        if (mma || mmx) {
            long expectedType = mmx?MATRIX:STRING;
            coordMx = f->Compute(last0);
            if (!coordMx || coordMx->ObjectClass() != expectedType) {
                if (mmx) {
                    WarnError (_String("Matrix expected but not supplied."));
                } else {
                    WarnError (_String("String key expected but not supplied."));
                }

                return 0;
            }
        } else {
            WarnError ("Matrix/List LHS expected but not supplied.");
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
                        warnError (_String("Can't assign non-polynomial entries to polynomial matrices."));
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

//__________________________________________________________________________________
long        HandleFormulaParsingError (_String errMsg, _String* saveError, _String& s, long index) {
    if (index >= 0) {
        errMsg = errMsg & " in the following context: '"&s.Cut(MAX(0,index-CONTEXT_TRUNCATION),index)&"<ERROR HERE>"&s.Cut(index+1,MIN (index+CONTEXT_TRUNCATION, s.sLength-1)) & "'";
    }
    if (saveError) {
        *saveError = errMsg;
    } else {
        WarnError(errMsg);
    }
    return HY_FORMULA_FAILED;
}

//__________________________________________________________________________________
bool        checkLHS (_List* levelOps, _List* levelData, _String& errMsg, char & deref, _Formula * f, _Variable*& lhs) {
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
                  
    deref = HY_STRING_DIRECT_REFERENCE;
    if (levelOps->lLength > 0) { // this is where 'f is non-empty' cases will go
        check = false;
        if (levelOps->lLength == 1) {
            char buffered_op = ((_Operation*)((*levelOps)(0)))->TheCode(); 
            if (buffered_op == HY_OP_CODE_MUL) {
                check = true; deref = HY_STRING_LOCAL_DEREFERENCE;
            } else {
                if (buffered_op == HY_OP_CODE_POWER) {
                    check = true; deref = HY_STRING_GLOBAL_DEREFERENCE;
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
        _Operation * theOp = (_Operation*)(*levelData)(0);
        if (!theOp->IsAVariable(false)) {
            errMsg = "The left-hand side of an assignment must be a variable (not a constant)";
            return false;
        }        
        lhs = LocateVar(theOp->GetAVariable());
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
        
    if (s.getChar(i)=='>') {
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
    
    if (s.getChar(i-1) != ':') {
        _PMathObj varObj = newF.Compute();
        if (!varObj) {
            return HandleFormulaParsingError ("Invalid RHS in an assignment ", parsingContext.errMsg(), s, i);
        }
        if (twoToken && s.getChar(i-1) == '+') {
            theV->SetValue(theV->Compute()->Execute(HY_OP_CODE_ADD,varObj));
        } else {
            theV->SetValue(varObj);
        }
    } else {
        theV->SetFormula (newF);
    }
    return HY_FORMULA_EXPRESSION;
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

 z op x/y        |   HY_STRING_DIRECT_REFERENCE
 *z op x/y       |   HY_STRING_LOCAL_DEREFERENCE
 ^z op x/y       |   HY_STRING_GLOBAL_DEREFERENCE


*/

{
    static bool inAssignment = false;

    expressionsParsed++;

    _List           operations,
                    operands,
                    *levelOps,
                    *levelData;

    /*04252006*/

    _SimpleList     squareBrackets,
                    mergeMAccess,
                    mergeMAccessLevel;
                    


    long            level                 = 0;
    /* 04252006 mlevel = -1, */
    /* mcount = 0 ; */

    _String curOp;

    bool            impliedMult = false,
                    globalKey   = false,
                    twoToken    = false;

    char            storage     = 0;


    operations.AppendNewInstance (new _List);
    operands.AppendNewInstance (new _List);

    levelOps  = (_List*)(operations(level));
    levelData = (_List*)(operands(level));

    for (long i = 0; i<=s.sLength; i++) {
        storage = 0; // no implied ops by default

        if (isspace(s.getChar(i))) { // skip spaces and tabs
            continue;
        }

        char     lookAtMe = s.getChar(i);

        if (i==s.sLength || lookAtMe == ')' || lookAtMe == ']' || lookAtMe == ',') // closing ) or ]
            // or a parameter list
        {

            /* 04252006 if (level == mlevel && s.getChar(i)!=']')*/
            if (squareBrackets.lLength && squareBrackets.lData[squareBrackets.lLength-1] == level && lookAtMe != ']') {
                return HandleFormulaParsingError ("Missing or unbalanced '[]' ", parsingContext.errMsg(), s, i);
             }

            /* 04252006 if (s.getChar(i)==']' && s.getChar(i+1)!='[')
                mlevel = -1; */

            if (lookAtMe != ',')
                if (i != s.sLength) {
                    level--;
                } else if (level) {
                    level = -1;
                }

            /* 04252206 if (i!=s.sLength && lookAtMe != ',')
                level --;
            else
                if (lookAtMe !=',' && level)
                    level = -1; */

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

                for (int i = 0; i<levelData->countitems(); i++) {
                    f->theFormula << (*levelData)(i);    // mod 07072006 to not duplicate
                }

                levelData->Clear();

                for (int k = levelOps->countitems()-1; k>=0; k--) {
                    f->theFormula << (*levelOps)(k);    // mod 07072006 to not duplicate
                }

                levelOps->Clear();
            } else {
                if (levelData->lLength>1) {
                    return HandleFormulaParsingError ("Syntax error ", parsingContext.errMsg(), s, i);
                } else if (levelData->lLength) {
                    f->theFormula << (*levelData)(0);    // mod 07072006 to not duplicate
                }

                levelData->Clear();
            }

            if (i<s.sLength && lookAtMe !=',' ) {
                operations.Delete(level+1);
                operands.Delete(level+1);
                levelOps    = (_List*)(operations(level));
                levelData   = (_List*)(operands(level));

                if (lookAtMe !=']')
                    if ( BinOps.Find(s.getChar(i+1))==-1 && i<s.sLength-1 && s.getChar(i+1)!=')' && s.getChar(i+1)!=']' && s.getChar(i+1)!='[' && HalfOps.Find(s.getChar(i+1))==-1 && s.getChar(i+1)!=',') {
                        storage = s.getChar(i);
                        s.setChar(i,'*');
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
                        previousMaccess->nInstances++;
                        f->theFormula.Delete (mergeIndex);
                        f->theFormula << previousMaccess;
                        previousMaccess->nInstances--;
                    }
                } else {
                    f->theFormula.AppendNewInstance(new _Operation (curOp ,2));
                }
            }

            if (!storage) {
                continue;
            }
        }


        if (s.getChar(i) == '=' && s.getChar(i+1) != '=' && (!twoToken || s.getChar(i-1)==':' || s.getChar (i-1) == '+')) { // assignment operator
            _String  errMsg;

            bool check               = !inAssignment,
                 is_array_assignment = f->IsArrayAccess();
                 
            char deref = 0; 
            
            _Variable *lhs_variable = nil;

 
            if (check) {
                if (is_array_assignment) {
                    (((_Operation*)((f->theFormula)(f->theFormula.lLength-1)))->TheCode()) = HY_OP_CODE_MCOORD;
                } else {
                    check = checkLHS (levelOps, levelData, errMsg, deref, f, lhs_variable);
                }
            } else {
                errMsg = "Can't assign within another assignment";
            }
            
            if (!check) {
                return HandleFormulaParsingError (errMsg, parsingContext.errMsg(), s, i);
            }

            inAssignment = true;
            _String ss (s,i+1,-1); // this is the RHS
            _Formula  newF;
           
            if (Parse(&newF,ss,parsingContext, f2) != HY_FORMULA_EXPRESSION) {
                inAssignment = false;
                return HY_FORMULA_FAILED;
            }
            inAssignment = false;
            if (!is_array_assignment && lhs_variable)
                // normal variable assignment
            {
                if (!f2) { // immediate execution
                    if (_parserHelperHandleInlineAssignmentCases (s,parsingContext, i, lhs_variable, f,  deref, newF, twoToken) == HY_FORMULA_FAILED) {
                        return HY_FORMULA_FAILED;
                    }
                } else { // this gets called from ExecuteCase0...
                    if (twoToken && s.getChar(i-1) == '+') { // += gets handled here
                    
                        _Operation* self = new _Operation ();
                        self->SetAVariable(lhs_variable->GetAVariable());
                        newF.theFormula.InsertElement (self,0,false);
                        DeleteObject (self);
                        if (deref != HY_STRING_DIRECT_REFERENCE) {
                             _Operation* ref = new _Operation (*(_String*)BuiltInFunctions(deref == HY_STRING_GLOBAL_DEREFERENCE ? HY_OP_CODE_POWER : HY_OP_CODE_MUL),1);
                             newF.theFormula.InsertElement (ref,1,false);
                             DeleteObject (ref);
                      }
                      newF.theFormula.AppendNewInstance (new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                    }
                    f->Duplicate((BaseRef)&newF);
                }
                twoToken     = false;

                parsingContext.assignmentRefID()   = lhs_variable->GetAVariable();
                parsingContext.assignmentRefType() = deref;

                return (s.getChar(i-1)==':')?HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT:HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT;
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

                        if (newF.IsAConstant() || s.getChar(i-1) !=':') {
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

                            if (mmo)
                                if (mmo->ObjectClass () == MATRIX) {
                                    mmx = (_Matrix*)(mmo->GetValue());
                                    ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                                } else if (mmo->ObjectClass () == ASSOCIATIVE_LIST) {
                                    mma = (_AssociativeList*)(mmo->GetValue());
                                    ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
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

                                mmx->MStore (&hC, &vC, newF, (twoToken && s.getChar(i-1) =='+')?HY_OP_CODE_ADD:HY_OP_CODE_NONE);
                            }
                        } else if (mma) {
                            _PMathObj coordIdx = f->Compute(last0);

                            if (!coordIdx|| coordIdx->ObjectClass() != STRING ) {
                                anError = true;
                            } else {
                                mma->MStore (coordIdx, newF.Compute(),true, (twoToken && s.getChar(i-1) =='+')?HY_OP_CODE_ADD:HY_OP_CODE_NONE);
                            }
                        } else {
                            anError = true;
                        }


                        if (anError) {
                            return HandleFormulaParsingError ("Invalid matrix/associative list ident supplied ", parsingContext.errMsg(), s, i);
                        }

                        return HY_FORMULA_EXPRESSION;
                    } else {
                        bool isSimple = (s.getChar(i-1) != ':');
                        f2->Duplicate   ((BaseRef)&newF);
                        if (last0 == 0) {
                            ((_Operation*)f->theFormula(0))->SetAVariable(-((_Operation*)f->theFormula(0))->GetAVariable()-3);
                        }
                        return isSimple?((s.getChar(i-1) == '+')?HY_FORMULA_FORMULA_VALUE_INCREMENT:HY_FORMULA_FORMULA_VALUE_ASSIGNMENT):HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT;
                    }
                } else {
                // *(expression) reference
                    if (f2) {
                        bool isSimple = (s.getChar(i-1) != ':');

                        if (twoToken && s.getChar(i-1) == '+') { // += gets handled here
                            newF.theFormula.InsertElement (new _Operation (*(_String*)BuiltInFunctions(deref == HY_STRING_GLOBAL_DEREFERENCE ? HY_OP_CODE_POWER : HY_OP_CODE_MUL),1), 0, false);
                            for (long lhs_ops = f->theFormula.lLength-1; lhs_ops >= 0; lhs_ops --) {  
                                newF.theFormula.InsertElement (f->theFormula(lhs_ops), 0, true);
                            }
                            newF.theFormula.AppendNewInstance (new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                        }
                        f2->Duplicate   ((BaseRef)&newF);
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

        if ( s.getChar(i-1)==':' && (s.getChar(i)=='<' || s.getChar(i)=='>')) { // variable bounds
            _Variable * lhs = nil;
            _String errMsg;
            char    deref;
            
            if (inAssignment||f->IsArrayAccess()||! checkLHS (levelOps, levelData, errMsg, deref, f, lhs)) {
               return HandleFormulaParsingError ("Can't set bounds like this ", parsingContext.errMsg(), s, i);
            }

            inAssignment = true;

            _String ss (s,i+1,-1);
            _Formula newF;

            if (Parse(&newF,ss,parsingContext,f2) != HY_FORMULA_EXPRESSION) {
                inAssignment = false;
                return HY_FORMULA_FAILED;
            }


            inAssignment = false;
            twoToken     = false;


    
            if (lhs) {
                if (!f2) {
                    if (_parserHelperHandleInlineBoundCases (s,parsingContext,i,lhs,f, deref, newF) == HY_FORMULA_FAILED) {
                        return HY_FORMULA_FAILED;
                    } 
                } else { // BOUND ASSIGNMENTS
                    f2->Duplicate   ((BaseRef)&newF);

                    parsingContext.assignmentRefID()   = lhs->GetAVariable();
                    parsingContext.assignmentRefType() = deref;

                    return (s.getChar(i)=='>')?HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT:HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT;
                }
            } else {
                if (!f2) {
                    if (_parserHelperHandleInlineBoundCases (s,parsingContext,i,lhs,f, deref, newF) == HY_FORMULA_FAILED) {
                        return HY_FORMULA_FAILED;
                    } 
                } else { // BOUND ASSIGNMENTS
                    f2->Duplicate   ((BaseRef)&newF);
                    parsingContext.assignmentRefType() = deref;
                    return (s.getChar(i)=='>')?HY_FORMULA_REFERENCE_UPPER_BOUND_ASSIGNMENT:HY_FORMULA_REFERENCE_LOWER_BOUND_ASSIGNMENT;
                }
            
            }

            return HY_FORMULA_EXPRESSION;
        }


        if (s.getChar(i) == '{') // a matrix
            /* 20090803 SLKP:
                 fixed the code to deal with
            */
        {
        
            parsingContext.isVolatile() = true;
            
            int     j       = s.ExtractEnclosedExpression (i,'{','}',true,true);

            if (j<0) {
                return HandleFormulaParsingError ("Poorly formed matrix/associative array construct ", parsingContext.errMsg(), s, i);
            }

            _String matrixDef   (s,i,j);

            if (matrixDef.sLength == 2 || matrixDef.sData[1] == '"') {
                _AssociativeList *theList = new _AssociativeList ();
                if (!theList) {
                    checkPointer (theList);
                }
                if (matrixDef.sLength > 2) {
                    matrixDef.Trim (1,matrixDef.sLength-2);
                    if (!theList->ParseStringRepresentation (matrixDef,parsingContext.errMsg() == nil, parsingContext.formulaScope())) {
                        return HandleFormulaParsingError ("Poorly formed associative array construct ", parsingContext.errMsg(), s, i);
                    }
                }

                levelData->AppendNewInstance (new _Operation (theList));
            } else {
                _Matrix *theMatrix = new _Matrix (matrixDef,false,parsingContext.formulaScope());
                if (!theMatrix) {
                    checkPointer (theMatrix);
                }
                levelData->AppendNewInstance (new _Operation (theMatrix));
            }

            i = j;
            continue;
            
        }

        if (s.getChar(i) == '[') { // opening [
            long  lastCode = -1;

            if (!f->IsEmpty()) {
                lastCode = ((_Operation*)((f->theFormula)(f->theFormula.lLength-1)))->TheCode();
            }

            if (lastCode == HY_OP_CODE_MACCESS && s.getChar(i-1) == ']') {
                mergeMAccess << f->theFormula.lLength-1;
                mergeMAccessLevel << level;
            } else {
                if (levelData->lLength == 0 && f->IsEmpty()) {
                   return HandleFormulaParsingError ("[..] must be preceded by an object to index ", parsingContext.errMsg(), s, i);
                }

                if (levelData->lLength) {
                    f->theFormula.AppendNewInstance((*levelData)[levelData->lLength-1]);
                    levelData->Delete(levelData->lLength-1,false);
                }
            }

            squareBrackets << ++level;

            curOp       = empty;
            operations.AppendNewInstance(new _List);
            operands.AppendNewInstance(new _List);
            levelOps   = (_List*) (operations(level));
            levelData  = (_List*) (operands  (level));

            continue;
        }


        if (s.getChar(i) == '(') { // opening (
            level++;
            operations.AppendNewInstance (new _List);
            operands.AppendNewInstance   (new _List);
            levelOps    =   (_List*)(operations(level));
            levelData   =   (_List*)(operands(level));
            curOp       =   empty;
            continue;
        }

        if (s.getChar(i)=='"') { // a string literal
            long j             = 1,
                 inPlaceID     = -1;

            _String * literal = (_String*)checkPointer(new _String (16,true));

            while (i+j<s.sLength) {
                char aChar = s.sData[i+j];
                if (aChar =='\\') {
                    if (i+j+1<s.sLength) {
                        if (s.sData[i+j+1]=='"' ||s.sData[i+j+1]=='`' ) {
                            j++;
                            (*literal)<<s.sData[i+j++];
                        } else {
                            (*literal)<<s.sData[i+j++];
                            (*literal)<<s.sData[i+j++];
                        }
                    }
                    continue;
                }

                if (aChar =='"' && inPlaceID < 0) {
                    break;
                }

                if (aChar == '`') {
                    if (inPlaceID < 0) {
                        inPlaceID = ++j;
                    } else if (j == inPlaceID) {
                        return HandleFormulaParsingError ("Attempted to string substitute an empty quotation ", parsingContext.errMsg(), s, i);
                    } else {
                        _String     inPlaceVID (s,i+inPlaceID,i+j-1),
                                    inPlaceValue = ProcessLiteralArgument(&inPlaceVID, parsingContext.formulaScope());

                        /*if (!inPlaceValue) {
                            inPlaceValue = (_FString*)ProcessLiteralArgument(&inPlaceVID, theParent);
                            if (!inPlaceValue) {
                                return HandleFormulaParsingError ("Attempted to string substitute something other that a string variable/expression ", saveError, s, i);
                            }
                        }*/

                        (*literal) << inPlaceValue;
                        inPlaceID = -1;
                        parsingContext.isVolatile() = true;
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
            if (inPlaceID >= 0) {
                return HandleFormulaParsingError ("Unterminated string substitution inside a literal ", parsingContext.errMsg(), s, i);
            }
            levelData->AppendNewInstance (new _Operation (new _FString(*literal)));
            DeleteObject(literal);

            i += j;
            continue;

        }

        if (alpha.isAllowed [(unsigned char)s.getChar(i)]) { // an identifier
            bool takeVarReference = false;
            
            if (twoToken) {
                char opChar = s.getChar(i-1);
                if (((_String*)BuiltInFunctions(HY_OP_CODE_REF))->Equal(opChar)) {
                    takeVarReference = true;
                } else {
                    _String thisOp (opChar);
                    levelOps->AppendNewInstance (new _Operation (thisOp,1L));
                }
            }
            
            impliedMult = (i && numeric.isAllowed [(unsigned char)s.getChar(i-1)]);

            long j = 1;
            while ( i+j<s.sLength && (alpha.isAllowed [(unsigned char)s.getChar(i+j)]|| numeric.isAllowed [(unsigned char)s.getChar(i+j)]) ) {
                j++;
            }

            curOp =  (s.Cut(i,i+j-1));
            i+=j-1;

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
                
            if (UnOps.Find(curOp)>=0) { // a standard function
                if (takeVarReference) {
                    return HandleFormulaParsingError ("Cannot make a reference from a built-in function", parsingContext.errMsg(), s, i);
                }
                                
                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else { // a variable
                // check if this is a function defined  in the list of "standard functions"
                long bLang = noneObject?-1:FunctionNameList.BinaryFind(&curOp);
                if (bLang>=0) {
                    if (takeVarReference) {
                        return HandleFormulaParsingError ("Cannot make a reference from a built-in function", parsingContext.errMsg(), s, i);
                    }
                    levelOps->AppendNewInstance (new _Operation (curOp,FunctionArgumentCount(bLang)));
                    continue;
                }

                // check if this is a function defined in the batch language

                if ((bLang =  noneObject?-1:FindBFFunctionName (curOp, parsingContext.formulaScope()))>=0) {
                    if (takeVarReference) {
                        return HandleFormulaParsingError ("Cannot make a reference from user-defined function", parsingContext.errMsg(), s, i);
                    }
                    levelOps->AppendNewInstance (new _Operation (curOp,-bLang-1));
                    continue;
                }

                long curOpl = curOp.sLength;
                if (curOpl>2 && curOp[curOpl-1]=='_' && curOp[curOpl-2]=='_') { // instant variable refrence
                    _String realVarName (curOp,0,curOpl-3);

                    realVarName = parsingContext.contextualizeRef (realVarName);
                   

                    long realVarLoc = LocateVarByName (realVarName);
                    if (realVarLoc<0) { // bad instant variable reference
                        return HandleFormulaParsingError ("Attempted to take value of undeclared variable ", parsingContext.errMsg(), s, i);
                     }
                    if (!f2) { // 03/25/2004 ? Confused why the else
                        levelData->AppendNewInstance(new _Operation((_MathObject*)FetchVar (realVarLoc)->Compute()->makeDynamic()));
                    } else {
                        _Operation theVar (true, realVarName, globalKey, parsingContext.formulaScope());
                        theVar.SetTerms(-variableNames.GetXtra (realVarLoc)-1);
                        theVar.SetAVariable(-2);
                        (*levelData) && (&theVar);
                    }
                } else {
                    if (noneObject)
                        levelData->AppendNewInstance (new _Operation (false, curOp));
                    else
                        if (parsingContext.formulaScope() && _hyApplicationGlobals.Find(&curOp) >= 0) {
                            levelData->AppendNewInstance (new _Operation(true, curOp, globalKey, nil, takeVarReference));
                        } else {
                            levelData->AppendNewInstance (new _Operation(true, curOp, globalKey, parsingContext.formulaScope(), takeVarReference));
                        }
                }
                globalKey = false;
                if (impliedMult) {
                    storage = s.getChar(i);
                    s.setChar(i,((_String*)BuiltInFunctions(HY_OP_CODE_MUL))->getChar(0));
                } else if (s.getChar(i+1)=='(') {
                    if (!storage) {
                        storage = s.getChar(i);
                        s.setChar(i,((_String*)BuiltInFunctions(HY_OP_CODE_MUL))->getChar(0));
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

        if (numeric.isAllowed [(unsigned char)s.getChar(i)]) {
            if (twoToken) {
                _String thisOp (s.getChar(i-1));
                levelOps->AppendNewInstance (new _Operation (thisOp,1L));
            }
            long j = 1;

            while ( i+j<s.sLength && (numeric.isAllowed [(unsigned char)s.getChar(i+j)] || ((s.getChar(i+j)=='-' || s.getChar(i+j)=='+' )&& tolower(s.getChar(i+j-1))=='e')) ) {
                j++;
            }

            curOp =  (s.Cut(i,i+j-1));
            i+=j-1;
            levelData->AppendNewInstance (new _Operation (false, curOp));
            if (i<s.sLength-1 && s.getChar(i+1)=='(') {
                storage = s.getChar(i);
                s.setChar(i,((_String*)BuiltInFunctions(HY_OP_CODE_MUL))->getChar(0));
            } else {
                continue;
            }
        }
        
        if ( BinOps.Find (s.getChar(i))!=-1 || (twoToken&& (BinOps.Find(s.getChar(i-1)*(long)256+s.getChar(i))!=-1)) ) {
            if (!twoToken && BinOps.Find(s.getChar(i)*(long)256+s.getChar(i+1)) != -1) {
                twoToken = true;
                continue;
            }

            if (twoToken||(BinOps.Find(s.getChar(i)*256+s.getChar(i+1))!=-1)) {
                if (!twoToken) {
                    i++;
                }
                curOp = _String(s.getChar(i-1))&(_String)(s.getChar(i));
            } else {
                curOp = s.getChar(i);
            }

            long twoOrOne = 2;

            if (storage) {
                s.setChar(i,storage);
            }

            if (levelData->countitems()==0) {
                if (s[i-curOp.sLength]!=')' && storage!=')' && s[i-curOp.sLength] !=']') {
                    if (!twoToken && UnOps.Find (s.getChar(i)) >= 0) {
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
                        f->theFormula&&((*levelData)(k));
                    }

                    levelData->Clear();
                    levelData->AppendNewInstance (newS);
                } else {
                    for (unsigned long k = 0; k<levelData->countitems(); k++) {
                        f->theFormula << ((*levelData)(k));
                    }
                    levelData->Clear();
                }
            }

            if (!levelOps->countitems()) {
                levelOps->AppendNewInstance (new _Operation (curOp,twoOrOne));
                if (terminateExecution) {
                    return HY_FORMULA_FAILED;
                }
                continue;
            }

            // check operation precedence

            long h,g;
            _String prevOp = *((((_Operation*)((*levelOps)(levelOps->countitems()-1)))->GetCode()));

            h = BinOps.Find(prevOp.sLength==2?prevOp.getChar(0)*256+prevOp.getChar(1):prevOp.getChar(0));
            g = BinOps.Find(curOp.sLength==2?curOp.getChar(0)*256+curOp.getChar(1):curOp.getChar(0));

            if (h!=-1) {
                h = opPrecedence (h);
            }

            g = opPrecedence (g);



            if (g>h && h!=-1) { // store the op, don't do it yet!
                levelOps->AppendNewInstance (new _Operation (curOp,twoOrOne));
                if (terminateExecution) {
                    return HY_FORMULA_FAILED;
                }
                continue;
            }

            // do the stored operations now

            for (int j = levelOps->countitems()-1; j>=0; j--) {
                _String  sss = (((_Operation*)((*levelOps)(levelOps->countitems()-1)))->GetCode());
                h = BinOps.Find(sss.sLength==2?sss.getChar(0)*256+sss.getChar(1):sss.getChar(0));
                if (h==-1) {
                    h=100;
                } else {
                    h = opPrecedence (h);
                }

                if (h<g) {
                    break;
                }
                f->theFormula&&((*levelOps)(j));
                levelOps->Delete((*levelOps).lLength-1);
            }
            levelOps->AppendNewInstance (new _Operation (curOp,twoOrOne));
            if (terminateExecution) {
                return HY_FORMULA_FAILED;
            }
            continue;
        } else if (UnOps.Find(s.getChar(i)) >= 0) {
            if ((s.getChar(i)=='-' || s.getChar(i)=='+') && (!i|| s.getChar(i-1)=='(')) { // unary minus?
                curOp   = s.getChar(i);
                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else {
                if (HalfOps.contains(s.getChar(i))) {
                    twoToken = true;
                    continue;
                }
                return HandleFormulaParsingError ("Bad binary operator placement ", parsingContext.errMsg(), s, i);
            }
        } else {
            if (!HalfOps.contains(s.getChar(i))) {
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
    checkParameter (VerbosityLevelString, verbosityLevel, -1.0);
    return verbosityLevel;
}

//__________________________________________________________________________________
void  checkParameter (_String& name, _Parameter& dest, _Parameter def, _VariableContainer* theP)
{
    long f;
    if (theP) {
        _String ppn = *theP->GetName() & '.' & name;
        f = LocateVarByName(ppn);
    } else {
        f = LocateVarByName (name);
    }
    if (f<0) {
        dest = def;
    } else {
        dest = FetchVar(f)->Value();
    }
}

//__________________________________________________________________________________
void  stashParameter (_String& name, _Parameter v, bool set)
{
    static  _Parameter stash = 0.0;

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
void  setParameter (_String& name, _Parameter def, _String* namespc)
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

void  setParameter (_String& name, _PMathObj def, bool dup, _String* namespc)
{
    if (namespc) {
        _String namespcd = AppendContainerName(name,namespc);
        setParameter (namespcd,def,dup);
    } else {
        long f = LocateVarByName (name);
        if (f<0) {
            _Variable cornholio(name);
            setParameter (name,def,dup);
        } else {
            FetchVar(f)->SetValue(def,dup);
        }
    }
}

//__________________________________________________________________________________

void ExportIndVariables (_String& glVars, _String& locVars, _SimpleList* indepVarList)
{
    _String * stIn,
              str;
             
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

void ExportDepVariables (_String& glVars, _String& locVars, _SimpleList* depVarList)
{
    if (depVarList->lLength) {
        _String * stIn,
                  str;

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

void ExportCatVariables (_String& rec, _SimpleList* catVarList)
{
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
