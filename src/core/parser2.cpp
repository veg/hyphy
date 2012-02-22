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
    for (long k2 = 0; k2 < vars.lLength; k2++) {
        _PMathObj varValue = LocateVar (vars.lData[k2])->Compute();
        if (varValue->ObjectClass() == NUMBER) {
            values[k2].value = varValue->Value();
        } else {
            values[k2].reference = (Ptr)((_Matrix*)varValue)->theData;
        }
    }
}


//__________________________________________________________________________________

void        WarnNotDefined (_PMathObj p, long opCode)
{
    _FString * t = (_FString*)p->Type();
    WarnError (_String("Operation '")&*(_String*)BuiltInFunctions(opCode)&"' is not implemented/defined for a " & *t->theString);
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

long       ExecuteFormula (_Formula*f , _Formula* f2, long code, long reference, _VariableContainer* nameSpace)
{
    if (code == HY_FORMULA_EXPRESSION || code == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
        _PMathObj  formulaValue = f->Compute(0, nameSpace);
        if (!formulaValue) {
            return 0;
        }

        if (code == HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT) {
            LocateVar (reference)->SetValue (formulaValue);
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

    if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT || code == HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT) {
        if (f2->IsEmpty()) {
            WarnError ("Empty RHS in a constraint assignment.");
            return 0;
        }

        _PMathObj varObj = f2->Compute();
        if (varObj->ObjectClass()!=NUMBER) {
            WarnError ("Not a numeric RHS in a constraint assignment.");
            return 0;
        }

        _Variable * theV = LocateVar (reference);

        if (code == HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT ) {
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

//__________________________________________________________________________________
long        Parse (_Formula* f, _String& s, long& variableReference, _VariableContainer* theParent, _Formula* f2, bool flagErrors, bool* isVolatile)
/* SLKP 20110908: added the concept of a 'volatile' formula, i.e. something that should be reparsed every time in ExecuteCase0
                : currently those include
                :    inline constructors (matrices, dictionaries)
                :    `` substitutions in strings
                
    
   SLKP 20100817: decoupled return code from variable reference return
*/

// returns:

/*

 case                       | return value                                  | variableReference value

 parse failed               | HY_FORMULA_FAILED                             | undefined
 expresion (no LHS)         | HY_FORMULA_EXPRESSION                         | undefined
 z = x/y                    | HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT          | index of the LHS
 z := x/y                   | HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT        | index of the LHS
 object[a] := x/y           | HY_FORMULA_FORMULA_FORMULA_ASSIGNMENT         | undefined
 object[a] = x/y            | HY_FORMULA_FORMULA_VALUE_ASSIGNMENT           | undefined
 z :< expr                  | HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT    | index of the LHS
 z :> expr                  | HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT    | index of the LHS


 -1     - valid formula - no variable assignment, e.g. x + y
 -3     - regular '=' with formula receptacle, e.g. matrix [a] = 1
 -4     = ":=" with formula receptacle, e.g. matrix [1] := x + y
 <-4        = ":=" with variable receptacle, e.g. z := x / y
 >=0        "=" with variable receptacle, e.g. z = x/y

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
                if (flagErrors) {
                    WarnError (_String("Missing or unbalanced '[]' in Expression:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return -2;
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
                if (flagErrors) {
                    WarnError (_String("Unbalanced Parentheses in Expression:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return HY_FORMULA_FAILED;
            }

            if (lookAtMe ==',' && (level<1 || (squareBrackets.lLength && squareBrackets.lData[squareBrackets.lLength-1] == level))) {
                if (flagErrors) {
                    WarnError (_String("Parameter list is out of context:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return HY_FORMULA_FAILED;
            }

            if (levelOps->lLength) { // there are some buffered operations left
                if (levelOps->lLength > 3 || levelData->lLength > 2) {
                    if (flagErrors) {
                        WarnError (_String ("Syntax error in expression:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                    }
                    return HY_FORMULA_FAILED;
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
                    if (flagErrors) {
                        WarnError (_String ("Syntax error in expression:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                    }
                    return HY_FORMULA_FAILED;
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
                    if (flagErrors) {
                        WarnError (_String ("Unexpected ']': ")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                    }
                    return HY_FORMULA_FAILED;
                }
                squareBrackets.Delete(squareBrackets.lLength-1);
                curOp = *(_String*)BuiltInFunctions(HY_OP_CODE_MACCESS);
                if (mergeMAccess.lLength && mergeMAccess.lData[mergeMAccess.lLength-1] >= 0 && mergeMAccessLevel.lData[mergeMAccessLevel.lLength-1] == level) {
                    long mergeIndex              = mergeMAccess.lData[mergeMAccess.lLength-1];
                    _Operation * previousMaccess = (_Operation*) f->theFormula (mergeIndex);
                    if (previousMaccess->GetCode () != curOp) {
                        if (flagErrors) {
                            WarnError (_String ("Internal error in Parse. Incorrect matrix access token code: ")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                        }
                        return HY_FORMULA_FAILED;
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
                    _Operation theVar (curOp ,2);
                    f->theFormula&&(&theVar);
                }
            }

            if (!storage) {
                continue;
            }
        }


        if (s.getChar(i) == '=' && s.getChar(i+1) != '=' && (!twoToken || s.getChar(i-1)==':' || s.getChar (i-1) == '+')) { // assignment operator
            _String*  sss = nil;

            if (f->IsEmpty() == false) { // have buffered operations
                sss = &(((_Operation*)((f->theFormula)(f->theFormula.lLength-1)))->GetCode());
            }

            bool check = !inAssignment;

            if (check) {
                if (sss) {
                    if (!sss->Equal((_String*)BuiltInFunctions(HY_OP_CODE_MACCESS))) {
                        check = false;
                    } else {
                        (((_Operation*)((f->theFormula)(f->theFormula.lLength-1)))->TheCode()) = HY_OP_CODE_MCOORD;
                    }
                    /* this will break if another operation is inserted between MAccess and MCoord */
                    //*sss = "MCoord";
                } else if (!f->IsEmpty() ||
                           levelData->countitems()!=1 ||
                           !((_Operation*)(*levelData)(0))->IsAVariable()) {
                    check = false;
                }
            }
            if (!check) {
                if (flagErrors) {
                    WarnError (_String("Can't assign like this:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return HY_FORMULA_FAILED;
            }

            inAssignment = true;
            _String ss (s,i+1,-1);
            _Formula  newF;
            long      refV;
            if (Parse(&newF,ss,refV,theParent,f2,flagErrors,isVolatile) != HY_FORMULA_EXPRESSION) {
                inAssignment = false;
                return HY_FORMULA_FAILED;
            }
            inAssignment = false;
            if (!sss)
                // normal variable assignment
            {
                _Variable * theV = (_Variable*)LocateVar((((_Operation*)(*levelData)(0))->GetAVariable()));
                if (!f2) {
                    if (s.getChar(i-1) != ':') {
                        _PMathObj varObj = newF.Compute();
                        if (!varObj) {
                            if (flagErrors) {
                                WarnError (_String("Invalid RHS in an assignment:")&s.Cut(i+1,-1));
                            }
                            return HY_FORMULA_FAILED;
                        }
                        if (twoToken && s.getChar(i-1) == '+') {
                            theV->SetValue(theV->Compute()->Execute(HY_OP_CODE_ADD,varObj));
                        } else {
                            theV->SetValue(varObj);
                        }
                    } else {
                        theV->SetFormula (newF);
                    }
                } else { // this gets called from ExecuteCase0...
                    if (twoToken && s.getChar(i-1) == '+') {
                        _Operation* self = new _Operation ();
                        self->SetAVariable(theV->GetAVariable());
                        newF.theFormula.InsertElement (self,0,false);
                        DeleteObject (self);
                        newF.theFormula.AppendNewInstance (new _Operation (*(_String*)BuiltInFunctions(HY_OP_CODE_ADD),2));
                    }
                    f->Duplicate((BaseRef)&newF);
                }
                twoToken     = false;

                variableReference = theV->GetAVariable();

                return (s.getChar(i-1)==':')?HY_FORMULA_VARIABLE_FORMULA_ASSIGNMENT:HY_FORMULA_VARIABLE_VALUE_ASSIGNMENT;
            } else
                // matrix/associative array element assignment
            {
                long stackD = -1,
                     last0  = 0;

                for (long opID = 0; opID < f->theFormula.lLength - 1; opID ++) {
                    ((_Operation*)f->theFormula(opID)) -> StackDepth (stackD);
                    if (stackD == 0) {
                        last0 = opID;
                    }
                }

                if (!f2) {
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
                        if (flagErrors) {
                            WarnError (_String("Invalid Matrix/Associative List Ident Supplied:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                        }
                        return HY_FORMULA_FAILED;
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
            }
        }

        if ( s.getChar(i-1)==':' && (s.getChar(i)=='<' || s.getChar(i)=='>')) { // variable bounds
            if (inAssignment||(!f->IsEmpty())||(levelData->countitems()!=1)||!(((_Operation*)(*levelData)(0))->IsAVariable())) {
                if (flagErrors) {
                    WarnError (_String("Can't set bounds like this: ")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return HY_FORMULA_FAILED;
            }

            inAssignment = true;

            _String ss (s,i+1,-1);
            _Formula newF;

            long     refV;

            if (Parse(&newF,ss,refV,theParent,f2,flagErrors,isVolatile) != HY_FORMULA_EXPRESSION) {
                inAssignment = false;
                return HY_FORMULA_FAILED;
            }


            inAssignment = false;
            twoToken = false;
            _Variable * theV = (_Variable*)LocateVar((((_Operation*)(*levelData)(0))->GetAVariable()));

            if (!f2) {
                _PMathObj varObj = newF.Compute();
                if (varObj->ObjectClass()!=NUMBER) {
                    if (flagErrors) {
                        WarnError (_String("Variable bound must evaluate to a number: ")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                    }
                    return HY_FORMULA_FAILED;
                }

                if (s.getChar(i)=='>') {
                    theV->SetBounds(varObj->Value(),theV->GetUpperBound());
                } else {
                    theV->SetBounds(theV->GetLowerBound(),varObj->Value());
                }
            } else {
                f2->Duplicate   ((BaseRef)&newF);
                variableReference = theV->GetAVariable();
                return (s.getChar(i)=='>')?HY_FORMULA_VARIABLE_UPPER_BOUND_ASSIGNMENT:HY_FORMULA_VARIABLE_LOWER_BOUND_ASSIGNMENT;
            }

            return HY_FORMULA_EXPRESSION;
        }


        if (s.getChar(i) == '{') // a matrix
            /* 20090803 SLKP:
                 fixed the code to deal with
            */
        {
        
            if (isVolatile) *isVolatile = true;
            
            if (flagErrors) {
                int     j       = s.ExtractEnclosedExpression (i,'{','}',true,true);

                if (j<0) {
                    if (flagErrors) {
                        WarnError (_String("Poorly formed matrix/associative array construct:")&s.Cut(i,-1));
                    }
                    return HY_FORMULA_FAILED;
                }

                _String matrixDef   (s,i,j);

                if (matrixDef.sLength == 2 || matrixDef.sData[1] == '"') {
                    _AssociativeList *theList = new _AssociativeList ();
                    if (!theList) {
                        checkPointer (theList);
                    }
                    if (matrixDef.sLength > 2) {
                        matrixDef.Trim (1,matrixDef.sLength-2);
                        if (!theList->ParseStringRepresentation (matrixDef,flagErrors, theParent)) {
                            if (flagErrors) {
                                WarnError (_String("Poorly formed associative array construct:")&s.Cut(i,-1));
                            }
                            return HY_FORMULA_FAILED;
                        }
                    }

                    levelData->AppendNewInstance (new _Operation (theList));
                } else {
                    _Matrix *theMatrix = new _Matrix (matrixDef,false,theParent);
                    if (!theMatrix) {
                        checkPointer (theMatrix);
                    }
                    levelData->AppendNewInstance (new _Operation (theMatrix));
                }

                i = j;
                continue;
            } else {
                return HY_FORMULA_FAILED;
            }
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
                    if (flagErrors) {
                        WarnError (_String("[..] must be preceded by an object to index:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                    }
                    return HY_FORMULA_FAILED;
                }

                if (levelData->lLength) {
                    f->theFormula&&((*levelData)[levelData->lLength-1]);
                    levelData->Delete(levelData->lLength-1);
                }
            }

            squareBrackets << ++level;

            curOp       = empty;
            _List          blank;
            operations && &blank;
            operands   && &blank;
            levelOps   = (_List*) (operations(level));
            levelData  = (_List*) (operands  (level));

            /*if ( mcount == 0 && (levelData->lLength == 0 || !((_Operation*)(*levelData)(levelData->lLength-1))->IsAVariable()))
                {
                if (flagErrors) WarnError (_String("[..] must be preceded by a matrix variable:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                return -2;
            }
            else
                if (mcount == 0)
                {
                    f->theFormula&&((*levelData)[levelData->lLength-1]);
                    levelData->Delete(levelData->lLength-1);
                }

            mcount++;
            if (mcount > 2)
            {
                if (flagErrors) WarnError (_String("Only single or double indexing of arrays is allowed: ")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                return -2;
            }
            level++;
            curOp       = empty;
            _List          blank;
            operations && &blank;
            operands   && &blank;
            levelOps   = (_List*) (operations(level));
            levelData  = (_List*) (operands  (level));
            mlevel     = level;*/
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

                if (aChar =='"') {
                    break;
                }

                if (aChar == '`') {
                    if (inPlaceID < 0) {
                        inPlaceID = ++j;
                    } else if (j == inPlaceID) {
                        if (flagErrors) {
                            WarnError (_String("Attempted to string substitute an empty quotation:")&s.Cut(0,i+j)&"?"&s.Cut(i+j+1,-1));
                        }
                        return HY_FORMULA_FAILED;
                    } else {
                        _String     inPlaceVID (s,i+inPlaceID,i+j-1);
                        _FString    *inPlaceValue = (_FString*)FetchObjectFromVariableByType (&inPlaceVID, STRING);

                        if (!inPlaceValue) {
                            if (flagErrors) {
                                WarnError (_String("Attempted to string substitute something other that a string variable:")&s.Cut(0,i+j)&"?"&s.Cut(i+j+1,-1));
                            }
                            return HY_FORMULA_FAILED;
                        }

                        (*literal) << inPlaceValue->theString;
                        inPlaceID = -1;
                        if (isVolatile) *isVolatile = true;
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
                if (flagErrors) {
                    WarnError (_String("Unterminated string substitution inside a literal:")&s.Cut(0,i+j)&"?"&s.Cut(i+j+1,-1));
                }
                return HY_FORMULA_FAILED;
            }
            _FString   *fString = new _FString (*literal);
            DeleteObject (literal);
            checkPointer (fString);

            levelData->AppendNewInstance (new _Operation (fString));

            i += j;
            continue;

        }

        if (alpha.isAllowed [s.getChar(i)]) { // an identifier
            if (twoToken) {
                _String thisOp (s.getChar(i-1));
                levelOps->AppendNewInstance (new _Operation (thisOp,1L));
            }
            impliedMult = (i && numeric.isAllowed [s.getChar(i-1)]);

            long j = 1;
            while ( i+j<s.sLength && (alpha.isAllowed [s.getChar(i+j)]|| numeric.isAllowed [s.getChar(i+j)]) ) {
                j++;
            }

            curOp =  (s.Cut(i,i+j-1));
            i+=j-1;

            if (curOp.Equal(&globalToken)) {
                globalKey = true;
                continue;
            }
            
            bool noneObject = false;
            if (curOp.Equal(&noneToken)) {
                noneObject = true;
                globalKey  = true;
            }
                
            if (UnOps.contains (_String(',')&curOp&',')) { // a standard function
                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else { // a variable
                // check if this is a function defined  in the list of "standard functions"
                long bLang = noneObject?-1:FunctionNameList.BinaryFind(&curOp);
                if (bLang>=0) {
                    levelOps->AppendNewInstance (new _Operation (curOp,FunctionArgumentCount(bLang)));
                    continue;
                }

                // check if this is a function defined in the batch language

                if ((bLang =  noneObject?-1:FindBFFunctionName (curOp, theParent))>=0) {
                    levelOps->AppendNewInstance (new _Operation (curOp,-bLang-1));
                    continue;
                }

                long curOpl = curOp.sLength;
                if (curOpl>2 && curOp[curOpl-1]=='_' && curOp[curOpl-2]=='_') { // instant variable refrence
                    _String realVarName (curOp,0,curOpl-3);

                    if (theParent) {
                        realVarName = *(theParent->GetName()) & '.' & realVarName;
                    }

                    long realVarLoc = LocateVarByName (realVarName);
                    if (realVarLoc<0) { // bad instant variable reference
                        if (flagErrors) {
                            WarnError (_String("Attempted to take value of undeclared variable:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                        }
                        return HY_FORMULA_FAILED;
                    }
                    if (!f2) { // 03/25/2004 ? Confused why the else
                        levelData->AppendNewInstance(new _Operation((_MathObject*)FetchVar (realVarLoc)->Compute()->makeDynamic()));
                    } else {
                        _Operation theVar (true, realVarName, globalKey, theParent);
                        theVar.SetTerms(-variableNames.GetXtra (realVarLoc)-1);
                        theVar.SetAVariable(-2);
                        (*levelData) && (&theVar);
                    }
                } else {
                    if (noneObject)
                        levelData->AppendNewInstance (new _Operation (false, curOp));
                    else
                        if (theParent && _hyApplicationGlobals.Find(&curOp) >= 0) {
                            levelData->AppendNewInstance (new _Operation(true, curOp, globalKey, nil));
                        } else {
                            levelData->AppendNewInstance (new _Operation(true, curOp, globalKey, theParent));
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
        if (numeric.isAllowed [s.getChar(i)]) {
            if (twoToken) {
                _String thisOp (s.getChar(i-1));
                levelOps->AppendNewInstance (new _Operation (thisOp,1L));
            }
            
            long j = 1;

            while ( i+j<s.sLength && (numeric.isAllowed [s.getChar(i+j)] || ((s.getChar(i+j)=='-' || s.getChar(i+j)=='+' )&& tolower(s.getChar(i+j-1))=='e')) ) {
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


            /*if (UnOps.contains(s.getChar(i)) && !twoToken)
            {
                char cim1 = s.getChar(i-1);

                if ( i == 0 || cim1=='(' || cim1=='[' || cim1==',')
                {
                    curOp = s.getChar(i);
                    levelOps->AppendNewInstance(new _Operation (curOp,1));
                    continue;
                }
            }*/

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
                    if (!twoToken && UnOps.contains(s.getChar(i))) {
                        twoOrOne = 1;
                    } else {
                        if (flagErrors) {
                            WarnError (_String("Bad Binary Operator Placement:")&s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                        }
                        return HY_FORMULA_FAILED;
                    }
                }
            }

            twoToken = false;

            if (levelData->countitems()) {
                int k;
                if (storage) {
                    BaseRef newS = (*levelData)(levelData->countitems()-1)->makeDynamic();
                    for (k = 0; k<levelData->countitems()-1; k++) {
                        f->theFormula&&((*levelData)(k));
                    }

                    levelData->Clear();
                    levelData->AppendNewInstance (newS);
                } else {
                    for (k = 0; k<levelData->countitems(); k++) {
                        f->theFormula&&((*levelData)(k));
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
        } else if (UnOps.contains(s.getChar(i))) {
            if ((s.getChar(i)=='-' || s.getChar(i)=='+') && (!i|| s.getChar(i-1)=='(')) { // unary minus?
                curOp   = s.getChar(i);
                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else {
                if (HalfOps.contains(s.getChar(i))) {
                    twoToken = true;
                    continue;
                }
                if (flagErrors) {
                    WarnError ((_String)"Bad unary operator placement " &s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return HY_FORMULA_FAILED;
            }
        } else {
            if (!HalfOps.contains(s.getChar(i))) {
                if (flagErrors) {
                    WarnError ((_String)"Bad symbols in expression " &s.Cut(0,i)&"?"&s.Cut(i+1,-1));
                }
                return HY_FORMULA_FAILED;
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
        _String namespcd = *namespc & '.' & name;
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
        _String namespcd = *namespc & '.' & name;
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
    _String * stIn;
    char    str[4096];

    for (long   i=0; i<indepVarList->lLength; i++) {
        _Variable *thisVar = LocateVar(indepVarList->lData[i]);
        if (thisVar->IsGlobal()) {
            sprintf (str, "\nglobal %s=%.16g;", thisVar->GetName()->getStr(),(double)thisVar->Compute()->Value());
            stIn = &glVars;
        } else {
            sprintf (str, "\n%s=%.16g;", thisVar->GetName()->getStr(),(double)thisVar->Compute()->Value());
            stIn = &locVars;
        }
        *stIn << str;
        if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
            sprintf (str, "\n%s:>%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetLowerBound());
            *stIn << str;
        }
        if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
            sprintf (str, "\n%s:<%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetUpperBound());
            *stIn << str;
        }
    }
}

//__________________________________________________________________________________

void ExportDepVariables (_String& glVars, _String& locVars, _SimpleList* depVarList)
{
    if (depVarList->lLength) {
        _String * stIn;
        char    str[4096];

        /* first we have to reorder global variables, so that dependent global variables which depend
           on other dependent global variables are written afterwards (lest they be implicitly declared 
           as local).
           The algorithm is very ugly, but since there are only a few global dependent variables (in general...) */

        _SimpleList     _globalVariablesList,
                        lfDepGlobs,
                        tl1;

        _List           dependancyLists;
        {
            for (long i=0; i<depVarList->lLength; i++)
                if (LocateVar(depVarList->lData[i])->IsGlobal()) {
                    lfDepGlobs << depVarList->lData[i];
                }
        }
        lfDepGlobs.Sort();

        for (long i=0; i<depVarList->lLength; i++) {
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
                sprintf (str, "\nglobal %s", thisVar->GetName()->getStr());
                stIn = &glVars;
            } else {
                sprintf (str, "\n%s", thisVar->GetName()->getStr());
                stIn = &locVars;
            }
            (*stIn)<<str;
            (*stIn)<<":=";
            _String* s = thisVar->GetFormulaString();
            (*stIn)<<s;
            DeleteObject(s);
            (*stIn)<<';';
            if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
                sprintf (str, "\n%s:>%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetLowerBound());
                (*stIn)<<str;
            }
            if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
                sprintf (str, "\n%s:<%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetUpperBound());
                (*stIn)<<str;
            }
        }

        if (_globalVariablesList.lLength)
            // check internal dependancies
        {
            _SimpleList writeOrder (_globalVariablesList.lLength,0,1),
                        indexList  (_globalVariablesList.lLength,0,1);


            for (long i2 = 0; i2 < _globalVariablesList.lLength; i2++) {
                long updatedIndex = writeOrder.lData[i2];
                _SimpleList * depList = (_SimpleList*)dependancyLists(i2);
                for (long i3 = 0; i3 < depList->lLength; i3 ++) {
                    long i4 = _globalVariablesList.Find (depList->lData[i3]);
                    if (i4 >= 0 && updatedIndex < writeOrder.lData[i4]) {
                        updatedIndex = writeOrder.lData[i4] + 1;
                    }
                }
                writeOrder.lData[i2] = updatedIndex;
            }

            SortLists (&writeOrder, &indexList);

            for (long i=0; i<_globalVariablesList.lLength; i++) {
                _Variable * thisVar = LocateVar(depVarList->lData[_globalVariablesList.lData[indexList.lData[i]]]);
                sprintf (str, "\nglobal %s", thisVar->GetName()->getStr());
                glVars<<str;
                glVars<<":=";
                _String* s = thisVar->GetFormulaString();
                glVars<<s;
                DeleteObject(s);
                glVars<<';';
                if (!CheckEqual(thisVar->GetLowerBound(),DEFAULTPARAMETERLBOUND)) {
                    sprintf (str, "\n%s:>%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetLowerBound());
                    glVars<<str;
                }
                if (!CheckEqual(thisVar->GetUpperBound(),DEFAULTPARAMETERUBOUND)) {
                    sprintf (str, "\n%s:<%.16g;", thisVar->GetName()->getStr(),(double)thisVar->GetUpperBound());
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
