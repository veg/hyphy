/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-2009
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon              (apoon@cfenet.ubc.ca)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

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

extern   _SimpleList    BinOps,
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

_Parameter  verbosityLevel = 0.0,
            twoOverSqrtPi   = 2./sqrtPi;


_SimpleList simpleOperationCodes,
            simpleOperationFunctions;

_String     internalRerootTreeID ("_INTERNAL_REROOT_TREE_");

long        subNumericValues = 0;

//__________________________________________________________________________________

_Parameter  InterpolateValue        (_Parameter*, _Parameter*, long, _Parameter*, _Parameter*, _Parameter, _Parameter&);
_Parameter  TrapezoidLevelK         (_Formula&, _Variable*, _Parameter, _Parameter, long);
_Parameter  TrapezoidLevelKSimple   (_Formula&, _Variable*, _Parameter, _Parameter, long, _SimpleFormulaDatum*, _SimpleFormulaDatum*, _SimpleList&, _SimpleList&);

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
//__________________________________________________________________________________

_Formula::_Formula (void)
{
    theTree     = nil;
    resultCache = nil;
}

//__________________________________________________________________________________

_Formula::_Formula (_PMathObj p, bool isAVar)
{
    theTree     = nil;
    resultCache = nil;
    if (!isAVar) {
        theFormula.AppendNewInstance (new _Operation (p));
    } else {
        _Variable* v = (_Variable*)p;
        theFormula.AppendNewInstance (new _Operation (true,*v->GetName(),v->IsGlobal(), nil));
    }
}

//__________________________________________________________________________________

void _Formula::Initialize (void) {}

//__________________________________________________________________________________

void _Formula::Duplicate  (BaseRef f)
{
    _Formula * f_cast = (_Formula*) f;

    theFormula.Duplicate       (& f_cast->theFormula);
    theStack.theStack.Duplicate(& f_cast->theStack.theStack);

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

void _Formula::DuplicateReference  (_Formula* f)
{
    for (int i=0; i<f->theFormula.lLength; i++) {
        _Operation *theO = ((_Operation*)f->theFormula(i));
        if (theO->GetAVariable()==-2) {
            theFormula.AppendNewInstance(new _Operation ((_PMathObj)LocateVar (-theO->GetNoTerms()-1)->Compute()->makeDynamic()));
        } else {
            theFormula&& theO;
        }
    }
}

//__________________________________________________________________________________
BaseRef _Formula::makeDynamic (void)
{
    _Formula * res = new _Formula;
    checkPointer(res);

    res->Duplicate((BaseRef)this);

    return (BaseRef)res;
}
//__________________________________________________________________________________

_Formula::~_Formula (void)
{

    Clear();
}

//__________________________________________________________________________________

void    _Formula::Clear (void)
{
    if (theTree) {
        theTree->delete_tree();
        delete theTree;
    }
    theTree = nil;
    if (resultCache) {
        DeleteObject (resultCache);
    }

    theFormula.Clear();
//  theStack.Clear();
}

//__________________________________________________________________________________

BaseRef _Formula::toStr (_List* matchedNames, bool dropTree)
{
    ConvertToTree();

    _String * result = (_String*)checkPointer(new _String((unsigned long)16,true));

    long          savepd = printDigits;
    printDigits          = 0;

    if (theTree) { // there is something to do
        internalToStr (*result, theTree, -1, matchedNames);
    } else {
        if (theFormula.lLength) {
            (*result) << "RPN:";
            internalToStr (*result,nil,0,nil,(_Operation*)(theFormula(0)));
            for (long k=1; k<theFormula.lLength; k++) {
                (*result)<<',';
                internalToStr (*result,nil,0,nil,(_Operation*)(theFormula(k)));
            }
        }
    }

    printDigits = savepd;
    result->Finalize ();
    if (theTree && dropTree) {
        theTree->delete_tree();
        delete theTree;
        theTree = nil;
    }
    return result;
}
//__________________________________________________________________________________

node<long>* _Formula::DuplicateFormula (node<long>* src, _Formula& tgt)
{
    node<long>* resNode = new node<long>;
    checkPointer (resNode);

    tgt.theFormula && (_Operation*) theFormula (src->in_object);

    resNode->in_object = tgt.theFormula.lLength-1;

    for (long k=1; k<=src->get_num_nodes(); k++) {
        resNode->add_node (*DuplicateFormula (src->go_down (k), tgt));
    }

    return     resNode;
}

//__________________________________________________________________________________

_Formula* _Formula::Differentiate (_String varName, bool bail)
{
    long          varID = LocateVarByName (varName),
                  k;

    if (varID<0) {
        return new _Formula (new _Constant (0.0));
    }

    varID = variableNames.GetXtra (varID);

    _Formula*     res = new _Formula ();
    checkPointer  (res);

    ConvertToTree    ();

    _SimpleList  varRefs,
                 dydx;


    {
        _AVLList al (&varRefs);
        ScanFForVariables (al, true, true, true);
        al.ReorderList ();
    }

    for (k=0; k < varRefs.lLength; k++) {
        _Variable* thisVar = LocateVar (varRefs.lData[k]);
        _Formula * dYdX;
        if (thisVar->IsIndependent()) {
            dYdX = new _Formula ((thisVar->GetName()->Equal (&varName))?new _Constant (1.0):new _Constant (0.0));
            checkPointer (dYdX);
            dYdX->ConvertToTree();
            dydx << (long)dYdX;
        } else {
            dYdX = thisVar->varFormula->Differentiate (varName);
            if (dYdX->theFormula.lLength == 0) {
                delete (dYdX);
                return res;
            }
            dydx << (long)dYdX;
        }
    }

    SortLists             (&varRefs, &dydx);
    node<long>*           dTree = nil;

    if (!(dTree = InternalDifferentiate (theTree, varID, varRefs, dydx, *res))) {
        for (k=0; k<dydx.lLength; k++) {
            delete ((_Formula*)dydx.lData[k]);
        }

        if (bail) {
            WarnError    (_String ("Differentiation of ") & _String((_String*)toStr()) & " failed.");
            res->Clear();
            return       res;
        } else {
            delete res;
            return nil;
        }
    }

    for (k=0; k<dydx.lLength; k++) {
        delete ((_Formula*)dydx.lData[k]);
    }

    res->theFormula.AppendNewInstance (new _Operation(new _Constant (0.0))) ;
    res->theTree         = dTree;
    res->InternalSimplify (dTree);
    res->ConvertFromTree  ();
    return res;

}

//__________________________________________________________________________________

node<long>* _Formula::InternalDifferentiate (node<long>* currentSubExpression, long varID, _SimpleList& varRefs, _SimpleList& dydx, _Formula& tgt)
{
    _Operation * op = (_Operation*)theFormula (currentSubExpression->in_object);

    if (op->theData!=-1) {
        long k     = varRefs.BinaryFind (op->GetAVariable());
        if (k<0) {
            return nil;
        }

        _Formula* dYdX = (_Formula*)dydx(k);
        return dYdX->DuplicateFormula (dYdX->theTree, tgt);
    }

    if (op->theNumber) {
        _Formula src (new _Constant (0.0));
        src.ConvertToTree ();
        return   src.DuplicateFormula (src.theTree, tgt);
    }

    node<long>* newNode = (node<long>*)checkPointer (new node<long>);


    switch (op->opCode) {
    case HY_OP_CODE_MUL: {
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                    * b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);

        if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
                b1->delete_tree (true);
            }
            if (b2) {
                b2->delete_tree (true);
            }
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('+');

        _Operation*       newOp  = new _Operation (opC2,2),
        *         newOp2 = new _Operation (opC ,2),
        *         newOp3 = new _Operation (opC ,2);

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);

        node<long>*       newNode2 = (node<long>*)checkPointer(new node<long>);
        node<long>*       newNode3 = (node<long>*)checkPointer(new node<long>);

        newNode2->add_node (*b1);
        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode3->add_node (*b2);
        newNode3->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode->add_node  (*newNode2);
        newNode->add_node  (*newNode3);

        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance (newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula. AppendNewInstance (newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance (newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_ADD: // +
    case HY_OP_CODE_SUB: { // -
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                    * b2 = nil;

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        long      isUnary = (currentSubExpression->get_num_nodes()==1);

        if (!isUnary) {
            b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);
            if (!b2) {
                b1->delete_tree      (true);
                newNode->delete_tree (true);
                return nil;
            }
        }


        _Operation*   newOp = new _Operation ();
        checkPointer  (newOp);
        newOp->Duplicate (op);
        newNode->add_node (*b1);
        if (!isUnary) {
            newNode->add_node (*b2);
        }
        newNode->in_object = tgt.theFormula.lLength;

        tgt.theFormula.AppendNewInstance(newOp);
        return          newNode;
    }
    break;

    case HY_OP_CODE_DIV: { // /
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                    * b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);

        if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
                b1->delete_tree (true);
            }
            if (b2) {
                b2->delete_tree (true);
            }
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('-'),
                          opC3 ('/'),
                          opC4 ('^');

        _Operation*       newOp  = new _Operation (opC3 ,2),
        *         newOp2 = new _Operation (opC4 ,2),
        *         newOp3 = new _Operation (opC2 ,2),
        *         newOp4 = new _Operation (opC  ,2),
        *         newOp5 = new _Operation (opC  ,2),
        *         newOp6 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);
        checkPointer      (newOp5);
        checkPointer      (newOp6);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;
        node<long>*       newNode6 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);
        checkPointer      (newNode5);
        checkPointer      (newNode6);

        newNode6->add_node (*b1);
        newNode6->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode5->add_node (*b2);
        newNode5->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode4->add_node  (*newNode6);
        newNode4->add_node  (*newNode5);

        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));
        newNode2->add_node (*newNode3);

        newNode->add_node  (*newNode4);
        newNode->add_node  (*newNode2);

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


        return          newNode;
    }
    break;

    case HY_OP_CODE_ARCTAN: { // Arctan
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('/'),
                          opC2 ('+'),
                          opC3 ('^');

        _Operation*       newOp  = new _Operation (opC ,2),
        *         newOp2 = new _Operation (opC2 ,2),
        *         newOp3 = new _Operation (new _Constant (1.0)),
        *         newOp4 = new _Operation (opC3 ,2),
        *         newOp5 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);
        checkPointer      (newNode5);

        newNode4->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
        newNode4->add_node (*newNode5);

        newNode2->add_node (*newNode3);
        newNode2->add_node (*newNode4);

        newNode->add_node (*b1);
        newNode->add_node (*newNode2);

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

        return          newNode;
    }
    break;

    case HY_OP_CODE_COS: { // Cos
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('-'),
                          opC3 ("Sin");

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,1),
        *         newOp3 = new _Operation (opC3 ,1);

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);

        newNode3->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode2->add_node (*newNode3);

        newNode->add_node  (*newNode2);
        newNode->add_node  (*b1);

        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula << newOp3;
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula << newOp2;
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula << newOp;

        DeleteObject    (newOp);
        DeleteObject    (newOp2);
        DeleteObject    (newOp3);

        return          newNode;
    }
    break;

    case HY_OP_CODE_ERF: { // Erf
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('/'),
                          opC3 ("Exp"),
                          opC4 ("-"),
                          opC5 ("^");

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,2),
        *         newOp3 = new _Operation (new _Constant (twoOverSqrtPi)),
        *         newOp4 = new _Operation (opC3 ,1),
        *         newOp5 = new _Operation (opC4 ,1),
        *         newOp6 = new _Operation (opC5 ,2),
        *         newOp7 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);
        checkPointer      (newOp5);
        checkPointer      (newOp6);
        checkPointer      (newOp7);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;
        node<long>*       newNode6 = new node<long>;
        node<long>*       newNode7 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);
        checkPointer      (newNode5);
        checkPointer      (newNode6);
        checkPointer      (newNode7);

        newNode6->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
        newNode6->add_node (*newNode7);

        newNode5->add_node (*newNode6);
        newNode4->add_node (*newNode5);
        newNode2->add_node (*newNode4);
        newNode2->add_node (*newNode3);

        newNode->add_node  (*b1);
        newNode->add_node  (*newNode2);

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

        return          newNode;
    }
    break;

    case HY_OP_CODE_EXP: // HY_OP_CODE_EXP
    case HY_OP_CODE_SIN: { // HY_OP_CODE_SIN
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('*'),
                          opC2;

        if (op->opCode==HY_OP_CODE_SIN) {
            opC2 = *(_String*)BuiltInFunctions(HY_OP_CODE_COS);
        } else {
            opC2 = *(_String*)BuiltInFunctions(HY_OP_CODE_EXP);
        }

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,1);

        checkPointer      (newOp);
        checkPointer      (newOp2);

        node<long>*       newNode2 = new node<long>;

        checkPointer      (newNode2);

        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode->add_node  (*newNode2);
        newNode->add_node  (*b1);

        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return           newNode;
    }
    break;

    case HY_OP_CODE_LOG: { // Log
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }
        _String           opC  ('/');

        _Operation*       newOp  = new _Operation (opC  ,2);

        checkPointer      (newOp);

        newNode->add_node  (*b1);
        newNode->add_node  (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return           newNode;
    }
    break;

    case HY_OP_CODE_SQRT: { // Sqrt
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('/'),
                          opC2 ('*'),
                          opC3 (*(_String*)BuiltInFunctions(HY_OP_CODE_SQRT));

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC2 ,2),
        *         newOp3 = new _Operation (opC3 ,1),
        *         newOp4 = new _Operation (new _Constant (2.0));

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);

        newNode3->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode2->add_node (*newNode4);
        newNode2->add_node (*newNode3);

        newNode->add_node  (*b1);
        newNode->add_node  (*newNode2);

        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_TAN: { // Tan
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt);

        if (!b1) {
            newNode->delete_tree (true);
            return nil;
        }

        _String           opC  ('/'),
                          opC2 ('^'),
                          opC3 (*(_String*)BuiltInFunctions(HY_OP_CODE_COS));

        _Operation*       newOp  = new _Operation (opC ,2),
        *         newOp2 = new _Operation (opC2,2),
        *         newOp3 = new _Operation (new _Constant (2.0)),
        *         newOp4 = new _Operation (opC3,1);

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);

        newNode4->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode2->add_node (*newNode4);
        newNode2->add_node (*newNode3);

        newNode->add_node (*b1);
        newNode->add_node (*newNode2);

        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula.AppendNewInstance(newOp);

        return          newNode;
    }
    break;

    case HY_OP_CODE_POWER: // ^
        // f[x]^g[x] (g'[x] Log[f[x]] + f'[x]g[x]/f[x])
    {
        node<long>* b1 = InternalDifferentiate (currentSubExpression->go_down(1), varID, varRefs, dydx, tgt),
                    * b2 = InternalDifferentiate (currentSubExpression->go_down(2), varID, varRefs, dydx, tgt);

        if (!b1 || !b2) {
            newNode->delete_tree(true);
            if (b1) {
                b1->delete_tree (true);
            }
            if (b2) {
                b2->delete_tree (true);
            }
            return nil;
        }

        _String           opC  ('*'),
                          opC2 ('+'),
                          opC3 ('/'),
                          opC4 ('^'),
                          opC5 (*(_String*)BuiltInFunctions(HY_OP_CODE_LOG));

        _Operation*       newOp  = new _Operation (opC  ,2),
        *         newOp2 = new _Operation (opC4 ,2),
        *         newOp3 = new _Operation (opC2 ,2),
        *         newOp4 = new _Operation (opC3 ,2),
        *         newOp5 = new _Operation (opC  ,2),
        *         newOp6 = new _Operation (opC  ,2),
        *         newOp7 = new _Operation (opC5, 1);

        checkPointer      (newOp);
        checkPointer      (newOp2);
        checkPointer      (newOp3);
        checkPointer      (newOp4);
        checkPointer      (newOp5);
        checkPointer      (newOp6);
        checkPointer      (newOp7);

        node<long>*       newNode2 = new node<long>;
        node<long>*       newNode3 = new node<long>;
        node<long>*       newNode4 = new node<long>;
        node<long>*       newNode5 = new node<long>;
        node<long>*       newNode6 = new node<long>;
        node<long>*       newNode7 = new node<long>;

        checkPointer      (newNode2);
        checkPointer      (newNode3);
        checkPointer      (newNode4);
        checkPointer      (newNode5);
        checkPointer      (newNode6);
        checkPointer      (newNode7);


        newNode7->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode6->add_node (*b2);
        newNode6->add_node (*newNode7);

        newNode5->add_node (*b1);
        newNode5->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode4->add_node (*newNode5);
        newNode4->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));

        newNode3->add_node (*newNode6);
        newNode3->add_node (*newNode4);

        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(1),tgt));
        newNode2->add_node (*DuplicateFormula (currentSubExpression->go_down(2),tgt));

        newNode->add_node  (*newNode2);
        newNode->add_node  (*newNode3);

        newNode7->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp7);
        newNode6->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp6);
        newNode5->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp5);
        newNode4->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp4);
        newNode3->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp3);
        newNode2->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp2);
        newNode->in_object = tgt.theFormula.lLength;
        tgt.theFormula .AppendNewInstance(newOp);

        return          newNode;
    }
    }

    delete (newNode);
    return nil;
}

//__________________________________________________________________________________

bool _Formula::InternalSimplify (node<long>* startNode)
// returns true if the subexpression at
// and below startnode is constant
{
    long        numChildren = startNode->get_num_nodes(),
                k,
                collapse2 = -1;

    bool        isConstant  = true,
                firstConst  = true,
                secondConst = (numChildren>1);

    _Parameter  theVal      = 0.0;

    _PMathObj   newVal      = nil;

    _Operation* op = (_Operation*)theFormula (startNode->get_data());

    if  (numChildren == 0) {
        return !op->IsAVariable();
    }

    for  (k=1; k<=numChildren; k++) {
        InternalSimplify (startNode->go_down(k));
        if (k==1) {
            firstConst = InternalSimplify (startNode->go_down(k));
        } else if (k==2) {
            secondConst = InternalSimplify (startNode->go_down(k));
        } else {
            isConstant = isConstant && InternalSimplify (startNode->go_down(k));
        }
    }

    isConstant = isConstant && firstConst && ((numChildren==1)||secondConst);

    if (op->opCode > HY_OP_CODE_NONE) {
        if (isConstant) {
            _Stack scrap;
            for  (k=1; k<=numChildren; k++) {
                ((_Operation*)theFormula (startNode->go_down(k)->get_data()))->Execute (scrap);
            }
            op->Execute (scrap);
            newVal = (_PMathObj)scrap.Pop()->makeDynamic();
        } else {
            if (firstConst||secondConst) {
                theVal  =((_Operation*)theFormula (startNode->go_down(firstConst?1:2)->get_data()))->GetANumber()->Value();

                switch (op->opCode) {
                case HY_OP_CODE_MUL: { // *
                    if (CheckEqual (theVal,0.0)) { // *0 => 0
                        newVal = new _Constant (0.0);
                        break;
                    }
                    if (CheckEqual (theVal,1.0)) { // ?*1 => ?
                        if (firstConst) {
                            collapse2 = 2;
                        } else {
                            collapse2 = 1;
                        }
                        break;
                    }
                }
                break;

                case HY_OP_CODE_ADD: { // +
                    if (CheckEqual (theVal,0.0)) { // ?+0 => ?
                        if (firstConst) {
                            collapse2 = 2;
                        } else {
                            collapse2 = 1;
                        }
                        break;
                    }
                }

                case HY_OP_CODE_SUB: { // -
                    if (CheckEqual (theVal,0.0)) {
                        if (firstConst) {
                            collapse2 = -2;    // 0-? => -?
                        } else {
                            collapse2 = 1;    // ?-0 => ?
                        }
                        break;
                    }
                }

                case HY_OP_CODE_DIV: { // /
                    if (firstConst&&CheckEqual (theVal,0.0)) { // 0/? => 0
                        newVal = new _Constant (0.0);
                        break;
                    }
                    if (secondConst&&CheckEqual (theVal,1.0)) { // ?/1 => ?
                        collapse2 = 1;
                        break;
                    }
                }
                break;

                case HY_OP_CODE_POWER: { // ^
                    if (firstConst&&CheckEqual (theVal,1.0)) { // 1^? => 1
                        newVal = new _Constant (1.0);
                        break;
                    }
                    if (secondConst&&CheckEqual (theVal,1.0)) { // ?^1 => ?
                        collapse2 = 1;
                        break;
                    }
                }
                break;
                }
            }
        }
    }

    if (newVal) {
        for  (k=numChildren; k; k--) {
            node<long>* theChild = startNode->go_down(k);
            delete (theChild);
            startNode->kill_node (k);
        }
        startNode->in_object = theFormula.lLength;
        theFormula.AppendNewInstance (new _Operation(newVal));
    }

    if (collapse2!=-1) {
        if (collapse2>0) {
            k = 3-collapse2;

            startNode->go_down(k)->delete_tree();
            delete (startNode->go_down(k));

            node <long>*    replaceWith = startNode->go_down(collapse2);

            startNode->kill_node(1);
            startNode->kill_node(2);

            for (k=1; k<=replaceWith->get_num_nodes(); k++) {
                startNode->add_node (*replaceWith->go_down(k));
            }

            startNode->in_object = replaceWith->in_object;

            delete (replaceWith);

        } else { // 0-? => -?
            delete   (startNode->go_down(1));
            startNode->kill_node(1);
        }
    }
    return false;
}


//__________________________________________________________________________________

void _Formula::internalToStr (_String& result, node<long>* currentNode, char opLevel, _List* matchNames, _Operation* thisNodeOperation)
{
    if (!thisNodeOperation) {
        thisNodeOperation = (_Operation*)theFormula (currentNode->get_data());
    }

    // decide what to do about this operation

    if (thisNodeOperation->IsAVariable(false))
        // this operation is just a variable - add ident to string and return
    {
        if (subNumericValues) {
            if (subNumericValues == 2) {
                _Variable* theV = LocateVar(thisNodeOperation->GetAVariable());
                if  (_x_&&(theV->GetAVariable()==_x_->GetAVariable())) {
                    result << _x_->GetName();
                    return;
                }
            }

            _Variable *thisVariable = LocateVar(thisNodeOperation->GetAVariable());
            _PMathObj subThisValue = thisVariable->Compute();

            if (subThisValue->ObjectClass () == NUMBER) {
                if (subNumericValues == 3) {
                    result << LocateVar(thisNodeOperation->GetAVariable())->GetName();
                    result << '[';
                    result.AppendNewInstance(new _String (subThisValue->Value()));
                    result << ':';
                    result.AppendNewInstance(new _String (thisVariable->GetLowerBound()));
                    result << '-';
                    result.AppendNewInstance(new _String (thisVariable->GetUpperBound()));
                    result << ']';

                } else {
                    result.AppendNewInstance(new _String (subThisValue->Value()));
                }
            } else if (subThisValue->ObjectClass () == STRING) {
                result.AppendNewInstance((_String*)subThisValue->toStr());
            } else {
                result << LocateVar(thisNodeOperation->GetAVariable())->GetName();
            }
        } else {
            long variableIDX = thisNodeOperation->GetAVariable();
            if (variableIDX>=0) {
                _String * vName = LocateVar(variableIDX)->GetName();
                if (matchNames) {
                    _List * p1 = (_List*)(*matchNames)(0),
                            * p2 = (_List*)(*matchNames)(1);

                    long  f = p1->Find (vName);

                    if (f<0) {
                        result<<vName;
                    } else {
                        result<<(_String*)(*p2)(f);
                    }
                } else {
                    result<<vName;
                }
            }
        }
        return;
    }

    long nOps = thisNodeOperation->GetNoTerms();
    if (nOps>0)
        // a built-in operation or a function call
    {
        // check if it's a built-in binary operation
        _String  opString (thisNodeOperation->GetCode());
        long f = BinOps.Find ((opString.sLength>1)?(opString.sData[0]*256+opString.sData[1]):opString.sData[0]);
        if (f!=-1)
            // indeed - a binary operation is what we have. check if need to wrap the return in parentheses
        {
            if ((!currentNode)||(currentNode->get_num_nodes()==2)) {
                char tOpLevel  = opPrecedence(f),
                     tOpLevel2 = tOpLevel;

                if (associativeOps.Find(f)<0) {
                    tOpLevel2 ++;
                }
                if (opLevel>=0) { // need to worry about op's precedence
                    {
                        bool parens = opPrecedence(f)<opLevel;

                        if (parens&&currentNode) { // put parentheses around the return of this expression
                            result<<'(';
                            internalToStr (result, currentNode->go_down(1),tOpLevel,matchNames);
                            result<<&thisNodeOperation->GetCode();
                            internalToStr (result, currentNode->go_down(2),tOpLevel2,matchNames);
                            result<<')';
                            return;
                        }
                    }
                }
                if (currentNode) {
                    internalToStr (result, currentNode->go_down(1),tOpLevel,matchNames);
                }
                result<<&thisNodeOperation->GetCode();
                if (currentNode) {
                    internalToStr (result, currentNode->go_down(2),tOpLevel2,matchNames);
                }
                return;
            } else { // mixed binary-unary operation
                result<<&thisNodeOperation->GetCode();
                if (currentNode) {
                    result<<'(';
                    internalToStr (result, currentNode->go_down(1),opPrecedence(f),matchNames);
                    result<<')';
                }
                return;
            }
        } else {
            _String mac ("MAccess");
            if (!thisNodeOperation->GetCode().Equal(&mac)) {
                result<<&thisNodeOperation->GetCode();
                if (currentNode) {
                    result<<'(';
                    for (long k=1; k<=nOps; k++) {
                        if (k>1) {
                            result<<',';
                        }
                        internalToStr (result, currentNode->go_down(k),-1,matchNames);
                    }
                    result<<')';
                }
            } else { // matrix element access - treat specially
                if (currentNode) {
                    internalToStr (result, currentNode->go_down(1),-1,matchNames);
                    for (long k=2; k<=nOps; k++) {
                        result<<'[';
                        internalToStr (result, currentNode->go_down(k),-1,matchNames);
                        result<<']';
                    }
                }
            }
        }
        return;
    }
    if (nOps<0)
        // a user-defined function
    {
        result<<((_String*)batchLanguageFunctionNames(-nOps-1));
        if (currentNode) {
            result<<'(';
            for (long k=1; k<=batchLanguageFunctionParameters(-nOps-1); k++) {
                if(k>1) {
                    result<<',';
                }
                internalToStr (result, currentNode->go_down(k),-1,matchNames);
            }
            result<<')';
        }
        return;
    }
    _PMathObj opValue = thisNodeOperation->GetANumber();
    _String* conv = (_String*)opValue->toStr();
    if (opValue->ObjectClass()==STRING) {
        result<<'"';
        result<<conv;
        result<<'"';
    } else {
        result<<conv;
    }
    DeleteObject(conv);
}
//__________________________________________________________________________________

bool     _Formula::IsEmpty(void)
// is there anything in the formula
{
    return bool(!theFormula.lLength);
}

//__________________________________________________________________________________

_Parameter   _Formula::Newton(_Formula& derivative, _Variable* unknown, _Parameter targetValue, _Parameter left, _Parameter right)
// find a root of the formulaic expression, using Newton's method, given the derivative and a bracketed root.
// will alternate between bisections and Newton iterations based on what is fatser
{
    // check that there is indeed a sign change on the interval
    _Constant   dummy;
    _Parameter  t1,t2,t3,t4,t5,lastCorrection = 100., newCorrection;
    _String     msg;
    long        iterCount = 0;
    dummy.SetValue (left);
    unknown->SetValue(&dummy);
    t1 = Compute()->Value()-targetValue;
    if (t1==0.0) {
        return left;
    }
    dummy.SetValue (right);
    unknown->SetValue(&dummy);
    t2 = Compute()->Value()-targetValue;
    if (t2==0) {
        return right;
    }
    if (t1*t2>0.0) {
        subNumericValues = 2;
        _String *s = (_String*)toStr();
        subNumericValues = 0;
        _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&",Inf)";
        ReportWarning (msg);
        DeleteObject (s);
        return    left;
    }
    // else all is good we can start the machine
    bool useNewton = false;
    while ((fabs((right-left)/MAX(left,right))>machineEps*10.)&&(iterCount<200)) { // stuff to do
        iterCount ++;
        if (!useNewton) {
            t3 = (right+left)/2;
        }
        dummy.SetValue(t3);
        unknown->SetValue(&dummy);
        t4 = Compute()->Value()-targetValue;
        // get the correction term from the derivative
        dummy.SetValue(t3);
        unknown->SetValue(&dummy);
        t5 = derivative.Compute()->Value();
        useNewton = true;
        if (t5==0.0) {
            useNewton = false;
        } else {
            newCorrection = -t4/t5;
            if (t3) {
                if (fabs(newCorrection/t3)<machineEps*2.) { // correction too small - the root has been found
                    return t3;
                }
            } else if (fabs(newCorrection)<machineEps*2.) { // correction too small - the root has been found
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
            dummy.SetValue(t3);
            unknown->SetValue(&dummy);
            t4 = Compute()->Value()-targetValue;
            if (t4==0.0) {
                return t3;
            }
            if (t4*t1 >0) {
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

#define MAX_BRENT_ITERATES 100L

//__________________________________________________________________________________

_Parameter   _Formula::Brent(_Variable* unknown, _Parameter a, _Parameter b, _Parameter tol, _List* storeEvals, _Parameter rhs)
// find a root of the formulaic expression, using Brent's method and a bracketed root.
{
    // check that there is indeed a sign change on the interval
    _Constant   dummy;

    _Parameter  fa = 0.0,fb = 0.0,fc,d,e,min1,min2,xm,p,q,r,s,tol1,
                c = b;

    min1 = unknown->GetLowerBound();
    min2 = unknown->GetUpperBound();

    long        it = 0;

    if (a>b) { // autobracket to the left
        dummy.SetValue (b);
        unknown->SetValue(&dummy);
        fb = Compute()->Value();
        if (storeEvals) {
            (*storeEvals) && & dummy;
            dummy.SetValue (fb);
            (*storeEvals) && & dummy;
        }

        if (b<0.00001 && b>-0.00001) {
            a = b-0.0001;
        } else {
            a = b-fabs(b)*0.1;
        }

        if (a<min1) {
            a = min1+0.5*(b-min1);
        }

        dummy.SetValue (a);
        unknown->SetValue(&dummy);
        fa = Compute()->Value()-rhs;
        if (storeEvals) {
            (*storeEvals) && & dummy;
            dummy.SetValue (fa);
            (*storeEvals) && & dummy;
        }

        for (long k=0; k<50; k++) {
            if (fb*fa<0.0) {
                break;
            }

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

            dummy.SetValue (a);
            unknown->SetValue(&dummy);
            fa = Compute()->Value()-rhs;
            if (storeEvals) {
                (*storeEvals) && & dummy;
                dummy.SetValue (fa);
                (*storeEvals) && & dummy;
            }
        }
    } else if (CheckEqual (a,b)) { // autobracket to the right
        dummy.SetValue (a);
        unknown->SetValue(&dummy);
        fa = Compute()->Value()-rhs;

        if (storeEvals) {
            (*storeEvals) && & dummy;
            dummy.SetValue (fa);
            (*storeEvals) && & dummy;
        }
        a = b;

        if ((b<0.00001)&&(b>-0.00001)) {
            b = b+0.0001;
        } else {
            b = b+fabs(b)*0.1;
        }

        if (b>min2) {
            b = a+0.5*(min2-a);
        }

        dummy.SetValue (b);
        unknown->SetValue(&dummy);
        fb = Compute()->Value()-rhs;

        if (storeEvals) {
            (*storeEvals) && & dummy;
            dummy.SetValue (fb);
            (*storeEvals) && & dummy;
        }

        for (long k=0; k<50; k++) {
            if (fb*fa<0.0) {
                break;
            }

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

            dummy.SetValue (b);
            unknown->SetValue(&dummy);
            fb = Compute()->Value()-rhs;
            if (storeEvals) {
                (*storeEvals) && & dummy;
                dummy.SetValue (fb);
                (*storeEvals) && & dummy;
            }
        }
    }


    if (fa == 0.0) {
        dummy.SetValue (a);
        unknown->SetValue(&dummy);
        fa = Compute()->Value()-rhs;
        if (storeEvals) {
            (*storeEvals) && & dummy;
            dummy.SetValue (fa);
            (*storeEvals) && & dummy;
        }
        if (fa == 0.0) {
            return a;
        }
    }

    if (fb == 0.0) {
        dummy.SetValue (b);
        unknown->SetValue(&dummy);
        fb = Compute()->Value()-rhs;
        if (storeEvals) {
            (*storeEvals) && & dummy;
            dummy.SetValue (fb);
            (*storeEvals) && & dummy;
        }
        if (fb==0) {
            return b;
        }
    }

    if (fa*fb<0.0) {
        fc = fb;

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

            tol1 = 2.*fabs(b)*machineEps+.5*tol;

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
                b+=fabs(tol1)*(1.*2.*(xm>=.0));
            }

            dummy.SetValue (b);
            unknown->SetValue(&dummy);
            fb = Compute()->Value()-rhs;
            if (storeEvals) {
                (*storeEvals) && & dummy;
                dummy.SetValue (fb);
                (*storeEvals) && & dummy;
            }
        }
    }

    subNumericValues = 2;
    _String msg ((_String*)toStr());
    subNumericValues = 0;
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

_Parameter   _Formula::Newton(_Formula& derivative, _Parameter targetValue, _Parameter left, _Parameter max_right, _Variable* unknown)
// given a monotone function and a left bracket bound, found the right bracket bound and solve
{
    // check that there is indeed a sign change on the interval
    _Constant   dummy;
    dummy.SetValue (left);
    unknown->SetValue(&dummy);
    _Parameter  t1 = Compute()->Value(), right = left, t2, step = 1.0;

    if (max_right-left < step * 100) {
        step = (max_right-left)/100;
    }
    if  (step==0.0) {
        return left;
    }
    do {
        right += step;
        if (right>max_right) { // function doesn't seem to have a root
            subNumericValues = 2;
            _String *s = (_String*)toStr();
            subNumericValues = 0;
            _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&","&_String(right)&")";
            ReportWarning (msg);
            DeleteObject (s);
            return    left;
        }
        dummy.SetValue (right);
        unknown->SetValue(&dummy);
        t2 = Compute()->Value();
        step*=2;
        if (right+step>max_right)
            if (right<max_right) {
                step = max_right-right;
            }
    } while ((targetValue-t1)*(targetValue-t2)>0);
    return Newton (derivative, unknown, targetValue, left, right);

}

//__________________________________________________________________________________

_Parameter   _Formula::Newton(_Variable* unknown, _Parameter targetValue, _Parameter x_min, _Parameter left, _Parameter right)
{
    // check that there is indeed a sign change on the interval
    _Constant   dummy;
    _Parameter  t1,t2,t3,t4,t5,lastCorrection = 100, newCorrection;
    _String     msg;
    t1 =Integral(unknown, x_min, left)-targetValue;
    if (t1==0.0) {
        return left;
    }
    t2 = t1+Integral(unknown, left, right);
    if (t2==0) {
        return right;
    }
    if (t1*t2>0.0) {
        subNumericValues = 2;
        _String *s = (_String*)toStr();
        subNumericValues = 0;
        _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&","
                      &_String (right)&"]";
        ReportWarning (msg);
        DeleteObject (s);
        return    left;
    }
    // else all is good we can start the machine
    bool useNewton = false;
    while (right-left>1e-6) { // stuff to do
        if (!useNewton) {
            t3 = (right+left)/2;
        }
        dummy.SetValue(t3);
        unknown->SetValue(&dummy);
        t4 = Integral(unknown, x_min, t3)-targetValue;
        // get the correction term from the derivative
        dummy.SetValue(t3);
        unknown->SetValue(&dummy);
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
            t4 = Integral(unknown, x_min, t3)-targetValue;
            if (t4==0.0) {
                return t3;
            }
            if (t4*t1 >0) {
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

//__________________________________________________________________________________

_Parameter   _Formula::Newton( _Variable* unknown, _Parameter targetValue,_Parameter x_min, _Parameter left)
// given a monotone function and a left bracket bound, found the right bracket bound and solve
{
    // check that there is indeed a sign change on the interval
    _Parameter  t1 = Integral(unknown, x_min, left), right = left, t2, step = 1.0;
    do {
        right += step;
        t2 = Integral(unknown, right-step, right);
        step*=2;
        if (right>=1e10) { // function doesn't seem to have a root
            subNumericValues = 2;
            _String *s = (_String*)toStr();
            subNumericValues = 0;
            _String msg = *s&"="&_String(targetValue)&" has no (or multiple) roots in ["&_String(left)&",Inf)";
            WarnError (msg);
            DeleteObject (s);
            return    0.0;
        }
    } while ((targetValue-t1)*(targetValue-t2-t1)>=0);
    return Newton ( unknown, targetValue,x_min, left, right);

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


_Parameter   _Formula::Integral(_Variable* dx, _Parameter left, _Parameter right, bool infinite)
// uses Romberg's intergation method
{
    if (infinite) { // tweak "right" here
        _Parameter value = 1.0, step = 1.0, right1=-1;
        right = left;
        while (value>machineEps) {
            right+=step;
            _Constant dummy (right);
            dx->SetValue(&dummy);
            value = fabs(Compute()->Value());
            if ((value<0.001)&&(right1<0)) {
                right1 = right;
            }
            step *= 2;
            if (step>100000) { // integrand decreasing too slowly
                _String msg, *s = (_String*)toStr();
                msg = *s & " decreases too slowly to be integrated over an infinite interval";
                DeleteObject(s);
                WarnError (msg);
                return 0.0;
            }
        }
        if (right1<right-machineEps) {
            return Integral(dx,left,right1,false)+Integral(dx,right1,right,false);
        } else {
            return Integral(dx,left,right1,false);
        }
    }

    checkParameter(intPrecFact,integrationPrecisionFactor,1e-4);
    checkParameter(intMaxIter,maxRombergSteps,8);

    _Parameter ss,
               dss,
               *s,
               *h;

    s = new _Parameter [(long)maxRombergSteps];
    h = new _Parameter [(long)(maxRombergSteps+1)];
    checkPointer(s);
    checkPointer(h);

    h[0]=1.0;

    long         interpolateSteps = 5,
                 stackDepth = 0;

    _SimpleList  fvidx,
                 changingVars,
                 idxMap;

    _Parameter   * ic = new _Parameter[interpolateSteps],
    * id = new _Parameter[interpolateSteps];

    _SimpleFormulaDatum
    * stack = nil,
      * vvals = nil;

    checkPointer (ic);
    checkPointer (id);

    if (AmISimple (stackDepth,fvidx)) {
        stack = new _SimpleFormulaDatum [stackDepth];
        checkPointer (stack);
        vvals = new _SimpleFormulaDatum [fvidx.lLength];
        checkPointer (vvals);
        ConvertToSimple (fvidx);
        stackDepth = dx->GetAVariable();
        for (long vi = 0; vi < fvidx.lLength; vi++) {
            _Variable * checkvar = LocateVar (fvidx.lData[vi]);
            if (checkvar->CheckFForDependence (stackDepth,true)) {
                changingVars << fvidx.lData[vi];
                idxMap << vi;
            }
            vvals[vi].value = checkvar->Compute()->Value();
        }
        changingVars.InsertElement ((BaseRef)stackDepth,0,false,false);
        idxMap.InsertElement ((BaseRef)fvidx.Find(stackDepth),0,false,false);
    } else {
        stackDepth = -1;
    }

    for (long j=0; j<(long)maxRombergSteps; j++) {
        if (stackDepth>=0) {
            s[j] = TrapezoidLevelKSimple(*this, dx, left, right, j+1, stack, vvals,changingVars,idxMap);
        } else {
            s[j] = TrapezoidLevelK(*this, dx, left, right, j+1);
        }
        if (j>=4) {
            ss = InterpolateValue(&h[j-4],&s[j-4],interpolateSteps,ic,id,0.0, dss);
            if (fabs(dss)<= integrationPrecisionFactor * fabs(ss)) {
                delete s;
                delete h;
                delete ic;
                delete id;
                if (stackDepth>=0) {
                    ConvertFromSimple(fvidx);
                    delete (stack);
                    delete (vvals);
                }
                return ss;
            }
        }
        h[j+1] =  h[j]/9.0;
    }

    if (stackDepth>=0) {
        ConvertFromSimple(fvidx);
        delete (stack);
        delete (vvals);
    }
    _String *str = (_String*)toStr(),
             msg = _String("Integral of ")&*str & " over ["&_String(left)&","&_String(right)&"] converges slowly, loss of precision may occur. Change either INTEGRATION_PRECISION_FACTOR or INTEGRATION_MAX_ITERATES";

    DeleteObject (str);
    ReportWarning (msg);

    delete s;
    delete h;
    delete ic;
    delete id;
    return ss;
}

//__________________________________________________________________________________

_Parameter   _Formula::MeanIntegral(_Variable* dx, _Parameter left, _Parameter right, bool infinite)
{
    _Formula newF;
    _String tim ("*");
    _Operation times (tim,2);
    _Operation term (true, *(dx->GetName()));
    newF.Duplicate((BaseRef)this);
    newF.theFormula&& (&term);
    newF.theFormula&& (& times);
    return newF.Integral (dx,left,right, infinite);
}

//__________________________________________________________________________________

long     _Formula::NumberOperations(void)
// number of operations in the formula
{
    return theFormula.lLength;
}

//__________________________________________________________________________________

_PMathObj _Formula::Compute (long startAt, _VariableContainer * nameSpace) // compute the value of the formula
{
    if (theFormula.lLength == 0) {
        return nil;
    }

    bool wellDone = true;

    if (startAt == 0) {
        theStack.theStack.Clear();
    }

    if (startAt == 0 && resultCache && resultCache->lLength) {
        long cacheID     = 0;
        bool cacheResult = false;

        for (long i=0; i<theFormula.lLength; i++) {
            _Operation* thisOp ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
            if (i < theFormula.lLength-1) {
                _Operation* nextOp  ((_Operation*)(((BaseRef**)theFormula.lData)[i+1]));

                if (! cacheResult && nextOp->CanResultsBeCached(thisOp)) {
                    if (!thisOp->Execute(theStack,nameSpace)) {
                        wellDone = false;
                        break;
                    }

                    _Matrix *currentArg = (_Matrix*)theStack.Pop(false),
                             *cachedArg  = (_Matrix*)((_PMathObj)(*resultCache)(cacheID)),
                              *diff     = nil;

                    if (cachedArg->ObjectClass() == MATRIX) {
                        diff =  (_Matrix*)cachedArg->SubObj(currentArg);
                    }

                    if (diff && diff->MaxElement() <= 1e-12) {
                        DeleteObject  (theStack.Pop  ());
                        theStack.Push ((_PMathObj)(*resultCache)(cacheID+1));
                        cacheID += 2;
                        i ++;
                    } else {
                        cacheResult = true;
                        resultCache->Replace(cacheID++,theStack.Pop(false),true);
                    }
                    DeleteObject (diff);
                    continue;
                }
            }
            if (!thisOp->Execute(theStack,nameSpace)) {
                wellDone = false;
                break;
            }
            if (cacheResult) {
                resultCache->Replace(cacheID++,theStack.Pop(false),true);
                cacheResult = false;
            }
        }
    } else {
        for (long i=startAt; i<theFormula.lLength; i++)
            if (!((_Operation*)(((BaseRef**)theFormula.lData)[i]))->Execute(theStack, nameSpace)) {
                wellDone = false;
                break;
            }
    }
    if (theStack.theStack.lLength != 1 || !wellDone) {
        WarnError (_String((_String*)toStr()) & _String(" contains errors."));
        return    new _Constant (0.0);
    }


    return theStack.Pop(false);
}

//__________________________________________________________________________________

bool _Formula::CheckSimpleTerm (_PMathObj thisObj)
{
    if (thisObj) {
        long oc = thisObj->ObjectClass();
        if (oc !=NUMBER) {
            if (oc ==MATRIX) {
                _Matrix * mv = (_Matrix*)thisObj->Compute();
                if (mv->IsIndependent () && !mv->SparseDataStructure()) {
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

void _Formula::ConvertMatrixArgumentsToSimpleOrComplexForm (bool makeComplex)
{
    if (makeComplex) {
        if (resultCache) {
            DeleteObject (resultCache);
            resultCache = nil;
        }
    } else {
        if (!resultCache) {
            resultCache = new _List();
            for (int i=1; i<theFormula.lLength; i++) {
                _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
                if (thisOp->CanResultsBeCached(((_Operation*)(((BaseRef**)theFormula.lData)[i-1])))) {
                    resultCache->AppendNewInstance(new _MathObject());
                    resultCache->AppendNewInstance(new _MathObject());
                }
            }
        }
    }

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));

        _Matrix   * thisMatrix = nil;

        if (thisOp->theNumber) {
            if (thisOp->theNumber->ObjectClass() == MATRIX) {
                thisMatrix = (_Matrix*)thisOp->theNumber;
            }
        } else {
            if (thisOp->theData>-1) {
                _Variable* thisVar = LocateVar (thisOp->theData);
                if (thisVar->ObjectClass() == MATRIX) {
                    thisMatrix = (_Matrix*)thisVar->GetValue();
                }
            }
        }

        if (thisMatrix)
            if (makeComplex) {
                thisMatrix->MakeMeGeneral();
            } else {
                thisMatrix->MakeMeSimple();
            }
    }
}

//__________________________________________________________________________________

bool _Formula::AmISimple (long& stackDepth, _SimpleList& variableIndex)
{
    if (!theFormula.lLength) {
        return true;
    }

    long locDepth = 0;

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
        locDepth++;
        if ( thisOp->theData<-2 || thisOp->numberOfTerms<0) {
            return false;
        }

        if (thisOp->theNumber) {
            if (thisOp->theNumber->ObjectClass() != NUMBER) {
                return false;
            }
        } else {
            if (thisOp->theData>-1) {
                _Variable* thisVar = LocateVar (thisOp->theData);
                if (thisVar->ObjectClass()!=NUMBER) {
                    _PMathObj cv = thisVar->GetValue();
                    if (!CheckSimpleTerm (cv)) {
                        return false;
                    }
                }
                if (variableIndex.Find (thisOp->theData)==-1) {
                    variableIndex<<thisOp->theData;
                }
            } else {
                if (simpleOperationCodes.Find(thisOp->opCode)==-1) {
                    return false;
                } else if (thisOp->opCode == HY_OP_CODE_MACCESS && thisOp->numberOfTerms != 2) {
                    return false;
                }

                locDepth-=thisOp->numberOfTerms;
            }
        }
        if (locDepth>stackDepth) {
            stackDepth = locDepth;
        } else if (locDepth==0) {
            _String errStr = _String("Invalid formula passed to _Formula::AmISimple") & _String ((_String*)toStr());
            WarnError (errStr);
            return false;
        }
    }
    return true;
}

//__________________________________________________________________________________

void _Formula::ConvertToSimple (_SimpleList& variableIndex)
{
    if (theFormula.lLength)
        for (int i=0; i<theFormula.lLength; i++) {
            _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
            if (thisOp->theNumber) {
                continue;
            } else if (thisOp->theData >= 0) {
                thisOp->theData = variableIndex.Find (thisOp->theData);
            } else if (thisOp->opCode == HY_OP_CODE_SUB && thisOp->numberOfTerms == 1) {
                thisOp->opCode = (long)MinusNumber;
            } else {
                if (thisOp->opCode == HY_OP_CODE_MACCESS) {
                    thisOp->numberOfTerms = -2;
                }
                thisOp->opCode = simpleOperationFunctions(simpleOperationCodes.Find(thisOp->opCode));
            }

        }
}

//__________________________________________________________________________________

void _Formula::ConvertFromSimple (_SimpleList& variableIndex)
{
    if (!theFormula.lLength) {
        return;
    }

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
        if (thisOp->theNumber) {
            continue;
        } else {
            if (thisOp->theData>-1) {
                thisOp->theData = variableIndex[thisOp->theData];
            } else if (thisOp->opCode == (long)MinusNumber) {
                thisOp->opCode = HY_OP_CODE_SUB;
            } else {
                if (thisOp->opCode == (long)FastMxAccess) {
                    thisOp->numberOfTerms = 2;
                }
                thisOp->opCode = simpleOperationCodes(simpleOperationFunctions.Find(thisOp->opCode));
            }
        }
    }
}

//__________________________________________________________________________________

_Parameter _Formula::ComputeSimple (_SimpleFormulaDatum* stack, _SimpleFormulaDatum* varValues)
{
    if (!theFormula.lLength) {
        return 0.0;
    }

    long stackTop = 0;

    for (int i=0; i<theFormula.lLength; i++) {
        _Operation* thisOp = ((_Operation*)(((BaseRef**)theFormula.lData)[i]));
        if (thisOp->theNumber) {
            stack[stackTop++].value = thisOp->theNumber->Value();
            continue;
        } else {
            if (thisOp->theData>-1) {
                stack[stackTop++] = varValues[thisOp->theData];
            } else {
                stackTop--;
                if (thisOp->numberOfTerms==2) {
                    _Parameter  (*theFunc) (_Parameter, _Parameter);
                    theFunc = (_Parameter(*)(_Parameter,_Parameter))thisOp->opCode;
                    if (stackTop<1) {
                        _String errMsg = "Internal error in _Formula::ComputeSimple - stack underflow.)";
                        WarnError (errMsg);
                        return 0.0;
                    }
                    stack[stackTop-1].value = (*theFunc)(stack[stackTop-1].value,stack[stackTop].value);
                } else {
                    if (thisOp->numberOfTerms==-2) {
                        _Parameter  (*theFunc) (Ptr,_Parameter);
                        theFunc = (_Parameter(*)(Ptr,_Parameter))thisOp->opCode;
                        if (stackTop<1) {
                            _String errMsg = "Internal error in _Formula::ComputeSimple - stack underflow.)";
                            WarnError (errMsg);
                            return 0.0;
                        }
                        stack[stackTop-1].value = (*theFunc)(stack[stackTop-1].reference,stack[stackTop].value);
                    } else {
                        _Parameter  (*theFunc) (_Parameter);
                        theFunc = (_Parameter(*)(_Parameter))thisOp->opCode;
                        stack[stackTop++].value = (*theFunc)(stack[stackTop].value);
                    }
                }
            }
        }
    }

    return stack->value;
}

//__________________________________________________________________________________

bool _Formula::EqualFormula (_Formula* f)
{
    if (theFormula.lLength == f->theFormula.lLength) {
        for (int i=0; i<theFormula.lLength; i++) {
            if (!((_Operation*)(((BaseRef**)theFormula.lData)[i]))->EqualOp (((_Operation*)(((BaseRef**)f->theFormula.lData)[i])))) {
                return false;
            }
        }
        return true;
    }
    return false;
}

//__________________________________________________________________________________

_PMathObj _Formula::ConstructPolynomial (void) // compute the value of the formula
{
    theStack.Reset();
    bool wellDone = true;

    for (long i=0; i<theFormula.lLength; i++)
        if (!((_Operation*)((BaseRef**)theFormula.lData)[i])->ExecutePolynomial(theStack)) {
            wellDone = false;
            break;
        }

    if (theStack.theStack.lLength !=1 || !wellDone) {
        return    nil;
    }

    return theStack.Pop(false);
}
//__________________________________________________________________________________

bool _Formula::HasChanged (bool ingoreCats)
{
    _Operation *thisOp;
    long dataID;
    for (int i = 0; i<theFormula.lLength; i++) {
        thisOp = (_Operation*)((BaseRef**)theFormula.lData)[i];
        if (thisOp->IsAVariable()) {
            dataID = thisOp->GetAVariable();
            if (dataID>=0) {
                if (((_Variable*)(((BaseRef*)(variablePtrs.lData))[dataID]))->HasChanged(ingoreCats)) {
                    return true;
                }
            } else if (thisOp->theNumber->HasChanged()) {
                return true;
            }
        } else if (thisOp->opCode == HY_OP_CODE_BRANCHLENGTH||thisOp->opCode == HY_OP_CODE_RANDOM||thisOp->opCode == HY_OP_CODE_TIME)
            // time, random or branch length
        {
            return true;
        } else if (thisOp->numberOfTerms<0) {
            dataID = -thisOp->numberOfTerms-2;
            if (dataID<batchLanguageFunctionClassification.lLength) {
                if (batchLanguageFunctionClassification.lData[dataID] == BL_FUNCTION_NORMAL_UPDATE) {
                    continue;
                }
            }
            return true;
        }
    }
    return false;
}

//__________________________________________________________________________________

bool _Formula::HasChangedSimple (_SimpleList& variableIndex)
{
    _Operation *thisOp;
    for (long i = 0; i<theFormula.lLength; i++) {
        thisOp = (_Operation*)((BaseRef**)theFormula.lData)[i];
        if (thisOp->theNumber) {
            continue;
        } else if (thisOp->theData >= 0) {
            if (((_Variable*)(((BaseRef*)(variablePtrs.lData))[variableIndex.lData[thisOp->theData]]))->HasChanged(false)) {
                return true;
            }
        } else {
            long codeOp = simpleOperationFunctions(simpleOperationCodes.Find(thisOp->opCode));
            if (codeOp == HY_OP_CODE_RANDOM || codeOp == HY_OP_CODE_TIME) {
                return true;
            }
        }
    }
    return false;
}


//__________________________________________________________________________________

void _Formula::ScanFForVariables (_AVLList&l, bool includeGlobals, bool includeAll, bool includeCategs, bool skipMatrixAssignments)
{
    for (long i = 0; i<theFormula.lLength; i++) {
        _Operation* theObj = ((_Operation**)theFormula.lData)[i];
        if (theObj->IsAVariable()) {
            if (!includeGlobals)
                if ((((_Variable*)LocateVar(theObj->GetAVariable()))->IsGlobal())||
                        (((_Variable*)LocateVar(theObj->GetAVariable()))->ObjectClass()!=NUMBER)) {
                    continue;
                }

            long f = theObj->GetAVariable();

            if (f>=0) {
                _Variable * v = LocateVar(f);

                if (v->IsCategory()&&includeCategs) {
                    v->ScanForVariables (l,includeGlobals);
                }

                if(includeAll || v->ObjectClass()==NUMBER) {
                    l.Insert ((BaseRef)f);
                }

                if (skipMatrixAssignments) {
                    if (v->ObjectClass()!=MATRIX || !theObj->AssignmentVariable()) {
                        v->ScanForVariables(l,includeGlobals);
                    }
                } else if (!v->IsIndependent()) {
                    v->ScanForVariables(l,includeGlobals);
                }
            } else if (theObj->theNumber)
                if (theObj->theNumber->ObjectClass()==MATRIX) {
                    ((_Matrix*)theObj->theNumber)->ScanForVariables(l,includeGlobals);
                }
        }
    }
}

//__________________________________________________________________________________

void _Formula::ScanFForType (_SimpleList &l, int type)
{
    for (long i = 0; i<theFormula.lLength; i++) {
        _Operation* theObj = ((_Operation**)theFormula.lData)[i];
        if (theObj->IsAVariable()) {
            long f = theObj->GetAVariable();

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

bool _Formula::CheckFForDependence (long varID, bool checkAll)
{
    for (int i = 0; i<theFormula.lLength; i++) {
        _Operation* theObj = (_Operation*)theFormula(i);
        long f;
        if (theObj->IsAVariable()) {
            f = theObj->GetAVariable();
            if (f>=0) {
                if (f == varID) {
                    return true;
                }
                if (checkAll) {
                    _Variable * v = LocateVar(f);
                    if (!v->IsIndependent())
                        if (v->CheckFForDependence(varID)) {
                            return true;
                        }
                }
            }
        }
    }
    return false;
}

//__________________________________________________________________________________

void  _Formula::LocalizeFormula (_Formula& ref, _String& parentName, _SimpleList& iv, _SimpleList& iiv, _SimpleList& dv, _SimpleList& idv)
{
    for (int i = 0; i<ref.theFormula.lLength; i++) {
        if (((_Operation*)ref.theFormula(i))->IsAVariable()) {
            long     vIndex = ((_Operation*)ref.theFormula(i))->GetAVariable();
            _Variable* theV = LocateVar (vIndex);
            if (theV->IsGlobal()) {
                theFormula&& ref.theFormula(i);
                continue;
            }
            if (theV->IsContainer()) {
                continue;
            }
            _String  fullName = parentName&"."&*(theV->GetName());
            long lvIndex = LocateVarByName(fullName);
            if (lvIndex==-1)
                // local variable doesn't yet exist - create
            {
                _Variable dummy (fullName);
                lvIndex = LocateVarByName(fullName);
                if (theV->IsIndependent()) {
                    iv<<lvIndex;
                    iiv<<vIndex;
                } else {
                    dv<<lvIndex;
                    idv<<vIndex;
                }
            }
            _Operation newVar (true, fullName);
            theFormula && &newVar;
        } else {
            theFormula&& ref.theFormula(i);
        }
    }
}

//__________________________________________________________________________________

bool _Formula::DependsOnVariable (long idx)
{
    for (long f = 0; f<theFormula.lLength; f++) {
        _Operation * anOp = ((_Operation**)theFormula.lData)[f];
        if (anOp->IsAVariable() && anOp->GetAVariable() == idx) {
            return true;
        }
    }

    return false;
}

//__________________________________________________________________________________

_Operation* _Formula::GetIthTerm (long idx)
{
    if (idx >= 0 && idx < theFormula.lLength) {
        return ((_Operation**)theFormula.lData)[idx];
    }

    return nil;
}


//__________________________________________________________________________________

bool _Formula::IsAConstant (void)
{
    for (int i = 0; i<theFormula.lLength; i++)
        if (((_Operation*)((BaseRef**)theFormula.lData)[i])->IsAVariable()) {
            return false;
        }

    return true;
}

//__________________________________________________________________________________

bool _Formula::IsConstant (void)
{
    for (int i = 0; i<theFormula.lLength; i++)
        if (((_Operation*)((BaseRef**)theFormula.lData)[i])->IsConstant() == false) {
            return false;
        }

    return true;
}

//__________________________________________________________________________________

void _Formula::SimplifyConstants (void)
{
    theStack.theStack.Clear();
    for (long i = 0; i<theFormula.countitems(); i++) {
        long j;
        _Operation* thisOp = ((_Operation*)((BaseRef**)theFormula.lData)[i]);
        if ((thisOp->theData==-1)&&(thisOp->opCode>=0)&&(thisOp->numberOfTerms)) {
            long nt = thisOp->numberOfTerms;
            if (nt<0) {
                nt = batchLanguageFunctionParameters.lData[-nt-1];
            }

            for (j = 1; j<=nt; j++) {
                _Operation*  aTerm = ((_Operation*)((BaseRef**)theFormula.lData)[i-j]);
                if ((aTerm->IsAVariable())||(aTerm->opCode>=0)) {
                    break;
                }
            }

            if (j>nt)
                // all terms are constant -> Evaluate
            {
                for (j=i-thisOp->numberOfTerms; j<=i; j++) {
                    ((_Operation*)((BaseRef**)theFormula.lData)[j])->Execute(theStack);
                }
                long n = i-thisOp->numberOfTerms;
                thisOp = new _Operation (theStack.Pop());
                for (j=n; j<=i; j++) {
                    theFormula.Delete (n);
                }
                theFormula.InsertElement (thisOp,n,false);
                i = n+1;
                theStack.theStack.Clear();
                thisOp->nInstances--;
            } else {
                if (thisOp->numberOfTerms > 0 &&
                        (thisOp->opCode==HY_OP_CODE_MUL||thisOp->opCode==HY_OP_CODE_DIV||thisOp->opCode==HY_OP_CODE_POWER))
                    // *,/,^ 1 can be removed
                {
                    _Operation*  aTerm = ((_Operation*)((BaseRef**)theFormula.lData)[i-1]);
                    if (!((aTerm->IsAVariable())||(aTerm->opCode>=0))) {
                        if (aTerm->theNumber->ObjectClass()==NUMBER && aTerm->theNumber->Value() == 1.) {
                            theFormula.Delete (i);
                            theFormula.Delete (i-1);
                            i--;
                        }
                    }
                }
            }
        }
    }
}

//__________________________________________________________________________________

_PMathObj _Formula::GetTheMatrix (void)
{
    if (theFormula.lLength==1) {
        _Operation* firstOp = (_Operation*)theFormula(0);
        _PMathObj   ret = firstOp->GetANumber();
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

long _Formula::ObjectClass (void)
{
    if (theStack.theStack.lLength) {
        return ((_PMathObj)theStack.theStack.lData[0])->ObjectClass();
    }

    _PMathObj res =   Compute();

    if (res) {
        return res->ObjectClass();
    }

    return HY_UNDEFINED;
}

//__________________________________________________________________________________

_FString::_FString (void)
{
    theString = new _String;
}

//__________________________________________________________________________________

_FString::~_FString ()
{
    if (nInstances<=1) {
        DeleteObject (theString);
    } else {
        nInstances--;
    }

}

//__________________________________________________________________________________

_FString::_FString (_String* data)
{
    theString = data;
}

//__________________________________________________________________________________

_FString::_FString (long inData)
{
    checkPointer (theString = new _String (inData));
}

//__________________________________________________________________________________

_FString::_FString (_String& data, bool meta)
{
    if (meta) {
        unsigned long ssi = _String::storageIncrement;

        if (data.sLength>ssi) {
            _String::storageIncrement = data.sLength;
        }

        theString = new _String (data.sLength,true);
        for (long k=0; k<data.sLength; k++) {
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

BaseRef _FString::makeDynamic (void)
{
    //_FString* res = new _FString;
    //res->theString->Duplicate(theString);
    return new _FString(*theString, false);
}

//__________________________________________________________________________________

void _FString::Duplicate (BaseRef o)
{
    DeleteObject (theString);
    _FString* m = (_FString*)o;
    theString = (_String*)m->theString->makeDynamic();
}

//__________________________________________________________________________________

_PMathObj _FString::Add (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        _String * res    = new _String (theString->sLength+theStr->theString->sLength,true);
        if (res) {
            (*res) << theString;
            (*res) << theStr->theString;
            res->Finalize();
            return new _FString (res);
        } else {
            checkPointer (nil);
        }
    }
    _String* convStr = (_String*)p->toStr();
    _String res (*theString& *((_String*)convStr));
    DeleteObject (convStr);
    return new _FString (res, false);
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
            checkPointer (theString = new _String (s, true));
            return s;
        } else {
            theString->Finalize();
        }
    } else {
        WarnError ("Invalid 2nd argument in call to string*number");
    }
    return 0;
}



//__________________________________________________________________________________

_PMathObj _FString::AreEqual (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
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
        return new _Constant ((_Parameter)equal);
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
        return new _Constant ((_Parameter)equal);
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
            Ptr regex = PrepRegExp (theStr->theString, errNo, true);
            if (regex) {
                theString->RegExpMatchAll(regex, matches);
                FlushRegExp (regex);
            } else {
                WarnError (GetRegExpError (errNo));
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
        WarnError ("Invalid 2nd argument in call to string$reg.exp.");
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
                Ptr regex = PrepRegExp (theStr->theString, errNo, true);

                if (!regex) {
                    WarnError (GetRegExpError (errNo));
                    return new _FString (empty);
                }

                theString->RegExpMatchAll(regex, matches);
                _FString * res;
                if (matches.lLength) {
                    _String * newString = new _String (theString->sLength+1,true);
                    checkPointer (newString);
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

        WarnError ("Invalid 2nd argument in call to string^{{pattern,replacement}}");
    }
    return new _FString (empty,false);
}

//__________________________________________________________________________________

_PMathObj _FString::NotEqual (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)!equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)!equal);
    }
}

//__________________________________________________________________________________

_PMathObj _FString::Less (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Less(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Less(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________

_PMathObj _FString::LessEq (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Less(theStr->theString)||theString->Equal(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Less(convStr)||theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
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
    long    varRef;
    if (Parse (&F,copyMe, varRef) == HY_FORMULA_EXPRESSION) {
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
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Greater(convStr)||theString->Equal(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________

_PMathObj _FString::Greater (_PMathObj p)
{
    if (p->ObjectClass()==STRING) {
        _FString* theStr = (_FString*)p;
        bool     equal = theString->Greater(theStr->theString);
        return new _Constant ((_Parameter)equal);
    } else {
        _String* convStr = (_String*)p->toStr();
        bool     equal = theString->Greater(convStr);
        DeleteObject (convStr);
        return new _Constant ((_Parameter)equal);
    }
}

//__________________________________________________________________________________

BaseRef  _FString::toStr()
{
    theString->nInstances++;
    return theString;
}

//__________________________________________________________________________________


_PMathObj _FString::RerootTree (void)
{
    long     stashedModelID = lastMatrixDeclared,
             totalNodeCount = 0;

    lastMatrixDeclared      = HY_NO_MODEL;
    /* unset current model; do not want the internal tree to have an attached model */

    _TheTree    rTree (internalRerootTreeID,*theString);

    if (rTree.IsDegenerate()) { // no need to reroot two-sequence trees
        lastMatrixDeclared = stashedModelID;
        DeleteVariable  (internalRerootTreeID);
        return new _FString (*theString, false);
    }

    if (terminateExecution) {
        lastMatrixDeclared = stashedModelID;
        DeleteVariable  (internalRerootTreeID);
        return new _FString;
    }

    _CalcNode   *iterator = rTree.DepthWiseTraversal (true),
                 *rerootAt;

    node<long>  *cNode;

    _GrowingVector  valueCache;

    while       (iterator)
        // count the number of descendants of a given node, store as numeric value of the CalcNode
    {
        cNode    = &rTree.GetCurrentNode();
        valueCache.Store(iterator->Value());
        if (long myNodeCount = cNode->get_num_nodes()) {
            _Parameter tNodeCount = 0.0;

            for (long k = 1; k <= myNodeCount; k++) {
                tNodeCount += ((_CalcNode*)LocateVar(cNode->go_down(k)->in_object))->Value();
            }

            iterator->SetNumericValue(tNodeCount+1.0);
        } else {
            iterator->SetNumericValue(1.0);
        }

        iterator = rTree.DepthWiseTraversal (false);
        totalNodeCount ++;
    }

    iterator = rTree.DepthWiseTraversal (true);

    long        maxMin = 0;
    _Parameter  bRatio  = 0.0;

    while       (iterator) {
        _Parameter      nodeMin   = totalNodeCount-iterator->Value(),
                        thisRatio = nodeMin/(_Parameter)iterator->Value();

        if (thisRatio>1.0) {
            thisRatio = 1./thisRatio;
        }

        cNode    = &rTree.GetCurrentNode();
        if (cNode->get_num_nodes()) {
            for (long k = cNode->get_num_nodes(); k; k--) {
                long tt = ((_CalcNode*)LocateVar(cNode->go_down(k)->in_object))->Value();
                if (tt<nodeMin) {
                    nodeMin = tt;
                }
            }
        } else {
            nodeMin = 1;
        }

        if ((nodeMin>maxMin)||((nodeMin==maxMin)&&(thisRatio>bRatio))) {
            bRatio = thisRatio;
            maxMin = nodeMin;
            rerootAt = iterator;
            if (!cNode->get_parent()) {
                rerootAt = nil;
            }
        }
        iterator = rTree.DepthWiseTraversal (false);
    }

    iterator        = rTree.DepthWiseTraversal (true);
    totalNodeCount  = 0;
    while       (iterator)
        // restore branch lengths
    {
        iterator->SetNumericValue(valueCache.theData[totalNodeCount]);
        iterator = rTree.DepthWiseTraversal (false);
        totalNodeCount ++;
    }

    _FString* res;
    if (rerootAt) {
        _String stringCopy = *rerootAt->GetName();
        stringCopy.Trim (stringCopy.FindBackwards ('.',0,-1)+1,-1);
        _FString    rAt  (stringCopy);
        res = (_FString*)rTree.RerootTree (&rAt);
    } else {
        res = new _FString (*theString, false);
    }


    DeleteVariable  (internalRerootTreeID);

    lastMatrixDeclared = stashedModelID;

    return      res;
}

//__________________________________________________________________________________


_PMathObj _FString::Evaluate ()
{
    if (theString && theString->sLength) {
        _String     s (*theString);
        _Formula    evaluator (s);
        _PMathObj   evalTo = evaluator.Compute();

        if (evalTo && !terminateExecution) {
            evalTo->AddAReference();
            return evalTo;
        }
    }
    return new _Constant (.0);
}



//__________________________________________________________________________________


_PMathObj _FString::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{
    switch (opCode) {
    case HY_OP_CODE_NOT: // !
        return FileExists();
    case HY_OP_CODE_NEQ: // !=
        return NotEqual(p);
        break;
    case HY_OP_CODE_IDIV: // $ match regexp
        return EqualRegExp(p);
        break;
    case HY_OP_CODE_MOD: // % equal case insenstive
        return AreEqualCIS(p);
        break;
    case HY_OP_CODE_AND: { // && upcase or lowercase
        _Parameter pVal = 0.0;
        if (p->ObjectClass () == NUMBER) {
            pVal = p->Value();
        }

        if (pVal < 0.0) {
            return (_PMathObj)makeDynamic();
        } else {
            _String * t = nil;


            if (CheckEqual(pVal,2.0) || CheckEqual(pVal,3.0) || CheckEqual(pVal,4.0) || CheckEqual(pVal,5.0) || CheckEqual(pVal,6.0)) {
                checkPointer (t = new _String (theString->sLength+1,true));
                t->EscapeAndAppend (*theString, CheckEqual(pVal,3.0) + 2*CheckEqual(pVal,4.0) + 4*CheckEqual(pVal,5.0) + 5*CheckEqual(pVal,6.0));
                t->Finalize();
            } else {
                t = new _String (*theString);
                checkPointer (t);
                if (CheckEqual(pVal,1.0)) {
                    t->UpCase();
                } else {
                    t->LoCase();
                }
            }

            return new _FString (t);
        }
    }
    break;
    case HY_OP_CODE_MUL: // *
        if (p->ObjectClass() == MATRIX) {
            return      MapStringToVector (p);
        } else {
            return new _Constant(AddOn(p));
        }
        break;
    case HY_OP_CODE_ADD: // +
        if (p) {
            return Add(p);
        } else {
            return Sum();
        }
        break;
    case HY_OP_CODE_DIV: // /
        return EqualAmb(p);
        break;
    case HY_OP_CODE_LESS: // <
        return Less(p);
        break;
    case HY_OP_CODE_LEQ: // <=
        return LessEq(p);
        break;
    case HY_OP_CODE_EQ: // ==
        return AreEqual(p);
        break;
    case HY_OP_CODE_GREATER: // >
        return Greater(p);
        break;
    case HY_OP_CODE_GEQ: // >=
        return GreaterEq(p);
        break;
    case HY_OP_CODE_ABS: // Abs
        return new _Constant (theString->sLength);
        break;
    case HY_OP_CODE_DIFF: // Differentiate
        return Differentiate(p);
        break;
    case HY_OP_CODE_EVAL: // Eval
        return Evaluate();
        break;
    case HY_OP_CODE_EXP: // Exp
        return new _Constant (theString->LempelZivProductionHistory(nil));
        break;
    case HY_OP_CODE_FORMAT: { // Format
        _String cpyString (*theString);
        _Formula f (cpyString);
        _PMathObj fv = f.Compute();
        if (fv && fv->ObjectClass () == NUMBER) {
            return ((_Constant*)fv)->FormatNumberString (p,p2);
        } else {
            ReportWarning (_String("Failed to evaluate ")& *theString & " to a number in call to Format (string...)");
            return new _FString();
        }
    }
    break;
    case HY_OP_CODE_INVERSE: { // Inverse
        _FString * res = new _FString (*theString, false);
        checkPointer (res);
        for (long i1 = 0, i2 = theString->sLength-1; i1<theString->sLength; i1++, i2--) {
            res->theString->sData[i1] = theString->sData[i2];
        }

        return res;
    }
    break;
    case HY_OP_CODE_JOIN: // Inverse
        return Join (p);

    case HY_OP_CODE_LOG: // Log - check sum
        return new _Constant (theString->Adler32());
    case HY_OP_CODE_MACCESS: // MAccess
        return CharAccess(p,p2);
        break;
    case HY_OP_CODE_REROOTTREE: // RerootTree
        return RerootTree ();
        break;
    case HY_OP_CODE_ROWS: // Count Objects of given type
        return CountGlobalObjects();
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    case HY_OP_CODE_POWER: // Replace (^)
        return ReplaceReqExp (p);
        break;
    case HY_OP_CODE_OR: // Match all instances of the reg.ex (||)
        return EqualRegExp (p, true);
        break;
    }

    WarnNotDefined (this, opCode);
    return new _FString;

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
                    char thisChar = aKey->theString->sData[0];
                    if (mapper[thisChar] < 0) {
                        mapper[thisChar] = r;
                    }
                }
            }

            _SimpleList mapped;
            for (long s = 0; s < theString->sLength; s++) {
                mapped << mapper[theString->sData[s]];
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
_PMathObj   _FString::FileExists (void)
{
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
_PMathObj   _FString::CountGlobalObjects (void)
{
    _Parameter res = 0.0;

    long      standardType = _HY_GetStringGlobalTypes.Find(theString);
    if (standardType >=0 ) {
        standardType = _HY_GetStringGlobalTypes.GetXtra (standardType);
    }

    switch (standardType) {
    case 0:
        return new _Constant (likeFuncList.lLength);
    case 1:
        return new _Constant (dataSetList.lLength);
    case 2:
        return new _Constant (dataSetFilterList.lLength);
    case 3:
        return new _Constant (batchLanguageFunctionNames.lLength);
    case 4: {
        _SimpleList tc;
        long        si,
                    vi = variableNames.Traverser (tc,si,variableNames.GetRoot());

        for (; vi >= 0; vi = variableNames.Traverser (tc,si))
            if (((_Variable*)FetchVar(vi))->ObjectClass () == TREE) {
                res += 1.;
            }

        break;
    }

    case 5:
        return new _Constant (scfgList.lLength);
    case 6:
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
            isAllowed [s.sData[r2]] = true;
        }
    }
    bool     isAllowed [256];
}
alpha       ("ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz"),
            numeric     (".0123456789eE");

_String     globalToken ("global");

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

            if (lookAtMe ==',' && (level<1 || squareBrackets.lLength && squareBrackets.lData[squareBrackets.lLength-1] == level)) {
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

            if (UnOps.contains (_String(',')&curOp&',')) { // a standard function
                levelOps->AppendNewInstance (new _Operation (curOp,1));
                continue;
            } else { // a variable
                // check if this is a function defined  in the list of "standard functions"
                long bLang = FunctionNameList.BinaryFind(&curOp);
                if (bLang>=0) {
                    levelOps->AppendNewInstance (new _Operation (curOp,FunctionArgumentCount(bLang)));
                    continue;
                }

                // check if this is a function defined in the batch language



                if ((bLang = FindBFFunctionName (curOp, theParent))>=0) {
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
        if ( BinOps.Find (s.getChar(i))!=-1 || twoToken&& (BinOps.Find(s.getChar(i-1)*(long)256+s.getChar(i))!=-1) ) {
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
_Formula::_Formula (_String&s, _VariableContainer* theParent, bool errors)
// the parser itself
{
    theTree     = nil;
    resultCache = nil;

    long      refV;
    if (Parse (this, s, refV, theParent,nil,errors)== HY_FORMULA_FAILED) {
        Clear();
    }
}

//__________________________________________________________________________________
void    _Formula::ConvertToTree (void)
{
    if (!theTree&&theFormula.lLength) { // work to do
        _SimpleList nodeStack;
        _Operation* currentOp;
        for (long i=0; i<theFormula.lLength; i++) {
            currentOp = (_Operation*)theFormula(i);
            if (currentOp->TheCode()<0) { // a data bit
                node<long>* leafNode = new node<long>;
                checkPointer(leafNode);
                leafNode->init(i);
                nodeStack<<(long)leafNode;
            } else { // an operation
                long nTerms = currentOp->GetNoTerms();
                if (nTerms<0) {
                    nTerms = batchLanguageFunctionParameters(-nTerms-1);
                }

                if (nTerms>nodeStack.lLength) {
                    _String errMsg ("Expression syntax/semantics error.");
                    WarnError (errMsg);
                    theTree = nil;
                    return;
                }

                node<long>* operationNode = new node<long>;
                checkPointer(operationNode);
                operationNode->init(i);
                for (long j=0; j<nTerms; j++) {
                    operationNode->prepend_node(*((node<long>*)nodeStack(nodeStack.lLength-1)));
                    nodeStack.Delete(nodeStack.lLength-1);
                }
                nodeStack<<(long)operationNode;
            }
        }
        if (nodeStack.lLength!=1) {
            WarnError ((_String)"Expression syntax/sematics error.");
            theTree = nil;
            return;
        } else {
            theTree = (node<long>*)nodeStack(0);
        }
    }
}

//__________________________________________________________________________________
void    _Formula::ConvertFromTree (void)
{
    if (theTree) { // work to do
        _SimpleList termOrder;
        node<long>* currentNode = DepthWiseStepTraverser (theTree);
        while (currentNode) {
            termOrder<<currentNode->get_data();
            currentNode = DepthWiseStepTraverser ((node<long>*)nil);
        }
        if (termOrder.lLength!=theFormula.lLength) { // something has changed
            _List newFormula;
            for (long i=0; i<termOrder.lLength; i++) {
                newFormula<<theFormula(termOrder(i));
            }
            theFormula.Clear();
            theFormula.Duplicate(&newFormula);
            theTree->delete_tree();
            delete (theTree);
            theTree = nil;
            ConvertToTree();
        }
    }
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
