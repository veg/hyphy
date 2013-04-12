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

#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "polynoml.h"
#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif


//_String polyBOperations ("+-*^"), polyUOperations ("-");

_Parameter dropPrecision                = log(1e-9),
           drop2Precision               = log(1e-18),
           topPolyCap                   = 5.0,
           dropTerms                    = 1.0,
           enforcePolyCap               = 1.0,
           *varCheckArray               = nil,
            dropThreshold               = 0.0,
            maximumPolyTermsPerVariable = 500,
            maxPolynomialExpIterates    = 20,
            polynomialExpPrecision      = 1e-10;

long       varCheckAllocated            = 0,
           polyTermCap                  = 0x1fffffff;

bool       checkReset                   = true;


//__________________________________________________________________________________

_Polynomial::_Polynomial (void)
{
    theTerms = new _PolynomialData;
    checkPointer(theTerms);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_SimpleList& vars)
{
    variableIndex.Duplicate (&vars);
    theTerms = new _PolynomialData (vars.countitems());
    checkPointer(theTerms);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_Parameter value)
// a constant polynomial
{
    theTerms = new _PolynomialData;
    checkPointer(theTerms);
    theTerms->AddTerm (nil, value);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_Variable& v)
// a monic monomial
{
    variableIndex<<v.GetIndex();
    theTerms = new _PolynomialData (1);
    checkPointer(theTerms);
    long  vIndex = 1;
    theTerms->AddTerm (&vIndex,1.0);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_Polynomial& p)
{
    variableIndex.Duplicate (&p.variableIndex);
    theTerms = new _PolynomialData;
    checkPointer(theTerms);
    if (p.theTerms) {
        theTerms->Duplicate (p.theTerms);
    } else {
        theTerms->numberVars = variableIndex.countitems();
    }
    compList1.Duplicate (&p.compList1);
    compList2.Duplicate (&p.compList2);
}

//__________________________________________________________________________________

_Polynomial::~_Polynomial ()
{
    if (theTerms) {
        DeleteObject(theTerms);
    }
}

//__________________________________________________________________________________

BaseObj*    _Polynomial::makeDynamic(void)
{
    _Polynomial* res = new _Polynomial;
    checkPointer(res);

    res->variableIndex.Duplicate (&variableIndex);
    res->compList1.Duplicate (&compList1);
    res->compList2.Duplicate (&compList2);

    if (theTerms) {
        res->theTerms->Duplicate (theTerms);
    } else {
        DeleteObject (res->theTerms);
        res->theTerms = nil;
    }
    return res;
}

//__________________________________________________________________________________

void    _Polynomial::Duplicate  (BaseRef tp)
{
    _Polynomial* p = (_Polynomial*)tp;
    variableIndex.Clear();
    variableIndex.Duplicate (&p->variableIndex);
    compList1.Duplicate (&p->compList1);
    compList2.Duplicate (&p->compList2);
    DeleteObject(theTerms);
    if (p->theTerms) {
        theTerms =  new _PolynomialData (*(p->theTerms));
        checkPointer(theTerms);
    }
}
//__________________________________________________________________________________


_PMathObj _Polynomial::Execute (long opCode, _PMathObj p, _PMathObj, _hyExecutionContext* context)   // execute this operation with the second arg if necessary
{
    switch (opCode) {
    case HY_OP_CODE_MUL: //*
        if (p)
            return Mult(p);
        break;
    case HY_OP_CODE_ADD: // +
        if (p) {
            return Add(p);
        } else {
            return Sum ();
        }
        break;
    case HY_OP_CODE_SUB: // -
        if (p) {
            return Sub(p);
        } else {
            return Minus();
        }
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    case HY_OP_CODE_POWER: // ^
        if (p)
            Raise(p);
        break;
    }

    WarnNotDefined (this, opCode, context);
    return nil;

}

//__________________________________________________________________________________
_MathObject*    _Polynomial::IsANumber (bool returnLeading)
{
    long nV = variableIndex.countitems();
    if (!nV) {
        if (theTerms->NumberOfTerms()>0) {
            return new _Constant(theTerms->theCoeff[0]);
        } else {
            return new _Constant(0.0);
        }

    }
    if (theTerms->NumberOfTerms()<=1) {
        if (theTerms->NumberOfTerms()==0) {
            return new _Constant(0.0);
        }
        if (theTerms->IsFirstANumber() || returnLeading) {
            return new _Constant(theTerms->theCoeff[0]);
        }
    }
    return nil;
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Add (_MathObject* m)
{
    return Plus(m);
}
//__________________________________________________________________________________

void    ResetPolynomialCheck (_Polynomial* p)
{
    if (dropTerms) {
        if(!enforcePolyCap) {
            if (varCheckAllocated!=p->variableIndex.countitems()) {
                if (varCheckArray) {
                    free (varCheckArray);
                }
                varCheckAllocated=p->variableIndex.countitems();
                varCheckArray = (_Parameter*)MemAllocate (varCheckAllocated*sizeof(_Parameter));
                _Parameter lb, ub;
                for (long j=varCheckAllocated-1; j>=0; j--) {
                    _Variable* theV = LocateVar(p->variableIndex(j));
                    lb = fabs(theV->GetLowerBound());
                    ub = fabs(theV->GetUpperBound());
                    varCheckArray[j]=log(lb>ub?lb:ub);
                }
            }
        }
    }
    checkReset = true;

}

//__________________________________________________________________________________
void _Polynomial::CheckTerm (void)
{
    long myIndex = theTerms->actTerms-1;
    _Parameter myCoeff = theTerms->GetCoeff(myIndex);
    if (!theTerms->checkTerm (myCoeff,myIndex)) {
        theTerms->DeleteTerm(myIndex);
    }

}

//__________________________________________________________________________________

bool         _Polynomial::Equal(_MathObject* m)
{
    bool result = false;
    if (m->ObjectClass() == POLYNOMIAL || m->ObjectClass() == NUMBER) {
        _Polynomial * diff = (_Polynomial *)Sub(m);
        if (diff) {
            _Constant * v = (_Constant*)diff->IsANumber(true);
            if (v!=nil) {
                result = fabs (v->Value()) < 1.e-6;
                DeleteObject (v);
            }
            //_String * diffS = (_String*)diff->toStr();
            //printf ("%s\n", diffS->getStr());
            DeleteObject (diff);
        }

    }
    return result;
}


//__________________________________________________________________________________
_MathObject* _Polynomial::Plus (_MathObject* m, bool subtract)
{
    long objectT = m->ObjectClass();

    if (objectT==1) { // a number
        Convert2OperationForm();
        _Polynomial* result = new _Polynomial (*this);
        checkPointer(result);
        _Parameter mV = subtract?-m->Value():m->Value();

        if (!variableIndex.lLength) { // constant poly
            if (theTerms->NumberOfTerms()) {
                result->theTerms->GetCoeff(0)+=mV;
                if (result->theTerms->GetCoeff(0)==0.0) {
                    result->theTerms->DeleteTerm(0);
                }
            } else {
                if (mV) {
                    result->theTerms->AddTerm(nil,mV);
                }
            }
        } else { // variable containing poly
            if (theTerms->NumberOfTerms()) {
                if (theTerms->IsFirstANumber()) {
                    result->theTerms->GetCoeff(0)+=mV;
                    if (result->theTerms->GetCoeff(0)==0.0) {
                        result->theTerms->DeleteTerm(0);
                    }
                } else {
                    result->theTerms->AddTerm(mV);
                }
            } else {
                result->theTerms->AddTerm (mV);
            }
        }
        return result;
    }

    if (objectT == POLYNOMIAL) { // another polynomial
        Convert2OperationForm();
        _Polynomial* p2 = (_Polynomial*)m;

        if (variableIndex.lLength == 0) {
            if (theTerms->NumberOfTerms()) {
                _Constant coef1 (theTerms->GetCoeff(0));
                if (!subtract) {
                    return p2->Plus (&coef1,false);
                } else {
                    _Polynomial *ppp = (_Polynomial*)p2->Plus (&coef1,true);
                    _Parameter  *invC = ppp->theTerms->theCoeff;
                    for (long inv = 0; inv<ppp->theTerms->actTerms; inv++, invC++) {
                        (*invC) *= -1.0;
                    }
                    return ppp;
                }
            } else {
                if (!subtract) {
                    return new _Polynomial(*p2);
                } else {
                    _Polynomial* ppp = new _Polynomial (*p2);
                    checkPointer(ppp);
                    _Parameter *invC = ppp->theTerms->theCoeff;
                    for (long inv = 0; inv<ppp->theTerms->actTerms; inv++, invC++) {
                        (*invC)*=(-1.0);
                    }

                    return ppp;
                }
            }

        }
        if (p2->variableIndex.lLength == 0) {
            if (p2->theTerms->NumberOfTerms()) {
                _Constant coef2 (p2->theTerms->GetCoeff(0));
                return    Plus (&coef2,subtract);
            } else {
                if (!subtract) {
                    return new _Polynomial(*this);
                } else {
                    _Polynomial* ppp = new _Polynomial (*this);
                    checkPointer(ppp);
                    _Parameter *invC = ppp->theTerms->theCoeff;
                    for (long inv = 0; inv<ppp->theTerms->actTerms; inv++, invC++) {
                        (*invC)*=(-1.0);
                    }

                    return ppp;
                }
            }
        }

        p2->Convert2OperationForm();

        long nt2        = p2->theTerms->NumberOfTerms(),
             nt1      = theTerms->NumberOfTerms(),
             pos1      = 0,
             pos2        = 0,
             * term1,
             * term2;

        char c          = 0,
             advancing  = -1;

        _Parameter * coeff1 = theTerms->GetCoeff(),
                     *coeff2 = p2->theTerms->GetCoeff();

        _Polynomial* res;

        if (variableIndex.Equal(p2->variableIndex))
            // same variable arrays - proceed to the operation directly
        {
            res = new _Polynomial (variableIndex); // create a blank new result holder
            checkPointer        (res);
            ResetPolynomialCheck(res);
            while (1) { // stuff left to do
                if (advancing == 0) { // advancing in the 1st polynomial
                    pos1++;
                    if (pos1>=nt1) {
                        advancing = 2;
                        continue;
                    }
                    coeff1++;
                    term1 = theTerms->GetTerm (pos1);
                    c = res->theTerms->CompareTerms (term1, term2);
                    if (c<=0) {
                        res->theTerms->AddTerm (term1,*coeff1);
                        if (c<0) {
                            continue;
                        }
                    }
                    if (c>0) {
                        advancing = 1;
                        res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2);
                        continue;
                    }
                } else if (advancing == 1) { // advancing in the 2nd polynomial
                    pos2++;
                    if (pos2>=nt2) {
                        advancing = 3;
                        continue;
                    }
                    coeff2++;
                    term2 = p2->theTerms->GetTerm (pos2);
                    c = res->theTerms->CompareTerms (term2, term1);
                    if (c<=0) {
                        res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2);
                        if (c<0) {
                            continue;
                        }
                    }
                    if (c>0) {
                        advancing = 0;
                        res->theTerms->AddTerm (term1,*coeff1);
                        continue;
                    }
                } else if (advancing == 2) { // flush out the 2nd poly
                    while (pos2<nt2) {
                        res->theTerms->AddTerm (p2->theTerms->GetTerm(pos2),subtract?-*coeff2:*coeff2);
                        pos2++;
                        coeff2++;
                    }
                    break;
                } else if (advancing == 3) { // flush out the 2nd poly
                    while (pos1<nt1) {
                        res->theTerms->AddTerm (theTerms->GetTerm(pos1),*coeff1);
                        pos1++;
                        coeff1++;
                    }
                    break;
                } else if (advancing == -1) { // just starting
                    if (!nt1) { // first poly is empty!
                        advancing = 3;
                        continue;
                    }
                    if (!nt2) { // second poly is empty!
                        advancing = 2;
                        continue;
                    }
                    term1 = theTerms->GetTerm(0);
                    term2 = p2->theTerms->GetTerm(0);
                    c = res->theTerms->CompareTerms (term1, term2);
                    if (c<=0) { // begin with the first poly
                        res->theTerms->AddTerm (term1, *coeff1);
                        advancing = 0;
                        if (c) {
                            continue;
                        } else {
                            if (subtract) {
                                res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff2;
                            } else {
                                res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff2;
                            }
                            pos2++;
                            if (pos2==nt2) {
                                pos1++;
                                coeff1++;
                                advancing = 3;
                            } else {
                                term2 = p2->theTerms->GetTerm (pos2);
                                coeff2++;
                            }
                            res->CheckTerm();
                            continue;
                        }
                    } else {
                        res->theTerms->AddTerm (term2, subtract?-*coeff2:*coeff2);
                        advancing = 1;
                        continue;
                    }

                }

                // will only get here when there are like terms
                // see which one is being added
                if (advancing == 0) { // moving up in the first term
                    if (subtract) {
                        res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff2;
                    } else {
                        res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff2;
                    }
                    pos1++;
                    if (pos1>=nt1) {
                        advancing = 2;
                        pos2++;
                        term2 = p2->theTerms->GetTerm (pos2);
                        coeff2++;
                    } else {
                        term1 = theTerms->GetTerm (pos1);
                        coeff1++;
                        advancing = 1;
                    }

                } else {
                    //if (subtract)
                    //  res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff1;
                    //else
                    res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff1;
                    pos2++;
                    if (pos2>=nt2) {
                        advancing = 3;
                        pos1++;
                        term1 = theTerms->GetTerm (pos1);
                        coeff1++;
                    } else {
                        term2 = p2->theTerms->GetTerm (pos2);
                        coeff2++;
                        advancing = 0;
                    }
                }

                // check to see if the term is significant
                res->CheckTerm(); // by default will check the last term
            }


        } else {
            _SimpleList mergedVariables, merge1, merge2;
            mergedVariables.Merge(variableIndex, p2->variableIndex, &merge1, &merge2);
            if ((merge1.countitems()==0)||(merge2.countitems()==0))
                // one of the poly's variables are a superset for the other - treat accordingly
            {
                bool  firstBigger = merge1.countitems()==0; // is the first poly bigger than the second
                long  *reindexList = firstBigger?merge2.quickArrayAccess():merge1.quickArrayAccess(),
                       reindexLength = firstBigger?merge2.countitems():merge1.countitems();
                res = new _Polynomial (mergedVariables); // create a blank new result holder
                checkPointer(res);
                ResetPolynomialCheck(res);
                while (1) { // stuff left to do
                    if (advancing == 0) { // advancing in the 1st polynomial
                        pos1++;
                        if (pos1>=nt1) {
                            advancing = 2;
                            continue;
                        }
                        coeff1++;
                        term1 = theTerms->GetTerm (pos1);
                        c = firstBigger?res->theTerms->CompareTerms (term1, term2,reindexList,reindexLength):
                            -res->theTerms->CompareTerms (term2, term1,reindexList,reindexLength);
                        if (c<=0) {
                            if (firstBigger) {
                                res->theTerms->AddTerm (term1,*coeff1);
                            } else {
                                res->theTerms->AddTerm (term1,*coeff1,reindexList, reindexLength);
                            }
                            if (c<0) {
                                continue;
                            }
                        }
                        if (c>0) {
                            advancing = 1;
                            if (firstBigger) {
                                res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2,reindexList,reindexLength);
                            } else {
                                res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2);
                            }
                            continue;
                        }

                    } else if (advancing == 1) { // advancing in the 2nd polynomial
                        pos2++;
                        if (pos2>=nt2) {
                            advancing = 3;
                            continue;
                        }
                        coeff2++;
                        term2 = p2->theTerms->GetTerm (pos2);
                        c = firstBigger?-res->theTerms->CompareTerms (term1, term2,reindexList,reindexLength):
                            res->theTerms->CompareTerms (term2, term1,reindexList,reindexLength);
                        if (c<=0) {
                            if (firstBigger) {
                                res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2,reindexList, reindexLength);
                            } else {
                                res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2);
                            }
                            if (c<0) {
                                continue;
                            }
                        }
                        if (c>0) {
                            advancing = 0;
                            if (firstBigger) {
                                res->theTerms->AddTerm (term1,*coeff1);
                            } else {
                                res->theTerms->AddTerm (term1,*coeff1,reindexList, reindexLength);
                            }
                            continue;
                        }
                    } else if (advancing == 2) { // flush out the 2nd poly
                        if (firstBigger)
                            while (pos2<nt2) {
                                res->theTerms->AddTerm (p2->theTerms->GetTerm(pos2),subtract?-*coeff2:*coeff2,reindexList,reindexLength);
                                pos2++;
                                coeff2++;
                            }
                        else
                            while (pos2<nt2) {
                                res->theTerms->AddTerm (p2->theTerms->GetTerm(pos2),subtract?-*coeff2:*coeff2);
                                pos2++;
                                coeff2++;
                            }
                        break;
                    } else if (advancing == 3) { // flush out the 2nd poly
                        if (firstBigger)
                            while (pos1<nt1) {
                                res->theTerms->AddTerm (theTerms->GetTerm(pos1),*coeff1);
                                pos1++;
                                coeff1++;
                            }
                        else
                            while (pos1<nt1) {
                                res->theTerms->AddTerm (theTerms->GetTerm(pos1),*coeff1,reindexList,reindexLength);
                                pos1++;
                                coeff1++;
                            }
                        break;
                    } else if (advancing == -1) { // just starting
                        if (!nt1) { // first poly is empty!
                            advancing = 3;
                            continue;
                        }
                        if (!nt2) { // second poly is empty!
                            advancing = 2;
                            continue;
                        }
                        term1 = theTerms->GetTerm(0);
                        term2 = p2->theTerms->GetTerm(0);
                        c = firstBigger?res->theTerms->CompareTerms (term1, term2,reindexList,reindexLength):
                            -res->theTerms->CompareTerms (term2, term1,reindexList,reindexLength);
                        if (c<=0) { // begin with the first poly
                            if (firstBigger) {
                                res->theTerms->AddTerm (term1, *coeff1);
                            } else {
                                res->theTerms->AddTerm (term1, *coeff1,reindexList, reindexLength);
                            }
                            advancing = 0;
                            if (c) {
                                continue;
                            } else {
                                if (subtract) {
                                    res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff2;
                                } else {
                                    res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff2;
                                }
                                pos2++;
                                if (pos2>=nt2) {
                                    advancing = 3;
                                    coeff1++;
                                    pos1++;
                                } else {
                                    term2 = p2->theTerms->GetTerm (pos2);
                                    coeff2++;
                                }
                                res->CheckTerm();
                                continue;
                            }
                        } else {
                            if (firstBigger) {
                                res->theTerms->AddTerm (term2, subtract?-*coeff2:*coeff2,reindexList, reindexLength);
                            } else {
                                res->theTerms->AddTerm (term2, subtract?-*coeff2:*coeff2);
                            }
                            advancing = 1;
                            continue;
                        }

                    }

                    // will only get here when there are like terms
                    // see which one is being added
                    if (advancing == 0) { // moving up in the first term
                        if (subtract) {
                            res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff2;
                        } else {
                            res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff2;
                        }
                        pos1++;
                        if (pos1>=nt1) {
                            advancing = 2;
                            pos2++;
                            term2 = p2->theTerms->GetTerm (pos2);
                            coeff2++;
                        } else {
                            term1 = theTerms->GetTerm (pos1);
                            coeff1++;
                            advancing = 1;
                        }

                    } else {
                        //if (subtract)
                        //  res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff1;
                        //else
                        res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff1;
                        pos2++;
                        if (pos2>=nt2) {
                            advancing = 3;
                            pos1++;
                            coeff1++;
                            term1 = theTerms->GetTerm (pos1);
                        } else {
                            term2 = p2->theTerms->GetTerm (pos2);
                            coeff2++;
                            advancing = 0;
                        }
                    }

                    // check to see if the term is significant
                    res->CheckTerm(); // by default will check the last term
                }
            } else
                // both variable indices must be reindexed
            {
                long  *ri1 = merge1.quickArrayAccess(),
                       *ri2 = merge2.quickArrayAccess(),
                        rl1 = merge1.countitems(),
                        rl2 = merge2.countitems() ;

                res = new _Polynomial (mergedVariables); // create a blank new result holder
                checkPointer(res);
                ResetPolynomialCheck(res);

                while (1) { // stuff left to do
                    if (advancing == 0) { // advancing in the 1st polynomial
                        pos1++;
                        if (pos1>=nt1) {
                            advancing = 2;
                            continue;
                        }
                        coeff1++;
                        term1 = theTerms->GetTerm (pos1);
                        c = res->theTerms->CompareTerms (term1, term2,ri1,ri2,rl1,rl2);
                        if (c<=0) {
                            res->theTerms->AddTerm (term1,*coeff1,ri1, rl1);
                            if (c<0) {
                                continue;
                            }
                        }
                        if (c>0) {
                            advancing = 1;
                            res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2,ri2,rl2);
                            continue;
                        }

                    } else if (advancing == 1) { // advancing in the 2nd polynomial
                        pos2++;
                        if (pos2>=nt2) {
                            advancing = 3;
                            continue;
                        }
                        coeff2++;
                        term2 = p2->theTerms->GetTerm (pos2);
                        c = res->theTerms->CompareTerms (term2, term1,ri2,ri1,rl2,rl1);
                        if (c<=0) {
                            res->theTerms->AddTerm (term2,subtract?-*coeff2:*coeff2,ri2, rl2);
                            if (c<0) {
                                continue;
                            }
                        }
                        if (c>0) {
                            advancing = 0;
                            res->theTerms->AddTerm (term1,*coeff1,ri1, rl1);
                            continue;
                        }
                    } else if (advancing == 2) { // flush out the 2nd poly
                        while (pos2<nt2) {
                            res->theTerms->AddTerm (p2->theTerms->GetTerm(pos2),subtract?-*coeff2:*coeff2,ri2,rl2);
                            pos2++;
                            coeff2++;
                        }
                        break;
                    } else if (advancing == 3) { // flush out the 2nd poly
                        while (pos1<nt1) {
                            res->theTerms->AddTerm (theTerms->GetTerm(pos1),*coeff1,ri1,rl1);
                            pos1++;
                            coeff1++;
                        }
                        break;
                    } else if (advancing == -1) { // just starting
                        if (!nt1) { // first poly is empty!
                            advancing = 3;
                            continue;
                        }
                        if (!nt2) { // second poly is empty!
                            advancing = 2;
                            continue;
                        }
                        term1 = theTerms->GetTerm(0);
                        term2 = p2->theTerms->GetTerm(0);
                        c = res->theTerms->CompareTerms (term1, term2,ri1,ri2,rl1,rl2);
                        if (c<=0) { // begin with the first poly
                            res->theTerms->AddTerm (term1, *coeff1,ri1, rl1);
                            advancing = 0;
                            if (c) {
                                continue;
                            } else {
                                if (subtract) {
                                    res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff2;
                                } else {
                                    res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff2;
                                }
                                pos2++;
                                if (pos2==nt2) {
                                    pos1++;
                                    coeff1++;
                                    advancing = 3;
                                } else {
                                    term2 = p2->theTerms->GetTerm (pos2);
                                    coeff2++;
                                }
                                res->CheckTerm();
                                continue;
                            }
                        } else {
                            res->theTerms->AddTerm (term2, subtract?-*coeff2:*coeff2,ri2, rl2);
                            advancing = 1;
                            continue;
                        }

                    }

                    // will only get here when there are like terms
                    // see which one is being added
                    if (advancing == 0) { // moving up in the first term
                        if (subtract) {
                            res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff2;
                        } else {
                            res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff2;
                        }
                        pos1++;
                        if (pos1>=nt1) {
                            advancing = 2;
                            pos2++;
                            term2 = p2->theTerms->GetTerm (pos2);
                            coeff2++;
                        } else {
                            term1 = theTerms->GetTerm (pos1);
                            coeff1++;
                            advancing = 1;
                        }

                    } else {
                        //if (subtract)
                        //  res->theTerms->GetCoeff (res->theTerms->actTerms-1)-= *coeff1;
                        //else
                        res->theTerms->GetCoeff (res->theTerms->actTerms-1)+= *coeff1;
                        pos2++;
                        if (pos2>=nt2) {
                            advancing = 3;
                            pos1++;
                            term1 = theTerms->GetTerm (pos1);
                            coeff1++;
                        } else {
                            term2 = p2->theTerms->GetTerm (pos2);
                            coeff2++;
                            advancing = 0;
                        }
                    }

                    // check to see if the term is significant
                    res->CheckTerm(); // by default will check the last term
                }
            }
        }

        if (!res->theTerms->checkMe()) {
            //BufferToConsole (_String((_String*)toStr()));
            //NLToConsole();
            //BufferToConsole (_String((_String*)m->toStr()));
            //NLToConsole();
            //BufferToConsole (_String((_String*)res->toStr()));
            //NLToConsole();
            return nil;
        }
//      res->theTerms->ChopTerms();
        if (res->theTerms->GetNoTerms()==0) {
            DeleteObject (res);
            return new _Polynomial (0.0);
        }
        return res;
    }

    _String errMsg ("An incompatible operand was supplied to polynomial addition/subtraction");
    FlagError (&errMsg);
    return nil;
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Sub (_MathObject* m)
{
    return Plus (m,true);
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Minus (void)
{
    _Constant min (-1.0);
    return Mult (&min);
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Mult (_MathObject* m)
{
    long objectT = m->ObjectClass();

    if (objectT==NUMBER) { // a number or a monomial
        Convert2OperationForm();
        _Parameter nb = ((_Constant*)m)->Value();
        if (nb==0.0) {
            return new _Polynomial;
        }

        _Polynomial* result = new _Polynomial (*this);
        checkPointer(result);
        for (long i = theTerms->NumberOfTerms()-1; i>=0; i--) {
            result->theTerms->GetCoeff(i)*=nb;
        }
        return result;
    }

    if (objectT == POLYNOMIAL) { // another polynomial
        /*      Convert2OperationForm();
                _Polynomial* p2 = (_Polynomial*)m, *res;
                p2->Convert2OperationForm();
                if (!variableIndex.countitems())
                {
                    if (theTerms->NumberOfTerms())
                    {
                        _Constant coef1 (theTerms->GetCoeff(0));
                        return p2->Mult (&coef1);
                    }
                    else
                        return new _Polynomial();
                }
                if (!p2->variableIndex.countitems())
                {
                    if (p2->theTerms->NumberOfTerms())
                    {
                        _Constant coef2 (p2->theTerms->GetCoeff(0));
                        return Mult (&coef2);
                    }
                    else
                        return new _Polynomial();
                }

                long nt1 = theTerms->NumberOfTerms ()-1, nt2 = p2->theTerms->NumberOfTerms ()-1,
                     csp = 0, cwp = 0;

                long *positionArray = (long*)MemAllocate ((nt2+1)*sizeof(long)), temp1, temp2, resVars;
                char compRes, compRes2, compRes3, status = 0;
                if ((nt1<0) || (nt2<0)) return new _Polynomial();
                memset (positionArray,0,(nt2+1)*sizeof(long));
                long* terms1, *terms2, *holder1, *holder2, *holder3;

                //set up proper variable/term agreement
                if (variableIndex.Equal(p2->variableIndex))
                // same variables - the simplest case
                {
                    terms1 = theTerms->GetTerm(0);
                    terms2 = p2->theTerms->GetTerm(0);
                    resVars = variableIndex.countitems();
                    res = new _Polynomial (variableIndex);
                }
                else
                {
                    _SimpleList mergedVars, merge1, merge2;
                    mergedVars.Merge (variableIndex,p2->variableIndex, &merge1, &merge2);
                    if (!merge1.countitems()) // 1 is a superset of 2
                    {
                        res = new _Polynomial (variableIndex);
                        temp2 = p2->variableIndex.countitems();
                        resVars = temp1 = variableIndex.countitems();
                        terms1 = theTerms->GetTerm(0);
                        holder2 = terms2 = (long*)MemAllocate (nt2*sizeof(long)*temp1);
                        holder3 = p2->theTerms->GetTerm(0);
                        holder1 = merge2.quickArrayAccess();
                        for (;csp<nt2;csp++,holder2+=temp1, holder3+= temp2)
                        {
                            _PolynomialData::RearrangeTerm (holder2, holder3, holder1, temp2);
                        }
                    }
                    else
                        if (!merge2.countitems()) // 2 is a superset of 1
                        {
                            res = new _Polynomial (p2->variableIndex);
                            resVars = temp2 = p2->variableIndex.countitems();
                            temp1 = variableIndex.countitems();
                            terms2 = p2->theTerms->GetTerm(0);
                            holder1 = terms1 = (long*)MemAllocate (nt1*sizeof(long)*temp2);
                            holder3 = theTerms->GetTerm(0);
                            holder2 = merge1.quickArrayAccess();
                            for (;csp<nt1;csp++,holder1+=temp2, holder3+= temp1)
                            {
                                _PolynomialData::RearrangeTerm (holder1, holder3, holder2, temp1);
                            }
                        }
                        else // no inclusions
                        {
                            res = new _Polynomial (mergedVars);
                            resVars = mergedVars.countitems();

                            temp2 = p2->variableIndex.countitems();
                            temp1 = variableIndex.countitems();
                            holder2 = terms2 = (long*)MemAllocate (nt2*sizeof(long)*resVars);
                            holder3 = p2->theTerms->GetTerm(0);
                            holder1 = merge2.quickArrayAccess();
                            for (;csp<nt2;csp++,holder2+=resVars, holder3+= temp2)
                            {
                                _PolynomialData::RearrangeTerm (holder2, holder3, holder1, temp2);
                            }
                            holder1 = terms1 = (long*)MemAllocate (nt1*sizeof(long)*resVars);
                            holder3 = theTerms->GetTerm(0);
                            holder2 = merge1.quickArrayAccess();
                            for (;csp<nt1;csp++,holder1+=resVars, holder3+= temp1)
                            {
                                _PolynomialData::RearrangeTerm (holder1, holder3, holder2, temp1);
                            }
                        }

                    } // done arranging variables

                    ResetPolynomialCheck(res);
                    holder1 = (long*)MemAllocate (resVars*sizeof(long));
                    holder2 = (long*)MemAllocate (resVars*sizeof(long));
                    holder3 = (long*)MemAllocate (resVars*sizeof(long));

                    _Parameter *tempCoeff = (_Parameter*)MemAllocate (theTerms->NumberOfTerms()*sizeof(_Parameter));
                    memcpy (tempCoeff, theTerms->GetCoeff(), theTerms->NumberOfTerms()*sizeof(_Parameter));


                    res->theTerms->MultiplyTerms    (holder1, terms1, terms2); // add the first term

                    res->theTerms->AddTerm (holder1, *tempCoeff * *p2->theTerms->GetCoeff());


                    while (1)
                    {
                        if ((status==0)&&(cwp==csp))
                        {
                            temp1 = positionArray[csp]+1;
                            if (temp1>nt2)
                            {
                                csp++;
                                cwp = csp;
                                if (csp>nt1) break;
                                if (csp==nt1)
                                {
                                    long i = positionArray[csp]+1;
                                    _Parameter coe = theTerms->GetCoeff (csp);
                                    while (i<=nt2)
                                    {
                                        res->theTerms->MultiplyTerms (holder1, theTerms->GetTerm(csp), p2->theTerms->GetTerm(i));
                                        res->theTerms->AddTerm (holder1, coe*p2->theTerms->GetCoeff (i));
                                        i++;
                                    }
                                    break;
                                }
                                continue;
                            }
                            compRes = 0;
                        }
                        else
                        {
                            if (positionArray[cwp]+1>nt2)
                            {
                                csp = cwp+1;
                                if (csp==nt1)
                                {
                                    long i = positionArray[csp]+1;
                                    _Parameter coe = theTerms->GetCoeff (csp);
                                    while (i<=nt2)
                                    {
                                        res->theTerms->MultiplyTerms (holder1, theTerms->GetTerm(csp), p2->theTerms->GetTerm(i));
                                        res->theTerms->AddTerm (holder1, coe*p2->theTerms->GetCoeff (i));
                                        i++;
                                    }
                                    break;
                                }
                                cwp = csp;
                                status = 0;
                                continue;
                            }
                        }
                        res->theTerms->MultiplyTerms (holder1, terms1+resVars*cwp, terms2+resVars*(positionArray[cwp]+1));

                        if ((compRes>=0)||status)
                            res->theTerms->MultiplyTerms    (holder2, terms1+resVars*(cwp+1), terms2+resVars*positionArray[cwp+1]);

                        if (status==0)
                        {
                            compRes = res->theTerms->CompareTerms (holder1, holder2);
                        }
                        else
                            compRes = 1;

                        if (compRes<0)
                        {
                            if (status==1)
                            {
                                cwp --;
                                if (cwp==csp) status = 0;
                                continue;
                            }
                            res->theTerms->AddTerm(holder1,tempCoeff[cwp]*p2->theTerms->GetCoeff (positionArray[cwp]));
                            if (cwp>csp)
                            {
                                tempCoeff[cwp]=theTerms->GetCoeff(cwp);
                            }
                            res->CheckTerm();
                            positionArray[cwp]++;
                            status = 0;
                            continue;
                        }
                        if (compRes==0)
                        {
                            positionArray[cwp]++;
                            tempCoeff[cwp+1]+= tempCoeff[cwp]*p2->theTerms->GetCoeff (positionArray[cwp])/
                                                         p2->theTerms->GetCoeff (positionArray[cwp+1]);
                            status = 0;
                            continue;
                        }

                        if ((cwp-csp>1)||(cwp==nt1-1))
                        {
                            compRes3 = 1;
                            for (long k=csp; k<cwp-1;k++)
                            {
                                res->theTerms->MultiplyTerms    (holder3, terms1+resVars*(k), terms2+resVars*(positionArray[k]+1));
                                compRes3 = res->theTerms->CompareTerms (holder3, holder2);
                                if (!compRes3)
                                {
                                    tempCoeff[cwp+1]+= tempCoeff[k]*p2->theTerms->GetCoeff (positionArray[k])/
                                                         p2->theTerms->GetCoeff (positionArray[k+1]);
                                    positionArray[csp]++;
                                    if (k==csp)
                                        if (positionArray[csp]>=nt2)
                                            csp++;
                                    break;
                                }
                                if (compRes3<0) break;
                            }

                        }
                        res->theTerms->AddTerm(holder2,tempCoeff[cwp+1]*p2->theTerms->GetCoeff (positionArray[cwp+1]));
                        tempCoeff[cwp+1]=theTerms->GetCoeff(cwp+1);
                        res->CheckTerm();

                        temp1 = (positionArray[cwp+1]+1);
                        res->theTerms->MultiplyTerms    (holder2, terms1+resVars*(cwp+1), terms2+resVars*(temp1));
                        compRes2 = res->theTerms->CompareTerms (holder1, holder2);
                        if (compRes2<=0)
                        {
                            if (!compRes2)
                            {
                                positionArray[cwp]++;
                                tempCoeff[cwp+1]+= tempCoeff[cwp]*p2->theTerms->GetCoeff (positionArray[cwp])/
                                                         p2->theTerms->GetCoeff (temp1);
                                if (cwp>csp)
                                {
                                    tempCoeff[cwp]=theTerms->GetCoeff(cwp);
                                }
                            }

                            if (cwp<=nt1-2)
                            {
                                res->theTerms->MultiplyTerms    (holder3,  terms1+resVars*(cwp+2),  terms2+resVars*(positionArray[cwp+2]));
                                compRes3 = res->theTerms->CompareTerms (holder2, holder3);
                                compRes = res->theTerms->CompareTerms (holder1, holder3);
                                if (!compRes)
                                {
                                    positionArray[csp]++;
                                    tempCoeff[cwp+2]+= tempCoeff[csp]*p2->theTerms->GetCoeff (positionArray[csp])/
                                                             p2->theTerms->GetCoeff (positionArray[cwp+2]);
                                    if (positionArray[csp]>nt2)
                                        csp++;
                                }
                            }
                            else
                            {
                                compRes3 = -1;
                            }

                            if (compRes3<0)
                            {
                                positionArray[cwp+1]=temp1;
                                res->theTerms->AddTerm (holder1, tempCoeff[cwp+1]*p2->theTerms->GetCoeff (positionArray[cwp+1]));
                                tempCoeff[cwp+1]=theTerms->GetCoeff(cwp+1);
                                res->CheckTerm();
                                compRes = 1;
                                if (status==1)
                                {
                                    if (compRes2==0)
                                    {
                                        positionArray[cwp+1]++;
                                        if (positionArray[cwp+1]>nt2)
                                        {
                                            csp = cwp+1;
                                            if (cwp<csp)
                                            {
                                                cwp = csp;
                                                status = 0;
                                            }
                                        }
                                        else
                                            if (positionArray[cwp+1]==nt2)
                                            {
                                                positionArray[cwp+1]--;
                                            }
                                    }
                                }
                                continue;
                            }

                            if (!compRes3)
                            {
                                tempCoeff[cwp+2]+= tempCoeff[cwp+1]*p2->theTerms->GetCoeff (positionArray[cwp+1])/
                                                         p2->theTerms->GetCoeff (positionArray[cwp+2]);
                                tempCoeff[cwp+1]=theTerms->GetCoeff(cwp+1);
                                status = 0;
                                continue;
                            }

                            cwp++;
                            status = 1;
                            continue;
                        }
                        else
                        {
                            if (cwp<=nt1-2)
                            {
                                res->theTerms->MultiplyTerms    (holder3,  terms1+resVars*(cwp+2),  terms2+resVars*(positionArray[cwp+2]));
                                compRes3 = res->theTerms->CompareTerms (holder2, holder3);
                            }
                            else
                            {
                                compRes3 = -1;
                            }
                            if (compRes3 < 0)
                            {
                                positionArray[cwp+1]=temp1;
                                tempCoeff[cwp+1]=theTerms->GetCoeff(cwp+1);
                                status = 1;
                                continue;
                            }

                            cwp++;
                            if (compRes3 == 0)
                            {
                                tempCoeff[cwp+1]+= tempCoeff[cwp]*p2->theTerms->GetCoeff (positionArray[cwp])/p2->theTerms->GetCoeff (positionArray[cwp+1]);
                                tempCoeff[cwp]=theTerms->GetCoeff(cwp);
                                continue;
                            }

                            status = 1;
                            continue;

                        }

                    }

                    free (positionArray);
                    free (holder1);
                    free (holder2);
                    free (holder3);
                    free (tempCoeff);
                    if (terms1!=theTerms->GetTerm(0))
                        free (terms1);
                    if (terms2!=p2->theTerms->GetTerm(0))
                        free (terms2);
                    res->theTerms->checkMe();
                    return res;*/


        Convert2OperationForm();
        _Polynomial* p2 = (_Polynomial*)m;
        p2->Convert2OperationForm();
        if (!variableIndex.countitems()) {
            if (theTerms->NumberOfTerms()) {
                _Constant coef1 (theTerms->GetCoeff(0));
                return p2->Mult (&coef1);
            } else {
                return new _Polynomial();
            }

        }
        if (!p2->variableIndex.countitems()) {
            if (p2->theTerms->NumberOfTerms()) {
                _Constant coef2 (p2->theTerms->GetCoeff(0));
                return Mult (&coef2);
            } else {
                return new _Polynomial();
            }
        }
        long nt2 = p2->theTerms->NumberOfTerms(), nt1 = theTerms->NumberOfTerms(),
             ref = (variableIndex.countitems()+p2->variableIndex.countitems())*maximumPolyTermsPerVariable,
             *idx1;
        _Polynomial* res;
        _SimpleList  reIndex, terms1, terms2;

        if (nt1*nt2>ref)
            // too many terms expected - chop b4 multiplying
        {
            _SimpleList  index;
            long i;
            _Parameter  logTop = log(topPolyCap), *theCoeff = theTerms->GetCoeff();
            for (i = 0; i<nt1; i++,theCoeff++) {
                index<<i;
                terms1<<(long)(log(fabs(*theCoeff))+logTop*theTerms->SumOfPowers(i));
            }
            theCoeff = p2->theTerms->GetCoeff();
            for (i = 0; i<nt2; i++,theCoeff++) {
                index<<i+nt1;
                terms1<<(long)(log(fabs(*theCoeff))+logTop*p2->theTerms->SumOfPowers(i));
            }
            SortLists (&terms1,&index);
            terms1.Clear();
            terms2.Clear();
            idx1 = index.quickArrayAccess()+nt1+nt2-1;
            long temp1 = sqrt((_Parameter)ref*nt1/nt2)+1, temp2 = sqrt((_Parameter)ref*nt2/nt1)+1;
            if (temp1>nt1) {
                temp1 = nt1;
            }
            if (temp2>nt2) {
                temp2 = nt2;
            }
            for (i = 0; i<ref; i++, idx1--) {
                if (*idx1>=nt1) {
                    terms2<<*idx1;
                } else {
                    terms1<<*idx1;
                }
            }
            if (terms1.countitems()==0) {
                terms1<<0;
            }
            if (terms2.countitems()==0) {
                terms2<<0;
            }

        }

        if (!terms1.countitems()) {
            if (variableIndex.Equal(p2->variableIndex)) { // same variable arrays - proceed to the operation directly
                res = new _Polynomial (variableIndex);
                checkPointer(res);
                ResetPolynomialCheck(res);
                long* newTerm = new long[res->variableIndex.countitems()],*term1,f;
                checkPointer(newTerm);
                for (long i=0; i<nt1; i++) {
                    ref = 0;
                    _Parameter c1 = theTerms->GetCoeff(i);
                    term1 = theTerms->GetTerm(i);
                    for (long j=0; j<nt2; j++) {
                        res->theTerms->MultiplyTerms (newTerm, term1, p2->theTerms->GetTerm(j));
                        f = res->theTerms->FindTerm (newTerm,reIndex.quickArrayAccess(), ref);
                        if (f>=0) { // existing term - add or sub
                            res->theTerms->GetCoeff(reIndex(f))+=c1*p2->theTerms->GetCoeff(j);
                            ref = f+1;
                        } else { // new term - insert
                            reIndex.InsertElement ((BaseRef)res->theTerms->actTerms,-f-2, false,false);
                            res->theTerms->AddTerm (newTerm, c1*p2->theTerms->GetCoeff(j));
                            ref = -f-1;
                        }
                    }
                }

                delete newTerm;
            } else {
                _SimpleList merge1, merge2, joint;
                joint.Merge (variableIndex, p2->variableIndex,&merge1,&merge2);
                res = new _Polynomial (joint);
                checkPointer(res);
                ResetPolynomialCheck(res);
                long  f,nv = res->variableIndex.countitems(),* scratch1 = new long[nv],
                * scratch3,k,* scratch2;
                checkPointer(scratch1);
                bool do1 = merge1.countitems()!=0, do2 = merge2.countitems()!=0;
                if (do1) {
                    scratch2 = new long[nv];
                    checkPointer(scratch2);
                }
                if (do2) {
                    scratch3 = new long[nv];
                    checkPointer(scratch3);
                }
                for (long i=0; i<nt1; i++) {
                    _Parameter c1 = theTerms->GetCoeff(i);
                    ref = 0;
                    if (do1) {
                        for (k=0; k<nv; k++,scratch2++) {
                            *scratch2 = 0;
                        }
                        scratch2 -= nv;
                        theTerms->RearrangeTerm (scratch2, theTerms->GetTerm(i), merge1.quickArrayAccess(), merge1.countitems());
                    } else {
                        scratch2 = theTerms->GetTerm(i);
                    }
                    for (long j=0; j<nt2; j++) {
                        if (do2) {
                            for (k=0; k<nv; k++,scratch3++) {
                                *scratch3 = 0;
                            }
                            scratch3 -= nv;
                            theTerms->RearrangeTerm (scratch3, p2->theTerms->GetTerm(j), merge2.quickArrayAccess(),merge2.countitems());
                        } else {
                            scratch3 = p2->theTerms->GetTerm(j);
                        }
                        res->theTerms->MultiplyTerms (scratch1,scratch2,scratch3);
                        f = res->theTerms->FindTerm (scratch1,reIndex.quickArrayAccess(), ref);
                        if (f>=0) { // existing term - add or sub
                            res->theTerms->GetCoeff(reIndex(f))+=c1*p2->theTerms->GetCoeff(j);
                            ref = f+1;
                        } else { // new term - insert
                            reIndex.InsertElement ((BaseRef)res->theTerms->actTerms,-f-2, false,false);
                            res->theTerms->AddTerm (scratch1, c1*p2->theTerms->GetCoeff(j));
                            ref = -f-1;
                        }
                    }
                }
                delete [] scratch1;
                if (do1) {
                    delete [] scratch2;
                }
                if (do2) {
                    delete [] scratch3;
                }
            }

        } else {
            terms1.Sort();
            terms2.Sort();
            long* onOff = (long*)MemAllocate (nt2*sizeof(long));
            onOff+=(nt2-1);
            for (long i=nt2-1; i>=0; i--,onOff--)
                if (terms2.BinaryFind(i)>=0) {
                    *onOff = 1;
                } else {
                    *onOff = 0;
                }
            onOff++;
            if (variableIndex.Equal(p2->variableIndex)) { // same variable arrays - proceed to the operation directly
                res = new _Polynomial (variableIndex);
                checkPointer(res);
                ResetPolynomialCheck(res);
                long* newTerm = new long[res->variableIndex.countitems()],*term1,f;
                checkPointer(newTerm);
                for (long i=0; i<nt1; i++) {
                    if (terms1.BinaryFind(i)<0) {
                        continue;
                    }
                    ref = 0;
                    _Parameter c1 = theTerms->GetCoeff(i);
                    term1 = theTerms->GetTerm(i);
                    for (long j=0; j<nt2; j++) {
                        if (onOff[j]==0) {
                            continue;
                        }
                        res->theTerms->MultiplyTerms (newTerm, term1, p2->theTerms->GetTerm(j));
                        f = res->theTerms->FindTerm (newTerm,reIndex.quickArrayAccess(), ref);
                        if (f>=0) { // existing term - add or sub
                            res->theTerms->GetCoeff(reIndex(f))+=c1*p2->theTerms->GetCoeff(j);
                            ref = f+1;
                        } else { // new term - insert
                            reIndex.InsertElement ((BaseRef)res->theTerms->actTerms,-f-2, false,false);
                            res->theTerms->AddTerm (newTerm, c1*p2->theTerms->GetCoeff(j));
                            ref = -f-1;
                        }
                    }
                }

                delete []newTerm;
            } else {
                _SimpleList merge1, merge2, joint;
                joint.Merge (variableIndex, p2->variableIndex,&merge1,&merge2);
                res = new _Polynomial (joint);
                checkPointer(res);
                ResetPolynomialCheck(res);
                long  f,nv = res->variableIndex.countitems(),* scratch1 = new long[nv],
                * scratch3,k,* scratch2;
                checkPointer(scratch1);
                bool do1 = merge1.countitems()!=0, do2 = merge2.countitems()!=0;
                if (do1) {
                    scratch2 = new long[nv];
                    checkPointer(scratch2);
                }
                if (do2) {
                    scratch3 = new long[nv];
                    checkPointer(scratch3);
                }
                for (long i=0; i<nt1; i++) {
                    _Parameter c1 = theTerms->GetCoeff(i);
                    ref = 0;
                    if (do1) {
                        for (k=0; k<nv; k++,scratch2++) {
                            *scratch2 = 0;
                        }
                        scratch2 -= nv;
                        theTerms->RearrangeTerm (scratch2, theTerms->GetTerm(i), merge1.quickArrayAccess(), merge1.countitems());
                    } else {
                        scratch2 = theTerms->GetTerm(i);
                    }
                    for (long j=0; j<nt2; j++) {
                        if (do2) {
                            for (k=0; k<nv; k++,scratch3++) {
                                *scratch3 = 0;
                            }
                            scratch3 -= nv;
                            theTerms->RearrangeTerm (scratch3, p2->theTerms->GetTerm(j), merge2.quickArrayAccess(),merge2.countitems());
                        } else {
                            scratch3 = p2->theTerms->GetTerm(j);
                        }
                        res->theTerms->MultiplyTerms (scratch1,scratch2,scratch3);
                        f = res->theTerms->FindTerm (scratch1,reIndex.quickArrayAccess(), ref);
                        if (f>=0) { // existing term - add or sub
                            res->theTerms->GetCoeff(reIndex(f))+=c1*p2->theTerms->GetCoeff(j);
                            ref = f+1;
                        } else { // new term - insert
                            reIndex.InsertElement ((BaseRef)res->theTerms->actTerms,-f-2, false,false);
                            res->theTerms->AddTerm (scratch1, c1*p2->theTerms->GetCoeff(j));
                            ref = -f-1;
                        }
                    }
                }
                delete []scratch1;
                if (do1) {
                    delete []scratch2;
                }
                if (do2) {
                    delete []scratch3;
                }
            }
            free (onOff);
        }

        res->theTerms->ResortTerms(reIndex.quickArrayAccess());
        res->theTerms->checkMe();
//      res->theTerms->ChopTerms();
        return res;

    }

    _String errMsg ("An incompatible operand was supplied to polynomial multiplication");
    FlagError (&errMsg);
    return nil;
}


//__________________________________________________________________________________
_MathObject* _Polynomial::Raise (_MathObject* m)
{
    long objectT = m->ObjectClass();
    bool del = false;
    if (objectT==POLYNOMIAL) {
        m = ((_Polynomial*)m)->IsANumber();
        if (!m) {
            return nil;
        }
        del = true;
        objectT = m->ObjectClass();
    }
    if (objectT==1) { // a number
        Convert2OperationForm();
        _Polynomial* result;
        if (theTerms->NumberOfTerms()==1) { // just a monomial
            long power = (long)m->Value();
            result = new _Polynomial (*this);
            checkPointer(result);
            result->theTerms->RaiseTerm (result->theTerms->GetTerm(0),power);
            result->theTerms->GetCoeff(0)=_PolynomialData::BinaryRaise(result->theTerms->GetCoeff(0),power);

        } else { // binary raise
            result = new _Polynomial(1.0);
            checkPointer(result);
            _Polynomial* oldR;
            char bits[sizeof(long)*8], nLength = 0;
            long pwr = m->Value();
            while (pwr) {
                bits[nLength]=pwr%2;
                pwr/=2;
                nLength++;
            }
            while (nLength) {
                oldR = result;
                nLength--;
                if (bits[nLength]) {
                    result=(_Polynomial*)result->Mult(this);
                    DeleteObject( oldR);
                }
                oldR = result;
                if (nLength) {
                    result=(_Polynomial*)result->Mult(result);
                    DeleteObject( oldR);
                }
            }
        }
        if (del) {
            DeleteObject(m);
        }
        return result;
    }
    _String errMsg ("An incompatible operand was supplied to polynomial raise to power");
    FlagError (&errMsg);
    return nil;
}
//__________________________________________________________________________________
void    _Polynomial::Convert2ComputationForm(_SimpleList* c1, _SimpleList* c2, _SimpleList* termsToInclude)
{
    // compList1 and two one will contain pairs of numbers to be interpreted as follows
    // for pair (m,n)
    // m>=0, n>=0, m!=last var
    // m<0, n<0 - mult by -n-th power of var -m-1, do not add yet
    // m<0, n>0 - mult by -n-th power of var -m-1, do add
    // m>=0, n<0 - reset all vars after m and set the power of m-the var to -n, do not compute yet
    // m>=0, n>0 - reset all vars after m and set the power of m-the var to -n, do add
    // last var, n>0 - add n successive powers of last var
    if ((theTerms->NumberOfTerms()!=0)&&(compList1.countitems()==0)) { // stuff to do
        _SimpleList* cL1, *cL2, ti;

        long i = 0,n=variableIndex.countitems()-1,p,l, s=0,j,*vi;
        if (c1&&c2&&termsToInclude) {
            cL1 = c1;
            cL2 = c2;
            ti.Duplicate (termsToInclude);
        } else {
            cL1 = &compList1;
            cL2 = &compList2;
            for (i=0; i<theTerms->actTerms; i++) {
                ti<<i;
            }
            i = 0;
        }
        l = ti.countitems();
        vi = ti.quickArrayAccess();
        (*cL1).Clear();
        (*cL2).Clear();

        if (!theTerms->IsFirstANumber()) {
            long* fst = theTerms->GetTerm(vi[0]);
            p = fst[n];
            // set up the first term
            (*cL1)<<-n-1;
            (*cL2)<<-p;

            for (j=variableIndex.countitems()-2; j>=0; j--) {
                p = fst[j];
                if (p) {
                    (*cL1)<<-j-1;
                    (*cL2)<<-p;
                }
            }
            (*cL2)[(*cL2).countitems()-1]*=-1;
            if ((*cL2).countitems()>1) {
                if (fst[n]==0) {
                    (*cL2).Delete(0);
                    (*cL1).Delete(0);
                }
            }
        } else {
            (*cL1)<<n;
            (*cL2)<<0;
        }
        i++;
        long* powerDiff = new long[n+2];
        checkPointer(powerDiff);
        for (; i<l; i++) { // even more stuff to do

            long* cM = theTerms->GetTerm(vi[i]),
                  * prevM = theTerms->GetTerm(vi[i-1]);
            long k = 0, ch = -1;
            bool    reset = false;
            for (j=0; j<n; j++,cM++,prevM++) {
                powerDiff[j]=*cM-*prevM;
                if (powerDiff[j]) {
                    if (ch<0) {
                        ch = j;
                    }
                    k--;
                    if (!reset) {
                        reset = powerDiff[j]<0;
                    }
                }
            }

            powerDiff[j]=*cM-*prevM;

            if (!reset) {
                reset = powerDiff[n]<0;
            }
            if (!k) {
                k = powerDiff[n];
            } else if (powerDiff[n]) {
                k--;
            }
            // analyze the difference
            if (k==1) {
                s++;
                continue;
            } else {
                if (s>0) {
                    (*cL1)<<n;
                    (*cL2)<<s;
                    s=0;
                }
                if (k>0) {
                    (*cL1)<<n;
                    (*cL2)<<-k;
                    continue;
                }
                if (k<0) {
                    // difference in more than one power
                    if (k==-1) { // change in exactly one variable
                        (*cL1)<<-ch-1;
                        (*cL2)<<powerDiff[ch];

                        continue;
                    }
                    // change in more than one variable
                    (*cL1)<<(reset?ch:-ch-1);
                    (*cL2)<<-powerDiff[ch];
                    long c = ch+1;
                    prevM = theTerms->GetTerm(vi[i-1])+c;
                    for (; c<n; c++,prevM++) {
                        if (powerDiff[c]>0) {
                            (*cL1)<<-c-1;
                            (*cL2)<<(reset?-(*prevM+powerDiff[c]):-powerDiff[c]);
                        } else if (powerDiff[c]<0) {
                            long tp = -(*prevM+powerDiff[c]);
                            if (tp) {
                                (*cL1)<<-c-1;
                                (*cL2)<<tp;
                            }
                        } else if (reset) {
                            if (*prevM) {
                                (*cL1)<<-c-1;
                                (*cL2)<<-*prevM;
                            }
                        }

                    }
                    if (powerDiff[n]>0) {
                        (*cL1)<<-n-1;
                        (*cL2)<<(reset?(-(*prevM+powerDiff[n])):-powerDiff[n]);
                    } else if (powerDiff[n]<0) {
                        long tp = -(*prevM+powerDiff[n]);
                        if (tp) {
                            (*cL1)<<-n-1;
                            (*cL2)<<tp;
                        }
                    } else if (reset) {
                        if (*prevM) {
                            (*cL1)<<-n-1;
                            (*cL2)<<-*prevM;

                        }
                    }
                    (*cL2)[(*cL2).countitems()-1]*=-1;

                }
            }
        }
        if (s>0) {
            (*cL1)<<n;
            (*cL2)<<s;
        }

        delete (powerDiff);
        if (!(c1&&c2)) {
            free(theTerms->thePowers);
            theTerms->thePowers = nil;
        }
    }
}

//__________________________________________________________________________________
void _Polynomial::Convert2OperationForm (void)
{
    // see if anything needs to be done at all
    if (compList1.countitems()&&!theTerms->thePowers) {
        long  n = variableIndex.countitems()+1, m = compList1.countitems();
        long  i1,i2,i,lp=n-2,*scratch=nil,index=0;
        if (n>1) {
            theTerms->thePowers = (long*)MemAllocate (theTerms->allocTerms*sizeof(long)*(n-1));
            scratch = new long [n-1];
            checkPointer(scratch);
            memset (scratch,0,sizeof(long)*(n-1));
            memset (theTerms->thePowers,0,theTerms->allocTerms*sizeof(long)*(n-1));
        }

        for (i=0; i<m; i++) { // loop over all commands
            i1=compList1(i);
            i2=compList2(i);
            if (i1==lp) { // operations with the last var
                if (i2>0) { // do several iterations
                    for (long k=0; k<i2; k++,index++) {
                        scratch[lp]++;
                        theTerms->WriteTerm(scratch,index);
                    }
                } else {
                    if (!i2) {
                        theTerms->WriteTerm(scratch,index);
                    } else {
                        scratch[lp]-=i2;
                        theTerms->WriteTerm(scratch,index);
                    }
                    index++;
                }
            } else {
                if (i1<0) {
                    bool compute = i2<0;
                    i1=-i1-1;
                    if (i2<0) {
                        i2=-i2;
                    }
                    if (i2==1) {
                        scratch[i1]++;
                    } else {
                        scratch[i1]+=i2;
                    }
                    if (compute) {
                        continue;
                    }
                    theTerms->WriteTerm(scratch,index);
                    index++;
                } else {
                    bool compute = i2<0;
                    for (long k=i1+1; k<=lp; k++) {
                        scratch[k]=0;
                    }
                    if (i2<0) {
                        i2=-i2;
                    }
                    if (i2==1) {
                        scratch[i1]++;
                    } else {
                        scratch[i1]+=i2;
                    }

                    if (compute) {
                        continue;
                    }
                    theTerms->WriteTerm(scratch,index);
                    index++;
                }


            }
        }
        if (scratch) {
            delete scratch;
        }
        compList1.Clear();
        compList2.Clear();
    }

}
//__________________________________________________________________________________
_MathObject* _Polynomial::Compute (void)
{
    return new _Constant (ComputePolynomial());
}

//__________________________________________________________________________________
_Parameter _Polynomial::ComputePolynomial (void)
//assumed that the poly is already in the comp form
{
    Convert2ComputationForm();
    long n = variableIndex.countitems()+1;
    _Parameter * varValues = new _Parameter[n];
    checkPointer(varValues);
    for (long i=0; i<n-1; i++) {
        varValues[i]=LocateVar(variableIndex(i))->Compute()->Value();
    }
    _Parameter result = ComputeP (varValues, theTerms->GetCoeff(),n,compList1.countitems(),
                                  compList1.quickArrayAccess(), compList2.quickArrayAccess());
    delete varValues;
    return result;
}
//__________________________________________________________________________________

_Parameter      _Polynomial::ComputeP (_Parameter* varValues, _Parameter* compCoeff, long n, long m, long* c1, long* c2)
{
    _Parameter * holder = new _Parameter[n],term=1,result = 0,lv;
    checkPointer(holder);
    long  i1,i2,i,lp=n-2;
    for (i=0; i<n-1; i++) {
        holder[i]=1;
    }
    lv = n>1?varValues[n-2]:1;
    for (i=0; i<m; i++,c1++,c2++) { // loop over all commands
        i1=*c1;
        i2=*c2;
        if (i1==lp) { // operations with the last var
            if (i2>0) { // do several iterations
                for (long k=0; k<i2; k++,compCoeff++) {
                    term*=lv;
                    result+=term**compCoeff;
                }
            } else {
                if (!i2) {
                    result+=*compCoeff;
                } else {
                    term*=_PolynomialData::BinaryRaise(lv,-i2);
                    result+=term**compCoeff;
                }
                compCoeff++;
            }
        } else {
            if (i1<0) {
                bool compute = i2<0;
                i1=-i1-1;
                if (compute) {
                    i2=-i2;
                }
                if (i2==1) {
                    holder[i1]*=varValues[i1];
                    term*=varValues[i1];
                } else {
                    _Parameter p2 = _PolynomialData::BinaryRaise(varValues[i1],i2);
                    holder[i1]*=p2;
                    term*=p2;
                }
                if (compute) {
                    continue;
                }
                result+=term**compCoeff;
                compCoeff++;
            } else {
                bool compute = i2<0;
                long k;
                for (k=i1+1; k<=lp; k++) {
                    holder[k]=1;
                }
                if (i2<0) {
                    i2=-i2;
                }
                if (i2==1) {
                    holder[i1]*=varValues[i1];
                } else {
                    holder[i1]*=_PolynomialData::BinaryRaise(varValues[i1],i2);
                }
                term = 1;
                for (k=0; k<=i1; k++) {
                    term*=holder[k];
                }

                if (compute) {
                    continue;
                }
                result+=term**compCoeff;
                compCoeff++;
            }


        }
    }
    delete holder;
    return result;
}
//__________________________________________________________________________________

void    _Polynomial::RankTerms(_SimpleList* receptacle)
{
    receptacle->Clear();
    _Parameter  logTop = log(topPolyCap), *theCoeff = theTerms->GetCoeff();

    for (long i = 0; i<theTerms->actTerms; i++,theCoeff++) {
        (*receptacle)<<(long)(log(fabs(*theCoeff))+logTop*theTerms->SumOfPowers(i));
    }
}

//__________________________________________________________________________________

BaseObj* _Polynomial::toStr (void)
{
    _String result (10, true);
    if (theTerms->NumberOfTerms()) {
        long i;
        _List _varNames;
        for (i=0; i<variableIndex.countitems(); i++) {
            _varNames<<LocateVar(variableIndex(i))->GetName();
        }

        for (i=0; i<theTerms->NumberOfTerms(); i++) {
            char        number [100];
            snprintf (number, sizeof(number),PRINTF_FORMAT_STRING,theTerms->GetCoeff(i));
            if (i>0 && number[0]!='-') {
                result<<'+';
            }

            result<<number;
            bool       firstN = theTerms->IsFirstANumber();
            if (i>0 || !firstN) {
                result<<'*';
                long *cT = theTerms->GetTerm(i);
                bool      printedFirst = false;
                for (long k=0; k<variableIndex.countitems(); k++,cT++) {
                    if (*cT>0) {
                        if (printedFirst) {
                            result<<'*';
                        } else {
                            printedFirst = true;
                        }

                        result<<(_String*)_varNames(k);
                        if (*cT>1) {
                            result<<'^';
                            _String st (*cT);
                            result<< *cT;
                        }
                    }
                }
            }
        }
    } else {
        _String*s = (_String*)compList1.toStr();
        result<<s;
        result<<'\n';
        DeleteObject(s);
        s = (_String*)compList2.toStr();
        result<<s;
        result<<'\n';
        DeleteObject(s);
    }
    result.Finalize();
    return result.makeDynamic();
}

//__________________________________________________________________________________

void _Polynomial::toFileStr (FILE*f)
{
    if (theTerms->NumberOfTerms()&&theTerms->thePowers) {
        fprintf(f,"p(");
        _List _varNames;
        long i;
        for (i=0; i<variableIndex.countitems(); i++) {
            _varNames<<LocateVar(variableIndex(i))->GetName();
            fprintf(f,"%s",((_String*)_varNames(i))->sData);
            if (i<variableIndex.countitems()-1) {
                fprintf(f,",");
            }
        }
        fprintf(f,")=");
        for (i=0; i<theTerms->NumberOfTerms(); i++) {
            char number [100];
            snprintf (number, sizeof(number),PRINTF_FORMAT_STRING,theTerms->GetCoeff(i));
            if ((i>0)&&(number[0]!='-')) {
                fprintf(f,"+");
            }
            fprintf(f,"%s",number);
            if ((i>0)||!theTerms->IsFirstANumber()) {
                fprintf(f,"*");
                long *cT = theTerms->GetTerm(i);
                for (long k=0; k<variableIndex.countitems(); k++,cT++) {
                    if (*cT>0) {
                        fprintf(f,"%s",((_String*)_varNames(k))->sData);
                        if (*cT>1) {
                            fprintf(f,"^");;
                            fprintf(f,"%ld",*cT);
                        }
                    }
                }
            }
        }
    } else {
        compList1.toFileStr(f);
        compList2.toFileStr(f);
    }
}

//__________________________________________________________________________________
void    _Polynomial::ScanForVariables (_AVLList&l, bool globals, _AVLListX* tagger, long weight)
{
    for (long i = 0; i<variableIndex.lLength; i++) {
        long vi = variableIndex(i);

        _Variable* v = LocateVar (vi);
        if (v->IsGlobal()) {
            if (globals) {
                l.Insert ((BaseRef)vi);
                if (tagger) {
                    tagger->UpdateValue((BaseRef)vi, weight, 0);
                }
            }
        } else {
            l.Insert ((BaseRef)vi);
            if (tagger) {
                tagger->UpdateValue((BaseRef)vi, weight, 0);
            }
       }
    }
}

//__________________________________________________________________________________
bool    _Polynomial::IsObjectEmpty (void)
{
    if (compList1.countitems()) {
        return false;
    }
    if (theTerms->NumberOfTerms()) {
        if (theTerms->NumberOfTerms()==1) {
            if(theTerms->IsFirstANumber()) {
                return theTerms->theCoeff[0]==0.0;
            }
        }
        return false;
    }
    return true;
}
//__________________________________________________________________________________
bool    _Polynomial::HasChanged (void)
{
    for (long k=variableIndex.countitems()-1; k>=0; k--) {
        if (LocateVar(variableIndex(k))->HasChanged()) {
            return true;
        }
    }
    return false;
}

//__________________________________________________________________________________
bool    _Polynomial::IsMaxElement (_Parameter bench)
{
    _Parameter* tc = theTerms->GetCoeff();
    for (long k=0; k<theTerms->actTerms; k++, tc++) {
        if (fabs(*tc)>=bench) {
            return true;
        }
    }
    return false;
}

//__________________________________________________________________________________
void    SetPolyTermCap (long t)
{
    polyTermCap = t;
}
