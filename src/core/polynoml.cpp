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

#include <string.h>
#include <math.h>

#include "polynoml.h"
#include "global_things.h"

using namespace hy_global;


hyFloat dropPrecision        = log(1e-9),
           drop2Precision      = log(1e-18),
           topPolyCap           = 5.0,
           dropTerms          = 1.0,
           enforcePolyCap      = 1.0,
           *varCheckArray       = nil,
            dropThreshold        = 0.0,
            maximumPolyTermsPerVariable  = 500,
            maxPolynomialExpIterates     = 20,
            polynomialExpPrecision         = 1e-10;

long       varCheckAllocated  = 0,
           polyTermCap        = 0x1fffffff;

bool       checkReset         = true;


//__________________________________________________________________________________
_PolynomialData::_PolynomialData (void) : BaseObj(){
    theCoeff = nil;
    thePowers = nil;
    numberVars = 0;
    actTerms = 0;
    allocTerms = 0;
}

//__________________________________________________________________________________
_PolynomialData::_PolynomialData (long vars)
{
    numberVars = vars>=0?vars:0;
    theCoeff = (hyFloat*)MemAllocate (sizeof(hyFloat)*POLY_DATA_INCREMENT);
    if (numberVars) {
        thePowers =  (long*)MemAllocate (sizeof(long)*POLY_DATA_INCREMENT*vars);
    } else {
        thePowers = nil;
    }
    allocTerms = POLY_DATA_INCREMENT;
    actTerms = 0;
}

//__________________________________________________________________________________
_PolynomialData::_PolynomialData (_PolynomialData const& source) {
    Duplicate(&source);
}

//__________________________________________________________________________________
_PolynomialData const & _PolynomialData::operator = (_PolynomialData const& source) {
    if (this != &source) {
        Duplicate(&source);
    }
    return *this;
}

//__________________________________________________________________________________
_PolynomialData::_PolynomialData (long vars, long terms, hyFloat* theCoeffs)
{
    numberVars = vars>=0?vars:0;
    allocTerms = (terms/POLY_DATA_INCREMENT+1)*POLY_DATA_INCREMENT;
    actTerms = terms;
    theCoeff = (hyFloat*)MemAllocate (sizeof(hyFloat)*allocTerms);
    memcpy (theCoeff,theCoeffs,sizeof(hyFloat)*terms);
    thePowers = nil;
}

//__________________________________________________________________________________
_PolynomialData::~_PolynomialData (void)
{
    if (CanFreeMe()) {
        if (theCoeff) {
            free (theCoeff);
        }
        if (thePowers) {
            free (thePowers);
        }
        allocTerms = 0;
        actTerms = 0;
    }
}

//__________________________________________________________________________________
BaseRef _PolynomialData::makeDynamic (void) const {
    _PolynomialData *ret = new _PolynomialData;
    ret->allocTerms = allocTerms;
    ret->actTerms = actTerms;
    ret->numberVars = numberVars;
    ret->theCoeff = theCoeff;
    ret->thePowers = thePowers;
    return ret;
}

//__________________________________________________________________________________
void _PolynomialData::Duplicate (BaseRefConst source)
{
    _PolynomialData const *s = (_PolynomialData const*)source;
    allocTerms = s->allocTerms;
    actTerms = s->actTerms;
    numberVars = s->numberVars;
    if (actTerms) {
        theCoeff = (hyFloat*)MemAllocate (sizeof(hyFloat)*allocTerms);
        memcpy (theCoeff, s->theCoeff, sizeof (hyFloat)*actTerms);
        if (numberVars) {
            thePowers = (long*)MemAllocate (sizeof(long)*allocTerms*numberVars);
            memcpy (thePowers, s->thePowers, sizeof(long)*actTerms*numberVars);
        } else {
            thePowers = nil;
        }

    } else {
        thePowers = nil;
        theCoeff = nil;
    }
}

//__________________________________________________________________________________
bool    _PolynomialData::checkMe (void)
{
    if (actTerms>1) {
        long *t1 = GetTerm(0),
              *t2;
        for (long i=1; i<actTerms; i++) {
            t2 = GetTerm(i);
            if (CompareTerms (t1,t2)>=0) {
                HandleApplicationError ("Internal polynomial error!", false);
                return false;
            }
            t1 = t2;
        }
    }
    return true;
}
//__________________________________________________________________________________
long*   _PolynomialData::GetTerm (long index)
{
    return thePowers+index*numberVars;
}

//__________________________________________________________________________________
void    _PolynomialData::AddTerm (long * theTerm, hyFloat theC)
{
    if (!(allocTerms-actTerms)) { // no space left - reallocate
        allocTerms+=POLY_DATA_INCREMENT;
        if (theCoeff) {
            theCoeff = (hyFloat*)MemReallocate ((char*)theCoeff, allocTerms*sizeof(hyFloat));
        } else {
            theCoeff = (hyFloat*)MemAllocate (allocTerms*sizeof(hyFloat));
        }
        if (numberVars) {
            if (thePowers) {
                thePowers = (long*)MemReallocate ((char*)thePowers, allocTerms*sizeof(long)*numberVars);
            } else {
                thePowers = (long*)MemAllocate (allocTerms*sizeof(long)*numberVars);
            }
        }
    }
    theCoeff[actTerms] = theC;
    if (numberVars>2) {
        long *stTerm = thePowers+actTerms*numberVars;
        for (long i=0; i<numberVars; i++,theTerm++,stTerm++) {
            *stTerm=*theTerm;
        }
//      memcpy (thePowers+actTerms*numberVars, theTerm, numberVars*sizeof(long));
    } else {
        if (numberVars == 2) {
            *(thePowers+actTerms*numberVars) = *theTerm;
            *(thePowers+actTerms*numberVars+1) = *(theTerm+1);
        } else if (numberVars == 1) {
            *(thePowers+actTerms*numberVars) = *theTerm;
        }
    }
    actTerms++;

}

//__________________________________________________________________________________
void    _PolynomialData::WriteTerm (long * theTerm, long index)
{
    if (numberVars>2) {
        long *stTerm = thePowers+index*numberVars;
        for (long i=0; i<numberVars; i++,theTerm++,stTerm++) {
            *stTerm=*theTerm;
        }
//      memcpy (thePowers+index*numberVars, theTerm, numberVars*sizeof(long));
    } else {
        if (numberVars == 2) {
            *(thePowers+index*numberVars) = *theTerm;
            *(thePowers+index*numberVars+1) = *(theTerm+1);
        } else if (numberVars == 1) {
            *(thePowers+index*numberVars) = *theTerm;
        }
    }
    actTerms++;

}

//__________________________________________________________________________________
void    _PolynomialData::AddTerm (long * theTerm, hyFloat theC, long* reindexer, long actLength)
{
    if (!(allocTerms-actTerms)) { // no space left - reallocate
        allocTerms+=POLY_DATA_INCREMENT;
        if (theCoeff) {
            theCoeff = (hyFloat*)MemReallocate ((char*)theCoeff, allocTerms*sizeof(hyFloat));
        } else {
            theCoeff = (hyFloat*)MemAllocate (allocTerms*sizeof(hyFloat));
        }
        if (numberVars) {
            if (thePowers) {
                thePowers = (long*)MemReallocate ((char*)thePowers, allocTerms*sizeof(long)*numberVars);
            } else {
                thePowers = (long*)MemAllocate (allocTerms*sizeof(long)*numberVars);
            }
        }
    }
    theCoeff[actTerms] = theC;
    if (numberVars>2) {
        long *newTerm = thePowers+actTerms*numberVars;
        for (long i=0; i<numberVars; i++) {
            newTerm[i] = 0L;
        }
 
        for (long i=0; i<actLength; i++) {
            newTerm[reindexer[i]] = theTerm[i];
        }

    } else {
        *(thePowers+actTerms*numberVars) = 0;
        *(thePowers+actTerms*numberVars+1) = 0;
        *(thePowers+actTerms*numberVars+*reindexer) = *theTerm;
    }
    actTerms++;

}

//__________________________________________________________________________________
void    _PolynomialData::AddTerm (hyFloat theC)
{
    if (numberVars==0) {
        AddTerm (nil,0);
    } else {
        if (!(allocTerms-actTerms)) { // no space left - reallocate
            allocTerms+=POLY_DATA_INCREMENT;
            if (theCoeff) {
                theCoeff = (hyFloat*)MemReallocate ((char*)theCoeff, allocTerms*sizeof(hyFloat));
            } else {
                theCoeff = (hyFloat*)MemAllocate (allocTerms*sizeof(hyFloat));
            }
            if (numberVars) {
                if (thePowers) {
                    thePowers = (long*)MemReallocate ((char*)thePowers, allocTerms*sizeof(long)*numberVars);
                } else {
                    thePowers = (long*)MemAllocate (allocTerms*sizeof(long)*numberVars);
                }
            }
        }
        memmove (thePowers+numberVars,thePowers, numberVars*actTerms*sizeof(long));
        memmove (theCoeff+1,theCoeff,sizeof(hyFloat)*actTerms);
        *theCoeff = theC;
        for (long i=0; i<numberVars; *thePowers=0,i++,thePowers++) {}
        thePowers-=numberVars;
//      memset (thePowers,0,numberVars*sizeof(long));
        actTerms++;
    }

}


//__________________________________________________________________________________
void    _PolynomialData::DeleteTerm (long index)
{
    actTerms--;
    if (index!=actTerms) { // shift stuff
        memmove ((void*)(theCoeff+index), (void*)(theCoeff+index+1), (actTerms-index)*sizeof(hyFloat));
        if (numberVars)
            memmove ((void*)(thePowers+index*numberVars), (void*)(thePowers+(1+index)*numberVars),
                     (actTerms-index)*numberVars*sizeof(long));
    }
    if (allocTerms-actTerms>POLY_DATA_INCREMENT) {
        allocTerms-=POLY_DATA_INCREMENT;
        theCoeff = (hyFloat*)MemReallocate ((char*)theCoeff, allocTerms*sizeof(hyFloat));
        if (numberVars) {
            thePowers = (long*)MemReallocate ((char*)thePowers, allocTerms*sizeof(long)*numberVars);
        }
    }
}

//__________________________________________________________________________________
long    _PolynomialData::SumOfPowers (long index)
{
    if (numberVars) {
        long* theTerm = GetTerm (index), res = 0;
        for (long i=0; i<numberVars; i++, theTerm++) {
            res+=*theTerm;
        }
        return res;
    }
    return 0;
}

//__________________________________________________________________________________
long    _PolynomialData::WeightedSumOfPowers (long index, hyFloat* w)
{
    if (numberVars) {
        long* theTerm = GetTerm (index), res = 0;
        for (long i=0; i<numberVars; i++, theTerm++, w++) {
            res+=((hyFloat)*theTerm)**w;
        }
        return res;
    }
    return 0;
}

//__________________________________________________________________________________
bool    _PolynomialData::IsFirstANumber (void)
{
    long* fst = GetTerm(0);
    for (long i=0; i<numberVars; i++)
        if (fst[i]) {
            return false;
        }

    return true;
}


//__________________________________________________________________________________
void    _PolynomialData::MultiplyTerms (long* target, long* s1, long* s2)
{
    for (long i=0; i<numberVars; i++,target++,s1++,s2++) {
        *target = *s1+*s2;
    }
}

//__________________________________________________________________________________
void    _PolynomialData::RaiseTerm (long* target, long power) {
    for (long i=0L; i<numberVars; i++) {
        target[i] *= power;
    }
}

//__________________________________________________________________________________
hyFloat  _PolynomialData::BinaryRaise (hyFloat base, long pwr) {
    hyFloat result = 1;
    char bits[sizeof(long)*8];
    
    bool invert = false;
    if (pwr < 0) {
        invert = true;
        pwr = -pwr;
    }
    
    unsigned char nLength = 0;
    while (pwr) {
        bits[nLength]=pwr%2;
        pwr = pwr >> 1;
        nLength++;
    }
    while (nLength) {
        nLength--;
        if (bits[nLength]) {
            result*=base;
        }
        if (nLength) {
            result*=result;
        }
    }
    return invert ? 1./result : result;
}

//__________________________________________________________________________________
long    _PolynomialData::FindTerm (long* theTerm,long* reIndex, long start)
{
    long top=actTerms-1, bottom=start, middle;
    char comp;

    if (top==-1) {
        return -2;
    }

    while (top>bottom) {
        middle = (top+bottom)/2;
        // compare the two
        comp = CompareTerms (GetTerm (reIndex[middle]), theTerm);
        if (comp==1) {
            top = middle==top?top-1:middle;
        } else if (comp==-1) {
            bottom = middle==bottom?bottom+1:middle;
        } else {
            return middle;
        }


    }
    middle = top;
    comp = CompareTerms (GetTerm (reIndex[middle]), theTerm);
    if (!comp) {
        return middle;
    }
    return comp<0?-middle-3:-middle-2;
}

//__________________________________________________________________________________
void    _PolynomialData::ResortTerms (long* reIndex)
{
    long i,j,*source,*target, deleted = 0;

    hyFloat* newCoeff = (hyFloat*)MemAllocate (allocTerms*sizeof(hyFloat));
    long*       newPowers = (long*)MemAllocate(allocTerms*numberVars*sizeof(long));

    // pass 1 to check for deletions
    for (i=0; i<actTerms; i++, reIndex++) {
        if (checkTerm (theCoeff[*reIndex], *reIndex)) {
            newCoeff[i] = theCoeff[*reIndex];
        } else {
            newCoeff[i]=0.0;
        }
    }

    reIndex-=actTerms;
    // pass 2 copy
    for (i=0; i<actTerms; i++, newCoeff++, reIndex++) {
        if (*newCoeff!=0.0) {
            if (deleted) {
                *(newCoeff-deleted)=*newCoeff;
            }
            target = newPowers+numberVars*(i-deleted);
            source = thePowers+numberVars*(*reIndex);
            for (j=0; j<numberVars; j++,source++,target++) {
                *target = *source;
            }
        } else {
            deleted++;
        }
    }
    free (theCoeff);
    free (thePowers);
    theCoeff = newCoeff-actTerms;
    thePowers = newPowers;
    actTerms -= deleted;
    if (allocTerms-actTerms>POLY_DATA_INCREMENT) {
        long theCut = ((allocTerms-actTerms)/POLY_DATA_INCREMENT)*POLY_DATA_INCREMENT;
        allocTerms-=theCut;
        theCoeff = (hyFloat*)MemReallocate ((char*)theCoeff, allocTerms*sizeof(hyFloat));
        if (numberVars) {
            thePowers = (long*)MemReallocate ((char*)thePowers, allocTerms*sizeof(long)*numberVars);
        }
    }

}
//__________________________________________________________________________________
void    _PolynomialData::ChopTerms (void)
{
    long maxAllowedTerms = maximumPolyTermsPerVariable*numberVars;
    if (actTerms<=maxAllowedTerms) {
        return;    // nothing to do
    } else {
        _SimpleList terms, index;
        long i,*qa,k,j;
        hyFloat  logTop = log(topPolyCap),*qc;
        for (i = 0; i<actTerms; i++,theCoeff++) {
            index<<i;
            terms<<(long)(log(fabs(*theCoeff))+logTop*SumOfPowers(i));
        }
        SortLists (&terms,&index);
        terms.Clear();
        theCoeff -= actTerms;
        qa = index.quickArrayAccess()+maxAllowedTerms;
        for (i=maxAllowedTerms; i<actTerms; i++,qa++) {
            theCoeff[*qa] = 0.0;
        }

        allocTerms = (maxAllowedTerms/POLY_DATA_INCREMENT+1)*POLY_DATA_INCREMENT;
        hyFloat* newCoeff = (hyFloat*)MemAllocate (allocTerms*sizeof(hyFloat)),*nP = newCoeff;
        long*       newPowers = (long*)MemAllocate(allocTerms*numberVars*sizeof(long)), *target, *source;

        target = newPowers;
        source = thePowers;
        k = 0;
        qc = GetCoeff();
        for (i=0; i<actTerms; i++, qc++, source+=numberVars) {
            if (*qc!=0.0) {
                *nP = *qc;
                nP++;
                for (j=0; j<numberVars; j++,source++,target++) {
                    *target = *source;
                }
            } else {
                k++;
            }
        }
        free (theCoeff);
        free (thePowers);
        theCoeff = newCoeff;
        thePowers = newPowers;
        actTerms -= k;
    }

}

//__________________________________________________________________________________
void  _PolynomialData::RearrangeTerm (long* target, long* source, long* markup, long items)
{
    for (long i=0; i<items; i++, source++, markup++) {
        target[*markup] = *source;
    }
}
//__________________________________________________________________________________
char    _PolynomialData::CompareTerms (long* s1, long* s2)
{
    for (long i=0; i<numberVars; i++) {
        long comp = s1[i]-s2[i];
        if (comp>0L) {
            return 1;
        }
        if (comp<0L) {
            return -1;
        }
    }
    return 0;
}

//__________________________________________________________________________________
char    _PolynomialData::CompareTerms (long* s1, long* s2, long* secondReindex, long actLength)
{
    long second_index = 0;
    for (long k = 0L; k < numberVars; k++) {
      long v1 = s1[k],
            v2 = 0;
      if (second_index < actLength && k == secondReindex[second_index]) {
        v2 = s2[second_index];
        second_index++;
      }
      
      if (v1 != v2) {
        return v1 > v2 ? 1 : -1;
      }
    }
  
    return 0;
  
    /*long comp,i;
  
    for (i=0; i<actLength; i++,s1++,s2++,secondReindex++) {
        comp = *secondReindex-i;
        if (comp>0) {
            if (*s1) {
                return 1;
            } else {
                s2--;
                secondReindex--;
                actLength++;
                continue;
            }
        }
        comp = *s1-*s2;
        if (comp>0) {
            return 1;
        }
        if (comp<0) {
            return -1;
        }
    }
    for (; i<numberVars; i++,s1++) {
        if (*s1) {
            return 1;
        }
    }
    return 0;*/
}

//__________________________________________________________________________________
char    _PolynomialData::CompareTerms (long* s1, long* s2, long* firstReindex, long* secondReindex,long actLength1,  long actLength2)
{
    long first_index = 0L,
         second_index = 0L,
         k;
  
    for (k = 0; k < numberVars; k ++) {
      
      long v1 = 0,
           v2 = 0;
      
      if (first_index < actLength1) {
        if (k == firstReindex[first_index]) {
          v1 = s1 [first_index];
          first_index ++;
        }
      }

      if (second_index < actLength2) {
        if (k == secondReindex[second_index]) {
          v2 = s2 [second_index];
          second_index ++;
        }
      }
      
      if (v1 != v2) {
        return v1 < v2 ? -1 : 1;
      }
      
    }
      
    return 0;
  
    /*bool secondLonger = actLength1<actLength2,
         firstLonger  = actLength1>actLength2;

    long minLength = actLength1<actLength2?actLength1:actLength2,
         i;

    for (i=0; i<minLength; i++,s1++,s2++,firstReindex++,secondReindex++) {
        long comp = *firstReindex-*secondReindex;
        if (comp<0) {
            if (*s1) {
                return 1;
            }
            secondReindex--;
            s2--;

            if (firstLonger && minLength!=actLength1) {
                minLength++;
            }

            continue;
        }
        if (comp>0) {
            if (*s2) {
                return -1;
            }
            firstReindex--;
            s1--;

            if (minLength!=actLength2 && secondLonger) {
                minLength++;
            }
            continue;
        }

        comp = *s1-*s2;
        if (comp>0) {
            return 1;
        }
        if (comp<0) {
            return -1;
        }
    }
    if (actLength1>minLength)
        for (; i<actLength1; i++,s1++) {
            if (*s1) {
                return 1;
            }
        }
    else
        for (; i<actLength2; i++,s2++) {
            if (*s2) {
                return -1;
            }
        }
    return 0;*/
}

//__________________________________________________________________________________
bool        _PolynomialData::checkTerm (hyFloat myCoeff, long myIndex)
{
    if (myCoeff==0.0) { // delete in any case
        return false;
    }
    if (checkReset) {
        checkReset = false;
        dropThreshold = dropPrecision+log(fabs(myCoeff));
        if (dropThreshold<drop2Precision) {
            dropThreshold = drop2Precision;
            if (enforcePolyCap) {
                dropThreshold += SumOfPowers(myIndex)*log(topPolyCap);
            } else {
                dropThreshold += WeightedSumOfPowers(myIndex,varCheckArray);
            }
            return false;
        }
        if (enforcePolyCap) {
            dropThreshold += SumOfPowers(myIndex)*log(topPolyCap);
        } else {
            dropThreshold += WeightedSumOfPowers(myIndex,varCheckArray);
        }
        return true;

    }
    if (dropTerms) {
        if (enforcePolyCap) {
            if (log(fabs(myCoeff))+log(topPolyCap)*(hyFloat)(SumOfPowers(myIndex))<dropThreshold) {
                return false;
            }
        } else if (enforcePolyCap) {
            if (log(fabs(myCoeff))+WeightedSumOfPowers(myIndex,varCheckArray)<dropThreshold) {
                return false;
            }
        }
    }
    return true;
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (void) {
    theTerms = new _PolynomialData;
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_SimpleList& vars) {
    variableIndex.Duplicate (&vars);
    theTerms = new _PolynomialData (vars.countitems());
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (hyFloat value) {
// a constant polynomial
    theTerms = new _PolynomialData;
    theTerms->AddTerm (nil, value);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_Variable& v)
// a monic monomial
{
    variableIndex<<v.get_index();
    theTerms = new _PolynomialData (1);
    long  vIndex = 1;
    theTerms->AddTerm (&vIndex,1.0);
}

//__________________________________________________________________________________

_Polynomial::_Polynomial (_Polynomial const& p) {
    *this = p;
}

//__________________________________________________________________________________

_Polynomial const & _Polynomial::operator = (_Polynomial const& p) {
    if (this != &p) {
        variableIndex.Duplicate (&p.variableIndex);
        theTerms = new _PolynomialData;
        if (p.theTerms) {
            theTerms->Duplicate (p.theTerms);
        } else {
            theTerms->numberVars = variableIndex.countitems();
        }
        compList1.Duplicate (&p.compList1);
        compList2.Duplicate (&p.compList2);
    }
    return *this;
}
//__________________________________________________________________________________

_Polynomial::~_Polynomial ()
{
    if (theTerms) {
        DeleteObject(theTerms);
    }
}

//__________________________________________________________________________________

BaseObj*    _Polynomial::makeDynamic(void) const{
    _Polynomial* res = new _Polynomial;
 
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

void    _Polynomial::Duplicate  (BaseRefConst tp) {
    _Polynomial const* p = (_Polynomial const*)tp;
    variableIndex.Clear();
    variableIndex.Duplicate (&p->variableIndex);
    compList1.Duplicate (&p->compList1);
    compList2.Duplicate (&p->compList2);
    DeleteObject(theTerms);
    if (p->theTerms) {
        theTerms =  new _PolynomialData (*(p->theTerms));
    }
}
//__________________________________________________________________________________


HBLObjectRef _Polynomial::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context,HBLObjectRef cache)
{
  switch (opCode) { // first check operations without arguments
    case HY_OP_CODE_TYPE: // Type
      return Type(cache);
  }
  
  _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
  switch (opCode) { // next check operations without arguments or with one argument
    case HY_OP_CODE_SUB: // -
      if (arg0) {
        return Sub(arg0,cache);
      } else {
        return Minus(cache);
      }
    case HY_OP_CODE_ADD: // +
      if (arg0) {
        return Add(arg0,cache);
      } else {
        return Sum (cache);
      }
  }
  
  if (arg0) {
    switch (opCode) {
      case HY_OP_CODE_MUL: //*
        return Mult(arg0,cache);
      case HY_OP_CODE_POWER: // ^
        return Raise(arg0,cache);
    }
    
  }
  
  switch (opCode) {
    case HY_OP_CODE_MUL:
    case HY_OP_CODE_POWER:
      WarnWrongNumberOfArguments (this, opCode,context, arguments);
      break;
    default:
      WarnNotDefined (this, opCode,context);
  }
  
  return new _MathObject;
  
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
_MathObject* _Polynomial::Add (_MathObject* m, HBLObjectRef cache) {
    return Plus(m,false,cache);
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
                varCheckArray = (hyFloat*)MemAllocate (varCheckAllocated*sizeof(hyFloat));
                hyFloat lb, ub;
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
    hyFloat myCoeff = theTerms->GetCoeff(myIndex);
    if (!theTerms->checkTerm (myCoeff,myIndex)) {
        theTerms->DeleteTerm(myIndex);
    }

}

//__________________________________________________________________________________

bool         _Polynomial::Equal(_MathObject* m)
{
    bool result = false;
    if (m->ObjectClass() == POLYNOMIAL || m->ObjectClass() == NUMBER) {
        _Polynomial * diff = (_Polynomial *)Sub(m, nil);
        if (diff) {
            _Constant * v = (_Constant*)diff->IsANumber(true);
            if (v) {
                result = CheckEqual (v->Value(),0.0);
                DeleteObject (v);
            }
            DeleteObject (diff);
        }

    }
    return result;
}


//__________________________________________________________________________________
_MathObject* _Polynomial::Plus (_MathObject* m, bool subtract, HBLObjectRef cache) {
    long objectT = m->ObjectClass();

    if (objectT==1) { // a number
        Convert2OperationForm();
        _Polynomial* result = new _Polynomial (*this);
        hyFloat mV = subtract?-m->Value():m->Value();

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
                    return p2->Plus (&coef1,false,nil);
                } else {
                    _Polynomial *ppp = (_Polynomial*)p2->Plus (&coef1,true,nil);
                    hyFloat  *invC = ppp->theTerms->theCoeff;
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
                    hyFloat *invC = ppp->theTerms->theCoeff;
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
                return    Plus (&coef2,subtract,nil);
            } else {
                if (!subtract) {
                    return new _Polynomial(*this);
                } else {
                    _Polynomial* ppp = new _Polynomial (*this);
                    hyFloat *invC = ppp->theTerms->theCoeff;
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

        hyFloat * coeff1 = theTerms->GetCoeff(),
                     *coeff2 = p2->theTerms->GetCoeff();

        _Polynomial* res;

        if (variableIndex.Equal(p2->variableIndex))
            // same variable arrays - proceed to the operation directly
        {
            res = new _Polynomial (variableIndex); // create a blank new result holder
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
            if ( merge1.countitems()==mergedVariables.countitems() || merge2.countitems()==mergedVariables.countitems()) {
                bool  firstBigger = merge1.countitems()==mergedVariables.countitems(); // is the first poly bigger than the second
              
                long  *reindexList = firstBigger?merge2.quickArrayAccess():merge1.quickArrayAccess(),
                       reindexLength = firstBigger?merge2.countitems():merge1.countitems();
              
              
                res = new _Polynomial (mergedVariables); // create a blank new result holder

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

 
      
      /*
      NLToConsole();
      NLToConsole();
      BufferToConsole (_String((_String*)this->toStr()).getStr());
      BufferToConsole("\n+\n");
      BufferToConsole (_String((_String*)p2->toStr()).getStr());
      BufferToConsole("\n=\n");
      BufferToConsole (_String((_String*)res->toStr()).getStr());
      NLToConsole();
       */
      

      if (!res->theTerms->checkMe()) {
            return nil;
      }
//      res->theTerms->ChopTerms();
        if (res->theTerms->GetNoTerms()==0) {
            DeleteObject (res);
            return new _Polynomial (0.0);
        }
        return res;
    }

    HandleApplicationError ("An incompatible operand was supplied to polynomial addition/subtraction");
    return nil;
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Sub (_MathObject* m,HBLObjectRef cache) {
    return Plus (m,true,cache);
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Minus (HBLObjectRef cache) {
    _Constant min (-1.0);
    return Mult (&min, cache);
}

//__________________________________________________________________________________
_MathObject* _Polynomial::Mult (_MathObject* m, HBLObjectRef cache) {
    long objectT = m->ObjectClass();

    if (objectT==NUMBER) { // a number or a monomial
        Convert2OperationForm();
        hyFloat nb = ((_Constant*)m)->Value();
        if (nb==0.0) {
            return new _Polynomial (0.);
        }

        _Polynomial* result = new _Polynomial (*this);
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

                    hyFloat *tempCoeff = (hyFloat*)MemAllocate (theTerms->NumberOfTerms()*sizeof(hyFloat));
                    memcpy (tempCoeff, theTerms->GetCoeff(), theTerms->NumberOfTerms()*sizeof(hyFloat));


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
                                    hyFloat coe = theTerms->GetCoeff (csp);
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
                                    hyFloat coe = theTerms->GetCoeff (csp);
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
                return p2->Mult (&coef1, nil);
            } else {
                return new _Polynomial();
            }

        }
        if (!p2->variableIndex.countitems()) {
            if (p2->theTerms->NumberOfTerms()) {
                _Constant coef2 (p2->theTerms->GetCoeff(0));
                return Mult (&coef2,nil);
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
            hyFloat  logTop = log(topPolyCap), *theCoeff = theTerms->GetCoeff();
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
                ResetPolynomialCheck(res);
                long* newTerm = new long[res->variableIndex.countitems()],*term1,f;
                for (long i=0; i<nt1; i++) {
                    ref = 0;
                    hyFloat c1 = theTerms->GetCoeff(i);
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

                delete []newTerm;
            } else {
                _SimpleList merge1, merge2, joint;
                joint.Merge (variableIndex, p2->variableIndex,&merge1,&merge2);
                res = new _Polynomial (joint);
                ResetPolynomialCheck(res);
                long  f,nv = res->variableIndex.countitems(),* scratch1 = new long[nv],
                * scratch3,k,* scratch2;
                bool do1 = merge1.countitems()!=0, do2 = merge2.countitems()!=0;
                if (do1) {
                    scratch2 = new long[nv];
                }
                if (do2) {
                    scratch3 = new long[nv];
                }
                for (long i=0; i<nt1; i++) {
                    hyFloat c1 = theTerms->GetCoeff(i);
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
                ResetPolynomialCheck(res);
                long* newTerm = new long[res->variableIndex.countitems()],*term1,f;
                for (long i=0; i<nt1; i++) {
                    if (terms1.BinaryFind(i)<0) {
                        continue;
                    }
                    ref = 0;
                    hyFloat c1 = theTerms->GetCoeff(i);
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
                ResetPolynomialCheck(res);
                long  f,nv = res->variableIndex.countitems(),* scratch1 = new long[nv],
                * scratch3,k,* scratch2;
                bool do1 = merge1.countitems()!=0, do2 = merge2.countitems()!=0;
                if (do1) {
                    scratch2 = new long[nv];
                }
                if (do2) {
                    scratch3 = new long[nv];
                }
                for (long i=0; i<nt1; i++) {
                    hyFloat c1 = theTerms->GetCoeff(i);
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

    HandleApplicationError ("An incompatible operand was supplied to polynomial multiplication");
    return nil;
}


//__________________________________________________________________________________
_MathObject* _Polynomial::Raise (_MathObject* m, HBLObjectRef cache) {
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
    if (objectT==NUMBER) { // a number
        Convert2OperationForm();
        _Polynomial* result;
        if (theTerms->NumberOfTerms()==1) { // just a monomial
            long power = (long)m->Value();
            result = new _Polynomial (*this);
            result->theTerms->RaiseTerm (result->theTerms->GetTerm(0),power);
            result->theTerms->GetCoeff(0)=_PolynomialData::BinaryRaise(result->theTerms->GetCoeff(0),power);

        } else { // binary raise
            result = new _Polynomial(1.0);
          
            _Polynomial* oldR;
            unsigned char bits[sizeof(long)*8];
            long pwr = m->Value(), nLength = 0L;
            while (pwr) {
                bits[nLength]=pwr%2;
                pwr/=2;
                nLength++;
            }
            while (nLength) {
                oldR = result;
                nLength--;
                if (bits[nLength]) {
                    result=(_Polynomial*)result->Mult(this,nil);
                    DeleteObject( oldR);
                }
                oldR = result;
                if (nLength) {
                    result=(_Polynomial*)result->Mult(result,nil);
                    DeleteObject( oldR);
                }
            }
        }
        if (del) {
            DeleteObject(m);
        }
        return result;
    }
    HandleApplicationError ("An incompatible operand was supplied to polynomial raise to power");
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

        delete [] powerDiff;
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
            delete [] scratch;
        }
        compList1.Clear();
        compList2.Clear();
    }

}
//__________________________________________________________________________________
_MathObject* _Polynomial::Compute (void) {
    return new _Constant (ComputePolynomial());
}

//__________________________________________________________________________________
hyFloat _Polynomial::ComputePolynomial (void)
//assumed that the poly is already in the comp form
{
    Convert2ComputationForm();
    long n = variableIndex.countitems()+1;
    hyFloat * varValues = new hyFloat[n];
    for (long i=0; i<n-1; i++) {
        varValues[i]=LocateVar(variableIndex(i))->Compute()->Value();
    }
    hyFloat result = ComputeP (varValues, theTerms->GetCoeff(),n,compList1.countitems(),
                                  compList1.quickArrayAccess(), compList2.quickArrayAccess());
    delete [] varValues;
    return result;
}
//__________________________________________________________________________________

hyFloat      _Polynomial::ComputeP (hyFloat* varValues, hyFloat* compCoeff, long n, long m, long* c1, long* c2)
{
    hyFloat * holder = new hyFloat[n],term=1,result = 0,lv;
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
                    hyFloat p2 = _PolynomialData::BinaryRaise(varValues[i1],i2);
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
    delete [] holder;
    return result;
}
//__________________________________________________________________________________

void    _Polynomial::RankTerms(_SimpleList* receptacle)
{
    receptacle->Clear();
    hyFloat  logTop = log(topPolyCap), *theCoeff = theTerms->GetCoeff();

    for (long i = 0; i<theTerms->actTerms; i++,theCoeff++) {
        (*receptacle)<<(long)(log(fabs(*theCoeff))+logTop*theTerms->SumOfPowers(i));
    }
}



//__________________________________________________________________________________

BaseObj* _Polynomial::toStr (unsigned long padding) {
  _StringBuffer * result = new _StringBuffer (32UL);
  if (theTerms->NumberOfTerms()) {
    long i;
    _List _varNames;
    for (i=0; i<variableIndex.countitems(); i++) {
      _varNames<<LocateVar(variableIndex(i))->GetName();
    }
    
    bool       firstN = theTerms->IsFirstANumber();
    for (i=0; i<theTerms->NumberOfTerms(); i++) {
      char        number [100];
      
      hyFloat  coeff = theTerms->GetCoeff(i);
      
      bool print_coeff = false;
      
      if ((i == 0 && firstN) || !CheckEqual(coeff, 1)) {
        print_coeff = true;
      }
      
      snprintf (number, sizeof(number),PRINTF_FORMAT_STRING,coeff);
      if (i>0 && coeff >= 0.) {
        *result<<'+';
      }
      
      if (print_coeff) {
        *result<<number;
      }
      
      if (i>0 || !firstN) {
        if (print_coeff) {
          *result<<'*';
        }
        long *cT = theTerms->GetTerm(i);
        bool      printedFirst = false;
        for (long k=0; k<variableIndex.countitems(); k++,cT++) {
          if (*cT != 0) {
            if (printedFirst) {
              *result<<'*';
            } else {
              printedFirst = true;
            }
            
            *result<<(_String*)_varNames(k);
            if (*cT>1) {
              *result<<'^'<< _String(*cT);
            } else {
                if (*cT < 0) {
                    *result<<'^'<< _String(*cT).Enquote('(',')');
                }
            }
          }
        }
      }
    }
  } else {
    *result <<_String((_String*)compList1.toStr(padding))
            <<'\n'
            <<_String((_String*)compList2.toStr(padding))
            <<'\n';
  }
  result->TrimSpace();
  return result;
}

//__________________________________________________________________________________

void _Polynomial::toFileStr (FILE*f, unsigned long padding)
{
    if (theTerms->NumberOfTerms()&&theTerms->thePowers) {
        fprintf(f,"p(");
        _List _varNames;
        long i;
        for (i=0; i<variableIndex.countitems(); i++) {
            _varNames<<LocateVar(variableIndex(i))->GetName();
            fprintf(f,"%s",((_String*)_varNames(i))->get_str());
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
                        fprintf(f,"%s",((_String*)_varNames(k))->get_str());
                        if (*cT>1) {
                            fprintf(f,"^");;
                            fprintf(f,"%ld",*cT);
                        }
                    }
                }
            }
        }
    } else {
        compList1.toFileStr(f, padding);
        compList2.toFileStr(f, padding);
    }
}

//__________________________________________________________________________________
void    _Polynomial::ScanForVariables (_AVLList&l, bool globals, _AVLListX* tagger, long weight) const {
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
bool    _Polynomial::HasChanged (bool)
{
    for (long k=variableIndex.countitems()-1; k>=0; k--) {
        if (LocateVar(variableIndex(k))->HasChanged()) {
            return true;
        }
    }
    return false;
}

//__________________________________________________________________________________
bool    _Polynomial::IsMaxElement (hyFloat bench)
{
    hyFloat* tc = theTerms->GetCoeff();
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
