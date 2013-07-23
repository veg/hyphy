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

#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include "string.h"
#include "stdlib.h"
#include "time.h"

#include "constant.h"

#include "legacy_parser.h"

_Formula *chi2 = nil,
         *derchi2 = nil;

long randomCount = 0;

extern _Parameter machineEps;
extern _Parameter tolerance;

long            lastMatrixDeclared = -1,
                dummyVariable1,
                dummyVariable2,
                expressionsParsed = 0;

_Parameter gammaCoeff [7] = {
    2.50662827463100050,
    190.9551718944012,
    -216.8366818451899,
    60.19441758801798,
    -3.087513097785903,
    0.003029460875352382,
    -0.00001345152485367085
};

_Parameter lngammaCoeff [6] = {
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
};

//______________________________________________________________________________
_Constant::_Constant (_Parameter value)
{
    theValue = value;
}

//______________________________________________________________________________
void _Constant::Initialize (void)
{
    BaseObj::Initialize();
    theValue = 0;
}

//______________________________________________________________________________
void _Constant::Duplicate (BaseRef c)
{
    BaseObj::Initialize();
    theValue = dynamic_cast<_Constant*>(c)->theValue;
}

//______________________________________________________________________________
BaseRef _Constant::makeDynamic (void)
{
    _Constant * res = (_Constant*)checkPointer(new _Constant);
    res->Duplicate(this);
    return res;
}

//______________________________________________________________________________
_Constant::_Constant (_String& s)
{
    theValue = atof (s.sData);
}

//______________________________________________________________________________
_Constant::_Constant (void)
{
    theValue = 0;
}

//______________________________________________________________________________
_Parameter    _Constant::Value (void)
{
    return theValue;
}

//______________________________________________________________________________
BaseRef _Constant::toStr(void)
{
    return parameterToString(Value());
}

//______________________________________________________________________________
_PMathObj _Constant::Add (_PMathObj theObj) {
    if (theObj->ObjectClass() == STRING) {
        return new _Constant ((theValue+((_FString*)theObj)->theString->toNum()));
    } else {
        return new _Constant (theValue+theObj->Value());
    }
}

//______________________________________________________________________________
_PMathObj _Constant::Sub (_PMathObj theObj) {
    return new _Constant ((theValue-theObj->Value()));
}

//______________________________________________________________________________
_PMathObj _Constant::Minus (void)
{
    return     new  _Constant (-Value());
}

//______________________________________________________________________________
_PMathObj _Constant::Sum (void)
{
    return     new  _Constant (Value());
}

//______________________________________________________________________________
_PMathObj _Constant::Mult (_PMathObj theObj)
{
    return new _Constant (theValue*theObj->Value());
}
//______________________________________________________________________________
_PMathObj _Constant::Div (_PMathObj theObj)
{
    return new _Constant (theValue/theObj->Value());
}

//______________________________________________________________________________
_PMathObj _Constant::lDiv (_PMathObj theObj) // %
{
    _Parameter       denom = floor (theObj->Value());
    return     denom != 0.0 ? new _Constant  (fmod (floor(theValue), denom)):new _Constant  (floor (theValue));
}
//______________________________________________________________________________
_PMathObj _Constant::longDiv (_PMathObj theObj) // $
{
    _Parameter       denom = floor (theObj->Value());
    return     denom != 0.0 ? new _Constant  (floor (floor (theValue) / denom)):new _Constant  (0.0);
}
//______________________________________________________________________________
_PMathObj _Constant::Raise (_PMathObj theObj)
{
    if (!theObj) {
        return nil;
    }

    _Parameter    base  = Value(),
                  expon = theObj->Value();

    if (base>0.0) {
        return    new  _Constant (exp (log(base)*(expon)));;
    } else {
        if (base<0.0)
            if (CheckEqual (expon, (long)expon)) {
                return new _Constant (((((long)expon)%2)?-1:1)*exp (log(-base)*(expon)));
            } else {
                _String errMsg ("An invalid base/exponent pair passed to ^");
                WarnError (errMsg.sData);
            }

        return     new _Constant (0.0);
    }
}

//______________________________________________________________________________
_PMathObj _Constant::Random (_PMathObj upperB)
{
    if (randomCount == 0) {
        randomCount++;
    }
    _Parameter l = theValue, 
               u=upperB->Value(),
               r = l;
  
    if (u>l) {
        r=genrand_int32();
        r/=RAND_MAX_32;
        r =l+(u-l)*r;
    }
    return new _Constant (r);

}

//______________________________________________________________________________
void     _Constant::Assign (_PMathObj theObj)
{
    this->~_Constant ();
    theValue = theObj->Value();
}

//______________________________________________________________________________
bool     _Constant::Equal (_PMathObj theObj)
{
    return theValue==theObj->Value();
}

//______________________________________________________________________________
_PMathObj _Constant::Abs (void)
{
  return new _Constant(fabs(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Sin (void)
{
  return new _Constant(sin(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Cos (void)
{
    return     new _Constant  (cos(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Tan (void)
{
    return     new _Constant  (tan(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Exp (void)
{
    return     new _Constant  (exp(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::FormatNumberString (_PMathObj p, _PMathObj p2)
{
    long       a1 = p->Value(),
               a2 = p2->Value();

    char       format[32],
               buffer[256];

#ifdef     __USE_LONG_DOUBLE__
    if (a1>=0 && a2>=0) {
        if (a1>0) {
            snprintf    (format,32, "%%%ld.%ldLf",(long)a1,(long)a2);
        } else {
            snprintf    (format,32,"%%.%ldLf",(long)a2);
        }
    } else if (a1>=0) {
        snprintf    (format,32,"%%%ldLf",(long)a1);
    } else if (a2>=0) {
        snprintf    (format,32,"%%.%ldLf",(long)a2);
    } else {
        snprintf    (format,32,"%%Lg");
    }
#else
    if (a1>=0 && a2>=0) {
        if (a1>0) {
            snprintf    (format,32, "%%%ld.%ldf",(long)a1,(long)a2);
        } else {
            snprintf    (format,32, "%%.%ldf",(long)a2);
        }
    } else if (a1>=0) {
        snprintf    (format,32, "%%%ldf",(long)a1);
    } else if (a2>=0) {
        snprintf    (format,32, "%%.%ldf",(long)a2);
    } else {
        snprintf    (format,32, "%%g");
    }

#endif
    a1 = snprintf    (buffer,256, format,Value());
    _String    t (buffer);
    return     new _FString (t);
}

//______________________________________________________________________________
_PMathObj _Constant::Log (void)
{
    return     new _Constant  (log(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Sqrt (void)
{
    return     new _Constant  (sqrt(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Arctan (void)
{
    return     new _Constant  (atan(theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Gamma (void) {
    return new _Constant (_Gamma (theValue));
}

//______________________________________________________________________________
_Parameter _Constant::_Gamma (_Parameter val)
{
    _Parameter theV = val>=1.0?val:2.-val, result = gammaCoeff[0], temp = theV;

    for (long i=1; i<7; i++, temp+=1.0) {
        result+=gammaCoeff[i]/temp;
    }

    temp = theV+4.5;
    result *= exp(-temp+log(temp)*(theV-.5));

    if (val>=1.0) {
        return    result;
    }
    temp = pi_const*(1.-val);
    return     temp/result/sin(temp);
}

//______________________________________________________________________________
_Parameter _Constant::_LnGamma (_Parameter val)
{
    // obtained from Numerical Recipes in C, p. 214 by afyp, February 7, 2007
    _Parameter  tmp, ser;

    tmp = val + 5.5;
    tmp -= (val+0.5) * log(tmp);
    ser = 1.000000000190015 +
          lngammaCoeff [0] / (val + 1.0) +
          lngammaCoeff [1] / (val + 2.0) +
          lngammaCoeff [2] / (val + 3.0) +
          lngammaCoeff [3] / (val + 4.0) +
          lngammaCoeff [4] / (val + 5.0) +
          lngammaCoeff [5] / (val + 6.0);
          

    /*for (long j = 0; j <= 5; j++) {
        ser += lngammaCoeff[j] / ++y;
    }*/

    return -tmp + log(2.506628274631005*ser/val);
}

//______________________________________________________________________________
_PMathObj _Constant::LnGamma (void) {
    return new _Constant (_LnGamma (theValue));
}

//______________________________________________________________________________
_PMathObj _Constant::Beta (_PMathObj arg)
{
    if (arg->ObjectClass()!=NUMBER) {
        WarnError ("A non-numerical argument passed to Beta(x,y)");
        return    nil;
    }
    
    _Parameter  y  = arg->Value(),
                lgX = _LnGamma (theValue),
                lgY = _LnGamma (y),
                lgXY = _LnGamma (theValue + y);
    
    
    return new _Constant (exp (lgX + lgY - lgXY));
}

//______________________________________________________________________________
_PMathObj _Constant::IBeta (_PMathObj arg1, _PMathObj arg2) {
  if (theValue<=0.0) {
    if (theValue < 0.0) {
      ReportWarning   (_String ("IBeta is defined for x betweeen 0 and 1. Had: ") & theValue);
    }
    return new _Constant (0.0);
  }
  
  if (theValue>=1.0) {
    if (theValue>1.0) {
      ReportWarning   (_String ("IBeta is defined for x betweeen 0 and 1. Had: ") & theValue);
    }
    return new _Constant (1.0);
  }
  
  
  if ((arg1->ObjectClass()!=NUMBER)||(arg2->ObjectClass()!=NUMBER)) {
    WarnError   ("IBeta called with a non-scalar argument.");
    return      nil;
  }
  
  _Parameter       a = arg1->Value(),
  b = arg2->Value(),
  ga = _LnGamma(a),
  gb = _LnGamma(b),
  x = theValue,
  aa,
  c,
  d,
  del,
  h,
  qab,
  qam,
  qap,
  FPMIN = 1e-100;
  
  bool        swap = false;
  
  long        m,
  m2;
  
  if (x >= (a+1.)/(a+b+2.)) {
    swap = true;
    c = b;
    b = a;
    a = c;
    x = 1. - x;
  }
  
  qab = a+b;
  qap = a+1.;
  qam = a-1.;
  c   = 1.;
  d   = 1. - qab*x/qap;
  if  ((d<FPMIN)&&(d>-FPMIN)) {
    d = FPMIN;
  }
  d   = 1./d;
  h   = d;
  
  for (m=1; m<100; m++) {
    m2 = 2*m;
    aa = m*(b-m)*x / ((qam+m2)*(a+m2));
    d = 1.+aa*d;
    if  ((d<FPMIN)&&(d>-FPMIN)) {
      d = FPMIN;
    }
    c = 1.+aa/c;
    if  ((c<FPMIN)&&(c>-FPMIN)) {
      c = FPMIN;
    }
    d = 1./d;
    h*= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.+aa*d;
    if  ((d<FPMIN)&&(d>-FPMIN)) {
      d = FPMIN;
    }
    c = 1.+aa/c;
    if  ((c<FPMIN)&&(c>-FPMIN)) {
      c = FPMIN;
    }
    d = 1./d;
    del = d*c;
    h*= del;
    del -= 1.;
    if  ((del<1.e-14)&&(del>-1.e-14))   {
      break;
    }
  }
  
  d  = _LnGamma(a+b);
  c   = exp (a*log(x)+b*log(1-x)+d-ga-gb);
  
  if (swap) {
    return new _Constant (1.-c*h/a);
  } else {
    return new _Constant (c*h/a);
  }
  
  return nil;
}

//______________________________________________________________________________
_Parameter _Constant::_IGamma (_Parameter v, _Parameter x) {
   _Parameter sum = 0.;
   if (x<=v+1) // use the series representation
        // IGamma (a,x)=exp(-x) x^a \sum_{n=0}^{\infty} \frac{\Gamma((a)}{\Gamma(a+1+n)} x^n
    {
        _Parameter term = 1.0/v, den = v+1.;
        long count = 0;
        while (fabs(term)>=fabs(sum)*machineEps && count<500) {
            sum+=term;
            term*=x/den;
            den += 1.0;
            count++;
        }
    } else // use the continue fraction representation
        // IGamma (a,x)=exp(-x) x^a 1/x+/1-a/1+/1/x+/2-a/1+/2/x+...
    {
        _Parameter lastTerm = 0, a0 = 1.0, a1 = x, b0 = 0.0, b1 = 1.0, factor = 1.0, an, ana, anf;
        for (long count = 1; count<500; count++) {
            an = count;
            ana = an - v;
            a0 = (a1+a0*ana)*factor;
            b0 = (b1+b0*ana)*factor;
            anf = an*factor;
            a1  = x*a0+anf*a1;
            b1  = x*b0+anf*b1;
            if (a1!=0.0) {
                factor=1.0/a1;
                sum = b1*factor;
                if (fabs(sum-lastTerm)/sum<machineEps) {
                    break;
                }
                lastTerm = sum;
            }

        }
    }
    sum = exp (-x + v * log (x)) / _Gamma (v);
    if (x > v + 1.) {
      sum = 1. - sum;
    }
    return sum;
}


//______________________________________________________________________________
_PMathObj _Constant::IGamma (_PMathObj arg)
{
    if (arg->ObjectClass()!=NUMBER) {
        WarnError ("A non-numerical argument passed to IGamma(a,x)");
        return new _Constant (0.0);
    }
    _Parameter x = arg->Value();
    if (x>1e25) {
        x=1e25;
    } else if (x<0) {
        WarnError ("The domain of x is {x>0} for IGamma (a,x)");
        return new _Constant (0.0);
    } else if (x==0.0) {
        return new _Constant (0.0);
    }


    
    return new _Constant (_IGamma (theValue, x));
}

//______________________________________________________________________________
_PMathObj _Constant::Erf (void)
{
    _Parameter IG = _IGamma(0.5, theValue * theValue);
               
    if (theValue<0) {
        return new _Constant (-IG);
    }
    return new _Constant (IG);
}

//______________________________________________________________________________
_PMathObj _Constant::ZCDF (void) {
  _Parameter IG = _IGamma (0.5, theValue*theValue*0.5) * 0.5;
 
  if (theValue>0.) {
    return new _Constant (IG + 0.5);
  } 
  return new _Constant (0.5-IG);
  
}

//______________________________________________________________________________
_PMathObj _Constant::Time (void) {
    _Parameter tval;
    if (theValue<1.0) {
        tval = ((_Parameter)clock()/CLOCKS_PER_SEC);
    } else {
        time_t tt;
        tval = ((_Parameter)time(&tt));
    }
    return     new _Constant (tval);
}



//______________________________________________________________________________
_PMathObj _Constant::GammaDist (_PMathObj alpha, _PMathObj beta) {

    _Parameter x = theValue, 
               a = alpha->Value(),
               b = beta->Value(), 
               gd = exp(a * log(b) -b*x +(a-1)*log(x));
               
    return new _Constant (gd / _Gamma (a));
}

//______________________________________________________________________________
_PMathObj _Constant::CGammaDist (_PMathObj alpha, _PMathObj beta) {
    return new _Constant (_IGamma (alpha->Value(), theValue*beta->Value()));
}

//______________________________________________________________________________
_PMathObj _Constant::CChi2 (_PMathObj n)
// chi^2 n d.f. probability up to x
{
    _Parameter halfn = n->Value() * 0.5;

    if (theValue < 0. || halfn <= 0.) {
        ReportWarning ("CChi2(x,n) only makes sense for both arguments positive");
        return new _Constant (0.0);
    }
    return new _Constant (_IGamma (halfn, theValue * 0.5));
}

//______________________________________________________________________________
_PMathObj _Constant::InvChi2 (_PMathObj n)
// chi^2 n d.f. probability up to x
{
    if (!chi2) {
        _String fla ("IGamma(_n_,_x_)");
        chi2 = new _Formula (fla, nil);
        fla = "_x_^(_n_-1)/Gamma(_n_)/Exp(_x_)";
        derchi2 = new _Formula (fla,nil);
    }
    _Constant * halfn = new _Constant (n->Value()*.5);
    if ( theValue<0. || halfn->theValue<0. || theValue>1.) {
        ReportWarning ("InvChi2(x,n) only makes sense for n positive, and x in [0,1]");
        DeleteObject (halfn);
        return new _Constant (0.0);
    }
    LocateVar(dummyVariable2)->SetValue (halfn);
    halfn->SetValue(chi2->Newton(*derchi2,theValue,1e-25,1.e100,LocateVar(dummyVariable1))*2);
    return halfn;
}


//______________________________________________________________________________
_PMathObj _Constant::Less (_PMathObj theObj){
    return new _Constant (theValue<theObj->Value());
}

//______________________________________________________________________________
_PMathObj _Constant::Greater (_PMathObj theObj) {
    return new _Constant (theValue>theObj->Value());
}
//______________________________________________________________________________
_PMathObj _Constant::LessEq (_PMathObj theObj) {
   return new _Constant (theValue<=theObj->Value());
}

//______________________________________________________________________________
_PMathObj _Constant::GreaterEq (_PMathObj theObj) {
   return new _Constant (theValue>=theObj->Value());
}

  //______________________________________________________________________________
_PMathObj _Constant::AreEqual (_PMathObj theObj)
{
  if (!theObj) {
    return nil;
  }
  
  _Parameter a = theValue,
  b = theObj->Value();
  
  if (a==0.0) {
    return new _Constant (b==0.0);
  }
  
  return new _Constant(fabs ((a-b)/a)<tolerance);
}

  //______________________________________________________________________________
_PMathObj _Constant::NotEqual (_PMathObj theObj)
{
  if (!theObj) {
    return nil;
  }
  _Parameter a = theValue,
  b = theObj->Value();
  
  if (a==0.0) {
    return new _Constant (b!=0.0);
  }
  
  return new _Constant(fabs ((a-b)/a)>=tolerance);
}

//______________________________________________________________________________
_PMathObj _Constant::LAnd (_PMathObj theObj) {
  return new _Constant ( (!CheckEqual (theValue, 0.0) && !CheckEqual (theObj->Value(), 0.0)) ? 1.0 : 0.0);
}

//______________________________________________________________________________
_PMathObj _Constant::LOr (_PMathObj theObj) {
  return new _Constant ( (!CheckEqual (theValue, 0.0) || !CheckEqual (theObj->Value(), 0.0)) ? 1.0 : 0.0);
}

//______________________________________________________________________________
_PMathObj _Constant::LNot () {
    return new _Constant (CheckEqual(theValue, 0.0) ? 1.0 : 0.0);
}

//______________________________________________________________________________
_PMathObj _Constant::Min (_PMathObj theObj)
{
   _Parameter ov = theObj->Value();
  return new _Constant (theValue < ov? theValue : ov);
}

//______________________________________________________________________________
_PMathObj _Constant::Max (_PMathObj theObj)
{
   _Parameter ov = theObj->Value();
  return new _Constant (theValue > ov? theValue : ov);
}

//______________________________________________________________________________
_PMathObj _Constant::Execute(
    long opCode, _PMathObj p, _PMathObj p2,
    _hyExecutionContext *
        context, _PMathObj p3) // execute this operation with the second arg if necessary
    {
  
  switch (opCode) {
    case HY_OP_CODE_REF:
      return LocateVar (Value()) -> ComputeReference (context->GetContext());
      
  }
  return _MathObject::Execute (opCode, p, p2, context, p3);
}
