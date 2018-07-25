/*
 
 HyPhy - Hypothesis Testing Using Phylogenies.
 
 Copyright (C) 1997-now
 Core Developers:
 Sergei L Kosakovsky Pond (sergeilkp@icloud.com)
 Art FY Poon    (apoon@cfenet.ubc.ca)
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
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include "string.h"
#include "stdlib.h"
#include "time.h"

#include "constant.h"

//SW: This should be just helper functions
#include "parser.h"

_Formula *chi2 = nil,
         *derchi2 = nil;


extern _Parameter machineEps;
extern _Parameter tolerance;

long            lastMatrixDeclared = -1,
                dummyVariable1,
                dummyVariable2,
                expressionsParsed = 0;





//__________________________________________________________________________________


_Parameter _gamma (_Parameter alpha) {
  _Parameter static gammaCoeff [7] = {
    2.50662827463100050,
    190.9551718944012,
    -216.8366818451899,
    60.19441758801798,
    -3.087513097785903,
    0.003029460875352382,
    -0.00001345152485367085
  };
  
  _Parameter theV = alpha >=1.0? alpha : 2.-alpha,
  result = gammaCoeff[0],
  temp = theV;
  
  for (int i = 1; i < 7; ++i , temp += 1.) {
    result += gammaCoeff[i] / temp;
  }
  
  temp = theV + 4.5;
  result *= exp(-temp+log(temp)*(theV-.5));
  
  if (alpha >= 1.0) {
    return result;
  }
  temp = pi_const * (1-alpha);
  return temp / result / sin (temp);
}

  //__________________________________________________________________________________


_Parameter _ln_gamma (_Parameter alpha) {
  // obtained from Numerical Recipes in C, p. 214 by afyp, February 7, 2007

  _Parameter static lngammaCoeff [6] = {
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  };
  
  
  _Parameter  x, y, tmp, ser;
  
  y = x = alpha;
  tmp = x + 5.5;
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015;
  
  for (int j = 0; j < 6 ; ++j ) {
    ser += lngammaCoeff[j] / ( y += 1. );
  }
  
  return -tmp + log(2.506628274631005*ser/x);

 
}

//__________________________________________________________________________________

_Parameter _ibeta (_Parameter x, _Parameter a, _Parameter b) {
    // check ranges
  if (x > 0. && x < 1.) { // in range
    
 
    _Parameter  aa,
                c,
                d,
                del,
                h,
                qab,
                qam,
                qap,
                FPMIN = 1e-100;

    
    bool        swap = false;
    
    
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
    
    for (int m=1; m<100; m++) {
      _Parameter m2 = 2*m;
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
    
    c = exp (a*log(x)+b*log(1.-x)+_ln_gamma(a+b)-_ln_gamma(a)-_ln_gamma(b));
    
    if (swap) {
      return 1.-c*h/a;
    } else {
      return c*h/a;
    }
  }
  
  
  if (x <= 0.) {
    if ( x < 0.) {
      ReportWarning (_String ("IBeta is defined for x in [0,1]. Had x = ") & x);
    }
    return 0.;
  }
  
  if ( x > 1.) {
    ReportWarning (_String ("IBeta is defined for x in [0,1]. Had x = ") & x);
  }
  
  return 1.;
}

//__________________________________________________________________________________

_Parameter _igamma (_Parameter a, _Parameter x) {
  _Parameter sum = 0.;
  if (x>1e25) {
    x=1.e25;
  } else if (x<0.) {
    WarnError ("The domain of x is {x>0} for IGamma (a,x)");
    return 0.0;
  } else if (x==0.0) {
    return 0.0;
  }
  
  
  _Parameter gamma = _gamma (a);
  
  if (x <= a + 1.) {
      // use the series representation
      // IGamma (a,x)=exp(-x) x^a \sum_{n=0}^{\infty} \frac{\Gamma((a)}{\Gamma(a+1+n)} x^n
    
    _Parameter term = 1.0/a, den = a+1.;
    
    for (int count = 0; fabs (term) >= fabs (sum) * machineEps && count < 500; ++ count) {
      sum+=term;
      term*=x/den;
      den += 1.0;
    }
    
    return sum * exp (-x + a * log (x)) / gamma;
    
  }
    // use the continue fraction representation
    // IGamma (a,x)=exp(-x) x^a 1/x+/1-a/1+/1/x+/2-a/1+/2/x+...
  
  _Parameter lastTerm = 0., a0 = 1.0, a1 = x, b0 = 0.0, b1 = 1.0, factor = 1.0, an, ana, anf;
  for (int count = 1; count<500; ++count) {
    an = count;
    ana = an - a;
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
  
  return 1.0 - sum * exp (-x + a * log (x)) / gamma;
  
}

//__________________________________________________________________________________
_Constant::_Constant (_Parameter value)
{
    theValue = value;
}
//__________________________________________________________________________________

void _Constant::Initialize (bool)
{
    BaseObj::Initialize();
    theValue = 0;
}
//__________________________________________________________________________________

void _Constant::Duplicate (BaseRef c)
{
    BaseObj::Initialize();
    theValue = ((_Constant*)c)->theValue;
}

//__________________________________________________________________________________

BaseRef _Constant::makeDynamic (void)
{
    _Constant * res = (_Constant*)checkPointer(new _Constant);
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________

_Constant::_Constant (_String& s)
{
    theValue = atof (s.sData);
}

//__________________________________________________________________________________
_Constant::_Constant (void)
{
    theValue = 0;
}

//__________________________________________________________________________________
//_Constant::~_Constant (void)  {
//}

//__________________________________________________________________________________
_Parameter    _Constant::Value (void) {
    return theValue;
}

//__________________________________________________________________________________
BaseRef _Constant::toStr(unsigned long)
{
    return parameterToString(Value());
}

//__________________________________________________________________________________
_PMathObj _Constant::Add (_PMathObj theObj) {
  
    if (theObj->ObjectClass() == STRING) {
        return new _Constant ((theValue+((_FString*)theObj)->theString->toNum()));
    }

    return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {return a + b;});
}

//__________________________________________________________________________________
_PMathObj _Constant::Sub (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {return a - b;});
}

//__________________________________________________________________________________
_PMathObj _Constant::Minus (void) {
    return     new  _Constant (-Value());
}

//__________________________________________________________________________________
_PMathObj _Constant::Sum (void) {
    return     new  _Constant (Value());
}

//__________________________________________________________________________________
_PMathObj _Constant::Mult (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {return a * b;});

}
//__________________________________________________________________________________
_PMathObj _Constant::Div (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {return a / b;});
}
//__________________________________________________________________________________
_PMathObj _Constant::lDiv (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {
    long       denom = b;
    return     denom != 0L ? (long(a) % denom): a;
  });
}
//__________________________________________________________________________________
_PMathObj _Constant::longDiv (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {
    long       denom = b;
    return     denom != 0L ? (long(a) / denom): 0.0;
  });
}

//__________________________________________________________________________________
_PMathObj _Constant::Raise (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter base, _Parameter expon) -> _Parameter {
    if (base>0.0) {
      return    exp (log(base)*(expon));
    } else {
      if (base<0.0) {
        if (CheckEqual (expon, (long)expon)) {
          return ((((long)expon)%2)?-1:1)*exp (log(-base)*(expon));
        } else {
          _String errMsg ("An invalid base/exponent pair passed to ^");
          WarnError (errMsg.sData);
        }
      }
      
      if (expon != 0.0)
        return     0.0;
      else
        return     1.0;
    }
  });
}

//__________________________________________________________________________________
_PMathObj _Constant::Random (_PMathObj theObj) {
 
  return _check_type_and_compute (theObj, [] (_Parameter l, _Parameter u) -> _Parameter {
      _Parameter r = l;
      if (u>l) {
          r=genrand_int32();
          r/=RAND_MAX_32;
          r =l+(u-l)*r;
      }
      return r;
  });
}

//__________________________________________________________________________________
void     _Constant::Assign (_PMathObj theObj)
{
    this->~_Constant ();
    theValue = ((_Constant*)theObj)->theValue;
}

//__________________________________________________________________________________
bool     _Constant::Equal (_PMathObj theObj) {
    return theValue==((_Constant*)theObj)->theValue;
}

//__________________________________________________________________________________
_PMathObj _Constant::Abs (void) {
    return     new _Constant (fabs(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Sin (void) {
    return     new  _Constant (sin(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Cos (void){
    return     new _Constant  (cos(theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Tan (void)
{
    return     new _Constant  (tan(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Exp (void)
{
    return     new _Constant  (exp(theValue));
}
//__________________________________________________________________________________
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
    snprintf    (buffer,256, format,Value());
    return     new _FString (new _String (buffer));
}
//__________________________________________________________________________________
_PMathObj _Constant::Log (void) {
    return     new _Constant  (log(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Sqrt (void) {
    return     new _Constant  (sqrt(theValue));
}
//__________________________________________________________________________________
_PMathObj _Constant::Arctan (void) {
    return     new _Constant  (atan(theValue));
}


//__________________________________________________________________________________

_PMathObj _Constant::Gamma (void) {
  return new _Constant (_gamma (theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::LnGamma (void) {
  return new _Constant (_ln_gamma (theValue));
}

//__________________________________________________________________________________
_PMathObj _Constant::Beta (_PMathObj arg) {
    return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
      return exp (_ln_gamma (a) + _ln_gamma (b) - _ln_gamma (a + b));
    });
}

//__________________________________________________________________________________
_PMathObj _Constant::IBeta (_PMathObj arg1, _PMathObj arg2) {
  return _check_type_and_compute_3 (arg1, arg2, [] (_Parameter a, _Parameter b, _Parameter c) -> _Parameter {
    return _ibeta(a,b,c);
  });
}


//__________________________________________________________________________________
_PMathObj _Constant::IGamma (_PMathObj arg) {
  return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter { return _igamma(a, b);});
}

//__________________________________________________________________________________
_PMathObj _Constant::Erf (void) {
  
    _Parameter ig = _igamma (0.5, theValue * theValue);
    if (theValue < 0.) {
      ig = -ig;
    }
    return new _Constant (ig);
}

//__________________________________________________________________________________
_PMathObj _Constant::ZCDF (void) {
    _Parameter ig = _igamma (0.5, theValue * theValue * 0.5);
  
    if (theValue > 0) {
      return new _Constant (0.5 * (ig + 1.));
    }
    return new _Constant (0.5 * ( 1. - ig));
 }

//__________________________________________________________________________________
_PMathObj _Constant::Time (void) {
    _Constant * result = new _Constant;
    if (theValue<1.0) {
        result->theValue = ((_Parameter)clock()/CLOCKS_PER_SEC);
    } else {
        time_t tt;
        result->theValue = ((_Parameter)time(&tt));
    }
    return     result;
}

//__________________________________________________________________________________
_PMathObj _Constant::Less (_PMathObj theObj) {
    return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {return a < b;});
}

//__________________________________________________________________________________
_PMathObj _Constant::Greater (_PMathObj theObj) {
  return _check_type_and_compute (theObj, [] (_Parameter a, _Parameter b) -> _Parameter {return a > b;});
}

//__________________________________________________________________________________
_PMathObj _Constant::GammaDist (_PMathObj alpha, _PMathObj beta) {
  return _check_type_and_compute_3 (alpha, beta, [] (_Parameter x, _Parameter a, _Parameter b) -> _Parameter {
    _Parameter gamma_dist = exp(a * log(b) -b*x +(a-1.)*log(x));
    return gamma_dist / _gamma (a);
  });
}

//__________________________________________________________________________________
_PMathObj _Constant::CGammaDist (_PMathObj alpha, _PMathObj beta) {
    return _check_type_and_compute_3 (alpha, beta, [] (_Parameter x, _Parameter a, _Parameter b) -> _Parameter {
      return _igamma (a, b * x);
    });
 }

//__________________________________________________________________________________
_PMathObj _Constant::CChi2 (_PMathObj theObj) {
    return _check_type_and_compute (theObj, [] (_Parameter x, _Parameter b) -> _Parameter {
      if (x < 0.0 || b <= 0.) {
        ReportWarning ("CChi2(x,n) only makes sense for both arguments positive");
        return 0.0;
      }
      return _igamma( b*0.5 , x * 0.5);
    });
}

//__________________________________________________________________________________
_PMathObj _Constant::InvChi2 (_PMathObj n) {
// chi^2 n d.f. probability up to x
    if (!chi2) {
        _String fla ("IGamma(_n_,_x_)");
        chi2 = new _Formula (fla, nil);
        fla = "_x_^(_n_-1)/Gamma(_n_)/Exp(_x_)";
        derchi2 = new _Formula (fla,nil);
    }
    _Constant halfn (((_Constant*)n)->theValue*.5);
    if ((theValue<0)||(halfn.theValue<0)||(theValue>1.0)) {
        _String warnMsg ("InvChi2(x,n) only makes sense for n positive, and x in [0,1]");
        ReportWarning (warnMsg);
        return new _Constant (0.0);
    }
    LocateVar(dummyVariable2)->SetValue (&halfn);
    halfn.SetValue(chi2->Newton(*derchi2,theValue,1e-25,1.e100,LocateVar(dummyVariable1))*2);
    return (_PMathObj)halfn.makeDynamic();
}

//__________________________________________________________________________________
_PMathObj _Constant::LessEq (_PMathObj arg) {
    return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter { return a <= b;});
}

//__________________________________________________________________________________
_PMathObj _Constant::GreaterEq (_PMathObj arg) {
  return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter { return a >= b;});
}
//__________________________________________________________________________________
_PMathObj _Constant::AreEqual (_PMathObj arg) {
  return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
    if (a==0.0) {
      return b==0.0;
    }
    
    return fabs ((a-b)/a)<tolerance;
  
  });
  
}

//__________________________________________________________________________________
_PMathObj _Constant::NotEqual (_PMathObj arg) {
  return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
    if (a==0.0) {
      return b!=0.0;
    }
    
    return fabs ((a-b)/a)>=tolerance;
  });
}

//__________________________________________________________________________________
_PMathObj _Constant::LAnd (_PMathObj arg) {
    return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
      return long (a) && long (b);
    });
}

  //__________________________________________________________________________________
_PMathObj _Constant::LOr (_PMathObj arg) {
  return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
    return long (a) || long (b);
  });
}

//__________________________________________________________________________________
_PMathObj _Constant::LNot () {
    return new _Constant (CheckEqual(theValue, 0.0));
}

//__________________________________________________________________________________
_PMathObj _Constant::Min (_PMathObj arg) {
    return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
      return a < b ? a : b;
    });
}

//__________________________________________________________________________________
_PMathObj _Constant::Max (_PMathObj arg) {
  return _check_type_and_compute (arg, [] (_Parameter a, _Parameter b) -> _Parameter {
    return a > b ? a : b;
  });
}
