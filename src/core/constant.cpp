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
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include "string.h"
#include "time.h"

#include "mersenne_twister.h"
#include "constant.h"

//SW: This should be just helper functions
#include "parser.h"
#include "global_things.h"

using namespace hy_global;

_Formula *chi2 = nil,
         *derchi2 = nil;


extern hyFloat tolerance;

long                        lastMatrixDeclared = -1,
                            expressionsParsed = 0;

unsigned char               _Constant::preallocated_buffer [_HY_CONSTANT_PREALLOCATE_SLOTS*sizeof (_Constant)];
_SimpleList                 _Constant::free_slots;

//___________________________________________________________________________________________
hyFloat  gaussDeviate (void) {
  /*
   Use Box-Muller transform to generate random deviates from Gaussian distribution with
   zero mean and unit variance (Numerical Recipes).
   */
  
  static int      iset = 0;
  static double   gset;
  double          fac, rsq, v1, v2;
  
  if (iset == 0) {
    do {
      v1 = 2.0 * genrand_real2() - 1.0;   // uniform random number on (0,1), i.e. no endpoints
      v2 = 2.0 * genrand_real2() - 1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    
    fac = sqrt(-2.0 * log(rsq)/rsq);
    gset = v1 * fac;
    iset = 1;           // set flag to indicate second deviate available
    
    return (hyFloat) (v2 * fac);
  } else {
    iset = 0;
    return (hyFloat) gset;       // use second deviate
  }
}


//__________________________________________________________________________________________________________
hyFloat  exponDeviate (void) {
  return -log(1.0-genrand_real2());   // uniform random number on interval (0,1]
}


//__________________________________________________________________________________________________________
hyFloat  gammaDeviate (hyFloat a, hyFloat scale) {
  
  /* -----------------------------------------
   GS algorithm from GNU GPL rgamma.c
   *  Mathlib : A C Library of Special Functions
   *  Copyright (C) 1998 Ross Ihaka
   *  Copyright (C) 2000--2008 The R Development Core Team
   
   * Ziggurat algorithm from Marsaglia and Tsang (2000) "A Simple Method
   *  for Generating Gamma Variables" ACM Trans Math Soft 26(3) 363-372
   ----------------------------------------- */
  
  const static double exp_m1 = 0.36787944117144232159;    /* exp(-1) = 1/e */
  
  double              e, x, p;
  
  if (a < 0.0) {
    ReportWarning ("NaN in gammaDeviate()");
    return 0.;
  } else if (a == 0.0) {
    return 0.;
  } else if (a < 1.0) {   // GS algorithm for parameters 0 < a < 1
    if(a == 0) {
      return 0.;
    }
    
    e = 1.0 + exp_m1 * a;
    
    while (1) {
      p = e * genrand_real2();    // should be uniform random number on open interval (0,1)
                                  // but genrand_real3 scoping (baseobj.cpp) is insufficient
      if (p >= 1.0) {
        x = -log((e - p) / a);
        if (exponDeviate() >= (1.0 - a) * log(x)) {
          break;
        }
      } else {
        x = exp(log(p) / a);
        if (exponDeviate() >= x) {
          break;
        }
      }
    }
    
    return x*scale;
  }
  
  else if (a == 1.0) {
    return exponDeviate() * scale;
  }
  
  else {  // a > 1., Ziggurat algorithm
    double  x, v, u,
    d   = a - 1./3.,
    c = 1. / sqrt(9.*d);
    
    for (;;) {
      do {
        x = gaussDeviate();
        v = 1. + c * x;
      } while (v <= 0.);
      
      v = v * v * v;
      u = genrand_real2();
      
      if (u < 1. - 0.0331 * (x*x)*(x*x) ) {
        return (d * v * scale);
      }
      
      if ( log(u) < 0.5*x*x + d*(1.-v+log(v)) ) {
        return (d * v * scale);
      }
    }
  }
}



//__________________________________________________________________________________
hyFloat  chisqDeviate (double df) {
  if (df < 0.0) {
    HandleApplicationError (_String("ERROR in chisqDeviate(): require positive degrees of freedom"));
    return HY_INVALID_RETURN_VALUE;
  }
  
  return gammaDeviate(df/2.0, 2.0);   // chi-square distribution is special case of gamma
}



//__________________________________________________________________________________


hyFloat _gamma (hyFloat alpha) {
    hyFloat static gammaCoeff [7] = {
        2.50662827463100050,
        190.9551718944012,
        -216.8366818451899,
        60.19441758801798,
        -3.087513097785903,
        0.003029460875352382,
        -0.00001345152485367085
    };
    
    hyFloat theV = alpha >=1.0? alpha : 2.-alpha,
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


hyFloat _ln_gamma (hyFloat alpha) {
    // obtained from Numerical Recipes in C, p. 214 by afyp, February 7, 2007
    
    if (alpha <= 0.) {
        return NAN;
    }
    
    hyFloat static lngammaCoeff [6] = {
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };
  
    static hyFloat lookUpTable [20] = {   0.       ,  0.       ,  0.6931472,  1.7917595,  3.1780538,
      4.7874917,  6.5792512,  8.5251614, 10.6046029, 12.8018275,
      15.1044126, 17.5023078, 19.9872145, 22.5521639, 25.1912212,
      27.8992714, 30.6718601, 33.5050735, 36.3954452, 39.3398842
    };
  
    if (alpha <= 20. && floor (alpha) == alpha) {
      return lookUpTable [(long) alpha - 1];
    }

    
    hyFloat  x, y, tmp, ser;
    
    y = x = alpha;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
  
    ser = 1.000000000190015 +
          lngammaCoeff[0] / (y + 1.) +
          lngammaCoeff[1] / (y + 2.) +
          lngammaCoeff[2] / (y + 3.) +
          lngammaCoeff[3] / (y + 4.) +
          lngammaCoeff[4] / (y + 5.) +
          lngammaCoeff[5] / (y + 6.);

    /*
    for (int j = 0; j < 6 ; ++j ) {
        ser += lngammaCoeff[j] / ( y += 1. );
    }
    */
    
    return -tmp + log(2.506628274631005*ser/x);
    
    
}

//__________________________________________________________________________________

hyFloat _ibeta (hyFloat x, hyFloat a, hyFloat b) {
    // check ranges
    if (x > 0. && x < 1.) { // in range
        
        
        hyFloat  aa,
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
            hyFloat m2 = 2*m;
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

hyFloat _igamma (hyFloat a, hyFloat x) {
    hyFloat sum = 0.;
    if (x>1e25) {
        x=1.e25;
    } else if (x<0.) {
        HandleApplicationError ("The domain of x is {x>0} for IGamma (a,x)");
        return 0.0;
    } else if (x==0.0) {
        return 0.0;
    }
    
    
    hyFloat gamma = _gamma (a);
    
    if (x <= a + 1.) {
        // use the series representation
        // IGamma (a,x)=exp(-x) x^a \sum_{n=0}^{\infty} \frac{\Gamma((a)}{\Gamma(a+1+n)} x^n
        
        hyFloat term = 1.0/a, den = a+1.;
        
        for (int count = 0; fabs (term) >= fabs (sum) * kMachineEpsilon && count < 500; ++ count) {
            sum+=term;
            term*=x/den;
            den += 1.0;
        }
        
        return sum * exp (-x + a * log (x)) / gamma;
        
    }
    // use the continue fraction representation
    // IGamma (a,x)=exp(-x) x^a 1/x+/1-a/1+/1/x+/2-a/1+/2/x+...
    
    hyFloat lastTerm = 0., a0 = 1.0, a1 = x, b0 = 0.0, b1 = 1.0, factor = 1.0, an, ana, anf;
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
            if (fabs(sum-lastTerm)/sum < kMachineEpsilon) {
                break;
            }
            lastTerm = sum;
        }
    }
    
    return 1.0 - sum * exp (-x + a * log (x)) / gamma;
    
}

//__________________________________________________________________________________
_Constant::_Constant (hyFloat value) {
    theValue = value;
}
//__________________________________________________________________________________

void _Constant::Initialize (bool) {
    BaseObj::Initialize();
    theValue = 0.;
}
//__________________________________________________________________________________

void _Constant::Duplicate (BaseRefConst c) {
    BaseObj::Initialize();
    theValue = ((_Constant const*)c)->theValue;
}

//__________________________________________________________________________________

BaseRef _Constant::makeDynamic (void) const{
    _Constant * res = new _Constant;
    res->Duplicate(this);
    return res;
}

//__________________________________________________________________________________
void * _Constant::operator new (size_t size) {
    if (_Constant::free_slots.nonempty()) {
        _Constant * result = ((_Constant*)_Constant::preallocated_buffer)+_Constant::free_slots.Pop();
        return result;
    }
    return ::operator new (size);
}

//__________________________________________________________________________________
void  _Constant::operator delete (void * p) {
    _Constant * cp = (_Constant*)p,
              * cpp = (_Constant*)_Constant::preallocated_buffer;
    if (cp >= cpp &&  cp < (cpp + _HY_CONSTANT_PREALLOCATE_SLOTS)) {
        free_slots << long (cp - cpp);
    } else {
        ::operator delete (p);
    }
}

//__________________________________________________________________________________

_Constant::_Constant (_String& s) {
    theValue = s.to_float();
}

//__________________________________________________________________________________
_Constant::_Constant (void) : theValue(0.0) {
}

//__________________________________________________________________________________
//_Constant::~_Constant (void)  {
//}

//__________________________________________________________________________________
hyFloat    _Constant::Value (void) {
    return theValue;
}
//__________________________________________________________________________________
BaseRef _Constant::toStr(unsigned long) {
    return parameterToString(Value());
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Add (HBLObjectRef theObj, HBLObjectRef cache) {
    if (theObj->ObjectClass() == STRING) {
        return _returnConstantOrUseCache (theValue+((_FString*)theObj)->get_str().to_float(), cache);
    }
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {return a + b;}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Sub (HBLObjectRef theObj, HBLObjectRef cache) {
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {return a - b;}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Minus (HBLObjectRef cache) {
    return     _returnConstantOrUseCache (-Value(), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Sum (HBLObjectRef cache) {
    return     _returnConstantOrUseCache (Value(), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Mult (HBLObjectRef theObj, HBLObjectRef cache) {
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {return a * b;}, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::Div (HBLObjectRef theObj,  HBLObjectRef cache) {
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {return a / b;}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::lDiv (HBLObjectRef theObj, HBLObjectRef cache) { // %
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {
            long       denom = b;
            return     denom != 0L ? (long(a) % denom): a;
    }, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::longDiv (HBLObjectRef theObj, HBLObjectRef cache)  {// div
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {
        long       denom = b;
        return     denom != 0L ? (long(a) / denom): 0.0;
    }, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::Raise (HBLObjectRef theObj, HBLObjectRef cache) {
    return _check_type_and_compute (theObj, [] (hyFloat base, hyFloat expon) -> hyFloat {
        if (base>0.0) {
            if (expon == 1.) {
                return base;
            }
            return    exp (log(base)*(expon));
        } else {
            if (base<0.0) {
                if (CheckEqual (expon, (long)expon)) {
                    return ((((long)expon)%2)?-1:1)*exp (log(-base)*(expon));
                } else {
                    HandleApplicationError("An invalid base/exponent pair passed to ^");
                }
            }
            
            if (expon != 0.0)
                return     0.0;
            else
                return     1.0;
        }
    }, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Random (HBLObjectRef upperB, HBLObjectRef cache) {
    return _check_type_and_compute (upperB, [] (hyFloat l, hyFloat u) -> hyFloat {
        hyFloat r = l;
        if (u>l) {
            r = l + (u-l) * genrand_real1();
        }
        return r;
    }, cache);
}

//__________________________________________________________________________________
bool     _Constant::Equal (HBLObjectRef theObj) {
    return theValue==((_Constant*)theObj)->theValue;
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Abs (HBLObjectRef cache) {
    return     _returnConstantOrUseCache(fabs(theValue), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Sin (HBLObjectRef cache){
    return     _returnConstantOrUseCache(sin(theValue), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Cos (HBLObjectRef cache){
    return     _returnConstantOrUseCache(cos(theValue) , cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Tan (HBLObjectRef cache){
    return     _returnConstantOrUseCache(tan(theValue) , cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Exp (HBLObjectRef cache){
    return     _returnConstantOrUseCache(exp(theValue) , cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::FormatNumberString (HBLObjectRef p, HBLObjectRef p2, HBLObjectRef cache) {
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
    return _returnStringOrUseCache(buffer, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::Log (HBLObjectRef cache) {
    return     _returnConstantOrUseCache(log(theValue), cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::Sqrt (HBLObjectRef cache) {
    return     _returnConstantOrUseCache(sqrt (theValue), cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::Arctan (HBLObjectRef cache) {
    return     _returnConstantOrUseCache(atan (theValue), cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::Gamma (HBLObjectRef cache) {
    return     _returnConstantOrUseCache(_gamma (theValue), cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::LnGamma (HBLObjectRef cache) {
    return     _returnConstantOrUseCache(_ln_gamma (theValue), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Beta (HBLObjectRef arg,HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
        return exp (_ln_gamma (a) + _ln_gamma (b) - _ln_gamma (a + b));
    }, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::IBeta (HBLObjectRef arg1, HBLObjectRef arg2,HBLObjectRef cache) {
    return _check_type_and_compute_3 (arg1, arg2, [] (hyFloat a, hyFloat b, hyFloat c) -> hyFloat {
        return _ibeta(a,b,c);
    },cache);
}


//__________________________________________________________________________________
HBLObjectRef _Constant::IGamma (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat { return _igamma(a, b);}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Erf (HBLObjectRef cache) {
    hyFloat ig = _igamma (0.5, theValue * theValue);
    if (theValue < 0.) {
        ig = -ig;
    }
    return _returnConstantOrUseCache (ig, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::ZCDF (HBLObjectRef cache) {
    hyFloat ig = _igamma (0.5, theValue * theValue * 0.5);
    
    if (theValue > 0.) {
        return _returnConstantOrUseCache (0.5 * (ig + 1.), cache);
    }
    return _returnConstantOrUseCache (0.5 * ( 1. - ig), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Time (HBLObjectRef cache) {
    //_Constant * result = new _Constant;
    if (theValue<1.0) {
        return  _returnConstantOrUseCache (((hyFloat)clock()/CLOCKS_PER_SEC), cache);
    } else {
        time_t tt;
        return  _returnConstantOrUseCache  ((hyFloat)time(&tt), cache);
    }
    //return     result;
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Less (HBLObjectRef theObj, HBLObjectRef cache) {
   return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {return a < b;}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Greater (HBLObjectRef theObj, HBLObjectRef cache) {
    return _check_type_and_compute (theObj, [] (hyFloat a, hyFloat b) -> hyFloat {return a > b;}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::GammaDist (HBLObjectRef alpha, HBLObjectRef beta, HBLObjectRef cache) {
    return _check_type_and_compute_3 (alpha, beta, [] (hyFloat x, hyFloat a, hyFloat b) -> hyFloat {
        hyFloat gamma_dist = exp(a * log(b) -b*x +(a-1.)*log(x));
        return gamma_dist / _gamma (a);
    }, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::CGammaDist (HBLObjectRef alpha, HBLObjectRef beta, HBLObjectRef cache) {
    return _check_type_and_compute_3 (alpha, beta, [] (hyFloat x, hyFloat a, hyFloat b) -> hyFloat {
        return _igamma (a, b * x);
    }, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::CChi2 (HBLObjectRef n, HBLObjectRef cache) {
// chi^2 n d.f. probability up to x
    return _check_type_and_compute (n, [] (hyFloat x, hyFloat b) -> hyFloat {
        if (x < 0.0 || b <= 0.) {
            ReportWarning ("CChi2(x,n) only makes sense for both arguments positive");
            return 0.0;
        }
        return _igamma( b*0.5 , x * 0.5);
    }, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::InvChi2 (HBLObjectRef n, HBLObjectRef cache) {
// chi^2 n d.f. probability up to x
    if (!chi2) {
        chi2 = new _Formula (_String ("IGamma(") &  kNVariableName & "," & kXVariableName & ")", nil);
        derchi2 = new _Formula (kXVariableName & "^(" & kNVariableName & "-1)/Gamma(" & kNVariableName & ")/Exp(" & kXVariableName & ")",nil);
    }
    hyFloat half_n = ((_Constant*)n)->theValue*.5;
    if (theValue<0. || half_n < 0.|| theValue> 1.0) {
        ReportWarning ("InvChi2(x,n) is defined for n > 0, and x in [0,1]");
        return _returnConstantOrUseCache(0., cache);
    }
    
    _Constant* result = new _Constant (half_n);
    hy_n_variable->SetValue (result, true, true, NULL);
    result->SetValue(chi2->Newton(*derchi2,theValue,1e-25,1.e100,hy_x_variable)*2);
    return result;
}

//__________________________________________________________________________________
HBLObjectRef _Constant::LessEq (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat { return a <= b;}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::GreaterEq (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat { return a >= b;}, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::AreEqual (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
        if (a==0.0) {
            return b==0.0;
        }
        return fabs ((a-b)/a)<tolerance;
    }, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::NotEqual (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
        if (a==0.0) {
            return b!=0.0;
        }
        
        return fabs ((a-b)/a)>=tolerance;
    }, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::LAnd (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
        return long (a) && long (b);
    }, cache);
}
//__________________________________________________________________________________
HBLObjectRef _Constant::LOr (HBLObjectRef arg,HBLObjectRef cache) {
return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
    return long (a) || long (b);
}, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::LNot (HBLObjectRef cache) {
    return _returnConstantOrUseCache(CheckEqual(theValue, 0.0), cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Min (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
        return a < b ? a : b;
    }, cache);
}

//__________________________________________________________________________________
HBLObjectRef _Constant::Max (HBLObjectRef arg, HBLObjectRef cache) {
    return _check_type_and_compute (arg, [] (hyFloat a, hyFloat b) -> hyFloat {
        return a > b ? a : b;
    }, cache);
}
