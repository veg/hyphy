
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

#include "batchlan.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#include "polynoml.h"
#include "likefunc.h"

#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

/*SLKP 20110209; include progress report updates */
#if !defined __UNIX__ && !defined __HEADLESS__
#include "HYConsoleWindow.h"
#endif
/*SLKP*/

//#include "profiler.h"

#define MEMORYERROR "Out of Memory"
#define ZEROOBJECT  0.0
#define ZEROPOINTER nil


_String     MATRIX_AGREEMENT            = "CONVERT_TO_POLYNOMIALS",
            ANAL_COMP_FLAG              = "ANALYTIC_COMPUTATIONS",
            ANAL_MATRIX_TOLERANCE      = "ANAL_MATRIX_TOLERANCE",
            PROFILE_MEAN_VAR_MULT      = "PROFILE_MEAN_VAR_MULT",
            CACHE_FORMULA_DEPENDANCY  = "CACHE_FORMULA_DEPENDANCY",
            BRANCH_LENGTH_STENCIL      = "BRANCH_LENGTH_STENCIL",
            AVL_ITERATOR_ORDER          = "INDEXORDER",
            AVL_ITERATOR_ORDER_VALUE    = "VALUEINDEXORDER";

int _Matrix::precisionArg = 0;
int _Matrix::storageIncrement = 16;
//  percent of total size (reasonable values divide 100)
int _Matrix::switchThreshold = 50;

#ifndef     __HYALTIVEC__
_Parameter  _Matrix::truncPrecision = 1e-13;
#define     MatrixMemAllocate(X) MemAllocate(X)
#define     MatrixMemFree(X)     free(X)
#else
#define     MatrixMemAllocate(X) VecMemAllocate(X)
#define     MatrixMemFree(X)     vec_free(X)
extern      vector float VECTOR_ZERO;
_Parameter  _Matrix::truncPrecision = 1e-8;
#endif

_Parameter  analMatrixTolerance = 1e-6,
            zero = 0,
            ANALYTIC_COMPUTATION_FLAG = 0,
            AUTO_PAD_DIAGONAL = 1,
            toPolyOrNot=0.0,
            toMorNot2M=1.0;


_List       builtInMatrixFunctions;

_Matrix     *GlobalFrequenciesMatrix;

long        matrixExpCount = 0,
            taylorTermsCount = 0,
            squaringsCount = 0,
            non0count = 0;

extern      _String         printDigitsSpec;

//__________________________________________________________________________________

int         fexact_                 (long , long , double *, double , double , double , double *, double *);
void        MatrixIndexError        (long, long, long, long);


// function prototypes
_Parameter  lnGamma (_Parameter),
            gammaDeviate (double, double = 1.);


//__________________________________________________________________________________________________________
_Parameter  lnGamma(_Parameter theValue)
{
    //  Returns Log(gamma(x))
    if (theValue <= 0) {
        _String oops ("ERROR (matrix.cpp): Requested lnGamma(x) for x <= 0.");
        WarnError (oops);

        return 0.;
    }

    static _Parameter lngammaCoeff [6] = {   76.18009172947146,
                                         -86.50532032941677,
                                         24.01409824083091,
                                         -1.231739572450155,
                                         0.1208650973866179e-2,
                                         -0.5395239384953e-5
                                         };

    static _Parameter lookUpTable [20] = {   0.       ,  0.       ,  0.6931472,  1.7917595,  3.1780538,
                                         4.7874917,  6.5792512,  8.5251614, 10.6046029, 12.8018275,
                                         15.1044126, 17.5023078, 19.9872145, 22.5521639, 25.1912212,
                                         27.8992714, 30.6718601, 33.5050735, 36.3954452, 39.3398842
                                         };

    // use look-up table for small integer values
    if (theValue <= 20 && (theValue - (long)theValue) == 0.) {
        return (lookUpTable [(long) theValue - 1]);
    }

    // else do it the hard way
    _Parameter  x, y, tmp, ser;

    y = x = theValue;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
    ser = 1.000000000190015;

    for (long j = 0; j <= 5; j++) {
        ser += lngammaCoeff[j] / ++y;
    }

    return (-tmp + log(2.506628274631005*ser/x));
}


//___________________________________________________________________________________________
_Parameter  gaussDeviate (void)
{
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

        return (_Parameter) (v2 * fac);
    } else {
        iset = 0;
        return (_Parameter) gset;       // use second deviate
    }
}


//__________________________________________________________________________________________________________
_Parameter  exponDeviate (void)
{
    return -log(1.0-genrand_real2());   // uniform random number on interval (0,1]
}



//__________________________________________________________________________________________________________
_Parameter  gammaDeviate (double a, double scale)
{
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
_Parameter  chisqDeviate (double df)
{
    if (df < 0.0) {
        WarnError (_String("ERROR in chisqDeviate(): require positive degrees of freedom"));
        return 0;
    }

    return gammaDeviate(df/2.0, 2.0);   // chi-square distribution is special case of gamma
}



//__________________________________________________________________________________

void    MatrixIndexError (long hPos, long vPos, long hDim, long vDim)
{
    _String errMsg ("Invalid Matrix Index [");
    errMsg = errMsg & _String ((long)hPos) & "][" & _String ((long)vPos) &
             "] in a " &_String (hDim) & " by " &_String (vDim) & " matrix.";
    WarnError (errMsg);
}



//_____________________________________________________________________________________________

inline  bool    _Matrix::IsNonEmpty  (long logicalIndex)
{
    return  (theIndex?theIndex [logicalIndex]!=-1:(storageType!=1?GetMatrixObject(logicalIndex)!=ZEROPOINTER:true));
}

//__________________________________________________________________________________

bool        _Matrix::HasChanged(void)
{
    if (storageType == 2) {
        _Formula* theF, **theFormulae = (_Formula**)theData;
        if (theIndex) {
            for (long i=0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    theF = theFormulae[i];
                    if (theF->HasChanged()) {
                        return true;
                    }
                }
            }
        } else {
            for (long i=0; i<lDim; i++) {
                theF = theFormulae[i];
                if ((theF!=(_Formula*)ZEROPOINTER)&&(theF->HasChanged())) {
                    return true;
                }
            }
        }

    } else if (storageType == 0) {
        _MathObject *theO, **objData=(_MathObject**)theData;
        if (theIndex) {
            for (long i=0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    if (objData[i]->HasChanged()) {
                        return true;
                    }
                }
            }
        } else {
            for (long i=0; i<lDim; i++) {
                theO = objData[i];
                if (theO&&theO->HasChanged()) {
                    return true;
                }
            }
        }

    } else if (storageType == 3) {
        for (long fid = 0; fid < cmd->formulasToEval.lLength; fid++)
            if (((_Formula*)cmd->formulasToEval.lData[fid])->HasChangedSimple(cmd->varIndex)) {
                return true;
            }
    }
    return false;
}
//__________________________________________________________________________________


inline static void ROTATE(_Parameter * a, long i, long j, long k, long l, _Parameter & g, _Parameter & h, _Parameter s, _Parameter tau, long hDim)
{
    g = a[i*hDim + j];
    h = a[k*hDim + l];
    a[i*hDim + j] = g - s*(h + g*tau);
    a[k*hDim + l] = h + s*(g - h*tau);
}

//__________________________________________________________________________________
void        _Matrix::Balance (void)
{
    if (storageType!=1 || hDim!=vDim || hDim==0) { // only works for numerical matrices at this stage
        _String errorMsg ("Balance only works with numerical non-empty square dense matrices");
        WarnError (errorMsg);
        return;
    }

    _Parameter       Squared_Radix = 2.0 * 2.0;

    bool             done = false;

    while (!done) {
        done = true;

        for (long i = 0; i < hDim; i++) {
            _Parameter r = 0.0,
                       c = 0.0;

            for (long j = 0; j < vDim; j++)
                if (i!=j) {
                    r += fabs (theData[i*vDim+j]);
                    c += fabs (theData[j*vDim+i]);
                }

            if (r > 0.0 && c > 0.0) {
                _Parameter g = r / Squared_Radix,
                           f = 1.,
                           s = c+r;

                while (c<g) {
                    f *= 2.0;
                    c *= Squared_Radix;
                }

                g = r * 2.0;

                while (c>g) {
                    f /= 2.0;
                    c /= Squared_Radix;
                }

                if ((c+r)/f < 0.95*s) {
                    done = false;
                    g = 1. / f;
                    for (long j = 0; j < vDim; j++) {
                        theData[i*vDim+j] *= g;
                        theData[j*vDim+i] *= f;
                    }
                }

            }
        }
    }
}

//__________________________________________________________________________________
void        _Matrix::Schur (void)
{
    if (storageType!=1 || hDim!=vDim || hDim==0) { // only works for numerical matrices at this stage
        _String errorMsg ("Hessenberg only works with numerical non-empty square dense matrices");
        WarnError (errorMsg);
        return;
    }

    for (long m = 1; m < hDim-1; m++) {
        _Parameter x = 0.0;
        long       i = m;

        for (long j = m; j < hDim; j++)
            if (fabs (theData[j*vDim + m-1]) > x) {
                x = theData[j*vDim + m-1];
                i = j;
            }

        if (i!=m) {
            for (long j=m-1; j<hDim; j++) {
                _Parameter t = theData[i*vDim + j];
                theData[i*vDim + j] = theData[m*vDim + j];
                theData[m*vDim + j] = t;
            }
            {
                for (long j=0; j<hDim; j++) {
                    _Parameter t = theData[j*vDim + i];
                    theData[j*vDim + i] = theData[j*vDim + m];
                    theData[j*vDim + m] = t;
                }
            }
        }

        if (x)
            for (long i = m+1; i < hDim; i++) {
                _Parameter y = theData[i*vDim + m -1];
                if (y != 0.0) {
                    y /= x;
                    theData[i*vDim + m -1] = y;
                    for (long j = m; j < hDim; j++) {
                        theData[i*vDim+j] -= y*theData[m*vDim+j];
                    }
                    {
                        for (long j = 0; j < hDim; j++) {
                            theData[j*vDim+m] += y*theData[j*vDim+i];
                        }
                    }
                }
            }
    }

    for (long r = 2; r < hDim; r++)
        for (long c = 0; c<r-1; c++) {
            theData[r*hDim+c] = 0.0;
        }

}

//__________________________________________________________________________________

#define MX_ACCESS(a,b) theData[(a)*hDim+(b)]

//__________________________________________________________________________________
void        _Matrix::EigenDecomp (_Matrix& real, _Matrix & imag)
{
    if (storageType!=1 || hDim!=vDim || hDim==0) { // only works for numerical matrices at this stage
        _String errorMsg ("EigenDecomp only works with numerical non-empty square dense matrices");
        WarnError (errorMsg);
        return;
    }

    _Parameter anorm = 0.0;

    for (long k = 0; k < hDim; k++)
        for (long k2 = k?k-1:0; k2 < hDim; k2++) {
            anorm += fabs (MX_ACCESS(k,k2));
        }

    long        nn = hDim - 1;
    _Parameter  t  = 0;

    CreateMatrix (&real, hDim, 1, false, true, false);
    CreateMatrix (&imag, hDim, 1, false, true, false);

    while (nn >= 0) {
        long its = 0,
             l   = 0;
        do {
            for (l = nn; l>=1; l--) {
                _Parameter s = fabs (MX_ACCESS(l-1,l-1)) + fabs (MX_ACCESS(l,l));
                if (s == 0.0) {
                    s = anorm;
                }

                if (fabs (MX_ACCESS(l,l-1)) + s == s) {
                    break;
                }
            }

            _Parameter x = MX_ACCESS(nn,nn);
            if (l  == nn) { // one root
                real.theData[nn]   = x + t;
                imag.theData[nn--] = 0.0;
            } else {
                _Parameter y = MX_ACCESS(nn-1,nn-1),
                           w = MX_ACCESS(nn,nn-1)*MX_ACCESS(nn-1,nn);

                if ( l == nn - 1) { // two roots
                    _Parameter p = 0.5 * (y-x),
                               q = p*p + w,
                               z = sqrt (fabs(q));

                    x += t;

                    if (q >= 0.0) { // real pair
                        z = p + (p>0.0?z:-z);
                        real.theData[nn] = real.theData[nn-1] = x+z;
                        if (z) {
                            real.theData[nn] = x-w/z;
                        }
                        imag.theData[nn] = imag.theData[nn-1] = 0.0;
                    } else { // complex pair
                        real.theData[nn]   = real.theData [nn-1] = x+p;
                        imag.theData[nn-1] = -(imag.theData[nn] = z);
                    }
                    nn -= 2;
                } else { // no roots; continue iteration

                    _Parameter p,q,r,z,s;

                    if (its == 30) {
                        _String errMsg ("Too many QR iterations in EigenDecomp");
                        WarnError (errMsg);
                        return;
                    }

                    if (its == 10 || its == 20) {
                        t += x;
                        for (long i=0; i<hDim; i++) {
                            MX_ACCESS(i,i) -= x;
                        }
                        _Parameter s = fabs(MX_ACCESS(nn,nn-1)) + fabs (MX_ACCESS(nn-1,nn-2));
                        y = x = 0.75 * s;
                        w = -0.4375*s*s;
                    }
                    its++;

                    long m = nn-2;

                    for (; m>=l; m--) {
                        z = MX_ACCESS(m,m);
                        r = x-z;
                        s = y-z;
                        p = (r*s - w)/MX_ACCESS(m+1,m) + MX_ACCESS(m,m+1);
                        q = MX_ACCESS(m+1,m+1)-z-r-s;
                        r = MX_ACCESS(m+2,m+1);
                        s = fabs (p) + fabs (q) + fabs (r);

                        p/=s;
                        q/=s;
                        r/=s;

                        if (m == l) {
                            break;
                        }

                        _Parameter u = fabs (MX_ACCESS(m,m-1)) * (fabs (q) + fabs (r)),
                                   v = fabs (p) * (fabs (MX_ACCESS(m-1,m-1)) + fabs(z) + fabs(MX_ACCESS(m+1,m+1)));

                        if (u+v == v) {
                            break;
                        }

                    }

                    for (long i = m+2; i<hDim; i++) {
                        MX_ACCESS(i,i-2) = 0.0;
                        if (i!=m+2) {
                            MX_ACCESS(i,i-3) = 0.0;
                        }
                    }

                    for (long k = m; k <= nn-1; k++) {
                        if (k!=m) {
                            p = MX_ACCESS(k,k-1),
                            q = MX_ACCESS(k+1,k-1),
                            r = 0.0;
                            if (k != nn-1) {
                                r = MX_ACCESS(k+2,k-1);
                            }

                            if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }

                        s = sqrt (p*p+q*q+r*r);

                        if (s != 0.0) {
                            if (p < 0.0) {
                                s = -s;
                            }

                            if (k == m) {
                                if (l!=m) {
                                    MX_ACCESS(k,k-1) = -MX_ACCESS(k,k-1);
                                }
                            } else {
                                MX_ACCESS(k,k-1) = -s*x;
                            }

                            p += s;
                            x = p/s;
                            y = q/s;
                            z=r/s;
                            q/=p;
                            r/=p;
                            for (long j=k; j<=nn; j++) {
                                p = MX_ACCESS(k,j)+q*MX_ACCESS(k+1,j);
                                if (k!=nn-1) {
                                    p += r*MX_ACCESS(k+2,j);
                                    MX_ACCESS(k+2,j) -= p*z;
                                }
                                MX_ACCESS(k+1,j) -= p*y;
                                MX_ACCESS(k,j)   -= p*x;
                            }

                            long mmin = nn < k+3 ? nn: k+3;
                            for (long i = 0; i<=mmin; i++) {
                                p = x*MX_ACCESS(i,k) + y*MX_ACCESS(i,k+1);
                                if (k!=nn-1) {
                                    p += z*MX_ACCESS(i,k+2);
                                    MX_ACCESS(i,k+2) -= p*r;
                                }
                                MX_ACCESS(i,k+1) -= p*q;
                                MX_ACCESS(i,k) -= p;
                            }
                        }
                    }
                }
            }

        } while (l < nn - 1);
    }

}

//__________________________________________________________________________________
_PMathObj   _Matrix::Eigensystem (void)
{
    // find the eigenvectors of a symmetric matrix using Jacobi rotations
    // The original matrix is preserved.
    // returns an associative list with a sorted vector of eigenvalues and
    // a square matrix where columns are the corresponding eigenvalues

    if ((storageType!=1)||(hDim!=vDim)||(hDim==0)) { // only works for numerical matrices at this stage
        _String errorMsg ("Eigensystem only works with numerical non-empty square matrices");
        WarnError (errorMsg);
        return    new _AssociativeList();
    }

    // check for symmetry

    for (long k=0; k<hDim; k++) {
        for (long v=k+1; v<hDim; v++) {
            if (!CheckEqual((*this)(k,v), (*this)(v,k))) {
                //_String errorMsg ("Eigensystem presently only works on symmetric matrices");
                //WarnError (errorMsg);
                //return      new _AssociativeList();

                //_String nonSym = _String("Failed symmetry check:" ) & k & ":" & v;
                //WarnError (nonSym);

                _AssociativeList * res = new _AssociativeList;

                _Matrix            *cpy = new _Matrix (*this),
                *rl  = new _Matrix,
                *im  = new _Matrix;

                _String            key ("0");

                cpy->CheckIfSparseEnough(true);

                checkPointer (cpy);
                cpy->Balance ();
                cpy->Schur   ();
                cpy->EigenDecomp (*rl,*im);


                key = "0";
                {
                    _FString fk (key, false);
                    res->MStore (&fk, rl, false);
                }
                key = "1";
                {
                    _FString fk (key, false);
                    res->MStore (&fk, im, false);
                }
                key = "2";
                DeleteObject (cpy);
                return res;
            }
        }
    }
    _Matrix a (*this);
    a.CheckIfSparseEnough (true);

    _Parameter* b = new _Parameter[hDim],
    *   z = new _Parameter[hDim];

    _Matrix * d = new _Matrix(hDim,1,false, true),
    * v = new _Matrix(hDim,hDim,false,true);

    checkPointer ((Ptr)(b && z && d && v));

    for (long cnt = 0, diagIndex=0; cnt < hDim; cnt ++, diagIndex+=hDim+1) {
        v->theData[diagIndex] = 1.;
        b[cnt] = (d->theData[cnt] = a.theData [diagIndex]);
        z[cnt] = 0.0;
    }

    for (long pass = 0; pass < 50; pass ++) {
        _Parameter sm = 0.,
                   tresh = 0.;

        {
            for (long ec = 0; ec < hDim-1; ec ++)
                for (long ec2 = ec+1; ec2 < hDim; ec2++) {
                    sm += fabs(a.theData[ec*hDim+ec2]);
                }
        }

        if (sm == 0.0) {
            break;
        }

        if (pass < 3) {
            tresh = 0.2 * sm / (hDim*hDim);
        }

        for (long ec=0; ec < hDim-1; ec++) {
            for (long ec2=ec+1; ec2 < hDim; ec2++) {
                long       midx = ec*hDim+ec2;

                _Parameter mel = a.theData[midx],
                           g   = 100. * fabs (mel),
                           t   = fabs (d->theData[ec]),
                           c   = fabs (d->theData[ec2]);

                if (pass>3 && t+g == t && c+g == c) {
                    a.theData[midx] = 0.;
                } else if (fabs(mel) > tresh) {
                    _Parameter h = d->theData[ec2]-d->theData[ec];
                    if (fabs (h) + g == fabs (h)) {
                        t = mel/h;
                    } else {
                        _Parameter theta = 0.5*h/mel;
                        t = 1./(fabs(theta)+sqrt(1.+theta*theta));
                        if (theta<0.0) {
                            t = -t;
                        }
                    }

                    c = 1.0/sqrt(1.0+t*t);

                    _Parameter s    = t*c;
                    _Parameter tau  = s/(1.0+c);

                    h = t*mel;

                    z[ec]           -= h;
                    z[ec2]          += h;
                    d->theData[ec]  -= h;
                    d->theData[ec2] += h;

                    /*_String checkV = (c*c-s*s)*mel + s*c*(a.theData[ec*hDim+ec]-a.theData[ec2*hDim+ec2]);
                    StringToConsole (checkV);
                    BufferToConsole ("\n");*/

                    a.theData[midx] = 0.;

                    for (long j=0; j<ec; j++) {
                        ROTATE (a.theData, j, ec, j, ec2, g, h, s, tau, hDim);
                    }

                    {
                        for (long j=ec+1; j<ec2; j++) {
                            ROTATE (a.theData, ec, j, j, ec2, g, h, s, tau, hDim);
                        }
                    }
                    {
                        for (long j=ec2+1; j<hDim; j++) {
                            ROTATE (a.theData, ec, j, ec2, j, g, h, s, tau, hDim);
                        }
                    }
                    {
                        for (long j=0; j<hDim; j++) {
                            ROTATE (v->theData, j, ec, j, ec2, g, h, s, tau, hDim);
                        }
                    }
                }
            }
        }
        {
            for (long ec=0; ec<hDim; ec++) {
                b[ec] += z[ec];
                d->theData[ec] = b[ec];
                z[ec] = 0.;
            }
        }
    }



    _Matrix * ds = new _Matrix(hDim,2,false, true),
    * vs = new _Matrix(hDim,hDim,false,true),
    * dss;

    for (long r=0; r<hDim; r++) {
        ds->theData[2*r]   = -d->theData[r];
        ds->theData[2*r+1] = r;
    }

    _Constant sc (0.0);
    dss = (_Matrix*)ds->SortMatrixOnColumn (&sc);
    DeleteObject (ds);
    {
        for (long r=0; r<hDim; r++) {
            d->theData[r] = -dss->theData[2*r];
            for (long c1 = r, c2 = dss->theData[2*r+1]; c1<hDim*hDim; c1+=hDim, c2+=hDim) {
                vs->theData[c1] = v->theData[c2];
            }
        }
    }

    _AssociativeList * res = new _AssociativeList ();
    checkPointer (res);

    _String            key ("0");

    {
        _FString fk (key, false);
        res->MStore (&fk, d, false);
    }

    key = "1";

    {
        _FString fk (key, false);
        res->MStore (&fk, vs, false);
    }

    DeleteObject (v);
    //DeleteObject (d);
    DeleteObject (dss);
    //DeleteObject (vs);

    delete [] b;
    delete [] z;


    return res;
}

//__________________________________________________________________________________
_PMathObj   _Matrix::LUDecompose (void)
{
    // perform the LU decomposition using Crout's algorithm with partial pivoting
    // The original matrix is preserved.
    // after performing this decomposition, the routine LUSolve can be called with an arbitrary vector
    // the return object is an nx(n+1) matrix which contains the LU decomposition followed
    // by a vector of row interchanges
    if (storageType!=1 || hDim!=vDim || hDim==0) { // only works for numerical matrices at this stage
        _String errorMsg ("LUDecompose only works with numerical non-empty square matrices");
        WarnError (errorMsg);
        return    new _Matrix();
    }

    _Parameter *        scalings = new _Parameter[hDim];
    checkPointer        (scalings);

    long perRow = vDim+1;
    _Matrix * result = new _Matrix (hDim,perRow,false,true);

    checkPointer (result);
    // duplicate the original matrix into result
    if (theIndex) //matrix is sparse
        for (long i=0; i<lDim; i++) {
            if (IsNonEmpty(i)) {
                long maxRow = theIndex[i];
                result->Store(maxRow/vDim,maxRow%vDim,theData[i]);
            }
        }
    else
        for (long i=0; i<hDim; i++) {
            long maxRow = i*vDim;
            for (long j=0; j<vDim; j++) {
                result->theData[maxRow+i+j]=theData[maxRow+j];
            }
        }

    // produce the scaling vector used in interchanging the rows
    for (long i=0; i<vDim; i++) {
        _Parameter rowMax = 0.0,
                   cell;

        for (long j=i*perRow; j<(i+1)*perRow-1; j++)
            if ((cell=fabs(result->theData[j]))>rowMax) {
                rowMax = cell;
            }

        if (rowMax==0.0) {
            _String errorMsg = _String("LUDecompose doesn't work on singular matrices (row ") & i & ')';
            WarnError (errorMsg);
            return    nil;
        }
        scalings[i]=1.0/rowMax;
    }
    // main loop for finding L and U

    for (long j=0; j<vDim; j++) {
        {
            for (long i=0; i<j; i++) {
                // fill in superdiagonal elements (U) in column j
                _Parameter sum = result->theData[i*perRow+j];
                for (long k=0; k<i; k++) {
                    sum -= result->theData[i*perRow+k]*result->theData[k*perRow+j];
                }
                result->theData[i*perRow+j] = sum;
            }
        }
        long       maxRow = 0;
        _Parameter rowMax = 0.0,
                   cell;

        for (long i=j; i<hDim; i++) {
            // calculate the unscaled version of elements of L and the diagonal
            _Parameter sum = result->theData[i*perRow+j];

            for (long k=0; k<j; k++) {
                sum -= result->theData[i*perRow+k]*result->theData[k*perRow+j];
            }

            result->theData[i*perRow+j]=sum;

            if ((cell=scalings[i]*fabs(sum))>= rowMax) { // find max under the diagonal in column j
                rowMax = cell;
                maxRow = i;
            }
        }

        if (j!=maxRow) { // interchange rows
            for (long k=0; k<hDim; k++) {
                cell                                    =   result->theData[maxRow*perRow+k];
                result->theData[maxRow*perRow+k]        =   result->theData[j*perRow+k];
                result->theData[j*perRow+k]             =   cell;
            }
            scalings[maxRow]=scalings[j];
        }
        // store the index permutation
        result->theData[j*perRow+vDim] = maxRow;

        if (result->theData[j*perRow+j]==0.0) {
            result->theData[j*perRow+j] = 1.0e-25;
        }

        // divide by the pivoting element

        if (j!=hDim-1) {
            cell = 1.0/result->theData[j*perRow+j];
            for (long i=j+1; i<hDim; i++) {
                result->theData[i*perRow+j]*=cell;
            }
        }
    }
    delete [] scalings;
    return result;
}
//__________________________________________________________________________________
_PMathObj   _Matrix::LUSolve (_PMathObj p)
// takes a matrix in LU decomposed state and a vector of row permutation returned by LU
// returns a vector of solutions
{

    if ((storageType!=1)||(hDim+1!=vDim)||(vDim<=0)) { // only works for numerical matrices at this stage
        _String errorMsg ("LUSolve only works with numerical non-empty matrices of dimension nx(n+1) returned by LUDecompose.");
        WarnError (errorMsg);
        return    nil;
    }
    if (p->ObjectClass()==MATRIX) {
        _Matrix *b=(_Matrix*)p;
        if (!((b->hDim!=hDim)||(b->vDim!=1)||(b->storageType!=1))) {
            _Parameter sum;
            _Matrix result (*b);
            result.CheckIfSparseEnough(true);
            long i,j,trueI, firstI = -1;
            for (i=0; i<hDim; i++) {
                trueI = (*this)(i,vDim-1);
                if ((trueI<0)||(trueI>=hDim)) {
                    break;
                }
                sum = result.theData[trueI];
                result.theData[trueI]=result.theData[i];
                if (firstI>=0)
                    for (j=firstI; j<i; j++) {
                        sum -= theData[i*vDim+j]*result.theData[j];
                    }
                else if (sum) {
                    firstI = i;
                }
                result.theData[i]=sum;
            }
            if (i==hDim) {
                for (i=hDim-1; i>-1; i--) {
                    sum = result.theData[i];
                    for (j=i+1; j<hDim; j++) {
                        sum-=theData[i*vDim+j]*result.theData[j];
                    }
                    result.theData[i]=sum/theData[i*vDim+i];
                }
                return (_PMathObj)result.makeDynamic();
            }
        }
    }
    _String errorMsg ("LUSolve expects the 2nd parameter to be a column vector defining the right hand side of LUx=b");
    WarnError (errorMsg);
    return new _Matrix(1,1,false,true);

}



//__________________________________________________________________________________
_PMathObj   _Matrix::CholeskyDecompose (void)
{
    /* ---------------------------------------------------
        CholeskyDecompose()
            Constrcts lower triangular matrix L such that
            its own transpose can serve as upper part in
            LU decomposition.
            Requires that matrix is symmetric and positive
            definite.
        * algorithm based on Numerical Recipes
       --------------------------------------------------- */

    if (storageType!=1 || hDim!=vDim || hDim==0) { // only works for numerical square matrices at this stage
        _String errorMsg ("CholeskyDecompose only works with numerical non-empty square matrices");
        WarnError (errorMsg);
        return    new _Matrix();
    }

    long        n           = GetHDim();
    _Parameter  sum;
    _Matrix *   lowerTri    = new _Matrix ((_Matrix &)*this);   // duplication constructor

    checkPointer (lowerTri);

    for (long i = 0; i < n; i++) {
        for (long j = i; j < n; j++) {
            sum = (*lowerTri)(i,j);

            for (long k = i-1; k >= 0; k--) {
                sum -= (*lowerTri)(i,k) * (*lowerTri)(j,k);
            }

            if (i==j) {
                if (sum <= 0.0) {   // matrix is not positive-definite
                    WarnError(_String("In CholeskyDecompose(): matrix not positive definite, (row ") & i & ')');
                    return nil;
                }

                lowerTri->Store (i, i, sqrt(sum));
            }

            else {
                lowerTri->Store (j, i, sum / ((*lowerTri)(i,i)));
            }
        }
    }

    /* zero upper triagonal entries */
    for (long i = 0; i < n; i++) {
        for (long j = i+1; j < n; j++) {
            lowerTri->Store (i, j, 0.);
        }
    }

    ReportWarning (_String("_Matrix::CholeskyDecompose returning with ") & (_String *) lowerTri->toStr());

    return lowerTri;
}



//__________________________________________________________________________________
_PMathObj   _Matrix::Log (void)
{
    if (storageType==1) {
        _Matrix* res = new _Matrix;
        checkPointer (res);
        res->Duplicate (this);
        if (theIndex) {
            for (long k=0; k<lDim; k++)
                if (theIndex[k]>=0) {
                    res->theData[k] = log(theData[k]);
                }
        } else {
            for (long k=0; k<lDim; k++) {
                res->theData[k] = log(theData[k]);
            }
        }
        return res;
    }
    _String errorMsg ("Can't apply logs to non-numeric matrices.");
    WarnError (errorMsg);
    return new _Matrix(1,1,false,true);
}

//__________________________________________________________________________________
_PMathObj   _Matrix::Inverse (void)
{
    if ((storageType!=1)||(hDim!=vDim)||(hDim==0)) {
        _String errorMsg ("Inverse only works with numerical non-empty square matrices.");
        WarnError (errorMsg);
        return    nil;
    }
    _Matrix * LUdec = (_Matrix*)LUDecompose();
    if (LUdec) {
        _Matrix b      (hDim,1,false,true),
                result (hDim,vDim,false,true);
        b.theData[0]=1.0;
        for (long i=0; i<hDim; i++) {
            if (i) {
                b.theData[i]=1.0;
                b.theData[i-1]=0.0;
            }
            _Matrix* invVector = (_Matrix*)LUdec->LUSolve(&b);
            _Matrix* corrTerm = (_Matrix*)(*this*(*invVector)-b).makeDynamic();
            _Matrix* corrX =  (_Matrix*)LUdec->LUSolve(corrTerm);
            *invVector-=*corrX;
            DeleteObject (corrX);
            DeleteObject (corrTerm);
            for (long j=0; j<hDim; j++) {
                result.theData[j*vDim+i]=invVector->theData[j];
            }
            DeleteObject (invVector);
        }
        DeleteObject (LUdec);
        return (_PMathObj)result.makeDynamic();
    }
    return new _Matrix (1,1,false,true);

}

//__________________________________________________________________________________
_PMathObj   _Matrix::MultByFreqs (long freqID)
// multiply this transition probs matrix by frequencies
{
    _PMathObj value = ComputeNumeric(true);

    if (freqID>=0) {
        _Matrix* freqMatrix = nil;
        freqID = modelFrequenciesIndices.lData[freqID];
        if (freqID>=0) {
            freqMatrix = (_Matrix*)LocateVar(freqID)->GetValue();
            if (freqMatrix->storageType != 1) {
                if (freqMatrix->theValue) {
                    freqMatrix = (_Matrix*)freqMatrix->theValue;
                } else {
                    freqMatrix = (_Matrix*)freqMatrix->ComputeNumeric();
                }
            }
        }

        if (theIndex) {
            _Matrix*    vm = (_Matrix*) value;
            _Parameter *dp = vm ->theData;

            _Parameter *tempDiags = new _Parameter [hDim];

            for (long i=0; i<hDim; i++) {
                tempDiags[i] = 0.0;
            }

            if (freqMatrix)
                for (long i=0; i<lDim; i++) {
                    long p;
                    if ((p = theIndex[i])!=-1) {
                        long h = p/vDim;
                        p %= vDim;
                        if (h!=p) {
                            tempDiags[h] += (dp[i] *= freqMatrix->theData[p]);
                        }
                    }
                }
            else
                for (long i=0; i<lDim; i++) {
                    long p;
                    if ((p = theIndex[i])!=-1) {
                        long h = p/vDim;
                        p %= vDim;
                        if (h!=p) {
                            tempDiags[h] += dp[i];
                        }
                    }
                }

            for (long j=0; j<hDim; j++) {
                vm->Store (j,j,-tempDiags[j]);
            }

            delete [] tempDiags;
        } else {
            _Parameter * theMatrix = ((_Matrix*)value)->theData;

            if (freqMatrix) {
                if (freqMatrix->theIndex) {
                    for (long i=0; i<lDim; i++) {
                        theMatrix[i] *= (*freqMatrix)[i%vDim];
                    }
                } else {
                    for (long i=0; i<lDim; i++) {
                        theMatrix[i] *= freqMatrix->theData[i%vDim];
                    }
                }
            }
            {
                for (long i=0; i<lDim; i+=(vDim+1)) {
                    theMatrix[i] = 0.;
                }
            }
            for (long i=0; i<lDim; i++) {
                long h = i/vDim,v = i%vDim;
                if (h!=v) {
                    theMatrix[h*vDim+h]-=theMatrix[h*vDim+v];
                }
            }
        }

    }
    return value;
}


//__________________________________________________________________________________
_PMathObj   _Matrix::Compute (void)
{
    //if ((storageType != 1)&&(storageType != 2))
    if (storageType != 1) {
        if (storageType == 0) {
            if (ANALYTIC_COMPUTATION_FLAG) {
                return this;
            }
        }
        if (IsAStringMatrix()) {
            return this;
        }

        if (theValue) {
            DeleteObject (theValue);
        }
        if (storageType != 3) {
            theValue = Evaluate(false);
        } else {
            theValue  = EvaluateSimple ();
        }
        return theValue;
    }
    return this;
}

//__________________________________________________________________________________
_PMathObj   _Matrix::ComputeNumeric (bool copy)
{
    if (storageType != 1) {
        if (storageType == 0 && ANALYTIC_COMPUTATION_FLAG) {
            return this;
        }

        if (theValue) {
            DeleteObject (theValue);
        }

        if (storageType != 3) {
            theValue  = Evaluate(false);
        } else {
            theValue  = EvaluateSimple ();
        }
        return theValue;
    }
    if (copy) {
        if (theValue) {
            DeleteObject (theValue);
        }

        theValue = (_Matrix*)makeDynamic();
        return theValue;
    }
    return this;
}

//__________________________________________________________________________________
_PMathObj   _Matrix::RetrieveNumeric (void)
{
    if (storageType != 1) {
        if (theValue) {
            return theValue;
        }

        return ComputeNumeric();
    }
    return this;
}

//__________________________________________________________________________________
_PMathObj   _Matrix::Sum (void)
{
    return new _Constant (MaxElement (1));
}

//__________________________________________________________________________________


_PMathObj _Matrix::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{
    //_Constant res;
    // why was static?? mod 07/21/2003


    switch (opCode) {
    case HY_OP_CODE_IDIV: // $
        return MultElements(p);
        break;
    case HY_OP_CODE_MOD: // %
        return SortMatrixOnColumn (p);
        break;
    case HY_OP_CODE_AND: // &&
        return pFDR (p);
        break;
    case HY_OP_CODE_MUL: // *
        return MultObj(p);
        break;
    case HY_OP_CODE_ADD: // +
        if (p) {
            return AddObj (p);
        } else {
            return Sum ();
        }
        break;
    case HY_OP_CODE_SUB: // -
        if (p) {
            return SubObj(p);
        } else {
            return (_PMathObj)((*this)*(-1.0)).makeDynamic();
        }
        break;
    case HY_OP_CODE_LESS: // <
        return PathLogLikelihood(p);
        break;
    case HY_OP_CODE_LEQ: // <=
        return K_Means(p);
        break;
    case HY_OP_CODE_EQ: // ==
        return ProfileMeanFit(p);
        break;
    case HY_OP_CODE_GREATER: // >
        return NeighborJoin (!CheckEqual(p->Value(),0.0));
        break;
    case HY_OP_CODE_GEQ: // >=
        return MakeTreeFromParent (p->Value());
        break;
    case HY_OP_CODE_ABS: // Abs
        return Abs();
        break;
    case HY_OP_CODE_CCHI2: //CChi2
        if (p->ObjectClass()==NUMBER && p->Value()>0.999 ) {
            return new _Constant (FisherExact(5.,80.,1.));
        } else {
            return new _Constant (FisherExact(0.,0.,0.));
        }
        break;
    case HY_OP_CODE_COLUMNS:  //Columns
        return new _Constant (vDim);
        break;
    case HY_OP_CODE_EIGENSYSTEM: //Eigensystem
        return Eigensystem();
        break;
    case HY_OP_CODE_EXP: //Exp
        return Exponentiate();
        break;
    case HY_OP_CODE_INVERSE: //Inverse
        return Inverse();
        break;
    case HY_OP_CODE_LUDECOMPOSE: // LUDecompose
        return LUDecompose();
        break;
    case HY_OP_CODE_LUSOLVE: // LUSolve
        return LUSolve (p);
        break;
    case HY_OP_CODE_LOG: // Log
        return Log();
        break;
    case HY_OP_CODE_MACCESS: // MAccess
        return MAccess (p,p2);
        break;
    case HY_OP_CODE_MAX: // Max
    case HY_OP_CODE_MIN: // Max
        if (p->ObjectClass()==NUMBER) {
            if (CheckEqual (p->Value(), 1)) {
                long index;
                _Parameter v[2] = {opCode == HY_OP_CODE_MAX?MaxElement (0,&index):MinElement(0,&index),0.0};
                v[1] = index;
                return new _Matrix (v,1,2);
            }
        }
        return new _Constant (opCode == HY_OP_CODE_MAX?MaxElement (0):MinElement (0));
        break;

    case HY_OP_CODE_MCOORD: // MCoord
        return MCoord (p,p2);
        break;
    case HY_OP_CODE_RANDOM: // Random
        return Random (p);
        break;
    case HY_OP_CODE_ROWS: // Rows
        return new _Constant (hDim);
        break;
    case HY_OP_CODE_SIMPLEX: // Simplex
        return SimplexSolve();
        break;
    case HY_OP_CODE_TRANSPOSE: { // Transpose
        _Matrix* result = (_Matrix*)makeDynamic();
        result->Transpose ();
        return result;
    }
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    case HY_OP_CODE_POWER: // ^ (Poisson log-likelihood)
        return  PoissonLL (p);
    }

    WarnNotDefined (this, opCode);
    return nil;
}
//_____________________________________________________________________________________________

_Matrix::_Matrix ()                             // default constructor, doesn't do much
{
    Initialize();
}
//_____________________________________________________________________________________________

void _Matrix::Initialize ()                             // default constructor, doesn't do much
{
    theData         = nil;
    theIndex        = nil;
    vDim = hDim = lDim  = bufferPerRow = overflowBuffer = 0;
    storageType     = 1;
    allocationBlock = 1;
    theValue        = nil;

}

//_____________________________________________________________________________________________
void    CreateMatrix    (_Matrix* theMatrix, long theHDim, long theVDim,  bool sparse = false, bool allocateStorage = false, bool isFla = false)
{
    long i;

    theMatrix->theValue     = nil;
    theMatrix->storageType  = allocateStorage;
    if (theHDim && theVDim) {
        if (sparse) { // store matrix as sparse
            theMatrix->lDim = theHDim*theVDim*theMatrix->storageIncrement/100+1; // size of storage in elements
            if (theMatrix->lDim-1<theHDim)
                // either the matrix or the allocation block are too small
                // to sensibly store the matrix as sparse.
            {
                CreateMatrix (theMatrix, theHDim, theVDim, false, allocateStorage, isFla);
                return;
            }
            if (!(theMatrix->theIndex = (long*)MatrixMemAllocate(sizeof(long)*theMatrix->lDim))) { // allocate element index storage
                warnError(-108);
                return;
            } else {
                for (i = 0; i<theMatrix->lDim; i++) {
                    theMatrix->theIndex[i] = -1;
                }
            }

        } else {
            theMatrix->lDim = theHDim*theVDim;
            theMatrix->theIndex = nil; // no index storage needed
        }

        if (!allocateStorage)
            // matrix will store pointers to elements
        {
            if (!(theMatrix->theData =(_Parameter*)MatrixMemAllocate(theMatrix->lDim*sizeof(void*)))) { // allocate element index storage
                warnError(-108);
                return;
            } else {
                // populate with zero data pointers
                if (!isFla)
                    for (i = 0; i<theMatrix->lDim; i++) {
                        ((_MathObject**)theMatrix->theData)[i] = ZEROPOINTER;
                    }
                else
                    for (i = 0; i<theMatrix->lDim; i++) {
                        ((_Formula**)theMatrix->theData)[i]    = ZEROPOINTER;
                    }
            }

        } else {
            if (!(theMatrix->theData =(_Parameter*)MatrixMemAllocate (sizeof(_Parameter)*theMatrix->lDim))) { // allocate element index storage
                warnError(-108);
                return;
            }

            else {
                // populate with zero data objects
                for (i = 0; i<theMatrix->lDim; i++) {
                    ((_Parameter*)theMatrix->theData)[i] = ZEROOBJECT;
                }
            }
        }
    } else {
        theMatrix->lDim      = 0;
        theMatrix->theIndex  = nil;
        theMatrix->theData   = nil;
    }

    theMatrix->hDim = theHDim;
    theMatrix->vDim = theVDim;
    theMatrix->bufferPerRow = theMatrix->overflowBuffer =  theMatrix->allocationBlock = 0;

}
//_____________________________________________________________________________________________
bool    _Matrix::AmISparse(void)
{
    if (theIndex) {
        return true;    // duh!
    }
    if (storageType==2) {
        return false;
    }

    long k=0;
    if (storageType==1) {
        for (long i=0; i<lDim; i++)
            if (theData[i]!=ZEROOBJECT) {
                k++;
            }
    } else {
        for (long i=0; i<lDim; i++)
            if (IsNonEmpty(i) && !GetMatrixObject(i)->IsObjectEmpty()) {
                k++;
            }
    }


    if ((_Parameter(k)/lDim*100.)<=_Matrix::switchThreshold) {
        // we indeed are sparse enough
        _Matrix sparseMe (hDim,vDim,true,storageType==1);
        if (storageType==1) {
            for (long i=0; i<lDim; i++) {
                if (theData[i]!=ZEROOBJECT) {
                    sparseMe[i]=theData[i];
                }
            }
        } else if (storageType==0) {
            for (long i=0; i<lDim; i++) {
                if ((GetMatrixObject(i)!=ZEROPOINTER)&&(!GetMatrixObject(i)->IsObjectEmpty())) {
                    sparseMe.StoreObject(i,GetMatrixObject(i));
                }
                GetMatrixObject(i)->nInstances++;
            }
        }

        Clear();
        DuplicateMatrix (this, &sparseMe);
        return true;
    }
    return false;
}

//_____________________________________________________________________________________________
bool    _Matrix::AmISparseFast (_Matrix& whereTo)
{
    if (theIndex) {
        return true;    // duh!
    }

    long k=0,
         i;

    for (i=0; i<lDim; i++)
        if (theData[i]!=ZEROOBJECT) {
            k++;
        }

    if ((k*100)/lDim<=_Matrix::switchThreshold) {
        // we indeed are sparse enough
        if (k == 0) {
            k = 1;
        }

        _Parameter *          newData  = (_Parameter*)MatrixMemAllocate (k*sizeof(_Parameter));
        if (whereTo.theIndex) {
            free (whereTo.theIndex);
        }
        whereTo.theIndex               = (long*)MemAllocate (k*sizeof(long));

        if (!(newData&&whereTo.theIndex)) {
            warnError (-108);
        }

        long p = 0;

        whereTo.theIndex[0] = -1;

        for (i=0; i<lDim; i++)
            if (theData[i]!=ZEROOBJECT) {
                whereTo.theIndex[p] = i;
                newData[p++]=theData[i];
            }

        whereTo.lDim     = k;
        free     (whereTo.theData);
        whereTo.theData = newData;
        return true;
    }

    return false;
}

//_____________________________________________________________________________________________

bool    _Matrix::IsReversible(_Matrix* freqs)
{
    if (hDim != vDim || (freqs && freqs->hDim * freqs->vDim != hDim)
            || (storageType != 1 && storageType != 2) ||
            (freqs && freqs->storageType != 1 && freqs->storageType != 2)) {
        return false;
    }

    bool   needAnalytics = storageType == 2 || (freqs && freqs->storageType == 2);
    if (needAnalytics) {
        if (freqs) {
            for (long r = 0; r < hDim; r++)
                for (long c = r+1; c < hDim; c++) {
                    bool compResult = true;
                    if (storageType == 2) {
                        _Formula* rc = GetFormula(r,c),
                                  * cr = GetFormula(c,r);

                        if (rc && cr) {
                            _Polynomial *rcp = (_Polynomial *)rc->ConstructPolynomial(),
                                         *crp = (_Polynomial *)cr->ConstructPolynomial();

                            if (rcp && crp) {
                                _PMathObj     tr = nil,
                                              tc = nil;

                                if (freqs->storageType == 2) {
                                    if (freqs->GetFormula(r,0)) {
                                        tr = freqs->GetFormula(r,0)->ConstructPolynomial();
                                        if (tr) {
                                            tr->nInstances++;
                                        } else {
                                            return false;
                                        }
                                    }
                                    if (freqs->GetFormula(c,0)) {
                                        tc = freqs->GetFormula(c,0)->ConstructPolynomial();
                                        if (tc) {
                                            tc->nInstances++;
                                        } else {
                                            DeleteObject (tr);
                                            return false;
                                        }
                                    }

                                } else {
                                    tr = new _Constant ((*freqs)[r]);
                                    tc = new _Constant ((*freqs)[c]);
                                }
                                if (tr && tc) {
                                    _Polynomial        * rcpF = (_Polynomial*)rcp->Mult(tr),
                                                         * crpF = (_Polynomial*)crp->Mult(tc);

                                    compResult         = rcpF->Equal(crpF);
                                    DeleteObject (rcpF);
                                    DeleteObject (crpF);
                                } else {
                                    compResult = !(tr||tc);
                                }

                                DeleteObject (tr);
                                DeleteObject (tc);
                            } else {
                                compResult = false;
                            }

                            //DeleteObject (rcp); DeleteObject (crp);
                        } else {
                            compResult = !(rc || cr);
                        }
                    }
                    if (!compResult) {
                        return false;
                    }
                }
        } else {
            for (long r = 0; r < hDim; r++)
                for (long c = r+1; c < hDim; c++) {
                    bool compResult = true;
                    _Formula* rc = GetFormula(r,c),
                              * cr = GetFormula(c,r);

                    if (rc && cr) {
                        _Polynomial *rcp = (_Polynomial *)rc->ConstructPolynomial(),
                                     *crp = (_Polynomial *)cr->ConstructPolynomial();

                        if (rcp && crp) {
                            compResult = rcp->Equal(crp);
                        } else {
                            compResult = rc->EqualFormula(cr);
                        }

                        //DeleteObject (rcp); DeleteObject (crp);
                    } else {
                        compResult = !(rc || cr);
                    }

                    if (!compResult) {
                        return false;
                    }
                }
        }
        return true;
    } else {
        if (freqs) {
            for (long r = 0; r < hDim; r++)
                for (long c = r+1; c < hDim; c++)
                    if (! CheckEqual ((*this)(r,c)*(*freqs)[r], (*this)(c,r)*(*freqs)[c])) {
                        return false;
                    }
        } else {
            for (long r = 0; r < hDim; r++)
                for (long c = r+1; c < hDim; c++)
                    if (! CheckEqual ((*this)(r,c), (*this)(c,r))) {
                        return false;
                    }
        }
        return true;
    }
    return false;
}

//_____________________________________________________________________________________________

void    _Matrix::CheckIfSparseEnough(bool force)

// check if matrix is sparse enough to justify compressed storage
{
    long i;

    if (theIndex && (force || lDim>hDim*vDim*switchThreshold/100))
        // switch to normal matrix storage - more than half elements are non-zero
    {
        // -= allocationBlock;
        if (storageType!=1)
            // pointers
        {
            long* tempData;

            if (!(tempData = (long*) MatrixMemAllocate(hDim*vDim*sizeof(Ptr)))) {
                warnError(-108);
            } else {
                for (i = 0; i<hDim*vDim; i++) {
                    tempData[i]=0;
                }
                for (i = 0; i<lDim; i++) {
                    if (IsNonEmpty(i)) {
                        tempData[theIndex[i]]=((long*)theData)[i];
                    }
                }
                MatrixMemFree (theData);
                theData = (_Parameter*)tempData;
            }
        } else
            //objects
        {
            _Parameter* tempData;
            if (!(tempData =  (_Parameter*)MatrixMemAllocate (sizeof(_Parameter)*hDim*vDim))) {
                warnError(-108);
            } else {
                for (i = 0; i<hDim*vDim; i++) {
                    tempData [i] = ZEROOBJECT;
                }
                for (i = 0; i<lDim; i++) {
                    long k = theIndex[i];
                    if (k!=-1) {
                        tempData [k] = ((_Parameter*)theData) [i];
                    }
                }
                MatrixMemFree( theData);
                theData = tempData;
            }

        }
        MatrixMemFree (theIndex);
        theIndex = nil;
        bufferPerRow = overflowBuffer = allocationBlock = 0;
        lDim = vDim*hDim;
    }
}

//_____________________________________________________________________________________________
bool    _Matrix::IncreaseStorage    (void)
{
    lDim += allocationBlock;

    long* tempIndex, i;

    if (!(tempIndex = (long*)MatrixMemAllocate(lDim*sizeof(long)))) {
        warnError(-108);
    } else {
        memcpy (tempIndex, theIndex, (lDim-allocationBlock)*sizeof(long));
        MatrixMemFree( theIndex);

        for (i = lDim-1; i>=lDim-allocationBlock; i--) {
            tempIndex [i] = -1;
        }
        theIndex = tempIndex;
    }

    if (storageType != 1)
        // pointers or formulas
    {
        _MathObject** tempData;
        if (!(tempData = (_MathObject**) MatrixMemAllocate(sizeof( char)* lDim*sizeof(void*)))) {
            warnError(-108);
        } else {
            memcpy (tempData, theData, (lDim-allocationBlock)*sizeof(void*));
            MatrixMemFree (theData);
            for (i = lDim-1; i>=lDim-allocationBlock; i--) {
                tempData [i] = ZEROPOINTER;
            }
            theData = (_Parameter*)tempData;
        }
    } else
        //objects
    {
        _Parameter* tempData;
        if (!(tempData =  (_Parameter*)MatrixMemAllocate(sizeof(_Parameter)* lDim))) {
            warnError(-108);
        } else {
            for (i = lDim-1; i>=lDim-allocationBlock; i--) {
                tempData [i] = ZEROOBJECT;
            }
            for (; i>=0; i--) {
                tempData [i] = ((_Parameter*)theData) [i];
            }
            MatrixMemFree( theData);
            theData = tempData;
        }
    }
    return TRUE;

}

//_____________________________________________________________________________________________
void    DuplicateMatrix (_Matrix* targetMatrix, _Matrix* sourceMatrix)
{
    if (targetMatrix==sourceMatrix) {
        return;
    }
    targetMatrix->lDim = sourceMatrix->lDim;
    targetMatrix->hDim = sourceMatrix->hDim;
    targetMatrix->vDim = sourceMatrix->vDim;
    targetMatrix->storageType = sourceMatrix->storageType;
    targetMatrix->bufferPerRow =sourceMatrix->bufferPerRow;
    targetMatrix->overflowBuffer = sourceMatrix->overflowBuffer;
    targetMatrix->allocationBlock = sourceMatrix->allocationBlock;
    targetMatrix->theValue = nil;

    if (sourceMatrix->theIndex) {
        if (!(targetMatrix->theIndex = (long*)MatrixMemAllocate(sizeof(long) *sourceMatrix->lDim))) { // allocate element index storage
            warnError(-108);
        } else {
            memcpy ((void*)targetMatrix->theIndex,(void*)sourceMatrix->theIndex,sourceMatrix->lDim*sizeof(long));
        }
    } else {
        targetMatrix->theIndex = nil;
    }


    targetMatrix->theData = nil;

    if (sourceMatrix->lDim) {
        if (sourceMatrix->storageType==0)
            // matrix will store pointers to elements
        {
            if (targetMatrix->lDim) {
                if (!(targetMatrix->theData = (_Parameter*)MatrixMemAllocate(sizeof( char)*sourceMatrix->lDim*sizeof(void*)))) { // allocate element index storage
                    warnError(-108);
                } else {
                    memcpy ((void*)targetMatrix->theData,(void*)sourceMatrix->theData,sourceMatrix->lDim*sizeof(void*));
                    if (!sourceMatrix->theIndex) { // non-sparse matrix
                        for (long i=0; i<sourceMatrix->lDim; i++)
                            if (sourceMatrix->GetMatrixObject(i)) {
                                (sourceMatrix->GetMatrixObject(i))->nInstances++;
                            }
                    } else
                        for (long i=0; i<sourceMatrix->lDim; i++) {
                            _MathObject* theO = (sourceMatrix->GetMatrixObject(i));
                            if (theO!=ZEROPOINTER) {
                                theO->nInstances++;
                            }
                        }

                }
            }
        } else if (sourceMatrix->storageType==2) {
            if (targetMatrix->lDim) {
                targetMatrix->theData = (_Parameter*)MatrixMemAllocate(sourceMatrix->lDim*sizeof(void*));
                _Formula ** theFormulas = (_Formula**)(sourceMatrix->theData), **newFormulas =
                                              (_Formula**)(targetMatrix->theData);
                if (sourceMatrix->theIndex) {
                    for (long i = 0; i<sourceMatrix->lDim; i++)
                        if (sourceMatrix->IsNonEmpty(i)) {
                            newFormulas[i] = (_Formula*)theFormulas[i]->makeDynamic();
                        }
                } else
                    for (long i = 0; i<sourceMatrix->lDim; i++)
                        if(theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                            newFormulas[i] = (_Formula*)theFormulas[i]->makeDynamic();
                        } else {
                            newFormulas[i]=ZEROPOINTER;
                        }
            }
        } else {
            if (targetMatrix->lDim) {
                if (!(targetMatrix->theData =(_Parameter*)MatrixMemAllocate(sizeof( _Parameter)*targetMatrix->lDim))) { // allocate element index storage
                    warnError(-108);
                } else {
                    memcpy ((Ptr)targetMatrix->theData,(Ptr)sourceMatrix->theData,sizeof(_Parameter)*sourceMatrix->lDim);
                }
            }
        }
    } else {
        targetMatrix->theData = nil;
        targetMatrix->lDim    = 0;
    }

}
//_____________________________________________________________________________________________
BaseRef _Matrix::makeDynamic (void)
{
    _Matrix * result = new _Matrix;
    checkPointer    (result);
    DuplicateMatrix (result, this);

    return result;
}

//_____________________________________________________________________________________________
void _Matrix::Duplicate (BaseRef obj)
{

    _Matrix* m = (_Matrix*)obj;
    Clear();
    DuplicateMatrix (this,m);
}


//_____________________________________________________________________________________________

_Matrix::_Matrix (long theHDim, long theVDim, bool sparse, bool allocateStorage)    // create an empty matrix of given dimensions;
// the flag specifies whether it is sparse or not

{
    CreateMatrix (this, theHDim, theVDim, sparse, allocateStorage);
}


//_____________________________________________________________________________________________

void    _Matrix::Convert2Formulas (void)
{
    if (storageType == 1) {
        storageType = 2;
        _Formula** tempData = (_Formula**)MatrixMemAllocate (sizeof(void*)*lDim);
        if (!theIndex) {
            for (long i = 0; i<lDim; i++) {
                tempData[i] = new _Formula (new _Constant (((_Parameter*)theData)[i]));
            }
        } else
            for (long i = 0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    //_Constant c (((_Parameter*)theData)[i]);
                    //_Formula f((_PMathObj)c.makeDynamic());
                    //tempData[i] = (_Formula*)f.makeDynamic();
                    tempData[i] = new _Formula (new _Constant (((_Parameter*)theData)[i]));
                } else {
                    tempData[i]=nil;
                }
            }

        MatrixMemFree (theData);
        theData = (_Parameter*)tempData;
    }
}



//_____________________________________________________________________________________________

_Matrix::_Matrix (_String& s, bool isNumeric, _VariableContainer* theP)
// takes two separate formats
// 1st : {{i11,...,i1n}{i21,...,i2n}....{in1,...,inn}} // all elements must be explicitly specified
// 2st : {hor dim, <vert dim>,{hor index, vert index, value or formula}{...}...}
{
    // decide which input type is being presented
    Initialize();

    bool    isAConstant = true; // is this a matrix of numbers, or formulas
    char    cc;


    long    i=s.FirstNonSpaceIndex(),
            j=s.FirstNonSpaceIndex(i+1),
            k=0,
            hPos = 0,
            vPos = 0;

    _String terminators (",}");

    if (j>i && s.sLength>4) { // non-empty string
        _String term;
        if (s.sData[i]=='{' && s.sData[j]=='{') { // first type
            i = j+1;
            // read the dimensions first

            while (i<s.sLength) {
                i = s.FindTerminator (i, terminators);
                if (i < 0) {
                    WarnError ("Unterminated matrix definition");
                    return;
                }
                cc = s.sData[i];

                if (cc=='}') {
                    break;
                }

                if (cc==',') {
                    vDim++;
                }
                i++;
            }

            vDim++;
            hDim = 1;

            for (i = i + 1; i<s.sLength-1; i++) {
                i = s.ExtractEnclosedExpression (i,'{','}',true, true);
                if (i < 0) {
                    break;
                }
                hDim ++;
            }

            if ( hDim<=0 || vDim<=0) {
                return;
            }

            if (isNumeric) {
                CreateMatrix (this, hDim, vDim, false, true, false);
            } else {
                CreateMatrix (this, hDim, vDim, false, false, true);
            }

            // scan the elements one-by-one

            for (i=1; i<s.sLength-1; i++) {
                if (s.sData[i]=='{') {
                    while (s.sData[i]!='}') {
                        i++;
                        j = s.FindTerminator (i, terminators);

                        if (j<0) {
                            WarnError ("Unterminated matrix definition");
                            return;
                        }

                        _String lterm (s,s.FirstNonSpaceIndex(i,j-1,1),j-1); // store the term in a string

                        //printf ("%s\n", lterm.sData);

                        if (isNumeric) {
                            if (lterm.sLength == 1 && lterm.sData[0]=='*') {
                                lterm = empty;    // dummy element in probability matrix
                            }

                            theData[vDim*hPos+vPos] = lterm.toNum();
                        } else {
                            if (lterm.sLength == 1 && lterm.sData[0]=='*') {
                                lterm = empty;    // dummy element in probability matrix
                            }

                            _Formula*  theTerm = (_Formula*)checkPointer(new _Formula (lterm, theP));

                            if (isAConstant) // there is hope that this matrix is of numbers
                                if (theTerm->ObjectClass() == NUMBER) {
                                    isAConstant = theTerm->IsAConstant();
                                } else {
                                    isAConstant = false;
                                }

                            ((_Formula**)theData)[vDim*hPos+vPos] = theTerm;
                        }

                        vPos++;
                        if (vPos>vDim) {
                            WarnError ("Rows of unequal lengths in matrix definition");
                            return;
                        }

                        i=j;
                    }
                }
                if (s[i]=='}') {
                    if (vPos!=vDim) {
                        warnError(-104);
                        return;
                    }
                    hPos++;
                    vPos = 0;
                    if (hPos>hDim) {
                        warnError(-104);
                        return;
                    }
                }
            }
            if (hPos!=hDim) {
                warnError(-104);
                return;
            }
        } else { // second type of input
            for (i=j,j=0; s.sData[i]!='{' && s.sData[i]!='}' && i<s.sLength; i++) {
                if (s.sData[i]==',') { // neither hDim nore vDim have been specified
                    if (j > 0) {
                        break;
                    }
                    term = s.Cut(1,i-1);
                    hDim = round(ProcessNumericArgument (&term,theP));
                    j    = i+1;
                }
            }

            if (j) { // both hDim and vDim specified
                term = s.Cut(j,i-1);
                vDim = ProcessNumericArgument (&term,theP);
            } else { // only one dim specified, matrix assumed to be square
                term = s.Cut(1,i-1);
                hDim = ProcessNumericArgument (&term,theP);
                vDim = hDim;
            }

            if (hDim<=0 || vDim<=0) {
                return;
            }

            if (isNumeric) {
                CreateMatrix (this, hDim, vDim, true, true, false);
            } else {
                CreateMatrix (this, hDim, vDim, true, false, true);
            }

            // read the terms now

            for (; i<s.sLength; i++) {
                if (s.sData[i]=='{') {
                    hPos = -1;
                    vPos = -1;
                    k    = i+1;

                    for (j=i+1; j<s.sLength && s.sData[j]!='}'; j++) {
                        j = s.FindTerminator (j, terminators);

                        if (j<0) {
                            WarnError ("Unterminated matrix definition");
                            return;
                        }

                        if (s.sData[j]==',') {
                            term = s.Cut (s.FirstNonSpaceIndex(k,j-1,1),j-1);
                            _Formula coordF (term,theP);
                            _Parameter coordV = coordF.Compute()->Value();
                            if (hPos == -1) {
                                hPos = coordV;
                            } else {
                                vPos = coordV;
                            }
                            k = j+1;
                        } else {
                            j--;
                        }
                    }

                    if (hPos <0 || vPos<0 || hPos>=hDim || vPos>=vDim)
                        // bad index
                    {
                        MatrixIndexError (hPos,vPos,hDim,vDim);
                        return;
                    }

                    term = s.Cut(k,j-1); // read the element

                    if (isNumeric) {
                        if ((term.sLength == 1) && (term.sData[0]=='*')) {
                            term = empty;    // dummy element in probability matrix
                        }

                        (*this)[vDim*hPos+vPos];
                        k = Hash (hPos,vPos);
                        theData[k]=term.toNum();
                    } else {
                        if ((term.sLength == 1) && (term.sData[0]=='*')) {
                            term = empty;    // dummy element in probability matrix
                        }

                        _Formula theTerm (term,theP);
                        if (isAConstant) { // there is hope that this matrix is of numbers
                            isAConstant = theTerm.IsAConstant();
                        }

                        (*this)[vDim*hPos+vPos];
                        k = Hash (hPos,vPos);
                        ((_Formula**)theData)[k]=(_Formula*)theTerm.makeDynamic();
                    }
                    i = j;
                }
            }
        } // end else

        if (!isNumeric) {
            storageType = 2; // formula elements
            checkParameter (ANAL_COMP_FLAG, ANALYTIC_COMPUTATION_FLAG, 0);
            if ((ANALYTIC_COMPUTATION_FLAG)&&!isAConstant) {
                ConvertFormulas2Poly (false);
            }

            if (isAConstant) { // a matrix of numbers - store as such
                Evaluate ();
            }
            AmISparse();
        }
    }
}


//_____________________________________________________________________________________________

_Matrix::_Matrix (_Matrix& m)
{
    DuplicateMatrix (this, &m);
}

//_____________________________________________________________________________________________

_Matrix::_Matrix (_SimpleList& sl, long colArg)
{
    if (sl.lLength) {
        if (colArg > 0 && colArg < sl.lLength) {
            CreateMatrix (this, sl.lLength/colArg + colArg*(sl.lLength%colArg > 0), colArg,     false, true, false);
        } else {
            CreateMatrix (this, 1, sl.lLength,  false, true, false);
        }
        for (long k=0; k<sl.lLength; k++) {
            theData[k] = sl.lData[k];
        }
    } else {
        Initialize();
    }
}

//_____________________________________________________________________________________________

_Matrix::_Matrix (_Parameter* inList, unsigned long rows, unsigned long columns)
{
    CreateMatrix (this, rows, columns, false, true, false);
    for (long k = 0; k < rows*columns; k++) {
        theData[k] = inList[k];
    }
}


//_____________________________________________________________________________________________

_Matrix::_Matrix (_List& sl)
// list of strings
{
    if (sl.lLength) {
        CreateMatrix     (this, 1, sl.lLength,  false, true, false);
        _Constant        hi (0.),
                         vi;

        for (long k=0; k<sl.lLength; k++) {
            _FString  *choiceString = new _FString (*(_String*) sl(k));
            _Formula  sf (choiceString);
            vi.SetValue (k);
            MStore(&hi,&vi,sf);
            //DeleteObject (choiceString);
        }
    } else {
        Initialize();
    }
}

//_____________________________________________________________________________________________

void    _Matrix:: ScanForVariables(_AVLList& theReceptacle, bool inclG)
{
    ScanForVariables2 (theReceptacle, inclG, -1);
}
//_____________________________________________________________________________________________

void    _Matrix:: ScanForVariables2(_AVLList& theReceptacle, bool inclG, long modelID, bool inclCat)
{
    if (storageType == 2) { // a formula based matrix, there is stuff to do
        if (modelID >= 0) {
            _AssociativeList*      definedCache = nil;
            _Variable*             cachedDeps = FetchVar(LocateVarByName (CACHE_FORMULA_DEPENDANCY));

            if (cachedDeps && cachedDeps->ObjectClass () == ASSOCIATIVE_LIST)
                // 20100316 SLKP: I am pretty sure this is broken...
            {
                definedCache = (_AssociativeList*)cachedDeps->GetValue();
                _String     matrixKey (modelID);
                _Matrix*    cachedValues = (_Matrix*)definedCache->GetByKey (matrixKey,MATRIX);

                if (cachedValues == nil) {
                    _Formula ** theFormulas = (_Formula**)theData;

                    _SimpleList sl1,
                                sl2;
                    _AVLList    a1 (&sl1),
                                a2 (&sl2);

                    if (theIndex) {
                        for (long i = 0; i<lDim; i++)
                            if (IsNonEmpty(i)) {
                                theFormulas[i]->ScanFForVariables(a1,false);
                                theFormulas[i]->ScanFForVariables(a2,true);
                            }
                    } else
                        for (long i = 0; i<lDim; i++)
                            if (theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                                theFormulas[i]->ScanFForVariables(a1,false);
                                theFormulas[i]->ScanFForVariables(a2,true);
                            }

                    a1.ReorderList();
                    a2.ReorderList();

                    cachedValues = new _Matrix (2,sl2.lLength,false,true);
                    checkPointer (cachedValues);

                    for (long k=0; k<sl1.lLength; k++) {
                        cachedValues->theData[k] = sl1.lData[k];
                    }
                    {
                        for (long k=sl1.lLength; k<sl2.lLength; k++) {
                            cachedValues->theData[k] = -1.;
                        }
                    }
                    {
                        for (long k=0; k<sl2.lLength; k++) {
                            cachedValues->theData[k+sl2.lLength] = sl2.lData[k];
                        }
                    }

                    _FString aKey (matrixKey,false);

                    definedCache->MStore (&aKey, cachedValues, false);

                }

                long colCount = cachedValues->GetVDim(),
                     rowIndex = inclG?colCount:0;

                for (long k=0; k<colCount; k++,rowIndex++) {
                    long vI = cachedValues->theData[rowIndex];
                    if (vI >= 0) {
                        theReceptacle.Insert ((BaseRef)vI);
                    } else {
                        break;
                    }
                }

                return;
            }

        }

        _Formula ** theFormulas = (_Formula**)theData;

        if (theIndex) {
            for (long i = 0; i<lDim; i++)
                if (IsNonEmpty(i)) {
                    theFormulas[i]->ScanFForVariables(theReceptacle,inclG,false,inclCat);
                }
        } else
            for (long i = 0; i<lDim; i++) {
                if (theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                    theFormulas[i]->ScanFForVariables (theReceptacle,inclG,false,inclCat);
                }
            }
    } else if (storageType == 0) { // a polynomial based matrix, there is stuff to do
        _MathObject ** thePoly = (_MathObject**)theData;
        if (theIndex)
            for (long i = 0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    thePoly[i]->ScanForVariables(theReceptacle,inclG);
                }
            }
        else
            for (long i = 0; i<lDim; i++) {
                if (thePoly[i]!=ZEROPOINTER) {
                    thePoly[i]->ScanForVariables (theReceptacle,inclG);
                }
            }
    }

}

//_____________________________________________________________________________________________

bool    _Matrix::IsConstant(void)
{
    if (storageType == 1) {
        return true;
    }

    if (storageType == 2) { // a formula based matrix, there is stuff to do
        _Formula ** theFormulas = (_Formula**)theData;
        if (theIndex) {
            for (long i = 0; i<lDim; i++)
                if (IsNonEmpty(i) && !theFormulas[i]->IsConstant()) {
                    return false;
                }
        } else
            for (long i = 0; i<lDim; i++)
                if (theFormulas[i]!=(_Formula*)ZEROPOINTER && !theFormulas[i]->IsConstant()) {
                    return false;
                }

        return true;

    }
    return false;
}

//_____________________________________________________________________________________________

bool        _Matrix::ProcessFormulas (long& stackLength, _SimpleList& varList,   _SimpleList& newFormulas,
                                      _SimpleList& references, _AVLListX& flaStrings,
                                      bool runAll, _Matrix * stencil)
{
    _Formula *      thisFormula = nil;
    _Formula **     theFormulas = (_Formula**)theData;

    bool isGood = true;

    if (theIndex) {
        for (long i = 0; i<lDim; i++) {
            long cellIndex = theIndex [i];
            if (cellIndex>-1) {
                if (stencil && CheckEqual(stencil->theData[cellIndex],0.0)) {
                    references << -1;
                    continue;
                }
                thisFormula = theFormulas[i];

                if (runAll || thisFormula->AmISimple(stackLength,varList)) {
                    _String * flaString = (_String*)thisFormula->toStr(nil,true);
                    long      fref = flaStrings.Insert(flaString,newFormulas.lLength);
                    if (fref < 0) {
                        references << flaStrings.GetXtra (-fref-1);
                        DeleteObject (flaString);
                    } else {
                        newFormulas << (long)thisFormula;
                        references << fref;
                    }

                } else {
                    isGood = false;
                    break;
                }
            } else {
                references << -1;
            }
        }
    } else {
        for (long i = 0; i<lDim; i++) {
            if ((theFormulas[i]!=(_Formula*)ZEROPOINTER)&&(!theFormulas[i]->IsEmpty())) {
                thisFormula = theFormulas[i];

                if (stencil && CheckEqual(stencil->theData[i],0.0)) {
                    references << -1;
                    continue;
                }

                if (runAll || thisFormula->AmISimple(stackLength,varList)) {
                    _String * flaString = (_String*)thisFormula->toStr(nil,true);
                    long      fref = flaStrings.Insert(flaString,newFormulas.lLength);
                    if (fref < 0) {
                        references << flaStrings.GetXtra (-fref-1);
                        DeleteObject (flaString);
                    } else {
                        newFormulas << (long)thisFormula;
                        references << fref;
                    }
                } else {
                    isGood = false;
                    break;
                }
            } else {
                references << -1;
            }
        }
    }
    return isGood;
}

//_____________________________________________________________________________________________
_Matrix*        _Matrix::branchLengthStencil (void)
{
    _Matrix * stencil = (_Matrix*)FetchObjectFromVariableByType (&BRANCH_LENGTH_STENCIL,MATRIX);
    if (stencil) {
        if (stencil->storageType==1 && stencil->hDim==stencil->vDim && stencil->hDim == hDim) {
            stencil->CheckIfSparseEnough (true);
        } else {
            stencil = nil;
        }
    }

    return stencil;
}

//_____________________________________________________________________________________________
_String*        _Matrix::BranchLengthExpression (_Matrix* baseFreqs, bool mbf)
{
    if (storageType == 2) {

        long            stackLength = 0;

        _SimpleList     varList,
                        newFormulas,
                        references;

        _List           flaStringsL;
        _AVLListX       flaStrings(&flaStringsL);
        _Matrix*        stencil = branchLengthStencil();

        ProcessFormulas (stackLength,varList,newFormulas,references,flaStrings,true,stencil);

        _String * sendMeBack = new _String(128L, true);
        if (baseFreqs->storageType == 1)
            // numerical base frequencies
        {
            _Matrix   multipliersByRate (newFormulas.lLength,1,false,true);
            for (long i = 0; i<lDim; i++) {
                long thisRef = references.lData[i];
                if (thisRef>=0) {
                    long cellIndex = i;
                    if (theIndex) {
                        cellIndex = theIndex[i];
                    }

                    multipliersByRate.theData[thisRef] += (*baseFreqs)(cellIndex/vDim,0) *
                                                          (mbf?(*baseFreqs)(cellIndex%vDim,0):1.0);
                }
            }


            bool    firstDone = false;
            for (long k=0; k<newFormulas.lLength; k++) {
                if (!CheckEqual(multipliersByRate.theData[k],0.0)) {
                    if (firstDone) {
                        (*sendMeBack) << '+';
                    }
                    _String * fStr = (_String*)flaStringsL(k);
                    (*sendMeBack) << '(';
                    (*sendMeBack) << fStr;
                    (*sendMeBack) << ")*";
                    (*sendMeBack) << _String(multipliersByRate.theData[k]);
                    firstDone = true;
                }
            }
        } else if (baseFreqs->storageType == 2)
            // formula-based equilibrium frequencies
        {
            _List   freqFla,
                    multipliersByRate;

            for (long k=0; k<newFormulas.lLength; k++) {
                multipliersByRate.AppendNewInstance(new _String (128L,true));
            }

            for (long i = 0; i<hDim; i++) {
                freqFla.AppendNewInstance ((_String*)baseFreqs->GetFormula(i,0)->toStr(nil,true));
            }

            for (long i = 0; i<lDim; i++) {
                long thisRef = references.lData[i];
                if (thisRef>=0) {
                    _String * thisAdder = (_String*)multipliersByRate(thisRef);
                    if (thisAdder->sLength) {
                        (*thisAdder) << '+';
                    }

                    long cellIndex = i;
                    if (theIndex) {
                        cellIndex = theIndex[i];
                    }

                    (*thisAdder) << '(';
                    if (mbf) {

                        (*thisAdder) << (_String*)freqFla(cellIndex%vDim);
                        (*thisAdder) << ")*(";
                    }

                    (*thisAdder) << (_String*)freqFla(cellIndex/vDim);
                    (*thisAdder) << ')';
                }
            }

            for (long k=0; k<newFormulas.lLength; k++) {
                ((_String*)multipliersByRate(k))->Finalize();
            }

            for (long k=0; k<newFormulas.lLength; k++) {
                if (k) {
                    (*sendMeBack) << '+';
                }
                (*sendMeBack) << '(';
                (*sendMeBack) << (_String*)flaStringsL(k);
                (*sendMeBack) << ")*(";
                (*sendMeBack) << (_String*)multipliersByRate(k);
                (*sendMeBack) << ')';
            }
        }
        sendMeBack->Finalize();
        if (sendMeBack->sLength) {
            _Formula        blF (*sendMeBack);
            _Polynomial*    isPoly = (_Polynomial*)blF.ConstructPolynomial();
            if (isPoly) {
                DeleteObject (sendMeBack);
                sendMeBack = (_String*)isPoly->toStr();
            }
        }
        return sendMeBack;
    }
    return new _String;
}

//_____________________________________________________________________________________________
void        _Matrix::MakeMeSimple (void)
{
    if (storageType == 2) {
        long            stackLength = 0;
        bool            isGood      = true;

        _SimpleList     varList,
                        newFormulas,
                        references;

        _List           flaStringsL;
        _AVLListX       flaStrings(&flaStringsL);

        isGood      =   ProcessFormulas (stackLength,varList,newFormulas,references,flaStrings);

        if (isGood) {
            storageType = 3;

            for (long k = 0; k < newFormulas.lLength; k++) {
                ((_Formula*)newFormulas.lData[k])->ConvertToSimple(varList);
            }

            cmd                         = new _CompiledMatrixData;
            cmd->varIndex.Duplicate (&varList);
            cmd->theStack               = (_SimpleFormulaDatum*)MatrixMemAllocate (stackLength*sizeof(_SimpleFormulaDatum));
            cmd->varValues              = (_SimpleFormulaDatum*)MatrixMemAllocate ((cmd->varIndex.lLength>0?varList.lLength:1)*sizeof(_SimpleFormulaDatum));
            cmd->formulaRefs            = references.lData;
            references.lData            = nil;
            cmd->formulaValues          = new _Parameter [newFormulas.lLength];
            checkPointer (cmd->formulaValues);
            cmd->formulasToEval.Duplicate (&newFormulas);
        }

    }
}
//_____________________________________________________________________________________________
void        _Matrix::MakeMeGeneral (void)
{
    if (storageType == 3) {
        for (long k = 0; k < cmd->formulasToEval.lLength; k++) {
            ((_Formula*)cmd->formulasToEval.lData[k])->ConvertFromSimple(cmd->varIndex);
        }

        delete [] cmd->formulaValues;
        free   (cmd->formulaRefs);

        MatrixMemFree   (cmd->theStack);
        MatrixMemFree   (cmd->varValues);
        delete          (cmd);
        cmd             = nil;
        storageType     = 2;
    }
}
//_____________________________________________________________________________________________
_PMathObj   _Matrix::Evaluate (bool replace)
// evaluate the matrix  overwriting (or not) the old one
{
    _Matrix result (hDim, vDim, bool (theIndex), true);

    if (storageType == 2) {
        _PMathObj formValue = nil;
        _Formula ** theFormulas = (_Formula**)theData;
        if (theIndex) {
            for (long i = 0; i<lDim; i++) {
                //long k =
                if (theIndex[i]!=-1) {
                    formValue = theFormulas[i]->Compute();
                    if (formValue) {
                        result[HashBack(i)] = formValue->Value();
                        //DeleteObject (formValue);
                    } else {
                        result[HashBack(i)] = 0;
                    }
                }
            }
            // check for probablilty matrices * fillers
            if ((hDim==vDim)&&(!replace))
                for (long i = 0; i<hDim; i++) {
                    long k = Hash(i,i);
                    if ((k>=0)&&theFormulas[k]->IsEmpty()) {
                        _Parameter *st = &result[k];
                        *st=0;
                        for (long j = 0; j<vDim; j++) {
                            if (j==i) {
                                continue;
                            }
                            *st-=result(i,j);
                        }
                    } else if (k<0) {
                        _Parameter *st = &result[i*vDim+i];
                        *st=0;
                        for (long j = 0; j<vDim; j++) {
                            if (j==i) {
                                continue;
                            }
                            *st-=result(i,j);
                        }
                    }
                }
        } else {
            for (long i = 0; i<lDim; i++) {
                if (theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                    formValue = theFormulas[i]->Compute();
                    if (formValue && formValue->ObjectClass() == NUMBER) {
                        result.theData[i] = formValue->Value();
                        //DeleteObject (formValue);
                    } else {
                        result.theData[i] = 0;
                    }
                }
            }
            // check for probablilty matrices * fillers

            if ((hDim==vDim)&&(!replace))
                for (long i = 0; i<lDim; i+=vDim+1) {
                    if (theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                        if (theFormulas[i]->IsEmpty()) {
                            _Parameter st = 0;
                            long k = i/vDim,j;
                            for (j = k*vDim; j<k*vDim+k; j++) {
                                st-=result.theData[j];
                            }
                            for (j = k*vDim+k+1; j<(k+1)*vDim; j++) {
                                st-=result.theData[j];
                            }
                            result.theData[i] = st;
                        }
                    }
                }
        }
    }
    if (storageType == 0) {
        _PMathObj polValue = nil;
        _MathObject ** thePoly = (_MathObject**)theData;
        if (theIndex) {
            for (long i = 0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    polValue = thePoly[i]->Compute();
                    if (polValue) {
                        result[HashBack(i)] = polValue->Value();
                        DeleteObject (polValue);
                    } else {
                        result[i] = 0;
                    }
                }
            }

        } else {
            for (long i = 0; i<lDim; i++) {
                if (thePoly[i]!=(_MathObject*)ZEROPOINTER) {
                    polValue = thePoly[i]->Compute();
                    if (polValue) {
                        result[i] = polValue->Value();
                        DeleteObject (polValue);
                    } else {
                        result[i] = 0;
                    }
                }
            }
        }
    }
    if (replace) {
        *this = result;
    } else {
        return (_PMathObj)result.makeDynamic();
    }
    return nil;
}

//_____________________________________________________________________________________________
void        _Matrix::ConvertToSimpleList (_SimpleList & sl)
{
    sl.Clear();
    if (storageType == 1) {
        sl.RequestSpace (hDim*vDim+1);

        for (long i=0; i<hDim; i++)
            for (long j=0; j<vDim; j++) {
                sl << (*this)(i,j);
            }
    }
}

//_____________________________________________________________________________________________
bool        _Matrix::IsAStringMatrix (void)
// check if a formula matrix contains strings
{
    if (storageType == 2) {
        _PMathObj   formValue = nil;
        _Formula ** theFormulas = (_Formula**)theData;
        if (theIndex)
            for (long i = 0; i<lDim; i++) {
                if (theIndex[i]!=-1 && !(theFormulas[i]->IsEmpty())) {
                    formValue = theFormulas[i]->Compute();
                    if (formValue) {
                        return formValue->ObjectClass() == STRING;
                    }
                }
            }
        else
            for (long i = 0; i<lDim; i++)
                if (theFormulas[i]!=(_Formula*)ZEROPOINTER && !(theFormulas[i]->IsEmpty())) {
                    formValue = theFormulas[i]->Compute();
                    if (formValue) {
                        return formValue->ObjectClass() == STRING;
                    }
                }

    }
    return false;
}

//_____________________________________________________________________________________________
void        _Matrix::FillInList (_List& fillMe, bool doNumeric)
// check if a formula matrix contains strings
{
    if (storageType == _FORMULA_TYPE)
        for (long r=0; r<hDim; r++)
            for (long c=0; c<vDim; c++) {
                _Formula * entryFla = GetFormula(r,c);
                if (entryFla) {
                    _PMathObj computedValue = entryFla->Compute();
                    if (computedValue)
                        if (computedValue->ObjectClass() == STRING) {
                            fillMe && ((_FString*)computedValue)->theString;
                        } else {
                            fillMe.Clear();
                            return;
                        }
                }
            }
    else {
        if (doNumeric && storageType == _NUMERICAL_TYPE) {
            for (long r=0; r<hDim; r++)
                for (long c=0; c<vDim; c++) {
                    fillMe.AppendNewInstance (new _String ((*this)(r,c)));
                }
        }
    }
}

//_____________________________________________________________________________________________
_PMathObj   _Matrix::EvaluateSimple (void)
// evaluate the matrix  overwriting the old one
{
    _Matrix * result = new _Matrix (hDim, vDim, bool (theIndex), true);
    checkPointer (result);


    if (cmd->varIndex.lLength) {
        for (long i=0; i<cmd->varIndex.lLength; i++) {
            _Variable* curVar = LocateVar(cmd->varIndex.lData[i]);
            if (curVar->ObjectClass () != MATRIX) {
                if (curVar->IsIndependent()) {
                    cmd->varValues[i].value = LocateVar (cmd->varIndex.lData[i])->Value();
                } else {
                    cmd->varValues[i].value = LocateVar (cmd->varIndex.lData[i])->Compute()->Value();
                }
            } else {
                cmd->varValues[i].reference = (Ptr)((_Matrix*)LocateVar (cmd->varIndex.lData[i])->Compute())->theData;
            }
        }
    }


    for (long f = 0; f < cmd->formulasToEval.lLength; f++) {
        cmd->formulaValues [f] = ((_Formula*)cmd->formulasToEval.lData[f])->ComputeSimple(cmd->theStack, cmd->varValues);
        /*if (terminateExecution)
        {
            ((_Formula*)cmd->formulasToEval.lData[f])->ConvertFromSimple(cmd->varIndex);
            _String * s = (_String*)((_Formula*)cmd->formulasToEval.lData[f])->toStr();
            WarnError (*s);
            DeleteObject (s);
            return result;
        }*/
    }

    long * fidx = cmd->formulaRefs;

    if (theIndex) {
        result->lDim = lDim;
        result->bufferPerRow = bufferPerRow;
        result->overflowBuffer = overflowBuffer;
        result->allocationBlock = allocationBlock;
        result->theIndex = (long*)MemReallocate((Ptr)result->theIndex,sizeof(long)*lDim);
        result->theData = (_Parameter*)MemReallocate ((Ptr)result->theData,sizeof(_Parameter)*lDim);

        /*memcpy (result->theIndex,theIndex,sizeof(long)*lDim);*/




        for (long i = 0; i<lDim; i++) {
            long idx = theIndex[i];

            if (idx != -1) {
                result->theData[i] = cmd->formulaValues[fidx[i]];
            }

            result->theIndex[i] = idx;
        }

        /*for (long i = 0; i<lDim; i++)
        {
            if (theIndex[i]!=-1)
            {
                formValue = theFormulas[i]->ComputeSimple(cmd->theStack, cmd->varValues);
                result.theData[i] = formValue;
            }
        } */

        if (hDim==vDim) {
            _Parameter* diagStorage = new _Parameter [hDim];
            checkPointer ((Ptr)diagStorage);
            {
                for (long i = 0; i<hDim; i++) {
                    diagStorage[i] = 0.0;
                }
            }
            for (long i = 0; i<lDim; i++) {
                long k = result->theIndex[i];
                if (k!=-1) {
                    diagStorage[k/hDim] -= result->theData[i];
                }
            }
            {
                for (long i = 0; i<hDim; i++) {
                    (*result)[i*hDim+i] = diagStorage[i];
                }
            }
            delete [] diagStorage;
        }
    } else {
        /*long i;
        for (i = 0; i<lDim; i++)
        {
            if (theFormulas[i]!=(_Formula*)ZEROPOINTER)
            {
                formValue = theFormulas[i]->ComputeSimple(cmd->theStack,cmd->varValues);
                result.theData[i] = formValue;
                //break;
            }
        }       */

        for (long i = 0; i<lDim; i++) {
            if (fidx[i]>= 0) {
                result->theData[i] = cmd->formulaValues[fidx[i]];
            }
        }

        if (hDim==vDim)
            for (long i = 0; i<lDim; i+=vDim+1) {
                if (fidx[i] < 0) { // mod Aug 2 2005
                    //if (theFormulas[i]->IsEmpty())
                    //{
                    _Parameter st = 0;
                    long k = i/vDim,j;
                    for (j = k*vDim; j<k*vDim+k; j++) {
                        st-=result->theData[j];
                    }

                    for (j = k*vDim+k+1; j<(k+1)*vDim; j++) {
                        st-=result->theData[j];
                    }

                    result->theData[i] = st;
                    //}
                }
            }
    }
    //return (_PMathObj)result.makeDynamic();
    return result;
}
//_____________________________________________________________________________________________
void    _Matrix::ClearFormulae (void)
{
    _Formula ** theFormulas = (_Formula**)theData;
    if (theIndex) {
        for (long i = 0; i<lDim; i++) {
            if (IsNonEmpty(i)) {
                delete (theFormulas[i]);
            }
        }
    } else
        for (long i = 0; i<lDim; i++) {
            if (theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                delete (theFormulas[i]);
            }
        }
}

//_____________________________________________________________________________________________
void    _Matrix::ClearObjects (void)
{
    _MathObject ** thePolys = (_MathObject**)theData;
    if (theIndex) {
        for (long i = 0; i<lDim; i++) {
            if (IsNonEmpty(i)) {
                DeleteObject (thePolys[i]);
            }
        }
    } else
        for (long i = 0; i<lDim; i++) {
            if (thePolys[i]!=(_MathObject*)ZEROPOINTER) {
                DeleteObject (thePolys[i]);
            }
        }
}

//_____________________________________________________________________________________________

void    _Matrix::Clear (void)
{
    DeleteObject (theValue);
    if (storageType == 2) { // has formulas in it - must delete
        ClearFormulae();
    }
    if (storageType == 0) { // has objects in it - must delete
        ClearObjects();
    }
    if (theIndex) {
        MatrixMemFree (theIndex);
        theIndex = nil;
    }
    if (theData) {
        MatrixMemFree (theData);
        hDim = vDim = 0;
        theData = nil;
    }

}

//_____________________________________________________________________________________________

void    _Matrix::Resize (long newH)
{
    if (newH >= 0 && newH != hDim && storageType == 1 && theIndex == nil) {
        hDim = newH;
        lDim = newH*vDim;

        if (theData) {
            theData = (_Parameter*) MemReallocate ((Ptr)theData,sizeof (_Parameter)*lDim);
        } else {
            theData = (_Parameter*) MemAllocate (sizeof (_Parameter)*lDim);
        }

    }
}

//_____________________________________________________________________________________________

_Matrix::~_Matrix (void)
{
    Clear();
}

//_____________________________________________________________________________________________

void    _Matrix::operator = (_Matrix& m)
{
    Clear();
    DuplicateMatrix (this, &m);
}

//_____________________________________________________________________________________________

void    _Matrix::operator = (_Matrix* m)
{
    Clear();
    DuplicateMatrix (this, m);
}


//_____________________________________________________________________________________________
_Parameter _Matrix::AbsValue (void)
{
    if (storageType == 1 && (hDim==1 || vDim == 1)) {
        _Parameter norm = 0.;
        if (theIndex) {
            for (long k = 0; k<lDim; k++)
                if (long i = theIndex[k] >= 0) {
                    norm += theData[i]*theData[i];
                }

            norm = sqrt(norm);
        } else {
            for (long k = 0; k<lDim; k++) {
                norm += theData[k]*theData[k];
            }
            norm = sqrt(norm);
        }
        return norm;
    }

    return 0.;
}

//_____________________________________________________________________________________________
_PMathObj _Matrix::Abs (void)
{
    if (storageType == 1 && (hDim==1 || vDim == 1)) {
        return new _Constant (AbsValue());
    }
    return new _Constant(MaxElement());

}

//_____________________________________________________________________________________________

void    _Matrix::Add  (_Matrix& storage, _Matrix& secondArg, bool subtract)
// addition operation on matrices
// internal function

{

    // check matrix dimensions to ensure that they are addable
    if (!((hDim==secondArg.hDim)&&(storage.hDim==secondArg.hDim)&&(vDim==secondArg.vDim)&&(storage.vDim==secondArg.vDim))) {
        _String  errMsg = _String ("Incompatible dimensions when trying to add or subtract matrices: first argument was a ") & _String (hDim) & 'x'
                          & _String (vDim) & " matrix and the second was a "& _String (secondArg.hDim) & 'x'  & _String (secondArg.vDim) & " matrix.";
        WarnError (errMsg);
        /*warnError( -103);*/
        return;
    }

    if (storageType == 1) {
        if (&storage != this) { // not an add&store operation
            // copy *this to storage
            if (theIndex) { //sparse matrix
                for (long i = 0; i<lDim; i++) {
                    long k = theIndex[i];
                    if (k!=-1) {
                        storage[k] = theData[i];
                    }
                }
            } else { // dense matrix
                _Parameter * _hprestrict_ stData = storage.fastIndex();
                for (long i = 0; i<lDim; i++) {
                    stData[i] = theData[i];
                }
            }
        }

        if (secondArg.theIndex) { //sparse matrix
            if (storage.theIndex) {
                if (subtract) {
                    for (long i = 0; i<secondArg.lDim; i++) {
                        long k = secondArg.theIndex[i];
                        if (k!=-1) {
                            storage[k]-=secondArg.theData[i];
                        }
                    }
                } else {
                    for (long i = 0; i<secondArg.lDim; i++) {
                        long k = secondArg.theIndex[i];
                        if (k!=-1) {
                            storage[k]+=secondArg.theData[i];
                        }
                    }
                }
            } else {
                if (subtract) {
                    for (long i = 0; i<secondArg.lDim; i++) {
                        long k = secondArg.theIndex[i];
                        if (k!=-1) {
                            storage.theData[k]-=secondArg.theData[i];
                        }
                    }
                } else {
                    for (long i = 0; i<secondArg.lDim; i++) {
                        long k = secondArg.theIndex[i];
                        if (k!=-1) {
                            storage.theData[k]+=secondArg.theData[i];
                        }
                    }
                }
            }

        } else {
            _Parameter _hprestrict_ * argData = secondArg.theData;
            _Parameter _hprestrict_ * stData = storage.theData;
            if (subtract)
                for (long idx = 0; idx < secondArg.lDim; idx++) {
                    stData[idx]-=argData[idx];
                }
            else
                for (long idx = 0; idx < secondArg.lDim; idx++) {
                    stData[idx]+=argData[idx];
                }
        }
    } else

        if (storageType == 0) {
            long i;
            if (&storage != this) { // not an add&store operation
                /*              if (theIndex) //sparse matrix
                                {
                                    for (i = 0; i<lDim; i++)
                                        if (IsNonEmpty(i))
                                            storage.StoreObject(theIndex[i],GetMatrixObject(i),true);
                                }
                                else // normal matrix
                                {
                                    for (i = 0; i<lDim; i++)
                                        storage.StoreObject(i,GetMatrixObject(i),true);
                                }*/
            }

            if (secondArg.theIndex) { //sparse matrix
                if (theIndex) { // both matrices are sparse
                    if (subtract) {
                        for (i = 0; i<secondArg.lDim; i++)
                            if (secondArg.IsNonEmpty(i)) {
                                long hb =secondArg.HashBack (i), h = Hash (hb/vDim, hb%vDim);
                                if (h<0) { // empty slot in matrix 1
                                    storage.StoreObject (hb,secondArg.GetMatrixObject(i)->Minus());
                                } else {
                                    storage.StoreObject (hb, GetMatrixObject(h)->Sub (secondArg.GetMatrixObject(i)));
                                }
                            }
                    } else {
                        for (i = 0; i<secondArg.lDim; i++)
                            if (secondArg.IsNonEmpty(i)) {
                                long hb =secondArg.HashBack (i), h = Hash (hb/vDim, hb%vDim);
                                if (h<0) { // empty slot in matrix 1
                                    storage.StoreObject (hb,secondArg.GetMatrixObject(i),true);
                                } else {
                                    storage.StoreObject (hb,GetMatrixObject(h)->Add (secondArg.GetMatrixObject(i)));
                                }
                            }
                    }
                } else { // *this is not sparse
                    DuplicateMatrix(&storage,this);
                    if (subtract) {
                        for (i = 0; i<secondArg.lDim; i++)
                            if (secondArg.IsNonEmpty(i)) {
                                long p = secondArg.HashBack (i);
                                if (CheckObject(p)) {
                                    storage.StoreObject (p,GetMatrixObject(p)->Sub(secondArg.GetMatrixObject(i)));
                                } else {
                                    storage.StoreObject (p,secondArg.GetMatrixObject(i)->Minus());
                                }
                            }
                    } else {
                        for (i = 0; i<secondArg.lDim; i++)
                            if (secondArg.IsNonEmpty(i)) {
                                long p = secondArg.HashBack (i);
                                if (CheckObject(p)) {
                                    storage.StoreObject (p,GetMatrixObject(p)->Add(secondArg.GetMatrixObject(i)));
                                } else {
                                    storage.StoreObject (p,secondArg.GetMatrixObject(i),true);
                                }
                            }
                    }
                }
            } else { // secondarg isn't sparse - storage must also be non-sparse
                if (storage.theIndex) { // storage is sparse - oops
                    storage.CheckIfSparseEnough(true);    // force to non-sparse storage
                }
                if (!theIndex) { // * this is not sparse
                    _PMathObj tempP;
                    DuplicateMatrix(&storage,this);
                    if (subtract) {
                        for (i = 0; i<secondArg.lDim; i++) {
                            tempP = secondArg.GetMatrixObject(i);
                            if (tempP) {
                                if (CheckObject(i)) {
                                    storage.StoreObject(i,GetMatrixObject(i)->Sub(tempP));
                                } else {
                                    storage.StoreObject(i,tempP->Minus());
                                }
                            }
                        }
                    } else {
                        for (i = 0; i<secondArg.lDim; i++) {
                            tempP = secondArg.GetMatrixObject(i);
                            if (tempP) {
                                if (CheckObject(i)) {
                                    storage.StoreObject(i,GetMatrixObject(i)->Add(tempP));
                                } else {
                                    storage.StoreObject(i,tempP,true);
                                }
                            }
                        }
                    }
                } else { // *this is sparse
                    _PMathObj tempP;
                    long h;
                    if (subtract) {
                        for (i = 0; i<secondArg.lDim; i++) {
                            tempP = secondArg.GetMatrixObject(i);
                            if (tempP) {
                                h = Hash (i/hDim,i%hDim);
                                if (h>=0) {
                                    storage.StoreObject(i,GetMatrixObject(h)->Sub(tempP));
                                } else {
                                    storage.StoreObject(i,tempP->Minus());
                                }
                            }
                        }
                    } else {
                        for (i = 0; i<secondArg.lDim; i++) {
                            tempP = secondArg.GetMatrixObject(i);
                            if (tempP) {
                                h = Hash (i/hDim,i%hDim);
                                if (h>=0) {
                                    storage.StoreObject(i,GetMatrixObject(h)->Add(tempP));
                                } else {
                                    storage.StoreObject(i,tempP,true);
                                }
                            }
                        }
                    }
                }
            }
        }
    if (storage.theIndex) {
        storage.CheckIfSparseEnough();
    }
}

//_____________________________________________________________________________________________

bool    _Matrix::AddWithThreshold  (_Matrix& secondArg, _Parameter prec)
{
    bool res = true;
    if (secondArg.theIndex) { //sparse matrix
        long i,k;
        for (i = 0; res&&(i<secondArg.lDim); i++) {
            k = secondArg.theIndex[i];
            if (k!=-1) {
                if (secondArg.theData[i]/theData[k] > prec) {
                    res = false;
                }
                theData[k]+=secondArg.theData[i];
            }
        }
        for (; i<secondArg.lDim; i++) {
            k = secondArg.theIndex[i];
            if (k!=-1) {
                theData[k]+=secondArg.theData[i];
            }
        }
    } else {
        _Parameter* argData = secondArg.theData, *stData = theData,
                    *bound = theData+lDim;
        for (; res&&(stData!=bound); argData++, stData++) {
            if (*argData/ *stData> prec) {
                res = false;
            }
            *stData+=*argData;
        }
        for (; stData!=bound; argData++, stData++) {
            *stData+=*argData;
        }
    }
    return !res;
}

//_____________________________________________________________________________________________

void    _Matrix::Subtract  (_Matrix& storage, _Matrix& secondArg)
// subtraction operation on matrices
// internal function

{
    Add (storage,secondArg,true);
}

//_____________________________________________________________________________________________

void    _Matrix::Multiply  (_Matrix& storage, _Parameter c)
// multiply a matrix by a scalar
// internal function

{
    if (storageType == 1) { // numbers
        if (theIndex) {
            for (long k = 0; k < lDim; k++)
                if (storage.theIndex[k] != -1) {
                    storage.theData[k] = theData[k]*c;
                }
        } else
            for (long k = 0; k < lDim; k++) {
                storage.theData[k] = theData[k]*c;
            }
    } else {
        _Constant * cc = new _Constant (c);
        checkPointer (cc);

        if (storageType == 2) {
            _String    star ('*');
            _Operation * cOp = new _Operation (cc),
            * mOp = new _Operation (star,2);

            checkPointer (cOp);
            checkPointer (mOp);
            for (long i=0; i<lDim; i++)
                if (IsNonEmpty (i)) {
                    long h       = HashBack (i);
                    _Formula * f = GetFormula (h/vDim,h%vDim);
                    f->GetList() && cOp;
                    f->GetList() && mOp;
                }
        } else {
            if (storageType != 3) {
                if (theIndex)
                    //sparse matrix
                {
                    for (long i=0; i<lDim; i++)
                        if (IsNonEmpty (i)) {
                            storage.StoreObject (HashBack(i),GetMatrixObject(i)->Mult (cc));
                        }
                } else {
                    for (long i=0; i<lDim; i++)
                        if (IsNonEmpty (i)) {
                            storage.StoreObject (i,GetMatrixObject(i)->Mult (cc));
                        }
                }
            }
            DeleteObject (cc);
        }

    }
}


//_____________________________________________________________________________________________

void    _Matrix::Multiply  (_Matrix& storage, _Matrix& secondArg )
// multiplication operation on matrices
// internal function
// storage is assumed to NOT be *this

{
    _PMathObj tempP, tempP2;

    if ( !theIndex && !secondArg.theIndex)
        // simplest case of two non-sparse matrices - multiply in a straightforward way
    {
        if ( storageType == 0 && secondArg.storageType ==0) { // both matrices are polynomial in nature
            for (long i=0; i<hDim; i++)
                for (long j=i*secondArg.vDim; j<(i+1)*secondArg.vDim; j++) {
                    _MathObject* secTerm = secondArg.GetMatrixObject(j%secondArg.vDim), *firstTerm = GetMatrixObject (i*vDim);
                    if (firstTerm&&secTerm) {
                        storage.StoreObject (j,firstTerm->Mult (secTerm));
                    } else {
                        storage.StoreObject (j,new _Polynomial(0.0));
                    }
                    for (long k=i*vDim+1, l=j%secondArg.vDim+secondArg.vDim; k<(i+1)*vDim; k++, l+=secondArg.vDim) {
                        tempP = GetMatrixObject (k), tempP2 = secondArg.GetMatrixObject(l);
                        if (tempP&&tempP2) {
                            _MathObject* temp = tempP->Mult (tempP2);
                            storage.StoreObject (j,temp->Add(storage.GetMatrixObject(j)));
                            DeleteObject (temp);
                        }
                    }
                }
        } else {
            if ( hDim == vDim && secondArg.hDim == secondArg.vDim)
                /* two square dense matrices */
            {
                long cumulativeIndex = 0;

                _Parameter * row = theData;

#ifndef _SLKP_SSE_VECTORIZATION_
                if (vDim % 4 == 0)
                    // manual loop unroll
                {
                    for (long i=0; i<hDim; i++, row += vDim) {
                        for (long j=0; j<secondArg.vDim; j++) {
                            _Parameter resCell  = 0.0;

                            for (long k = 0, column = j; k < vDim; k+=4, column += (secondArg.vDim<<2))
                                resCell += row[k]   * secondArg.theData [column]                         +
                                           row[k+1] * secondArg.theData [column + secondArg.vDim ]       +
                                           row[k+2] * secondArg.theData [column + (secondArg.vDim << 1)] +
                                           row[k+3] * secondArg.theData [column + secondArg.vDim * 3];

                            storage.theData[cumulativeIndex++] = resCell;

                        }
                    }
                } else {
                    for (long i=0; i<hDim; i++, row += vDim) {
                        for (long j=0; j<secondArg.vDim; j++) {
                            _Parameter resCell = 0.0;

                            for (long k = 0, column = j; k < vDim; k++, column += secondArg.vDim) {
                                resCell += row[k] * secondArg.theData[column];
                            }

                            storage.theData[cumulativeIndex++] = resCell;
                        }
                    }
                }
#else
                secondArg.Transpose();
                for (long i=0; i<hDim; i++, row += vDim) {
                    for (long j=0; j<hDim; j++) {
                        _Parameter resCell  = 0.0;
                        for (long k = 0, column = j*hDim; k < vDim; k++, column ++) {
                            resCell += row[k] * secondArg.theData [column];
                        }

                        storage.theData[cumulativeIndex++] = resCell;
                    }
                }
                secondArg.Transpose();
#endif
            } else
                /* rectangular matrices */
            {
                /*long  off1 = secondArg.vDim,
                        off2 = 2*secondArg.vDim,
                        off3 = 3*secondArg.vDim,
                        off4 = 4*secondArg.vDim;*/

                _Parameter *rD,
                           *fD,
                           *sD,
                           *stop,
                           resCell;

                for (long i=0; i<hDim; i++) {
                    rD = storage.theData+i*secondArg.vDim;
                    for (long j=0; j<secondArg.vDim; j++,rD++) {
                        resCell = 0.0;
                        fD = theData+i*vDim;
                        stop = fD+vDim;
                        sD = secondArg.theData+j;
                        //for (k=0; k<vDim; k++, fD++, sD+=secondArg.vDim)
                        //for (; stop!=fD; fD++, sD+=secondArg.vDim)
                        for (; fD!=stop;)
                            //resCell+=*fD**sD;
                        {
                            resCell += *(fD++)**sD;
                            sD += secondArg.vDim;
                        }
                        *rD = resCell;
                    }
                }
            }
        }

    } else if (theIndex && !secondArg.theIndex) { // sparse multiplied by non-sparse
        if (storageType == 1 && secondArg.storageType ==1) { // both numeric
            if ( vDim == hDim && secondArg.vDim==secondArg.hDim) { // both square and same dimension
                long loopBound = vDim - vDim%4;

                for (long k=0; k<lDim; k++) { // loop over entries in the 1st matrix
                    long m = theIndex[k];
                    if  (m!=-1) { // non-zero
                        long i = m%vDim;

                        // this element will contribute to (r, c' = [0..vDim-1]) entries in the result matrix
                        // in the form of A_rc * B_cc'

                        _Parameter  value                           = theData[k];
                        _Parameter  _hprestrict_ *res               = storage.theData    + (m-i);
                        _Parameter  _hprestrict_ *secArg            = secondArg.theData  + i*vDim;

#ifndef _SLKP_SSE_VECTORIZATION_
                        for (long i = 0; i < loopBound; i+=4) {
                            res[i] += value * secArg[i];
                            res[i+1] += value * secArg[i+1];
                            res[i+2] += value * secArg[i+2];
                            res[i+3] += value * secArg[i+3];
                        }
                        for (long j = loopBound; j < vDim; j++) {
                            res[j]   += value * secArg[j];
                        }
#else
                        for (long i = 0; i < vDim; i++) {
                            res[i]   += value * secArg[i];
                        }
#endif

                    }
                }
            } else {
                for (long k=0; k<lDim; k++) {
                    long m = theIndex[k];
                    if (m!=-1) {
                        long i = m/vDim;
                        long j = m%vDim;
                        _Parameter c = theData[k];
                        _Parameter* stData = storage.theData+i*secondArg.vDim,
                                    * secArgData=secondArg.theData+j*secondArg.vDim,
                                      * stopper = secArgData+secondArg.vDim;
                        for (; secArgData!=stopper; stData++,secArgData++) {
                            *stData+=c**secArgData;
                        }
                    }
                }
            }
        } else { // polynomial entries
            for (long k=0; k<lDim; k++) {
                if (IsNonEmpty(k)) {
                    long i = theIndex[k]/vDim;
                    long j = theIndex[k]%vDim;
                    _MathObject* p = GetMatrixObject(k);
                    for (long l=j*secondArg.vDim, m=i*secondArg.vDim; l<(j+1)*secondArg.vDim; l++,m++) {
                        tempP = secondArg.GetMatrixObject (l);
                        if (!tempP) {
                            continue;
                        }
                        _MathObject* temp = p->Mult (secondArg.GetMatrixObject (l));
                        tempP = storage.GetMatrixObject(m);
                        if (tempP) {
                            storage.StoreObject (m, tempP->Add (temp));
                        } else {
                            storage.StoreObject (m, temp, true);
                        }
                        DeleteObject (temp);
                    }
                }
            }
        }

    } else if ( !theIndex && secondArg.theIndex)
        // non-sparse multiplied by sparse
    {
        if ( storageType == 1 && secondArg.storageType ==1) {
            if (vDim == hDim && secondArg.vDim==secondArg.hDim)
                // both are square matrices
            {
                for (long k=0; k<secondArg.lDim; k++) {
                    long m = secondArg.theIndex[k];
                    if (m!=-1) { // a non-zero value
                        // because r_ij = sum_k a_ik * b_kj
                        // a non-zero b_kj will contribute a_ik * b_kj to the a_ij cell of the result
                        // loop over i...

                        _Parameter c = secondArg.theData[k];

                        for (long cell = m%secondArg.vDim, secondCell = m/secondArg.vDim; cell < lDim; cell += vDim, secondCell += vDim) {
                            storage.theData[cell] += c * theData[secondCell];
                        }

                    }
                }
            } else {
                for (long k=0; k<secondArg.lDim; k++) {
                    long m = secondArg.theIndex[k];
                    if (m!=-1) {
                        long i = m/secondArg.vDim;
                        long j = m%secondArg.vDim;
                        _Parameter c = secondArg.theData[k];
                        _Parameter *stData = storage.theData+j,
                                    *secData = theData+i,
                                     *stopper = theData+lDim;
                        for (; secData<stopper; secData+=vDim, stData+=secondArg.vDim) {
                            *stData += c**secData;
                        }
                    }
                }
            }
        } else { // polynomial entries
            for (long k=0; k<secondArg.lDim; k++) {

                if (secondArg.IsNonEmpty(k)) {
                    long i = secondArg.theIndex[k]/secondArg.vDim;
                    long j = secondArg.theIndex[k]%secondArg.vDim;
                    _MathObject* p = secondArg.GetMatrixObject(k);
                    for (long l=i, m=j; l<lDim; l+=vDim,m+=secondArg.vDim) {
                        tempP = GetMatrixObject (l);
                        if (!tempP) {
                            continue;
                        }
                        _MathObject* temp = p->Mult (tempP);
                        tempP = storage.GetMatrixObject(m);
                        if (tempP) {
                            storage.StoreObject (m, tempP->Add (temp));
                        } else {
                            storage.StoreObject (m, temp, true);
                        }
                        DeleteObject (temp);
                    }
                }
            }

        }
    } else
        //sparse by sparse
    {
        long *indexTable,
             *indexTable2,
             *indexVector,
             //indexTableDim = secondArg.hDim*(secondArg.vDim+1),
             indexTableDim = secondArg.hDim*secondArg.vDim,
             t,
             //dd = secondArg.vDim+1 ;
             dd = secondArg.vDim;
        indexTable  = (long*)MatrixMemAllocate( sizeof(long)*indexTableDim);
        indexTable2 = (long*)MatrixMemAllocate( sizeof(long)*indexTableDim);
        indexVector = (long*)MatrixMemAllocate( sizeof(long)*secondArg.hDim);

        if (!(indexTable&&indexTable2&&indexVector)) {
            warnError (-108);
        }

        memset (indexTable,0,indexTableDim*sizeof(long));
        memset (indexTable2,0,indexTableDim*sizeof(long));
        memset (indexVector,0,secondArg.hDim*sizeof(long));
        if (storageType == 1)
            // numeric
        {
            for (long i=0; i<secondArg.lDim; i++) {
                if ((t=secondArg.theIndex[i])!=-1) {
                    long k = t/secondArg.vDim;
                    long j = k*dd+(indexVector[k]++);
                    indexTable [j] = t%secondArg.vDim;
                    indexTable2[j] = i;
                }
            }
            for (long k=0; k<lDim; k++) {
                if ((t=theIndex[k])!=-1) {
                    long i = t/vDim;
                    long j = t%vDim;
                    _Parameter c = theData[k];
                    long n = j*dd;
                    long m = i*secondArg.vDim;
                    for (long l=n; l<n+indexVector[j]; l++) {
                        storage.theData[m+indexTable[l]]+= c*secondArg.theData[indexTable2[l]];
                    }
                }
            }
        } else { // polynomial entries
            for (long i=0; i<secondArg.lDim; i++) {
                t=secondArg.theIndex[i];
                if (IsNonEmpty(i)) {
                    long k = t/secondArg.vDim;
                    long j = k*dd+(indexVector[k]++);
                    indexTable [j] = t%secondArg.vDim;
                    indexTable2[j] = i;
                }
            }
            for (long k=0; k<lDim; k++) {
                if (IsNonEmpty(k)) {
                    long i = theIndex[k]/vDim;
                    long j = theIndex[k]%vDim;
                    _MathObject* p = GetMatrixObject(k);
                    long n = j*dd;
                    long m = i*secondArg.vDim;
                    for (long l=n; l<n+indexVector[j]; l++) {
                        _MathObject* temp = p->Mult (secondArg.GetMatrixObject (indexTable2[l]));
                        tempP = storage.GetMatrixObject(m+indexTable[l]%secondArg.vDim);
                        if (tempP) {
                            storage.StoreObject (m+indexTable[l]%secondArg.vDim, tempP->Add (temp));
                        } else {
                            storage.StoreObject (m+indexTable[l]%secondArg.vDim, temp, true);
                        }

                        DeleteObject (temp);
                    }
                }
            }
        }

        MatrixMemFree( indexTable);
        MatrixMemFree( indexTable2);
        MatrixMemFree( indexVector);

    }


}

//_____________________________________________________________________________________________

long    _Matrix::HashBack  (long logicalIndex)
// returns element's matrix index in the form vDim*(i-1)+j, where i, j are the matrix coordinates
// given a buffer index
{
    return theIndex?theIndex [logicalIndex]:logicalIndex;
}

//_____________________________________________________________________________________________

_Parameter  _Matrix::MaxElement  (char runMode, long* indexStore)
// returns matrix's largest abs value element
{
    if (storageType == 1) {
        _Parameter max  = 0.0,
                   temp;

        bool doAbsValue = runMode != 1 && runMode != 3,
             doMaxElement = runMode == 0 || runMode == 3;

        if (doMaxElement) {
            max = -A_LARGE_NUMBER;
        }

        if (theIndex) {
            for (long i = 0; i<lDim; i++) {
                long k = theIndex[i];
                if  (k != -1) {
                    temp = theData[i];
                    if (doAbsValue && temp<0.0) {
                        temp = -temp;
                    }

                    if (doMaxElement) {
                        if (temp>max) {
                            max = temp;
                            if (indexStore) {
                                *indexStore = k;
                            }
                        }
                    } else {
                        max += temp;
                    }
                }
            }
            return max;
        } else {
            for (long i = 0; i<lDim; i++) {
                temp = theData[i];
                if (doAbsValue && temp<0.0) {
                    temp = -temp;
                }

                if (doMaxElement) {
                    if (temp>max) {
                        max = temp;
                        if (indexStore) {
                            *indexStore = i;
                        }
                    }
                } else {
                    max += temp;
                }
            }
            return max;
        }
    }
    if (runMode) {
        return 0;
    }

    return 10.0;
}

//_____________________________________________________________________________________________

void    _Matrix::RowAndColumnMax  (_Parameter& r, _Parameter &c, _Parameter * cache)

// returns the maximum row sum / column sum
// the cache must be big enough to hold hDim + vDim
// leave as nil to allocate cache run time

{
    r = c = 10.;

    if (storageType == 1) { // numeric matrix
        _Parameter  *maxScratch = cache;
        r = c = 0.;

        if (maxScratch == nil) {
            checkPointer (maxScratch = (_Parameter*)calloc (hDim+vDim,sizeof(_Parameter)));
        } else
            for (long k = 0; k < hDim + vDim; k++) {
                maxScratch[k] = .0;
            }

        _Parameter * rowMax = maxScratch,
                     * colMax = maxScratch + hDim;

        if (theIndex)
            // sparse matrix
            for (long i = 0; i<lDim; i++) {
                long k = theIndex[i];
                if  (k!=-1) {
                    _Parameter temp = theData[i];

                    if (temp<0.0) {
                        rowMax[k/vDim] -= temp;
                        colMax[k%vDim] -= temp;
                    } else {
                        rowMax[k/vDim] += temp;
                        colMax[k%vDim] += temp;
                    }
                }
            }
        else
            // dense matrix
            for (long i = 0, k=0; i<hDim; i++) {
                for (long j=0; j<vDim; j++, k++) {
                    _Parameter temp = theData[k];
                    if (temp<0.0) {
                        rowMax[i] -= temp;
                        colMax[j] -= temp;
                    } else {
                        rowMax[i] += temp;
                        colMax[j] += temp;
                    }
                }
            }

        for (long i=0; i<hDim; i++) if (rowMax[i]>r)    {
                r = rowMax [i];
            }
        for (long j=0; j<vDim; j++) if (colMax[j]>c)    {
                c = colMax [j];
            }

        if (!cache) {
            free(maxScratch);
        }
    }
}

//_____________________________________________________________________________________________

bool    _Matrix::IsMaxElement  (_Parameter bench)
// returns matrix's largest abs value element
{
    if (storageType == 1) {
        _Parameter t,
                   mBench = -bench;
        for (long i = 0; i<lDim; i++) {
            t = theData[i];
            if ((t<mBench)||(t>bench)) {
                return true;
            }
        }
        return false;
    } else if (storageType == 0) {
        _Polynomial ** pData = (_Polynomial **)theData;
        for (long i = 0; i<lDim; i++, pData++) {
            if ((*pData)->IsMaxElement(bench)) {
                return true;
            }
        }
        return false;
    }
    return true;
}

//_____________________________________________________________________________________________

_Parameter  _Matrix::MaxRelError  (_Matrix& compMx)
// returns matrix's largest abs value element
{
    if (storageType == 1) {
        _Parameter max = 0, temp;
        for (long i = 0; i<lDim; i++) {
            temp = theData[i]/compMx.theData[i];
            if (temp<0.0) {
                temp*=-1.0;
            }
            if (temp>max) {
                max = temp;
            }
        }
        return max;
    }
    return 10.0;
}

//_____________________________________________________________________________________________

bool    _Matrix::IsAVector  (char type)
{
    if (GetHDim() == 1) {
        return type != HY_MATRIX_COLUMN_VECTOR;
    }
    if (GetVDim() == 1) {
        return (type != HY_MATRIX_ROW_VECTOR);
    }
    return false;
}

//_____________________________________________________________________________________________

_Parameter  _Matrix::MinElement  (char doAbsValue, long* storeIndex)
// returns matrix's smalles non-zero abs value element
{
    if (storageType == 1) {
        _Parameter min = DBL_MAX;

        if (theIndex)
            for (long i = 0; i<lDim; i++) {
                if (theIndex[i] < 0) {
                    continue;
                }

                _Parameter temp = theData[i];

                if (temp < 0.0 && doAbsValue) {
                    temp = -temp;
                }

                if (temp<min) {
                    if (storeIndex) {
                        *storeIndex = theIndex[i];
                    }
                    min = temp;
                }

            }
        else
            for (long i = 0; i<lDim; i++) {
                _Parameter temp = theData[i];

                if (temp < 0.0 && doAbsValue) {
                    temp = -temp;
                }

                if (temp<min) {
                    if (storeIndex) {
                        *storeIndex = i;
                    }
                    min = temp;
                }
            }

        return min;
    } else {
        return 1.0;
    }
}
//_____________________________________________________________________________________________
void    _Matrix::Transpose (void)
// transpose a matrix
{
    if (storageType == 1) {
        if (hDim == vDim) { // do an in place swap
            if (!theIndex) { // non-sparse
                for (long i = 0; i<hDim; i++)
                    for (long j = i+1; j<vDim; j++) {
                        _Parameter z      = theData[i*vDim+j];
                        theData[i*vDim+j] = theData[j*vDim+i];
                        theData[j*vDim+i] = z;
                    }
            } else { // sparse
                for (long i = 0; i<lDim; i++) {
                    long p = theIndex[i];
                    if (p!=-1) {
                        long k      = p/vDim;
                        long l      = p%vDim;

                        if (l!=k) { // off - diag
                            p            = Hash (l,k);
                            _Parameter z = theData[i];
                            if (p>=0) {
                                theData[i] = theData[p];
                                theData[p] = z;
                            } else {
                                theIndex[i]=-1;
                                (*this)[l*vDim+k]=z;
                            }
                        }
                    }
                }
            }
        } else {
            _Matrix result (vDim, hDim, bool(theIndex), true);
            if (!theIndex) { // dense
                for (long i = 0; i<hDim; i++)
                    for (long j = 0; j<vDim; j++) {
                        result.theData[j*hDim+i]=theData[i*vDim+j];
                    }
            } else {
                for (long i = 0; i<lDim; i++)
                    if (IsNonEmpty(i)) {
                        long r = theIndex[i]/vDim,
                             c = theIndex[i]%vDim;

                        result[c*hDim+r]=theData[i];
                    }
            }
            *this = result;
        }
    } else { // polynomial entries
        Ptr z;
        if (hDim == vDim) {
            if (!theIndex) { // non-sparse
                for (long i = 0; i<hDim; i++)
                    for (long j = i+1; j<vDim; j++) {
                        if (storageType==2) {
                            z = (Ptr)GetFormula(i,j);
                        } else {
                            z = (Ptr)GetMatrixObject(i*vDim+j);
                        }

                        if (storageType==2) {
                            ((_Formula**)theData)[i*vDim+j] = GetFormula(j,i);
                        } else {
                            ((_PMathObj*)theData)[i*vDim+j] = GetMatrixObject(j*vDim+i);
                        }

                        ((Ptr*)theData)[j*vDim+i] = z;
                    }
            } else { // sparse
                long i,k,l,p;
                for (i = 0; i<lDim; i++) {
                    if (IsNonEmpty(i)) {
                        k = theIndex[i]/vDim;
                        l = theIndex[i]%vDim;
                        if (l!=k) {
                            p = Hash (l,k);

                            if (storageType==2) {
                                z = (Ptr)GetFormula(k,l);
                            } else {
                                z = (Ptr)GetMatrixObject(i);
                            }

                            if (p>=0) {
                                if (storageType==2) {
                                    ((_Formula**)theData)[i]    = GetFormula(l,k);
                                } else {
                                    ((_MathObject**)theData)[i] = GetMatrixObject(p);
                                }

                                ((Ptr*)theData)[p] = z;
                            } else {
                                theIndex[i]=-1;
                                if (storageType==2) {
                                    StoreFormula(l,k,*(_Formula*)z,false,false);
                                } else {
                                    StoreObject(l*vDim+k,(_PMathObj)z);
                                }
                            }
                        }
                    }
                }
            }
        } else {

            _Matrix result;
            CreateMatrix (&result,vDim, hDim, bool(theIndex),false,storageType==2);
            result.storageType = storageType;
            if (!theIndex) {
                for (long i = 0; i<hDim; i++)
                    for (long j = 0; j<vDim; j++) {
                        if (storageType == 2) {
                            result.StoreFormula(j,i,*GetFormula(i,j),true,false);
                        } else {
                            z =   (Ptr)GetMatrixObject (i*vDim+j);
                            result.StoreObject(j*hDim+i,(_PMathObj)z);
                            ((_PMathObj)z)->nInstances++;
                        }
                    }
            } else {
                long r,c;
                for (long i = 0; i<lDim; i++)
                    if (IsNonEmpty(i)) {
                        r = theIndex[i]/vDim;
                        c = theIndex[i]%vDim;
                        if (storageType == 2) {
                            result.StoreFormula(c,r,*GetFormula(r,c),true,false);
                        } else {
                            z =   (Ptr)GetMatrixObject (i);
                            result.StoreObject(c*hDim+r,(_PMathObj)z);
                            ((_PMathObj)z)->nInstances++;
                        }
                    }
            }
            Swap (result);
            //*this = result;
        }
    }
}

//_____________________________________________________________________________________________

void    _Matrix::CompressSparseMatrix (bool transpose, _Parameter * stash)
{
    if (theIndex) {
        _SimpleList sortedIndex  ((unsigned long)lDim)
        ,sortedIndex3 ((unsigned long)lDim)
        ,sortedIndex2
        ;


        long blockChunk = 32,
             blockShift = hDim / blockChunk + 1,
             max        = 0;


        for (long i2=0; i2<lDim; i2++) {
            long k = theIndex[i2];
            if  (k!=-1) {
                long r = transpose?(k/vDim):(k%vDim),
                     c = transpose?(k%vDim):(k/vDim),
                     r2 = c / blockChunk * blockShift + r / blockChunk,
                     r3 = r2 * lDim + r * vDim + c;

                sortedIndex  << (c*vDim + r);
                sortedIndex3 << r3;
                stash[sortedIndex.lLength-1] = theData[i2];
                if (r3 > max) {
                    max = r3;
                }
            }
        }

        if (max > (lDim<<4)) {
            sortedIndex2. Populate(sortedIndex.lLength,0,1);
            SortLists(&sortedIndex3,&sortedIndex2);
        } else {
            DeleteObject (sortedIndex3.CountingSort(-1, &sortedIndex2));
        }

        for (long i=0; i<sortedIndex.lLength; i++) {
            theIndex[i] = sortedIndex.lData[sortedIndex2.lData[i]];
            theData[i]  = stash[sortedIndex2.lData[i]];
        }

        lDim = sortedIndex.lLength;
    }
}


//_____________________________________________________________________________________________

_Matrix*    _Matrix::Exponentiate (void)
{
    // find the maximal elements of the matrix
    long i,
         power2 = 0;

#ifndef _OPENMP
    matrixExpCount++;
#endif

    _Parameter max     = 1.0,
               *stash  = new _Parameter[hDim*(1+vDim)];

    if (storageType) {
        _Parameter t;
        RowAndColumnMax (max, t, stash);
        max *= t;
        if (max > .1) {
            max             = sqrt (10.*max);
            power2          = (long)((log (max)/log ((_Parameter)2.0)))+1;
            max             = exp (power2 * log ((_Parameter)2.0));
            (*this)         *= 1.0/max;
        } else {
            power2 = 0;
        }

        if (theIndex)
            // transpose sparse matrix
        {
            CompressSparseMatrix (true,stash);
        }

    } else {
        max = 1.;
    }

    _Matrix *result = new _Matrix(hDim, vDim , !storageType, storageType),
    temp    (*this);

    checkPointer (result);
    // put ones on the diagonal

    if (storageType) {
        for (i=0; i<result->lDim; i+=vDim+1) {
            result->theData[i]=1.0;
        }
    } else {
        _Polynomial one (1.0);
        for (i=0; i<(*result).hDim*(*result).vDim; i+=vDim+1) {
            (*result).StoreObject(i,&one,true);
        }
    }

    if (max == 0.0) {
        delete [] stash;
        return result;
    }

    (*result) += (*this);

    i = 2;

    if (precisionArg||!storageType) {
        if (storageType)
            for (; i<=precisionArg; i++) {
                temp      *= (*this);
                temp      *= 1.0/i;
                (*result) += temp;
            }
        else {
            while (temp.IsMaxElement (polynomialExpPrecision)) {
                if (i>maxPolynomialExpIterates) {
                    break;
                }
                temp        *= (*this);
                temp        *= 1.0/i;
                (*result)   += temp;
                i++;
            }
            if (i>maxPolynomialExpIterates) {
                _String   wM ("Polynomial Matrix Exponential Failed to achieve accuracy POLYNOMIAL_EXP_PRECISION in under MAX_POLYNOMIAL_EXP_ITERATES. Either decrease the precision, or increase the maximum number of iterates.");
                ReportWarning (wM);
            }
        }
    } else {
        _Parameter tMax = MAX(MinElement()*sqrt ((_Parameter)hDim),truncPrecision);

        i=2;

        _Matrix tempS (hDim, vDim, false, temp.storageType);
        do {
            temp.MultbyS        (*this,theIndex!=nil, &tempS, stash);
            temp      *= 1.0/i;
            (*result) += temp;
            i         ++;
#ifndef _OPENMP
            taylorTermsCount++;
#endif
        } while (temp.IsMaxElement(tMax*truncPrecision*i));

        // use Pade (4,4) here

        /*_Matrix temp (*this), top (result) , bottom (result);
        temp *= .5;
        top+=temp;
        bottom-=temp;
        temp *= *this;
        temp *= 3.0/14.0;
        top+=temp;
        bottom+=temp;
        temp *= *this;
        temp *= 1.0/9.0;
        top+=temp;
        bottom-=temp;
        temp *= 1.0/20.0;
        top+=temp;
        bottom+=temp;
        _Matrix* inv = (_Matrix*)bottom.Inverse();
        top *= *inv;
        DeleteObject (inv);
        result = top;*/
    }

    if (power2) {
        (*this)*=max;
    }

    if (theIndex)
        // transpose back
    {
        for (i=0; i<lDim; i++) {
            long k = theIndex[i];
            if  (k!=-1) {
                theIndex[i] = (k%vDim)*vDim + k/vDim;
            }
        }
        result->Transpose();
    }


    for (long s = 0; s<power2; s++) {
#ifndef _OPENMP
        squaringsCount++;
#endif
        result->Sqr(stash);
    }
    delete [] stash;

    return result;
}

//_____________________________________________________________________________________________

long    _Matrix::Hash  (long i, long j)
// returns element's position in the buffer (-1 if not found)
{
    if (!bufferPerRow)
        // first call, need to set up hashing parameters
    {
        bufferPerRow = lDim/hDim;
        overflowBuffer = hDim*storageIncrement/100;
        if (!(bufferPerRow = (lDim-overflowBuffer)/hDim)) {
            bufferPerRow = 1;
        }
        overflowBuffer = lDim-bufferPerRow*hDim;
        allocationBlock = hDim*vDim*storageIncrement/100+1;
    }

    if (!theIndex) {
        return i*vDim+j;
    }
    // ordinary matrix

    long elementIndex = i*vDim+j, k, l, m=i*bufferPerRow,n,p;

    for (k = 0; k<lDim/allocationBlock; k++,m+=allocationBlock) {
        for (l=m; l<m+bufferPerRow; l++) {
            p = theIndex[l];
            if (p!=elementIndex) {
                if (p==-1) {
                    return -l-2;
                }
            } else {
                return l;
            }
        }
        n = (k+1)*allocationBlock-1;
        for (l = n; l>n-overflowBuffer; l--) {
            p = theIndex[l];
            if (p!=elementIndex) {
                if (p==-1) {
                    return -l-2;
                }
            } else {
                return l;
            }
        }
    }
    return -1;
}

//_____________________________________________________________________________________________
_Parameter      _Matrix::operator () (long i, long j)
{
    long lIndex = Hash (i,j);
    if (lIndex<0) {
        return ZEROOBJECT;
    } else {
        return theData[lIndex];
    }
}

//_____________________________________________________________________________________________
_Matrix*        _Matrix::ExtractElementsByEnumeration (_SimpleList*h, _SimpleList*v, bool column) // extract by row
{
    if (storageType && h->lLength == v->lLength && h->lLength > 0) {
        _Matrix * result = new _Matrix (column?h->lLength:1,column?1:h->lLength,false,true);
        checkPointer (result);

        if (storageType == 2) // formulae
            for (long k=0; k<h->lLength; k++) {
                result->StoreFormula(column?k:0,column?0:k,*GetFormula(h->lData[k],v->lData[k]));
            }
        else
            for (long k=0; k<h->lLength; k++) {
                result->theData[k] = (*this)(h->lData[k],v->lData[k]);
            }

        return result;
    }
    return new _Matrix;
}



//_____________________________________________________________________________________________
_PMathObj _Matrix::MAccess (_PMathObj p, _PMathObj p2)
{
    if (!p) {
        warnError(-106);
        return new _Constant (0.0);
    }

    if (hDim == 0 || vDim == 0) {
        return new _Constant (0.0);
    }

    if (p->ObjectClass() == MATRIX) {
        if (p2 == nil) {
            _Matrix * nn = (_Matrix*)p;
            if (nn->storageType == 1)
                if (nn->hDim == hDim && nn->vDim == vDim) {
                    _SimpleList hL,
                                vL;

                    for (long r=0; r<hDim; r++)
                        for (long c=0; c<vDim; c++)
                            if ((*nn)(r,c) > 0.0) {
                                hL << r;
                                vL << c;
                            }

                    return ExtractElementsByEnumeration (&hL,&vL);
                } else {
                    if (nn->hDim > 0 && nn->vDim == 1) { // extract by row
                        _SimpleList hL;

                        for (long r=0; r<nn->hDim; r++) {
                            long v = (*nn)(r,0);
                            if (v>=0 && v<hDim) {
                                hL<<v;
                            }
                        }

                        if (hL.lLength) {
                            _Matrix * result = new _Matrix (hL.lLength,vDim,false,true);
                            checkPointer (result);
                            long k = 0;
                            for (long r=0; r<hL.lLength; r++) {
                                long ri = hL.lData[r];
                                for (long c=0; c<vDim; c++,k++) {
                                    result->theData[k] = (*this)(ri,c);
                                }
                            }
                            return result;
                        }

                        return new _Matrix;
                    } else if (nn->vDim > 0 && nn->hDim == 1) { // extract by column
                        _SimpleList hL;

                        for (long r=0; r<nn->vDim; r++) {
                            long v = (*nn)(0,r);
                            if (v>=0 && v<vDim) {
                                hL<<v;
                            }
                        }

                        if (hL.lLength) {
                            _Matrix * result = new _Matrix (hDim,hL.lLength,false,true);
                            checkPointer (result);
                            long k = 0;
                            for (long c=0; c<hDim; c++)
                                for (long r=0; r<hL.lLength; r++,k++) {
                                    result->theData[k] = (*this)(c,hL.lData[r]);
                                }
                            return result;
                        }

                        return new _Matrix;
                    }
                }

            ReportWarning ("Incorrect dimensions or matrix type (must be numeric) for an indexing matrix in call to []");
        } else {
            if (p2->ObjectClass() == MATRIX) {
                _Matrix * nn =  (_Matrix*)((_Matrix*)p)->ComputeNumeric();
                _Matrix * nn2 = (_Matrix*)((_Matrix*)p2)->ComputeNumeric();

                if (nn->hDim == 1 && nn->vDim == 2 && nn->storageType == 1 && nn2->hDim == 1 && nn2->vDim == 2 && nn2->storageType == 1) {
                    long left   = (*nn)(0,0),
                         top    = (*nn)(0,1),
                         bottom = (*nn2)(0,1),
                         right  = (*nn2)(0,0);

                    if (left >= 0 && left < hDim && right >= 0 && right < hDim && left <=right &&
                            top >= 0 && top < vDim && bottom >=0 && bottom < vDim && top <= bottom) {
                        _SimpleList hL,
                                    vL;

                        for (long r=left; r<=right; r++)
                            for (long c=top; c<=bottom; c++) {
                                hL << r;
                                vL << c;
                            }

                        _Matrix * subM = ExtractElementsByEnumeration (&hL,&vL);
                        subM->hDim = right-left+1;
                        subM->vDim = bottom-top+1;

                        return subM;
                    }
                }
                ReportWarning ("Incorrect dimensions or matrix type (must be numeric 2x1 matrices) for an rectangular extract in call to []");
            }

        }
        return new _Constant (0.0);
    } else {
        if (p->ObjectClass() == STRING) {


            _String aFormulaString = *((_FString*)p)->theString;
            _Formula f (aFormulaString);

            if (!f.IsEmpty()) {
                /* check formula validity */

                _String cell_value ("_MATRIX_ELEMENT_VALUE_"),
                        cell_row   ("_MATRIX_ELEMENT_ROW_"),
                        cell_column("_MATRIX_ELEMENT_COLUMN_");

                _Variable * cv = CheckReceptacle(&cell_value, empty, false),
                            * cr = CheckReceptacle(&cell_row, empty, false),
                              * cc = CheckReceptacle(&cell_column, empty, false);

                cv->CheckAndSet (0.0);
                cr->CheckAndSet (0.0);
                cc->CheckAndSet (0.0);

                f.Compute();
                if (terminateExecution) {
                    return new _Matrix ();
                } else {

                    _Formula * conditionalCheck = nil;

                    if (p2 && p2->ObjectClass() == STRING) {
                        conditionalCheck = new _Formula (*((_FString*)p2)->theString);
                        if (conditionalCheck->IsEmpty()) {
                            delete conditionalCheck;
                            conditionalCheck = nil;
                        }

                        conditionalCheck->Compute();
                        if (terminateExecution) {
                            delete conditionalCheck;
                            return new _Matrix ();
                        }
                    }

                    _Matrix   * retMatrix = new _Matrix (hDim,vDim,false,true);

                    long          stackDepth = 0;
                    _SimpleList   vIndex;

                    if (f.AmISimple (stackDepth,vIndex) && (!conditionalCheck || conditionalCheck->AmISimple(stackDepth,vIndex))) {
                        _SimpleFormulaDatum * stack     = new _SimpleFormulaDatum [stackDepth+1],
                        * varValues = new _SimpleFormulaDatum [vIndex.lLength];

                        bool                constantValue = false;
                        _Parameter          constantV     = f.Compute()->Value();

                        if (f.IsConstant()) {
                            constantValue = true;
                            constantV     = f.Compute()->Value();
                        } else {
                            f.ConvertToSimple (vIndex);
                        }

                        if (conditionalCheck) {
                            conditionalCheck->ConvertToSimple(vIndex);
                        }

                        if (constantValue && !conditionalCheck) {
                            for (long r=0; r<hDim; r++)
                                for (long c=0; c<vDim; c++) {
                                    retMatrix->Store (r,c,constantV);
                                }
                        } else {

                            long rid []= {cr->GetAVariable(),cc->GetAVariable(),cv->GetAVariable()};

                            for (long k=0; k<3; k++) {
                                rid[k] = vIndex.Find(rid[k]);
                            }

                            PopulateArraysForASimpleFormula(vIndex, varValues);

                            for (long r=0; r<hDim; r++) {

                                if (rid[0]>=0) {
                                    varValues[rid[0]].value = r;
                                }

                                for (long c=0; c<vDim; c++) {
                                    if (rid[1]>=0) {
                                        varValues[rid[1]].value = c;
                                    }

                                    if (rid[2]>=0) {
                                        varValues[rid[2]].value = (*this)(r,c);
                                    }

                                    if (conditionalCheck && CheckEqual(conditionalCheck->ComputeSimple(stack,varValues),0.0)) {
                                        if (rid[2]>=0) {
                                            retMatrix->Store (r,c,varValues[rid[2]].value);
                                        } else {
                                            retMatrix->Store (r,c, (*this)(r,c));
                                        }
                                        continue;
                                    }

                                    if (constantValue) {
                                        retMatrix->Store (r,c,constantV);
                                    } else {
                                        retMatrix->Store (r,c,f.ComputeSimple(stack,varValues));
                                    }
                                }
                            }

                            f.ConvertFromSimple (vIndex);
                        }
                        if (conditionalCheck) {
                            conditionalCheck->ConvertFromSimple(vIndex);
                        }

                        delete  [] stack;
                        delete  [] varValues;
                    } else {
                        for (long r=0; r<hDim; r++) {
                            cr->CheckAndSet (r);
                            for (long c=0; c<vDim; c++) {
                                cc->CheckAndSet (c);
                                cv->CheckAndSet ((*this)(r,c));
                                _PMathObj fv;

                                if (conditionalCheck) {
                                    fv = conditionalCheck->Compute();
                                    if (fv->ObjectClass() == NUMBER)
                                        if (CheckEqual (fv->Value(), 0.0)) {
                                            retMatrix->Store (r,c,cv->Value());
                                            continue;
                                        }
                                }

                                fv = f.Compute();
                                if (fv->ObjectClass()==NUMBER) {
                                    retMatrix->Store (r,c,fv->Value());
                                }
                            }
                        }
                    }
                    retMatrix->AmISparse();
                    if (conditionalCheck) {
                        delete conditionalCheck;
                    }
                    return retMatrix;
                }
            }
            ReportWarning (_String("Invalid formula expression for element-wise matrix operations: ") & *((_FString*)p)->theString);
            return new _Matrix;
        }
    }

    long    ind1 = p->Value(),
            ind2 = -1;

    if (p2) {
        ind2 = p2->Value();
        // handle the row/column access operations here i.e. [R][-1] or [-1][R]

        if (ind1 == -1 && ind2 >=0 && ind2 <vDim) { // valid column access
            _SimpleList hL (hDim,0,1),
                        vL (hDim,ind2,0);
            return ExtractElementsByEnumeration (&hL,&vL,true);
        }

        if (ind2 == -1 && ind1 >=0 && ind1 <hDim) { // valid row access
            _SimpleList hL (vDim,ind1,0),
                        vL (vDim,0,1);
            return ExtractElementsByEnumeration (&hL,&vL);
        }
    }

    if (hDim == 1) {
        if (ind2<0) {
            ind2 = ind1;
        }
        ind1=0;
    }

    if (vDim == 1) {
        ind2 = 0;
    }

    if (ind2<0) { // allow direct vectorlike indexing, i.e m[21] = m[3][3] (if the dim is *x6)
        ind2  = ind1%vDim;
        ind1 /=vDim;
    }

    if (ind1<0 || ind1>=hDim || ind2>=vDim) {
        MatrixIndexError     (ind1,ind2,hDim,vDim);
        return new _Constant (0.0);
    }

    if (ind2>=0) { // element access
        if (storageType == 2) { // formulas
            if (!theIndex) {
                _Formula * entryFla = (((_Formula**)theData)[ind1*vDim+ind2]);
                if (entryFla) {
                    return (_PMathObj)entryFla->Compute()->makeDynamic();
                } else {
                    return new _Constant (0.0);
                }
            } else {
                long p = Hash (ind1, ind2);
                if (p<0) {
                    return new _Constant (0.0);
                } else {
                    return (_PMathObj)(((_Formula**)theData)[p])->Compute()->makeDynamic();
                }
            }
        } else {
            if (storageType == 1) {
                if (theIndex) {
                    return new _Constant ((*this)(ind1,ind2));
                } else {
                    return new _Constant (theData[ind1*vDim+ind2]);
                }

            } else {
                _MathObject* cell;
                if (!theIndex) {
                    cell = (_MathObject*)GetMatrixObject (ind1*vDim+ind2)->makeDynamic();
                } else {
                    long p = Hash (ind1, ind2);
                    if (p<0) {
                        cell = new _Constant (0.0);
                    } else {
                        cell = (_MathObject*)GetMatrixObject (p)->makeDynamic();
                    }
                }
                return cell;
            }
        }
    }

    return new _Constant (0.0);
}

//_____________________________________________________________________________________________
_Formula* _Matrix::GetFormula (long ind1, long ind2)
{

    if (hDim == 1) {
        if (ind2<0) {
            ind2 = ind1;
        }
        ind1=0;
    }

    if (vDim == 1) {
        ind2 = 0;
    }

    if (ind2<0) {
        ind2 = ind1%vDim;
        ind1/=vDim;
    }

    if ( ind1<0 || ind1>=hDim || ind2>=vDim) {
        MatrixIndexError (ind1,ind2,hDim,vDim);
        return nil;
    }


    if (ind2>=0) { // element access
        if (storageType == 2) { // formulas
            if (!theIndex) {
                return (((_Formula**)theData)[ind1*vDim+ind2]);
            } else {
                long p = Hash (ind1, ind2);
                if (p<0) {
                    return nil;
                } else {
                    return (((_Formula**)theData)[p]);
                }
            }
        }
    }

    return nil;
}

//_____________________________________________________________________________________________
_PMathObj _Matrix::MCoord (_PMathObj p, _PMathObj p2)
{
    long ind1 = -1,
         ind2 = -1;

    if (!p) {
        warnError( -106);
    }

    ind1 = p->Value();
    if (p2) {
        ind2 = p2->Value();
    }


    if (hDim == 1) {
        if (ind2<0) {
            ind2 = ind1;
        }
        ind1=0;
    }

    if (vDim == 1) {
        ind2 = 0;
    }

    if (ind2<0) { // allow direct vectorlike indexing, i.e m[21] = m[3][3] (if the dim is *x6)
        ind2 = ind1%vDim;
    }

    _Matrix * res = new _Matrix (1,2,false,true);
    if (res) {
        res->theData[0]=ind1;
        res->theData[1]=ind2;
    } else {
        checkPointer (res);
    }
    return res;

}

//_____________________________________________________________________________________________
bool _Matrix::MResolve (_PMathObj p, _PMathObj p2, long& ind1, long& ind2)
{
    ind1 = -1;
    ind2 = -1;

    if (!p) {
        warnError(-106);
        return false;
    }

    ind1 = p->Value();
    if (p2) {
        ind2 = p2->Value();
    }

    return CheckCoordinates (ind1,ind2);
}

//_____________________________________________________________________________________________

bool _Matrix::CheckCoordinates (long& ind1, long& ind2)
{
    if (hDim == 1) {
        if (ind2<0) {
            ind2 = ind1;
        }
        ind1=0;
    }

    if (vDim == 1) {
        ind2 = 0;
    }

    if (ind2<0) { // allow direct vectorlike indexing, i.e m[21] = m[3][3] (if the dim is *x6)
        if (vDim > 1) {
            ind2 = ind1%vDim;
            ind1/= vDim;
        } else {
            ind2 = 0;
        }
    }

    if (ind1<0 || ind1>=hDim || ind2>=vDim) {
        MatrixIndexError (ind1,ind2, hDim, vDim);
        return false;
    }
    return true;
}


//_____________________________________________________________________________________________
void _Matrix::MStore (long ind1, long ind2, _Formula& f, long opCode)
{
    if (ind2>=0) { // element storage
        if (storageType == 2) { // formulas
            StoreFormula (ind1,ind2,f);
        } else {
            if (!f.IsAConstant()) {
                Convert2Formulas();
                StoreFormula (ind1,ind2,f);
            } else {
                _PMathObj res = f.Compute();
                _Parameter toStore = res->Value();
                if (opCode == HY_OP_CODE_ADD) {
                    toStore += (*this)(ind1,ind2);
                }
                Store(ind1,ind2,toStore);
            }
        }
    }
}

//_____________________________________________________________________________________________
void _Matrix::MStore (_PMathObj p, _PMathObj p2, _Formula& f, long opCode)
{
    long      ind1, ind2;
    if (MResolve (p,p2, ind1,ind2)) {
        MStore   (ind1,ind2,f, opCode);
    }
}

//_____________________________________________________________________________________________
void _Matrix::MStore (_PMathObj p, _PMathObj p2, _PMathObj poly)
{
    long      ind1, ind2;
    if (MResolve (p,p2, ind1,ind2)) {
        MStore   (ind1,ind2,poly);
    }

}
//_____________________________________________________________________________________________
void _Matrix::MStore (long ind1, long ind2, _PMathObj poly)
{
    if (ind2>=0) { // element storage
        if (storageType == 0) { // formulas
            StoreObject (ind1,ind2,poly,true);
            if (AUTO_PAD_DIAGONAL) {
                UpdateDiag (ind1,ind2,poly);
            }
        } else {
            _Polynomial* pp = (_Polynomial*)poly;
            poly = pp->IsANumber();
            if (!poly) { // just a number
                storageType==1?ConvertNumbers2Poly():ConvertFormulas2Poly();
                StoreObject (ind1,ind2,pp,true);
            } else {
                (*this)[Hash(ind1,ind2)] = poly->Value();
            }
        }
    }
}


//_____________________________________________________________________________________________
_Parameter&     _Matrix::operator [] (long i)
{
    long lIndex = Hash (i/vDim, i%vDim);
    if (lIndex == -1) {
        IncreaseStorage();
        lIndex = Hash (i/vDim, i%vDim);
    }
    if (lIndex<0) {
        theIndex[-lIndex-2] = i;
        return ((_Parameter*)theData)[-lIndex-2];
    } else {
        return ((_Parameter*)theData)[lIndex];
    }
}

//_____________________________________________________________________________________________
void        _Matrix::Store (long i, long j, _Parameter value)
{
    if (storageType!=1) {
        return;
    }

    long lIndex;

    if (theIndex) {
        lIndex = Hash (i, j);

        if (lIndex == -1) {
            IncreaseStorage();
            lIndex = Hash (i, j);
        }
    } else {
        lIndex = i*vDim + j;
    }

    if (lIndex<0) {
        theIndex[-lIndex-2] = i*vDim+j;
        ((_Parameter*)theData)[-lIndex-2] = value;
    } else {
        ((_Parameter*)theData)[lIndex] = value;
    }

}

//_____________________________________________________________________________________________
void        _Matrix::StoreObject (long i, long j, _MathObject* value, bool dup)
{
    if (storageType) {
        return;
    }

    long lIndex = Hash (i, j);
    if (lIndex == -1) {
        IncreaseStorage();
        lIndex = Hash (i, j);
    }

    if (dup) {
        value = (_MathObject*) value->makeDynamic();
    }
    if (lIndex<0) {
        theIndex[-lIndex-2] = i*vDim+j;
        ((_MathObject**)theData)[-lIndex-2] = value;
    } else {
        DeleteObject (GetMatrixObject(lIndex));
        ((_MathObject**)theData)[lIndex] = value;
    }
    if (AUTO_PAD_DIAGONAL) { // correct the diagonal entry
    }

}
//_____________________________________________________________________________________________

void        _Matrix::UpdateDiag  (long i,long j, _MathObject* value)
{
    if (i!=j) {
        _MathObject * diagCell = nil, *newCell;
        if (!theIndex) {
            diagCell = GetMatrixObject(i*hDim+i);
        } else {
            long lIndex = Hash (i,i);
            if (lIndex>=0) {
                diagCell = GetMatrixObject(lIndex);
            }
        }
        if (!diagCell) {
            newCell = value->Minus();
        } else {
            newCell = diagCell->Sub(value);
        }
        StoreObject(i,i,newCell,false);
    }
}
//_____________________________________________________________________________________________
void        _Matrix::StoreObject (long k, _MathObject* value, bool dup)
{
    StoreObject (k/vDim, k%vDim, value, dup);
}

//_____________________________________________________________________________________________
void        _Matrix::StoreFormula (long i, long j, _Formula& f, bool copyF, bool simplify)
{
    if (storageType!=2) {
        return;
    }

    long lIndex = Hash (i, j);
    if (lIndex == -1) {
        IncreaseStorage();
        lIndex = Hash (i, j);
    }

    if (lIndex<0) {
        theIndex[-lIndex-2] = i*vDim+j;
        ((_Formula**)theData)[-lIndex-2] = copyF?(_Formula*)f.makeDynamic():&f;
        if (simplify) {
            ((_Formula**)theData)[-lIndex-2]->SimplifyConstants();
        }
    } else {
        if (copyF && ((_Formula**)theData)[lIndex]!=(_Formula*)ZEROPOINTER) {
            delete ((_Formula**)theData)[lIndex];
        }
        ((_Formula**)theData)[lIndex] = copyF?(_Formula*)f.makeDynamic():&f;
        if (simplify) {
            ((_Formula**)theData)[lIndex]->SimplifyConstants();
        }
    }

    CheckIfSparseEnough();
}

//_____________________________________________________________________________________________


void        _Matrix::Swap (_Matrix& m)
{
    long         *tIndex,
                 t;

    _Parameter   *tData;
    _PMathObj    tObj;
    _CompiledMatrixData *tCmd;

    SWAP(theData,m.theData,tData);
    SWAP(hDim,m.hDim,t);
    SWAP(vDim,m.vDim,t);
    SWAP(lDim,m.lDim,t);
    SWAP(theIndex,m.theIndex,tIndex);
    SWAP(storageType,m.storageType,t);
    SWAP(bufferPerRow,m.bufferPerRow,t);
    SWAP(overflowBuffer,m.overflowBuffer,t);
    SWAP(allocationBlock,m.allocationBlock,t);
    SWAP(theValue,m.theValue,tObj);
    SWAP(cmd,m.cmd,tCmd);
}

//_____________________________________________________________________________________________
/*_Matrix       IterateStrassen (_Matrix& source1, _Matrix& source2)
// square the matrix
{

    if ((source1.storageType!=1)||(source2.storageType!=1)) warnError (-111);
    long iterationDim = source1.hDim, i;

    if (iterationDim>2)
    {
        // create four quadrants of each matrix and temporary storage
        _Matrix a11 (iterationDim/2, iterationDim/2, (bool)false, true),
                a12 (iterationDim/2, iterationDim/2, (bool)false, true),
                a21 (iterationDim/2, iterationDim/2, (bool)false, true),
                a22 (iterationDim/2, iterationDim/2, (bool)false, true),
                b11 (iterationDim/2, iterationDim/2, (bool)false, true),
                b12 (iterationDim/2, iterationDim/2, (bool)false, true),
                b21 (iterationDim/2, iterationDim/2, (bool)false, true),
                b22 (iterationDim/2, iterationDim/2, (bool)false, true),
                s1 (iterationDim/2, iterationDim/2, (bool)false, true),
                s2 (iterationDim/2, iterationDim/2, (bool)false, true),
                s3 (iterationDim/2, iterationDim/2, (bool)false, true),
                s4 (iterationDim/2, iterationDim/2, (bool)false, true),
                s5 (iterationDim/2, iterationDim/2, (bool)false, true),
                s6 (iterationDim/2, iterationDim/2, (bool)false, true),
                t1 (iterationDim/2, iterationDim/2, (bool)false, true),
                t2 (iterationDim/2, iterationDim/2, (bool)false, true),
                c11 (iterationDim/2, iterationDim/2, (bool)false, true),
                c12 (iterationDim/2, iterationDim/2, (bool)false, true),
                c21 (iterationDim/2, iterationDim/2, (bool)false, true),
                c22 (iterationDim/2, iterationDim/2, (bool)false, true);

        // copy data to quadrants

        iterationDim/=2;
        for (i=0;i<iterationDim;i++)
        {
            memcpy ((_Parameter*)a11.theData+i*iterationDim, (_Parameter*)source1.theData+i*source1.hDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)a12.theData+i*iterationDim, (_Parameter*)source1.theData+i*source1.hDim+iterationDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)a21.theData+i*iterationDim, (_Parameter*)source1.theData+(i+iterationDim)*source1.hDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)a22.theData+i*iterationDim, (_Parameter*)source1.theData+(i+iterationDim)*source1.hDim+iterationDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)b11.theData+i*iterationDim, (_Parameter*)source2.theData+i*source2.hDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)b12.theData+i*iterationDim, (_Parameter*)source2.theData+i*source2.hDim+iterationDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)b21.theData+i*iterationDim, (_Parameter*)source2.theData+(i+iterationDim)*source2.hDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)b22.theData+i*iterationDim, (_Parameter*)source2.theData+(i+iterationDim)*source2.hDim+iterationDim, iterationDim*sizeof (_Parameter));
        }

        // compute strassen blocks

        t1=a12;
        t1-=a22;
        t2=b21;
        t2+=b22;
        s1=IterateStrassen (t1,t2);

        t1=a11;
        t1+=a22;
        t2=b11;
        t2+=b22;
        s2=IterateStrassen (t1,t2);

        t1=a11;
        t1-=a21;
        t2=b11;
        t2+=b12;
        s3=IterateStrassen (t1,t2);

        t1=a11;
        t1+=a12;
        s4=IterateStrassen (t1,b22);

        t1=b12;
        t1-=b22;
        s5=IterateStrassen (a11,t1);

        t1=b21;
        t1-=b11;
        s6=IterateStrassen (a22,t1);

        c11 = s1;
        c11 += s2;
        c11 -= s4;
        c11 += s6;

        c12 = s4;
        c12 += s5;

        t1=a21;
        t1+=a22;
        s1=IterateStrassen (t1,b11);

        c21 = s6;
        c21 += s1;

        c22 = s2;
        c22 -= s3;
        c22 += s5;
        c22 -= s1;

        _Matrix result (source1.hDim, source1.hDim, (bool)false, true);

        // copy the results from four resultant quadrants c11..c22


        for (i=0;i<iterationDim;i++)
        {
            memcpy ((_Parameter*)result.theData+i*source1.hDim,(_Parameter*)c11.theData+i*iterationDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)result.theData+i*source1.hDim+iterationDim,(_Parameter*)c12.theData+i*iterationDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)result.theData+(i+iterationDim)*source1.hDim,(_Parameter*)c21.theData+i*iterationDim, iterationDim*sizeof (_Parameter));
            memcpy ((_Parameter*)result.theData+(i+iterationDim)*source1.hDim+iterationDim,(_Parameter*)c22.theData+i*iterationDim, iterationDim*sizeof (_Parameter));
        }

        return result;
    }
    else
        return source1*source2;
}


//_____________________________________________________________________________________________
void        _Matrix::SqrStrassen (void)
// square the matrix
// the matrix is assumed to be a non-pointer matrix
{
    if (hDim!=vDim) return;
    // a non-square matrix

    if (theIndex)
    // sparse matrix - multiply directly at better speed
    {
        _Matrix temp (hDim, vDim, (bool)false, true);
        Multiply (temp, *this);
        Swap(temp);
    }
    else
    {
        // pad the matrix with zeros to make the dimension a power of two;
        long newDim = 2,i;
        while (newDim<hDim) newDim*=2;


        // copies the entries of the original to the padded matrix

        if (newDim != hDim)
        {
            _Matrix paddedMatrix (newDim, newDim, (bool)false, true);
            for (i=0; i<hDim; i++)
                memcpy ((_Parameter*)paddedMatrix.theData+i*newDim, (_Parameter*)theData+i*hDim, hDim*sizeof (_Parameter));

            // call Strassen iteration function


            paddedMatrix .Swap( IterateStrassen (paddedMatrix,paddedMatrix));

            // strip unneeded elements and copy them back into the original;

            for (i=0; i<hDim; i++)
                memcpy ((_Parameter*)theData+i*hDim,(_Parameter*)paddedMatrix.theData+i*newDim, hDim*sizeof (_Parameter));
        }
        else
            Swap( IterateStrassen (*this,*this));

    }
}
*/

//_____________________________________________________________________________________________

void        _Matrix::AplusBx (_Matrix& B, _Parameter x)
{
    _Matrix temp (B);
    temp *= x;
    *this+=temp;
}

//_____________________________________________________________________________________________
void        _Matrix::Sqr (_Parameter* _hprestrict_ stash)
{
    if (hDim!=vDim) {
        return;
    }
    // not a square matrix

    if (theIndex|| storageType!=1 )
        // sparse or non-numeric matrix
    {
        _Matrix temp (hDim, vDim, storageType==0?theIndex!=nil:false, storageType);
        Multiply (temp, *this);
        Swap(temp);
    } else {
        if (hDim==4)
            // special case for nucleotides
        {
            for (long i=0, k = 0; i<16; i+=4) {
                for (long j=0; j<4; j++, k++) {
                    stash[k] = theData[i]   * theData [j]
                               + theData[i+1] * theData [j+4]
                               + theData[i+2] * theData [j+8]
                               + theData[i+3] * theData [j+12];
                }
            }
        } else {
            long loopBound = vDim - vDim % 4;

            // loop interchange rocks!

            _Parameter  * column = stash+lDim;

            for (long j = 0; j < vDim; j++) {
                for (long c = 0; c < vDim; c++) {
                    column[c] = theData[j + c * vDim];
                }

                for (long i = 0; i < lDim; i += vDim) {
                    _Parameter * row    = theData + i,
                                 buffer = 0.0;


#ifndef _SLKP_SSE_VECTORIZATION_
                    long        k = 0;

                    for (; k < loopBound; k+=4)
                        buffer += row[k]   * column [k] +
                                  row[k+1] * column [k+1] +
                                  row[k+2] * column [k+2] +
                                  row[k+3] * column [k+3];

                    for (; k < vDim; k++) {
                        buffer += row[k] * column [k];
                    }
#else
                    for (long k = 0; k < vDim; k++) {
                        buffer += row[k]   * column [k];
                    }

#endif
                    stash[i+j] = buffer;
                }
            }
        }

        for (long s = 0; s < lDim; s++) {
            theData[s] = stash[s];
        }
    }
}
//_____________________________________________________________________________________________
void        _Matrix::AgreeObjects (_Matrix& m)
{
    if (storageType==2)
        if (toPolyOrNot!=0.0) {
            ConvertFormulas2Poly ();
        } else {
            Evaluate(true);
        }

    if (m.storageType==2)
        if (toPolyOrNot!=0.0) {
            m.ConvertFormulas2Poly ();
        } else {
            m.Evaluate(true);
        }

    if (storageType!=m.storageType) {
        if (toPolyOrNot) {
            if (storageType == 1) {
                ConvertNumbers2Poly ();
            } else {
                m.ConvertNumbers2Poly ();
            }
        } else {
            if (storageType == 1) {
                m.Evaluate (true);
            } else {
                Evaluate ();
            }
        }
    }
}
//_____________________________________________________________________________________________
void        _Matrix::ConvertFormulas2Poly (bool force2numbers)
{
    bool conversionFlag = true;
    _MathObject** tempStorage = (_MathObject**)MatrixMemAllocate(sizeof(void*)*lDim);

    long i;

    for (i=0; i<lDim; i++) {
        tempStorage[i]=ZEROPOINTER;
    }


    if (theIndex) { // sparse
        for (i=0; i<lDim; i++) {
            if (IsNonEmpty(i)) {
                _PMathObj polyCell = ((_Formula**)theData)[i]->ConstructPolynomial();
                if (polyCell) { // valid polynomial conversion
                    tempStorage[i] = (_PMathObj)polyCell;
                    polyCell->nInstances++;
                } else {
                    conversionFlag = false;
                    break;
                }
            }
        }
        if (conversionFlag) {
            // check for "*" entries
            for (i=0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    if (((_Formula**)theData)[i]->IsEmpty()) { // "*" entry
                        long r = theIndex[i]/vDim, c = theIndex[i]%vDim;
                        _Polynomial diag;
                        for (long j=0; j<vDim; j++) {
                            if (j==c) {
                                continue;
                            }
                            long h = Hash (r,j);
                            if (h>=0) {
                                _Polynomial * temp = (_Polynomial *)diag.Sub(tempStorage[h]);
                                diag.Duplicate (temp);
                                DeleteObject (temp);
                            }
                        }
                        DeleteObject(tempStorage[i]);
                        tempStorage[i]=(_Polynomial*)diag.makeDynamic();
                    }
                }
            }
        }
    } else {
        for (long i=0; i<lDim; i++) {
            _Formula* f = ((_Formula**)theData)[i];
            if (f->IsEmpty()) {
                continue;
            }
            _PMathObj polyCell = f->ConstructPolynomial();
            if (polyCell) { // valid polynomial conversion
                tempStorage[i] = (_PMathObj)polyCell;
                polyCell->nInstances++;
            } else {
                conversionFlag = false;
                break;
            }
        }
        if (conversionFlag) {
            // check for "*" entries
            for (long i=0; i<lDim; i++) {
                if (((_Formula**)theData)[i]->IsEmpty()) { // "*" entry
                    long r = i/vDim;
                    _Polynomial diag;
                    for (long j=vDim*r; j<vDim*(r+1); j++) {
                        if (j==i) {
                            continue;
                        }
                        _Polynomial * temp = (_Polynomial *)diag.Sub(tempStorage[j]);
                        diag.Duplicate (temp);
                        DeleteObject (temp);
                    }
                    DeleteObject(tempStorage[i]);
                    tempStorage[i]=(_Polynomial*)diag.makeDynamic();
                }
            }
        }
    }

    if (conversionFlag) { // successful conversion
        ClearFormulae();
        MatrixMemFree (theData);
        theData = (_Parameter*) tempStorage;
        storageType = 0;
        if (!theIndex) {
            _Polynomial zero;
            for (i=0; i<lDim; i++)
                if (!GetMatrixObject (i)) {
                    StoreObject (i,&zero,true);
                }
        }
    } else {
        for (long i=0; i<lDim; i++) {
            DeleteObject (tempStorage[i]);
        }
        MatrixMemFree (tempStorage);
        if (force2numbers) {
            Evaluate(true);
        }
    }

}

//_____________________________________________________________________________________________
void        _Matrix::ConvertNumbers2Poly (void)
{
    _MathObject ** tempStorage = (_MathObject**)MatrixMemAllocate (lDim*sizeof (void*));
    if (!theIndex) {
        for (long i=0; i<lDim; i++) {
            tempStorage[i]=new _Polynomial (theData[i]);
        }
    } else {
        for (long i=0; i<lDim; i++)
            if (IsNonEmpty (i)) {
                tempStorage[i]=new _Polynomial (theData[i]);
            } else {
                tempStorage[i] = nil;
            }
    }
    MatrixMemFree (theData);
    theData = (_Parameter*) tempStorage;
    storageType = 0;
}


//_____________________________________________________________________________________________
void        _Matrix::operator += (_Matrix& m)
{
    AgreeObjects (m);
    if ((!m.theIndex) && theIndex) {
        CheckIfSparseEnough(true);
    }
    Add (*this,m);
}

//______________________________________________________________

void    _Matrix::RecursiveIndexSort (long from, long to, _SimpleList* index)
{
    long            middle          = (from+to)/2,
                    bottommove        = 1,
                    topmove        = 1,
                    imiddleV         = index->lData[middle];

    _Parameter middleV = theData[middle];

    if (middle)
        while ((middle-bottommove>=from)&&(theData[middle-bottommove]>=theData[middle])) {
            bottommove++;
        }
    if (from<to)
        while ((middle+topmove<=to)&&(theData[middle+topmove]<=theData[middle])) {
            topmove++;
        }

    for (long i=from; i<middle-bottommove; i++)
        if (theData[i]>=theData[middle]) {
            _Parameter temp = theData[middle-bottommove];
            theData[middle-bottommove] = theData[i];
            theData[i]=temp;
            long       ltemp = index->lData[middle-bottommove];
            index->lData[middle-bottommove] = index->lData[i];
            index->lData[i]=ltemp;
            bottommove++;
            while ((middle-bottommove>=from)&&(theData[middle-bottommove]>=theData[middle])) {
                bottommove++;
            }
        }

    {
        for (long i=middle+topmove+1; i<=to; i++)
            if (theData[i]<=theData[middle]) {
                _Parameter temp = theData[middle+topmove];
                theData[middle+topmove] = theData[i];
                theData[i]=temp;
                long ltemp = index->lData[middle+topmove];
                index->lData[middle+topmove] = index->lData[i];
                index->lData[i]=ltemp;
                topmove++;
                while ((middle+topmove<=to)&&(theData[middle+topmove]<=theData[middle])) {
                    topmove++;
                }
            }
    }

    if (topmove==bottommove) {
        for (long i=1; i<bottommove; i++) {
            _Parameter temp = theData[middle+i];
            theData[middle+i] = theData[middle-i];
            theData[middle-i]=temp;
            long ltemp = index->lData[middle+i];
            index->lData[middle+i] = index->lData[middle-i];
            index->lData[middle-i]=ltemp;
        }
    } else if (topmove>bottommove) {
        long shift = topmove-bottommove;
        for (long i=1; i<bottommove; i++) {
            _Parameter temp = theData[middle+i+shift];
            theData[middle+i+shift] = theData[middle-i];
            theData[middle-i]=temp;
            long ltemp = index->lData[middle+i+shift];
            index->lData[middle+i+shift] = index->lData[middle-i];
            index->lData[middle-i]=ltemp;
        }
        {
            for (long i=0; i<shift; i++) {
                theData[middle+i]       =   theData[middle+i+1];
                index->lData[middle+i]  =   index->lData[middle+i+1];
            }
        }
        middle+=shift;
        theData[middle]=middleV;
        index->lData[middle]=imiddleV;
    } else {
        long shift = bottommove-topmove;
        for (long i=1; i<topmove; i++) {
            _Parameter temp = theData[middle-i-shift];
            theData[middle-i-shift] = theData[middle+i];
            theData[middle+i]=temp;
            long ltemp = index->lData[middle-i-shift];
            index->lData[middle-i-shift] = index->lData[middle+i];
            index->lData[middle+i]=ltemp;
        }
        {
            for (long i=0; i<shift; i++) {
                theData[middle-i]=theData[middle-i-1];
                index->lData[middle-i]=index->lData[middle-i-1];
            }
        }
        middle-=shift;
        theData[middle]=middleV;
        index->lData[middle]=imiddleV;
    }

    if (to>middle+1) {
        RecursiveIndexSort (middle+1,to, index);
    }
    if (from<middle-1) {
        RecursiveIndexSort (from,middle-1, index);
    }
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::SortMatrixOnColumn (_PMathObj mp)
{
    if (storageType!=1) {
        _String    errMsg ("Only numeric matrices can be sorted");
        WarnError  (errMsg);
        return new _Matrix (1,1);
    }

    _SimpleList sortOn;

    if (mp->ObjectClass () != NUMBER || mp->Value() < 0.0 || mp->Value () > GetVDim()-1) {
        bool goodMe = false;
        if (mp->ObjectClass () == MATRIX) {
            _Matrix * sortOnM = (_Matrix*)((_Matrix*)mp)->ComputeNumeric();
            long sortBy = sortOnM->GetHDim()*sortOnM->GetVDim();
            for (long k=0; k<sortBy; k=k+1) {
                long idx = (*sortOnM)[k];
                if (idx>=0) {
                    sortOn << idx;
                }
            }
            goodMe = sortOn.lLength;
        }
        if (!goodMe) {
            _String    errMsg ("Invalid column index to sort the matrix on:");
            errMsg = errMsg & _String((_String*)mp->toStr());
            WarnError  (errMsg);
            return new _Matrix (1,1);
        }
    } else {
        sortOn << mp->Value();
    }

    if (theData == nil) {
        return new _Matrix (1,1);
    }

    _SimpleList             idx (hDim,0,1);
    for (long col2Sort = 0; col2Sort < sortOn.lLength; col2Sort++) {
        long colIdx = sortOn.lData[col2Sort];

        _Matrix theColumn   (1,hDim,false,true);

        if (theIndex)
            for (long k=0; k<hDim; k++) {
                theColumn.theData[k] = (*this)(k, colIdx);
            }
        else
            for (long k=0, j = colIdx; k<hDim; k++, j+=vDim) {
                theColumn.theData[k] = theData[j];
            }

        if (col2Sort == 0) {
            theColumn.RecursiveIndexSort (0, hDim-1, &idx);
        } else {
            long currentFrom = 0,
                 currentTo   = 1;

            bool didSort     = false;

            _SimpleList    revIdx (hDim,0,1),
                           cpyIdx (idx);

            SortLists      (&cpyIdx, &revIdx);

            while (currentFrom < vDim-1) {
                currentTo = currentFrom + 1;
                while (currentTo < vDim) {
                    if ((*this)(revIdx.lData[currentTo], colIdx) > (*this)(revIdx.lData[currentFrom], colIdx)) {
                        break;
                    } else {
                        currentTo ++;
                    }
                }
                if (currentTo-currentFrom>2) {
                    theColumn.RecursiveIndexSort (currentFrom, currentTo-1, &idx);
                    didSort = true;
                }

                currentFrom = currentTo;
            }

            if (!didSort) {
                break;
            }
        }
    }

    _Matrix                 *result     = new _Matrix (hDim, vDim, theIndex, 1);

    if (theIndex) {
        _SimpleList    revIdx (hDim,0,1);
        SortLists (&idx, &revIdx);
        for (long r=0; r<lDim; r++) {
            long oi = theIndex[r];

            if (oi >= 0) {
                long     v  = oi%vDim,
                         h  = oi/vDim,
                         ni = revIdx.lData[h]*vDim+v;

                (*result)[ni] = theData[r];
            }
        }
    } else
        for (long r=0; r<hDim; r++) {
            long remapped = idx.lData[r];
            remapped *= vDim;
            for (long c=r*vDim; c<r*vDim+vDim; c++, remapped++) {
                result->theData[c] = theData[remapped];
            }
        }

    if (!result) {
        checkPointer (result);
    }

    return result;
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::PoissonLL (_PMathObj mp)
{
    if (storageType!=1) {
        _String    errMsg ("Only numeric matrices can be passed to Poisson Log-Likelihood");
        WarnError  (errMsg);
        return new _Constant (0.0);
    }

    if (mp->ObjectClass () != NUMBER || mp->Value() < 0.0) {
        _String    errMsg ("Invalid Poisson distribution parameter");
        errMsg = errMsg & _String((_String*)mp->toStr());
        WarnError  (errMsg);
        return new _Constant (0.0);
    }

    _Parameter     loglik = 0.0,
                   *logFactorials = new _Parameter [101],
    lambda        = mp->Value(),
    logLambda      = log (lambda),
    log2p         = log (sqrt(8.*atan(1.)));

    checkPointer (logFactorials);

    logFactorials[0] = 0.;
    logFactorials[1] = 0.;

    long           maxFactorialDone = 1;

    for (long idx = 0; idx < lDim; idx++) {
        long  cellValue = 0;
        if (theIndex) {
            cellValue = theIndex[idx];
            if (cellValue<0) {
                continue;
            }

            cellValue = theData[cellValue];
        } else {
            cellValue = theData[idx];
        }

        if (cellValue>=0) {
            if (maxFactorialDone>=cellValue) {
                loglik += logLambda * cellValue - lambda - logFactorials [cellValue];
            } else {
                if (cellValue<=100) {
                    for (long idx2 = maxFactorialDone+1; idx2 <= cellValue; idx2++) {
                        logFactorials[idx2] = logFactorials[idx2-1]+log((_Parameter)idx2);
                    }
                    loglik += logLambda * cellValue - lambda - logFactorials [cellValue];
                    maxFactorialDone = cellValue;
                } else
                    // use Stirling's formula
                {
                    loglik += logLambda * cellValue - lambda + cellValue - (cellValue+0.5)*log((_Parameter)cellValue)-log2p;
                }
            }
        }
    }

    delete      [] logFactorials;

    return new _Constant (loglik);
}


//_____________________________________________________________________________________________
_PMathObj       _Matrix::PathLogLikelihood (_PMathObj mp)
{
    _Matrix                 *m          = nil;

    _String                 errMsg;

    if (storageType!=1 || hDim != 3) {
        errMsg = ("First argument in call to < (PathLogLikelihood) must be a numeric 3xN matrix");
    } else {
        errMsg = ("Second argument in call to < (PathLogLikelihood) must be a square matrix");
        if (mp->ObjectClass () == MATRIX) {
            m = (_Matrix*)mp->Compute();
            if (m->GetHDim() == m->GetVDim()) {
                errMsg = empty;
            }
        }
    }

    if (errMsg.sLength) {
        WarnError  (errMsg);
        return new _Constant (0.);
    }

    CheckIfSparseEnough     (true);

    _Parameter              res     = 0.0;
    long                    maxDim  = m->GetHDim();

    for (long step = 0; step < vDim; step++) {
        long       i1 = theData[step],
                   i2 = theData[vDim+step];

        _Parameter t  = theData[2*vDim+step];

        if (i1<0 || i2 < 0 || i1 >= maxDim || i2 >= maxDim || t<0.0) {
            errMsg = _String ("An invalid transition in step ") & (step+1) & " of the chain: " & i1 & " to " & i2 & " in time " & t;
            WarnError  (errMsg);
            return new _Constant (0.);
        }

        _Matrix         rateMx (*m);
        rateMx *= t;
        _Matrix   * tMatrix = rateMx.Exponentiate ();

        t = tMatrix->theData[maxDim*i1+i2];
        DeleteObject (tMatrix);

        if (t>0.0) {
            res += log (t);
        } else {
            return new _Constant (-1.e300);
        }
    }

    return new _Constant (res);
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::pFDR (_PMathObj classes)
{
    _String         errMsg;
    long            k,
                    steps     = 20,
                    iterCount = 500;

    _Parameter      pVal = 0.0,
                    maxLambda = 0.0;


    if (theIndex) {
        CheckIfSparseEnough (true);
    }

    if (storageType!=1) {
        errMsg = "Only numeric matrices can be passed to && (pFDR)";
    } else {
        if ((GetVDim () != 1 && GetHDim() != 1) || GetVDim()*GetHDim() < 1) {
            errMsg = "The first argument of && (pFDR) must be an Nx1 matrix.";
        } else if (classes->ObjectClass () != NUMBER || classes->Value() > 1. || (pVal = classes->Value()) < 0.0) {
            errMsg = _String ("Invalid baseline p-value (must be in (0,1)):") & _String((_String*)classes->toStr());
        } else {
            for (long i=1; i<lDim; i++) {
                _Parameter pCount = theData[i];
                if (pCount < 0.0 || pCount > 1.0) {
                    errMsg = _String ("Invalid p-value entry in matrix passed to pFDR (must be a positive integer):");
                }
                if (pCount > maxLambda) {
                    maxLambda = pCount;
                }
            }
        }
    }


    if (errMsg.sLength) {
        WarnError  (errMsg);
        return new _Constant (0.0);
    }

    _Matrix        lamdbaRange (steps,1,false,true),
                   pFDRs       (steps,1,false,true);

    _Parameter     anLamdba           = 0.0,
                   minPFDR          = 5.0,
                   uberPFDR        = 0.0,
                   uberPFDRUpperLimit = 0.0,
                   minMSE             = 1.e100,
                   aStep            = 1.0/steps;


    k = 0;
    while (anLamdba<1.0) {
        lamdbaRange.theData[k] = anLamdba;

        if ((pFDRs.theData[k] = computePFDR (anLamdba, pVal))<minPFDR) {
            minPFDR = pFDRs.theData[k];
        }

        k++;
        anLamdba += aStep;
    }

    for (k=0; k<steps; k++) {
        _Parameter mse    = 0.0;
        _Matrix    ITpDFR (iterCount,1,false,true);

        for (long it = 0; it < iterCount; it = it+1) {
            _Matrix         sampledPs (lDim,1,false,true);
            _SimpleList     sample    (lDim,0,1);
            sample.PermuteWithReplacement (1);

            for (long el = 0; el < lDim; el++) {
                sampledPs.theData[el] = theData[sample.lData[el]];
            }

            ITpDFR.theData[it] = sampledPs.computePFDR (lamdbaRange.theData[k], pVal);

            mse += (ITpDFR.theData[it]-minPFDR)*(ITpDFR.theData[it]-minPFDR);
        }

        mse /= iterCount;

        if (mse < minMSE) {
            minMSE = mse;
            uberPFDR = pFDRs.theData[k];
            _Constant  zer (0.0);
            _Matrix* sorted = (_Matrix*)ITpDFR.SortMatrixOnColumn (&zer);
            uberPFDRUpperLimit = sorted->theData[((long)(0.95*iterCount))];
            DeleteObject (sorted);
        }
    }

    _Matrix * resMx = new _Matrix(2,1,false,true);
    checkPointer (resMx);
    resMx->theData[0] = uberPFDR;
    resMx->theData[1] = uberPFDRUpperLimit;

    return resMx;
}

//_____________________________________________________________________________________________
_Parameter      _Matrix::computePFDR (_Parameter lambda, _Parameter gamma)
// assumes a non-sparse row/column matrix
{
    long        rejected    = 0,
                null         = 0;

    for (long idx = 0; idx < lDim; idx++) {
        if (theData[idx] <= gamma) {
            rejected++;
        }
        if (theData[idx] > lambda) {
            null++;
        }
    }

    if (null) {
        _Parameter pi_0 = null/(lDim*(1.-lambda)),
                   pr_p = 0;


        if (rejected) {
            pr_p = rejected/(_Parameter)lDim;
        } else {
            pr_p = 1./(_Parameter)lDim;
        }

        return     pi_0 * gamma / (pr_p /** (1.-exp(log(1.-gamma)*lDim))*/);
    } else {
        return 1;
    }
}

//_____________________________________________________________________________________________

_PMathObj       _Matrix::Random (_PMathObj kind)
{
    _String     errMsg;

    long myVDim = GetVDim(),
         myHDim = GetHDim();

    if (kind->ObjectClass() == NUMBER) {
        bool    resample = (kind->Compute()->Value()>0);


        _SimpleList     remapped (myVDim,0,1);

        if (resample) {
            remapped.PermuteWithReplacement(1);
        } else {
            remapped.Permute(1);
        }


        if (storageType==1) {   // numeric matrix
            _Matrix * res = new _Matrix (GetHDim(), GetVDim(),theIndex,true);
            checkPointer (res);

            if (!theIndex)
                for (long vv = 0; vv<lDim; vv+=myVDim)
                    for (long k2=0; k2<remapped.lLength; k2++) {
                        res->theData[vv+k2] = theData[vv+remapped.lData[k2]];
                    }
            else {
                for (long vv = 0; vv<myHDim; vv++)
                    for (long k=0; k<remapped.lLength; k++) {
                        long ki = remapped.lData[k];
                        if ((ki = Hash (vv,ki)) >= 0) {
                            res->Store (vv,k,theData[ki]);
                        }
                    }
            }
            return res;
        } else {            // formula matrix
            if (storageType==2) {
                _Matrix * res = new _Matrix (GetHDim(), GetVDim(),theIndex,false);
                checkPointer (res);

                for (long vv = 0; vv<myHDim; vv++)
                    for (long k=0; k<remapped.lLength; k++) {
                        long ki = remapped.lData[k];
                        _Formula * ff = GetFormula (vv,ki);
                        if (ff) {
                            res->StoreFormula (vv, k, *ff);
                        }
                    }
                return res;
            }
        }
    }

    else if (kind->ObjectClass() == ASSOCIATIVE_LIST) {
        ReportWarning (_String("_Matrix::Random() with associative list as first argument."));

        // Associative list should contain following arguments:
        //  "PDF" - string corresponding to p.d.f. ("Gamma", "Normal")
        //  "ARG0" ... "ARGn" - whatever parameter arguments (matrices) are required for the p.d.f.
        _AssociativeList    * pdfArgs   = (_AssociativeList *)kind;
        _List               * keys      = pdfArgs->GetKeys();
        _String             pdfkey      ("PDF"),
                            * arg0      = (_String *)(*keys)(0);

        if (arg0->Equal(&pdfkey)) {
            _String     pdf ((_String *) (pdfArgs->GetByKey(pdfkey,STRING))->toStr()),
                        arg ("ARG0");

            if (pdf == _String("Dirichlet")) {
                return (_Matrix *) DirichletDeviate();
            } else if (pdf == _String("Gaussian")) {
                return (_Matrix *) GaussianDeviate (*(_Matrix *) pdfArgs->GetByKey (arg, MATRIX));
            } else if (pdf == _String("Wishart")) {
                return (_Matrix *) WishartDeviate (*(_Matrix *) pdfArgs->GetByKey (arg, MATRIX));
            } else if (pdf == _String("InverseWishart")) {
                return (_Matrix *) InverseWishartDeviate (*(_Matrix *) pdfArgs->GetByKey (arg, MATRIX));
            } else if (pdf == _String("Multinomial")) {
                return (_Matrix *) MultinomialSample ((_Constant *) pdfArgs->GetByKey (arg, NUMBER));
            } else {
                errMsg = _String("String argument passed to Random not a supported PDF: '") & pdf & "'";
            }
        } else {
            errMsg = _String("Expecting \"PDF\" key in associative list argument passed to Random(), received: ") & ((_String *)(*keys)(0))->getStr();
        }

    } else if (kind->ObjectClass () == STRING) {
        _String key = *((_FString*)kind->Compute())->theString;
        if (key == _String("LHS"))
            // latin hypercube sampling: samples are in ROWS
        {
            _Matrix * lhc = new _Matrix (myHDim, myVDim, false, true);

            _SimpleList permutation (myHDim,0,1);

            for (long c = 0; c < myVDim; c++) {
                permutation.Permute (1);
                for (long r = 0; r < myHDim; r++) {
                    lhc->theData[r*myVDim + c] = theData[permutation.lData[r]*myVDim + c];
                }
            }

            return lhc;
        }
        errMsg = _String ("Invalid string argument passed to matrix Random :") & key;
    } else {
        errMsg = _String ("Invalid argument passes to matrix Random (should be a number, an associative list or a string):") & _String((_String*)kind->toStr());
    }

    // error handling
    WarnError (errMsg);
    return new _Matrix (1,1);
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::K_Means (_PMathObj classes)
{
    _String         errMsg;
    _Matrix     *   arg;
    long            clusterCount,
                    iterCount,
                    dataPoints = 0;

    if (theIndex) {
        CheckIfSparseEnough (true);
    }

    if (storageType!=1) {
        errMsg = "Only numeric matrices can be passed to <= (K-means)";
    } else {
        if (GetVDim () != 2) {
            errMsg = "The first argument of <= (K-means) must be an Nx2 matrix, with samples in the first columns, and counts in the 2nd.";
        } else if (classes->ObjectClass () != MATRIX) {
            errMsg = _String ("Invalid number of clusters is call to K-means (must be >=1):") & _String((_String*)classes->toStr());
        } else {
            arg = (_Matrix*)classes->Compute();
            if (arg->GetVDim () != 1 || arg->GetHDim () != 2 || (clusterCount=arg->theData[0]) < 1 || (iterCount = arg->theData[1]) < 1) {
                errMsg = _String ("Invalid second argument is call to K-means (must be a 2x1 matrix of positive integers):") & _String((_String*)classes->toStr());
            } else {
                for (long i=1; i<lDim; i+=2) {
                    long pCount = theData[i];
                    if (pCount <= 0) {
                        errMsg = _String ("Invalid count entry in matrix passed to K-means (must be a positive integer):");
                    }
                    dataPoints += pCount;
                }
            }
        }
    }


    if (errMsg.sLength) {
        WarnError  (errMsg);
        return new _Matrix (1,1);
    }

    _Matrix * res = new _Matrix (2, clusterCount, false, true);
    checkPointer (res);

    if (clusterCount == 1) {
        _Parameter sampleMean    = 0.,
                   errorEstimate = 0.;

        for (long c1=0, c2=1; c1 < 2*hDim; c1+=2, c2+=2) {
            sampleMean +=  theData[c1] * theData[c2];
        }

        sampleMean /= dataPoints;

        {
            for (long c1=0, c2=1; c1 < 2*hDim; c1+=2, c2+=2) {
                _Parameter locErr = theData[c1] - sampleMean;
                //if (locErr > 0.0)
                errorEstimate += locErr*locErr*theData[c2];
                //else
                //errorEstimate -= locErr*theData[c2];
            }
        }

        res->theData[0] = sampleMean;
        res->theData[1] = errorEstimate;
    } else {
        _Matrix     clusterBreaks       (clusterCount,1,false,true),
                    lastClusterBreaks   (clusterCount,1,false,true),
                    clusterMeans        (clusterCount,2,false,true);

        _Parameter  minError    = 1.e100;
        long        hitMinError,
                    toggle      = 0;

        _SimpleList randomMapper    ((unsigned long)dataPoints,0,0),
                    allowedPoints (hDim+1,1,0);


        for (long c1 = 0; c1 < hDim; c1++) {
            hitMinError = theData[c1*2+1];

            for (long c2 = 0; c2 < hitMinError; c2++) {
                randomMapper.lData[toggle+c2] = c1;
            }

            toggle += hitMinError;

        }

        allowedPoints.lData[hDim] = 0;

        hitMinError = 0;
        toggle      = 0;

        for (long sampleCount = 0; sampleCount < iterCount; sampleCount ++) {
            _SimpleList       chosenMeans (clusterCount,0,0);

            for (long cc = 0; cc < clusterCount; cc = cc+1) {
                long sfrom = hDim;

                while (allowedPoints.lData[sfrom] == 0) {
                    sfrom = randomMapper.lData[(long)(genrand_real2()*dataPoints)];
                }

                chosenMeans.lData[cc] = sfrom;
                allowedPoints.lData[sfrom] = 0;
            }

            chosenMeans.Sort ();
            {
                for (long cc = 0; cc < clusterCount; cc = cc+1) {
                    clusterMeans.theData[toggle+cc] = theData[2*chosenMeans.lData[cc]];
                    allowedPoints.lData[chosenMeans.lData[cc]] = 1;
                }
            }

            _Parameter            lastErrorEstimate = 1.e100,
                                  errorEstimate = 0.;

            for (long cIters = 0; cIters < 250; cIters = cIters + 1) {
                if (fabs (lastErrorEstimate-errorEstimate) > 0.0001 || (clusterBreaks-lastClusterBreaks).MaxElement() > 0.0001 ) {
                    lastClusterBreaks = clusterBreaks;
                    lastErrorEstimate = errorEstimate;

                    {
                        long    breakPoint = 0,
                                k = 0,
                                kp1 = 1;

                        for (; k<clusterCount-1 && breakPoint < hDim; k=k+1, kp1++) {
                            _Parameter cm4 = clusterMeans.theData[toggle+k],
                                       cm2 = clusterMeans.theData[toggle+kp1],
                                       cm3 = theData[2*breakPoint];

                            while (fabs(cm3-cm4)<fabs(cm3-cm2)) {
                                breakPoint ++;
                                if (breakPoint == hDim) {
                                    break;
                                }
                                cm3 = theData[2*breakPoint];
                            }
                            clusterBreaks.theData[k] = breakPoint;
                        }
                        clusterBreaks.theData[k] = hDim;
                    }

                    errorEstimate = 0.0;

                    long        currentStart = 0;

                    for (long cm1 = 0; cm1 < clusterCount && currentStart < hDim-1; cm1 = cm1+1) {
                        _Parameter  cm2 = 0.;


                        long        k,
                                    cm4 = 0;

                        for (k=currentStart; k<clusterBreaks.theData[cm1]; k=k+1) {
                            cm2 += theData[k*2] * theData[k*2+1];
                            cm4 += theData[k*2+1];
                        }

                        long  cm3 = k-currentStart;

                        if (cm3>0) {
                            cm2 /= cm4;

                            clusterMeans.theData[cm1+toggle] = cm2;

                            for (long k=currentStart; k<clusterBreaks.theData[cm1]; k=k+1) {
                                _Parameter locError = theData[2*k] - cm2;
                                //if (locError < 0.0)
                                //  errorEstimate -= locError*theData[2*k+1];
                                //else
                                errorEstimate += locError*locError*theData[2*k+1];
                            }

                            currentStart = clusterBreaks.theData[cm1];
                        }
                    }
                } else {
                    break;
                }
            }

            if (minError == 0.0 || fabs((minError-errorEstimate)/minError) < 0.001) {
                hitMinError ++;
            } else if (errorEstimate < minError) {
                hitMinError = 1;
                toggle      = (toggle+clusterCount)%(2*clusterCount);
                minError    = errorEstimate;
            }
        }

        toggle  = (toggle+clusterCount)%(2*clusterCount);
        for (long k2 = 0; k2 < clusterCount; k2++) {
            res->theData[k2] = clusterMeans.theData[toggle+k2];
        }
        res->theData[clusterCount]   = minError;
        res->theData[clusterCount+1] = hitMinError;
    }

    return res;
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::ProfileMeanFit (_PMathObj classes)
{
    _String         errMsg;
    _Matrix     *   arg;
    long            weightClasses;

    _Parameter      dataPoints = 0.;

    if (theIndex) {
        CheckIfSparseEnough (true);
    }

    if (storageType!=1) {
        errMsg = "Only numeric matrices can be passed to <= (K-means)";
    } else {
        if (GetHDim () != 2) {
            errMsg = "The first argument of ProfileMeanFit must be an 2xN matrix, with samples in the first row, and counts in the 2nd.";
        } else if (classes->ObjectClass () != MATRIX) {
            errMsg = _String ("Invalid second argument for ProfileMeanFit (must be a column vector):") & _String((_String*)classes->toStr());
        } else {
            arg = (_Matrix*)classes->Compute();
            if (arg->GetVDim () != 1) {
                errMsg = _String ("Invalid second argument is call to ProfileMeanFit (must be a column vector):") & _String((_String*)classes->toStr());
            } else {
                weightClasses = arg->GetHDim ();

                for (long i=vDim; i<lDim; i++) {
                    long pCount = theData[i];
                    if (pCount <= 0) {
                        errMsg = _String ("Invalid count entry in matrix passed to ProfileMeanFit (must be a positive integer):");
                    }
                    dataPoints += pCount;
                }
            }
        }
    }


    if (errMsg.sLength) {
        WarnError  (errMsg);
        return new _Matrix (1,1);
    }

    _Matrix * res           = new _Matrix (4, weightClasses, false, true);
    //_SimpleList               splitRuns (weightClasses,0,0);
    checkPointer (res);

    _Parameter      runningSum      = 0.,
                    targetSum        = arg->theData[0],
                    valueSum      = 0.,
                    logLikelihood   = 0.,
                    varMult;

    checkParameter  (PROFILE_MEAN_VAR_MULT, varMult, 1.);



    long            currentIndex    = 0,
                    currentSlider = 0,
                    runningSize     = 1,
                    currentSpan      = theData[vDim+currentIndex],
                    runningOffset = 0;

    while (currentIndex < vDim - 1) {
        runningSum = runningSum + theData[vDim+currentIndex]/dataPoints;

        if ((runningSum >= targetSum)||(vDim-currentIndex <= weightClasses - currentSlider)) {
            res->theData[currentSlider]                 = currentIndex;
            res->theData[weightClasses+currentSlider]   = runningSize;
            res->theData[weightClasses*2+currentSlider] =
                (theData[currentIndex]*theData[vDim+currentIndex]+valueSum)/(currentSpan+theData[vDim+currentIndex]);

            //splitRuns.lData[currentSlider] = currentSpan;
            runningSize   = 1;
            valueSum      = 0.;
            currentSlider ++;
            targetSum     = targetSum + arg->theData[currentSlider];
            currentSpan   = 0;
        } else {
            valueSum    +=  theData[currentIndex]*theData[vDim+currentIndex];
            runningSize ++;
            currentSpan +=  theData[vDim+currentIndex];
        }
        currentIndex++;
    }

    currentSpan += theData[vDim+currentIndex];
    valueSum    += theData[currentIndex]*theData[vDim+currentIndex];

    res->theData[currentSlider]                 = currentIndex;
    res->theData[weightClasses+currentSlider]   = runningSize;
    res->theData[2*weightClasses+currentSlider] = valueSum/currentSpan;
    //splitRuns.lData[currentSlider] = currentSpan;

    currentIndex    = 0;
    currentSlider   = 0;
    runningOffset   = 0;

    _Matrix          REWEIGHTED_MATRIX (vDim,1,false,true);

    while (currentSlider < weightClasses) {
        long        classSize   = res->theData[weightClasses+currentSlider];
        _Parameter  classWeight = arg->theData[currentSlider];//splitRuns.lData[currentSlider]/dataPoints;

        if (classWeight > 0.0) {
            if (classSize == 1) {
                logLikelihood += theData[vDim+runningOffset] * log (classWeight);
            } else {
                _Parameter      classMean       = res->theData[2*weightClasses+currentSlider],
                                //classNorm         = 0.,
                                classVar        = (fabs(classMean)>0.05)?0.5/(varMult*fabs(classMean)):0.5/(varMult*0.025);

                currentIndex    = runningOffset+classSize;

                for (long   reslider = runningOffset; reslider < currentIndex; reslider = reslider+1) {
                    targetSum = theData[reslider]-classMean;
                    targetSum = -targetSum*targetSum*classVar;
                    REWEIGHTED_MATRIX.theData[reslider] = targetSum;
                    //classNorm += exp(targetSum);
                }

                classWeight = log (classWeight);
                {
                    for (long reslider = runningOffset; reslider < currentIndex; reslider = reslider+1) {
                        logLikelihood += (REWEIGHTED_MATRIX.theData[reslider]+classWeight)*theData[vDim+reslider];
                    }
                }
            }
        } else {
            if (classSize>0) {
                logLikelihood = -1e100;
                break;
            }
        }
        runningOffset += classSize;
        currentSlider++;
    }

    res->theData[3*weightClasses] = logLikelihood;
    return res;
}

//_____________________________________________________________________________________________
void            _Matrix::PopulateConstantMatrix (const _Parameter v)
{
    if (storageType == 1)
        for (long r=0; r<lDim; r++) {
            theData[r] =v;
        }
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::AddObj (_PMathObj mp)
{
    if (_Matrix::ObjectClass()!=mp->ObjectClass()) {
        if (mp->ObjectClass () == STRING) {
            _Matrix * convMatrix = new _Matrix (*((_FString*)mp)->theString),
            * res;
            checkPointer (convMatrix);
            res = (_Matrix*)AddObj (convMatrix);
            DeleteObject (convMatrix);
            return res;
        }
        if (mp->ObjectClass () == NUMBER) {
            _Matrix* aNum = (_Matrix*)ComputeNumeric ();
            if (aNum->storageType == 1) {
                _Matrix * plusStuff = new _Matrix (hDim,vDim,false,true);
                checkPointer (plusStuff);
                _Parameter plusValue = mp->Value();

                if (theIndex) {
                    for (long k=0; k<hDim*vDim; k++) {
                        plusStuff->theData[k] = plusValue;
                    }

                    for (long l=0; l<lDim; l++) {
                        long rI = theIndex[l];
                        if (rI>0) {
                            plusStuff->theData[rI] += theData[l];
                        }
                    }
                } else
                    for (long r=0; r<lDim; r++) {
                        plusStuff->theData[r] = theData[r] + plusValue;
                    }

                return plusStuff;
            }
        }

        warnError( -101);
        return new _Matrix (1,1);
    }

    _Matrix * m = (_Matrix*)mp;
    AgreeObjects (*m);
    _Matrix * result = new _Matrix (hDim, vDim, bool((theIndex!=nil)&&(m->theIndex!=nil)), storageType);
    if (!result) {
        checkPointer (result);
    }
    Add (*result,*m);
    return result;
}

//_____________________________________________________________________________________________
void        _Matrix::operator -= (_Matrix& m)
{
    AgreeObjects (m);
    if ((!m.theIndex)&&theIndex) {
        CheckIfSparseEnough(true);
    }
    Subtract (*this,m);
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::SubObj (_PMathObj mp)
{
    if (mp->ObjectClass()!=ObjectClass()) {
        warnError( -101);
        return new _Matrix (1,1);
    }

    _Matrix * m = (_Matrix*)mp;
    AgreeObjects (*m);
    _Matrix * result = new _Matrix (hDim, vDim, bool((theIndex!=nil)&&(m->theIndex!=nil)), storageType);
    if (!result) {
        checkPointer (result);
    }
    Subtract (*result,*m);
    return result;
}

//_____________________________________________________________________________________________
void        _Matrix::operator *= (_Parameter c)
{
    Multiply (*this,c);
}

//_____________________________________________________________________________________________
_Matrix     _Matrix::operator * (_Parameter c)
{
    _Matrix result (*this);
    Multiply (result,c);
    return result;
}

//_____________________________________________________________________________________________
void        _Matrix::operator *= (_Matrix& m)
{
    CheckDimensions     (m);
    AgreeObjects        (m);
    _Matrix   result    (hDim, m.vDim, false, storageType);
    Multiply            (result,m);
    //if ((theIndex!=nil)||(m.theIndex!=nil)) result.AmISparse();
    if (theIndex!=nil && m.theIndex!=nil) {
        result.AmISparse();
    }
    Swap                (result);
}

//_____________________________________________________________________________________________
void        _Matrix::MultbyS (_Matrix& m, bool leftMultiply, _Matrix* externalStorage, _Parameter* stash)
{
    _Matrix * result = nil;
    if (!externalStorage) {
        result = new _Matrix (hDim, m.vDim, false, storageType);
    }

    _Matrix * receptacle = (externalStorage?externalStorage:result);

    if (leftMultiply) {
        m.Multiply (*receptacle,*this);
    } else {
        Multiply   (*receptacle,m);
    }

    if (theIndex&&m.theIndex) {
        if (receptacle->AmISparseFast(*this) == false) {
            Swap            (*receptacle);
        } else {
            CompressSparseMatrix(false,stash);
        }
    } else { // both dense
        Swap            (*receptacle);
    }

    if (!externalStorage) {
        DeleteObject (result);
    } else {
        externalStorage->CheckIfSparseEnough (true);
        memset (externalStorage->theData, 0, sizeof (_Parameter)*externalStorage->lDim);
        //for (long s = 0; s < externalStorage->lDim; s++) externalStorage->theData[s] = 0.0;
    }
}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::MultObj (_PMathObj mp)
{

    if (mp->ObjectClass()!=ObjectClass())
        if (mp->ObjectClass()!=NUMBER) {
            warnError(-101);
            return new _Matrix (1,1);
        } else {
            _Parameter theV = mp->Value();
            return (_PMathObj)((*this)*theV).makeDynamic();
        }

    _Matrix*        m = (_Matrix*)mp;
    CheckDimensions (*m);
    AgreeObjects    (*m);

    _Matrix*      result = new _Matrix (hDim, m->vDim, false, storageType);
    checkPointer  (result);
    Multiply      (*result,*m);
    return        result;

}

//_____________________________________________________________________________________________
_PMathObj       _Matrix::MultElements (_PMathObj mp)
{

    if (mp->ObjectClass()!=ObjectClass()) {
        warnError(-101);
        return new _Matrix (1,1);
    }

    _Matrix* m = (_Matrix*)mp;

    if ((GetHDim()!=m->GetHDim()) || (GetVDim()!=m->GetVDim())) {
        WarnError ("Element-wise multiplication requires matrixes of the same dimension.");
        return new _Matrix (1,1);
    }

    if ((storageType!=1)||(m->storageType != 1)) {
        WarnError ("Element-wise multiplication only works on numeric matrices");
        return new _Matrix (1,1);
    }

    _Matrix*      result = new _Matrix (hDim, vDim, false, true);
    checkPointer  (result);

    if (theIndex) {
        if (m->theIndex) {
            for (long k=0; k<lDim; k++) {
                long    i = theIndex[k];
                if (i>=0) {
                    result->theData [i] = theData[k] * (*m)(i/vDim, i%vDim);
                }
            }
        } else {
            for (long k=0; k<lDim; k++) {
                long    i = theIndex[k];
                if (i>=0) {
                    result->theData [i] = theData[k] * m->theData[i];
                }
            }
        }
    } else {
        if (m->theIndex) {
            for (long k=0; k<m->lDim; k++) {
                long    i = m->theIndex[k];
                if (i>=0) {
                    result->theData [i] = theData[i] * m->theData[k];
                }
            }
        } else
            for (long k=0; k<lDim; k++) {
                result->theData [k] = theData[k] * m->theData[k];
            }
    }

    if (theIndex||m->theIndex) {
        result->AmISparse();
    }

    return  result;
}

//_____________________________________________________________________________________________
void    _Matrix::CheckDimensions (_Matrix& secondArg)
// check matrix dimensions to ensure that they are multipliable
{
    if (vDim!=secondArg.hDim) {
        if (hDim == 1 && secondArg.hDim==1 && vDim == secondArg.vDim) { // handle scalar product separately
            secondArg.Transpose();
        } else {
            char str[255];
            sprintf (str,"Incompatible matrix dimensions in call to CheckDimension: %ldx%ld and %ldx%ld\n",hDim,vDim,secondArg.hDim,secondArg.vDim);
            ReportWarning (str);
            return;
        }
    }
}

//_____________________________________________________________________________________________
_Matrix     _Matrix::operator * (_Matrix& m)
{
    CheckDimensions (m);
    AgreeObjects (m);
    _Matrix result (hDim, m.vDim, false, storageType);
    Multiply (result,m);
    if ((theIndex!=nil)||(m.theIndex!=nil)) {
        result.AmISparse();
    }
    return result;

}
//_____________________________________________________________________________________________
_Matrix     _Matrix::operator + (_Matrix& m)
{
    AgreeObjects (m);
    _Matrix result (hDim, vDim, bool((theIndex!=nil)&&(m.theIndex!=nil)), storageType);
    Add (result,m);
    return result;

}
//_____________________________________________________________________________________________
_Matrix     _Matrix::operator - (_Matrix& m)
{
    AgreeObjects (m);
    _Matrix result (hDim, vDim, bool((theIndex!=nil)&&(m.theIndex!=nil)), storageType);
    Subtract (result,m);
    return result;
}


//_____________________________________________________________________________________________

BaseRef _Matrix::toStr(void)
{
    _String result(2048L,true);
    checkParameter (printDigitsSpec,printDigits,0);
    long digs = printDigits;

    //if (vDim<500)
    {
        if (storageType == 1 || (storageType == 2 && IsAStringMatrix())) {
            bool printStrings = storageType != 1;

            result << '{';
            result << '\n';
            for (long i = 0; i<hDim; i++) {
                result<<'{';
                char str[100];
                for (long j = 0; j<vDim; j++) {
                    if (printStrings) {
                        result << '"';
                        _Formula * f = GetFormula (i,j);
                        _PMathObj fv;
                        if (f && (fv=f->Compute())) {
                            result << ((_FString*)fv)->theString;
                        }
                        result << '"';
                    } else {
                        if (digs >= 0)
#ifdef __USE_LONG_DOUBLE__
                            sprintf(str, "%18.12Lg", (*this)(i,j));
#else
                            sprintf(str, "%18.12g", (*this)(i,j));
#endif
                        else
#ifdef __USE_LONG_DOUBLE__
                            sprintf(str, "%Lg", (*this)(i,j));
#else
                            sprintf(str, "%g", (*this)(i,j));
#endif
                        _String cell (str);
                        result<<&cell;

                    }
                    if (j<vDim-1) {
                        result<<',';
                    }
                }
                result<<'}';
                result<<'\n';
            }
            result<<'}';
            result<<'\n';
            result.Finalize();
        } else if (storageType == 0) {
            checkParameter (ANAL_COMP_FLAG, ANALYTIC_COMPUTATION_FLAG, 0);
            if (!ANALYTIC_COMPUTATION_FLAG) {
                return Compute()->toStr();
            }
            for (long i = 0; i<hDim; i++) {
                result<<'\n';
                result<<'[';
                for (long j = 0; j<vDim; j++) {
                    long p = Hash (i,j);
                    if (GetMatrixObject(p)) {
                        if (p>=0) {
                            _String *sp = (_String*) GetMatrixObject (p)->toStr();
                            result<<sp;
                            if (j<vDim-1) {
                                result<<',';
                            }
                            result<<' ';
                            DeleteObject (sp);
                            continue;
                        }
                    }
                    result<<'0';
                }
                result<<']';
            }
            result<<'\n';
            result<<'\n';
            result.Finalize();
        } else {
            _Matrix* eval = (_Matrix*)ComputeNumeric();
            result.Finalize();
            _String* ss = (_String*)eval->toStr();
            return ss;
        }


    }

    return result.makeDynamic();
}

//_____________________________________________________________________________________________

void     _Matrix::Serialize (_String& res, _String& myID)
{
    if (storageType) {
        res << '\n';
        res <<  myID;
        if (storageType == 1) {
            _String * mStr = (_String*)toStr();
            res << '=';
            res << *mStr;
            res << ';';
            DeleteObject (mStr);
        } else if (storageType == 2) {
            _String mHeader = _String ("={") & hDim & ',' & vDim & "};\n";
            res << mHeader;
            for (long h=0; h<hDim; h++)
                for (long v=0; v<vDim; v++) {
                    _Formula *theCell = GetFormula (h,v);
                    if (theCell&&(!theCell->IsEmpty())) {
                        _String * fStr = (_String*)theCell->toStr();
                        res << myID;
                        res << '[';
                        res << _String(h);
                        res << "][";
                        res << _String(v);
                        res << "]:=";
                        res << *fStr;
                        res << ";\n";
                        DeleteObject (fStr);
                    }
                }
        }
    }
}

//_________________________________________________________
void    _Matrix::toFileStr (FILE*dest)
{
    if (storageType == 1 || (storageType == 2 && IsAStringMatrix())) {
        bool printStrings = storageType != 1;
        long digs         = -1;

        if (!printStrings) {
            checkParameter (printDigitsSpec,printDigits,0);
            digs =  printDigits;
        }

        if (!printStrings && digs != -1) {
            _String formatStr;
            if (digs<=0 || digs>15) {
                digs = 8;
            }
            formatStr = "%";
            formatStr = formatStr&_String(digs+6)&'.'&_String(digs)&'g';
            char *fs = formatStr.getStr();
            fprintf (dest, "\n{");
            for (long i = 0; i<hDim; i++) {
                fprintf (dest, "{");
                for (long j = 0; j<vDim; j++) {
                    fprintf(dest, fs, (*this)(i,j));
                    if (j<vDim-1) {
                        fprintf (dest, ",");
                    }
                    if (j%100==0) {
                        fflush(dest);
                    }
                }
                fprintf (dest, "}\n");
            }
            fprintf (dest, "}\n");
        } else {
            fprintf (dest, "\n{");
            for (long i = 0; i<hDim; i++) {
                fprintf (dest, "{");
                for (long j = 0; j<vDim; j++) {
                    if (j) {
                        fprintf (dest,",");
                    }

                    if (printStrings) {
                        fprintf (dest,"\"");;
                        _Formula * f = GetFormula (i,j);
                        _PMathObj fv;
                        if (f && (fv=f->Compute())) {
                            fprintf (dest,"%s",((_FString*)fv)->theString->sData);
                        }
                        fprintf (dest,"\"");
                    } else {
                        fprintf(dest, "%g", (*this)(i,j));
                    }
                }
                fprintf (dest, "}\n");
            }
            fprintf (dest, "}\n");
        }
    } else if (storageType==0) {
        checkParameter (ANAL_COMP_FLAG, ANALYTIC_COMPUTATION_FLAG, 0);
        if (!ANALYTIC_COMPUTATION_FLAG) {
            Compute()->toFileStr(dest);
            return;
        }
        for (long i = 0; i<hDim; i++) {
            fprintf (dest, "\n[");
            for (long j = 0; j<vDim; j++) {
                long p = Hash (i,j);
                if (p>=0) {
                    _String *sp = (_String*) GetMatrixObject (p)->toStr();
                    fprintf(dest, "%s", sp->sData);
                    fprintf (dest, ",");
                    DeleteObject (sp);
                } else {
                    fprintf(dest, "%g", 0.0);
                }
            }
            fprintf (dest, "]");
        }
    } else {
        _Matrix* eval = (_Matrix*)(storageType==3?EvaluateSimple():Evaluate(false));
        eval->toFileStr(dest);
        DeleteObject (eval);
    }
}

//_____________________________________________________________________________________________

void    SetIncrement (int m)
{
    _Matrix::storageIncrement = m;
}
//_____________________________________________________________________________________________
void    _Matrix::InitMxVar (_SimpleList& mxVariables, _Parameter glValue)
{
    _Constant gv (glValue);
    for (long k=0; k<mxVariables.countitems(); k++) {
        _Variable* theV = LocateVar(mxVariables(k));
        theV->SetValue (&gv);
    }

}
//_____________________________________________________________________________________________
bool    _Matrix::ImportMatrixExp (FILE* theSource)
{
    long mDim=0,i,k=0,j,m;
    char buffer[255],fc=0;
    buffer[0]=0;
    while(1) {
        buffer[mDim]=fgetc(theSource);
        if (feof(theSource)) {
            return false;
        }
        if (buffer[mDim]==',') {
            break;
        }
        mDim++;
    }
    buffer[mDim] = 0;
    mDim = atol (buffer); // matrix dimension
    Clear();
    CreateMatrix (this,mDim,mDim,false,false);
    // read in the variables
    i = 0;
    _SimpleList varList,c1,c2;
    while (fc!=';') {
        fc = fgetc (theSource);
        if ((fc==',')||(fc==';')) {
            buffer [i] = 0;
            _String varName (buffer);

            _Variable * ppv = CheckReceptacle (&varName, empty, true);
            varList << ppv->GetAVariable();
            i = 0;
        } else {
            buffer[i]=fc;
            i++;
        }
        if (feof(theSource)) {
            return false;
        }
    }
    do {
        fc = fgetc (theSource);
        if (feof(theSource)) {
            return false;
        }
    } while (fc!=';');

    k = 0; // term counter

    while (k<mDim*mDim) {
        i = 0;
        _Polynomial* thisCell = new _Polynomial (varList);
        checkPointer(thisCell);
        while (fc!='{') {
            fc = fgetc (theSource);
            buffer[i] = fc;
            i++;
            if (feof(theSource)) {
                return false;
            }
        }
        m = atol (buffer);
        _Parameter* theCoeffs = (_Parameter*)MatrixMemAllocate(m*sizeof(_Parameter));
        j = 0;
        while (fc!='}') {
            i = 0;
            do {
                buffer[i] = fc = fgetc (theSource);
                i++;
                if (feof(theSource)) {
                    return false;
                }
            } while ((fc!=',')&&(fc!='}'));
            buffer[i]=0;
            theCoeffs[j]=atof (buffer);
            j++;
            if (j>m) {
                return false;
            }
        }
        _PolynomialData *pd = new _PolynomialData (varList.countitems(),j,theCoeffs);
        checkPointer(pd);
        MatrixMemFree (theCoeffs);
        fc = fgetc(theSource);
        if (fc != '{') {
            return false;
        }
        c1.Clear();
        while (fc!='}') {
            i = 0;
            do {
                buffer[i] = fc = fgetc (theSource);
                i++;
                if (feof(theSource)) {
                    return false;
                }
            } while ((fc!=',')&&(fc!='}'));
            buffer[i]=0;
            c1<<atol (buffer);
        }
        fc = fgetc(theSource);
        if (fc != '{') {
            return false;
        }
        c2.Clear();
        while (fc!='}') {
            i = 0;
            do {
                buffer[i] = fc = fgetc (theSource);
                i++;
                if (feof(theSource)) {
                    return false;
                }
            } while ((fc!=',')&&(fc!='}'));
            buffer[i]=0;
            c2<<atol (buffer);
        }
        thisCell->SetTheTerms(pd);
        thisCell->SetCLists (c1,c2);
        StoreObject(k,thisCell);
        k++;
    }

    return true;
}

//_____________________________________________________________________________________________

void    _Matrix::ExportMatrixExp (_Matrix* theBase, FILE* theDump)
// export the matrix's computational form in the following format
// matrix dimension followed by a comma
// a comma separated list of variable names followed by a semicolon
// a list of (precision, maxcap) followed by a semicolon
// for each matrix entry
// number of coeffs,
// a {} enclosed list of coefficients
// a {} enclosed first computational list
// a {} enclosed second computational list
// followed by a comma
{
    // write out the preliminaries
    if (storageType!=0) {
        warnError(-200);
        return;
    }
    fprintf (theDump,"%ld,",hDim);
    _SimpleList mxVariables;
    {
        _AVLList        mxA (&mxVariables);
        ScanForVariables(mxA,true);
        mxA.ReorderList();
    }


    long k, i=0;
    _Parameter* varPool = (_Parameter*)MatrixMemAllocate (mxVariables.countitems()*sizeof(_Parameter));
    for (k=0; k<mxVariables.countitems(); k++) {
        fprintf (theDump,"%s",LocateVar(mxVariables(k))->GetName()->sData);
        if (k<mxVariables.countitems()-1) {
            fprintf (theDump,"%c",',');
        } else {
            fprintf (theDump,"%c",';');
        }
        varPool[k]=topPolyCap;
    }

    // begin by computing the actual "numerical exponential"
    // initialize all the variables to the polycap value

    InitMxVar   (mxVariables, topPolyCap);

    _Matrix     *dummy = (_Matrix*)theBase->Evaluate(false);
    _Matrix     *numExp = (_Matrix*)(dummy->Exponentiate());

    DeleteObject(dummy);
    checkParameter (ANAL_MATRIX_TOLERANCE,analMatrixTolerance,1e-6);
    fprintf (theDump,"%g,%g;",analMatrixTolerance,topPolyCap);

    // now loop thru the cells and check the precision term by term
    for (k=0; k<lDim; k++) {
        _SimpleList termRank, termIndex,c1,c2;
        _Polynomial* thisCell = ((_Polynomial**)theData)[k];
        long nTerms = thisCell->GetTheTerms()->NumberOfTerms(),
             step = nTerms/10+1, upTo = step, tup,j;
        _Parameter* coeffHolder =  (_Parameter*)MatrixMemAllocate (nTerms*sizeof(_Parameter)), error, bestError = 1;

        thisCell->RankTerms(&termRank);
        for (i=0; i<nTerms; i++) {
            termIndex<<i;
        }
        SortLists (&termRank,&termIndex);
        termRank.Clear();
        for (i=0; i<nTerms; i++) {
            termRank<<(nTerms-termIndex.Find(i)-1);
        }
        bestError = 1;
        while(upTo<nTerms+step) {
            if (upTo<nTerms) {
                tup = upTo;
            } else {
                tup = nTerms-1;
            }
            termIndex.Clear();
            for (i=0,j=0; (i<nTerms)&&(j<=tup); i++) {
                if (termRank.lData[i]<=tup) {
                    coeffHolder[j]=thisCell->GetTheTerms()->GetCoeff(i);
                    j++;
                    termIndex<<i;
                }
            }
            thisCell->Convert2ComputationForm(&c1,&c2,&termIndex);
            error = fabs(thisCell->ComputeP(varPool,coeffHolder,thisCell->GetNoVariables()+1,c1.countitems(),c1.quickArrayAccess(),
                                            c2.quickArrayAccess())-numExp->directIndex(k));
            if (error<bestError) {
                bestError = error;
            }
            if (bestError<=analMatrixTolerance) {
                break;
            }
            upTo+=step;
        }

        if (bestError>analMatrixTolerance) {
            char be[100];
            sprintf (be,"%g",bestError);
            _String wm ("Polynomial Matrix Exp approximation failed tolerance test in cell (");
            wm = wm&_String(k/hDim)&","&_String(k%hDim)&"). Tolerance achieved is:"&be;
            ReportWarning (wm);
        }
        fprintf(theDump,"%ld{",tup+1);
        for (i=0; i<=tup; i++) {
            if (i) {
                fprintf(theDump,",%18.16g",coeffHolder[i]);
            } else {
                fprintf(theDump,"%18.16g",coeffHolder[i]);
            }
        }
        fprintf(theDump,"}%ld",tup);
        c1.toFileStr(theDump);
        c2.toFileStr(theDump);
        MatrixMemFree (coeffHolder);

    }
    MatrixMemFree (varPool);
    DeleteObject (numExp);

    /*  // scan the maximum number of terms in the polynomial cells
        for (k=0;k<lDim;k++)
        {
            _Polynomial* currentCell = ((_Polynomial**)theData)[k];
            currentCell->Compute();
            if (currentCell->ComputationalSize()>maxterms)
            {
                maxterms = currentCell->ComputationalSize();
            }
        }

        // divide the number of terms in 10 subgroups and compute the absolute error, then try for the maxcap
        _Parameter  stdCap = 0.0;
        if (enforcePolyCap)
            stdCap = topPolyCap;
        else
        {
            _Parameter ub,lb;
            for (k=0;k<mxVariables.countitems(); k++)
            {
                _Variable* theV = LocateVar(mxVariables(k));
                lb = fabs(theV->GetLowerBound());
                ub = fabs(theV->GetUpperBound());
                if (ub>stdCap)
                    stdCap = ub;
                if (lb>stdCap)
                    stdCap = lb;
            }
        }
        if (stdCap>10.0)
            stdCap = 10.0;

        long  partitionSize = maxterms/10+1, currentSize = maxterms-partitionSize;
        _Matrix     *polym, *numm, *numm1, *errMx;
        _Parameter  error =0.0 ,step = stdCap/10, tryCap = stdCap+step, stdError;
        InitMxVar (mxVariables,stdCap);
        numm = (_Matrix*)theBase->Evaluate(false);
        // first run the precision analysis on the max terms setting

        InitMxVar (mxVariables,stdCap);
        polym = (_Matrix*)Evaluate(false);

        errMx = (_Matrix*)polym->SubObj(&numm->Exponentiate());
        stdError = errMx->MaxElement();
        DeleteObject (errMx);
        DeleteObject (polym);
        stdError = exp(log(10.0)*long(log(stdError)/log(10.0)));

        while (error<=stdError)
        {
            InitMxVar (mxVariables,tryCap);
            polym = (_Matrix*)Evaluate(false);
            numm1 = (_Matrix*)theBase->Evaluate(false);
            _Matrix *errMx = (_Matrix*)polym->SubObj(&numm1->Exponentiate());
            error = errMx->MaxElement();
            DeleteObject (errMx);
            DeleteObject (polym);
            DeleteObject (numm1);
            tryCap+=step;
            if (tryCap>10.0) break;
        }

        fprintf(theDump,"%d,%g,%g;", maxterms, stdError,tryCap-step);


        while (currentSize>0)
        {
            SetPolyTermCap (currentSize);
            InitMxVar (mxVariables,stdCap);
            polym = (_Matrix*)Evaluate(false);
            _Matrix *errMx = (_Matrix*)polym->SubObj(&numm->Exponentiate());
            error = errMx->MaxElement();
            DeleteObject (errMx);
            DeleteObject (polym);
            error = exp(log(10.0)*long(log(error)/log(10.0)));
            tryCap = stdCap;
            while (error>stdError)
            {
                InitMxVar (mxVariables,tryCap);
                polym = (_Matrix*)Evaluate(false);
                numm1 = (_Matrix*)theBase->Evaluate(false);
                errMx = (_Matrix*)polym->SubObj(&numm1->Exponentiate());
                error = errMx->MaxElement();
                DeleteObject (errMx);
                DeleteObject (polym);
                DeleteObject (numm1);
                tryCap-=step;
                if (tryCap<0) break;
            }
            if ((tryCap<0)||(error>stdError)) break;
            fprintf(theDump,"%d,%g,%g;", currentSize, stdError,tryCap+step);
            currentSize-=partitionSize;
        }

        DeleteObject (numm);
        SetPolyTermCap (0x1fffffff);

        for (k=0;k<lDim;k++)
        {
            if (k)
                fprintf(theDump,",");
            _Polynomial* currentCell = ((_Polynomial**)theData)[k];
            fprintf(theDump,"{");
            _PolynomialData* cD = currentCell->GetTheTerms();
            for (long l=0; l<cD->GetNoTerms();l++)
            {
                if (l)
                    fprintf(theDump,",%18.16g",cD->GetCoeff(l));
                else
                    fprintf(theDump,"%18.16g",cD->GetCoeff(l));
            }
            fprintf(theDump,"}");
            currentCell->toFileStr(theDump);
        }*/

}

//_____________________________________________________________________________________________

_Parameter  _Matrix::ExpNumberOfSubs  (_Matrix* freqs, bool mbf)
{
    if (storageType!=1 || freqs->storageType!=1 || hDim!=vDim) {
        return 0.0;
    }

    _Parameter      result      =   0.0;
    _Matrix         *nf         =   nil,
                     *stencil =   branchLengthStencil();

    if (freqs->theIndex) {
        nf = new _Matrix (*freqs);
        checkPointer (nf);
        nf->CheckIfSparseEnough (true);
    } else {
        nf = freqs;
    }


    if (theIndex) {
        _Parameter*   diags = new _Parameter[hDim];
        checkPointer (diags);
        for (long k=0; k<hDim; k++) {
            diags[k] = 0.0;
        }

        if (mbf) {
            for (long k=0; k<lDim; k++) {
                long l=theIndex[k];
                if (l>=0) {
                    long i = l%vDim;
                    long j = l/vDim;

                    if (i!=j && (!stencil || stencil->theData[l] > 0.0)) {
                        diags[j] += theData[k]*nf->theData[i];
                    }
                }
            }
        } else {
            for (long k=0; k<lDim; k++) {
                long l=theIndex[k];
                if (l>=0) {
                    long i = l%vDim;
                    long j = l/vDim;

                    if (i!=j && (!stencil || stencil->theData[l] > 0.0)) {
                        diags[j] += theData[k];
                    }
                }
            }
        }
        {
            for (long k=0; k<hDim; k++) {
                result += diags[k]*nf->theData[k];
            }
        }
        delete [] diags;
    } else {
        if (mbf) {
            for (long k=0; k<hDim; k++) {
                long i,j;

                _Parameter diag = 0.;
                for (j=k*vDim,i=0; j<k*(vDim+1); j++,i++)
                    if (!stencil || stencil->theData[j] > 0.0) {
                        diag += theData[j]*nf->theData[i];
                    }

                j++;
                i++;
                for (; j<vDim*(k+1); j++,i++)
                    if (!stencil || stencil->theData[j] > 0.0) {
                        diag += theData[j]*nf->theData[i];
                    }

                result += diag * nf->theData[k];
            }
        } else {
            for (long k=0; k<hDim; k++) {
                long j;
                _Parameter diag = 0.;

                for (j=k*vDim; j<k*(vDim+1); j++)
                    if (!stencil || stencil->theData[j] > 0.0) {
                        diag += theData[j];
                    }

                j++;
                for (; j<vDim*(k+1); j++)
                    if (!stencil || stencil->theData[j] > 0.0) {
                        diag += theData[j];
                    }

                result += diag * nf->theData[k];
            }
        }
    }

    if (nf!=freqs) {
        DeleteObject (nf);
    }

    return result;
}

//_____________________________________________________________________________________________
_List*      _Matrix::ComputeRowAndColSums (void)
// the first entry is the matrix with row sums
// the second - the entry with column sums
// the third  - a constant with the total sum
{
    if ((storageType == 1) && (hDim >= 1) && (vDim >= 1)) {
        _List*      resList = new _List;
        _Matrix     *rowSums     = new _Matrix (hDim,1,false,true),
        *columnSums  = new _Matrix (vDim,1,false,true);

        if (!(rowSums && columnSums && resList)) {
            checkPointer (nil);
        }

        _Parameter totals = 0.0;

        if (theIndex) {
            for (long item = 0; item < lDim; item ++) {
                long idx = theIndex[item];
                if (idx>=0) {

                    _Parameter      v = theData[idx];

                    rowSums->theData[idx/vDim] += v;
                    columnSums->theData[idx%vDim] += v;
                    totals += v;
                }
            }
        } else {
            for (long rows = 0; rows < hDim; rows++) {
                _Parameter rowSum = 0.;

                for (long columns = 0; columns < vDim; columns ++) {
                    rowSum += theData[rows*vDim+columns];
                }

                rowSums->theData[rows] = rowSum;
                totals += rowSum;
            }

            for (long columns = 0; columns < vDim; columns++) {
                _Parameter colSum = 0.;

                for (long rows = 0; rows < hDim; rows ++) {
                    colSum += theData[rows*vDim+columns];
                }

                columnSums->theData[columns] = colSum;
            }
        }

        (*resList) << rowSums;
        (*resList) << columnSums;

        DeleteObject (rowSums);
        DeleteObject (columnSums);

        {
            _Constant tempC (totals);
            (*resList) && & tempC;
        }

        return resList;

    }
    return nil;
}

//_____________________________________________________________________________________________

_Matrix* _Matrix::NeighborJoin (bool methodIndex)
{
    long          specCount = GetHDim();

    if (storageType != 1 ||  specCount!= GetVDim() || specCount < 4) {
        WarnError ("NeigborJoin needs a square numeric matrix of dimension >= 4");
        return    new _Matrix;
    }

    CheckIfSparseEnough (true);

    _Matrix              netDivergence (specCount,1,false,true);
    _SimpleList          useColumn     (specCount,0,1),
                         columnIndex   (specCount,0,1);

    _Matrix*             res = new _Matrix         ((specCount+1)*2,3,false,true);

    checkPointer         (res);


    for (long k=0; k<specCount ; k=k+1) {
        for (long j=0; j<k; j=j+1) {
            _Parameter d = theData[j*specCount+k];

            netDivergence.theData[k] += d;
            netDivergence.theData[j] += d;

        }
        res->theData[k*3+2] = 1;
    }

    long   cladesMade = 1;

    while (cladesMade < specCount) {
        _Parameter      min = 1.e100;

        long            minIndex  = -1,
                        minIndex2 = -1,
                        minIndexR = -1,
                        minIndexC = -1,
                        k = specCount-1-cladesMade;

        _Parameter      recRemaining = 1./k;

        if (cladesMade == specCount-1) {
            minIndex = useColumn.lData[1];

            _Parameter d = theData[minIndex];

            if ((d<0)&&methodIndex) {
                d = 0;
            }

            k = columnIndex.lData[1];

            if (k>=specCount+cladesMade-2) {
                k = columnIndex[0];
            }

            long    m = specCount+cladesMade-2;

            res->theData[k*3+1]  = d;
            res->theData[k*3]    = m;
            res->theData[3*m+2] += res->theData[3*k+2];
            res->theData[3*m]    = -1;

            break;
        }

        for (long i=1; i<useColumn.lLength; i=i+1) {
            long c1 = useColumn.lData[i];

            for (long j=0; j<i; j=j+1) {
                long c2 = useColumn.lData[j];

                //if (c2>=c1)
                //break;

                _Parameter d = theData[c2*specCount+c1]-(netDivergence.theData[c1]+netDivergence.theData[c2])*recRemaining;

                if (d<min) {
                    min         = d;
                    minIndex    = c2;
                    minIndex2   = c1;
                    minIndexR   = j;
                    minIndexC   = i;
                }
            }
        }

        if (minIndex < 0 || minIndex2 < 0 || minIndexR < 0 || minIndexC < 0) {
            _String err ("Invalid distance matrix passed to NeighborJoin. Matrices written onto messages.log"),
                    invalidMx ((_String*)toStr());
            ReportWarning (invalidMx);
            ReportWarning (_String((_String*)netDivergence.toStr()));
            ReportWarning (_String((_String*)useColumn.toStr()));
            WarnError (err);
            DeleteObject (res);
            return new _Matrix;
        }

        _Parameter      D  = theData[minIndex*specCount+minIndex2],
                        d  = (D - (netDivergence.theData[minIndex2]-netDivergence.theData[minIndex])*recRemaining)*0.5,
                        d2 = D - d;

        if (methodIndex) {
            if (d<0) {
                d = 0.0;
                d2 = D;
            }
            if (d2<0) {
                d2 = 0.0;
                d = D;
                if (d<0) {
                    d = 0;
                }
            }
        }

        long    m = columnIndex.lData [minIndexC],
                n = columnIndex.lData [minIndexR];

        k       = specCount+cladesMade-1;

        res->theData[n*3]       =   k;
        res->theData[n*3+1]     =   d;

        res->theData[m*3]       =   k;
        res->theData[m*3+1]     =   d2;

        res->theData[k*3+2] = res->theData[n*3+2]+res->theData[m*3+2]+1;

        d = theData[minIndex*specCount+minIndex2];

        netDivergence.theData[minIndex]  = 0;
        netDivergence.theData[minIndex2] = 0;

        useColumn.Delete(minIndexC);
        columnIndex.Delete(minIndexC);

        for (k=0; k<useColumn.lLength; k++) {
            long  k2 = useColumn.lData[k];

            if (k2>=minIndex) {
                if (k2 == minIndex) {
                    k++;
                }
                break;
            }

            _Parameter d2 = theData[k2*specCount+minIndex]+theData[k2*specCount+minIndex2],
                       t  =  (d2-d)*.5;

            netDivergence.theData  [k2]               += t-d2;
            theData [k2*specCount+minIndex]            = t;
            netDivergence.theData[minIndex]           += t;

        }

        for (; k<useColumn.lLength; k++) {
            long  k2 = useColumn.lData[k];
            if (k2 >= minIndex2) {
                if (k2 == minIndex2) {
                    k++;
                }
                break;
            }

            _Parameter  d2 = theData[minIndex*specCount+k2]+theData[k2*specCount+minIndex2],
                        t =  (d2-d)*.5;

            netDivergence.theData [k2]                  += t-d2;
            theData[minIndex*specCount+k2]               = t;
            netDivergence.theData[minIndex]             += t;

        }

        //for (k=minIndex2+1;k<ds.species; k=k+1)
        for (; k<useColumn.lLength; k++) {
            long  k2 = useColumn.lData[k];

            _Parameter  d2 = theData[minIndex*specCount+k2]+theData[minIndex2*specCount+k2],
                        t =  (d2-d)*.5;

            netDivergence.theData [k2]                   += t-d2;
            theData[minIndex*specCount+k2]                = t;
            netDivergence.theData[minIndex]              += t;
        }

        columnIndex.lData[minIndexR] = specCount+cladesMade-1;
        {
            for (long i=0; i<minIndex2; i++) {
                theData[i*specCount+minIndex2] = 0;
            }
        }
        {
            for (long i=minIndex2+1; i<specCount; i++) {
                theData[minIndex2*specCount+i]=0;
            }
        }

        cladesMade ++;
    }


    //_Matrix    *tree  = res->MakeTreeFromParent (specCount);
    //DeleteObject (res);
    //return tree;
    return res;
}

//_____________________________________________________________________________________________
_Matrix*        _Matrix::MakeTreeFromParent (long specCount)
{
    if (hDim == 0 || vDim == 0) {
        return new _Matrix;
    }

    _Matrix     *tree = new _Matrix (2*(specCount+1),5,false,true),
    CI  (2*(specCount+1),1,false,true);

    checkPointer (tree);

    for (long kk = 0; kk < specCount-1; kk++) {
        tree->theData[kk*5+4] = -1;
    }

    long cladesMade = 0;

    for (long nodeID2 = 0; nodeID2 < specCount; nodeID2 ++) {
        long        nodeID       = nodeID2,
                    nodeDepth    = 0,
                    saveNodeID   = nodeID,
                    parentID     = theData[nodeID*3],
                    layoutOffset = cladesMade,
                    m,
                    n;

        while (parentID>=0) {
            n = tree->theData[(parentID-specCount)*5+4];
            if (n >= 0) {
                layoutOffset = n+tree->theData[(parentID-specCount)*5+3];
                break;
            }
            parentID  = theData[parentID*3];
        }

        parentID   = theData[nodeID*3];

        while (parentID>=0) {
            n = parentID-specCount;
            m = theData[nodeID*3+2];

            if (tree->theData[n*5+4] < 0)
                /* this node hasn't been laid out yet */
            {
                if (theData[parentID*3]>=0) {
                    tree->theData[n*5+4] = layoutOffset; /* where the layout for the clade begins */
                    tree->theData[n*5+3]   = m; /* offset for that layout */
                }

                m += layoutOffset - 1;

                tree->theData[m*5]   = nodeID;
                tree->theData[m*5+2] = theData[nodeID*3+1];

                CI.theData[nodeID] = m;
            } else
                /* it has been laid out */
            {
                m += tree->theData[n*5+3]+tree->theData[n*5+4] - 1;

                tree->theData[m*5]   = nodeID;
                tree->theData[m*5+2] = theData[nodeID*3+1];

                tree->theData[n*5+3] = m + theData[nodeID*3+2];

                CI.theData[nodeID]   = m;
                nodeDepth ++;

                break;
            }
            nodeDepth++;
            nodeID    = parentID;
            parentID  = theData[nodeID*3];
        }

        /* update levels of nodes */

        if (parentID<0) {
            nodeID   = saveNodeID;
            parentID = theData[nodeID*3];

            while (parentID>=0) {
                m = CI.theData[nodeID];
                tree->theData[m*5+1] = nodeDepth;
                nodeDepth --;
                saveNodeID = nodeID;
                nodeID     = parentID;
                parentID   = theData[nodeID*3];
            }

            cladesMade += theData[3*saveNodeID+2];
        } else {
            m = CI.theData[parentID];

            n = tree->theData[m*5+1];/* depth of the parent */

            nodeID   = saveNodeID;

            while (nodeDepth >= 0) {
                m = CI.theData[nodeID];

                tree->theData[m*5+1] = nodeDepth+n;

                nodeDepth --;
                nodeID  = theData[nodeID*3];
            }
        }
    }
    tree->theData[cladesMade*5]      = 2*specCount-2;
    tree->theData[cladesMade*5+1]    = 0;
    tree->theData[(specCount-2)*5+4] = 0;
    return tree;
}


//_____________________________________________________________________________________________
_Parameter      _Matrix::FisherExact (_Parameter p1, _Parameter p2, _Parameter p3)
{
    if ((hDim>=1)&&(vDim>=1)&&(hDim+vDim>2)) {
        if (vDim<hDim) {
            _Matrix temp (*this);
            temp.Transpose();
            return  temp.FisherExact (p1,p2,p3);
        }
        _Matrix *  numericMx = (_Matrix*)ComputeNumeric();

        double     prob,
                   pval;

        numericMx->CheckIfSparseEnough (true);

        double        *tempArray = new double [numericMx->lDim];
        checkPointer (tempArray);

        for (long i=0; i<hDim; i++)
            for (long j=0; j<vDim; j++) {
                tempArray[j*hDim+i] = numericMx->theData[i*vDim+j];
            }

        fexact_ (hDim,vDim,tempArray,p1,p2,p3,&prob,&pval);
        delete  []  tempArray;
        return pval;

    }
    return 1.;
}

//_____________________________________________________________________________________________

void        _Matrix::SimplexHelper1 (long rowIndex, _SimpleList& columnList, long columnCount, bool useAbsValue, long& maxIndex, _Parameter& maxValue)
// find the maximum element (using absolute value of not) in row rowIndex+1 of this matrix,
// over first columnCount columns indexed by columnList
//  column indexing is offset by + 1 to account for the first column not being eligible for pivoting
{
    if (columnCount <= 0) {
        maxValue = 0.0;
    } else {
        rowIndex = (rowIndex+1)*vDim;
        maxIndex = columnList.lData[0];
        maxValue = theData[rowIndex+maxIndex+1];
        for (long k=1; k<columnCount; k++) {
            _Parameter t = useAbsValue?
                           (fabs(theData[rowIndex+columnList.lData[k]+1])-fabs(maxValue))
                           :(theData[rowIndex+columnList.lData[k]+1]-maxValue);
            if (t>0.) {
                maxValue = theData[rowIndex+columnList.lData[k]+1];
                maxIndex = columnList.lData[k];
            }
        }
    }
}

//_____________________________________________________________________________________________

void        _Matrix::SimplexHelper2 (long& pivotIndex, long columnToExamine, _Parameter eps)
{
    long            m = hDim-2,
                    n = vDim-1,
                    i = 0;

    _Parameter      q1,
                    q;

    pivotIndex = -1;
    for (; i<m; i++)
        if (theData[(i+1)*vDim+columnToExamine+1] < -eps) {
            break;
        }
    if (i>=m) {
        return;    // function is unbounded
    }
    q1              = -theData[(i+1)*vDim]/theData[(i+1)*vDim+columnToExamine+1];
    pivotIndex      = i;
    for (i=pivotIndex+1; i<m; i++) {
        if (theData[(i+1)*vDim+columnToExamine+1] < -eps) {
            q = -theData[(i+1)*vDim]/theData[(i+1)*vDim+columnToExamine+1];
            if (q<q1) {
                pivotIndex = i;
                q1         = q;
            } else {
                _Parameter q0, qp;
                if (q==q1) { // degeneracy
                    for (long k=0; k<n; k++) {
                        qp = -theData[(pivotIndex+1)*vDim + k + 1]/theData[(pivotIndex+1)*vDim+columnToExamine+1];
                        q0 = -theData[(i+1)*vDim + k + 1]/theData[(i+1)*vDim+columnToExamine+1];
                        if (q0!=qp) {
                            break;
                        }
                    }
                    if (q0 < qp) {
                        pivotIndex = i;
                    }
                }
            }
        }
    }

}

//_____________________________________________________________________________________________

void        _Matrix::SimplexHelper3 (long i1, long k1, long ip, long kp)
{
    _Parameter piv = 1./theData[(ip+1)*vDim+kp+1];
    for (long i=0; i<=i1+1; i++)
        if (i-1 != ip) { // not the pivot row
            theData[i*vDim+kp+1] *= piv;
            for (long k=0; k<=k1+1; k++)
                if (k-1 != kp) {
                    theData[i*vDim+k] -= theData[(ip+1)*vDim+k] * theData[i*vDim+kp+1];
                }
        }
    for (long k=0; k<=k1+1; k++)
        if (k-1 != kp)  {
            theData[(ip+1)*vDim+k] *= -piv;
        }
    theData[(ip+1)*vDim+kp+1] = piv;
}

//_____________________________________________________________________________________________
_Matrix*    _Matrix::SimplexSolve (_Parameter desiredPrecision )
// this function is adapted from the Num. Recipes in C version; but with 0 indexing
// hyphy primitives
// and without goto labels

// the of dimension RxC is interpreted as follows
// R-1 constaints
// C-2 variables

// the first row:
//      cell   0      - current value of the objective function
//      cells  1-C-2  - coefficient of variable x_k in the objective function
//      cell   C-1    - if >=0. then maximize the function
//                    - if <0   then minimize the function

// other rows
// constraint j written in the form
// b_j - a_j1 x_1 - a_j2 x_2 - ... -a_jk x_k
// last cell is what type of constraint it is:
// < 0: <= inequality
// > 0: >= inequality
// = 0 equality

// upon return, will contain a row matrix of either C-1 cells:
// extreme value of the objective function in the first cell
// variable values in the same order as originally supplied
// if an empty matrix is returned - no feasible solution could be found
// if a 1x1 matrix is returned - the objective function is unbounded

{
    _String errMsg;

    long n = vDim-2, // number of variables
         m = hDim-1; // number of constraints

    if (storageType == 1 && n>0 && m>0) {
        while (1) // artificial construct used to break out on error
            // an to avoid goto statements
        {
            bool        doMaximize = (*this)(0,n+1) >= 0.0;

            // allocate temporary storage
            _Matrix     tempMatrix (m+2,n+1,false,true);
            // first, copy the objective function row
            for (long i=0; i<=n; i++) {
                tempMatrix.Store(0,i,doMaximize?(*this)(0,i):-(*this)(0,i));
            }

            // now, count the number of constraints of each type and reorder things

            long    m1 = 0, // <= constraints
                    m2 = 0, // >= constraints
                    m3 = 0; // == constraints

            {
                for (long i=1; i<=m; i++) {
                    _Parameter t = (*this)(i,n+1);
                    if (t<0.0) {
                        m1++;
                    } else if (t>0.0) {
                        m2++;
                    } else {
                        m3++;
                    }
                    if ((*this)(i,0) < 0.0) {
                        errMsg = "Negative values are not allowed in the first column of the simplex tableau";
                        break;
                    }
                }
            }

            if (errMsg.sLength) {
                break;
            }

            // copy coefficients into the temp matrix, sorting the constraints in the <=, >= and == order
            {
                for (long i=1, t1=0, t2=0, t3=0; i<=m; i++) {
                    _Parameter t   = (*this)(i,n+1);
                    long       idx;
                    if (t<0.0) {
                        idx=1+t1++;
                    } else if (t>0.0) {
                        idx=1+m1+t2++;
                    } else {
                        idx=1+m1+m2+t3++;
                    }
                    for (long j=0; j<=n; j++) {
                        tempMatrix.Store(idx,j,(*this)(i,j));
                    }

                }
            }
            // allocate temporary storage

            _SimpleList l1      (n+1,0,1),
                        l3       (m,0,0),
                        izrov    (n,0,1),
                        iposv    (m,n,1);

            long        nl1     = n;

            if (m2+m3) { // >= and == constraints exist; origin is not a feasible solution
                for (long i=0; i<m2; l3.lData[i] = 1,i++) ; // slack variables in the list of 'basis' variables
                for (long k=0; k<=n; k++) { // compute the auxiliary objective function
                    _Parameter q = 0.;
                    for (long k2 = m1+1; k2<=m; k2++) {
                        q += tempMatrix(k2,k);
                    }
                    tempMatrix.Store (m+1,k,-q);
                }
                while (1) { // initial artifical construct
                    long        pivotColumn,
                                ip;
                    _Parameter  pivotValue;

                    tempMatrix.SimplexHelper1 (m,l1,nl1,false,pivotColumn, pivotValue);
                    if (pivotValue <= desiredPrecision && tempMatrix(m+1,0) < -desiredPrecision)
                        // aux objective is still negative and can't be improved
                        // no feasible solution
                    {
                        return new _Matrix;
                    }
                    if (pivotValue <= desiredPrecision && tempMatrix(m+1,0) <= desiredPrecision)
                        // aux objective is zero and can't be improved
                        // found a feasible solution; clean up artificial variables and move on to phase 2
                    {
                        for (ip = m1+m2; ip < m; ip++) {
                            if (iposv.lData[ip] == ip + n) {
                                tempMatrix.SimplexHelper1 (ip,l1,nl1,true,pivotColumn, pivotValue);
                                if (pivotValue > desiredPrecision) {
                                    goto one;
                                }
                            }
                        }
                        for (long i = m1; i<m1+m2; i++)
                            if (l3.lData[i-m1] == 1)
                                for (long k=0; k<=n; k++) {
                                    tempMatrix.Store (i+1,k,-tempMatrix(i+1,k));
                                }

                        break;
                    }

                    tempMatrix.SimplexHelper2 (ip,pivotColumn,desiredPrecision);
                    if (ip<0) {
                        return new _Matrix (1,1,false,true);    // unbounded function
                    }

one:
                    tempMatrix.SimplexHelper3 (m,n-1,ip,pivotColumn);
                    if (iposv.lData[ip] >= n+m1+m2) {
                        long k = 0;
                        for (k=0; k<nl1; k++)
                            if (l1.lData[k] == pivotColumn) {
                                break;
                            }
                        nl1--;
                        for (long i2=k; i2<nl1; i2++) {
                            l1.lData[i2] = l1.lData[i2+1];
                        }
                    } else {
                        long k2 = iposv.lData[ip] - m1 - n;
                        if (k2 >= 0 && l3.lData[k2]) {
                            l3.lData[k2] = 0;
                            tempMatrix.theData[(m+1)*tempMatrix.vDim + pivotColumn + 1] ++;
                            for (long i=0; i<m+2; i++) {
                                tempMatrix.theData[i*tempMatrix.vDim + pivotColumn + 1] *= -1.0;
                            }
                        }
                    }
                    long s = izrov.lData[pivotColumn];
                    izrov.lData[pivotColumn] = iposv.lData[ip];
                    iposv.lData[ip] = s;
                }// end of phase 1
            }

            while (1) {
                long            pivotColumn,
                                pivotRow;
                _Parameter      pivotValue;

                tempMatrix.SimplexHelper1 (-1,l1,nl1,false,pivotColumn,pivotValue);
                if (pivotValue < desiredPrecision) { // done!
                    // produce the final solution
                    _Matrix * resMatrix = new _Matrix (1,n+1,false,true);
                    checkPointer (resMatrix);
                    resMatrix->Store(0,0,doMaximize?tempMatrix(0,0):-tempMatrix(0,0));
                    for (long k=0; k<iposv.lLength; k++)
                        if (iposv.lData[k]<n) {
                            resMatrix->Store(0,iposv.lData[k]+1,tempMatrix(k+1,0));
                        }
                    return resMatrix;
                }
                tempMatrix.SimplexHelper2 (pivotRow,pivotColumn,desiredPrecision);
                if (pivotRow<0) {
                    return new _Matrix (1,1,false,true);
                }
                tempMatrix.SimplexHelper3 (m-1,n-1,pivotRow,pivotColumn);
                long s = izrov.lData[pivotColumn];
                izrov.lData[pivotColumn] = iposv.lData[pivotRow];
                iposv.lData[pivotRow] = s;
            }

        }
    } else {
        errMsg = "SimplexSolve requires a numeric matrix with > 1 row and > 2 columns";
    }

    WarnError (errMsg);
    return new _Matrix;
}

//_____________________________________________________________________________________________

void    _Matrix::CopyABlock (_Matrix * source, long startRow, long startColumn, long rowSpan, long colSpan)
{
    long indexTarget = startRow*vDim + startColumn,
         indexSource = 0,
         sourceHDim  = rowSpan<=0?source->hDim:rowSpan,
         sourceVDim  = colSpan<=0?source->vDim:colSpan,
         maxRow         = MIN (hDim, startRow    + sourceHDim),
         maxColumn   = MIN (vDim, startColumn + sourceVDim);

    for  (long r = startRow; r < maxRow; r++) {
        for (long c = startColumn, c2 = 0; c < maxColumn; c++, c2++) {
            theData[indexTarget+c2] = source->theData[indexSource+c2];
        }

        indexSource += sourceVDim;
        indexTarget += vDim;
    }
}


//_____________________________________________________________________________________________
_PMathObj   _Matrix::DirichletDeviate (void)
{
    /* -----------------------------------------------------------
        DirichletDeviate()
            Generate vector of random deviates from the Dirichlet
            distribution defined by contents of this matrix as
            hyperparameters (a > 0).
       ----------------------------------------------------------- */

    _String     errMsg;

    long        dim;

    _Parameter  denom   = 0.;

    _Matrix     res (1, dim = GetHDim()*GetVDim(), false, true);    // row vector


    if (storageType != 1) {
        errMsg = "Only numeric vectors can be passed to <= (DirichletDeviate)";
    }

    if (IsAVector()) {
        // generate a random deviate from gamma distribution for each hyperparameter
        for (long i = 0; i < dim; i++) {
            if (theData[i] < 0) {
                WarnError (_String("Dirichlet not defined for negative parameter values."));
                return new _Matrix (1,1,false,true);
            }

            res.Store (0, i, gammaDeviate(theData[i]));
            denom += res(0,i);
        }

        // normalize by sum
        for (long i = 0; i < dim; i++) {
            res.Store (0, i, res(0,i)/denom);
        }

        return (_PMathObj) res.makeDynamic();
    } else {
        errMsg = "Argument must be a row- or column-vector.";
    }

    WarnError (errMsg);
    return new _Matrix (1,1,false,true);
}



//_____________________________________________________________________________________________
_PMathObj   _Matrix::GaussianDeviate (_Matrix & cov)
{
    /* ------------------------------------------------------
        GaussianDeviate()
            Generate vector of random deviates from k-
            dimensional Gaussian distribution given contents
            of this matrix as mean parameters, and argument
            as covariance matrix.

            Use algorithm described in Numerical Recipes
            3rd ed., p.379
       ------------------------------------------------------ */

    ReportWarning (_String("Entered _Matrix::GaussianDeviate() with cov = ") & (_String *)(cov.toStr()));

    _String     errMsg;

    if (storageType != 1 || GetHDim() > 1) {
        WarnError (_String("ERROR in _Matrix::GaussianDeviate(), expecting to be called on numeric row vector matrix, current dimensions: ") & GetHDim() & "x" & GetVDim());
        return new _Matrix;
    }

    long        kdim        = GetVDim();    // number of entries in this _Matrix object as vector of means

    if (cov.GetHDim() == kdim && cov.GetVDim() == kdim) {
        _Matrix     * cov_cd    = (_Matrix *) cov.CholeskyDecompose();
        _Matrix     gaussvec (1, kdim, false, true);

        ReportWarning (_String("\nCholesky decomposition of cov = ") & (_String *) cov_cd->toStr());

        // fill column vector with independent standard normal deviates
        for (long i = 0; i < kdim; i++) {
            gaussvec.Store (0, i, gaussDeviate());
        }

        ReportWarning (_String ("\nvector of gaussian deviates = ") & (_String *) gaussvec.toStr());

        // left multiply vector by Cholesky decomposition of covariance matrix
        gaussvec *= (_Matrix &) (*cov_cd);

        // shift mean
        for (long i = 0; i < kdim; i++) {
            gaussvec.Store (0, i, gaussvec(0,i) + theData[i]);
        }

        DeleteObject (cov_cd);
        return (_PMathObj) gaussvec.makeDynamic();
    }

    WarnError (_String("Error in _Matrix::GaussianDeviate(), incompatible dimensions in covariance matrix: ") & cov.GetHDim() & "x" & cov.GetVDim());
    return new _Matrix;
}


//_____________________________________________________________________________________________
_PMathObj   _Matrix::MultinomialSample (_Constant *replicates)
{
    _String       errMsg;
    long          values      = GetHDim();
    unsigned long samples     = replicates?replicates->Value ():0;

    _Matrix     *eval    = (_Matrix*)Compute (),
                 * sorted = nil,
                   * result = nil;

    if (samples < 1) {
        errMsg = "Expected a numerical (>=1) value for the number of replicates";
    } else if (eval->storageType != 1 || GetVDim() != 2 || values < 2) {
        errMsg = "Expecting numerical Nx2 (with N>=1) matrix.";
    } else {
        _Constant one (1.);
        sorted = (_Matrix*) eval->SortMatrixOnColumn(&one);


        _Parameter      sum = 0.;

        for (long n = 1; n < 2*values; n+=2) {
            _Parameter v = sorted->theData[n];
            if (v < 0.) {
                sum = 0.;
                break;
            }
            sum += v;
        }



        if (CheckEqual (sum, 0.)) {
            errMsg = "The probabilities (second column) cannot add to 0 or be negative";
        } else {
            sum = 1./sum;

            _Matrix     *raw_result  = new _Matrix (1, values, false, true),
            *normalized  = new _Matrix (1, values, false, true);

            for (long v = 0; v < values; v++) {
                normalized->theData[values-1-v] = sorted->theData[1+2*v] * sum;
            }

            //BufferToConsole (_String(*(_String*)normalized->toStr()).sData);
            //NLToConsole ();

            _String      _HYMultinomialStatus       ("Generating multinomial samples");


#if !defined __UNIX__ || defined __HEADLESS__
            TimerDifferenceFunction(false); // save initial timer; will only update every 1 second
#if !defined __HEADLESS__
            SetStatusLine     (empty,_HYMultinomialStatus, empty, 0, HY_SL_TASK|HY_SL_PERCENT);
#else
            SetStatusLine     (_HYMultinomialStatus);
#endif
            _Parameter  seconds_accumulator = .0,
                        temp;
#endif

            for (unsigned long it = 0; it < samples; it++) {
                _Parameter randomValue = genrand_real2(),
                           sum   = normalized->theData[0];
                long       index = 0;

                while (sum < randomValue) {
                    index++;
                    sum += normalized->theData[index];
                }

                raw_result->theData[index] += 1.;
#if !defined __UNIX__ || defined __HEADLESS__
                if ((it % 1000 == 0) && (temp=TimerDifferenceFunction(true))>1.0) { // time to update
                    seconds_accumulator += temp;

                    _String statusLine = _HYMultinomialStatus & " " & (_Parameter)(it+1) & "/" & (_Parameter)samples
                                         & " samples drawn (" & (1.0+it)/seconds_accumulator & "/second)";

#if defined __HEADLESS__
                    SetStatusLine (statusLine);
#else
                    SetStatusLine (empty,statusLine,empty,100*(float)it/(samples),HY_SL_TASK|HY_SL_PERCENT);
#endif
                    TimerDifferenceFunction (false); // reset timer for the next second
                    yieldCPUTime (); // let the GUI handle user actions

                    if (terminateExecution) { // user wants to cancel the analysis
                        break;
                    }
                }
#endif
            }

            result = new _Matrix (1, values, false, true);

            for (long v = 0; v < values; v++) {
                result->theData[(long)sorted->theData[2*(values-1-v)]] = raw_result->theData[v];
            }

            DeleteObject (raw_result);
            DeleteObject (sorted);
            sorted = normalized;
        }
    }


    DeleteObject (sorted);
    if (errMsg.sLength) {
        WarnError (_String("Error in _Matrix::MultinomialSample(). ") & errMsg);
        DeleteObject (result);
        return new _Matrix;
    }
    return result;
}


//_____________________________________________________________________________________________
_PMathObj   _Matrix::InverseWishartDeviate (_Matrix & df)
{
    /* ---------------------------------------------------
        InverseWishartDeviate()
            Generates a random matrix whose inverse
            has the Wishart distribution with this matrix
            supplying the covariance matrix parameter and
            a degrees of freedom vector argument.
       --------------------------------------------------- */

    _String     errMsg;
    long        n       = GetHDim();


    if (storageType != 1 || GetHDim() != GetVDim()) {
        errMsg = "expecting numerical symmetric matrix.";
    }

    else if (df.MatrixType() != 1 || df.GetHDim() != n || df.GetVDim() > 1) {
        errMsg = "expecting numerical row vector for second argument (degrees of freedom).";
    } else {
        // compute Cholesky factor for this matrix inverse, extract the diagonal
        _Matrix * inv       = (_Matrix *) Inverse();
        _Matrix * invCD     = (_Matrix *) (inv->CholeskyDecompose());

        _Matrix decomp ((_Matrix &) *invCD);    // duplication constructor

        DeleteObject (invCD);

        return WishartDeviate (df, decomp);
    }



    WarnError (_String ("ERROR in _Matrix::InverseWishartDeviate, ") & errMsg);
    return new _Matrix;
}

//_____________________________________________________________________________________________
_PMathObj   _Matrix::WishartDeviate (_Matrix & df)
{
    _Matrix     diag;   // calls default constructor
    return WishartDeviate (df, diag);
}


_PMathObj   _Matrix::WishartDeviate (_Matrix & df, _Matrix & decomp)
{
    /* ---------------------------------------------------
     WishartDeviate()
        Generates a random matrix following the Wishart
        distribution with this matrix supplying the
        covariance matrix parameter.

        First argument: degrees of freedom vector.
        Second argument (optional):
            Diagonal of Cholesky decomposition of
            covariance matrix, overrides this matrix.
     --------------------------------------------------- */


    // debugging
    ReportWarning (_String("Entered _Matrix::WishartDeviate() with this matrix: ") & (_String *) this->toStr() & " and df vector " & (_String *) df.toStr());


    long        n   = GetHDim();

    _Matrix     rdeviates (n, n, false, true),
                rd_transpose;


    if (!df.IsAVector(0)) {
        WarnError (_String ("ERROR in _Matrix::WishartDeviate(), expecting row vector for degrees of freedom argument."));
        return new _Matrix (1,1,false,true);
    } else if (df.IsAVector(1)) {
        df.Transpose(); // convert column vector to row vector
    }


    if (decomp.GetHDim() == 0) {    // no second argument, perform Cholesky decomposition
        if (storageType != 1 || GetHDim() != GetVDim()) {
            WarnError (_String ("ERROR in _Matrix::WishartDeviate(), expecting square numeric matrix."));
            return new _Matrix (1,1,false,true);
        } else {
            _Matrix     * cholesky = (_Matrix *) CholeskyDecompose();

            if (cholesky->GetHDim() > 0) {
                // duplicate
                CreateMatrix (&decomp, cholesky->GetHDim(), cholesky->GetVDim(), false, true, false);

                for (long i = 0; i < cholesky->GetHDim(); i++) {
                    for (long j = 0; j < cholesky->GetVDim(); j++) {
                        decomp.Store(i, j, (*cholesky)(i,j));
                    }
                }

                DeleteObject (cholesky);
            } else {
                return (cholesky);  // empty _Matrix from error in CholeskyDecompose()
            }
        }
    }

    ReportWarning (_String("diag=") & (_String *)decomp.toStr());   // column vector


    // populate diagonal with square root of i.i.d. chi-square random deviates
    for (long i = 0; i < n; i++) {
        rdeviates.Store (i, i, sqrt(chisqDeviate(df(0,i)-i+1)) );

        // populate upper triagonal with i.i.d. standard normal N(0,1) deviates
        for (long j = i+1; j < n; j++) {
            rdeviates.Store (i, j, gaussDeviate());
        }
    }

    ReportWarning (_String("rdeviates(A)=") & (_String *)rdeviates.toStr());


    // result is obtained from D^T B D, where B = A^T A, ^T is matrix transpose
    rd_transpose = (_Matrix &) rdeviates;
    rd_transpose.Transpose();
    ReportWarning (_String("transpose(A)=") & (_String *)rd_transpose.toStr());
    rd_transpose *= (_Matrix &) rdeviates;  // A^T A
    ReportWarning (_String("A^T A=") & (_String *)rd_transpose.toStr());
    rd_transpose *= (_Matrix &) decomp; // A^T A D
    ReportWarning (_String("A^T A D=") & (_String *)rd_transpose.toStr());

    decomp.Transpose();
    decomp *= (_Matrix &) rd_transpose; // D^T A^T A D
    ReportWarning (_String("D^T A^T A D=") & (_String *)decomp.toStr());

    return (_PMathObj) decomp.makeDynamic();
}




//_____________________________________________________________________________________________
// AssociativeList
//_____________________________________________________________________________________________

_AssociativeList::_AssociativeList (void):avl(&theData)
{
}

//_____________________________________________________________________________________________

BaseRef _AssociativeList::makeDynamic (void)
{
    _AssociativeList * newAL = new _AssociativeList ();
    newAL->Duplicate (this);
    return newAL;
}

//_____________________________________________________________________________________________

bool _AssociativeList::ParseStringRepresentation (_String& serializedForm, bool doErrors, _VariableContainer* theP)
{
    _List               splitKeys;
    _ElementaryCommand::ExtractConditions (serializedForm, 0, splitKeys, ',' , false);
    for (long k = 0; k < splitKeys.lLength; k = k + 1) {
        _List aPair;
        _ElementaryCommand::ExtractConditions (*(_String*)splitKeys(k), 0, aPair, ':' , false);
        if (aPair.lLength == 2) {
            _String  key        (ProcessLiteralArgument((_String*)aPair(0),theP));
            _Formula value      (*(_String*)aPair(1),theP, doErrors);

            _PMathObj   valueC  = value.Compute();
            if (valueC) {
                MStore (key, valueC, true);
            } else {
                if (doErrors) {
                    WarnError (*(_String*)aPair(1) & " could not be evaluated");
                }
                return false;

            }
        } else {
            if (doErrors) {
                WarnError (*(_String*)splitKeys(k) & " does not appear to specify a valid key:value pair");
            }
            return false;
        }
    }
    return true;
}

//_____________________________________________________________________________________________

BaseRef _AssociativeList::toStr (void)
{
    _String defName   ("_hyphyAssociativeArray");
    return Serialize  (defName);
}

//_____________________________________________________________________________________________

void _AssociativeList::Duplicate (BaseRef br)
{
    nInstances = 1;
    _AssociativeList * copyMe = (_AssociativeList*)br;
    theData.Duplicate (&copyMe->theData);
    avl.leftChild.Duplicate (&copyMe->avl.leftChild);
    avl.rightChild.Duplicate (&copyMe->avl.rightChild);
    avl.balanceFactor.Duplicate (&copyMe->avl.balanceFactor);
    avl.emptySlots.Duplicate (&copyMe->avl.emptySlots);
    avl.xtraD.Duplicate (&copyMe->avl.xtraD);
    avl.root = copyMe->avl.root;
}

//_____________________________________________________________________________________________

_PMathObj _AssociativeList::MCoord (_PMathObj p)
{
    return new _FString ((_String*)p->toStr());
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::MAccess (_PMathObj p)
{
    long        f;

    if (p->ObjectClass() == STRING) {
        f = avl.Find (((_FString*)p)->theString);
    } else {
        _String s ((_String*)p->toStr());
        f = avl.Find (&s);
    }
    if (f>=0) {
        _PMathObj res = (_PMathObj)avl.GetXtra (f);
        res->nInstances++;
        return res;
    } else {
        return new _Constant (0.0);
    }
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::MIterator (_PMathObj p, _PMathObj p2)
{
    long done = 0;


    if (p->ObjectClass() == STRING && p2->ObjectClass() == STRING) {

        long avlRoot = avl.GetRoot();

        if (avlRoot >= 0) {
            _String * s  = (_String*)p->toStr(),
                      * s2 = (_String*)p2->toStr();

            long    fID  = FindBFFunctionName (*s),
                    fID2 = FindBFFunctionName (*s2);

            if (fID < 0 || batchLanguageFunctionParameters.lData[fID] != 2) {
                WarnError ("The first argument in an iterator call for Associative Arrays must be a valid identifier of a function taking two arguments (key, value)");
            } else {
                if (fID2 >= 0 && batchLanguageFunctionParameters.lData[fID2] != 1) {
                    WarnError ("The second argument in an iterator call for Associative Arrays must be either empty or a valid identifier of a function taking a single argument");
                }

                _Formula      testFormula,
                              actionFormula;

                actionFormula.GetList().AppendNewInstance(new _Operation());
                actionFormula.GetList().AppendNewInstance(new _Operation());
                actionFormula.GetList().AppendNewInstance(new _Operation(empty,-fID-1));

                if (fID2 >= 0) {
                    testFormula.GetList().AppendNewInstance(new _Operation());
                    testFormula.GetList().AppendNewInstance(new _Operation(empty,-fID2-1));
                }


                _SimpleList  hist;
                long         ls,
                             cn = avl.Traverser (hist,ls,avlRoot);

                _FString * fKey = new _FString;
                while (cn >= 0) {
                    _String* aKey = ((_String**)avl.dataList->lData)[cn];
                    if (aKey) {
                        DeleteObject (fKey->theString);
                        fKey->theString = (_String*)aKey->toStr();
                        if (fID2 >= 0) {
                            ((_Operation**)testFormula.GetList().lData)[0]->SetNumber(fKey);
                            if (CheckEqual(testFormula.Compute()->Value(),0.0)) {
                                cn = avl.Traverser (hist,ls);
                                continue;
                            }
                        }
                        ((_Operation**)actionFormula.GetList().lData)[0]->SetNumber(fKey);
                        ((_Operation**)actionFormula.GetList().lData)[1]->SetNumber((_PMathObj)avl.GetXtra (cn));
                        actionFormula.Compute();
                        done ++;
                    }
                    cn = avl.Traverser (hist,ls);
                }
                DeleteObject (fKey);

                ((_Operation**)actionFormula.GetList().lData)[0]->SetNumber(nil);
                ((_Operation**)actionFormula.GetList().lData)[1]->SetNumber(nil);
                if (fID2 >= 0) {
                    ((_Operation**)testFormula.GetList().lData)[0]->SetNumber(nil);
                }

            }
            DeleteObject (s) ;
            DeleteObject (s2);
        }
    } else if (p->ObjectClass () == STRING && p2->ObjectClass () == NUMBER) {
        _String * s  = (_String*)p->toStr();

        if (s->Equal (&AVL_ITERATOR_ORDER) || s->Equal (&AVL_ITERATOR_ORDER_VALUE)) {
            long index = avl.GetByIndex(p2->Compute()->Value());
            if (index >= 0) {
                return s->Equal (&AVL_ITERATOR_ORDER)? (new _FString(*((_String**)avl.dataList->lData)[index],false)): ((_PMathObj)avl.GetXtra (index)->makeDynamic());
            } else {
                WarnError ("Index out of bounds in call to AVL iterator (by index)");
            }
        }

        DeleteObject (s);
    } else {
        WarnError ("Both arguments must be Strings (or a String Literal and a number) in an iterator call for Associative Arrays");
    }
    return new _Constant (done);
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::GetByKey (_String& key, long objType)
{
    long        f = avl.Find (&key);

    if (f>=0) {
        _PMathObj res = (_PMathObj)avl.GetXtra (f);
        if (res->ObjectClass () == objType) {
            return res;
        }
    }

    return nil;
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::GetByKey (_String& key)
{
    long        f = avl.Find (&key);

    if (f>=0) {
        return (_PMathObj)avl.GetXtra (f);
    }

    return nil;
}

//_____________________________________________________________________________________________
_PMathObj _AssociativeList::GetByKey (long nKey, long objType)
{
    _String key (nKey);
    return GetByKey (key, objType);
}

//_____________________________________________________________________________________________
void _AssociativeList::DeleteByKey (_PMathObj p)
{
    if (p->ObjectClass() == STRING) {
        avl.Delete (((_FString*)p)->theString,true);
    } else {
        if (p->ObjectClass() == ASSOCIATIVE_LIST) {
            _List * keys2remove = ((_AssociativeList*)p)->GetKeys();
            for (long ki = 0; ki < keys2remove->lLength; ki++) {
                avl.Delete ((_String*)(*keys2remove)(ki),true);
            }
        } else {
            _String * s = (_String*)p->toStr();
            avl.Delete (s,true);
            DeleteObject (s);
        }
    }
}


//_____________________________________________________________________________________________
void _AssociativeList::MStore (_PMathObj p, _PMathObj inObject, bool repl, long opCode)
{
    if (!p) {
        return;
    }

    _FString * index = (_FString*)p;
    long       f     = avl.Find (index->theString);

    if (f>=0) { // already exists - replace
        if (opCode == HY_OP_CODE_ADD) {
            _PMathObj newObject = ((_PMathObj)avl.GetXtra(f))->Execute (HY_OP_CODE_ADD,inObject);
            if (repl == false) {
                DeleteObject (inObject);
            } else {
                repl = false;
            }
            inObject = newObject;
        }
        avl.xtraD.Replace (f, inObject, repl);
    } else { // insert new
        if (repl) {
            BaseRef br = inObject->makeDynamic();
            avl.Insert (index->theString->makeDynamic(),(long)br,false);
            //br->nInstances--;
        } else {
            avl.Insert (index->theString->makeDynamic(),(long)inObject,false);
        }
    }
}



//_____________________________________________________________________________________________
void _AssociativeList::MStore (_String obj, _PMathObj inObject, bool repl)
{
    _FString f (obj);
    MStore (&f,inObject, repl);
}

//_____________________________________________________________________________________________
void _AssociativeList::MStore (_String obj, _String info)
{
    _FString inf (info);
    MStore (obj, &inf, true);
}


//_____________________________________________________________________________________________
_String* _AssociativeList::Serialize (_String& avlName)
{
    _String * outString = new _String (1024L,true);
    checkPointer (outString);

    (*outString) << "{";
    bool        doComma = false;
    _List * meKeys = GetKeys();
    for (long k = 0; k < meKeys->lLength; k=k+1) {
        _String   *thisKey  = (_String*)(*meKeys)(k);
        if (thisKey) {
            if (doComma) {
                (*outString) << ',';
                (*outString) << '\n';
            }

            (*outString) << '"';
            outString->EscapeAndAppend(*thisKey, false);
            (*outString) << '"';

            _PMathObj anObject = GetByKey (*thisKey);

            (*outString) << ':';
            if (anObject->ObjectClass() == STRING) {
                (*outString) << '"';
                outString->EscapeAndAppend(_String ((_String*)anObject->toStr()),0);
                (*outString) << '"';
            } else {
                (*outString) << _String ((_String*)anObject->toStr());
            }
            doComma = true;
        }
    }
    (*outString) << "}";
    outString->Finalize ();
    return outString;
}


//__________________________________________________________________________________
_PMathObj _AssociativeList::Compute (void)
{
    return this;
}

//__________________________________________________________________________________
_List* _AssociativeList::GetKeys (void)
{
    return (_List*)avl.dataList;
}

//_____________________________________________________________________________________________
void        _AssociativeList::FillInList (_List& fillMe)
{
    _SimpleList  hist;
    long         ls,
                 cn = avl.Traverser (hist,ls,avl.GetRoot());

    while (cn >= 0) {
        _String* aKey = ((_String**)avl.dataList->lData)[cn];
        if (aKey) {
            fillMe.AppendNewInstance(avl.GetXtra (cn)->toStr());
        }
        cn = avl.Traverser (hist,ls);
    }
}


//_____________________________________________________________________________________________
void        _AssociativeList::Merge (_PMathObj p)
{
    if (p && p->ObjectClass() == ASSOCIATIVE_LIST) {
        _AssociativeList    * rhs = (_AssociativeList*) p;

        _SimpleList  hist;
        long         ls,
                     cn = rhs->avl.Traverser (hist,ls,rhs->avl.GetRoot());

        while (cn >= 0) {
            _String* aKey = new _String(*((_String**)rhs->avl.dataList->lData)[cn]);
            long insAt = 0;
            if ((insAt = avl.Insert (aKey, (long)rhs->avl.GetXtra(cn)->makeDynamic(), false)) < 0) { // already exists
                avl.SetXtra (-insAt-1,new _MathObject(), false);
                DeleteObject (aKey);
            }

            cn = rhs->avl.Traverser (hist,ls);
        }
    } else {
        WarnError ("Associative list merge operation requires an associative list argument.");
    }
}

//_____________________________________________________________________________________________
_PMathObj        _AssociativeList::Sum (void)
{
    _Parameter sum = 0.;
        
    _SimpleList  hist;
    long         ls,
    cn = avl.Traverser (hist,ls,avl.GetRoot());
        
    while (cn >= 0) {
        _PMathObj value = (_PMathObj)avl.GetXtra (cn);
        switch (value->ObjectClass()){
            case NUMBER:
                sum += ((_Constant*)value)->Value();
                break;
            case STRING:
                sum += ((_FString*)value)->theString->toNum();
                break;
            case MATRIX: {
                _Constant * sumOfValue = (_Constant*) ((_Matrix*)value->Compute())->Sum();
                sum += sumOfValue->Value();
                DeleteObject (sumOfValue);
                break;
            }
            case ASSOCIATIVE_LIST: {
                _Constant * sumOfValue = (_Constant*) ((_AssociativeList*)value->Compute())->Sum();
                sum += sumOfValue->Value();
                DeleteObject (sumOfValue);
                break;
            }
        }
        cn = avl.Traverser (hist,ls);
    }
    
    return new _Constant (sum);
}

//__________________________________________________________________________________


_PMathObj _AssociativeList::Execute (long opCode, _PMathObj p, _PMathObj p2)   // execute this operation with the second arg if necessary
{

    switch (opCode) {
    case HY_OP_CODE_ADD: // +
        if (p){
            MStore (_String((long)avl.countitems()), p, true);
            return new _Constant (avl.countitems());
        } else {
            return Sum();
        }

    case HY_OP_CODE_MUL: // merge
        Merge (p);
        return new _Constant (avl.countitems());

    case HY_OP_CODE_SUB:
    case HY_OP_CODE_ABS:
        if (opCode == HY_OP_CODE_SUB) {
            DeleteByKey (p);
        }
        return new _Constant (avl.countitems());
        break;
    case HY_OP_CODE_MACCESS: // MAccess
        if (p2) {
            return MIterator (p,p2);
        } else {
            return MAccess   (p);
        }
        break;
    case HY_OP_CODE_MCOORD: // MCoord
        return MCoord (p);
        break;
    case HY_OP_CODE_ROWS: // Rows - get keys
        if (avl.emptySlots.lLength) {
            _List  dataListCompact;
            for (long k=0; k<avl.dataList->lLength; k++) {
                BaseRef anItem = ((BaseRef*)avl.dataList->lData)[k];
                if (anItem) {
                    dataListCompact << anItem;
                }
            }
            return new _Matrix (dataListCompact);
        } else {
            return new _Matrix (*(_List*)avl.dataList);
        }
        break;
    case HY_OP_CODE_TYPE: // Type
        return Type();
        break;
    }


    WarnNotDefined (this, opCode);
    return nil;

}

/*--------------------------------------------------------------------------------------------------------------------------------*/
// Growing Vector
/*--------------------------------------------------------------------------------------------------------------------------------*/

_GrowingVector::_GrowingVector (bool iscol) : _Matrix (64,1,false,true)
{
    used = 0;
    isColumn = iscol;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

BaseRef     _GrowingVector::makeDynamic (void)
{
    _GrowingVector * result = (_GrowingVector*)checkPointer(new _GrowingVector);
    result->_Matrix::Duplicate (this);
    result->used = used;
    result->vDim = 1;
    result->isColumn = isColumn;
    return result;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

long        _GrowingVector::Store (_Parameter toStore)
{
    if (used < hDim) {
        theData[used++] = toStore;  // increment AFTER argument is sent to function
        return used-1;
    } else {
        Resize (used + MAX (used/8,64));    // allocate another block of 64
        return Store (toStore);
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        _GrowingVector::operator << (const _SimpleList& theSource)
{
    for (long k = 0; k < theSource.lLength; k++) {
        Store (theSource.lData[k]);
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void        _GrowingVector::Clear (void)
{
    _Matrix::Clear();
    ZeroUsed();
    vDim = 1;
}

//_____________________________________________________________________________________________
void _GrowingVector::Duplicate (BaseRef obj)
{
    _Matrix::Duplicate (obj);
    used = ((_GrowingVector*)obj)->used;
    isColumn = ((_GrowingVector*)obj)->isColumn;
}


/*--------------------------------------------------------------------------------------------------------------------------------*/
// _NTupleStorage
/*--------------------------------------------------------------------------------------------------------------------------------*/

_NTupleStorage::_NTupleStorage (unsigned long N, unsigned long K)
{
    storageN        = N;
    storageK        = (K>N)?MIN(N,1):K;
    unsigned long   matrixDimension = 1;

    // now compute what dimension of the matrix will be required
    // handle the special cases first

    if (storageK) { // not just the empty set
        // populate C_NK_Lookup
        for (long i2=0; i2<=storageN; i2++) {
            C_NK_Lookup << 1;    // N choose 0 is always 1 for every N
        }
        for (long i=1; i<=storageK; i++) {
            for (long filler = 0; filler < i; filler++, C_NK_Lookup<<0); // N choose K where K>N is invalid
            C_NK_Lookup << 1;
            // K choose K is 1
            for (long j=i+1; j<=storageN; j=j+1) {
                // N choose K = N/(N-K) times (N-1) choose K
                C_NK_Lookup << C_NK_Lookup.lData[C_NK_Lookup.lLength-1] * j/(j-i);
            }
        }
    }

    /*  for (long i=0; i<=storageK; i++)
        {
            long offset = i*(storageN+1);
            for (long j=0; j<=storageN; j++)
                printf ("(%d,%d) = %d\n", j,i, C_NK_Lookup.lData[offset+j]);
        }
    */

    matrixDimension = C_NK_Lookup.lData[C_NK_Lookup.lLength-1];
    CreateMatrix (this, 1, matrixDimension, false, true);
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

BaseRef _NTupleStorage::makeDynamic (void)
{
    _NTupleStorage* copy = new _NTupleStorage;
    checkPointer (copy);
    copy->_Matrix::Duplicate (this);
    copy->storageN = storageN;
    copy->storageK = storageK;
    copy->C_NK_Lookup.Duplicate (&C_NK_Lookup);
    return copy;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

bool    _NTupleStorage::CheckKTuple (_SimpleList& kTuple)
{
    if (kTuple.lLength == storageK) {
        if (storageK) {
            kTuple.Sort();
            for (long k=0; k<kTuple.lLength; k++)
                if (kTuple.lData[k] < 0 || kTuple.lData[k] >= storageN || (k && kTuple.lData[k] == kTuple.lData[k-1])) {
                    return false;
                }
        }
        return true;
    }
    return false;
}


/*--------------------------------------------------------------------------------------------------------------------------------*/

unsigned long   _NTupleStorage::Index (_SimpleList& kTuple)
{
    unsigned long myIndex = 0;
    if (storageK)
        for (long k=kTuple.lLength-1; k >=0; k--) {
            myIndex += C_NK_Lookup.lData[(k+1)*(1+storageN) + kTuple.lData[k]];
        }
    return myIndex;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_Parameter  _NTupleStorage::DirectIndex (unsigned long directIndex)
{
    return theData[directIndex];
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

unsigned long   _NTupleStorage::Store (_Parameter value, _SimpleList& kTuple)
{
    unsigned long myIndex = Index (kTuple);
    theData[myIndex] = value;
    return myIndex;
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

_Parameter  _NTupleStorage::Retrieve (_SimpleList& kTuple)
{
    unsigned long myIndex = Index (kTuple);
    return theData[myIndex];
}

/*--------------------------------------------------------------------------------------------------------------------------------*/

void    _NTupleStorage::IndexToTuple (unsigned long directIndex, _SimpleList& kTuple)
{
    kTuple.Clear();
    if (storageK && directIndex < C_NK_Lookup.lData[C_NK_Lookup.lLength-1]) {
        long currentN = storageN-1;
        for (long k=storageK; k > 0; k--) {
            long  i             = currentN,
                  lookup_offset = k*(storageN+1);

            while (C_NK_Lookup.lData[lookup_offset+i] > directIndex) {
                //printf ("(%d, %d) -> %d\n", k, i, C_NK_Lookup.lData[lookup_offset+i]);
                i--;
            }
            //printf ("(%d, %d) -> %d\n", k, i, C_NK_Lookup.lData[lookup_offset+i]);
            kTuple << i;
            //printf ("Stored %d\n", i);

            currentN     = i-1;
            directIndex -= C_NK_Lookup.lData[lookup_offset+i];
        }
    }
    kTuple.Flip();
}

