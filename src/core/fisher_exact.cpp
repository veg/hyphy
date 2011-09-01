/*

This file is a C++ translation of the Mehta-Patel network
algorithm for computing association p-values on NxM contigency
table

Original source from

    http://portal.acm.org/citation.cfm?id=214326&jmp=indexterms&dl=GUIDE&dl=ACM

Cleaned up and dapted for C++. Included a bug fix for negative keys
in hash tables in f3xact_.

Sergei L. Kosakovsky Pond,
October 2003

Function prototype is declared in matrix.h

*/

#include "hy_strings.h"
#include <stdlib.h>
#include <math.h>
#ifdef    __HYPHYDMALLOC__
#include "dmalloc.h"
#endif

/* ----------------------------------------------------------------------- */
// Function Prototypes

void        allocate_fexact_keys (long, long);
void        free_fexact_keys     (void);

int         fexact_(long,    long, double *, double , double , double , double *, double *);
int         isort_ (long*,   long *);
double      alogam_(double*, long *);
double      gammds_(double*, double *, long *);
int         f2xact_(long *, long *, double *,
                    double *, double *, double *,
                    double *, double *, double *, long *, long *, long *, long *, long *,
                    long *, long *, long *, double *, long *,
                    long *, double *, double *, double *,
                    long *, long *, double *);


int         f3xact_(long *, long *, long *, long *, double *, long *, double *,
                    long *, long *, long *, long *, long *, long *, long *, long *, long *, double *, double *, double *);
int         f4xact_(long *, long *, long *, long *, double *, double *, long *,
                    long *, long *, long *, long *, long *, long *, double *, double *);
int         f5xact_(double *, double *, long *, long *, long *, long *,
                    double *, long *, long *, long *, long *, long *, long *, long *, bool *);
int         f6xact_(long *,  long *, long *, long *, long *, long *, long *, long *);
int         f7xact_(long *,  long *, long *, long *, long *, long *);
int         f8xact_(long *,  long *, long *, long *, long *);
double      f9xact_(long *,  long *, long *, double *);
int         f10act_(long,  long *, long , long * , double * , bool * , double * , long * , long * , long *);
void        f11act_(long *,  long, long, long *);
long        i_dnnt (double *);
/* ----------------------------------------------------------------------- */

#define Minimum(a,b) (((a)<(b))?(a):(b))
#define Maximum(a,b) (((a)>(b))?(a):(b))

/* ----------------------------------------------------------------------- */
inline long i_dnnt(double *x)
{
    return( (*x)>=0 ? (long) (*x + .5) : (long) (*x - .5) );
}

long    *fexact_i4 = nil,
         *fexact_i5 = nil,
          *fexact_i7 = nil,
           *fexact_i10 = nil;

double  *fexact_i6 = nil,
         *fexact_i8 = nil,
          *fexact_i9 = nil,
           *fexact_i9a = nil;

long    fexact_ldkey = 0,
        fexact_ldstp = 0;

/* ----------------------------------------------------------------------- */

void    allocate_fexact_keys (long ldkey, long mult)
{
    fexact_ldkey = ldkey;
    fexact_ldstp = ldkey * mult;
    long    i__1    = ldkey << 1;

    fexact_i4  = (long*) MemAllocate(i__1*sizeof (long));
    fexact_i5 = (long*) MemAllocate(i__1*sizeof (long));
    i__1 = fexact_ldstp << 1;
    fexact_i6 = (double*) MemAllocate(i__1*sizeof (double));
    i__1 = fexact_ldstp * 6;
    fexact_i7 = (long*) MemAllocate(i__1*sizeof (long));
    i__1 = ldkey << 1;
    fexact_i8 = (double*) MemAllocate(i__1*sizeof (double));
    fexact_i9  = (double*) MemAllocate(i__1*sizeof (double));
    fexact_i9a = (double*) MemAllocate(i__1*sizeof (double));
    fexact_i10  = (long*) MemAllocate(i__1*sizeof (long));
}

/* ----------------------------------------------------------------------- */

void    free_fexact_keys (void)
{
    free (fexact_i4);
    free (fexact_i5);
    free (fexact_i6);
    free (fexact_i7);
    free (fexact_i8);
    free (fexact_i9);
    free (fexact_i9a);
    free (fexact_i10);

    fexact_i4 = nil;
}



/* ----------------------------------------------------------------------- */
/*  Name:       isort_ */

/*  Purpose:    Shell sort for an long vector. */

/*  Usage:      isort_ (N, IX) */

/*  Arguments: */
/*     n      - Lenth of vector IX.  (Input) */
/*     ix     - Vector to be sorted.  (Input/output) */
/* ----------------------------------------------------------------------- */

int isort_(long *n, long *ix)
{
    long i__,
         j,
         m,
         il[10],
         kl,
         it,
         iu[10],
         ku,
         ikey;
    /* Parameter adjustments */
    --ix;

    /* Function Body */
    m = 1;
    i__ = 1;
    j = *n;
L10:
    if (i__ >= j) {
        goto L40;
    }

    kl = i__;
    ku = j;
    ikey = i__;
    ++j;
    /*                                  Find element in first half */
L20:
    ++i__;
    if (i__ < j) {
        if (ix[ikey] > ix[i__]) {
            goto L20;
        }
    }
    /*                                  Find element in second half */
L30:
    --j;
    if (ix[j] > ix[ikey]) {
        goto L30;
    }
    /*                                  Interchange */
    if (i__ < j) {
        it = ix[i__];
        ix[i__] = ix[j];
        ix[j] = it;
        goto L20;
    }
    it = ix[ikey];
    ix[ikey] = ix[j];
    ix[j] = it;
    /*                                  Save upper and lower subscripts of */
    /*                                  the array yet to be sorted */
    if (m < 11) {
        if (j - kl < ku - j) {
            il[m - 1] = j + 1;
            iu[m - 1] = ku;
            i__ = kl;
            --j;
        } else {
            il[m - 1] = kl;
            iu[m - 1] = j - 1;
            i__ = j + 1;
            j = ku;
        }
        ++m;
        goto L10;
    } else {
        _String errorMsg ("Internal error in shell sort");
        WarnError (errorMsg);
    }
    /*                                  Use another segment */
L40:
    --m;
    if (m == 0) {
        goto L9000;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
    goto L10;

L9000:
    return 0;
} /* isort_ */

/* ----------------------------------------------------------------------- */
/*  Name:       GAMMDS */

/*  Purpose:    Cumulative distribution for the gamma distribution. */

/*  Usage:      PGAMMA (Q, ALPHA,IFAULT) */

/*  Arguments: */
/*     Q      - Value at which the distribution is desired.  (Input) */
/*     ALPHA  - Parameter in the gamma distribution.  (Input) */
/*     IFAULT - Error indicator.  (Output) */
/*               IFAULT  DEFINITION */
/*                 0     No error */
/*                 1     An argument is misspecified. */
/*                 2     A numerical error has occurred. */
/*     PGAMMA - The cdf for the gamma distribution with parameter alpha */
/*              evaluated at Q.  (Output) */
/* ----------------------------------------------------------------------- */

double gammds_(double *y, double *p, long *ifault)
{
    /* Initialized data */

    double e = 1e-6;
    double zero = 0.;
    double one = 1.;

    /* System generated locals */
    double ret_val, d__1, d__2;


    /* Local variables */
    double a, c__, f;
    long ifail;

    *ifault = 1;
    ret_val = zero;
    if (*y <= zero || *p <= zero) {
        return ret_val;
    }
    *ifault = 2;

    d__2 = *p + one;
    d__1 = *p * log(*y) - alogam_(&d__2, &ifail) - *y;
    f = exp(d__1);
    if (f == zero) {
        return ret_val;
    }
    *ifault = 0;

    /*       Series begins */

    c__ = one;
    ret_val = one;
    a = *p;
L10:
    a += one;
    c__ = c__ * *y / a;
    ret_val += c__;
    if (c__ / ret_val > e) {
        goto L10;
    }
    ret_val *= f;
    return ret_val;
} /* gammds_ */

/* ----------------------------------------------------------------------- */
/*  Name:       ALOGAM */

/*  Purpose:    Value of the log-gamma function. */

/*  Usage:      ALOGAM (X, IFAULT) */

/*  Arguments: */
/*     X      - Value at which the log-gamma function is to be evaluated. */
/*              (Input) */
/*     IFAULT  - Error indicator.  (Output) */
/*               IFAULT  DEFINITION */
/*                 0     No error */
/*                 1     X .LT. 0 */
/*     ALGAMA - The value of the log-gamma function at XX.  (Output) */
/* ----------------------------------------------------------------------- */

double alogam_(double *x, long *ifault)
{
    /* Initialized data */

    static double a1 = .918938533204673;
    static double a2 = 5.95238095238e-4;
    static double a3 = 7.93650793651e-4;
    static double a4 = .002777777777778;
    static double a5 = .083333333333333;
    static double half = .5;
    static double zero = 0.;
    static double one = 1.;
    static double seven = 7.;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double f, y, z__;

    ret_val = zero;
    *ifault = 1;
    if (*x < zero) {
        return ret_val;
    }
    *ifault = 0;
    y = *x;
    f = zero;
    if (y >= seven) {
        goto L30;
    }
    f = y;
L10:
    y += one;
    if (y >= seven) {
        goto L20;
    }
    f *= y;
    goto L10;
L20:
    f = -log(f);
L30:
    z__ = one / (y * y);
    ret_val = f + (y - half) * log(y) - y + a1 + (((-a2 * z__ + a3) * z__ -
              a4) * z__ + a5) / y;
    return ret_val;
} /* alogam_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F11ACT */

/*  Purpose:    Routine for revising row totals. */

/*  Arguments: */
/*     IROW   - Vector containing the row totals.  (Input) */
/*     I1     - Indicator.  (Input) */
/*     I2     - Indicator.  (Input) */
/*     N      - Vector containing the row totals.  (Input) */
/* ----------------------------------------------------------------------- */
void f11act_(long *irow, long i1, long i2, long *n)
{
    for (long i = 0; i < i1-1; i++) {
        n[i] = irow[i];
    }
    {
        for (long i = i1-1; i < i2; i++) {
            n[i] = irow[i + 1];
        }
    }
} /* f11act_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F10ACT */

/*  Purpose:    Computes the shortest path length for special tables. */

/*  Usage:      CALL F10ACT (NROW, IROW, NCOL, ICOL, VAL, XMIN, FACT, ND, */
/*                          NE, M) */

/*  Arguments: */
/*     NROW   - The number of rows in the table.  (Input) */
/*     IROW   - Vector of length NROW containing the row totals.  (Input) */
/*     NCOL   - The number of columns in the table.  (Input) */
/*     ICO    - Vector of length NCOL containing the column totals. */
/*              (Input) */
/*     VAL    - The shortest path.  (Output) */
/*     XMIN   - Set to true if shortest path obtained.  (Output) */
/*     FACT   - Vector containing the logarithms of factorials. */
/*              (Input) */
/*     ND     - Workspace vector of length NROW. */
/*     NE     - Workspace vector of length NCOL. */
/*     M      - Workspace vector of length NCOL. */

/*  Chapter:    STAT/LIBRARY Categorical and Discrete Data Analysis */
/* ----------------------------------------------------------------------- */
int f10act_(long nrow, long *irow, long ncol, long *icol, double *val, bool *xmin, double *fact, long *nd, long *ne, long *m)
{
    /* System generated locals */
    long is, ix;

    /* Parameter adjustments */

    long i__,
         i__1;

    --m;
    --ne;
    --nd;
    --icol;
    --irow;

    /* Function Body */

    //for (long i1 = 0; i1 < nrow; i1++)
    //  nd[i1] = 0;

    i__1 = nrow - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        nd[i__] = 0;
    }

    is = icol[1] / nrow;
    ne[1] = is;

    //is = icol[0]/nrow;
    //ne[0] = is;

    ix = icol[1] - nrow * is;
    m[1] = ix;

    //ix = icol[0] - nrow*is;
    //m[0] = ix;

    //if (ix != 0)
    //++nd[ix-1];

    if (ix != 0) {
        ++nd[ix];
    }

    i__1 = ncol;
    for (i__ = 2; i__ <= i__1; ++i__) {
        ix = icol[i__] / nrow;
        ne[i__] = ix;
        is += ix;
        ix = icol[i__] - nrow * ix;
        m[i__] = ix;
        if (ix != 0) {
            ++nd[ix];
        }
    }

    /*for (long i2 = 1; i2 < ncol; i2++)
    {
        long t = icol[i2]/nrow;
        ne[i2] = t;
        is += t;
        t = icol[i2] - nrow*t;
        m[i2] = t;
        if (t != 0)
            nd[t-1]++;
    }*/

    for (i__ = nrow - 2; i__ >= 1; --i__) {
        nd[i__] += nd[i__ + 1];
    }

    //for (long i3 = nrow-3; i3>=0; i3--)
    //  nd[i3] += nd [i3+1];

    ix = 0;
    long nrw1 = nrow + 1;
    for (i__ = nrow; i__ >= 2; --i__) {
        ix = ix + is + nd[nrw1 - i__] - irow[i__];
        if (ix < 0) {
            return 0;
        }
    }

    /*ix = 0;
    for (long i4 = nrow; i4>=1; i4--)
    {
        ix += is + nd[nrow-i4+1] - irow[i4];
        if (ix<0)
            return 0;
    }*/

    i__1 = ncol;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ix = ne[i__];
        is = m[i__];
        *val = *val + is * fact[ix + 1] + (nrow - is) * fact[ix];
    }
    *xmin = true;

    /* double   temp = 0.0;

     for (long i5 = 0; i5 < ncol; i5++)
     {
        ix = ne[i5];
        is = m[i5];

        temp += is *fact[ix+1] + (nrow-is)*fact[ix];
     }
     *val += temp;
     *xmin = true;*/

    return 0;
}

/* ----------------------------------------------------------------------- */
/*  Name:       F9XACT */

/*  Purpose:    Computes the log of a multinomial coefficient. */

/*  Usage:      F9XACT(N, MM, IR, FACT) */

/*  Arguments: */
/*     N      - Length of IR.  (Input) */
/*     MM     - Number for factorial in numerator.  (Input) */
/*     IR     - Vector of length N containing the numebers for the */
/*              denominator of the factorial.  (Input) */
/*     FACT   - Table of log factorials.  (Input) */
/*     F9XACT  - The log of the multinomal coefficient.  (Output) */
/* ----------------------------------------------------------------------- */
double f9xact_(long *n, long *mm, long *ir, double *fact)
{
    /* System generated locals */
    long    i__1,
            k;

    double  ret_val;

    /* Parameter adjustments */
    --ir;

    /* Function Body */
    ret_val = fact[*mm];
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
        ret_val -= fact[ir[k]];
    }

    return ret_val;
} /* f9xact_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F8XACT */

/*  Purpose:    Routine for reducing a vector when there is a zero */
/*              element. */

/*  Usage:      CALL F8XACT (IROW, IS, I1, IZERO, NEW) */

/*  Arguments: */
/*     IROW   - Vector containing the row counts.  (Input) */
/*     IS     - Indicator.  (Input) */
/*     I1     - Indicator.  (Input) */
/*     IZERO  - Position of the zero.  (Input) */
/*     NEW    - Vector of new row counts.  (Output) */
/* ----------------------------------------------------------------------- */
int f8xact_(long *irow, long *is, long *i1, long *izero, long *new__)
{
    /* System generated locals */
    long i__1,
         i__;

    /* Parameter adjustments */
    --new__;
    --irow;

    /* Function Body */
    i__1 = *i1 - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
        new__[i__] = irow[i__];
        /* L10: */
    }

    i__1 = *izero - 1;
    for (i__ = *i1; i__ <= i__1; ++i__) {
        if (*is >= irow[i__ + 1]) {
            goto L30;
        }
        new__[i__] = irow[i__ + 1];
        /* L20: */
    }

    i__ = *izero;
L30:
    new__[i__] = *is;
L40:
    ++i__;
    if (i__ > *izero) {
        return 0;
    }
    new__[i__] = irow[i__];
    goto L40;
} /* f8xact_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F7XACT */

/*  Purpose:    Generate the new nodes for given marinal totals. */

/*  Usage:      CALL F7XACT (NROW, IMAX, IDIF, K, KS, IFLAG) */

/*  Arguments: */
/*     NROW   - The number of rows in the table.  (Input) */
/*     IMAX   - The row marginal totals.  (Input) */
/*     IDIF   - The column counts for the new column.  (Input/output) */
/*     K      - Indicator for the row to decrement.  (Input/output) */
/*     KS     - Indicator for the row to increment.  (Input/output) */
/*     IFLAG  - Status indicator.  (Output) */
/*              If IFLAG is zero, a new table was generated.  For */
/*              IFLAG = 1, no additional tables could be generated. */
/* ----------------------------------------------------------------------- */
int f7xact_(long *nrow, long *imax, long *idif, long *k, long *ks, long *iflag)
{
    /* System generated locals */
    long i__1, i__2, i__, m, k1, mm;

    /* Parameter adjustments */
    --idif;
    --imax;

    /* Function Body */
    *iflag = 0;
    /*                                  Find node which can be */
    /*                                  incremented, ks */
    if (*ks == 0) {
L10:
        ++(*ks);
        if (idif[*ks] == imax[*ks]) {
            goto L10;
        }
    }
    /*                                 Find node to decrement (>ks) */
    /* L20: */
    if (idif[*k] > 0 && *k > *ks) {
        --idif[*k];
L30:
        --(*k);
        if (imax[*k] == 0) {
            goto L30;
        }
        m = *k;
        /*                                 Find node to increment (>=ks) */
L40:
        if (idif[m] >= imax[m]) {
            --m;
            goto L40;
        }
        ++idif[m];
        /*                                 Change ks */
        if (m == *ks) {
            if (idif[m] == imax[m]) {
                *ks = *k;
            }
        }
    } else {
        /*                                 Check for finish */
L50:
        i__1 = *nrow;
        for (k1 = *k + 1; k1 <= i__1; ++k1) {
            if (idif[k1] > 0) {
                goto L70;
            }
            /* L60: */
        }
        *iflag = 1;
        goto L9000;
        /*                                 Reallocate counts */
L70:
        mm = 1;
        i__1 = *k;
        for (i__ = 1; i__ <= i__1; ++i__) {
            mm += idif[i__];
            idif[i__] = 0;
            /* L80: */
        }
        *k = k1;
L90:
        --(*k);
        /* Computing MIN */
        i__1 = mm, i__2 = imax[*k];
        m = Minimum(i__1,i__2);
        idif[*k] = m;
        mm -= m;
        if (mm > 0 && *k != 1) {
            goto L90;
        }
        /*                                 Check that all counts */
        /*                                 reallocated */
        if (mm > 0) {
            if (k1 != *nrow) {
                *k = k1;
                goto L50;
            }
            *iflag = 1;
            goto L9000;
        }
        /*                                 Get ks */
        --idif[k1];
        *ks = 0;
L100:
        ++(*ks);
        if (*ks > *k) {
            goto L9000;
        }
        if (idif[*ks] >= imax[*ks]) {
            goto L100;
        }
    }

L9000:
    return 0;
} /* f7xact_ */


/* ----------------------------------------------------------------------- */
/*  Name:       F6XACT */

/*  Purpose:    Pop a node off the stack. */

/*  Usage:      CALL F6XACT (NROW, IROW, IFLAG, KYY, KEY, LDKEY, LAST, */
/*                          IPN) */

/*  Arguments: */
/*     NROW   - The number of rows in the table.  (Input) */
/*     IROW   - Vector of length nrow containing the row sums on output. */
/*              (Output) */
/*     IFLAG  - Set to 3 if there are no additional nodes to process. */
/*              (Output) */
/*     KYY    - Constant mutlipliers used in forming the hash table key. */
/*              (Input) */
/*     KEY    - Vector of length LDKEY containing the hash table keys. */
/*              (Input/output) */
/*     LDKEY  - Length of vector KEY.  (Input) */
/*     LAST   - Index of the last key popped off the stack. */
/*              (Input/output) */
/*     IPN    - Pointer to the linked list of past path lengths. */
/*              (Output) */
/* ----------------------------------------------------------------------- */
int f6xact_(long *nrow, long *irow, long *iflag, long *kyy, long *key, long *ldkey, long *last, long *ipn)
{
    long j,
         kval;


    /* Parameter adjustments */
    --key;
    --kyy;
    --irow;

    /* Function Body */
L10:
    ++(*last);
    if (*last <= *ldkey) {
        if (key[*last] < 0) {
            goto L10;
        }
        /*                                  Get KVAL from the stack */
        kval = key[*last];
        key[*last] = -9999;
        for (j = *nrow; j >= 2; --j) {
            irow[j] = kval / kyy[j];
            kval -= irow[j] * kyy[j];
            /* L20: */
        }
        irow[1] = kval;
        *ipn = *last;
    } else {
        *last = 0;
        *iflag = 3;
    }
    return 0;
} /* f6xact_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F5XACT */

/*  Purpose:    Put node on stack in network algorithm. */

/*  Usage:      CALL F5XACT (PASTP, TOL, KVAL, KEY, LDKEY, IPOIN, STP, */
/*                          LDSTP, IFRQ, NPOIN, NR, NL, IFREQ, ITOP, */
/*                          IPSH) */

/*  Arguments: */
/*     PASTP  - The past path length.  (Input) */
/*     TOL    - Tolerance for equivalence of past path lengths.  (Input) */
/*     KVAL   - Key value.  (Input) */
/*     KEY    - Vector of length LDKEY containing the key values. */
/*              (Input/output) */
/*     LDKEY  - Length of vector KEY.  (Input) */
/*     IPOIN  - Vector of length LDKEY pointing to the linked list */
/*              of past path lengths.  (Input/output) */
/*     STP    - Vector of length LSDTP containing the linked lists */
/*              of past path lengths.  (Input/output) */
/*     LDSTP  - Length of vector STP.  (Input) */
/*     IFRQ   - Vector of length LDSTP containing the past path */
/*              frequencies.  (Input/output) */
/*     NPOIN  - Vector of length LDSTP containing the pointers to */
/*              the next past path length.  (Input/output) */
/*     NR     - Vector of length LDSTP containing the right object */
/*              pointers in the tree of past path lengths. */
/*              (Input/output) */
/*     NL     - Vector of length LDSTP containing the left object */
/*              pointers in the tree of past path lengths. */
/*              (Input/output) */
/*     IFREQ  - Frequency of the current path length.  (Input) */
/*     ITOP   - Pointer to the top of STP.  (Input) */
/*     IPSH   - Option parameter.  (Input) */
/*              If IPSH is true, the past path length is found in the */
/*              table KEY.  Otherwise the location of the past path */
/*              length is assumed known and to have been found in */
/*              a previous call. */
/* ----------------------------------------------------------------------- */
int f5xact_(double *pastp, double *tol, long *kval, long *key, long *ldkey, long *ipoin,
            double *stp, long *ldstp, long *ifrq, long *npoin, long *nr, long *nl, long *ifreq, long *itop, bool *ipsh)
{
    /* System generated locals */
    static  long    itp;

    long    i__1,
            ird,
            ipn,
            itmp;

    double  test1,
            test2;

    _String errMsg;

    /* Parameter adjustments */
    --nl;
    --nr;
    --npoin;
    --ifrq;
    --stp;
    --ipoin;
    --key;

    /* Function Body */
    if (*ipsh) {
        /*                                  Convert KVAL to long in range */
        /*                                  1, ..., LDKEY. */
        ird = *kval % *ldkey + 1;
        /*                                  Search for an unused location */
        i__1 = *ldkey;
        for (itp = ird; itp <= i__1; ++itp) {
            if (key[itp] == *kval) {
                goto L40;
            }
            if (key[itp] < 0) {
                goto L30;
            }
            /* L10: */
        }

        i__1 = ird - 1;
        for (itp = 1; itp <= i__1; ++itp) {
            if (key[itp] == *kval) {
                goto L40;
            }
            if (key[itp] < 0) {
                goto L30;
            }
            /* L20: */
        }
        /*                                  Return if KEY array is full */
        errMsg = "Fisher Exact:LDKEY is too small for this problem.  It is not possible to estimate the value of LDKEY required, but twice the current value may be sufficient.";
        WarnError (errMsg);
        return 0;
        /*                                  Update KEY */
L30:
        key[itp] = *kval;
        ++(*itop);
        ipoin[itp] = *itop;
        /*                                  Return if STP array full */
        if (*itop > *ldstp) {
            errMsg = "Fisher Exact: LDSTP is too small for this problem.  It is not possible to estimate the value of LDSTP required, but twice the current value may be sufficient.";
            WarnError (errMsg);
            return 0;
        }
        /*                                  Update STP, etc. */
        npoin[*itop] = -1;
        nr[*itop] = -1;
        nl[*itop] = -1;
        stp[*itop] = *pastp;
        ifrq[*itop] = *ifreq;
        goto L9000;
    }
    /*                                  Find location, if any, of pastp */
L40:
    ipn = ipoin[itp];
    test1 = *pastp - *tol;
    test2 = *pastp + *tol;

L50:
    if (stp[ipn] < test1) {
        ipn = nl[ipn];
        if (ipn > 0) {
            goto L50;
        }
    } else if (stp[ipn] > test2) {
        ipn = nr[ipn];
        if (ipn > 0) {
            goto L50;
        }
    } else {
        ifrq[ipn] += *ifreq;
        goto L9000;
    }
    /*                                  Return if STP array full */
    ++(*itop);
    if (*itop > *ldstp) {
        errMsg = "Fisher Exact: LDSTP is too small for this problem.  It is not possible to estimate the value of LDSTP required, but twice the current value may be sufficient.";
        WarnError (errMsg);
        goto L9000;
    }
    /*                                  Find location to add value */
    ipn = ipoin[itp];
    itmp = ipn;
L60:
    if (stp[ipn] < test1) {
        itmp = ipn;
        ipn = nl[ipn];
        if (ipn > 0) {
            goto L60;
        } else {
            nl[itmp] = *itop;
        }
    } else if (stp[ipn] > test2) {
        itmp = ipn;
        ipn = nr[ipn];
        if (ipn > 0) {
            goto L60;
        } else {
            nr[itmp] = *itop;
        }
    }
    /*                                  Update STP, etc. */
    npoin[*itop] = npoin[itmp];
    npoin[itmp] = *itop;
    stp[*itop] = *pastp;
    ifrq[*itop] = *ifreq;
    nl[*itop] = -1;
    nr[*itop] = -1;

L9000:
    return 0;
} /* f5xact_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F4XACT */

/*  Purpose:    Computes the longest path length for a given table. */

/*  Usage:      CALL F4XACT (NROW, IROW, NCOL, ICOL, DSP, FACT, ICSTK, */
/*                          NCSTK, LSTK, MSTK, NSTK, NRSTK, IRSTK, YSTK, */
/*                          TOL) */

/*  Arguments: */
/*     NROW   - The number of rows in the table.  (Input) */
/*     IROW   - Vector of length NROW containing the row sums for the */
/*              table.  (Input) */
/*     NCOL   - The number of columns in the table.  (Input) */
/*     ICOL   - Vector of length K containing the column sums for the */
/*              table.  (Input) */
/*     DSP    - The shortest path for the table.  (Output) */
/*     FACT   - Vector containing the logarithms of factorials.  (Input) */
/*     ICSTK  - NCOL by NROW+NCOL+1 work array. */
/*     NCSTK  - Work vector of length NROW+NCOL+1. */
/*     LSTK   - Work vector of length NROW+NCOL+1. */
/*     MSTK   - Work vector of length NROW+NCOL+1. */
/*     NSTK   - Work vector of length NROW+NCOL+1. */
/*     NRSTK  - Work vector of length NROW+NCOL+1. */
/*     IRSTK  - NROW by MAX(NROW,NCOL) work array. */
/*     YSTK   - Work vector of length NROW+NCOL+1. */
/*     TOL    - Tolerance.  (Input) */
/* ----------------------------------------------------------------------- */
int f4xact_(long *nrow, long *irow, long *ncol,
            long *icol, double *dsp, double *fact, long *icstk,
            long *ncstk, long *lstk, long *mstk, long *nstk,
            long *nrstk, long *irstk, double *ystk, double *tol)
{
    /* System generated locals */
    long icstk_dim1,
         icstk_offset,
         irstk_dim1,
         irstk_offset,
         i__1,
         i__,
         j,
         k,
         l,
         m,
         n,
         mn,
         ic1,
         ir1,
         ict,
         nco,
         irt,
         nro,
         istk;

    double y,
           amx;


    /* Parameter adjustments */
    irstk_dim1 = *nrow;
    irstk_offset = 1 + irstk_dim1;
    irstk -= irstk_offset;
    --irow;
    icstk_dim1 = *ncol;
    icstk_offset = 1 + icstk_dim1;
    icstk -= icstk_offset;
    --icol;
    --ncstk;
    --lstk;
    --mstk;
    --nstk;
    --nrstk;
    --ystk;

    /* Function Body */
    if (*nrow == 1) {
        i__1 = *ncol;
        for (i__ = 1; i__ <= i__1; ++i__) {
            *dsp -= fact[icol[i__]];
            /* L10: */
        }
        goto L9000;
    }

    if (*ncol == 1) {
        i__1 = *nrow;
        for (i__ = 1; i__ <= i__1; ++i__) {
            *dsp -= fact[irow[i__]];
            /* L20: */
        }
        goto L9000;
    }

    if (*nrow * *ncol == 4) {
        if (irow[2] <= icol[2]) {
            *dsp = *dsp - fact[irow[2]] - fact[icol[1]] - fact[icol[2] - irow[
                        2]];
        } else {
            *dsp = *dsp - fact[icol[2]] - fact[irow[1]] - fact[irow[2] - icol[
                        2]];
        }
        goto L9000;
    }
    /*                                  initialization before loop */
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
        irstk[i__ + irstk_dim1] = irow[*nrow - i__ + 1];
        /* L30: */
    }

    i__1 = *ncol;
    for (j = 1; j <= i__1; ++j) {
        icstk[j + icstk_dim1] = icol[*ncol - j + 1];
        /* L40: */
    }

    nro = *nrow;
    nco = *ncol;
    nrstk[1] = nro;
    ncstk[1] = nco;
    ystk[1] = 0.f;
    y = 0.f;
    istk = 1;
    l = 1;
    amx = 0.f;

L50:
    ir1 = irstk[istk * irstk_dim1 + 1];
    ic1 = icstk[istk * icstk_dim1 + 1];
    if (ir1 > ic1) {
        if (nro >= nco) {
            m = nco - 1;
            n = 2;
        } else {
            m = nro;
            n = 1;
        }
    } else if (ir1 < ic1) {
        if (nro <= nco) {
            m = nro - 1;
            n = 1;
        } else {
            m = nco;
            n = 2;
        }
    } else {
        if (nro <= nco) {
            m = nro - 1;
            n = 1;
        } else {
            m = nco - 1;
            n = 2;
        }
    }

L60:
    if (n == 1) {
        i__ = l;
        j = 1;
    } else {
        i__ = 1;
        j = l;
    }

    irt = irstk[i__ + istk * irstk_dim1];
    ict = icstk[j + istk * icstk_dim1];
    mn = irt;
    if (mn > ict) {
        mn = ict;
    }
    y += fact[mn];
    if (irt == ict) {
        --nro;
        --nco;
        f11act_(&irstk[istk * irstk_dim1 + 1], i__, nro, &irstk[(istk + 1) *
                irstk_dim1 + 1]);
        f11act_(&icstk[istk * icstk_dim1 + 1], j, nco, &icstk[(istk + 1) *
                icstk_dim1 + 1]);
    } else if (irt > ict) {
        --nco;
        f11act_(&icstk[istk * icstk_dim1 + 1], j, nco, &icstk[(istk + 1) *
                icstk_dim1 + 1]);
        i__1 = irt - ict;
        f8xact_(&irstk[istk * irstk_dim1 + 1], &i__1, &i__, &nro, &irstk[(
                    istk + 1) * irstk_dim1 + 1]);
    } else {
        --nro;
        f11act_(&irstk[istk * irstk_dim1 + 1], i__, nro, &irstk[(istk + 1) *
                irstk_dim1 + 1]);
        i__1 = ict - irt;
        f8xact_(&icstk[istk * icstk_dim1 + 1], &i__1, &j, &nco, &icstk[(istk
                + 1) * icstk_dim1 + 1]);
    }

    if (nro == 1) {
        i__1 = nco;
        for (k = 1; k <= i__1; ++k) {
            y += fact[icstk[k + (istk + 1) * icstk_dim1]];
            /* L70: */
        }
        goto L90;
    }

    if (nco == 1) {
        i__1 = nro;
        for (k = 1; k <= i__1; ++k) {
            y += fact[irstk[k + (istk + 1) * irstk_dim1]];
            /* L80: */
        }
        goto L90;
    }

    lstk[istk] = l;
    mstk[istk] = m;
    nstk[istk] = n;
    ++istk;
    nrstk[istk] = nro;
    ncstk[istk] = nco;
    ystk[istk] = y;
    l = 1;
    goto L50;

L90:
    if (y > amx) {
        amx = y;
        if (*dsp - amx <= *tol) {
            *dsp = 0.f;
            goto L9000;
        }
    }

L100:
    --istk;
    if (istk == 0) {
        *dsp -= amx;
        if (*dsp - amx <= *tol) {
            *dsp = 0.f;
        }
        goto L9000;
    }
    l = lstk[istk] + 1;

L110:
    if (l > mstk[istk]) {
        goto L100;
    }
    n = nstk[istk];
    nro = nrstk[istk];
    nco = ncstk[istk];
    y = ystk[istk];
    if (n == 1) {
        if (irstk[l + istk * irstk_dim1] < irstk[l - 1 + istk * irstk_dim1]) {
            goto L60;
        }
    } else if (n == 2) {
        if (icstk[l + istk * icstk_dim1] < icstk[l - 1 + istk * icstk_dim1]) {
            goto L60;
        }
    }

    ++l;
    goto L110;
L9000:
    return 0;
} /* f4xact_ */


/* ----------------------------------------------------------------------- */
/*  Name:       F3XACT */

/*  Purpose:    Computes the shortest path length for a given table. */

/*  Usage:      CALL F3XACT (NROW, IROW, NCOL, ICOL, DLP, MM, FACT, ICO, */
/*                          IRO, IT, LB, NR, NT, NU, ITC, IST, STV, ALEN, */
/*                          TOL) */

/*  Arguments: */
/*     NROW   - The number of rows in the table.  (Input) */
/*     IROW   - Vector of length NROW containing the row sums for the */
/*              table.  (Input) */
/*     NCOL   - The number of columns in the table.  (Input) */
/*     ICOL   - Vector of length K containing the column sums for the */
/*              table.  (Input) */
/*     DLP    - The longest path for the table.  (Output) */
/*     MM     - The total count in the table.  (Output) */
/*     FACT   - Vector containing the logarithms of factorials.  (Input) */
/*     ICO    - Work vector of length MAX(NROW,NCOL). */
/*     IRO    - Work vector of length MAX(NROW,NCOL). */
/*     IT     - Work vector of length MAX(NROW,NCOL). */
/*     LB     - Work vector of length MAX(NROW,NCOL). */
/*     NR     - Work vector of length MAX(NROW,NCOL). */
/*     NT     - Work vector of length MAX(NROW,NCOL). */
/*     NU     - Work vector of length MAX(NROW,NCOL). */
/*     ITC    - Work vector of length 400. */
/*     IST    - Work vector of length 400. */
/*     STV    - Work vector of length 400. */
/*     ALEN   - Work vector of length MAX(NROW,NCOL). */
/*     TOL    - Tolerance.  (Input) */
/* ----------------------------------------------------------------------- */
int f3xact_(long *nrow, long *irow, long *ncol,
            long *icol, double *dlp, long *mm, double *fact,
            long *ico, long *iro, long *it, long *lb, long *nr,
            long *nt, long *nu, long *itc, long *ist, double *stv,
            double *alen, double *tol)
{
    /* Initialized data */

    long ldst = 200;
    long nst = 0;
    long nitc = 0;

    /* System generated locals */
    long i__1;
    double d__1, d__2;

    /* Local variables */
    long i__, k;
    double v;
    long n11, n12, ii, nn, ks, ic1, ic2, nc1, nn1, nr1, nco;
    double val;
    long nct, ipn, irl, key, lev, itp, nro;
    double vmn;
    long nrt, kyy, nc1s;
    bool xmin;
    _String errMsg;

    /* Parameter adjustments */
    --stv;
    --ist;
    --itc;
    --nu;
    --nt;
    --nr;
    --lb;
    --it;
    --iro;
    --ico;
    --icol;
    --irow;

    /* Function Body */

    i__1 = *ncol;
    for (i__ = 0; i__ <= i__1; ++i__) {
        alen[i__] = 0.f;
        /* L10: */
    }
    for (i__ = 1; i__ <= 400; ++i__) {
        ist[i__] = -1;
        /* L20: */
    }
    /*                                  nrow is 1 */
    if (*nrow <= 1) {
        if (*nrow > 0) {
            *dlp -= fact[icol[1]];
            i__1 = *ncol;
            for (i__ = 2; i__ <= i__1; ++i__) {
                *dlp -= fact[icol[i__]];
                /* L30: */
            }
        }
        goto L9000;
    }
    /*                                  ncol is 1 */
    if (*ncol <= 1) {
        if (*ncol > 0) {
            *dlp = *dlp - fact[irow[1]] - fact[irow[2]];
            i__1 = *nrow;
            for (i__ = 3; i__ <= i__1; ++i__) {
                *dlp -= fact[irow[i__]];
                /* L40: */
            }
        }
        goto L9000;
    }
    /*                                  2 by 2 table */
    if (*nrow * *ncol == 4) {
        n11 = (irow[1] + 1) * (icol[1] + 1) / (*mm + 2);
        n12 = irow[1] - n11;
        *dlp = *dlp - fact[n11] - fact[n12] - fact[icol[1] - n11] - fact[icol[
                    2] - n12];
        goto L9000;
    }
    /*                                  Test for optimal table */
    val = 0.f;
    xmin = false;
    if (irow[*nrow] <= irow[1] + *ncol) {
        f10act_(*nrow, &irow[1], *ncol, &icol[1], &val, &xmin, fact, &lb[1], &
                nu[1], &nr[1]);
    }
    if (! xmin) {
        if (icol[*ncol] <= icol[1] + *nrow) {
            f10act_(*ncol, &icol[1], *nrow, &irow[1], &val, &xmin, fact, &lb[1],
                    &nu[1], &nr[1]);
        }
    }

    if (xmin) {
        *dlp -= val;
        goto L9000;
    }
    /*                                  Setup for dynamic programming */
    nn = *mm;
    /*                                  Minimize ncol */
    if (*nrow >= *ncol) {
        nro = *nrow;
        nco = *ncol;

        i__1 = *nrow;
        for (i__ = 1; i__ <= i__1; ++i__) {
            iro[i__] = irow[i__];
            /* L50: */
        }

        ico[1] = icol[1];
        nt[1] = nn - ico[1];
        i__1 = *ncol;
        for (i__ = 2; i__ <= i__1; ++i__) {
            ico[i__] = icol[i__];
            nt[i__] = nt[i__ - 1] - ico[i__];
            /* L60: */
        }
    } else {
        nro = *ncol;
        nco = *nrow;

        ico[1] = irow[1];
        nt[1] = nn - ico[1];
        i__1 = *nrow;
        for (i__ = 2; i__ <= i__1; ++i__) {
            ico[i__] = irow[i__];
            nt[i__] = nt[i__ - 1] - ico[i__];
            /* L70: */
        }

        i__1 = *ncol;
        for (i__ = 1; i__ <= i__1; ++i__) {
            iro[i__] = icol[i__];
            /* L80: */
        }
    }
    /*                                  Initialize pointers */
    vmn = 1e10;
    nc1s = nco - 1;
    irl = 1;
    ks = 0;
    k = ldst;
    kyy = ico[nco] + 1;
    goto L100;
    /*                                  Test for optimality */
L90:
    xmin = false;
    if (iro[nro] <= iro[irl] + nco) {
        f10act_(nro, &iro[irl], nco, &ico[1], &val, &xmin, fact, &lb[1], &
                nu[1], &nr[1]);
    }
    if (! xmin) {
        if (ico[nco] <= ico[1] + nro) {
            f10act_(nco, &ico[1], nro, &iro[irl], &val, &xmin, fact, &lb[1],
                    &nu[1], &nr[1]);
        }
    }

    if (xmin) {
        if (val < vmn) {
            vmn = val;
        }
        goto L200;
    }
    /*                                  Setup to generate new node */
L100:
    lev = 1;
    nr1 = nro - 1;
    nrt = iro[irl];
    nct = ico[1];
    lb[1] = (long) ((double) ((nrt + 1) * (nct + 1)) / (double) (
                        nn + nr1 * nc1s + 1) - *tol) - 1;
    nu[1] = (long) ((double) ((nrt + nc1s) * (nct + nr1)) / (
                        double) (nn + nr1 + nc1s)) - lb[1] + 1;
    nr[1] = nrt - lb[1];
    /*                                  Generate a node */
L110:
    --nu[lev];
    if (nu[lev] == 0) {
        if (lev == 1) {
            goto L200;
        }
        --lev;
        goto L110;
    }
    ++lb[lev];
    --nr[lev];
L120:
    alen[lev] = alen[lev - 1] + fact[lb[lev]];
    if (lev < nc1s) {
        nn1 = nt[lev];
        nrt = nr[lev];
        ++lev;
        nc1 = nco - lev;
        nct = ico[lev];
        lb[lev] = (long) ((double) ((nrt + 1) * (nct + 1)) / (
                              double) (nn1 + nr1 * nc1 + 1) - *tol);
        nu[lev] = (long) ((double) ((nrt + nc1) * (nct + nr1)) / (
                              double) (nn1 + nr1 + nc1) - lb[lev] + 1);
        nr[lev] = nrt - lb[lev];
        goto L120;
    }
    alen[nco] = alen[lev] + fact[nr[lev]];
    lb[nco] = nr[lev];

    v = val + alen[nco];
    if (nro == 2) {
        /*                                  Only 1 row left */
        v = v + fact[ico[1] - lb[1]] + fact[ico[2] - lb[2]];
        i__1 = nco;
        for (i__ = 3; i__ <= i__1; ++i__) {
            v += fact[ico[i__] - lb[i__]];
            /* L130: */
        }
        if (v < vmn) {
            vmn = v;
        }
    } else if (nro == 3 && nco == 2) {
        /*                                  3 rows and 2 columns */
        nn1 = nn - iro[irl] + 2;
        ic1 = ico[1] - lb[1];
        ic2 = ico[2] - lb[2];
        n11 = (iro[irl + 1] + 1) * (ic1 + 1) / nn1;
        n12 = iro[irl + 1] - n11;
        v = v + fact[n11] + fact[n12] + fact[ic1 - n11] + fact[ic2 - n12];
        if (v < vmn) {
            vmn = v;
        }
    } else {
        /*                                  Column marginals are new node */
        i__1 = nco;
        for (i__ = 1; i__ <= i__1; ++i__) {
            it[i__] = ico[i__] - lb[i__];
            /* L140: */
        }
        /*                                  Sort column marginals */
        if (nco == 2) {
            if (it[1] > it[2]) {
                ii = it[1];
                it[1] = it[2];
                it[2] = ii;
            }
        } else if (nco == 3) {
            ii = it[1];
            if (ii > it[3]) {
                if (ii > it[2]) {
                    if (it[2] > it[3]) {
                        it[1] = it[3];
                        it[3] = ii;
                    } else {
                        it[1] = it[2];
                        it[2] = it[3];
                        it[3] = ii;
                    }
                } else {
                    it[1] = it[3];
                    it[3] = it[2];
                    it[2] = ii;
                }
            } else if (ii > it[2]) {
                it[1] = it[2];
                it[2] = ii;
            } else if (it[2] > it[3]) {
                ii = it[2];
                it[2] = it[3];
                it[3] = ii;
            }
        } else {
            isort_(&nco, &it[1]);
        }
        /*                                  Compute hash value */
        key = it[1] * kyy + it[2];
        for (i__ = 3; i__ <= nco; ++i__) {
            key = it[i__] + key * kyy;
            /* L150: */
        }
        /*                                  Table index */

        /*if(key < 0)
        {
            errMsg = "Negative key in f3xact";
            WarnError (errMsg);
            //return 0;
        }*/

        if ((ipn = key % ldst + 1) < 1) {
            ipn += ldst;
        }
        /*                                  Find empty position */
        ii = ks + ipn;

        for (itp = ipn; itp <= ldst; ++itp) {
            if (ist[ii] < 0) {
                goto L180;
            } else if (ist[ii] == key) {
                goto L190;
            }
            ++ii;
            /* L160: */
        }

        ii = ks + 1;
        i__1 = ipn - 1;
        for (itp = 1; itp <= i__1; ++itp) {
            if (ist[ii] < 0) {
                goto L180;
            } else if (ist[ii] == key) {
                goto L190;
            }
            ++ii;
            /* L170: */
        }

        errMsg = ("Fisher Exact: Stack length exceeded in f3xact.  This problem should not occur.");
        ReportWarning (errMsg);
        return 0;
        /*                                  Push onto stack */
L180:
        ist[ii] = key;
        stv[ii] = v;
        ++nst;
        ii = nst + ks;
        itc[ii] = itp;
        goto L110;
        /*                                  Marginals already on stack */
L190:
        /* Computing MIN */
        d__1 = v, d__2 = stv[ii];
        stv[ii] = Minimum(d__1,d__2);
    }
    goto L110;
    /*                                  Pop item from stack */
L200:
    if (nitc > 0) {
        /*                                  Stack index */
        itp = itc[nitc + k] + k;
        --nitc;
        val = stv[itp];
        key = ist[itp];
        ist[itp] = -1;
        /*                                  Compute marginals */
        for (i__ = nco; i__ >= 2; --i__) {
            ico[i__] = key % kyy;
            key /= kyy;
            /* L210: */
        }
        ico[1] = key;
        /*                                  Set up nt array */
        nt[1] = nn - ico[1];
        i__1 = nco;
        for (i__ = 2; i__ <= i__1; ++i__) {
            nt[i__] = nt[i__ - 1] - ico[i__];
            /* L220: */
        }
        goto L90;

    } else if (nro > 2 && nst > 0) {
        /*                                  Go to next level */
        nitc = nst;
        nst = 0;
        k = ks;
        ks = ldst - ks;
        nn -= iro[irl];
        ++irl;
        --nro;
        goto L200;
    }

    *dlp -= vmn;
L9000:
    return 0;
} /* f3xact_ */

/* ----------------------------------------------------------------------- */
/*  Name:       F2XACT */

/*  Purpose:    Computes Fisher's exact test for a contingency table, */
/*              routine with workspace variables specified. */

/*  Usage:      CALL F2XACT (NROW, NCOL, TABLE, EXPECT, PERCNT, */
/*                          EMIN, PRT, PRE, FACT, ICO, IRO, KYY, IDIF, */
/*                          IRN, KEY, LDKEY, IPOIN, STP, LDSTP, IFRQ, */
/*                          DLP, DSP, TM, KEY2, IWK, RWK) */
/* ----------------------------------------------------------------------- */
int f2xact_(long *nrow, long *ncol, double *table,
            double *expect, double *percnt, double *
            emin, double *prt, double *pre, double *fact, long *
            ico, long *iro, long *kyy, long *idif, long *irn, long
            *key, long *ldkey, long *ipoin, double *stp, long *ldstp,
            long *ifrq, double *dlp, double *dsp, double *tm,
            long *key2, long *iwk, double *rwk)
{
    /* Initialized data */

    static long   imax = 2147483647;
    static double amiss = -12345.f;
    static double tol = 3.45254e-7;
    static double emx = 1e30f;

    /* System generated locals */
    long table_dim1, table_offset, i__1, i__2;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    long i__, j, k, n, k1;
    double dd, df;
    long i31, i32, i33, i34, i35, i36, i37, i38, i39, i41, i42, i43,
         i44, i45, i46, i47, i48, ii, kb, kd, ks;
    double pv;
    long i310, i311;
    double ddf;
    long nco, nrb;
    double emn, drn, dro, obs;
    long ipn, ipo, itp, nro;
    double tmp, obs2, obs3;
    long nro2, kval, kmax, jkey, last;
    bool ipsh;
    long itmp;
    double dspt;
    long itop, jstp, ntot, jstp2, jstp3, jstp4, iflag, ncell, ifreq;
    bool chisq;
    long ikkey;
    double pastp;
    long ikstp;
    long ikstp2;
    long ifault;

    _String errMsg ("Fisher Exact:");

    /*                                  SPECIFICATIONS FOR ARGUMENTS */
    /*                                  SPECIFICATIONS FOR LOCAL VARIABLES */
    /*                                  SPECIFICATIONS FOR INTRINSICS */
    /*                                  SPECIFICATIONS FOR SUBROUTINES */
    /*                                  SPECIFICATIONS FOR FUNCTIONS */
    /* *********************************************************************** */
    /*                                  IMAX is the largest representable */
    /*                                  long on the machine */
    /* *********************************************************************** */
    /* Parameter adjustments */
    table_dim1 = *nrow;
    table_offset = 1 + table_dim1;
    table -= table_offset;
    --ico;
    --iro;
    --kyy;
    --idif;
    --irn;
    --key;
    --ipoin;
    --stp;
    --ifrq;
    --dlp;
    --dsp;
    --tm;
    --key2;
    --iwk;
    --rwk;

    /* Function Body */
    /* *********************************************************************** */
    /*                                  AMISS is a missing value indicator */
    /*                                  which is returned when the */
    /*                                  probability is not defined. */
    /* *********************************************************************** */
    /* *********************************************************************** */
    /*                                  TOL is chosen as the square root of */
    /*                                  the smallest relative spacing */
    /* *********************************************************************** */
    /* *********************************************************************** */
    /*                                  EMX is a large positive value used */
    /*                                  in comparing expected values */
    /* *********************************************************************** */
    /*                                  Initialize KEY array */
    i__1 = 2 * *ldkey;
    for (i__ = 1; i__ <= i__1; ++i__) {
        key[i__] = -9999;
        key2[i__] = -9999;
        /* L10: */
    }
    /*                                  Initialize parameters */
    *pre = 0.f;
    itop = 0;
    if (*expect > 0.) {
        emn = *emin;
    } else {
        emn = emx;
    }
    /*                                  Initialize pointers for workspace */
    k = Maximum(*nrow,*ncol);
    /*                                  f3xact */
    i31 = 1;
    i32 = i31 + k;
    i33 = i32 + k;
    i34 = i33 + k;
    i35 = i34 + k;
    i36 = i35 + k;
    i37 = i36 + k;
    i38 = i37 + k;
    i39 = i38 + 400;
    i310 = 1;
    i311 = 401;
    /*                                  f4xact */
    k = *nrow + *ncol + 1;
    i41 = 1;
    i42 = i41 + k;
    i43 = i42 + k;
    i44 = i43 + k;
    i45 = i44 + k;
    i46 = i45 + k;
    i47 = i46 + k * Maximum(*nrow,*ncol);
    i48 = 1;
    /*                                  Check table dimensions */
    /* if (*nrow > *ldtabl) {
        errMsg = errMsg & "NROW must be less than or equal to LDTABL.";
        WarnError (errMsg);
        return 0;
     }*/
    if (*ncol <= 1) {
        errMsg = errMsg & "NCOL must be greater than 1.0";
        WarnError (errMsg);
        return 0;
    }
    /*                                  Compute row marginals and total */
    ntot = 0;
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
        iro[i__] = 0;
        i__2 = *ncol;
        for (j = 1; j <= i__2; ++j) {
            if (table[i__ + j * table_dim1] < -1e-4) {
                errMsg = errMsg & "All elements of TABLE must be positive.";
                WarnError (errMsg);
                return 0;
            }
            iro[i__] += i_dnnt(&table[i__ + j * table_dim1]);
            ntot += i_dnnt(&table[i__ + j * table_dim1]);
            /* L20: */
        }
        /* L30: */
    }

    if (ntot == 0) {
        errMsg = errMsg & "All elements of TABLE are zero.  PRT and PRE are set to missing values (NaN, not a number)";
        ReportWarning (errMsg);
        return 0;
        *prt = amiss;
        *pre = amiss;
        goto L9000;
    }
    /*                                  Column marginals */
    i__1 = *ncol;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ico[i__] = 0;
        i__2 = *nrow;
        for (j = 1; j <= i__2; ++j) {
            ico[i__] += i_dnnt(&table[j + i__ * table_dim1]);
            /* L40: */
        }
        /* L50: */
    }
    /*                                  sort */
    isort_(nrow, &iro[1]);
    isort_(ncol, &ico[1]);
    /*                                  Determine row and column marginals */

    if (*nrow > *ncol) {
        nro = *ncol;
        nco = *nrow;
        /*                                  Interchange row and column marginals */
        i__1 = *nrow;
        for (i__ = 1; i__ <= i__1; ++i__) {
            itmp = iro[i__];
            if (i__ <= *ncol) {
                iro[i__] = ico[i__];
            }
            ico[i__] = itmp;
            /* L60: */
        }
    } else {
        nro = *nrow;
        nco = *ncol;
    }

    /*                                  Get multiplers for stack */
    kyy[1] = 1;
    i__1 = nro;
    for (i__ = 2; i__ <= i__1; ++i__) {
        /*                                  Hash table multipliers */
        if (iro[i__ - 1] + 1 <= imax / kyy[i__ - 1]) {
            kyy[i__] = kyy[i__ - 1] * (iro[i__ - 1] + 1);
            j /= kyy[i__ - 1];
        } else {
            errMsg = errMsg & "The hash table key cannot be computed because the largest key is larger than the largest representable integer.  The algorithm cannot proceed.";
            ReportWarning (errMsg);
            return 0;
        }
        /* L70: */
    }
    /*                                  Maximum product */
    if (iro[nro - 1] + 1 <= imax / kyy[nro - 1]) {
        kmax = (iro[nro] + 1) * kyy[nro - 1];
    } else {
        errMsg = errMsg & "The hash table key cannot be computed because the largest key is larger than the largest representable integer.  The algorithm cannot proceed.";
        ReportWarning (errMsg);
        goto L9000;
    }
    /*                                  Compute log factorials */
    fact[0] = 0.;
    fact[1] = 0.;
    fact[2] = log(2.);
    i__1 = ntot;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
        fact[i__] = fact[i__ - 1] + log((double) i__);
        j = i__ + 1;
        if (j <= ntot) {
            fact[j] = fact[i__] + fact[2] + fact[j / 2] - fact[j / 2 - 1];
        }
        /* L80: */
    }
    /*                                  Compute observed path length: OBS */
    obs = tol;
    ntot = 0;
    i__1 = nco;
    for (j = 1; j <= i__1; ++j) {
        dd = 0.f;
        i__2 = nro;
        for (i__ = 1; i__ <= i__2; ++i__) {
            if (*nrow <= *ncol) {
                dd += fact[i_dnnt(&table[i__ + j * table_dim1])];
                ntot += i_dnnt(&table[i__ + j * table_dim1]);
            } else {
                dd += fact[i_dnnt(&table[j + i__ * table_dim1])];
                ntot += i_dnnt(&table[j + i__ * table_dim1]);
            }
            /* L90: */
        }
        obs = obs + fact[ico[j]] - dd;
        /* L100: */
    }
    /*                                  Denominator of observed table: DRO */
    dro = f9xact_(&nro, &ntot, &iro[1], fact);
    *prt = exp(obs - dro);
    /*                                  Initialize pointers */
    k = nco;
    last = *ldkey + 1;
    jkey = *ldkey + 1;
    jstp = *ldstp + 1;
    jstp2 = *ldstp * 3 + 1;
    jstp3 = (*ldstp << 2) + 1;
    jstp4 = *ldstp * 5 + 1;
    ikkey = 0;
    ikstp = 0;
    ikstp2 = *ldstp << 1;
    ipo = 1;
    ipoin[1] = 1;
    stp[1] = 0.f;
    ifrq[1] = 1;
    ifrq[ikstp2 + 1] = -1;

L110:
    kb = nco - k + 1;
    ks = 0;
    n = ico[kb];
    kd = nro + 1;
    kmax = nro;
    /*                                  IDIF is the difference in going to th */
    /*                                  daughter */
    i__1 = nro;
    for (i__ = 1; i__ <= i__1; ++i__) {
        idif[i__] = 0;
        /* L120: */
    }
    /*                                  Generate the first daughter */
L130:
    --kd;
    /* Computing MIN */
    i__1 = n, i__2 = iro[kd];
    ntot = Minimum(i__1,i__2);
    idif[kd] = ntot;
    if (idif[kmax] == 0) {
        --kmax;
    }
    n -= ntot;
    if (n > 0 && kd != 1) {
        goto L130;
    }
    if (n != 0) {
        goto L310;
    }

    k1 = k - 1;
    n = ico[kb];
    ntot = 0;
    i__1 = nco;
    for (i__ = kb + 1; i__ <= i__1; ++i__) {
        ntot += ico[i__];
        /* L140: */
    }
    /*                                  Arc to daughter length=ICO(KB) */
L150:
    i__1 = nro;
    for (i__ = 1; i__ <= i__1; ++i__) {
        irn[i__] = iro[i__] - idif[i__];
        /* L160: */
    }
    /*                                  Sort irn */
    if (k1 > 1) {
        if (nro == 2) {
            if (irn[1] > irn[2]) {
                ii = irn[1];
                irn[1] = irn[2];
                irn[2] = ii;
            }
        } else if (nro == 3) {
            ii = irn[1];
            if (ii > irn[3]) {
                if (ii > irn[2]) {
                    if (irn[2] > irn[3]) {
                        irn[1] = irn[3];
                        irn[3] = ii;
                    } else {
                        irn[1] = irn[2];
                        irn[2] = irn[3];
                        irn[3] = ii;
                    }
                } else {
                    irn[1] = irn[3];
                    irn[3] = irn[2];
                    irn[2] = ii;
                }
            } else if (ii > irn[2]) {
                irn[1] = irn[2];
                irn[2] = ii;
            } else if (irn[2] > irn[3]) {
                ii = irn[2];
                irn[2] = irn[3];
                irn[3] = ii;
            }
        } else {
            i__1 = nro;
            for (j = 2; j <= i__1; ++j) {
                i__ = j - 1;
                ii = irn[j];
L170:
                if (ii < irn[i__]) {
                    irn[i__ + 1] = irn[i__];
                    --i__;
                    if (i__ > 0) {
                        goto L170;
                    }
                }
                irn[i__ + 1] = ii;
                /* L180: */
            }
        }
        /*                                  Adjust start for zero */
        i__1 = nro;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if (irn[i__] != 0) {
                goto L200;
            }
            /* L190: */
        }
L200:
        nrb = i__;
        nro2 = nro - i__ + 1;
    } else {
        nrb = 1;
        nro2 = nro;
    }
    /*                                  Some table values */
    ddf = f9xact_(&nro, &n, &idif[1], fact);
    drn = f9xact_(&nro2, &ntot, &irn[nrb], fact) - dro + ddf;
    /*                                  Get hash value */
    if (k1 > 1) {
        kval = irn[1] + irn[2] * kyy[2];
        i__1 = nro;
        for (i__ = 3; i__ <= i__1; ++i__) {
            kval += irn[i__] * kyy[i__];
            /* L210: */
        }
        /*                                  Get hash table entry */
        i__ = kval % (*ldkey << 1) + 1;
        /*                                  Search for unused location */
        i__1 = *ldkey << 1;
        for (itp = i__; itp <= i__1; ++itp) {
            ii = key2[itp];
            if (ii == kval) {
                goto L240;
            } else if (ii < 0) {
                key2[itp] = kval;
                dlp[itp] = 1.;
                dsp[itp] = 1.;
                goto L240;
            }
            /* L220: */
        }

        i__1 = i__ - 1;
        for (itp = 1; itp <= i__1; ++itp) {
            ii = key2[itp];
            if (ii == kval) {
                goto L240;
            } else if (ii < 0) {
                key2[itp] = kval;
                dlp[itp] = 1.f;
                goto L240;
            }
            /* L230: */
        }

        errMsg = errMsg & "LDKEY is too small.  It is not possible to give thevalue of LDKEY required, but you could try doubling LDKEY (and possibly LDSTP).";
        ReportWarning (errMsg);
        return 0;
    }

L240:
    ipsh = true;
    /*                                  Recover pastp */
    ipn = ipoin[ipo + ikkey];
    pastp = stp[ipn + ikstp];
    ifreq = ifrq[ipn + ikstp];
    /*                                  Compute shortest and longest path */
    if (k1 > 1) {
        obs2 = obs - fact[ico[kb + 1]] - fact[ico[kb + 2]] - ddf;
        i__1 = k1;
        for (i__ = 3; i__ <= i__1; ++i__) {
            obs2 -= fact[ico[kb + i__]];
            /* L250: */
        }

        if (dlp[itp] > 0.) {
            dspt = obs - obs2 - ddf;
            /*                                  Compute longest path */
            dlp[itp] = 0.;
            f3xact_(&nro2, &irn[nrb], &k1, &ico[kb + 1], &dlp[itp], &ntot,
                    fact, &iwk[i31], &iwk[i32], &iwk[i33], &iwk[i34], &iwk[
                        i35], &iwk[i36], &iwk[i37], &iwk[i38], &iwk[i39], &rwk[
                        i310], &rwk[i311], &tol);
            /* Computing Minimum */
            d__1 = 0., d__2 = dlp[itp];
            dlp[itp] = Minimum(d__1,d__2);
            /*                                  Compute shortest path */
            dsp[itp] = dspt;
            f4xact_(&nro2, &irn[nrb], &k1, &ico[kb + 1], &dsp[itp], fact, &
                    iwk[i47], &iwk[i41], &iwk[i42], &iwk[i43], &iwk[i44], &
                    iwk[i45], &iwk[i46], &rwk[i48], &tol);
            /* Computing Minimum */
            d__1 = 0., d__2 = dsp[itp] - dspt;
            dsp[itp] = Minimum(d__1,d__2);
            /*                                  Use chi-squared approximation? */
            if ((double) (irn[nrb] * ico[kb + 1]) / (double) ntot >
                    emn) {
                ncell = 0.f;
                i__1 = nro2;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    i__2 = k1;
                    for (j = 1; j <= i__2; ++j) {
                        if ((double) (irn[nrb + i__ - 1] * ico[kb + j]) >=
                                ntot * *expect) {
                            ++ncell;
                        }
                        /* L260: */
                    }
                    /* L270: */
                }
                if ((double) (ncell * 100) >= k1 * nro2 * *percnt) {
                    tmp = 0.f;
                    i__1 = nro2;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        tmp = tmp + fact[irn[nrb + i__ - 1]] - fact[irn[nrb +
                                i__ - 1] - 1];
                        /* L280: */
                    }
                    tmp *= k1 - 1;
                    i__1 = k1;
                    for (j = 1; j <= i__1; ++j) {
                        tmp += (nro2 - 1) * (fact[ico[kb + j]] - fact[ico[kb
                                             + j] - 1]);
                        /* L290: */
                    }
                    df = (double) ((nro2 - 1) * (k1 - 1));
                    tmp += df * 1.83787706640934548356065947281;
                    tmp -= (nro2 * k1 - 1) * (fact[ntot] - fact[ntot - 1]);
                    tm[itp] = (obs - dro) * -2. - tmp;
                } else {
                    /*                                  tm(itp) set to a flag value */
                    tm[itp] = -9876.;
                }
            } else {
                tm[itp] = -9876.;
            }
        }
        obs3 = obs2 - dlp[itp];
        obs2 -= dsp[itp];
        if (tm[itp] == -9876.) {
            chisq = false;
        } else {
            chisq = true;
            tmp = tm[itp];
        }
    } else {
        obs2 = obs - drn - dro;
        obs3 = obs2;
    }
    /*                                  Process node with new PASTP */
L300:
    if (pastp <= obs3) {
        /*                                  Update pre */
        *pre += (double) ifreq * exp(pastp + drn);

    } else if (pastp < obs2) {
        if (chisq) {
            df = (double) ((nro2 - 1) * (k1 - 1));
            /* Computing MAX */
            d__2 = 0., d__3 = tmp + (pastp + drn) * 2.;
            d__1 = Maximum(d__2,d__3) / 2.;
            d__4 = df / 2.;
            pv = 1.f - gammds_(&d__1, &d__4, &ifault);
            *pre += (double) ifreq * exp(pastp + drn) * pv;
        } else {
            /*                                  Put daughter on queue */
            d__1 = pastp + ddf;
            f5xact_(&d__1, &tol, &kval, &key[jkey], ldkey, &ipoin[jkey], &stp[
                        jstp], ldstp, &ifrq[jstp], &ifrq[jstp2], &ifrq[jstp3], &
                    ifrq[jstp4], &ifreq, &itop, &ipsh);
            ipsh = false;
        }
    }
    /*                                  Get next PASTP on chain */
    ipn = ifrq[ipn + ikstp2];
    if (ipn > 0) {
        pastp = stp[ipn + ikstp];
        ifreq = ifrq[ipn + ikstp];
        goto L300;
    }
    /*                                  Generate a new daughter node */
    f7xact_(&kmax, &iro[1], &idif[1], &kd, &ks, &iflag);
    if (iflag != 1) {
        goto L150;
    }
    /*                                  Go get a new mother from stage K */
L310:
    iflag = 1;
    f6xact_(&nro, &iro[1], &iflag, &kyy[1], &key[ikkey + 1], ldkey, &last, &
            ipo);
    /*                                  Update pointers */
    if (iflag == 3) {
        --k;
        itop = 0;
        ikkey = jkey - 1;
        ikstp = jstp - 1;
        ikstp2 = jstp2 - 1;
        jkey = *ldkey - jkey + 2;
        jstp = *ldstp - jstp + 2;
        jstp2 = (*ldstp << 1) + jstp;
        i__1 = *ldkey << 1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            key2[i__] = -9999;
            /* L320: */
        }
        if (k >= 2) {
            goto L310;
        }
    } else {
        goto L110;
    }

L9000:
    return 0;
} /* f2xact_ */

/* ----------------------------------------------------------------------- */
/*  Name:       FEXACT */

/*  Purpose:    Computes Fisher's exact test probabilities and a hybrid */
/*              approximation to Fisher exact test probabilities for a */
/*              contingency table using the network algorithm. */

/*  Usage:      CALL FEXACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT, */
/*                          EMIN, PRT, PRE) */

/*  Arguments: */
/*     NROW   - The number of rows in the table.  (Input) */
/*     NCOL   - The number of columns in the table.  (Input) */
/*     TABLE  - NROW by NCOL matrix containing the contingency table. */
/*              (Input) */
/*     LDTABL - Leading dimension of TABLE exactly as specified in the */
/*              dimension statement in the calling program.  (Input) */
/*     EXPECT - Expected value used in the hybrid algorithm for */
/*              deciding when to use asymptotic theory probabilities. */
/*              (Input) */
/*              If EXPECT .LE. 0.0 then asymptotic theory probabilities */
/*              are not used and Fisher exact test probabilities are */
/*              computed.  Otherwise, if PERCNT or more of the cells in */
/*              the remaining table have estimated expected values of */
/*              EXPECT or more, with no remaining cell having expected */
/*              value less than EMIN, then asymptotic chi-squared */
/*              probabilities are used.  See the algorithm section of the */
/*              manual document for details.  Use EXPECT = 5.0 to obtain */
/*              the 'Cochran' condition. */
/*     PERCNT - Percentage of remaining cells that must have estimated */
/*              expected  values greater than EXPECT before asymptotic */
/*              probabilities can be used.  (Input) */
/*              See argument EXPECT for details.  Use PERCNT = 80.0 to */
/*              obtain the 'Cochran' condition. */
/*     EMIN   - Minimum cell estimated expected value allowed for */
/*              asymptotic chi-squared probabilities to be used.  (Input) */
/*              See argument EXPECT for details.  Use EMIN = 1.0 to */
/*              obtain the 'Cochran' condition. */
/*     PRT    - Probability of the observed table for fixed marginal */
/*              totals.  (Output) */
/*     PRE    - Table p-value.  (Output) */
/*              PRE is the probability of a more extreme table, where */
/*              'extreme' is in a probabilistic sense. */
/*              If EXPECT .LT. 0 then the Fisher exact probability */
/*              is returned.  Otherwise, an approximation to the */
/*              Fisher exact probability is computed based upon */
/*              asymptotic chi-squared probabilities for ``large'' */
/*              table expected values.  The user defines ``large'' */
/*              through the arguments EXPECT, PERCNT, and EMIN. */

/*  Remarks: */
/*  1. For many problems one megabyte or more of workspace can be */
/*     required.  If the environment supports it, the user should begin */
/*     by increasing the workspace used to 200,000 units. */

/*  2. In FEXACT, LDSTP = 30*LDKEY.  The proportion of table space used */
/*     by STP may be changed by changing the line MULT = 30 below to */
/*     another value. */

/*  3. FEXACT may be converted to single precision by setting IREAL = 3, */
/*     and converting all DOUBLE PRECISION specifications (except the */
/*     specifications for RWRK, IWRK, and DWRK) to REAL.  This will */
/*     require changing the names and specifications of the intrinsic */
/*     functions ALOG, AMAX1, AMIN1, EXP, and REAL.  In addition, the */
/*     machine specific constants will need to be changed, and the name */
/*     DWRK will need to be changed to RWRK in the call to F2XACT. */

/*  4. Machine specific constants are specified and documented in F2XACT. */
/*     A missing value code is specified in both FEXACT and F2XACT. */

/*  5. Although not a restriction, is is not generally practical to call */
/*     this routine with large tables which are not sparse and in */
/*     which the 'hybrid' algorithm has little effect.  For example, */
/*     although it is feasible to compute exact probabilities for the */
/*     table */
/*            1 8 5 4 4 2 2 */
/*            5 3 3 4 3 1 0 */
/*           10 1 4 0 0 0 0, */
/*     computing exact probabilities for a similar table which has been */
/*     enlarged by the addition of an extra row (or column) may not be */
/*     feasible. */
/* ----------------------------------------------------------------------- */
int fexact_(long nrow, long ncol, double *table, double expect, double percnt, double emin, double *prt, double *pre)
{
    long ntot;

    _String errMsg ("Fisher Exact:");

    /* *********************************************************************** */
    /*                                  Set IREAL = 4 for DOUBLE PRECISION */
    /*                                  Set IREAL = 3 for SINGLE PRECISION */
    /* *********************************************************************** */

    /* if (*nrow > *ldtabl) {
        errMsg = errMsg &  "NROW must be less than or equal to LDTABL.";
        WarnError(errMsg);
        free (equiv_1);
        return 0;
     }*/

    ntot = 0;
    for (long i = 0; i<ncol*nrow; i++) {
        if (table[i] < 0.) {
            errMsg = errMsg &  "All elements of TABLE must be non-negative.";
            WarnError(errMsg);
            return 0;
        }
        ntot += (long)(table[i]+0.5);
    }


    if (ntot == 0) {
        errMsg = errMsg & "All elements of TABLE are zero.  PRT and PRE are set to missing values (NaN, not a number).";
        ReportWarning(errMsg);
        *prt = -1.;
        *pre = -1.;
        return 0;
    }

    long    k       = nrow + ncol + 1,
            kk   = k * ncol;


    double * i1   = (double*)MemAllocate((ntot + 1)*sizeof(double)),
             * irwk = (double*)MemAllocate(Maximum(ncol + 401,k)*sizeof(double));

    bool    didAlloc = false;

    long    *i2     = (long*) MemAllocate  (ncol*sizeof (long)),
             *i3      = (long*) MemAllocate  (ncol*sizeof (long)),
              *i3a  = (long*) MemAllocate  (ncol*sizeof (long)),
               *i3b   = (long*) MemAllocate  (nrow*sizeof (long)),
                *i3c    = (long*) MemAllocate  (nrow*sizeof (long)),
                 *iiwk    = (long*) MemAllocate  (Maximum(k * 5 + (kk << 1),ncol * 7 + 800)*sizeof (long));


    if (!fexact_i4) {
        didAlloc = true;
        allocate_fexact_keys(4096,30);
    }


    f2xact_ (&nrow, &ncol, table, &expect, &percnt, &emin, prt, pre, i1, i2, i3, i3a,
             i3b, i3c, fexact_i4, &fexact_ldkey, fexact_i5, fexact_i6, &fexact_ldstp, fexact_i7, fexact_i8,
             fexact_i9, fexact_i9a, fexact_i10, iiwk, irwk);

    free (i1);
    free (i2);
    free (i3);
    free (i3a);
    free (i3b);
    free (i3c);

    free (irwk);
    free (iiwk);

    if   (didAlloc) {
        free_fexact_keys();
    }
    return 0;
} /* fexact_ */



