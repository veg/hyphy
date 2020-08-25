/*

HyPhy - Hypothesis Testing Using Phylogenies.

Copyright (C) 1997-now
Core Developers:
  Sergei L Kosakovsky Pond (spond@ucsd.edu)
  Art FY Poon    (apoon42@uwo.ca)
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

#ifndef     __POLY__
#define     __POLY__

#include "parser.h"

#define   POLY_DATA_INCREMENT 10
#define   GLOBAL_VARIABLE   1
#define   CATEGORY_VARIABLE 2
#define   RANDOM_VARIABLE   3


//__________________________________________________________________________________

class _PolynomialData : public BaseObj
{

public:

    _PolynomialData (void);
    _PolynomialData (long);
    _PolynomialData (_PolynomialData const&);
    _PolynomialData (long,long, hyFloat*);
    _PolynomialData const & operator = (_PolynomialData const&);

    virtual ~_PolynomialData ();

    virtual BaseObj*    makeDynamic(void) const;
    virtual void        Duplicate  (BaseRefConst);

    inline  hyFloat*         GetCoeff (void) {
        return theCoeff;
    }
    inline  hyFloat&         GetCoeff (long index) {
        return theCoeff[index];
    }

    long    *           GetTerm (long);
    long                GetNoTerms (void) {
        return actTerms;
    }
    void                AddTerm (long*, hyFloat);
    void                AddTerm (long*, hyFloat, long*, long);
    void                AddTerm (hyFloat);
    void                WriteTerm (long*,long);
    void                DeleteTerm (long);
    bool                IsFirstANumber (void);
    inline long         NumberOfTerms (void) {
        return actTerms;
    }
    long                SumOfPowers (long);
    long                WeightedSumOfPowers (long,hyFloat*);

    // temp!

    bool                checkMe (void);

    friend class _Polynomial;

    void                MultiplyTerms (long*, long*, long*);
    void                RaiseTerm     (long*, long);
    static  hyFloat  BinaryRaise   (hyFloat, long);
    static  void        RearrangeTerm (long*, long*, long*,long);
    char                CompareTerms  (long*, long*);
    char                CompareTerms  (long*, long*, long*, long);
    char                CompareTerms  (long*, long*, long*, long*, long, long);
    long                FindTerm      (long*, long*, long start = 0);
    void                ResortTerms   (long*);
    void                ChopTerms     (void);
    bool                checkTerm     (hyFloat, long);


protected:

    hyFloat*     theCoeff;
    long*           thePowers;
    long            numberVars, actTerms, allocTerms;

};

//__________________________________________________________________________________

class _Polynomial : public _MathObject {

public:

    _Polynomial             (void);
    _Polynomial             (_SimpleList&);
    _Polynomial             (_Polynomial const&);
    _Polynomial const &     operator = (_Polynomial const&);
    _Polynomial             (hyFloat);
    _Polynomial             (_Variable&);
    virtual                 ~_Polynomial ();
    virtual                 _MathObject* ExecuteSingleOp (long opCode, _List *arguments = nil, _hyExecutionContext* context = _hyDefaultExecutionContext, HBLObjectRef cache = nil);   // execute this operation with the list of Args

    virtual BaseObj*        makeDynamic(void) const;
    virtual void            Duplicate  (BaseRefConst);

    virtual _MathObject*    Add                 (_MathObject*, HBLObjectRef cache);
    virtual _MathObject*    Plus                (_MathObject*, bool subtract, HBLObjectRef cache);
    virtual _MathObject*    Sub                 (_MathObject*, HBLObjectRef cache);
    virtual _MathObject*    Raise               (_MathObject*,  HBLObjectRef cache);
    virtual _MathObject*    Minus               (HBLObjectRef cache);
    virtual _MathObject*    Mult                (_MathObject*, HBLObjectRef cache);
    virtual _MathObject*    Compute             (void);
    virtual bool            Equal               (_MathObject*);
    hyFloat                 ComputePolynomial   (void);

    hyFloat                 ComputeP            (hyFloat* , hyFloat* , long , long, long*, long*);
    _MathObject*            IsANumber           (bool = false);
    virtual  bool           IsObjectEmpty       (void);

    virtual unsigned long            ObjectClass (void) const {
        return POLYNOMIAL;
    }
    virtual hyFloat      Value (void) {
        return ComputePolynomial();
    }

    virtual BaseObj*        toStr (unsigned long = 0UL);
    void                    CheckTerm(void);

    virtual void            toFileStr (FILE*, unsigned long = 0UL);

    long                    GetNoVariables(void) const {
        return variableIndex.countitems();
    }
    
    _Variable*              GetIthVariable (unsigned long i) const {
        return LocateVar (variableIndex.get (i));
    }
    
    _PolynomialData*        GetTheTerms(void) const {
        return theTerms;
    }
    void                    SetTheTerms(_PolynomialData* td) {
        theTerms = td;
    }
    void                    SetCLists(_SimpleList& c1,_SimpleList& c2) {
        compList1.Duplicate(&c1);
        compList2.Duplicate(&c2);
    }
    virtual void            ScanForVariables (_AVLList &l, bool globals = false, _AVLListX* tagger = nil, long weight = 0) const;
    virtual bool            HasChanged (bool = false);
    friend  void            ResetPolynomialCheck
    (_Polynomial*);
    long                    ComputationalSize (void) {
        return compList1.countitems();
    }
    bool                    IsMaxElement    (hyFloat);
    void                    Convert2ComputationForm
    (_SimpleList *c1 = nil, _SimpleList *c2 = nil, _SimpleList* termsToInclude = nil);
    void                    RankTerms       (_SimpleList*);
        
protected:

    void                DropSmallTerms(void);
    void                Convert2OperationForm(void);

    _SimpleList         variableIndex,
                        compList1,
                        compList2;

    _PolynomialData *   theTerms;


};

extern hyFloat dropPrecision, topPolyCap, dropTerms, enforcePolyCap,
       maximumPolyTermsPerVariable, maxPolynomialExpIterates,polynomialExpPrecision;
void    SetPolyTermCap (long);
#endif
