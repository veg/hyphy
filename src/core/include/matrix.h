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

#ifndef     __MATRIX__
#define     __MATRIX__

#include "hy_strings.h"
#include "avllistx.h"
#include "variablecontainer.h"
#include "trie.h"

#define     _POLYNOMIAL_TYPE 0
#define     _NUMERICAL_TYPE  1
#define     _FORMULA_TYPE 2
#define     _SIMPLE_FORMULA_TYPE 3



//_____________________________________________________________________________________________

#define      _HY_MATRIX_RANDOM_DIRICHLET         01L
#define      _HY_MATRIX_RANDOM_GAUSSIAN          02L
#define      _HY_MATRIX_RANDOM_WISHART           03L
#define      _HY_MATRIX_RANDOM_INVERSE_WISHART   04L
#define      _HY_MATRIX_RANDOM_MULTINOMIAL       05L

extern        _Trie        _HY_MatrixRandomValidPDFs;

//_____________________________________________________________________________________________

class _Formula;
/*__________________________________________________________________________________________________________________________________________ */

struct      _CompiledMatrixData {

    _SimpleFormulaDatum * theStack,
                        * varValues;

    hyFloat         * formulaValues;

    long      * formulaRefs;
    bool        has_volatile_entries;

    _SimpleList varIndex,
                formulasToEval;

};

/*__________________________________________________________________________________________________________________________________________ */

class       _Matrix: public _MathObject {
    
// data members
public:
    hyFloat   *theData;                            // matrix elements
    
protected:
    
    // data
    
    long        hDim, vDim, lDim;               // matrix physical dimensions; lDim - length of
    // actual storage allocated
    
    long*       theIndex;                       // indices of matrix elements in logical storage
    long*       compressedIndex;
        /**
            20200821 SLKP to speed sparse caclulations CompressSparseMatrix will create this DENSE index, which has the following structure
                - Entries in theIndex are expected to be compressed (no -1) and arranged BY ROW
                - First hDim values: the INDEX of the non-first non zero index for this ROW in theData
                - Next lDim values: the COLUMN index for the corresponding entry
         
                [x,1,2,x]
                [1,2,x,x]
                [x,x,2,x]
                [x,x,x,3]
         
                theIndex : [1,2,4,5,10,15]
                compressedIndex: [0,2,4,5,1,2,0,1,2,3]
                        
        */
    
private:
    
    long        storageType,                    // true element matrix(1) or a matrix of pointers(0) which do not need to be deleted
    // 2, if elements of the matrix are actually formulas which must be initialized to numerics before use
    // matrices of type two are merely storage tables and can not be operated on directly, i.e their
    // numerical values are computed first
                bufferPerRow,                   // values reflecting internal storage structure for
                overflowBuffer,                 // sparse matrices
                allocationBlock;
    
    _CompiledMatrixData*
    cmd;
    
    HBLObjectRef   theValue;                       // stores evaluated values of the matrix

    static      int     storageIncrement,       // how many percent of full matrix size
    // to allocate to the matrix storage per increment
    
                        precisionArg,                    // how many elements in exp series to truncate after
    
                        switchThreshold;                 // maximum percentage of non-zero elements
    // to keep the matrix sparse
    
    static      hyFloat truncPrecision;
    
    // matrix exp truncation precision
    
 

public:

    // constructors

    _Matrix ();                                 // default constructor, doesn't do much

    _Matrix (_String const&, bool, _FormulaParsingContext&, bool use_square_brackets = false);
    // matrix from a string of the form
    // {{i11,i12,...}{i21,i22,..}{in1,in2..)})
    // or {# rows,<# cols>{i1,j1,expr}{i2,j2,expr}..}
    // elements may be arbitrary formulas

    _Matrix (long theHDim, long theVDim, bool sparse = false, bool allocateStorage = false);    // create an empty matrix of given dimensions;



    // creates an empty matrix of given dimensions;
    // the first flag specifies whether it is sparse or not
    // the second is the storage type -- see below

//  _Matrix (long, long, int, bool);
    // a test function which generates a random matrix of given dimensions
    // where the third parameter specifies the percentage of 0 entries and
    // the first flag indicates how to store the matrix: as spars or usual

    _Matrix ( _Matrix const &);                       //duplicator

    _Matrix ( _SimpleList const &, long = -1);        // make matrix from simple list
    // the optional argument C (if > 0) tells HyPhy
    // to make a matrix with C elements per row
    // if <= 0 - a row matrix is returned

    _Matrix ( _List const &, bool parse_escapes = true);                         //make string matrix from a list

    _Matrix (const hyFloat *, unsigned long, unsigned long);
    /*
        20110616 SLKP
            added a simple constructor from a list of floating point values
            first argument: the values (must be at least 2nd arg x 3rd arg long)
            second argument: how many rows
            second argument: how many columns

    */
    _Matrix (hyFloat, unsigned long, unsigned long);
    /**
      20180919 SLKP
          make an rxc matrix that is constant (each cell is the same)
    */

    ~_Matrix (void);                            //destructor

    virtual void    Clear (bool complete = true);               //deletes all the entries w/o destroying the matrix
    virtual void    ZeroNumericMatrix (void);               //deletes all the entries w/o destroying the matrix

    void    Initialize (bool = true);                  // zeros all matrix structures

    virtual void        Serialize (_StringBuffer&,_String&);
    // write the matrix definition in HBL

    //_____________________________________________________________________________________________
    
    

    inline bool        is_empty (void) const {return GetVDim () == 0UL || GetHDim () == 0UL;}
    inline bool        is_row (void) const { return GetHDim() == 1UL;}
    inline bool        is_column (void) const { return GetVDim() == 1UL;}
    inline bool        is_square (void) const { return GetVDim() == GetHDim();}
    inline bool        is_dense (void) const {return theIndex == nil;}
    inline bool        is_expression_based (void) const {return storageType == _FORMULA_TYPE;}
    inline bool        is_numeric (void) const {return storageType == _NUMERICAL_TYPE;}
    inline bool        is_polynomial (void) const {return storageType == _POLYNOMIAL_TYPE;}
    inline bool        has_type (int t) const {return storageType == t;}

    HBLObjectRef           Evaluate (bool replace = true); // evaluates the matrix if contains formulas
    // if replace is true, overwrites the original

    virtual HBLObjectRef   ExecuteSingleOp (long opCode, _List* arguments = nil, _hyExecutionContext* context = _hyDefaultExecutionContext, HBLObjectRef cache = nil);
    // execute this operation with the list of Args

    HBLObjectRef   MAccess (HBLObjectRef, HBLObjectRef, HBLObjectRef cache = nil);
    // implements the M[i][j] operation for formulas
    HBLObjectRef   MCoord (HBLObjectRef, HBLObjectRef, HBLObjectRef cache = nil);
    // implements the M[i][j] operation for formulas

    void        MStore (long, long, _Formula&, long = -1);
    void        MStore (long, long, HBLObjectRef, long);
    bool        MResolve (HBLObjectRef, HBLObjectRef, long&, long&);
    // resolve coordiates from two Number arguments

    bool        CheckCoordinates ( long&, long&);
    // validate matrix coordinates

    bool        ValidateFormulaEntries (bool (long, long, _Formula*));
    // validate matrix coordinates

    void        MStore (HBLObjectRef, HBLObjectRef, _Formula&, long = HY_OP_CODE_NONE);
    // implements the M[i][j]= operation for formulas
    /*
        20100811: the last argument provides an op code (-1 = none)
        to perform on the _Formula argument and the current value in the matrix;
        this only applies to constant _Formula arguments

        e.g. passing HY_OP_CODE_ADD implements +=
     */

    void        MStore (long, long, HBLObjectRef);
    void        MStore (HBLObjectRef, HBLObjectRef, HBLObjectRef);
    // implements the M[i][j]= operation for objects
    virtual HBLObjectRef   Compute (void);         // returns the numeric value of this matrix

    virtual HBLObjectRef   ComputeNumeric (bool = false);  // returns the numeric value of this matrix
    virtual HBLObjectRef   RetrieveNumeric (void); // returns the numeric value of this matrix

    virtual void        ScanForVariables  (_AVLList&, bool inclG = false, _AVLListX* tagger = nil,long weight = 0) const;
    virtual void        ScanForVariables2 (_AVLList&, bool inclG = false, long modelID = -1, bool inclCat = true, _AVLListX* tagger = nil,long weight = 0) const;
    // scans for all local independent variables on which the matrix depends
    // and stores them in the list passed as the parameter

    virtual bool        IsIndependent (void)   {
        return storageType!=_FORMULA_TYPE;
    }
    // used to determine whether the matrix contains references
    // to other unknowns

    virtual unsigned long        ObjectClass (void) const      {
        return MATRIX;
    }

    _Matrix const&     operator = (_Matrix const&);             // assignment operation on matrices
    _Matrix const&     operator = (_Matrix const*);             // assignment operation on matrices with temp results

    virtual HBLObjectRef    Random (HBLObjectRef, HBLObjectRef cache);    // reshuffle the matrix

    virtual HBLObjectRef    AddObj (HBLObjectRef, HBLObjectRef cache);    // addition operation on matrices

    virtual HBLObjectRef    SubObj (HBLObjectRef, HBLObjectRef cache);    // subtraction operation on matrices

    virtual HBLObjectRef    MultObj (HBLObjectRef, HBLObjectRef cache);   // multiplication operation on matrices

    virtual HBLObjectRef    MultElements (HBLObjectRef, bool elementWiseDivide, HBLObjectRef cache);  // element wise multiplication/division operation on matrices

    virtual HBLObjectRef    Sum          (HBLObjectRef cache);

    _Matrix     operator + (_Matrix&);          // addition operation on matrices

    _Matrix     operator - (_Matrix&);          // subtraction operation on matrices

    _Matrix     operator * (_Matrix&);          // multiplication operation on matrices

    _Matrix     operator * (hyFloat);        // multiplication of a matrix by a number

    void        operator += (_Matrix&);         // add/store operation on matrices

    void        operator -= (_Matrix&);         // subtract/store operation on matrices

    void        operator *= (_Matrix&);         // multiply/store operation on matrices

    void        operator *= (hyFloat);       // multiply by a #/store operation on matrices

    void        AplusBx  (_Matrix&, hyFloat); // A = A + B*x (scalar)

    hyFloat        Sqr         (hyFloat* _hprestrict_);
    // square the matrix; takes a scratch vector
    // of at least lDim doubles
    // return the maximum absolute element-wise difference between X and X^2

    _List*      ComputeRowAndColSums            (void);
    _Matrix*    MutualInformation               (void);
    void        FillInList                      (_List&, bool convert_numbers = false) const;
    // SLKP 20101108:
    //               added a boolean flag to allow numeric matrices
    //               to be implicitly converted to strings

    _Matrix*    NeighborJoin                    (bool, HBLObjectRef cache);
    _Matrix*    MakeTreeFromParent              (long, HBLObjectRef cache);
    _Matrix*    ExtractElementsByEnumeration    (_SimpleList*,_SimpleList*,bool=false);
    _Matrix*    SimplexSolve                    (hyFloat = 1.e-6);


//  void        SqrStrassen (void);
// square the matrix by Strassen's Multiplication


    _Matrix*    Exponentiate (hyFloat scale_to = 1.0, bool check_transition = false, _Matrix * write_here = nil);                // exponent of a matrix
    void        Transpose (void);                   // transpose a matrix
    _Matrix     Gauss   (void);                     // Gaussian Triangularization process
    HBLObjectRef   LUDecompose (void) const;
    HBLObjectRef   CholeskyDecompose (void) const;
    // added by afyp July 6, 2009
    HBLObjectRef   Eigensystem (HBLObjectRef) const;
    HBLObjectRef   LUSolve (HBLObjectRef) const;
    HBLObjectRef   Inverse (HBLObjectRef cache) const;
    HBLObjectRef   Abs (HBLObjectRef cache);                     // returns the norm of a matrix
    // if it is a vector - returns the Euclidean length
    // otherwise returns the largest element

    hyFloat  AbsValue                        (void) const;
    
    template <typename CALLBACK>  HBLObjectRef ApplyScalarOperation (CALLBACK && functor, HBLObjectRef cache) const;
    
    // return the matrix of logs of every matrix element
    
    void        SwapRows (const long, const long);
    long        CompareRows (const long, const long);

    hyFloat  operator () (long, long) const;       // read access to an element in a matrix
    hyFloat& operator [] (long);             // read/write access to an element in a matrix
  
    hyFloat& get_dense_numeric_cell (unsigned long r, unsigned long c) {
        return theData[r*vDim + c];
    }

    void        Store               (long, long, hyFloat);                       // write access to an element in a matrix
    void        StoreObject         (long, long, _MathObject*, bool dup = false);
    void        StoreObject         (long,  _MathObject*,bool dup = false);
    void        StoreFormula        (long, long, _Formula&, bool = true, bool = true);
    void        NonZeroEntries      (_SimpleList&);

    void        UpdateDiag          (long ,long , _MathObject*);

    void        Swap                (_Matrix&);         // fast swap matrix data
    friend      void                SetIncrement (int); // storage parameter access
    // an auxiliary function which creates an empty
    // matrix of given dimensions and storage class (normal/sparse)
    // and storage type (pointer/array)

    friend      void                DuplicateMatrix (_Matrix*,  _Matrix const*);
    // an auxiliary function which duplicates a matrix


    hyFloat          MaxElement      (char doSum = 0, long * = nil) const;
    // SLKP 20101022: added an option to return the sum of all elements as an option (doSum = 1) or
    // the sum of all absolute values (doSum == 2)
    // returns the largest element's abs value for given matrix
    // SLKP 20110523: added an option (doSum == 3) to return the largest element (no abs value)
    // for run modes 0 and 3, if the 2nd argument is non-nil, the index of the winning element will be stored

    hyFloat          MinElement      (char doAbs = 1, long * = nil);

    // SLKP 20110620: added the 2nd argument to optionally store the index of the smallest element
    //              : added the option to NOT do absolute values
    // returns the smallest, non-zero element value for given matrix

    bool                IsMaxElement    (hyFloat);
    // true if there is an elem abs val of which is greater than the arg
    // false otherwise


    hyFloat              MaxRelError(_Matrix&);

//  friend      _Matrix     IterateStrassen (_Matrix&, _Matrix&);
    // used in Strassen Squaring

    virtual     BaseRef     makeDynamic (void) const; // duplicate this object into a dynamic copy

    virtual     void        Duplicate   (BaseRefConst obj); // duplicate this object into a dynamic copy

    virtual     BaseRef     toStr       (unsigned long = 0UL);       // convert this matrix to a string

    virtual     void        toFileStr   (FILE*dest, unsigned long = 0UL);

    bool        AmISparse               (void);

    hyFloat  ExpNumberOfSubs         (_Matrix*,bool);

    virtual     bool        IsVariable  (void) {
        return storageType != _NUMERICAL_TYPE;
    }
    // is this matrix a constant or a variable quantity?

    virtual     bool        IsConstant  (void);

    virtual     bool        IsPrintable (void) {
        return storageType != _FORMULA_TYPE;
    }
    
    virtual     bool        Equal       (HBLObjectRef);

    void        ExportMatrixExp         (_Matrix*, FILE*);
    bool        ImportMatrixExp         (FILE*);

    hyFloat  FisherExact             (hyFloat, hyFloat, hyFloat);

    virtual     bool        HasChanged  (bool = false);
    // have any variables which are referenced by the elements changed?

    virtual     unsigned long
    GetHDim                     (void) const{
        return hDim;
    }
    
    bool     check_dimension                         (unsigned long rows, unsigned long columns) const {
        return hDim == rows && vDim == columns;
    }
    
    unsigned long        GetVDim                     (void) const {
        return vDim;
    }
    unsigned long        GetSize                     (void) const {
        return lDim;
    }
    long        GetMySize                   (void) {
        return sizeof(_Matrix)+lDim*(storageType==_NUMERICAL_TYPE?sizeof(hyFloat):sizeof(hyPointer));
    }

    void        PopulateConstantMatrix      (hyFloat);
    /* SLKP 20090818
            fill out a numeric matrix with a fixed value
            if the matrix is sparse, only will out the non-void entries
     */

    _Formula*      GetFormula                  (long, long) const;
    HBLObjectRef   GetMatrixCell               (long, long, HBLObjectRef cache = nil) const;
    HBLObjectRef   MultByFreqs                 (long, bool = false);
    HBLObjectRef   EvaluateSimple              (_Matrix* existing_receptacle = nil);
    HBLObjectRef   SortMatrixOnColumn          (HBLObjectRef, HBLObjectRef cache);
    HBLObjectRef   K_Means                     (HBLObjectRef, HBLObjectRef cache);
    HBLObjectRef   pFDR                        (HBLObjectRef, HBLObjectRef cache);    // positive false discovery rate
    HBLObjectRef   PoissonLL                   (HBLObjectRef, HBLObjectRef cache);    // log likelihood of a vector of poisson samples given a parameter value


    // added by afyp, July 1, 2009
    HBLObjectRef   DirichletDeviate            (void);         // this matrix used for alpha hyperparameters
    HBLObjectRef   GaussianDeviate             (_Matrix &);    //  "   "   "   "       mean hyperparameter, additional argument for variance
    HBLObjectRef   InverseWishartDeviate       (_Matrix &);    //  "   "   "   "       rho hyperparameter, additional for phi matrix
    HBLObjectRef   WishartDeviate              (_Matrix &, _Matrix &),
                WishartDeviate                (_Matrix &);

    HBLObjectRef   MultinomialSample           (_Constant*);
    /* SLKP 20110208: an internal function to draw the multinomial sample

        the matrix _base_ must be 2xN, where each _row_ lists

            value (integer 0 to N-1, but not enforced), probability

        the function will normalize by sum of all the values in the second column

        the constant argument is the number of replicates (M) to draw

        returns an 1xN matrix with counts of how often each value has been drawn

     */


    bool        IsValidTransitionMatrix     () const;

    bool        IsReversible                (_Matrix* = nil);
    // check if the matrix is reversible
    // if given a base frequencies assumes that rate matrix entries will not be multiplied by freq terms

    bool        IsAStringMatrix             (void) const;
    void        MakeMeSimple                (void);
    void        MakeMeGeneral               (void);
    void        ConvertToSimpleList         (_SimpleList&);
    void        CompressSparseMatrix        (bool, hyFloat*);
    //prepare the transition probs matrix for exponentiation

    long        Hash (long, long) const;                  // hashing function, which finds matrix
    // physical element in local storage buffer
    // returns -1 if insufficient storage
    // returns a negative number
    // if element was not found, the number returned
    // indicates the first available slot

    hyFloat*       fastIndex(void)  const {
        return (!theIndex)&&(storageType==_NUMERICAL_TYPE)?(hyFloat*)theData:nil;
    }
    inline            hyFloat&         directIndex(long k)   {
        return theData[k];
    }
    long              MatrixType (void) {
        return storageType;
    }
    bool              SparseDataStructure (void) {
        return theIndex;
    }
    
    template <typename CALLBACK, typename EXTRACTOR>  void ForEach (CALLBACK&& cbv, EXTRACTOR&& accessor) const {
        if (theIndex) {
            for (unsigned long i=0UL; i<(unsigned long)lDim; i++) {
                if (theIndex[i] >= 0L) {
                    cbv (accessor (i), theIndex[i], i);
                }
            }
        } else {
            for (unsigned long i=0UL; i<(unsigned long)lDim; i++) {
                cbv (accessor (i), i, i);
            }
        }
    }

    template <typename CALLBACK> void ForEachCellNumeric (CALLBACK&& cbv) const {
        if (theIndex) {
            for (unsigned long i=0UL; i<(unsigned long)lDim; i++) {
                long idx = theIndex[i];
                if (idx >= 0L) {
                    long row = idx / vDim;
                    cbv (theData[i], idx, row, idx - row*vDim);
                }
            }
        } else {
            for (unsigned long i=0UL, c = 0UL; i<(unsigned long)hDim; i++) {
                for (unsigned long j=0UL; j<(unsigned long)vDim; j++, c++) {
                    cbv (theData[c], c, i, j);
                }
            }
        }
    }

    template <typename CALLBACK, typename EXTRACTOR>  bool Any (CALLBACK&& cbv, EXTRACTOR&& accessor) const {
        if (theIndex) {
            for (unsigned long i=0UL; i<(unsigned long)lDim; i++) {
                if (theIndex[i] >= 0L) {
                    if (cbv (accessor (i), theIndex[i])) {
                        return true;
                    }
                }
            }
        } else {
            for (unsigned long i=0UL; i<(unsigned long)lDim; i++) {
                if (cbv (accessor (i), i)) {
                    return true;
                }
            }
        }
        return true;
    }

    void              CheckIfSparseEnough (bool = false);       // check if matrix is sparse enough to justify
    void              Convert2Formulas      (void);     // converts a numeric matrix to formula-based mx
    // sparse storage

    void              Resize                (long);     // resize a dense numeric matrix to have more rows
  

    inline            hyFloat    get_direct        (long const index) const {
        return theData [index];
    }

    inline            hyFloat    get        (long const row, long const column) const {
      return theData [row * vDim + column];
    }
  
    inline            hyFloat&    set        (long const row, long const column)  {
      return theData [row * vDim + column];
    }

    _String*          BranchLengthExpression(_Matrix*, bool);

    void              CopyABlock                        (_Matrix*, long, long, long = 0, long = 0);
    /* starting at element (row -- 2nd argument, column -- 3rd argument)
       copy the source matrix (1st argument) row by row

       e.g. if this matrix is 3x4 and the source matrix is 2x2
       then copying from element 2,2 (0 - based as always)
       will result in

    [[ x x x x]
       [ x x x x]
     [ x x y y]]

       where y is used to denote an element copied from the source
     4th and 5th arguments override source.hDim and source.vDim,
       respectively, if they are positive

     Note that both matrices are ASSUMED to be numeric and dense

       NO ERROR CHECKING IS DONE!

     */
    /*---------------------------------------------------*/

    static void    CreateMatrix    (_Matrix* theMatrix, long theHDim, long theVDim,  bool sparse = false, bool allocateStorage = false, bool isFla = false);

    void        RecursiveIndexSort      (long, long, _SimpleList*);



private:

    void     internal_to_str (_StringBuffer*, FILE*, unsigned long padding);
    void     SetupSparseMatrixAllocations (void);
    bool     is_square_numeric   (bool dense = true) const;
    
    hyFloat  computePFDR         (hyFloat, hyFloat);
    void        InitMxVar           (_SimpleList&   , hyFloat);
    bool        ProcessFormulas     (long&, _AVLList&, _SimpleList&, _SimpleList&, _AVLListX&, bool = false, _Matrix* = nil);

    HBLObjectRef   PathLogLikelihood   (HBLObjectRef, HBLObjectRef cache);
    /* SLKP: 20100812

     This function assumes that 'this' an 3xK matrix, where each column is of the form
     A: integer in [0,N-1], B: integer in [0,N-1], T: real >= 0

     and the argument is an NxN RATE matrix for a Markov chain

     The return value is the \sum_ j = 0 ^ {K-1} Prob {A_j -> B_j | T_j}
     */

    _Matrix*    BranchLengthStencil (void) const;

    //bool      IsAStringMatrix     (void);
    void        AddMatrix           (_Matrix&, _Matrix&, bool sub = false);
    // aux arithmetic rountines
    bool        AddWithThreshold    (_Matrix&, hyFloat);
    void        RowAndColumnMax     (hyFloat&, hyFloat&, hyFloat* = nil);
    void        Subtract            (_Matrix&, _Matrix&);
    void        Multiply            (_Matrix&, hyFloat);
    void        Multiply            (_Matrix&, _Matrix const &) const;
    bool        IsNonEmpty          (long) const;
    // checks to see if the i-th position in the storage is non-empty
    bool        CheckDimensions     (_Matrix&) const;
    // compare dims of 2 matrices to see if they can be multiplied
    long        HashBack            (long) const;
    // hashing function, which finds matrix
    // physical element given local storage
    void        MultbyS             (_Matrix&,bool,_Matrix* = nil, hyFloat* = nil);
    // internal function used in exponentiating sparse matrices

    void        Balance             (void);  // perform matrix balancing; i.e. a norm reduction which preserves the eigenvalues
    // lifted from balanc function in NR

    void        Schur               (void);  // reduce the matrix to Hessenberg form preserving eigenvalues
    // lifted from elmhes function in NR

    void        EigenDecomp         (_Matrix&,_Matrix&) const;  // find the eigenvalues of a real matrix
    // return real and imaginary parts
    // lifted from hqr in NR


    bool        AmISparseFast       (_Matrix&);

    bool        IncreaseStorage     (void);
    // expand the buffer, where elements are held
    // returns TRUE/FALSE for success/failure


    void        BreakPoints             (long, long, _SimpleList*);
    void        ConvertFormulas2Poly    (bool force2numbers=true);
    void        ConvertNumbers2Poly     (void);
    void        AgreeObjects            (_Matrix&);

    void        ClearFormulae           (void);             // internal reuseable purger
    void        ClearObjects            (void);             // internal reuseable purger
    inline
    _MathObject*&
    GetMatrixObject         (long ind) const {
        return ((_MathObject**)theData)[ind];
    }
    inline
    bool        CheckObject             (long ind) const{
        return ((_MathObject**)theData)[ind]!=nil;
    }

    void        SimplexHelper1          (long, _SimpleList&, long, bool, long&, hyFloat&);
    void        SimplexHelper2          (long&, long, hyFloat);
    void        SimplexHelper3          (long,long,long,long);
    //  helper functions for SimplexSolver

    // if nil - matrix stored conventionally

};

/*__________________________________________________________________________________________________________________________________________ */


/*__________________________________________________________________________________________________________________________________________ */

extern  _Matrix *GlobalFrequenciesMatrix;
// the matrix of frequencies for the trees to be set by block likelihood evaluator
extern  long  ANALYTIC_COMPUTATION_FLAG;

HBLObjectRef _returnMatrixOrUseCache (long nrow, long ncol, long type, bool is_sparse, HBLObjectRef cache);

#ifdef  _SLKP_USE_AVX_INTRINSICS
    inline double _avx_sum_4 (__m256d const & x) {
      /*__m256d t = _mm256_add_pd (_mm256_shuffle_pd (x, x, 0x0),
                                 // (x3,x3,x1,x1)
                                 _mm256_shuffle_pd (x, x, 0xf)
                                 // (x2,x2,x0,x0);
                                 );
      return _mm_cvtsd_f64 (_mm_add_pd(
                                       _mm256_castpd256_pd128 (t), // (x3+x2,x3+x2)
                                       _mm256_extractf128_pd(t,1)  // (x1+x0,x0+x1);
                                       ));
       */
        __m256d sum      = _mm256_hadd_pd(x, x);
        // sum now has (x[0]+x[1],x[0]+x[1],x[2]+x[3],x[2]+x[3])
        return _mm_cvtsd_f64(_mm_add_pd(_mm256_extractf128_pd(sum, 1), _mm256_castpd256_pd128(sum)));
        /*double  __attribute__ ((aligned (32))) array[4];
        _mm256_store_pd (array, x);
        return (array[0]+array[1])+(array[2]+array[3])  ;*/
      
    }
#endif

#endif
