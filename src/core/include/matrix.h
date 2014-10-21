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

#define     HY_MATRIX_COLUMN_VECTOR     1
#define     HY_MATRIX_ROW_VECTOR        2

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

    _Parameter         * formulaValues;

    long      * formulaRefs;
    bool        has_volatile_entries;

    _SimpleList varIndex,
                formulasToEval;

};

/*__________________________________________________________________________________________________________________________________________ */

class       _Matrix: public _MathObject
{

public:

    // constructors

    _Matrix ();                                 // default constructor, doesn't do much

    _Matrix (_String&, bool = false, _VariableContainer* = nil);
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

    _Matrix ( _Matrix &);                       //duplicator

    _Matrix ( _SimpleList &, long = -1);        // make matrix from simple list
    // the optional argument C (if > 0) tells HyPhy
    // to make a matrix with C elements per row
    // if <= 0 - a row matrix is returned

    _Matrix ( _List &);                         //make string matrix from a list

    _Matrix (_Parameter *, unsigned long, unsigned long);
    /*
        20110616 SLKP
            added a simple constructor from a list of floating point values
            first argument: the values (must be at least 2nd arg x 3rd arg long)
            second argument: how many rows
            second argument: how many columns

    */

    ~_Matrix (void);                            //destructor

    virtual void    Clear (void);               //deletes all the entries w/o destroying the matrix

    void    Initialize (void);                  // zeros all matrix structures

    virtual void        Serialize (_String&,_String&);
    // write the matrix definition in HBL

    virtual bool        IsAVector (char = 0);   // is a vector? 0 - either row or column; 1 column; 2 row

    _PMathObj           Evaluate (bool replace = true); // evaluates the matrix if contains formulas
    // if replace is true, overwrites the original

    virtual _PMathObj   Execute (long opCode, _PMathObj p = nil, _PMathObj p2 = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    // execute this operation with the list of Args

    _PMathObj   MAccess (_PMathObj, _PMathObj);
    // implements the M[i][j] operation for formulas
    _PMathObj   MCoord (_PMathObj, _PMathObj);
    // implements the M[i][j] operation for formulas

    void        MStore (long, long, _Formula&, long = -1);
    bool        MResolve (_PMathObj, _PMathObj, long&, long&);
    // resolve coordiates from two Number arguments

    bool        CheckCoordinates ( long&, long&);
    // validate matrix coordinates

    void        MStore (_PMathObj, _PMathObj, _Formula&, long = HY_OP_CODE_NONE);
    // implements the M[i][j]= operation for formulas
    /*
        20100811: the last argument provides an op code (-1 = none)
        to perform on the _Formula argument and the current value in the matrix;
        this only applies to constant _Formula arguments

        e.g. passing HY_OP_CODE_ADD implements +=
     */

    void        MStore (long, long, _PMathObj);
    void        MStore (_PMathObj, _PMathObj, _PMathObj);
    // implements the M[i][j]= operation for objects
    virtual _PMathObj   Compute (void);         // returns the numeric value of this matrix

    virtual _PMathObj   ComputeNumeric (bool = false);  // returns the numeric value of this matrix
    virtual _PMathObj   RetrieveNumeric (void); // returns the numeric value of this matrix

    virtual void        ScanForVariables  (_AVLList&, bool inclG = false, _AVLListX* tagger = nil,long weight = 0);
    virtual void        ScanForVariables2 (_AVLList&, bool inclG = false, long modelID = -1, bool inclCat = true, _AVLListX* tagger = nil,long weight = 0);
    // scans for all local independent variables on which the matrix depends
    // and stores them in the list passed as the parameter

    virtual bool        IsIndependent (void)   {
        return storageType!=_FORMULA_TYPE;
    }
    // used to determine whether the matrix contains references
    // to other unknowns

    virtual unsigned long        ObjectClass (void)      {
        return MATRIX;
    }

    void     operator = (_Matrix&);             // assignment operation on matrices
    void     operator = (_Matrix*);             // assignment operation on matrices with temp results

    virtual _PMathObj    Random (_PMathObj);    // reshuffle the matrix

    virtual _PMathObj    AddObj (_PMathObj);    // addition operation on matrices

    virtual _PMathObj    SubObj (_PMathObj);    // subtraction operation on matrices

    virtual _PMathObj    MultObj (_PMathObj);   // multiplication operation on matrices

    virtual _PMathObj    MultElements (_PMathObj, bool elementWiseDivide = false);  // element wise multiplication/division operation on matrices

    virtual _PMathObj    Sum          (void);

    _Matrix     operator + (_Matrix&);          // addition operation on matrices

    _Matrix     operator - (_Matrix&);          // subtraction operation on matrices

    _Matrix     operator * (_Matrix&);          // multiplication operation on matrices

    _Matrix     operator * (_Parameter);        // multiplication of a matrix by a number

    void        operator += (_Matrix&);         // add/store operation on matrices

    void        operator -= (_Matrix&);         // subtract/store operation on matrices

    void        operator *= (_Matrix&);         // multiply/store operation on matrices

    void        operator *= (_Parameter);       // multiply by a #/store operation on matrices

    void        AplusBx  (_Matrix&, _Parameter); // A = A + B*x (scalar)

    void        Sqr         (_Parameter* _hprestrict_);
    // square the matrix; takes a scratch vector
    // of at least lDim doubles

    _List*      ComputeRowAndColSums            (void);
    _Matrix*    MutualInformation               (void);
    void        FillInList                      (_List&, bool = false);
    // SLKP 20101108:
    //               added a boolean flag to allow numeric matrices
    //               to be implicitly converted to strings

    _Matrix*    NeighborJoin                    (bool);
    _Matrix*    MakeTreeFromParent              (long);
    _Matrix*    ExtractElementsByEnumeration    (_SimpleList*,_SimpleList*,bool=false);
    _Matrix*    SimplexSolve                    (_Parameter = 1.e-6);


//  void        SqrStrassen (void);
// square the matrix by Strassen's Multiplication


    _Matrix*    Exponentiate (void);                // exponent of a matrix
    void        Transpose (void);                   // transpose a matrix
    _Matrix     Gauss   (void);                     // Gaussian Triangularization process
    _PMathObj   LUDecompose (void);
    _PMathObj   CholeskyDecompose (void);
    // added by afyp July 6, 2009
    _PMathObj   Eigensystem (void);
    _PMathObj   LUSolve (_PMathObj);
    _PMathObj   Inverse (void);
    _PMathObj   Abs (void);                     // returns the norm of a matrix
    // if it is a vector - returns the Euclidean length
    // otherwise returns the largest element

    _Parameter  AbsValue                        (void);
    virtual     _PMathObj Log                   (void);
    // return the matrix of logs of every matrix element
    
    void        SwapRows (const long, const long);
    long        CompareRows (const long, const long);

    _Parameter  operator () (long, long);       // read access to an element in a matrix
    _Parameter& operator [] (long);             // read/write access to an element in a matrix

    void        Store               (long, long, _Parameter);                       // write access to an element in a matrix
    void        StoreObject         (long, long, _MathObject*, bool dup = false);
    void        StoreObject         (long,  _MathObject*,bool dup = false);
    void        StoreFormula        (long, long, _Formula&, bool = true, bool = true);
    void        NonZeroEntries      (_SimpleList&);

    void        UpdateDiag          (long ,long , _MathObject*);

    void        Swap                (_Matrix&);         // fast swap matrix data
    friend      void                SetIncrement (int); // storage parameter access
    friend      void                CreateMatrix    (_Matrix*, long, long, bool, bool, bool);
    // an auxiliary function which creates an empty
    // matrix of given dimensions and storage class (normal/sparse)
    // and storage type (pointer/array)

    friend      void                DuplicateMatrix (_Matrix*,  _Matrix*);
    // an auxiliary function which duplicates a matrix


    _Parameter          MaxElement      (char doSum = 0, long * = nil);
    // SLKP 20101022: added an option to return the sum of all elements as an option (doSum = 1) or
    // the sum of all absolute values (doSum == 2)
    // returns the largest element's abs value for given matrix
    // SLKP 20110523: added an option (doSum == 3) to return the largest element (no abs value)
    // for run modes 0 and 3, if the 2nd argument is non-nil, the index of the winning element will be stored

    _Parameter          MinElement      (char doAbs = 1, long * = nil);

    // SLKP 20110620: added the 2nd argument to optionally store the index of the smallest element
    //              : added the option to NOT do absolute values
    // returns the smallest, non-zero element value for given matrix

    bool                IsMaxElement    (_Parameter);
    // true if there is an elem abs val of which is greater than the arg
    // false otherwise


    _Parameter              MaxRelError(_Matrix&);

//  friend      _Matrix     IterateStrassen (_Matrix&, _Matrix&);
    // used in Strassen Squaring

    virtual     BaseRef     makeDynamic (void); // duplicate this object into a dynamic copy

    virtual     void        Duplicate   (BaseRef obj); // duplicate this object into a dynamic copy

    virtual     BaseRef     toStr       (void);       // convert this matrix to a string

    virtual     void        toFileStr   (FILE*dest);

    bool        AmISparse               (void);

    _Parameter  ExpNumberOfSubs         (_Matrix*,bool);

    virtual     bool        IsVariable  (void) {
        return storageType != _NUMERICAL_TYPE;
    }
    // is this matrix a constant or a variable quantity?

    virtual     bool        IsConstant  (void);

    virtual     bool        IsPrintable (void) {
        return storageType != _FORMULA_TYPE;
    }
    
    virtual     bool        Equal       (_PMathObj);

    void        ExportMatrixExp         (_Matrix*, FILE*);
    bool        ImportMatrixExp         (FILE*);

    _Parameter  FisherExact             (_Parameter, _Parameter, _Parameter);

    virtual     bool        HasChanged  (void);
    // have any variables which are referenced by the elements changed?

    virtual     long
    GetHDim                     (void) {
        return hDim;
    }
    long        GetVDim                     (void) {
        return vDim;
    }
    long        GetSize                     (void) {
        return lDim;
    }
    long        GetMySize                   (void) {
        return sizeof(_Matrix)+lDim*(storageType==_NUMERICAL_TYPE?sizeof(_Parameter):sizeof(Ptr));
    }

    void        PopulateConstantMatrix      (const _Parameter);
    /* SLKP 20090818
            fill out a numeric matrix with a fixed value
            if the matrix is sparse, only will out the non-void entries
     */

    _Formula*   GetFormula                  (long, long);
    _PMathObj   MultByFreqs                 (long);
    _PMathObj   EvaluateSimple              (void);
    _PMathObj   SortMatrixOnColumn          (_PMathObj);
    _PMathObj   K_Means                     (_PMathObj);
    _PMathObj   pFDR                        (_PMathObj);    // positive false discovery rate
    _PMathObj   PoissonLL                   (_PMathObj);    // log likelihood of a vector of poisson samples given a parameter value


    // added by afyp, July 1, 2009
    _PMathObj   DirichletDeviate            (void);         // this matrix used for alpha hyperparameters
    _PMathObj   GaussianDeviate             (_Matrix &);    //  "   "   "   "       mean hyperparameter, additional argument for variance
    _PMathObj   InverseWishartDeviate       (_Matrix &);    //  "   "   "   "       rho hyperparameter, additional for phi matrix
    _PMathObj   WishartDeviate              (_Matrix &, _Matrix &),
                WishartDeviate                (_Matrix &);

    _PMathObj   MultinomialSample           (_Constant*);
    /* SLKP 20110208: an internal function to draw the multinomial sample

        the matrix _base_ must be 2xN, where each _row_ lists

            value (integer 0 to N-1, but not enforced), probability

        the function will normalize by sum of all the values in the second column

        the constant argument is the number of replicates (M) to draw

        returns an 1xN matrix with counts of how often each value has been drawn

     */



    bool        IsReversible                (_Matrix* = nil);
    // check if the matrix is reversible
    // if given a base frequencies assumes that rate matrix entries will not be multiplied by freq terms

    bool        IsAStringMatrix             (void);
    void        MakeMeSimple                (void);
    void        MakeMeGeneral               (void);
    void        ConvertToSimpleList         (_SimpleList&);
    void        CompressSparseMatrix        (bool, _Parameter*);
    //prepare the transition probs matrix for exponentiation

    long        Hash (long, long);                  // hashing function, which finds matrix
    // physical element in local storage buffer
    // returns -1 if insufficient storage
    // returns a negative number
    // if element was not found, the number returned
    // indicates the first available slot

    _Parameter*       fastIndex(void)   {
        return (!theIndex)&&(storageType==_NUMERICAL_TYPE)?(_Parameter*)theData:nil;
    }
    inline            _Parameter&         directIndex(long k)   {
        return theData[k];
    }
    long              MatrixType (void) {
        return storageType;
    }
    bool              SparseDataStructure (void) {
        return theIndex;
    }
    void              CheckIfSparseEnough (bool = false);       // check if matrix is sparse enough to justify
    void              Convert2Formulas      (void);     // converts a numeric matrix to formula-based mx
    // sparse storage

    void              Resize                (long);     // resize a dense numeric matrix to have more rows

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

    _Parameter   *theData;                            // matrix elements


protected:

    // data

    long        hDim, vDim, lDim;               // matrix physical dimensions; lDim - length of
    // actual storage allocated

    long*       theIndex;                       // indices of matrix elements in logical storage

private:

    _Parameter  computePFDR         (_Parameter, _Parameter);
    void        InitMxVar           (_SimpleList&   , _Parameter);
    bool        ProcessFormulas     (long&, _SimpleList&, _SimpleList&, _SimpleList&, _AVLListX&, bool = false, _Matrix* = nil);

    _PMathObj   PathLogLikelihood   (_PMathObj);
    /* SLKP: 20100812

     This function assumes that 'this' an 3xK matrix, where each column is of the form
     A: integer in [0,N-1], B: integer in [0,N-1], T: real >= 0

     and the argument is an NxN RATE matrix for a Markov chain

     The return value is the \sum_ j = 0 ^ {K-1} Prob {A_j -> B_j | T_j}
     */

    _PMathObj   ProfileMeanFit      (_PMathObj);

    _Matrix*    branchLengthStencil (void);

    //bool      IsAStringMatrix     (void);
    void        Add                 (_Matrix&, _Matrix&, bool sub = false);
    // aux arithmetic rountines
    bool        AddWithThreshold    (_Matrix&, _Parameter);
    void        RowAndColumnMax     (_Parameter&, _Parameter&, _Parameter* = nil);
    void        Subtract            (_Matrix&, _Matrix&);
    void        Multiply            (_Matrix&, _Parameter);
    void        Multiply            (_Matrix&, _Matrix&);
    bool        IsNonEmpty          (long);
    // checks to see if the i-th position in the storage is non-empty
    bool        CheckDimensions     (_Matrix&);
    // compare dims of 2 matrices to see if they can be multiplied
    long        HashBack            (long);
    // hashing function, which finds matrix
    // physical element given local storage
    void        MultbyS             (_Matrix&,bool,_Matrix* = nil, _Parameter* = nil);
    // internal function used in exponentiating sparse matrices

    void        Balance             (void);  // perform matrix balancing; i.e. a norm reduction which preserves the eigenvalues
    // lifted from balanc function in NR

    void        Schur               (void);  // reduce the matrix to Hessenberg form preserving eigenvalues
    // lifted from elmhes function in NR

    void        EigenDecomp         (_Matrix&,_Matrix&);  // find the eigenvalues of a real matrix
    // return real and imaginary parts
    // lifted from hqr in NR


    bool        AmISparseFast       (_Matrix&);

    bool        IncreaseStorage     (void);
    // expand the buffer, where elements are held
    // returns TRUE/FALSE for success/failure


    void        BreakPoints             (long, long, _SimpleList*);
    void        RecursiveIndexSort      (long, long, _SimpleList*);
    void        ConvertFormulas2Poly    (bool force2numbers=true);
    void        ConvertNumbers2Poly     (void);
    void        AgreeObjects            (_Matrix&);

    void        ClearFormulae           (void);             // internal reuseable purger
    void        ClearObjects            (void);             // internal reuseable purger
    inline
    _MathObject*&
    GetMatrixObject         (long ind) {
        return ((_MathObject**)theData)[ind];
    }
    inline
    bool        CheckObject             (long ind) {
        return ((_MathObject**)theData)[ind]!=nil;
    }

    void        SimplexHelper1          (long, _SimpleList&, long, bool, long&, _Parameter&);
    void        SimplexHelper2          (long&, long, _Parameter);
    void        SimplexHelper3          (long,long,long,long);
    //  helper functions for SimplexSolver

    // if nil - matrix stored conventionally

    static      int     storageIncrement,       // how many percent of full matrix size
                // to allocate to the matrix storage per increment

                precisionArg,                    // how many elements in exp series to truncate after

                switchThreshold;                 // maximum percentage of non-zero elements
    // to keep the matrix sparse

    static      _Parameter truncPrecision;

    // matrix exp truncation precision

    long        storageType,                    // true element matrix(1) or a matrix of pointers(0) which do not need to be deleted
                // 2, if elements of the matrix are actually formulas which must be initialized to numerics before use
                // matrices of type two are merely storage tables and can not be operated on directly, i.e their
                // numerical values are computed first
                bufferPerRow,                   // values reflecting internal storage structure for
                overflowBuffer,                 // sparse matrices
                allocationBlock;

    _CompiledMatrixData*
    cmd;

    _PMathObj   theValue;                       // stores evaluated values of the matrix
};

/*__________________________________________________________________________________________________________________________________________ */

class   _GrowingVector: public _Matrix
// automatically growing matrix class
{

public:
    _GrowingVector  (bool = true);
    virtual     ~_GrowingVector (void) {};

    virtual     BaseRef     makeDynamic (void); // duplicate this object into a dynamic copy
    virtual     void        Duplicate   (BaseRef); // duplicate an object from reference

    virtual     void        Clear (void);

    virtual     long        GetHDim                     (void) {
        if (isColumn) {
            return GetUsed();
        }
        return 1;
    }
    virtual     long        GetVDim                     (void) {
        if (!isColumn) {
            return GetUsed();
        }
        return 1;
    }
    long   Store            (_Parameter);
    long   GetUsed          (void) {
        return used;
    }
    void     ZeroUsed       (void) {
        used = 0;
    }

    void    operator <<     (const _SimpleList&);

    long   used;
    bool   isColumn;
};

/*__________________________________________________________________________________________________________________________________________ */

class   _NTupleStorage: public _Matrix
// a way to associate floating points
// with ordered K-tuples indices of the form 0<=i_1 < i_2 < ... < i_K < N, where N>=K
// see http://en.wikipedia.org/wiki/Combinadic
{

public:
    _NTupleStorage  (void) {
        storageN = 0;
        storageK = 0;
    }
    // null constructor; does next to nothing
    _NTupleStorage  (unsigned long,unsigned long);
    // [long] - what is the maximum value for an index (N>=0; N=0 - empty set)
    // [long] - how many indices in an ordered K-tuple (K), must be <=N; if K>N, then it will be reset to Min(1,N)
    // will allocate an array of the appropriate dimension and zero it
    // _NTupleStorage will set N=K=0 if the memory required to store the array exceeds MAX_LONG bytes (a good proxy for
    // the maximum addressable space in the system)

    virtual     ~_NTupleStorage (void) {};

    virtual     BaseRef     makeDynamic (void);
    // create a copy of the object dynamic copy

    // getters for N and K
    unsigned long GetN      (void) {
        return storageN;
    }
    unsigned long GetK      (void) {
        return storageK;
    }

    bool            CheckKTuple     (_SimpleList&);
    // modify (in-place) the argument to make it an ordered K-tuple; will also check dimensions and such
    // returns true if the k-tuple if valid; false otherwise

    unsigned long   Index           (_SimpleList&);
    // return an index into the linear array pointed to by the K-tuple in the argument (ASSUMED to be valid here!)

    _Parameter      DirectIndex     (unsigned long);
    // retrieve a value using a direct index for the K-tuple (computed by Index)

    unsigned long   Store           (_Parameter, _SimpleList&);
    // associate a value with the K-tuple; returns the direct index for the value

    _Parameter      Retrieve        (_SimpleList&);
    // return the value associated with the K-tuple

    void            IndexToTuple    (unsigned long, _SimpleList&);
    // given a direct index in 0..(NcK-1); return a K-tuple corresponding to that index


private:

    unsigned long   storageN,
             storageK;

    _SimpleList     C_NK_Lookup; // an (K+1)x(N+1) linear array which stores I choose J in element (J,I)

};



/*__________________________________________________________________________________________________________________________________________ */

class           _AssociativeList: public _MathObject
{
public:
    _AssociativeList                    (void);
    virtual ~_AssociativeList           (void) {}

    bool    ParseStringRepresentation   (_String&, bool = true, _VariableContainer* = nil);
    /* SLKP 20090803

        Parse the list represented as
            {"key": value, ...}

        the boolean argument is supplied to request reporting/suppression of error messages

        returns true on successful parse

     */

    virtual BaseRef     toStr           (void);
    virtual _PMathObj   Execute         (long opCode, _PMathObj = nil, _PMathObj = nil, _hyExecutionContext* context = _hyDefaultExecutionContext);
    virtual BaseRef     makeDynamic     (void);
    virtual _PMathObj   Compute         (void);
    virtual void        Merge           (_PMathObj);
    /* 20100907: SLKP
            A simple function to merge two lists;
            the combined list will have the key set equal to the union of the two input key sets
            if there are conflicting values for a given key, an undefined value will be stored in
            for the corresponding key


     */

    virtual void        Duplicate       (BaseRef);
    _PMathObj           MAccess         (_PMathObj);

    _PMathObj           MIterator       (_PMathObj, _PMathObj);
    /* perform a function call (ID stored in the first argument)
       having performed [an optional] conditional check on the associated key (either empty for noop or a function ID)
       Both functional IDs MUST be defined and take TWO and ONE argumens respectively

       returns the number of items processed
    */

    _PMathObj           GetByKey        (_String&, long);
    _PMathObj           GetByKey        (_String&);
    _PMathObj           GetByKey        (long, long);
    void                DeleteByKey     (_PMathObj);
    _PMathObj           MCoord          (_PMathObj);
    void                MStore          (_PMathObj, _PMathObj, bool = true, long = HY_OP_CODE_NONE);
    // SLKP 20100811: see the comment for _Matrix::MStore

    void                MStore          (_String  , _PMathObj, bool = true);
    void                MStore          (_String  , _String);
    virtual unsigned long        ObjectClass     (void)      {
        return ASSOCIATIVE_LIST;
    }
    _List*              GetKeys         (void);
    void                FillInList      (_List&);
    _String*            Serialize       (_String&);
    
    /**
     * Traverse the dictionary, cast each value into a float and return their sum.
     * Note that matrices and dictionary values will be processed recursively, i.e. "Sum" will be called on them.
     * All values that cannot be cast to a float will be treated as 0.
     * @return The sum of all dictionary elements.
     */
    _PMathObj           Sum             (void);

    _AVLListXL          avl;

private:

    _List           theData;
};

/*__________________________________________________________________________________________________________________________________________ */

extern  _Matrix *GlobalFrequenciesMatrix;
// the matrix of frequencies for the trees to be set by block likelihood evaluator
extern  _Parameter  ANALYTIC_COMPUTATION_FLAG;

void       InsertStringListIntoAVL  (_AssociativeList* , _String, _SimpleList&, _List&);
void       InsertVarIDsInList       (_AssociativeList* , _String, _SimpleList&);
#endif
