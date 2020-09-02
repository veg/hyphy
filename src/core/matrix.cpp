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
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#include "matrix.h"
#include "polynoml.h"
#include "associative_list.h"
#include "batchlan.h"
#include "likefunc.h"

#include "function_templates.h"
#include "mersenne_twister.h"
#include "global_things.h"
#include "string_file_wrapper.h"


//#include "profiler.h"

using namespace hy_global;

#define MEMORYERROR "Out of Memory"
#define ZEROOBJECT  0.0
#define ZEROPOINTER nil


_String     MATRIX_AGREEMENT            = "CONVERT_TO_POLYNOMIALS",
            ANAL_COMP_FLAG              = "ANALYTIC_COMPUTATIONS",
            ANAL_MATRIX_TOLERANCE      = "ANAL_MATRIX_TOLERANCE",
            CACHE_FORMULA_DEPENDANCY  = "CACHE_FORMULA_DEPENDANCY";


int _Matrix::precisionArg = 0;
int _Matrix::storageIncrement = 16;
int _Matrix::switchThreshold = 25;

hyFloat  _Matrix::truncPrecision = 1e-13;
#define     MatrixMemAllocate(X) MemAllocate(X, false, 64)
#define     MatrixMemFree(X)     free(X)
#define     MX_ACCESS(a,b) theData[(a)*hDim+(b)]


hyFloat     analMatrixTolerance = 1e-6,
            zero = 0,
            AUTO_PAD_DIAGONAL = 1,
            toPolyOrNot=0.0,
            toMorNot2M=1.0,
            _log2 = log (2.0);

long        ANALYTIC_COMPUTATION_FLAG = 0;

_Trie       _HY_MatrixRandomValidPDFs;



//__________________________________________________________________________________

int         fexact_                 (long , long , double *, double , double , double , double *, double *);
void        MatrixIndexError        (long, long, long, long);


// function prototypes


#ifdef _SLKP_USE_AVX_INTRINSICS
  void echo_avx_sum_4 (__m256d const x) {
    double a[4];
    _mm256_storeu_pd(a, x);
    printf ("%g|%g|%g|%g\n", a[0], a[1], a[2], a[3]);
  }
#endif




//__________________________________________________________________________________

void    MatrixIndexError (long hPos, long vPos, long hDim, long vDim) {
    HandleApplicationError (
                            kErrorStringInvalidMatrixIndex & _String ((long)hPos).Enquote('[',']') & _String ((long)vPos).Enquote('[',']') &
                            " in an " &_String (hDim) & " by " &_String (vDim) & " matrix.");
}

//_____________________________________________________________________________________________

_Matrix::_Matrix () {                             // default constructor, doesn't do much
  Initialize();
}
//_____________________________________________________________________________________________

void _Matrix::Initialize (bool) {                            // default constructor, doesn't do much
  theData         = nil;
  theIndex        = nil;
  vDim = hDim = lDim  = bufferPerRow = overflowBuffer = 0L;
  storageType     = _NUMERICAL_TYPE;
  allocationBlock = 1;
  theValue        = nil;
  compressedIndex = nil;
  
}

//_____________________________________________________________________________________________

_Matrix::_Matrix (_String const& s, bool isNumeric, _FormulaParsingContext & fpc, bool use_square_brackets) {
  // takes two separate formats
  // 1st : {{i11,...,i1n}{i21,...,i2n}....{in1,...,inn}} // all elements must be explicitly specified
  // 2st : {hor dim, <vert dim>,{hor index, vert index, value or formula}{...}...}
  
  bool                doErrors            = fpc.errMsg() == nil,
                      compute_keys_values = fpc.buildComplexObjects();
  
  _VariableContainer const* theP = fpc.formulaScope();
  
  Initialize();
  
  bool    isAConstant = true; // is this a matrix of numbers, or formulas
  char    cc;
  
  
  long    i=s.FirstNonSpaceIndex(),
          j=s.FirstNonSpaceIndex(i+1),
  k=0,
  hPos = 0,
  vPos = 0;
  
  char open_terminator  = use_square_brackets ? '[' : '{',
       close_terminator = use_square_brackets ? ']' : '}';
    
  bool    terminators [256] {false};
  terminators [(unsigned char)','] = true;
  terminators [(unsigned char)close_terminator] = true;
  
  
  try {
     auto  handle_numeric_parameter = [&] (_String& term) -> long {
       if (compute_keys_values) {
        return round (ProcessNumericArgument(&term, theP));
       }
       else {
         _String   err_msg;
         _Formula  exp (term,theP, &err_msg);
         if (exp.IsConstant()) {
           return round(exp.Compute()->Value());
         }
         return 1;
       }
     };

    if (j>i && s.length()>4) { // non-empty string
      _String term;
      if (s.char_at (i) == open_terminator && s.char_at (j) == open_terminator) { // first type
        i = j+1;
        // read the dimensions first
        
        while (i<s.length()) {
          long i2 = s.FindTerminator (i, terminators);
          if (i2 == kNotFound) {
            HandleApplicationError (kErrorStringUnterminatedMatrix & PrepareErrorContext(s, i));
          }
          i = i2;
          cc = s.char_at (i);
          
          if (cc==close_terminator) {
            break;
          }
          
          if (cc==',') {
            vDim++;
          }
          i++;
        }
        
        vDim++;
        hDim = 1;
        
        for (i = i + 1L; i<s.length()-1; i++) {
          i = s.ExtractEnclosedExpression (i,open_terminator,close_terminator,fExtractRespectQuote | fExtractRespectEscape);
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
        
        for (i=1; i<s.length()-1; i++) {
          if (s.char_at(i) == open_terminator) {
            while (s.char_at(i) != close_terminator) {
              i++;
              j = s.FindTerminator (i, terminators);
              
              if (j<0) {
                if (doErrors) {
                  HandleApplicationError (kErrorStringUnterminatedMatrix & PrepareErrorContext(s, i));
                }
                return;
              }
              
              _String lterm (s,s.FirstNonSpaceIndex(i,j-1,kStringDirectionForward),j-1); // store the term in a string
              
              //printf ("%s\n", lterm.sData);
              
              if (isNumeric) {
                if (lterm.length() == 1 && lterm.char_at(0) =='*') {
                  lterm = kEmptyString;    // dummy element in probability matrix
                }
                
                theData[vDim*hPos+vPos] = lterm.to_float ();
              } else {
                if (lterm.length() == 1 && lterm.char_at(0) =='*') {
                  lterm = kEmptyString;    // dummy element in probability matrix
                }
                
                _String errMsg;
                _Formula*  theTerm;
                
                 if (doErrors) {
                   theTerm = new _Formula (lterm, theP);
                 } else {
                   theTerm = new _Formula (lterm, theP, fpc.errMsg());
                 }
                   
                 isAConstant = isAConstant && theTerm->IsAConstant() && theTerm->ObjectClass() == NUMBER;
                ((_Formula**)theData)[vDim*hPos+vPos] = theTerm;
              }
              
              vPos++;
              if (vPos>vDim) {
                throw _String ("Rows of unequal lengths in matrix definition in ") & PrepareErrorContext(lterm, 0);
              }
              
              i=j;
            }
          }
          if (s[i]==close_terminator) {
            if (vPos!=vDim) {
              throw  kErrorStringBadMatrixDefinition & PrepareErrorContext (s,i-16);
            }
            hPos++;
            vPos = 0;
            if (hPos>hDim) {
              throw kErrorStringBadMatrixDefinition & PrepareErrorContext (s,i-16);
            }
          }
        }
        if (hPos!=hDim) {
          throw kErrorStringBadMatrixDefinition & PrepareErrorContext (s,i-16);
        }
      } else { // second type of input
        for (i=j,j=0; s.char_at (i) !=open_terminator && s.char_at (i) !=close_terminator && i<s.length(); i++) {
          if (s.char_at(i)==',') { // neither hDim nore vDim have been specified
            if (j > 0) {
              break;
            }
            term = s.Cut(1,i-1);
            hDim = handle_numeric_parameter (term);
            j    = i+1;
          }
        }
        
        if (j) { // both hDim and vDim specified
          term = s.Cut(j,i-1);
          vDim = handle_numeric_parameter (term);
        } else { // only one dim specified, matrix assumed to be square
          term = s.Cut(1,i-1);
          hDim = handle_numeric_parameter (term);
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
        
        for (; i<s.length(); i++) {
          if (s.char_at (i) ==open_terminator) {
            hPos = -1;
            vPos = -1;
            k    = i+1;
            
            for (j=i+1; j<s.length () && s.char_at (j) !=close_terminator; j++) {
              long j2 = s.FindTerminator (j, terminators);
              
              if (j2<0) {
                throw (kErrorStringUnterminatedMatrix & PrepareErrorContext(s,j));
              }
              j = j2;
              
              if (s.char_at (j) ==',') {
                term = s.Cut (s.FirstNonSpaceIndex(k,j-1,kStringDirectionForward),j-1);
                _Formula coordF (term,theP);
                hyFloat coordV = coordF.Compute()->Value();
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
              if (term.length() == 1UL && term.get_char (0)=='*') {
                term = kEmptyString;    // dummy element in probability matrix
              }
              
              (*this)[vDim*hPos+vPos];
              k = Hash (hPos,vPos);
              theData[k]=term.to_float ();
            } else {
              if (term.length() == 1UL && term.get_char (0)=='*') {
                term = kEmptyString;    // dummy element in probability matrix
              }
              
              _Formula * theTerm = new _Formula (term,theP);
              isAConstant = isAConstant && theTerm->IsAConstant();
              
              (*this)[vDim*hPos+vPos];
              k = Hash (hPos,vPos);
              ((_Formula**)theData)[k]=theTerm;
            }
            i = j;
          }
        }
      } // end else
      
      if (!isNumeric) {
        storageType = 2; // formula elements
        checkParameter (ANAL_COMP_FLAG, ANALYTIC_COMPUTATION_FLAG, 0L);
        if ((ANALYTIC_COMPUTATION_FLAG)&&!isAConstant) {
          _Matrix::ConvertFormulas2Poly (false);
        }
        
        if (isAConstant) { // a matrix of numbers - store as such
          _Matrix::Evaluate ();
        }
        _Matrix::AmISparse();
      }
    }
  } catch (const _String& err) {
    if (doErrors) {
      HandleApplicationError(err);
    }
  }
}


//_____________________________________________________________________________________________

_Matrix::_Matrix (_Matrix const& m) {
  DuplicateMatrix (this, &m);
}

//_____________________________________________________________________________________________

_Matrix::_Matrix (_SimpleList const& sl, long colArg) {
  if (sl.lLength) {
    if (colArg > 0 && colArg < sl.lLength) {
      CreateMatrix (this, sl.lLength/colArg + colArg*(sl.lLength%colArg > 0), colArg,     false, true, false);
    } else {
      CreateMatrix (this, 1, sl.lLength,  false, true, false);
    }
    for (long k=0; k<sl.lLength; k++) {
      theData[k] = sl.list_data[k];
    }
  } else {
    Initialize();
  }
}

//_____________________________________________________________________________________________

_Matrix::_Matrix (hyFloat const* inList, unsigned long rows, unsigned long columns) {
  CreateMatrix (this, rows, columns, false, true, false);
  for (unsigned long k = 0; k < rows*columns; k++) {
    theData[k] = inList[k];
  }
}

//_____________________________________________________________________________________________

_Matrix::_Matrix (hyFloat constant, unsigned long rows, unsigned long columns) {
  CreateMatrix (this, rows, columns, false, true, false);
  for (unsigned long k = 0; k < rows*columns; k++) {
    theData[k] = constant;
  }
}


//_____________________________________________________________________________________________

_Matrix::_Matrix (_List const& sl, bool parse_escapes)
// list of strings
{
  if (sl.nonempty()) {
    CreateMatrix     (this, 1, sl.lLength,  false, false, true);
    this->storageType = _FORMULA_TYPE;
    
    if (parse_escapes) {
      for (unsigned long k=0UL; k<sl.lLength; k++) {
        StoreFormula (0L,k,*new _Formula (new _FString (*(_String*) sl.GetItem(k))), false, false);
      }
    } else {
      for (unsigned long k=0UL; k<sl.lLength; k++) {
        _String* entry_k = (_String*) sl.GetItem(k);
        entry_k -> AddAReference();
        StoreFormula (0L,k,*new _Formula (new _FString (entry_k)), false, false);
      }
      
    }
  } else {
    Initialize();
  }
}

//_____________________________________________________________________________________________
void    _Matrix::CreateMatrix    (_Matrix* populate_me, long rows, long columns,  bool sparse, bool allocateStorage, bool expression_matrix) {
    
    populate_me->theValue     = nil;
    populate_me->storageType  = allocateStorage ? _NUMERICAL_TYPE : _POLYNOMIAL_TYPE;
    
    if (rows && columns) {
        if (sparse) { // store matrix as sparse
            populate_me->lDim = rows*columns*populate_me->storageIncrement/100+1; // size of storage in elements
            if (populate_me->lDim-1L<rows) {
                // either the matrix or the allocation block are too small
                // to sensibly store the matrix as sparse.
                CreateMatrix (populate_me, rows, columns, false, allocateStorage, expression_matrix);
                return;
            }
            populate_me->theIndex = (long*)MatrixMemAllocate(sizeof(long)*populate_me->lDim);
            InitializeArray(populate_me->theIndex, populate_me->lDim, -1L);
            
        } else {
            populate_me->lDim = rows*columns;
            populate_me->theIndex = nil; // no index storage needed
        }
        
        if (!allocateStorage) {
            // matrix will store pointers to elements
            populate_me->theData =(hyFloat*)MatrixMemAllocate(populate_me->lDim*sizeof(void*));
            if (expression_matrix) {
                InitializeArray ((_Formula**)populate_me->theData,    populate_me->lDim, (_Formula*)ZEROPOINTER);
            } else {
                InitializeArray ((_MathObject**)populate_me->theData, populate_me->lDim, (HBLObjectRef)ZEROPOINTER);
            }
            
        } else {
            populate_me->theData =(hyFloat*)MatrixMemAllocate (sizeof(hyFloat)*populate_me->lDim);
            memset (populate_me->theData, 0, populate_me->lDim*sizeof(hyFloat));
        }
        populate_me->hDim = rows;
        populate_me->vDim = columns;
        populate_me->SetupSparseMatrixAllocations ();
        populate_me->compressedIndex = nil;
    } else {
        populate_me->lDim      = 0L;
        populate_me->theIndex  = nil;
        populate_me->theData   = nil;
        populate_me->compressedIndex = nil;
        populate_me->hDim = 0UL;
        populate_me->vDim = 0UL;
    }
    
}


//_____________________________________________________________________________________________
void    DuplicateMatrix (_Matrix* targetMatrix, _Matrix const* sourceMatrix) {
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
    
  targetMatrix->compressedIndex = nil;
  
  if (sourceMatrix->theIndex) {
    if (!(targetMatrix->theIndex = (long*)MatrixMemAllocate(sizeof(long) *sourceMatrix->lDim))) { // allocate element index storage
      HandleApplicationError ( kErrorStringMemoryFail );
    } else {
      memcpy ((void*)targetMatrix->theIndex,(void*)sourceMatrix->theIndex,sourceMatrix->lDim*sizeof(long));
    }
    if (sourceMatrix->compressedIndex) {
        targetMatrix->compressedIndex = (long*)MatrixMemAllocate(sizeof(long) *(sourceMatrix->lDim+sourceMatrix->hDim));
        memcpy ((void*)targetMatrix->compressedIndex,(void*)sourceMatrix->compressedIndex,sizeof(long) *(sourceMatrix->lDim+sourceMatrix->hDim));
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
        if (!(targetMatrix->theData = (hyFloat*)MatrixMemAllocate(sizeof( char)*sourceMatrix->lDim*sizeof(void*)))) { // allocate element index storage
          HandleApplicationError ( kErrorStringMemoryFail );
        } else {
          memcpy ((void*)targetMatrix->theData,(void*)sourceMatrix->theData,sourceMatrix->lDim*sizeof(void*));
          if (!sourceMatrix->theIndex) { // non-sparse matrix
            for (long i=0; i<sourceMatrix->lDim; i++)
              if (sourceMatrix->GetMatrixObject(i)) {
                (sourceMatrix->GetMatrixObject(i))->AddAReference();
              }
          } else
            for (long i=0; i<sourceMatrix->lDim; i++) {
              _MathObject* theO = (sourceMatrix->GetMatrixObject(i));
              if (theO!=ZEROPOINTER) {
                theO->AddAReference();
              }
            }
          
        }
      }
    } else if (sourceMatrix->storageType==2) {
      if (targetMatrix->lDim) {
        targetMatrix->theData = (hyFloat*)MatrixMemAllocate(sourceMatrix->lDim*sizeof(void*));
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
        if (!(targetMatrix->theData =(hyFloat*)MatrixMemAllocate(sizeof( hyFloat)*targetMatrix->lDim))) { // allocate element index storage
          HandleApplicationError ( kErrorStringMemoryFail );
        } else {
          memcpy ((hyPointer)targetMatrix->theData,(hyPointer)sourceMatrix->theData,sizeof(hyFloat)*sourceMatrix->lDim);
        }
      }
    }
  } else {
    targetMatrix->theData = nil;
    targetMatrix->lDim    = 0;
  }
  
}
//_____________________________________________________________________________________________
BaseRef _Matrix::makeDynamic (void) const {
  _Matrix * result = new _Matrix;
  DuplicateMatrix (result, this);
  
  return result;
}

//_____________________________________________________________________________________________
void _Matrix::Duplicate (BaseRefConst obj) {
  Clear();
  DuplicateMatrix (this,(_Matrix const*)obj);
}


//_____________________________________________________________________________________________

_Matrix::_Matrix (long theHDim, long theVDim, bool sparse, bool allocateStorage)    // create an kEmptyString matrix of given dimensions;
                                                                                    // the flag specifies whether it is sparse or not

{
  CreateMatrix (this, theHDim, theVDim, sparse, allocateStorage);
}



//_____________________________________________________________________________________________

inline  bool    _Matrix::IsNonEmpty  (long logicalIndex) const {
    if (theIndex) {
        return theIndex[logicalIndex] != -1;
    }
    if (storageType == _NUMERICAL_TYPE) {
        return true;
    }
    return GetMatrixObject(logicalIndex)!=ZEROPOINTER;
}

//__________________________________________________________________________________

bool        _Matrix::HasChanged(bool) {
    
    switch (storageType) {
        case _POLYNOMIAL_TYPE: {
          return Any ([&] (_MathObject * f, unsigned long) -> bool {if (f) return f->HasChanged(); return false;},
                      [&] (unsigned long i) -> _MathObject * {return ((_MathObject**)theData)[i];});
        }
        break;
        
        case _FORMULA_TYPE: {
            return Any ([&] (_Formula * f, unsigned long) -> bool {if (f) return f->HasChanged(); return false;},
                        [&] (unsigned long i) -> _Formula * {return ((_Formula**)theData)[i];});
        }
        break;

        case _SIMPLE_FORMULA_TYPE: {
            if (cmd->has_volatile_entries) return true;
            return cmd->varIndex.Any ([&] (long value, unsigned long) -> bool {
                return LocateVar (value)->HasChanged();
            });
        }
        break;
    }

    return false;
}
//__________________________________________________________________________________


inline static void ROTATE(hyFloat * a, long i, long j, long k, long l, hyFloat & g, hyFloat & h, hyFloat s, hyFloat tau, long hDim) {
    // this is from NR
    g = a[i*hDim + j];
    h = a[k*hDim + l];
    a[i*hDim + j] = g - s*(h + g*tau);
    a[k*hDim + l] = h + s*(g - h*tau);
}

//__________________________________________________________________________________

bool        _Matrix::is_square_numeric(bool dense) const {
    if (storageType!=_NUMERICAL_TYPE || hDim != vDim || hDim==0L){
        HandleApplicationError ("Square numerical matrix required");
    }
    if (dense && theIndex) {
        HandleApplicationError ("Square dense numerical matrix required");
    }
    return true;
}

//__________________________________________________________________________________
void        _Matrix::Balance (void) {
    if (!is_square_numeric (true)) {
        return;
    }

    hyFloat       Squared_Radix = 2.0 * 2.0;
    bool          done = false;

    while (!done) {
        done = true;

        for (long i = 0L; i < hDim; i++) {
            hyFloat r = 0.0,
                       c = 0.0;

            for (long j = 0L; j < vDim; j++)
                if (i!=j) {
                    r += fabs (theData[i*vDim+j]);
                    c += fabs (theData[j*vDim+i]);
                }

            if (r > 0.0 && c > 0.0) {
                hyFloat g = r / Squared_Radix,
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
void        _Matrix::Schur (void) {
    if (!is_square_numeric (true)) {
        return;
    }

    for (long m = 1L; m < hDim-1; m++) {
        hyFloat x = 0.0;
        long       i = m;

        for (long j = m; j < hDim; j++)
            if (fabs (theData[j*vDim + m-1]) > x) {
                x = theData[j*vDim + m-1];
                i = j;
            }

        if (i!=m) {
            for (long j=m-1; j<hDim; j++) {
                hyFloat t = theData[i*vDim + j];
                theData[i*vDim + j] = theData[m*vDim + j];
                theData[m*vDim + j] = t;
            }
            {
                for (long j=0; j<hDim; j++) {
                    hyFloat t = theData[j*vDim + i];
                    theData[j*vDim + i] = theData[j*vDim + m];
                    theData[j*vDim + m] = t;
                }
            }
        }

        if (x)
            for (long i = m+1; i < hDim; i++) {
                hyFloat y = theData[i*vDim + m -1];
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

    for (long r = 2L; r < hDim; r++)
        InitializeArray(theData + r*hDim, r-1, 0.0);
}


//__________________________________________________________________________________
void        _Matrix::EigenDecomp (_Matrix& real, _Matrix & imag) const {
    if (!is_square_numeric()) {
        return;
    }

    hyFloat anorm = 0.0;

    for (long k = 0L; k < hDim; k++) {
        for (long k2 = k?k-1:0; k2 < hDim; k2++) {
            anorm += fabs (MX_ACCESS(k,k2));
        }
    }

    long        nn = hDim - 1L;
    hyFloat  t  = 0;

    CreateMatrix (&real, hDim, 1, false, true, false);
    CreateMatrix (&imag, hDim, 1, false, true, false);

    while (nn >= 0) {
        long its = 0,
             l   = 0;
        do {
            for (l = nn; l>=1; l--) {
                hyFloat s = fabs (MX_ACCESS(l-1,l-1)) + fabs (MX_ACCESS(l,l));
                if (s == 0.0) {
                    s = anorm;
                }

                if (fabs (MX_ACCESS(l,l-1)) + s == s) {
                    break;
                }
            }

            hyFloat x = MX_ACCESS(nn,nn);
            if (l  == nn) { // one root
                real.theData[nn]   = x + t;
                imag.theData[nn--] = 0.0;
            } else {
                hyFloat y = MX_ACCESS(nn-1,nn-1),
                           w = MX_ACCESS(nn,nn-1)*MX_ACCESS(nn-1,nn);

                if ( l == nn - 1) { // two roots
                    hyFloat p = 0.5 * (y-x),
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

                    hyFloat p,q,r,z,s;

                    if (its == 30) {
                        HandleApplicationError ("Too many QR iterations in EigenDecomp");
                        return;
                    }

                    if (its == 10 || its == 20) {
                        t += x;
                        for (long i=0; i<hDim; i++) {
                            MX_ACCESS(i,i) -= x;
                        }
                        hyFloat s = fabs(MX_ACCESS(nn,nn-1)) + fabs (MX_ACCESS(nn-1,nn-2));
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

                        hyFloat u = fabs (MX_ACCESS(m,m-1)) * (fabs (q) + fabs (r)),
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

//_____________________________________________________________________________________________
bool        _Matrix::ValidateFormulaEntries (bool callback (long, long, _Formula*)) {
    if (storageType == _FORMULA_TYPE) {
        _Formula ** formula_entires = (_Formula**)theData;
        
        long direct_index = 0L;
        for (unsigned long row = 0UL; row < hDim ; row++) {
            for (unsigned long col = 0UL; col < vDim ; col++) {
                _Formula * this_cell;
                if (is_dense()) {
                    this_cell = formula_entires[direct_index++];
                } else {
                    direct_index = Hash (row,col);
                    if (direct_index >= 0) {
                        this_cell = formula_entires[direct_index];
                    } else {
                        this_cell = nil;
                    }
                }
                if (! callback (row, col, this_cell)) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}


//__________________________________________________________________________________
HBLObjectRef   _Matrix::Eigensystem (HBLObjectRef cache) const {
    // find the eigenvectors of a symmetric matrix using Jacobi rotations
    // The original matrix is preserved.
    // returns an associative list with a sorted vector of eigenvalues and
    // a square matrix where columns are the corresponding eigenvalues
    if (!is_square_numeric()) {
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

 
                _Matrix            cpy (*this),
                *rl  = new _Matrix,
                *im  = new _Matrix;

                cpy.CheckIfSparseEnough(true);
                cpy.Balance ();
                cpy.Schur   ();
                cpy.EigenDecomp (*rl,*im);

                return & ((*new _AssociativeList) < _associative_list_key_value {"0", rl}
                       < _associative_list_key_value {"1", im});
            }
        }
    }
    _Matrix a (*this);
    a.CheckIfSparseEnough (true);

    hyFloat* b = new hyFloat[hDim],
    *   z = new hyFloat[hDim];

    _Matrix * d = new _Matrix(hDim,1,false, true),
    * v = new _Matrix(hDim,hDim,false,true);

    for (long cnt = 0, diagIndex=0; cnt < hDim; cnt ++, diagIndex+=hDim+1) {
        v->theData[diagIndex] = 1.;
        b[cnt] = (d->theData[cnt] = a.theData [diagIndex]);
        z[cnt] = 0.0;
    }

    for (int pass = 0; pass < 50; pass ++) {
        hyFloat sm = 0.,
                   tresh = 0.;

        for (long ec = 0; ec < hDim-1; ec ++)
            for (long ec2 = ec+1; ec2 < hDim; ec2++) {
                sm += fabs(a.theData[ec*hDim+ec2]);
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

                hyFloat mel = a.theData[midx],
                           g   = 100. * fabs (mel),
                           t   = fabs (d->theData[ec]),
                           c   = fabs (d->theData[ec2]);

                if (pass>3 && t+g == t && c+g == c) {
                    a.theData[midx] = 0.;
                } else if (fabs(mel) > tresh) {
                    hyFloat h = d->theData[ec2]-d->theData[ec];
                    if (fabs (h) + g == fabs (h)) {
                        t = mel/h;
                    } else {
                        hyFloat theta = 0.5*h/mel;
                        t = 1./(fabs(theta)+sqrt(1.+theta*theta));
                        if (theta<0.0) {
                            t = -t;
                        }
                    }

                    c = 1.0/sqrt(1.0+t*t);

                    hyFloat s    = t*c;
                    hyFloat tau  = s/(1.0+c);

                    h = t*mel;

                    z[ec]           -= h;
                    z[ec2]          += h;
                    d->theData[ec]  -= h;
                    d->theData[ec2] += h;


                    a.theData[midx] = 0.;

                    for (long j=0; j<ec; j++) {
                        ROTATE (a.theData, j, ec, j, ec2, g, h, s, tau, hDim);
                    }

                    for (long j=ec+1; j<ec2; j++) {
                        ROTATE (a.theData, ec, j, j, ec2, g, h, s, tau, hDim);

                    }
                    for (long j=ec2+1; j<hDim; j++) {
                        ROTATE (a.theData, ec, j, ec2, j, g, h, s, tau, hDim);
                    }
                  
                    for (long j=0; j<hDim; j++) {
                        ROTATE (v->theData, j, ec, j, ec2, g, h, s, tau, hDim);
                    }
                }
            }
        }
        for (long ec=0; ec<hDim; ec++) {
            b[ec] += z[ec];
            d->theData[ec] = b[ec];
            z[ec] = 0.;
        }
    }



    _Matrix ds (hDim,2,false, true),
    * vs = new _Matrix(hDim,hDim,false,true),
    * dss;

    for (long r=0; r<hDim; r++) {
        ds.theData[2*r]   = -d->theData[r];
        ds.theData[2*r+1] = r;
    }

    _Constant sc (0.0);
    dss = (_Matrix*)ds.SortMatrixOnColumn (&sc, nil);

    for (long r=0; r<hDim; r++) {
        d->theData[r] = -dss->theData[2*r];
        for (long c1 = r, c2 = dss->theData[2*r+1]; c1<hDim*hDim; c1+=hDim, c2+=hDim) {
            vs->theData[c1] = v->theData[c2];
        }
    }

    DeleteObject (v);
    DeleteObject (dss);

    delete [] b;
    delete [] z;

    return & ((*new _AssociativeList) << _associative_list_key_value {"0", d}
            << _associative_list_key_value {"1", vs});

}

//__________________________________________________________________________________
HBLObjectRef   _Matrix::LUDecompose (void) const {
    // perform the LU decomposition using Crout's algorithm with partial pivoting
    // The original matrix is preserved.
    // after performing this decomposition, the routine LUSolve can be called with an arbitrary vector
    // the return object is an nx(n+1) matrix which contains the LU decomposition followed
    // by a vector of row interchanges
    if (!is_square_numeric(false)) { // only works for numerical matrices at this stage
        return    new _Matrix();
    }

    hyFloat *        scalings = new hyFloat[hDim];

    long per_row = vDim+1;
    _Matrix * result = new _Matrix (hDim,per_row,false,true);
    // result is a dense matrix

    // duplicate the original matrix into result
  
  
    if (is_dense()) {//matrix is sparse
      for (long i=0L; i<hDim; i++) {
        long row_start = i*vDim;
        for (long j=0; j<vDim; j++) {
          result->theData[row_start+i+j]=theData[row_start+j];
        }
      }
    }
    else {
      for (long i=0; i<lDim; i++) {
        if (IsNonEmpty(i)) {
          long cell_coord = theIndex[i];
          result->Store(cell_coord/vDim,cell_coord%vDim,theData[i]);
        }
      }
    }

    // produce the scaling vector used in interchanging the rows
    for (long i=0L; i<vDim; i++) {
      
        hyFloat row_max = 0.0;

        for (long j=i*per_row; j<(i+1)*per_row-1; j++) {
            StoreIfGreater(row_max, fabs(result->theData[j]));
        }

        if (row_max==0.0) {
            HandleApplicationError (_String("LUDecompose doesn't work on singular matrices (row ") & i & ')');
            return    new _MathObject;
        }
        scalings[i]=1.0/row_max;
    }
    // main loop for finding L and U

    for (long j=0L; j<vDim; j++) {
        for (long i=0L; i<j; i++) {
            // fill in superdiagonal elements (U) in column j
            hyFloat sum = result->get(i,j);
            for (long k=0L; k<i; k++) {
                sum -= result->get(i,k) * result->get (k,j);
            }
            result->set (i,j) = sum;
        }
        long       max_row_index   = 0;
        hyFloat    max_row_value     = 0.0;

        for (long i=j; i<hDim; i++) {
            // calculate the unscaled version of elements of L and the diagonal
            hyFloat sum = result->get(i,j);

            for (long k=0L; k<j; k++) {
               sum -= result->get(i,k) * result->get (k,j);
            }
            result->set (i,j) = sum;

            if (StoreIfGreater(max_row_value, scalings[i]*fabs(sum))) { // find max under the diagonal in column j
                max_row_index = i;
            }
        }

        if (j!=max_row_index) { // interchange rows
            for (long k=0L; k<hDim; k++) {
                Exchange(result->set(max_row_index,k), result->set(j,k));
            }
            scalings[max_row_index]=scalings[j];
        }
        // store the index permutation
        result->theData[j*per_row+vDim] = max_row_index;

        if (result->get(j,j) == 0.0) {
            result->set(j,j) = 1.0e-25;
        }

        // divide by the pivoting element

        if (j!=hDim-1) {
            hyFloat scaler = 1.0/result->get(j,j);
            for (long i=j+1L; i<hDim; i++) {
              result->set(i,j) *= scaler;
            }
        }
    }
    delete [] scalings;
    return result;
}
//__________________________________________________________________________________
HBLObjectRef   _Matrix::LUSolve (HBLObjectRef p) const {
// takes a matrix in LU decomposed state and a vector of row permutation returned by LU
// returns a vector of solutions

    if (storageType != _NUMERICAL_TYPE || hDim+1!=vDim || vDim<=0 ) { // only works for numerical matrices at this stage
        HandleApplicationError ("LUSolve only works with numerical non-empty matrices of dimension nx(n+1) returned by LUDecompose.");
        return  nil;
    }
    if (p->ObjectClass()==MATRIX) {
      _Matrix *b=(_Matrix*)p;
      if (!((b->hDim!=hDim)||(b->vDim!=1)||(b->storageType!=1))) {
        hyFloat sum;
        _Matrix * result = new _Matrix (*b);
        result->CheckIfSparseEnough(true);
        long i = 0L,
             first_index = -1L;

        for (; i<hDim; i++) {
          long row_index = get (i, vDim - 1L);
          if (row_index<0 || row_index>=hDim) {
            break;
          }
          hyFloat sum = result->theData[row_index];
          result->theData[row_index]=result->theData[i];
          if (first_index>=0)
            for (long j=first_index; j<i; j++) {
              sum -= get (i,j) *result->theData[j];
            }
          else if (sum != 0.0) {
            first_index = i;
          }
          result->theData[i]=sum;
        }
        if (i==hDim) {
          for (i=hDim-1; i>-1; i--) {
            hyFloat sum = result->theData[i];
            for (long j=i+1L; j<hDim; j++) {
              sum -= get (i,j) *result->theData[j];
            }
            result->theData[i]=sum/get(i,i);
          }
          return result;
        }
      }
    }
    HandleApplicationError ("LUSolve expects the 2nd parameter to be a column vector defining the right hand side of LUx=b");
    return new _Matrix(1,1,false,true);
}



//__________________________________________________________________________________
HBLObjectRef   _Matrix::CholeskyDecompose (void) const
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

    if (!is_square_numeric()) { // only works for numerical square matrices at this stage
        return    new _Matrix();
    }

    long        n           = GetHDim();
    hyFloat  sum;
    _Matrix *   lower_triangular    = new _Matrix (*this);   // duplication constructor

    for (long i = 0; i < n; i++) {
        for (long j = i; j < n; j++) {
            sum = lower_triangular->get (i,j);

            for (long k = i-1L; k >= 0L; k--) {
                sum -= lower_triangular->get (i,k) * lower_triangular->get (j,k);
            }

            if (i==j) {
                if (sum <= 0.0) {   // matrix is not positive-definite
                    HandleApplicationError (_String("In CholeskyDecompose(): matrix not positive definite, (row ") & i & ')');
                    return new _MathObject;
                }

                lower_triangular->set (i, i) = sqrt(sum);
            }
            else {
                lower_triangular->set (j, i) =  sum / lower_triangular->get(i,i);
            }
        }
    }

    /* zero upper triagonal entries */
    for (long i = 0L; i < n; i++) {
        for (long j = i+1L; j < n; j++) {
            lower_triangular->set (i, j) = 0.;
        }
    }

    return lower_triangular;
}



//__________________________________________________________________________________
template <typename CALLBACK> HBLObjectRef   _Matrix::ApplyScalarOperation (CALLBACK && functor, HBLObjectRef cache) const {
    if (storageType==_NUMERICAL_TYPE) {
        _Matrix* res;
        
        if (cache && cache->ObjectClass() == MATRIX) {
            res = (_Matrix*)cache;
            *res = *this;
            res->AddAReference();
        } else {
            res = new _Matrix (*this);
        }
      
        res->ForEach ([&] (hyFloat&& value, unsigned long index, long hashed) -> void {res->theData[hashed] = functor(value);},
                      [&] (unsigned long index) -> hyFloat {return theData[index];});
      
        return res;
    }
    HandleApplicationError ("Can't apply scalar opetarations to non-numeric matrices.");
    return new _Matrix(1,1,false,true);
}

//__________________________________________________________________________________
HBLObjectRef   _Matrix::Inverse (HBLObjectRef cache) const {
    if (!is_square_numeric(false)) {
        return    new _MathObject;
    }
  
    _Matrix * LUdec = (_Matrix*)LUDecompose();
    if (LUdec) {
        _Matrix b      (hDim,1,false,true),
                * result = (_Matrix*)_returnMatrixOrUseCache(hDim,vDim,_NUMERICAL_TYPE,false, cache);
        b.theData[0]=1.0;
      for (long i=0L; i<hDim; i++) {
            if (i) {
                b.theData[i]=1.0;
                b.theData[i-1L]=0.0;
            }
            _Matrix* invVector = (_Matrix*)LUdec->LUSolve(&b);
            _Matrix corrTerm (GetHDim(),1, false, true);
             Multiply(corrTerm, *invVector);
            corrTerm -= b;
            //_Matrix* corrTerm = (_Matrix*)(*this*(*invVector)-b).makeDynamic();
            _Matrix* corrX =  (_Matrix*)LUdec->LUSolve(&corrTerm);
            *invVector-=*corrX;
            DeleteObject (corrX);
            for (long j=0; j<hDim; j++) {
                result->set (j,i) = invVector->theData[j];
            }
            DeleteObject (invVector);
        }
        DeleteObject (LUdec);
        return result;
    }
    return new _Matrix (1,1,false,true);

}

//__________________________________________________________________________________
HBLObjectRef   _Matrix::MultByFreqs (long freqID, bool reuse_value_object) {
// multiply this transition probs matrix by frequencies
    HBLObjectRef value = ComputeNumeric(true);//!reuse_value_object);
    
    //printf ("\n%s\n", _String ((_String*)toStr()).get_str());

    if (freqID>=0) {
        _Matrix* freq_matrix = nil;
        freqID = modelFrequenciesIndices.list_data[freqID];
        if (freqID>=0) {
            freq_matrix = (_Matrix*)LocateVar(freqID)->GetValue();
            if (freq_matrix->storageType != _NUMERICAL_TYPE) {
                if (freq_matrix->theValue) {
                    freq_matrix = (_Matrix*)freq_matrix->theValue;
                } else {
                    freq_matrix = (_Matrix*)freq_matrix->ComputeNumeric();
                }
            }
        }

        if (theIndex) {
            _Matrix*    vm = (_Matrix*) value;
            hyFloat *dp = vm ->theData;
            hyFloat *tempDiags = (hyFloat*) alloca (sizeof(hyFloat) * hDim);
            InitializeArray(tempDiags, hDim, 0.0);

            if (freq_matrix) {
                  for (long i=0; i<lDim; i++) {
                      long p = theIndex[i];
                      if (p != -1) {
                          long h = p/vDim;
                          p -= h*vDim;
                          if (h!=p) {
                              tempDiags[h] += (dp[i] *= freq_matrix->theData[p]);
                          }
                      }
                  }
            }
            else {
                  for (long i=0; i<lDim; i++) {
                      long p = theIndex[i];
                      if (p != -1) {
                          long h = p/vDim;
                          p -= h*vDim;
                          if (h!=p) {
                              tempDiags[h] += dp[i];
                          }
                      }
                  }
            }
          
            for (long j=0L; j<hDim; j++) {
                vm->Store (j,j,-tempDiags[j]);
            }

        } else {
            hyFloat * theMatrix = ((_Matrix*)value)->theData;

            if (freq_matrix) {
                if (freq_matrix->theIndex) {
                    for (long i=0; i<lDim; i++) {
                        theMatrix[i] *= (*freq_matrix)[i%vDim];
                    }
                } else {
                    for (unsigned long column=0UL; column<vDim; column++) {
                      const hyFloat freq_i = freq_matrix->theData[column];
                      unsigned long entry = column;
                      for (;entry < lDim - vDim; entry += vDim) {
                        theMatrix[entry] *= freq_i;
                        theMatrix[entry+=vDim] *= freq_i;
                      }
                      if (entry < lDim) {
                        theMatrix[entry] *= freq_i;
                      }
                    }
                }
            }
          
            for (unsigned long row_start = 0UL, row = 0UL; row_start < lDim; row_start+=vDim, row++) {
              unsigned long diag = row_start + row;
              theMatrix [diag] = 0.;
              for (unsigned long col = 0UL; col < row; col++) {
                theMatrix[diag] -= theMatrix[row_start + col];
              }
              for (unsigned long col = row+1; col < vDim; col++) {
                theMatrix[diag] -= theMatrix[row_start + col];
              }
            }
        }

    }
    return value;
}


//__________________________________________________________________________________
HBLObjectRef   _Matrix::Compute (void) {
  //if ((storageType != 1)&&(storageType != 2))
  if (storageType != _NUMERICAL_TYPE) {
    if (storageType == _POLYNOMIAL_TYPE) {
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
    
    if (storageType != _SIMPLE_FORMULA_TYPE) {
      theValue  = Evaluate(false);
    } else {
      theValue  = EvaluateSimple ();
    }
    return theValue;
  }
  return this;
}

//__________________________________________________________________________________
HBLObjectRef   _Matrix::ComputeNumeric (bool copy) {
    if (storageType != _NUMERICAL_TYPE) {
        if (storageType == 0 && ANALYTIC_COMPUTATION_FLAG) {
            return this;
        }
        if (storageType != _SIMPLE_FORMULA_TYPE) {
            if (theValue) {
                DeleteObject (theValue);
            }
            theValue  = Evaluate(false);
        } else {
            if (copy) {
                if (theValue) {
                    DeleteObject (theValue);
                }
                theValue = EvaluateSimple();
            } else {
                theValue = EvaluateSimple((_Matrix*)theValue);
            }
        }
        return theValue;
    }
    if (copy) {
        if (theValue) {
            DeleteObject (theValue);
        }
        return (theValue = (_Matrix*)makeDynamic());
    }
    return this;
}

//__________________________________________________________________________________
HBLObjectRef   _Matrix::RetrieveNumeric (void) {
    if (storageType != _NUMERICAL_TYPE) {
        if (theValue) {
            return theValue;
        }

        return ComputeNumeric();
    }
    return this;
}

//__________________________________________________________________________________
HBLObjectRef   _Matrix::Sum (HBLObjectRef cache) {
    return _returnConstantOrUseCache(MaxElement (1), cache);
}

//__________________________________________________________________________________


HBLObjectRef _Matrix::ExecuteSingleOp (long opCode, _List* arguments, _hyExecutionContext* context, HBLObjectRef cache)  {

  
    switch (opCode) { // first check operations without arguments
      case HY_OP_CODE_ABS: // Abs
        return Abs(cache);
      case HY_OP_CODE_COLUMNS:  //Columns
        return _returnConstantOrUseCache(vDim, cache);
      case HY_OP_CODE_INVERSE: //Inverse
        return Inverse(cache);
      case HY_OP_CODE_EIGENSYSTEM: //Eigensystem
        return Eigensystem(cache);
      case HY_OP_CODE_EVAL: //Eval
        return (HBLObjectRef)ComputeNumeric()->makeDynamic();
      case HY_OP_CODE_EXP: //Exp
        return Exponentiate();
      case HY_OP_CODE_LUDECOMPOSE: // LUDecompose
        return LUDecompose();
      case HY_OP_CODE_LOG: // Log
        return ApplyScalarOperation ([] (hyFloat h) -> hyFloat {return log (h);}, cache);
      case HY_OP_CODE_ROWS: // Rows
        return _returnConstantOrUseCache (hDim, cache);
      case HY_OP_CODE_SIMPLEX: // Simplex
        return SimplexSolve();
      case HY_OP_CODE_TRANSPOSE: { // Transpose
        _Matrix* result = new _Matrix (*this);
        result->Transpose();
        return result;
      }
      case HY_OP_CODE_TYPE: // Type
        return Type(cache);
   }
  
  _MathObject * arg0 = _extract_argument (arguments, 0UL, false);
  
  switch (opCode) { // next check operations without arguments or with one argument
    case HY_OP_CODE_ADD: // +
      if (arg0) {
        return AddObj (arg0, cache);
      } else {
        return Sum (cache);
      }
      break;
    case HY_OP_CODE_SUB: // -
      if (arg0) {
        return SubObj(arg0, cache);
      } else {
        return ApplyScalarOperation ([] (hyFloat h) -> hyFloat {return -h;}, cache);

        //return (HBLObjectRef)((*this)*(-1.0)).makeDynamic();
      }
      break;
  }
  
  if (arg0) {
    switch (opCode) { // operations that require exactly one argument
      case HY_OP_CODE_IDIV: // $
      case HY_OP_CODE_DIV:  // /
        return MultElements(arg0,opCode == HY_OP_CODE_DIV, cache);
      case HY_OP_CODE_MOD: // %
        return SortMatrixOnColumn (arg0, cache);
      case HY_OP_CODE_AND: // &&
        return pFDR (arg0, cache);
      case HY_OP_CODE_MUL: // *
        return MultObj(arg0, cache);
      case HY_OP_CODE_LESS: // <
        return PathLogLikelihood(arg0, cache);
      case HY_OP_CODE_LEQ: // <=
        return K_Means(arg0, cache);
      case HY_OP_CODE_EQ: // ==
        return _returnConstantOrUseCache(Equal (arg0), cache);
        //return ProfileMeanFit(arg0);
      case HY_OP_CODE_GREATER: // >
        return NeighborJoin (!CheckEqual(arg0->Value(),0.0), cache);
      case HY_OP_CODE_GEQ: // >=
        return MakeTreeFromParent (arg0->Value(), cache);
      case HY_OP_CODE_CCHI2: //CChi2
        if (arg0->ObjectClass()==NUMBER && arg0->Value()>0.999 ) {
          return _returnConstantOrUseCache (FisherExact(5.,80.,1.), cache);
        } else {
          return _returnConstantOrUseCache (FisherExact(0.,0.,0.), cache);
        }
      case HY_OP_CODE_LUSOLVE: // LUSolve
        return LUSolve (arg0);
      case HY_OP_CODE_RANDOM: // Random
        return Random (arg0, cache);
      case HY_OP_CODE_POWER: // ^ (Poisson log-likelihood)
          return  PoissonLL (arg0, cache);
      case HY_OP_CODE_MAX: // Max
      case HY_OP_CODE_MIN: // Max
        if (arg0->ObjectClass()==NUMBER) {
          if (CheckEqual (arg0->Value(), 1)) {
            long index = 0L;
            hyFloat v[2] = {opCode == HY_OP_CODE_MAX?MaxElement (0,&index):MinElement(0,&index),0.0};
            v[1] = index;
            return new _Matrix (v,1,2);
          }
        }
        return _returnConstantOrUseCache (opCode == HY_OP_CODE_MAX?MaxElement (0):MinElement (0), cache);
   }
    _MathObject * arg1 = _extract_argument (arguments, 1UL, false);
    
     switch (opCode) {
        
      case HY_OP_CODE_MACCESS: // MAccess
        return MAccess (arg0,arg1, cache);
        
      case HY_OP_CODE_MCOORD: // MCoord
        return MCoord (arg0, arg1, cache);
    }
    
  }
  
  switch (opCode) {
    case HY_OP_CODE_ADD: // +
    case HY_OP_CODE_SUB: // -
    case HY_OP_CODE_IDIV: // $
    case HY_OP_CODE_DIV:  // /
    case HY_OP_CODE_MOD: // %
    case HY_OP_CODE_AND: // &&
    case HY_OP_CODE_MUL: // *
    case HY_OP_CODE_LESS: // <
    case HY_OP_CODE_LEQ: // <=
    case HY_OP_CODE_EQ: // ==
    case HY_OP_CODE_GREATER: // >
    case HY_OP_CODE_GEQ: // >=
    case HY_OP_CODE_CCHI2: //CChi2
    case HY_OP_CODE_LUSOLVE: // LUSolve
    case HY_OP_CODE_RANDOM: // Random
    case HY_OP_CODE_POWER: // ^ (Poisson log-likelihood)
    case HY_OP_CODE_MAX: // Max
    case HY_OP_CODE_MIN: // Max
    case HY_OP_CODE_MACCESS: // MAccess
    case HY_OP_CODE_MCOORD: // MCoord
      WarnWrongNumberOfArguments (this, opCode,context, arguments);
      break;
    default:
      WarnNotDefined (this, opCode,context);
  }

   return new _MathObject;
}

//_____________________________________________________________________________________________
bool    _Matrix::AmISparse(void)
{
    if (theIndex) {
        return true;    // duh!
    }
  
    if (storageType== _FORMULA_TYPE || _SIMPLE_FORMULA_TYPE) {
        return false;
    }

    long k=0L;
    if (storageType==_NUMERICAL_TYPE) {
      for (long i=0; i<lDim; i++) {
          if (theData[i]!=ZEROOBJECT) {
              k++;
          }
      }
    } else {
      for (long i=0; i<lDim; i++) {
          if (IsNonEmpty(i) && !GetMatrixObject(i)->IsObjectEmpty()) {
              k++;
          }
      }
    }


    if ((hyFloat(k)/lDim*100.)<=_Matrix::switchThreshold) {
        // we indeed are sparse enough
        _Matrix sparseMe (hDim,vDim,true,storageType==_NUMERICAL_TYPE);
        if (storageType==_NUMERICAL_TYPE) {
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
                GetMatrixObject(i)->AddAReference();
            }
        }

        Clear();
        DuplicateMatrix (this, &sparseMe);
        return true;
    }
    return false;
}

//_____________________________________________________________________________________________
bool    _Matrix::AmISparseFast (_Matrix& whereTo) {
    if (theIndex) {
        return true;    // duh!
    }

    long k = 0L,
         threshold = lDim*_Matrix::switchThreshold/100;
    
    for (unsigned long i=0UL; i<lDim ; i++) {
          if (__builtin_expect(theData[i] != ZEROOBJECT, 1L)) {
              if (++k >= threshold) {
                  break;
              }
          }
    }

    if (k < threshold) {
        // we indeed are sparse enough
        
        if (k == 0L) {
            k = 1L;
        }

        hyFloat *          newData  = (hyFloat*)MatrixMemAllocate (k*sizeof(hyFloat));
        
        if (whereTo.theIndex) {
            free (whereTo.theIndex);
        }
        whereTo.theIndex               = (long*)MemAllocate (k*sizeof(long));

        long p = 0;
        whereTo.theIndex[0] = -1;

        for (unsigned long i=0; i < lDim; i++)
            if (theData[i]!=ZEROOBJECT) {
                whereTo.theIndex[p] = i;
                newData[p++] = theData[i];
            }

        whereTo.lDim     = k;
        free     (whereTo.theData);
        whereTo.theData = newData;
        return true;
    }

    return false;
}

//_____________________________________________________________________________________________

bool    _Matrix::IsValidTransitionMatrix() const {
    if (is_square() && is_numeric()) {
        long d = GetHDim();
        hyFloat * sums = new hyFloat [d] {0.0};
        long idx = 0L;
        for (long r = 0L; r < d; r++) {
            for (long c = 0L; c < d; c++, idx++) {
                hyFloat term = theData[idx];
                if (term < 0.0 || term > 1.0) {
                    delete [] sums;
                    return false;
                }
                sums[r] += term;
            }
        }
        for (long r = 0L; r < d; r++) {
            if (!CheckEqual(1.0, sums[r])) {
                delete [] sums;
                return false;
            }
        }
        delete [] sums;
        return true;
    }
    return false;
}

    
//_____________________________________________________________________________________________

bool    _Matrix::IsReversible(_Matrix* freqs) {
    
    try {
        
        if (!is_square()) {
            throw _String ("Not a square matrix in _Matrix::IsReversible");
        }
        
        if (freqs && freqs->GetHDim () * freqs->GetVDim () != GetHDim()) {
            throw _String ("Incompatible frequency and rate matrix dimensions in _Matrix::IsReversible");
        }
        
        if (!is_numeric() && !is_expression_based()) {
            throw _String ("Unsupported rate matrix type in _Matrix::IsReversible");
        }

        if (freqs && !freqs->is_numeric() && !freqs->is_expression_based()) {
            throw _String ("Unsupported frequency matrix type in _Matrix::IsReversible");
        }


        bool   needAnalytics = is_expression_based() || (freqs && freqs->is_expression_based());
        if (needAnalytics) {
            if (freqs) {
                for (long r = 0; r < hDim; r++)
                    for (long c = r+1; c < hDim; c++) {
                        bool compResult = true;
                        if (is_expression_based()) {
                            _Formula* rc = GetFormula(r,c),
                                      * cr = GetFormula(c,r);

                            if (rc && cr) {
                                _Polynomial *rcp = (_Polynomial *)rc->ConstructPolynomial(),
                                             *crp = (_Polynomial *)cr->ConstructPolynomial();

                                if (rcp && crp) {
                                    HBLObjectRef     tr = nil,
                                                     tc = nil;

                                    if (freqs->is_expression_based()) {
                                        if (freqs->GetFormula(r,0)) {
                                            tr = freqs->GetFormula(r,0)->ConstructPolynomial();
                                            if (tr) {
                                                tr->AddAReference();
                                            } else {
                                                throw _String ("Could not convert matrix cell (") & r & ',' & c & ") to a polynomial";
                                            }
                                        }
                                        if (freqs->GetFormula(c,0)) {
                                            tc = freqs->GetFormula(c,0)->ConstructPolynomial();
                                            if (tc) {
                                                tc->AddAReference();
                                            } else {
                                                DeleteObject (tr);
                                                throw _String ("Could not convert frequency cell (") & c & ") to a polynomial";
                                             }
                                        }
                                    } else {
                                        tr = new _Constant ((*freqs)[r]);
                                        tc = new _Constant ((*freqs)[c]);
                                    }
                                    if (tr && tc) {
                                        _Polynomial        * rcpF = (_Polynomial*)rcp->Mult(tr, nil),
                                                           * crpF = (_Polynomial*)crp->Mult(tc, nil);

                                        compResult         = rcpF->Equal(crpF);
                                        DeleteObject (rcpF);
                                        DeleteObject (crpF);
                                        
                                    } else {
                                        compResult = !(tr||tc);
                                    }
                                    
                                    DeleteObject (tr);
                                    DeleteObject (tc);
                                    
                                    if (!compResult) {
                                         throw _String ("Unequal cells at (") & r & ',' & c & ")";
                                    }
                                } else {
                                    throw _String ("Could not convert a matrix cell at (") & r & ',' & c & ") to a polynomial: " & _StringBuffer ((_StringBuffer*)rc->toStr(kFormulaStringConversionNormal));
                                }
                             } else {
                                compResult = !(rc || cr);
                            }
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
    } catch (const _String& reason) {
        ReportWarning (_String ("Reversibility checks failed: ") & reason);
        return false;
    }
    return false;
}

//_____________________________________________________________________________________________

void    _Matrix::CheckIfSparseEnough(bool force) {

// check if matrix is sparse enough to justify compressed storage

    if (theIndex && (force || lDim>hDim*vDim*::_Matrix::switchThreshold/100)) {
        // switch to normal matrix storage - more than half elements are non-zero
        // -= allocationBlock;

        long square_dimension = vDim*hDim;

        if (!is_numeric()) {
            // pointers
            hyPointer* tempData = (hyPointer*) MemAllocate (square_dimension*sizeof(hyPointer));
            InitializeArray(tempData, square_dimension, (hyPointer)nil);
            
            for (unsigned long i = 0UL; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    tempData[theIndex[i]]=((hyPointer*)theData)[i];
                }
            }
            MatrixMemFree( theData);
            theData = (hyFloat*)tempData;
       } else {
            //objects
            hyFloat* tempData = (hyFloat*) MemAllocate (square_dimension*sizeof(hyFloat));
            InitializeArray(tempData, square_dimension, 0.0);

            for (unsigned long i = 0UL; i<lDim; i++) {
                long k = theIndex[i];
                if (k >= 0) {
                    tempData [k] = ((hyFloat*)theData) [i];
                }
            }
            MatrixMemFree( theData);
            theData = (hyFloat*)tempData;

        }
        lDim = square_dimension;
        MatrixMemFree( theIndex);
        theIndex = nil;
    }
}

//_____________________________________________________________________________________________
bool    _Matrix::IncreaseStorage    (void)
{
    lDim += allocationBlock;

    long* tempIndex, i;

    if (!(tempIndex = (long*)MatrixMemAllocate(lDim*sizeof(long)))) {
        HandleApplicationError ( kErrorStringMemoryFail );
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
            HandleApplicationError ( kErrorStringMemoryFail );
        } else {
            memcpy (tempData, theData, (lDim-allocationBlock)*sizeof(void*));
            MatrixMemFree (theData);
            for (i = lDim-1; i>=lDim-allocationBlock; i--) {
                tempData [i] = ZEROPOINTER;
            }
            theData = (hyFloat*)tempData;
        }
    } else
        //objects
    {
        hyFloat* tempData;
        if (!(tempData =  (hyFloat*)MatrixMemAllocate(sizeof(hyFloat)* lDim))) {
            HandleApplicationError ( kErrorStringMemoryFail );
        } else {
            for (i = lDim-1; i>=lDim-allocationBlock; i--) {
                tempData [i] = ZEROOBJECT;
            }
            for (; i>=0; i--) {
                tempData [i] = ((hyFloat*)theData) [i];
            }
            MatrixMemFree( theData);
            theData = tempData;
        }
    }
    return TRUE;

}


//_____________________________________________________________________________________________

void    _Matrix::Convert2Formulas (void)
{
    if (storageType == 1) {
        storageType = 2;
        _Formula** tempData = (_Formula**)MatrixMemAllocate (sizeof(void*)*lDim);
        if (!theIndex) {
            for (long i = 0; i<lDim; i++) {
                tempData[i] = new _Formula (new _Constant (((hyFloat*)theData)[i]));
            }
        } else
            for (long i = 0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    //_Constant c (((hyFloat*)theData)[i]);
                    //_Formula f((_PMathObj)c.makeDynamic());
                    //tempData[i] = (_Formula*)f.makeDynamic();
                    tempData[i] = new _Formula (new _Constant (((hyFloat*)theData)[i]));
                } else {
                    tempData[i]=nil;
                }
            }

        MatrixMemFree (theData);
        theData = (hyFloat*)tempData;
    }
}




//_____________________________________________________________________________________________

void    _Matrix:: ScanForVariables(_AVLList& theReceptacle, bool inclG, _AVLListX* tagger, long weights) const {
    ScanForVariables2 (theReceptacle, inclG, -1, true, tagger, weights);
}
//_____________________________________________________________________________________________

void    _Matrix:: ScanForVariables2(_AVLList& theReceptacle, bool inclG, long modelID, bool inclCat, _AVLListX* tagger, long weights) const {
    if (is_expression_based()) { // a formula based matrix, there is stuff to do
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

                    for (unsigned long k=0; k<sl1.lLength; k++) {
                        cachedValues->theData[k] = sl1.list_data[k];
                    }
                    {
                        for (unsigned long k=sl1.lLength; k<sl2.lLength; k++) {
                            cachedValues->theData[k] = -1.;
                        }
                    }
                    {
                        for (unsigned long k=0; k<sl2.lLength; k++) {
                            cachedValues->theData[k+sl2.lLength] = sl2.list_data[k];
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
                        if (tagger) {
                            tagger->UpdateValue((BaseRef)vI, weights, 0);
                        }
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
                    theFormulas[i]->ScanFForVariables(theReceptacle,inclG,false,inclCat, false, tagger, weights);
                }
        } else
            for (long i = 0; i<lDim; i++) {
                if (theFormulas[i]!=(_Formula*)ZEROPOINTER) {
                    theFormulas[i]->ScanFForVariables (theReceptacle,inclG,false,inclCat, false, tagger, weights);
                }
            }
    } else if (storageType == 0) { // a polynomial based matrix, there is stuff to do
        _MathObject ** thePoly = (_MathObject**)theData;
        if (theIndex)
            for (long i = 0; i<lDim; i++) {
                if (IsNonEmpty(i)) {
                    thePoly[i]->ScanForVariables(theReceptacle,inclG,tagger, weights);
                }
            }
        else
            for (long i = 0; i<lDim; i++) {
                if (thePoly[i]!=ZEROPOINTER) {
                    thePoly[i]->ScanForVariables (theReceptacle,inclG,tagger, weights);
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

bool        _Matrix::ProcessFormulas (long& stackLength, _AVLList& varList,   _SimpleList& newFormulas,
                                      _SimpleList& references, _AVLListX& flaStrings,
                                      bool runAll, _Matrix * stencil) {
    _Formula *      thisFormula = nil;
    _Formula **     theFormulas = (_Formula**)theData;

    bool isGood = true;

    if (theIndex) {
        for (long i = 0L; i<lDim; i++) {
            long cellIndex = theIndex [i];
            if (cellIndex>-1) {
                if (stencil && CheckEqual(stencil->theData[cellIndex],0.0)) {
                    references << -1;
                    continue;
                }
                thisFormula = theFormulas[i];

                if (runAll || thisFormula->AmISimple(stackLength,varList)) {
                    _String * flaString = (_String*)thisFormula->toStr(kFormulaStringConversionNormal, nil,true);
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
        for (long i = 0L; i<lDim; i++) {
            if ((theFormulas[i]!=(_Formula*)ZEROPOINTER)&&(!theFormulas[i]->IsEmpty())) {
                thisFormula = theFormulas[i];

                if (stencil && CheckEqual(stencil->theData[i],0.0)) {
                    references << -1;
                    continue;
                }

                if (runAll || thisFormula->AmISimple(stackLength,varList)) {
                    _String * flaString = (_String*)thisFormula->toStr(kFormulaStringConversionNormal, nil,true);
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
_Matrix*        _Matrix::BranchLengthStencil (void) const {
    
    _Matrix * stencil = (_Matrix*)hy_env::EnvVariableGet(hy_env::branch_length_stencil, MATRIX);
    if (stencil) {
        if (stencil->storageType == _NUMERICAL_TYPE && stencil->hDim==stencil->vDim && stencil->hDim == hDim) {
            stencil->CheckIfSparseEnough (true);
        } else {
            stencil = nil;
        }
    }

    return stencil;
}

//_____________________________________________________________________________________________
_String*        _Matrix::BranchLengthExpression (_Matrix* baseFreqs, bool mbf) {
    if (storageType == _FORMULA_TYPE) {

        long            stack_length = 0L;

        _SimpleList     new_formulas,
                        references;

        _List           converted_expressions;
        _AVLListX       converted_expressions_avl(&converted_expressions);
        _Matrix*        stencil = BranchLengthStencil();

        print_digit_specification = hy_env::EnvVariableGetDefaultNumber(hy_env::print_float_digits);
      
       _SimpleList varList = PopulateAndSort([&] (_AVLList & list) -> void {
           ProcessFormulas (stack_length,list,new_formulas,references,converted_expressions_avl,true,stencil);
       });

        _StringBuffer * sendMeBack = new _StringBuffer(256L);
      
        if (baseFreqs->is_numeric()) {
            // numerical base frequencies
            _Matrix   multipliersByRate (new_formulas.countitems(),1,false,true);
            
            ForEach([this, &multipliersByRate, &references, baseFreqs, mbf] (BaseRef object, unsigned long direct_index, unsigned long index) -> void {
                long this_ref = references.get (index);
                if (this_ref >= 0) {
                    multipliersByRate.set(0,this_ref) += (*baseFreqs)(direct_index/vDim,0) *
                                                         (mbf?(*baseFreqs)(direct_index%vDim,0):1.0);
                }
                
            }, [] (unsigned long) -> BaseRef {return nil;});
            

            for (unsigned long k=0UL; k<new_formulas.countitems(); k++) {
                hyFloat this_multiplier = multipliersByRate (k, 0);
                
                if (!CheckEqual(this_multiplier,0.0)) {
                    if (sendMeBack->nonempty()) {
                        (*sendMeBack) << '+';
                    }
                    
                    (*sendMeBack) << '('
                                  << (_String*)converted_expressions.GetItem(k)
                                  << ")*"
                                  << _String(this_multiplier, print_digit_specification);
                  
                }
            }
        } else if (baseFreqs->is_expression_based()) {
            // formula-based equilibrium frequencies
            
            _List   freqFla,
                    multipliersByRate;

            for (long k=0L; k<new_formulas.countitems(); k++) {
                multipliersByRate.AppendNewInstance(new _StringBuffer (128L));
            }

            for (long k=0L; k<hDim; k++) {
                freqFla.AppendNewInstance ((_String*)baseFreqs->GetFormula(k,0)->toStr(kFormulaStringConversionNormal, nil,true));
            }

            ForEach([this, &multipliersByRate, &references, baseFreqs, mbf, &freqFla] (BaseRef object, unsigned long direct_index, unsigned long index) -> void {
                long this_ref = references.get (index);
                if (this_ref >= 0L) {
                    _StringBuffer * thisAdder = (_StringBuffer*)multipliersByRate(this_ref);
                    if (thisAdder->nonempty()) {
                        (*thisAdder) << '+';
                    }
                    (*thisAdder) << '(';
                    if (mbf) {
                        (*thisAdder) << (_String*)freqFla(direct_index%vDim)
                        << ")*(";
                    }
                    (*thisAdder) << (_String*)freqFla(direct_index/vDim) << ')';
                }
            }, [] (unsigned long) -> BaseRef {return nil;});

            for (long k=0L; k<new_formulas.countitems(); k++) {
                ((_StringBuffer*)multipliersByRate(k))->TrimSpace();
                if (k) {
                    (*sendMeBack) << '+';
                }
                
                (*sendMeBack) << '('
                    << (_String*)converted_expressions.GetItem(k)
                    << ")*("
                    << (_String*)multipliersByRate(k)
                    << ')';
            }
        }
        sendMeBack->TrimSpace();
        if (sendMeBack->nonempty()) {
            _Formula        blF (*sendMeBack);
            _Polynomial*    isPoly = (_Polynomial*)blF.ConstructPolynomial();
            if (isPoly) {
                DeleteObject (sendMeBack);
                sendMeBack = (_StringBuffer*)isPoly->toStr();
            }
        }
        return sendMeBack;
    }
    return new _String;
}

//_____________________________________________________________________________________________
void        _Matrix::MakeMeSimple (void) {
    if (storageType == _FORMULA_TYPE) {
        long            stackLength = 0L;

        _SimpleList     newFormulas,
                        references;

        _List           flaStringsL;
        _AVLListX       flaStrings(&flaStringsL);

 
        _SimpleList varListAux;
        _AVLList    varList (&varListAux);
      
        if (ProcessFormulas (stackLength,varList,newFormulas,references,flaStrings)) {
            storageType = _SIMPLE_FORMULA_TYPE;

            cmd                         = new _CompiledMatrixData;
            cmd->has_volatile_entries   = false;
          
            for (unsigned long k = 0; k < newFormulas.lLength; k++) {
                cmd->has_volatile_entries = ((_Formula*)newFormulas.get(k))->ConvertToSimple(varList) || cmd->has_volatile_entries;
            }

            cmd->varIndex.Duplicate     (&varListAux);
            cmd->theStack               = (_SimpleFormulaDatum*)MatrixMemAllocate (stackLength*sizeof(_SimpleFormulaDatum));
            cmd->varValues              = (_SimpleFormulaDatum*)MatrixMemAllocate ((cmd->varIndex.countitems()>0?varList.countitems():1)*sizeof(_SimpleFormulaDatum));
            long allocation_size = MAX (references.lLength, 1) * sizeof (long);
            cmd->formulaRefs            = (long*)MemAllocate (allocation_size);
            memcpy (cmd->formulaRefs, references.list_data, allocation_size);
            cmd->formulaValues          = new hyFloat [newFormulas.lLength];
            cmd->formulasToEval.Duplicate (&newFormulas);
        }

    }
}
//_____________________________________________________________________________________________
void        _Matrix::MakeMeGeneral (void) {
    if (storageType == _SIMPLE_FORMULA_TYPE) {
        for (long k = 0L; k < cmd->formulasToEval.lLength; k++) {
            ((_Formula*)cmd->formulasToEval.list_data[k])->ConvertFromSimpleList(cmd->varIndex);
        }

        delete [] cmd->formulaValues;
        free   (cmd->formulaRefs);

        MatrixMemFree   (cmd->theStack);
        MatrixMemFree   (cmd->varValues);
        delete          (cmd);
        cmd             = nil;
        storageType     = _FORMULA_TYPE;
    }
}
//_____________________________________________________________________________________________
HBLObjectRef   _Matrix::Evaluate (bool replace)
// evaluate the matrix  overwriting (or not) the old one
{
    _Matrix result (hDim, vDim, bool (theIndex), true);

    if (storageType == 2) {
        HBLObjectRef formValue = nil;
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
                        hyFloat *st = &result[k];
                        *st=0;
                        for (long j = 0; j<vDim; j++) {
                            if (j==i) {
                                continue;
                            }
                            *st-=result(i,j);
                        }
                    } else if (k<0) {
                        hyFloat *st = &result[i*vDim+i];
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
                            hyFloat st = 0;
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
        HBLObjectRef polValue = nil;
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
        return (HBLObjectRef)result.makeDynamic();
    }
    return nil;
}

//_____________________________________________________________________________________________
void        _Matrix::ConvertToSimpleList (_SimpleList & sl)
{
    sl.Clear();
    if (storageType == _NUMERICAL_TYPE) {
        sl.RequestSpace (hDim*vDim+1);

        for (long i=0; i<hDim; i++)
            for (long j=0; j<vDim; j++) {
                sl << (*this)(i,j);
            }
    } else {
      if (storageType == _FORMULA_TYPE) {
        _Matrix * c = ((_Matrix*)Compute());
        if (c->storageType == _NUMERICAL_TYPE) {
          c -> ConvertToSimpleList (sl);
        }
      }
    }
}

//_____________________________________________________________________________________________
bool        _Matrix::IsAStringMatrix (void) const
// check if a formula matrix contains strings
{
    if (is_expression_based()) {
        try {
            return Any ([&] (_Formula * f, unsigned long) -> bool {
                            if (f) {
                                if (f->ObjectClass() == STRING)
                                    return true;
                                throw (0);
                            }
                            return false;
                        },
                        [&] (unsigned long i) -> _Formula * {return ((_Formula**)theData)[i];});
            
        } catch (int ) {
            return false;
        }

    }
    return false;
}

//_____________________________________________________________________________________________
void        _Matrix::FillInList (_List& fillMe, bool convert_numbers) const {
// check if a formula matrix contains strings
    if (is_expression_based()) {
          for (unsigned long r=0UL; r<hDim; r++)
              for (unsigned long c=0UL; c<vDim; c++) {
                  _Formula * entryFla = GetFormula(r,c);
                  if (entryFla) {
                      HBLObjectRef computedValue = FetchObjectFromFormulaByType (*entryFla, STRING);
                      if (computedValue) {
                          fillMe < new _StringBuffer (((_FString*)computedValue)->get_str());
                      } else {
                        fillMe.Clear();
                        return;
                      }
                  }
              }
    } else {
        if (convert_numbers && is_numeric()) {
            for (unsigned long r=0UL; r<hDim; r++) {
                for (unsigned long c=0UL; c<vDim; c++) {
                    fillMe.AppendNewInstance (new _String ((*this)(r,c)));
                }
            }
        }
    }
}

//_____________________________________________________________________________________________
HBLObjectRef   _Matrix::EvaluateSimple (_Matrix* existing_storage) {
// evaluate the matrix  overwriting the old one
    _Matrix * result;
    
    if (existing_storage && existing_storage->hDim == hDim && existing_storage->vDim == vDim && existing_storage->is_numeric() && ((bool)existing_storage->theIndex == (bool)theIndex)) {
        existing_storage->ZeroNumericMatrix();
        result = existing_storage;
    } else {
        if (existing_storage) {
            DeleteObject (existing_storage);
        }
        result = new _Matrix (hDim, vDim, bool (theIndex), true);
    }


    if (cmd->varIndex.lLength) {
        for (long i=0; i<cmd->varIndex.lLength; i++) {
            _Variable* curVar = LocateVar(cmd->varIndex.list_data[i]);
            if (curVar->ObjectClass () != MATRIX) {
                if (curVar->IsIndependent()) {
                    cmd->varValues[i].value = LocateVar (cmd->varIndex.list_data[i])->Value();
                } else {
                    cmd->varValues[i].value = LocateVar (cmd->varIndex.list_data[i])->Compute()->Value();
                }
            } else {
                cmd->varValues[i].reference = (hyPointer)((_Matrix*)LocateVar (cmd->varIndex.list_data[i])->Compute())->theData;
            }
        }
    }


    for (long f = 0; f < cmd->formulasToEval.lLength; f++) {
        cmd->formulaValues [f] = ((_Formula*)cmd->formulasToEval.list_data[f])->ComputeSimple(cmd->theStack, cmd->varValues);
    }

    long * fidx = cmd->formulaRefs;

    if (theIndex) {
        if (result->lDim != lDim) {
            result->lDim = lDim;
            result->theIndex = (long*)MemReallocate((hyPointer)result->theIndex,sizeof(long)*lDim);
            result->theData = (hyFloat*)MemReallocate ((hyPointer)result->theData,sizeof(hyFloat)*lDim);
        }
        result->bufferPerRow = bufferPerRow;
        result->overflowBuffer = overflowBuffer;
        result->allocationBlock = allocationBlock;


        for (long i = 0; i<lDim; i++) {
            long idx = theIndex[i];

            if (idx != -1) {
                result->theData[i] = cmd->formulaValues[fidx[i]];
            }

            result->theIndex[i] = idx;
        }


        if (hDim==vDim) {
            hyFloat* diagStorage = (hyFloat*)alloca (sizeof(hyFloat) * hDim);
            InitializeArray(diagStorage, hDim, 0.0);
            for (long i = 0; i<lDim; i++) {
                long k = result->theIndex[i];
                if (k!=-1) {
                    diagStorage[k/hDim] += result->theData[i];
                }
            }
            for (long i = 0; i<hDim; i++) {
                (*result)[i*hDim+i] = -diagStorage[i];
            }
        }
    } else {

        for (long i = 0; i<lDim; i++) {
            if (fidx[i]>= 0) {
                result->theData[i] = cmd->formulaValues[fidx[i]];
            }
        }

        if (hDim==vDim)
            for (long i = 0L, r = 0L; i<lDim; i+=vDim+1L, r++) {
                if (fidx[i] < 0) { // mod Aug 2 2005
                    //if (theFormulas[i]->IsEmpty())
                    //{
                    
                    hyFloat st = 0.;
                    long j;
                    
                    for (j = r*vDim; j<r*vDim+r; j++) {
                        st += result->theData[j];
                    }

                    for (j = r*vDim+r+1; j<(r+1)*vDim; j++) {
                        st += result->theData[j];
                    }

                    result->theData[i] = -st;
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

void    _Matrix::Clear (bool complete) {
    DeleteObject (theValue);
    if (is_expression_based()) { // has formulas in it - must delete
        ClearFormulae();
    }
    if (is_polynomial()) { // has objects in it - must delete
        ClearObjects();
    }
        
    if (theIndex) {
        if (complete) {
            MatrixMemFree (theIndex);
            theIndex = nil;
        } else {
            InitializeArray(theIndex, lDim, -1L);
        }
    }
    if (theData) {
        if (complete) {
            MatrixMemFree (theData);
            hDim = vDim = 0;
            theData = nil;
        } else {
            memset (theData, 0, lDim * (is_numeric() ? sizeof (hyFloat): sizeof (void*)));
        }
    }
    if (compressedIndex) {
        MatrixMemFree (compressedIndex);
        compressedIndex = nil;
    }

}

//_____________________________________________________________________________________________

void    _Matrix::ZeroNumericMatrix (void) {
    if (is_numeric()) {
        InitializeArray (theData, lDim, 0.0);
        if (!is_dense()) {
            InitializeArray (theIndex, lDim, -1L);
        }
    }
}

//_____________________________________________________________________________________________

void    _Matrix::Resize (long newH)
{
    if (newH >= 0 && newH != hDim && storageType == 1 && theIndex == nil) {
        hDim = newH;
        lDim = newH*vDim;

        if (theData) {
            theData = (hyFloat*) MemReallocate ((hyPointer)theData,sizeof (hyFloat)*lDim);
        } else {
            theData = (hyFloat*) MemAllocate (sizeof (hyFloat)*lDim);
        }
    }
}

//_____________________________________________________________________________________________

_Matrix::~_Matrix (void) {
    _Matrix::Clear();
}

//_____________________________________________________________________________________________

_Matrix const&    _Matrix::operator = (_Matrix const& m) {
    // SLKP 20180917 : reuse memory if copying dense numeric matrices of the same dimension
    if (m.is_numeric() && is_numeric() && CanFreeMe() && m.theIndex == nil && theIndex == nil && m.GetHDim () == GetHDim () && GetVDim () == m.GetVDim()) {
      unsigned long i = 0UL;
      for (unsigned long r = 0UL; r < hDim; r++) {
        for (unsigned long c = 0UL; c < vDim; c++, i++) {
          theData[i] = m.theData[i];
        }
      }
    } else {
      Clear();
      DuplicateMatrix (this, &m);
    }
    return *this;
}

//_____________________________________________________________________________________________

_Matrix const&    _Matrix::operator = (_Matrix const* m) {
    //Clear();
    //DuplicateMatrix (this, m);
    *this = *m;
    return *this;
}


//_____________________________________________________________________________________________
hyFloat _Matrix::AbsValue (void) const{
    if (is_numeric() && (is_row() || is_column())) {
        hyFloat norm = 0.;

        this->ForEach ([&] (hyFloat&& value, unsigned long, long) -> void {norm += value * value;},
                 [&] (unsigned long index) -> hyFloat {return theData[index];});

        return sqrt(norm);
    }

    return 0.;
}

//_____________________________________________________________________________________________
HBLObjectRef _Matrix::Abs (HBLObjectRef cache)
{
    if (storageType == 1 && (hDim==1 || vDim == 1)) {
        return _returnConstantOrUseCache(AbsValue(), cache);
    }
    return _returnConstantOrUseCache(MaxElement(), cache);

}

//_____________________________________________________________________________________________

void    _Matrix::AddMatrix  (_Matrix& storage, _Matrix& secondArg, bool subtract)
// addition operation on matrices
// internal function

{

    // check matrix dimensions to ensure that they are addable
    if (!((hDim==secondArg.hDim)&&(storage.hDim==secondArg.hDim)&&(vDim==secondArg.vDim)&&(storage.vDim==secondArg.vDim))) {
        HandleApplicationError  (_String ("Incompatible dimensions when trying to add or subtract matrices: first argument was a ") & _String (hDim) & 'x'
                          & _String (vDim) & " matrix and the second was a "& _String (secondArg.hDim) & 'x'  & _String (secondArg.vDim) & " matrix.");
        return;
    }

    if (is_numeric()) {
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
                memcpy (storage.theData, theData, sizeof (hyFloat)*lDim);
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
            hyFloat * _hprestrict_ argData = secondArg.theData;
            hyFloat * _hprestrict_ stData  = storage.theData;
            
            long    upto = secondArg.lDim >> 4 << 4;
                       
            if (subtract) {
#ifdef  _SLKP_USE_AVX_INTRINSICS
            #define     CELL_OP(x) _mm256_storeu_pd (stData+x, _mm256_sub_pd (_mm256_loadu_pd (stData+x), _mm256_loadu_pd (argData+x)))
                
                for (long idx = 0; idx < upto; idx+=16) {
                    CELL_OP (idx);
                    CELL_OP (idx+4);
                    CELL_OP (idx+8);
                    CELL_OP (idx+12);
                }
#else
                for (long idx = 0; idx < upto; idx+=4) {
                    stData[idx]-=argData[idx];
                    stData[idx+1]-=argData[idx+1];
                    stData[idx+2]-=argData[idx+2];
                    stData[idx+3]-=argData[idx+3];
                }
#endif
            } else {
#ifdef  _SLKP_USE_AVX_INTRINSICS
            #define     CELL_OP(x) _mm256_storeu_pd (stData+x, _mm256_add_pd (_mm256_loadu_pd (stData+x), _mm256_loadu_pd (argData+x)))
                 
                 for (long idx = 0; idx < upto; idx+=16) {
                     CELL_OP (idx);
                     CELL_OP (idx+4);
                     CELL_OP (idx+8);
                     CELL_OP (idx+12);
                 }

 #else
                for (long idx = 0; idx < upto; idx+=4) {
                    stData[idx]+=argData[idx];
                    stData[idx+1]+=argData[idx+1];
                    stData[idx+2]+=argData[idx+2];
                    stData[idx+3]+=argData[idx+3];
                }
#endif
            }
            if (subtract)
                for (long idx = upto; idx < secondArg.lDim; idx++) {
                    stData[idx]-=argData[idx];
                 }
            else
                for (long idx = upto; idx < secondArg.lDim; idx++) {
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
                                if (h<0) { // kEmptyString slot in matrix 1
                                    storage.StoreObject (hb,secondArg.GetMatrixObject(i)->Minus());
                                } else {
                                    storage.StoreObject (hb, GetMatrixObject(h)->Sub (secondArg.GetMatrixObject(i)));
                                }
                            }
                    } else {
                        for (i = 0; i<secondArg.lDim; i++)
                            if (secondArg.IsNonEmpty(i)) {
                                long hb =secondArg.HashBack (i), h = Hash (hb/vDim, hb%vDim);
                                if (h<0) { // kEmptyString slot in matrix 1
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
                    HBLObjectRef tempP;
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
                    HBLObjectRef tempP;
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

bool    _Matrix::AddWithThreshold  (_Matrix& secondArg, hyFloat prec)
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
        hyFloat* argData = secondArg.theData, *stData = theData,
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
    AddMatrix (storage,secondArg,true);
}

//_____________________________________________________________________________________________

void    _Matrix::Multiply  (_Matrix& storage, hyFloat c)
// multiply a matrix by a scalar
// internal function

{
    if (is_numeric()) { // numbers
        hyFloat * _hprestrict_  destination = storage.theData;
        hyFloat const *  source      = theData;
            
        if (theIndex) {
            for (long k = 0L; k < lDim; k++)
                if (storage.theIndex[k] != -1) {
                    destination[k] = source[k]*c;
                }
        } else {
  #ifdef  _SLKP_USE_AVX_INTRINSICS
  #define                 CELL_OP(k) _mm256_storeu_pd (destination + k, _mm256_mul_pd(value_op, _mm256_loadu_pd (source+k)))
        long lDimM4 = lDim >> 4 << 4,
             k = 0;
            
        __m256d  value_op = _mm256_set1_pd (c);
         for (k = 0L; k < lDimM4; k+=16) {
             CELL_OP (k);
             CELL_OP (k+4);
             CELL_OP (k+8);
             CELL_OP (k+12);
        }
        for (; k < lDim; k++) {
            destination[k] = source[k]*c;
        }
            
  #else
            for (long k = 0L; k < lDim; k++) {
                destination[k] = source[k]*c;
            }
  #endif
        }
            
    } else {
        _Constant * cc = new _Constant (c);

        if (storageType == 2) {
            _String const    star ('*');

            for (long i=0; i<lDim; i++)
                if (IsNonEmpty (i)) {
                    long h       = HashBack (i);
                    _Formula * f = GetFormula (h/vDim,h%vDim);
                    f->GetList().AppendNewInstance (new _Operation (cc));
                    f->GetList().AppendNewInstance (new _Operation (star,2));
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

void    _Matrix::Multiply  (_Matrix& storage, _Matrix const& secondArg) const
// multiplication operation on matrices
// internal function
// storage is assumed to NOT be *this

{
    HBLObjectRef tempP, tempP2;

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
                unsigned long cumulativeIndex = 0UL;
                const unsigned long dimm4 = (vDim >> 2) << 2;

                const hyFloat * row = theData;
                hyFloat  * dest = storage.theData;
              

#ifndef _SLKP_SSE_VECTORIZATION_
            
              
              if (dimm4 == vDim) {
                InitializeArray (dest, lDim, 0.0);
                for (unsigned long c = 0UL; c < secondArg.vDim; c ++) {
                  
#ifdef  _SLKP_USE_AVX_INTRINSICS

                  if (vDim == 20UL) { // special case for amino-acids
                    
                    __m256d __attribute__ ((aligned (32))) col_buffer[5];
                    
                      hyFloat   quad1[4] __attribute__ ((aligned (32))),
                                quad2[4] __attribute__ ((aligned (32))),
                                quad3[4] __attribute__ ((aligned (32))),
                                quad4[4] __attribute__ ((aligned (32))),
                                quad5[4] __attribute__ ((aligned (32)));
                      

                                quad1 [0] = secondArg.theData[c];
                                quad1 [1] = secondArg.theData[c + 20UL];
                                quad1 [2] = secondArg.theData[c + 40UL];
                                quad1 [3] = secondArg.theData[c + 60UL];

                                quad2 [0] = secondArg.theData[c + 80UL];
                                quad2 [1] = secondArg.theData[c + 100UL];
                                quad2 [2] = secondArg.theData[c + 120UL];
                                quad2 [3] = secondArg.theData[c + 140UL];

                                quad3 [0] = secondArg.theData[c + 160UL];
                                quad3 [1] = secondArg.theData[c + 180UL];
                                quad3 [2] = secondArg.theData[c + 200UL];
                                quad3 [3] = secondArg.theData[c + 220UL];

                                quad4 [0] = secondArg.theData[c + 240UL];
                                quad4 [1] = secondArg.theData[c + 260UL];
                                quad4 [2] = secondArg.theData[c + 280UL];
                                quad4 [3] = secondArg.theData[c + 300UL];

                                quad5 [0] = secondArg.theData[c + 320UL];
                                quad5 [1] = secondArg.theData[c + 340UL];
                                quad5 [2] = secondArg.theData[c + 360UL];
                                quad5 [3] = secondArg.theData[c + 380UL];

                              col_buffer[0] = _mm256_load_pd (quad1);
                              col_buffer[1] = _mm256_load_pd (quad2);
                              col_buffer[2] = _mm256_load_pd (quad3);
                              col_buffer[3] = _mm256_load_pd (quad4);
                              col_buffer[4] = _mm256_load_pd (quad5);
                      
                    hyFloat const * p = theData;
                    hyFloat row [20];
                    for (unsigned long r = 0UL; r < 20UL; r ++, p += 20UL) {
#ifdef _SLKP_USE_FMA3_INTRINSICS
                        __m256d sum1 = _mm256_fmadd_pd (_mm256_loadu_pd(p), col_buffer[0], _mm256_mul_pd(_mm256_loadu_pd(p+4UL), col_buffer[1]));
                        __m256d sum2 = _mm256_fmadd_pd (_mm256_loadu_pd(p+8UL), col_buffer[2],
                                                        _mm256_fmadd_pd(_mm256_loadu_pd(p+12UL), col_buffer[3],
                                                        _mm256_mul_pd(_mm256_loadu_pd(p+16UL), col_buffer[4])));
                        
                        //dest[r*vDim + c] = _avx_sum_4  (_mm256_add_pd (sum1, sum2));
                        row[r] = _avx_sum_4  (_mm256_add_pd (sum1, sum2));
#else
                      __m256d r0 = _mm256_mul_pd(_mm256_loadu_pd(p), col_buffer[0]);
                      __m256d r1 = _mm256_mul_pd(_mm256_loadu_pd(p+4UL), col_buffer[1]);
                      __m256d r2 = _mm256_mul_pd(_mm256_loadu_pd(p+8UL), col_buffer[2]);
                      __m256d r3 = _mm256_mul_pd(_mm256_loadu_pd(p+12UL), col_buffer[3]);
                      __m256d r4 = _mm256_mul_pd(_mm256_loadu_pd(p+16UL), col_buffer[4]);
                      
                      __m256d s01 = _mm256_add_pd(r0, r1);
                      __m256d s23 = _mm256_add_pd(r2, r3);
                      __m256d s234 = _mm256_add_pd(s23, r4);
                      row[r] = _avx_sum_4 (_mm256_add_pd(s01, s234));
                       //dest[r*vDim + c] = _avx_sum_4 (_mm256_add_pd(s01, s234));
#endif
                    }
                      
                      
                    for (unsigned long r = 0UL, idx = c; r < 20UL; r++, idx += 20UL) {
                        dest[idx] = row[r];
                    }
                    continue;
                  }

#endif
                  /*
                   load a series of 4 consecutive elements from a column in the second matrix,
                   say c [] = [i,i+1,i+2,i+3: c]
                   
                   next, iterate over all rows in the first matrix, looking for matched consecutive
                   elements, e.g.
                   
                   r [] = [r: i,i+1,i+2,i+3]
                   
                   compute sum_{t=0..3} c[t] * r[t]
                   
                   add to the element (r,c) in the destination matrix
                   
                   */
                  
                    const unsigned long
                                            column_shift2 = secondArg.vDim << 1,
                                            column_shift3 = (secondArg.vDim << 1) + secondArg.vDim,
                                            column_shift4 = secondArg.vDim << 2;
                    
                    for (unsigned long i = 0UL, vector_index = c; i < secondArg.hDim; i += 4UL, vector_index += column_shift4) {
                      hyFloat c0 = secondArg.theData[vector_index],
                                 c1 = secondArg.theData[vector_index+secondArg.vDim],
                                 c2 = secondArg.theData[vector_index+column_shift2],
                                 c3 = secondArg.theData[vector_index+column_shift3];
                
                      for (unsigned long r = 0UL; r < hDim; r ++) {
                        
                        unsigned long element = r*vDim + i;
                        
                        hyFloat r0 = theData[element]   * c0,
                                   r1 = theData[element+1] * c1,
                                   r2 = theData[element+2] * c2,
                                   r3 = theData[element+3] * c3;
                        
                        r0 += r1;
                        r2 += r3;
                        dest[r*vDim + c] += r0 + r2;
                  
                      }
                   }
                }
              } else {
                  const unsigned long
                          column_shift2 = secondArg.vDim << 1,
                          column_shift3 = (secondArg.vDim << 1) + secondArg.vDim,
                          column_shift4 = secondArg.vDim << 2;

                  for (unsigned long i=0UL; i<hDim; i++, row += vDim) {
                      for (unsigned long j=0UL; j<secondArg.vDim; j++) {
                          hyFloat resCell  = 0.0;

                          unsigned long k = 0UL,
                                       column = j;
                        
                          
                          for (; k < dimm4; k+=4, column += column_shift4) {
                              hyFloat pr1 = row[k]   * secondArg.theData [column],                         
                                         pr2 = row[k+1] * secondArg.theData [column + secondArg.vDim ],      
                                         pr3 = row[k+2] * secondArg.theData [column + column_shift2],
                                         pr4 = row[k+3] * secondArg.theData [column + column_shift3];
                            
                              pr1 += pr2;
                              pr3 += pr4;
                            
                              resCell += pr1 + pr3;
                          }
                          
                          for (; k < vDim; k++, column += secondArg.vDim) {
                              resCell += row[k] * secondArg.theData[column];
                          }
                        
                          dest[cumulativeIndex++] = resCell;
                     }
                  }
              }
              
#else
                secondArg.Transpose();
                for (long i=0; i<hDim; i++, row += vDim) {
                    for (long j=0; j<hDim; j++) {
                        hyFloat resCell  = 0.0;
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
#define _HY_MATRIX_CACHE_BLOCK 128
                 if (vDim >= 256) {
                     long nt = 1;
#ifdef _OPENMP
                      #define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)
#ifdef __HYPHYMPI__
                     if (hy_mpi_node_rank == 0)
                         
#endif
                     nt           = MIN(omp_get_max_threads(),secondArg.vDim / _HY_MATRIX_CACHE_BLOCK + 1);
#endif
                     for (long r = 0; r < hDim; r ++) {
#ifdef _OPENMP
  #if _OPENMP>=201511
    #pragma omp parallel for default(none) shared(r,secondArg,storage) schedule(monotonic:guided) proc_bind(spread) if (nt>1)  num_threads (nt)
  #else
    #if _OPENMP>=200803
      #pragma omp parallel for default(none) shared(r,secondArg,storage) schedule(guided) proc_bind(spread) if (nt>1)  num_threads (nt)
    #endif
  #endif
#endif
                         for (long c = 0; c < secondArg.vDim; c+= _HY_MATRIX_CACHE_BLOCK) {
                             hyFloat cacheBlockInMatrix2 [_HY_MATRIX_CACHE_BLOCK][_HY_MATRIX_CACHE_BLOCK];
                             const long upto_p = (secondArg.vDim-c>=_HY_MATRIX_CACHE_BLOCK)?_HY_MATRIX_CACHE_BLOCK:(secondArg.vDim-c);
                             for (long r2 = 0; r2 < secondArg.hDim; r2+= _HY_MATRIX_CACHE_BLOCK) {
                                 const long upto_p2 = (secondArg.hDim-r2)>=_HY_MATRIX_CACHE_BLOCK?_HY_MATRIX_CACHE_BLOCK:(secondArg.hDim-r2);
                                 for (long p = 0; p < upto_p; p++) {
                                     for (long p2 = 0; p2 < upto_p2; p2++) {
                                         cacheBlockInMatrix2[p][p2] = secondArg.theData [(r2+p2)*secondArg.vDim+c+p];
                                     }
                                 }
                                 if (upto_p2 % 4 == 0) {
                                     for (long p = 0; p < upto_p; p++) {
                                         hyFloat updater = 0.;
                                         for (long p2 = 0; p2 < upto_p2; p2+=4) {
                                             hyFloat pr1 = theData[r*vDim + r2 + p2]*cacheBlockInMatrix2[p][p2],
                                                        pr2 = theData[r*vDim + r2 + p2+1]*cacheBlockInMatrix2[p][p2+1],
                                                        pr3 = theData[r*vDim + r2 + p2+2]*cacheBlockInMatrix2[p][p2+2],
                                                        pr4 = theData[r*vDim + r2 + p2+3]*cacheBlockInMatrix2[p][p2+3];
                                             pr1 += pr2;
                                             pr3 += pr4;
                                             updater += pr1 + pr3;
                                         }
                                         storage.theData[r*secondArg.vDim + c + p] += updater;
                                     } 
                                 } else
                                     for (long p = 0; p < upto_p; p++) {
                                         hyFloat updater = 0.;
                                         for (long p2 = 0; p2 < upto_p2; p2++) {
                                             updater += theData[r*vDim + r2 + p2]*cacheBlockInMatrix2[p][p2];
                                         }
                                         storage.theData[r*secondArg.vDim + c + p] += updater;
                                     } 
                             }
                         }
                     }
                     
                 } else {
                     
                     
                    if (vDim % 4) {
                        long mod4 = vDim-vDim%4;
                        for (long i=0; i<hDim; i++) {
                            for (long j=0; j<secondArg.vDim; j++) {
                                hyFloat resCell = 0.0;
                                long k = 0;
                                for (; k < mod4; k+=4) {
                                    resCell += theData[i*vDim + k] * secondArg.theData[k*secondArg.vDim + j] +
                                    theData[i*vDim + k + 1] * secondArg.theData[(k+1)*secondArg.vDim + j] +
                                    theData[i*vDim + k + 2] * secondArg.theData[(k+2)*secondArg.vDim + j] +
                                    theData[i*vDim + k + 3] * secondArg.theData[(k+3)*secondArg.vDim + j];
                                }
                                for (; k < vDim; k++) {
                                    resCell += theData[i*vDim + k] * secondArg.theData[k*secondArg.vDim + j];
                                }
                                
                                storage.theData[i*secondArg.vDim + j] = resCell;
                            }
                        }
                    } else {
                        for (long i=0; i<hDim; i++) {
                            for (long j=0; j<secondArg.vDim; j++) {
                                hyFloat resCell = 0.0;
                                for (long k = 0; k < vDim; k+=4) {
                                    resCell += theData[i*vDim + k] * secondArg.theData[k*secondArg.vDim + j] +
                                    theData[i*vDim + k + 1] * secondArg.theData[(k+1)*secondArg.vDim + j] +
                                    theData[i*vDim + k + 2] * secondArg.theData[(k+2)*secondArg.vDim + j] +
                                    theData[i*vDim + k + 3] * secondArg.theData[(k+3)*secondArg.vDim + j];
                                }
                                
                                
                                storage.theData[i*secondArg.vDim + j] = resCell;
                            }
                        }
                    }
                }
                
            }
        }

    } else if (theIndex && !secondArg.theIndex) { // sparse multiplied by non-sparse
        if (storageType == 1 && secondArg.storageType ==1) { // both numeric
            if ( vDim == hDim && secondArg.vDim==secondArg.hDim) { // both square and same dimension
                /*
                  break out a special case for universal code
                */
              
              if (vDim == 61) {
                
                if (compressedIndex) {
                    long currentXIndex = 0L;
                    hyFloat  * _hprestrict_ res               = storage.theData;
                    
                    for (long i = 0; i < hDim; i++) { // row in source
                      while (currentXIndex < compressedIndex[i]) {
                            long currentXColumn = compressedIndex[currentXIndex + hDim];
                            // go into the second matrix and look up all the non-zero entries in the currentXColumn row
                          
                            hyFloat value = theData[currentXIndex];
                            hyFloat  * _hprestrict_ secArg            = secondArg.theData  + currentXColumn*61;
                            #ifdef  _SLKP_USE_AVX_INTRINSICS
                                __m256d  value_op = _mm256_set1_pd (value);
                                #ifdef _SLKP_USE_FMA3_INTRINSICS
                                    #define                 CELL_OP(x) _mm256_storeu_pd (res+x, _mm256_fmadd_pd (value_op, _mm256_loadu_pd (secArg+x),_mm256_loadu_pd(res+x)))
                                #else
                                    #define                 CELL_OP(x) _mm256_storeu_pd (res+x,   _mm256_add_pd (_mm256_loadu_pd(res+x),    _mm256_mul_pd(value_op, _mm256_loadu_pd (secArg+x))))
                                #endif
                              
                                CELL_OP(0);CELL_OP(4);CELL_OP(8);CELL_OP(12);
                                CELL_OP(16);CELL_OP(20);CELL_OP(24);CELL_OP(28);
                                CELL_OP(32);CELL_OP(36);CELL_OP(40);CELL_OP(44);
                                CELL_OP(48);CELL_OP(52);CELL_OP(56);
                          
                            #else
                                for (unsigned long i = 0UL; i < 60UL; i+=4UL) {
                                    res[i]   += value * secArg[i];
                                    res[i+1] += value * secArg[i+1];
                                    res[i+2] += value * secArg[i+2];
                                    res[i+3] += value * secArg[i+3];
                                }
                            #endif
                            res[60]   += value * secArg[60];
                            currentXIndex ++;
                      }
                      res += vDim;

                  }
                    
                  
                    
                } else {
                     for (unsigned long k=0UL; k<lDim; k++) { // loop over entries in the sparse matrix
                      long m = theIndex[k];
                      if (m >= 0L) {
                        long i = ((unsigned long)m)%61;
                      
                        hyFloat  value                            = theData[k];
                        hyFloat  * _hprestrict_ res               = storage.theData    + (m-i);
                        hyFloat  * _hprestrict_ secArg            = secondArg.theData  + i*61;
                        
                        #ifdef  _SLKP_USE_AVX_INTRINSICS
                          __m256d  value_op = _mm256_set1_pd (value);
                        
                         #ifdef _SLKP_USE_FMA3_INTRINSICS
                            #define                 CELL_OP(x) _mm256_storeu_pd (res+x, _mm256_fmadd_pd (value_op, _mm256_loadu_pd (secArg+x),_mm256_loadu_pd(res+x)))
                          #else
                            #define                 CELL_OP(x) _mm256_storeu_pd (res+x,   _mm256_add_pd (_mm256_loadu_pd(res+x),    _mm256_mul_pd(value_op, _mm256_loadu_pd (secArg+x))))
                          #endif
                          CELL_OP(0);CELL_OP(4);CELL_OP(8);CELL_OP(12);
                          CELL_OP(16);CELL_OP(20);CELL_OP(24);CELL_OP(28);
                          CELL_OP(32);CELL_OP(36);CELL_OP(40);CELL_OP(44);
                          CELL_OP(48);CELL_OP(52);CELL_OP(56);
                        #else
                          for (unsigned long i = 0UL; i < 60UL; i+=4UL) {
                            res[i]   += value * secArg[i];
                            res[i+1] += value * secArg[i+1];
                            res[i+2] += value * secArg[i+2];
                            res[i+3] += value * secArg[i+3];
                           }
                        #endif
                          res[60]   += value * secArg[60];
                      }
                     }
                }
                
              } else {
                  long loopBound = (vDim >> 2) << 2;
                  
                  if (compressedIndex) {
                    long currentXIndex = 0L;
                    hyFloat  * _hprestrict_ res               = storage.theData;
                    
                    for (long i = 0; i < hDim; i++) { // row in source
                      while (currentXIndex < compressedIndex[i]) {
                            long currentXColumn = compressedIndex[currentXIndex + hDim];
                            // go into the second matrix and look up all the non-zero entries in the currentXColumn row
                          
                            hyFloat value = theData[currentXIndex];
                            hyFloat  * _hprestrict_ secArg            = secondArg.theData  + currentXColumn*vDim;
                            #ifdef  _SLKP_USE_AVX_INTRINSICS
                                __m256d  value_op = _mm256_set1_pd (value);
                            #endif
                            for (unsigned long i = 0UL; i < loopBound; i+=4) {
                                  #ifdef  _SLKP_USE_AVX_INTRINSICS
                                  #ifdef _SLKP_USE_FMA3_INTRINSICS
                                        _mm256_storeu_pd (res+i, _mm256_fmadd_pd (value_op, _mm256_loadu_pd (secArg+i),_mm256_loadu_pd(res+i)));
                                  #else
                                        _mm256_storeu_pd (res+i, _mm256_add_pd (_mm256_loadu_pd(res+i),  _mm256_mul_pd(value_op, _mm256_loadu_pd (secArg+i))));
                                  #endif
                                  #else
                                    res[i]   += value * secArg[i];
                                    res[i+1] += value * secArg[i+1];
                                    res[i+2] += value * secArg[i+2];
                                    res[i+3] += value * secArg[i+3];
                                  #endif
                            }
                            for (unsigned long i = loopBound; i < vDim; i++) {
                                res[i]   += value * secArg[i];
                            }
                            currentXIndex ++;
                      }
                      res += vDim;

                 }
             } else {
                      for (unsigned long k=0UL; k<lDim; k++) { // loop over entries in the sparse matrix
                          long m = theIndex[k];
                          if  (m != -1L ) { // non-zero
                              long i = ((unsigned long)m)%vDim;
                              // this element will contribute to (r, c' = [0..vDim-1]) entries in the result matrix
                              // in the form of A_rc * B_cc'

                              hyFloat  value                           = theData[k];
                              hyFloat  *_hprestrict_ res               = storage.theData    + (m-i);
                              hyFloat  *_hprestrict_ secArg            = secondArg.theData  + i*vDim;
                              #ifdef  _SLKP_USE_AVX_INTRINSICS
                                    __m256d  value_op = _mm256_set1_pd (value);
                              #endif
                              
                              for (unsigned long i = 0UL; i < loopBound; i+=4) {
                              #ifdef  _SLKP_USE_AVX_INTRINSICS
                                #ifdef _SLKP_USE_FMA3_INTRINSICS
                                      _mm256_storeu_pd (res+i, _mm256_fmadd_pd (value_op, _mm256_loadu_pd (secArg+i),_mm256_loadu_pd(res+i)));
                                #else
                                      _mm256_storeu_pd (res+i, _mm256_add_pd (_mm256_loadu_pd(res+i),  _mm256_mul_pd(value_op, _mm256_loadu_pd (secArg+i))));
                                #endif
                                #else
                                  res[i]   += value * secArg[i];
                                  res[i+1] += value * secArg[i+1];
                                  res[i+2] += value * secArg[i+2];
                                  res[i+3] += value * secArg[i+3];
                                #endif
                              }
                               for (unsigned long i = loopBound; i < vDim; i++) {
                                  res[i]   += value * secArg[i];
                              }

                          }
                      }
                 }
              } // special codon case
            } else {
                for (long k=0; k<lDim; k++) {
                    long m = theIndex[k];
                    if (m!=-1) {
                        long i = ((unsigned long)m)/vDim;
                        long j = m - i*vDim;
                        hyFloat c = theData[k];
                        hyFloat* stData = storage.theData+i*secondArg.vDim,
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

                        hyFloat c = secondArg.theData[k];

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
                        hyFloat c = secondArg.theData[k];
                        hyFloat *stData = storage.theData+j,
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
    } else {
        //sparse by sparse
        
        /**
             X Y where both X and Y are sparse can be multipled more efficiently than O (N^3)
             if cell (i,j) is non-zero in X, it will contribute to cells (i,k) in the result matrix,
             it will contribute to the cells of the result matrix (i,k) where k is such that there are non-zero entries in thh k-th row of matrix Y
        */
        
        if (compressedIndex && secondArg.compressedIndex) {
            long currentXIndex = 0,
                 storageIndex  = 0;
            
            for (long i = 0; i < hDim; i++) { // row in source
                
                while (currentXIndex < compressedIndex[i]) {
                    long currentXColumn = compressedIndex[currentXIndex + hDim];
                    // go into the second matrix and look up all the non-zero entries in the currentXColumn row
                    hyFloat c = theData[currentXIndex];
                    long from = currentXColumn ? secondArg.compressedIndex[currentXColumn-1] : 0,
                          to   = secondArg.compressedIndex[currentXColumn];
                    for (long secondIndex = from; secondIndex < to; secondIndex ++) {
                        storage.theData[storageIndex + secondArg.compressedIndex[secondIndex + secondArg.hDim]] += c*secondArg.theData[secondIndex];
                    }
                    currentXIndex ++;
                }
                
                storageIndex += secondArg.vDim;

            }
        } else {
        
            long * indexVector = (long*)alloca ( sizeof(long)*secondArg.hDim);
                // how many non-zero elements are there in the i-th row of matrix Y
            memset (indexVector,0,secondArg.hDim*sizeof(long));
            
            long *indexTable,
                 *indexTable2,
                 indexTableColumnWidth = secondArg.vDim,
                 indexTableDim = secondArg.hDim*indexTableColumnWidth;

            indexTable  = (long*)MatrixMemAllocate( sizeof(long)*indexTableDim);
              // element (i,j) of this matrix is the COLUMN index in which the j-th non-zero entry appears in row i of matrix Y
            indexTable2 = (long*)MatrixMemAllocate( sizeof(long)*indexTableDim);
              // element (i,j) of this matrix is the index (in theData) for the j-th non-zero entry appears in row i of matrix Y

            memset (indexTable,0,indexTableDim*sizeof(long));
            memset (indexTable2,0,indexTableDim*sizeof(long));
             
            if (is_numeric()) {
                // numeric
                
                for (long i=0; i<secondArg.lDim; i++) {
                    long elementIndex = secondArg.theIndex[i];
                    if (__builtin_expect (elementIndex >= 0L,1)) {
                        long k = elementIndex/secondArg.vDim;
                        long j = k*secondArg.vDim+(indexVector[k]++);
                        indexTable [j] = elementIndex % secondArg.vDim;
                        indexTable2[j] = i;
                    }
                }
                for (long k=0; k<lDim; k++) {
                    long elementIndex = theIndex[k];
                    if (__builtin_expect (elementIndex >= 0L,1)) {
                        hyFloat c = theData[k];
                        long i = elementIndex / vDim,
                             j = elementIndex - i*vDim,
                             n = j*indexTableColumnWidth,
                             m = i*secondArg.vDim,
                             nonZeroCount = indexVector[j];
                        
                        for (long l=n; l<n+nonZeroCount; l++) {
                            storage.theData[m+indexTable[l]] += c*secondArg.theData[indexTable2[l]];
                        }
                    }
                }
            } else { // polynomial entries
                 
                for (long i=0; i<secondArg.lDim; i++) {
                    long elementIndex  = secondArg.theIndex[i];
                    if (IsNonEmpty(i)) {
                        long k = elementIndex/secondArg.vDim;
                        long j = k*secondArg.vDim+(indexVector[k]++);
                        indexTable [j] = elementIndex%secondArg.vDim;
                        indexTable2[j] = i;
                    }
                }
                for (long k=0; k<lDim; k++) {
                    if (IsNonEmpty(k)) {
                        long i = theIndex[k]/vDim;
                        long j = theIndex[k]%vDim;
                        _MathObject* p = GetMatrixObject(k);
                        long n = j*secondArg.vDim;
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
        }
    }


}

//_____________________________________________________________________________________________

long    _Matrix::HashBack  (long logicalIndex) const {
// returns element's matrix index in the form vDim*(i-1)+j, where i, j are the matrix coordinates
// given a buffer index
    return theIndex?theIndex [logicalIndex]:logicalIndex;
}

//_____________________________________________________________________________________________

hyFloat  _Matrix::MaxElement  (char runMode, long* indexStore) const {
// returns matrix's largest abs value element
    if (storageType == _NUMERICAL_TYPE) {
        hyFloat max  = 0.0,
                   temp;

        bool doAbsValue = runMode != 1 && runMode != 3,
             doMaxElement = runMode == 0 || runMode == 3;

        if (doMaxElement) {
            max = -INFINITY;
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

void    _Matrix::RowAndColumnMax  (hyFloat& r, hyFloat &c, hyFloat * cache)

// returns the maximum row sum / column sum
// the cache must be big enough to hold hDim + vDim
// leave as nil to allocate cache run time

{
    r = c = 10.;

    if (is_numeric()) { // numeric matrix
        hyFloat  *maxScratch = cache;
        r = c = 0.;

        if (maxScratch == nil) {
            maxScratch = (hyFloat*)MemAllocate ((hDim+vDim)*sizeof(hyFloat), true);
        } else
            InitializeArray(maxScratch, hDim + vDim, 0.0);

        hyFloat * rowMax = maxScratch,
                     * colMax = maxScratch + hDim;

        if (theIndex)
            // sparse matrix
            for (long i = 0; i<lDim; i++) {
                long k = theIndex[i];
                if  (k!=-1) {
                    hyFloat temp = theData[i];

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
                    hyFloat temp = theData[k];
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

bool    _Matrix::IsMaxElement  (hyFloat bench) {
// returns matrix's largest abs value element
    if (is_numeric()) {
        
        hyFloat mBench = -bench;
        if (!theIndex || compressedIndex) {
            for (long i = 0; i<lDim; i++) {
                hyFloat t = theData[i];
                if ( t<mBench || t>bench ) {
                    return true;
                }
            }
        } else {
            for (long i = 0; i<lDim; i++) {
                if (theIndex [i] >= 0) {
                    hyFloat t = theData[i];
                    if ( t<mBench || t>bench ) {
                        return true;
                    }
                }
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

hyFloat  _Matrix::MaxRelError  (_Matrix& compMx)
// returns matrix's largest abs value element
{
    if (storageType == 1) {
        hyFloat max = 0, temp;
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

hyFloat  _Matrix::MinElement  (char doAbsValue, long* storeIndex)
// returns matrix's smalles non-zero abs value element
{
    if (storageType == 1) {
        hyFloat min = DBL_MAX;

        if (theIndex)
            for (long i = 0; i<lDim; i++) {
                if (theIndex[i] < 0) {
                    continue;
                }

                hyFloat temp = theData[i];

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
                hyFloat temp = theData[i];

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
                        hyFloat z      = theData[i*vDim+j];
                        theData[i*vDim+j] = theData[j*vDim+i];
                        theData[j*vDim+i] = z;
                    }
            } else { // sparse
                for (long i = 0; i<lDim; i++) {
                    long p = theIndex[i];
                    if (p!=-1) {
                        long k      = p/vDim;
                        long l      = p - k*vDim;

                        if (l!=k) { // off - diag
                            p            = Hash (l,k);
                            hyFloat z = theData[i];
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
        hyPointer z;
        if (hDim == vDim) {
            if (!theIndex) { // non-sparse
                for (long i = 0; i<hDim; i++)
                    for (long j = i+1; j<vDim; j++) {
                        if (storageType==2) {
                            z = (hyPointer)GetFormula(i,j);
                        } else {
                            z = (hyPointer)GetMatrixObject(i*vDim+j);
                        }

                        if (storageType==2) {
                            ((_Formula**)theData)[i*vDim+j] = GetFormula(j,i);
                        } else {
                            ((HBLObjectRef*)theData)[i*vDim+j] = GetMatrixObject(j*vDim+i);
                        }

                        ((hyPointer*)theData)[j*vDim+i] = z;
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
                                z = (hyPointer)GetFormula(k,l);
                            } else {
                                z = (hyPointer)GetMatrixObject(i);
                            }

                            if (p>=0) {
                                if (storageType==2) {
                                    ((_Formula**)theData)[i]    = GetFormula(l,k);
                                } else {
                                    ((_MathObject**)theData)[i] = GetMatrixObject(p);
                                }

                                ((hyPointer*)theData)[p] = z;
                            } else {
                                theIndex[i]=-1;
                                if (storageType==2) {
                                    StoreFormula(l,k,*(_Formula*)z,false,false);
                                } else {
                                    StoreObject(l*vDim+k,(HBLObjectRef)z);
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
                            z =   (hyPointer)GetMatrixObject (i*vDim+j);
                            result.StoreObject(j*hDim+i,(HBLObjectRef)z);
                            ((HBLObjectRef)z)->AddAReference();
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
                            z =   (hyPointer)GetMatrixObject (i);
                            result.StoreObject(c*hDim+r,(HBLObjectRef)z);
                            ((HBLObjectRef)z)->AddAReference();
                        }
                    }
            }
            Swap (result);
            //*this = result;
        }
    }
}

//_____________________________________________________________________________________________

void    _Matrix::CompressSparseMatrix (bool transpose, hyFloat * stash) {
    /**
        The purpose of this function is to "pack" the elements of a sparse matrix into an index array (theIndex) so that
                    (1) there are no "unused" elements in the middle
                    (2) the elements in theIndex are in the same order as they are in the matrix (row by row)
     
        Optonally, the packe matrix can be transposed.
    */
    
    
    if (theIndex) {
        _SimpleList sortedIndex  ((unsigned long)lDim, (long*)alloca (lDim * sizeof (long))),
                    sortedIndex3 ((unsigned long)lDim, (long*)alloca (lDim * sizeof (long))),
                    sortedIndex2;


        /*const long blockChunk = 32,
                   blockShift = hDim / blockChunk + 1*/
        
        long  max        = 0,
              secondDim = transpose ? hDim : vDim,
              firstDim  = transpose ? vDim : hDim;


        if (transpose) {
            for (long i2=0; i2<lDim; i2++) {
                long k = theIndex[i2];
                if  (__builtin_expect (k!=-1,1)) {
                    /*
                    long r = k/vDim,
                         c = k%vDim,
                         r2 = c / blockChunk * blockShift + r / blockChunk,
                         r3 = r2 * lDim + r * vDim + c;

                    sortedIndex  << (c*vDim + r);
                    sortedIndex3 << r3;
                    stash[sortedIndex.lLength-1] = theData[i2];
                    if (r3 > max) {
                        max = r3;
                    }*/
                    long transposed_index = (k%vDim) * vDim + k/vDim;
                    stash[sortedIndex.lLength] = theData[i2];
                    sortedIndex  << transposed_index;
                    sortedIndex3 << transposed_index;
                    if (max < transposed_index) max = transposed_index;
                    
                }
            }
        } else {
            for (long i2=0; i2<lDim; i2++) {
                long k = theIndex[i2];
                if  (__builtin_expect (k!=-1,1)) {
                    /*long r = k%vDim,
                         c = k/vDim,
                         r2 = c / blockChunk * blockShift + r / blockChunk,
                         r3 = r2 * lDim + r * vDim + c;

                    sortedIndex  << (c*vDim + r);
                    sortedIndex3 << r3;
                    stash[sortedIndex.lLength-1] = theData[i2];
                    if (r3 > max) {
                        max = r3;
                    }*/
                    stash[sortedIndex.lLength] = theData[i2];
                    sortedIndex  << k;
                    sortedIndex3 << k;
                    if (max < k) max = k;
                }
            }
        }
        

        //printf ("_Matrix::CompressSparseMatrix %d\n%s\n%ld\n", transpose, _String ((_String*)sortedIndex3.toStr()).get_str(),max);
        
        /*if (max > (lDim<<4)) {
            sortedIndex2. Populate(sortedIndex.lLength,0,1);
            SortLists(&sortedIndex3,&sortedIndex2);
        } else {*/
            sortedIndex3.CountingSort(max+1, &sortedIndex2,false);
        //}
        
       lDim = sortedIndex.lLength;

       /*for (long i=0; i<sortedIndex.lLength; i++) {
           //printf ("%ld %ld", i, theIndex[i]);
           theIndex[i] = sortedIndex.list_data[sortedIndex2.list_data[i]];
           theData[i]  = stash[sortedIndex2.list_data[i]];
       }*/

        if (compressedIndex) {
            MatrixMemFree(compressedIndex);
        }
        
        compressedIndex = (long*) MatrixMemAllocate((lDim + firstDim) * sizeof (long));
        
        long              currentRow = 0;
        lDim = sortedIndex.lLength;

        for (long i=0; i<sortedIndex.lLength; i++) {
            //printf ("%ld %ld", i, theIndex[i]);
            long entryIndex = sortedIndex.list_data[sortedIndex2.list_data[i]];
            theIndex[i] = entryIndex;
            
            long indexRow = entryIndex / secondDim,
                 indexColumn = entryIndex - indexRow * secondDim;
 
            //printf ("[%ld] %ld %ld %ld\n", i, indexRow, indexColumn, currentRow);

            compressedIndex[i + firstDim] = indexColumn;
            if (indexRow > currentRow) {
                for (long l = currentRow; l < indexRow; l++) {
                    compressedIndex[l] = i;
                }
                //printf (">[%ld] %ld\n", currentRow, compressedIndex[currentRow]);
                currentRow = indexRow;
            } /*else {
                if (currentRow > indexRow) {
                    printf ("\n\n\nBARF\n\n\n");
                }
            }*/
            
            //printf (" %ld\n", theIndex[i]);
            theData[i]  = stash[sortedIndex2.list_data[i]];
        }
        compressedIndex[firstDim-1] = lDim;
        //printf (">[%ld] %ld\n", firstDim-1, compressedIndex[firstDim-1]);
        //exit (1);
        
        
        // this code checks conversion consistency
         
        /*long ei = 0;
        long from = 0;
        for (long i=0; i<hDim; i++) {
            for (long k = from; k < compressedIndex[i]; k++, ei++) {
                long my_row = theIndex[ei] / vDim,
                     my_column = theIndex[ei] % vDim;
                
                if (my_row != i || my_column != compressedIndex[hDim+k]) {
                    printf ("BARF\n");
                }
            }
            from = compressedIndex[i];
        }
        
        if (ei != lDim) {
            printf ("BARF\n");
        }*/
        

    }
}


//_____________________________________________________________________________________________

_Matrix*    _Matrix::Exponentiate (hyFloat scale_to, bool check_transition, _Matrix * existing_storage) {
    // find the maximal elements of the matrix
    
    
    try {
        if (!is_square()) {
            throw _String ("Exponentiate is not defined for non-square matrices");
        }
        
        long i,
             power2 = 0L;
        
#ifndef _OPENMP
        matrix_exp_count++;
#endif
        
        hyFloat max     = 1.0,
                *stash,
                *stash2 = 0;
        //  = new hyFloat[hDim*(1+vDim)];
        
        if (!is_polynomial()) {
            stash = (hyFloat*)alloca(sizeof (hyFloat) * hDim*(1+vDim));
            if (theIndex) {
                // transpose sparse matrix
                CompressSparseMatrix (true,stash);
            }

            hyFloat t;
            //bool    censor = false;
            RowAndColumnMax (max, t, stash);
            /*if (t > 10000. || max > 10000.) {
                censor = true;
                t   = MIN (t,10000.);
                max = MIN (max, 10000.);
            }*/
            max *= t;
            //max = MaxElement();
            //max = max * max;
            if (max > .1) {
                max             = scale_to*sqrt (10.*max);
                power2          = (long)((log (max)/_log2))+1L;
                max             = exp (power2 * _log2);
                stash2 = (hyFloat*)alloca(sizeof (hyFloat) * lDim);
                memcpy   (stash2, theData, sizeof (hyFloat) * lDim);
                /*if (censor) {
                    for (long i = 0; i < lDim; i++) {
                        if (theData[i] < -10000.) {
                            theData[i] = -10000.;
                        } else {
                            if (theData[i] > 10000.) {
                                theData[i] = 10000;
                            }
                        }
                    }
                }*/
                (*this)         *= 1.0/max;
            } else {
                power2 = 0;
            }
            
            
        } else {
            max = 1.;
        }
       
        _Matrix *result;
        
        if (!is_polynomial() && existing_storage && existing_storage->hDim == hDim && existing_storage->vDim == vDim && existing_storage->is_numeric() && existing_storage->is_dense()) {
            result = existing_storage;
            InitializeArray(result->theData, result->lDim, 0.0);
        } else {
            result = new _Matrix(hDim, vDim , is_polynomial(), !is_polynomial());
        }
    
        
        // put ones on the diagonal
        
        if (!is_polynomial()) {
            long step = vDim + 1;
            for (long diag = 0; diag < result->lDim; diag += step) {
                result->theData[diag] = 1.;
            }
        } else {
            for (i=0; i<(*result).hDim*(*result).vDim; i+=vDim+1) {
                (*result).StoreObject(i,new _Polynomial (1.),false);
            }
        }
        
        if (max == 0.0) {
            return result;
        }
        
        (*result) += (*this);
        
        i = 2;
        
        if (precisionArg || is_polynomial()) {
            _Matrix temp    (*this);

            if (!is_polynomial()) {
                for (; i<=precisionArg; i++) {
                    temp      *= (*this);
                    temp      *= 1.0/i;
                    (*result) += temp;
                }
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
            hyFloat tMax = MAX(MinElement()*sqrt ((hyFloat)hDim),truncPrecision);
            
            i=2;
            
            
            if (is_dense()) { // avoid matrix allocation
                _Matrix temp ;
                temp.hDim = hDim;
                temp.vDim = vDim;
                temp.lDim = lDim;
                temp.theData = (hyFloat*)alloca(sizeof (hyFloat) * hDim*vDim);
                memcpy (temp.theData, theData, sizeof (hyFloat) * hDim*vDim);

                _Matrix tempS;
                tempS.hDim = hDim;
                tempS.vDim = vDim;
                tempS.lDim = lDim;
                tempS.theData = stash;
                // zero out the stash
                memset (stash, 0, sizeof (hyFloat)*lDim);
                do {
                    temp.MultbyS        (*this,false, &tempS, nil);
                    // after this call, temp and tempS are gonna swap pointers to theData
                    temp      *= 1.0/i;
                    (*result) += temp;
                    i         ++;
#ifndef _OPENMP
                    taylor_terms_count++;
#endif
                } while (temp.IsMaxElement(tMax*truncPrecision*i));
                
                
                if (tempS.theData != stash) { // need to copy data to tempS
                    Exchange (tempS.theData, temp.theData);
                    memcpy (temp.theData, tempS.theData, lDim * sizeof (hyFloat));
                }
                tempS.theData = nil;
                temp.theData = nil;
                
            } else  {
                _Matrix temp    (*this);
                _Matrix tempS (hDim, vDim, false, temp.storageType);
                do {
                    temp.MultbyS        (*this,theIndex!=nil, &tempS, stash);
                    temp      *= 1.0/i;
                    (*result) += temp;
                    i         ++;
    #ifndef _OPENMP
                    taylor_terms_count++;
    #endif
                } while (temp.IsMaxElement(tMax*truncPrecision*i));
            }
            
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
            //(*this)*=max;
            memcpy(theData, stash2, sizeof (hyFloat) * lDim);
            
        }
        
        if (theIndex) {
            // transpose back
            for (i=0; i<lDim; i++) {
                long k = theIndex[i];
                if  (k!=-1) {
                    long div = k / vDim;
                    theIndex[i] = (k - div * vDim)*vDim + div;
                }
            }
            result->Transpose();
        }
        
       
        //_Matrix stash_mx (*result);
        
        for (long s = 0; s<power2; s++) {
#ifndef _OPENMP
            squarings_count++;
#endif
           /* for (i = 0; i < hDim; i++) {
                if ((*result)(i,i) > 1.) {
                    printf ("\n%ld\n", s);
                    ObjectToConsole(result);
                    NLToConsole();
                    ObjectToConsole(&stash_mx);
                    NLToConsole();
                    _Matrix cp = *this;
                    cp.CheckIfSparseEnough(true);
                    ObjectToConsole(&cp);
                    abort();
                }
            }*/
            
            if (result->Sqr(stash) < DBL_EPSILON * 1.e3) {
                break;
            }
        }
        
        
        if (check_transition) {
            bool pass = true;
            if (result->is_dense()) {
                for (unsigned long r = 0L; r < result->lDim; r += result->vDim) {
                    if (result->theData[r] > 1.) {
                        pass = false;
                        break;
                    }
                }
            } else {
                for (unsigned long r = 0L; r < result->hDim; r ++) {
                    if ((*result)(r,r) > 1.) {
                        pass = false;
                        break;
                    }
                }
            }
            if (!pass) {
                if (scale_to < 1.e100) {
                    DeleteObject (result);
                    return this->Exponentiate(scale_to * 100, true);
                }
                
                /*printf ("SCALE %lg : \n", scale_to);
                
                for (unsigned long r = 0L; r < hDim; r ++) {
                    hyFloat sum = 0.;
                    printf ("%ld %18.16lg %18.16lg\n", r, (*result)(r,r), (*result)(r,r) - 1.);
                    for (unsigned long c = 0L; c < vDim; c++) {
                        sum += (*this)(r,c);
                    }
                    printf ("%ld %g\n", r, sum);
                }

                
                ObjectToConsole(this);
                ObjectToConsole(result);*/
                
                throw _String ("Failed to compute a valid transition matrix; this is usually caused by ill-conditioned rate matrices (e.g. very large rate values)");
            }
        }
        
        return result;
    }
    catch (const _String& e) {
        HandleApplicationError(e);
    }
    
    return new _Matrix;
    
}

//_____________________________________________________________________________________________

void     _Matrix::SetupSparseMatrixAllocations (void) {
    overflowBuffer = hDim*storageIncrement/100;
    bufferPerRow = MAX (1, (lDim-overflowBuffer)/hDim);
    overflowBuffer = lDim-bufferPerRow*hDim;
    allocationBlock = hDim*vDim*storageIncrement/100+1;
}

//_____________________________________________________________________________________________

long    _Matrix::Hash  (long i, long j) const {
// returns element's position in the buffer (-1 if not found)

    if (is_dense()) {
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
hyFloat      _Matrix::operator () (long i, long j) const {
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

        if (storageType == 2) // formulae
            for (long k=0; k<h->lLength; k++) {
                result->StoreFormula(column?k:0,column?0:k,*GetFormula(h->list_data[k],v->list_data[k]));
            }
        else
            for (long k=0; k<h->lLength; k++) {
                result->theData[k] = (*this)(h->list_data[k],v->list_data[k]);
            }

        return result;
    }
    return new _Matrix;
}



//_____________________________________________________________________________________________
HBLObjectRef _Matrix::MAccess (HBLObjectRef p, HBLObjectRef p2, HBLObjectRef cache) {
  if (!p) {
    HandleApplicationError ( kErrorStringInvalidMatrixIndex );
    return _returnConstantOrUseCache (0., cache);
  }
  
  if (hDim <= 0L || vDim <= 0L) {
    return _returnConstantOrUseCache (0., cache);
  }
  
  if (p->ObjectClass() == MATRIX) {
    if (p2 == nil) {
      _Matrix * nn = (_Matrix*)p;
      if (nn->storageType == 1) {
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
            
            for (unsigned long r=0UL; r<nn->hDim; r++) {
              long v = floor((*nn)(r,0L));
              if (v>=0L && v<hDim) {
                hL<<v;
              }
            }
            
            if (hL.lLength) {
              _Matrix * result = new _Matrix (hL.lLength,vDim,false,true);
              unsigned long k = 0UL;
              for (unsigned long r=0UL; r<hL.lLength; r++) {
                unsigned long ri = hL.list_data[r];
                for (unsigned long c=0UL; c<vDim; c++,k++) {
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
              long k = 0;
              for (long c=0; c<hDim; c++)
                for (long r=0; r<hL.lLength; r++,k++) {
                  result->theData[k] = (*this)(c,hL.list_data[r]);
                }
              return result;
            }
            
            return new _Matrix;
          }
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
       _String aFormulaString (((_FString*)p)->get_str());
      _Formula f (aFormulaString, currentExecutionList ? currentExecutionList->nameSpacePrefix : nil);
      
      if (!f.IsEmpty()) {
        /* check formula validity */
        
        
          _Variable * cv = CheckReceptacle(&hy_env::matrix_element_value, kEmptyString, false),
                    * cr = CheckReceptacle(&hy_env::matrix_element_row, kEmptyString, false),
                    * cc = CheckReceptacle(&hy_env::matrix_element_column, kEmptyString, false);
        
        cv->CheckAndSet (0.0, false, NULL);
        cr->CheckAndSet (0.0, false, NULL);
        cc->CheckAndSet (0.0, false, NULL);
        
        f.Compute();
        if (terminate_execution) {
          return new _Matrix ();
        } else {
          
          _Formula * conditionalCheck = nil;
          
          if (p2 && p2->ObjectClass() == STRING) {
            conditionalCheck = new _Formula (((_FString*)p2)->get_str(), currentExecutionList ? currentExecutionList->nameSpacePrefix : nil);
            if (conditionalCheck->IsEmpty()) {
              delete conditionalCheck;
              conditionalCheck = nil;
            }
            
            conditionalCheck->Compute();
            if (terminate_execution) {
              delete conditionalCheck;
              return new _Matrix ();
            }
          }
          
          _Matrix   * retMatrix = new _Matrix (hDim,vDim,false,true);
          
          long          stackDepth = 0;
          _SimpleList   vIndexAux;
          _AVLList      vIndex (&vIndexAux);
          
          if (f.AmISimple (stackDepth,vIndex) && (!conditionalCheck || conditionalCheck->AmISimple(stackDepth,vIndex))) {
            _SimpleFormulaDatum * stack     = new _SimpleFormulaDatum [stackDepth+1],
            * varValues = new _SimpleFormulaDatum [vIndex.countitems()];
            
            bool                constantValue = false;
            hyFloat          constantV     = f.Compute()->Value();
            
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
              
              long rid []= {cr->get_index(),cc->get_index(),cv->get_index()};
              
              for (long k=0; k<3; k++) {
                rid[k] = vIndexAux.Find(rid[k]);
              }
              
              PopulateArraysForASimpleFormula(vIndexAux, varValues);
              
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
                    //printf ("Formula eval (stack depth= %d) (%d, %g, %g) %g\n", stackDepth, rid[2], varValues[rid[2]], f.ComputeSimple(stack,varValues));
                    
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
              cr->CheckAndSet (r,false, NULL);
              for (long c=0; c<vDim; c++) {
                cc->CheckAndSet (c,false, NULL);
                  cv->CheckAndSet ((*this)(r,c),false, NULL);
                HBLObjectRef fv;
                
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
      ReportWarning (_String("Invalid formula expression for element-wise matrix operations: ") & ((_FString*)p)->get_str());
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
    return _returnConstantOrUseCache (0., cache);
  }
  
  if (ind2>=0) { // element access
    return GetMatrixCell (ind1, ind2, cache);
  }
  
  return _returnConstantOrUseCache (0., cache);
}

//_____________________________________________________________________________________________
HBLObjectRef _Matrix::GetMatrixCell (long ind1, long ind2, HBLObjectRef cache) const {
    if (is_expression_based()) { // formulas
      if (!theIndex) {
        _Formula * entryFla = (((_Formula**)theData)[ind1*vDim+ind2]);
        if (entryFla) {
          return (HBLObjectRef)entryFla->Compute()->makeDynamic();
        } else {
          return _returnConstantOrUseCache (0., cache);
        }
      } else {
        long p = Hash (ind1, ind2);
        if (p<0) {
          return _returnConstantOrUseCache (0., cache);
        } else {
          return (HBLObjectRef)(((_Formula**)theData)[p])->Compute()->makeDynamic();
        }
      }
    } else {
      if (is_numeric()) {
        if (theIndex) {
          return _returnConstantOrUseCache ((*this)(ind1,ind2), cache);
        } else {
          return _returnConstantOrUseCache (theData[ind1*vDim+ind2], cache);
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

    
//_____________________________________________________________________________________________
_Formula* _Matrix::GetFormula (long ind1, long ind2) const {

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
HBLObjectRef _Matrix::MCoord (HBLObjectRef p, HBLObjectRef p2, HBLObjectRef cachedResult) {
    long ind1 = -1L,
         ind2 = -1L;

    if (!p) {
        HandleApplicationError ( kErrorStringInvalidMatrixIndex );
        return new _MathObject;
    }

    ind1 = p->Value();
    if (p2) {
        ind2 = p2->Value();
    }


    if (hDim == 1L) {
        if (ind2<0L) {
            ind2 = ind1;
        }
        ind1=0L;
    }

    if (vDim == 1L) {
        ind2 = 0L;
    }

    if (ind2<0L) { // allow direct vectorlike indexing, i.e m[21] = m[3][3] (if the dim is *x6)
        ind2 = ind1%vDim;
        ind1 = ind1/vDim;
    }
    _Matrix * res = nil;
    if (cachedResult && cachedResult->ObjectClass() == MATRIX) {
        res = (_Matrix*)cachedResult;
        if (!(res->is_numeric() && res->check_dimension(1,2))){
            res = nil;
        }
    }
    if (!res)
        res = new _Matrix (1L,2L,false,true);

    res->theData[0]=ind1;
    res->theData[1]=ind2;
    return res;

}

//_____________________________________________________________________________________________
bool _Matrix::MResolve (HBLObjectRef p, HBLObjectRef p2, long& ind1, long& ind2)
{
    ind1 = -1;
    ind2 = -1;

    if (!p) {
        HandleApplicationError ( kErrorStringInvalidMatrixIndex );
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
void _Matrix::MStore (long ind1, long ind2, _Formula& f, long opCode) {
    if (ind2>=0) { // element storage
        if (is_expression_based()) { // formulas
            if (opCode == HY_OP_CODE_ADD) {
                _Formula * addOn = GetFormula(ind1,ind2);
                if (addOn) {
                    StoreFormula (ind1,ind2,*_Formula::PatchFormulasTogether(*addOn, f, HY_OP_CODE_ADD),false);
                    return;
                }
            } 
            StoreFormula (ind1,ind2,f);
        } else {
            if (!f.IsAConstant()) {
                Convert2Formulas();
                StoreFormula (ind1,ind2,f);
            } else {
                HBLObjectRef res = f.Compute();
                hyFloat toStore = res->Value();
                if (opCode == HY_OP_CODE_ADD) {
                    toStore += (*this)(ind1,ind2);
                }
                Store(ind1,ind2,toStore);
            }
        }
    }
}

//_____________________________________________________________________________________________
void _Matrix::MStore (long ind1, long ind2, HBLObjectRef value, long opCode) {
    if (ind2>=0) { // element storage
        if (is_expression_based()) { // formulas
            value->AddAReference();
            if (opCode == HY_OP_CODE_ADD) {
                _Formula * addOn = GetFormula(ind1,ind2);
                if (addOn) {
                    StoreFormula (ind1,ind2,*_Formula::PatchFormulasTogether(*addOn, value, HY_OP_CODE_ADD),false);
                    return;
                }
            }
            _Formula * f = new _Formula (value);
            StoreFormula (ind1,ind2,*f,false);
        } else {
            if (value->ObjectClass() != NUMBER) {
                Convert2Formulas();
                value->AddAReference();
                _Formula * f = new _Formula (value);
                StoreFormula (ind1,ind2,*f,false);
            } else {
                hyFloat toStore = value->Value();
                if (opCode == HY_OP_CODE_ADD) {
                    toStore += (*this)(ind1,ind2);
                }
                Store(ind1,ind2,toStore);
            }
        }
    }
}

//_____________________________________________________________________________________________
void _Matrix::MStore (HBLObjectRef p, HBLObjectRef p2, _Formula& f, long opCode)
{
    long      ind1, ind2;
    if (MResolve (p,p2, ind1,ind2)) {
        MStore   (ind1,ind2,f, opCode);
    }
}

//_____________________________________________________________________________________________
void _Matrix::MStore (HBLObjectRef p, HBLObjectRef p2, HBLObjectRef poly)
{
    long      ind1, ind2;
    if (MResolve (p,p2, ind1,ind2)) {
        MStore   (ind1,ind2,poly);
    }

}
//_____________________________________________________________________________________________
void _Matrix::MStore (long ind1, long ind2, HBLObjectRef poly) {
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
hyFloat&     _Matrix::operator [] (long i) {
    if (is_dense()) {
      return theData [i];
    }
    
    unsigned long r = (unsigned long)i / vDim,
                  c = i - vDim * r;
    
    long lIndex = Hash (r, c);
    if (lIndex == -1) {
        IncreaseStorage();
        lIndex = Hash (r, c);
    }
    if (lIndex<0) {
        theIndex[-lIndex-2] = i;
        return ((hyFloat*)theData)[-lIndex-2];
    } else {
        return ((hyFloat*)theData)[lIndex];
    }
}

//_____________________________________________________________________________________________
void        _Matrix::Store (long i, long j, hyFloat value) {
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
        ((hyFloat*)theData)[-lIndex-2] = value;
    } else {
        ((hyFloat*)theData)[lIndex] = value;
    }

}

//_____________________________________________________________________________________________
void        _Matrix::StoreObject (long i, long j, _MathObject* value, bool dup) {
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
    if (is_expression_based()) {
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
            if (((_Formula**)theData)[lIndex]!=(_Formula*)ZEROPOINTER) {
                delete ((_Formula**)theData)[lIndex];
            }
            ((_Formula**)theData)[lIndex] = copyF?(_Formula*)f.makeDynamic():&f;
            if (simplify) {
                ((_Formula**)theData)[lIndex]->SimplifyConstants();
            }
        }

        CheckIfSparseEnough();
    }
}

//_____________________________________________________________________________________________


void        _Matrix::Swap (_Matrix& m){
    Exchange(theData,m.theData);
    Exchange(hDim,m.hDim);
    Exchange(vDim,m.vDim);
    Exchange(lDim,m.lDim);
    Exchange(theIndex,m.theIndex);
    Exchange(storageType,m.storageType);
    Exchange(bufferPerRow,m.bufferPerRow);
    Exchange(overflowBuffer,m.overflowBuffer);
    Exchange(allocationBlock,m.allocationBlock);
    Exchange(theValue,m.theValue);
    Exchange(cmd,m.cmd);
}

//_____________________________________________________________________________________________

void        _Matrix::AplusBx (_Matrix& B, hyFloat x)
{
    _Matrix temp (B);
    temp *= x;
    *this+=temp;
}

//#define _SLKP_USE_SSE_INTRINSICS

//_____________________________________________________________________________________________
hyFloat        _Matrix::Sqr (hyFloat* _hprestrict_ stash) {
    
    hyFloat diff = 0.;
    if (hDim!=vDim) {
        return diff;
    }
    // not a square matrix

    if (theIndex|| storageType!=1 )
        // sparse or non-numeric matrix
    {
        _Matrix temp (hDim, vDim, storageType==0?theIndex!=nil:false, storageType);
        Multiply (temp, *this);
        Swap(temp);
        return DBL_EPSILON * 1.e4;
    } else {
        if (hDim==4)
            // special case for nucleotides
        {
            for (unsigned long i=0UL, k = 0UL; i<16; i+=4) {
                for (unsigned long j=0UL; j<4UL; j++, k++) {
                  hyFloat p1 = theData[i]   * theData [j];
                  hyFloat p2 = theData[i+1] * theData [j+4];
                   p1 += theData[i+2] * theData [j+8];
                   p2 += theData[i+3] * theData [j+12];
                  
                   stash[k] = p1+p2;
                }
            }
        } else {
            long loopBound = (vDim >> 2) << 2;
            //vDim - vDim % 4;

            // loop interchange rocks!

          
            hyFloat  * _hprestrict_ column = stash+lDim;
            hyFloat const  * source = theData;

            for (long j = 0; j < vDim; j++) {
                for (long c = 0; c < vDim; c++) {
                    column[c] = source[j + c * vDim];
                }

#ifdef _SLKP_USE_AVX_INTRINSICS
                if (vDim == 61UL) {
                  for (unsigned long i = 0; i < lDim; i += 61) {
                    hyFloat * _hprestrict_ row = theData + i;
                    
                    
                     __m256d   sum256;
                    
#ifdef _SLKP_USE_FMA3_INTRINSICS
                      
                    
                        sum256 = _mm256_fmadd_pd (_mm256_loadu_pd (row), _mm256_loadu_pd (column),
                                                 _mm256_fmadd_pd (_mm256_loadu_pd (row+4), _mm256_loadu_pd (column+4),
                                         _mm256_mul_pd (_mm256_loadu_pd (row+8), _mm256_loadu_pd (column+8))));
                      
                        for (unsigned long k = 12UL; k < 60UL; k += 12UL) {
                            sum256 =  _mm256_fmadd_pd (_mm256_loadu_pd (row+k), _mm256_loadu_pd (column+k),
                                                        _mm256_fmadd_pd (_mm256_loadu_pd (row+k+4), _mm256_loadu_pd (column+k+4),
                                                        _mm256_fmadd_pd (_mm256_loadu_pd (row+k+8), _mm256_loadu_pd (column+k+8), sum256))
                                                       );
                        }
#else
                      
                    sum256 = _mm256_setzero_pd();
                    for (unsigned long k = 0UL; k < 60UL; k += 12UL) {
                      __m256d term0 = _mm256_mul_pd (_mm256_loadu_pd (row+k), _mm256_loadu_pd (column+k));
                      __m256d term1 = _mm256_mul_pd (_mm256_loadu_pd (row+k+4), _mm256_loadu_pd (column+k+4));
                      __m256d term2 = _mm256_mul_pd (_mm256_loadu_pd (row+k+8), _mm256_loadu_pd (column+k+8));
                      
                      __m256d sum01 = _mm256_add_pd(term0,term1);
                      __m256d plus2 = _mm256_add_pd(term2, sum256);
                      
                      sum256 = _mm256_add_pd (sum01, plus2);
                    }
#endif
                    
                    stash[i+j] = _avx_sum_4(sum256) + row[60] * column [60];
                    
                  }
                  
                } else {
                  for (unsigned long i = 0; i < lDim; i += vDim) {
                      hyFloat * row = theData + i;
                      
                      __m256d   sum256 = _mm256_setzero_pd();
                    
                      long k;
                      
                      for (k = 0; k < loopBound; k += 4) {
#ifdef _SLKP_USE_FMA3_INTRINSICS
                          sum256 = _mm256_fmadd_pd (_mm256_loadu_pd (row+k), _mm256_loadu_pd (column+k), sum256);
#else
                          sum256 = _mm256_add_pd (_mm256_mul_pd (_mm256_loadu_pd (row+k), _mm256_loadu_pd (column+k)), sum256);
#endif
                      }
                    
                      hyFloat result = _avx_sum_4(sum256);
                    
                      for (; k < vDim; k++) {
                          result += row[k] * column [k];
                      }
                      
                      stash[i+j] = result;
                      
                  }
                }

#else
                for (long i = 0; i < lDim; i += vDim) {
                    hyFloat * row    = theData + i,
                                 buffer [4] = {0.,0.,0.,0.};


                    unsigned long        k;

                    for (k = 0UL; k < loopBound; k += 4UL) {
                        buffer [0] += row[k] * column [k];
                        buffer [1] += row[k+1] * column [k+1];
                        buffer [2] += row[k+2] * column [k+2];
                        buffer [3] += row[k+3] * column [k+3];
                    }

                    for (; k < vDim; k++) {
                        buffer[0] += row[k] * column [k];
                    }

                    stash[i+j] = (buffer[0] + buffer[1]) + (buffer[2] + buffer[3]);
                }
#endif
           }
        }
        
        //memcpy (theData, stash, lDim * sizeof (hyFloat));

        for (long s = 0; s < lDim; s++) {
            StoreIfGreater(diff, fabs (theData[s] - stash[s]));
            theData[s] = stash[s];
        }
    }
    return diff;
}
//_____________________________________________________________________________________________
void        _Matrix::AgreeObjects (_Matrix& m)
{
  if (storageType==2) {
    if (toPolyOrNot!=0.0) {
      ConvertFormulas2Poly ();
    } else {
      Evaluate(true);
    }
  }
  
  if (m.storageType==2) {
    if (toPolyOrNot!=0.0) {
      m.ConvertFormulas2Poly ();
    } else {
      m.Evaluate(true);
    }
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
                HBLObjectRef polyCell = ((_Formula**)theData)[i]->ConstructPolynomial();
                if (polyCell) { // valid polynomial conversion
                    tempStorage[i] = (HBLObjectRef)polyCell;
                    polyCell->AddAReference();
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
                                _Polynomial * temp = (_Polynomial *)diag.Sub(tempStorage[h], nil);
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
            HBLObjectRef polyCell = f->ConstructPolynomial();
            if (polyCell) { // valid polynomial conversion
                tempStorage[i] = (HBLObjectRef)polyCell;
                polyCell->AddAReference();
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
                        _Polynomial * temp = (_Polynomial *)diag.Sub(tempStorage[j], nil);
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
        theData = (hyFloat*) tempStorage;
        storageType = 0;
        if (!theIndex) {
            _Polynomial zero_polynomial;
            for (i=0; i<lDim; i++)
                if (!GetMatrixObject (i)) {
                    StoreObject (i,&zero_polynomial,true);
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
    theData = (hyFloat*) tempStorage;
    storageType = 0;
}


//_____________________________________________________________________________________________
void        _Matrix::operator += (_Matrix& m)
{
    AgreeObjects (m);
    if ((!m.theIndex) && theIndex) {
        CheckIfSparseEnough(true);
    }
    AddMatrix (*this,m);
}

//______________________________________________________________

long    _Matrix::CompareRows (const long row1, const long row2) {
    for (long column_id = 0; column_id < vDim; column_id ++) {
        hyFloat v1 = theData[row1*vDim+column_id],
                   v2 = theData[row2*vDim+column_id];
        if (!CheckEqual (v1,v2)) {
            return (v1 < v2)?-1L:1L;
        }
    }
    return 0L;
}

//______________________________________________________________

void    _Matrix::SwapRows (const long row1, const long row2) {
    long idx1 = row1*vDim,
         idx2 = row2*vDim;
    for (long column_id = 0; column_id < vDim; column_id ++) {
        hyFloat t = theData[idx1];
        theData[idx1++] = theData[idx2];
        theData[idx2++] = t;
    }
}
//______________________________________________________________

void    _Matrix::RecursiveIndexSort (long from, long to, _SimpleList* index) {
    long            middle          = (from+to) >> 1,
                    bottommove      = 1L,
                    topmove         = 1L;

    /*
        Use '+' to denote an element that is gretae than 'M' (the 'middle' element)
        and '-' to denote an element than is less than 'M'
     
        Initially we may have something like
     
        --++--+M--+++--++-
        and we want to end up with
        ---------M+++++++
     
        Initially, we arrange the elements as
     
        ----+++M-----++++++, and then swap 'bottommove' pluses (of which there are 3 in this case)
                            with 'topmove' minuses (of which there are 5)
    
     */
    

    if (middle)
        while (middle-bottommove>=from && CompareRows (middle-bottommove, middle) >= 0L) {
            bottommove++;
        }
    if (from<to)
        while (middle+topmove<=to && CompareRows (middle+topmove,middle) <= 0L) {
            topmove++;
        }

    for (long i=from; i<middle-bottommove; i++)
        if (CompareRows (i, middle) >= 0L) {
            SwapRows (middle-bottommove, i);
            index->Swap(middle-bottommove,i);
            bottommove++;

            while (middle-bottommove>=from && CompareRows (middle-bottommove, middle) >= 0L) {
                bottommove++;
            }
        }

    {
        for (long i=middle+topmove+1; i<=to; i++)
            if (CompareRows(i,middle) <= 0L) {
                SwapRows   (i, middle+topmove);
                index->Swap(i, middle+topmove);
                
                topmove++;
                while (middle+topmove<=to && CompareRows (middle+topmove,middle) <= 0L) {
                   topmove++;
                }
            }
    }

    if (topmove==bottommove) {
        for (long i=1; i<bottommove; i++) {
            SwapRows(middle+i, middle-i);
            index->Swap (middle+i, middle-i);
        }
    } else if (topmove>bottommove) {
        long shift = topmove-bottommove;
        // in the example above, shift = 2
        
        for (long i=1; i<bottommove; i++) {
             SwapRows (middle-i, middle+i+shift);
             index->Swap(middle-i, middle+i+shift);
        }
        // at the end of this loop, the example above will look like 
        // -------M--+++++++++, so now if we swap 'M' with the last '-', we'll arrive at the desired configuration
        
        SwapRows    (middle, middle+shift);
        index->Swap (middle, middle+shift);
        middle+=shift;
        
    } else {
        long shift = bottommove-topmove;
        for (long i=1; i<topmove; i++) {
            SwapRows (middle+i, middle-i-shift);
            index->Swap (middle+i, middle-i-shift);
        }

        SwapRows    (middle, middle-shift);
        index->Swap (middle, middle-shift);
        middle-=shift;
    }

    if (to>middle+1) {
        RecursiveIndexSort (middle+1,to, index);
    }
    if (from<middle-1) {
        RecursiveIndexSort (from,middle-1, index);
    }
}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::SortMatrixOnColumn (HBLObjectRef mp, HBLObjectRef cache)
{
    if (storageType!=1) {
        HandleApplicationError  ("Only numeric matrices can be sorted");
        return new _MathObject();
    }

    if (theData == nil) {
        return new _Matrix (0,0);
    }

    _SimpleList sortOn;

    if (mp->ObjectClass () != NUMBER || mp->Value() < 0.0 || mp->Value () > GetVDim()-1) {
        bool goodMe = false;
        if (mp->ObjectClass () == MATRIX) {
            _Matrix * sortOnM = (_Matrix*)((_Matrix*)mp)->ComputeNumeric();
            long sortBy      = sortOnM->GetHDim()*sortOnM->GetVDim(),
                 maxColumnID = GetVDim();
                  
            for (long k=0; k<sortBy; k=k+1) {
                long idx = (*sortOnM)[k];
                if (idx < 0 || idx >= maxColumnID) {
                    HandleApplicationError (_String("Invalid column index to sort on in call to ") & __PRETTY_FUNCTION__ & " : " & idx);
                    return new _MathObject();               
                }
                sortOn << idx;
            }
            goodMe = sortOn.lLength;
        }
        if (!goodMe) {
            HandleApplicationError  (_String ("Invalid column index to sort the matrix on:") & _String((_String*)mp->toStr()).Enquote());
            return new _MathObject;
        }
    } else {
        sortOn << mp->Value();
    }

    // TODO SLKP 20111109 -- replace with a generic sort function
                     // the code below is BROKEN
    
    _SimpleList             idx (hDim,0,1);
    _Matrix theColumn   (hDim,sortOn.lLength,false,true);

    for (unsigned long col2Sort = 0; col2Sort < sortOn.lLength; col2Sort++) {
        long colIdx = sortOn.list_data[col2Sort];

        if (theIndex)
            for (long k=0; k<hDim; k++) {
                theColumn.theData[col2Sort+k*sortOn.lLength] = (*this)(k, colIdx);
            }
        else
            for (long k=0, j = colIdx; k<hDim; k++, j+=vDim) {
                theColumn.theData[col2Sort+k*sortOn.lLength] = theData[j];
            }

    }

    theColumn.RecursiveIndexSort (0, hDim-1, &idx);
    _Matrix                 *result     = (_Matrix*)_returnMatrixOrUseCache(hDim, vDim, _NUMERICAL_TYPE, theIndex, cache);

    if (theIndex) {
        _SimpleList    revIdx (hDim,0,1);
        SortLists (&idx, &revIdx);
        for (long r=0; r<lDim; r++) {
            long oi = theIndex[r];

            if (oi >= 0) {
                long     v  = oi%vDim,
                         h  = oi/vDim,
                         ni = revIdx.list_data[h]*vDim+v;

                (*result)[ni] = theData[r];
            }
        }
    } else
        for (long r=0; r<hDim; r++) {
            long remapped = idx.list_data[r];
            remapped *= vDim;
            for (long c=r*vDim; c<r*vDim+vDim; c++, remapped++) {
                result->theData[c] = theData[remapped];
            }
        }


    return result;
}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::PoissonLL (HBLObjectRef mp, HBLObjectRef cache)
{
    if (!is_numeric()) {
        HandleApplicationError ("Only numeric matrices can be passed to Poisson Log-Likelihood");
        return new _MathObject;
    }

    if (mp->ObjectClass () != NUMBER || mp->Value() < 0.0) {
        HandleApplicationError  (_String ("Invalid Poisson distribution parameter") & (_String((_String*)mp->toStr())).Enquote());
        return new _MathObject;
    }

    hyFloat     loglik = 0.0,
                   *logFactorials = new hyFloat [101],
    lambda        = mp->Value(),
    logLambda     = log (lambda),
    log2p         = log (sqrt(8.*atan(1.)));


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
                        logFactorials[idx2] = logFactorials[idx2-1]+log((hyFloat)idx2);
                    }
                    loglik += logLambda * cellValue - lambda - logFactorials [cellValue];
                    maxFactorialDone = cellValue;
                } else
                    // use Stirling's formula
                {
                    loglik += logLambda * cellValue - lambda + cellValue - (cellValue+0.5)*log((hyFloat)cellValue)-log2p;
                }
            }
        }
    }

    delete      [] logFactorials;

    return _returnConstantOrUseCache(loglik, cache);
}


//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::PathLogLikelihood (HBLObjectRef mp, HBLObjectRef cache) {
    try {
        _Matrix                 *m          = nil;

        if (! is_numeric() || hDim != 3) {
            throw (_String("First argument must be a numeric 3xN matrix"));
        } else {
            //errMsg = "Second argument in call to < (PathLogLikelihood) must be a square matrix";
            if (mp->ObjectClass () == MATRIX) {
                m = (_Matrix*)mp->Compute();
                if (m->GetHDim() != m->GetVDim()) {
                    throw (_String("Second argument must be a square matrix"));
                }
            } else {
                throw (_String("Second argument must be a matrix"));
            }
        }


        CheckIfSparseEnough     (true);

        hyFloat              res     = 0.0;
        long                    maxDim  = m->GetHDim();

        for (unsigned long step = 0UL; step < vDim; step++) {
            
            long        i1 = get (0,step),
                        i2 = get (1,step);
            hyFloat     t  = get (2,step);

            if (i1<0 || i2 < 0 || i1 >= maxDim || i2 >= maxDim || t<0.0) {
                throw (_String ("An invalid transition in step ") & _String ((long)(step+1L)) & " of the chain: " & i1 & " to " & i2 & " in time " & t);
            }

            _Matrix         rateMx (*m);
            rateMx *= t;
            _Matrix   * tMatrix = rateMx.Exponentiate ();
            t = tMatrix->theData[maxDim*i1+i2];
            DeleteObject (tMatrix);

            if (t>0.0) {
                res += log (t);
            } else {
                return _returnConstantOrUseCache(-1.e300, cache);
            }
        }
        return _returnConstantOrUseCache(res, cache);
    } catch (const _String& err) {
        HandleApplicationError  (err);
        return new _MathObject;
    }

}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::pFDR (HBLObjectRef classes, HBLObjectRef cache) {
    try {
        long            steps     = 20,
                        iter_count = 500;

        hyFloat         p_value = 0.0,
                        max_lambda = 0.0;


        if (theIndex) {
            CheckIfSparseEnough (true);
        }

        if (!is_numeric()) {
            throw _String("Only numeric matrices can be passed to && (pFDR)");
        } else {
            if (!(is_column() || is_row()) || is_empty())   {
                throw _String("The first argument of && (pFDR) must be an Nx1/1xN matrix.");
            } else if (classes->ObjectClass () != NUMBER || classes->Value() > 1. || (p_value = classes->Value()) < 0.0) {
                throw _String ("Invalid baseline p-value (must be in (0,1)):") & _String((_String*)classes->toStr());
            } else {
                for (unsigned long i=0UL; i<lDim; i++) {
                    hyFloat p_count = theData[i];
                    if (p_count < 0.0 || p_count > 1.0) {
                        throw _String ("Invalid p-value entry in matrix passed to pFDR (must be a positive integer):") & p_count;
                    }
                    StoreIfGreater(max_lambda, p_count);
                }
            }
        }


        _Matrix        lamdbaRange (steps,1,false,true),
                       pFDRs       (steps,1,false,true);

        hyFloat     anLamdba           = 0.0,
                       minPFDR          = 5.0,
                       uberPFDR        = 0.0,
                       uberPFDRUpperLimit = 0.0,
                       minMSE             = 1.e100,
                       aStep            = 1.0/steps;


        unsigned long k = 0;
        while (anLamdba<1.0) {
            lamdbaRange.theData[k] = anLamdba;

            if ((pFDRs.theData[k] = computePFDR (anLamdba, p_value))<minPFDR) {
                minPFDR = pFDRs.theData[k];
            }

            k++;
            anLamdba += aStep;
        }

        for (unsigned long k=0UL; k<steps; k++) {
            hyFloat mse    = 0.0;
            _Matrix    ITpDFR (iter_count,1,false,true);

            for (unsigned long it = 0; it < iter_count; it++) {
                _Matrix         sampledPs (lDim,1,false,true);
                _SimpleList     sample    (lDim,0,1);
                sample.PermuteWithReplacement (1);

                for (long el = 0; el < lDim; el++) {
                    sampledPs.theData[el] = theData[sample.list_data[el]];
                }

                ITpDFR.theData[it] = sampledPs.computePFDR (lamdbaRange.theData[k], p_value);
                mse += (ITpDFR.theData[it]-minPFDR)*(ITpDFR.theData[it]-minPFDR);
            }

            mse /= iter_count;

            if (mse < minMSE) {
                minMSE = mse;
                uberPFDR = pFDRs.theData[k];
                _Constant  zer (0.0);
                _Matrix* sorted = (_Matrix*)ITpDFR.SortMatrixOnColumn (&zer, nil);
                uberPFDRUpperLimit = sorted->theData[((long)(0.95*iter_count))];
                DeleteObject (sorted);
            }
        }

        _Matrix * resMx = (_Matrix *) _returnMatrixOrUseCache (2,1,_NUMERICAL_TYPE,false, cache);
        resMx->theData[0] = uberPFDR;
        resMx->theData[1] = uberPFDRUpperLimit;
        return resMx;
    } catch (const _String& err) {
        HandleApplicationError  (err);
        return new _MathObject;
    }
}

//_____________________________________________________________________________________________
hyFloat      _Matrix::computePFDR (hyFloat lambda, hyFloat gamma)
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
        hyFloat pi_0 = null/(lDim*(1.-lambda)),
                   pr_p = 0;

        if (rejected) {
            pr_p = rejected/(hyFloat)lDim;
        } else {
            pr_p = 1./(hyFloat)lDim;
        }

        return pi_0 * gamma / (pr_p /** (1.-exp(log(1.-gamma)*lDim))*/);

    } else {
        return 1;
    }
}

//_____________________________________________________________________________________________

HBLObjectRef _Matrix::Random (HBLObjectRef kind, HBLObjectRef cache) {

    try {
        long columns = GetVDim(),
             rows    = GetHDim();

        if (kind->ObjectClass() == NUMBER) {
            bool    resample = (kind->Compute()->Value()>0);

            _SimpleList     remapped (columns,0,1);

            if (resample) {
                remapped.PermuteWithReplacement(1);
            } else {
                remapped.Permute(1);
            }


            if (is_numeric()) {   // numeric matrix
                _Matrix * res = (_Matrix *)_returnMatrixOrUseCache(rows, columns,_NUMERICAL_TYPE, theIndex != nil,cache);

                if (is_dense())
                    for (unsigned long vv = 0; vv<lDim; vv+=columns)
                        for (unsigned long k2=0; k2<remapped.lLength; k2++) {
                            res->theData[vv+k2] = theData[vv+remapped.list_data[k2]];
                        }
                else {
                    for (unsigned long vv = 0; vv< rows; vv++)
                        for (unsigned long k=0; k<remapped.lLength; k++) {
                            long ki = remapped.list_data[k];
                            if ((ki = Hash (vv,ki)) >= 0) {
                                res->Store (vv,k,theData[ki]);
                            }
                        }
                }
                return res;
            } else {            // formula matrix
                if (is_expression_based()) {
                    _Matrix * res = new _Matrix (rows, columns,theIndex != nil,false);

                    for (unsigned long vv = 0UL; vv< rows; vv++)
                        for (unsigned long k=0UL; k<remapped.lLength; k++) {
                            _Formula * ff = GetFormula (vv,remapped.get (k));
                            if (ff) {
                                res->StoreFormula (vv, k, *ff);
                            }
                        }
                    return res;
                }
            }
        }

        else if (kind->ObjectClass() == ASSOCIATIVE_LIST) {
            //ReportWarning (_String("_Matrix::Random() with associative list as first argument."));

            // Associative list should contain following arguments:
            //  "PDF" - string corresponding to p.d.f. ("Gamma", "Normal")
            //  "ARG0" ... "ARGn" - whatever parameter arguments (matrices) are required for the p.d.f.
            
            _AssociativeList    * pdfArgs   = (_AssociativeList *)kind;
            _List               * keys      = pdfArgs->GetKeys();
            _String             pdfkey      ("PDF"),
                                * arg0      = (_String *)pdfArgs->GetByKey(pdfkey,STRING);
            DeleteObject (keys);
            
            
            if (arg0) {
                _String     pdf ((_String*)arg0->toStr()),
                            arg ("ARG0");
                
                long        pdfCode = _HY_MatrixRandomValidPDFs.GetValueFromString (pdf);
                
                 switch (pdfCode) {
                    case _HY_MATRIX_RANDOM_DIRICHLET:
                        return (_Matrix *) DirichletDeviate();
                    case _HY_MATRIX_RANDOM_GAUSSIAN:
                        return (_Matrix *) GaussianDeviate (*(_Matrix *) pdfArgs->GetByKey (arg, MATRIX));
                    case _HY_MATRIX_RANDOM_WISHART:
                        return (_Matrix *) WishartDeviate (*(_Matrix *) pdfArgs->GetByKey (arg, MATRIX));
                    case _HY_MATRIX_RANDOM_INVERSE_WISHART:
                        return (_Matrix *) InverseWishartDeviate (*(_Matrix *) pdfArgs->GetByKey (arg, MATRIX));
                    case _HY_MATRIX_RANDOM_MULTINOMIAL:
                        return (_Matrix *) MultinomialSample ((_Constant *) pdfArgs->GetByKey (arg, NUMBER));
                    default:
                        throw _String("String argument passed to Random not a supported PDF: ") & pdf.Enquote();
                }
            } else {
                throw _String("Expecting 'PDF' key in associative list argument passed to Random(), received: ") & *arg0;
            }

        } else if (kind->ObjectClass () == STRING) {
            _String key = ((_FString*)kind->Compute())->get_str();
            if (key == _String("LHS")) {
                // latin hypercube sampling: samples are in ROWS
                _Matrix * lhc = new _Matrix ( rows, columns, false, true);

                _SimpleList permutation ( rows ,0,1);

                for (unsigned long c = 0; c < columns; c++) {
                    permutation.Permute (1);
                    for (long r = 0; r < rows ; r++) {
                        lhc->set(r,c) = get (permutation.get(r),c);
                    }
                }

                return lhc;
            }
            throw _String ("Invalid string argument passed to matrix Random :") & key;
        } else {
            throw _String ("Invalid argument passes to matrix Random (should be a number, an associative list or a string):") & _String((_String*)kind->toStr());
        }
    } catch (_String const& err) {
        HandleApplicationError (err);
    }
    return new _Matrix (1,1);
}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::K_Means (HBLObjectRef classes, HBLObjectRef cache) {
    // K-means clustering on scalar data
    /*
     
     this    : Nx2 matrix {{value1, count1}{value2, count 2}...}}
     classes : 2x1 matrix {{cluster count}{number of random restarts}}
     
     reutn   : 2 x Max (cluster count, 2)  {{cluster mean 1, cluster mean 2, ... , cluster mean N}{total L^2 error, how many restarts hit the min}}
     
     */
    
    // TODO: 20171026 SLKP revised. check correctness
    
     try {
        _Matrix     *   arg;
        long            cluster_count,
                        iter_count,
                        sample_count = 0L;

        if (theIndex) {
            CheckIfSparseEnough (true);
        }

        if (!is_numeric()) {
            throw _String("Only numeric matrices can be passed to <= (K-means)");
        } else {
            if (GetVDim () != 2) {
                throw _String("The first argument of <= (K-means) must be an Nx2 matrix, with samples in the first column, and counts in the 2nd.");
            } else if (classes->ObjectClass () != MATRIX) {
                throw _String ("Invalid number of clusters is call to K-means (must be >=1):") & _String((_String*)classes->toStr());
            } else {
                arg = (_Matrix*)classes->Compute();
                if (!arg->check_dimension (1,2) || (cluster_count=arg->theData[0]) < 1 || (iter_count = arg->theData[1]) < 1) {
                    throw _String ("Invalid second argument is call to K-means (must be a 2x1 matrix of positive integers specifying cluster_count and maximum number of random restarts. Had ") & _String((_String*)classes->toStr()).Enquote();
                } else {
                    for (unsigned long i=1UL; i<lDim; i+=2UL) {
                        long pCount = theData[i];
                        if (pCount <= 0L) {
                             throw _String ("Invalid count entry in matrix passed to K-means (must be a positive integer):") & pCount;
                        }
                        sample_count += pCount;
                    }
                }
            }
        }
        
        if (cluster_count > sample_count) {
            throw _String ("More clusters requested than available data points");
        }

        _Matrix * res = (_Matrix*)_returnMatrixOrUseCache(2, cluster_count, _NUMERICAL_TYPE, false, cache);

        if (cluster_count == 1L) {
            hyFloat sampleMean    = 0.,
                    errorEstimate = 0.;

            for (unsigned long c1=0UL, c2=1UL; c1 < 2*hDim; c1+=2UL, c2+=2UL) {
                sampleMean +=  theData[c1] * theData[c2];
            }

            sampleMean /= sample_count;

            for (unsigned long c1=0UL, c2=1UL; c1 < 2*hDim; c1+=2UL, c2+=2UL) {
                hyFloat locErr = theData[c1] - sampleMean;
                errorEstimate += locErr*locErr*theData[c2];
            }

            res->theData[0] = sampleMean;
            res->theData[1] = errorEstimate;
        } else {

            hyFloat  minError    = 1.e100;
            

            _SimpleList full_copy_list    ((unsigned long)sample_count,0,0);

            for (unsigned long c1 = 0UL, overall = 0UL; c1 < hDim; c1++) {
                unsigned long copies = get (c1, 1);

                for (unsigned long c2 = 0; c2 < copies; c2++, overall++) {
                    full_copy_list.list_data[overall] = c1;
                }
            }


            long        hit_min_error;
            _Matrix     best_cluster_means;

            for (unsigned long sampleCount = 0UL; sampleCount < iter_count; sampleCount ++) {
                // choose N random cluster centers to start with
                _SimpleList       chosen_means        = full_copy_list.Sample(cluster_count),
                                  cluster_assignments (hDim, 0, 0);
                
                _Matrix           cluster_means        (cluster_count,2,false,true);
                
                

                for (long cc = 0; cc < cluster_count; cc = cc+1) {
                    cluster_means.set (cc, 0) = get(chosen_means.get(cc), 0);
                }

                hyFloat            last_error_estimate = 1.e100,
                                   error_estimate      = 0.;
                
                
                for (unsigned long cIters = 0UL; cIters < hDim * 5; cIters ++) {
                    
                        bool moved = false;
                        // assign each point to the nearest cluster centroid
                        for (unsigned long data_point = 0UL; data_point < hDim; data_point++) {
                            
                            unsigned long best_cluster = cluster_assignments.get (data_point);
                            
                            hyFloat         this_value              = get  (data_point, 0),
                                            current_min_distance    = fabs (this_value-cluster_means.get (cluster_assignments.get (data_point),0));
                            
                            for (unsigned long cluster_id = 0UL; cluster_id < cluster_count; cluster_id++) {
                                if (StoreIfLess (current_min_distance, fabs (this_value-cluster_means.get (cluster_id,0)))) {
                                    best_cluster = cluster_id;
                                }
                            }
                                    
                            if (best_cluster != cluster_assignments.get (data_point)) {
                                moved = true;
                                cluster_assignments[data_point] = best_cluster;
                            }
                        }
                    
                        if (moved) {
                            for (long cc = 0; cc < cluster_count; cc = cc+1) {
                                cluster_means.set (cc,0) = 0.;
                                cluster_means.set (cc,1) = 0.;
                            }
                            for (unsigned long data_point = 0UL; data_point < hDim; data_point++) {
                                cluster_means.set (cluster_assignments.get (data_point), 0) += get (data_point,0);
                                cluster_means.set (cluster_assignments.get (data_point), 1) += get (data_point,1);
                           }
                           for (long cc = 0; cc < cluster_count; cc = cc+1) {
                               if (cluster_means.get (cc, 1) != .0) {
                                   cluster_means.set (cc,0) /= cluster_means.get (cc,1);
                               }
                           }
                        } else {
                            break;
                        }
                        
                        
                }

                if (minError == 0.0 || fabs((minError-error_estimate)/minError) < 0.001) {
                    hit_min_error ++;
                } else if (error_estimate < minError) {
                    hit_min_error = 1;
                    minError    = error_estimate;
                    best_cluster_means = cluster_means;
                }
            }

              for (long k2 = 0; k2 < cluster_count; k2++) {
                res->theData[k2] = best_cluster_means.get (k2, 0);
            }
            
            res->theData[cluster_count]   = minError;
            res->theData[cluster_count+1] = hit_min_error;
            return res;
        }
    } catch (_String const& err) {
        HandleApplicationError (err);
    }

   return new _Matrix;
}


//_____________________________________________________________________________________________
void            _Matrix::PopulateConstantMatrix (hyFloat v) {
    if (is_numeric()) {
        InitializeArray(theData, lDim, (hyFloat&&)v);
    }
}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::AddObj (HBLObjectRef mp, HBLObjectRef cache)
{
    if (_Matrix::ObjectClass()!=mp->ObjectClass()) {
        if (mp->ObjectClass () == STRING) {
            _FormulaParsingContext def;
            _Matrix * convMatrix = new _Matrix (((_FString*)mp)->get_str(), false, def),
            * res;
            res = (_Matrix*)AddObj (convMatrix, cache);
            DeleteObject (convMatrix);
            return res;
        }
        if (mp->ObjectClass () == NUMBER) {
            _Matrix* aNum = (_Matrix*)ComputeNumeric ();
            
            hyFloat pValue = mp->Value();
            
            if (aNum->is_numeric()) {
                return ApplyScalarOperation ([=] (hyFloat h) -> hyFloat {return h + pValue;}, cache);
            }
        }

        HandleApplicationError ( kErrorStringIncompatibleOperands );
        return new _Matrix (1,1);
    }

    _Matrix * m = (_Matrix*)mp;
    AgreeObjects (*m);
    _Matrix * result = (_Matrix *)_returnMatrixOrUseCache( hDim, vDim, storageType, theIndex && m->theIndex, cache);
    AddMatrix (*result,*m);
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
void       _Matrix::NonZeroEntries (_SimpleList& target) {
    if (theIndex && storageType == 1) {
        target.Clear();
        target.RequestSpace(lDim);
        for (long elementID = 0; elementID < lDim; elementID ++) {
            if (theIndex[elementID] >= 0) {
                target << theIndex[elementID];
            }
        }
        target.Sort();
    }
}

//_____________________________________________________________________________________________
bool       _Matrix::Equal(HBLObjectRef mp)
{
    if (mp->ObjectClass()!=ObjectClass()) {
        return false;
    }

    _Matrix * m = (_Matrix*)mp;
    
    if (m->storageType == storageType && m->hDim == hDim && m->vDim == vDim) {
        if (is_numeric()) {
            if (theIndex || m->theIndex) {
                for (long r = 0L; r < hDim; r ++) {
                    for (long c = 0L; c < vDim; c++) {
                        if (!CheckEqual((*this)(r,c), (*m)(r,c))) {
                            return false;
                        }
                    }
                }
                
            } else {
                for (long elementID = 0; elementID < lDim; elementID ++) {
                    if (!CheckEqual(theData[elementID], m->theData[elementID])) {
                        return false;
                    }
                }
            }
            
            return true;
        } else {
            if (IsAStringMatrix() && m->IsAStringMatrix()) {
                for (long r = 0L; r < hDim; r ++) {
                    for (long c = 0L; c < vDim; c++) {
                        _Formula * f1 = GetFormula(r,c),
                                 * f2 = GetFormula(r,c);
                        
                        if (f1 && f2) {
                            if (((_FString*)f1->Compute())->get_str() != ((_FString*)f2->Compute())->get_str()) {
                                return false;
                            }
                        } else {
                            if (f1 || f2) {
                                return false;
                            }
                        }
                    }
                }
                return true;
            }
        }
    }
    
    return false;
}


//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::SubObj (HBLObjectRef mp, HBLObjectRef cache)
{
    if (mp->ObjectClass()!=ObjectClass()) {
        HandleApplicationError ( kErrorStringIncompatibleOperands );
        return new _Matrix (1,1);
    }

    _Matrix * m = (_Matrix*)mp;
    AgreeObjects (*m);
    _Matrix * result = (_Matrix*) _returnMatrixOrUseCache(hDim, vDim, storageType,theIndex && m->theIndex, cache);
    Subtract (*result,*m);
    return result;
}

//_____________________________________________________________________________________________
void        _Matrix::operator *= (hyFloat c) {
    Multiply (*this,c);
}

//_____________________________________________________________________________________________
_Matrix     _Matrix::operator * (hyFloat c) {
    _Matrix result (*this);
    Multiply (result,c);
    return result;
}

//_____________________________________________________________________________________________
void        _Matrix::operator *= (_Matrix& m) {
    if (CheckDimensions     (m)) {
        AgreeObjects        (m);
        _Matrix   result    (hDim, m.vDim, false, storageType);
        Multiply            (result,m);
        //if ((theIndex!=nil)||(m.theIndex!=nil)) result.AmISparse();
        if (theIndex!=nil && m.theIndex!=nil) {
            result.AmISparse();
        }
        Swap                (result);
    }
}

//_____________________________________________________________________________________________
void        _Matrix::MultbyS (_Matrix& m, bool leftMultiply, _Matrix* externalStorage, hyFloat* stash) {
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
        memset (externalStorage->theData, 0, sizeof (hyFloat)*externalStorage->lDim);
        //for (long s = 0; s < externalStorage->lDim; s++) externalStorage->theData[s] = 0.0;
    }
}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::MultObj (HBLObjectRef mp, HBLObjectRef cache) {
  
  if (mp->ObjectClass()!=ObjectClass()) {
    if (mp->ObjectClass()!=NUMBER) {
      HandleApplicationError ( kErrorStringIncompatibleOperands );
      return new _Matrix (1,1);
    } else {
      hyFloat theV = mp->Value();
      return (HBLObjectRef)((*this)*theV).makeDynamic();
    }
  }
  
  _Matrix*        m = (_Matrix*)mp;
  if (!CheckDimensions (*m)) return new _MathObject;
  AgreeObjects    (*m);
  
  _Matrix*      result = (_Matrix*) _returnMatrixOrUseCache(hDim, m->vDim, storageType, false, cache);
  Multiply      (*result,*m);
  return        result;
  
}

//_____________________________________________________________________________________________
HBLObjectRef       _Matrix::MultElements (HBLObjectRef mp, bool elementWiseDivide, HBLObjectRef cache) {
    
    if (mp->ObjectClass()!=ObjectClass()) {
        HandleApplicationError ( kErrorStringIncompatibleOperands );
        return new _Matrix (1,1);
    }
    
    _Matrix* m = (_Matrix*)mp;

    bool by_column = false;
    // if the second argument has dimension 1xcolumns of the first matrix, then
    // result [i][j] is assigned this [i][j] * / argument [0][j]
    // in other words, divide or multiply each column
    
    bool by_row    = false;
    // if the first argument has dimension rows of the second matrix x 1 then
    // result [i][j] is assigned argument [i][j] * / this [i][0]
    // in other words, divide or multiply each row
    
    
    if ( GetHDim()!=m->GetHDim()  || GetVDim()!=m->GetVDim()) {
        if (GetVDim() == m->GetVDim() && m->GetHDim () == 1) {
            by_column = true;
        } else {
            if (GetHDim() == m->GetHDim() && GetVDim () == 1) {
                by_row = true;
            } else {
                HandleApplicationError ("Element-wise multiplication/division requires matrixes of the same dimension, or (NxM) $ (1xM) or (Nx1) $ (NxM) matrices ");
                return new _Matrix (1,1);
            }
        }
    }
    
    if (! is_numeric() || ! m->is_numeric() ) {
        HandleApplicationError ("Element-wise multiplication/division only works on numeric matrices");
        return new _Matrix (1,1);
    }
    
    _Matrix*      result = (_Matrix*) _returnMatrixOrUseCache(GetHDim(), m->GetVDim(), _NUMERICAL_TYPE, false, cache);
    
    if (theIndex || m->theIndex) {
        auto operation = elementWiseDivide ? DivNumbers : MultNumbers;
        
        long index = 0L;
        if (by_row) {
            for (long row = 0; row < hDim; row++) {
                for (long column = 0; column < m->vDim; column++, index++) {
                    result->theData[index] = operation ( (*this)(row,0), (*m)(row,column));
                }
            }
        } else {
            if (by_column) {
                for (long row = 0; row < hDim; row++) {
                    for (long column = 0; column < m->vDim; column++, index++) {
                        result->theData[index] = operation ( (*this)(row,column), (*m)(0,column));
                    }
                }
            }
            else {
                for (long row = 0; row < hDim; row++) {
                    for (long column = 0; column < m->vDim; column++, index++) {
                        result->theData[index] = operation ( (*this)(row,column), (*m)(row,column));
                    }
                }
            }
        }
    } else {
        if (elementWiseDivide) {
            long index = 0L;
            if (by_row) {
                for (long row = 0; row < hDim; row++) {
                    for (long column = 0; column < m->vDim; column++, index++) {
                        result->theData[index] = theData[row] / m->theData [index];
                    }
                }
            } else {
                if (by_column) {
                    for (long row = 0; row < hDim; row++) {
                        for (long column = 0; column < m->vDim; column++, index++) {
                            result->theData[index] = theData[index] / m->theData [column];
                        }
                    }
                }
                else {
                    for (long row = 0; row < hDim; row++) {
                        for (long column = 0; column < m->vDim; column++, index++) {
                            result->theData[index] = theData[index] / m->theData [index];
                        }
                    }
                }
            }
        } else {
            long index = 0L;
            if (by_row) {
                for (long row = 0; row < hDim; row++) {
                    for (long column = 0; column < m->vDim; column++, index++) {
                        result->theData[index] = theData[row] * m->theData [index];
                    }
                }
            } else {
                if (by_column) {
                    for (long row = 0; row < hDim; row++) {
                        for (long column = 0; column < m->vDim; column++, index++) {
                            result->theData[index] = theData[index] * m->theData [column];
                        }
                    }
                }
                else {
                    for (long row = 0; row < hDim; row++) {
                        for (long column = 0; column < m->vDim; column++, index++) {
                            result->theData[index] =theData[index] * m->theData [index];
                        }
                    }
                }
            }
        }
    }
    
    if (theIndex||m->theIndex) {
        result->AmISparse();
    }
    
    return  result;
}

//_____________________________________________________________________________________________
bool    _Matrix::CheckDimensions (_Matrix& secondArg) const {
// check matrix dimensions to ensure that they are multipliable
    if (vDim!=secondArg.hDim) {
        if (hDim == 1 && secondArg.hDim==1 && vDim == secondArg.vDim) { // handle scalar product separately
            secondArg.Transpose();
        } else {
            char str[255];
            snprintf (str, sizeof(str),"Incompatible matrix dimensions in call to CheckDimension: %ldx%ld and %ldx%ld\n",hDim,vDim,secondArg.hDim,secondArg.vDim);
            HandleApplicationError (str);
            return false;
        }
    }
    return true;
}

//_____________________________________________________________________________________________
_Matrix     _Matrix::operator * (_Matrix& m)
{
    if (!CheckDimensions (m)) {   
        _Matrix d;
        return d;
    }
    
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
    AddMatrix (result,m);
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

//_________________________________________________________
void    _Matrix::internal_to_str (_StringBuffer* string, FILE * file, unsigned long padding) {
    
    StringFileWrapper res (string, file);
   _String padder (" ", padding);
    
    static const _String kUseJSONForMatrix ("USE_JSON_FOR_MATRIX");
    
    bool is_numeric_mx = is_numeric ();
    bool directly_printable  = is_numeric_mx || IsAStringMatrix ();
    
    if (directly_printable) {
        
        long digs         = -1L;
        
        bool doJSON = hy_env::EnvVariableTrue(kUseJSONForMatrix);
        
        char openBracket  = doJSON ? '[' : '{',
             closeBracket = doJSON ? ']' : '}';
        
        if (is_numeric_mx) {
            digs = MIN (print_digit_specification = hy_env::EnvVariableGetDefaultNumber(hy_env::print_float_digits), 15);
        }
        res << padder << openBracket << kStringFileWrapperNewLine;

        if (is_numeric_mx) {
            
            _String formatStr = _String("%") &_String(digs+6)&'.'&_String(digs)&'g';
            
             char  number_buffer [256];
 
            for (long i = 0L; i<hDim; i++) {
                if (i) {
                    res << padder;
                }
                res << openBracket;

                for (long j = 0L; j<vDim; j++) {
                    if (j) {
                        res << ", ";
                    }
                    parameterToCharBuffer ((*this)(i,j), number_buffer, 255, doJSON);
                    res << number_buffer;
                }
                res << closeBracket << (doJSON && i != hDim -1 ? ',' : ' ') << kStringFileWrapperNewLine;
             }
        } else {
            for (long i = 0L; i<hDim; i++) {
                if (i) {
                    res << padder;
                }
                res << openBracket;

                for (long j = 0L; j<vDim; j++) {
                    if (j) {
                        res << ", ";
                    }
                    res << '"';
                    
                    _Formula * f = GetFormula (i,j);
                    if (f) {
                        HBLObjectRef fv = f->Compute();
                        if (fv) {
                          res << _String ((_String*)fv->toStr());
                          //;((_FString*)fv)->get_str();
                        }
                    }
                    res << '"';

                }
                res << closeBracket << (doJSON && i != hDim -1  ? ',' : ' ') << kStringFileWrapperNewLine;
            }
        }
        res << padder << closeBracket;
    } else if (storageType==_POLYNOMIAL_TYPE) {
        ANALYTIC_COMPUTATION_FLAG  = hy_env::EnvVariableTrue (ANAL_COMP_FLAG);
        if (!ANALYTIC_COMPUTATION_FLAG) {
            ((_Matrix*)Compute())->internal_to_str (string, file, padding);
            return;
        }
        for (long i = 0; i<hDim; i++) {
            res << "\n[";
            for (long j = 0; j<vDim; j++) {
                long p = Hash (i,j);
                if (j) {
                    res << ",";
                }
                if (p>=0) {
                    _String *sp = (_String*) GetMatrixObject (p)->toStr();
                     res << *sp;
                     DeleteObject (sp);
                } else {
                    res << "0.";
                }
            }
            res << "]\n";
        }
    } else {
        _Matrix* eval = (_Matrix*)(storageType==3?EvaluateSimple():Evaluate(false));
        eval->internal_to_str(string, file, padding);
        DeleteObject (eval);
    }
}
//_________________________________________________________
void    _Matrix::toFileStr (FILE*dest, unsigned long padding){
    internal_to_str(nil, dest, padding);
}
//_____________________________________________________________________________________________

BaseRef _Matrix::toStr(unsigned long padding) {
    _StringBuffer * serialized = new _StringBuffer (2048L);
    internal_to_str (serialized, nil, padding);
    return serialized;
}

//_____________________________________________________________________________________________

void     _Matrix::Serialize (_StringBuffer& res, _String& myID) {
    if (storageType != _POLYNOMIAL_TYPE) {
        res << '\n';
        res <<  myID;
        if (is_numeric()) {
            res << '=';
            res.AppendNewInstance((_String*)toStr());
            res << ';';
        } else if (is_expression_based()) {
            res << (_String ("={") & hDim & ',' & vDim & "};\n");
            for (long h=0L; h<hDim; h++) {
                for (long v=0L; v<vDim; v++) {
                    _Formula *theCell = GetFormula (h,v);
                    if (theCell&& !theCell->IsEmpty()) {
                        res << myID << '[' << _String(h) << "][" << _String(v) << "]:=";
                        res.AppendNewInstance((_String*)theCell->toStr(kFormulaStringConversionNormal));
                        res << ";\n";
                    }
                }
            }
        }
    }
}


//_____________________________________________________________________________________________

void    SetIncrement (int m) {
    _Matrix::storageIncrement = m;
}
//_____________________________________________________________________________________________
void    _Matrix::InitMxVar (_SimpleList& mxVariables, hyFloat glValue) {
    mxVariables.Each ([&] (long value, unsigned long) -> void {
        LocateVar(value)->SetValue (new _Constant (glValue), false,true, NULL);
    });
}
//_____________________________________________________________________________________________
bool    _Matrix::ImportMatrixExp (FILE* theSource) {
    // TODO: SLKP 20171027, need to review and possibly deprecate
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

            _Variable * ppv = CheckReceptacle (&varName, kEmptyString, true);
            varList << ppv->get_index();
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
        while (fc!='{') {
            fc = fgetc (theSource);
            buffer[i] = fc;
            i++;
            if (feof(theSource)) {
                return false;
            }
        }
        _Polynomial* thisCell = new _Polynomial (varList);
        m = atol (buffer);
        hyFloat* theCoeffs = (hyFloat*)MatrixMemAllocate(m*sizeof(hyFloat));
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
        fc = fgetc(theSource);
        if (fc != '{') {
            return false;
        }
        _PolynomialData *pd = new _PolynomialData (varList.countitems(),j,theCoeffs);
        MatrixMemFree (theCoeffs);
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
// TODO: SLKP 20171027, need to review and possibly deprecate

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
        HandleApplicationError ( kErrorStringMatrixExportError );
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
    hyFloat* varPool = (hyFloat*)MatrixMemAllocate (mxVariables.countitems()*sizeof(hyFloat));
    for (k=0; k<mxVariables.countitems(); k++) {
        fprintf (theDump,"%s",LocateVar(mxVariables(k))->GetName()->get_str());
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
        hyFloat* coeffHolder =  (hyFloat*)MatrixMemAllocate (nTerms*sizeof(hyFloat)), error, bestError = 1;

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
                if (termRank.list_data[i]<=tup) {
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
            snprintf (be, sizeof(be),"%g",bestError);
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
}

//_____________________________________________________________________________________________

hyFloat  _Matrix::ExpNumberOfSubs  (_Matrix* freqs, bool mbf) {
    // TODO SLKP 20171027 SLKP reviewed and edited; check correctness
    
    if (!is_square_numeric(false) || !freqs->is_numeric()) {
        return 0.0;
    }

    hyFloat      result      =   0.0;
    _Matrix      *stencil    =   BranchLengthStencil();

    if ( freqs->is_dense() == false ) {
        freqs->CheckIfSparseEnough(true);
    }
    
    if (stencil) {
        if (mbf) {
            ForEachCellNumeric ([&] (hyFloat value, unsigned long index, unsigned long row, unsigned long column) -> void {
                if (row != column && stencil->theData[index]) {
                    result += value * freqs->theData[row] * freqs->theData[column];
                }
            });
        } else {
            ForEachCellNumeric ([&] (hyFloat value, unsigned long index, unsigned long row, unsigned long column) -> void {
                if (row != column && stencil->theData[index]) {
                    result += value * freqs->theData[row];
                }
            });
        }
    } else {
        if (mbf) {
            ForEachCellNumeric ([&] (hyFloat value, unsigned long index, unsigned long row, unsigned long column) -> void {
                if (row != column) {
                    result += value * freqs->theData[row] * freqs->theData[column];
                }
            });
        } else {
            ForEachCellNumeric ([&] (hyFloat value, unsigned long index, unsigned long row, unsigned long column) -> void {
                if (row != column) {
                    result += value * freqs->theData[row];
                }
            });
        }
    }
    return result;
}

//_____________________________________________________________________________________________
_List*      _Matrix::ComputeRowAndColSums (void) {
// the first entry is the matrix with row sums
// the second - the entry with column sums
// the third  - a constant with the total sum
    if ((storageType == 1) && (hDim >= 1) && (vDim >= 1)) {
        _List*      resList = new _List;
        _Matrix     *rowSums     = new _Matrix (hDim,1,false,true),
        *columnSums  = new _Matrix (vDim,1,false,true);

       
        hyFloat totals = 0.0;

        if (theIndex) {
            for (long item = 0; item < lDim; item ++) {
                long idx = theIndex[item];
                if (idx>=0) {

                    hyFloat      v = theData[idx];

                    rowSums->theData[idx/vDim] += v;
                    columnSums->theData[idx%vDim] += v;
                    totals += v;
                }
            }
        } else {
            for (long rows = 0; rows < hDim; rows++) {
                hyFloat rowSum = 0.;

                for (long columns = 0; columns < vDim; columns ++) {
                    rowSum += theData[rows*vDim+columns];
                }

                rowSums->theData[rows] = rowSum;
                totals += rowSum;
            }

            for (long columns = 0; columns < vDim; columns++) {
                hyFloat colSum = 0.;

                for (long rows = 0; rows < hDim; rows ++) {
                    colSum += theData[rows*vDim+columns];
                }

                columnSums->theData[columns] = colSum;
            }
        }

        (*resList) < rowSums
                   < columnSums
                   < new _Constant (totals);

        return resList;

    }
    return nil;
}

//_____________________________________________________________________________________________

_Matrix* _Matrix::NeighborJoin (bool methodIndex, HBLObjectRef cache)
{
    long          specCount = GetHDim();

    if (storageType != 1 ||  specCount!= GetVDim() || specCount < 4) {
        HandleApplicationError ("NeigborJoin needs a square numeric matrix of dimension >= 4");
        return    new _Matrix;
    }

    CheckIfSparseEnough (true);

    _Matrix              netDivergence (specCount,1,false,true);
    _SimpleList          useColumn     (specCount,0,1),
                         columnIndex   (specCount,0,1);

    _Matrix*             res = (_Matrix* )_returnMatrixOrUseCache((specCount+1)*2,3,_NUMERICAL_TYPE,false,cache);

    for (long k=0; k<specCount ; k++) {
        for (long j=0; j<k; j++) {
            hyFloat d = theData[j*specCount+k];

            netDivergence.theData[k] += d;
            netDivergence.theData[j] += d;

        }
        res->theData[k*3+2] = 1;
    }

    long   cladesMade = 1;

    while (cladesMade < specCount) {
        hyFloat      min = 1.e100;

        long            minIndex  = -1,
                        minIndex2 = -1,
                        minIndexR = -1,
                        minIndexC = -1,
                        k = specCount-1-cladesMade;

        hyFloat      recRemaining = 1./k;

        if (cladesMade == specCount-1) {
            minIndex = useColumn.list_data[1];

            hyFloat d = theData[minIndex];

            if ((d<0)&&methodIndex) {
                d = 0;
            }

            k = columnIndex.list_data[1];

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
            long c1 = useColumn.list_data[i];

            for (long j=0; j<i; j=j+1) {
                long c2 = useColumn.list_data[j];

                //if (c2>=c1)
                //break;

                hyFloat d = theData[c2*specCount+c1]-(netDivergence.theData[c1]+netDivergence.theData[c2])*recRemaining;

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
            _String err = _String ("Invalid distance matrix passed to NeighborJoin. Matrices written onto ") & hy_messages_log_name;
            ReportWarning ((_String*)toStr());
            ReportWarning (_String((_String*)netDivergence.toStr()));
            ReportWarning (_String((_String*)useColumn.toStr()));
            HandleApplicationError (err);
            DeleteObject (res);
            return new _Matrix;
        }

        hyFloat      D  = theData[minIndex*specCount+minIndex2],
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

        long    m = columnIndex.list_data [minIndexC],
                n = columnIndex.list_data [minIndexR];

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
            long  k2 = useColumn.list_data[k];

            if (k2>=minIndex) {
                if (k2 == minIndex) {
                    k++;
                }
                break;
            }

            hyFloat d2 = theData[k2*specCount+minIndex]+theData[k2*specCount+minIndex2],
                       t  =  (d2-d)*.5;

            netDivergence.theData  [k2]               += t-d2;
            theData [k2*specCount+minIndex]            = t;
            netDivergence.theData[minIndex]           += t;

        }

        for (; k<useColumn.lLength; k++) {
            long  k2 = useColumn.list_data[k];
            if (k2 >= minIndex2) {
                if (k2 == minIndex2) {
                    k++;
                }
                break;
            }

            hyFloat  d2 = theData[minIndex*specCount+k2]+theData[k2*specCount+minIndex2],
                        t =  (d2-d)*.5;

            netDivergence.theData [k2]                  += t-d2;
            theData[minIndex*specCount+k2]               = t;
            netDivergence.theData[minIndex]             += t;

        }

        //for (k=minIndex2+1;k<ds.species; k=k+1)
        for (; k<useColumn.lLength; k++) {
            long  k2 = useColumn.list_data[k];

            hyFloat  d2 = theData[minIndex*specCount+k2]+theData[minIndex2*specCount+k2],
                        t =  (d2-d)*.5;

            netDivergence.theData [k2]                   += t-d2;
            theData[minIndex*specCount+k2]                = t;
            netDivergence.theData[minIndex]              += t;
        }

        columnIndex.list_data[minIndexR] = specCount+cladesMade-1;
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
_Matrix*        _Matrix::MakeTreeFromParent (long specCount, HBLObjectRef cache) {
    if (is_empty()) {
        return new _Matrix;
    }
    
    try {

        if (specCount<0L ) {
            throw (_String ("Parameter to ") & __PRETTY_FUNCTION__ & " must be greater than or equal to 0");
        }
        
        if (GetVDim () != 3) {
            throw (_String ("Expected a matrix with 3 columns"));
        }
        if (GetHDim () <= 2*specCount + 1) {
            throw (_String ("Expected a matrix with at least ") & (2*specCount + 1) & " rows");
        }

        const long result_rows = 2*(specCount+1);
        _Matrix     *tree = (_Matrix* )_returnMatrixOrUseCache(result_rows,5,_NUMERICAL_TYPE,false,cache),
                    CI  (2*(specCount+1),1,false,true);


        for (long kk = 0; kk < specCount-1; kk++) {
            tree->theData[kk*5+4] = -1; // set parent records to
        }

        long cladesMade = 0L;

        for (long nodeID2 = 0L; nodeID2 < specCount; nodeID2 ++) {
            long        nodeID       = nodeID2,
                        nodeDepth    = 0,
                        saveNodeID   = nodeID,
                        parentID     = theData[nodeID*3],
                        layoutOffset = cladesMade,
                        m,
                        n;

            while (parentID>=0) {
                long idx = parentID-specCount;
                if (idx < 0 || idx >= result_rows) {
                    throw (_String ("Invalid parent index in row ") & nodeID2);
                }
                n = tree->theData[idx*5+4];
                if (n >= 0) {
                    layoutOffset = n+tree->theData[idx*5+3];
                    break;
                }
                parentID  = theData[parentID*3];
            }

            parentID   = theData[nodeID*3];

            while (parentID>=0) {
                n = parentID-specCount;
                if (n < 0 || n >= result_rows) {
                    throw (_String ("Invalid parent index in row ") & nodeID);
                }
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
                    if (m < 0 || m >= result_rows) {
                        throw (_String ("Invalid parent index in row ") & nodeID);
                    }
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
    } catch (const _String& error) {
        HandleApplicationError(error);
        return new _Matrix (1,1,false,true);
    }
}


//_____________________________________________________________________________________________
hyFloat      _Matrix::FisherExact (hyFloat p1, hyFloat p2, hyFloat p3)
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

void        _Matrix::SimplexHelper1 (long rowIndex, _SimpleList& columnList, long columnCount, bool useAbsValue, long& maxIndex, hyFloat& maxValue)
// find the maximum element (using absolute value of not) in row rowIndex+1 of this matrix,
// over first columnCount columns indexed by columnList
//  column indexing is offset by + 1 to account for the first column not being eligible for pivoting
{
    if (columnCount <= 0) {
        maxValue = 0.0;
    } else {
        rowIndex = (rowIndex+1)*vDim;
        maxIndex = columnList.list_data[0];
        maxValue = theData[rowIndex+maxIndex+1];
        for (long k=1; k<columnCount; k++) {
            hyFloat t = useAbsValue?
                           (fabs(theData[rowIndex+columnList.list_data[k]+1])-fabs(maxValue))
                           :(theData[rowIndex+columnList.list_data[k]+1]-maxValue);
            if (t>0.) {
                maxValue = theData[rowIndex+columnList.list_data[k]+1];
                maxIndex = columnList.list_data[k];
            }
        }
    }
}

//_____________________________________________________________________________________________

void        _Matrix::SimplexHelper2 (long& pivotIndex, long columnToExamine, hyFloat eps)
{
    long            m = hDim-2,
                    n = vDim-1,
                    i = 0;

    hyFloat      q1,
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
                hyFloat q0, qp;
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
    hyFloat piv = 1./theData[(ip+1)*vDim+kp+1];
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
_Matrix*    _Matrix::SimplexSolve (hyFloat desiredPrecision ) {
// this function is adapted from the Num. Recipes in C version; but with 0 indexing
// hyphy primitives
// and without goto labels

// the of dimension RxC is interpreted as follows
// R-1 constraints
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
// if an kEmptyString matrix is returned - no feasible solution could be found
// if a 1x1 matrix is returned - the objective function is unbounded

    try {

        long n = vDim-2, // number of variables
             m = hDim-1; // number of constraints

        if (is_numeric() && n>0 && m>0) {
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
                        hyFloat t = (*this)(i,n+1);
                        if (t<0.0) {
                            m1++;
                        } else if (t>0.0) {
                            m2++;
                        } else {
                            m3++;
                        }
                        if ((*this)(i,0) < 0.0) {
                            throw _String("Negative values are not allowed in the first column of the simplex tableau");
                         }
                    }
                }


                // copy coefficients into the temp matrix, sorting the constraints in the <=, >= and == order
                {
                    for (long i=1, t1=0, t2=0, t3=0; i<=m; i++) {
                        hyFloat t   = (*this)(i,n+1);
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
                    for (long i=0; i<m2; l3.list_data[i] = 1,i++) ; // slack variables in the list of 'basis' variables
                    for (long k=0; k<=n; k++) { // compute the auxiliary objective function
                        hyFloat q = 0.;
                        for (long k2 = m1+1; k2<=m; k2++) {
                            q += tempMatrix(k2,k);
                        }
                        tempMatrix.Store (m+1,k,-q);
                    }
                    while (1) { // initial artifical construct
                        long        pivotColumn,
                                    ip;
                        hyFloat  pivotValue;

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
                                if (iposv.list_data[ip] == ip + n) {
                                    tempMatrix.SimplexHelper1 (ip,l1,nl1,true,pivotColumn, pivotValue);
                                    if (pivotValue > desiredPrecision) {
                                        goto one;
                                    }
                                }
                            }
                            for (long i = m1; i<m1+m2; i++)
                                if (l3.list_data[i-m1] == 1)
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
                        if (iposv.list_data[ip] >= n+m1+m2) {
                            long k = 0;
                            for (k=0; k<nl1; k++)
                                if (l1.list_data[k] == pivotColumn) {
                                    break;
                                }
                            nl1--;
                            for (long i2=k; i2<nl1; i2++) {
                                l1.list_data[i2] = l1.list_data[i2+1];
                            }
                        } else {
                            long k2 = iposv.list_data[ip] - m1 - n;
                            if (k2 >= 0 && l3.list_data[k2]) {
                                l3.list_data[k2] = 0;
                                tempMatrix.theData[(m+1)*tempMatrix.vDim + pivotColumn + 1] ++;
                                for (long i=0; i<m+2; i++) {
                                    tempMatrix.theData[i*tempMatrix.vDim + pivotColumn + 1] *= -1.0;
                                }
                            }
                        }
                        long s = izrov.list_data[pivotColumn];
                        izrov.list_data[pivotColumn] = iposv.list_data[ip];
                        iposv.list_data[ip] = s;
                    }// end of phase 1
                }

                while (1) {
                    long            pivotColumn,
                                    pivotRow;
                    hyFloat      pivotValue;

                    tempMatrix.SimplexHelper1 (-1,l1,nl1,false,pivotColumn,pivotValue);
                    if (pivotValue < desiredPrecision) { // done!
                        // produce the final solution
                        _Matrix * resMatrix = new _Matrix (1,n+1,false,true);
                        resMatrix->Store(0,0,doMaximize?tempMatrix(0,0):-tempMatrix(0,0));
                        for (long k=0; k<iposv.lLength; k++)
                            if (iposv.list_data[k]<n) {
                                resMatrix->Store(0,iposv.list_data[k]+1,tempMatrix(k+1,0));
                            }
                        return resMatrix;
                    }
                    tempMatrix.SimplexHelper2 (pivotRow,pivotColumn,desiredPrecision);
                    if (pivotRow<0) {
                        return new _Matrix (1,1,false,true);
                    }
                    tempMatrix.SimplexHelper3 (m-1,n-1,pivotRow,pivotColumn);
                    long s = izrov.list_data[pivotColumn];
                    izrov.list_data[pivotColumn] = iposv.list_data[pivotRow];
                    iposv.list_data[pivotRow] = s;
                }

            }
        } else {
            throw _String("SimplexSolve requires a numeric matrix with > 1 row and > 2 columns");
        }


    } catch (_String const& err) {
        HandleApplicationError (err);
    }
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
HBLObjectRef   _Matrix::DirichletDeviate (void)
{
    /* -----------------------------------------------------------
        DirichletDeviate()
            Generate vector of random deviates from the Dirichlet
            distribution defined by contents of this matrix as
            hyperparameters (a > 0).
       ----------------------------------------------------------- */
    try {

        long        dim;

        hyFloat  denom   = 0.;

        _Matrix     res (1, dim = GetHDim()*GetVDim(), false, true);    // row vector


        if (!is_numeric()) {
            throw _String("Only numeric vectors can be passed to DirichletDeviate");
        }

        if (is_row() || is_column ()) {
            // generate a random deviate from gamma distribution for each hyperparameter
            
            for (long i = 0L; i < dim; i++) {
                if (theData[i] < 0.) {
                    throw _String("Dirichlet not defined for negative parameter values.");
                }

                res.Store (0, i, gammaDeviate(theData[i]));
                denom += res(0,i);
            }

            // normalize by sum
            for (long i = 0; i < dim; i++) {
                res.Store (0, i, res(0,i)/denom);
            }

            return (HBLObjectRef) res.makeDynamic();
        } else {
            throw _String("Argument must be a row- or column-vector.");
        }
    } catch (_String const& err) {
        HandleApplicationError (err);
    }
    return new _Matrix (1,1,false,true);
}



//_____________________________________________________________________________________________
HBLObjectRef   _Matrix::GaussianDeviate (_Matrix & cov)
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

    //ReportWarning (_String("Entered _Matrix::GaussianDeviate() with cov = ") & (_String *)(cov.toStr()));

    try {

        if (storageType != 1 || GetHDim() > 1) {
            HandleApplicationError (_String("ERROR in _Matrix::GaussianDeviate(), expecting to be called on numeric row vector matrix, current dimensions: ") & GetHDim() & "x" & GetVDim());
            return new _Matrix;
        }

        long kdim = GetVDim();    // number of entries in this _Matrix object as vector of means

        if (cov.check_dimension(kdim, kdim)) {
            _Matrix* cov_cd = (_Matrix *) cov.CholeskyDecompose(),
                    * gaussvec = new _Matrix (1, kdim, false, true);

            //ReportWarning (_String("\nCholesky decomposition of cov = ") & (_String *) cov_cd->toStr());

            // fill column vector with independent standard normal deviates
            for (long i = 0L; i < kdim; i++) {
                gaussvec->Store (0, i, gaussDeviate());
            }

            //ReportWarning (_String ("\nvector of gaussian deviates = ") & (_String *) gaussvec.toStr());

            // left multiply vector by Cholesky decomposition of covariance matrix
            *gaussvec *= *cov_cd;

            // shift mean
            for (long i = 0L; i < kdim; i++) {
                gaussvec->Store (0, i, (*gaussvec)(0,i) + theData[i]);
            }

            DeleteObject (cov_cd);
            return gaussvec;
        } else {
            throw (_String("Error in _Matrix::GaussianDeviate(), incompatible dimensions in covariance matrix: ") & cov.GetHDim() & "x" & cov.GetVDim());

        }
    } catch (const _String& err) {
        HandleApplicationError (err);
    }

    return new _Matrix;
}


//_____________________________________________________________________________________________
HBLObjectRef   _Matrix::MultinomialSample (_Constant *replicates) {
    
    try {
        _List         reference_manager;
        
        long          values      = GetHDim();
        unsigned long samples     = replicates?replicates->Value ():0;

        _Matrix     *eval    = (_Matrix*)Compute (),
                    * sorted = nil,
                    * result = nil;

        if (samples == 0UL) {
            throw _String ("Expected a numerical (>=1) value for the number of replicates");
        } else if ( ! eval->is_numeric() || GetVDim() != 2 || values < 2) {
            throw _String ("Expecting numerical Nx2 (with N>=1) matrix.");
        } else {
            _Constant one (1.);
            sorted = (_Matrix*) eval->SortMatrixOnColumn(&one, nil);
            reference_manager < sorted;
            hyFloat      sum = 0.;

            for (long n = 1L; n < 2*values; n+=2L) {
                hyFloat v = sorted->theData[n];
                if (v < 0.) {
                    sum = 0.;
                    break;
                }
                sum += v;
            }
            if (CheckEqual (sum, 0.)) {
                throw _String ("The probabilities (second column) cannot add to 0 or be negative");
            } else {
                sum = 1./sum;

                _Matrix     *raw_result  = new _Matrix (1, values, false, true),
                *normalized  = new _Matrix (1, values, false, true);

                reference_manager <raw_result;
                reference_manager <normalized;

                
                for (long v = 0; v < values; v++) {
                    normalized->theData[values-1-v] = sorted->theData[1+2*v] * sum;
                }

 
                 hyFloat  seconds_accumulator = .0,
                            temp;

                for (unsigned long it = 0UL; it < samples; it++) {
                     raw_result->theData[DrawFromDiscrete(normalized->theData, values)] += 1.;
                }

                result = new _Matrix (values, 2, false, true);

                for (long v = 0; v < values; v++) {
                    result->theData[2*v]   = (long)sorted->theData[2*(values-1-v)];
                    result->theData[2*v+1] = raw_result->theData[v];
                }

                
                return result;
            }
        }
    }
    catch (_String const& err) {
        HandleApplicationError (err);
    }
    return new _Matrix;
}

//_____________________________________________________________________________________________
HBLObjectRef   _Matrix::InverseWishartDeviate (_Matrix & df)
{
    /* ---------------------------------------------------
        InverseWishartDeviate()
            Generates a random matrix whose inverse
            has the Wishart distribution with this matrix
            supplying the covariance matrix parameter and
            a degrees of freedom vector argument.
       --------------------------------------------------- */

    try {
        long        n       = GetHDim();


        if (!is_square_numeric()) {
            throw _String("Expecting a numerical square matrix.");
        }

        else if (!df.is_numeric() || !df.check_dimension(n,1)) {
            throw _String("Expecting numerical column vector for second argument (degrees of freedom).");
        } else {
            // compute Cholesky factor for this matrix inverse, extract the diagonal
            _List   reference_manager;
            
            _Matrix * inv       = (_Matrix *) Inverse(nil);
            _Matrix * invCD     = (_Matrix *) (inv->CholeskyDecompose());
            
            DeleteObject (inv);
            reference_manager < invCD;
            
            return WishartDeviate (df, *invCD);
        }
    } catch (const _String& err) {
        HandleApplicationError (err);
    }
    return new _Matrix;
}

//_____________________________________________________________________________________________
HBLObjectRef   _Matrix::WishartDeviate (_Matrix & df) {
    _Matrix     diag;   // calls default constructor
    return WishartDeviate (df, diag);
}

//_____________________________________________________________________________________________

HBLObjectRef   _Matrix::WishartDeviate (_Matrix & df, _Matrix & decomp) {
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


    try {

        long        n   = GetHDim();

        _Matrix     rdeviates (n, n, false, true),
                    rd_transpose;


        if (!(df.is_row () || df.is_column())) {
            throw _String("Expecting row vector for degrees of freedom argument.");
        } else if (df.is_column()) {
            df.Transpose(); // convert column vector to row vector
        }

        if (decomp.is_empty()) {    // no second argument, perform Cholesky decomposition
            if (!is_square_numeric()) {
                throw _String("Expecting square numeric matrix.");
            } else {
                _Matrix     * cholesky = (_Matrix *) CholeskyDecompose();

                if (cholesky->GetHDim() > 0) {
                    decomp = *cholesky;
                    DeleteObject (cholesky);
                } else {
                    return cholesky;  // empty _Matrix from error in CholeskyDecompose()
                }
            }
        }


        // populate diagonal with square root of i.i.d. chi-square random deviates
        for (unsigned long i = 0UL; i < n; i++) {
            rdeviates.Store (i, i, sqrt(chisqDeviate(df(0,i)-i+1)) );

            // populate upper triagonal with i.i.d. standard normal N(0,1) deviates
            for (unsigned long j = i+1UL; j < n; j++) {
                rdeviates.Store (i, j, gaussDeviate());
            }
        }

  
        // result is obtained from D^T B D, where B = A^T A, ^T is matrix transpose
        rd_transpose = rdeviates;
        rd_transpose.Transpose();
        rd_transpose *= rdeviates;  // A^T A
        rd_transpose *= decomp; // A^T A D

        decomp.Transpose();
        decomp *= rd_transpose; // D^T A^T A D
 
        return (HBLObjectRef) decomp.makeDynamic();
    } catch (const _String& err) {
        HandleApplicationError(err);
    }
    return new _Matrix;
}

//-----------------------------------------------------------------------------------------------------------------

HBLObjectRef _returnMatrixOrUseCache (long nrow, long ncol, long type, bool is_sparse, HBLObjectRef cache) {
    if (cache && cache->ObjectClass() == MATRIX) {
        _Matrix *cached_mx = (_Matrix*)cache;
        if (cached_mx->check_dimension(nrow, ncol) && cached_mx->has_type (type) && cached_mx->is_dense() == !is_sparse) {
            cached_mx->Clear(false);
        } else {
            cached_mx->Clear();
            _Matrix::CreateMatrix (cached_mx, nrow, ncol, is_sparse, type == _NUMERICAL_TYPE ? true : false);
        }
        //cached_mx->AddAReference();
        return cached_mx;
    }
    return new _Matrix (nrow, ncol, is_sparse, type == _NUMERICAL_TYPE ? true : false);
}

