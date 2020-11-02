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
#include <ctype.h>
#include <float.h>


#include "global_things.h"
#include "likefunc.h"

extern  long likeFuncEvalCallCount;


using namespace hy_global;
using namespace hy_env;



#define __ll_loop_preamble if (tcc) { \
    if (__builtin_expect (parentTCCIBit==_HY_BITMASK_WIDTH_,0)) {\
        parentTCCIBit   = 0;\
        parentTCCIIndex ++;\
    }\
    if (__builtin_expect(siteID > siteFrom && (tcc->list_data[parentTCCIIndex] & bitMaskArray.masks[parentTCCIBit]) > 0,0)) {\
        if (!isLeaf) {\
            childVector     += alphabetDimension;\
            if (__builtin_expect(++currentTCCBit == _HY_BITMASK_WIDTH_,0)) {\
                currentTCCBit   = 0;\
                currentTCCIndex ++;\
            }\
        }\
        parentTCCIBit++;\
        continue;\
    }\
    parentTCCIBit++;\
    }\
    hyFloat  const * _hprestrict_ tMatrix = transitionMatrix;\
    hyFloat  sum         = 0.0;\
    long        didScale = 0;\

#define __ll_loop_epilogue if (didScale) {\
if (siteCorrectionCounts) {\
    siteCorrectionCounts [siteOrdering.list_data[siteID]] += didScale;\
}\
if (tcc) {\
    long cparentTCCIIndex   =   parentTCCIIndex,\
    cparentTCCIBit   =   parentTCCIBit;\
    hyFloat              scM;\
    if (didScale < 0) {\
        scM = _lfScalingFactorThreshold;\
        if (didScale < -1) {\
            for (long k = 0; k < -didScale-1; k++) {\
                scM *= _lfScalingFactorThreshold;\
            }\
        }\
    } else {\
        scM = _lfScalerUpwards;\
        if (didScale > 1) {\
            for (long k = 1; k < didScale; k++) {\
                scM *= _lfScalerUpwards;\
            }\
        }\
    }\
    for (long sid = siteID + 1; sid < siteTo; sid++,cparentTCCIBit++) {\
        if (cparentTCCIBit == _HY_BITMASK_WIDTH_) {\
            cparentTCCIBit   = 0;\
            cparentTCCIIndex ++;\
        }\
        if ((tcc->list_data[cparentTCCIIndex] & bitMaskArray.masks[cparentTCCIBit]) > 0) {\
            if (siteCorrectionCounts) {\
                siteCorrectionCounts [siteOrdering.list_data[sid]] += didScale;\
            }\
            scalingAdjustments   [parentCode*siteCount + sid] *= scM;\
            localScalerChange                               += didScale * theFilter->theFrequencies.get(siteOrdering.list_data[sid]);\
        } else {\
            break;\
        }\
    }\
}\
}\

#ifdef _SLKP_USE_SSE_INTRINSICS
inline double _sse_sum_2 (__m128d const & x) {
    return _mm_cvtsd_f64(_mm_hadd_pd(x, x));
}
#endif


template<long D> inline void __ll_handle_matrix_transpose (hyFloat const * __restrict transitionMatrix, hyFloat * __restrict tMatrixT) {
    long i = 0L;
    for (long r = 0L; r < D; r++) {
        #pragma unroll(4)
        #pragma GCC unroll 4
        for (long c = 0L; c < D; c++, i++) {
            tMatrixT[c*D+r] = transitionMatrix[i];
        }
    }
}

template<long D> inline bool __ll_handle_conditional_array_initialization ( long * __restrict lNodeFlags, bool isLeaf, long nodeCode, long setBranch, long iNodes, long siteID, long siteFrom, long siteCount, _SimpleList&            siteOrdering, hyFloat * __restrict parentConditionals, hyFloat const * __restrict tMatrix, _Vector const*     lNodeResolutions, hyFloat *& childVector, _SimpleList const* tcc, long &currentTCCBit, long& currentTCCIndex, hyFloat * & lastUpdatedSite, long* __restrict  setBranchTo) {
    

    if (isLeaf) {
        long siteState;
        if (setBranch != nodeCode + iNodes) {
            siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.list_data[siteID]] ;
        } else {
            siteState = setBranchTo[siteOrdering.list_data[siteID]] ;
        }
        if (__builtin_expect(siteState >= 0L,1)) {
            // a single character state; sweep down the appropriate column
            /*if (likeFuncEvalCallCount == 15098 && nodeCode == 3706 && siteID == 91) {
                fprintf (stderr, "\nSITE CHECK: ID %ld, STATE %ld\n", nodeCode, siteState);
                for (long e = 0; e < 4; e++) {
                    fprintf (stderr, "%ld => %lg\n", e, parentConditionals[e]);
                }
            }*/
            
            #pragma unroll(4)
            #pragma GCC unroll 4
            for (long k = 0L; k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState+D*k];
            }
            /*if (likeFuncEvalCallCount == 15098 && nodeCode == 3706 && siteID == 91) {
                fprintf (stderr, "\nSITE CHECK: ID %ld, STATE %ld\n", nodeCode, siteState);
                for (long e = 0; e < 4; e++) {
                    fprintf (stderr, "%ld => %lg\n", e, parentConditionals[e]);
                }
            }*/
            return true;
        } else {
            childVector = lNodeResolutions->theData + (-siteState-1) * D;
        }
    } else {
        if (tcc) {
            if (__builtin_expect((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom,0)) {
                #pragma unroll(4)
                #pragma GCC unroll 4
                for (long k = 0L; k < D; k++) {
                    childVector[k] = lastUpdatedSite[k];
                }
            }
            if (__builtin_expect(++currentTCCBit == _HY_BITMASK_WIDTH_,0)) {
                currentTCCBit   = 0;
                currentTCCIndex ++;
            }
            lastUpdatedSite = childVector;
        }
    }
    return false;
}

inline bool __ll_handle_conditional_array_initialization_generic ( long * __restrict lNodeFlags, bool isLeaf, long nodeCode, long setBranch, long iNodes, long siteID, long siteFrom, long siteCount, _SimpleList&            siteOrdering, hyFloat * __restrict parentConditionals, hyFloat const * __restrict tMatrix, _Vector const*     lNodeResolutions, hyFloat *& childVector, _SimpleList const* tcc, long &currentTCCBit, long& currentTCCIndex, hyFloat * & lastUpdatedSite, long* __restrict  setBranchTo, long D) {
    

    if (isLeaf) {
        long siteState;
        if (setBranch != nodeCode + iNodes) {
            siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.list_data[siteID]] ;
        } else {
            siteState = setBranchTo[siteOrdering.list_data[siteID]] ;
        }
        if (__builtin_expect(siteState >= 0L,1)) {
            // a single character state; sweep down the appropriate column
            #pragma unroll(4)
            #pragma GCC unroll 4
            for (long k = 0L; k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState+D*k];
            }
            return true;
        } else {
            childVector = lNodeResolutions->theData + (-siteState-1) * D;
        }
    } else {
        if (tcc) {
            if (__builtin_expect((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom,0)) {
                #pragma unroll(4)
                #pragma GCC unroll 4
                for (long k = 0L; k < D; k++) {
                    childVector[k] = lastUpdatedSite[k];
                }
            }
            if (__builtin_expect(++currentTCCBit == _HY_BITMASK_WIDTH_,0)) {
                currentTCCBit   = 0;
                currentTCCIndex ++;
            }
            lastUpdatedSite = childVector;
        }
    }
    return false;
}

#ifdef _SLKP_USE_SSE_INTRINSICS
template<long D, long S> inline void __ll_handle_block10_product_sum_linear (hyFloat const* _hprestrict_ tT, hyFloat * _hprestrict_ childVector, __m128d * lastTotal) {
        
    __m128d T [5] = {_mm_loadu_pd (tT+S),
                     _mm_loadu_pd (tT+S+2),
                     _mm_loadu_pd (tT+S+4),
                     _mm_loadu_pd (tT+S+6),
                     _mm_loadu_pd (tT+S+8)};
    
    __m128d C [5] =  {  _mm_loadu_pd (childVector+S),
                        _mm_loadu_pd (childVector+S+2),
                        _mm_loadu_pd (childVector+S+4),
                        _mm_loadu_pd (childVector+S+6),
                        _mm_loadu_pd (childVector+S+8)
                     };

     if (S) {
        T[0] = _mm_mul_pd (T[0], C[0]);
        T[1] = _mm_mul_pd (T[1], C[1]);
        T[2] = _mm_mul_pd (T[2], C[2]);
        T[3] = _mm_mul_pd (T[3], C[3]);
        T[4] = _mm_mul_pd (T[4], C[4]);
        T[0] = _mm_add_pd (T[0],T[1]);
        T[2] = _mm_add_pd (T[2],T[3]);
        *lastTotal = _mm_add_pd (*lastTotal, T[4]);
        T[0] = _mm_add_pd (T[0],T[2]);
        *lastTotal = _mm_add_pd (T[0], *lastTotal);
    } else {
        *lastTotal = _mm_mul_pd (T[0], C[0]);
        T[1] = _mm_mul_pd (T[1], C[1]);
        T[2] = _mm_mul_pd (T[2], C[2]);
        T[3] = _mm_mul_pd (T[3], C[3]);
        T[4] = _mm_mul_pd (T[4], C[4]);
        T[0] = _mm_add_pd (T[1],T[3]);
        T[2] = _mm_add_pd (T[2],T[4]);
        T[0] = _mm_add_pd (T[0],T[2]);
        *lastTotal = _mm_add_pd (T[0], *lastTotal);
    }
}
#endif

#ifdef _SLKP_USE_AVX_INTRINSICS
template<long D, long S> inline void __ll_handle_block20_product_sum_linear (hyFloat const* _hprestrict_ tT, hyFloat * _hprestrict_ childVector, __m256d * lastTotal) {
        
    __m256d T [5] = {_mm256_loadu_pd (tT+S),
                     _mm256_loadu_pd (tT+S+4),
                     _mm256_loadu_pd (tT+S+8),
                     _mm256_loadu_pd (tT+S+12),
                     _mm256_loadu_pd (tT+S+16)};
    __m256d C [5] =  {  _mm256_loadu_pd (childVector+S),
                        _mm256_loadu_pd (childVector+S+4),
                        _mm256_loadu_pd (childVector+S+8),
                        _mm256_loadu_pd (childVector+S+12),
                        _mm256_loadu_pd (childVector+S+16)
                     };

    #ifdef _SLKP_USE_FMA3_INTRINSICS
        if (S) {
            T[0] = _mm256_mul_pd (T[0], C[0]);
            T[1] = _mm256_mul_pd (T[1], C[1]);
            T[2] = _mm256_fmadd_pd (T[2], C[2], T[0]); // 0+2
            T[3] = _mm256_fmadd_pd (T[3], C[3], T[1]); // 1+3
            T[4] = _mm256_fmadd_pd (T[4], C[4], *lastTotal); // 4 + last
            T[2] = _mm256_add_pd (T[2],T[3]); // 0+1+2+3
            *lastTotal = _mm256_add_pd(T[4], T[2]);
        } else {
            T[0] = _mm256_mul_pd (T[0], C[0]);
            T[1] = _mm256_mul_pd (T[1], C[1]);
            T[2] = _mm256_fmadd_pd (T[2], C[2], T[0]); // 0+2
            *lastTotal = _mm256_fmadd_pd (T[3], C[3], T[1]); // 1+3
            T[0] = _mm256_fmadd_pd(T[4], C[4],T[2]); // 0+2+4
            *lastTotal = _mm256_add_pd(T[0], *lastTotal);
        }
    #else
        if (S) {
            T[0] = _mm256_mul_pd (T[0], C[0]);
            T[1] = _mm256_mul_pd (T[1], C[1]);
            T[2] = _mm256_mul_pd (T[2], C[2]);
            T[3] = _mm256_mul_pd (T[3], C[3]);
            T[4] = _mm256_mul_pd (T[4], C[4]);
            T[0] = _mm256_add_pd (T[0],T[1]);
            T[2] = _mm256_add_pd (T[2],T[3]);
            *lastTotal = _mm256_add_pd (*lastTotal, T[4]);
            T[0] = _mm256_add_pd (T[0],T[2]);
            *lastTotal = _mm256_add_pd (T[0], *lastTotal);
        } else {
            *lastTotal = _mm256_mul_pd (T[0], C[0]);
            T[1] = _mm256_mul_pd (T[1], C[1]);
            T[2] = _mm256_mul_pd (T[2], C[2]);
            T[3] = _mm256_mul_pd (T[3], C[3]);
            T[4] = _mm256_mul_pd (T[4], C[4]);
            T[0] = _mm256_add_pd (T[1],T[3]);
            T[2] = _mm256_add_pd (T[2],T[4]);
            T[0] = _mm256_add_pd (T[0],T[2]);
            *lastTotal = _mm256_add_pd (T[0], *lastTotal);
        }
    #endif

}

#endif

#ifdef _SLKP_USE_SSE_INTRINSICS
template<long D, long S> inline __m128d __ll_handle_block10_product_sum (hyFloat const* _hprestrict_ transposedMatrix, hyFloat * _hprestrict_ childVector, hyFloat * _hprestrict_ parentConditionals, __m128d * grandTotal) {
    
        __m128d P [5] = {
                    _mm_set1_pd(0.),
                    _mm_set1_pd(0.),
                    _mm_set1_pd(0.),
                    _mm_set1_pd(0.),
                    _mm_set1_pd(0.)
                    };
                                    
    hyFloat const * __restrict tT = transposedMatrix;

    for (long col_idx = 0; col_idx < D; col_idx++, tT+=D) {
        if (childVector[col_idx] > 0.) {
            __m128d C = _mm_set1_pd(childVector[col_idx]);
                                        
            __m128d T [5] = {_mm_loadu_pd (tT+S),
                             _mm_loadu_pd (tT+S+2),
                             _mm_loadu_pd (tT+S+4),
                             _mm_loadu_pd (tT+S+6),
                             _mm_loadu_pd (tT+S+8)};
            
            P[0] = _mm_add_pd (_mm_mul_pd (T[0], C), P[0]);
            P[1] = _mm_add_pd (_mm_mul_pd (T[1], C), P[1]);
            P[2] = _mm_add_pd (_mm_mul_pd (T[2], C), P[2]);
            P[3] = _mm_add_pd (_mm_mul_pd (T[3], C), P[3]);
            P[4] = _mm_add_pd (_mm_mul_pd (T[4], C), P[4]);
        }
    }

    
    P[0] = _mm_mul_pd(_mm_loadu_pd(parentConditionals+S),     P[0]);
    P[1] = _mm_mul_pd(_mm_loadu_pd(parentConditionals+S+2),   P[1]);
    P[2] = _mm_mul_pd(_mm_loadu_pd(parentConditionals+S+4),   P[2]);
    P[3] = _mm_mul_pd(_mm_loadu_pd(parentConditionals+S+6),   P[3]);
    P[4] = _mm_mul_pd(_mm_loadu_pd(parentConditionals+S+8),   P[4]);
    
    _mm_storeu_pd(parentConditionals+S, P[0]);
    _mm_storeu_pd(parentConditionals+S+2, P[1]);
    _mm_storeu_pd(parentConditionals+S+4, P[2]);
    _mm_storeu_pd(parentConditionals+S+6, P[3]);
    _mm_storeu_pd(parentConditionals+S+8, P[4]);

    P[0] =_mm_add_pd (P[0],P[1]);
    P[1] =_mm_add_pd (P[2],P[3]);
    P[0] =_mm_add_pd (P[0],P[4]);
    P[2] =_mm_add_pd (P[0],P[1]);
    if (grandTotal) {
        *grandTotal = _mm_add_pd (P[2],*grandTotal);
        return *grandTotal;
    }
    return P[2];
    
}
#endif

#ifdef _SLKP_USE_AVX_INTRINSICS
template<long D, long S> inline __m256d __ll_handle_block20_product_sum (hyFloat const* _hprestrict_ transposedMatrix, hyFloat * _hprestrict_ childVector, hyFloat * _hprestrict_ parentConditionals, __m256d * grandTotal) {
    __m256d P [5] = {
                      _mm256_set1_pd(0.),
                      _mm256_set1_pd(0.),
                      _mm256_set1_pd(0.),
                      _mm256_set1_pd(0.),
                      _mm256_set1_pd(0.)
                    };
                                    
    hyFloat const * __restrict tT = transposedMatrix;

    for (long col_idx = 0; col_idx < D; col_idx++, tT+=D) {
        if (childVector[col_idx] > 0.) {
            __m256d C = _mm256_set1_pd(childVector[col_idx]);
                                        
            __m256d T [5] = {_mm256_loadu_pd (tT+S),
                             _mm256_loadu_pd (tT+S+4),
                             _mm256_loadu_pd (tT+S+8),
                             _mm256_loadu_pd (tT+S+12),
                             _mm256_loadu_pd (tT+S+16)};
            
            #ifdef _SLKP_USE_FMA3_INTRINSICS
                P[0] = _mm256_fmadd_pd (T[0], C, P[0]);
                P[1] = _mm256_fmadd_pd (T[1], C, P[1]);
                P[2] = _mm256_fmadd_pd (T[2], C, P[2]);
                P[3] = _mm256_fmadd_pd (T[3], C, P[3]);
                P[4] = _mm256_fmadd_pd (T[4], C, P[4]);
            #else
                P[0] = _mm256_add_pd (_mm256_mul_pd (T[0], C), P[0]);
                P[1] = _mm256_add_pd (_mm256_mul_pd (T[1], C), P[1]);
                P[2] = _mm256_add_pd (_mm256_mul_pd (T[2], C), P[2]);
                P[3] = _mm256_add_pd (_mm256_mul_pd (T[3], C), P[3]);
                P[4] = _mm256_add_pd (_mm256_mul_pd (T[4], C), P[4]);
            #endif
        }
    }

    
    P[0] = _mm256_mul_pd(_mm256_loadu_pd(parentConditionals+S),      P[0]);
    P[1] = _mm256_mul_pd(_mm256_loadu_pd(parentConditionals+S+4),    P[1]);
    P[2] = _mm256_mul_pd(_mm256_loadu_pd(parentConditionals+S+8),    P[2]);
    P[3] = _mm256_mul_pd(_mm256_loadu_pd(parentConditionals+S+12),   P[3]);
    P[4] = _mm256_mul_pd(_mm256_loadu_pd(parentConditionals+S+16),   P[4]);
    _mm256_storeu_pd(parentConditionals+S, P[0]);
    _mm256_storeu_pd(parentConditionals+S+4, P[1]);
    _mm256_storeu_pd(parentConditionals+S+8, P[2]);
    _mm256_storeu_pd(parentConditionals+S+12, P[3]);
    _mm256_storeu_pd(parentConditionals+S+16, P[4]);

    P[0] =_mm256_add_pd (P[0],P[1]);
    P[1] =_mm256_add_pd (P[2],P[3]);
    P[0] =_mm256_add_pd (P[0],P[4]);
    P[2] = _mm256_add_pd (P[0],P[1]);
    if (grandTotal) {
        *grandTotal = _mm256_add_pd (P[2],*grandTotal);
        return *grandTotal;
    }
    return P[2];
    
}
#endif

template<long D> inline void __ll_product_sum_loop (hyFloat const* _hprestrict_ tMatrix, hyFloat* _hprestrict_ childVector, hyFloat* _hprestrict_ parentConditionals, hyFloat& sum) {
    for (long p = 0; p < D; p++) {
        hyFloat      accumulator = 0.0;
        
        #pragma GCC unroll 8
        #pragma clang loop vectorize(enable)
        #pragma clang loop interleave(enable)
        #pragma clang loop unroll(enable)
        for (long c = 0; c < D; c++)
            accumulator +=  tMatrix[c]   * childVector[c];
        
        tMatrix               += D;
        sum += (parentConditionals[p] *= accumulator);
    }
}

inline void __ll_product_sum_loop_generic (hyFloat const* _hprestrict_ tMatrix, hyFloat* _hprestrict_ childVector, hyFloat* _hprestrict_ parentConditionals, hyFloat& sum, long const D) {
    for (long p = 0; p < D; p++) {
        hyFloat      accumulator = 0.0;
        
        #pragma GCC unroll 8
        #pragma clang loop vectorize(enable)
        #pragma clang loop interleave(enable)
        #pragma clang loop unroll(enable)
        for (long c = 0; c < D; c++)
            accumulator +=  tMatrix[c]   * childVector[c];
        
        tMatrix               += D;
        sum += (parentConditionals[p] *= accumulator);
    }
}

template<long D, bool ADJUST> inline void __ll_loop_handle_scaling (hyFloat& sum, hyFloat* _hprestrict_ parentConditionals, hyFloat* _hprestrict_ scalingAdjustments, long& didScale, long parentCode, long siteCount, long siteID, long& localScalerChange, long siteFrequency) {
    
    /*if (sum == 0.) {
        fprintf (stderr, "THE SUM IS EXACTLY ZERO parent code %ld\n", parentCode);
    }*/
    
    /*
     if (isnan (sum)) {
        HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch at branch " & parentCode & "; this is not a recoverable error and indicates some serious COVFEFE taking place.");
    }*/
    
    /*if (sum == 0.) {
        printf ("Exactly 0 in __ll_loop_handle_scaling: site = %ld, branch = %ld\n", siteID, parentCode);
    }
    if (isinf(sum)) {
        printf ("Infinity at __ll_loop_handle_scaling: site = %ld, branch = %ld\n", siteID, parentCode);
    }*/
    
    if (__builtin_expect(sum < _lfScalingFactorThreshold && sum > 0.0,0)) {
        
        hyFloat scaler = _computeBoostScaler(scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards, sum, didScale);
        
        //fprintf (stderr, "UP %ld (%ld) %lg\n", didScale, parentCode, scaler);
        
        /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
            fprintf (stderr, "UP %ld (%ld) %lg\n", didScale, parentCode, scaler);
        }*/
        if (didScale) {
            #pragma unroll(4)
            #pragma GCC unroll 4
            for (long c = 0; c < D; c++) {
                parentConditionals [c] *= scaler;
                //if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                //fprintf (stderr, "%ld=>%g\n", c, parentConditionals [c]);
                // }
                //if (isnan (parentConditionals [c])) {
                //    HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch at branch " & parentCode & "; this is not a recoverable error and indicates some serious COVFEFE taking place.");
                //}
            }
            
            if (siteFrequency == 1L) {
                localScalerChange += didScale;
            } else {
                localScalerChange += didScale * siteFrequency;
            }
            if (ADJUST)
                scalingAdjustments [parentCode*siteCount + siteID]  *= scaler;
        }
        
    } else {
        if (__builtin_expect(sum > _lfScalerUpwards,0)) {
            if (sum < HUGE_VAL) { // no point scaling an infinity
                
                hyFloat scaler = _computeReductionScaler (scalingAdjustments [parentCode*siteCount + siteID] * _lfScalingFactorThreshold, sum, didScale);
                /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                    fprintf (stderr, "DOWN %ld (%ld) %lg\n", didScale, parentCode, scaler);
                }*/
                
                if (didScale) {
                    #pragma unroll(4)
                    #pragma GCC unroll 4
                    for (long c = 0; c < D; c++) {
                        parentConditionals [c] *= scaler;
                        //if (isnan (parentConditionals [c])) {
                        //    HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch at branch " & parentCode & "; this is not a recoverable error and indicates some serious COVFEFE taking place.");
                        //}
                        //if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                        //   fprintf (stderr, "%ld=>%g\n", c, parentConditionals [c]);
                        // }
                    }
                    
                    if (siteFrequency == 1L) {
                        localScalerChange += didScale;
                    } else {
                        localScalerChange += didScale * siteFrequency;
                    }
                    if (ADJUST)
                        scalingAdjustments [parentCode*siteCount + siteID] *= scaler;
                }
            }
        }
    }
}

template<bool ADJUST> inline void __ll_loop_handle_scaling_generic (hyFloat& sum, hyFloat* _hprestrict_ parentConditionals, hyFloat* _hprestrict_ scalingAdjustments, long& didScale, long parentCode, long siteCount, long siteID, long& localScalerChange, long siteFrequency, long D) {
    if (__builtin_expect(sum < _lfScalingFactorThreshold && sum > 0.0,0)) {
        
        hyFloat scaler = _computeBoostScaler(scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards, sum, didScale);
        
        if (didScale) {
            #pragma unroll(8)
            #pragma GCC unroll 8
            for (long c = 0; c < D; c++) {
                parentConditionals [c] *= scaler;
            }
            
            if (siteFrequency == 1L) {
                localScalerChange += didScale;
            } else {
                localScalerChange += didScale * siteFrequency;
            }
            if (ADJUST)
                scalingAdjustments [parentCode*siteCount + siteID]  *= scaler;
        }
        
    } else {
        if (__builtin_expect(sum > _lfScalerUpwards,0)) {
            if (sum < HUGE_VAL) { // no point scaling an infinity
                
                hyFloat scaler = _computeReductionScaler (scalingAdjustments [parentCode*siteCount + siteID] * _lfScalingFactorThreshold, sum, didScale);
                
                if (didScale) {
                    #pragma unroll(8)
                    #pragma GCC unroll 8
                    for (long c = 0; c < D; c++) {
                         parentConditionals [c] *= scaler;
                     }
                    
                    if (siteFrequency == 1L) {
                        localScalerChange += didScale;
                    } else {
                        localScalerChange += didScale * siteFrequency;
                    }
                    if (ADJUST)
                        scalingAdjustments [parentCode*siteCount + siteID] *= scaler;
                }
            }
        }
    }
}

template<long D> inline void __ll_loop_handle_leaf_case (hyFloat* _hprestrict_ pp, hyFloat *  _hprestrict_ localScalingFactor , long siteFrom, long siteTo, _SimpleList&        siteOrdering, bool matchSet, long * _hprestrict_ setBranchTo) {
    if (matchSet) {
        memset (pp, 0, (siteTo-siteFrom) * sizeof (hyFloat));
        for (long k = siteFrom; k < siteTo; k++, pp += D) {
             pp[setBranchTo[siteOrdering.list_data[k]]] = localScalingFactor[k];
        }
    } else {
        for (long k = siteFrom; k < siteTo; k++, pp += D) {
            hyFloat lsf = localScalingFactor[k];
#pragma unroll(4)
#pragma GCC unroll 4
            for (long s = 0; s < D; s++) {
                pp[s] = lsf;
            }
        }
    }
}


inline void __ll_loop_handle_leaf_generic (hyFloat* _hprestrict_ pp, hyFloat *  _hprestrict_ localScalingFactor , long siteFrom, long siteTo, _SimpleList&        siteOrdering, bool matchSet, long * _hprestrict_ setBranchTo, long D) {
    
    if (matchSet) {
        memset (pp, 0, (siteTo-siteFrom) * sizeof (hyFloat));
        for (long k = siteFrom; k < siteTo; k++, pp += D) {
             pp[setBranchTo[siteOrdering.list_data[k]]] = localScalingFactor[k];
        }
    } else {
        if (D >= 4) {
            long DM4 = D-4;
            for (long k = siteFrom; k < siteTo; k++, pp += D) {
                hyFloat lsf = localScalingFactor[k];
                long s;
                
                for (s = 0; s < DM4; s+=4) {
                    pp[s] = lsf;
                    pp[s+1] = lsf;
                    pp[s+2] = lsf;
                    pp[s+3] = lsf;
                }
                
                for (; s < D; s++) {
                    pp[s] = lsf;
                }
            }
        } else {
            for (long k = siteFrom; k < siteTo; k++, pp += D) {
                hyFloat lsf = localScalingFactor[k];
                
                for (long s = 0; s < D; s++) {
                    pp[s] = lsf;
                }
            }
        }
    }
}


inline void _handle4x4_pruning_case (double const* childVector, double const* tMatrix, double* parentConditionals, void* transposed_mx) {
#ifdef _SLKP_USE_SSE_INTRINSICS
    double tv [4]     __attribute__ ((aligned (16))) = {childVector[0],
        childVector[1],
        childVector[2],
        childVector[3]};
    
    __m128d buffer0 = _mm_loadu_pd (tv),
    buffer1 = _mm_loadu_pd (tv+2),
    matrix01 = _mm_loadu_pd (tMatrix),
    matrix12 = _mm_loadu_pd (tMatrix+2),
    matrix34 = _mm_loadu_pd (tMatrix+4),
    matrix56 = _mm_loadu_pd (tMatrix+6),
    reg_storage  = _mm_mul_pd (buffer0, matrix01),
    reg_storage2 = _mm_mul_pd (buffer0, matrix34);
    
    matrix34     = _mm_mul_pd(buffer1, matrix12);
    matrix56     = _mm_mul_pd(buffer1, matrix56);
    reg_storage  = _mm_add_pd (reg_storage, matrix34);
    reg_storage2 = _mm_add_pd (reg_storage2, matrix56);
    reg_storage  = _mm_hadd_pd (reg_storage,reg_storage2);
    matrix01 = _mm_loadu_pd (parentConditionals);
    matrix01 = _mm_mul_pd (reg_storage, matrix01);
    _mm_storeu_pd (parentConditionals, matrix01);
    
    
    
    matrix01 = _mm_loadu_pd (tMatrix+8);
    matrix12 = _mm_loadu_pd (tMatrix+10);
    matrix34 = _mm_loadu_pd (tMatrix+12);
    matrix56 = _mm_loadu_pd (tMatrix+14);
    reg_storage  = _mm_mul_pd (buffer0, matrix01);
    reg_storage2 = _mm_mul_pd (buffer0, matrix34);
    
    matrix34     = _mm_mul_pd(buffer1, matrix12);
    matrix56     = _mm_mul_pd(buffer1, matrix56);
    reg_storage  = _mm_add_pd (reg_storage, matrix34);
    reg_storage2 = _mm_add_pd (reg_storage2, matrix56);
    reg_storage  = _mm_hadd_pd (reg_storage,reg_storage2);
    
    matrix01 = _mm_loadu_pd (parentConditionals+2);
    matrix01 = _mm_mul_pd (reg_storage, matrix01);
    _mm_storeu_pd (parentConditionals+2, matrix01);
    
    
    /*
     A1*B1 + A2*B2 + A3*B3 + A4*B4, where A4 = 1-A1-A2-A3 can be done with three multipications
     and 3 extra additions, like
     
     A1*(B1-B4) + A2*(B2-B4) + A3*(B3-B4) + B4

     20180914: SLKP, turning this off because of unstable numerical behavior if B1, B2, B3, B4 have
     very different magnitudes (that occurs for poorly initialized trees that require
     massive scaling)

     */
    
#elif defined _SLKP_USE_AVX_INTRINSICS
    
    __m256d c3     = _mm256_set1_pd(childVector[3]),
    c0     = _mm256_set1_pd(childVector[0]),
    c1     = _mm256_set1_pd(childVector[1]),
    c2     = _mm256_set1_pd(childVector[2]),
    t0,t1,t2,t3;
  
    if (transposed_mx) {
      t0    = ((__m256d*)transposed_mx)[0];
      t1    = ((__m256d*)transposed_mx)[1];
      t2    = ((__m256d*)transposed_mx)[2];
      t3    = ((__m256d*)transposed_mx)[3];
    } else {
      t0     = (__m256d) {tMatrix[0],tMatrix[4],tMatrix[8],tMatrix[12]};
      t1     = (__m256d) {tMatrix[1],tMatrix[5],tMatrix[9],tMatrix[13]};
      t2     = (__m256d) {tMatrix[2],tMatrix[6],tMatrix[10],tMatrix[14]};
      t3     = (__m256d) {tMatrix[3],tMatrix[7],tMatrix[11],tMatrix[15]};
    }
  
  // load transition matrix by column
#ifdef _SLKP_USE_FMA3_INTRINSICS
    __m256d sum01 = _mm256_fmadd_pd (c0, t0,_mm256_mul_pd(c1,t1)),
            sum23 = _mm256_fmadd_pd (c2,t2, _mm256_mul_pd(c3,t3));
    
    _mm256_storeu_pd(parentConditionals, _mm256_mul_pd (_mm256_loadu_pd (parentConditionals), _mm256_add_pd (sum01, sum23)));

#else
    __m256d sum01 = _mm256_add_pd (_mm256_mul_pd(c0,t0),_mm256_mul_pd(c1,t1)),
    sum23 = _mm256_add_pd (_mm256_mul_pd(c2,t2), _mm256_mul_pd(c3,t3));
    
    _mm256_storeu_pd(parentConditionals, _mm256_mul_pd (_mm256_loadu_pd (parentConditionals), _mm256_add_pd (sum01, sum23)));

#endif
    
    
    
#else
    // 12 multiplications, 16 additions, 3 subtractions
    
    /*hyFloat t1 = childVector[0] - childVector[3],
    t2 = childVector[1] - childVector[3],
    t3 = childVector[2] - childVector[3],
    t4 = childVector[3];
    
    parentConditionals [0] *= (tMatrix[0]  * t1 + tMatrix[1] * t2) + (tMatrix[2] * t3 + t4);
    parentConditionals [1] *= (tMatrix[4]  * t1 + tMatrix[5] * t2) + (tMatrix[6] * t3 + t4);
    parentConditionals [2] *= (tMatrix[8]  * t1 + tMatrix[9] * t2) + (tMatrix[10] * t3 + t4);
    parentConditionals [3] *= (tMatrix[12] * t1 + tMatrix[13] * t2) + (tMatrix[14] * t3 + t4);*/
  
    parentConditionals [0] *= tMatrix[0] * childVector[0] + tMatrix[1] * childVector[1] + tMatrix[2] * childVector[2] + tMatrix[3] * childVector[3];
    parentConditionals [1] *= tMatrix[4] * childVector[0] + tMatrix[5] * childVector[1] + tMatrix[6] * childVector[2] + tMatrix[7] * childVector[3];
    parentConditionals [2] *= tMatrix[8] * childVector[0] + tMatrix[9] * childVector[1] + tMatrix[10] * childVector[2] + tMatrix[11] * childVector[3];
    parentConditionals [3] *= tMatrix[12] * childVector[0] + tMatrix[13] * childVector[1] + tMatrix[14] * childVector[2] + tMatrix[15] * childVector[3];

#endif
    
}

inline void __ll_handle_tcc_init (_SimpleList const* __restrict tcc, bool isLeaf, long siteCount, long siteFrom, long nodeCode, long parentCode, long& parentTCCIBit, long& parentTCCIIndex, long &currentTCCBit, long& currentTCCIndex) {
    if (tcc) {
        parentTCCIIndex = siteCount * parentCode + siteFrom;
        parentTCCIBit   = parentTCCIIndex % _HY_BITMASK_WIDTH_;
        parentTCCIIndex = parentTCCIIndex / _HY_BITMASK_WIDTH_;
        if (! isLeaf) {
            currentTCCIndex = siteCount * nodeCode + siteFrom;
            currentTCCBit   = currentTCCIndex % _HY_BITMASK_WIDTH_;
            currentTCCIndex /= _HY_BITMASK_WIDTH_;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------*/

hyFloat      _TheTree::ComputeTreeBlockByBranch  (                   _SimpleList&        siteOrdering,
                                                  _SimpleList&        updateNodes,
                                                  _SimpleList* __restrict      tcc,
                                                  _DataSetFilter const*     theFilter,
                                                  hyFloat* __restrict         iNodeCache,
                                                  long      * __restrict        lNodeFlags,
                                                  hyFloat* __restrict        scalingAdjustments,
                                                  _Vector* __restrict    lNodeResolutions,
                                                  long&               overallScaler,
                                                  long                siteFrom,
                                                  long                siteTo,
                                                  long                catID,
                                                  hyFloat* __restrict         storageVec,
                                                  long* __restrict              siteCorrectionCounts,
                                                  long                setBranch,
                                                  long* __restrict              setBranchTo
                                                  )
// the updateNodes flags the nodes (leaves followed by inodes in the same order as flatLeaves and flatNodes)
// that must be recomputed
{
    // process the leaves first
    
    
    long * storage = (long*)alloca (sizeof (long) * flatNodes.lLength);
    InitializeArray(storage, flatNodes.lLength, 0L);
    _SimpleList     taggedInternals                 (flatNodes.lLength, storage);
    unsigned long   const alphabetDimension     =         theFilter->GetDimension(),
    siteCount           =         theFilter->GetPatternCount();
    
    _CalcNode       *currentTreeNode;
    long            localScalerChange     =         0;
    
    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }
    
    #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
        hyFloat * tMatrixT = nil;
        switch (alphabetDimension) {
            case 20:
                tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 20*20);
                break;
            case 60:
                tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 60*60);
                break;
            case 61:
                tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 61*61);
                break;
            case 62:
                tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 62*62);
                break;
            case 63:
                tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 63*63);
                break;
        }
    #endif
    
    for  (unsigned long nodeID = 0; nodeID < updateNodes.lLength; nodeID++) {
        long    nodeCode   = updateNodes.list_data [nodeID],
        parentCode = flatParents.list_data [nodeCode];
        hyFloat  *       childVector,
                 *       lastUpdatedSite;

        bool    isLeaf;
        if (nodeCode < flatLeaves.lLength) {
            isLeaf = true;
            currentTreeNode = ((_CalcNode*) flatCLeaves (nodeCode));
        } else {
            isLeaf = false;
            nodeCode -=  flatLeaves.lLength;
            childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
            currentTreeNode = ((_CalcNode*) flatTree    (nodeCode));
        }
        
        hyFloat  *  _hprestrict_ parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedInternals.list_data[parentCode] == 0) {
            // mark the parent for update and clear its conditionals if needed
            taggedInternals.list_data[parentCode]     = 1;
            hyFloat    *  _hprestrict_ localScalingFactor      = scalingAdjustments + parentCode*siteCount;
            
            bool    matchSet   = (parentCode == setBranch);
            
            switch (alphabetDimension) {
                case 4L:
                    __ll_loop_handle_leaf_case<4> (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo);
                    break;
                case 20L:
                    __ll_loop_handle_leaf_case<20> (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo);
                    break;
                case 60L:
                    __ll_loop_handle_leaf_case<60> (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo);
                    break;
                case 61L:
                    __ll_loop_handle_leaf_case<61> (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo);
                    break;
                case 62L:
                    __ll_loop_handle_leaf_case<62> (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo);
                    break;
                case 63L:
                    __ll_loop_handle_leaf_case<63> (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo);
                    break;
                default:
                    __ll_loop_handle_leaf_generic (parentConditionals, localScalingFactor , siteFrom, siteTo, siteOrdering, matchSet, setBranchTo, alphabetDimension);
            }
        }
        
        
        
        hyFloat  const * _hprestrict_ transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
        
        /*
         if (likeFuncEvalCallCount == 15098 && parentCode == 3688) {
            //if (currentTreeNode->GetName()->EndsWith("mt811400_SARS2_orf1ab_usa__4")) {
            fprintf (stderr, "\nBRANCH ID %ld (%ld = parent ID, parent name = %s) (%s)\n", nodeCode, parentCode, ((_CalcNode*)flatTree (parentCode))->GetName()->get_str(), currentTreeNode->GetName()->get_str());
                for (long e = 0; e < 16; e++) {
                    fprintf (stderr, "%ld => %lg\n", transitionMatrix[e]);
                }
            //}
        }
         */
        
        long currentTCCIndex,currentTCCBit,parentTCCIIndex,parentTCCIBit;
        __ll_handle_tcc_init (tcc, isLeaf, siteCount, siteFrom, nodeCode, parentCode, parentTCCIBit, parentTCCIIndex, currentTCCBit, currentTCCIndex);
        
        //long successiveSkips = 0;
        /**
          20200929: SLKP trying to refactor this code.
          handle specialized used cases
         */
        
// BEGIN NUCLEOTIDE CASE
        if (alphabetDimension == 4UL) {
            #ifdef _SLKP_USE_AVX_INTRINSICS
            __m256d tmatrix_transpose [4] = {
                (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
            };
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 4UL) {
                __ll_loop_preamble
                if (__ll_handle_conditional_array_initialization<4> (
                                                                      lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                    continue;
                }
                
                
                     
                #ifdef _SLKP_USE_AVX_INTRINSICS
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, tmatrix_transpose);
                #else
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, nil);
                #endif
                 
                sum     = (parentConditionals [0] + parentConditionals [1]) + (parentConditionals [2] + parentConditionals [3]);
                    
                __ll_loop_handle_scaling<4L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                
                childVector += 4L;
                __ll_loop_epilogue
            }
            
        }
// END NUCLEOTIDE CASE
        
// START AMINO-ACID CASE
        else if (alphabetDimension == 20UL) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                __ll_handle_matrix_transpose<20> (transitionMatrix, tMatrixT);
            #endif
            
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 20L) {
                __ll_loop_preamble
                if (__ll_handle_conditional_array_initialization<20> (
                    lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                    continue;
                }
                #if defined _SLKP_USE_AVX_INTRINSICS
                    sum = _avx_sum_4(__ll_handle_block20_product_sum<20,0> (tMatrixT, childVector, parentConditionals, nil));
                #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal = _mm_set1_pd(0.);
                    grandTotal = __ll_handle_block10_product_sum<20,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<20,10> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal;
                    sum = _sse_sum_2 (grandTotal);
                #else
                    __ll_product_sum_loop<20L> (tMatrix, childVector, parentConditionals, sum);
                #endif

                __ll_loop_handle_scaling<20L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                
                childVector += 20;
                __ll_loop_epilogue
            }

    }
    // END AMINO-ACID CASE
    // START UNIVERSAL CODE
    else if (alphabetDimension == 60UL) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                __ll_handle_matrix_transpose<60> (transitionMatrix, tMatrixT);
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 60L) {
                __ll_loop_preamble
                if (__ll_handle_conditional_array_initialization<60> (
                                                                      lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                    continue;
                }
                        
                
                #ifdef _SLKP_USE_AVX_INTRINSICS
                    __m256d grandTotal;
                    grandTotal = __ll_handle_block20_product_sum<60,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block20_product_sum<60,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block20_product_sum<60,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    sum += _avx_sum_4(grandTotal);
                #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal = _mm_set1_pd(0.);
                    grandTotal = __ll_handle_block10_product_sum<60,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<60,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,50> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    sum += _sse_sum_2(grandTotal);
                #else
                    __ll_product_sum_loop<60L> (tMatrix, childVector, parentConditionals, sum);
                #endif

                __ll_loop_handle_scaling<60L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

                childVector += 60;
                __ll_loop_epilogue
            }
    }
    else if (alphabetDimension == 61UL) {
        #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
            __ll_handle_matrix_transpose<61> (transitionMatrix, tMatrixT);
        #endif
        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 61L) {
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization<61> (
                                                                  lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                continue;
            }
                    
            
            #ifdef _SLKP_USE_AVX_INTRINSICS
            //Site 297 evaluated to a NaN probability in ComputeTreeBlockByBranch at branch 9698; this is not a recoverable error and indicates some serious COVFEFE taking place.
            
                /*if (siteID == 297 && parentCode == 9507) {
                    printf ("\nCONDITIONAL DUMP @ %s %ld (%ld parent)\n", currentTreeNode->GetName()->get_str(), nodeCode, parentCode);
                    for (int i = 0; i < 61; i++) {
                        printf ("%d %lg %lg\n", i, childVector[i], parentConditionals[i]);
                    }
                }*/

                __m256d grandTotal = __ll_handle_block20_product_sum<61,0> (tMatrixT, childVector, parentConditionals, nil);
                /*if (siteID ==  297 && parentCode == 9698) {
                    double checkGT [4];
                    _mm256_storeu_pd (checkGT, grandTotal);
                    printf ("\n%g %g %g %g", checkGT [0], checkGT [1], checkGT [2], checkGT [3]);
                }*/
                grandTotal = __ll_handle_block20_product_sum<61,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                /*if (siteID ==  297 && parentCode == 9698) {
                    double checkGT [4];
                    _mm256_storeu_pd (checkGT, grandTotal);
                    printf ("\n%g %g %g %g", checkGT [0], checkGT [1], checkGT [2], checkGT [3]);
                }*/
                grandTotal = __ll_handle_block20_product_sum<61,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                /*if (siteID ==  297 && parentCode == 9698) {
                    double checkGT [4];
                    _mm256_storeu_pd (checkGT, grandTotal);
                    printf ("\n%g %g %g %g", checkGT [0], checkGT [1], checkGT [2], checkGT [3]);
                }*/
                
                hyFloat const * __restrict tT = tMatrix + 61*60;
                __m256d lastTotal;
                __ll_handle_block20_product_sum_linear<61,0> (tT, childVector, &lastTotal);
                __ll_handle_block20_product_sum_linear<61,20> (tT, childVector, &lastTotal);
                __ll_handle_block20_product_sum_linear<61,40> (tT, childVector, &lastTotal);

                hyFloat s60 = _avx_sum_4(lastTotal) + childVector[60] * tT[60];
                parentConditionals[60] *= s60;
                sum += _avx_sum_4(grandTotal) + s60;
                /*if (siteID ==  297 && parentCode == 9698) {
                    printf ("\n%g => %g\n", s60, sum);
                }*/
            
                /*if (siteID == 297 && parentCode == 9507) {
                printf ("\nPRE-SCALE @ %s %ld\n", currentTreeNode->GetName()->get_str(), nodeCode);
                    for (int i = 0; i < 61; i++) {
                        printf ("%d %lg %lg\n", i, childVector[i], parentConditionals[i]);
                    }
                }
                for (int i = 0; i < 61; i++) {
                    if (isnan (parentConditionals[i])) {
                        HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch; this is not a recoverable error and indicates some serious COVFEFE taking place, node " & long(nodeCode) & "\n\n");
                    }
                }*/
                
            
            #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal;
                    grandTotal = __ll_handle_block10_product_sum<61,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<61,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,50> (tMatrixT, childVector, parentConditionals, &grandTotal);

                    hyFloat const * __restrict tT = tMatrix + 61*60;
                    __m128d lastTotal;
                    __ll_handle_block10_product_sum_linear<61,0> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,10> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,20> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,30> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,40> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,50> (tT, childVector, &lastTotal);
                    hyFloat s60 = _sse_sum_2(lastTotal) + childVector[60] * tT[60];
                    parentConditionals[60] *= s60;
                    sum += _sse_sum_2(grandTotal) + s60;
            #else
                __ll_product_sum_loop<61L> (tMatrix, childVector, parentConditionals, sum);
            #endif
            
            __ll_loop_handle_scaling<61L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
            /*if (siteID == 297 && parentCode == 9507) {
                printf ("\nPOST-SCALE @ %s %ld\n", currentTreeNode->GetName()->get_str(), nodeCode);
                for (int i = 0; i < 61; i++) {
                    printf ("%d %lg %lg\n", i, childVector[i], parentConditionals[i]);
                }
            }*/

            childVector += 61;
            __ll_loop_epilogue
        }
    } else if (alphabetDimension == 62UL) {
        #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
            __ll_handle_matrix_transpose<62> (transitionMatrix, tMatrixT);
        #endif
        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 62L) {
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization<62> (
                                                                  lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                continue;
            }
                    
            
            #if defined _SLKP_USE_AVX_INTRINSICS
                __m256d grandTotal = _mm256_set1_pd(0.);
                grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<62,0> (tMatrixT, childVector, parentConditionals, nil),grandTotal);
                grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<62,20> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<62,40> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                    
                
                hyFloat const * __restrict tT = tMatrix + 62*60;
                __m256d lastTotal;
                __ll_handle_block20_product_sum_linear<62,0> (tT, childVector, &lastTotal);
                __ll_handle_block20_product_sum_linear<62,20> (tT, childVector, &lastTotal);
                __ll_handle_block20_product_sum_linear<62,40> (tT, childVector, &lastTotal);
                hyFloat s60 = _avx_sum_4(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61];
                parentConditionals[60] *= s60;

                tT +=62;
                __m256d lastTotal2;
                __ll_handle_block20_product_sum_linear<62,0> (tT, childVector, &lastTotal2);
                __ll_handle_block20_product_sum_linear<62,20> (tT, childVector, &lastTotal2);
                __ll_handle_block20_product_sum_linear<62,40> (tT, childVector, &lastTotal2);

                hyFloat s61 = _avx_sum_4(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61];
                parentConditionals[61] *= s61;
                sum += _avx_sum_4(grandTotal) + s60 + s61;
            #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal;
                    grandTotal = __ll_handle_block10_product_sum<62,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<62,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,50> (tMatrixT, childVector, parentConditionals, &grandTotal);

                    hyFloat const * __restrict tT = tMatrix + 62*60;
                    __m128d lastTotal;
                    __ll_handle_block10_product_sum_linear<62,0> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,10> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,20> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,30> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,40> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,50> (tT, childVector, &lastTotal);
                    hyFloat s60 = _sse_sum_2(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61];
                    parentConditionals[60] *= s60;

                    tT += 62;
                    __m128d lastTotal2;
                    __ll_handle_block10_product_sum_linear<62,0> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,10> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,20> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,30> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,40> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,50> (tT, childVector, &lastTotal2);
            
                    hyFloat s61 = _sse_sum_2(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61];
                    parentConditionals[61] *= s61;

                    sum += _sse_sum_2(grandTotal) + s60 + s61;
            #else
                __ll_product_sum_loop<62L> (tMatrix, childVector, parentConditionals, sum);
            #endif

            __ll_loop_handle_scaling<62L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

            childVector += 62;
            __ll_loop_epilogue
        }
    } else if (alphabetDimension == 63UL) {
        #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
            __ll_handle_matrix_transpose<63> (transitionMatrix, tMatrixT);
        #endif
        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 63L) {
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization<63> (
                                                                  lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                continue;
            }
                    
            
            #if defined _SLKP_USE_AVX_INTRINSICS
                __m256d grandTotal = _mm256_set1_pd(0.);
                grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<63,0> (tMatrixT, childVector, parentConditionals, nil),grandTotal);
                grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<63,20> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<63,40> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                    
                
                hyFloat const * __restrict tT = tMatrix + 63*60;
                __m256d lastTotal;
                __ll_handle_block20_product_sum_linear<63,0> (tT, childVector, &lastTotal);
                __ll_handle_block20_product_sum_linear<63,20> (tT, childVector, &lastTotal);
                __ll_handle_block20_product_sum_linear<63,40> (tT, childVector, &lastTotal);
                hyFloat s60 = _avx_sum_4(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                parentConditionals[60] *= s60;

                tT += 63;
                __m256d lastTotal2;
                __ll_handle_block20_product_sum_linear<63,0> (tT, childVector, &lastTotal2);
                __ll_handle_block20_product_sum_linear<63,20> (tT, childVector, &lastTotal2);
                __ll_handle_block20_product_sum_linear<63,40> (tT, childVector, &lastTotal2);

                hyFloat s61 = _avx_sum_4(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                parentConditionals[61] *= s61;
 
                tT += 63;
                __m256d lastTotal3;
                __ll_handle_block20_product_sum_linear<63,0> (tT, childVector, &lastTotal3);
                __ll_handle_block20_product_sum_linear<63,20> (tT, childVector, &lastTotal3);
                __ll_handle_block20_product_sum_linear<63,40> (tT, childVector, &lastTotal3);

                hyFloat s62 = _avx_sum_4(lastTotal3) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                parentConditionals[62] *= s62;
                sum += _avx_sum_4(grandTotal) + s60 + s61 + s62;
            #elif defined _SLKP_USE_SSE_INTRINSICS
                __m128d grandTotal;
                grandTotal = __ll_handle_block10_product_sum<63,0> (tMatrixT, childVector, parentConditionals, nil);
                grandTotal = __ll_handle_block10_product_sum<63,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                grandTotal = __ll_handle_block10_product_sum<63,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                grandTotal = __ll_handle_block10_product_sum<63,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                grandTotal = __ll_handle_block10_product_sum<63,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                grandTotal = __ll_handle_block10_product_sum<63,50> (tMatrixT, childVector, parentConditionals, &grandTotal);

                hyFloat const * __restrict tT = tMatrix + 63*60;
                __m128d lastTotal;
                __ll_handle_block10_product_sum_linear<63,0> (tT, childVector, &lastTotal);
                __ll_handle_block10_product_sum_linear<63,10> (tT, childVector, &lastTotal);
                __ll_handle_block10_product_sum_linear<63,20> (tT, childVector, &lastTotal);
                __ll_handle_block10_product_sum_linear<63,30> (tT, childVector, &lastTotal);
                __ll_handle_block10_product_sum_linear<63,40> (tT, childVector, &lastTotal);
                __ll_handle_block10_product_sum_linear<63,50> (tT, childVector, &lastTotal);
                hyFloat s60 = _sse_sum_2(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                parentConditionals[60] *= s60;

                tT += 63;
                __m128d lastTotal2;
                __ll_handle_block10_product_sum_linear<63,0> (tT, childVector, &lastTotal2);
                __ll_handle_block10_product_sum_linear<63,10> (tT, childVector, &lastTotal2);
                __ll_handle_block10_product_sum_linear<63,20> (tT, childVector, &lastTotal2);
                __ll_handle_block10_product_sum_linear<63,30> (tT, childVector, &lastTotal2);
                __ll_handle_block10_product_sum_linear<63,40> (tT, childVector, &lastTotal2);
                __ll_handle_block10_product_sum_linear<63,50> (tT, childVector, &lastTotal2);
        
                hyFloat s61 = _sse_sum_2(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];;
                parentConditionals[61] *= s61;

                tT += 63;
                __m128d lastTotal3;
                __ll_handle_block10_product_sum_linear<63,0> (tT, childVector, &lastTotal3);
                __ll_handle_block10_product_sum_linear<63,10> (tT, childVector, &lastTotal3);
                __ll_handle_block10_product_sum_linear<63,20> (tT, childVector, &lastTotal3);
                __ll_handle_block10_product_sum_linear<63,30> (tT, childVector, &lastTotal3);
                __ll_handle_block10_product_sum_linear<63,40> (tT, childVector, &lastTotal3);
                __ll_handle_block10_product_sum_linear<63,50> (tT, childVector, &lastTotal3);
        
                hyFloat s62 = _sse_sum_2(lastTotal3) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];;
                parentConditionals[62] *= s62;

                sum += _sse_sum_2(grandTotal) + s60 + s61 + s62;
            #else
                __ll_product_sum_loop<63L> (tMatrix, childVector, parentConditionals, sum);
            #endif

            __ll_loop_handle_scaling<63L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

            childVector += 63;
            __ll_loop_epilogue
        }
    } else
    // END CODON CASES; GENERIC CASE
         for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += alphabetDimension) {
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization_generic (
                lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo, alphabetDimension)) {
                continue;
            }
            __ll_product_sum_loop_generic (tMatrix, childVector, parentConditionals, sum, alphabetDimension);

            __ll_loop_handle_scaling_generic <true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]),alphabetDimension);

            childVector += alphabetDimension;
            __ll_loop_epilogue
         }
    }
    
    // assemble the entire likelihood
    
    hyFloat * _hprestrict_ rootConditionals = iNodeCache + alphabetDimension * (siteFrom + (flatTree.lLength-1)  * siteCount);
    hyFloat                result = 0.0,
    correction = 0.0;
    
    
    for (long siteID = siteFrom, rootIndex = 0L; siteID < siteTo; siteID++) {
        hyFloat accumulator = 0.;
        
        if (setBranch == flatTree.lLength-1) {
            long                rootState = setBranchTo[siteOrdering.list_data[siteID]];
            accumulator         = rootConditionals[rootIndex + rootState] * theProbs[rootState];
            rootIndex           += alphabetDimension;
        } else {
            #pragma unroll(4)
            #pragma GCC unroll 4
            for (long p = 0; p < alphabetDimension; p++,rootIndex++) {
                accumulator += rootConditionals[rootIndex] * theProbs[p];
            }
        }
        /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
            fprintf (stderr, "\nREGULAR COMPUTE %lg (%ld) (%lg %lg %lg %lg)\n", accumulator, setBranch, rootConditionals[rootIndex-4],rootConditionals[rootIndex-3],rootConditionals[rootIndex-2],rootConditionals[rootIndex-1]);
        }*/
        if (storageVec) {
            storageVec [siteOrdering.list_data[siteID]] = accumulator;
        } else {
            if (accumulator <= 0.0) {
                result = -INFINITY;
#pragma omp critical
                {
                    //printf ("BAILING WITH INFINITY %ld\n", siteID);
                    hy_global::ReportWarning (_String("Site ") & (1L+siteOrdering.list_data[siteID]) & " evaluated to a 0 probability in ComputeTreeBlockByBranch");
                }
                break;
            }
            
            if (!isnan (accumulator)) {
                hyFloat term;
                long   const    site_frequency = theFilter->theFrequencies (siteOrdering.list_data[siteID]);
                
                if (site_frequency > 1L) {
                    term = log(accumulator) * site_frequency - correction;
                } else {
                    term = log(accumulator) - correction;
                }
                // Kahan sum
                

                hyFloat temp_sum = result + term;
                correction = (temp_sum - result) - term;
                result = temp_sum;
            } else {
                HandleApplicationError(_String("Site ") & (1L+siteOrdering.list_data[siteID]) & " evaluated to a NaN probability in ComputeTreeBlockByBranch; this is not a recoverable error and indicates some serious COVFEFE taking place.");
            }

        }
    }
    
    if (!storageVec && localScalerChange) {
#pragma omp atomic
        overallScaler += localScalerChange;
    }
    
    return result;
}

/*---------------------------------------------------------------------------------------------------*/
 
template<long D> inline bool __lcache_loop_preface (bool isLeaf, long* __restrict lNodeFlags, long siteID, _SimpleList const& siteOrdering, long nodeCode, long siteCount, long siteFrom, hyFloat* __restrict parentConditionals, hyFloat const* __restrict tMatrix, bool& canScale, hyFloat *& childVector, hyFloat *& lastUpdatedSite, _SimpleList const*  tcc, long&currentTCCBit, long& currentTCCIndex, long&parentTCCIBit, long& parentTCCIIndex, bool notPassedRoot, _Vector const*     lNodeResolutions) {
    if (isLeaf) {
        long siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.list_data[siteID]] ;
        if (siteState >= 0L) {
            unsigned long target_index = siteState;
            #pragma unroll(4)
            #pragma GCC unroll 4
            for (long k = 0L; k < D; k++, target_index+=D) {
                parentConditionals[k]   *= tMatrix[target_index];
            }
            return true;
        } else {
            childVector = lNodeResolutions->theData + (-siteState-1) * D;
        }
        canScale = false;
    } else {
        if (tcc&&notPassedRoot) {
            if ((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom)
                // the value of this conditional vector needs to be copied from a previously stored site
                // subtree duplication
                #pragma unroll(4)
                #pragma GCC unroll 4
                for (long k = 0UL; k < D; k++) {
                    childVector[k] = lastUpdatedSite[k];
                }
            else {
                lastUpdatedSite = childVector;
            }
            
            if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                currentTCCBit   = 0;
                currentTCCIndex ++;
            }
            if (++parentTCCIBit == _HY_BITMASK_WIDTH_) {
                parentTCCIBit   = 0;
                parentTCCIIndex ++;
            }
        }
    }
    return false;
}

/*---------------------------------------------------------------------------------------------------*/
 
inline bool __lcache_loop_preface_generic (bool isLeaf, long* __restrict lNodeFlags, long siteID, _SimpleList const& siteOrdering, long nodeCode, long siteCount, long siteFrom, hyFloat* __restrict parentConditionals, hyFloat const* __restrict tMatrix, bool& canScale, hyFloat *& childVector, hyFloat *& lastUpdatedSite, _SimpleList const*  tcc, long&currentTCCBit, long& currentTCCIndex, long&parentTCCIBit, long& parentTCCIIndex, bool notPassedRoot, _Vector const*     lNodeResolutions, long D) {
    if (isLeaf) {
        long siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.list_data[siteID]] ;
        if (siteState >= 0L) {
            unsigned long target_index = siteState;
            #pragma unroll(4)
            #pragma GCC unroll 4
            for (long k = 0L; k < D; k++, target_index+=D) {
                parentConditionals[k]   *= tMatrix[target_index];
            }
            
            return true;
        } else {
            childVector = lNodeResolutions->theData + (-siteState-1) * D;
        }
        canScale = false;
    } else {
        if (tcc&&notPassedRoot) {
            if ((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom)
                // the value of this conditional vector needs to be copied from a previously stored site
                // subtree duplication
                #pragma unroll(4)
                #pragma GCC unroll 4
                for (long k = 0UL; k < D; k++) {
                    childVector[k] = lastUpdatedSite[k];
                }
            else {
                lastUpdatedSite = childVector;
            }
            
            if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                currentTCCBit   = 0;
                currentTCCIndex ++;
            }
            if (++parentTCCIBit == _HY_BITMASK_WIDTH_) {
                parentTCCIBit   = 0;
                parentTCCIIndex ++;
            }
        }
    }
    return false;
}
 
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/

void            _TheTree::ComputeBranchCache    (
                                                 _SimpleList&            siteOrdering,
                                                 long                    brID,
                                                 hyFloat*   __restrict      cache,
                                                 hyFloat*   __restrict      iNodeCache,
                                                 _DataSetFilter const*     theFilter,
                                                 long           * __restrict       lNodeFlags,
                                                 hyFloat*  __restrict       scalingAdjustments,
                                                 long        *  __restrict         siteCorrectionCounts,
                                                 _Vector const*     lNodeResolutions,
                                                 long&                   overallScaler,
                                                 long const                   siteFrom,
                                                 long                    siteTo,
                                                 long const                  catID,
                                                 _SimpleList const*            tcc,
                                                 hyFloat* __restrict        siteRes
                                                 )
{
    
    
    /*
     the cache matrix (linearized into a vector) will have TWO rows with siteCount blocks of alphabetDimension doubles, storing the conditional likelihoods of individual sites at a given branch
     in the virtually rerooted tree
     
     cache ->
     Row 0 [brID node -- the branch that is being rerooted on]
     Row 1 [conditional likelihoods for the new root]
     */
    
    
    auto __handle_site_corrections = [&] (long didScale, long siteID)->void {
        if (didScale&&siteCorrectionCounts) {
            //
            siteCorrectionCounts [siteOrdering.list_data[siteID]] += didScale;
            if (didScale == 1L) {
                siteRes[siteOrdering.list_data[siteID]] *= _lfScalerUpwards;
            } else {
                if (didScale == -1L) {
                    siteRes[siteOrdering.list_data[siteID]] *= _lfScalingFactorThreshold;
                } else {
                    if (didScale > 0) {
                        for (long k = 0; k < didScale; k++) {
                            siteRes[siteOrdering.list_data[siteID]] *= _lfScalerUpwards;
                        }
                    } else{
                        for (long k = 0; k < -didScale; k++) {
                             siteRes[siteOrdering.list_data[siteID]] *= _lfScalingFactorThreshold;
                         }
                    }
                }
            }
        }
    };
    
    long *tagged_node_cache = (long*)alloca ((flatLeaves.lLength + flatNodes.lLength)*sizeof (long));
    
    _SimpleList taggedNodes (flatLeaves.lLength + flatNodes.lLength, tagged_node_cache),
                nodesToProcess,
                rootPath;
    
    taggedNodes.Populate (flatLeaves.lLength + flatNodes.lLength, 0, 0);
    
     
    long        myParent               = brID       -flatLeaves.lLength;
    
    const long  alphabetDimension     =            theFilter->GetDimension(),
    alphabetDimensionmod4  =         alphabetDimension - alphabetDimension % 4,
    siteCount               =            theFilter->GetPatternCount();

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }
    
    do {
        taggedNodes.list_data[myParent+flatLeaves.lLength] = 1;
        myParent = flatParents.list_data[myParent+flatLeaves.lLength];
    } while (myParent >= 0);

#if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
    hyFloat * tMatrixT = nil;
    switch (alphabetDimension) {
        case 20:
            tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 20*20);
            break;
        case 60:
            tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 60*60);
            break;
        case 61:
            tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 61*61);
            break;
        case 62:
            tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 62*62);
            break;
        case 63:
            tMatrixT = (hyFloat*)alloca (sizeof (hyFloat) * 63*63);
            break;
    }
#endif

    
    for (unsigned long k = 0UL; k <  flatLeaves.lLength+flatNodes.lLength; k++) {
        myParent = flatParents.list_data[k];
        if (taggedNodes.list_data[myParent+flatLeaves.lLength] == 1 && taggedNodes.list_data[k] == 0) {
            if (myParent != brID - flatLeaves.lLength) {
                nodesToProcess << k;
            }
        }
        if (taggedNodes.list_data[k]) {
            rootPath << k;
        }
    }
    
    
    hyFloat * state = cache + alphabetDimension * siteFrom,
    * childVector;
    
    long        localScalerChange = 0;
    
    // first populate the downward looking vector of conditionals
    
    if (brID < flatLeaves.lLength) { // a leaf
        if (alphabetDimension == 4) {
            for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += 4) {
                long siteState = lNodeFlags[brID*siteCount + siteOrdering.list_data[siteID]] ;
                if (siteState >= 0) {
                    state[0] = 0.; state [1] = 0.; state [2] = 0.; state [3] = 0.;
                    state[siteState] = 1.;
                } else {
                    childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                    state [0] = childVector [0]; state [1] = childVector [1]; state [2] = childVector [2]; state [3] = childVector [3];
                }
            }

        } else {
            for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension) {
                long siteState = lNodeFlags[brID*siteCount + siteOrdering.list_data[siteID]] ;
                if (siteState >= 0) {
                    // a single character state; sweep down the appropriate column
                    memset (state, 0, sizeof (hyFloat) * alphabetDimension);
                    state[siteState] = 1.;
                } else {
                    childVector = lNodeResolutions->theData + (-siteState-1) * alphabetDimension;
                    memcpy (state, childVector, sizeof (hyFloat) * alphabetDimension);
                }
            }
        }
        
    } else { // an internal branch
        long        nodeCode = brID - flatLeaves.lLength;
        hyFloat *lastUpdated = iNodeCache + (nodeCode * siteCount + siteFrom) * alphabetDimension;
        
        long currentTCCIndex        ,
        currentTCCBit            ;
        
        if (tcc) {
            currentTCCIndex = siteCount * nodeCode + siteFrom;
            currentTCCBit   = currentTCCIndex % _HY_BITMASK_WIDTH_;
            currentTCCIndex /= _HY_BITMASK_WIDTH_;
        }
        
        for (long siteID = siteFrom; siteID < siteTo; siteID ++, state += alphabetDimension) {
            if (tcc) {
                if ((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) == 0) {
                    lastUpdated = iNodeCache + (nodeCode * siteCount + siteID) * alphabetDimension;
                }
            }
            
            for (long s = 0; s < alphabetDimension; s++) {
                state[s] = lastUpdated[s];
            }
            
            if (tcc) {
                if (++currentTCCBit == _HY_BITMASK_WIDTH_) {
                    currentTCCBit   = 0;
                    currentTCCIndex ++;
                }
            } else {
                lastUpdated += alphabetDimension;
            }
        }
    }
    
    taggedNodes.Populate (flatTree.lLength, 0, 0);
    rootPath.Flip ();
    
    /*if (likeFuncEvalCallCount == 15098) {
        fprintf (stderr, "\nCaching branch %ld\n", brID);
        rootPath.Each ([&](long n, unsigned long i) -> void {
            _CalcNode * current_node = (_CalcNode*) (n < flatLeaves.lLength ?  flatCLeaves (n):
                                                                              flatTree    (n-flatLeaves.lLength));
            fprintf (stderr, "[%d: %d] %s\n", i, n,current_node->GetName()->get_str());
        });
        ObjectToConsole(&nodesToProcess);NLToConsole();
    }*/
    
    long const node_count = nodesToProcess.lLength + rootPath.lLength - 2L;
    
    for  (long nodeID = 0; nodeID < node_count; nodeID++) {
        bool    notPassedRoot = nodeID<nodesToProcess.lLength;
        
        long nodeCode, parentCode;
        
        if (notPassedRoot) {
            nodeCode = nodesToProcess.list_data [nodeID];
            parentCode = flatParents.list_data [nodeCode];
        } else {
            nodeCode = rootPath.list_data[nodeID-nodesToProcess.lLength];
            parentCode = (rootPath.list_data[nodeID-nodesToProcess.lLength+1] - flatLeaves.lLength);
        }
        
        bool    isLeaf     = nodeCode < flatLeaves.lLength;
        
        if (!isLeaf) {
            nodeCode -=  flatLeaves.lLength;
        }
        
        hyFloat * parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedNodes.list_data[parentCode] == 0L) {
            // mark the parent for update and clear its conditionals if needed
            //printf ("Resetting parentCode = %ld\n", parentCode);
            taggedNodes.list_data[parentCode]     = 1L;
            hyFloat     const *localScalingFactor      = scalingAdjustments + parentCode*siteCount;
            if (alphabetDimension == 4L) {
                unsigned long k3     = 0UL;
                for (long k = siteFrom; k < siteTo; k++, k3+=4) {
                    hyFloat scaler = localScalingFactor[k];
                    parentConditionals [k3]   = scaler;
                    parentConditionals [k3+1UL] = scaler;
                    parentConditionals [k3+2UL] = scaler;
                    parentConditionals [k3+3UL] = scaler;
                }
            } else {
                unsigned long k3     = 0UL;
                for (unsigned long k = siteFrom; k < siteTo; k++) {
                    hyFloat scaler = localScalingFactor[k];
                    #pragma unroll(4)
                    #pragma GCC unroll 4
                    for (unsigned long k2 = 0UL; k2 < alphabetDimension; k2++, k3++) {
                        parentConditionals [k3] = scaler;
                    }
                }
            }
        }
        
        _CalcNode    * currentTreeNode = (_CalcNode*) (isLeaf?  flatCLeaves (nodeCode):
                                                       flatTree    (notPassedRoot?nodeCode:parentCode));
        
        /*if (likeFuncEvalCallCount == 15098) {
            fprintf (stderr, "%ld/%ld (%ld/%ld) => %s\n", nodeID, nodeCode, isLeaf, notPassedRoot, currentTreeNode->GetName()->get_str());
        }*/
        
        hyFloat  const *  transitionMatrix = currentTreeNode->GetCompExp(catID)->theData;
        hyFloat  *       childVector,*     lastUpdatedSite;
        
        if (!isLeaf) {
            lastUpdatedSite = childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
        }
        
        
        long currentTCCIndex,currentTCCBit,parentTCCIIndex,parentTCCIBit;
        __ll_handle_tcc_init (tcc, isLeaf, siteCount, siteFrom, nodeCode, parentCode, parentTCCIBit, parentTCCIIndex, currentTCCBit, currentTCCIndex);

 
        if (alphabetDimension == 4L) {
            #ifdef _SLKP_USE_AVX_INTRINSICS
                __m256d tmatrix_transpose [4] = {
                    (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                    (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                    (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                    (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                };
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 4L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<4>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                        fprintf (stderr, "__lcache_loop_preface (%ld) %g %g %g %g\n", nodeCode, parentConditionals[0], parentConditionals[1], parentConditionals[2], parentConditionals[3]);
                    }*/
                    continue;
                }
                long     didScale =  0;
                #ifdef _SLKP_USE_AVX_INTRINSICS
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, tmatrix_transpose);
                #else
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, nil);
                #endif
                if (canScale) {
                    hyFloat sum     = (parentConditionals [0] + parentConditionals [1]) + (parentConditionals [2] + parentConditionals [3]);
                    
                    __ll_loop_handle_scaling<4L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

                }
                /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                    fprintf (stderr, "NODE = %ld, PARENT = %ld (%ld), P(G) = %lg, P(T) = %lg, scale = %ld\n", nodeCode, parentCode, canScale, parentConditionals[2], parentConditionals[3], didScale);
                }*/
                childVector += 4L;
                __handle_site_corrections(didScale, siteID);
            }
        } else if (alphabetDimension == 20L) {
                #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                    __ll_handle_matrix_transpose<20>(transitionMatrix, tMatrixT);
                #endif
                for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 20L) {
                    bool canScale = !notPassedRoot;
                    hyFloat  const *tMatrix = transitionMatrix;
                    if (__lcache_loop_preface<20>(
                    isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                        continue;
                    }
                    long     didScale =  0;
                    hyFloat sum     = 0.;
                    #if defined _SLKP_USE_AVX_INTRINSICS
                        __m256d s256 = __ll_handle_block20_product_sum<20,0> (tMatrixT, childVector, parentConditionals, nil);
                    #elif  defined _SLKP_USE_SSE_INTRINSICS
                        __m128d s128;
                        s128 =  __ll_handle_block10_product_sum<20,0> (tMatrixT, childVector, parentConditionals, nil);
                        s128 =  __ll_handle_block10_product_sum<20,10> (tMatrixT, childVector, parentConditionals, &s128);
                    #else
                        __ll_product_sum_loop<20L> (tMatrix, childVector, parentConditionals, sum);
                    #endif
                    if (canScale) {
                        #if defined _SLKP_USE_AVX_INTRINSICS
                            sum = _avx_sum_4(s256);
                        #elif defined _SLKP_USE_SSE_INTRINSICS
                            sum = _sse_sum_2(s128);
                        #endif
                        __ll_loop_handle_scaling<20L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                    }
                    childVector += 20L;
                    __handle_site_corrections(didScale, siteID);
                }
        } else if (alphabetDimension == 60L) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                __ll_handle_matrix_transpose<60L>(transitionMatrix, tMatrixT);
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 60L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<60>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    continue;
                }
                long     didScale =  0;
                hyFloat sum     = 0.;
                #if defined _SLKP_USE_AVX_INTRINSICS
                    __m256d grandTotal = _mm256_set1_pd(0.);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<60,0> (tMatrixT, childVector, parentConditionals, nil),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<60,20> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<60,40> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal;
                    grandTotal = __ll_handle_block10_product_sum<60,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<60,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<60,50> (tMatrixT, childVector, parentConditionals, &grandTotal);
                #else
                    __ll_product_sum_loop<60L> (tMatrix, childVector, parentConditionals, sum);
                #endif
                if (canScale) {
                    #if defined _SLKP_USE_AVX_INTRINSICS
                        sum = _avx_sum_4(grandTotal);
                    #elif defined _SLKP_USE_SSE_INTRINSICS
                        sum = _sse_sum_2(grandTotal);
                    #endif
                    __ll_loop_handle_scaling<60L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 60L;
                __handle_site_corrections(didScale, siteID);
            }
        } else if (alphabetDimension == 61L) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                __ll_handle_matrix_transpose<61L>(transitionMatrix, tMatrixT);
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 61L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<61>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    continue;
                }
                long     didScale =  0;
                hyFloat sum     = 0.;
                #if defined _SLKP_USE_AVX_INTRINSICS
                    __m256d grandTotal = _mm256_set1_pd(0.);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<61,0> (tMatrixT, childVector, parentConditionals, nil),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<61,20> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<61,40> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                        
                    hyFloat const * __restrict tT = tMatrix + 61*60;
                    __m256d lastTotal;
                    __ll_handle_block20_product_sum_linear<61,0> (tT, childVector, &lastTotal);
                    __ll_handle_block20_product_sum_linear<61,20> (tT, childVector, &lastTotal);
                    __ll_handle_block20_product_sum_linear<61,40> (tT, childVector, &lastTotal);

                    hyFloat s60 = _avx_sum_4(lastTotal) + childVector[60] * tT[60];
                    parentConditionals[60] *= s60;
                #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal;
                    grandTotal = __ll_handle_block10_product_sum<61,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<61,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<61,50> (tMatrixT, childVector, parentConditionals, &grandTotal);

                    hyFloat const * __restrict tT = tMatrix + 61*60;
                    __m128d lastTotal;
                    __ll_handle_block10_product_sum_linear<61,0> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,10> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,20> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,30> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,40> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<61,50> (tT, childVector, &lastTotal);
                    hyFloat s60 = _sse_sum_2(lastTotal) + childVector[60] * tT[60];
                    parentConditionals[60] *= s60;
                #else
                    __ll_product_sum_loop<61L> (tMatrix, childVector, parentConditionals, sum);
                #endif
                if (canScale) {
                    bool canScale = !notPassedRoot;
                    #if defined _SLKP_USE_AVX_INTRINSICS
                        sum = _avx_sum_4(grandTotal) + s60;
                    #elif defined _SLKP_USE_SSE_INTRINSICS
                        sum = _sse_sum_2(grandTotal) + s60;
                    #endif
                    __ll_loop_handle_scaling<61L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 61L;
                __handle_site_corrections(didScale, siteID);
            }
        }  else if (alphabetDimension == 62L) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                __ll_handle_matrix_transpose<62L>(transitionMatrix, tMatrixT);
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 62L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<62>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    continue;
                }
                long     didScale =  0;
                hyFloat sum     = 0.;
                #if defined _SLKP_USE_AVX_INTRINSICS
                    __m256d grandTotal = _mm256_set1_pd(0.);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<62,0> (tMatrixT, childVector, parentConditionals, nil),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<62,20> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<62,40> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                        
                    hyFloat const * __restrict tT = tMatrix + 62*60;
                    __m256d lastTotal;
                    __ll_handle_block20_product_sum_linear<62,0> (tT, childVector, &lastTotal);
                    __ll_handle_block20_product_sum_linear<62,20> (tT, childVector, &lastTotal);
                    __ll_handle_block20_product_sum_linear<62,40> (tT, childVector, &lastTotal);
                    hyFloat s60 = _avx_sum_4(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61];
                    parentConditionals[60] *= s60;

                    tT += 62;
                    __m256d lastTotal2;
                    __ll_handle_block20_product_sum_linear<62,0> (tT, childVector, &lastTotal2);
                    __ll_handle_block20_product_sum_linear<62,20> (tT, childVector, &lastTotal2);
                    __ll_handle_block20_product_sum_linear<62,40> (tT, childVector, &lastTotal2);

                    hyFloat s61 = _avx_sum_4(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61];
                    parentConditionals[61] *= s61;
                
                #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal;
                    grandTotal = __ll_handle_block10_product_sum<62,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<62,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<62,50> (tMatrixT, childVector, parentConditionals, &grandTotal);

                    hyFloat const * __restrict tT = tMatrix + 62*60;
                    __m128d lastTotal;
                    __ll_handle_block10_product_sum_linear<62,0> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,10> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,20> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,30> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,40> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<62,50> (tT, childVector, &lastTotal);
                    hyFloat s60 = _sse_sum_2(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61];
                    parentConditionals[60] *= s60;
                
                    tT += 62;
                    __m128d lastTotal2;
                    __ll_handle_block10_product_sum_linear<62,0> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,10> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,20> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,30> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,40> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<62,50> (tT, childVector, &lastTotal2);
                    hyFloat s61 = _sse_sum_2(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61];
                    parentConditionals[61] *= s61;
                #else
                    __ll_product_sum_loop<62L> (tMatrix, childVector, parentConditionals, sum);
                #endif
                if (canScale) {
                    #if defined _SLKP_USE_AVX_INTRINSICS
                        sum = _avx_sum_4(grandTotal) + s60 + s61;
                    #elif defined _SLKP_USE_SSE_INTRINSICS
                        sum = _sse_sum_2(grandTotal) + s60 + s61;
                    #endif
                    __ll_loop_handle_scaling<62L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 62L;
                __handle_site_corrections(didScale, siteID);
            }
        } else if (alphabetDimension == 63L) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS
                __ll_handle_matrix_transpose<63L>(transitionMatrix, tMatrixT);
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 63L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<63>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    continue;
                }
                long     didScale =  0;
                hyFloat sum     = 0.;
                #if defined _SLKP_USE_AVX_INTRINSICS
                    __m256d grandTotal = _mm256_set1_pd(0.);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<63,0> (tMatrixT, childVector, parentConditionals, nil),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<63,20> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                    grandTotal = _mm256_add_pd (__ll_handle_block20_product_sum<63,40> (tMatrixT, childVector, parentConditionals, &grandTotal),grandTotal);
                        
                    hyFloat const * __restrict tT = tMatrix + 63*60;
                    __m256d lastTotal;
                    __ll_handle_block20_product_sum_linear<63,0> (tT, childVector, &lastTotal);
                    __ll_handle_block20_product_sum_linear<63,20> (tT, childVector, &lastTotal);
                    __ll_handle_block20_product_sum_linear<63,40> (tT, childVector, &lastTotal);
                    hyFloat s60 = childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62] + _avx_sum_4(lastTotal);
                    parentConditionals[60] *= s60;

                    tT += 63;
                    __m256d lastTotal2;
                    __ll_handle_block20_product_sum_linear<63,0> (tT, childVector, &lastTotal2);
                    __ll_handle_block20_product_sum_linear<63,20> (tT, childVector, &lastTotal2);
                    __ll_handle_block20_product_sum_linear<63,40> (tT, childVector, &lastTotal2);

                    hyFloat s61 = childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62] + _avx_sum_4(lastTotal2);
                    parentConditionals[61] *= s61;
     
                    tT += 63;
                    __m256d lastTotal3;
                    __ll_handle_block20_product_sum_linear<63,0> (tT, childVector, &lastTotal3);
                    __ll_handle_block20_product_sum_linear<63,20> (tT, childVector, &lastTotal3);
                    __ll_handle_block20_product_sum_linear<63,40> (tT, childVector, &lastTotal3);

                    hyFloat s62 = childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62] + _avx_sum_4(lastTotal3);
                    parentConditionals[62] *= s62;
                #elif defined _SLKP_USE_SSE_INTRINSICS
                    __m128d grandTotal;
                    grandTotal = __ll_handle_block10_product_sum<63,0> (tMatrixT, childVector, parentConditionals, nil);
                    grandTotal = __ll_handle_block10_product_sum<63,10> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<63,20> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<63,30> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<63,40> (tMatrixT, childVector, parentConditionals, &grandTotal);
                    grandTotal = __ll_handle_block10_product_sum<63,50> (tMatrixT, childVector, parentConditionals, &grandTotal);

                    hyFloat const * __restrict tT = tMatrix + 63*60;
                    __m128d lastTotal;
                    __ll_handle_block10_product_sum_linear<63,0> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<63,10> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<63,20> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<63,30> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<63,40> (tT, childVector, &lastTotal);
                    __ll_handle_block10_product_sum_linear<63,50> (tT, childVector, &lastTotal);
                    hyFloat s60 = _sse_sum_2(lastTotal) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                    parentConditionals[60] *= s60;
                
                    tT += 63;
                    __m128d lastTotal2;
                    __ll_handle_block10_product_sum_linear<63,0> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<63,10> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<63,20> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<63,30> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<63,40> (tT, childVector, &lastTotal2);
                    __ll_handle_block10_product_sum_linear<63,50> (tT, childVector, &lastTotal2);
                    hyFloat s61 = _sse_sum_2(lastTotal2) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                    parentConditionals[61] *= s61;

                    tT += 63;
                    __m128d lastTotal3;
                    __ll_handle_block10_product_sum_linear<63,0> (tT, childVector, &lastTotal3);
                    __ll_handle_block10_product_sum_linear<63,10> (tT, childVector, &lastTotal3);
                    __ll_handle_block10_product_sum_linear<63,20> (tT, childVector, &lastTotal3);
                    __ll_handle_block10_product_sum_linear<63,30> (tT, childVector, &lastTotal3);
                    __ll_handle_block10_product_sum_linear<63,40> (tT, childVector, &lastTotal3);
                    __ll_handle_block10_product_sum_linear<63,50> (tT, childVector, &lastTotal3);
                    hyFloat s62 = _sse_sum_2(lastTotal3) + childVector[60] * tT[60] + childVector[61] * tT[61] + childVector[62] * tT[62];
                    parentConditionals[62] *= s62;
                #else
                    __ll_product_sum_loop<63L> (tMatrix, childVector, parentConditionals, sum);
                #endif
                if (canScale) {
                    #if defined _SLKP_USE_AVX_INTRINSICS
                        sum = _avx_sum_4(grandTotal) + s60 + s61 + s62;
                    #elif defined _SLKP_USE_SSE_INTRINSICS
                        sum = _sse_sum_2(grandTotal) + s60 + s61 + s62;
                    #endif
                    __ll_loop_handle_scaling<63L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 63L;
                __handle_site_corrections(didScale, siteID);
            }
        } else  {
              
              for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += alphabetDimension) {
                bool canScale = !notPassedRoot;

                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface_generic(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions, alphabetDimension)) {
                    continue;
                }
                long     didScale =  0;
                hyFloat sum     = 0.;
                __ll_product_sum_loop_generic (tMatrix, childVector, parentConditionals, sum, alphabetDimension);
                if (canScale) {
                    //printf ("scale generic\n");
                    #pragma GCC unroll 8
                    #pragma clang loop vectorize(enable)
                    #pragma clang loop interleave(enable)
                    #pragma clang loop unroll(enable)
                    for (long k = 0; k < alphabetDimension; k++) {
                        sum += parentConditionals[k];
                    }
                    __ll_loop_handle_scaling_generic<false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]), alphabetDimension);
                }
                childVector += alphabetDimension;
                __handle_site_corrections(didScale, siteID);
            }
        }
    }
    
    
    
    //printf ("root name %s\n", ((_CalcNode    *)flatTree(rootPath.list_data[rootPath.lLength-2] - flatLeaves.lLength))->GetName()->sData);
    
    hyFloat const *rootConditionals   = iNodeCache +  (rootPath.list_data[rootPath.lLength-2] - flatLeaves.lLength)  * siteCount * alphabetDimension;
    
    state = cache + alphabetDimension * siteCount;
    const unsigned long site_bound = alphabetDimension*siteTo;
    for (unsigned long ii = siteFrom * alphabetDimension; ii < site_bound; ii++) {
        state[ii] = rootConditionals[ii];
        /*if (likeFuncEvalCallCount == 15098 && ii / alphabetDimension == 91) {
            printf ("Site %ld, Root conditional [%ld] = %g, node state [%ld] = %g\n", ii/alphabetDimension, ii, state[ii], ii, cache[ii]);
        }*/
    }
    
    if (!siteCorrectionCounts && localScalerChange) {
#pragma omp atomic
        overallScaler += localScalerChange;
        
        //#pragma omp atomic
        // printf ("Rescale in ComputeBranchCache at branch %ld %ld\n", brID, localScalerChange);
    }
}
