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


#if not defined _SLKP_USE_ARM_NEON && not defined _SLKP_USE_SSE_INTRINSICS && not defined _SLKP_USE_AVX_INTRINSICS
    #define _HY_NO_INTRINSICS_MARKER
#endif

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

#if defined _SLKP_USE_SSE_INTRINSICS  or defined _SLKP_USE_AVX_INTRINSICS
inline double _sse_sum_2 (__m128d const & x) {
    return _mm_cvtsd_f64(_mm_hadd_pd(x, x));
}
#endif

#ifdef _SLKP_USE_APPLE_BLAS_2
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
  AtlasConj=114};
extern "C" void cblas_dgemv(const enum CBLAS_ORDER __Order,
                 const enum CBLAS_TRANSPOSE __TransA,
                 const int __M, const int __N,
                 const double __alpha, const double *__A,
                 const int __lda, const double *__X, const int __incX,
                 const double __beta, double *__Y, const int __incY);
#endif

/*
template<long D> inline void __ll_handle_matrix_transpose (hyFloat const * __restrict transitionMatrix, hyFloat * __restrict tMatrixT) {
    long i = 0L;
    for (long r = 0L; r < D; r++) {
        //#pragma unroll(4)
        #pragma GCC unroll 4
        for (long c = 0L; c < D; c++, i++) {
            tMatrixT[c*D+r] = transitionMatrix[i];
        }
    }
}
*/

template<long D> inline bool __ll_handle_conditional_array_initialization ( long * __restrict lNodeFlags, bool isLeaf, long nodeCode, long setBranch, long iNodes, long siteID, long siteFrom, long siteCount, _SimpleList&            siteOrdering, hyFloat * __restrict parentConditionals, hyFloat const * __restrict tMatrix, _Vector const*     lNodeResolutions, hyFloat *& childVector, _SimpleList const* tcc, long &currentTCCBit, long& currentTCCIndex, hyFloat * & lastUpdatedSite, long* __restrict  setBranchTo) {
    

    if (isLeaf) {
        long siteState;
        if (setBranch != nodeCode + iNodes) {
            siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.list_data[siteID]] ;
            //if (siteID == 524) {
            //    fprintf (stderr, "Site-state @ %ld: %d\n", nodeCode, siteState);
            //}
        } else {
            siteState = setBranchTo[siteOrdering.list_data[siteID]] ;
        }
        if (__builtin_expect(siteState >= 0L,1)) {
            // a single character state; sweep down the appropriate column
#ifdef _SLKP_USE_ARM_NEON
            
            const long ub = (D>>2<<2);
            
            for (long k = 0; k < ub; k+=4) {
                float64x2x2_t PC = vld1q_f64_x2 (parentConditionals+k),
                              CC;
                
                CC.val[0] = vld1q_lane_f64 (tMatrix + siteState + D*k,      CC.val[0],0);
                CC.val[1] = vld1q_lane_f64 (tMatrix + siteState + D*(k+2),  CC.val[1],0);
                
                CC.val[0] = vld1q_lane_f64 (tMatrix + siteState + D*(k+1),  CC.val[0],1);
                CC.val[1] = vld1q_lane_f64 (tMatrix + siteState + D*(k+3),  CC.val[1],1);

                PC.val[0] = vmulq_f64(PC.val[0], CC.val[0]);
                PC.val[1] = vmulq_f64(PC.val[1], CC.val[1]);

                vst1q_f64 (parentConditionals+k, PC.val[0]);
                vst1q_f64 (parentConditionals+k+2, PC.val[1]);

                
            }
            for (long k = ub; k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState+D*k];
            }
            
            /*
            if (parentConditionals[siteState] == 0.) {
                if (siteID == 524) {
                    fprintf (stderr, "%ld/%ld / %g\n", siteState, siteID, parentConditionals[siteState]);
                }
            }
            */
#else
            #pragma GCC unroll 4
            for (long k = 0L; k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState+D*k];
            }
#endif
            return true;
        } else {
            childVector = lNodeResolutions->theData + (-siteState-1) * D;
        }
    } else {
        if (tcc) {
            if (__builtin_expect((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom,0)) {
                #pragma GCC unroll 4
                for (long k = 0L; k < D; k++) {
                    childVector[k] = lastUpdatedSite[k];
                }
                //memcpy (childVector, lastUpdatedSite, sizeof (hyFloat)*D);
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

template<long D> inline bool __ll_handle_conditional_array_initialization_transposed ( long * __restrict lNodeFlags, bool isLeaf, long nodeCode, long setBranch, long iNodes, long siteID, long siteFrom, long siteCount, _SimpleList&            siteOrdering, hyFloat * __restrict parentConditionals, hyFloat const * __restrict tMatrix, _Vector const*     lNodeResolutions, hyFloat *& childVector, _SimpleList const* tcc, long &currentTCCBit, long& currentTCCIndex, hyFloat * & lastUpdatedSite, long* __restrict  setBranchTo) {
    

    if (isLeaf) {
        long siteState;
        if (setBranch != nodeCode + iNodes) {
            siteState = lNodeFlags[nodeCode*siteCount + siteOrdering.list_data[siteID]] ;
        } else {
            siteState = setBranchTo[siteOrdering.list_data[siteID]] ;
        }
        if (__builtin_expect(siteState >= 0L,1)) {
            // a single character state; sweep down the appropriate column
#ifdef _SLKP_USE_ARM_NEON
            for (long k = 0; k < (D>>2<<2); k+=4) {
                float64x2x2_t PC = vld1q_f64_x2 (parentConditionals+k),
                              CC = vld1q_f64_x2 (tMatrix+siteState*D+k);
                
  
                PC.val[0] = vmulq_f64(PC.val[0], CC.val[0]);
                PC.val[1] = vmulq_f64(PC.val[1], CC.val[1]);

                vst1q_f64_x2 (parentConditionals+k, PC);
                
                
            }
            for (long k = (D>>2<<2); k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState*D+k];
            }
#endif
            
#ifdef _SLKP_USE_AVX_INTRINSICS
            for (long k = 0; k < (D>>2<<2); k+=4) {
                __m256d PC = _mm256_loadu_pd (parentConditionals+k),
                        CC = _mm256_loadu_pd (tMatrix+siteState*D+k);
                
                _mm256_storeu_pd (parentConditionals+k, _mm256_mul_pd (_mm256_loadu_pd (parentConditionals+k),_mm256_loadu_pd (tMatrix+siteState*D+k)));
  
            }
            for (long k = (D>>2<<2); k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState*D+k];
            }
#endif
            
#ifdef _SLKP_USE_SSE_INTRINSICS
            for (long k = 0; k < (D>>2<<2); k+=4) {
                
                _mm_storeu_pd (parentConditionals+k, _mm_mul_pd (_mm_loadu_pd (parentConditionals+k),_mm_loadu_pd (tMatrix+siteState*D+k)));
                _mm_storeu_pd (parentConditionals+k + 2, _mm_mul_pd (_mm_loadu_pd (parentConditionals+k+2),_mm_loadu_pd (tMatrix+siteState*D+k+2)));
                
            }
            for (long k = (D>>2<<2); k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState*D+k];
            }
#endif
            
#if not defined _SLKP_USE_AVX_INTRINSICS && not defined _SLKP_USE_SSE_INTRINSICS && not defined _SLKP_USE_ARM_NEON
            #pragma GCC unroll 4
            for (long k = 0L; k < D; k++) {
                parentConditionals[k] *= tMatrix[siteState*D+k];
            }
#endif
            return true;
        } else {
            childVector = lNodeResolutions->theData + (-siteState-1) * D;
        }
    } else {
        if (tcc) {
            if (__builtin_expect((tcc->list_data[currentTCCIndex] & bitMaskArray.masks[currentTCCBit]) > 0 && siteID > siteFrom,0)) {
                #pragma GCC unroll 4
                for (long k = 0L; k < D; k++) {
                    childVector[k] = lastUpdatedSite[k];
                }
                //memcpy (childVector, lastUpdatedSite, sizeof (hyFloat)*D);
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
            //#pragma unroll(4)
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
                //#pragma unroll(4)
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
    }*/
    /*if (isinf(sum)) {
        printf ("Infinity at __ll_loop_handle_scaling: site = %ld, branch = %ld @ eval %ld %g\n", siteID, parentCode,  likeFuncEvalCallCount, scalingAdjustments [parentCode*siteCount + siteID]);
        for (int k = 0; k < D; k++)
            printf ("\t %d %g\n", k, parentConditionals[k]);
    }*/
    
    if (sum < _lfScalingFactorThreshold && sum > 0.0) {
        
        hyFloat scaler = _computeBoostScaler(scalingAdjustments [parentCode*siteCount + siteID] * _lfScalerUpwards, sum, didScale);
        
        //fprintf (stderr, "UP %ld (%ld) %lg\n", didScale, parentCode, scaler);
        
        /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
            fprintf (stderr, "UP %ld (%ld) %lg\n", didScale, parentCode, scaler);
        }*/
        if (didScale) {
           // #pragma unroll(4)
            #pragma GCC unroll 4
            for (long c = 0; c < D; c++) {
                parentConditionals [c] *= scaler;
                //if (!isfinite (parentConditionals [c])) {
                //    HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch at branch " & parentCode & "; this is not a recoverable error and indicates some serious COVFEFE taking place.");
                //}
                //if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                //fprintf (stderr, "%ld=>%g\n", c, parentConditionals [c]);
                // }
                //if (isnan (parentConditionals [c])) {
                //    HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch at branch " & parentCode & "; this is not a recoverable error and indicates some serious COVFEFE taking place.");
                //}
                //if (likeFuncEvalCallCount > 1261 && siteID == 1 && parentCode == 117) {
                //   fprintf (stderr, "%ld=>%g\n", c, parentConditionals [c]);
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
        if (sum > _lfScalerUpwards) {
            if (sum < HUGE_VAL) { // no point scaling an infinity
                
                hyFloat scaler = _computeReductionScaler (scalingAdjustments [parentCode*siteCount + siteID] * _lfScalingFactorThreshold, sum, didScale);
                /*if (likeFuncEvalCallCount == 15098 && siteID == 91) {
                    fprintf (stderr, "DOWN %ld (%ld) %lg\n", didScale, parentCode, scaler);
                }*/
                
                if (didScale) {
                    //#pragma unroll(4)
                    #pragma GCC unroll 4
                    for (long c = 0; c < D; c++) {
                        parentConditionals [c] *= scaler;
                        //if (!isfinite (parentConditionals [c])) {
                        //    HandleApplicationError(_String("Site ") & siteID & " evaluated to a NaN probability in ComputeTreeBlockByBranch at branch " & parentCode & "; this is not a recoverable error and indicates some serious COVFEFE taking place.");
                        //}
                        //if (likeFuncEvalCallCount > 1261 && siteID == 1 && parentCode == 117) {
                        //   fprintf (stderr, "%ld=>%g\n", c, parentConditionals [c]);
                        //}
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
            //#pragma unroll(8)
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
                    //#pragma unroll(8)
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
        memset (pp, 0, (siteTo-siteFrom) * D * sizeof (hyFloat));
        for (long k = siteFrom; k < siteTo; k++, pp += D) {
            pp[setBranchTo[siteOrdering.list_data[k]]] = localScalingFactor[k];
        }
    } else {
        for (long k = siteFrom; k < siteTo; k++, pp += D) {
            hyFloat lsf = localScalingFactor[k];
//#pragma unroll(4)
#pragma GCC unroll 4
            for (long s = 0; s < D; s++) {
                pp[s] = lsf;
            }
        }
    }
}

#ifdef _HY_NO_INTRINSICS_MARKER

    void _mx_vect_4x4 (double *cv, double const *M, double const *V, int stride) {
        
        int O2 = stride << 1,
            O3 = O2 + stride;
        
        cv[0] = M[0]*V[0] + M[1]*V[1] + M[2]*V[2] + M[3]*V[3];
        
        cv[1] = M[stride]*V[0] + M[stride+1]*V[1] + M[stride+2]*V[2] + M[stride+3]*V[3];
        
        cv[2] = M[O2]*V[0] + M[O2+1]*V[1] + M[O2+2]*V[2] + M[O2+3]*V[3];
        
        cv[3] = M[O3]*V[0] + M[O3+1]*V[1] + M[O3+2]*V[2] + M[O3+3]*V[3];
    }
    
    void _mx_vect_4x4_add (double *cv, double const *M, double const *V, int stride) {
        int O2 = stride << 1,
            O3 = O2 + stride;
        
        cv[0] += M[0]*V[0] + M[1]*V[1] + M[2]*V[2] + M[3]*V[3];
        
        cv[1] += M[stride]*V[0] + M[stride+1]*V[1] + M[stride+2]*V[2] + M[stride+3]*V[3];
        
        cv[2] += M[O2]*V[0] + M[O2+1]*V[1] + M[O2+2]*V[2] + M[O2+3]*V[3];
        
        cv[3] += M[O3]*V[0] + M[O3+1]*V[1] + M[O3+2]*V[2] + M[O3+3]*V[3];
    }

    void _hy_matrix_vector_product_blocked_4x4 (double * C, double const *M, double const *V, int D) {
        auto offset = [D](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
 

        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    double accumulator[4];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    C[off]   = accumulator[0];
                    C[off+1] = accumulator[1];
                    C[off+2] = accumulator[2];
                    C[off+3] = accumulator[3];

                }
                break;
            }
                
            case 1: {
                for (int i = 0; i < blocks; i++) {
                    
                    double accumulator[4];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    double last_v = V [offset (0,blocks)];

                    C[off]     = accumulator [0] + M[moffset]     * last_v;
                    C[off + 1] = accumulator [1] + M[moffset+D]   * last_v;
                    C[off + 2] = accumulator [2] + M[moffset+2*D] * last_v;
                    C[off + 3] = accumulator [3] + M[moffset+3*D] * last_v;
                }
                
                M += offset (blocks,0);
                
                double accumulator[4];

                accumulator[0] = M[0]*V[0];
                accumulator[1] = M[1]*V[1];
                accumulator[2] = M[2]*V[2];
                accumulator[3] = M[3]*V[3];

                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    accumulator[0] += M[off]*V[off];
                    accumulator[1] += M[off+1]*V[off+1];
                    accumulator[2] += M[off+2]*V[off+2];
                    accumulator[3] += M[off+3]*V[off+3];
                }
                
                C[D-1] = M[D-1] * V[D-1] + (accumulator[0] + accumulator[1]) + (accumulator[2] + accumulator[3]);
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    
                    double accumulator[4];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    double  last_v_m1 = V [offset (0,blocks)],
                            last_v    = V [offset (0,blocks) + 1];

                    C[off]     = accumulator [0] + M[moffset]     * last_v_m1 + M[moffset + 1] * last_v;
                    C[off + 1] = accumulator [1] + M[moffset+D]   * last_v_m1 + M[moffset + 1 + D] * last_v;
                    C[off + 2] = accumulator [2] + M[moffset+2*D] * last_v_m1 + M[moffset + 1 + 2*D] * last_v;
                    C[off + 3] = accumulator [3] + M[moffset+3*D] * last_v_m1 + M[moffset + 1 + 3*D] * last_v;
                }
                
                
                M += offset (blocks,0);
                
                double accumulator  [4],
                       accumulator2 [4];

                accumulator[0] = M[0]*V[0];
                accumulator[1] = M[1]*V[1];
                accumulator[2] = M[2]*V[2];
                accumulator[3] = M[3]*V[3];
                
                accumulator2[0] = M[D]*V[0];
                accumulator2[1] = M[D+1]*V[1];
                accumulator2[2] = M[D+2]*V[2];
                accumulator2[3] = M[D+3]*V[3];

                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    accumulator[0] += M[off]*V[off];
                    accumulator[1] += M[off+1]*V[off+1];
                    accumulator[2] += M[off+2]*V[off+2];
                    accumulator[3] += M[off+3]*V[off+3];
                    accumulator2[0] += M[D+off]*V[off];
                    accumulator2[1] += M[D+off+1]*V[off+1];
                    accumulator2[2] += M[D+off+2]*V[off+2];
                    accumulator2[3] += M[D+off+3]*V[off+3];
                }
                
                C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1] + (accumulator[0] + accumulator[1]) + (accumulator[2] + accumulator[3]);
                C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + (accumulator2[0] + accumulator2[1]) + (accumulator2[2] + accumulator2[3]);
     
                break;
            }
                
            case 3: {
                for (int i = 0; i < blocks; i++) {
                    
                    double accumulator[4];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    double  last_v_m2 = V [offset (0,blocks)],
                    last_v_m1 = V [offset (0,blocks) + 1],
                    last_v    = V [offset (0,blocks) + 2];
                    
                    C[off]     = accumulator [0] + M[moffset]     * last_v_m2 + M[moffset + 1] * last_v_m1 + M[moffset + 2] * last_v;
                    C[off + 1] = accumulator [1] + M[moffset+D]   * last_v_m2 + M[moffset + 1 + D] * last_v_m1 + M[moffset + 2 + D] * last_v;
                    C[off + 2] = accumulator [2] + M[moffset+2*D] * last_v_m2 + M[moffset + 1 + 2*D] * last_v_m1 + M[moffset + 2 + 2*D] * last_v;
                    C[off + 3] = accumulator [3] + M[moffset+3*D] * last_v_m2 + M[moffset + 1 + 3*D] * last_v_m1 + M[moffset + 2 + 3*D] * last_v;
                }
                
                
                M += offset (blocks,0);
                
                double accumulator  [4],
                       accumulator2 [4],
                       accumulator3 [4];
                
                accumulator[0] = M[0]*V[0];
                accumulator[1] = M[1]*V[1];
                accumulator[2] = M[2]*V[2];
                accumulator[3] = M[3]*V[3];
                
                accumulator2[0] = M[D]*V[0];
                accumulator2[1] = M[D+1]*V[1];
                accumulator2[2] = M[D+2]*V[2];
                accumulator2[3] = M[D+3]*V[3];
                
                accumulator3[0] = M[2*D]*V[0];
                accumulator3[1] = M[2*D+1]*V[1];
                accumulator3[2] = M[2*D+2]*V[2];
                accumulator3[3] = M[2*D+3]*V[3];
                

                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    accumulator[0] += M[off]*V[off];
                    accumulator[1] += M[off+1]*V[off+1];
                    accumulator[2] += M[off+2]*V[off+2];
                    accumulator[3] += M[off+3]*V[off+3];
                    
                    accumulator2[0] += M[D+off]*V[off];
                    accumulator2[1] += M[D+off+1]*V[off+1];
                    accumulator2[2] += M[D+off+2]*V[off+2];
                    accumulator2[3] += M[D+off+3]*V[off+3];

                    accumulator3[0] += M[2*D+off]*V[off];
                    accumulator3[1] += M[2*D+off+1]*V[off+1];
                    accumulator3[2] += M[2*D+off+2]*V[off+2];
                    accumulator3[3] += M[2*D+off+3]*V[off+3];
                }
                
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + (accumulator[0] + accumulator[1]) + (accumulator[2] + accumulator[3]);
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + (accumulator2[0] + accumulator2[1]) + (accumulator2[2] + accumulator2[3]);
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + (accumulator3[0] + accumulator3[1]) + (accumulator3[2] + accumulator3[3]);
            }

        }
    }

    template<int D> void _hy_mvp_blocked_4x4 (double * C, double const *M, double const *V) {
        _hy_matrix_vector_product_blocked_4x4 (C,M,V,D);
    }

    inline double _handle4x4_pruning_case_direct (double const* childVector, void* T, double* parentConditionals) {
        
        hyFloat * tMatrix = (hyFloat*)T;
        
        hyFloat t1 = childVector[0] - childVector[3],
        t2 = childVector[1] - childVector[3],
        t3 = childVector[2] - childVector[3],
        t4 = childVector[3];
        
        parentConditionals [0] *= (tMatrix[0]  * t1 + tMatrix[1] * t2) + (tMatrix[2] * t3 + t4);
        parentConditionals [1] *= (tMatrix[4]  * t1 + tMatrix[5] * t2) + (tMatrix[6] * t3 + t4);
        parentConditionals [2] *= (tMatrix[8]  * t1 + tMatrix[9] * t2) + (tMatrix[10] * t3 + t4);
        parentConditionals [3] *= (tMatrix[12] * t1 + tMatrix[13] * t2) + (tMatrix[14] * t3 + t4);
        
        return (parentConditionals[0] + parentConditionals[1]) + (parentConditionals[2] + parentConditionals[3]);
    }




    double _hy_vvmult_sum_generic (double * C, double const *M, int D) {

        double sum = 0.;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            double accumulator[8];
            
            accumulator [0] = (C[0] *= M[0]);
            accumulator [1] = (C[1] *= M[1]);
            accumulator [2] = (C[2] *= M[2]);
            accumulator [3] = (C[3] *= M[3]);
            accumulator [4] = (C[4] *= M[4]);
            accumulator [5] = (C[5] *= M[5]);
            accumulator [6] = (C[6] *= M[6]);
            accumulator [7] = (C[7] *= M[7]);

            for (int i = 1; i < blocks; i++) {
                
                int off = i<<3;
                
                accumulator [0] += (C[off    ] *= M[off    ]);
                accumulator [1] += (C[off + 1] *= M[off + 1]);
                accumulator [2] += (C[off + 2] *= M[off + 2]);
                accumulator [3] += (C[off + 3] *= M[off + 3]);
                accumulator [4] += (C[off + 4] *= M[off + 4]);
                accumulator [5] += (C[off + 5] *= M[off + 5]);
                accumulator [6] += (C[off + 6] *= M[off + 6]);
                accumulator [7] += (C[off + 7] *= M[off + 7]);

            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = (accumulator[0] + accumulator[1]) + (accumulator[2] + accumulator[3]) +
                  (accumulator[4] + accumulator[5]) + (accumulator[6] + accumulator[7]);
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 4:
                    sum += (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 5:
                    sum += (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 6:
                    sum += (C[D-6] *= M[D-6]) +
                           (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 7:
                    sum += (C[D-7] *= M[D-7])
                         + (C[D-6] *= M[D-6])
                         + (C[D-5] *= M[D-5])
                         + (C[D-4] *= M[D-4])
                         + (C[D-3] *= M[D-3])
                         + (C[D-2] *= M[D-2])
                         + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

    template<int D> double _hy_vvmult_sum (double * C, double const *M) {
        return _hy_vvmult_sum_generic(C,M,D);
    }

#else

#if defined  _SLKP_USE_SSE_INTRINSICS

void _mx_vect_4x4 (__m128d *cv, double const *M, double const *V, int stride) {
        
        __m128d       col[2],
                      accumulator[2],
                      accumulator2[2],
                      accumulator3[2],
                      accumulator4[2];
                      
        col[0] = _mm_loadu_pd  (V);
        col[1] = _mm_loadu_pd  (V+2);

 
        accumulator[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M));
        accumulator[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+2));

        accumulator2[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M+stride));
        accumulator2[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+stride+2));

        accumulator3[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M+stride*2));
        accumulator3[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+stride*2+2));

        accumulator4[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M+stride*3));
        accumulator4[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+stride*3+2));

        accumulator[0]  = _mm_add_pd (accumulator[0],accumulator[1]);
        accumulator2[0] = _mm_add_pd (accumulator2[0],accumulator2[1]);
        accumulator3[0] = _mm_add_pd (accumulator3[0],accumulator3[1]);
        accumulator4[0] = _mm_add_pd (accumulator4[0],accumulator4[1]);
        
        cv[0] = _mm_add_pd (_mm_shuffle_pd(accumulator[0],accumulator2[0],0),_mm_shuffle_pd(accumulator[0],accumulator2[0],3));
        cv[1] = _mm_add_pd (_mm_shuffle_pd(accumulator3[0],accumulator4[0],0),_mm_shuffle_pd(accumulator3[0],accumulator4[0],3));
  
         
    }
    
    void _mx_vect_4x4_add (__m128d *cv, double const *M, double const *V, int stride) {
         __m128d       col[2],
                      accumulator[2],
                      accumulator2[2],
                      accumulator3[2],
                      accumulator4[2];
                      
        col[0] = _mm_loadu_pd  (V);
        col[1] = _mm_loadu_pd  (V+2);

 
        accumulator[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M));
        accumulator[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+2));

        accumulator2[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M+stride));
        accumulator2[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+stride+2));

        accumulator3[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M+stride*2));
        accumulator3[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+stride*2+2));

        accumulator4[0] =  _mm_mul_pd(col[0],_mm_loadu_pd(M+stride*3));
        accumulator4[1] =  _mm_mul_pd(col[1],_mm_loadu_pd(M+stride*3+2));

        accumulator[0]  = _mm_add_pd (accumulator[0],accumulator[1]);
        accumulator2[0] = _mm_add_pd (accumulator2[0],accumulator2[1]);
        accumulator3[0] = _mm_add_pd (accumulator3[0],accumulator3[1]);
        accumulator4[0] = _mm_add_pd (accumulator4[0],accumulator4[1]);
        
        cv[0] = _mm_add_pd(cv[0],_mm_add_pd (_mm_shuffle_pd(accumulator[0],accumulator2[0],0),_mm_shuffle_pd(accumulator[0],accumulator2[0],3)));
        cv[1] = _mm_add_pd(cv[1],_mm_add_pd (_mm_shuffle_pd(accumulator3[0],accumulator4[0],0),_mm_shuffle_pd(accumulator3[0],accumulator4[0],3)));
        
    }

    void _hy_matrix_vector_product_blocked_4x4 (double * C, double const *M, double const *V, int D) {
        auto offset = [D](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
 
        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);

                }
                break;
            }
                
            case 1: {
                for (int i = 0; i < blocks; i++) {
                    
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    __m128d last_v = _mm_loaddup_pd (V + offset (0,blocks));
                    
                    accumulator [0] = _mm_add_pd (accumulator [0], _mm_mul_pd (last_v, _mm_setr_pd (M[moffset],M[moffset+D])));
                    accumulator [1] = _mm_add_pd (accumulator [1], _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2*D],M[moffset+3*D])));

                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);
                }
                
                M += offset (blocks,0);
                
                __m128d accumulator[2];
                accumulator[0] = _mm_mul_pd (_mm_loadu_pd(M), _mm_loadu_pd(V));
                accumulator[1] = _mm_mul_pd (_mm_loadu_pd(M+2), _mm_loadu_pd(V+2));
                
  
 
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
  
                    accumulator[0] = _mm_add_pd (accumulator[0], _mm_mul_pd (_mm_loadu_pd(M+off), _mm_loadu_pd(V+off)));
                    accumulator[1] = _mm_add_pd (accumulator[1], _mm_mul_pd (_mm_loadu_pd(M+off+2), _mm_loadu_pd(V+off+2)));
               }
                
                C[D-1] = M[D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator[0],accumulator[1]));
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    __m128d last_vm1 = _mm_loaddup_pd (V + offset (0,blocks)),
                            last_v  =  _mm_loaddup_pd (V + offset (0,blocks) + 1);
                    
                    accumulator [0] = _mm_add_pd (accumulator [0], 
                        _mm_add_pd(
                            _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+1],M[moffset+D+1])), 
                            _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset],M[moffset+D]))));
                            
                    accumulator [1] = _mm_add_pd (accumulator [1], _mm_add_pd(_mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2*D+1],M[moffset+3*D+1])), _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset+2*D],M[moffset+3*D]))));

  
                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);
                }
                
                M += offset (blocks,0);
                
                __m128d accumulator  [2],
                        accumulator2 [2], 
                        v1 = _mm_loadu_pd(V),
                        v2 = _mm_loadu_pd(V+2);
                        
                        
                accumulator[0] = _mm_mul_pd (_mm_loadu_pd(M), v1);
                accumulator[1] = _mm_mul_pd (_mm_loadu_pd(M+2), v2);
                accumulator2[0] = _mm_mul_pd (_mm_loadu_pd(M+D), v1);
                accumulator2[1] = _mm_mul_pd (_mm_loadu_pd(M+D+2), v2);
                
  
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v1 = _mm_loadu_pd(V+off);
                    v2 = _mm_loadu_pd(V+2+off);
                    accumulator[0] = _mm_add_pd (accumulator[0], _mm_mul_pd (_mm_loadu_pd(M+off), v1));
                    accumulator[1] = _mm_add_pd (accumulator[1], _mm_mul_pd (_mm_loadu_pd(M+off+2), v2));
                    accumulator2[0] = _mm_add_pd (accumulator2[0], _mm_mul_pd (_mm_loadu_pd(M+off+D), v1));
                    accumulator2[1] = _mm_add_pd (accumulator2[1], _mm_mul_pd (_mm_loadu_pd(M+off+D+2), v2));
               }
                
               C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1]+ _sse_sum_2 (_mm_add_pd (accumulator[0],accumulator[1]));
               C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator2[0],accumulator2[1]));

                break;
            }
                
            case 3: {
            
              for (int i = 0; i < blocks; i++) {
                    
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    __m128d last_vm2 = _mm_loaddup_pd (V + offset (0,blocks)),
                            last_vm1 = _mm_loaddup_pd (V + offset (0,blocks) + 1),
                            last_v  =  _mm_loaddup_pd (V + offset (0,blocks) + 2);
                    
                    accumulator [0] = _mm_add_pd (
                        _mm_add_pd(
                            accumulator [0], 
                            _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2],M[moffset+D+2]))
                        ), 
                        _mm_add_pd(
                            _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset+1],M[moffset+D+1])), 
                            _mm_mul_pd (last_vm2, _mm_setr_pd (M[moffset],M[moffset+D]))
                        )
                    );

                    accumulator [1] = _mm_add_pd (
                        _mm_add_pd(
                            accumulator [1], 
                            _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2*D+2],M[moffset+3*D+2]))
                        ), 
                        _mm_add_pd(
                            _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset+2*D+1],M[moffset+3*D+1])), 
                            _mm_mul_pd (last_vm2, _mm_setr_pd (M[moffset+2*D],M[moffset+3*D]))
                        )
                    );

                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);
                }
                
                
                M += offset (blocks,0);
                
                __m128d accumulator   [2],
                        accumulator2 [2], 
                        accumulator3 [2],
                        v1 = _mm_loadu_pd(V),
                        v2 = _mm_loadu_pd(V+2);
                        
                        
                accumulator[0] = _mm_mul_pd (_mm_loadu_pd(M), v1);
                accumulator[1] = _mm_mul_pd (_mm_loadu_pd(M+2), v2);
                accumulator2[0] = _mm_mul_pd (_mm_loadu_pd(M+D), v1);
                accumulator2[1] = _mm_mul_pd (_mm_loadu_pd(M+D+2), v2);
                accumulator3[0] = _mm_mul_pd (_mm_loadu_pd(M+2*D), v1);
                accumulator3[1] = _mm_mul_pd (_mm_loadu_pd(M+2*D+2), v2);
                
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v1 = _mm_loadu_pd(V+off);
                    v2 = _mm_loadu_pd(V+2+off);
                    
                    accumulator[0]  = _mm_add_pd (accumulator[0], _mm_mul_pd (_mm_loadu_pd(M+off), v1));
                    accumulator[1]  = _mm_add_pd (accumulator[1], _mm_mul_pd (_mm_loadu_pd(M+off+2), v2));
                    accumulator2[0] = _mm_add_pd (accumulator2[0], _mm_mul_pd (_mm_loadu_pd(M+off+D), v1));
                    accumulator2[1] = _mm_add_pd (accumulator2[1], _mm_mul_pd (_mm_loadu_pd(M+off+D+2), v2));
                    accumulator3[0] = _mm_add_pd (accumulator3[0], _mm_mul_pd (_mm_loadu_pd(M+off+2*D), v1));
                    accumulator3[1] = _mm_add_pd (accumulator3[1], _mm_mul_pd (_mm_loadu_pd(M+off+2*D+2), v2));
               }
                
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator[0],accumulator[1]));
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator2[0],accumulator2[1]));
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator3[0],accumulator3[1]));


            }

        }
    }

    template<int D> void _hy_mvp_blocked_4x4 (double * C, double const *M, double const *V) {
        auto offset = [](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
 
        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);

                }
                break;
            }
                
            case 1: {
                for (int i = 0; i < blocks; i++) {
                    
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    __m128d last_v = _mm_loaddup_pd (V + offset (0,blocks));
                    
                    accumulator [0] = _mm_add_pd (accumulator [0], _mm_mul_pd (last_v, _mm_setr_pd (M[moffset],M[moffset+D])));
                    accumulator [1] = _mm_add_pd (accumulator [1], _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2*D],M[moffset+3*D])));

                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);
                }
                
                M += offset (blocks,0);
                
                __m128d accumulator[2];
                accumulator[0] = _mm_mul_pd (_mm_loadu_pd(M), _mm_loadu_pd(V));
                accumulator[1] = _mm_mul_pd (_mm_loadu_pd(M+2), _mm_loadu_pd(V+2));
                
  
 
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
  
                    accumulator[0] = _mm_add_pd (accumulator[0], _mm_mul_pd (_mm_loadu_pd(M+off), _mm_loadu_pd(V+off)));
                    accumulator[1] = _mm_add_pd (accumulator[1], _mm_mul_pd (_mm_loadu_pd(M+off+2), _mm_loadu_pd(V+off+2)));
               }
                
                C[D-1] = M[D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator[0],accumulator[1]));
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    __m128d last_vm1 = _mm_loaddup_pd (V + offset (0,blocks)),
                            last_v  =  _mm_loaddup_pd (V + offset (0,blocks) + 1);
                    
                    accumulator [0] = _mm_add_pd (accumulator [0], 
                        _mm_add_pd(
                            _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+1],M[moffset+D+1])), 
                            _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset],M[moffset+D]))));
                            
                    accumulator [1] = _mm_add_pd (accumulator [1], _mm_add_pd(_mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2*D+1],M[moffset+3*D+1])), _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset+2*D],M[moffset+3*D]))));

  
                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);
                }
                
                M += offset (blocks,0);
                
                __m128d accumulator  [2],
                        accumulator2 [2], 
                        v1 = _mm_loadu_pd(V),
                        v2 = _mm_loadu_pd(V+2);
                        
                        
                accumulator[0] = _mm_mul_pd (_mm_loadu_pd(M), v1);
                accumulator[1] = _mm_mul_pd (_mm_loadu_pd(M+2), v2);
                accumulator2[0] = _mm_mul_pd (_mm_loadu_pd(M+D), v1);
                accumulator2[1] = _mm_mul_pd (_mm_loadu_pd(M+D+2), v2);
                
  
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v1 = _mm_loadu_pd(V+off);
                    v2 = _mm_loadu_pd(V+2+off);
                    accumulator[0] = _mm_add_pd (accumulator[0], _mm_mul_pd (_mm_loadu_pd(M+off), v1));
                    accumulator[1] = _mm_add_pd (accumulator[1], _mm_mul_pd (_mm_loadu_pd(M+off+2), v2));
                    accumulator2[0] = _mm_add_pd (accumulator2[0], _mm_mul_pd (_mm_loadu_pd(M+off+D), v1));
                    accumulator2[1] = _mm_add_pd (accumulator2[1], _mm_mul_pd (_mm_loadu_pd(M+off+D+2), v2));
               }
                
               C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1]+ _sse_sum_2 (_mm_add_pd (accumulator[0],accumulator[1]));
               C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator2[0],accumulator2[1]));

                break;
            }
                
            case 3: {
            
              for (int i = 0; i < blocks; i++) {
                    
                    __m128d accumulator[2];
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                    
                    __m128d last_vm2 = _mm_loaddup_pd (V + offset (0,blocks)),
                            last_vm1 = _mm_loaddup_pd (V + offset (0,blocks) + 1),
                            last_v  =  _mm_loaddup_pd (V + offset (0,blocks) + 2);
                    
                    accumulator [0] = _mm_add_pd (
                        _mm_add_pd(
                            accumulator [0], 
                            _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2],M[moffset+D+2]))
                        ), 
                        _mm_add_pd(
                            _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset+1],M[moffset+D+1])), 
                            _mm_mul_pd (last_vm2, _mm_setr_pd (M[moffset],M[moffset+D]))
                        )
                    );

                    accumulator [1] = _mm_add_pd (
                        _mm_add_pd(
                            accumulator [1], 
                            _mm_mul_pd (last_v, _mm_setr_pd (M[moffset+2*D+2],M[moffset+3*D+2]))
                        ), 
                        _mm_add_pd(
                            _mm_mul_pd (last_vm1, _mm_setr_pd (M[moffset+2*D+1],M[moffset+3*D+1])), 
                            _mm_mul_pd (last_vm2, _mm_setr_pd (M[moffset+2*D],M[moffset+3*D]))
                        )
                    );

                    _mm_storeu_pd (C+off, accumulator[0]);
                    _mm_storeu_pd (C+off+2, accumulator[1]);
                }
                
                
                M += offset (blocks,0);
                
                __m128d accumulator   [2],
                        accumulator2 [2], 
                        accumulator3 [2],
                        v1 = _mm_loadu_pd(V),
                        v2 = _mm_loadu_pd(V+2);
                        
                        
                accumulator[0] = _mm_mul_pd (_mm_loadu_pd(M), v1);
                accumulator[1] = _mm_mul_pd (_mm_loadu_pd(M+2), v2);
                accumulator2[0] = _mm_mul_pd (_mm_loadu_pd(M+D), v1);
                accumulator2[1] = _mm_mul_pd (_mm_loadu_pd(M+D+2), v2);
                accumulator3[0] = _mm_mul_pd (_mm_loadu_pd(M+2*D), v1);
                accumulator3[1] = _mm_mul_pd (_mm_loadu_pd(M+2*D+2), v2);
                
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v1 = _mm_loadu_pd(V+off);
                    v2 = _mm_loadu_pd(V+2+off);
                    
                    accumulator[0]  = _mm_add_pd (accumulator[0], _mm_mul_pd (_mm_loadu_pd(M+off), v1));
                    accumulator[1]  = _mm_add_pd (accumulator[1], _mm_mul_pd (_mm_loadu_pd(M+off+2), v2));
                    accumulator2[0] = _mm_add_pd (accumulator2[0], _mm_mul_pd (_mm_loadu_pd(M+off+D), v1));
                    accumulator2[1] = _mm_add_pd (accumulator2[1], _mm_mul_pd (_mm_loadu_pd(M+off+D+2), v2));
                    accumulator3[0] = _mm_add_pd (accumulator3[0], _mm_mul_pd (_mm_loadu_pd(M+off+2*D), v1));
                    accumulator3[1] = _mm_add_pd (accumulator3[1], _mm_mul_pd (_mm_loadu_pd(M+off+2*D+2), v2));
               }
                
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator[0],accumulator[1]));
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator2[0],accumulator2[1]));
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + _sse_sum_2 (_mm_add_pd (accumulator3[0],accumulator3[1]));


            }

        }
    }

    inline double _handle4x4_pruning_case_direct (double const* childVector, void* T, double* parentConditionals) {
    
        __m128d  *  TM = (__m128d*)T;
        
        __m128d     c0     = _mm_set1_pd (childVector[0]-childVector[3]),
                    c1     = _mm_set1_pd (childVector[1]-childVector[3]),
                    c2     = _mm_set1_pd (childVector[2]-childVector[3]),
                    c3     = _mm_set1_pd (childVector[3]);
        
        __m128d     t[4];
      
        
        t[0] = _mm_add_pd (_mm_mul_pd (c1, TM[2]), _mm_mul_pd (c0, TM[0]));
        t[1] = _mm_add_pd (_mm_mul_pd (c1, TM[3]), _mm_mul_pd (c0, TM[1]));
        t[2] = _mm_add_pd (_mm_mul_pd (c2, TM[4]), c3);
        t[3] = _mm_add_pd (_mm_mul_pd (c2, TM[5]), c3);
       
        
        t[0] = _mm_mul_pd  (_mm_loadu_pd (parentConditionals), _mm_add_pd (t[0], t[2]));
        _mm_storeu_pd (parentConditionals,t[0]);
        t[1] = _mm_mul_pd  (_mm_loadu_pd (parentConditionals+2), _mm_add_pd (t[1], t[3]));
        _mm_storeu_pd (parentConditionals+2,t[1]);
        
        return _sse_sum_2 (_mm_add_pd (t[0],t[1]));
        
    }



    double _hy_vvmult_sum_generic (double * C, double const *M, int D) {

        double sum = 0.;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            __m128d  accumulator [4];
            
            accumulator[0] = _mm_mul_pd (_mm_loadu_pd (M), _mm_loadu_pd (C));
            accumulator[1] = _mm_mul_pd (_mm_loadu_pd (M+2), _mm_loadu_pd (C+2));
            accumulator[2] = _mm_mul_pd (_mm_loadu_pd (M+4), _mm_loadu_pd (C+4));
            accumulator[3] = _mm_mul_pd (_mm_loadu_pd (M+6), _mm_loadu_pd (C+6));
            
            _mm_storeu_pd (C, accumulator[0]);
            _mm_storeu_pd (C+2, accumulator[1]);
            _mm_storeu_pd (C+4, accumulator[2]);
            _mm_storeu_pd (C+6, accumulator[3]);

            for (int i = 1; i < blocks; i++) {
                
                int off = i<<3;
                
                __m128d CA[4];
                
                CA[0] =  _mm_mul_pd (_mm_loadu_pd (M+off), _mm_loadu_pd (C+off));
                CA[1] =  _mm_mul_pd (_mm_loadu_pd (M+off+2), _mm_loadu_pd (C+off+2));
                CA[2] =  _mm_mul_pd (_mm_loadu_pd (M+off+4), _mm_loadu_pd (C+off+4));
                CA[3] =  _mm_mul_pd (_mm_loadu_pd (M+off+6), _mm_loadu_pd (C+off+6));
                
                
                accumulator[0] = _mm_add_pd (accumulator[0], CA[0]);
                accumulator[1] = _mm_add_pd (accumulator[1], CA[1]);
                accumulator[2] = _mm_add_pd (accumulator[2], CA[2]);
                accumulator[3] = _mm_add_pd (accumulator[3], CA[3]);

                _mm_storeu_pd (C+off, CA[0]);
                _mm_storeu_pd (C+off+2, CA[1]);
                _mm_storeu_pd (C+off+4, CA[2]);
                _mm_storeu_pd (C+off+6, CA[3]);

            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = _sse_sum_2(_mm_add_pd (_mm_add_pd (accumulator[0], accumulator[1]), _mm_add_pd (accumulator[2], accumulator[3])));
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 4:
                    sum += (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 5:
                    sum += (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 6:
                    sum += (C[D-6] *= M[D-6]) +
                           (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 7:
                    sum += (C[D-7] *= M[D-7])
                         + (C[D-6] *= M[D-6])
                         + (C[D-5] *= M[D-5])
                         + (C[D-4] *= M[D-4])
                         + (C[D-3] *= M[D-3])
                         + (C[D-2] *= M[D-2])
                         + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

    template<int D> double _hy_vvmult_sum (double * C, double const *M) {
    
        double sum = 0.;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            __m128d  accumulator [4];
            
            accumulator[0] = _mm_mul_pd (_mm_loadu_pd (M), _mm_loadu_pd (C));
            accumulator[1] = _mm_mul_pd (_mm_loadu_pd (M+2), _mm_loadu_pd (C+2));
            accumulator[2] = _mm_mul_pd (_mm_loadu_pd (M+4), _mm_loadu_pd (C+4));
            accumulator[3] = _mm_mul_pd (_mm_loadu_pd (M+6), _mm_loadu_pd (C+6));
            
            _mm_storeu_pd (C, accumulator[0]);
            _mm_storeu_pd (C+2, accumulator[1]);
            _mm_storeu_pd (C+4, accumulator[2]);
            _mm_storeu_pd (C+6, accumulator[3]);

            for (int i = 1; i < blocks; i++) {
                
                int off = i<<3;
                
                __m128d CA[4];
                
                CA[0] =  _mm_mul_pd (_mm_loadu_pd (M+off), _mm_loadu_pd (C+off));
                CA[1] =  _mm_mul_pd (_mm_loadu_pd (M+off+2), _mm_loadu_pd (C+off+2));
                CA[2] =  _mm_mul_pd (_mm_loadu_pd (M+off+4), _mm_loadu_pd (C+off+4));
                CA[3] =  _mm_mul_pd (_mm_loadu_pd (M+off+6), _mm_loadu_pd (C+off+6));
                
                
                accumulator[0] = _mm_add_pd (accumulator[0], CA[0]);
                accumulator[1] = _mm_add_pd (accumulator[1], CA[1]);
                accumulator[2] = _mm_add_pd (accumulator[2], CA[2]);
                accumulator[3] = _mm_add_pd (accumulator[3], CA[3]);

                _mm_storeu_pd (C+off, CA[0]);
                _mm_storeu_pd (C+off+2, CA[1]);
                _mm_storeu_pd (C+off+4, CA[2]);
                _mm_storeu_pd (C+off+6, CA[3]);

            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = _sse_sum_2(_mm_add_pd (_mm_add_pd (accumulator[0], accumulator[1]), _mm_add_pd (accumulator[2], accumulator[3])));
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 4:
                    sum += (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 5:
                    sum += (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 6:
                    sum += (C[D-6] *= M[D-6]) +
                           (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 7:
                    sum += (C[D-7] *= M[D-7])
                         + (C[D-6] *= M[D-6])
                         + (C[D-5] *= M[D-5])
                         + (C[D-4] *= M[D-4])
                         + (C[D-3] *= M[D-3])
                         + (C[D-2] *= M[D-2])
                         + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

#endif  // SSE

#ifdef  _SLKP_USE_AVX_INTRINSICS

inline __m256d _hy_matrix_handle_axv_mfma (__m256d c, __m256d a, __m256d b) {
    #if defined _SLKP_USE_FMA3_INTRINSICS
        return _mm256_fmadd_pd (a,b,c);
    #else
        return _mm256_add_pd (c, _mm256_mul_pd (a,b));
    #endif
}



/*
bool echo_avx_sum_4_non_zero (__m256d const x) {
    double a[4];
    _mm256_storeu_pd(a, x);
    if (a[0] > 0 && a[1] > 0 && a[2] > 0 && a[3] > 0) {
        printf ("%g|%g|%g|%g\n", a[0], a[1], a[2], a[3]);
        return true;
    }
    return false;
}
*/
    
void _mx_vect_4x4 (__m256d &cv, double const *M, double const *V, int stride) {
        
        __m256d       col,
                      accumulator,
                      accumulator2,
                      accumulator3,
                      accumulator4;
                      
        col = _mm256_loadu_pd  (V);
        
        accumulator  = _mm256_mul_pd (col, _mm256_loadu_pd (M));            // 0,1,2,3
        accumulator2 = _mm256_mul_pd (col, _mm256_loadu_pd (M+stride));     // 4,5,6,7
        accumulator3 = _mm256_mul_pd (col, _mm256_loadu_pd (M+2*stride));   // 8,9,10,11
        accumulator4 = _mm256_mul_pd (col, _mm256_loadu_pd (M+3*stride));   // 12,13,14,15
        
        //  Need to transpose the 4x4 then add them to storage
        
        __m256d    t0 = _mm256_unpacklo_pd (accumulator, accumulator2),  // 0,4,2,6
                   t1 = _mm256_unpackhi_pd (accumulator, accumulator2),  // 1,5,3,7
                   t2 = _mm256_unpacklo_pd (accumulator3, accumulator4), // 8,12,10,14
                   t3 = _mm256_unpackhi_pd (accumulator3, accumulator4); // 9,13,11,15
                 
         
        accumulator  = _mm256_permute2f128_pd (t0,t2, 2*16); // 0,4,8,12
        accumulator2 = _mm256_permute2f128_pd (t1,t3, 2*16); // 1,5,9,13
        accumulator3 = _mm256_permute2f128_pd (t0,t2 ,1+3*16);   // 2,6,19,14
        accumulator4 = _mm256_permute2f128_pd (t1,t3, 1+3*16);    // 3,7,11,15
        
       
        cv = _mm256_add_pd (_mm256_add_pd (accumulator,accumulator2), _mm256_add_pd (accumulator3, accumulator4));
        
         
    }
    
    void _mx_vect_4x4_add (__m256d & cv, double const *M, double const *V, int stride) {
         __m256d       col,
                      accumulator,
                      accumulator2,
                      accumulator3,
                      accumulator4;
                      
        col = _mm256_loadu_pd  (V);
        
        accumulator  = _mm256_mul_pd (col, _mm256_loadu_pd (M));
        accumulator2 = _mm256_mul_pd (col, _mm256_loadu_pd (M+stride));
        accumulator3 = _mm256_mul_pd (col, _mm256_loadu_pd (M+2*stride));
        accumulator4 = _mm256_mul_pd (col, _mm256_loadu_pd (M+3*stride));
 
         __m256d    t0 = _mm256_unpacklo_pd (accumulator, accumulator2),  // 0,4,2,6
                   t1 = _mm256_unpackhi_pd (accumulator, accumulator2),  // 1,5,3,7
                   t2 = _mm256_unpacklo_pd (accumulator3, accumulator4), // 8,12,10,14
                   t3 = _mm256_unpackhi_pd (accumulator3, accumulator4); // 9,13,11,15
        
        accumulator  = _mm256_permute2f128_pd (t0,t2, 2*16); // 0,4,8,12
        accumulator2 = _mm256_permute2f128_pd (t1,t3, 2*16); // 1,5,9,13
        accumulator3 = _mm256_permute2f128_pd (t0,t2 ,1+3*16);   // 2,6,19,14
        accumulator4 = _mm256_permute2f128_pd (t1,t3, 1+3*16);    // 3,7,11,15

        cv = _mm256_add_pd (cv, _mm256_add_pd (_mm256_add_pd (accumulator,accumulator2), _mm256_add_pd (accumulator3, accumulator4)));
        
    }

    void _hy_matrix_vector_product_blocked_4x4 (double * C, double const *M, double const *V, int D) {
        auto offset = [D](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
 
        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                        
                                        
                    _mm256_storeu_pd (C+off, accumulator);

                }
                break;
            }
                
            case 1: {
                
                for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                     
                    __m256d last_v = _mm256_broadcast_sd (V + offset (0,blocks));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_v, _mm256_setr_pd (M[moffset], M[moffset + D], M[moffset + 2*D], M[moffset + 3*D]));
                    _mm256_storeu_pd (C+off, accumulator);

                }
                
                
                M += offset (blocks,0);
                
                __m256d accumulator;
                accumulator = _mm256_mul_pd (_mm256_loadu_pd(M), _mm256_loadu_pd(V));
                
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;                   
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, _mm256_loadu_pd (M+off), _mm256_loadu_pd(V+off));
               }
                
                C[D-1] = M[D-1] * V[D-1] + _avx_sum_4 (accumulator);
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                     
                    __m256d last_vm1 = _mm256_broadcast_sd (V + offset (0,blocks)),
                            last_vm  = _mm256_broadcast_sd (V + offset (0,blocks)+1);
                            
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm1, _mm256_setr_pd (M[moffset], M[moffset + D], M[moffset + 2*D], M[moffset + 3*D]));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm, _mm256_setr_pd (M[moffset + 1], M[moffset + D + 1], M[moffset + 2*D + 1], M[moffset + 3*D + 1]));
                    _mm256_storeu_pd (C+off, accumulator);

                }
                
                
                M += offset (blocks,0);
                
                 __m256d accumulator[2], 
                         v =  _mm256_loadu_pd(V);
                         
                 accumulator[0] = _mm256_mul_pd (_mm256_loadu_pd(M), v);
                 accumulator[1] = _mm256_mul_pd (_mm256_loadu_pd(M+D), v);
 
   
                 for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v =  _mm256_loadu_pd(V + off);
                    accumulator[0] = _hy_matrix_handle_axv_mfma (accumulator[0], _mm256_loadu_pd(M+off), v); 
                    accumulator[1] = _hy_matrix_handle_axv_mfma (accumulator[1], _mm256_loadu_pd(M+D+off), v); 
               }
                
                C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1]+ _avx_sum_4 (accumulator[0]);
                C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _avx_sum_4 (accumulator[1]);

                break;
            }
                
            case 3: {
            
              for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                     
                    __m256d last_vm2 = _mm256_broadcast_sd (V + offset (0,blocks)),
                            last_vm1 = _mm256_broadcast_sd (V + offset (0,blocks)+1),
                            last_vm  = _mm256_broadcast_sd (V + offset (0,blocks)+2);
                            
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm2, _mm256_setr_pd (M[moffset], M[moffset + D], M[moffset + 2*D], M[moffset + 3*D]));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm1, _mm256_setr_pd (M[moffset + 1], M[moffset + D + 1], M[moffset + 2*D + 1], M[moffset + 3*D + 1]));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm,  _mm256_setr_pd (M[moffset + 2], M[moffset + D + 2], M[moffset + 2*D + 2], M[moffset + 3*D + 2]));
                    _mm256_storeu_pd (C+off, accumulator);

                }
                 
                 M += offset (blocks,0);
                
                 __m256d accumulator[3], 
                         v =  _mm256_loadu_pd(V);
                         
                 accumulator[0] = _mm256_mul_pd (_mm256_loadu_pd(M), v);
                 accumulator[1] = _mm256_mul_pd (_mm256_loadu_pd(M+D), v);
                 accumulator[2] = _mm256_mul_pd (_mm256_loadu_pd(M+2*D), v);
 
   
                 for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v =  _mm256_loadu_pd(V + off);
                    accumulator[0] = _hy_matrix_handle_axv_mfma (accumulator[0], _mm256_loadu_pd(M+off), v); 
                    accumulator[1] = _hy_matrix_handle_axv_mfma (accumulator[1], _mm256_loadu_pd(M+D+off), v); 
                    accumulator[2] = _hy_matrix_handle_axv_mfma (accumulator[2], _mm256_loadu_pd(M+2*D+off), v); 
                }
                 
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + _avx_sum_4 (accumulator[0]);
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _avx_sum_4 (accumulator[1]);
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + _avx_sum_4 (accumulator[2]);

 
                break;

            }

        }
    }

    template<int D> void _hy_mvp_blocked_4x4 (double * C, double const *M, double const *V) {
        auto offset = [](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
 
        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                        
                                        
                    _mm256_storeu_pd (C+off, accumulator);

                }
                break;
            }
                
            case 1: {
                
                for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                     
                    __m256d last_v = _mm256_broadcast_sd (V + offset (0,blocks));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_v, _mm256_setr_pd (M[moffset], M[moffset + D], M[moffset + 2*D], M[moffset + 3*D]));
                    _mm256_storeu_pd (C+off, accumulator);

                }
                
                
                M += offset (blocks,0);
                
                __m256d accumulator;
                accumulator = _mm256_mul_pd (_mm256_loadu_pd(M), _mm256_loadu_pd(V));
                
                for (int j = 1; j < blocks; j++) {
                    int off = j << 2;                   
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, _mm256_loadu_pd (M+off), _mm256_loadu_pd(V+off));
               }
                
                C[D-1] = M[D-1] * V[D-1] + _avx_sum_4 (accumulator);
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                     
                    __m256d last_vm1 = _mm256_broadcast_sd (V + offset (0,blocks)),
                            last_vm  = _mm256_broadcast_sd (V + offset (0,blocks)+1);
                            
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm1, _mm256_setr_pd (M[moffset], M[moffset + D], M[moffset + 2*D], M[moffset + 3*D]));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm, _mm256_setr_pd (M[moffset + 1], M[moffset + D + 1], M[moffset + 2*D + 1], M[moffset + 3*D + 1]));
                    _mm256_storeu_pd (C+off, accumulator);

                }
                
                
                M += offset (blocks,0);
                
                 __m256d accumulator[2], 
                         v =  _mm256_loadu_pd(V);
                         
                 accumulator[0] = _mm256_mul_pd (_mm256_loadu_pd(M), v);
                 accumulator[1] = _mm256_mul_pd (_mm256_loadu_pd(M+D), v);
 
   
                 for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v =  _mm256_loadu_pd(V + off);
                    accumulator[0] = _hy_matrix_handle_axv_mfma (accumulator[0], _mm256_loadu_pd(M+off), v); 
                    accumulator[1] = _hy_matrix_handle_axv_mfma (accumulator[1], _mm256_loadu_pd(M+D+off), v); 
               }
                
                C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1]+ _avx_sum_4 (accumulator[0]);
                C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _avx_sum_4 (accumulator[1]);

                break;
            }
                
            case 3: {
            
              for (int i = 0; i < blocks; i++) {
                    __m256d accumulator;
                    int    off = offset (0,i);
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int    moffset = offset (i,blocks);
                     
                    __m256d last_vm2 = _mm256_broadcast_sd (V + offset (0,blocks)),
                            last_vm1 = _mm256_broadcast_sd (V + offset (0,blocks)+1),
                            last_vm  = _mm256_broadcast_sd (V + offset (0,blocks)+2);
                            
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm2, _mm256_setr_pd (M[moffset], M[moffset + D], M[moffset + 2*D], M[moffset + 3*D]));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm1, _mm256_setr_pd (M[moffset + 1], M[moffset + D + 1], M[moffset + 2*D + 1], M[moffset + 3*D + 1]));
                    accumulator = _hy_matrix_handle_axv_mfma (accumulator, last_vm,  _mm256_setr_pd (M[moffset + 2], M[moffset + D + 2], M[moffset + 2*D + 2], M[moffset + 3*D + 2]));
                    _mm256_storeu_pd (C+off, accumulator);

                }
                 
                 M += offset (blocks,0);
                
                 __m256d accumulator[3], 
                         v =  _mm256_loadu_pd(V);
                         
                 accumulator[0] = _mm256_mul_pd (_mm256_loadu_pd(M), v);
                 accumulator[1] = _mm256_mul_pd (_mm256_loadu_pd(M+D), v);
                 accumulator[2] = _mm256_mul_pd (_mm256_loadu_pd(M+2*D), v);
 
   
                 for (int j = 1; j < blocks; j++) {
                    int off = j << 2;
                    v =  _mm256_loadu_pd(V + off);
                    accumulator[0] = _hy_matrix_handle_axv_mfma (accumulator[0], _mm256_loadu_pd(M+off), v); 
                    accumulator[1] = _hy_matrix_handle_axv_mfma (accumulator[1], _mm256_loadu_pd(M+D+off), v); 
                    accumulator[2] = _hy_matrix_handle_axv_mfma (accumulator[2], _mm256_loadu_pd(M+2*D+off), v); 
                }
                 
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + _avx_sum_4 (accumulator[0]);
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + _avx_sum_4 (accumulator[1]);
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + _avx_sum_4 (accumulator[2]);

 
                break;

            }
        }

    }

    inline double _handle4x4_pruning_case_direct (double const* childVector, void* T, double* parentConditionals) {
    
    
        __m256d  *  TM = (__m256d*)T;
        
        __m256d     c0     = _mm256_set1_pd  (childVector[0]-childVector[3]),
                    c1     = _mm256_set1_pd  (childVector[1]-childVector[3]),
                    c2     = _mm256_set1_pd  (childVector[2]-childVector[3]),
                    c3     = _mm256_set1_pd  (childVector[3]);
        
        __m256d     t[2];
      
        
        t[0] = _hy_matrix_handle_axv_mfma (_mm256_mul_pd (c0, TM[0]), c1, TM[1]);
        t[1] = _hy_matrix_handle_axv_mfma (c3,c2, TM[2]);
        
        t[0] = _mm256_mul_pd (_mm256_loadu_pd (parentConditionals), _mm256_add_pd (t[0], t[1]));
        
        _mm256_storeu_pd (parentConditionals, t[0]);        
        return _avx_sum_4 (t[0]);
        
    }



    double _hy_vvmult_sum_generic (double * C, double const *M, int D) {

        double sum = 0.;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            __m256d  accumulator [2];
            
            accumulator[0] = _mm256_mul_pd (_mm256_loadu_pd (M), _mm256_loadu_pd (C));
            accumulator[1] = _mm256_mul_pd (_mm256_loadu_pd (M+4), _mm256_loadu_pd (C+4));
            
            _mm256_storeu_pd (C, accumulator[0]);
            _mm256_storeu_pd (C+4, accumulator[1]);

            for (int i = 1; i < blocks; i++) {
                
                int off = i<<3;
                
                __m256d CA[2];
                
                CA[0] =  _mm256_mul_pd (_mm256_loadu_pd (M+off), _mm256_loadu_pd (C+off));
                CA[1] =  _mm256_mul_pd (_mm256_loadu_pd (M+off+4), _mm256_loadu_pd (C+off+4));
                
                
                accumulator[0] = _mm256_add_pd (accumulator[0], CA[0]);
                accumulator[1] = _mm256_add_pd (accumulator[1], CA[1]);

                _mm256_storeu_pd (C+off, CA[0]);
                _mm256_storeu_pd (C+off+4, CA[1]);
  
            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = _avx_sum_4(_mm256_add_pd (accumulator[0], accumulator[1]));
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 4:
                    sum += (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 5:
                    sum += (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 6:
                    sum += (C[D-6] *= M[D-6]) +
                           (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 7:
                    sum += (C[D-7] *= M[D-7])
                         + (C[D-6] *= M[D-6])
                         + (C[D-5] *= M[D-5])
                         + (C[D-4] *= M[D-4])
                         + (C[D-3] *= M[D-3])
                         + (C[D-2] *= M[D-2])
                         + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

    template<int D> double _hy_vvmult_sum (double * C, double const *M) {
       double sum = 0.;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            __m256d  accumulator [2];
            
            accumulator[0] = _mm256_mul_pd (_mm256_loadu_pd (M), _mm256_loadu_pd (C));
            accumulator[1] = _mm256_mul_pd (_mm256_loadu_pd (M+4), _mm256_loadu_pd (C+4));
            
            _mm256_storeu_pd (C, accumulator[0]);
            _mm256_storeu_pd (C+4, accumulator[1]);

            for (int i = 1; i < blocks; i++) {
                
                int off = i<<3;
                
                __m256d CA[2];
                
                CA[0] =  _mm256_mul_pd (_mm256_loadu_pd (M+off), _mm256_loadu_pd (C+off));
                CA[1] =  _mm256_mul_pd (_mm256_loadu_pd (M+off+4), _mm256_loadu_pd (C+off+4));
                
                
                accumulator[0] = _mm256_add_pd (accumulator[0], CA[0]);
                accumulator[1] = _mm256_add_pd (accumulator[1], CA[1]);

                _mm256_storeu_pd (C+off, CA[0]);
                _mm256_storeu_pd (C+off+4, CA[1]);
  
            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = _avx_sum_4(_mm256_add_pd (accumulator[0], accumulator[1]));
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 4:
                    sum += (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 5:
                    sum += (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 6:
                    sum += (C[D-6] *= M[D-6]) +
                           (C[D-5] *= M[D-5]) +
                           (C[D-4] *= M[D-4]) +
                           (C[D-3] *= M[D-3]) +
                           (C[D-2] *= M[D-2]) +
                           (C[D-1] *= M[D-1]);
                    break;
                case 7:
                    sum += (C[D-7] *= M[D-7])
                         + (C[D-6] *= M[D-6])
                         + (C[D-5] *= M[D-5])
                         + (C[D-4] *= M[D-4])
                         + (C[D-3] *= M[D-3])
                         + (C[D-2] *= M[D-2])
                         + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

#endif // AVX

#ifdef _SLKP_USE_ARM_NEON

    void _mx_vect_4x4 (float64x2x2_t& cv, double const *M, double const *V, int stride) {
        float64x2x2_t col,
                      accumulator,
                      accumulator2,
                      accumulator3,
                      accumulator4,
                      row,
                      row2,
                      row3,
                      row4;

        col  = vld1q_f64_x2 (V);
        row  = vld1q_f64_x2 (M);

        accumulator.val[0] = vmulq_f64 (col.val[0], row.val[0]);
        accumulator.val[1] = vmulq_f64 (col.val[1], row.val[1]);
                             
        row2 = vld1q_f64_x2 (M+stride);
        accumulator2.val[0] = vmulq_f64 (col.val[0], row2.val[0]);
        accumulator2.val[1] = vmulq_f64 (col.val[1], row2.val[1]);
        
        row3 = vld1q_f64_x2 (M+2*stride);
        accumulator3.val[0] = vmulq_f64 (col.val[0], row3.val[0]);
        accumulator3.val[1] = vmulq_f64 (col.val[1], row3.val[1]);
                             
        row4 = vld1q_f64_x2 (M+3*stride);
        accumulator4.val[0] = vmulq_f64 (col.val[0], row4.val[0]);
        accumulator4.val[1] = vmulq_f64 (col.val[1], row4.val[1]);
        
        //float64x2x2_t cv = vld1q_f64_x2 (C);
        accumulator.val[0]  = vaddq_f64 (accumulator.val[0],accumulator.val[1]);
        accumulator2.val[0] = vaddq_f64 (accumulator2.val[0],accumulator2.val[1]);
        accumulator3.val[0] = vaddq_f64 (accumulator3.val[0],accumulator3.val[1]);
        accumulator4.val[0] = vaddq_f64 (accumulator4.val[0],accumulator4.val[1]);

        cv.val[0] = vaddq_f64(vzip1q_f64 (accumulator.val[0],accumulator2.val[0]),vzip2q_f64 (accumulator.val[0],accumulator2.val[0]));
        cv.val[1] = vaddq_f64(vzip1q_f64 (accumulator3.val[0],accumulator4.val[0]),vzip2q_f64 (accumulator3.val[0],accumulator4.val[0]));
    }


    void _mx_vect_4x4_add (float64x2x2_t &cv, double const *M, double const *V, int stride) {
        
        /*float64x2x2_t col;
        float64x2x2_t
                      row,
                      row2,
                      row3,
                      row4;
        
        float64x2_t accumulator,
                    accumulator2,
                    accumulator3,
                    accumulator4;
        
                      
        col   = vld1q_f64_x2 (V);
        row   = vld1q_f64_x2 (M);
        row2  = vld1q_f64_x2 (M+stride);

        accumulator  = vfmaq_f64 (vmulq_f64 (col.val[0], row.val[0]),  col.val[1], row.val[1]);
        accumulator2 = vfmaq_f64 (vmulq_f64 (col.val[0], row2.val[0]), col.val[1], row2.val[1]);

        row   = vld1q_f64_x2 (M+stride*2);
        row2  = vld1q_f64_x2 (M+stride*3);
        accumulator3 = vfmaq_f64 (vmulq_f64 (col.val[0], row.val[0]), col.val[1], row.val[1]);
        accumulator4 = vfmaq_f64 (vmulq_f64 (col.val[0], row2.val[0]), col.val[1], row2.val[1]);

        cv.val[0] = vaddq_f64(cv.val[0],vaddq_f64(vzip1q_f64 (accumulator,accumulator2),vzip2q_f64 (accumulator,accumulator2)));
        cv.val[1] = vaddq_f64(cv.val[1],vaddq_f64(vzip1q_f64 (accumulator3,accumulator4),vzip2q_f64 (accumulator3,accumulator4)));*/
        
        float64x2x2_t col,
                      accumulator,
                      accumulator2,
                      accumulator3,
                      accumulator4,
                      row,
                      row2,
                      row3,
                      row4;

        col  = vld1q_f64_x2 (V);

        row  = vld1q_f64_x2 (M);
        accumulator.val[0] = vmulq_f64 (col.val[0], row.val[0]);
        accumulator.val[1] = vmulq_f64 (col.val[1], row.val[1]);

        row2 = vld1q_f64_x2 (M+stride);

        accumulator2.val[0] = vmulq_f64 (col.val[0], row2.val[0]);
        accumulator2.val[1] = vmulq_f64 (col.val[1], row2.val[1]);
        
        row3 = vld1q_f64_x2 (M+2*stride);
        accumulator3.val[0] = vmulq_f64 (col.val[0], row3.val[0]);
        accumulator3.val[1] = vmulq_f64 (col.val[1], row3.val[1]);
                             
        row4 = vld1q_f64_x2 (M+3*stride);
        accumulator4.val[0] = vmulq_f64 (col.val[0], row4.val[0]);
        accumulator4.val[1] = vmulq_f64 (col.val[1], row4.val[1]);
        
        accumulator.val[0]  = vaddq_f64 (accumulator.val[0],accumulator.val[1]);
        accumulator2.val[0] = vaddq_f64 (accumulator2.val[0],accumulator2.val[1]);
        accumulator3.val[0] = vaddq_f64 (accumulator3.val[0],accumulator3.val[1]);
        accumulator4.val[0] = vaddq_f64 (accumulator4.val[0],accumulator4.val[1]);


        cv.val[0] = vaddq_f64(cv.val[0], vaddq_f64(vzip1q_f64 (accumulator.val[0],accumulator2.val[0]),vzip2q_f64 (accumulator.val[0],accumulator2.val[0])));
        cv.val[1] = vaddq_f64(cv.val[1], vaddq_f64(vzip1q_f64 (accumulator3.val[0],accumulator4.val[0]),vzip2q_f64 (accumulator3.val[0],accumulator4.val[0])));
        
    }

    inline double _handle4x4_pruning_case_direct (double const* childVector, void* tMatrix, double* parentConditionals) {
        
        
        float64x2x2_t * TM = (float64x2x2_t*)tMatrix;
        
        float64x2_t c0     = vdupq_n_f64(childVector[0]-childVector[3]),
                    c1     = vdupq_n_f64(childVector[1]-childVector[3]),
                    c2     = vdupq_n_f64(childVector[2]-childVector[3]),
                    c3     = vdupq_n_f64(childVector[3]);
        
        float64x2_t t[4];
      
        t[0] = vfmaq_f64 (vmulq_f64  (c1,TM[1].val[0]),c0,TM[0].val[0]);
        t[1] = vfmaq_f64 (vmulq_f64  (c1,TM[1].val[1]),c0,TM[0].val[1]);
        t[2] = vaddq_f64 (vmulq_f64  (c2,TM[2].val[0]),c3);
        t[3] = vaddq_f64 (vmulq_f64  (c2,TM[2].val[1]),c3);
        
        t[0] = vmulq_f64  (vld1q_f64 (parentConditionals), vaddq_f64 (t[0], t[2]));
        vst1q_f64 (parentConditionals,t[0]);
        t[1] = vmulq_f64  (vld1q_f64 (parentConditionals+2), vaddq_f64 (t[1],t[3]));
        vst1q_f64 (parentConditionals+2,t[1]);
        return vaddvq_f64(vaddq_f64(t[0],t[1]));
        //return (parentConditionals[0] + parentConditionals[1]) + (parentConditionals[2] + parentConditionals[3]);
    }



    void _hy_matrix_vector_product_blocked_4x4 (double * C, double const *M, double const *V, int D) {
        auto offset = [D](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
        
        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                break;
            }
                
            case 1: {
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int moffset = offset (i,blocks);
                    float64x2_t mv1,
                                mv2;
                    mv1 = vsetq_lane_f64 ( M[moffset],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D],mv2,1);
                    
                    float64x2_t v = vdupq_n_f64 (V [offset (0,blocks)]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                M += offset (blocks,0);
                float64x2x2_t accumulator,
                              row = vld1q_f64_x2 (M),
                              col = vld1q_f64_x2 (V);

                accumulator.val[0] =  vmulq_f64 (row.val[0], col.val[0]);
                accumulator.val[1] =  vmulq_f64 (row.val[1], col.val[1]);
                
                for (int j = 1; j < blocks; j++) {
                    row = vld1q_f64_x2 (M + (j<<2));
                    col = vld1q_f64_x2 (V + (j<<2));
                    
                    accumulator.val[0] =  vfmaq_f64 (accumulator.val[0] ,row.val[0], col.val[0]);
                    accumulator.val[1] =  vfmaq_f64 (accumulator.val[1] ,row.val[1], col.val[1]);
                }
                
                C[D-1] = M[D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]));
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    
                    int moffset = offset (i,blocks);
                    float64x2_t mv1,
                                mv2;
                    
                    mv1 = vsetq_lane_f64 ( M[moffset],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D],mv2,1);
                    
                    float64x2_t v = vdupq_n_f64 (V [offset (0,blocks)]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    
                    mv1 = vsetq_lane_f64 ( M[moffset+1],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D+1],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D+1],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D+1],mv2,1);
                    
                    v = vdupq_n_f64 (V [offset (0,blocks)+1]);
                    
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                
                M += offset (blocks,0);
                float64x2x2_t accumulator,
                              accumulator2,
                              row  = vld1q_f64_x2 (M),
                              row2 = vld1q_f64_x2 (M+D),
                              col = vld1q_f64_x2 (V);

                accumulator.val[0] =  vmulq_f64 (row.val[0], col.val[0]);
                accumulator.val[1] =  vmulq_f64 (row.val[1], col.val[1]);
                accumulator2.val[0] =  vmulq_f64 (row2.val[0], col.val[0]);
                accumulator2.val[1] =  vmulq_f64 (row2.val[1], col.val[1]);

                for (int j = 1; j < blocks; j++) {
                    row = vld1q_f64_x2 (M + (j<<2));
                    row2 = vld1q_f64_x2 (M + D+ (j<<2));
                    col = vld1q_f64_x2 (V + (j<<2));
                    
                    accumulator.val[0] =  vfmaq_f64 (accumulator.val[0] ,row.val[0], col.val[0]);
                    accumulator.val[1] =  vfmaq_f64 (accumulator.val[1] ,row.val[1], col.val[1]);
                    accumulator2.val[0] =  vfmaq_f64 (accumulator2.val[0] ,row2.val[0], col.val[0]);
                    accumulator2.val[1] =  vfmaq_f64 (accumulator2.val[1] ,row2.val[1], col.val[1]);
                }
                C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]));
                C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator2.val[0],accumulator2.val[1]));
     
                break;
            }
                
            case 3: {
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    
                    int moffset = offset (i,blocks);
                    float64x2_t mv1,
                                mv2;
                    
                    mv1 = vsetq_lane_f64 ( M[moffset],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D],mv2,1);
                    
                    float64x2_t v = vdupq_n_f64 (V [offset (0,blocks)]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    
                    mv1 = vsetq_lane_f64 ( M[moffset+1],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D+1],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D+1],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D+1],mv2,1);
                    
                    v = vdupq_n_f64 (V [offset (0,blocks)+1]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);

                    mv1 = vsetq_lane_f64 ( M[moffset+2],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D+2],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D+2],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D+2],mv2,1);
                    
                    v = vdupq_n_f64 (V [offset (0,blocks)+2]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                
                M += offset (blocks,0);
                float64x2x2_t accumulator,
                              accumulator2,
                              accumulator3,
                              row  = vld1q_f64_x2 (M),
                              row2 = vld1q_f64_x2 (M+D),
                              row3 = vld1q_f64_x2 (M+2*D),
                              col = vld1q_f64_x2 (V);

                accumulator.val[0] =  vmulq_f64 (row.val[0], col.val[0]);
                accumulator.val[1] =  vmulq_f64 (row.val[1], col.val[1]);
                accumulator2.val[0] =  vmulq_f64 (row2.val[0], col.val[0]);
                accumulator2.val[1] =  vmulq_f64 (row2.val[1], col.val[1]);
                accumulator3.val[0] =  vmulq_f64 (row3.val[0], col.val[0]);
                accumulator3.val[1] =  vmulq_f64 (row3.val[1], col.val[1]);

                for (int j = 1; j < blocks; j++) {
                    row = vld1q_f64_x2 (M + (j<<2));
                    row2 = vld1q_f64_x2 (M + D+ (j<<2));
                    row3 = vld1q_f64_x2 (M + 2*D+ (j<<2));
                    col = vld1q_f64_x2 (V + (j<<2));
                    
                    accumulator.val[0] =  vfmaq_f64 (accumulator.val[0] ,row.val[0], col.val[0]);
                    accumulator.val[1] =  vfmaq_f64 (accumulator.val[1] ,row.val[1], col.val[1]);
                    accumulator2.val[0] =  vfmaq_f64 (accumulator2.val[0] ,row2.val[0], col.val[0]);
                    accumulator2.val[1] =  vfmaq_f64 (accumulator2.val[1] ,row2.val[1], col.val[1]);
                    accumulator3.val[0] =  vfmaq_f64 (accumulator3.val[0] ,row3.val[0], col.val[0]);
                    accumulator3.val[1] =  vfmaq_f64 (accumulator3.val[1] ,row3.val[1], col.val[1]);
                }
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]));;
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator2.val[0],accumulator2.val[1]));
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator3.val[0],accumulator3.val[1]));

                break;
            }

        }
    }

    template<int D> void _hy_mvp_blocked_4x4 (double * C, double const *M, double const *V) {
        auto offset = [](int i, int j) -> int {
            return (i<<2)*D + (j << 2);
        };
        
        int blocks = D>>2;
        int remainder = D - (blocks<<2);
        
        switch (remainder) {
                
            case 0: {
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                break;
            }
                
            case 1: {
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    int moffset = offset (i,blocks);
                    float64x2_t mv1,
                                mv2;
                    mv1 = vsetq_lane_f64 ( M[moffset],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D],mv2,1);
                    
                    float64x2_t v = vdupq_n_f64 (V [offset (0,blocks)]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                M += offset (blocks,0);
                float64x2x2_t accumulator,
                              row = vld1q_f64_x2 (M),
                              col = vld1q_f64_x2 (V);

                accumulator.val[0] =  vmulq_f64 (row.val[0], col.val[0]);
                accumulator.val[1] =  vmulq_f64 (row.val[1], col.val[1]);
                
                for (int j = 1; j < blocks; j++) {
                    row = vld1q_f64_x2 (M + (j<<2));
                    col = vld1q_f64_x2 (V + (j<<2));
                    
                    accumulator.val[0] =  vfmaq_f64 (accumulator.val[0] ,row.val[0], col.val[0]);
                    accumulator.val[1] =  vfmaq_f64 (accumulator.val[1] ,row.val[1], col.val[1]);
                }
                
                C[D-1] = M[D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]));
                break;
            }
                
            case 2: {
                
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    
                    int moffset = offset (i,blocks);
                    float64x2_t mv1,
                                mv2;
                    
                    mv1 = vsetq_lane_f64 ( M[moffset],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D],mv2,1);
                    
                    float64x2_t v = vdupq_n_f64 (V [offset (0,blocks)]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    
                    mv1 = vsetq_lane_f64 ( M[moffset+1],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D+1],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D+1],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D+1],mv2,1);
                    
                    v = vdupq_n_f64 (V [offset (0,blocks)+1]);
                    
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                
                M += offset (blocks,0);
                float64x2x2_t accumulator,
                              accumulator2,
                              row  = vld1q_f64_x2 (M),
                              row2 = vld1q_f64_x2 (M+D),
                              col = vld1q_f64_x2 (V);

                accumulator.val[0] =  vmulq_f64 (row.val[0], col.val[0]);
                accumulator.val[1] =  vmulq_f64 (row.val[1], col.val[1]);
                accumulator2.val[0] =  vmulq_f64 (row2.val[0], col.val[0]);
                accumulator2.val[1] =  vmulq_f64 (row2.val[1], col.val[1]);

                for (int j = 1; j < blocks; j++) {
                    row = vld1q_f64_x2 (M + (j<<2));
                    row2 = vld1q_f64_x2 (M + D+ (j<<2));
                    col = vld1q_f64_x2 (V + (j<<2));
                    
                    accumulator.val[0] =  vfmaq_f64 (accumulator.val[0] ,row.val[0], col.val[0]);
                    accumulator.val[1] =  vfmaq_f64 (accumulator.val[1] ,row.val[1], col.val[1]);
                    accumulator2.val[0] =  vfmaq_f64 (accumulator2.val[0] ,row2.val[0], col.val[0]);
                    accumulator2.val[1] =  vfmaq_f64 (accumulator2.val[1] ,row2.val[1], col.val[1]);
                }
                C[D-2] = M[D-2] * V[D-2] + M[D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]));
                C[D-1] = M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator2.val[0],accumulator2.val[1]));
     
                break;
            }
                
            case 3: {
                for (int i = 0; i < blocks; i++) {
                    float64x2x2_t accumulator;
                    _mx_vect_4x4 (accumulator, M + offset (i,0), V, D);
                    for (int j=1; j < blocks; j++) {
                        _mx_vect_4x4_add (accumulator, M + offset (i,j), V + offset (0,j), D);
                    }
                    
                    
                    int moffset = offset (i,blocks);
                    float64x2_t mv1,
                                mv2;
                    
                    mv1 = vsetq_lane_f64 ( M[moffset],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D],mv2,1);
                    
                    float64x2_t v = vdupq_n_f64 (V [offset (0,blocks)]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    
                    mv1 = vsetq_lane_f64 ( M[moffset+1],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D+1],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D+1],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D+1],mv2,1);
                    
                    v = vdupq_n_f64 (V [offset (0,blocks)+1]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);

                    mv1 = vsetq_lane_f64 ( M[moffset+2],mv1,0);
                    mv1 = vsetq_lane_f64 ( M[moffset+D+2],mv1,1);
                    mv2 = vsetq_lane_f64 ( M[moffset+2*D+2],mv2,0);
                    mv2 = vsetq_lane_f64 ( M[moffset+3*D+2],mv2,1);
                    
                    v = vdupq_n_f64 (V [offset (0,blocks)+2]);
                    accumulator.val[0] = vfmaq_f64 (accumulator.val[0], mv1, v);
                    vst1q_f64 (C + offset (0,i), accumulator.val[0]);
                    accumulator.val[1] = vfmaq_f64 (accumulator.val[1], mv2, v);
                    vst1q_f64 (C + offset (0,i) + 2, accumulator.val[1]);
                    //handle_small_product_add <4,1> (C + offset (0,i), M + offset (i,blocks), V + offset (0,blocks), D);
                }
                
                M += offset (blocks,0);
                float64x2x2_t accumulator,
                              accumulator2,
                              accumulator3,
                              row  = vld1q_f64_x2 (M),
                              row2 = vld1q_f64_x2 (M+D),
                              row3 = vld1q_f64_x2 (M+2*D),
                              col = vld1q_f64_x2 (V);

                accumulator.val[0] =  vmulq_f64 (row.val[0], col.val[0]);
                accumulator.val[1] =  vmulq_f64 (row.val[1], col.val[1]);
                accumulator2.val[0] =  vmulq_f64 (row2.val[0], col.val[0]);
                accumulator2.val[1] =  vmulq_f64 (row2.val[1], col.val[1]);
                accumulator3.val[0] =  vmulq_f64 (row3.val[0], col.val[0]);
                accumulator3.val[1] =  vmulq_f64 (row3.val[1], col.val[1]);

                for (int j = 1; j < blocks; j++) {
                    row = vld1q_f64_x2 (M + (j<<2));
                    row2 = vld1q_f64_x2 (M + D+ (j<<2));
                    row3 = vld1q_f64_x2 (M + 2*D+ (j<<2));
                    col = vld1q_f64_x2 (V + (j<<2));
                    
                    accumulator.val[0] =  vfmaq_f64 (accumulator.val[0] ,row.val[0], col.val[0]);
                    accumulator.val[1] =  vfmaq_f64 (accumulator.val[1] ,row.val[1], col.val[1]);
                    accumulator2.val[0] =  vfmaq_f64 (accumulator2.val[0] ,row2.val[0], col.val[0]);
                    accumulator2.val[1] =  vfmaq_f64 (accumulator2.val[1] ,row2.val[1], col.val[1]);
                    accumulator3.val[0] =  vfmaq_f64 (accumulator3.val[0] ,row3.val[0], col.val[0]);
                    accumulator3.val[1] =  vfmaq_f64 (accumulator3.val[1] ,row3.val[1], col.val[1]);
                }
                C[D-3] = M[D-3] * V[D-3]  + M[D-2] * V[D-2] + M[D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]));;
                C[D-2] = M[2*D-3] * V[D-3] + M[2*D-2] * V[D-2] + M[2*D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator2.val[0],accumulator2.val[1]));
                C[D-1] = M[3*D-3] * V[D-3] + M[3*D-2] * V[D-2] + M[3*D-1] * V[D-1] + vaddvq_f64(vpaddq_f64(accumulator3.val[0],accumulator3.val[1]));

                break;
            }

        }
    }

// ddot multiply in and accumulate


    template<int D> double _hy_vvmult_sum (double * C, double const *M) {
        double sum = 0.;
        float64x2x2_t accumulator;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            float64x2x4_t CA = vld1q_f64_x4 (C),
            MA = vld1q_f64_x4 (M),
            accumulator;
            
            accumulator.val[0] = vmulq_f64 (CA.val[0], MA.val[0]);
            accumulator.val[1] = vmulq_f64 (CA.val[1], MA.val[1]);
            accumulator.val[2] = vmulq_f64 (CA.val[2], MA.val[2]);
            accumulator.val[3] = vmulq_f64 (CA.val[3], MA.val[3]);

            vst1q_f64_x4 (C, accumulator);
            
            for (int i = 1; i < blocks; i++) {
                CA = vld1q_f64_x4 (C+(i<<3));
                MA = vld1q_f64_x4 (M+(i<<3));
                
                CA.val[0] = vmulq_f64 (CA.val[0], MA.val[0]);
                CA.val[1] = vmulq_f64 (CA.val[1], MA.val[1]);
                CA.val[2] = vmulq_f64 (CA.val[2], MA.val[2]);
                CA.val[3] = vmulq_f64 (CA.val[3], MA.val[3]);

                vst1q_f64_x4 (C+(i<<3), CA);
                
                accumulator.val[0] = vaddq_f64 (CA.val[0], accumulator.val[0]);
                accumulator.val[1] = vaddq_f64 (CA.val[1], accumulator.val[1]);
                accumulator.val[2] = vaddq_f64 (CA.val[2], accumulator.val[2]);
                accumulator.val[3] = vaddq_f64 (CA.val[3], accumulator.val[3]);
            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = vaddvq_f64(vpaddq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]),vpaddq_f64(accumulator.val[2],accumulator.val[3])));
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) + (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) + (C[D-2] *= M[D-2]) + (C[D-1] *= M[D-1]);
                    break;
                case 4: {
                    float64x2x2_t CA = vld1q_f64_x2 (C+D-4),
                    MA = vld1q_f64_x2 (M+D-4);
                    
                    CA.val[0] = vmulq_f64 (CA.val[0],MA.val[0]);
                    CA.val[1] = vmulq_f64 (CA.val[1],MA.val[1]);
                    
                    vst1q_f64_x2 (C+D-4,CA);
                    
                    sum += vaddvq_f64(vpaddq_f64(CA.val[0],CA.val[1]));
                    break;
                }
                case 5: {
                    
                    float64x2x2_t CA = vld1q_f64_x2 (C+D-5),
                    MA = vld1q_f64_x2 (M+D-5);
                    
                    CA.val[0] = vmulq_f64 (CA.val[0],MA.val[0]);
                    CA.val[1] = vmulq_f64 (CA.val[1],MA.val[1]);
                    
                    vst1q_f64_x2 (C+D-5,CA);
                    
                    sum += (C[D-1] *= M[D-1]) + vaddvq_f64(vpaddq_f64(CA.val[0],CA.val[1]));
                    break;
                }
                case 6: {
                    float64x2x2_t CA  = vld1q_f64_x2 (C+D-6),
                    MA  = vld1q_f64_x2 (M+D-6);
                    
                    float64x2_t  CA2 = vld1q_f64 (C+D-2),
                    MA2 = vld1q_f64 (M+D-2);
                    
                    CA.val[0] = vmulq_f64 (CA.val[0],MA.val[0]);
                    CA.val[1] = vmulq_f64 (CA.val[1],MA.val[1]);
                    CA2 = vmulq_f64 (CA2,MA2);
                    
                    vst1q_f64_x2 (C+D-6,CA);
                    vst1q_f64 (C+D-2,CA2);
                    
                    sum += vaddvq_f64(vpaddq_f64(CA.val[0],CA.val[1])) + vaddvq_f64 (CA2);
                    break;
                }
                case 7:
                    sum += (C[D-7] *= M[D-7])
                    + (C[D-6] *= M[D-6])
                    + (C[D-5] *= M[D-5])
                    + (C[D-4] *= M[D-4])
                    + (C[D-3] *= M[D-3])
                    + (C[D-2] *= M[D-2])
                    + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

    double _hy_vvmult_sum_generic (double * C, double const *M, int D) {
        double sum = 0.;
        float64x2x2_t accumulator;
        
        int blocks = D>>3;
        int remainder = D - (blocks<<3);
        
        if (blocks) {
            float64x2x4_t CA = vld1q_f64_x4 (C),
            MA = vld1q_f64_x4 (M),
            accumulator;
            
            accumulator.val[0] = vmulq_f64 (CA.val[0], MA.val[0]);
            accumulator.val[1] = vmulq_f64 (CA.val[1], MA.val[1]);
            accumulator.val[2] = vmulq_f64 (CA.val[2], MA.val[2]);
            accumulator.val[3] = vmulq_f64 (CA.val[3], MA.val[3]);

            vst1q_f64_x4 (C, accumulator);
            
            for (int i = 1; i < blocks; i++) {
                CA = vld1q_f64_x4 (C+(i<<3));
                MA = vld1q_f64_x4 (M+(i<<3));
                
                CA.val[0] = vmulq_f64 (CA.val[0], MA.val[0]);
                CA.val[1] = vmulq_f64 (CA.val[1], MA.val[1]);
                CA.val[2] = vmulq_f64 (CA.val[2], MA.val[2]);
                CA.val[3] = vmulq_f64 (CA.val[3], MA.val[3]);

                vst1q_f64_x4 (C+(i<<3), CA);
                
                accumulator.val[0] = vaddq_f64 (CA.val[0], accumulator.val[0]);
                accumulator.val[1] = vaddq_f64 (CA.val[1], accumulator.val[1]);
                accumulator.val[2] = vaddq_f64 (CA.val[2], accumulator.val[2]);
                accumulator.val[3] = vaddq_f64 (CA.val[3], accumulator.val[3]);
            }
            
            //accumulator.val[0] = vaddq_f64 (accumulator.val[0], accumulator.val[1]);
            sum = vaddvq_f64(vpaddq_f64(vpaddq_f64(accumulator.val[0],accumulator.val[1]),vpaddq_f64(accumulator.val[2],accumulator.val[3])));
        }
        
        if (remainder) {
            switch (remainder) {
                case 1:
                    sum += (C[D-1] *= M[D-1]);
                    break;
                case 2:
                    sum += (C[D-2] *= M[D-2]) + (C[D-1] *= M[D-1]);
                    break;
                case 3:
                    sum += (C[D-3] *= M[D-3]) + (C[D-2] *= M[D-2]) + (C[D-1] *= M[D-1]);
                    break;
                case 4: {
                    float64x2x2_t CA = vld1q_f64_x2 (C+D-4),
                    MA = vld1q_f64_x2 (M+D-4);
                    
                    CA.val[0] = vmulq_f64 (CA.val[0],MA.val[0]);
                    CA.val[1] = vmulq_f64 (CA.val[1],MA.val[1]);
                    
                    vst1q_f64_x2 (C+D-4,CA);
                    
                    sum += vaddvq_f64(vpaddq_f64(CA.val[0],CA.val[1]));
                    break;
                }
                case 5: {
                    
                    float64x2x2_t CA = vld1q_f64_x2 (C+D-5),
                    MA = vld1q_f64_x2 (M+D-5);
                    
                    CA.val[0] = vmulq_f64 (CA.val[0],MA.val[0]);
                    CA.val[1] = vmulq_f64 (CA.val[1],MA.val[1]);
                    
                    vst1q_f64_x2 (C+D-5,CA);
                    
                    sum += (C[D-1] *= M[D-1]) + vaddvq_f64(vpaddq_f64(CA.val[0],CA.val[1]));
                    break;
                }
                case 6: {
                    float64x2x2_t CA  = vld1q_f64_x2 (C+D-6),
                    MA  = vld1q_f64_x2 (M+D-6);
                    
                    float64x2_t  CA2 = vld1q_f64 (C+D-2),
                    MA2 = vld1q_f64 (M+D-2);
                    
                    CA.val[0] = vmulq_f64 (CA.val[0],MA.val[0]);
                    CA.val[1] = vmulq_f64 (CA.val[1],MA.val[1]);
                    CA2 = vmulq_f64 (CA2,MA2);
                    
                    vst1q_f64_x2 (C+D-6,CA);
                    vst1q_f64 (C+D-2,CA2);
                    
                    sum += vaddvq_f64(vpaddq_f64(CA.val[0],CA.val[1])) + vaddvq_f64 (CA2);
                    break;
                }
                case 7:
                    sum += (C[D-7] *= M[D-7])
                    + (C[D-6] *= M[D-6])
                    + (C[D-5] *= M[D-5])
                    + (C[D-4] *= M[D-4])
                    + (C[D-3] *= M[D-3])
                    + (C[D-2] *= M[D-2])
                    + (C[D-1] *= M[D-1]);
                    break;
            }
        }
        
        return sum;
    }

    
#endif

#endif




inline void __ll_loop_handle_leaf_generic (hyFloat* _hprestrict_ pp, hyFloat *  _hprestrict_ localScalingFactor , long siteFrom, long siteTo, _SimpleList&        siteOrdering, bool matchSet, long * _hprestrict_ setBranchTo, long D) {
    
    if (matchSet) {
        memset (pp, 0, (siteTo-siteFrom) * D * sizeof (hyFloat));
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
    //tcc = nil;
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
    
    hyFloat * mvs = (hyFloat*)alloca (sizeof (hyFloat) * alphabetDimension);
        
    //if (setBranch >=0 )
       // printf ("\nSet to %d (%s)\n", setBranch, ((_CalcNode*) flatTree    (setBranch))->GetName()->get_str());
    
    /*if (likeFuncEvalCallCount == 0) {
        fprintf (stderr, "\nSite ID: %ld\n%s\n%s\n%s\n", siteOrdering.get(0), theFilter->GetColumn(siteOrdering.get(0)), theFilter->GetColumn(siteOrdering.get(0)+1), theFilter->GetColumn(siteOrdering.get(0)+2));
    }*/
    
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
        
        //printf ("\n***\nNode %d = %s\n", nodeID, currentTreeNode->GetName()->get_str());
        
        hyFloat  *  _hprestrict_ parentConditionals = iNodeCache +            (siteFrom + parentCode  * siteCount) * alphabetDimension;
        if (taggedInternals.list_data[parentCode] == 0) {
            // mark the parent for update and clear its conditionals if needed
            taggedInternals.list_data[parentCode]     = 1;
            hyFloat    *  _hprestrict_ localScalingFactor      = scalingAdjustments + parentCode*siteCount;
            
      
            
            bool    matchSet   = (parentCode == setBranch);
            
            //printf ("At %s, set to %d (%d)\n", currentTreeNode->GetName()->get_str(), parentCode, matchSet);
            
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
            
            /*if (matchSet) {
                for (long i = 0; i < alphabetDimension; i++) {
                    printf ("%g\t", parentConditionals[i]);
                }
                printf ("\n");
            }*/
        }
        
        _Matrix  const * transitionMatrixObj           = currentTreeNode->GetCompExp(catID);
        hyFloat  const * _hprestrict_ transitionMatrix = transitionMatrixObj->theData;
        
        
         /*if (likeFuncEvalCallCount == 52) {
            //if (currentTreeNode->GetName()->EndsWith("mt811400_SARS2_orf1ab_usa__4")) {
            fprintf (stderr, "\nBRANCH ID %ld (%ld = parent ID, parent name = %s) (%s)\n", nodeCode, parentCode, ((_CalcNode*)flatTree (parentCode))->GetName()->get_str(), currentTreeNode->GetName()->get_str());
                for (long e = 0; e < 4; e++) {
                    fprintf (stderr, "\n");
                    for (long e2 = 0; e2 < 4; e2++) {
                    fprintf (stderr, "%18.12g\t", transitionMatrix[e*4+e2]);
                }
            }
            fprintf (stderr, "\n");
        }*/
         
        
        long currentTCCIndex,currentTCCBit,parentTCCIIndex,parentTCCIBit;
        __ll_handle_tcc_init (tcc, isLeaf, siteCount, siteFrom, nodeCode, parentCode, parentTCCIBit, parentTCCIIndex, currentTCCBit, currentTCCIndex);
        
        //long successiveSkips = 0;
        /**
          20200929: SLKP trying to refactor this code.
          handle specialized used cases
         */
        
// BEGIN NUCLEOTIDE CASE
        if (alphabetDimension == 4UL) {
            
            #ifdef _HY_NO_INTRINSICS_MARKER
                hyFloat * tmatrix_transpose = (hyFloat *)transitionMatrix;
            #else
            
                #ifdef _SLKP_USE_AVX_INTRINSICS
                    __m256d tmatrix_transpose [4] = {
                        (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                        (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                        (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                        (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                    };
                #endif
                
                #ifdef _SLKP_USE_ARM_NEON
                    float64x2x2_t tmatrix_transpose [4] = {
                        (float64x2x2_t) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                        (float64x2x2_t) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                        (float64x2x2_t) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                        (float64x2x2_t) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                    };
                #endif
                
                #ifdef _SLKP_USE_SSE_INTRINSICS
                    __m128d tmatrix_transpose [8] = {(__m128d) {transitionMatrix[0],transitionMatrix[4]}, (__m128d) {transitionMatrix[8],transitionMatrix[12]},
                    (__m128d) {transitionMatrix[1],transitionMatrix[5]}, (__m128d) {transitionMatrix[9],transitionMatrix[13]},
                    (__m128d) {transitionMatrix[2],transitionMatrix[6]}, (__m128d) {transitionMatrix[10],transitionMatrix[14]},
                    (__m128d) {transitionMatrix[3],transitionMatrix[7]}, (__m128d) {transitionMatrix[11],transitionMatrix[15]}
                    };
                #endif
        #endif
            
            
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 4UL) {
                __ll_loop_preamble
                if (__ll_handle_conditional_array_initialization<4> (
                                                                      lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                    continue;
                }
                
                /*
                #if defined _SLKP_USE_AVX_INTRINSICS or defined _SLKP_USE_ARM_NEON
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, tmatrix_transpose);
                #else
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, nil);
                #endif
                 
                sum     = (parentConditionals [0] + parentConditionals [1]) + (parentConditionals [2] + parentConditionals [3]);
                */
                
                /*if (likeFuncEvalCallCount == 52 && siteID == 7) {
                   //if (currentTreeNode->GetName()->EndsWith("mt811400_SARS2_orf1ab_usa__4")) {
                   fprintf (stderr, "\nBEFORE Parent conditionals (%ld = parent ID, parent name = %s) (%s)\n", parentCode, ((_CalcNode*)flatTree (parentCode))->GetName()->get_str(), currentTreeNode->GetName()->get_str());
                   for (long e = 0; e < 4; e++) {
                        fprintf (stderr, "%18.12g\t", parentConditionals[e]);
                   }
                   fprintf (stderr, "\n");
               }*/

                
                sum = _handle4x4_pruning_case_direct (childVector, tmatrix_transpose, parentConditionals);
                
                __ll_loop_handle_scaling<4L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                
                
                /*if (likeFuncEvalCallCount == 52 && siteID == 7) {
                   //if (currentTreeNode->GetName()->EndsWith("mt811400_SARS2_orf1ab_usa__4")) {
                   fprintf (stderr, "\nAFTER Parent conditionals (%ld = parent ID, parent name = %s) (%s), sum = %g\n", parentCode, ((_CalcNode*)flatTree (parentCode))->GetName()->get_str(), currentTreeNode->GetName()->get_str(), sum);
                   for (long e = 0; e < 4; e++) {
                        fprintf (stderr, "%18.12g/%18.12g\t", childVector[e], parentConditionals[e]);
                   }
                   fprintf (stderr, "\n");
               }*/

                
                childVector += 4L;
                __ll_loop_epilogue
            }
            
        }
// END NUCLEOTIDE CASE
        
// START AMINO-ACID CASE
        else if (alphabetDimension == 20UL) {
            
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 20L) {
                __ll_loop_preamble
                if (__ll_handle_conditional_array_initialization<20> (
                                                                      lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                    continue;
                }

                
                /*if (siteID == 1 && likeFuncEvalCallCount >= 1263 && parentCode == 105) {
                    printf ("\n@%ld (%d : %d)\n", likeFuncEvalCallCount, parentCode, nodeCode);
                    for (int k = 0; k < 20; k++)
                        printf ("\t %d %g %g\n", k, parentConditionals[k], childVector[k]);
                }*/
                
                #ifdef _SLKP_USE_APPLE_BLAS_2
                         cblas_dgemv(CblasRowMajor,
                                     CblasNoTrans,
                                     20,
                                     20,
                                     1.,
                                     tMatrix,
                                     20,
                                     childVector,
                                     1,
                                     0.,
                                     mvs,
                                     1);
                #else
                        _hy_mvp_blocked_4x4<20> (mvs, tMatrix, childVector);
                #endif
                sum += _hy_vvmult_sum<20> (parentConditionals, mvs);

                __ll_loop_handle_scaling<20L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

                childVector += 20;
                __ll_loop_epilogue
            }

    }
    // END AMINO-ACID CASE
    // START UNIVERSAL CODE
    else if (alphabetDimension == 60UL) {
             for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 60L) {
                __ll_loop_preamble
                if (__ll_handle_conditional_array_initialization<60> (
                                                                      lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                    continue;
                }

                
                _hy_mvp_blocked_4x4<60> (mvs, tMatrix, childVector);
                sum += _hy_vvmult_sum<60> (parentConditionals, mvs);
            
                __ll_loop_handle_scaling<60L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

                childVector += 60;
                __ll_loop_epilogue
            }
    }
    else if (alphabetDimension == 61UL) {
        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 61L) {

/*
#pragma omp critical
            if (siteID == 524) {
                printf ("(%ld/%ld)%s\n", likeFuncEvalCallCount, siteFrom, currentTreeNode->GetName()->get_str());
            }
*/

                    
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization<61> (
                                                                  lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                                                                      
                                                                      
                continue;
            }


/*
#pragma omp critical
            if (1 || hy_env::EnvVariableTrue("UBER_VERBOSE_DEBUG")) {
                if (siteID == 524) {
                //if (siteOrdering.list_data[siteID] == 2 && nodeID == 75) {
                    //printf ("%ld\t%ld\t%ld\t%ld\t%16.12g\n",nodeID, siteFrom, siteTo, siteID, sum);
                    for (int kk = 0; kk < 61; kk++) {
                        printf ("%s (%ld)\t%ld\t%ld\t%d\t%16.12g\t%16.12g\t%16.12g\n", currentTreeNode->GetName()->get_str(),nodeID, siteFrom, siteID, kk, parentConditionals[kk], childVector[kk], mvs[kk]);
                    }
                }
            }
*/

        #ifdef _SLKP_USE_APPLE_BLAS_2
                 cblas_dgemv(CblasRowMajor,
                             CblasNoTrans,
                             61,
                             61,
                             1.,
                             tMatrix,
                             61,
                             childVector,
                             1,
                             0.,
                             mvs,
                             1);
        #else
                _hy_mvp_blocked_4x4<61> (mvs, tMatrix, childVector);
        #endif
        sum += _hy_vvmult_sum<61> (parentConditionals, mvs);

/*
#pragma omp critical
            if (1 || hy_env::EnvVariableTrue("UBER_VERBOSE_DEBUG")) {
                if (siteID == 524) {
                //if (siteOrdering.list_data[siteID] == 2 && nodeID == 75) {
                    //printf ("%ld\t%ld\t%ld\t%ld\t%16.12g\n",nodeID, siteFrom, siteTo, siteID, sum);
                    for (int kk = 0; kk < 61; kk++) {
                        printf ("%s (%d)\t%ld\t%ld\t%d\t%16.12g\t%16.12g\t%16.12g\n", currentTreeNode->GetName()->get_str(), nodeID, siteFrom, siteID, kk, parentConditionals[kk], childVector[kk], mvs[kk]);
                    }
                }
            }
*/


            
            __ll_loop_handle_scaling<61L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

            childVector += 61;
            __ll_loop_epilogue
        }
    } else if (alphabetDimension == 62UL) {
        #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS || defined _SLKP_USE_ARM_NEON
            //__ll_handle_matrix_transpose<62> (transitionMatrix, tMatrixT);
            //transitionMatrixObj->TransposeIntoStorage (tMatrixT);
        #endif
        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 62L) {
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization<62> (
                                                                  lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                continue;
            }

            

            #ifdef _SLKP_USE_APPLE_BLAS_2
                     cblas_dgemv(CblasRowMajor,
                                 CblasNoTrans,
                                 62,
                                 62,
                                 1.,
                                 tMatrix,
                                 62,
                                 childVector,
                                 1,
                                 0.,
                                 mvs,
                                 1);
            #else
                    _hy_mvp_blocked_4x4<62> (mvs, tMatrix, childVector);
            #endif
            sum += _hy_vvmult_sum<62> (parentConditionals, mvs);
            __ll_loop_handle_scaling<62L, true> (sum, parentConditionals, scalingAdjustments, didScale, parentCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

            childVector += 62;
            __ll_loop_epilogue
        }
    } else if (alphabetDimension == 63UL) {
        for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 63L) {
            __ll_loop_preamble
            if (__ll_handle_conditional_array_initialization<63> (
                                                                  lNodeFlags, isLeaf, nodeCode, setBranch, flatTree.lLength, siteID, siteFrom, siteCount, siteOrdering, parentConditionals, tMatrix, lNodeResolutions, childVector, tcc, currentTCCBit, currentTCCIndex, lastUpdatedSite, setBranchTo)) {
                continue;
            }
            
           
            #ifdef _SLKP_USE_APPLE_BLAS_2
                     cblas_dgemv(CblasRowMajor,
                                 CblasNoTrans,
                                 63,
                                 63,
                                 1.,
                                 tMatrix,
                                 63,
                                 childVector,
                                 1,
                                 0.,
                                 mvs,
                                 1);
            #else
                    _hy_mvp_blocked_4x4<63> (mvs, tMatrix, childVector);
            #endif
            
            sum += _hy_vvmult_sum<63> (parentConditionals, mvs);

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
            #ifdef _SLKP_USE_APPLE_BLAS_2
                 cblas_dgemv(CblasRowMajor,
                             CblasNoTrans,
                             alphabetDimension,
                             alphabetDimension,
                             1.,
                             tMatrix,
                             alphabetDimension,
                             childVector,
                             1,
                             0.,
                             mvs,
                             1);
            #else
                _hy_matrix_vector_product_blocked_4x4 (mvs, tMatrix, childVector, alphabetDimension);
            #endif
             
             //_hy_mvp_blocked_4x4<16> (mvs, tMatrix, childVector);
             sum += _hy_vvmult_sum_generic (parentConditionals, mvs, alphabetDimension);
             //sum += _hy_vvmult_sum<16> (parentConditionals, mvs);
             
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
            //#pragma unroll(4)
            #pragma GCC unroll 4
            for (long p = 0; p < alphabetDimension; p++,rootIndex++) {
                accumulator += rootConditionals[rootIndex] * theProbs[p];
            }
        }
      
/*
#pragma omp critical
        if (hy_env::EnvVariableTrue("UBER_VERBOSE_DEBUG")) {
            if (siteOrdering.list_data[siteID] == 2) {
                printf ("\n");
                printf ("%ld\t%ld\t%ld\t%16.12g\n",siteFrom, siteTo, siteID, accumulator);
            }
        }
*/

        if (storageVec) {
            storageVec [siteOrdering.list_data[siteID]] = accumulator;
            if (accumulator <= 0.0)
#pragma omp critical
                {
                    //printf ("BAILING WITH INFINITY %ld / %ld\n", siteID, siteOrdering.list_data[siteID]);
                    hy_global::ReportWarning (_String("Site ") & (1L+siteOrdering.list_data[siteID]) & " evaluated to a 0 probability in ComputeTreeBlockByBranch with category " & catID);
                }
       } else {
            if (accumulator <= 0.0) {
                result = -INFINITY;
#pragma omp critical
                {
                    //printf ("BAILING WITH INFINITY %ld / %ld (%ld)\n", siteID, siteOrdering.list_data[siteID], siteOrdering.lLength);
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
                
                /*if (likeFuncEvalCallCount == 328 && siteID == 30) {
                    fprintf (stderr, "\n%ld\t%18.14lg\t%18.14lg\t%18.14lg\t%ld\n", siteID, accumulator, term, result,  site_frequency   );
                    for (long k = 0; k < alphabetDimension; k++) {
                        fprintf (stderr, "\t %d = %16.12g %16.12g\n", k, rootConditionals[rootIndex-alphabetDimension+k], theProbs[k]);
                    }
                }*/

                hyFloat temp_sum = result + term;
                correction = (temp_sum - result) - term;
                result = temp_sum;
            } else {
                
                long lfID = -1;
                IsLinkedToALF (lfID);
                ((_LikelihoodFunction*)likeFuncList (lfID))->
                _TerminateAndDump(_String("Site ") & (1L+siteOrdering.list_data[siteID]) & " evaluated to a NaN probability in ComputeTreeBlockByBranch; this is not a recoverable error and indicates some serious COVFEFE taking place.");
            }

        }
    }
     
    /*if (setBranch >= 0) {
        
        for (long nn = 0; nn < flatTree.lLength ; nn++) {
            hyFloat * _hprestrict_ nc = iNodeCache + alphabetDimension * (siteFrom + (nn)  * siteCount);
            printf ("\nNode %s\n", ((_CalcNode*)flatTree.GetItem(nn))->GetName()->get_str());
            for (long p = 0; p < alphabetDimension; p++) {
                printf ("%10.6g\t", nc[p]);
            }
            printf ("\n\n");
        }
        
        //exit (0);
    }*/
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
            //#pragma unroll(4)
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
                //#pragma unroll(4)
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
            //#pragma unroll(4)
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
                //#pragma unroll(4)
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
             siteCorrectionCounts [siteOrdering.list_data[siteID]] += didScale;
            if (didScale == 1L) {
                siteRes[siteOrdering.list_data[siteID]] *= _lfScalerUpwards;
            } else {
                if (didScale == -1L) {
                    siteRes[siteOrdering.list_data[siteID]] *= _lfScalingFactorThreshold;
                } else {
                    if (didScale > 0) {
                        siteRes[siteOrdering.list_data[siteID]] *= ComputePower (_lfScalerUpwards, didScale);
                        /*for (long k = 0; k < didScale; k++) {
                            siteRes[siteOrdering.list_data[siteID]] *= _lfScalerUpwards;
                        }*/
                    } else{
                        siteRes[siteOrdering.list_data[siteID]] *= ComputePower (_lfScalingFactorThreshold, -didScale);
                        /*for (long k = 0; k < -didScale; k++) {
                             siteRes[siteOrdering.list_data[siteID]] *= _lfScalingFactorThreshold;
                         }*/
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
    //alphabetDimensionmod4  =         alphabetDimension - alphabetDimension % 4,
    siteCount               =            theFilter->GetPatternCount();

    if (siteTo  > siteCount)    {
        siteTo = siteCount;
    }
    
    do {
        taggedNodes.list_data[myParent+flatLeaves.lLength] = 1;
        myParent = flatParents.list_data[myParent+flatLeaves.lLength];
    } while (myParent >= 0);
    
    hyFloat * mvs = (hyFloat*)alloca (sizeof (hyFloat) * alphabetDimension);

    
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
    /*
    if (likeFuncEvalCallCount == 3013) {
        fprintf (stderr, "\nCaching branch %ld\n", brID);
        rootPath.Each ([&](long n, unsigned long i) -> void {
            _CalcNode * current_node = (_CalcNode*) (n < flatLeaves.lLength ?  flatCLeaves (n):
                                                                              flatTree    (n-flatLeaves.lLength));
            fprintf (stderr, "[%d: %d] %s\n", i, n,current_node->GetName()->get_str());
        });
        ObjectToConsole(&nodesToProcess);NLToConsole();
    }
    */
    
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
                    //#pragma unroll(4)
                    #pragma GCC unroll 4
                    for (unsigned long k2 = 0UL; k2 < alphabetDimension; k2++, k3++) {
                        parentConditionals [k3] = scaler;
                    }
                }
            }
        }
        
        _CalcNode    * currentTreeNode = (_CalcNode*) (isLeaf?  flatCLeaves (nodeCode):
                                                       flatTree    (notPassedRoot?nodeCode:parentCode));
        
        //if (likeFuncEvalCallCount == 3013) {
        //    fprintf (stderr, "%ld/%ld (%ld/%ld) => %s\n", nodeID, nodeCode, isLeaf, notPassedRoot, currentTreeNode->GetName()->get_str());
        //}
        
        _Matrix  const * transitionMatrixObj           = currentTreeNode->GetCompExp(catID);
        hyFloat  const *  transitionMatrix = transitionMatrixObj->theData;
        hyFloat  *       childVector,*     lastUpdatedSite;
        
        if (!isLeaf) {
            lastUpdatedSite = childVector = iNodeCache + (siteFrom + nodeCode * siteCount) * alphabetDimension;
        }
        
        
        long currentTCCIndex,currentTCCBit,parentTCCIIndex,parentTCCIBit;
        __ll_handle_tcc_init (tcc, isLeaf, siteCount, siteFrom, nodeCode, parentCode, parentTCCIBit, parentTCCIIndex, currentTCCBit, currentTCCIndex);

 
        if (alphabetDimension == 4L) {
            #ifdef _HY_NO_INTRINSICS_MARKER
                hyFloat * tmatrix_transpose = (hyFloat *)transitionMatrix;
            #else

                #ifdef _SLKP_USE_AVX_INTRINSICS
                    __m256d tmatrix_transpose [4] = {
                        (__m256d) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                        (__m256d) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                        (__m256d) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                        (__m256d) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                    };
                #endif
                
                #ifdef _SLKP_USE_ARM_NEON
                    float64x2x2_t tmatrix_transpose [4] = {
                        (float64x2x2_t) {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
                        (float64x2x2_t) {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
                        (float64x2x2_t) {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
                        (float64x2x2_t) {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
                    };
                #endif
                
                #ifdef _SLKP_USE_SSE_INTRINSICS
                    __m128d tmatrix_transpose [8] = {(__m128d) {transitionMatrix[0],transitionMatrix[4]}, (__m128d) {transitionMatrix[8],transitionMatrix[12]},
                    (__m128d) {transitionMatrix[1],transitionMatrix[5]}, (__m128d) {transitionMatrix[9],transitionMatrix[13]},
                    (__m128d) {transitionMatrix[2],transitionMatrix[6]}, (__m128d) {transitionMatrix[10],transitionMatrix[14]},
                    (__m128d) {transitionMatrix[3],transitionMatrix[7]}, (__m128d) {transitionMatrix[11],transitionMatrix[15]}
                    };
                #endif
            #endif
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 4L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<4>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    /*if (likeFuncEvalCallCount == 328 && siteID == 30) {
                        fprintf (stderr, "__lcache_loop_preface (%ld isleaf = %d canscale = %d, %s) %g %g %g %g\n", nodeCode, isLeaf, canScale, currentTreeNode->GetName()->get_str(), parentConditionals[0], parentConditionals[1], parentConditionals[2], parentConditionals[3]);
                    }*/
                    continue;
                }
                long     didScale =  0;

                /*#if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_ARM_NEON
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, tmatrix_transpose);
                #else
                    _handle4x4_pruning_case (childVector, tMatrix, parentConditionals, nil);
                #endif
                if (canScale) {
                    hyFloat sum     = (parentConditionals [0] + parentConditionals [1]) + (parentConditionals [2] + parentConditionals [3]);
                    
                    __ll_loop_handle_scaling<4L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

                }*/
                
                hyFloat sum     = _handle4x4_pruning_case_direct (childVector, tmatrix_transpose, parentConditionals);

                if (canScale) {
                    
                    __ll_loop_handle_scaling<4L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));

                }
                
                /*if (likeFuncEvalCallCount == 328 && siteID == 30) {
                    fprintf (stderr, "NODE = %ld, PARENT = %ld (%ld), P(A) = %lg, P(C) = %lg, P(G) = %lg, P(T) = %lg, scale = %ld\n", nodeCode, parentCode, canScale, parentConditionals[0], parentConditionals[1], parentConditionals[2], parentConditionals[3], didScale);
                }*/
                
                childVector += 4L;
                __handle_site_corrections(didScale, siteID);
            }
        } else if (alphabetDimension == 20L) {
                for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 20L) {
                    bool canScale = !notPassedRoot;
                    hyFloat  const *tMatrix = transitionMatrix;
                    if (__lcache_loop_preface<20>(
                    isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                        continue;
                    }
                    long     didScale =  0;
                    hyFloat sum     = 0.;
                    _hy_mvp_blocked_4x4<20> (mvs, tMatrix, childVector);
                    sum = _hy_vvmult_sum<20> (parentConditionals, mvs);
                    if (canScale) {
                        __ll_loop_handle_scaling<20L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                    }
                    childVector += 20L;
                    __handle_site_corrections(didScale, siteID);
                }
        } else if (alphabetDimension == 60L) {
            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 60L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<60>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    continue;
                }
                long     didScale =  0;
                hyFloat sum     = 0.;
                _hy_mvp_blocked_4x4<60> (mvs, tMatrix, childVector);
                sum = _hy_vvmult_sum<60> (parentConditionals, mvs);
                if (canScale) {
                    __ll_loop_handle_scaling<60L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 60L;
                __handle_site_corrections(didScale, siteID);
            }
        } else if (alphabetDimension == 61L) {

            for (long siteID = siteFrom; siteID < siteTo; siteID++, parentConditionals += 61L) {
                bool canScale = !notPassedRoot;
                hyFloat  const *tMatrix = transitionMatrix;
                if (__lcache_loop_preface<61>(
                isLeaf, lNodeFlags, siteID, siteOrdering, nodeCode, siteCount, siteFrom, parentConditionals, tMatrix, canScale, childVector, lastUpdatedSite, tcc, currentTCCBit, currentTCCIndex, parentTCCIBit, parentTCCIIndex, notPassedRoot, lNodeResolutions)) {
                    continue;
                }
                long     didScale =  0;
  
                _hy_mvp_blocked_4x4<61> (mvs, tMatrix, childVector);
                hyFloat sum = _hy_vvmult_sum<61> (parentConditionals, mvs);
                
                if (canScale) {
                    __ll_loop_handle_scaling<61L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 61L;
                __handle_site_corrections(didScale, siteID);
            }
        }  else if (alphabetDimension == 62L) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS || defined _SLKP_USE_ARM_NEON
                //__ll_handle_matrix_transpose<62L>(transitionMatrix, tMatrixT);
                //transitionMatrixObj->TransposeIntoStorage (tMatrixT);
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

                _hy_mvp_blocked_4x4<62> (mvs, tMatrix, childVector);
                sum = _hy_vvmult_sum<62> (parentConditionals, mvs);

                if (canScale) {
                    __ll_loop_handle_scaling<62L, false> (sum, parentConditionals, scalingAdjustments, didScale, nodeCode, siteCount, siteID, localScalerChange, theFilter->theFrequencies.get (siteOrdering.list_data[siteID]));
                }
                childVector += 62L;
                __handle_site_corrections(didScale, siteID);
            }
        } else if (alphabetDimension == 63L) {
            #if defined _SLKP_USE_AVX_INTRINSICS || defined _SLKP_USE_SSE_INTRINSICS || defined _SLKP_USE_ARM_NEON
                //__ll_handle_matrix_transpose<63L>(transitionMatrix, tMatrixT);
                //transitionMatrixObj->TransposeIntoStorage (tMatrixT);
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
                _hy_mvp_blocked_4x4<63> (mvs, tMatrix, childVector);
                sum = _hy_vvmult_sum<63> (parentConditionals, mvs);
                if (canScale) {
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
                _hy_matrix_vector_product_blocked_4x4 (mvs, tMatrix, childVector, alphabetDimension);
                sum = _hy_vvmult_sum_generic (parentConditionals, mvs, alphabetDimension);
   
                if (canScale) {
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
        /*if (likeFuncEvalCallCount == 328) {
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
