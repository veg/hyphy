#include <Accelerate/Accelerate.h>
#include <arm_neon.h>
#include <benchmark/benchmark.h>
#include <cstdio>

typedef double hyFloat;
#define _hprestrict_ __restrict

#define N 4

// #define _SLKP_USE_ARM_NEON

//_____________________________________________________________________________________________
hyFloat HYPHY_SQR(hyFloat *_hprestrict_ theData, hyFloat *_hprestrict_ stash,
                  const long hDim, const long vDim, const long lDim) {

  hyFloat diff = 0.;
  if (hDim != vDim) {
    return diff;
  }
  // not a square matrix

  if (false) {
    return DBL_EPSILON * 1.e4;
  } else {
    if (hDim == 4L) {
      // special case for nucleotides
      for (unsigned long i = 0UL, k = 0UL; i < 16; i += 4) {
#pragma unroll 4
        for (unsigned long j = 0UL; j < 4UL; j++, k++) {
          hyFloat p1 = theData[i] * theData[j];
          hyFloat p2 = theData[i + 1] * theData[j + 4];
          p1 += theData[i + 2] * theData[j + 8];
          p2 += theData[i + 3] * theData[j + 12];
          stash[k] = p1 + p2;
        }
      }
    } else {
      long loopBound = (vDim >> 2) << 2;
      // vDim - vDim % 4;
      //  loop interchange rocks!

#ifdef _SLKP_USE_ARM_NEON
#define DO_GROUP_OP0(X, Y, k)                                                  \
  Y = vld1q_f64(theData + col_offset + k);                                     \
  X = vmulq_f64(A4, Y);
#define DO_GROUP_OP1(X, Y, k)                                                  \
  X = vld1q_f64(dest + row_offset + k);                                        \
  Y = vld1q_f64(theData + col_offset + k);                                     \
  X = vfmaq_f64(X, A4, Y);
#define DO_GROUP_OP2(X, k) vst1q_f64(dest + row_offset + k, X);

      if (true) {

        hyFloat *_hprestrict_ dest = stash;

        long ti = 0L, row_offset = 0L;

        if (loopBound == 60UL) {
          if (hDim == 61) { // special case for universal code

            for (long r = 0; r < 61; r++, row_offset += 61) {
              long col_offset = 0L;
              { // row 1
                float64x2_t A4 = vdupq_n_f64(theData[ti]);
                float64x2_t X1, X2, X3, X4, Y1, Y2, Y3, Y4;

                DO_GROUP_OP0(X1, Y1, 0);
                DO_GROUP_OP0(X2, Y2, 2);
                DO_GROUP_OP0(X3, Y3, 4);
                DO_GROUP_OP0(X4, Y4, 6);
                DO_GROUP_OP2(X1, 0);
                DO_GROUP_OP2(X2, 2);
                DO_GROUP_OP2(X3, 4);
                DO_GROUP_OP2(X4, 6);
                DO_GROUP_OP0(X1, Y1, 8);
                DO_GROUP_OP0(X2, Y2, 10);
                DO_GROUP_OP0(X3, Y3, 12);
                DO_GROUP_OP0(X4, Y4, 14);
                DO_GROUP_OP2(X1, 8);
                DO_GROUP_OP2(X2, 10);
                DO_GROUP_OP2(X3, 12);
                DO_GROUP_OP2(X4, 14);
                DO_GROUP_OP0(X1, Y1, 16);
                DO_GROUP_OP0(X2, Y2, 18);
                DO_GROUP_OP0(X3, Y3, 20);
                DO_GROUP_OP0(X4, Y4, 22);
                DO_GROUP_OP2(X1, 16);
                DO_GROUP_OP2(X2, 18);
                DO_GROUP_OP2(X3, 20);
                DO_GROUP_OP2(X4, 22);
                DO_GROUP_OP0(X1, Y1, 24);
                DO_GROUP_OP0(X2, Y2, 26);
                DO_GROUP_OP0(X3, Y3, 28);
                DO_GROUP_OP0(X4, Y4, 30);
                DO_GROUP_OP2(X1, 24);
                DO_GROUP_OP2(X2, 26);
                DO_GROUP_OP2(X3, 28);
                DO_GROUP_OP2(X4, 30);
                DO_GROUP_OP0(X1, Y1, 32);
                DO_GROUP_OP0(X2, Y2, 34);
                DO_GROUP_OP0(X3, Y3, 36);
                DO_GROUP_OP0(X4, Y4, 38);
                DO_GROUP_OP2(X1, 32);
                DO_GROUP_OP2(X2, 34);
                DO_GROUP_OP2(X3, 36);
                DO_GROUP_OP2(X4, 38);
                DO_GROUP_OP0(X1, Y1, 40);
                DO_GROUP_OP0(X2, Y2, 42);
                DO_GROUP_OP0(X3, Y3, 44);
                DO_GROUP_OP0(X4, Y4, 46);
                DO_GROUP_OP2(X1, 40);
                DO_GROUP_OP2(X2, 42);
                DO_GROUP_OP2(X3, 44);
                DO_GROUP_OP2(X4, 46);
                DO_GROUP_OP0(X1, Y1, 48);
                DO_GROUP_OP0(X2, Y2, 50);
                DO_GROUP_OP0(X3, Y3, 52);
                DO_GROUP_OP0(X4, Y4, 54);
                DO_GROUP_OP2(X1, 48);
                DO_GROUP_OP2(X2, 50);
                DO_GROUP_OP2(X3, 52);
                DO_GROUP_OP2(X4, 54);
                DO_GROUP_OP0(X1, Y1, 56);
                DO_GROUP_OP0(X2, Y2, 58);
                DO_GROUP_OP2(X1, 56);
                DO_GROUP_OP2(X2, 58);

                dest[row_offset + 60] = theData[ti] * theData[60];
                col_offset = 61;
                ti++;
              }

#define DO_GROUP_OP_4_1(X, Y, k)                                               \
  X = vld1q_f64_x4(dest + row_offset + k);                                     \
  Y = vld1q_f64_x4(theData + col_offset + k);                                  \
  X.val[0] = vfmaq_f64(X.val[0], A4, Y.val[0]);                                \
  X.val[1] = vfmaq_f64(X.val[1], A4, Y.val[1]);                                \
  X.val[2] = vfmaq_f64(X.val[2], A4, Y.val[2]);                                \
  X.val[3] = vfmaq_f64(X.val[3], A4, Y.val[3]);

#define DO_GROUP_OP_4_2(X, k) vst1q_f64_x4(dest + row_offset + k, X);

              for (long c = 1; c < 61; c++, ti++, col_offset += 61) {
                float64x2_t A4 = vdupq_n_f64(theData[ti]);
                /*float64x2x4_t X,Y;

                DO_GROUP_OP_4_1 (X,Y,0); DO_GROUP_OP_4_2 (X,0);
                DO_GROUP_OP_4_1 (X,Y,8); DO_GROUP_OP_4_2 (X,8);
                DO_GROUP_OP_4_1 (X,Y,16); DO_GROUP_OP_4_2 (X,16);
                DO_GROUP_OP_4_1 (X,Y,24); DO_GROUP_OP_4_2 (X,24);
                DO_GROUP_OP_4_1 (X,Y,32); DO_GROUP_OP_4_2 (X,32);
                DO_GROUP_OP_4_1 (X,Y,40); DO_GROUP_OP_4_2 (X,40);
                DO_GROUP_OP_4_1 (X,Y,48); DO_GROUP_OP_4_2 (X,48);
                DO_GROUP_OP1 (X.val[0], Y.val[0], 56); DO_GROUP_OP1 (X.val[1],
                Y.val[1], 58); DO_GROUP_OP2 (X.val[0],56); DO_GROUP_OP2
                (X.val[1],58);*/

                float64x2_t X1, X2, X3, X4, Y1, Y2, Y3, Y4;

                DO_GROUP_OP1(X1, Y1, 0);
                DO_GROUP_OP1(X2, Y2, 2);
                DO_GROUP_OP1(X3, Y3, 4);
                DO_GROUP_OP1(X4, Y4, 6);
                DO_GROUP_OP2(X1, 0);
                DO_GROUP_OP2(X2, 2);
                DO_GROUP_OP2(X3, 4);
                DO_GROUP_OP2(X4, 6);
                DO_GROUP_OP1(X1, Y1, 8);
                DO_GROUP_OP1(X2, Y2, 10);
                DO_GROUP_OP1(X3, Y3, 12);
                DO_GROUP_OP1(X4, Y4, 14);
                DO_GROUP_OP2(X1, 8);
                DO_GROUP_OP2(X2, 10);
                DO_GROUP_OP2(X3, 12);
                DO_GROUP_OP2(X4, 14);
                DO_GROUP_OP1(X1, Y1, 16);
                DO_GROUP_OP1(X2, Y2, 18);
                DO_GROUP_OP1(X3, Y3, 20);
                DO_GROUP_OP1(X4, Y4, 22);
                DO_GROUP_OP2(X1, 16);
                DO_GROUP_OP2(X2, 18);
                DO_GROUP_OP2(X3, 20);
                DO_GROUP_OP2(X4, 22);
                DO_GROUP_OP1(X1, Y1, 24);
                DO_GROUP_OP1(X2, Y2, 26);
                DO_GROUP_OP1(X3, Y3, 28);
                DO_GROUP_OP1(X4, Y4, 30);
                DO_GROUP_OP2(X1, 24);
                DO_GROUP_OP2(X2, 26);
                DO_GROUP_OP2(X3, 28);
                DO_GROUP_OP2(X4, 30);
                DO_GROUP_OP1(X1, Y1, 32);
                DO_GROUP_OP1(X2, Y2, 34);
                DO_GROUP_OP1(X3, Y3, 36);
                DO_GROUP_OP1(X4, Y4, 38);
                DO_GROUP_OP2(X1, 32);
                DO_GROUP_OP2(X2, 34);
                DO_GROUP_OP2(X3, 36);
                DO_GROUP_OP2(X4, 38);
                DO_GROUP_OP1(X1, Y1, 40);
                DO_GROUP_OP1(X2, Y2, 42);
                DO_GROUP_OP1(X3, Y3, 44);
                DO_GROUP_OP1(X4, Y4, 46);
                DO_GROUP_OP2(X1, 40);
                DO_GROUP_OP2(X2, 42);
                DO_GROUP_OP2(X3, 44);
                DO_GROUP_OP2(X4, 46);
                DO_GROUP_OP1(X1, Y1, 48);
                DO_GROUP_OP1(X2, Y2, 50);
                DO_GROUP_OP1(X3, Y3, 52);
                DO_GROUP_OP1(X4, Y4, 54);
                DO_GROUP_OP2(X1, 48);
                DO_GROUP_OP2(X2, 50);
                DO_GROUP_OP2(X3, 52);
                DO_GROUP_OP2(X4, 54);
                DO_GROUP_OP1(X1, Y1, 56);
                DO_GROUP_OP1(X2, Y2, 58);
                DO_GROUP_OP2(X1, 56);
                DO_GROUP_OP2(X2, 58);

                dest[row_offset + 60] += theData[ti] * theData[col_offset + 60];
              }
            }
          } else {
            for (long r = 0; r < hDim; r++, row_offset += vDim) {
              long col_offset = 0L;
              { // row 1
                float64x2_t A4 = vdupq_n_f64(theData[ti]);
                float64x2_t X1, X2, X3, X4, Y1, Y2, Y3, Y4;

                DO_GROUP_OP0(X1, Y1, 0);
                DO_GROUP_OP0(X2, Y2, 2);
                DO_GROUP_OP0(X3, Y3, 4);
                DO_GROUP_OP0(X4, Y4, 6);
                DO_GROUP_OP2(X1, 0);
                DO_GROUP_OP2(X2, 2);
                DO_GROUP_OP2(X3, 4);
                DO_GROUP_OP2(X4, 6);
                DO_GROUP_OP0(X1, Y1, 8);
                DO_GROUP_OP0(X2, Y2, 10);
                DO_GROUP_OP0(X3, Y3, 12);
                DO_GROUP_OP0(X4, Y4, 14);
                DO_GROUP_OP2(X1, 8);
                DO_GROUP_OP2(X2, 10);
                DO_GROUP_OP2(X3, 12);
                DO_GROUP_OP2(X4, 14);
                DO_GROUP_OP0(X1, Y1, 16);
                DO_GROUP_OP0(X2, Y2, 18);
                DO_GROUP_OP0(X3, Y3, 20);
                DO_GROUP_OP0(X4, Y4, 22);
                DO_GROUP_OP2(X1, 16);
                DO_GROUP_OP2(X2, 18);
                DO_GROUP_OP2(X3, 20);
                DO_GROUP_OP2(X4, 22);
                DO_GROUP_OP0(X1, Y1, 24);
                DO_GROUP_OP0(X2, Y2, 26);
                DO_GROUP_OP0(X3, Y3, 28);
                DO_GROUP_OP0(X4, Y4, 30);
                DO_GROUP_OP2(X1, 24);
                DO_GROUP_OP2(X2, 26);
                DO_GROUP_OP2(X3, 28);
                DO_GROUP_OP2(X4, 30);
                DO_GROUP_OP0(X1, Y1, 32);
                DO_GROUP_OP0(X2, Y2, 34);
                DO_GROUP_OP0(X3, Y3, 36);
                DO_GROUP_OP0(X4, Y4, 38);
                DO_GROUP_OP2(X1, 32);
                DO_GROUP_OP2(X2, 34);
                DO_GROUP_OP2(X3, 36);
                DO_GROUP_OP2(X4, 38);
                DO_GROUP_OP0(X1, Y1, 40);
                DO_GROUP_OP0(X2, Y2, 42);
                DO_GROUP_OP0(X3, Y3, 44);
                DO_GROUP_OP0(X4, Y4, 46);
                DO_GROUP_OP2(X1, 40);
                DO_GROUP_OP2(X2, 42);
                DO_GROUP_OP2(X3, 44);
                DO_GROUP_OP2(X4, 46);
                DO_GROUP_OP0(X1, Y1, 48);
                DO_GROUP_OP0(X2, Y2, 50);
                DO_GROUP_OP0(X3, Y3, 52);
                DO_GROUP_OP0(X4, Y4, 54);
                DO_GROUP_OP2(X1, 48);
                DO_GROUP_OP2(X2, 50);
                DO_GROUP_OP2(X3, 52);
                DO_GROUP_OP2(X4, 54);
                DO_GROUP_OP0(X1, Y1, 56);
                DO_GROUP_OP0(X2, Y2, 58);
                DO_GROUP_OP2(X1, 56);
                DO_GROUP_OP2(X2, 58);

                for (long k = 60; k < vDim; k++) {
                  dest[row_offset + k] = theData[ti] * theData[k];
                }
                col_offset = vDim;
                ti++;
              }

              for (long c = 1; c < hDim; c++, ti++, col_offset += vDim) {
                float64x2_t A4 = vdupq_n_f64(theData[ti]);
                float64x2_t X1, X2, X3, X4, Y1, Y2, Y3, Y4;

                DO_GROUP_OP1(X1, Y1, 0);
                DO_GROUP_OP1(X2, Y2, 2);
                DO_GROUP_OP1(X3, Y3, 4);
                DO_GROUP_OP1(X4, Y4, 6);
                DO_GROUP_OP2(X1, 0);
                DO_GROUP_OP2(X2, 2);
                DO_GROUP_OP2(X3, 4);
                DO_GROUP_OP2(X4, 6);
                DO_GROUP_OP1(X1, Y1, 8);
                DO_GROUP_OP1(X2, Y2, 10);
                DO_GROUP_OP1(X3, Y3, 12);
                DO_GROUP_OP1(X4, Y4, 14);
                DO_GROUP_OP2(X1, 8);
                DO_GROUP_OP2(X2, 10);
                DO_GROUP_OP2(X3, 12);
                DO_GROUP_OP2(X4, 14);
                DO_GROUP_OP1(X1, Y1, 16);
                DO_GROUP_OP1(X2, Y2, 18);
                DO_GROUP_OP1(X3, Y3, 20);
                DO_GROUP_OP1(X4, Y4, 22);
                DO_GROUP_OP2(X1, 16);
                DO_GROUP_OP2(X2, 18);
                DO_GROUP_OP2(X3, 20);
                DO_GROUP_OP2(X4, 22);
                DO_GROUP_OP1(X1, Y1, 24);
                DO_GROUP_OP1(X2, Y2, 26);
                DO_GROUP_OP1(X3, Y3, 28);
                DO_GROUP_OP1(X4, Y4, 30);
                DO_GROUP_OP2(X1, 24);
                DO_GROUP_OP2(X2, 26);
                DO_GROUP_OP2(X3, 28);
                DO_GROUP_OP2(X4, 30);
                DO_GROUP_OP1(X1, Y1, 32);
                DO_GROUP_OP1(X2, Y2, 34);
                DO_GROUP_OP1(X3, Y3, 36);
                DO_GROUP_OP1(X4, Y4, 38);
                DO_GROUP_OP2(X1, 32);
                DO_GROUP_OP2(X2, 34);
                DO_GROUP_OP2(X3, 36);
                DO_GROUP_OP2(X4, 38);
                DO_GROUP_OP1(X1, Y1, 40);
                DO_GROUP_OP1(X2, Y2, 42);
                DO_GROUP_OP1(X3, Y3, 44);
                DO_GROUP_OP1(X4, Y4, 46);
                DO_GROUP_OP2(X1, 40);
                DO_GROUP_OP2(X2, 42);
                DO_GROUP_OP2(X3, 44);
                DO_GROUP_OP2(X4, 46);
                DO_GROUP_OP1(X1, Y1, 48);
                DO_GROUP_OP1(X2, Y2, 50);
                DO_GROUP_OP1(X3, Y3, 52);
                DO_GROUP_OP1(X4, Y4, 54);
                DO_GROUP_OP2(X1, 48);
                DO_GROUP_OP2(X2, 50);
                DO_GROUP_OP2(X3, 52);
                DO_GROUP_OP2(X4, 54);
                DO_GROUP_OP1(X1, Y1, 56);
                DO_GROUP_OP1(X2, Y2, 58);
                DO_GROUP_OP2(X1, 56);
                DO_GROUP_OP2(X2, 58);

                for (long k = 60; k < vDim; k++) {
                  dest[row_offset + k] += theData[ti] * theData[col_offset + k];
                }
              }
            }
          }
        } else {
          for (long r = 0; r < hDim; r++, row_offset += vDim) {
            long col_offset = 0L;
            { // row 1
              float64x2_t A4 = vdupq_n_f64(theData[ti]);
#pragma GCC unroll 4
#pragma clang loop vectorize(enable)
#pragma clang loop interleave(enable)
// #pragma clang loop unroll(enable)
#pragma GCC ivdep
#pragma ivdep
              for (long k = 0; k < loopBound; k += 4) {
                float64x2_t D4, B4, X1, X2;
                DO_GROUP_OP0(D4, B4, k);
                DO_GROUP_OP0(X1, X2, k + 2);
                DO_GROUP_OP2(D4, k);
                DO_GROUP_OP2(X1, k + 2);
              }

              for (long k = loopBound; k < vDim; k++) {
                dest[row_offset + k] = theData[ti] * theData[k];
              }
              col_offset = vDim;
              ti++;
            }

            for (long c = 1; c < hDim; c++, ti++, col_offset += vDim) {
              float64x2_t A4 = vdupq_n_f64(theData[ti]);
              for (long k = 0; k < loopBound; k += 4) {
                float64x2_t D4, B4, X1, X2;
                DO_GROUP_OP1(D4, B4, k);
                DO_GROUP_OP1(X1, X2, k + 2);
                DO_GROUP_OP2(D4, k);
                DO_GROUP_OP2(X1, k + 2);
              }

              for (long k = loopBound; k < vDim; k++) {
                dest[row_offset + k] += theData[ti] * theData[col_offset + k];
              }
            }
          }
        }
      } else {

#endif

#ifdef _SLKP_USE_AVX_INTRINSICS
#define DO_GROUP_OP0(X, Y, k)                                                  \
  Y = _mm256_loadu_pd(theData + col_offset + k);                               \
  X = _mm256_mul_pd(A4, Y);
#ifdef _SLKP_USE_FMA3_INTRINSICS
#define DO_GROUP_OP1(X, Y, k)                                                  \
  X = _mm256_loadu_pd(dest + row_offset + k);                                  \
  Y = _mm256_loadu_pd(theData + col_offset + k);                               \
  X = _mm256_fmadd_pd(A4, Y, X);
#define DO_GROUP_OP2(X, k) _mm256_storeu_pd(dest + row_offset + k, X);
#else
#define DO_GROUP_OP1(X, Y, k)                                                  \
  X = _mm256_loadu_pd(dest + row_offset + k);                                  \
  Y = _mm256_loadu_pd(theData + col_offset + k);                               \
  X = _mm256_add_pd(X, _mm256_mul_pd(A4, Y));
#define DO_GROUP_OP2(X, k) _mm256_storeu_pd(dest + row_offset + k, X);
#endif
        if (true) {
          hyFloat *_hprestrict_ dest = stash;

          long ti = 0L, row_offset = 0L;

          if (loopBound == 60UL) { // codons
            if (hDim == 61) {      // special case universal genetic code
              for (long r = 0; r < 61; r++, row_offset += 61) {
                long col_offset = 0L;

                { // handle first row separately to zero out the dest entries
                  __m256d A4 = _mm256_set1_pd(theData[ti]);
                  __m256d D4_1, D4_2, D4_3, D4_4;
                  __m256d B4_1, B4_2, B4_3, B4_4;

                  DO_GROUP_OP0(D4_1, B4_1, 0);
                  DO_GROUP_OP0(D4_2, B4_2, 4);
                  DO_GROUP_OP0(D4_3, B4_3, 8);
                  DO_GROUP_OP0(D4_4, B4_4, 12);
                  DO_GROUP_OP2(D4_1, 0);
                  DO_GROUP_OP2(D4_2, 4);
                  DO_GROUP_OP2(D4_3, 8);
                  DO_GROUP_OP2(D4_4, 12);

                  DO_GROUP_OP0(D4_1, B4_1, 16);
                  DO_GROUP_OP0(D4_2, B4_2, 20);
                  DO_GROUP_OP0(D4_3, B4_3, 24);
                  DO_GROUP_OP0(D4_4, B4_4, 28);
                  DO_GROUP_OP2(D4_1, 16);
                  DO_GROUP_OP2(D4_2, 20);
                  DO_GROUP_OP2(D4_3, 24);
                  DO_GROUP_OP2(D4_4, 28);

                  DO_GROUP_OP0(D4_1, B4_1, 32);
                  DO_GROUP_OP0(D4_2, B4_2, 36);
                  DO_GROUP_OP0(D4_3, B4_3, 40);
                  DO_GROUP_OP0(D4_4, B4_4, 44);
                  DO_GROUP_OP2(D4_1, 32);
                  DO_GROUP_OP2(D4_2, 36);
                  DO_GROUP_OP2(D4_3, 40);
                  DO_GROUP_OP2(D4_4, 44);

                  DO_GROUP_OP0(D4_1, B4_1, 48);
                  DO_GROUP_OP0(D4_2, B4_2, 52);
                  DO_GROUP_OP0(D4_3, B4_3, 56);
                  DO_GROUP_OP2(D4_1, 48);
                  DO_GROUP_OP2(D4_2, 52);
                  DO_GROUP_OP2(D4_3, 56);

                  dest[row_offset + 60] = theData[ti] * theData[60];
                  ti++;
                  col_offset = vDim;
                }
                for (long c = 1; c < 61; c++, ti++, col_offset += 61) {
                  __m256d A4 = _mm256_set1_pd(theData[ti]);

                  __m256d D4_1, D4_2, D4_3, D4_4;
                  __m256d B4_1, B4_2, B4_3, B4_4;

                  DO_GROUP_OP1(D4_1, B4_1, 0);
                  DO_GROUP_OP1(D4_2, B4_2, 4);
                  DO_GROUP_OP1(D4_3, B4_3, 8);
                  DO_GROUP_OP1(D4_4, B4_4, 12);
                  DO_GROUP_OP2(D4_1, 0);
                  DO_GROUP_OP2(D4_2, 4);
                  DO_GROUP_OP2(D4_3, 8);
                  DO_GROUP_OP2(D4_4, 12);

                  DO_GROUP_OP1(D4_1, B4_1, 16);
                  DO_GROUP_OP1(D4_2, B4_2, 20);
                  DO_GROUP_OP1(D4_3, B4_3, 24);
                  DO_GROUP_OP1(D4_4, B4_4, 28);
                  DO_GROUP_OP2(D4_1, 16);
                  DO_GROUP_OP2(D4_2, 20);
                  DO_GROUP_OP2(D4_3, 24);
                  DO_GROUP_OP2(D4_4, 28);

                  DO_GROUP_OP1(D4_1, B4_1, 32);
                  DO_GROUP_OP1(D4_2, B4_2, 36);
                  DO_GROUP_OP1(D4_3, B4_3, 40);
                  DO_GROUP_OP1(D4_4, B4_4, 44);
                  DO_GROUP_OP2(D4_1, 32);
                  DO_GROUP_OP2(D4_2, 36);
                  DO_GROUP_OP2(D4_3, 40);
                  DO_GROUP_OP2(D4_4, 44);

                  DO_GROUP_OP1(D4_1, B4_1, 48);
                  DO_GROUP_OP1(D4_2, B4_2, 52);
                  DO_GROUP_OP1(D4_3, B4_3, 56);
                  DO_GROUP_OP2(D4_1, 48);
                  DO_GROUP_OP2(D4_2, 52);
                  DO_GROUP_OP2(D4_3, 56);

                  dest[row_offset + 60] +=
                      theData[ti] * theData[col_offset + 60];
                }
              }

            } else {
              for (long r = 0; r < hDim; r++, row_offset += vDim) {
                long col_offset = 0L;

                { // handle first row separately to zero out the dest entries
                  __m256d A4 = _mm256_set1_pd(theData[ti]);
                  __m256d D4_1, D4_2, D4_3, D4_4;
                  __m256d B4_1, B4_2, B4_3, B4_4;

                  DO_GROUP_OP0(D4_1, B4_1, 0);
                  DO_GROUP_OP0(D4_2, B4_2, 4);
                  DO_GROUP_OP0(D4_3, B4_3, 8);
                  DO_GROUP_OP0(D4_4, B4_4, 12);
                  DO_GROUP_OP2(D4_1, 0);
                  DO_GROUP_OP2(D4_2, 4);
                  DO_GROUP_OP2(D4_3, 8);
                  DO_GROUP_OP2(D4_4, 12);

                  DO_GROUP_OP0(D4_1, B4_1, 16);
                  DO_GROUP_OP0(D4_2, B4_2, 20);
                  DO_GROUP_OP0(D4_3, B4_3, 24);
                  DO_GROUP_OP0(D4_4, B4_4, 28);
                  DO_GROUP_OP2(D4_1, 16);
                  DO_GROUP_OP2(D4_2, 20);
                  DO_GROUP_OP2(D4_3, 24);
                  DO_GROUP_OP2(D4_4, 28);

                  DO_GROUP_OP0(D4_1, B4_1, 32);
                  DO_GROUP_OP0(D4_2, B4_2, 36);
                  DO_GROUP_OP0(D4_3, B4_3, 40);
                  DO_GROUP_OP0(D4_4, B4_4, 44);
                  DO_GROUP_OP2(D4_1, 32);
                  DO_GROUP_OP2(D4_2, 36);
                  DO_GROUP_OP2(D4_3, 40);
                  DO_GROUP_OP2(D4_4, 44);

                  DO_GROUP_OP0(D4_1, B4_1, 48);
                  DO_GROUP_OP0(D4_2, B4_2, 52);
                  DO_GROUP_OP0(D4_3, B4_3, 56);
                  DO_GROUP_OP2(D4_1, 48);
                  DO_GROUP_OP2(D4_2, 52);
                  DO_GROUP_OP2(D4_3, 56);

                  for (long k = loopBound; k < vDim; k++) {
                    dest[row_offset + k] = theData[ti] * theData[k];
                  }
                  ti++;
                  col_offset = vDim;
                }
                for (long c = 1; c < hDim; c++, ti++, col_offset += vDim) {
                  __m256d A4 = _mm256_set1_pd(theData[ti]);

                  __m256d D4_1, D4_2, D4_3, D4_4;
                  __m256d B4_1, B4_2, B4_3, B4_4;

                  DO_GROUP_OP1(D4_1, B4_1, 0);
                  DO_GROUP_OP1(D4_2, B4_2, 4);
                  DO_GROUP_OP1(D4_3, B4_3, 8);
                  DO_GROUP_OP1(D4_4, B4_4, 12);
                  DO_GROUP_OP2(D4_1, 0);
                  DO_GROUP_OP2(D4_2, 4);
                  DO_GROUP_OP2(D4_3, 8);
                  DO_GROUP_OP2(D4_4, 12);

                  DO_GROUP_OP1(D4_1, B4_1, 16);
                  DO_GROUP_OP1(D4_2, B4_2, 20);
                  DO_GROUP_OP1(D4_3, B4_3, 24);
                  DO_GROUP_OP1(D4_4, B4_4, 28);
                  DO_GROUP_OP2(D4_1, 16);
                  DO_GROUP_OP2(D4_2, 20);
                  DO_GROUP_OP2(D4_3, 24);
                  DO_GROUP_OP2(D4_4, 28);

                  DO_GROUP_OP1(D4_1, B4_1, 32);
                  DO_GROUP_OP1(D4_2, B4_2, 36);
                  DO_GROUP_OP1(D4_3, B4_3, 40);
                  DO_GROUP_OP1(D4_4, B4_4, 44);
                  DO_GROUP_OP2(D4_1, 32);
                  DO_GROUP_OP2(D4_2, 36);
                  DO_GROUP_OP2(D4_3, 40);
                  DO_GROUP_OP2(D4_4, 44);

                  DO_GROUP_OP1(D4_1, B4_1, 48);
                  DO_GROUP_OP1(D4_2, B4_2, 52);
                  DO_GROUP_OP1(D4_3, B4_3, 56);
                  DO_GROUP_OP2(D4_1, 48);
                  DO_GROUP_OP2(D4_2, 52);
                  DO_GROUP_OP2(D4_3, 56);

                  for (long k = loopBound; k < vDim; k++) {
                    dest[row_offset + k] +=
                        theData[ti] * theData[col_offset + k];
                  }
                }
              }
            }
          } else { // something else
            for (long r = 0; r < hDim; r++, row_offset += vDim) {
              long col_offset = 0L;

              {
                __m256d A4 = _mm256_set1_pd(theData[ti]);
                for (long k = 0; k < loopBound; k += 4) {
                  __m256d D4, B4;
                  DO_GROUP_OP0(D4, B4, k);
                  DO_GROUP_OP2(D4, k);
                }

                for (long k = loopBound; k < vDim; k++) {
                  dest[row_offset + k] = theData[ti] * theData[k];
                }
                col_offset = vDim;
                ti++;
              }

              for (long c = 1; c < hDim; c++, ti++, col_offset += vDim) {
                __m256d A4 = _mm256_set1_pd(theData[ti]);
                for (long k = 0; k < loopBound; k += 4) {
                  __m256d D4, B4;
                  DO_GROUP_OP1(D4, B4, k);
                  DO_GROUP_OP2(D4, k);
                }

                for (long k = loopBound; k < vDim; k++) {
                  dest[row_offset + k] += theData[ti] * theData[col_offset + k];
                }
              }
            }
          }
        } else {

#endif

          hyFloat *_hprestrict_ column = stash + lDim;
          hyFloat const *source = theData;

          for (long j = 0; j < vDim; j++) {
            for (long c = 0; c < vDim; c++) {
              column[c] = source[j + c * vDim];
            }

#ifdef _SLKP_USE_AVX_INTRINSICS
            if (vDim == 61UL) {
              for (unsigned long i = 0; i < lDim; i += 61) {
                hyFloat *_hprestrict_ row = theData + i;

                __m256d sum256;

#ifdef _SLKP_USE_FMA3_INTRINSICS

                sum256 = _mm256_fmadd_pd(
                    _mm256_loadu_pd(row), _mm256_loadu_pd(column),
                    _mm256_fmadd_pd(
                        _mm256_loadu_pd(row + 4), _mm256_loadu_pd(column + 4),
                        _mm256_mul_pd(_mm256_loadu_pd(row + 8),
                                      _mm256_loadu_pd(column + 8))));

                for (unsigned long k = 12UL; k < 60UL; k += 12UL) {
                  sum256 = _mm256_fmadd_pd(
                      _mm256_loadu_pd(row + k), _mm256_loadu_pd(column + k),
                      _mm256_fmadd_pd(
                          _mm256_loadu_pd(row + k + 4),
                          _mm256_loadu_pd(column + k + 4),
                          _mm256_fmadd_pd(_mm256_loadu_pd(row + k + 8),
                                          _mm256_loadu_pd(column + k + 8),
                                          sum256)));
                }
#else

                sum256 = _mm256_setzero_pd();
                for (unsigned long k = 0UL; k < 60UL; k += 12UL) {
                  __m256d term0 = _mm256_mul_pd(_mm256_loadu_pd(row + k),
                                                _mm256_loadu_pd(column + k));
                  __m256d term1 =
                      _mm256_mul_pd(_mm256_loadu_pd(row + k + 4),
                                    _mm256_loadu_pd(column + k + 4));
                  __m256d term2 =
                      _mm256_mul_pd(_mm256_loadu_pd(row + k + 8),
                                    _mm256_loadu_pd(column + k + 8));

                  __m256d sum01 = _mm256_add_pd(term0, term1);
                  __m256d plus2 = _mm256_add_pd(term2, sum256);

                  sum256 = _mm256_add_pd(sum01, plus2);
                }
#endif

                stash[i + j] = _avx_sum_4(sum256) + row[60] * column[60];
              }

            } else {
              for (unsigned long i = 0; i < lDim; i += vDim) {
                hyFloat *row = theData + i;

                __m256d sum256 = _mm256_setzero_pd();

                long k;

                for (k = 0; k < loopBound; k += 4) {
#ifdef _SLKP_USE_FMA3_INTRINSICS
                  sum256 = _mm256_fmadd_pd(_mm256_loadu_pd(row + k),
                                           _mm256_loadu_pd(column + k), sum256);
#else
                  sum256 =
                      _mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(row + k),
                                                  _mm256_loadu_pd(column + k)),
                                    sum256);
#endif
                }

                hyFloat result = _avx_sum_4(sum256);

                for (; k < vDim; k++) {
                  result += row[k] * column[k];
                }

                stash[i + j] = result;
              }
            }
#ifdef _SLKP_USE_AVX_INTRINSICS
          }
#endif
#else
        for (long i = 0; i < lDim; i += vDim) {
          hyFloat *row = theData + i, buffer[4] = {0., 0., 0., 0.};

          unsigned long k;

          for (k = 0UL; k < loopBound; k += 4UL) {
            buffer[0] += row[k] * column[k];
            buffer[1] += row[k + 1] * column[k + 1];
            buffer[2] += row[k + 2] * column[k + 2];
            buffer[3] += row[k + 3] * column[k + 3];
          }

          for (; k < vDim; k++) {
            buffer[0] += row[k] * column[k];
          }

          if (i + j >= lDim || i + j < 0) {
            printf("Well, shit\n");
          }
          stash[i + j] = (buffer[0] + buffer[1]) + (buffer[2] + buffer[3]);
        }
#endif
        }
      }
#ifdef _SLKP_USE_ARM_NEON
    }
#endif

    /*long lDimmod4 = (lDim >> 2) << 2;
    hyFloat diffs[4] = {0.0,0.0,0.0,0.0};

    for (long s = 0; s < lDimmod4; s+=4) {
        hyFloat d1 = fabs (theData[s  ] - stash[s  ]);
        hyFloat d2 = fabs (theData[s+1] - stash[s+1]);
        hyFloat d3 = fabs (theData[s+2] - stash[s+2]);
        hyFloat d4 = fabs (theData[s+3] - stash[s+3]);
        if (d1 > diffs[0]) diffs[0] = d1;
        if (d2 > diffs[1]) diffs[1] = d2;
        if (d3 > diffs[2]) diffs[2] = d3;
        if (d4 > diffs[3]) diffs[3] = d4;
    }

    for (long s = lDimmod4; s < lDim; s++) {
        hyFloat d1 = fabs (theData[s] - stash[s]);
        if (d1 > diffs[0]) diffs[0] = d1;
    }

    diff = MAX (MAX (diffs[0], diffs[1]), MAX (diffs[2], diffs[3]));

    memcpy (theData, stash, lDim * sizeof (hyFloat));*/

    /*for (long s = 0; s < lDim; s++) {
        StoreIfGreater(diff, fabs (theData[s] - stash[s]));
        theData[s] = stash[s];
    }*/
  }
  return diff;
}

double *src = new double[N * N];
double *rec = new double[2 * N * N];
double *check = new double[N * N];

bool already_init = false;

void init_mx(double *src) {
  if (!already_init) {
    for (int i = 0; i < N; i++) {
      double s = 0.;
      for (int j = 0; j < N; j++) {
        if (i != j) {
          double e = rand() / (double)RAND_MAX;
          src[i * N + j] = e;
          s += e;
        }
      }
      src[i * N + i] = -s;
    }
    already_init = true;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1., src, N,
                src, N, 0.0, check, N);
  }
}

void naive_multiply(double *C, double const *A, double const *B, int D) {
  const double *row = A;
  for (int r = 0; r < D; r++, row += D) {
    for (int c = 0; c < D; c++, C++) {
      double acc = 0.;
      const double *col = B + c;
      for (int k = 0; k < D; k++, col += D) {
        acc += row[k] * *col;
      }
      *C = acc;
    }
  }
}

void very_naive_multiply(double *C, double const *A, double const *B, int D) {
  for (int r = 0; r < D; r++) {
    for (int c = 0; c < D; c++) {
      double acc = 0.;
      for (int k = 0; k < D; k++) {
        acc += A[r * D + k] * B[c + D * k];
      }
      C[r * D + c] = acc;
    }
  }
}

void direct4x4_multiply(double *C, double const *A, double const *B) {
  for (int r = 0; r < 4; r++) {
    for (int c = 0; c < 4; c++) {
      double s1 = A[4 * r] * B[c] + A[4 * r + 1] * B[c + 4];
      double s2 = A[4 * r + 2] * B[c + 8] + A[4 * r + 3] * B[c + 12];
      C[r * 4 + c] = s1 + s2;
    }
  }
}

double check_err(const double *a, const double *b) {
  double max_d = 0.;
  for (int i = 0; i < N * N; i++) {
    double d = fabs(a[i] - b[i]);
    if (d > max_d)
      max_d = d;
  }
  return max_d;
}

// BENCHMARK(BM_LoadBy4)
//     ->RangeMultiplier(2)->Range(1<<4, 1<<12)->Complexity();

static void naive_SQR(benchmark::State &state) {
  init_mx(src);
  for (auto _ : state) {
    naive_multiply(rec, src, src, N);
  }
  printf("NAIVE err %g\n", check_err(rec, check));
}

static void vnaive_SQR(benchmark::State &state) {
  init_mx(src);
  for (auto _ : state) {
    very_naive_multiply(rec, src, src, N);
  }
  printf("V. NAIVE err %g\n", check_err(rec, check));
}

static void direct4x4_SQR(benchmark::State &state) {
  init_mx(src);
  for (auto _ : state) {
    direct4x4_multiply(rec, src, src);
  }
  printf("Direct 4x4 err %g\n", check_err(rec, check));
}

static void blas_SQR(benchmark::State &state) {
  init_mx(src);
  for (auto _ : state) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1., src, N,
                src, N, 0.0, rec, N);
  }
  printf("BLAS err %g\n", check_err(rec, check));
}

BENCHMARK(vnaive_SQR);
BENCHMARK(naive_SQR);
BENCHMARK(direct4x4_SQR);
BENCHMARK(blas_SQR);

BENCHMARK_MAIN();
