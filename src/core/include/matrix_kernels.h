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

#ifndef __MATRIX_KERNELS__
#define __MATRIX_KERNELS__

#include "matrix.h"

template <int N>
void _matrix_kernel_multiply_by_compressed_sparse_T(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr);

#if defined _SLKP_USE_ARM_SVE
#include <arm_sve.h>

template <int N>
void _matrix_kernel_multiply_by_compressed_sparse_T(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  uint64_t vl = svcntd();

  for (long i = 0; i < N; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {
      long j = 0;
      // Unrolled loop
      while (j + 4 * vl <= N) {
        svbool_t pgAll = svptrue_b64();
        svfloat64_t acc0 = svld1(pgAll, res + j);
        svfloat64_t acc1 = svld1(pgAll, res + j + vl);
        svfloat64_t acc2 = svld1(pgAll, res + j + 2 * vl);
        svfloat64_t acc3 = svld1(pgAll, res + j + 3 * vl);

        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;
          svfloat64_t valVec = svdup_f64(theDataPtr[cxi]);

          svfloat64_t b0 = svld1(pgAll, secArg + j);
          svfloat64_t b1 = svld1(pgAll, secArg + j + vl);
          svfloat64_t b2 = svld1(pgAll, secArg + j + 2 * vl);
          svfloat64_t b3 = svld1(pgAll, secArg + j + 3 * vl);

          acc0 = svmla_f64_x(pgAll, acc0, valVec, b0);
          acc1 = svmla_f64_x(pgAll, acc1, valVec, b1);
          acc2 = svmla_f64_x(pgAll, acc2, valVec, b2);
          acc3 = svmla_f64_x(pgAll, acc3, valVec, b3);
        }

        svst1(pgAll, res + j, acc0);
        svst1(pgAll, res + j + vl, acc1);
        svst1(pgAll, res + j + 2 * vl, acc2);
        svst1(pgAll, res + j + 3 * vl, acc3);

        j += 4 * vl;
      }

      // Tail loop
      while (j < N) {
        svbool_t pg = svwhilelt_b64(j, N);
        svfloat64_t acc = svld1(pg, res + j);

        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;
          svfloat64_t valVec = svdup_f64(theDataPtr[cxi]);
          svfloat64_t b = svld1(pg, secArg + j);
          acc = svmla_f64_m(pg, acc, valVec, b);
        }
        svst1(pg, res + j, acc);
        j += vl;
      }
    }
    res += N;
    currentXIndex = up;
  }
}

#elif defined _SLKP_USE_ARM_NEON

template <>
inline void _matrix_kernel_multiply_by_compressed_sparse_T<61>(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // N = 61
  // 61 doubles = 30.5 vectors of 2 doubles
  // We use 30 vector registers + 1 scalar
  // Regs: v0..v29 for vector data, v30 for scalar part (or scalar reg)
  // We need 1 register for 'val' splat, and 1 temp for loading 'secArg'.
  // ARM64 has 32 vector registers v0-v31.
  // So:
  // v0..v29 : Accumulators (60 doubles)
  // scalar_acc : Last double
  // v30 : temp load
  // v31 : val splat

  for (long i = 0; i < 61; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {
      // 1. Load all 61 elements into registers
      float64x2_t acc0 = vld1q_f64(res + 0);
      float64x2_t acc1 = vld1q_f64(res + 2);
      float64x2_t acc2 = vld1q_f64(res + 4);
      float64x2_t acc3 = vld1q_f64(res + 6);
      float64x2_t acc4 = vld1q_f64(res + 8);
      float64x2_t acc5 = vld1q_f64(res + 10);
      float64x2_t acc6 = vld1q_f64(res + 12);
      float64x2_t acc7 = vld1q_f64(res + 14);
      float64x2_t acc8 = vld1q_f64(res + 16);
      float64x2_t acc9 = vld1q_f64(res + 18);
      float64x2_t acc10 = vld1q_f64(res + 20);
      float64x2_t acc11 = vld1q_f64(res + 22);
      float64x2_t acc12 = vld1q_f64(res + 24);
      float64x2_t acc13 = vld1q_f64(res + 26);
      float64x2_t acc14 = vld1q_f64(res + 28);
      float64x2_t acc15 = vld1q_f64(res + 30);
      float64x2_t acc16 = vld1q_f64(res + 32);
      float64x2_t acc17 = vld1q_f64(res + 34);
      float64x2_t acc18 = vld1q_f64(res + 36);
      float64x2_t acc19 = vld1q_f64(res + 38);
      float64x2_t acc20 = vld1q_f64(res + 40);
      float64x2_t acc21 = vld1q_f64(res + 42);
      float64x2_t acc22 = vld1q_f64(res + 44);
      float64x2_t acc23 = vld1q_f64(res + 46);
      float64x2_t acc24 = vld1q_f64(res + 48);
      float64x2_t acc25 = vld1q_f64(res + 50);
      float64x2_t acc26 = vld1q_f64(res + 52);
      float64x2_t acc27 = vld1q_f64(res + 54);
      float64x2_t acc28 = vld1q_f64(res + 56);
      float64x2_t acc29 = vld1q_f64(res + 58);
      hyFloat acc30_scalar = res[60];

      // 2. Accumulate
      for (long cxi = currentXIndex; cxi < up; cxi++) {
        long currentXColumn = compIdx[cxi + 61];
        const hyFloat *secArg = secondArgBase + currentXColumn * 61;

        hyFloat val = theDataPtr[cxi];
        float64x2_t val_vec = vdupq_n_f64(val);

        float64x2_t C;
        C = vld1q_f64(secArg + 0);
        acc0 = vfmaq_f64(acc0, val_vec, C);
        C = vld1q_f64(secArg + 2);
        acc1 = vfmaq_f64(acc1, val_vec, C);
        C = vld1q_f64(secArg + 4);
        acc2 = vfmaq_f64(acc2, val_vec, C);
        C = vld1q_f64(secArg + 6);
        acc3 = vfmaq_f64(acc3, val_vec, C);
        C = vld1q_f64(secArg + 8);
        acc4 = vfmaq_f64(acc4, val_vec, C);
        C = vld1q_f64(secArg + 10);
        acc5 = vfmaq_f64(acc5, val_vec, C);
        C = vld1q_f64(secArg + 12);
        acc6 = vfmaq_f64(acc6, val_vec, C);
        C = vld1q_f64(secArg + 14);
        acc7 = vfmaq_f64(acc7, val_vec, C);
        C = vld1q_f64(secArg + 16);
        acc8 = vfmaq_f64(acc8, val_vec, C);
        C = vld1q_f64(secArg + 18);
        acc9 = vfmaq_f64(acc9, val_vec, C);
        C = vld1q_f64(secArg + 20);
        acc10 = vfmaq_f64(acc10, val_vec, C);
        C = vld1q_f64(secArg + 22);
        acc11 = vfmaq_f64(acc11, val_vec, C);
        C = vld1q_f64(secArg + 24);
        acc12 = vfmaq_f64(acc12, val_vec, C);
        C = vld1q_f64(secArg + 26);
        acc13 = vfmaq_f64(acc13, val_vec, C);
        C = vld1q_f64(secArg + 28);
        acc14 = vfmaq_f64(acc14, val_vec, C);
        C = vld1q_f64(secArg + 30);
        acc15 = vfmaq_f64(acc15, val_vec, C);
        C = vld1q_f64(secArg + 32);
        acc16 = vfmaq_f64(acc16, val_vec, C);
        C = vld1q_f64(secArg + 34);
        acc17 = vfmaq_f64(acc17, val_vec, C);
        C = vld1q_f64(secArg + 36);
        acc18 = vfmaq_f64(acc18, val_vec, C);
        C = vld1q_f64(secArg + 38);
        acc19 = vfmaq_f64(acc19, val_vec, C);
        C = vld1q_f64(secArg + 40);
        acc20 = vfmaq_f64(acc20, val_vec, C);
        C = vld1q_f64(secArg + 42);
        acc21 = vfmaq_f64(acc21, val_vec, C);
        C = vld1q_f64(secArg + 44);
        acc22 = vfmaq_f64(acc22, val_vec, C);
        C = vld1q_f64(secArg + 46);
        acc23 = vfmaq_f64(acc23, val_vec, C);
        C = vld1q_f64(secArg + 48);
        acc24 = vfmaq_f64(acc24, val_vec, C);
        C = vld1q_f64(secArg + 50);
        acc25 = vfmaq_f64(acc25, val_vec, C);
        C = vld1q_f64(secArg + 52);
        acc26 = vfmaq_f64(acc26, val_vec, C);
        C = vld1q_f64(secArg + 54);
        acc27 = vfmaq_f64(acc27, val_vec, C);
        C = vld1q_f64(secArg + 56);
        acc28 = vfmaq_f64(acc28, val_vec, C);
        C = vld1q_f64(secArg + 58);
        acc29 = vfmaq_f64(acc29, val_vec, C);
        acc30_scalar += val * secArg[60];
      }

      // 3. Store all 61 elements
      vst1q_f64(res + 0, acc0);
      vst1q_f64(res + 2, acc1);
      vst1q_f64(res + 4, acc2);
      vst1q_f64(res + 6, acc3);
      vst1q_f64(res + 8, acc4);
      vst1q_f64(res + 10, acc5);
      vst1q_f64(res + 12, acc6);
      vst1q_f64(res + 14, acc7);
      vst1q_f64(res + 16, acc8);
      vst1q_f64(res + 18, acc9);
      vst1q_f64(res + 20, acc10);
      vst1q_f64(res + 22, acc11);
      vst1q_f64(res + 24, acc12);
      vst1q_f64(res + 26, acc13);
      vst1q_f64(res + 28, acc14);
      vst1q_f64(res + 30, acc15);
      vst1q_f64(res + 32, acc16);
      vst1q_f64(res + 34, acc17);
      vst1q_f64(res + 36, acc18);
      vst1q_f64(res + 38, acc19);
      vst1q_f64(res + 40, acc20);
      vst1q_f64(res + 42, acc21);
      vst1q_f64(res + 44, acc22);
      vst1q_f64(res + 46, acc23);
      vst1q_f64(res + 48, acc24);
      vst1q_f64(res + 50, acc25);
      vst1q_f64(res + 52, acc26);
      vst1q_f64(res + 54, acc27);
      vst1q_f64(res + 56, acc28);
      vst1q_f64(res + 58, acc29);
      res[60] = acc30_scalar;
    }

    res += 61;
    currentXIndex = up;
  }
}

template <int N>
void _matrix_kernel_multiply_by_compressed_sparse_T(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // Compute layout constants at compile time
  constexpr int STRIP_SIZE = 32; // 16 vectors
  constexpr int NUM_FULL_STRIPS = N / STRIP_SIZE;
  constexpr int TAIL_SIZE = N % STRIP_SIZE;

  for (long i = 0; i < N; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {

// -------------------------------------------------
// PART A: Compile-time loop for Full Strips (0..31, 32..63, etc)
// -------------------------------------------------
#pragma unroll
      for (int s = 0; s < NUM_FULL_STRIPS; s++) {
        // We know exactly where we are: s * 32
        const int base = s * STRIP_SIZE;

        float64x2_t acc[16];

        // 1. Load Accumulators
        for (int k = 0; k < 16; k++) {
          acc[k] = vld1q_f64(res + base + (k * 2));
        }

        // 2. Accumulate
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N]; // constant N
          const hyFloat *secArg = secondArgBase + currentXColumn * N;

          float64x2_t val_vec = vdupq_n_f64(theDataPtr[cxi]);

          // Fully unrolled FMA sequence (4 blocks of 4 vectors)
          for (int k = 0; k < 4; k++) {
            float64x2x4_t C = vld1q_f64_x4(secArg + base + (k * 8));
            acc[k * 4 + 0] = vfmaq_f64(acc[k * 4 + 0], val_vec, C.val[0]);
            acc[k * 4 + 1] = vfmaq_f64(acc[k * 4 + 1], val_vec, C.val[1]);
            acc[k * 4 + 2] = vfmaq_f64(acc[k * 4 + 2], val_vec, C.val[2]);
            acc[k * 4 + 3] = vfmaq_f64(acc[k * 4 + 3], val_vec, C.val[3]);
          }
        }

        // 3. Store Accumulators
        for (int k = 0; k < 16; k++) {
          vst1q_f64(res + base + (k * 2), acc[k]);
        }
      }

      // -------------------------------------------------
      // PART B: The Tail (Specialized for 60, 61, 62, 63...)
      // -------------------------------------------------
      if constexpr (TAIL_SIZE > 0) {
        constexpr int base = NUM_FULL_STRIPS * STRIP_SIZE;
        constexpr int vecs_needed = TAIL_SIZE / 2;
        constexpr bool has_scalar = (TAIL_SIZE % 2 != 0);

        // These arrays are sized by constant expressions;
        // compiler allocates strictly registers, no malloc.
        float64x2_t acc[16];
        hyFloat scalar_acc = 0.0;

        // 1. Load Tail
        for (int k = 0; k < vecs_needed; k++) {
          acc[k] = vld1q_f64(res + base + (k * 2));
        }
        if constexpr (has_scalar) {
          scalar_acc = res[base + TAIL_SIZE - 1];
        }

        // 2. Accumulate Tail
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;

          hyFloat val = theDataPtr[cxi];
          float64x2_t val_vec = vdupq_n_f64(val);

          // Unrolling logic for the tail
          int vec_offset = 0;
          int d_offset = 0;

          // Chunk 8 doubles (4 vectors) - heavily unrolled
          // For N=61, Tail=29. This runs 3 times (24 elements).
          constexpr int chunks_of_8 = TAIL_SIZE / 8;
          for (int k = 0; k < chunks_of_8; ++k) {
            float64x2x4_t C = vld1q_f64_x4(secArg + base + d_offset);
            acc[vec_offset + 0] =
                vfmaq_f64(acc[vec_offset + 0], val_vec, C.val[0]);
            acc[vec_offset + 1] =
                vfmaq_f64(acc[vec_offset + 1], val_vec, C.val[1]);
            acc[vec_offset + 2] =
                vfmaq_f64(acc[vec_offset + 2], val_vec, C.val[2]);
            acc[vec_offset + 3] =
                vfmaq_f64(acc[vec_offset + 3], val_vec, C.val[3]);
            vec_offset += 4;
            d_offset += 8;
          }

          // Chunk 2 doubles (1 vector)
          // For N=61, remaining=5. This runs 2 times (4 elements).
          constexpr int remaining_vectors = (TAIL_SIZE % 8) / 2;
          for (int k = 0; k < remaining_vectors; ++k) {
            float64x2_t C = vld1q_f64(secArg + base + d_offset);
            acc[vec_offset] = vfmaq_f64(acc[vec_offset], val_vec, C);
            vec_offset++;
            d_offset += 2;
          }

          // Scalar
          if constexpr (has_scalar) {
            scalar_acc += val * secArg[base + TAIL_SIZE - 1];
          }
        }

        // 3. Store Tail
        for (int k = 0; k < vecs_needed; k++) {
          vst1q_f64(res + base + (k * 2), acc[k]);
        }
        if constexpr (has_scalar) {
          res[base + TAIL_SIZE - 1] = scalar_acc;
        }
      }
    }

    res += N;
    currentXIndex = up;
  }
}

#elif defined _SLKP_USE_SSE_INTRINSICS

template <>
inline void _matrix_kernel_multiply_by_compressed_sparse_T<61>(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // N=61, SSE (128-bit, 2 doubles).
  // 61 doubles = 30.5 vectors.
  // 16 XMM registers available. We can't hold all 31 vectors.
  // Split into 2 passes:
  // Pass 1: 0..31 (16 vectors = 32 doubles) -> Uses all 16 regs. Spills
  // val/temp. Pass 2: 32..60 (14 vectors + 1 scalar = 29 doubles).

  for (long i = 0; i < 61; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {
      // --- Pass 1: Elements 0..31 (16 vectors) ---
      __m128d acc0 = _mm_loadu_pd(res + 0);
      __m128d acc1 = _mm_loadu_pd(res + 2);
      __m128d acc2 = _mm_loadu_pd(res + 4);
      __m128d acc3 = _mm_loadu_pd(res + 6);
      __m128d acc4 = _mm_loadu_pd(res + 8);
      __m128d acc5 = _mm_loadu_pd(res + 10);
      __m128d acc6 = _mm_loadu_pd(res + 12);
      __m128d acc7 = _mm_loadu_pd(res + 14);
      __m128d acc8 = _mm_loadu_pd(res + 16);
      __m128d acc9 = _mm_loadu_pd(res + 18);
      __m128d acc10 = _mm_loadu_pd(res + 20);
      __m128d acc11 = _mm_loadu_pd(res + 22);
      __m128d acc12 = _mm_loadu_pd(res + 24);
      __m128d acc13 = _mm_loadu_pd(res + 26);
      __m128d acc14 = _mm_loadu_pd(res + 28);
      __m128d acc15 = _mm_loadu_pd(res + 30);

      for (long cxi = currentXIndex; cxi < up; cxi++) {
        long currentXColumn = compIdx[cxi + 61];
        const hyFloat *secArg = secondArgBase + currentXColumn * 61;
        __m128d val_vec = _mm_set1_pd(theDataPtr[cxi]);
        __m128d C;

        C = _mm_loadu_pd(secArg + 0);
        acc0 = _mm_add_pd(acc0, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 2);
        acc1 = _mm_add_pd(acc1, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 4);
        acc2 = _mm_add_pd(acc2, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 6);
        acc3 = _mm_add_pd(acc3, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 8);
        acc4 = _mm_add_pd(acc4, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 10);
        acc5 = _mm_add_pd(acc5, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 12);
        acc6 = _mm_add_pd(acc6, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 14);
        acc7 = _mm_add_pd(acc7, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 16);
        acc8 = _mm_add_pd(acc8, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 18);
        acc9 = _mm_add_pd(acc9, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 20);
        acc10 = _mm_add_pd(acc10, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 22);
        acc11 = _mm_add_pd(acc11, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 24);
        acc12 = _mm_add_pd(acc12, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 26);
        acc13 = _mm_add_pd(acc13, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 28);
        acc14 = _mm_add_pd(acc14, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 30);
        acc15 = _mm_add_pd(acc15, _mm_mul_pd(val_vec, C));
      }
      _mm_storeu_pd(res + 0, acc0);
      _mm_storeu_pd(res + 2, acc1);
      _mm_storeu_pd(res + 4, acc2);
      _mm_storeu_pd(res + 6, acc3);
      _mm_storeu_pd(res + 8, acc4);
      _mm_storeu_pd(res + 10, acc5);
      _mm_storeu_pd(res + 12, acc6);
      _mm_storeu_pd(res + 14, acc7);
      _mm_storeu_pd(res + 16, acc8);
      _mm_storeu_pd(res + 18, acc9);
      _mm_storeu_pd(res + 20, acc10);
      _mm_storeu_pd(res + 22, acc11);
      _mm_storeu_pd(res + 24, acc12);
      _mm_storeu_pd(res + 26, acc13);
      _mm_storeu_pd(res + 28, acc14);
      _mm_storeu_pd(res + 30, acc15);

      // --- Pass 2: Elements 32..60 (14 vectors + 1 scalar) ---
      // Re-use registers acc0..acc14
      acc0 = _mm_loadu_pd(res + 32);
      acc1 = _mm_loadu_pd(res + 34);
      acc2 = _mm_loadu_pd(res + 36);
      acc3 = _mm_loadu_pd(res + 38);
      acc4 = _mm_loadu_pd(res + 40);
      acc5 = _mm_loadu_pd(res + 42);
      acc6 = _mm_loadu_pd(res + 44);
      acc7 = _mm_loadu_pd(res + 46);
      acc8 = _mm_loadu_pd(res + 48);
      acc9 = _mm_loadu_pd(res + 50);
      acc10 = _mm_loadu_pd(res + 52);
      acc11 = _mm_loadu_pd(res + 54);
      acc12 = _mm_loadu_pd(res + 56);
      acc13 = _mm_loadu_pd(res + 58);
      hyFloat scalar_acc = res[60];

      for (long cxi = currentXIndex; cxi < up; cxi++) {
        long currentXColumn = compIdx[cxi + 61];
        const hyFloat *secArg = secondArgBase + currentXColumn * 61;
        hyFloat val = theDataPtr[cxi];
        __m128d val_vec = _mm_set1_pd(val);
        __m128d C;

        C = _mm_loadu_pd(secArg + 32);
        acc0 = _mm_add_pd(acc0, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 34);
        acc1 = _mm_add_pd(acc1, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 36);
        acc2 = _mm_add_pd(acc2, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 38);
        acc3 = _mm_add_pd(acc3, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 40);
        acc4 = _mm_add_pd(acc4, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 42);
        acc5 = _mm_add_pd(acc5, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 44);
        acc6 = _mm_add_pd(acc6, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 46);
        acc7 = _mm_add_pd(acc7, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 48);
        acc8 = _mm_add_pd(acc8, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 50);
        acc9 = _mm_add_pd(acc9, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 52);
        acc10 = _mm_add_pd(acc10, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 54);
        acc11 = _mm_add_pd(acc11, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 56);
        acc12 = _mm_add_pd(acc12, _mm_mul_pd(val_vec, C));
        C = _mm_loadu_pd(secArg + 58);
        acc13 = _mm_add_pd(acc13, _mm_mul_pd(val_vec, C));
        scalar_acc += val * secArg[60];
      }
      _mm_storeu_pd(res + 32, acc0);
      _mm_storeu_pd(res + 34, acc1);
      _mm_storeu_pd(res + 36, acc2);
      _mm_storeu_pd(res + 38, acc3);
      _mm_storeu_pd(res + 40, acc4);
      _mm_storeu_pd(res + 42, acc5);
      _mm_storeu_pd(res + 44, acc6);
      _mm_storeu_pd(res + 46, acc7);
      _mm_storeu_pd(res + 48, acc8);
      _mm_storeu_pd(res + 50, acc9);
      _mm_storeu_pd(res + 52, acc10);
      _mm_storeu_pd(res + 54, acc11);
      _mm_storeu_pd(res + 56, acc12);
      _mm_storeu_pd(res + 58, acc13);
      res[60] = scalar_acc;
    }

    res += 61;
    currentXIndex = up;
  }
}

template <int N>
void _matrix_kernel_multiply_by_compressed_sparse_T(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // Compute layout constants at compile time
  constexpr int STRIP_SIZE = 32; // 16 vectors
  constexpr int NUM_FULL_STRIPS = N / STRIP_SIZE;
  constexpr int TAIL_SIZE = N % STRIP_SIZE;

  for (long i = 0; i < N; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {

// -------------------------------------------------
// PART A: Compile-time loop for Full Strips (0..31, 32..63, etc)
// -------------------------------------------------
#pragma unroll
      for (int s = 0; s < NUM_FULL_STRIPS; s++) {
        // We know exactly where we are: s * 32
        const int base = s * STRIP_SIZE;

        __m128d acc[16];

        // 1. Load Accumulators
        for (int k = 0; k < 16; k++) {
          acc[k] = _mm_loadu_pd(res + base + (k * 2));
        }

        // 2. Accumulate
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N]; // constant N
          const hyFloat *secArg = secondArgBase + currentXColumn * N;

          __m128d val_vec = _mm_set1_pd(theDataPtr[cxi]);

          // Fully unrolled FMA sequence (4 blocks of 4 vectors)
          for (int k = 0; k < 4; k++) {
            __m128d C0 = _mm_loadu_pd(secArg + base + (k * 8) + 0);
            __m128d C1 = _mm_loadu_pd(secArg + base + (k * 8) + 2);
            __m128d C2 = _mm_loadu_pd(secArg + base + (k * 8) + 4);
            __m128d C3 = _mm_loadu_pd(secArg + base + (k * 8) + 6);

            acc[k * 4 + 0] =
                _mm_add_pd(acc[k * 4 + 0], _mm_mul_pd(val_vec, C0));
            acc[k * 4 + 1] =
                _mm_add_pd(acc[k * 4 + 1], _mm_mul_pd(val_vec, C1));
            acc[k * 4 + 2] =
                _mm_add_pd(acc[k * 4 + 2], _mm_mul_pd(val_vec, C2));
            acc[k * 4 + 3] =
                _mm_add_pd(acc[k * 4 + 3], _mm_mul_pd(val_vec, C3));
          }
        }

        // 3. Store Accumulators
        for (int k = 0; k < 16; k++) {
          _mm_storeu_pd(res + base + (k * 2), acc[k]);
        }
      }

      // -------------------------------------------------
      // PART B: The Tail (Specialized for 60, 61, 62, 63...)
      // -------------------------------------------------
      if constexpr (TAIL_SIZE > 0) {
        constexpr int base = NUM_FULL_STRIPS * STRIP_SIZE;
        constexpr int vecs_needed = TAIL_SIZE / 2;
        constexpr bool has_scalar = (TAIL_SIZE % 2 != 0);

        // These arrays are sized by constant expressions;
        // compiler allocates strictly registers, no malloc.
        __m128d acc[16];
        hyFloat scalar_acc = 0.0;

        // 1. Load Tail
        for (int k = 0; k < vecs_needed; k++) {
          acc[k] = _mm_loadu_pd(res + base + (k * 2));
        }
        if constexpr (has_scalar) {
          scalar_acc = res[base + TAIL_SIZE - 1];
        }

        // 2. Accumulate Tail
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;

          hyFloat val = theDataPtr[cxi];
          __m128d val_vec = _mm_set1_pd(val);

          // Unrolling logic for the tail
          int vec_offset = 0;
          int d_offset = 0;

          // Chunk 8 doubles (4 vectors) - heavily unrolled
          // For N=61, Tail=29. This runs 3 times (24 elements).
          constexpr int chunks_of_8 = TAIL_SIZE / 8;
          for (int k = 0; k < chunks_of_8; ++k) {
            __m128d C0 = _mm_loadu_pd(secArg + base + d_offset + 0);
            __m128d C1 = _mm_loadu_pd(secArg + base + d_offset + 2);
            __m128d C2 = _mm_loadu_pd(secArg + base + d_offset + 4);
            __m128d C3 = _mm_loadu_pd(secArg + base + d_offset + 6);

            acc[vec_offset + 0] =
                _mm_add_pd(acc[vec_offset + 0], _mm_mul_pd(val_vec, C0));
            acc[vec_offset + 1] =
                _mm_add_pd(acc[vec_offset + 1], _mm_mul_pd(val_vec, C1));
            acc[vec_offset + 2] =
                _mm_add_pd(acc[vec_offset + 2], _mm_mul_pd(val_vec, C2));
            acc[vec_offset + 3] =
                _mm_add_pd(acc[vec_offset + 3], _mm_mul_pd(val_vec, C3));

            vec_offset += 4;
            d_offset += 8;
          }

          // Chunk 2 doubles (1 vector)
          // For N=61, remaining=5. This runs 2 times (4 elements).
          constexpr int remaining_vectors = (TAIL_SIZE % 8) / 2;
          for (int k = 0; k < remaining_vectors; ++k) {
            __m128d C = _mm_loadu_pd(secArg + base + d_offset);
            acc[vec_offset] =
                _mm_add_pd(acc[vec_offset], _mm_mul_pd(val_vec, C));
            vec_offset++;
            d_offset += 2;
          }

          // Scalar
          if constexpr (has_scalar) {
            scalar_acc += val * secArg[base + TAIL_SIZE - 1];
          }
        }

        // 3. Store Tail
        for (int k = 0; k < vecs_needed; k++) {
          _mm_storeu_pd(res + base + (k * 2), acc[k]);
        }
        if constexpr (has_scalar) {
          res[base + TAIL_SIZE - 1] = scalar_acc;
        }
      }
    }

    res += N;
    currentXIndex = up;
  }
}

#elif defined _SLKP_USE_AVX_INTRINSICS

template <>
inline void _matrix_kernel_multiply_by_compressed_sparse_T<61>(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // N=61, AVX (256-bit, 4 doubles).
  // 61 doubles = 15.25 vectors.
  // We use 16 YMM registers available.
  // acc0..acc14 (15 vectors = 60 doubles).
  // scalar_acc (1 double).
  // Total 15 regs for vectors. leaves 1 for 'val' splat and temp loads.
  // This fits with minor register pressure.

  for (long i = 0; i < 61; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {
      __m256d acc0 = _mm256_loadu_pd(res + 0);
      __m256d acc1 = _mm256_loadu_pd(res + 4);
      __m256d acc2 = _mm256_loadu_pd(res + 8);
      __m256d acc3 = _mm256_loadu_pd(res + 12);
      __m256d acc4 = _mm256_loadu_pd(res + 16);
      __m256d acc5 = _mm256_loadu_pd(res + 20);
      __m256d acc6 = _mm256_loadu_pd(res + 24);
      __m256d acc7 = _mm256_loadu_pd(res + 28);
      __m256d acc8 = _mm256_loadu_pd(res + 32);
      __m256d acc9 = _mm256_loadu_pd(res + 36);
      __m256d acc10 = _mm256_loadu_pd(res + 40);
      __m256d acc11 = _mm256_loadu_pd(res + 44);
      __m256d acc12 = _mm256_loadu_pd(res + 48);
      __m256d acc13 = _mm256_loadu_pd(res + 52);
      __m256d acc14 = _mm256_loadu_pd(res + 56);
      hyFloat scalar_acc = res[60];

      for (long cxi = currentXIndex; cxi < up; cxi++) {
        long currentXColumn = compIdx[cxi + 61];
        const hyFloat *secArg = secondArgBase + currentXColumn * 61;
        hyFloat val = theDataPtr[cxi];
        __m256d val_vec = _mm256_set1_pd(val);
        __m256d C;

#ifdef _SLKP_USE_FMA3_INTRINSICS
        C = _mm256_loadu_pd(secArg + 0);
        acc0 = _mm256_fmadd_pd(val_vec, C, acc0);
        C = _mm256_loadu_pd(secArg + 4);
        acc1 = _mm256_fmadd_pd(val_vec, C, acc1);
        C = _mm256_loadu_pd(secArg + 8);
        acc2 = _mm256_fmadd_pd(val_vec, C, acc2);
        C = _mm256_loadu_pd(secArg + 12);
        acc3 = _mm256_fmadd_pd(val_vec, C, acc3);
        C = _mm256_loadu_pd(secArg + 16);
        acc4 = _mm256_fmadd_pd(val_vec, C, acc4);
        C = _mm256_loadu_pd(secArg + 20);
        acc5 = _mm256_fmadd_pd(val_vec, C, acc5);
        C = _mm256_loadu_pd(secArg + 24);
        acc6 = _mm256_fmadd_pd(val_vec, C, acc6);
        C = _mm256_loadu_pd(secArg + 28);
        acc7 = _mm256_fmadd_pd(val_vec, C, acc7);
        C = _mm256_loadu_pd(secArg + 32);
        acc8 = _mm256_fmadd_pd(val_vec, C, acc8);
        C = _mm256_loadu_pd(secArg + 36);
        acc9 = _mm256_fmadd_pd(val_vec, C, acc9);
        C = _mm256_loadu_pd(secArg + 40);
        acc10 = _mm256_fmadd_pd(val_vec, C, acc10);
        C = _mm256_loadu_pd(secArg + 44);
        acc11 = _mm256_fmadd_pd(val_vec, C, acc11);
        C = _mm256_loadu_pd(secArg + 48);
        acc12 = _mm256_fmadd_pd(val_vec, C, acc12);
        C = _mm256_loadu_pd(secArg + 52);
        acc13 = _mm256_fmadd_pd(val_vec, C, acc13);
        C = _mm256_loadu_pd(secArg + 56);
        acc14 = _mm256_fmadd_pd(val_vec, C, acc14);
#else
        C = _mm256_loadu_pd(secArg + 0);
        acc0 = _mm256_add_pd(acc0, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 4);
        acc1 = _mm256_add_pd(acc1, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 8);
        acc2 = _mm256_add_pd(acc2, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 12);
        acc3 = _mm256_add_pd(acc3, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 16);
        acc4 = _mm256_add_pd(acc4, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 20);
        acc5 = _mm256_add_pd(acc5, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 24);
        acc6 = _mm256_add_pd(acc6, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 28);
        acc7 = _mm256_add_pd(acc7, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 32);
        acc8 = _mm256_add_pd(acc8, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 36);
        acc9 = _mm256_add_pd(acc9, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 40);
        acc10 = _mm256_add_pd(acc10, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 44);
        acc11 = _mm256_add_pd(acc11, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 48);
        acc12 = _mm256_add_pd(acc12, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 52);
        acc13 = _mm256_add_pd(acc13, _mm256_mul_pd(val_vec, C));
        C = _mm256_loadu_pd(secArg + 56);
        acc14 = _mm256_add_pd(acc14, _mm256_mul_pd(val_vec, C));
#endif
        scalar_acc += val * secArg[60];
      }
      _mm256_storeu_pd(res + 0, acc0);
      _mm256_storeu_pd(res + 4, acc1);
      _mm256_storeu_pd(res + 8, acc2);
      _mm256_storeu_pd(res + 12, acc3);
      _mm256_storeu_pd(res + 16, acc4);
      _mm256_storeu_pd(res + 20, acc5);
      _mm256_storeu_pd(res + 24, acc6);
      _mm256_storeu_pd(res + 28, acc7);
      _mm256_storeu_pd(res + 32, acc8);
      _mm256_storeu_pd(res + 36, acc9);
      _mm256_storeu_pd(res + 40, acc10);
      _mm256_storeu_pd(res + 44, acc11);
      _mm256_storeu_pd(res + 48, acc12);
      _mm256_storeu_pd(res + 52, acc13);
      _mm256_storeu_pd(res + 56, acc14);
      res[60] = scalar_acc;
    }

    res += 61;
    currentXIndex = up;
  }
}

template <int N>
void _matrix_kernel_multiply_by_compressed_sparse_T(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // Compute layout constants at compile time
  constexpr int STRIP_SIZE = 32; // 8 vectors of 4 doubles
  constexpr int NUM_FULL_STRIPS = N / STRIP_SIZE;
  constexpr int TAIL_SIZE = N % STRIP_SIZE;

  for (long i = 0; i < N; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {

// -------------------------------------------------
// PART A: Compile-time loop for Full Strips
// -------------------------------------------------
#pragma unroll
      for (int s = 0; s < NUM_FULL_STRIPS; s++) {
        // We know exactly where we are: s * 32
        const int base = s * STRIP_SIZE;

        __m256d acc[8];

        // 1. Load Accumulators
        for (int k = 0; k < 8; k++) {
          acc[k] = _mm256_loadu_pd(res + base + (k * 4));
        }

        // 2. Accumulate
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N]; // constant N
          const hyFloat *secArg = secondArgBase + currentXColumn * N;

          __m256d val_vec = _mm256_set1_pd(theDataPtr[cxi]);

          for (int k = 0; k < 8; k++) {
            __m256d C = _mm256_loadu_pd(secArg + base + (k * 4));
#ifdef _SLKP_USE_FMA3_INTRINSICS
            acc[k] = _mm256_fmadd_pd(val_vec, C, acc[k]);
#else
            acc[k] = _mm256_add_pd(acc[k], _mm256_mul_pd(val_vec, C));
#endif
          }
        }

        // 3. Store Accumulators
        for (int k = 0; k < 8; k++) {
          _mm256_storeu_pd(res + base + (k * 4), acc[k]);
        }
      }

      // -------------------------------------------------
      // PART B: The Tail
      // -------------------------------------------------
      if constexpr (TAIL_SIZE > 0) {
        constexpr int base = NUM_FULL_STRIPS * STRIP_SIZE;
        constexpr int vecs_needed = TAIL_SIZE / 4;
        constexpr int remainder_scalars = TAIL_SIZE % 4;

        // These arrays are sized by constant expressions
        __m256d acc[8];
        hyFloat scalar_acc[3];

        // 1. Load Tail
        for (int k = 0; k < vecs_needed; k++) {
          acc[k] = _mm256_loadu_pd(res + base + (k * 4));
        }
        for (int k = 0; k < remainder_scalars; k++) {
          scalar_acc[k] = res[base + vecs_needed * 4 + k];
        }

        // 2. Accumulate Tail
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;

          hyFloat val = theDataPtr[cxi];
          __m256d val_vec = _mm256_set1_pd(val);

          // Vectors
          for (int k = 0; k < vecs_needed; ++k) {
            __m256d C = _mm256_loadu_pd(secArg + base + (k * 4));
#ifdef _SLKP_USE_FMA3_INTRINSICS
            acc[k] = _mm256_fmadd_pd(val_vec, C, acc[k]);
#else
            acc[k] = _mm256_add_pd(acc[k], _mm256_mul_pd(val_vec, C));
#endif
          }

          // Scalars
          for (int k = 0; k < remainder_scalars; ++k) {
            scalar_acc[k] += val * secArg[base + vecs_needed * 4 + k];
          }
        }

        // 3. Store Tail
        for (int k = 0; k < vecs_needed; k++) {
          _mm256_storeu_pd(res + base + (k * 4), acc[k]);
        }
        for (int k = 0; k < remainder_scalars; k++) {
          res[base + vecs_needed * 4 + k] = scalar_acc[k];
        }
      }
    }

    res += N;
    currentXIndex = up;
  }
}

#else

// Non-SIMD implementation of _matrix_kernel_multiply_by_compressed_sparse_T
template <int N>
void _matrix_kernel_multiply_by_compressed_sparse_T(
    hyFloat *__restrict__ resBase, const long *__restrict__ compIdx,
    const hyFloat *__restrict__ secondArgBase,
    const hyFloat *__restrict__ theDataPtr) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  constexpr int STRIP_SIZE = 32;
  constexpr int NUM_FULL_STRIPS = N / STRIP_SIZE;
  constexpr int TAIL_SIZE = N % STRIP_SIZE;

  for (long i = 0; i < N; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {
      // Full strips
      for (int s = 0; s < NUM_FULL_STRIPS; s++) {
        const int base = s * STRIP_SIZE;
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;
          hyFloat val = theDataPtr[cxi];

          for (int k = 0; k < STRIP_SIZE; k++) {
            res[base + k] += val * secArg[base + k];
          }
        }
      }

      // Tail
      if constexpr (TAIL_SIZE > 0) {
        constexpr int base = NUM_FULL_STRIPS * STRIP_SIZE;
        for (long cxi = currentXIndex; cxi < up; cxi++) {
          long currentXColumn = compIdx[cxi + N];
          const hyFloat *secArg = secondArgBase + currentXColumn * N;
          hyFloat val = theDataPtr[cxi];

          for (int k = 0; k < TAIL_SIZE; k++) {
            res[base + k] += val * secArg[base + k];
          }
        }
      }
    }
    res += N;
    currentXIndex = up;
  }
}

#endif

void _matrix_kernel_multiply_by_compressed_sparse(hyFloat *res,
                                                  const long *compIdx,
                                                  const hyFloat *secondArgBase,
                                                  const hyFloat *theDataPtr,
                                                  long N);

#endif