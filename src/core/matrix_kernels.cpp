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

#include "matrix_kernels.h"

/**** START: ARM NEON CODE *****/

#ifdef _SLKP_USE_ARM_NEON

void _matrix_kernel_multiply_by_compressed_sparse(hyFloat *resBase,
                                                  const long *compIdx,
                                                  const hyFloat *secondArgBase,
                                                  const hyFloat *theDataPtr,
                                                  long N) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // Compute layout constants at compile time
  constexpr long STRIP_SIZE = 32; // 16 vectors
  long NUM_FULL_STRIPS = N / STRIP_SIZE;
  long TAIL_SIZE = N % STRIP_SIZE;

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
      if (TAIL_SIZE > 0) {
        long base = NUM_FULL_STRIPS * STRIP_SIZE;
        long int vecs_needed = TAIL_SIZE >> 1;
        bool has_scalar = (TAIL_SIZE % 2 != 0);

        // These arrays are sized by constant expressions;
        // compiler allocates strictly registers, no malloc.
        float64x2_t acc[STRIP_SIZE >> 1]; // maximum needed
        hyFloat scalar_acc = 0.0;

        // 1. Load Tail
        for (int k = 0; k < vecs_needed; k++) {
          acc[k] = vld1q_f64(res + base + (k * 2));
        }
        if (has_scalar) {
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
          long chunks_of_8 = TAIL_SIZE / 8;
          for (long k = 0; k < chunks_of_8; ++k) {
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
          long remaining_vectors = (TAIL_SIZE % 8) / 2;
          for (long k = 0; k < remaining_vectors; ++k) {
            float64x2_t C = vld1q_f64(secArg + base + d_offset);
            acc[vec_offset] = vfmaq_f64(acc[vec_offset], val_vec, C);
            vec_offset++;
            d_offset += 2;
          }

          // Scalar
          if (has_scalar) {
            scalar_acc += val * secArg[base + TAIL_SIZE - 1];
          }
        }

        // 3. Store Tail
        for (int k = 0; k < vecs_needed; k++) {
          vst1q_f64(res + base + (k * 2), acc[k]);
        }
        if (has_scalar) {
          res[base + TAIL_SIZE - 1] = scalar_acc;
        }
      }
    }

    res += N;
    currentXIndex = up;
  }
}

#endif

/**** END: ARM NEON CODE *****/

/**** START: SSE CODE *****/

#ifdef _SLKP_USE_SSE_INTRINSICS

void _matrix_kernel_multiply_by_compressed_sparse(hyFloat *resBase,
                                                  const long *compIdx,
                                                  const hyFloat *secondArgBase,
                                                  const hyFloat *theDataPtr,
                                                  long N) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // Compute layout constants at compile time
  constexpr long STRIP_SIZE = 32; // 16 vectors
  long NUM_FULL_STRIPS = N / STRIP_SIZE;
  long TAIL_SIZE = N % STRIP_SIZE;

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
      if (TAIL_SIZE > 0) {
        long base = NUM_FULL_STRIPS * STRIP_SIZE;
        long int vecs_needed = TAIL_SIZE >> 1;
        bool has_scalar = (TAIL_SIZE % 2 != 0);

        // These arrays are sized by constant expressions;
        // compiler allocates strictly registers, no malloc.
        __m128d acc[STRIP_SIZE >> 1]; // maximum needed
        hyFloat scalar_acc = 0.0;

        // 1. Load Tail
        for (int k = 0; k < vecs_needed; k++) {
          acc[k] = _mm_loadu_pd(res + base + (k * 2));
        }
        if (has_scalar) {
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
          long chunks_of_8 = TAIL_SIZE / 8;
          for (long k = 0; k < chunks_of_8; ++k) {
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
          long remaining_vectors = (TAIL_SIZE % 8) / 2;
          for (long k = 0; k < remaining_vectors; ++k) {
            __m128d C = _mm_loadu_pd(secArg + base + d_offset);
            acc[vec_offset] =
                _mm_add_pd(acc[vec_offset], _mm_mul_pd(val_vec, C));
            vec_offset++;
            d_offset += 2;
          }

          // Scalar
          if (has_scalar) {
            scalar_acc += val * secArg[base + TAIL_SIZE - 1];
          }
        }

        // 3. Store Tail
        for (int k = 0; k < vecs_needed; k++) {
          _mm_storeu_pd(res + base + (k * 2), acc[k]);
        }
        if (has_scalar) {
          res[base + TAIL_SIZE - 1] = scalar_acc;
        }
      }
    }

    res += N;
    currentXIndex = up;
  }
}

#endif

/**** END: SSE CODE *****/

/**** START: AVX CODE *****/

#ifdef _SLKP_USE_AVX_INTRINSICS

void _matrix_kernel_multiply_by_compressed_sparse(hyFloat *resBase,
                                                  const long *compIdx,
                                                  const hyFloat *secondArgBase,
                                                  const hyFloat *theDataPtr,
                                                  long N) {
  long currentXIndex = 0L;
  hyFloat *res = resBase;

  // Compute layout constants
  constexpr long STRIP_SIZE = 32; // 8 vectors of 4 doubles
  long NUM_FULL_STRIPS = N / STRIP_SIZE;
  long TAIL_SIZE = N % STRIP_SIZE;

  for (long i = 0; i < N; i++) {
    long up = compIdx[i];

    if (currentXIndex < up) {

// -------------------------------------------------
// PART A: Full Strips
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
      if (TAIL_SIZE > 0) {
        long base = NUM_FULL_STRIPS * STRIP_SIZE;
        long vecs_needed = TAIL_SIZE / 4;
        long remainder_scalars = TAIL_SIZE % 4;

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

#endif

/**** END: AVX CODE *****/
