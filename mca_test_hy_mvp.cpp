#include <arm_neon.h>

#ifndef __restrict__
#define __restrict__
#endif

// Copied from tree_evaluator.cpp (NEON SIMD versions)

void _mx_vect_4x4(float64x2x2_t &cv, double const *__restrict M,
                  double const *__restrict V, int stride) {

  const float64x2x2_t v = vld1q_f64_x2(V);
  const float64x2x2_t m_row1 = vld1q_f64_x2(M);
  const float64x2x2_t m_row2 = vld1q_f64_x2(M + stride);
  const float64x2x2_t m_row3 = vld1q_f64_x2(M + (stride << 1));
  const float64x2x2_t m_row4 = vld1q_f64_x2(M + (stride << 1) + stride);

  float64x2_t dp1 = vmulq_f64(m_row1.val[0], v.val[0]);
  float64x2_t dp2 = vmulq_f64(m_row2.val[0], v.val[0]);

  dp1 = vfmaq_f64(dp1, m_row1.val[1], v.val[1]);
  dp2 = vfmaq_f64(dp2, m_row2.val[1], v.val[1]);

  float64x2_t dp3 = vmulq_f64(m_row3.val[0], v.val[0]);
  float64x2_t dp4 = vmulq_f64(m_row4.val[0], v.val[0]);

  dp3 = vfmaq_f64(dp3, m_row3.val[1], v.val[1]);
  dp4 = vfmaq_f64(dp4, m_row4.val[1], v.val[1]);

  cv.val[0] = vaddq_f64(vzip1q_f64(dp1, dp2), vzip2q_f64(dp1, dp2));
  cv.val[1] = vaddq_f64(vzip1q_f64(dp3, dp4), vzip2q_f64(dp3, dp4));
}

void _mx_vect_4x4_add(float64x2x2_t &cv, double const *__restrict M,
                      double const *__restrict V, int stride) {

  const float64x2x2_t v = vld1q_f64_x2(V);
  const float64x2x2_t m_row1 = vld1q_f64_x2(M);
  const float64x2x2_t m_row2 = vld1q_f64_x2(M + stride);
  const float64x2x2_t m_row3 = vld1q_f64_x2(M + (stride << 1));
  const float64x2x2_t m_row4 = vld1q_f64_x2(M + (stride << 1) + stride);

  float64x2_t dp1 = vmulq_f64(m_row1.val[0], v.val[0]);
  float64x2_t dp2 = vmulq_f64(m_row2.val[0], v.val[0]);

  dp1 = vfmaq_f64(dp1, m_row1.val[1], v.val[1]);
  dp2 = vfmaq_f64(dp2, m_row2.val[1], v.val[1]);

  float64x2_t dp3 = vmulq_f64(m_row3.val[0], v.val[0]);
  float64x2_t dp4 = vmulq_f64(m_row4.val[0], v.val[0]);

  dp3 = vfmaq_f64(dp3, m_row3.val[1], v.val[1]);
  dp4 = vfmaq_f64(dp4, m_row4.val[1], v.val[1]);

  const float64x2_t sum12 =
      vaddq_f64(vzip1q_f64(dp1, dp2), vzip2q_f64(dp1, dp2));
  const float64x2_t sum34 =
      vaddq_f64(vzip1q_f64(dp3, dp4), vzip2q_f64(dp3, dp4));

  cv.val[0] = vaddq_f64(cv.val[0], sum12);
  cv.val[1] = vaddq_f64(cv.val[1], sum34);
}

template <int D>
void _hy_mvp_blocked_4x4(double *C, double const *M, double const *V) {
  auto offset = [](int i, int j) -> int { return (i << 2) * D + (j << 2); };

  int blocks = D >> 2;
  int remainder = D - (blocks << 2);

  switch (remainder) {

  case 0: {
    for (int i = 0; i < blocks; i++) {
      float64x2x2_t accumulator;
      _mx_vect_4x4(accumulator, M + offset(i, 0), V, D);
      for (int j = 1; j < blocks; j++) {
        _mx_vect_4x4_add(accumulator, M + offset(i, j), V + offset(0, j), D);
      }

      vst1q_f64(C + offset(0, i), accumulator.val[0]);
      vst1q_f64(C + offset(0, i) + 2, accumulator.val[1]);
      // handle_small_product_add <4,1> (C + offset (0,i), M + offset
      // (i,blocks), V + offset (0,blocks), D);
    }
    break;
  }

  case 1: {
    for (int i = 0; i < blocks; i++) {
      float64x2x2_t accumulator;
      _mx_vect_4x4(accumulator, M + offset(i, 0), V, D);

      for (int j = 1; j < blocks; j++) {
        _mx_vect_4x4_add(accumulator, M + offset(i, j), V + offset(0, j), D);
      }

      int moffset = offset(i, blocks);

      float64x2_t v = vdupq_n_f64(V[offset(0, blocks)]),
                  m1 = vcombine_f64(vld1_f64(M + moffset),
                                    vld1_f64(M + moffset + D)),
                  m2 = vcombine_f64(vld1_f64(M + moffset + 2 * D),
                                    vld1_f64(M + moffset + 3 * D));

      accumulator.val[0] = vfmaq_f64(accumulator.val[0], m1, v);
      accumulator.val[1] = vfmaq_f64(accumulator.val[1], m2, v);

      vst1q_f64(C + offset(0, i), accumulator.val[0]);
      vst1q_f64(C + offset(0, i) + 2, accumulator.val[1]);

      // handle_small_product_add <4,1> (C + offset (0,i), M + offset
      // (i,blocks), V + offset (0,blocks), D);
    }
    M += offset(blocks, 0);
    float64x2x2_t accumulator, row = vld1q_f64_x2(M), col = vld1q_f64_x2(V);

    accumulator.val[0] = vmulq_f64(row.val[0], col.val[0]);
    accumulator.val[1] = vmulq_f64(row.val[1], col.val[1]);

    for (int j = 1; j < blocks; j++) {
      row = vld1q_f64_x2(M + (j << 2));
      col = vld1q_f64_x2(V + (j << 2));

      accumulator.val[0] =
          vfmaq_f64(accumulator.val[0], row.val[0], col.val[0]);
      accumulator.val[1] =
          vfmaq_f64(accumulator.val[1], row.val[1], col.val[1]);
    }

    C[D - 1] = M[D - 1] * V[D - 1] +
               vaddvq_f64(vpaddq_f64(accumulator.val[0], accumulator.val[1]));
    break;
  }

  case 2: {

    for (int i = 0; i < blocks; i++) {
      float64x2x2_t accumulator;
      _mx_vect_4x4(accumulator, M + offset(i, 0), V, D);
      for (int j = 1; j < blocks; j++) {
        _mx_vect_4x4_add(accumulator, M + offset(i, j), V + offset(0, j), D);
      }

      int moffset = offset(i, blocks);
      float64x2_t mv1 = vcombine_f64(vld1_f64(M + moffset),
                                     vld1_f64(M + moffset + D)),
                  mv2 = vcombine_f64(vld1_f64(M + moffset + 2 * D),
                                     vld1_f64(M + moffset + 3 * D));

      float64x2_t v = vdupq_n_f64(V[offset(0, blocks)]);

      accumulator.val[0] = vfmaq_f64(accumulator.val[0], mv1, v);
      accumulator.val[1] = vfmaq_f64(accumulator.val[1], mv2, v);

      mv1 = vcombine_f64(vld1_f64(M + moffset + 1),
                         vld1_f64(M + moffset + D + 1));

      mv2 = vcombine_f64(vld1_f64(M + moffset + 2 * D + 1),
                         vld1_f64(M + moffset + 3 * D + 1));

      v = vdupq_n_f64(V[offset(0, blocks) + 1]);

      accumulator.val[0] = vfmaq_f64(accumulator.val[0], mv1, v);
      vst1q_f64(C + offset(0, i), accumulator.val[0]);
      accumulator.val[1] = vfmaq_f64(accumulator.val[1], mv2, v);
      vst1q_f64(C + offset(0, i) + 2, accumulator.val[1]);
      // handle_small_product_add <4,1> (C + offset (0,i), M + offset
      // (i,blocks), V + offset (0,blocks), D);
    }

    M += offset(blocks, 0);
    float64x2x2_t accumulator, accumulator2, row = vld1q_f64_x2(M),
                                             row2 = vld1q_f64_x2(M + D),
                                             col = vld1q_f64_x2(V);

    accumulator.val[0] = vmulq_f64(row.val[0], col.val[0]);
    accumulator.val[1] = vmulq_f64(row.val[1], col.val[1]);
    accumulator2.val[0] = vmulq_f64(row2.val[0], col.val[0]);
    accumulator2.val[1] = vmulq_f64(row2.val[1], col.val[1]);

    for (int j = 1; j < blocks; j++) {
      row = vld1q_f64_x2(M + (j << 2));
      row2 = vld1q_f64_x2(M + D + (j << 2));
      col = vld1q_f64_x2(V + (j << 2));

      accumulator.val[0] =
          vfmaq_f64(accumulator.val[0], row.val[0], col.val[0]);
      accumulator.val[1] =
          vfmaq_f64(accumulator.val[1], row.val[1], col.val[1]);
      accumulator2.val[0] =
          vfmaq_f64(accumulator2.val[0], row2.val[0], col.val[0]);
      accumulator2.val[1] =
          vfmaq_f64(accumulator2.val[1], row2.val[1], col.val[1]);
    }
    C[D - 2] = M[D - 2] * V[D - 2] + M[D - 1] * V[D - 1] +
               vaddvq_f64(vpaddq_f64(accumulator.val[0], accumulator.val[1]));
    C[D - 1] = M[2 * D - 2] * V[D - 2] + M[2 * D - 1] * V[D - 1] +
               vaddvq_f64(vpaddq_f64(accumulator2.val[0], accumulator2.val[1]));

    break;
  }

  case 3: {
    for (int i = 0; i < blocks; i++) {
      float64x2x2_t accumulator;
      _mx_vect_4x4(accumulator, M + offset(i, 0), V, D);
      for (int j = 1; j < blocks; j++) {
        _mx_vect_4x4_add(accumulator, M + offset(i, j), V + offset(0, j), D);
      }

      int moffset = offset(i, blocks);
      float64x2_t mv1 = vcombine_f64(vld1_f64(M + moffset),
                                     vld1_f64(M + moffset + D)),
                  mv2 = vcombine_f64(vld1_f64(M + moffset + 2 * D),
                                     vld1_f64(M + moffset + 3 * D));

      float64x2_t v = vdupq_n_f64(V[offset(0, blocks)]);
      accumulator.val[0] = vfmaq_f64(accumulator.val[0], mv1, v);
      accumulator.val[1] = vfmaq_f64(accumulator.val[1], mv2, v);

      mv1 = vcombine_f64(vld1_f64(M + moffset + 1),
                         vld1_f64(M + moffset + D + 1));

      mv2 = vcombine_f64(vld1_f64(M + moffset + 2 * D + 1),
                         vld1_f64(M + moffset + 3 * D + 1));

      v = vdupq_n_f64(V[offset(0, blocks) + 1]);
      accumulator.val[0] = vfmaq_f64(accumulator.val[0], mv1, v);
      accumulator.val[1] = vfmaq_f64(accumulator.val[1], mv2, v);

      mv1 = vcombine_f64(vld1_f64(M + moffset + 2),
                         vld1_f64(M + moffset + D + 2));

      mv2 = vcombine_f64(vld1_f64(M + moffset + 2 * D + 2),
                         vld1_f64(M + moffset + 3 * D + 2));

      v = vdupq_n_f64(V[offset(0, blocks) + 2]);
      accumulator.val[0] = vfmaq_f64(accumulator.val[0], mv1, v);
      vst1q_f64(C + offset(0, i), accumulator.val[0]);
      accumulator.val[1] = vfmaq_f64(accumulator.val[1], mv2, v);
      vst1q_f64(C + offset(0, i) + 2, accumulator.val[1]);
      // handle_small_product_add <4,1> (C + offset (0,i), M + offset
      // (i,blocks), V + offset (0,blocks), D);
    }

    M += offset(blocks, 0);
    float64x2x2_t accumulator, accumulator2, accumulator3,
        row = vld1q_f64_x2(M), row2 = vld1q_f64_x2(M + D),
        row3 = vld1q_f64_x2(M + 2 * D), col = vld1q_f64_x2(V);

    accumulator.val[0] = vmulq_f64(row.val[0], col.val[0]);
    accumulator.val[1] = vmulq_f64(row.val[1], col.val[1]);
    accumulator2.val[0] = vmulq_f64(row2.val[0], col.val[0]);
    accumulator2.val[1] = vmulq_f64(row2.val[1], col.val[1]);
    accumulator3.val[0] = vmulq_f64(row3.val[0], col.val[0]);
    accumulator3.val[1] = vmulq_f64(row3.val[1], col.val[1]);

    for (int j = 1; j < blocks; j++) {
      row = vld1q_f64_x2(M + (j << 2));
      row2 = vld1q_f64_x2(M + D + (j << 2));
      row3 = vld1q_f64_x2(M + 2 * D + (j << 2));
      col = vld1q_f64_x2(V + (j << 2));

      accumulator.val[0] =
          vfmaq_f64(accumulator.val[0], row.val[0], col.val[0]);
      accumulator.val[1] =
          vfmaq_f64(accumulator.val[1], row.val[1], col.val[1]);
      accumulator2.val[0] =
          vfmaq_f64(accumulator2.val[0], row2.val[0], col.val[0]);
      accumulator2.val[1] =
          vfmaq_f64(accumulator2.val[1], row2.val[1], col.val[1]);
      accumulator3.val[0] =
          vfmaq_f64(accumulator3.val[0], row3.val[0], col.val[0]);
      accumulator3.val[1] =
          vfmaq_f64(accumulator3.val[1], row3.val[1], col.val[1]);
    }
    C[D - 3] = M[D - 3] * V[D - 3] + M[D - 2] * V[D - 2] + M[D - 1] * V[D - 1] +
               vaddvq_f64(vpaddq_f64(accumulator.val[0], accumulator.val[1]));
    ;
    C[D - 2] = M[2 * D - 3] * V[D - 3] + M[2 * D - 2] * V[D - 2] +
               M[2 * D - 1] * V[D - 1] +
               vaddvq_f64(vpaddq_f64(accumulator2.val[0], accumulator2.val[1]));
    C[D - 1] = M[3 * D - 3] * V[D - 3] + M[3 * D - 2] * V[D - 2] +
               M[3 * D - 1] * V[D - 1] +
               vaddvq_f64(vpaddq_f64(accumulator3.val[0], accumulator3.val[1]));

    break;
  }
  }
}

int main() {
  double C[60 * 60];
  double M[60 * 60];
  double V[60];
  _hy_mvp_blocked_4x4<60>(C, M, V);
  return 0;
}
