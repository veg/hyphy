#include <arm_neon.h>

#ifndef __restrict__
#define __restrict__
#endif

void _mx_vect_4x4_add(float64x2x2_t &cv, double const *__restrict M,
                      double const *__restrict V, int stride) {

  // Load the entire vector V into two registers {V0,V1} and {V2,V3}
  const float64x2x2_t v = vld1q_f64_x2(V);

  // Group loads to allow the CPU to hide memory latency
  const float64x2x2_t m_row1 = vld1q_f64_x2(M);
  const float64x2x2_t m_row2 = vld1q_f64_x2(M + stride);
  const float64x2x2_t m_row3 = vld1q_f64_x2(M + (stride << 1));
  const float64x2x2_t m_row4 = vld1q_f64_x2(M + (stride << 1) + stride);

  float64x2_t dp1 = vmulq_f64(m_row1.val[0], v.val[0]);
  float64x2_t dp2 = vmulq_f64(m_row2.val[0], v.val[0]);

  dp1 = vfmaq_f64(dp1, m_row1.val[1],
                  v.val[1]); // FMA: dp1 += m_row1.val[1] * v.val[1]

  dp2 = vfmaq_f64(dp2, m_row2.val[1],
                  v.val[1]); // FMA: dp2 += m_row2.val[1] * v.val[1]

  // --- Calculate dot product parts for each row using FMA ---
  // The pattern is: dp = (m_row.part1 * v.part1) + (m_row.part2 * v.part2)
  // We compute the first product, then use FMA for the second.

  // Process rows 1 & 2
  float64x2_t dp3 = vmulq_f64(m_row3.val[0], v.val[0]);
  float64x2_t dp4 = vmulq_f64(m_row4.val[0], v.val[0]);

  // Process rows 3 & 4
  dp3 = vfmaq_f64(dp3, m_row3.val[1],
                  v.val[1]); // FMA: dp3 += m_row3.val[1] * v.val[1]
  dp4 = vfmaq_f64(dp4, m_row4.val[1],
                  v.val[1]); // FMA: dp4 += m_row4.val[1] * v.val[1]

  // --- Horizontally sum parts and accumulate into the result vector cv ---
  // The `vadd(vzip1, vzip2)` pattern efficiently performs the horizontal sum.
  const float64x2_t sum12 =
      vaddq_f64(vzip1q_f64(dp1, dp2), vzip2q_f64(dp1, dp2));

  const float64x2_t sum34 =
      vaddq_f64(vzip1q_f64(dp3, dp4), vzip2q_f64(dp3, dp4));

  cv.val[0] = vaddq_f64(cv.val[0], sum12);
  cv.val[1] = vaddq_f64(cv.val[1], sum34);
}

int main() {
  double M[16], V[4];
  float64x2x2_t cv;
  cv.val[0] = vdupq_n_f64(0.0);
  cv.val[1] = vdupq_n_f64(0.0);
  _mx_vect_4x4_add(cv, M, V, 4);
  return 0;
}
