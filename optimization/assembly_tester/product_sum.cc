#include <arm_neon.h>

typedef double hyFloat;
#define _hprestrict_ __restrict

template <long D, long S>
inline float64x2_t
__ll_handle_block10_product_sum(hyFloat const *_hprestrict_ transposedMatrix,
                                hyFloat *_hprestrict_ childVector,
                                hyFloat *_hprestrict_ parentConditionals,
                                float64x2_t *grandTotal) {

  float64x2_t P[5] = {vdupq_n_f64(0.), vdupq_n_f64(0.), vdupq_n_f64(0.),
                      vdupq_n_f64(0.), vdupq_n_f64(0.)};

  hyFloat const *__restrict tT = transposedMatrix;

#pragma unroll 2
  for (long col_idx = 0; col_idx < D; col_idx++, tT += D) {
    if (childVector[col_idx] > 0.) {
      float64x2_t C = vdupq_n_f64(childVector[col_idx]);

      P[0] = vfmaq_f64(P[0], vld1q_f64(tT + S), C);
      P[1] = vfmaq_f64(P[1], vld1q_f64(tT + S + 2), C);
      P[2] = vfmaq_f64(P[2], vld1q_f64(tT + S + 4), C);
      P[3] = vfmaq_f64(P[3], vld1q_f64(tT + S + 6), C);
      P[4] = vfmaq_f64(P[4], vld1q_f64(tT + S + 8), C);
    }
  }

  P[0] = vmulq_f64(vld1q_f64(parentConditionals + S), P[0]);
  P[1] = vmulq_f64(vld1q_f64(parentConditionals + S + 2), P[1]);
  P[2] = vmulq_f64(vld1q_f64(parentConditionals + S + 4), P[2]);
  P[3] = vmulq_f64(vld1q_f64(parentConditionals + S + 6), P[3]);
  P[4] = vmulq_f64(vld1q_f64(parentConditionals + S + 8), P[4]);

  vst1q_f64(parentConditionals + S, P[0]);
  vst1q_f64(parentConditionals + S + 2, P[1]);
  vst1q_f64(parentConditionals + S + 4, P[2]);
  vst1q_f64(parentConditionals + S + 6, P[3]);
  vst1q_f64(parentConditionals + S + 8, P[4]);

  P[0] = vaddq_f64(P[0], P[1]);
  P[2] = vaddq_f64(P[2], P[3]);
  P[0] = vaddq_f64(P[0], P[4]);
  P[2] = vaddq_f64(P[0], P[2]);
  if (grandTotal) {
    *grandTotal = vaddq_f64(P[2], *grandTotal);
    return *grandTotal;
  }
  return P[2];
}

template float64x2_t __ll_handle_block10_product_sum<61L, 0L>(
    hyFloat const *_hprestrict_ transposedMatrix,
    hyFloat *_hprestrict_ childVector, hyFloat *_hprestrict_ parentConditionals,
    float64x2_t *grandTotal);