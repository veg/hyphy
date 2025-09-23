#include <arm_neon.h>
#include <stdio.h>

double _neon_sum_2(float64x2_t const x) {
  return vgetq_lane_f64(x, 0) + vgetq_lane_f64(x, 1);
}

int main(void) {
  double transitionMatrix[16];
  // void * transposed_mx = NULL;
  for (int i = 0; i < 16; i++) {
    transitionMatrix[i] = i * 0.1;
  }

  float64x2x2_t tmatrix_transpose[4] = {
      (float64x2x2_t){transitionMatrix[0], transitionMatrix[4],
                      transitionMatrix[8], transitionMatrix[12]},
      (float64x2x2_t){transitionMatrix[1], transitionMatrix[5],
                      transitionMatrix[9], transitionMatrix[13]},
      (float64x2x2_t){transitionMatrix[2], transitionMatrix[6],
                      transitionMatrix[10], transitionMatrix[14]},
      (float64x2x2_t){transitionMatrix[3], transitionMatrix[7],
                      transitionMatrix[11], transitionMatrix[15]}};

  double childVector[4] = {1. / 2., 1. / 4., 1. / 8., 1. / 16.};
  double parentConditionals[4] = {1.1, 1.2, 1.3, 1.4};

  float64x2_t c0 = vdupq_n_f64(childVector[0]),
              c1 = vdupq_n_f64(childVector[1]),
              c2 = vdupq_n_f64(childVector[2]),
              c3 = vdupq_n_f64(childVector[3]);

  float64x2x2_t t0, t1, t2, t3;

  void *transposed_mx = tmatrix_transpose;

  if (transposed_mx) {
    t0 = ((float64x2x2_t *)transposed_mx)[0];
    t1 = ((float64x2x2_t *)transposed_mx)[1];
    t2 = ((float64x2x2_t *)transposed_mx)[2];
    t3 = ((float64x2x2_t *)transposed_mx)[3];
  } else {
    t0 = (float64x2x2_t){transitionMatrix[0], transitionMatrix[4],
                         transitionMatrix[8], transitionMatrix[12]};
    t1 = (float64x2x2_t){transitionMatrix[1], transitionMatrix[5],
                         transitionMatrix[9], transitionMatrix[13]};
    t2 = (float64x2x2_t){transitionMatrix[2], transitionMatrix[6],
                         transitionMatrix[10], transitionMatrix[14]};
    t3 = (float64x2x2_t){transitionMatrix[3], transitionMatrix[7],
                         transitionMatrix[11], transitionMatrix[15]};
  }

  t0.val[0] = vfmaq_f64(vmulq_f64(c1, t1.val[0]), c0, t0.val[0]);
  t0.val[1] = vfmaq_f64(vmulq_f64(c1, t1.val[1]), c0, t0.val[1]);
  t2.val[0] = vfmaq_f64(vmulq_f64(c3, t3.val[0]), c2, t2.val[0]);
  t2.val[1] = vfmaq_f64(vmulq_f64(c3, t3.val[1]), c2, t2.val[1]);

  float64x2_t pv0 = vld1q_f64(parentConditionals),
              pv1 = vld1q_f64(parentConditionals + 2);

  pv0 = vmulq_f64(pv0, vaddq_f64(t0.val[0], t2.val[0]));
  pv1 = vmulq_f64(pv1, vaddq_f64(t0.val[1], t2.val[1]));
  // pv.val[0] = vaddq_f64 (t0.val[0], t2.val[0]);
  // pv.val[1] = vaddq_f64 (t0.val[1], t2.val[1]);

  double check[4];

  check[0] = parentConditionals[0] * (transitionMatrix[0] * childVector[0] +
                                      transitionMatrix[1] * childVector[1] +
                                      transitionMatrix[2] * childVector[2] +
                                      transitionMatrix[3] * childVector[3]);
  check[1] = parentConditionals[1] * (transitionMatrix[4] * childVector[0] +
                                      transitionMatrix[5] * childVector[1] +
                                      transitionMatrix[6] * childVector[2] +
                                      transitionMatrix[7] * childVector[3]);
  check[2] = parentConditionals[2] * (transitionMatrix[8] * childVector[0] +
                                      transitionMatrix[9] * childVector[1] +
                                      transitionMatrix[10] * childVector[2] +
                                      transitionMatrix[11] * childVector[3]);
  check[3] = parentConditionals[3] * (transitionMatrix[12] * childVector[0] +
                                      transitionMatrix[13] * childVector[1] +
                                      transitionMatrix[14] * childVector[2] +
                                      transitionMatrix[15] * childVector[3]);

  vst1q_f64(parentConditionals, pv0);
  vst1q_f64(parentConditionals + 2, pv1);

  for (int i = 0; i < 4; i++) {
    printf("%g => %g\n", parentConditionals[i], check[i]);
  }

  printf("\n%g %g\n", _neon_sum_2(pv1),
         parentConditionals[2] + parentConditionals[3]);

  /*

  float64x2x2_t tmatrix_transpose [4] = {
      (float64x2x2_t)
  {transitionMatrix[0],transitionMatrix[4],transitionMatrix[8],transitionMatrix[12]},
      (float64x2x2_t)
  {transitionMatrix[1],transitionMatrix[5],transitionMatrix[9],transitionMatrix[13]},
      (float64x2x2_t)
  {transitionMatrix[2],transitionMatrix[6],transitionMatrix[10],transitionMatrix[14]},
      (float64x2x2_t)
  {transitionMatrix[3],transitionMatrix[7],transitionMatrix[11],transitionMatrix[15]}
  };

  vst1q_f64_x2  (transitionMatrix, tmatrix_transpose[0]);
  vst1q_f64_x2  (transitionMatrix+4, tmatrix_transpose[1]);
  vst1q_f64_x2  (transitionMatrix+8, tmatrix_transpose[2]);
  vst1q_f64_x2  (transitionMatrix+12, tmatrix_transpose[3]);

  for (int i = 0; i < 16; i++) {
      printf ("%d %g\n", i, transitionMatrix[i]);
  }

  float64x2_t a = {1.,2.},
              b = {3.,4.},
              c = {1.,1.};

  double res[2];

  c = vfmaq_f64 (c, a, b);

  vst1q_f64 (res, c);

  printf ("\n%g %g\n", res[0], res[1]);

  */

  return 0;
}