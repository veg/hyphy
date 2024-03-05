//
//  matrix_mult.cpp
//  HyPhy
//
//  Created by Sergei Pond on 12/29/22.
//

#include <matrix.h>

#ifdef _SLKP_USE_APPLE_BLAS
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
  AtlasConj=114};
extern "C" void cblas_dgemm(const enum CBLAS_ORDER __Order,
                 const enum CBLAS_TRANSPOSE __TransA,
                 const enum CBLAS_TRANSPOSE __TransB, const int __M, const int __N,
                 const int __K, const double __alpha, const double *__A,
                 const int __lda, const double *__B, const int __ldb,
                 const double __beta, double *__C, const int __ldc);
#endif

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------


#ifdef _SLKP_USE_ARM_NEON
void _hy_matrix_multiply_4x4 (double * C, double *A, double *B, int stride, bool add) {
    
    float64x2x2_t A1, A2, A3, A4, A5, A6, A7, A8;
    float64x2x2_t B1;
    float64x2x2_t C1, C2, C3, C4;
    
    auto handle_block = [&] ()-> void {
        C1.val[0] = vmulq_f64 (A1.val[0],B1.val[0]);
        // C00 = A00 * B00; C01 = A00 * B01
        C1.val[1] = vmulq_f64 (A1.val[0],B1.val[1]);
        // C02 = A00 * B02; C01 = A00 * B03
        C2.val[0] = vmulq_f64 (A2.val[0],B1.val[0]);
        // C10 = A10 * B00; C01 = A00 * B01
        C2.val[1] = vmulq_f64 (A2.val[0],B1.val[1]);
        // C12 = A10 * B02; C13 = A10 * B03
        C3.val[0] = vmulq_f64 (A3.val[0],B1.val[0]);
        C3.val[1] = vmulq_f64 (A3.val[0],B1.val[1]);
        C4.val[0] = vmulq_f64 (A4.val[0],B1.val[0]);
        // C30 = A30 * B00; C30
        C4.val[1] = vmulq_f64 (A4.val[0],B1.val[1]);
    };
    
    auto handle_block2 = [&] (float64x2_t a1, float64x2_t a2, float64x2_t a3, float64x2_t a4) -> void {
        C1.val[0] = vfmaq_f64 (C1.val[0],a1,B1.val[0]);
        // C00 = A00 * B00; C01 = A00 * B01
        C2.val[0] = vfmaq_f64 (C2.val[0],a2,B1.val[0]);
        C3.val[0] = vfmaq_f64 (C3.val[0],a3,B1.val[0]);
        C4.val[0] = vfmaq_f64 (C4.val[0],a4,B1.val[0]);
        C1.val[1] = vfmaq_f64 (C1.val[1],a1,B1.val[1]);
        // C02 = A00 * B02; C01 = A00 * B03
        // C10 = A10 * B00; C01 = A00 * B01
        C2.val[1] = vfmaq_f64 (C2.val[1],a2,B1.val[1]);
        // C12 = A10 * B02; C13 = A10 * B03
        C3.val[1] = vfmaq_f64 (C3.val[1],a3,B1.val[1]);
        // C30 = A30 * B00; C30
        C4.val[1] = vfmaq_f64 (C4.val[1],a4,B1.val[1]);
    };
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;
    

    A1 = vld2q_dup_f64 (A);
    A2 = vld2q_dup_f64 (A+S1);
    A3 = vld2q_dup_f64 (A+S2);
    A4 = vld2q_dup_f64 (A+S3);

    A5 = vld2q_dup_f64 (A+2);
    A6 = vld2q_dup_f64 (A+S1+2);
    A7 = vld2q_dup_f64 (A+S2+2);
    A8 = vld2q_dup_f64 (A+S3+2);

    B1 = vld1q_f64_x2 (B);
    if (add) {
        C1 = vld1q_f64_x2 (C);
        C2 = vld1q_f64_x2 (C+S1);
        C3 = vld1q_f64_x2 (C+S2);
        C4 = vld1q_f64_x2 (C+S3);
        handle_block2 (A1.val[0],A2.val[0],A3.val[0],A4.val[0]);
    } else {
        handle_block ();
    }
    
    B1 = vld1q_f64_x2 (B+S1);
    handle_block2 (A1.val[1],A2.val[1],A3.val[1],A4.val[1]);

    B1 = vld1q_f64_x2 (B+S2);
    handle_block2 (A5.val[0],A6.val[0],A7.val[0],A8.val[0]);

    B1 = vld1q_f64_x2 (B+S3);
    handle_block2 (A5.val[1],A6.val[1],A7.val[1],A8.val[1]);

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
    vst1q_f64_x2 (C+S2, C3);
    vst1q_f64_x2 (C+S3, C4);
}



/**
        Nx1 cases
*/


void _hy_matrix_multiply_4x1x4 (double *C, double* A, double *B, int stride) {

    float64x2_t   A1,A2,A3,A4;
    float64x2x2_t B1;
    float64x2x2_t C1,C2,C3,C4;
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;
    
    A1  = vdupq_n_f64 (A[0]);
    A2  = vdupq_n_f64 (A[S1]);
    A3  = vdupq_n_f64 (A[S2]);
    A4  = vdupq_n_f64 (A[S3]);
    
    B1 = vld1q_f64_x2 (B);
    
    C1 = vld1q_f64_x2 (C);
    C2 = vld1q_f64_x2 (C + S1);
    C3 = vld1q_f64_x2 (C + S2);
    C4 = vld1q_f64_x2 (C + S3);

    C1.val[0] = vfmaq_f64 (C1.val[0],A1,B1.val[0]);
    C2.val[0] = vfmaq_f64 (C2.val[0],A2,B1.val[0]);
    C3.val[0] = vfmaq_f64 (C3.val[0],A3,B1.val[0]);
    C4.val[0] = vfmaq_f64 (C4.val[0],A4,B1.val[0]);
 
    C1.val[1] = vfmaq_f64 (C1.val[1],A1,B1.val[1]);
    C2.val[1] = vfmaq_f64 (C2.val[1],A2,B1.val[1]);
    C3.val[1] = vfmaq_f64 (C3.val[1],A3,B1.val[1]);
    C4.val[1] = vfmaq_f64 (C4.val[1],A4,B1.val[1]);
    
    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
    vst1q_f64_x2 (C+S2, C3);
    vst1q_f64_x2 (C+S3, C4);
}

void _hy_matrix_multiply_4x4x1 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;
    
    float64x2x2_t   A1,A2;
    float64x2x2_t   C1,C2;
    float64x2_t     B1,B2;
    
    A1.val[0] = vsetq_lane_f64 (A[0],A1.val[0],0);
    A1.val[0] = vsetq_lane_f64 (A[S1],A1.val[0],1);
    A1.val[1] = vsetq_lane_f64 (A[S2],A1.val[1],0);
    A1.val[1] = vsetq_lane_f64 (A[S3],A1.val[1],1);
    
    A2.val[0] = vsetq_lane_f64 (A[1],A2.val[0],0);
    A2.val[0] = vsetq_lane_f64 (A[S1+1],A2.val[0],1);
    A2.val[1] = vsetq_lane_f64 (A[S2+1],A2.val[1],0);
    A2.val[1] = vsetq_lane_f64 (A[S3+1],A2.val[1],1);

    B1 = vdupq_n_f64 (B[0]);
    B2 = vdupq_n_f64 (B[S1]);
    
    C1.val[0] = vmulq_f64(A1.val[0], B1);
    C1.val[1] = vmulq_f64(A1.val[1], B1);
    C2.val[0] = vmulq_f64(A2.val[0], B2);
    C2.val[1] = vmulq_f64(A2.val[1], B2);

    B1 = vdupq_n_f64 (B[S2]);
    B2 = vdupq_n_f64 (B[S3]);
    
    A1.val[0] = vsetq_lane_f64 (A[2],A1.val[0],0);
    A1.val[0] = vsetq_lane_f64 (A[S1+2],A1.val[0],1);
    A1.val[1] = vsetq_lane_f64 (A[S2+2],A1.val[1],0);
    A1.val[1] = vsetq_lane_f64 (A[S3+2],A1.val[1],1);
    
    A2.val[0] = vsetq_lane_f64 (A[3],A2.val[0],0);
    A2.val[0] = vsetq_lane_f64 (A[S1+3],A2.val[0],1);
    A2.val[1] = vsetq_lane_f64 (A[S2+3],A2.val[1],0);
    A2.val[1] = vsetq_lane_f64 (A[S3+3],A2.val[1],1);

    C1.val[0] = vfmaq_f64(C1.val[0], A1.val[0], B1);
    C1.val[1] = vfmaq_f64(C1.val[1], A1.val[1], B1);

    C2.val[0] = vfmaq_f64(C2.val[0], A2.val[0], B2);
    C2.val[1] = vfmaq_f64(C2.val[1], A2.val[1], B2);
    
    C1.val[0] = vaddq_f64 (C1.val[0],C2.val[0]);
    C1.val[1] = vaddq_f64 (C1.val[1],C2.val[1]);

    if (add) {
        C[0]  += vgetq_lane_f64 (C1.val[0],0);
        C[S1] += vgetq_lane_f64 (C1.val[0],1);
        C[S2] += vgetq_lane_f64 (C1.val[1],0);
        C[S3] += vgetq_lane_f64 (C1.val[1],1);
    } else {
        C[0]  = vgetq_lane_f64 (C1.val[0],0);
        C[S1] = vgetq_lane_f64 (C1.val[0],1);
        C[S2] = vgetq_lane_f64 (C1.val[1],0);
        C[S3] = vgetq_lane_f64 (C1.val[1],1);
    }
}


void _hy_matrix_multiply_1x4x4 (double *C, double* A, double *B, int stride, bool add) {
   
    int S1 = stride,
        S2 = stride << 1,
        S3 = stride + S2;

    float64x2_t     A1 = vdupq_n_f64 (A[0]),
                    A2 = vdupq_n_f64 (A[1]),
                    A3 = vdupq_n_f64 (A[2]),
                    A4 = vdupq_n_f64 (A[3]);
    
    float64x2x2_t   B1 = vld1q_f64_x2 (B),
                    B2 = vld1q_f64_x2 (B + S1),
                    B3 = vld1q_f64_x2 (B + S2),
                    B4 = vld1q_f64_x2 (B + S3);

    float64x2x2_t   C1,C2;
     
    C2.val[0] = vmulq_f64 (A1,B1.val[0]);
    C2.val[1] = vmulq_f64 (A1,B1.val[1]);
    if (add) {
        C1 = vld1q_f64_x2 (C);
        C1.val[0] = vfmaq_f64 (C1.val[0],A2,B2.val[0]);
        C1.val[1] = vfmaq_f64 (C1.val[1],A2,B2.val[1]);
    } else {
        C1.val[0] = vmulq_f64 (A2,B2.val[0]);
        C1.val[1] = vmulq_f64 (A2,B2.val[1]);
    }
    
    C2.val[0] = vfmaq_f64 (C2.val[0],A3,B3.val[0]);
    C2.val[1] = vfmaq_f64 (C2.val[1],A3,B3.val[1]);
    C1.val[0] = vfmaq_f64 (C1.val[0],A4,B4.val[0]);
    C1.val[1] = vfmaq_f64 (C1.val[1],A4,B4.val[1]);
    
    C1.val[0] = vaddq_f64 (C1.val[0],C2.val[0]);
    C1.val[1] = vaddq_f64 (C1.val[1],C2.val[1]);

    vst1q_f64_x2 (C, C1);
}

void _hy_matrix_multiply_1x4x1 (double *C, double* A, double *B, int stride, bool add) {
    
    double s1  = (A[0] * B[0]),
           s2  = A[1] * B[stride],
           s3  = A[2] * B [stride << 1],
           s4  = A[3] * B [(stride << 1) + stride];
    s1 += s3; s2 += s4;
    if (add) {
        C[0] += s1+s2;
    } else {
        C[0] = s1+s2;
    }
}


void _hy_matrix_multiply_4x1x1 (double *C, double* A, double *B, int stride) {

    int S1 = stride,
        S2 = stride << 1,
        S3 = stride + S2;

    float64x2x2_t   A1,
                    C1;

    float64x2_t     B1 = vdupq_n_f64 (B[0]);
    
    A1.val[0] = vsetq_lane_f64 (A[0],A1.val[0],0);
    A1.val[0] = vsetq_lane_f64 (A[S1],A1.val[0],1);
    A1.val[1] = vsetq_lane_f64 (A[S2],A1.val[1],0);
    A1.val[1] = vsetq_lane_f64 (A[S3],A1.val[1],1);

    C1.val[0] = vsetq_lane_f64 (C[0],C1.val[0],0);
    C1.val[0] = vsetq_lane_f64 (C[S1],C1.val[0],1);
    C1.val[1] = vsetq_lane_f64 (C[S2],C1.val[1],0);
    C1.val[1] = vsetq_lane_f64 (C[S3], C1.val[1],1);
 
    C1.val[0] = vfmaq_f64 (C1.val[0],A1.val[0],B1);
    C1.val[1] = vfmaq_f64 (C1.val[1],A1.val[1],B1);
    
    vst1q_lane_f64 (C, C1.val[0],0);
    vst1q_lane_f64 (C+S1, C1.val[0],1);
    vst1q_lane_f64 (C+S2, C1.val[1],0);
    vst1q_lane_f64 (C+S3, C1.val[1],1);
}

void _hy_matrix_multiply_1x1x4 (double *C, double* A, double *B, int stride) {
    float64x2_t     A1 = vdupq_n_f64 (A[0]);
    vst1q_f64 (C, vfmaq_f64 (vld1q_f64 (C), A1, vld1q_f64 (B)));
    vst1q_f64 (C+2, vfmaq_f64 (vld1q_f64 (C+2), A1, vld1q_f64 (B+2)));
 }

/**
        Nx2 cases
 */
void _hy_matrix_multiply_4x2x4 (double *C, double* A, double *B, int stride) {
    
    float64x2_t   A1,A2,A3,A4;
    float64x2x2_t B1;
    float64x2x2_t C1,C2,C3,C4;
    
    auto mult_block = [&]() -> void {
        C1.val[0] = vfmaq_f64 (C1.val[0],A1,B1.val[0]);
        C1.val[1] = vfmaq_f64 (C1.val[1],A1,B1.val[1]);
        C2.val[0] = vfmaq_f64 (C2.val[0],A2,B1.val[0]);
        C2.val[1] = vfmaq_f64 (C2.val[1],A2,B1.val[1]);
        C3.val[0] = vfmaq_f64 (C3.val[0],A3,B1.val[0]);
        C3.val[1] = vfmaq_f64 (C3.val[1],A3,B1.val[1]);
        C4.val[0] = vfmaq_f64 (C4.val[0],A4,B1.val[0]);
        C4.val[1] = vfmaq_f64 (C4.val[1],A4,B1.val[1]);
    };

    int S1 = stride, S2 = stride << 1, S3 = stride + S2;
    
    A1  = vdupq_n_f64 (A[0]);
    A2  = vdupq_n_f64 (A[S1]);
    A3  = vdupq_n_f64 (A[S2]);
    A4  = vdupq_n_f64 (A[S3]);
    
    B1 = vld1q_f64_x2 (B);
    
    C1 = vld1q_f64_x2 (C);
    C2 = vld1q_f64_x2 (C + S1);
    C3 = vld1q_f64_x2 (C + S2);
    C4 = vld1q_f64_x2 (C + S3);

    mult_block ();
    
    A1  = vdupq_n_f64 (A[1]);
    A2  = vdupq_n_f64 (A[S1+1]);
    A3  = vdupq_n_f64 (A[S2+1]);
    A4  = vdupq_n_f64 (A[S3+1]);
    B1 = vld1q_f64_x2 (B+S1);
    
    mult_block ();

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
    vst1q_f64_x2 (C+S2, C3);
    vst1q_f64_x2 (C+S3, C4);
}

void _hy_matrix_multiply_4x4x2 (double *C, double* A, double *B, int stride, bool add) {

    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t A1,
                  A2,
                  A3,
                  A4;
    
    float64x2_t B1 = vld1q_f64 (B),
                B2 = vld1q_f64 (B+S1);

    float64x2_t C1, C2, C3, C4;
    
    if (add) {
        C1 = vld1q_f64 (C),
        C2 = vld1q_f64 (C + S1),
        C3 = vld1q_f64 (C + S2),
        C4 = vld1q_f64 (C + S3);
    }
    
    A1 = vld2q_dup_f64 (A);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2); // (A[2][0]x2, A[2][1]x2)
    A4 = vld2q_dup_f64 (A+S3); // (A[3][0]x2, A[3][1]x2)
 
    if (add) {
        C1 = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2 = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3 = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C4 = vfmaq_f64 (C4,A4.val[0],B1); // 30 += 30*00, 31 += 30*01
    } else {
        C1 = vmulq_f64 (A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2 = vmulq_f64 (A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3 = vmulq_f64 (A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C4 = vmulq_f64 (A4.val[0],B1); // 30 += 30*00, 31 += 30*01
    }

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    C4 = vfmaq_f64 (C4,A4.val[1],B2); // 30 += 31*10, 31 += 31*11

    B1 = vld1q_f64 (B+S2),
    B2 = vld1q_f64 (B+S3);
    A1 = vld2q_dup_f64 (A+2);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1+2); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2+2); // (A[2][0]x2, A[2][1]x2)
    A4 = vld2q_dup_f64 (A+S3+2); // (A[3][0]x2, A[3][1]x2)

    C1 = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
    C2 = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
    C3 = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
    C4 = vfmaq_f64 (C4,A4.val[0],B1); // 30 += 30*00, 31 += 30*01

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    C4 = vfmaq_f64 (C4,A4.val[1],B2); // 30 += 31*10, 31 += 31*11

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
    vst1q_f64 (C+S2, C3);
    vst1q_f64 (C+S3, C4);
}

void _hy_matrix_multiply_2x4x2 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t   A1 = vld2q_dup_f64 (A),     // A00x2 A01x2
                    A2 = vld2q_dup_f64 (A+S1);   // A02x2 A03x2

                    
    float64x2_t     C1, C2, B1, B2;

    
    B1 = vld1q_f64 (B);
    B2 = vld1q_f64 (B+S1);

    
    if (add) {
        C1 = vld1q_f64 (C);
        C2 = vld1q_f64 (C+S1);
        C1 = vfmaq_f64 (C1, A1.val[0], B1); // 00 += 00 * 00 ; 01 += 00*01
        C2 = vfmaq_f64 (C2, A2.val[0], B1); // 10 += 10 * 00 ; 11 += 10*01
    } else {
        C1 = vmulq_f64 (A1.val[0], B1);
        C2 = vmulq_f64 (A2.val[0], B1);
    }

    C1 = vfmaq_f64 (C1, A1.val[1], B2); // 00 += 02*20, 01 += 02*21
    C2 = vfmaq_f64 (C2, A2.val[1], B2); // 10 += 12*20, 11 += 12*21

    B1 = vld1q_f64 (B+S2);
    B2 = vld1q_f64 (B+S3);
    A1 = vld2q_dup_f64 (A+2);
    A2 = vld2q_dup_f64 (A+S1+2);

    C1 = vfmaq_f64 (C1, A1.val[0], B1);
    C2 = vfmaq_f64 (C2, A2.val[0], B1);
    C1 = vfmaq_f64 (C1, A1.val[1], B2);
    C2 = vfmaq_f64 (C2, A2.val[1], B2);

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);


}

void _hy_matrix_multiply_2x4x4 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t   A1 = vld2q_dup_f64 (A),     // A00x2 A01x2
                    A2 = vld2q_dup_f64 (A+S1);   // A02x2 A03x2

                    
    float64x2x2_t   C1, C2,
                    B1, B2;

    
    B1 = vld1q_f64_x2 (B);
    B2 = vld1q_f64_x2 (B+S1);

    
    if (add) {
        C1 = vld1q_f64_x2 (C);
        C2 = vld1q_f64_x2 (C+S1);
        C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[0], B1.val[0]); // 00 += 00 * 00 ; 01 += 00*01
        C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[0], B1.val[1]); // 02 += 00 * 20 ; 03 += 00*30
        C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[0], B1.val[0]); // 10 += 10 * 00;  11 += 10*01
        C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[0], B1.val[1]); // 12 += 10 * 02;  13 += 10*03
    } else {
        C1.val[0] = vmulq_f64 (A1.val[0], B1.val[0]);
        C1.val[1] = vmulq_f64 (A1.val[0], B1.val[1]);
        C2.val[0] = vmulq_f64 (A2.val[0], B1.val[0]);
        C2.val[1] = vmulq_f64 (A2.val[0], B1.val[1]);
    }

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[1], B2.val[0]); // 00 += 10 * 10 +
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[1], B2.val[1]); //
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[1], B2.val[0]); //
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[1], B2.val[1]); //

    B1 = vld1q_f64_x2 (B+S2);
    B2 = vld1q_f64_x2 (B+S3);
    
    A1 = vld2q_dup_f64 (A+2);
    A2 = vld2q_dup_f64 (A+S1+2);

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[0], B1.val[0]); // 00 += 00 * 00 ; 01 += 00*01
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[0], B1.val[1]); // 02 += 00 * 20 ; 03 += 00*30
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[0], B1.val[0]); // 10 += 10 * 00;  11 += 10*01
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[0], B1.val[1]); // 12 += 10 * 02;  13 += 10*03

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[1], B2.val[0]); // 00 += 10 * 10 +
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[1], B2.val[1]); //
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[1], B2.val[0]); //
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[1], B2.val[1]); //

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);


}

void _hy_matrix_multiply_2x2x4 (double *C, double* A, double *B, int stride) {
    
    int S1 = stride;

    float64x2x2_t       A1 = vld2q_dup_f64 (A),      // A00x2 A01x2
                        A2 = vld2q_dup_f64 (A+S1);   // A10x2 A11x2

                    
    float64x2x2_t       C1, C2;
    float64x2x2_t       B1, B2;

    B1 = vld1q_f64_x2 (B);
    B2 = vld1q_f64_x2 (B+S1);
 
    C1 = vld1q_f64_x2 (C);
    C2 = vld1q_f64_x2 (C+S1);
    
    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[0], B1.val[0]); // 00 += 00 * 00 ; 01 += 00*01
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[0], B1.val[0]); // 10 += 10 * 00 ; 11 += 00*01
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[0], B1.val[1]); // 02 += 00 * 02 ; 03 += 00*03
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[0], B1.val[1]); // 12 += 10 * 02 ; 13 += 10*03

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[1], B2.val[0]); // 00 += 01 * 10 ; 01 += 00*01
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[1], B2.val[0]);
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[1], B2.val[1]);
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[1], B2.val[1]);

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
}



void _hy_matrix_multiply_4x2x2 (double *C, double* A, double *B, int stride) {
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t       A1 = vld2q_dup_f64 (A),      // A00x2 A01x2
                        A2 = vld2q_dup_f64 (A+S1),   // A10x2 A11x2
                        A3 = vld2q_dup_f64 (A+S2),   // A20x2 A21x2
                        A4 = vld2q_dup_f64 (A+S3);   // A30x2 A31x2

                    
    float64x2_t       C1, C2, C3, C4;
    float64x2_t       B1, B2;

    B1 = vld1q_f64 (B);
    B2 = vld1q_f64 (B+S1);
 
    C1 = vld1q_f64 (C);
    C2 = vld1q_f64 (C+S1);
    C3 = vld1q_f64 (C+S2);
    C4 = vld1q_f64 (C+S3);

    C1 = vfmaq_f64 (C1, A1.val[0], B1); // 00 += 00 * 00 ; 01 += 00*01
    C2 = vfmaq_f64 (C2, A2.val[0], B1); // 10 += 10 * 00 ; 11 += 10*01
    C3 = vfmaq_f64 (C3, A3.val[0], B1); // 20 += 20 * 00 ; 21 += 20*01
    C4 = vfmaq_f64 (C4, A4.val[0], B1); // 30 += 30 * 00 ; 31 += 30*01

    C1 = vfmaq_f64 (C1, A1.val[1], B2); // 00 += 01 * 10 ; 01 += 01*11
    C2 = vfmaq_f64 (C2, A2.val[1], B2); // 10 += 10 * 00 ; 11 += 10*01
    C3 = vfmaq_f64 (C3, A3.val[1], B2); // 20 += 20 * 00 ; 21 += 20*01
    C4 = vfmaq_f64 (C4, A4.val[1], B2); // 30 += 30 * 00 ; 31 += 30*01

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
    vst1q_f64 (C+S2, C3);
    vst1q_f64 (C+S3, C4);
}

void _hy_matrix_multiply_2x2x2 (double *C, double* A, double *B, int stride, bool add) {
   
    int S1 = stride;

    float64x2x2_t       A1 = vld2q_dup_f64 (A),      // A00x2 A01x2
                        A2 = vld2q_dup_f64 (A+S1);   // A10x2 A11x2

    float64x2_t       C1, C2;
    float64x2_t       B1, B2;

    B1 = vld1q_f64 (B);
    B2 = vld1q_f64 (B+S1);
 
    if (add) {
        C1 = vld1q_f64 (C);
        C2 = vld1q_f64 (C+S1);
        
        C1 = vfmaq_f64 (C1, A1.val[0], B1); // 00 += 00 * 00 ; 01 += 00*01
        C2 = vfmaq_f64 (C2, A2.val[0], B1); // 10 += 10 * 00 ; 11 += 10*01
    } else {
        C1 = vmulq_f64 (A1.val[0], B1);
        C2 = vmulq_f64 (A2.val[0], B1);
    }

    C1 = vfmaq_f64 (C1, A1.val[1], B2); // 00 += 01 * 10 ; 01 += 01*11
    C2 = vfmaq_f64 (C2, A2.val[1], B2); // 10 += 10 * 00 ; 11 += 10*01

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
}

/**
 Nx3 cases
*/

void _hy_matrix_multiply_4x3x4 (double *C, double* A, double *B, int stride) {
    
    float64x2x2_t A1,A2,A3,A4;
    float64x2x2_t B1;
    float64x2x2_t C1,C2,C3,C4;
    
    auto mult_block = [&](int i) -> void {
        C1.val[0] = vfmaq_f64 (C1.val[0],A1.val[i],B1.val[0]);
        C2.val[0] = vfmaq_f64 (C2.val[0],A2.val[i],B1.val[0]);
        C3.val[0] = vfmaq_f64 (C3.val[0],A3.val[i],B1.val[0]);
        C4.val[0] = vfmaq_f64 (C4.val[0],A4.val[i],B1.val[0]);
        C1.val[1] = vfmaq_f64 (C1.val[1],A1.val[i],B1.val[1]);
        C2.val[1] = vfmaq_f64 (C2.val[1],A2.val[i],B1.val[1]);
        C3.val[1] = vfmaq_f64 (C3.val[1],A3.val[i],B1.val[1]);
        C4.val[1] = vfmaq_f64 (C4.val[1],A4.val[i],B1.val[1]);
    };

    int S1 = stride, S2 = stride << 1, S3 = stride + S2;
    
    A1  = vld2q_dup_f64 (A);
    A2  = vld2q_dup_f64 (A+S1);
    A3  = vld2q_dup_f64 (A+S2);
    A4  = vld2q_dup_f64 (A+S3);
    
    B1 = vld1q_f64_x2 (B);
    
    C1 = vld1q_f64_x2 (C);
    C2 = vld1q_f64_x2 (C + S1);
    C3 = vld1q_f64_x2 (C + S2);
    C4 = vld1q_f64_x2 (C + S3);

    mult_block (0);
    
    B1 = vld1q_f64_x2 (B+S1);
    
    mult_block (1);

    A1.val[0]  = vdupq_n_f64 (A[2]);
    A2.val[0]  = vdupq_n_f64 (A[S1+2]);
    A3.val[0]  = vdupq_n_f64 (A[S2+2]);
    A4.val[0]  = vdupq_n_f64 (A[S3+2]);
    B1 = vld1q_f64_x2 (B+S2);
    
    mult_block (0);

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
    vst1q_f64_x2 (C+S2, C3);
    vst1q_f64_x2 (C+S3, C4);
}

void _hy_matrix_multiply_4x4x3 (double *C, double* A, double *B, int stride, bool add) {
 
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t A1,
                  A2,
                  A3,
                  A4;
    
    float64x2_t B1 = vld1q_f64 (B),
                B2 = vld1q_f64 (B+S1),
                B_LC; // last col of B

    float64x2_t C1,
                C2,
                C3,
                C4;
    
    float64x2x2_t C_LR, A_C;
    
    if (add) {
        C1 = vld1q_f64 (C),
        C2 = vld1q_f64 (C + S1),
        C3 = vld1q_f64 (C + S2),
        C4 = vld1q_f64 (C + S3);
        C_LR.val[0] = vsetq_lane_f64 (C[2],C_LR.val[0],0);
        C_LR.val[0] = vsetq_lane_f64 (C[S1+2],C_LR.val[0],1);
        C_LR.val[1] = vsetq_lane_f64 (C[S2+2],C_LR.val[1],0);
        C_LR.val[1] = vsetq_lane_f64 (C[S3+2],C_LR.val[1],1);

    }
    
    A1 = vld2q_dup_f64 (A);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2); // (A[2][0]x2, A[2][1]x2)
    A4 = vld2q_dup_f64 (A+S3); // (A[3][0]x2, A[3][1]x2)
 
    B_LC = vdupq_n_f64 (B[2]);
    A_C.val[0]  = vsetq_lane_f64 (A[0],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3],A_C.val[1],1);

    if (add) {
        C1   = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2   = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3   = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C4   = vfmaq_f64 (C4,A4.val[0],B1); // 30 += 30*00, 31 += 30*01
        
        C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
        C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);
    } else {
        C1 = vmulq_f64 (A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2 = vmulq_f64 (A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3 = vmulq_f64 (A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C4 = vmulq_f64 (A4.val[0],B1); // 30 += 30*00, 31 += 30*01
        C_LR.val[0] = vmulq_f64 (A_C.val[0],B_LC);
        C_LR.val[1] = vmulq_f64 (A_C.val[1],B_LC);
    }

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    C4 = vfmaq_f64 (C4,A4.val[1],B2); // 30 += 31*10, 31 += 31*11
    B_LC = vdupq_n_f64 (B[S1+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[1],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+1],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2+1],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3+1],A_C.val[1],1);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);


    B1 = vld1q_f64 (B+S2),
    B2 = vld1q_f64 (B+S3);
    A1 = vld2q_dup_f64 (A+2);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1+2); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2+2); // (A[2][0]x2, A[2][1]x2)
    A4 = vld2q_dup_f64 (A+S3+2); // (A[3][0]x2, A[3][1]x2)

    C1 = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
    C2 = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
    C3 = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
    C4 = vfmaq_f64 (C4,A4.val[0],B1); // 30 += 30*00, 31 += 30*01
    B_LC = vdupq_n_f64 (B[S2+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[2],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+2],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2+2],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3+2],A_C.val[1],1);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    C4 = vfmaq_f64 (C4,A4.val[1],B2); // 30 += 31*10, 31 += 31*11
    B_LC = vdupq_n_f64 (B[S3+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[3],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+3],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2+3],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3+3],A_C.val[1],1);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
    vst1q_f64 (C+S2, C3);
    vst1q_f64 (C+S3, C4);
    
    vst1q_lane_f64 (C+2,    C_LR.val[0], 0);
    vst1q_lane_f64 (C+S1+2, C_LR.val[0], 1);
    vst1q_lane_f64 (C+S2+2, C_LR.val[1], 0);
    vst1q_lane_f64 (C+S3+2, C_LR.val[1], 1);
}

void _hy_matrix_multiply_4x3x3 (double *C, double* A, double *B, int stride, bool add) {
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t A1,
                  A2,
                  A3,
                  A4;
    
    float64x2_t B1 = vld1q_f64 (B),
                B2 = vld1q_f64 (B+S1),
                B_LC; // last col of B

    float64x2_t C1,
                C2,
                C3,
                C4;
    
    float64x2x2_t C_LR, A_C;
    
    if (add) {
        C1 = vld1q_f64 (C),
        C2 = vld1q_f64 (C + S1),
        C3 = vld1q_f64 (C + S2),
        C4 = vld1q_f64 (C + S3);
        C_LR.val[0] = vsetq_lane_f64 (C[2],C_LR.val[0],0);
        C_LR.val[0] = vsetq_lane_f64 (C[S1+2],C_LR.val[0],1);
        C_LR.val[1] = vsetq_lane_f64 (C[S2+2],C_LR.val[1],0);
        C_LR.val[1] = vsetq_lane_f64 (C[S3+2],C_LR.val[1],1);
   }
    
    A1 = vld2q_dup_f64 (A);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2); // (A[2][0]x2, A[2][1]x2)
    A4 = vld2q_dup_f64 (A+S3); // (A[3][0]x2, A[3][1]x2)
 
    B_LC = vdupq_n_f64 (B[2]);
    A_C.val[0]  = vsetq_lane_f64 (A[0],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3],A_C.val[1],1);

    if (add) {
        C1   = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2   = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3   = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C4   = vfmaq_f64 (C4,A4.val[0],B1); // 30 += 30*00, 31 += 30*01
        C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
        C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);
    } else {
        C1 = vmulq_f64 (A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2 = vmulq_f64 (A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3 = vmulq_f64 (A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C4 = vmulq_f64 (A4.val[0],B1); // 30 += 30*00, 31 += 30*01
        C_LR.val[0] = vmulq_f64 (A_C.val[0],B_LC);
        C_LR.val[1] = vmulq_f64 (A_C.val[1],B_LC);
    }

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    C4 = vfmaq_f64 (C4,A4.val[1],B2); // 30 += 31*10, 31 += 31*11
    
    B_LC = vdupq_n_f64 (B[S1+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[1],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+1],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2+1],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3+1],A_C.val[1],1);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);


    B1 = vld1q_f64 (B+S2),
    A1 = vld2q_dup_f64 (A+2);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1+2); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2+2); // (A[2][0]x2, A[2][1]x2)
    A4 = vld2q_dup_f64 (A+S3+2); // (A[3][0]x2, A[3][1]x2)

    C1 = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
    C2 = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
    C3 = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
    C4 = vfmaq_f64 (C4,A4.val[0],B1); // 30 += 30*00, 31 += 30*01
    
    B_LC = vdupq_n_f64 (B[S2+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[2],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+2],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2+2],A_C.val[1],0);
    A_C.val[1]  = vsetq_lane_f64 (A[S3+2],A_C.val[1],1);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
    vst1q_f64 (C+S2, C3);
    vst1q_f64 (C+S3, C4);
    
    vst1q_lane_f64 (C+2,    C_LR.val[0], 0);
    vst1q_lane_f64 (C+S1+2, C_LR.val[0], 1);
    vst1q_lane_f64 (C+S2+2, C_LR.val[1], 0);
    vst1q_lane_f64 (C+S3+2, C_LR.val[1], 1);
}

void _hy_matrix_multiply_3x3x4 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride, S2 = stride << 1;

    float64x2x2_t A1,
                  A2,
                  A3,
                  B1 = vld1q_f64_x2 (B);
    
  
    float64x2x2_t C1,
                  C2,
                  C3;
    
    float64x2_t   A4, A5, A6;
    
    
    if (add) {
        C1 = vld1q_f64_x2 (C),
        C2 = vld1q_f64_x2 (C + S1),
        C3 = vld1q_f64_x2 (C + S2);
    }
    
    A1 = vld2q_dup_f64 (A);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2); // (A[2][0]x2, A[2][1]x2)
    A4 = vdupq_n_f64 (A[2]);   // (A[0][2]x2)
    A5 = vdupq_n_f64 (A[S1+2]);   // (A[1][2]x2)
    A6 = vdupq_n_f64 (A[S2+2]);   // (A[2][2]x2)

    if (add) {
        C1.val[0]   = vfmaq_f64 (C1.val[0],A1.val[0],B1.val[0]); // 00 += 00*00, 01 += 00*01
        C2.val[0]   = vfmaq_f64 (C2.val[0],A2.val[0],B1.val[0]); // 10 += 10*00, 01 += 10*01
        C3.val[0]   = vfmaq_f64 (C3.val[0],A3.val[0],B1.val[0]); // 20 += 20*00, 21 += 20*01
        C1.val[1]   = vfmaq_f64 (C1.val[1],A1.val[0],B1.val[1]); // 02 += 00*02, 03 += 00*03
        C2.val[1]   = vfmaq_f64 (C2.val[1],A2.val[0],B1.val[1]); // 12 += 10*02, 13 += 10*03
        C3.val[1]   = vfmaq_f64 (C3.val[1],A3.val[0],B1.val[1]); // 22 += 20*02, 23 += 20*03
    } else {
        C1.val[0]   = vmulq_f64 (A1.val[0],B1.val[0]); // 00 += 00*00, 01 += 00*01
        C2.val[0]   = vmulq_f64 (A2.val[0],B1.val[0]); // 10 += 10*00, 01 += 10*01
        C3.val[0]   = vmulq_f64 (A3.val[0],B1.val[0]); // 20 += 20*00, 21 += 20*01
        C1.val[1]   = vmulq_f64 (A1.val[0],B1.val[1]); // 02 += 00*02, 03 += 00*03
        C2.val[1]   = vmulq_f64 (A2.val[0],B1.val[1]); // 12 += 10*02, 13 += 10*03
        C3.val[1]   = vmulq_f64 (A3.val[0],B1.val[1]); // 22 += 20*02, 23 += 20*03
    }
    
    B1 = vld1q_f64_x2 (B+S1);
    C1.val[0]   = vfmaq_f64 (C1.val[0],A1.val[1],B1.val[0]); // 00 += 01*10, 01 += 01*11
    C2.val[0]   = vfmaq_f64 (C2.val[0],A2.val[1],B1.val[0]); // 10 += 11*10, 11 += 11*11
    C3.val[0]   = vfmaq_f64 (C3.val[0],A3.val[1],B1.val[0]); // 20 += 21*10, 21 += 21*11
    C1.val[1]   = vfmaq_f64 (C1.val[1],A1.val[1],B1.val[1]); // 02 += 02*22, 03 += 02*23
    C2.val[1]   = vfmaq_f64 (C2.val[1],A2.val[1],B1.val[1]); // 12 += 12*22, 13 += 12*23
    C3.val[1]   = vfmaq_f64 (C3.val[1],A3.val[1],B1.val[1]); // 22 += 02*22, 23 += 22*23

    B1 = vld1q_f64_x2 (B+S2);
    C1.val[0]   = vfmaq_f64 (C1.val[0],A4,B1.val[0]); // 00 += 02*20, 01 += 02*21
    C2.val[0]   = vfmaq_f64 (C2.val[0],A5,B1.val[0]); // 10 += 12*20, 11 += 12*21
    C3.val[0]   = vfmaq_f64 (C3.val[0],A6,B1.val[0]); // 20 += 22*20, 21 += 22*21
    C1.val[1]   = vfmaq_f64 (C1.val[1],A4,B1.val[1]); // 02 += 02*22, 03 += 02*23
    C2.val[1]   = vfmaq_f64 (C2.val[1],A5,B1.val[1]); // 12 += 12*22, 13 += 12*23
    C3.val[1]   = vfmaq_f64 (C3.val[1],A6,B1.val[1]); // 22 += 22*22, 23 += 22*23

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
    vst1q_f64_x2 (C+S2, C3);
}

void _hy_matrix_multiply_3x4x4 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t   A1 = vld2q_dup_f64 (A),      // A00x2 A01x2
                    A2 = vld2q_dup_f64 (A+S1),   // A10x2 A11x2
                    A3 = vld2q_dup_f64 (A+S2);   // A20x2 A21x2
                    
    float64x2x2_t   C1, C2, C3,
                    B1, B2;

    
    B1 = vld1q_f64_x2 (B);
    B2 = vld1q_f64_x2 (B+S1);

    
    if (add) {
        C1 = vld1q_f64_x2 (C);
        C2 = vld1q_f64_x2 (C+S1);
        C3 = vld1q_f64_x2 (C+S2);
        C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[0], B1.val[0]); // 00 += 00 * 00 ; 01 += 00*01
        C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[0], B1.val[1]); // 02 += 00 * 20 ; 03 += 00*30
        C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[0], B1.val[0]); // 10 += 10 * 00;  11 += 10*01
        C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[0], B1.val[1]); // 12 += 10 * 02;  13 += 10*03
        C3.val[0] = vfmaq_f64 (C3.val[0], A3.val[0], B1.val[0]); // 10 += 10 * 00;  11 += 10*01
        C3.val[1] = vfmaq_f64 (C3.val[1], A3.val[0], B1.val[1]); // 12 += 10 * 02;  13 += 10*03
    } else {
        C1.val[0] = vmulq_f64 (A1.val[0], B1.val[0]);
        C1.val[1] = vmulq_f64 (A1.val[0], B1.val[1]);
        C2.val[0] = vmulq_f64 (A2.val[0], B1.val[0]);
        C2.val[1] = vmulq_f64 (A2.val[0], B1.val[1]);
        C3.val[0] = vmulq_f64 (A3.val[0], B1.val[0]);
        C3.val[1] = vmulq_f64 (A3.val[0], B1.val[1]);
    }

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[1], B2.val[0]); // 00 += 10 * 10 +
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[1], B2.val[1]); //
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[1], B2.val[0]); //
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[1], B2.val[1]); //
    C3.val[0] = vfmaq_f64 (C3.val[0], A3.val[1], B2.val[0]); //
    C3.val[1] = vfmaq_f64 (C3.val[1], A3.val[1], B2.val[1]); //

    B1 = vld1q_f64_x2 (B+S2);
    B2 = vld1q_f64_x2 (B+S3);
    
    A1 = vld2q_dup_f64 (A+2);
    A2 = vld2q_dup_f64 (A+S1+2);
    A3 = vld2q_dup_f64 (A+S2+2);

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[0], B1.val[0]); // 00 += 00 * 00 ; 01 += 00*01
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[0], B1.val[1]); // 02 += 00 * 20 ; 03 += 00*30
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[0], B1.val[0]); // 10 += 10 * 00;  11 += 10*01
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[0], B1.val[1]); // 12 += 10 * 02;  13 += 10*03
    C3.val[0] = vfmaq_f64 (C3.val[0], A3.val[0], B1.val[0]); // 10 += 10 * 00;  11 += 10*01
    C3.val[1] = vfmaq_f64 (C3.val[1], A3.val[0], B1.val[1]); // 12 += 10 * 02;  13 += 10*03

    C1.val[0] = vfmaq_f64 (C1.val[0], A1.val[1], B2.val[0]); // 00 += 10 * 10 +
    C1.val[1] = vfmaq_f64 (C1.val[1], A1.val[1], B2.val[1]); //
    C2.val[0] = vfmaq_f64 (C2.val[0], A2.val[1], B2.val[0]); //
    C2.val[1] = vfmaq_f64 (C2.val[1], A2.val[1], B2.val[1]); //
    C3.val[0] = vfmaq_f64 (C3.val[0], A3.val[1], B2.val[0]); //
    C3.val[1] = vfmaq_f64 (C3.val[1], A3.val[1], B2.val[1]); //

    vst1q_f64_x2 (C,    C1);
    vst1q_f64_x2 (C+S1, C2);
    vst1q_f64_x2 (C+S2, C3);
}

void _hy_matrix_multiply_3x3x3 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride, S2 = stride << 1;

    float64x2x2_t A1,
                  A2,
                  A3;
    
    float64x2_t B1 = vld1q_f64 (B),
                B2 = vld1q_f64 (B+S1),
                B_LC; // last col of B

    float64x2_t C1,
                C2,
                C3;
    
    float64x2x2_t C_LR, A_C;
    
    if (add) {
        C1 = vld1q_f64 (C),
        C2 = vld1q_f64 (C + S1),
        C3 = vld1q_f64 (C + S2),
        C_LR.val[0] = vsetq_lane_f64 (C[2],C_LR.val[0],0);
        C_LR.val[0] = vsetq_lane_f64 (C[S1+2],C_LR.val[0],1);
        C_LR.val[1] = vdupq_n_f64 (C[S2+2]);
    }
    
    A1 = vld2q_dup_f64 (A);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2); // (A[2][0]x2, A[2][1]x2)
 
    B_LC = vdupq_n_f64 (B[2]);
    A_C.val[0]  = vsetq_lane_f64 (A[0],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2],A_C.val[1],0);

    if (add) {
        C1   = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2   = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3   = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        
        C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
        C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);
    } else {
        C1 = vmulq_f64 (A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2 = vmulq_f64 (A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3 = vmulq_f64 (A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C_LR.val[0] = vmulq_f64 (A_C.val[0],B_LC);
        C_LR.val[1] = vmulq_f64 (A_C.val[1],B_LC);
    }

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    B_LC = vdupq_n_f64 (B[S1+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[1],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+1],A_C.val[0],1);
    A_C.val[1]  = vdupq_n_f64 (A[S2+1]);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    B1 = vld1q_f64 (B+S2),
    A1 = vld2q_dup_f64 (A+2);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1+2); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2+2); // (A[2][0]x2, A[2][1]x2)
 
    C1 = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
    C2 = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
    C3 = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
    B_LC = vdupq_n_f64 (B[S2+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[2],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+2],A_C.val[0],1);
    A_C.val[1]  = vdupq_n_f64 (A[S2+2]);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
    vst1q_f64 (C+S2, C3);
    
    vst1q_lane_f64 (C+2,    C_LR.val[0], 0);
    vst1q_lane_f64 (C+S1+2, C_LR.val[0], 1);
    vst1q_lane_f64 (C+S2+2, C_LR.val[1], 0);
}

void _hy_matrix_multiply_3x4x3 (double *C, double* A, double *B, int stride, bool add) {
    

    int S1 = stride, S2 = stride << 1, S3 = stride + S2;

    float64x2x2_t A1,
                  A2,
                  A3;
    
    float64x2_t B1 = vld1q_f64 (B),
                B2 = vld1q_f64 (B+S1),
                B_LC; // last col of B

    float64x2_t C1,
                C2,
                C3;
    
    float64x2x2_t C_LR, A_C;
    
    if (add) {
        C1 = vld1q_f64 (C),
        C2 = vld1q_f64 (C + S1),
        C3 = vld1q_f64 (C + S2),
        C_LR.val[0] = vsetq_lane_f64 (C[2],C_LR.val[0],0);
        C_LR.val[0] = vsetq_lane_f64 (C[S1+2],C_LR.val[0],1);
        C_LR.val[1] = vdupq_n_f64 (C[S2+2]);
    }
    
    A1 = vld2q_dup_f64 (A);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2); // (A[2][0]x2, A[2][1]x2)
 
    B_LC = vdupq_n_f64 (B[2]);
    A_C.val[0]  = vsetq_lane_f64 (A[0],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1],A_C.val[0],1);
    A_C.val[1]  = vsetq_lane_f64 (A[S2],A_C.val[1],0);

    if (add) {
        C1   = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2   = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3   = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        
        C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
        C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);
    } else {
        C1 = vmulq_f64 (A1.val[0],B1); // 00 += 00*00, 01 += 00*01
        C2 = vmulq_f64 (A2.val[0],B1); // 10 += 10*00, 11 += 10*01
        C3 = vmulq_f64 (A3.val[0],B1); // 20 += 20*00, 21 += 20*01
        C_LR.val[0] = vmulq_f64 (A_C.val[0],B_LC);
        C_LR.val[1] = vmulq_f64 (A_C.val[1],B_LC);
    }

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    B_LC = vdupq_n_f64 (B[S1+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[1],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+1],A_C.val[0],1);
    A_C.val[1]  = vdupq_n_f64 (A[S2+1]);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    B1 = vld1q_f64 (B+S2),
    B2 = vld1q_f64 (B+S3);
    A1 = vld2q_dup_f64 (A+2);    // (A[0][0]x2, A[0][1]x2)
    A2 = vld2q_dup_f64 (A+S1+2); // (A[1][0]x2, A[1][1]x2)
    A3 = vld2q_dup_f64 (A+S2+2); // (A[2][0]x2, A[2][1]x2)
 
    C1 = vfmaq_f64 (C1,A1.val[0],B1); // 00 += 00*00, 01 += 00*01
    C2 = vfmaq_f64 (C2,A2.val[0],B1); // 10 += 10*00, 11 += 10*01
    C3 = vfmaq_f64 (C3,A3.val[0],B1); // 20 += 20*00, 21 += 20*01
    B_LC = vdupq_n_f64 (B[S2+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[2],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+2],A_C.val[0],1);
    A_C.val[1]  = vdupq_n_f64 (A[S2+2]);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    C1 = vfmaq_f64 (C1,A1.val[1],B2); // 00 += 01*10, 01 += 01*11
    C2 = vfmaq_f64 (C2,A2.val[1],B2); // 10 += 11*10, 11 += 11*11
    C3 = vfmaq_f64 (C3,A3.val[1],B2); // 20 += 21*10, 21 += 21*11
    B_LC = vdupq_n_f64 (B[S3+2]);
    A_C.val[0]  = vsetq_lane_f64 (A[3],A_C.val[0],0);
    A_C.val[0]  = vsetq_lane_f64 (A[S1+3],A_C.val[0],1);
    A_C.val[1]  = vdupq_n_f64 (A[S2+3]);
    C_LR.val[0] = vfmaq_f64 (C_LR.val[0],A_C.val[0],B_LC);
    C_LR.val[1] = vfmaq_f64 (C_LR.val[1],A_C.val[1],B_LC);

    vst1q_f64 (C,    C1);
    vst1q_f64 (C+S1, C2);
    vst1q_f64 (C+S2, C3);
    
    vst1q_lane_f64 (C+2,    C_LR.val[0], 0);
    vst1q_lane_f64 (C+S1+2, C_LR.val[0], 1);
    vst1q_lane_f64 (C+S2+2, C_LR.val[1], 0);
}

#endif

#if defined _SLKP_USE_SSE_INTRINSICS
#include <smmintrin.h>

void _hy_matrix_multiply_4x4 (double * C, double *A, double *B, int stride, bool add) {
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22, C31, C32, C41, C42;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2

    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03

    auto handle_block_mult = [&] () -> void {
        C11 = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C12 = _mm_mul_pd (A1,B2); // 00*02, 00*03
        C21 = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C22 = _mm_mul_pd (A2,B2); // 10*02, 10*03
        C31 = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C32 = _mm_mul_pd (A3,B2); // 20*02, 20*03
        C41 = _mm_mul_pd (A4,B1); // 20*00, 20*01
        C42 = _mm_mul_pd (A4,B2); // 20*02, 20*03
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C32 = _mm_add_pd (C32,_mm_mul_pd (A3,B2)); // 20*02, 20*03
        C41 = _mm_add_pd (C41,_mm_mul_pd (A4,B1)); // 20*00, 20*01
        C42 = _mm_add_pd (C42,_mm_mul_pd (A4,B2)); // 20*02, 20*03
    };
        
    if (add) {
        C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
        C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
        C31 = _mm_loadu_pd (C+S2); C32 = _mm_loadu_pd (C+S2+2);
        C41 = _mm_loadu_pd (C+S3); C42 = _mm_loadu_pd (C+S3+2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+2);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B2 = _mm_loadu_pd    (B+S2+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+3);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+3);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+3);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+3);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S3);      // 10,11
    B2 = _mm_loadu_pd    (B+S3+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S2+2, C32);
    _mm_storeu_pd (C+S3,   C41);
    _mm_storeu_pd (C+S3+2, C42);
}



/**
        Nx1 cases
*/


void _hy_matrix_multiply_4x1x4 (double *C, double* A, double *B, int stride) {

    __m128d   A1,A2,A3,A4;
    __m128d   B1,B2;
    __m128d   C1,C2,C3,C4;
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;
    
    A1  = _mm_loaddup_pd (A);     // A[0,0]
    A2  = _mm_loaddup_pd (A+S1);  // A[1,0]
    A3  = _mm_loaddup_pd (A+S2);  // A[2,0]
    A4  = _mm_loaddup_pd (A+S3);  // A[3,0]
    
    B1 = _mm_loadu_pd (B);
    B2 = _mm_loadu_pd (B+2);
    
    C1 = _mm_add_pd (_mm_loadu_pd (C),_mm_mul_pd (A1,B1)); // 00 += 00*00 01 += 00*01
    C2 = _mm_add_pd (_mm_loadu_pd (C+2),_mm_mul_pd (A1,B2)); // 02 += 00*02 03 += 00*03
    C3 = _mm_add_pd (_mm_loadu_pd (C+S1),_mm_mul_pd (A2,B1)); // 10 += 10*00 11 += 10*01
    C4 = _mm_add_pd (_mm_loadu_pd (C+S1+2),_mm_mul_pd (A2,B2)); // 00 += 00*00 01 += 00*01
    
    _mm_storeu_pd (C, C1);
    _mm_storeu_pd (C+2, C2);
    _mm_storeu_pd (C+S1, C3);
    _mm_storeu_pd (C+S1+2, C4);
    
    C1 = _mm_add_pd (_mm_loadu_pd (C+S2),_mm_mul_pd (A3,B1)); // 00 += 00*00 01 += 00*01
    C2 = _mm_add_pd (_mm_loadu_pd (C+S2+2),_mm_mul_pd (A3,B2)); // 02 += 00*02 03 += 00*03
    C3 = _mm_add_pd (_mm_loadu_pd (C+S3),_mm_mul_pd (A4,B1)); // 10 += 10*00 11 += 10*01
    C4 = _mm_add_pd (_mm_loadu_pd (C+S3+2),_mm_mul_pd (A4,B2)); // 00 += 00*00 01 += 00*01

    _mm_storeu_pd (C+S2, C1);
    _mm_storeu_pd (C+S2+2, C2);
    _mm_storeu_pd (C+S3, C3);
    _mm_storeu_pd (C+S3+2, C4);
}


void _hy_matrix_multiply_4x4x1 (double *C, double* A, double *B, int stride, bool add) {
        int  S1 = stride,
             S2 = stride << 1,
             S3 = S2 + stride;
         
        __m128d   A11,A12,A21,A22;
        __m128d   C11,C12,C21,C22;
        __m128d   B1,B2;
        
        // 00*00 + 01*10 + 02*20 + 03*30
    
        B1 = _mm_loadh_pd (_mm_load_sd (B), B+S1);
        B2 = _mm_loadh_pd (_mm_load_sd (B+S2), B+S3);

        A11 = _mm_loadu_pd (A);
        A12 = _mm_loadu_pd (A+2);
        A21 = _mm_loadu_pd (A+S1);
        A22 = _mm_loadu_pd (A+S1+2);
        
        C11 = _mm_mul_pd   (A11,B1); //00*00, 01*10
        C12 = _mm_mul_pd   (A12,B2); //02*20, 03*30
        C21 = _mm_mul_pd   (A21,B1);
        C22 = _mm_mul_pd   (A22,B2);
    
        C11 = _mm_add_pd   (C11,C12); // 00*00+02*20,01*10+03*30
        C21 = _mm_add_pd   (C21,C22);
        C11 = _mm_hadd_pd  (C11,C21);
        
        if (add) {
            C[0]  += _mm_cvtsd_f64 (C11);
            C[S1] += _mm_cvtsd_f64(_mm_unpackhi_pd (C11,C11));
        } else {
            _mm_storel_pd (C, C11);
            _mm_storeh_pd (C+S1,C11);
        }
        
        A11 = _mm_loadu_pd (A+S2);
        A12 = _mm_loadu_pd (A+S2+2);
        A21 = _mm_loadu_pd (A+S3);
        A22 = _mm_loadu_pd (A+S3+2);
        
        C11 = _mm_mul_pd   (A11,B1);
        C12 = _mm_mul_pd   (A12,B2);
        C21 = _mm_mul_pd   (A21,B1);
        C22 = _mm_mul_pd   (A22,B2);
    
        C11 = _mm_add_pd   (C11,C12); // 00*00+02*20,01*10+03*30
        C21 = _mm_add_pd   (C21,C22);
        C11 = _mm_hadd_pd  (C11,C21);
        
        if (add) {
            C[S2]  += _mm_cvtsd_f64 (C11);
            C[S3] += _mm_cvtsd_f64(_mm_unpackhi_pd (C11,C11));
        } else {
            _mm_storel_pd (C+S2, C11);
            _mm_storeh_pd (C+S3,C11);
        }
    
        
}


void _hy_matrix_multiply_1x4x4 (double *C, double* A, double *B, int stride, bool add) {

    int S1 = stride,
        S2 = stride << 1,
        S3 = stride + S2;

    __m128d       A1 = _mm_load1_pd (A),
                  A2 = _mm_load1_pd (A+1),
                  A3 = _mm_load1_pd (A+2),
                  A4 = _mm_load1_pd (A+3);
    
    __m128d       B11 = _mm_loadu_pd (B),
                  B12 = _mm_loadu_pd (B + 2),
                  B21 = _mm_loadu_pd (B + S1),
                  B22 = _mm_loadu_pd (B + S1 + 2);
                  
    // 00 = 00 * 00 + 01 * 10 + 02 * 20 + 03 * 30
    // 01 = 00 * 01 + 01 * 11 + 02 * 21 + 03 * 31
    // 02 = 00 * 02 + 01 * 12 + 02 * 22 + 03 * 32
    // 03 = 00 * 03 + 01 * 13 + 02 * 23 + 03 * 33

    __m128d       C11,C12;
     
    if (add) {
        C11 = _mm_add_pd (_mm_loadu_pd(C), _mm_mul_pd (A1,B11));
        C12 = _mm_add_pd (_mm_loadu_pd(C+2), _mm_mul_pd (A1,B12));
 
    } else {
        C11 = _mm_mul_pd (A1,B11);
        C12 = _mm_mul_pd (A1,B12);
    }
    
    C11 = _mm_add_pd (C11, _mm_mul_pd (A2,B21));
    C12 = _mm_add_pd (C12, _mm_mul_pd (A2,B22));
    
    B11 = _mm_loadu_pd (B + S2);
    B12 = _mm_loadu_pd (B + S2 + 2);
    B21 = _mm_loadu_pd (B + S3);
    B22 = _mm_loadu_pd (B + S3 + 2);
                  
    C11 = _mm_add_pd (C11, _mm_mul_pd (A3,B11));
    C12 = _mm_add_pd (C12, _mm_mul_pd (A3,B12));
    C11 = _mm_add_pd (C11, _mm_mul_pd (A4,B21));
    C12 = _mm_add_pd (C12, _mm_mul_pd (A4,B22));
    
    _mm_storeu_pd (C, C11);
    _mm_storeu_pd (C+2, C12);
}


void _hy_matrix_multiply_1x4x1 (double *C, double* A, double *B, int stride, bool add) {
    
    double s1  = (A[0] * B[0]),
           s2  = A[1] * B[stride],
           s3  = A[2] * B [stride << 1],
           s4  = A[3] * B [(stride << 1) + stride];
    
    s1 += s3; s2 += s4;
    if (add) {
        C[0] += s1+s2;
    } else {
        C[0] = s1+s2;
    }
}


void _hy_matrix_multiply_4x1x1 (double *C, double* A, double *B, int stride) {

    // C[0] = A[0]*B[0]
    // C[1] = A[1]*B[0]
    // C[2] = A[2]*B[0]
    // C[3] = A[3]*B[0]
    
    
    
    int S1 = stride,
        S2 = stride << 1,
        S3 = stride + S2;

    __m128d         B1 = _mm_load1_pd (B),
                    C1,C2;

    C1 = _mm_mul_pd (B1,_mm_loadh_pd (_mm_load_sd (A), A+S1));
    C2 = _mm_mul_pd (B1,_mm_loadh_pd (_mm_load_sd (A+S2), A+S3));
    
    C[0]  += _mm_cvtsd_f64 (C1);
    C[S1] += _mm_cvtsd_f64 ( _mm_unpackhi_pd (C1,C1));
    C[S2] += _mm_cvtsd_f64 (C2);
    C[S3] += _mm_cvtsd_f64 ( _mm_unpackhi_pd (C2,C2));
}

void _hy_matrix_multiply_1x1x4 (double *C, double* A, double *B, int stride) {

    __m128d         A1 = _mm_load1_pd (A),
                    C1 = _mm_loadu_pd (C),
                    C2 = _mm_loadu_pd (C+2);
                    
    _mm_storeu_pd (C, _mm_add_pd(_mm_mul_pd (A1,_mm_loadu_pd (B)), C1));
    _mm_storeu_pd (C+2, _mm_add_pd(_mm_mul_pd (A1,_mm_loadu_pd (B+2)), C2));

 }

/**
        Nx2 cases
 */

void _hy_matrix_multiply_4x2x4 (double *C, double* A, double *B, int stride) {

   int   S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22, C31, C32, C41, C42;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2

    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03


    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C32 = _mm_add_pd (C32,_mm_mul_pd (A3,B2)); // 20*02, 20*03
        C41 = _mm_add_pd (C41,_mm_mul_pd (A4,B1)); // 20*00, 20*01
        C42 = _mm_add_pd (C42,_mm_mul_pd (A4,B2)); // 20*02, 20*03
    };
        
    C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
    C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
    C31 = _mm_loadu_pd (C+S2); C32 = _mm_loadu_pd (C+S2+2);
    C41 = _mm_loadu_pd (C+S3); C42 = _mm_loadu_pd (C+S3+2);
   
    handle_block_madd ();
   
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S2+2, C32);
    _mm_storeu_pd (C+S3,   C41);
    _mm_storeu_pd (C+S3+2, C42);
}

void _hy_matrix_multiply_4x4x2 (double *C, double* A, double *B, int stride, bool add) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C21, C31, C41;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2

    B1 = _mm_loadu_pd (B);      // 00,01
 

    auto handle_block_mult = [&] () -> void {
        C11 = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C21 = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C31 = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C41 = _mm_mul_pd (A4,B1); // 20*00, 20*01
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C41 = _mm_add_pd (C41,_mm_mul_pd (A4,B1)); // 20*00, 20*01
    };
        
    if (add) {
        C11 = _mm_loadu_pd (C);
        C21 = _mm_loadu_pd (C+S1);
        C31 = _mm_loadu_pd (C+S2);
        C41 = _mm_loadu_pd (C+S3);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+2);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B2 = _mm_loadu_pd    (B+S2+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+3);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+3);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+3);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+3);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S3);      // 10,11
    B2 = _mm_loadu_pd    (B+S3+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S3,   C41);
}

void _hy_matrix_multiply_2x4x2 (double *C, double* A, double *B, int stride, bool add) {

    if (add) {
        for (int i = 0; i < 2; i++) {
#pragma unroll
            for (int k = 0; k < 4; k++) {
                for (int j=0; j < 2; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
#pragma unroll
            for (int j=0; j < 2; j++) {
                C[i*stride + j]  = A[i*stride] * B[j];
                for (int k = 1; k < 4; k++)  {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }

    }
}

void _hy_matrix_multiply_2x4x4 (double *C, double* A, double *B, int stride, bool add) {
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
 
    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03


    auto handle_block_mult = [&] () -> void {
        C11 = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C12 = _mm_mul_pd (A1,B2); // 00*02, 00*03
        C21 = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C22 = _mm_mul_pd (A2,B2); // 10*02, 10*03
     };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
     };
        
    if (add) {
        C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
        C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][1] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B2 = _mm_loadu_pd    (B+S2+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+3);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+3);    // A[1][1] x2
    B1 = _mm_loadu_pd    (B+S3);      // 10,11
    B2 = _mm_loadu_pd    (B+S3+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
}

void _hy_matrix_multiply_2x2x4 (double *C, double* A, double *B, int stride) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
 
    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03


     
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
     };
        
    C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
    C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
    handle_block_madd ();
    
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
}



void _hy_matrix_multiply_4x2x2 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 4; i++) {
        #pragma unroll
        for (int k = 0; k < 2; k++) {
            for (int j=0; j < 2; j++) {
                C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
            }
        }
    }
    
   /* int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4, A5, A6, A7, A8;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C21, C31, C41;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A5 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A6 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A7 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2
    A8 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2

    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+S1);
 
   
    C11 = _mm_add_pd ( _mm_loadu_pd (C),_mm_mul_pd (A1,B1)); // 00*00, 00*01
    C21 = _mm_add_pd ( _mm_loadu_pd (C+S1),_mm_mul_pd (A2,B1)); // 10*00, 10*01
    C31 = _mm_add_pd ( _mm_loadu_pd (C+S2),_mm_mul_pd (A3,B1)); // 20*00, 20*01
    C41 = _mm_add_pd ( _mm_loadu_pd (C+S3),_mm_mul_pd (A4,B1));
          

    C11 = _mm_add_pd ( C11,_mm_mul_pd (A5,B2)); // 00*00, 00*01
    C21 = _mm_add_pd ( C21,_mm_mul_pd (A6,B2)); // 10*00, 10*01
    C31 = _mm_add_pd ( C31,_mm_mul_pd (A7,B2)); // 20*00, 20*01
    C41 = _mm_add_pd ( C41,_mm_mul_pd (A8,B2));
 
    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S3,   C41);*/
}

void _hy_matrix_multiply_2x2x2 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                for (int j=0; j < 2; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
            for (int j=0; j < 2; j++) {
                C[i*stride + j] = A[i*stride] * B[j] + A[i*stride+1] * B[stride+j];
            }
        }
    }
}

/**
 Nx3 cases
*/

void _hy_matrix_multiply_4x3x4 (double *C, double* A, double *B, int stride) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22, C31, C32, C41, C42;
    
    // 00*00 + 01*10 + 02*20
    // 00*01 + 01*11 + 01*21
    // ....
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2

    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03

    auto handle_block_mult = [&] () -> void {
        C11 = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C12 = _mm_mul_pd (A1,B2); // 00*02, 00*03
        C21 = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C22 = _mm_mul_pd (A2,B2); // 10*02, 10*03
        C31 = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C32 = _mm_mul_pd (A3,B2); // 20*02, 20*03
        C41 = _mm_mul_pd (A4,B1); // 20*00, 20*01
        C42 = _mm_mul_pd (A4,B2); // 20*02, 20*03
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C32 = _mm_add_pd (C32,_mm_mul_pd (A3,B2)); // 20*02, 20*03
        C41 = _mm_add_pd (C41,_mm_mul_pd (A4,B1)); // 20*00, 20*01
        C42 = _mm_add_pd (C42,_mm_mul_pd (A4,B2)); // 20*02, 20*03
    };
        
    C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
    C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
    C31 = _mm_loadu_pd (C+S2); C32 = _mm_loadu_pd (C+S2+2);
    C41 = _mm_loadu_pd (C+S3); C42 = _mm_loadu_pd (C+S3+2);
    handle_block_madd ();
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+2);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B2 = _mm_loadu_pd    (B+S2+2);    // 12,13
      
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S2+2, C32);
    _mm_storeu_pd (C+S3,   C41);
    _mm_storeu_pd (C+S3+2, C42);
}

void _hy_matrix_multiply_4x4x3 (double *C, double* A, double *B, int stride, bool add) {
    // 00*00 + 01*10 + 02*20 + 03*30
    // 00*01 + 01*11 + 02*21 + 03*31
    // 02 = 00*02 + 01*12 + 02*22 + 03*32
    // 12 = 10*02 + 11*12 + 12*22 + 13*23
    // 22 = 20*02 + 21*12 + 22*22 + 23*23
    // 32 = 30*02 + 31*12 + 32*22 + 33*32
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B_LC1, B_LC2;       // current row in B (first two elements)
    __m128d     C11, C21, C31, C41, C_LC1, C_LC2;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2

    B1 = _mm_loadu_pd (B);      // 00,01
    
    B_LC1 = _mm_loaddup_pd (B+2);   // 02, 02
 
    auto handle_block_mult = [&] () -> void {
        C11    = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C21    = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C31    = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C41    = _mm_mul_pd (A4,B1); // 20*00, 20*01
        C_LC1  = _mm_mul_pd (B_LC1,  _mm_blend_pd (A1,A2,2));
        C_LC2  = _mm_mul_pd (B_LC1,  _mm_blend_pd (A3,A4,2));
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C41 = _mm_add_pd (C41,_mm_mul_pd (A4,B1)); // 20*00, 20*01
        C_LC1  = _mm_add_pd (C_LC1, _mm_mul_pd (B_LC1,  _mm_blend_pd (A1,A2,2)));
        C_LC2  = _mm_add_pd (C_LC2, _mm_mul_pd (B_LC1,  _mm_blend_pd (A3,A4,2)));
    };
        
    if (add) {
        C11  = _mm_loadu_pd (C);
        C21  = _mm_loadu_pd (C+S1);
        C31  = _mm_loadu_pd (C+S2);
        C41  = _mm_loadu_pd (C+S3);
        C_LC1 =  _mm_loadh_pd (_mm_load_sd (C+2), C+S1+2);
        C_LC2 =  _mm_loadh_pd (_mm_load_sd (C+S2+2), C+S3+2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B_LC1 = _mm_loaddup_pd (B+S1+2);   // 02, 02
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][2] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][3] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][2] x2
    A4 = _mm_loaddup_pd  (A+S3+2);    // A[3][2] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B_LC1 = _mm_loaddup_pd (B+S2+2);   // 02, 02
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+3);       // A[0][3] x2
    A2 = _mm_loaddup_pd  (A+S1+3);    // A[1][3] x2
    A3 = _mm_loaddup_pd  (A+S2+3);    // A[2][3] x2
    A4 = _mm_loaddup_pd  (A+S3+3);    // A[3][3] x2
    B1 = _mm_loadu_pd    (B+S3);      // 10,11
    B_LC1 = _mm_loaddup_pd (B+S3+2);   // 02, 02
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S3,   C41);
    _mm_storel_pd (C+2,    C_LC1);
    _mm_storeh_pd (C+S1+2, C_LC1);
    _mm_storel_pd (C+S2+2, C_LC2);
    _mm_storeh_pd (C+S3+2, C_LC2);

}

void _hy_matrix_multiply_4x3x3 (double *C, double* A, double *B, int stride, bool add) {
    // 00*00 + 01*10 + 02*20
    // 00*01 + 01*11 + 02*21
    // 02 = 00*02 + 01*12 + 02*22
    // 12 = 10*02 + 11*12 + 12*22
    // 22 = 20*02 + 21*12 + 22*22
    // 32 = 30*02 + 31*12 + 32*22
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3, A4;
    __m128d     B1, B_LC1, B_LC2;       // current row in B (first two elements)
    __m128d     C11, C21, C31, C41, C_LC1, C_LC2;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
    A4 = _mm_loaddup_pd  (A+S3);    // A[3][0] x2

    B1 = _mm_loadu_pd (B);      // 00,01
    
    B_LC1 = _mm_loaddup_pd (B+2);   // 02, 02
 
    auto handle_block_mult = [&] () -> void {
        C11    = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C21    = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C31    = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C41    = _mm_mul_pd (A4,B1); // 20*00, 20*01
        C_LC1  = _mm_mul_pd (B_LC1,  _mm_blend_pd (A1,A2,2));
        C_LC2  = _mm_mul_pd (B_LC1,  _mm_blend_pd (A3,A4,2));
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C41 = _mm_add_pd (C41,_mm_mul_pd (A4,B1)); // 20*00, 20*01
        C_LC1  = _mm_add_pd (C_LC1, _mm_mul_pd (B_LC1,  _mm_blend_pd (A1,A2,2)));
        C_LC2  = _mm_add_pd (C_LC2, _mm_mul_pd (B_LC1,  _mm_blend_pd (A3,A4,2)));
    };
        
    if (add) {
        C11  = _mm_loadu_pd (C);
        C21  = _mm_loadu_pd (C+S1);
        C31  = _mm_loadu_pd (C+S2);
        C41  = _mm_loadu_pd (C+S3);
        C_LC1 =  _mm_loadh_pd (_mm_load_sd (C+2), C+S1+2);
        C_LC2 =  _mm_loadh_pd (_mm_load_sd (C+S2+2), C+S3+2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    A4 = _mm_loaddup_pd  (A+S3+1);    // A[3][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B_LC1 = _mm_loaddup_pd (B+S1+2);   // 02, 02
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][2] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][3] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][2] x2
    A4 = _mm_loaddup_pd  (A+S3+2);    // A[3][2] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B_LC1 = _mm_loaddup_pd (B+S2+2);   // 02, 02
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S3,   C41);
    _mm_storel_pd (C+2,    C_LC1);
    _mm_storeh_pd (C+S1+2, C_LC1);
    _mm_storel_pd (C+S2+2, C_LC2);
    _mm_storeh_pd (C+S3+2, C_LC2);
    
}

void _hy_matrix_multiply_3x3x4 (double *C, double* A, double *B, int stride, bool add) {
    //      00*00 + 01*10 + 02*20
    //      00*01 + 01*11 + 02*21
    // 02 = 00*02 + 01*12 + 02*22
    // 03 = 00*03 + 01*13 + 02*23
    
    // 13 = 10*03 + 11*13 + 12*23
    // 23 = 20*03 + 21*13 + 22*23
    
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22, C31, C32;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
 
    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03

    auto handle_block_mult = [&] () -> void {
        C11 = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C12 = _mm_mul_pd (A1,B2); // 00*02, 00*03
        C21 = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C22 = _mm_mul_pd (A2,B2); // 10*02, 10*03
        C31 = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C32 = _mm_mul_pd (A3,B2); // 20*02, 20*03
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C32 = _mm_add_pd (C32,_mm_mul_pd (A3,B2)); // 20*02, 20*03
    };
        
    if (add) {
        C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
        C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
        C31 = _mm_loadu_pd (C+S2); C32 = _mm_loadu_pd (C+S2+2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][1] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B2 = _mm_loadu_pd    (B+S2+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S2+2, C32);
}

void _hy_matrix_multiply_3x4x4 (double *C, double* A, double *B, int stride, bool add) {
     int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m128d     A1, A2, A3;
    __m128d     B1, B2;       // current row in B
    __m128d     C11, C12, C21, C22, C31, C32;
    
    A1 = _mm_loaddup_pd  (A);       // A[0][0] x2
    A2 = _mm_loaddup_pd  (A+S1);    // A[1][0] x2
    A3 = _mm_loaddup_pd  (A+S2);    // A[2][0] x2
 
    B1 = _mm_loadu_pd (B);      // 00,01
    B2 = _mm_loadu_pd (B+2);    // 02,03

    auto handle_block_mult = [&] () -> void {
        C11 = _mm_mul_pd (A1,B1); // 00*00, 00*01
        C12 = _mm_mul_pd (A1,B2); // 00*02, 00*03
        C21 = _mm_mul_pd (A2,B1); // 10*00, 10*01
        C22 = _mm_mul_pd (A2,B2); // 10*02, 10*03
        C31 = _mm_mul_pd (A3,B1); // 20*00, 20*01
        C32 = _mm_mul_pd (A3,B2); // 20*02, 20*03
    };
    
    auto handle_block_madd = [&] () -> void {
        C11 = _mm_add_pd (C11,_mm_mul_pd (A1,B1)); // 00*00, 00*01
        C12 = _mm_add_pd (C12,_mm_mul_pd (A1,B2)); // 00*02, 00*03
        C21 = _mm_add_pd (C21,_mm_mul_pd (A2,B1)); // 10*00, 10*01
        C22 = _mm_add_pd (C22,_mm_mul_pd (A2,B2)); // 10*02, 10*03
        C31 = _mm_add_pd (C31,_mm_mul_pd (A3,B1)); // 20*00, 20*01
        C32 = _mm_add_pd (C32,_mm_mul_pd (A3,B2)); // 20*02, 20*03
    };
        
    if (add) {
        C11 = _mm_loadu_pd (C);    C12 = _mm_loadu_pd (C+2);
        C21 = _mm_loadu_pd (C+S1); C22 = _mm_loadu_pd (C+S1+2);
        C31 = _mm_loadu_pd (C+S2); C32 = _mm_loadu_pd (C+S2+2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm_loaddup_pd  (A+1);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+1);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+1);    // A[2][1] x2
    B1 = _mm_loadu_pd    (B+S1);      // 10,11
    B2 = _mm_loadu_pd    (B+S1+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+2);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+2);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+2);    // A[2][1] x2
    B1 = _mm_loadu_pd    (B+S2);      // 10,11
    B2 = _mm_loadu_pd    (B+S2+2);    // 12,13
   
    handle_block_madd ();

    A1 = _mm_loaddup_pd  (A+3);       // A[0][1] x2
    A2 = _mm_loaddup_pd  (A+S1+3);    // A[1][1] x2
    A3 = _mm_loaddup_pd  (A+S2+3);    // A[2][1] x2
    B1 = _mm_loadu_pd    (B+S3);      // 10,11
    B2 = _mm_loadu_pd    (B+S3+2);    // 12,13
   
    handle_block_madd ();

    _mm_storeu_pd (C,      C11);
    _mm_storeu_pd (C+2,    C12);
    _mm_storeu_pd (C+S1,   C21);
    _mm_storeu_pd (C+S1+2, C22);
    _mm_storeu_pd (C+S2,   C31);
    _mm_storeu_pd (C+S2+2, C32);
}

void _hy_matrix_multiply_3x3x3 (double *C, double* A, double *B, int stride, bool add) {
    
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                C[i*stride]      += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
             }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride]   * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]);
            }
        }
    }
     
}

void _hy_matrix_multiply_3x4x3 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 4; k++) {
                C[i*stride ]     += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j=0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride] * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]) +
                                   (A[i*stride+3] * B[(stride<<1) + stride + j]);
            }
        }
    }
}

#endif

//====================================================================================================================
//====================================================================================================================
//====================================================================================================================


#if defined _SLKP_USE_AVX_INTRINSICS

#include <immintrin.h>

inline __m256d _hy_matrix_handle_axv_mfma (__m256d c, __m256d a, __m256d b) {
    #if defined _SLKP_USE_FMA3_INTRINSICS
        return _mm256_fmadd_pd (a,b,c);
    #else
        return _mm256_add_pd (c, _mm256_mul_pd (a,b));
    #endif
}

void _hy_matrix_multiply_4x4 (double * C, double *A, double *B, int stride, bool add) {
    // 4x4 NEON multiplication using NEON intrinsics
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;


    __m256d     A1, A2, A3, A4;
    __m256d     B1;         // current row in B
    __m256d     C1, C2, C3, C4;
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4
    A4 = _mm256_broadcast_sd  (A+S3);    // A[3][0] x4

    B1 = _mm256_loadu_pd (B);      // 00,01
 
    auto handle_block_mult = [&] () -> void {
        C1 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C2 = _mm256_mul_pd (A2,B1); // 10*01, 01*11, 01*12, 01*13
        C3 = _mm256_mul_pd (A3,B1); // 10*00, 10*01
        C4 = _mm256_mul_pd (A4,B1); // 10*02, 10*03
     };
    
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
        C4 = _hy_matrix_handle_axv_mfma (C4,A4,B1); // 10*02, 10*03
    };
        
    if (add) {
        C1 = _mm256_loadu_pd (C);
        C2 = _mm256_loadu_pd (C+S1);
        C3 = _mm256_loadu_pd (C+S2);
        C4 = _mm256_loadu_pd (C+S3);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+1);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S1);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+2);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S2);      // 00,01

    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+3);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+3);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+3);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+3);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S3);      // 00,01

   
    handle_block_madd ();

    _mm256_storeu_pd  (C,      C1);
    _mm256_storeu_pd  (C+S1,   C2);
    _mm256_storeu_pd  (C+S2,   C3);
    _mm256_storeu_pd  (C+S3,   C4);
    
}



/**
        Nx1 cases
*/


void _hy_matrix_multiply_4x1x4 (double *C, double* A, double *B, int stride) {
    __m256d   A1,A2,A3,A4;
    __m256d   B1;
    __m256d   C1,C2,C3,C4;
    
    int S1 = stride, S2 = stride << 1, S3 = stride + S2;
    
    A1  = _mm256_broadcast_sd (A);     // A[0,0]
    A2  = _mm256_broadcast_sd (A+S1);  // A[1,0]
    A3  = _mm256_broadcast_sd (A+S2);  // A[2,0]
    A4  = _mm256_broadcast_sd (A+S3);  // A[3,0]
    
    B1 = _mm256_loadu_pd (B);
         
    _mm256_storeu_pd (C, _hy_matrix_handle_axv_mfma (_mm256_loadu_pd (C), A1, B1));
    _mm256_storeu_pd (C+S1, _hy_matrix_handle_axv_mfma (_mm256_loadu_pd (C+S1), A2, B1));
    _mm256_storeu_pd (C+S2, _hy_matrix_handle_axv_mfma (_mm256_loadu_pd (C+S2), A3, B1));
    _mm256_storeu_pd (C+S3, _hy_matrix_handle_axv_mfma (_mm256_loadu_pd (C+S3), A4, B1));
    
}


void _hy_matrix_multiply_4x4x1 (double *C, double* A, double *B, int stride, bool add) {
        int  S1 = stride,
             S2 = stride << 1,
             S3 = S2 + stride;
         
         // 00 = 00*00 + 01*10 + 02*20 + 03*30
         // 10 = 10*00 + 11*10 + 12*20 + 13*30
         // 20 = 20*00 + ... + 23*30
         // 30 = 30*00 + ... + 33*30
         
    
         
         
        __m256d   A1,A2,A3,A4;
        __m256d   B1;
        __m256i   LOAD_IDX = _mm256_set_epi64x (S3,S2,S1,0);
    
        B1 = _mm256_i64gather_pd (B,LOAD_IDX,8);
        
        A1 = _mm256_mul_pd(_mm256_loadu_pd (A),    B1);
        A2 = _mm256_mul_pd(_mm256_loadu_pd (A+S1), B1);
        A3 = _mm256_mul_pd(_mm256_loadu_pd (A+S2), B1);
        A4 = _mm256_mul_pd(_mm256_loadu_pd (A+S3), B1);
    

        A1 = _mm256_hadd_pd (A1,A2); // A1[0]+A1[1], A2[0]+A2[1], A1[2]+A1[3], A2[2]+A2[3]
        A3 = _mm256_hadd_pd (A3,A4);
         
        A2 = _mm256_permute4x64_pd (A1,0+4*2+16*1+64*3); // A1[0]+A1[1], A1[2]+A1[3], A2[0]+A2[1], A2[2]+A2[3]
        A4 = _mm256_permute4x64_pd (A3,0+4*2+16*1+64*3);
        
        A1 = _mm256_permute4x64_pd(_mm256_hadd_pd (A2,A4),0+4*2+16*1+64*3); // A1[0]+A1[1] + A1[2]+A1[3], A2[0]+A2[1] +  A2[2]+A2[3] , x2
        
        
        __m128d C11 = _mm256_extractf128_pd (A1,0);
        __m128d C21 = _mm256_extractf128_pd (A1,1);

        if (add) {
            C[0]  += _mm_cvtsd_f64 (C11);
            C[S1] += _mm_cvtsd_f64(_mm_unpackhi_pd (C11,C11));
            C[S2] += _mm_cvtsd_f64 (C21);
            C[S3] += _mm_cvtsd_f64(_mm_unpackhi_pd (C21,C21));
                
        } else {
            _mm_storel_pd (C,    C11);
            _mm_storeh_pd (C+S1, C11);
            _mm_storel_pd (C+S2, C21);
            _mm_storeh_pd (C+S3, C21);
        }
    
    
}


void _hy_matrix_multiply_1x4x4 (double *C, double* A, double *B, int stride, bool add) {
    
    int S1 = stride,
        S2 = stride << 1,
        S3 = stride + S2;

    __m256d       A1 = _mm256_broadcast_sd (A),
                  A2 = _mm256_broadcast_sd (A+1),
                  A3 = _mm256_broadcast_sd (A+2),
                  A4 = _mm256_broadcast_sd (A+3);
    
    __m256d       B1 = _mm256_loadu_pd (B);
                     
    // 00 = 00 * 00 + 01 * 10 + 02 * 20 + 03 * 30
    // 01 = 00 * 01 + 01 * 11 + 02 * 21 + 03 * 31
    // 02 = 00 * 02 + 01 * 12 + 02 * 22 + 03 * 32
    // 03 = 00 * 03 + 01 * 13 + 02 * 23 + 03 * 33

    __m256d       C1, C2;
     
    if (add) {
        C1= _hy_matrix_handle_axv_mfma (_mm256_loadu_pd (C), A1,B1);
    } else {
        C1 = _mm256_mul_pd (A1,B1);
    }
    
    C2 = _mm256_mul_pd (A2,_mm256_loadu_pd (B+S1));
    
    C1 = _hy_matrix_handle_axv_mfma (C1, A3,_mm256_loadu_pd (B+S2));
    C2 = _hy_matrix_handle_axv_mfma (C2, A4,_mm256_loadu_pd (B+S3));
    
    _mm256_storeu_pd (C, _mm256_add_pd (C1,C2));
    
}


void _hy_matrix_multiply_1x4x1 (double *C, double* A, double *B, int stride, bool add) {
    
    double s1  = (A[0] * B[0]),
           s2  = A[1] * B[stride],
           s3  = A[2] * B [stride << 1],
           s4  = A[3] * B [(stride << 1) + stride];
    
    s1 += s3; s2 += s4;
    if (add) {
        C[0] += s1+s2;
    } else {
        C[0] = s1+s2;
    }
}


void _hy_matrix_multiply_4x1x1 (double *C, double* A, double *B, int stride) {
    int S1 = stride,
        S2 = stride << 1,
        S3 = stride + S2;

    __m256d  C1,
             B1= _mm256_broadcast_sd (B);
             
    __m256i  LOAD_IDX = _mm256_set_epi64x (S3,S2,S1,0);
    
    C1 =  _hy_matrix_handle_axv_mfma (_mm256_i64gather_pd (C,LOAD_IDX,8), _mm256_i64gather_pd (A,LOAD_IDX,8),B1);
    
    __m128d C11 = _mm256_extractf128_pd (C1,0);
    __m128d C21 = _mm256_extractf128_pd (C1,1);
    _mm_storel_pd (C,    C11);
    _mm_storeh_pd (C+S1, C11);
    _mm_storel_pd (C+S2, C21);
    _mm_storeh_pd (C+S3, C21);
}

void _hy_matrix_multiply_1x1x4 (double *C, double* A, double *B, int stride) {
    _mm256_storeu_pd (C,_hy_matrix_handle_axv_mfma (_mm256_loadu_pd (C), _mm256_broadcast_sd (A), _mm256_loadu_pd (B)));
 }

/**
        Nx2 cases
 */

void _hy_matrix_multiply_4x2x4 (double *C, double* A, double *B, int stride) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;


    __m256d     A1, A2, A3, A4;
    __m256d     B1             ; // current row in B
    __m256d     C1, C2, C3, C4;
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4
    A4 = _mm256_broadcast_sd  (A+S3);    // A[3][0] x4

    B1 = _mm256_loadu_pd (B);      // 00,01
 
     
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
        C4 = _hy_matrix_handle_axv_mfma (C4,A4,B1); // 10*02, 10*03
    };
        
    C1 = _mm256_loadu_pd (C);
    C2 = _mm256_loadu_pd (C+S1);
    C3 = _mm256_loadu_pd (C+S2);
    C4 = _mm256_loadu_pd (C+S3);
    handle_block_madd ();
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+1);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S1);      // 00,01
   
    handle_block_madd ();

    _mm256_storeu_pd  (C,      C1);
    _mm256_storeu_pd  (C+S1,   C2);
    _mm256_storeu_pd  (C+S2,   C3);
    _mm256_storeu_pd  (C+S3,   C4);
}

void _hy_matrix_multiply_4x4x2 (double *C, double* A, double *B, int stride, bool add) {

// 00 = 00*00 + 01*10 + 02*20 + 03*30
// 01 = 00*01 + 01*11 + 02*21 + 03*30

// 10 = 10*00 + 11*10 + 12*20 + 13*30
// 11 = 10*01 + 11*11 + 12*21 + 13*31

// 20 = 20*00 + 21*10 + 22*20 + 23*30
// 21 = 20*01 + 21*11 + 22*21 + 23*31

// 30 = 30*00 + 31*10 + 32*20 + 33*30
// 31 = 30*01 + 31*11 + 32*21 + 33*31

    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;


    __m256d     A1, A2;
    __m256d     B1;                 // current row in B
    __m256d     C12,C34;
    
    A1 = _mm256_set_m128d  (_mm_loaddup_pd (A+S1), _mm_loaddup_pd (A));       // 00,00,10,10
    A2 = _mm256_set_m128d  (_mm_loaddup_pd (A+S3), _mm_loaddup_pd (A+S2));
    B1 = _mm256_broadcast_pd  ((__m128d const*)B);      // 00,01 x 2
     
    auto handle_block_madd = [&] () -> void {
        C12 = _hy_matrix_handle_axv_mfma (C12,A1,B1); // 00*00, 00*01
        C34 = _hy_matrix_handle_axv_mfma (C34,A2,B1); // 01*10, 01*11, 01*12, 01*13
     };
     auto handle_block_mult = [&] () -> void {
        C12 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C34 = _mm256_mul_pd (A2,B1); // 01*10, 01*11, 01*12, 01*13
     };
     
     if (add) {
         C12 =  _mm256_set_m128d  (_mm_loadu_pd (C+S1), _mm_loadu_pd (C));
         C34 =  _mm256_set_m128d  (_mm_loadu_pd (C+S3), _mm_loadu_pd (C+S2));
         handle_block_madd ();
    } else {
         handle_block_mult ();
    }
    
    A1 = _mm256_set_m128d  (_mm_loaddup_pd (A+S1+1), _mm_loaddup_pd (A+1));       // 00,00,10,10
    A2 = _mm256_set_m128d  (_mm_loaddup_pd (A+S3+1), _mm_loaddup_pd (A+S2+1));
    B1 = _mm256_broadcast_pd  ((__m128d const*)(B+S1));      // 00,01 x 2
    handle_block_madd ();

    A1 = _mm256_set_m128d  (_mm_loaddup_pd (A+S1+2), _mm_loaddup_pd (A+2));       // 00,00,10,10
    A2 = _mm256_set_m128d  (_mm_loaddup_pd (A+S3+2), _mm_loaddup_pd (A+S2+2));
    B1 = _mm256_broadcast_pd  ((__m128d const*)(B+S2));      // 00,01 x 2
    handle_block_madd ();

    A1 = _mm256_set_m128d  (_mm_loaddup_pd (A+S1+3), _mm_loaddup_pd (A+3));       // 00,00,10,10
    A2 = _mm256_set_m128d  (_mm_loaddup_pd (A+S3+3), _mm_loaddup_pd (A+S2+3));
    B1 = _mm256_broadcast_pd  ((__m128d const*)(B+S3));      // 00,01 x 2
    handle_block_madd ();

    _mm_store_pd     (C, _mm256_extractf128_pd (C12,0));
    _mm_store_pd     (C+S1, _mm256_extractf128_pd (C12,1));
    _mm_store_pd     (C+S2, _mm256_extractf128_pd (C34,0));
    _mm_store_pd     (C+S3, _mm256_extractf128_pd (C34,1));

    
}

void _hy_matrix_multiply_2x4x2 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 2; i++) {
#pragma unroll
            for (int k = 0; k < 4; k++) {
                for (int j=0; j < 2; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
#pragma unroll
            for (int j=0; j < 2; j++) {
                C[i*stride + j]  = A[i*stride] * B[j];
                for (int k = 1; k < 4; k++)  {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }

    }
}

void _hy_matrix_multiply_2x4x4 (double *C, double* A, double *B, int stride, bool add) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;


    __m256d     A1, A2, A3, A4;
    __m256d     B1;         // current row in B
    __m256d     C1, C2;
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4
    A4 = _mm256_broadcast_sd  (A+S3);    // A[3][0] x4

    B1 = _mm256_loadu_pd (B);      // 00,01
 
    auto handle_block_mult = [&] () -> void {
        C1 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C2 = _mm256_mul_pd (A2,B1); // 10*01, 01*11, 01*12, 01*13
     };
    
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
    };
        
    if (add) {
        C1 = _mm256_loadu_pd (C);
        C2 = _mm256_loadu_pd (C+S1);
         handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+1);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S1);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+2);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S2);      // 00,01

    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+3);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+3);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+3);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+3);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S3);      // 00,01

   
    handle_block_madd ();

    _mm256_storeu_pd  (C,      C1);
    _mm256_storeu_pd  (C+S1,   C2);
}

void _hy_matrix_multiply_2x2x4 (double *C, double* A, double *B, int stride) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m256d     A1, A2;
    __m256d     B1;
    __m256d     C1, C2;
    
    // 00 = 00 * 00 + 01 * 10
    // 01 = 00 * 01 + 01 * 11
    // 02 = 00 * 02 + 01 * 12
    // 03 = 00 * 03 + 01 * 13

    // 10 = 10 * 00 + 11 * 10
    // 11 = 10 * 01 + 11 * 11
    // 12 = 10 * 02 + 11 * 12
    // 13 = 10 * 03 + 11 * 13
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x2
 
    B1 = _mm256_loadu_pd (B);
 

     
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1, A1, B1);
        C2 = _hy_matrix_handle_axv_mfma (C2, A2, B1);
     };
        
    C1 = _mm256_loadu_pd (C);
    C2 = _mm256_loadu_pd (C+S1);
    
    handle_block_madd ();
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][0] x2
    B1 = _mm256_loadu_pd (B+S1);

    handle_block_madd ();

    _mm256_storeu_pd (C,      C1);
    _mm256_storeu_pd (C+S1,   C2);
}


void _hy_matrix_multiply_4x2x2 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 4; i++) {
#pragma unroll
        for (int k = 0; k < 2; k++) {
            for (int j=0; j < 2; j++) {
                C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
            }
        }
    }
    
    
}

void _hy_matrix_multiply_2x2x2 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                for (int j=0; j < 2; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
            for (int j=0; j < 2; j++) {
                C[i*stride + j] = A[i*stride] * B[j] + A[i*stride+1] * B[stride+j];
            }
        }
    }
}

/**
 Nx3 cases
*/

void _hy_matrix_multiply_4x3x4 (double *C, double* A, double *B, int stride) {
   
   int   S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m256d     A1, A2, A3, A4;
    __m256d     B1;         // current row in B
    __m256d     C1, C2, C3, C4;
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4
    A4 = _mm256_broadcast_sd  (A+S3);    // A[3][0] x4

    B1 = _mm256_loadu_pd (B);      // 00,01
     
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
        C4 = _hy_matrix_handle_axv_mfma (C4,A4,B1); // 10*02, 10*03
    };
        
    C1 = _mm256_loadu_pd (C);
    C2 = _mm256_loadu_pd (C+S1);
    C3 = _mm256_loadu_pd (C+S2);
    C4 = _mm256_loadu_pd (C+S3);
    handle_block_madd ();
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+1);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S1);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+2);    // A[3][1] x2
    B1 = _mm256_loadu_pd (B+S2);      // 00,01

    handle_block_madd ();

    _mm256_storeu_pd  (C,      C1);
    _mm256_storeu_pd  (C+S1,   C2);
    _mm256_storeu_pd  (C+S2,   C3);
    _mm256_storeu_pd  (C+S3,   C4);
    
}

void _hy_matrix_multiply_4x4x3 (double *C, double* A, double *B, int stride, bool add) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m256d     A1, A2, A3, A4;
    __m256d     B1;  // last value will be junk
    __m256d     C1, C2, C3, C4; // last value in each row will be junk
 
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4
    A4 = _mm256_broadcast_sd  (A+S3);    // A[3][0] x4

    __m256i   MASK = _mm256_set_epi64x (0,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF);

    B1 = _mm256_maskload_pd  (B,MASK);      // 00,01
 
    auto handle_block_mult = [&] () -> void {
        C1 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C2 = _mm256_mul_pd (A2,B1); // 10*01, 01*11, 01*12, 01*13
        C3 = _mm256_mul_pd (A3,B1); // 10*00, 10*01
        C4 = _mm256_mul_pd (A4,B1); // 10*02, 10*03
     };
    
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
        C4 = _hy_matrix_handle_axv_mfma (C4,A4,B1); // 10*02, 10*03
    };
        
    if (add) {
        C1 = _mm256_loadu_pd (C);
        C2 = _mm256_loadu_pd (C+S1);
        C3 = _mm256_loadu_pd (C+S2);
        C4 = _mm256_loadu_pd (C+S3);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+1);    // A[3][1] x2
    B1 = _mm256_maskload_pd (B+S1, MASK);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+2);    // A[3][1] x2
    B1 = _mm256_maskload_pd (B+S2, MASK);      // 00,01

    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+3);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+3);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+3);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+3);    // A[3][1] x2
    B1 = _mm256_maskload_pd (B+S3, MASK);      // 00,01

    handle_block_madd ();

    _mm256_maskstore_pd  (C,      MASK, C1);
    _mm256_maskstore_pd  (C+S1,   MASK, C2);
    _mm256_maskstore_pd  (C+S2,   MASK, C3);
    _mm256_maskstore_pd  (C+S3,   MASK, C4);
    
}

void _hy_matrix_multiply_4x3x3 (double *C, double* A, double *B, int stride, bool add) {
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    __m256d     A1, A2, A3, A4;
    __m256d     B1;  // last value will be junk
    __m256d     C1, C2, C3, C4; // last value in each row will be junk
 
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4
    A4 = _mm256_broadcast_sd  (A+S3);    // A[3][0] x4

    __m256i   MASK = _mm256_set_epi64x (0,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF);

    B1 = _mm256_maskload_pd  (B,MASK);      // 00,01
 
    auto handle_block_mult = [&] () -> void {
        C1 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C2 = _mm256_mul_pd (A2,B1); // 10*01, 01*11, 01*12, 01*13
        C3 = _mm256_mul_pd (A3,B1); // 10*00, 10*01
        C4 = _mm256_mul_pd (A4,B1); // 10*02, 10*03
     };
    
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
        C4 = _hy_matrix_handle_axv_mfma (C4,A4,B1); // 10*02, 10*03
    };
        
    if (add) {
        C1 = _mm256_loadu_pd (C);
        C2 = _mm256_loadu_pd (C+S1);
        C3 = _mm256_loadu_pd (C+S2);
        C4 = _mm256_loadu_pd (C+S3);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+1);    // A[3][1] x2
    B1 = _mm256_maskload_pd (B+S1, MASK);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    A4 = _mm256_broadcast_sd  (A+S3+2);    // A[3][1] x2
    B1 = _mm256_maskload_pd (B+S2, MASK);      // 00,01

    handle_block_madd ();


    _mm256_maskstore_pd  (C,      MASK, C1);
    _mm256_maskstore_pd  (C+S1,   MASK, C2);
    _mm256_maskstore_pd  (C+S2,   MASK, C3);
    _mm256_maskstore_pd  (C+S3,   MASK, C4);
    
}

void _hy_matrix_multiply_3x3x4 (double *C, double* A, double *B, int stride, bool add) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;


    __m256d     A1, A2, A3;
    __m256d     B1;         // current row in B
    __m256d     C1, C2, C3;
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4

    B1 = _mm256_loadu_pd (B);      // 00,01
 
    auto handle_block_mult = [&] () -> void {
        C1 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C2 = _mm256_mul_pd (A2,B1); // 10*01, 01*11, 01*12, 01*13
        C3 = _mm256_mul_pd (A3,B1); // 10*00, 10*01
     };
    
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
    };
        
    if (add) {
        C1 = _mm256_loadu_pd (C);
        C2 = _mm256_loadu_pd (C+S1);
        C3 = _mm256_loadu_pd (C+S2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    B1 = _mm256_loadu_pd (B+S1);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    B1 = _mm256_loadu_pd (B+S2);      // 00,01

    handle_block_madd ();

    _mm256_storeu_pd  (C,      C1);
    _mm256_storeu_pd  (C+S1,   C2);
    _mm256_storeu_pd  (C+S2,   C3);
}

void _hy_matrix_multiply_3x4x4 (double *C, double* A, double *B, int stride, bool add) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;


    __m256d     A1, A2, A3;
    __m256d     B1;         // current row in B
    __m256d     C1, C2, C3;
    
    A1 = _mm256_broadcast_sd  (A);       // A[0][0] x4
    A2 = _mm256_broadcast_sd  (A+S1);    // A[1][0] x4
    A3 = _mm256_broadcast_sd  (A+S2);    // A[2][0] x4

    B1 = _mm256_loadu_pd (B);      // 00,01
 
    auto handle_block_mult = [&] () -> void {
        C1 = _mm256_mul_pd (A1,B1); // 00*00, 00*01
        C2 = _mm256_mul_pd (A2,B1); // 10*01, 01*11, 01*12, 01*13
        C3 = _mm256_mul_pd (A3,B1); // 10*00, 10*01
     };
    
    auto handle_block_madd = [&] () -> void {
        C1 = _hy_matrix_handle_axv_mfma (C1,A1,B1); // 00*00, 00*01
        C2 = _hy_matrix_handle_axv_mfma (C2,A2,B1); // 01*10, 01*11, 01*12, 01*13
        C3 = _hy_matrix_handle_axv_mfma (C3,A3,B1); // 10*00, 10*01
    };
        
    if (add) {
        C1 = _mm256_loadu_pd (C);
        C2 = _mm256_loadu_pd (C+S1);
        C3 = _mm256_loadu_pd (C+S2);
        handle_block_madd ();
    } else {
        handle_block_mult ();
    }
    
    A1 = _mm256_broadcast_sd  (A+1);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+1);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+1);    // A[2][1] x2
    B1 = _mm256_loadu_pd (B+S1);      // 00,01
   
    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+2);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+2);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+2);    // A[2][1] x2
    B1 = _mm256_loadu_pd (B+S2);      // 00,01

    handle_block_madd ();

    A1 = _mm256_broadcast_sd  (A+3);       // A[0][1] x2
    A2 = _mm256_broadcast_sd  (A+S1+3);    // A[1][1] x2
    A3 = _mm256_broadcast_sd  (A+S2+3);    // A[2][1] x2
    B1 = _mm256_loadu_pd (B+S3);      // 00,01

   
    handle_block_madd ();

    _mm256_storeu_pd  (C,      C1);
    _mm256_storeu_pd  (C+S1,   C2);
    _mm256_storeu_pd  (C+S2,   C3);
}

void _hy_matrix_multiply_3x3x3 (double *C, double* A, double *B, int stride, bool add) {
    
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                C[i*stride]      += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
             }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride]   * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]);
            }
        }
    }
     
     
}

void _hy_matrix_multiply_3x4x3 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 4; k++) {
                C[i*stride ]     += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j=0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride] * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]) +
                                   (A[i*stride+3] * B[(stride<<1) + stride + j]);
            }
        }
    }
}

#endif

#if not defined _SLKP_USE_ARM_NEON && not defined _SLKP_USE_SSE_INTRINSICS && not defined _SLKP_USE_AVX_INTRINSICS

void _hy_matrix_multiply_4x4 (double * C, double *A, double *B, int stride, bool add) {
    // 4x4 NEON multiplication using NEON intrinsics
    
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    if (add) {
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                C[r*stride+c] += (A[r*stride]   * B[c]) +
                                 (A[r*stride+1] * B[c+S1]) +
                                 (A[r*stride+2] * B[c+S2]) +
                                 (A[r*stride+3] * B[c+S3]);
            }
        }
    } else {
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                C[r*stride+c] = (A[r*stride]   * B[c]) +
                                (A[r*stride+1] * B[c+S1]) +
                                (A[r*stride+2] * B[c+S2]) +
                                (A[r*stride+3] * B[c+S3]);
            }
        }

    }
}



/**
        Nx1 cases
*/


void _hy_matrix_multiply_4x1x4 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 4; i++) {
        C[i*stride]  += A[i*stride] * B[0];
        C[i*stride + 1]  += A[i*stride] * B[1];
        C[i*stride + 2]  += A[i*stride] * B[2];
        C[i*stride + 3]  += A[i*stride] * B[3];
   }
}


void _hy_matrix_multiply_4x4x1 (double *C, double* A, double *B, int stride, bool add) {
    int  S1 = stride,
         S2 = stride << 1,
         S3 = S2 + stride;

    if (add)
        for (int i = 0; i < 4; i++) {
            C[i*stride] += (A[i*stride] * B[0]) +
                       (A[i*stride+1] * B[S1]) +
                       (A[i*stride+2] * B[S2]) +
                       (A[i*stride+3] * B[S3]);
        }
    else {
        for (int i = 0; i < 4; i++) {
            C[i*stride] = (A[i*stride] * B[0]) +
                       (A[i*stride+1] * B[S1]) +
                       (A[i*stride+2] * B[S2]) +
                       (A[i*stride+3] * B[S3]);
        }

    }
    
}


void _hy_matrix_multiply_1x4x4 (double *C, double* A, double *B, int stride, bool add) {
    
    if (add) {
        for (int k = 0; k < 4; k++) {
#pragma unroll
            for (int j = 0; j < 4; j++) {
                C[j]  += A[k] * B[k*stride+j];
            }
        }
    } else {
         C[0] = A[0]*B[0];
         C[1] = A[0]*B[1];
         C[2] = A[0]*B[2];
         C[3] = A[0]*B[3];
         for (int k = 1; k < 4; k++) {
             for (int j = 0; j < 4; j++) {
                 C[j]  += A[k] * B[k*stride+j];
             }
         }

    }
    return;
    
    
    
}


void _hy_matrix_multiply_1x4x1 (double *C, double* A, double *B, int stride, bool add) {
    
    double s1  = (A[0] * B[0]),
           s2  = A[1] * B[stride],
           s3  = A[2] * B [stride << 1],
           s4  = A[3] * B [(stride << 1) + stride];
    
    s1 += s3; s2 += s4;
    if (add) {
        C[0] += s1+s2;
    } else {
        C[0] = s1+s2;
    }
}


void _hy_matrix_multiply_4x1x1 (double *C, double* A, double *B, int stride) {
#pragma unroll
    for (int i = 0; i < 4; i++) {
        C[i*stride]  += A[i*stride] * B[0];
    }

}

void _hy_matrix_multiply_1x1x4 (double *C, double* A, double *B, int stride) {
    C[0]  += A[0] * B[0];
    C[1]  += A[0] * B[1];
    C[2]  += A[0] * B[2];
    C[3]  += A[0] * B[3];
 }

/**
        Nx2 cases
 */

void _hy_matrix_multiply_4x2x4 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 2; k++) {
#pragma unroll
            for (int j=0; j < 4; j++) {
                C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
            }
        }
    }
}

void _hy_matrix_multiply_4x4x2 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 4; i++) {
#pragma unroll
                for (int k = 0; k < 4; k++) {
                    C[i*stride ]  += A[i*stride+k] * B[k*stride];
                    C[i*stride  + 1]  += A[i*stride+k] * B[k*stride + 1];
                }

        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int j=0; j < 2; j++){
                 C[i*stride + j]  = A[i*stride] * B[j];
                 for (int k = 1; k < 4; k++)  {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    }

    
}

void _hy_matrix_multiply_2x4x2 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 2; i++) {
#pragma unroll
            for (int k = 0; k < 4; k++) {
                for (int j=0; j < 2; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
#pragma unroll
            for (int j=0; j < 2; j++) {
                C[i*stride + j]  = A[i*stride] * B[j];
                for (int k = 1; k < 4; k++)  {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }

    }
}

void _hy_matrix_multiply_2x4x4 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 4; k++) {
#pragma unroll
                for (int j=0; j < 4; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++) {
                C[i*stride + j]  = A[i*stride] * B[j];
#pragma unroll
                for (int k=1; k < 4; k++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }

    }
}

void _hy_matrix_multiply_2x2x4 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 2; i++) {
        for (int k = 0; k < 2; k++) {
#pragma unroll
            for (int j=0; j < 4; j++) {
                C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
            }
        }
    }
}



void _hy_matrix_multiply_4x2x2 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 4; i++) {
#pragma unroll
        for (int k = 0; k < 2; k++) {
            for (int j=0; j < 2; j++) {
                C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
            }
        }
    }
    
    
}

void _hy_matrix_multiply_2x2x2 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                for (int j=0; j < 2; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
            for (int j=0; j < 2; j++) {
                C[i*stride + j] = A[i*stride] * B[j] + A[i*stride+1] * B[stride+j];
            }
        }
    }
}

/**
 Nx3 cases
*/

void _hy_matrix_multiply_4x3x4 (double *C, double* A, double *B, int stride) {
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 3; k++) {
#pragma unroll
            for (int j=0; j < 4; j++) {
                C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
            }
        }
    }
}

void _hy_matrix_multiply_4x4x3 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
#pragma unroll
                for (int j=0; j < 3; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                C[i*stride + j]  = A[i*stride] * B[j];
#pragma unroll
                for (int k=1; k < 4; k++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }

    }
     
}

void _hy_matrix_multiply_4x3x3 (double *C, double* A, double *B, int stride, bool add) {
    
    if (add) {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 3; k++) {
                for (int j=0; j < 3; j++) {
                    C[i*stride + j]  += A[i*stride+k] * B[k*stride+j];
                }
            }
        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride] * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]);
            }
        }
    }
}

void _hy_matrix_multiply_3x3x4 (double *C, double* A, double *B, int stride, bool add) {
    
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                C[i*stride]      += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
                C[i*stride + 3]  += A[i*stride+k] * B[k*stride+3];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                C[i*stride + j]  = (A[i*stride] * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]);
            }
        }

    }
     
     
}

void _hy_matrix_multiply_3x4x4 (double *C, double* A, double *B, int stride, bool add) {
    
    
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 4; k++) {
                C[i*stride]      += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
                C[i*stride + 3]  += A[i*stride+k] * B[k*stride+3];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j=0; j < 4; j++) {
                C[i*stride + j]  = (A[i*stride] * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]) +
                                   (A[i*stride+3] * B[(stride<<1) + stride + j]);
            }
        }
    }
}

void _hy_matrix_multiply_3x3x3 (double *C, double* A, double *B, int stride, bool add) {
    
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 3; k++) {
                C[i*stride]      += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
             }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride]   * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]);
            }
        }
    }
     
     
}

void _hy_matrix_multiply_3x4x3 (double *C, double* A, double *B, int stride, bool add) {
    if (add) {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < 4; k++) {
                C[i*stride ]     += A[i*stride+k] * B[k*stride];
                C[i*stride + 1]  += A[i*stride+k] * B[k*stride+1];
                C[i*stride + 2]  += A[i*stride+k] * B[k*stride+2];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int j=0; j < 3; j++) {
                C[i*stride + j]  = (A[i*stride] * B[j]) +
                                   (A[i*stride+1] * B[stride + j]) +
                                   (A[i*stride+2] * B[(stride<<1) + j]) +
                                   (A[i*stride+3] * B[(stride<<1) + stride + j]);
            }
        }
    }
}

#endif

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------

void _hy_matrix_multiply_NxN_blocked4 (double * C, double *A, double *B, int D) {
    
#ifdef _SLKP_USE_APPLE_BLAS
    if (D>=16) {
        //printf ("HERE\n");
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, D, D, D, 1., A, D, B, D, 0.0, C, D);
        return;
    }
#endif
    auto offset = [D](int i, int j) -> int {
        return (i<<2)*D + (j << 2);
    };
    
    int blocks = D>>2;
    int remainder = D - (blocks<<2);
    // fill in the first blocks in C
    
                              
    switch (remainder) {
        case 0: {
            for (int i = 0; i < blocks; i++) {
                for (int k = 0; k < blocks; k++) {
                    for (int j=0; j < blocks; j++) {
                        _hy_matrix_multiply_4x4 (C + offset (i,j), A + offset (i,k), B + offset (k,j), D, k > 0);
                    }
                }
            }
            break;
        }
        case 1: {
            
            if (blocks >= 1) {
                for (int i = 0; i < blocks; i++) {
                    for (int k = 0; k < blocks; k++) {
                        for (int j=0; j < blocks; j++) {
                            _hy_matrix_multiply_4x4 (C + offset (i,j), A + offset (i,k), B + offset (k,j), D, k > 0);
                            if (k == blocks - 1) {
                                _hy_matrix_multiply_4x1x4 (C + offset (i,j), A + offset (i,blocks), B + offset (blocks,j), D);
                            }
                        }
                    }
                }
            } else {
                C[0] = A[0] * B[0];
                return;
            }
            
            // last REMAINDER columns in the result matrix
            for (int i = 0; i < blocks; i++) {
                _hy_matrix_multiply_4x4x1 (C + offset (i,blocks), A + offset (i,0), B + offset (0,blocks), D, false);
                for (int k = 1; k < blocks; k++) {
                    _hy_matrix_multiply_4x4x1 (C + offset (i,blocks), A + offset (i,k), B + offset (k,blocks), D, true);
                }
                _hy_matrix_multiply_4x1x1 (C + offset (i,blocks), A + offset (i,blocks), B + offset (blocks,blocks), D);
            }
   
            _hy_matrix_multiply_1x4x1 (C + offset (blocks,blocks), A + offset (blocks,0), B + offset (0,blocks), D, false);
            for (int k = 1; k < blocks; k++) {
                _hy_matrix_multiply_1x4x1 (C + offset (blocks,blocks), A + offset (blocks,k), B + offset (k,blocks), D, true);
            }
            
            int last_element = offset (blocks,blocks);
            C [last_element] += A[last_element]*B[last_element];

            // last REMAINER ROWS in the result matrix
            for (int i = blocks-1; i >= 0; i--) {
                _hy_matrix_multiply_1x4x4 (C + offset (blocks,i), A + offset (blocks,0), B + offset (0,i), D, false);
                for (int k = 1; k < blocks; k++) {
                    _hy_matrix_multiply_1x4x4 (C + offset (blocks,i), A + offset (blocks,k), B + offset (k,i), D, true);
                }
                _hy_matrix_multiply_1x1x4 (C + offset (blocks,i), A + offset (blocks,blocks), B + offset (blocks,i), D);
            }
            
            break;
        }
        case 2: {
        // add overhanging pieces to the main blocks
            if (blocks >= 1) {
                for (int i = 0; i < blocks; i++) {
                    for (int k = 0; k < blocks; k++) {
                        for (int j=0; j < blocks; j++) {
                            _hy_matrix_multiply_4x4 (C + offset (i,j), A + offset (i,k), B + offset (k,j), D, k>0);
                            if (k == blocks - 1) {
                                _hy_matrix_multiply_4x2x4 (C + offset (i,j), A + offset (i,blocks), B + offset (blocks,j), D);
                            }
                        }
                    }
                }
            } else {
                _hy_matrix_multiply_2x2x2 (C, A, B, D, false);
                return;
            }
            
            // last REMAINDER columns in the result matrix
            for (int i = 0; i < blocks; i++) {
                _hy_matrix_multiply_4x4x2 (C + offset (i,blocks), A + offset (i,0), B + offset (0,blocks), D, false);
                for (int k = 1; k < blocks; k++) {
                    _hy_matrix_multiply_4x4x2 (C + offset (i,blocks), A + offset (i,k), B + offset (k,blocks), D, true);
                }
                _hy_matrix_multiply_4x2x2 (C + offset (i,blocks), A + offset (i,blocks), B + offset (blocks,blocks), D);
            }

            _hy_matrix_multiply_2x4x2(C + offset (blocks,blocks), A + offset (blocks,0), B + offset (0,blocks), D, false);
            for (int k = 1; k < blocks; k++) {
                _hy_matrix_multiply_2x4x2 (C + offset (blocks,blocks), A + offset (blocks,k), B + offset (k,blocks), D, true);
            }
            //handle_small_product_add<2,2,2> (C + offset (blocks,blocks), A + offset (blocks,blocks), B + offset (blocks,blocks), D);
            _hy_matrix_multiply_2x2x2 (C + offset (blocks,blocks), A + offset (blocks,blocks), B + offset (blocks,blocks), D, true);
            
            

            // last REMAINER ROWS in the result matrix
            for (int i = blocks-1; i >= 0; i--) {
                _hy_matrix_multiply_2x4x4(C + offset (blocks,i), A + offset (blocks,0), B + offset (0,i), D, false);
                for (int k = 1; k < blocks; k++) {
                    //handle_small_product_add <2,4,4>  (C + offset (blocks,i), A + offset (blocks,k), B + offset (k,i), D);
                    _hy_matrix_multiply_2x4x4 (C + offset (blocks,i), A + offset (blocks,k), B + offset (k,i), D, true);
                }
                _hy_matrix_multiply_2x2x4(C + offset (blocks,i), A + offset (blocks,blocks), B + offset (blocks,i), D);
            }
            
            break;
        }
        case 3: {
            // add overhanging pieces to the main blocks
            if (blocks >= 1) {
                for (int i = 0; i < blocks; i++) {
                    for (int k = 0; k < blocks; k++) {
                        for (int j=0; j < blocks; j++) {
                            _hy_matrix_multiply_4x4 (C + offset (i,j), A + offset (i,k), B + offset (k,j), D, k>0);
                            if (k == blocks - 1) {
                                _hy_matrix_multiply_4x3x4 (C + offset (i,j), A + offset (i,blocks), B + offset (blocks,j), D);
                            }
                        }
                    }
                }
            } else {
                _hy_matrix_multiply_3x3x3 (C, A, B, D, false);
                return;
            }
            
            // last REMAINDER columns in the result matrix
            for (int i = 0; i < blocks; i++) {
                _hy_matrix_multiply_4x4x3 (C + offset (i,blocks), A + offset (i,0), B + offset (0,blocks), D,false);
                for (int k = 1; k < blocks; k++) {
                    _hy_matrix_multiply_4x4x3 (C + offset (i,blocks), A + offset (i,k), B + offset (k,blocks), D, true);
                }
                _hy_matrix_multiply_4x3x3 (C + offset (i,blocks), A + offset (i,blocks), B + offset (blocks,blocks), D, true);
            }
            
            _hy_matrix_multiply_3x4x3 (C + offset (blocks,blocks), A + offset (blocks,0), B + offset (0,blocks), D, false);
            for (int k = 1; k < blocks; k++) {
                //handle_small_product_add <3,4,3> (C + offset (blocks,blocks), A + offset (blocks,k), B + offset (k,blocks), D);
                _hy_matrix_multiply_3x4x3 (C + offset (blocks,blocks), A + offset (blocks,k), B + offset (k,blocks), D, true);
            }
            //handle_small_product_add<3,3,3> (C + offset (blocks,blocks), A + offset (blocks,blocks), B + offset (blocks,blocks), D);
            _hy_matrix_multiply_3x3x3 (C + offset (blocks,blocks), A + offset (blocks,blocks), B + offset (blocks,blocks), D, true);
            
            // last REMAINER ROWS in the result matrix
            for (int i = blocks-1; i >= 0; i--) {
                _hy_matrix_multiply_3x4x4 (C + offset (blocks,i), A + offset (blocks,0), B + offset (0,i), D, false);
                for (int k = 1; k < blocks; k++) {
                    _hy_matrix_multiply_3x4x4  (C + offset (blocks,i), A + offset (blocks,k), B + offset (k,i), D,true);
                }
                //handle_small_product_add<3,3,4> (C + offset (blocks,i), A + offset (blocks,blocks), B + offset (blocks,i), D);
                _hy_matrix_multiply_3x3x4 (C + offset (blocks,i), A + offset (blocks,blocks), B + offset (blocks,i), D, true);
            }
            
            // HANDLE BLOCK/BLOCK
           

            break;
        }
    }
}

template<long B1, long B2> inline void mx_transpose_blocked (double  * __restrict C, double const  * __restrict A, int nrow, int ncol) {
    /* nrow and ncol is for matrix A
        "A" points to (r,c) in A (nrow X ncol dimension)
        "C" points to (c,r) in C (ncol x nrol dimension)
        
        A[r][c], A[r][c+1], A[r][c+2], A[r][c+3]
        ...
        A[r+3][c], .... A[r+3][c+3]
        
        This will get written to
        
        C[c][r], C[c][r+1], C[c][r+2], C[c][r+3]
        ...
        C[c+3][r], ... C[c+3][r+3]
    */
    
    
    for (int j = 0; j < B2; j++) {
        #pragma unroll
        for (int i = 0; i < B1; i++) {
            C[j*nrow+i] = A[i*ncol+j];
        }
    }
}


void _hy_matrix_transpose_blocked (double * __restrict C, double * __restrict A, int nrow, int ncol) {
    
    auto offset_A = [ncol](int i, int j) -> int {
        return (i<<2)*ncol + (j << 2);
    };
    
    auto offset_C = [nrow](int i, int j) -> int {
        return (i<<2)*nrow + (j << 2);
    };
    
    int row_blocks = nrow>>2;
    int col_blocks = ncol>>2;
    
    int remainder_rows    = nrow - (row_blocks<<2);
    int remainder_columns = ncol - (col_blocks<<2);
    
    for (int r = 0; r < row_blocks; r++) {
        for (int c = 0; c < col_blocks; c++) {
            mx_transpose_blocked<4,4> (C+ offset_C (c,r), A + offset_A(r,c), nrow, ncol);
        }
    }
    
    if (remainder_rows || remainder_columns) {
        switch (remainder_rows) {
            case 1:
                for (int c = 0; c < col_blocks; c++) {
                    mx_transpose_blocked<1,4> (C + offset_C (c,row_blocks), A + offset_A(row_blocks,c), nrow, ncol);
                }
                break;
            case 2:
                for (int c = 0; c < col_blocks; c++) {
                    mx_transpose_blocked<2,4> (C + offset_C (c,row_blocks), A + offset_A(row_blocks,c), nrow, ncol);
                }
                break;
            case 3:
                for (int c = 0; c < col_blocks; c++) {
                    mx_transpose_blocked<3,4> (C + offset_C (c,row_blocks), A + offset_A(row_blocks,c), nrow, ncol);
                }
                break;
        }

        switch (remainder_columns) {
            case 1:
                for (int r = 0; r < row_blocks; r++) {
                    mx_transpose_blocked<4,1> (C + offset_C (col_blocks,r), A + offset_A(r,col_blocks), nrow, ncol);
                }
                switch (remainder_rows) {
                    case 1:
                        mx_transpose_blocked<1,1> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                    case 2:
                        mx_transpose_blocked<2,1> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                    case 3:
                        mx_transpose_blocked<3,1> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                }
                break;
            case 2:
                for (int r = 0; r < row_blocks; r++) {
                    mx_transpose_blocked<4,2> (C + offset_C (col_blocks,r), A + offset_A(r,col_blocks), nrow, ncol);
                }
                switch (remainder_rows) {
                    case 1:
                        mx_transpose_blocked<1,2> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                    case 2:
                        mx_transpose_blocked<2,2> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                    case 3:
                        mx_transpose_blocked<3,2> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                }
                break;
            case 3:
                for (int r = 0; r < row_blocks; r++) {
                    mx_transpose_blocked<4,3> (C + offset_C (col_blocks,r), A + offset_A(r,col_blocks), nrow, ncol);
                }
                switch (remainder_rows) {
                    case 1:
                        mx_transpose_blocked<1,3> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                    case 2:
                        mx_transpose_blocked<2,3> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                    case 3:
                        mx_transpose_blocked<3,3> (C + offset_C (col_blocks,row_blocks), A + offset_A(row_blocks,col_blocks), nrow, ncol);
                        break;
                }
                break;
        }
    }
}

// Matrix-Vector blocked


