	.arch armv8-a
; GNU C++17 (Homebrew GCC 12.2.0) version 12.2.0 (aarch64-apple-darwin21)
;	compiled by GNU C version 12.2.0, GMP version 6.2.1, MPFR version 4.1.0, MPC version 1.2.1, isl version isl-0.25-GMP

; warning: MPFR header version 4.1.0 differs from library version 4.1.0-p13.
; GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
; options passed: -fPIC -mmacosx-version-min=13.0.0 -mlittle-endian -mabi=lp64 -O3
	.text
	.align	2
	.p2align 4,,11
	.globl __Z31__ll_handle_block10_product_sumILl61ELl0EE13__Float64x2_tPKdPdS3_PS0_
	.weak_definition __Z31__ll_handle_block10_product_sumILl61ELl0EE13__Float64x2_tPKdPdS3_PS0_
__Z31__ll_handle_block10_product_sumILl61ELl0EE13__Float64x2_tPKdPdS3_PS0_:
LFB4351:
; product_sum.cc:8:     float64x2_t P [5] = {
	movi	v1.2d, 0	; P$2
; product_sum.cc:18:     for (long col_idx = 0; col_idx < D; col_idx++, tT+=D) {
	mov	x4, 0	; col_idx,
; product_sum.cc:8:     float64x2_t P [5] = {
	mov	v2.16b, v1.16b	; P$3, P$2
	mov	v3.16b, v1.16b	; P$4, P$2
	mov	v5.16b, v1.16b	; P$0, P$2
	mov	v4.16b, v1.16b	; P$1, P$2
	.p2align 3,,7
L4:
; product_sum.cc:19:         if (childVector[col_idx] > 0.) {
	ldr	d0, [x1, x4, lsl 3]	; _5, MEM[(hyFloat *)childVector_46(D) + _101 * 8]
; product_sum.cc:19:         if (childVector[col_idx] > 0.) {
	fcmpe	d0, #0.0	; _5
	bgt	L6		;,
L2:
; product_sum.cc:18:     for (long col_idx = 0; col_idx < D; col_idx++, tT+=D) {
	add	x4, x4, 1	; col_idx, col_idx,
; product_sum.cc:18:     for (long col_idx = 0; col_idx < D; col_idx++, tT+=D) {
	add	x0, x0, 488	; transposedMatrix, transposedMatrix,
; product_sum.cc:18:     for (long col_idx = 0; col_idx < D; col_idx++, tT+=D) {
	cmp	x4, 61	; col_idx,
	bne	L4		;,
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:15269:   return __builtin_aarch64_ld1v2df ((const __builtin_aarch64_simd_df *) __a);
	ldp	q0, q17, [x2]	; MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D)], MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 16B], MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D)]
	ldp	q16, q7, [x2, 32]	; MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 32B], MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 48B], MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 32B]
	ldr	q6, [x2, 64]	; MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 64B], MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 64B]
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:1003:   return __a * __b;
	fmul	v0.2d, v0.2d, v5.2d	; _57, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D)], P$0
	fmul	v5.2d, v17.2d, v4.2d	; _55, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 16B], P$1
	fmul	v3.2d, v6.2d, v3.2d	; _47, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 64B], P$4
	fmul	v4.2d, v16.2d, v1.2d	; _27, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 32B], P$2
	fmul	v1.2d, v7.2d, v2.2d	; _29, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 48B], P$3
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:24715:   __builtin_aarch64_st1v2df ((__builtin_aarch64_simd_df *) __a, __b);
	stp	q0, q5, [x2]	; _57, _55, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D)]
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:345:   return __a + __b;
	fadd	v0.2d, v5.2d, v0.2d	; tmp137, _55, _57
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:24715:   __builtin_aarch64_st1v2df ((__builtin_aarch64_simd_df *) __a, __b);
	str	q3, [x2, 64]	; _47, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 64B]
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:345:   return __a + __b;
	fadd	v2.2d, v4.2d, v1.2d	; tmp139, _27, _29
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:24715:   __builtin_aarch64_st1v2df ((__builtin_aarch64_simd_df *) __a, __b);
	stp	q4, q1, [x2, 32]	; _27, _29, MEM <__Float64x2_t> [(double * {ref-all})parentConditionals_32(D) + 32B]
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:345:   return __a + __b;
	fadd	v0.2d, v0.2d, v3.2d	; tmp138, tmp137, _47
	fadd	v0.2d, v0.2d, v2.2d	; <retval>, tmp138, tmp139
; product_sum.cc:47:     if (grandTotal) {
	cbz	x3, L1	; grandTotal,
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:345:   return __a + __b;
	ldr	q1, [x3]	; *grandTotal_44(D), *grandTotal_44(D)
	fadd	v0.2d, v0.2d, v1.2d	; <retval>, <retval>, *grandTotal_44(D)
; product_sum.cc:48:         *grandTotal = vaddq_f64 (P[2],*grandTotal);
	str	q0, [x3]	; <retval>, *grandTotal_44(D)
L1:
; product_sum.cc:53: }
	ret	
	.p2align 2,,3
L6:
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:15269:   return __builtin_aarch64_ld1v2df ((const __builtin_aarch64_simd_df *) __a);
	ldp	q18, q17, [x0]	; MEM <__Float64x2_t> [(double * {ref-all})tT_84], MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 16B], MEM <__Float64x2_t> [(double * {ref-all})tT_84]
	ldp	q16, q7, [x0, 32]	; MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 32B], MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 48B], MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 32B]
	ldr	q6, [x0, 64]	; MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 64B], MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 64B]
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:14612:   return __builtin_aarch64_fmav2df (__b, __c, __a);
	fmla	v5.2d, v18.2d, v0.d[0]	; P$0, MEM <__Float64x2_t> [(double * {ref-all})tT_84], _5
	fmla	v4.2d, v17.2d, v0.d[0]	; P$1, MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 16B], _5
	fmla	v1.2d, v16.2d, v0.d[0]	; P$2, MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 32B], _5
	fmla	v2.2d, v7.2d, v0.d[0]	; P$3, MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 48B], _5
	fmla	v3.2d, v6.2d, v0.d[0]	; P$4, MEM <__Float64x2_t> [(double * {ref-all})tT_84 + 64B], _5
	b	L2		;
LFE4351:
	.section __TEXT,__eh_frame,coalesced,no_toc+strip_static_syms+live_support
EH_frame1:
	.set L$set$0,LECIE1-LSCIE1
	.long L$set$0
LSCIE1:
	.long	0
	.byte	0x1
	.ascii "zR\0"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x1e
	.uleb128 0x1
	.byte	0x10
	.byte	0xc
	.uleb128 0x1f
	.uleb128 0
	.align	3
LECIE1:
LSFDE1:
	.set L$set$1,LEFDE1-LASFDE1
	.long L$set$1
LASFDE1:
	.long	LASFDE1-EH_frame1
	.quad	LFB4351-.
	.set L$set$2,LFE4351-LFB4351
	.quad L$set$2
	.uleb128 0
	.align	3
LEFDE1:
	.ident	"GCC: (Homebrew GCC 12.2.0) 12.2.0"
	.subsections_via_symbols
