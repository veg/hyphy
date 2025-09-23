	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 13, 0	sdk_version 13, 0
	.globl	__Z31__ll_handle_block10_product_sumILl61ELl0EE19__simd128_float64_tPKdPdS3_PS0_ ; -- Begin function _Z31__ll_handle_block10_product_sumILl61ELl0EE19__simd128_float64_tPKdPdS3_PS0_
	.weak_definition	__Z31__ll_handle_block10_product_sumILl61ELl0EE19__simd128_float64_tPKdPdS3_PS0_
	.p2align	2
__Z31__ll_handle_block10_product_sumILl61ELl0EE19__simd128_float64_tPKdPdS3_PS0_: ; @_Z31__ll_handle_block10_product_sumILl61ELl0EE19__simd128_float64_tPKdPdS3_PS0_
	.cfi_startproc
; %bb.0:
	add	x8, x0, #552
	movi.2d	v0, #0000000000000000
	mov	w9, #8
	movi.2d	v1, #0000000000000000
	movi.2d	v2, #0000000000000000
	movi.2d	v3, #0000000000000000
	movi.2d	v4, #0000000000000000
	b	LBB0_2
LBB0_1:                                 ;   in Loop: Header=BB0_2 Depth=1
	add	x9, x9, #16
	add	x8, x8, #976
LBB0_2:                                 ; =>This Inner Loop Header: Depth=1
	add	x10, x1, x9
	ldur	d5, [x10, #-8]
	fcmp	d5, #0.0
	b.le	LBB0_4
; %bb.3:                                ;   in Loop: Header=BB0_2 Depth=1
	sub	x10, x8, #552
	ldr	q6, [x10]
	fmla.2d	v4, v6, v5[0]
	sub	x10, x8, #536
	ldr	q6, [x10]
	fmla.2d	v3, v6, v5[0]
	sub	x10, x8, #520
	ldr	q6, [x10]
	fmla.2d	v2, v6, v5[0]
	sub	x10, x8, #504
	ldr	q6, [x10]
	fmla.2d	v1, v6, v5[0]
	sub	x10, x8, #488
	ldr	q6, [x10]
	fmla.2d	v0, v6, v5[0]
LBB0_4:                                 ;   in Loop: Header=BB0_2 Depth=1
	cmp	x9, #488
	b.eq	LBB0_7
; %bb.5:                                ;   in Loop: Header=BB0_2 Depth=1
	ldr	d5, [x1, x9]
	fcmp	d5, #0.0
	b.le	LBB0_1
; %bb.6:                                ;   in Loop: Header=BB0_2 Depth=1
	ldp	q6, q7, [x8, #-64]
	fmla.2d	v4, v6, v5[0]
	fmla.2d	v3, v7, v5[0]
	ldp	q6, q7, [x8, #-32]
	fmla.2d	v2, v6, v5[0]
	fmla.2d	v1, v7, v5[0]
	ldr	q6, [x8]
	fmla.2d	v0, v6, v5[0]
	b	LBB0_1
LBB0_7:
	ldp	q5, q6, [x2]
	fmul.2d	v4, v4, v5
	fmul.2d	v3, v3, v6
	ldp	q5, q6, [x2, #32]
	fmul.2d	v2, v2, v5
	fmul.2d	v1, v1, v6
	ldr	q5, [x2, #64]
	fmul.2d	v0, v0, v5
	stp	q4, q3, [x2]
	stp	q2, q1, [x2, #32]
	str	q0, [x2, #64]
	fadd.2d	v3, v4, v3
	fadd.2d	v1, v2, v1
	fadd.2d	v0, v3, v0
	fadd.2d	v0, v1, v0
	cbz	x3, LBB0_9
; %bb.8:
	ldr	q1, [x3]
	fadd.2d	v0, v0, v1
	str	q0, [x3]
LBB0_9:
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
