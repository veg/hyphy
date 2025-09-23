	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 26, 0	sdk_version 26, 0
	.globl	__Z16_mx_vect_4x4_addR13float64x2x2_tPKdS2_i ; -- Begin function _Z16_mx_vect_4x4_addR13float64x2x2_tPKdS2_i
	.p2align	2
__Z16_mx_vect_4x4_addR13float64x2x2_tPKdS2_i: ; @_Z16_mx_vect_4x4_addR13float64x2x2_tPKdS2_i
	.cfi_startproc
; %bb.0:
                                        ; kill: def $w3 killed $w3 def $x3
	ld1.2d	{ v0, v1 }, [x2]
	lsl	w8, w3, #1
	sbfiz	x8, x8, #3, #32
	mov	x9, x1
	ld1.2d	{ v2, v3 }, [x9], x8
	sbfiz	x8, x3, #3, #32
	add	x10, x1, x8
	ld1.2d	{ v4, v5 }, [x10]
	ld1.2d	{ v6, v7 }, [x9], x8
	ld1.2d	{ v16, v17 }, [x9]
	fmul.2d	v18, v0, v2
	fmul.2d	v19, v0, v4
	fmla.2d	v18, v1, v3
	fmla.2d	v19, v1, v5
	fmul.2d	v2, v0, v6
	fmul.2d	v3, v0, v16
	fmla.2d	v2, v1, v7
	fmla.2d	v3, v1, v17
	faddp.2d	v0, v18, v19
	faddp.2d	v1, v2, v3
	ldp	q2, q3, [x0]
	fadd.2d	v0, v2, v0
	fadd.2d	v1, v3, v1
	stp	q0, q1, [x0]
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	_main                           ; -- Begin function main
	.p2align	2
_main:                                  ; @main
	.cfi_startproc
; %bb.0:
	mov	w0, #0                          ; =0x0
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
