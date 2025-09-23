	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 13, 0	sdk_version 13, 1
	.globl	__Z14naive_multiplyPdPKdS1_i    ; -- Begin function _Z14naive_multiplyPdPKdS1_i
	.p2align	2
__Z14naive_multiplyPdPKdS1_i:           ; @_Z14naive_multiplyPdPKdS1_i
	.cfi_startproc
; %bb.0:
                                        ; kill: def $w3 killed $w3 def $x3
	cmp	w3, #1
	b.lt	LBB0_7
; %bb.1:
	mov	w8, #0
	sxtw	x10, w3
	mov	w9, w3
	lsl	x10, x10, #3
	lsl	x11, x9, #3
LBB0_2:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB0_3 Depth 2
                                        ;       Child Loop BB0_4 Depth 3
	mov	x12, #0
	mov	x13, x2
LBB0_3:                                 ;   Parent Loop BB0_2 Depth=1
                                        ; =>  This Loop Header: Depth=2
                                        ;       Child Loop BB0_4 Depth 3
	mov	x14, #0
	movi	d0, #0000000000000000
	mov	x15, x13
LBB0_4:                                 ;   Parent Loop BB0_2 Depth=1
                                        ;     Parent Loop BB0_3 Depth=2
                                        ; =>    This Inner Loop Header: Depth=3
	ldr	d1, [x1, x14]
	ldr	d2, [x15]
	fmul	d1, d1, d2
	fadd	d0, d0, d1
	add	x14, x14, #8
	add	x15, x15, x10
	cmp	x11, x14
	b.ne	LBB0_4
; %bb.5:                                ;   in Loop: Header=BB0_3 Depth=2
	str	d0, [x0], #8
	add	x12, x12, #1
	add	x13, x13, #8
	cmp	x12, x9
	b.ne	LBB0_3
; %bb.6:                                ;   in Loop: Header=BB0_2 Depth=1
	add	w8, w8, #1
	add	x1, x1, x10
	cmp	w8, w3
	b.ne	LBB0_2
LBB0_7:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	__Z19very_naive_multiplyPdPKdS1_i ; -- Begin function _Z19very_naive_multiplyPdPKdS1_i
	.p2align	2
__Z19very_naive_multiplyPdPKdS1_i:      ; @_Z19very_naive_multiplyPdPKdS1_i
	.cfi_startproc
; %bb.0:
	cmp	w3, #1
	b.lt	LBB1_7
; %bb.1:
	mov	x8, #0
	mov	w9, w3
	lsl	x10, x9, #3
LBB1_2:                                 ; =>This Loop Header: Depth=1
                                        ;     Child Loop BB1_3 Depth 2
                                        ;       Child Loop BB1_4 Depth 3
	mov	x11, #0
	mul	x12, x8, x9
	mov	x13, x2
LBB1_3:                                 ;   Parent Loop BB1_2 Depth=1
                                        ; =>  This Loop Header: Depth=2
                                        ;       Child Loop BB1_4 Depth 3
	movi	d0, #0000000000000000
	mov	x14, x9
	mov	x15, x1
	mov	x16, x13
LBB1_4:                                 ;   Parent Loop BB1_2 Depth=1
                                        ;     Parent Loop BB1_3 Depth=2
                                        ; =>    This Inner Loop Header: Depth=3
	ldr	d1, [x15], #8
	ldr	d2, [x16]
	fmul	d1, d1, d2
	fadd	d0, d0, d1
	add	x16, x16, x10
	subs	x14, x14, #1
	b.ne	LBB1_4
; %bb.5:                                ;   in Loop: Header=BB1_3 Depth=2
	add	x14, x11, x12
	str	d0, [x0, x14, lsl #3]
	add	x11, x11, #1
	add	x13, x13, #8
	cmp	x11, x9
	b.ne	LBB1_3
; %bb.6:                                ;   in Loop: Header=BB1_2 Depth=1
	add	x8, x8, #1
	add	x1, x1, x10
	cmp	x8, x9
	b.ne	LBB1_2
LBB1_7:
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	__Z18direct4x4_multiplyPdPKdS1_ ; -- Begin function _Z18direct4x4_multiplyPdPKdS1_
	.p2align	2
__Z18direct4x4_multiplyPdPKdS1_:        ; @_Z18direct4x4_multiplyPdPKdS1_
	.cfi_startproc
; %bb.0:
	ldr	d0, [x2]
	ldp	d1, d2, [x1]
	fmul	d0, d1, d0
	ldr	d1, [x2, #32]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #64]
	ldp	d2, d3, [x1, #16]
	fmul	d1, d2, d1
	ldr	d2, [x2, #96]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0]
	ldr	d0, [x2, #8]
	ldp	d1, d2, [x1]
	fmul	d0, d1, d0
	ldr	d1, [x2, #40]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #72]
	ldp	d2, d3, [x1, #16]
	fmul	d1, d2, d1
	ldr	d2, [x2, #104]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #8]
	ldr	d0, [x2, #16]
	ldp	d1, d2, [x1]
	fmul	d0, d1, d0
	ldr	d1, [x2, #48]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #80]
	ldp	d2, d3, [x1, #16]
	fmul	d1, d2, d1
	ldr	d2, [x2, #112]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #16]
	ldr	d0, [x2, #24]
	ldp	d1, d2, [x1]
	fmul	d0, d1, d0
	ldr	d1, [x2, #56]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #88]
	ldp	d2, d3, [x1, #16]
	fmul	d1, d2, d1
	ldr	d2, [x2, #120]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #24]
	ldr	d0, [x2]
	ldp	d1, d2, [x1, #32]
	fmul	d0, d1, d0
	ldr	d1, [x2, #32]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #64]
	ldp	d2, d3, [x1, #48]
	fmul	d1, d2, d1
	ldr	d2, [x2, #96]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #32]
	ldr	d0, [x2, #8]
	ldp	d1, d2, [x1, #32]
	fmul	d0, d1, d0
	ldr	d1, [x2, #40]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #72]
	ldp	d2, d3, [x1, #48]
	fmul	d1, d2, d1
	ldr	d2, [x2, #104]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #40]
	ldr	d0, [x2, #16]
	ldp	d1, d2, [x1, #32]
	fmul	d0, d1, d0
	ldr	d1, [x2, #48]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #80]
	ldp	d2, d3, [x1, #48]
	fmul	d1, d2, d1
	ldr	d2, [x2, #112]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #48]
	ldr	d0, [x2, #24]
	ldp	d1, d2, [x1, #32]
	fmul	d0, d1, d0
	ldr	d1, [x2, #56]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #88]
	ldp	d2, d3, [x1, #48]
	fmul	d1, d2, d1
	ldr	d2, [x2, #120]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #56]
	ldr	d0, [x2]
	ldp	d1, d2, [x1, #64]
	fmul	d0, d1, d0
	ldr	d1, [x2, #32]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #64]
	ldp	d2, d3, [x1, #80]
	fmul	d1, d2, d1
	ldr	d2, [x2, #96]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #64]
	ldr	d0, [x2, #8]
	ldp	d1, d2, [x1, #64]
	fmul	d0, d1, d0
	ldr	d1, [x2, #40]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #72]
	ldp	d2, d3, [x1, #80]
	fmul	d1, d2, d1
	ldr	d2, [x2, #104]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #72]
	ldr	d0, [x2, #16]
	ldp	d1, d2, [x1, #64]
	fmul	d0, d1, d0
	ldr	d1, [x2, #48]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #80]
	ldp	d2, d3, [x1, #80]
	fmul	d1, d2, d1
	ldr	d2, [x2, #112]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #80]
	ldr	d0, [x2, #24]
	ldp	d1, d2, [x1, #64]
	fmul	d0, d1, d0
	ldr	d1, [x2, #56]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #88]
	ldp	d2, d3, [x1, #80]
	fmul	d1, d2, d1
	ldr	d2, [x2, #120]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #88]
	ldr	d0, [x2]
	ldp	d1, d2, [x1, #96]
	fmul	d0, d1, d0
	ldr	d1, [x2, #32]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #64]
	ldp	d2, d3, [x1, #112]
	fmul	d1, d2, d1
	ldr	d2, [x2, #96]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #96]
	ldr	d0, [x2, #8]
	ldp	d1, d2, [x1, #96]
	fmul	d0, d1, d0
	ldr	d1, [x2, #40]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #72]
	ldp	d2, d3, [x1, #112]
	fmul	d1, d2, d1
	ldr	d2, [x2, #104]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #104]
	ldr	d0, [x2, #16]
	ldp	d1, d2, [x1, #96]
	fmul	d0, d1, d0
	ldr	d1, [x2, #48]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #80]
	ldp	d2, d3, [x1, #112]
	fmul	d1, d2, d1
	ldr	d2, [x2, #112]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #112]
	ldr	d0, [x2, #24]
	ldp	d1, d2, [x1, #96]
	fmul	d0, d1, d0
	ldr	d1, [x2, #56]
	fmul	d1, d2, d1
	fadd	d0, d0, d1
	ldr	d1, [x2, #88]
	ldp	d2, d3, [x1, #112]
	fmul	d1, d2, d1
	ldr	d2, [x2, #120]
	fmul	d2, d3, d2
	fadd	d1, d1, d2
	fadd	d0, d0, d1
	str	d0, [x0, #120]
	ret
	.cfi_endproc
                                        ; -- End function
	.globl	__Z9check_errPKdS0_             ; -- Begin function _Z9check_errPKdS0_
	.p2align	2
__Z9check_errPKdS0_:                    ; @_Z9check_errPKdS0_
	.cfi_startproc
; %bb.0:
	ldp	d0, d1, [x0]
	ldp	d2, d3, [x1]
	fabd	d0, d0, d2
	movi	d2, #0000000000000000
	fmaxnm	d0, d0, d2
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #16]
	ldp	d3, d4, [x1, #16]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #32]
	ldp	d3, d4, [x1, #32]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #48]
	ldp	d3, d4, [x1, #48]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #64]
	ldp	d3, d4, [x1, #64]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #80]
	ldp	d3, d4, [x1, #80]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #96]
	ldp	d3, d4, [x1, #96]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ldp	d1, d2, [x0, #112]
	ldp	d3, d4, [x1, #112]
	fabd	d1, d1, d3
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	fabd	d1, d2, d4
	fcmp	d1, d0
	fcsel	d0, d1, d0, gt
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
