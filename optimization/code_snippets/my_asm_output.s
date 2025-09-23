	.arch armv8-a
; GNU C++17 (Homebrew GCC 12.2.0) version 12.2.0 (aarch64-apple-darwin21)
;	compiled by GNU C version 12.2.0, GMP version 6.2.1, MPFR version 4.1.0, MPC version 1.2.1, isl version isl-0.25-GMP

; warning: MPFR header version 4.1.0 differs from library version 4.1.0-p13.
; GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
; options passed: -fPIC -mmacosx-version-min=13.0.0 -mlittle-endian -mabi=lp64 -O3
	.text
	.cstring
	.align	3
lC1:
	.ascii "4 loads = %g, load 4 = %g\12\0"
	.section __TEXT,__text_startup,regular,pure_instructions
	.align	2
	.p2align 4,,11
	.globl _main
_main:
LFB4355:
	sub	sp, sp, #96	;,,
LCFI0:
; product_sum.cc:64:     const double D[8] = {1.,2.,3.,4.,5.,6.,7.,8.};
	adrp	x2, lC0@PAGE	; tmp107,
	add	x2, x2, lC0@PAGEOFF;momd	; tmp106, tmp107,
	add	x1, sp, 32	; tmp108,,
; product_sum.cc:83:     printf ("4 loads = %g, load 4 = %g\n", _neon_sum_2 (P[2]), _neon_sum_2 ( P4.val[2]));
	mov	x0, 4630263366890291200	; tmp124,
	fmov	d0, x0	; tmp122, tmp124
; product_sum.cc:63: int main (int argc, const char** argv) {
	stp	x29, x30, [sp, 16]	;,,
LCFI1:
	add	x29, sp, 16	;,,
; product_sum.cc:83:     printf ("4 loads = %g, load 4 = %g\n", _neon_sum_2 (P[2]), _neon_sum_2 ( P4.val[2]));
	adrp	x0, lC1@PAGE	; tmp118,
; product_sum.cc:63: int main (int argc, const char** argv) {
; product_sum.cc:64:     const double D[8] = {1.,2.,3.,4.,5.,6.,7.,8.};
	ldp	q3, q4, [x2]	; tmp110, tmp111,
; product_sum.cc:83:     printf ("4 loads = %g, load 4 = %g\n", _neon_sum_2 (P[2]), _neon_sum_2 ( P4.val[2]));
	str	d0, [sp]	; tmp122,
; product_sum.cc:64:     const double D[8] = {1.,2.,3.,4.,5.,6.,7.,8.};
	ldp	q1, q2, [x2, 32]	; tmp112, tmp113,
; product_sum.cc:83:     printf ("4 loads = %g, load 4 = %g\n", _neon_sum_2 (P[2]), _neon_sum_2 ( P4.val[2]));
	add	x0, x0, lC1@PAGEOFF;momd	;, tmp118,
; product_sum.cc:64:     const double D[8] = {1.,2.,3.,4.,5.,6.,7.,8.};
	stp	q3, q4, [x1]	; tmp110, tmp111, D
	stp	q1, q2, [x1, 32]	; tmp112, tmp113, D
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:16793:   return __builtin_aarch64_ld4v2df ((const __builtin_aarch64_simd_df *) __a);
	ld4	{v4.2d - v7.2d}, [x1]	; D.26810,
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:345:   return __a + __b;
	fadd	v1.2d, v7.2d, v6.2d	; tmp116, D.26810, D.26810
	fadd	v0.2d, v4.2d, v5.2d	; tmp115, D.26810, D.26810
	fadd	v0.2d, v0.2d, v1.2d	; _11, tmp115, tmp116
; /opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/current/gcc/aarch64-apple-darwin21/12/include/arm_neon.h:2750:   return __aarch64_vget_lane_any (__a, __b);
	dup	d1, v0.d[1]	; tmp120, _11,
; product_sum.cc:58:     return vgetq_lane_f64 (x,0) + vgetq_lane_f64 (x,1);
	fadd	d0, d1, d0	; tmp121, tmp120, _11
; product_sum.cc:83:     printf ("4 loads = %g, load 4 = %g\n", _neon_sum_2 (P[2]), _neon_sum_2 ( P4.val[2]));
	str	d0, [sp, 8]	; tmp121,
	bl	_printf		;
; product_sum.cc:86: }
	ldp	x29, x30, [sp, 16]	;,,
	mov	w0, 0	;,
	add	sp, sp, 96	;,,
LCFI2:
	ret	
LFE4355:
	.const
	.align	3
lC0:
	.word	0
	.word	1072693248
	.word	0
	.word	1073741824
	.word	0
	.word	1074266112
	.word	0
	.word	1074790400
	.word	0
	.word	1075052544
	.word	0
	.word	1075314688
	.word	0
	.word	1075576832
	.word	0
	.word	1075838976
	.section __TEXT,__text_startup,regular,pure_instructions
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
	.quad	LFB4355-.
	.set L$set$2,LFE4355-LFB4355
	.quad L$set$2
	.uleb128 0
	.byte	0x4
	.set L$set$3,LCFI0-LFB4355
	.long L$set$3
	.byte	0xe
	.uleb128 0x60
	.byte	0x4
	.set L$set$4,LCFI1-LCFI0
	.long L$set$4
	.byte	0x9d
	.uleb128 0xa
	.byte	0x9e
	.uleb128 0x9
	.byte	0x4
	.set L$set$5,LCFI2-LCFI1
	.long L$set$5
	.byte	0xdd
	.byte	0xde
	.byte	0xe
	.uleb128 0
	.align	3
LEFDE1:
	.ident	"GCC: (Homebrew GCC 12.2.0) 12.2.0"
	.subsections_via_symbols
