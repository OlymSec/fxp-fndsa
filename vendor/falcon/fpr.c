/*
 * Floating-point operations.
 *
 * This file implements the non-inline functions declared in
 * fpr.h, as well as the constants for FFT / iFFT.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include "inner.h"

#if FALCON_FPEMU // yyyFPEMU+1

/*
 * Normalize a provided unsigned integer to the 2^63..2^64-1 range by
 * left-shifting it if necessary. The exponent e is adjusted accordingly
 * (i.e. if the value was left-shifted by n bits, then n is subtracted
 * from e). If source m is 0, then it remains 0, but e is altered.
 * Both m and e must be simple variables (no expressions allowed).
 */
#define FPR_NORM64(m, e)   do { \
		uint32_t nt; \
 \
		(e) -= 63; \
 \
		nt = (uint32_t)((m) >> 32); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) << 32)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 5); \
 \
		nt = (uint32_t)((m) >> 48); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) << 16)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 4); \
 \
		nt = (uint32_t)((m) >> 56); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  8)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 3); \
 \
		nt = (uint32_t)((m) >> 60); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  4)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 2); \
 \
		nt = (uint32_t)((m) >> 62); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  2)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 1); \
 \
		nt = (uint32_t)((m) >> 63); \
		(m) ^= ((m) ^ ((m) <<  1)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt); \
	} while (0)

#if FALCON_ASM_CORTEXM4 // yyyASM_CORTEXM4+1

__attribute__((naked))
fpr
fpr_scaled(int64_t i __attribute__((unused)), int sc __attribute__((unused)))
{
	__asm__ (
	"push	{ r4, r5, r6, lr }\n\t"
	"\n\t"
	"@ Input i is in r0:r1, and sc in r2.\n\t"
	"@ Extract the sign bit, and compute the absolute value.\n\t"
	"@ -> sign bit in r3, with value 0 or -1\n\t"
	"asrs	r3, r1, #31\n\t"
	"eors	r0, r3\n\t"
	"eors	r1, r3\n\t"
	"subs	r0, r3\n\t"
	"sbcs	r1, r3\n\t"
	"\n\t"
	"@ Scale exponent to account for the encoding; if the source is\n\t"
	"@ zero or if the scaled exponent is negative, it is set to 32.\n\t"
	"addw	r2, r2, #1022\n\t"
	"orrs	r4, r0, r1\n\t"
	"bics	r4, r4, r2, asr #31\n\t"
	"rsbs	r5, r4, #0\n\t"
	"orrs	r4, r5\n\t"
	"ands	r2, r2, r4, asr #31\n\t"
	"adds	r2, #32\n\t"
	"\n\t"
	"@ Normalize value to a full 64-bit width, by shifting it left.\n\t"
	"@ The shift count is subtracted from the exponent (in r2).\n\t"
	"@ If the mantissa is 0, the exponent is set to 0.\n\t"
	"\n\t"
	"@ If top word is 0, replace with low word; otherwise, add 32 to\n\t"
	"@ the exponent.\n\t"
	"rsbs	r4, r1, #0\n\t"
	"orrs	r4, r1\n\t"
	"eors	r5, r0, r1\n\t"
	"bics	r5, r5, r4, asr #31\n\t"
	"eors	r1, r5\n\t"
	"ands	r0, r0, r4, asr #31\n\t"
	"lsrs	r4, r4, #31\n\t"
	"adds	r2, r2, r4, lsl #5\n\t"
	"\n\t"
	"@ Count leading zeros of r1 to finish the shift.\n\t"
	"clz	r4, r1\n\t"
	"subs	r2, r4\n\t"
	"rsbs	r5, r4, #32\n\t"
	"lsls	r1, r4\n\t"
	"lsrs	r5, r0, r5\n\t"
	"lsls	r0, r4\n\t"
	"orrs	r1, r5\n\t"
	"\n\t"
	"@ Clear the top bit; we know it's a 1 (unless the whole mantissa\n\t"
	"@ was zero, but then it's still OK to clear it)\n\t"
	"bfc	r1, #31, #1\n\t"
	"\n\t"
	"@ Now shift right the value by 11 bits; this puts the value in\n\t"
	"@ the 2^52..2^53-1 range. We also keep a copy of the pre-shift\n\t"
	"@ low bits in r5.\n\t"
	"movs	r5, r0\n\t"
	"lsrs	r0, #11\n\t"
	"orrs	r0, r0, r1, lsl #21\n\t"
	"lsrs	r1, #11\n\t"
	"\n\t"
	"@ Also plug the exponent at the right place. This must be done\n\t"
	"@ now so that, in case the rounding creates a carry, that carry\n\t"
	"@ adds to the exponent, which would be exactly what we want at\n\t"
	"@ that point.\n\t"
	"orrs	r1, r1, r2, lsl #20\n\t"
	"\n\t"
	"@ Rounding: we must add 1 to the mantissa in the following cases:\n\t"
	"@  - bits 11 to 9 of r5 are '011', '110' or '111'\n\t"
	"@  - bits 11 to 9 of r5 are '010' and one of the\n\t"
	"@    bits 0 to 8 is non-zero\n\t"
	"ubfx	r6, r5, #0, #9\n\t"
	"addw	r6, r6, #511\n\t"
	"orrs	r5, r6\n\t"
	"\n\t"
	"ubfx	r5, r5, #9, #3\n\t"
	"movs	r6, #0xC8\n\t"
	"lsrs	r6, r5\n\t"
	"ands	r6, #1\n\t"
	"adds	r0, r6\n\t"
	"adcs	r1, #0\n\t"
	"\n\t"
	"@ Put back the sign.\n\t"
	"orrs	r1, r1, r3, lsl #31\n\t"
	"\n\t"
	"pop	{ r4, r5, r6, pc}\n\t"
	);
}

#else // yyyASM_CORTEXM4+0

fpr
fpr_scaled(int64_t i, int sc)
{
	/*
	 * To convert from int to float, we have to do the following:
	 *  1. Get the absolute value of the input, and its sign
	 *  2. Shift right or left the value as appropriate
	 *  3. Pack the result
	 *
	 * We can assume that the source integer is not -2^63.
	 */
	int s, e;
	uint32_t t;
	uint64_t m;

	/*
	 * Extract sign bit.
	 * We have: -i = 1 + ~i
	 */
	s = (int)((uint64_t)i >> 63);
	i ^= -(int64_t)s;
	i += s;

	/*
	 * For now we suppose that i != 0.
	 * Otherwise, we set m to i and left-shift it as much as needed
	 * to get a 1 in the top bit. We can do that in a logarithmic
	 * number of conditional shifts.
	 */
	m = (uint64_t)i;
	e = 9 + sc;
	FPR_NORM64(m, e);

	/*
	 * Now m is in the 2^63..2^64-1 range. We must divide it by 512;
	 * if one of the dropped bits is a 1, this should go into the
	 * "sticky bit".
	 */
	m |= ((uint32_t)m & 0x1FF) + 0x1FF;
	m >>= 9;

	/*
	 * Corrective action: if i = 0 then all of the above was
	 * incorrect, and we clamp e and m down to zero.
	 */
	t = (uint32_t)((uint64_t)(i | -i) >> 63);
	m &= -(uint64_t)t;
	e &= -(int)t;

	/*
	 * Assemble back everything. The FPR() function will handle cases
	 * where e is too low.
	 */
	return FPR(s, e, m);
}

#endif // yyyASM_CORTEXM4-

#if FALCON_ASM_CORTEXM4 // yyyASM_CORTEXM4+1

// yyyPQCLEAN+0
#if 0
/* Debug code -- To get a printout of registers from a specific point
   in ARM Cortex M4 assembly code, uncomment this code and add a
   "bl DEBUG" call where wished for. */

void
print_regs(uint32_t *rr, uint32_t flags)
{
	int i;
	extern int printf(const char *fmt, ...);

	printf("\nRegs:\n");
	for (i = 0; i < 7; i ++) {
		int j;

		j = i + 7;
		printf("  %2d = %08X    %2d = %08X\n", i, rr[i], j, rr[j]);
	}
	printf("  flags = %08X  ", flags);
	if ((flags >> 31) & 1) {
		printf("N");
	}
	if ((flags >> 30) & 1) {
		printf("Z");
	}
	if ((flags >> 29) & 1) {
		printf("C");
	}
	if ((flags >> 28) & 1) {
		printf("V");
	}
	if ((flags >> 27) & 1) {
		printf("Q");
	}
	printf("\n");
}

__attribute__((naked))
void
DEBUG(void)
{
	__asm__ (
	"push	{ r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, lr }\n\t"
	"mov	r0, sp\n\t"
	"mrs	r1, apsr\n\t"
	"bl	print_regs\n\t"
	"pop	{ r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, pc }\n\t"
	);
}
#endif
// yyyPQCLEAN-

__attribute__((naked))
fpr
fpr_add(fpr x __attribute__((unused)), fpr y __attribute__((unused)))
{
	__asm__ (
	"push	{ r4, r5, r6, r7, r8, r10, r11, lr }\n\t"
	"\n\t"
	"@ Make sure that the first operand (x) has the larger absolute\n\t"
	"@ value. This guarantees that the exponent of y is less than\n\t"
	"@ or equal to the exponent of x, and, if they are equal, then\n\t"
	"@ the mantissa of y will not be greater than the mantissa of x.\n\t"
	"@ However, if absolute values are equal and the sign of x is 1,\n\t"
	"@ then we want to also swap the values.\n\t"
	"ubfx	r4, r1, #0, #31  @ top word without sign bit\n\t"
	"ubfx	r5, r3, #0, #31  @ top word without sign bit\n\t"
	"subs	r7, r0, r2       @ difference in r7:r4\n\t"
	"sbcs	r4, r5\n\t"
	"orrs	r7, r4\n\t"
	"rsbs	r5, r7, #0\n\t"
	"orrs	r7, r5      @ bit 31 of r7 is 0 iff difference is zero\n\t"
	"bics	r6, r1, r7\n\t"
	"orrs	r6, r4      @ bit 31 of r6 is 1 iff the swap must be done\n\t"
	"\n\t"
	"@ Conditional swap\n\t"
	"eors	r4, r0, r2\n\t"
	"eors	r5, r1, r3\n\t"
	"ands	r4, r4, r6, asr #31\n\t"
	"ands	r5, r5, r6, asr #31\n\t"
	"eors	r0, r4\n\t"
	"eors	r1, r5\n\t"
	"eors	r2, r4\n\t"
	"eors	r3, r5\n\t"
	"\n\t"
	"@ Extract mantissa of x into r0:r1, exponent in r4, sign in r5\n\t"
	"ubfx	r4, r1, #20, #11   @ Exponent in r4 (without sign)\n\t"
	"addw	r5, r4, #2047 @ Get a carry to test r4 for zero\n\t"
	"lsrs	r5, #11       @ r5 is the mantissa implicit high bit\n\t"
	"bfc	r1, #20, #11  @ Clear exponent bits (not the sign)\n\t"
	"orrs	r1, r1, r5, lsl #20  @ Set mantissa high bit\n\t"
	"asrs	r5, r1, #31   @ Get sign bit (sign-extended)\n\t"
	"bfc	r1, #31, #1   @ Clear the sign bit\n\t"
	"\n\t"
	"@ Extract mantissa of y into r2:r3, exponent in r6, sign in r7\n\t"
	"ubfx	r6, r3, #20, #11   @ Exponent in r6 (without sign)\n\t"
	"addw	r7, r6, #2047 @ Get a carry to test r6 for zero\n\t"
	"lsrs	r7, #11       @ r7 is the mantissa implicit high bit\n\t"
	"bfc	r3, #20, #11  @ Clear exponent bits (not the sign)\n\t"
	"orrs	r3, r3, r7, lsl #20  @ Set mantissa high bit\n\t"
	"asrs	r7, r3, #31   @ Get sign bit (sign-extended)\n\t"
	"bfc	r3, #31, #1   @ Clear the sign bit\n\t"
	"\n\t"
	"@ Scale mantissas up by three bits.\n\t"
	"lsls	r1, #3\n\t"
	"orrs	r1, r1, r0, lsr #29\n\t"
	"lsls	r0, #3\n\t"
	"lsls	r3, #3\n\t"
	"orrs	r3, r3, r2, lsr #29\n\t"
	"lsls	r2, #3\n\t"
	"\n\t"
	"@ x: exponent=r4, sign=r5, mantissa=r0:r1 (scaled up 3 bits)\n\t"
	"@ y: exponent=r6, sign=r7, mantissa=r2:r3 (scaled up 3 bits)\n\t"
	"\n\t"
	"@ At that point, the exponent of x (in r4) is larger than that\n\t"
	"@ of y (in r6). The difference is the amount of shifting that\n\t"
	"@ should be done on y. If that amount is larger than 59 then\n\t"
	"@ we clamp y to 0. We won't need y's exponent beyond that point,\n\t"
	"@ so we store that shift count in r6.\n\t"
	"subs	r6, r4, r6\n\t"
	"subs	r8, r6, #60\n\t"
	"ands	r2, r2, r8, asr #31\n\t"
	"ands	r3, r3, r8, asr #31\n\t"
	"\n\t"
	"@ Shift right r2:r3 by r6 bits. The shift count is in the 0..59\n\t"
	"@ range. r11 will be non-zero if and only if some non-zero bits\n\t"
	"@ were dropped.\n\t"
	"subs	r8, r6, #32\n\t"
	"bics	r11, r2, r8, asr #31\n\t"
	"ands	r2, r2, r8, asr #31\n\t"
	"bics	r10, r3, r8, asr #31\n\t"
	"orrs	r2, r2, r10\n\t"
	"ands	r3, r3, r8, asr #31\n\t"
	"ands	r6, r6, #31\n\t"
	"rsbs	r8, r6, #32\n\t"
	"lsls	r10, r2, r8\n\t"
	"orrs	r11, r11, r10\n\t"
	"lsrs	r2, r2, r6\n\t"
	"lsls	r10, r3, r8\n\t"
	"orrs	r2, r2, r10\n\t"
	"lsrs	r3, r3, r6\n\t"
	"\n\t"
	"@ If r11 is non-zero then some non-zero bit was dropped and the\n\t"
	"@ low bit of r2 must be forced to 1 ('sticky bit').\n\t"
	"rsbs	r6, r11, #0\n\t"
	"orrs	r6, r6, r11\n\t"
	"orrs	r2, r2, r6, lsr #31\n\t"
	"\n\t"
	"@ x: exponent=r4, sign=r5, mantissa=r0:r1 (scaled up 3 bits)\n\t"
	"@ y: sign=r7, value=r2:r3 (scaled to same exponent as x)\n\t"
	"\n\t"
	"@ If x and y don't have the same sign, then we should negate r2:r3\n\t"
	"@ (i.e. subtract the mantissa instead of adding it). Signs of x\n\t"
	"@ and y are in r5 and r7, as full-width words. We won't need r7\n\t"
	"@ afterwards.\n\t"
	"eors	r7, r5    @ r7 = -1 if y must be negated, 0 otherwise\n\t"
	"eors	r2, r7\n\t"
	"eors	r3, r7\n\t"
	"subs	r2, r7\n\t"
	"sbcs	r3, r7\n\t"
	"\n\t"
	"@ r2:r3 has been shifted, we can add to r0:r1.\n\t"
	"adds	r0, r2\n\t"
	"adcs	r1, r3\n\t"
	"\n\t"
	"@ result: exponent=r4, sign=r5, mantissa=r0:r1 (scaled up 3 bits)\n\t"
	"\n\t"
	"@ Normalize the result with some left-shifting to full 64-bit\n\t"
	"@ width. Shift count goes to r2, and exponent (r4) is adjusted.\n\t"
	"clz	r2, r0\n\t"
	"clz	r3, r1\n\t"
	"sbfx	r6, r3, #5, #1\n\t"
	"ands	r2, r6\n\t"
	"adds	r2, r2, r3\n\t"
	"subs	r4, r4, r2\n\t"
	"\n\t"
	"@ Shift r0:r1 to the left by r2 bits.\n\t"
	"subs	r7, r2, #32\n\t"
	"lsls	r7, r0, r7\n\t"
	"lsls	r1, r1, r2\n\t"
	"rsbs	r6, r2, #32\n\t"
	"orrs	r1, r1, r7\n\t"
	"lsrs	r6, r0, r6\n\t"
	"orrs	r1, r1, r6\n\t"
	"lsls	r0, r0, r2\n\t"
	"\n\t"
	"@ The exponent of x was in r4. The left-shift operation has\n\t"
	"@ subtracted some value from it, 8 in case the result has the\n\t"
	"@ same exponent as x. However, the high bit of the mantissa will\n\t"
	"@ add 1 to the exponent, so we only add back 7 (the exponent is\n\t"
	"@ added in because rounding might have produced a carry, which\n\t"
	"@ should then spill into the exponent).\n\t"
	"adds	r4, #7\n\t"
	"\n\t"
	"@ If the mantissa new mantissa is non-zero, then its bit 63 is\n\t"
	"@ non-zero (thanks to the normalizing shift). Otherwise, that bit\n\t"
	"@ is zero, and we should then set the exponent to zero as well.\n\t"
	"ands	r4, r4, r1, asr #31\n\t"
	"\n\t"
	"@ Shrink back the value to a 52-bit mantissa. This requires\n\t"
	"@ right-shifting by 11 bits; we keep a copy of the pre-shift\n\t"
	"@ low word in r3.\n\t"
	"movs	r3, r0\n\t"
	"lsrs	r0, #11\n\t"
	"orrs	r0, r0, r1, lsl #21\n\t"
	"lsrs	r1, #11\n\t"
	"\n\t"
	"@ Apply rounding.\n\t"
	"ubfx	r6, r3, #0, #9\n\t"
	"addw	r6, r6, #511\n\t"
	"orrs	r3, r6\n\t"
	"ubfx	r3, r3, #9, #3\n\t"
	"movs	r6, #0xC8\n\t"
	"lsrs	r6, r3\n\t"
	"ands	r6, #1\n\t"
	"adds	r0, r6\n\t"
	"adcs	r1, #0\n\t"
	"\n\t"
	"@Plug in the exponent with an addition.\n\t"
	"adds	r1, r1, r4, lsl #20\n\t"
	"\n\t"
	"@ If the new exponent is negative or zero, then it underflowed\n\t"
	"@ and we must clear the whole mantissa and exponent.\n\t"
	"rsbs	r4, r4, #0\n\t"
	"ands	r0, r0, r4, asr #31\n\t"
	"ands	r1, r1, r4, asr #31\n\t"
	"\n\t"
	"@ Put back the sign. This is the sign of x: thanks to the\n\t"
	"@ conditional swap at the start, this is always correct.\n\t"
	"bfi	r1, r5, #31, #1\n\t"
	"\n\t"
	"pop	{ r4, r5, r6, r7, r8, r10, r11, pc }\n\t"
	);
}

#else // yyyASM_CORTEXM4+0

fpr
fpr_add(fpr x, fpr y)
{
	uint64_t m, xu, yu, za;
	uint32_t cs;
	int ex, ey, sx, sy, cc;

	/*
	 * Make sure that the first operand (x) has the larger absolute
	 * value. This guarantees that the exponent of y is less than
	 * or equal to the exponent of x, and, if they are equal, then
	 * the mantissa of y will not be greater than the mantissa of x.
	 *
	 * After this swap, the result will have the sign x, except in
	 * the following edge case: abs(x) = abs(y), and x and y have
	 * opposite sign bits; in that case, the result shall be +0
	 * even if the sign bit of x is 1. To handle this case properly,
	 * we do the swap is abs(x) = abs(y) AND the sign of x is 1.
	 */
	m = ((uint64_t)1 << 63) - 1;
	za = (x & m) - (y & m);
	cs = (uint32_t)(za >> 63)
		| ((1U - (uint32_t)(-za >> 63)) & (uint32_t)(x >> 63));
	m = (x ^ y) & -(uint64_t)cs;
	x ^= m;
	y ^= m;

	/*
	 * Extract sign bits, exponents and mantissas. The mantissas are
	 * scaled up to 2^55..2^56-1, and the exponent is unbiased. If
	 * an operand is zero, its mantissa is set to 0 at this step, and
	 * its exponent will be -1078.
	 */
	ex = (int)(x >> 52);
	sx = ex >> 11;
	ex &= 0x7FF;
	m = (uint64_t)(uint32_t)((ex + 0x7FF) >> 11) << 52;
	xu = ((x & (((uint64_t)1 << 52) - 1)) | m) << 3;
	ex -= 1078;
	ey = (int)(y >> 52);
	sy = ey >> 11;
	ey &= 0x7FF;
	m = (uint64_t)(uint32_t)((ey + 0x7FF) >> 11) << 52;
	yu = ((y & (((uint64_t)1 << 52) - 1)) | m) << 3;
	ey -= 1078;

	/*
	 * x has the larger exponent; hence, we only need to right-shift y.
	 * If the shift count is larger than 59 bits then we clamp the
	 * value to zero.
	 */
	cc = ex - ey;
	yu &= -(uint64_t)((uint32_t)(cc - 60) >> 31);
	cc &= 63;

	/*
	 * The lowest bit of yu is "sticky".
	 */
	m = fpr_ulsh(1, cc) - 1;
	yu |= (yu & m) + m;
	yu = fpr_ursh(yu, cc);

	/*
	 * If the operands have the same sign, then we add the mantissas;
	 * otherwise, we subtract the mantissas.
	 */
	xu += yu - ((yu << 1) & -(uint64_t)(sx ^ sy));

	/*
	 * The result may be smaller, or slightly larger. We normalize
	 * it to the 2^63..2^64-1 range (if xu is zero, then it stays
	 * at zero).
	 */
	FPR_NORM64(xu, ex);

	/*
	 * Scale down the value to 2^54..s^55-1, handling the last bit
	 * as sticky.
	 */
	xu |= ((uint32_t)xu & 0x1FF) + 0x1FF;
	xu >>= 9;
	ex += 9;

	/*
	 * In general, the result has the sign of x. However, if the
	 * result is exactly zero, then the following situations may
	 * be encountered:
	 *   x > 0, y = -x   -> result should be +0
	 *   x < 0, y = -x   -> result should be +0
	 *   x = +0, y = +0  -> result should be +0
	 *   x = -0, y = +0  -> result should be +0
	 *   x = +0, y = -0  -> result should be +0
	 *   x = -0, y = -0  -> result should be -0
	 *
	 * But at the conditional swap step at the start of the
	 * function, we ensured that if abs(x) = abs(y) and the
	 * sign of x was 1, then x and y were swapped. Thus, the
	 * two following cases cannot actually happen:
	 *   x < 0, y = -x
	 *   x = -0, y = +0
	 * In all other cases, the sign bit of x is conserved, which
	 * is what the FPR() function does. The FPR() function also
	 * properly clamps values to zero when the exponent is too
	 * low, but does not alter the sign in that case.
	 */
	return FPR(sx, ex, xu);
}

#endif // yyyASM_CORTEXM4-

#if FALCON_ASM_CORTEXM4 // yyyASM_CORTEXM4+1

__attribute__((naked))
fpr
fpr_mul(fpr x __attribute__((unused)), fpr y __attribute__((unused)))
{
	__asm__ (
	"push	{ r4, r5, r6, r7, r8, r10, r11, lr }\n\t"
	"\n\t"
	"@ Extract mantissas: x.m = r4:r5, y.m = r6:r7\n\t"
	"@ r4 and r6 contain only 25 bits each.\n\t"
	"bics	r4, r0, #0xFE000000\n\t"
	"lsls	r5, r1, #7\n\t"
	"orrs	r5, r5, r0, lsr #25\n\t"
	"orrs	r5, r5, #0x08000000\n\t"
	"bics	r5, r5, #0xF0000000\n\t"
	"bics	r6, r2, #0xFE000000\n\t"
	"lsls	r7, r3, #7\n\t"
	"orrs	r7, r7, r2, lsr #25\n\t"
	"orrs	r7, r7, #0x08000000\n\t"
	"bics	r7, r7, #0xF0000000\n\t"
	"\n\t"
	"@ Perform product. Values are in the 2^52..2^53-1 range, so\n\t"
	"@ the product is at most 106-bit long. Of the low 50 bits,\n\t"
	"@ we only want to know if they are all zeros or not. Here,\n\t"
	"@ we get the top 56 bits in r10:r11, and r8 will be non-zero\n\t"
	"@ if and only if at least one of the low 50 bits is non-zero.\n\t"
	"umull	r8, r10, r4, r6      @ x0*y0\n\t"
	"lsls	r10, #7\n\t"
	"orrs	r10, r10, r8, lsr #25\n\t"
	"eors	r11, r11\n\t"
	"umlal	r10, r11, r4, r7     @ x0*y1\n\t"
	"umlal	r10, r11, r5, r6     @ x1*y0\n\t"
	"orrs	r8, r8, r10, lsl #7\n\t"
	"lsrs	r10, #25\n\t"
	"orrs	r10, r10, r11, lsl #7\n\t"
	"eors	r11, r11\n\t"
	"umlal	r10, r11, r5, r7     @ x1*y1\n\t"
	"\n\t"
	"@ Now r0, r2, r4, r5, r6 and r7 are free.\n\t"
	"@ If any of the low 50 bits was non-zero, then we force the\n\t"
	"@ low bit of r10 to 1.\n\t"
	"rsbs	r4, r8, #0\n\t"
	"orrs	r8, r8, r4\n\t"
	"orrs	r10, r10, r8, lsr #31\n\t"
	"\n\t"
	"@ r8 is free.\n\t"
	"@ r10:r11 contains the product in the 2^54..2^56-1 range. We\n\t"
	"@ normalize it to 2^54..2^55-1 (into r6:r7) with a conditional\n\t"
	"@ shift (low bit is sticky). r5 contains -1 if the shift was done,\n\t"
	"@ 0 otherwise.\n\t"
	"ands	r6, r10, #1\n\t"
	"lsrs	r5, r11, #23\n\t"
	"rsbs	r5, r5, #0\n\t"
	"orrs	r6, r6, r10, lsr #1\n\t"
	"orrs	r6, r6, r11, lsl #31\n\t"
	"lsrs	r7, r11, #1\n\t"
	"eors	r10, r10, r6\n\t"
	"eors	r11, r11, r7\n\t"
	"bics	r10, r10, r5\n\t"
	"bics	r11, r11, r5\n\t"
	"eors	r6, r6, r10\n\t"
	"eors	r7, r7, r11\n\t"
	"\n\t"
	"@ Compute aggregate exponent: ex + ey - 1023 + w\n\t"
	"@ (where w = 1 if the conditional shift was done, 0 otherwise)\n\t"
	"@ But we subtract 1 because the injection of the mantissa high\n\t"
	"@ bit will increment the exponent by 1.\n\t"
	"lsls	r0, r1, #1\n\t"
	"lsls	r2, r3, #1\n\t"
	"lsrs	r0, #21\n\t"
	"addw	r4, r0, #0x7FF   @ save ex + 2047 in r4\n\t"
	"lsrs	r2, #21\n\t"
	"addw	r8, r2, #0x7FF   @ save ey + 2047 in r8\n\t"
	"adds	r2, r0\n\t"
	"subw	r2, r2, #1024\n\t"
	"subs	r2, r5\n\t"
	"\n\t"
	"@ r5 is free.\n\t"
	"@ Also, if either of the source exponents is 0, or the result\n\t"
	"@ exponent is 0 or negative, then the result is zero and the\n\t"
	"@ mantissa and the exponent shall be clamped to zero. Since\n\t"
	"@ r2 contains the result exponent minus 1, we test on r2\n\t"
	"@ being strictly negative.\n\t"
	"ands	r4, r8    @ if bit 11 = 0 then one of the exponents was 0\n\t"
	"mvns	r5, r2\n\t"
	"ands	r5, r5, r4, lsl #20\n\t"
	"ands	r2, r2, r5, asr #31\n\t"
	"ands	r6, r6, r5, asr #31\n\t"
	"ands	r7, r7, r5, asr #31\n\t"
	"\n\t"
	"@ Sign is the XOR of the sign of the operands. This is true in\n\t"
	"@ all cases, including very small results (exponent underflow)\n\t"
	"@ and zeros.\n\t"
	"eors	r1, r3\n\t"
	"bfc	r1, #0, #31\n\t"
	"\n\t"
	"@ Plug in the exponent.\n\t"
	"bfi	r1, r2, #20, #11\n\t"
	"\n\t"
	"@ r2 and r3 are free.\n\t"
	"@ Shift back to the normal 53-bit mantissa, with rounding.\n\t"
	"@ Mantissa goes into r0:r1. For r1, we must use an addition\n\t"
	"@ because the rounding may have triggered a carry, that should\n\t"
	"@ be added to the exponent.\n\t"
	"movs	r4, r6\n\t"
	"lsrs	r0, r6, #2\n\t"
	"orrs	r0, r0, r7, lsl #30\n\t"
	"adds	r1, r1, r7, lsr #2\n\t"
	"ands	r4, #0x7\n\t"
	"movs	r3, #0xC8\n\t"
	"lsrs	r3, r4\n\t"
	"ands	r3, #1\n\t"
	"adds	r0, r3\n\t"
	"adcs	r1, #0\n\t"
	"\n\t"
	"pop	{ r4, r5, r6, r7, r8, r10, r11, pc }\n\t"
	);
}

#else // yyyASM_CORTEXM4+0

fpr
fpr_mul(fpr x, fpr y)
{
	uint64_t xu, yu, w, zu, zv;
	uint32_t x0, x1, y0, y1, z0, z1, z2;
	int ex, ey, d, e, s;

	/*
	 * Extract absolute values as scaled unsigned integers. We
	 * don't extract exponents yet.
	 */
	xu = (x & (((uint64_t)1 << 52) - 1)) | ((uint64_t)1 << 52);
	yu = (y & (((uint64_t)1 << 52) - 1)) | ((uint64_t)1 << 52);

	/*
	 * We have two 53-bit integers to multiply; we need to split
	 * each into a lower half and a upper half. Moreover, we
	 * prefer to have lower halves to be of 25 bits each, for
	 * reasons explained later on.
	 */
	x0 = (uint32_t)xu & 0x01FFFFFF;
	x1 = (uint32_t)(xu >> 25);
	y0 = (uint32_t)yu & 0x01FFFFFF;
	y1 = (uint32_t)(yu >> 25);
	w = (uint64_t)x0 * (uint64_t)y0;
	z0 = (uint32_t)w & 0x01FFFFFF;
	z1 = (uint32_t)(w >> 25);
	w = (uint64_t)x0 * (uint64_t)y1;
	z1 += (uint32_t)w & 0x01FFFFFF;
	z2 = (uint32_t)(w >> 25);
	w = (uint64_t)x1 * (uint64_t)y0;
	z1 += (uint32_t)w & 0x01FFFFFF;
	z2 += (uint32_t)(w >> 25);
	zu = (uint64_t)x1 * (uint64_t)y1;
	z2 += (z1 >> 25);
	z1 &= 0x01FFFFFF;
	zu += z2;

	/*
	 * Since xu and yu are both in the 2^52..2^53-1 range, the
	 * product is in the 2^104..2^106-1 range. We first reassemble
	 * it and round it into the 2^54..2^56-1 range; the bottom bit
	 * is made "sticky". Since the low limbs z0 and z1 are 25 bits
	 * each, we just take the upper part (zu), and consider z0 and
	 * z1 only for purposes of stickiness.
	 * (This is the reason why we chose 25-bit limbs above.)
	 */
	zu |= ((z0 | z1) + 0x01FFFFFF) >> 25;

	/*
	 * We normalize zu to the 2^54..s^55-1 range: it could be one
	 * bit too large at this point. This is done with a conditional
	 * right-shift that takes into account the sticky bit.
	 */
	zv = (zu >> 1) | (zu & 1);
	w = zu >> 55;
	zu ^= (zu ^ zv) & -w;

	/*
	 * Get the aggregate scaling factor:
	 *
	 *   - Each exponent is biased by 1023.
	 *
	 *   - Integral mantissas are scaled by 2^52, hence an
	 *     extra 52 bias for each exponent.
	 *
	 *   - However, we right-shifted z by 50 bits, and then
	 *     by 0 or 1 extra bit (depending on the value of w).
	 *
	 * In total, we must add the exponents, then subtract
	 * 2 * (1023 + 52), then add 50 + w.
	 */
	ex = (int)((x >> 52) & 0x7FF);
	ey = (int)((y >> 52) & 0x7FF);
	e = ex + ey - 2100 + (int)w;

	/*
	 * Sign bit is the XOR of the operand sign bits.
	 */
	s = (int)((x ^ y) >> 63);

	/*
	 * Corrective actions for zeros: if either of the operands is
	 * zero, then the computations above were wrong. Test for zero
	 * is whether ex or ey is zero. We just have to set the mantissa
	 * (zu) to zero, the FPR() function will normalize e.
	 */
	d = ((ex + 0x7FF) & (ey + 0x7FF)) >> 11;
	zu &= -(uint64_t)d;

	/*
	 * FPR() packs the result and applies proper rounding.
	 */
	return FPR(s, e, zu);
}

#endif // yyyASM_CORTEXM4-

#if FALCON_ASM_CORTEXM4 // yyyASM_CORTEXM4+1

__attribute__((naked))
fpr
fpr_div(fpr x __attribute__((unused)), fpr y __attribute__((unused)))
{
	__asm__ (
	"push	{ r4, r5, r6, r7, r8, r10, r11, lr }\n\t"

	"@ Extract mantissas of x and y, in r0:r4 and r2:r5, respectively.\n\t"
	"@ We don't touch r1 and r3 as they contain the exponents and\n\t"
	"@ signs, which we'll need later on.\n\t"
	"ubfx	r4, r1, #0, #20\n\t"
	"ubfx	r5, r3, #0, #20\n\t"
	"orrs	r4, r4, #0x00100000\n\t"
	"orrs	r5, r5, #0x00100000\n\t"
	"\n\t"
	"@ Perform bit-by-bit division. We want a 56-bit result in r8:r10\n\t"
	"@ (low bit is 0). Bits come from the carry flag and are\n\t"
	"@ injected with rrx, i.e. in position 31; we thus get bits in\n\t"
	"@ the reverse order. Bits accumulate in r8; after the first 24\n\t"
	"@ bits, we move the quotient bits to r10.\n\t"
	"eors	r8, r8\n\t"
	"\n\t"

#define DIVSTEP \
	"subs	r6, r0, r2\n\t" \
	"sbcs	r7, r4, r5\n\t" \
	"rrx	r8, r8\n\t" \
	"ands	r6, r2, r8, asr #31\n\t" \
	"ands	r7, r5, r8, asr #31\n\t" \
	"subs	r0, r6\n\t" \
	"sbcs	r4, r7\n\t" \
	"adds	r0, r0, r0\n\t" \
	"adcs	r4, r4, r4\n\t"

#define DIVSTEP4   DIVSTEP DIVSTEP DIVSTEP DIVSTEP
#define DIVSTEP8   DIVSTEP4 DIVSTEP4

	DIVSTEP8
	DIVSTEP8
	DIVSTEP8

	"\n\t"
	"@ We have the first 24 bits of the quotient, move them to r10.\n\t"
	"rbit	r10, r8\n\t"
	"\n\t"

	DIVSTEP8
	DIVSTEP8
	DIVSTEP8
	DIVSTEP4 DIVSTEP DIVSTEP DIVSTEP

#undef DIVSTEP
#undef DIVSTEP4
#undef DIVSTEP8

	"\n\t"
	"@ Lowest bit will be set if remainder is non-zero at this point\n\t"
	"@ (this is the 'sticky' bit).\n\t"
	"subs	r0, #1\n\t"
	"sbcs	r4, #0\n\t"
	"rrx	r8, r8\n\t"
	"\n\t"
	"@ We now have the next (low) 32 bits of the quotient.\n\t"
	"rbit	r8, r8\n\t"
	"\n\t"
	"@ Since both operands had their top bit set, we know that the\n\t"
	"@ result at this point is in 2^54..2^56-1. We scale it down\n\t"
	"@ to 2^54..2^55-1 with a conditional shift. We also write the\n\t"
	"@ result in r4:r5. If the shift is done, r6 will contain -1.\n\t"
	"ands	r4, r8, #1\n\t"
	"lsrs	r6, r10, #23\n\t"
	"rsbs	r6, r6, #0\n\t"
	"orrs	r4, r4, r8, lsr #1\n\t"
	"orrs	r4, r4, r10, lsl #31\n\t"
	"lsrs	r5, r10, #1\n\t"
	"eors	r8, r8, r4\n\t"
	"eors	r10, r10, r5\n\t"
	"bics	r8, r8, r6\n\t"
	"bics	r10, r10, r6\n\t"
	"eors	r4, r4, r8\n\t"
	"eors	r5, r5, r10\n\t"
	"\n\t"
	"@ Compute aggregate exponent: ex - ey + 1022 + w\n\t"
	"@ (where w = 1 if the conditional shift was done, 0 otherwise)\n\t"
	"@ But we subtract 1 because the injection of the mantissa high\n\t"
	"@ bit will increment the exponent by 1.\n\t"
	"lsls	r0, r1, #1\n\t"
	"lsls	r2, r3, #1\n\t"
	"lsrs	r0, r0, #21\n\t"
	"addw	r7, r0, #0x7FF  @ save ex + 2047 in r7\n\t"
	"subs	r0, r0, r2, lsr #21\n\t"
	"addw	r0, r0, #1021\n\t"
	"subs	r0, r6\n\t"
	"\n\t"
	"@ If the x operand was zero, then the computation was wrong and\n\t"
	"@ the result is zero. Also, if the result exponent is zero or\n\t"
	"@ negative, then the mantissa shall be clamped to zero. Since r0\n\t"
	"@ contains the result exponent minus 1, we test on r0 being\n\t"
	"@ strictly negative.\n\t"
	"mvns	r2, r0\n\t"
	"ands	r2, r2, r7, lsl #20\n\t"
	"ands	r0, r0, r2, asr #31\n\t"
	"ands	r4, r4, r2, asr #31\n\t"
	"ands	r5, r5, r2, asr #31\n\t"
	"\n\t"
	"@ Sign is the XOR of the sign of the operands. This is true in\n\t"
	"@ all cases, including very small results (exponent underflow)\n\t"
	"@ and zeros.\n\t"
	"eors	r1, r3\n\t"
	"bfc	r1, #0, #31\n\t"
	"\n\t"
	"@ Plug in the exponent.\n\t"
	"bfi	r1, r0, #20, #11\n\t"
	"\n\t"
	"@ Shift back to the normal 53-bit mantissa, with rounding.\n\t"
	"@ Mantissa goes into r0:r1. For r1, we must use an addition\n\t"
	"@ because the rounding may have triggered a carry, that should\n\t"
	"@ be added to the exponent.\n\t"
	"movs	r6, r4\n\t"
	"lsrs	r0, r4, #2\n\t"
	"orrs	r0, r0, r5, lsl #30\n\t"
	"adds	r1, r1, r5, lsr #2\n\t"
	"ands	r6, #0x7\n\t"
	"movs	r3, #0xC8\n\t"
	"lsrs	r3, r6\n\t"
	"ands	r3, #1\n\t"
	"adds	r0, r3\n\t"
	"adcs	r1, #0\n\t"
	"\n\t"
	"pop	{ r4, r5, r6, r7, r8, r10, r11, pc }\n\t"
	);
}

#else // yyyASM_CORTEXM4+0

fpr
fpr_div(fpr x, fpr y)
{
	uint64_t xu, yu, q, q2, w;
	int i, ex, ey, e, d, s;

	/*
	 * Extract mantissas of x and y (unsigned).
	 */
	xu = (x & (((uint64_t)1 << 52) - 1)) | ((uint64_t)1 << 52);
	yu = (y & (((uint64_t)1 << 52) - 1)) | ((uint64_t)1 << 52);

	/*
	 * Perform bit-by-bit division of xu by yu. We run it for 55 bits.
	 */
	q = 0;
	for (i = 0; i < 55; i ++) {
		/*
		 * If yu is less than or equal xu, then subtract it and
		 * push a 1 in the quotient; otherwise, leave xu unchanged
		 * and push a 0.
		 */
		uint64_t b;

		b = ((xu - yu) >> 63) - 1;
		xu -= b & yu;
		q |= b & 1;
		xu <<= 1;
		q <<= 1;
	}

	/*
	 * We got 55 bits in the quotient, followed by an extra zero. We
	 * want that 56th bit to be "sticky": it should be a 1 if and
	 * only if the remainder (xu) is non-zero.
	 */
	q |= (xu | -xu) >> 63;

	/*
	 * Quotient is at most 2^56-1. Its top bit may be zero, but in
	 * that case the next-to-top bit will be a one, since the
	 * initial xu and yu were both in the 2^52..2^53-1 range.
	 * We perform a conditional shift to normalize q to the
	 * 2^54..2^55-1 range (with the bottom bit being sticky).
	 */
	q2 = (q >> 1) | (q & 1);
	w = q >> 55;
	q ^= (q ^ q2) & -w;

	/*
	 * Extract exponents to compute the scaling factor:
	 *
	 *   - Each exponent is biased and we scaled them up by
	 *     52 bits; but these biases will cancel out.
	 *
	 *   - The division loop produced a 55-bit shifted result,
	 *     so we must scale it down by 55 bits.
	 *
	 *   - If w = 1, we right-shifted the integer by 1 bit,
	 *     hence we must add 1 to the scaling.
	 */
	ex = (int)((x >> 52) & 0x7FF);
	ey = (int)((y >> 52) & 0x7FF);
	e = ex - ey - 55 + (int)w;

	/*
	 * Sign is the XOR of the signs of the operands.
	 */
	s = (int)((x ^ y) >> 63);

	/*
	 * Corrective actions for zeros: if x = 0, then the computation
	 * is wrong, and we must clamp e and q to 0. We do not care
	 * about the case y = 0 (as per assumptions in this module,
	 * the caller does not perform divisions by zero).
	 */
	d = (ex + 0x7FF) >> 11;
	s &= d;
	e &= -d;
	q &= -(uint64_t)d;

	/*
	 * FPR() packs the result and applies proper rounding.
	 */
	return FPR(s, e, q);
}

#endif // yyyASM_CORTEXM4-

#if FALCON_ASM_CORTEXM4 // yyyASM_CORTEXM4+1

__attribute__((naked))
fpr
fpr_sqrt(fpr x __attribute__((unused)))
{
	__asm__ (
	"push	{ r4, r5, r6, r7, r8, r10, r11, lr }\n\t"
	"\n\t"
	"@ Extract mantissa (r0:r1) and exponent (r2). We assume that the\n\t"
	"@ sign is positive. If the source is zero, then the mantissa is\n\t"
	"@ set to 0.\n\t"
	"lsrs	r2, r1, #20\n\t"
	"bfc	r1, #20, #12\n\t"
	"addw	r3, r2, #0x7FF\n\t"
	"subw	r2, r2, #1023\n\t"
	"lsrs	r3, r3, #11\n\t"
	"orrs	r1, r1, r3, lsl #20\n\t"
	"\n\t"
	"@ If the exponent is odd, then multiply mantissa by 2 and subtract\n\t"
	"@ 1 from the exponent.\n\t"
	"ands	r3, r2, #1\n\t"
	"subs	r2, r2, r3\n\t"
	"rsbs	r3, r3, #0\n\t"
	"ands	r4, r1, r3\n\t"
	"ands	r3, r0\n\t"
	"adds	r0, r3\n\t"
	"adcs	r1, r4\n\t"
	"\n\t"
	"@ Left-shift the mantissa by 9 bits to put it in the\n\t"
	"@ 2^61..2^63-1 range (unless it is exactly 0).\n\t"
	"lsls	r1, r1, #9\n\t"
	"orrs	r1, r1, r0, lsr #23\n\t"
	"lsls	r0, r0, #9\n\t"
	"\n\t"
	"@ Compute the square root bit-by-bit.\n\t"
	"@ There are 54 iterations; first 30 can work on top word only.\n\t"
	"@   q = r3 (bit-reversed)\n\t"
	"@   s = r5\n\t"
	"eors	r3, r3\n\t"
	"eors	r5, r5\n\t"

#define SQRT_STEP_HI(bit) \
	"orrs	r6, r5, #(1 << (" #bit "))\n\t" \
	"subs	r7, r1, r6\n\t" \
	"rrx	r3, r3\n\t" \
	"ands	r6, r6, r3, asr #31\n\t" \
	"subs	r1, r1, r6\n\t" \
	"lsrs	r6, r3, #31\n\t" \
	"orrs	r5, r5, r6, lsl #((" #bit ") + 1)\n\t" \
	"adds	r0, r0\n\t" \
	"adcs	r1, r1\n\t"

#define SQRT_STEP_HIx5(b)  \
		SQRT_STEP_HI((b)+4) \
		SQRT_STEP_HI((b)+3) \
		SQRT_STEP_HI((b)+2) \
		SQRT_STEP_HI((b)+1) \
		SQRT_STEP_HI(b)

	SQRT_STEP_HIx5(25)
	SQRT_STEP_HIx5(20)
	SQRT_STEP_HIx5(15)
	SQRT_STEP_HIx5(10)
	SQRT_STEP_HIx5(5)
	SQRT_STEP_HIx5(0)

#undef SQRT_STEP_HI
#undef SQRT_STEP_HIx5

	"@ Top 30 bits of the result must be reversed: they were\n\t"
	"@ accumulated with rrx (hence from the top bit).\n\t"
	"rbit	r3, r3\n\t"
	"\n\t"
	"@ For the next 24 iterations, we must use two-word operations.\n\t"
	"@   bits of q now accumulate in r4\n\t"
	"@   s is in r6:r5\n\t"
	"eors	r4, r4\n\t"
	"eors	r6, r6\n\t"
	"\n\t"
	"@ First iteration is special because the potential bit goes into\n\t"
	"@ r5, not r6.\n\t"
	"orrs	r7, r6, #(1 << 31)\n\t"
	"subs	r8, r0, r7\n\t"
	"sbcs	r10, r1, r5\n\t"
	"rrx	r4, r4\n\t"
	"ands	r7, r7, r4, asr #31\n\t"
	"ands	r8, r5, r4, asr #31\n\t"
	"subs	r0, r0, r7\n\t"
	"sbcs	r1, r1, r8\n\t"
	"lsrs	r7, r4, #31\n\t"
	"orrs	r5, r5, r4, lsr #31\n\t"
	"adds	r0, r0\n\t"
	"adcs	r1, r1\n\t"

#define SQRT_STEP_LO(bit) \
	"orrs	r7, r6, #(1 << (" #bit "))\n\t" \
	"subs	r8, r0, r7\n\t" \
	"sbcs	r10, r1, r5\n\t" \
	"rrx	r4, r4\n\t" \
	"ands	r7, r7, r4, asr #31\n\t" \
	"ands	r8, r5, r4, asr #31\n\t" \
	"subs	r0, r0, r7\n\t" \
	"sbcs	r1, r1, r8\n\t" \
	"lsrs	r7, r4, #31\n\t" \
	"orrs	r6, r6, r7, lsl #((" #bit ") + 1)\n\t" \
	"adds	r0, r0\n\t" \
	"adcs	r1, r1\n\t"

#define SQRT_STEP_LOx4(b) \
		SQRT_STEP_LO((b)+3) \
		SQRT_STEP_LO((b)+2) \
		SQRT_STEP_LO((b)+1) \
		SQRT_STEP_LO(b)

	SQRT_STEP_LO(30)
	SQRT_STEP_LO(29)
	SQRT_STEP_LO(28)
	SQRT_STEP_LOx4(24)
	SQRT_STEP_LOx4(20)
	SQRT_STEP_LOx4(16)
	SQRT_STEP_LOx4(12)
	SQRT_STEP_LOx4(8)

#undef SQRT_STEP_LO
#undef SQRT_STEP_LOx4

	"@ Put low 24 bits in the right order.\n\t"
	"rbit	r4, r4\n\t"
	"\n\t"
	"@ We have a 54-bit result; compute the 55-th bit as the 'sticky'\n\t"
	"@ bit: it is non-zero if and only if r0:r1 is non-zero. We put the\n\t"
	"@ three low bits (including the sticky bit) in r5.\n\t"
	"orrs	r0, r1\n\t"
	"rsbs	r1, r0, #0\n\t"
	"orrs	r0, r1\n\t"
	"lsls	r5, r4, #1\n\t"
	"orrs	r5, r5, r0, lsr #31\n\t"
	"ands	r5, #0x7\n\t"
	"\n\t"
	"@ Compute the rounding: r6 is set to 0 or 1, and will be added\n\t"
	"@ to the mantissa.\n\t"
	"movs	r6, #0xC8\n\t"
	"lsrs	r6, r5\n\t"
	"ands	r6, #1\n\t"
	"\n\t"
	"@ Put the mantissa (53 bits, in the 2^52..2^53-1 range) in r0:r1\n\t"
	"@ (rounding not applied yet).\n\t"
	"lsrs	r0, r4, #1\n\t"
	"orrs	r0, r0, r3, lsl #23\n\t"
	"lsrs	r1, r3, #9\n\t"
	"\n\t"
	"@ Compute new exponent. This is half the old one (then reencoded\n\t"
	"@ by adding 1023). Exception: if the mantissa is zero, then the\n\t"
	"@ encoded exponent is set to 0. At that point, if the mantissa\n\t"
	"@ is non-zero, then its high bit (bit 52, i.e. bit 20 of r1) is\n\t"
	"@ non-zero. Note that the exponent cannot go out of range.\n\t"
	"lsrs	r2, r2, #1\n\t"
	"addw	r2, r2, #1023\n\t"
	"lsrs	r5, r1, #20\n\t"
	"rsbs	r5, r5, #0\n\t"
	"ands	r2, r5\n\t"
	"\n\t"
	"@ Place exponent. This overwrites the high bit of the mantissa.\n\t"
	"bfi	r1, r2, #20, #11\n\t"
	"\n\t"
	"@ Apply rounding. This may create a carry that will spill into\n\t"
	"@ the exponent, which is exactly what should be done in that case\n\t"
	"@ (i.e. increment the exponent).\n\t"
	"adds	r0, r0, r6\n\t"
	"adcs	r1, r1, #0\n\t"
	"\n\t"
	"pop	{ r4, r5, r6, r7, r8, r10, r11, pc }\n\t"
	);
}

#else // yyyASM_CORTEXM4+0

fpr
fpr_sqrt(fpr x)
{
	uint64_t xu, q, s, r;
	int ex, e;

	/*
	 * Extract the mantissa and the exponent. We don't care about
	 * the sign: by assumption, the operand is nonnegative.
	 * We want the "true" exponent corresponding to a mantissa
	 * in the 1..2 range.
	 */
	xu = (x & (((uint64_t)1 << 52) - 1)) | ((uint64_t)1 << 52);
	ex = (int)((x >> 52) & 0x7FF);
	e = ex - 1023;

	/*
	 * If the exponent is odd, double the mantissa and decrement
	 * the exponent. The exponent is then halved to account for
	 * the square root.
	 */
	xu += xu & -(uint64_t)(e & 1);
	e >>= 1;

	/*
	 * Double the mantissa.
	 */
	xu <<= 1;

	/*
	 * We now have a mantissa in the 2^53..2^55-1 range. It
	 * represents a value between 1 (inclusive) and 4 (exclusive)
	 * in fixed point notation (with 53 fractional bits). We
	 * compute the square root bit by bit.
	 */
	q = 0;
	s = 0;
	r = (uint64_t)1 << 53;
	for (int i = 0; i < 54; i ++) {
		uint64_t t, b;

		t = s + r;
		b = ((xu - t) >> 63) - 1;
		s += (r << 1) & b;
		xu -= t & b;
		q += r & b;
		xu <<= 1;
		r >>= 1;
	}

	/*
	 * Now, q is a rounded-low 54-bit value, with a leading 1,
	 * 52 fractional digits, and an additional guard bit. We add
	 * an extra sticky bit to account for what remains of the operand.
	 */
	q <<= 1;
	q |= (xu | -xu) >> 63;

	/*
	 * Result q is in the 2^54..2^55-1 range; we bias the exponent
	 * by 54 bits (the value e at that point contains the "true"
	 * exponent, but q is now considered an integer, i.e. scaled
	 * up.
	 */
	e -= 54;

	/*
	 * Corrective action for an operand of value zero.
	 */
	q &= -(uint64_t)((ex + 0x7FF) >> 11);

	/*
	 * Apply rounding and back result.
	 */
	return FPR(0, e, q);
}

#endif // yyyASM_CORTEXM4-

uint64_t
fpr_expm_p63(fpr x, fpr ccs)
{
	/*
	 * Polynomial approximation of exp(-x) is taken from FACCT:
	 *   https://eprint.iacr.org/2018/1234
	 * Specifically, values are extracted from the implementation
	 * referenced from the FACCT article, and available at:
	 *   https://github.com/raykzhao/gaussian
	 * Here, the coefficients have been scaled up by 2^63 and
	 * converted to integers.
	 *
	 * Tests over more than 24 billions of random inputs in the
	 * 0..log(2) range have never shown a deviation larger than
	 * 2^(-50) from the true mathematical value.
	 */
	static const uint64_t C[] = {
		0x00000004741183A3u,
		0x00000036548CFC06u,
		0x0000024FDCBF140Au,
		0x0000171D939DE045u,
		0x0000D00CF58F6F84u,
		0x000680681CF796E3u,
		0x002D82D8305B0FEAu,
		0x011111110E066FD0u,
		0x0555555555070F00u,
		0x155555555581FF00u,
		0x400000000002B400u,
		0x7FFFFFFFFFFF4800u,
		0x8000000000000000u
	};

	uint64_t z, y;
	unsigned u;
	uint32_t z0, z1, y0, y1;
	uint64_t a, b;

	y = C[0];
	z = (uint64_t)fpr_trunc(fpr_mul(x, fpr_ptwo63)) << 1;
	print_debug(("\t\t\t\tz = %016lx\n", z));
	for (u = 1; u < (sizeof C) / sizeof(C[0]); u ++) {
		/*
		 * Compute product z * y over 128 bits, but keep only
		 * the top 64 bits.
		 *
		 * TODO: On some architectures/compilers we could use
		 * some intrinsics (__umulh() on MSVC) or other compiler
		 * extensions (unsigned __int128 on GCC / Clang) for
		 * improved speed; however, most 64-bit architectures
		 * also have appropriate IEEE754 floating-point support,
		 * which is better.
		 */
		uint64_t c;

		z0 = (uint32_t)z;
		z1 = (uint32_t)(z >> 32);
		y0 = (uint32_t)y;
		y1 = (uint32_t)(y >> 32);
		a = ((uint64_t)z0 * (uint64_t)y1)
			+ (((uint64_t)z0 * (uint64_t)y0) >> 32);
		b = ((uint64_t)z1 * (uint64_t)y0);
		c = (a >> 32) + (b >> 32);
		c += (((uint64_t)(uint32_t)a + (uint64_t)(uint32_t)b) >> 32);
		c += (uint64_t)z1 * (uint64_t)y1;
		y = C[u] - c;
	}
	print_debug(("\t\t\t\ty = %016lx\n", y));

	/*
	 * The scaling factor must be applied at the end. Since y is now
	 * in fixed-point notation, we have to convert the factor to the
	 * same format, and do an extra integer multiplication.
	 */
	z = (uint64_t)fpr_trunc(fpr_mul(ccs, fpr_ptwo63)) << 1;
	print_debug(("\t\t\t\tz = %016lx\n", z));
	z0 = (uint32_t)z;
	z1 = (uint32_t)(z >> 32);
	y0 = (uint32_t)y;
	y1 = (uint32_t)(y >> 32);
	a = ((uint64_t)z0 * (uint64_t)y1)
		+ (((uint64_t)z0 * (uint64_t)y0) >> 32);
	b = ((uint64_t)z1 * (uint64_t)y0);
	y = (a >> 32) + (b >> 32);
	y += (((uint64_t)(uint32_t)a + (uint64_t)(uint32_t)b) >> 32);
	y += (uint64_t)z1 * (uint64_t)y1;

	return y;
}

const fpr fpr_gm_tab[] = {
	0, 0,
	 9223372036854775808U,  4607182418800017408U,
	 4604544271217802189U,  4604544271217802189U,
	13827916308072577997U,  4604544271217802189U,
	 4606496786581982534U,  4600565431771507043U,
	13823937468626282851U,  4606496786581982534U,
	 4600565431771507043U,  4606496786581982534U,
	13829868823436758342U,  4600565431771507043U,
	 4607009347991985328U,  4596196889902818827U,
	13819568926757594635U,  4607009347991985328U,
	 4603179351334086856U,  4605664432017547683U,
	13829036468872323491U,  4603179351334086856U,
	 4605664432017547683U,  4603179351334086856U,
	13826551388188862664U,  4605664432017547683U,
	 4596196889902818827U,  4607009347991985328U,
	13830381384846761136U,  4596196889902818827U,
	 4607139046673687846U,  4591727299969791020U,
	13815099336824566828U,  4607139046673687846U,
	 4603889326261607894U,  4605137878724712257U,
	13828509915579488065U,  4603889326261607894U,
	 4606118860100255153U,  4602163548591158843U,
	13825535585445934651U,  4606118860100255153U,
	 4598900923775164166U,  4606794571824115162U,
	13830166608678890970U,  4598900923775164166U,
	 4606794571824115162U,  4598900923775164166U,
	13822272960629939974U,  4606794571824115162U,
	 4602163548591158843U,  4606118860100255153U,
	13829490896955030961U,  4602163548591158843U,
	 4605137878724712257U,  4603889326261607894U,
	13827261363116383702U,  4605137878724712257U,
	 4591727299969791020U,  4607139046673687846U,
	13830511083528463654U,  4591727299969791020U,
	 4607171569234046334U,  4587232218149935124U,
	13810604255004710932U,  4607171569234046334U,
	 4604224084862889120U,  4604849113969373103U,
	13828221150824148911U,  4604224084862889120U,
	 4606317631232591731U,  4601373767755717824U,
	13824745804610493632U,  4606317631232591731U,
	 4599740487990714333U,  4606655894547498725U,
	13830027931402274533U,  4599740487990714333U,
	 4606912484326125783U,  4597922303871901467U,
	13821294340726677275U,  4606912484326125783U,
	 4602805845399633902U,  4605900952042040894U,
	13829272988896816702U,  4602805845399633902U,
	 4605409869824231233U,  4603540801876750389U,
	13826912838731526197U,  4605409869824231233U,
	 4594454542771183930U,  4607084929468638487U,
	13830456966323414295U,  4594454542771183930U,
	 4607084929468638487U,  4594454542771183930U,
	13817826579625959738U,  4607084929468638487U,
	 4603540801876750389U,  4605409869824231233U,
	13828781906679007041U,  4603540801876750389U,
	 4605900952042040894U,  4602805845399633902U,
	13826177882254409710U,  4605900952042040894U,
	 4597922303871901467U,  4606912484326125783U,
	13830284521180901591U,  4597922303871901467U,
	 4606655894547498725U,  4599740487990714333U,
	13823112524845490141U,  4606655894547498725U,
	 4601373767755717824U,  4606317631232591731U,
	13829689668087367539U,  4601373767755717824U,
	 4604849113969373103U,  4604224084862889120U,
	13827596121717664928U,  4604849113969373103U,
	 4587232218149935124U,  4607171569234046334U,
	13830543606088822142U,  4587232218149935124U,
	 4607179706000002317U,  4582730748936808062U,
	13806102785791583870U,  4607179706000002317U,
	 4604386048625945823U,  4604698657331085206U,
	13828070694185861014U,  4604386048625945823U,
	 4606409688975526202U,  4600971798440897930U,
	13824343835295673738U,  4606409688975526202U,
	 4600154912527631775U,  4606578871587619388U,
	13829950908442395196U,  4600154912527631775U,
	 4606963563043808649U,  4597061974398750563U,
	13820434011253526371U,  4606963563043808649U,
	 4602994049708411683U,  4605784983948558848U,
	13829157020803334656U,  4602994049708411683U,
	 4605539368864982914U,  4603361638657888991U,
	13826733675512664799U,  4605539368864982914U,
	 4595327571478659014U,  4607049811591515049U,
	13830421848446290857U,  4595327571478659014U,
	 4607114680469659603U,  4593485039402578702U,
	13816857076257354510U,  4607114680469659603U,
	 4603716733069447353U,  4605276012900672507U,
	13828648049755448315U,  4603716733069447353U,
	 4606012266443150634U,  4602550884377336506U,
	13825922921232112314U,  4606012266443150634U,
	 4598476289818621559U,  4606856142606846307U,
	13830228179461622115U,  4598476289818621559U,
	 4606727809065869586U,  4599322407794599425U,
	13822694444649375233U,  4606727809065869586U,
	 4601771097584682078U,  4606220668805321205U,
	13829592705660097013U,  4601771097584682078U,
	 4604995550503212910U,  4604058477489546729U,
	13827430514344322537U,  4604995550503212910U,
	 4589965306122607094U,  4607158013403433018U,
	13830530050258208826U,  4589965306122607094U,
	 4607158013403433018U,  4589965306122607094U,
	13813337342977382902U,  4607158013403433018U,
	 4604058477489546729U,  4604995550503212910U,
	13828367587357988718U,  4604058477489546729U,
	 4606220668805321205U,  4601771097584682078U,
	13825143134439457886U,  4606220668805321205U,
	 4599322407794599425U,  4606727809065869586U,
	13830099845920645394U,  4599322407794599425U,
	 4606856142606846307U,  4598476289818621559U,
	13821848326673397367U,  4606856142606846307U,
	 4602550884377336506U,  4606012266443150634U,
	13829384303297926442U,  4602550884377336506U,
	 4605276012900672507U,  4603716733069447353U,
	13827088769924223161U,  4605276012900672507U,
	 4593485039402578702U,  4607114680469659603U,
	13830486717324435411U,  4593485039402578702U,
	 4607049811591515049U,  4595327571478659014U,
	13818699608333434822U,  4607049811591515049U,
	 4603361638657888991U,  4605539368864982914U,
	13828911405719758722U,  4603361638657888991U,
	 4605784983948558848U,  4602994049708411683U,
	13826366086563187491U,  4605784983948558848U,
	 4597061974398750563U,  4606963563043808649U,
	13830335599898584457U,  4597061974398750563U,
	 4606578871587619388U,  4600154912527631775U,
	13823526949382407583U,  4606578871587619388U,
	 4600971798440897930U,  4606409688975526202U,
	13829781725830302010U,  4600971798440897930U,
	 4604698657331085206U,  4604386048625945823U,
	13827758085480721631U,  4604698657331085206U,
	 4582730748936808062U,  4607179706000002317U,
	13830551742854778125U,  4582730748936808062U,
	 4607181740574479067U,  4578227681973159812U,
	13801599718827935620U,  4607181740574479067U,
	 4604465633578481725U,  4604621949701367983U,
	13827993986556143791U,  4604465633578481725U,
	 4606453861145241227U,  4600769149537129431U,
	13824141186391905239U,  4606453861145241227U,
	 4600360675823176935U,  4606538458821337243U,
	13829910495676113051U,  4600360675823176935U,
	 4606987119037722413U,  4596629994023683153U,
	13820002030878458961U,  4606987119037722413U,
	 4603087070374583113U,  4605725276488455441U,
	13829097313343231249U,  4603087070374583113U,
	 4605602459698789090U,  4603270878689749849U,
	13826642915544525657U,  4605602459698789090U,
	 4595762727260045105U,  4607030246558998647U,
	13830402283413774455U,  4595762727260045105U,
	 4607127537664763515U,  4592606767730311893U,
	13815978804585087701U,  4607127537664763515U,
	 4603803453461190356U,  4605207475328619533U,
	13828579512183395341U,  4603803453461190356U,
	 4606066157444814153U,  4602357870542944470U,
	13825729907397720278U,  4606066157444814153U,
	 4598688984595225406U,  4606826008603986804U,
	13830198045458762612U,  4598688984595225406U,
	 4606761837001494797U,  4599112075441176914U,
	13822484112295952722U,  4606761837001494797U,
	 4601967947786150793U,  4606170366472647579U,
	13829542403327423387U,  4601967947786150793U,
	 4605067233569943231U,  4603974338538572089U,
	13827346375393347897U,  4605067233569943231U,
	 4590846768565625881U,  4607149205763218185U,
	13830521242617993993U,  4590846768565625881U,
	 4607165468267934125U,  4588998070480937184U,
	13812370107335712992U,  4607165468267934125U,
	 4604141730443515286U,  4604922840319727473U,
	13828294877174503281U,  4604141730443515286U,
	 4606269759522929756U,  4601573027631668967U,
	13824945064486444775U,  4606269759522929756U,
	 4599531889160152938U,  4606692493141721470U,
	13830064529996497278U,  4599531889160152938U,
	 4606884969294623682U,  4598262871476403630U,
	13821634908331179438U,  4606884969294623682U,
	 4602710690099904183U,  4605957195211051218U,
	13829329232065827026U,  4602710690099904183U,
	 4605343481119364930U,  4603629178146150899U,
	13827001215000926707U,  4605343481119364930U,
	 4594016801320007031U,  4607100477024622401U,
	13830472513879398209U,  4594016801320007031U,
	 4607068040143112603U,  4594891488091520602U,
	13818263524946296410U,  4607068040143112603U,
	 4603451617570386922U,  4605475169017376660U,
	13828847205872152468U,  4603451617570386922U,
	 4605843545406134034U,  4602900303344142735U,
	13826272340198918543U,  4605843545406134034U,
	 4597492765973365521U,  4606938683557690074U,
	13830310720412465882U,  4597492765973365521U,
	 4606618018794815019U,  4599948172872067014U,
	13823320209726842822U,  4606618018794815019U,
	 4601173347964633034U,  4606364276725003740U,
	13829736313579779548U,  4601173347964633034U,
	 4604774382555066977U,  4604305528345395596U,
	13827677565200171404U,  4604774382555066977U,
	 4585465300892538317U,  4607176315382986589U,
	13830548352237762397U,  4585465300892538317U,
	 4607176315382986589U,  4585465300892538317U,
	13808837337747314125U,  4607176315382986589U,
	 4604305528345395596U,  4604774382555066977U,
	13828146419409842785U,  4604305528345395596U,
	 4606364276725003740U,  4601173347964633034U,
	13824545384819408842U,  4606364276725003740U,
	 4599948172872067014U,  4606618018794815019U,
	13829990055649590827U,  4599948172872067014U,
	 4606938683557690074U,  4597492765973365521U,
	13820864802828141329U,  4606938683557690074U,
	 4602900303344142735U,  4605843545406134034U,
	13829215582260909842U,  4602900303344142735U,
	 4605475169017376660U,  4603451617570386922U,
	13826823654425162730U,  4605475169017376660U,
	 4594891488091520602U,  4607068040143112603U,
	13830440076997888411U,  4594891488091520602U,
	 4607100477024622401U,  4594016801320007031U,
	13817388838174782839U,  4607100477024622401U,
	 4603629178146150899U,  4605343481119364930U,
	13828715517974140738U,  4603629178146150899U,
	 4605957195211051218U,  4602710690099904183U,
	13826082726954679991U,  4605957195211051218U,
	 4598262871476403630U,  4606884969294623682U,
	13830257006149399490U,  4598262871476403630U,
	 4606692493141721470U,  4599531889160152938U,
	13822903926014928746U,  4606692493141721470U,
	 4601573027631668967U,  4606269759522929756U,
	13829641796377705564U,  4601573027631668967U,
	 4604922840319727473U,  4604141730443515286U,
	13827513767298291094U,  4604922840319727473U,
	 4588998070480937184U,  4607165468267934125U,
	13830537505122709933U,  4588998070480937184U,
	 4607149205763218185U,  4590846768565625881U,
	13814218805420401689U,  4607149205763218185U,
	 4603974338538572089U,  4605067233569943231U,
	13828439270424719039U,  4603974338538572089U,
	 4606170366472647579U,  4601967947786150793U,
	13825339984640926601U,  4606170366472647579U,
	 4599112075441176914U,  4606761837001494797U,
	13830133873856270605U,  4599112075441176914U,
	 4606826008603986804U,  4598688984595225406U,
	13822061021450001214U,  4606826008603986804U,
	 4602357870542944470U,  4606066157444814153U,
	13829438194299589961U,  4602357870542944470U,
	 4605207475328619533U,  4603803453461190356U,
	13827175490315966164U,  4605207475328619533U,
	 4592606767730311893U,  4607127537664763515U,
	13830499574519539323U,  4592606767730311893U,
	 4607030246558998647U,  4595762727260045105U,
	13819134764114820913U,  4607030246558998647U,
	 4603270878689749849U,  4605602459698789090U,
	13828974496553564898U,  4603270878689749849U,
	 4605725276488455441U,  4603087070374583113U,
	13826459107229358921U,  4605725276488455441U,
	 4596629994023683153U,  4606987119037722413U,
	13830359155892498221U,  4596629994023683153U,
	 4606538458821337243U,  4600360675823176935U,
	13823732712677952743U,  4606538458821337243U,
	 4600769149537129431U,  4606453861145241227U,
	13829825898000017035U,  4600769149537129431U,
	 4604621949701367983U,  4604465633578481725U,
	13827837670433257533U,  4604621949701367983U,
	 4578227681973159812U,  4607181740574479067U,
	13830553777429254875U,  4578227681973159812U,
	 4607182249242036882U,  4573724215515480177U,
	13797096252370255985U,  4607182249242036882U,
	 4604505071555817232U,  4604583231088591477U,
	13827955267943367285U,  4604505071555817232U,
	 4606475480113671417U,  4600667422348321968U,
	13824039459203097776U,  4606475480113671417U,
	 4600463181646572228U,  4606517779747998088U,
	13829889816602773896U,  4600463181646572228U,
	 4606998399608725124U,  4596413578358834022U,
	13819785615213609830U,  4606998399608725124U,
	 4603133304188877240U,  4605694995810664660U,
	13829067032665440468U,  4603133304188877240U,
	 4605633586259814045U,  4603225210076562971U,
	13826597246931338779U,  4605633586259814045U,
	 4595979936813835462U,  4607019963775302583U,
	13830392000630078391U,  4595979936813835462U,
	 4607133460805585796U,  4592167175087283203U,
	13815539211942059011U,  4607133460805585796U,
	 4603846496621587377U,  4605172808754305228U,
	13828544845609081036U,  4603846496621587377U,
	 4606092657816072624U,  4602260871257280788U,
	13825632908112056596U,  4606092657816072624U,
	 4598795050632330097U,  4606810452769876110U,
	13830182489624651918U,  4598795050632330097U,
	 4606778366364612594U,  4599006600037663623U,
	13822378636892439431U,  4606778366364612594U,
	 4602065906208722008U,  4606144763310860551U,
	13829516800165636359U,  4602065906208722008U,
	 4605102686554936490U,  4603931940768740167U,
	13827303977623515975U,  4605102686554936490U,
	 4591287158938884897U,  4607144295058764886U,
	13830516331913540694U,  4591287158938884897U,
	 4607168688050493276U,  4588115294056142819U,
	13811487330910918627U,  4607168688050493276U,
	 4604183020748362039U,  4604886103475043762U,
	13828258140329819570U,  4604183020748362039U,
	 4606293848208650998U,  4601473544562720001U,
	13824845581417495809U,  4606293848208650998U,
	 4599636300858866724U,  4606674353838411301U,
	13830046390693187109U,  4599636300858866724U,
	 4606898891031025132U,  4598136582470364665U,
	13821508619325140473U,  4606898891031025132U,
	 4602758354025980442U,  4605929219593405673U,
	13829301256448181481U,  4602758354025980442U,
	 4605376811039722786U,  4603585091850767959U,
	13826957128705543767U,  4605376811039722786U,
	 4594235767444503503U,  4607092871118901179U,
	13830464907973676987U,  4594235767444503503U,
	 4607076652372832968U,  4594673119063280916U,
	13818045155918056724U,  4607076652372832968U,
	 4603496309891590679U,  4605442656228245717U,
	13828814693083021525U,  4603496309891590679U,
	 4605872393621214213U,  4602853162432841185U,
	13826225199287616993U,  4605872393621214213U,
	 4597707695679609371U,  4606925748668145757U,
	13830297785522921565U,  4597707695679609371U,
	 4606637115963965612U,  4599844446633109139U,
	13823216483487884947U,  4606637115963965612U,
	 4601273700967202825U,  4606341107699334546U,
	13829713144554110354U,  4601273700967202825U,
	 4604811873195349477U,  4604264921241055824U,
	13827636958095831632U,  4604811873195349477U,
	 4586348876009622851U,  4607174111710118367U,
	13830546148564894175U,  4586348876009622851U,
	 4607178180169683960U,  4584498631466405633U,
	13807870668321181441U,  4607178180169683960U,
	 4604345904647073908U,  4604736643460027021U,
	13828108680314802829U,  4604345904647073908U,
	 4606387137437298591U,  4601072712526242277U,
	13824444749381018085U,  4606387137437298591U,
	 4600051662802353687U,  4606598603759044570U,
	13829970640613820378U,  4600051662802353687U,
	 4606951288507767453U,  4597277522845151878U,
	13820649559699927686U,  4606951288507767453U,
	 4602947266358709886U,  4605814408482919348U,
	13829186445337695156U,  4602947266358709886U,
	 4605507406967535927U,  4603406726595779752U,
	13826778763450555560U,  4605507406967535927U,
	 4595109641634432498U,  4607059093103722971U,
	13830431129958498779U,  4595109641634432498U,
	 4607107746899444102U,  4593797652641645341U,
	13817169689496421149U,  4607107746899444102U,
	 4603673059103075106U,  4605309881318010327U,
	13828681918172786135U,  4603673059103075106U,
	 4605984877841711338U,  4602646891659203088U,
	13826018928513978896U,  4605984877841711338U,
	 4598369669086960528U,  4606870719641066940U,
	13830242756495842748U,  4598369669086960528U,
	 4606710311774494716U,  4599427256825614420U,
	13822799293680390228U,  4606710311774494716U,
	 4601672213217083403U,  4606245366082353408U,
	13829617402937129216U,  4601672213217083403U,
	 4604959323120302796U,  4604100215502905499U,
	13827472252357681307U,  4604959323120302796U,
	 4589524267239410099U,  4607161910007591876U,
	13830533946862367684U,  4589524267239410099U,
	 4607153778602162496U,  4590406145430462614U,
	13813778182285238422U,  4607153778602162496U,
	 4604016517974851588U,  4605031521104517324U,
	13828403557959293132U,  4604016517974851588U,
	 4606195668621671667U,  4601869677011524443U,
	13825241713866300251U,  4606195668621671667U,
	 4599217346014614711U,  4606744984357082948U,
	13830117021211858756U,  4599217346014614711U,
	 4606841238740778884U,  4598582729657176439U,
	13821954766511952247U,  4606841238740778884U,
	 4602454542796181607U,  4606039359984203741U,
	13829411396838979549U,  4602454542796181607U,
	 4605241877142478242U,  4603760198400967492U,
	13827132235255743300U,  4605241877142478242U,
	 4593046061348462537U,  4607121277474223905U,
	13830493314328999713U,  4593046061348462537U,
	 4607040195955932526U,  4595545269419264690U,
	13818917306274040498U,  4607040195955932526U,
	 4603316355454250015U,  4605571053506370248U,
	13828943090361146056U,  4603316355454250015U,
	 4605755272910869620U,  4603040651631881451U,
	13826412688486657259U,  4605755272910869620U,
	 4596846128749438754U,  4606975506703684317U,
	13830347543558460125U,  4596846128749438754U,
	 4606558823023444576U,  4600257918160607478U,
	13823629955015383286U,  4606558823023444576U,
	 4600870609507958271U,  4606431930490633905U,
	13829803967345409713U,  4600870609507958271U,
	 4604660425598397818U,  4604425958770613225U,
	13827797995625389033U,  4604660425598397818U,
	 4580962600092897021U,  4607180892816495009U,
	13830552929671270817U,  4580962600092897021U,
	 4607180892816495009U,  4580962600092897021U,
	13804334636947672829U,  4607180892816495009U,
	 4604425958770613225U,  4604660425598397818U,
	13828032462453173626U,  4604425958770613225U,
	 4606431930490633905U,  4600870609507958271U,
	13824242646362734079U,  4606431930490633905U,
	 4600257918160607478U,  4606558823023444576U,
	13829930859878220384U,  4600257918160607478U,
	 4606975506703684317U,  4596846128749438754U,
	13820218165604214562U,  4606975506703684317U,
	 4603040651631881451U,  4605755272910869620U,
	13829127309765645428U,  4603040651631881451U,
	 4605571053506370248U,  4603316355454250015U,
	13826688392309025823U,  4605571053506370248U,
	 4595545269419264690U,  4607040195955932526U,
	13830412232810708334U,  4595545269419264690U,
	 4607121277474223905U,  4593046061348462537U,
	13816418098203238345U,  4607121277474223905U,
	 4603760198400967492U,  4605241877142478242U,
	13828613913997254050U,  4603760198400967492U,
	 4606039359984203741U,  4602454542796181607U,
	13825826579650957415U,  4606039359984203741U,
	 4598582729657176439U,  4606841238740778884U,
	13830213275595554692U,  4598582729657176439U,
	 4606744984357082948U,  4599217346014614711U,
	13822589382869390519U,  4606744984357082948U,
	 4601869677011524443U,  4606195668621671667U,
	13829567705476447475U,  4601869677011524443U,
	 4605031521104517324U,  4604016517974851588U,
	13827388554829627396U,  4605031521104517324U,
	 4590406145430462614U,  4607153778602162496U,
	13830525815456938304U,  4590406145430462614U,
	 4607161910007591876U,  4589524267239410099U,
	13812896304094185907U,  4607161910007591876U,
	 4604100215502905499U,  4604959323120302796U,
	13828331359975078604U,  4604100215502905499U,
	 4606245366082353408U,  4601672213217083403U,
	13825044250071859211U,  4606245366082353408U,
	 4599427256825614420U,  4606710311774494716U,
	13830082348629270524U,  4599427256825614420U,
	 4606870719641066940U,  4598369669086960528U,
	13821741705941736336U,  4606870719641066940U,
	 4602646891659203088U,  4605984877841711338U,
	13829356914696487146U,  4602646891659203088U,
	 4605309881318010327U,  4603673059103075106U,
	13827045095957850914U,  4605309881318010327U,
	 4593797652641645341U,  4607107746899444102U,
	13830479783754219910U,  4593797652641645341U,
	 4607059093103722971U,  4595109641634432498U,
	13818481678489208306U,  4607059093103722971U,
	 4603406726595779752U,  4605507406967535927U,
	13828879443822311735U,  4603406726595779752U,
	 4605814408482919348U,  4602947266358709886U,
	13826319303213485694U,  4605814408482919348U,
	 4597277522845151878U,  4606951288507767453U,
	13830323325362543261U,  4597277522845151878U,
	 4606598603759044570U,  4600051662802353687U,
	13823423699657129495U,  4606598603759044570U,
	 4601072712526242277U,  4606387137437298591U,
	13829759174292074399U,  4601072712526242277U,
	 4604736643460027021U,  4604345904647073908U,
	13827717941501849716U,  4604736643460027021U,
	 4584498631466405633U,  4607178180169683960U,
	13830550217024459768U,  4584498631466405633U,
	 4607174111710118367U,  4586348876009622851U,
	13809720912864398659U,  4607174111710118367U,
	 4604264921241055824U,  4604811873195349477U,
	13828183910050125285U,  4604264921241055824U,
	 4606341107699334546U,  4601273700967202825U,
	13824645737821978633U,  4606341107699334546U,
	 4599844446633109139U,  4606637115963965612U,
	13830009152818741420U,  4599844446633109139U,
	 4606925748668145757U,  4597707695679609371U,
	13821079732534385179U,  4606925748668145757U,
	 4602853162432841185U,  4605872393621214213U,
	13829244430475990021U,  4602853162432841185U,
	 4605442656228245717U,  4603496309891590679U,
	13826868346746366487U,  4605442656228245717U,
	 4594673119063280916U,  4607076652372832968U,
	13830448689227608776U,  4594673119063280916U,
	 4607092871118901179U,  4594235767444503503U,
	13817607804299279311U,  4607092871118901179U,
	 4603585091850767959U,  4605376811039722786U,
	13828748847894498594U,  4603585091850767959U,
	 4605929219593405673U,  4602758354025980442U,
	13826130390880756250U,  4605929219593405673U,
	 4598136582470364665U,  4606898891031025132U,
	13830270927885800940U,  4598136582470364665U,
	 4606674353838411301U,  4599636300858866724U,
	13823008337713642532U,  4606674353838411301U,
	 4601473544562720001U,  4606293848208650998U,
	13829665885063426806U,  4601473544562720001U,
	 4604886103475043762U,  4604183020748362039U,
	13827555057603137847U,  4604886103475043762U,
	 4588115294056142819U,  4607168688050493276U,
	13830540724905269084U,  4588115294056142819U,
	 4607144295058764886U,  4591287158938884897U,
	13814659195793660705U,  4607144295058764886U,
	 4603931940768740167U,  4605102686554936490U,
	13828474723409712298U,  4603931940768740167U,
	 4606144763310860551U,  4602065906208722008U,
	13825437943063497816U,  4606144763310860551U,
	 4599006600037663623U,  4606778366364612594U,
	13830150403219388402U,  4599006600037663623U,
	 4606810452769876110U,  4598795050632330097U,
	13822167087487105905U,  4606810452769876110U,
	 4602260871257280788U,  4606092657816072624U,
	13829464694670848432U,  4602260871257280788U,
	 4605172808754305228U,  4603846496621587377U,
	13827218533476363185U,  4605172808754305228U,
	 4592167175087283203U,  4607133460805585796U,
	13830505497660361604U,  4592167175087283203U,
	 4607019963775302583U,  4595979936813835462U,
	13819351973668611270U,  4607019963775302583U,
	 4603225210076562971U,  4605633586259814045U,
	13829005623114589853U,  4603225210076562971U,
	 4605694995810664660U,  4603133304188877240U,
	13826505341043653048U,  4605694995810664660U,
	 4596413578358834022U,  4606998399608725124U,
	13830370436463500932U,  4596413578358834022U,
	 4606517779747998088U,  4600463181646572228U,
	13823835218501348036U,  4606517779747998088U,
	 4600667422348321968U,  4606475480113671417U,
	13829847516968447225U,  4600667422348321968U,
	 4604583231088591477U,  4604505071555817232U,
	13827877108410593040U,  4604583231088591477U,
	 4573724215515480177U,  4607182249242036882U,
	13830554286096812690U,  4573724215515480177U,
	 4607182376410422530U,  4569220649180767418U,
	13792592686035543226U,  4607182376410422530U,
	 4604524701268679793U,  4604563781218984604U,
	13827935818073760412U,  4604524701268679793U,
	 4606486172460753999U,  4600616459743653188U,
	13823988496598428996U,  4606486172460753999U,
	 4600514338912178239U,  4606507322377452870U,
	13829879359232228678U,  4600514338912178239U,
	 4607003915349878877U,  4596305267720071930U,
	13819677304574847738U,  4607003915349878877U,
	 4603156351203636159U,  4605679749231851918U,
	13829051786086627726U,  4603156351203636159U,
	 4605649044311923410U,  4603202304363743346U,
	13826574341218519154U,  4605649044311923410U,
	 4596088445927168004U,  4607014697483910382U,
	13830386734338686190U,  4596088445927168004U,
	 4607136295912168606U,  4591947271803021404U,
	13815319308657797212U,  4607136295912168606U,
	 4603867938232615808U,  4605155376589456981U,
	13828527413444232789U,  4603867938232615808U,
	 4606105796280968177U,  4602212250118051877U,
	13825584286972827685U,  4606105796280968177U,
	 4598848011564831930U,  4606802552898869248U,
	13830174589753645056U,  4598848011564831930U,
	 4606786509620734768U,  4598953786765296928U,
	13822325823620072736U,  4606786509620734768U,
	 4602114767134999006U,  4606131849150971908U,
	13829503886005747716U,  4602114767134999006U,
	 4605120315324767624U,  4603910660507251362U,
	13827282697362027170U,  4605120315324767624U,
	 4591507261658050721U,  4607141713064252300U,
	13830513749919028108U,  4591507261658050721U,
	 4607170170974224083U,  4587673791460508439U,
	13811045828315284247U,  4607170170974224083U,
	 4604203581176243359U,  4604867640218014515U,
	13828239677072790323U,  4604203581176243359U,
	 4606305777984577632U,  4601423692641949331U,
	13824795729496725139U,  4606305777984577632U,
	 4599688422741010356U,  4606665164148251002U,
	13830037201003026810U,  4599688422741010356U,
	 4606905728766014348U,  4598029484874872834U,
	13821401521729648642U,  4606905728766014348U,
	 4602782121393764535U,  4605915122243179241U,
	13829287159097955049U,  4602782121393764535U,
	 4605393374401988274U,  4603562972219549215U,
	13826935009074325023U,  4605393374401988274U,
	 4594345179472540681U,  4607088942243446236U,
	13830460979098222044U,  4594345179472540681U,
	 4607080832832247697U,  4594563856311064231U,
	13817935893165840039U,  4607080832832247697U,
	 4603518581031047189U,  4605426297151190466U,
	13828798334005966274U,  4603518581031047189U,
	 4605886709123365959U,  4602829525820289164U,
	13826201562675064972U,  4605886709123365959U,
	 4597815040470278984U,  4606919157647773535U,
	13830291194502549343U,  4597815040470278984U,
	 4606646545123403481U,  4599792496117920694U,
	13823164532972696502U,  4606646545123403481U,
	 4601323770373937522U,  4606329407841126011U,
	13829701444695901819U,  4601323770373937522U,
	 4604830524903495634U,  4604244531615310815U,
	13827616568470086623U,  4604830524903495634U,
	 4586790578280679046U,  4607172882816799076U,
	13830544919671574884U,  4586790578280679046U,
	 4607178985458280057U,  4583614727651146525U,
	13806986764505922333U,  4607178985458280057U,
	 4604366005771528720U,  4604717681185626434U,
	13828089718040402242U,  4604366005771528720U,
	 4606398451906509788U,  4601022290077223616U,
	13824394326931999424U,  4606398451906509788U,
	 4600103317933788342U,  4606588777269136769U,
	13829960814123912577U,  4600103317933788342U,
	 4606957467106717424U,  4597169786279785693U,
	13820541823134561501U,  4606957467106717424U,
	 4602970680601913687U,  4605799732098147061U,
	13829171768952922869U,  4602970680601913687U,
	 4605523422498301790U,  4603384207141321914U,
	13826756243996097722U,  4605523422498301790U,
	 4595218635031890910U,  4607054494135176056U,
	13830426530989951864U,  4595218635031890910U,
	 4607111255739239816U,  4593688012422887515U,
	13817060049277663323U,  4607111255739239816U,
	 4603694922063032361U,  4605292980606880364U,
	13828665017461656172U,  4603694922063032361U,
	 4605998608960791335U,  4602598930031891166U,
	13825970966886666974U,  4605998608960791335U,
	 4598423001813699022U,  4606863472012527185U,
	13830235508867302993U,  4598423001813699022U,
	 4606719100629313491U,  4599374859150636784U,
	13822746896005412592U,  4606719100629313491U,
	 4601721693286060937U,  4606233055365547081U,
	13829605092220322889U,  4601721693286060937U,
	 4604977468824438271U,  4604079374282302598U,
	13827451411137078406U,  4604977468824438271U,
	 4589744810590291021U,  4607160003989618959U,
	13830532040844394767U,  4589744810590291021U,
	 4607155938267770208U,  4590185751760970393U,
	13813557788615746201U,  4607155938267770208U,
	 4604037525321326463U,  4605013567986435066U,
	13828385604841210874U,  4604037525321326463U,
	 4606208206518262803U,  4601820425647934753U,
	13825192462502710561U,  4606208206518262803U,
	 4599269903251194481U,  4606736437002195879U,
	13830108473856971687U,  4599269903251194481U,
	 4606848731493011465U,  4598529532600161144U,
	13821901569454936952U,  4606848731493011465U,
	 4602502755147763107U,  4606025850160239809U,
	13829397887015015617U,  4602502755147763107U,
	 4605258978359093269U,  4603738491917026584U,
	13827110528771802392U,  4605258978359093269U,
	 4593265590854265407U,  4607118021058468598U,
	13830490057913244406U,  4593265590854265407U,
	 4607045045516813836U,  4595436449949385485U,
	13818808486804161293U,  4607045045516813836U,
	 4603339021357904144U,  4605555245917486022U,
	13828927282772261830U,  4603339021357904144U,
	 4605770164172969910U,  4603017373458244943U,
	13826389410313020751U,  4605770164172969910U,
	 4596954088216812973U,  4606969576261663845U,
	13830341613116439653U,  4596954088216812973U,
	 4606568886807728474U,  4600206446098256018U,
	13823578482953031826U,  4606568886807728474U,
	 4600921238092511730U,  4606420848538580260U,
	13829792885393356068U,  4600921238092511730U,
	 4604679572075463103U,  4604406033021674239U,
	13827778069876450047U,  4604679572075463103U,
	 4581846703643734566U,  4607180341788068727U,
	13830552378642844535U,  4581846703643734566U,
	 4607181359080094673U,  4579996072175835083U,
	13803368109030610891U,  4607181359080094673U,
	 4604445825685214043U,  4604641218080103285U,
	13828013254934879093U,  4604445825685214043U,
	 4606442934727379583U,  4600819913163773071U,
	13824191950018548879U,  4606442934727379583U,
	 4600309328230211502U,  4606548680329491866U,
	13829920717184267674U,  4600309328230211502U,
	 4606981354314050484U,  4596738097012783531U,
	13820110133867559339U,  4606981354314050484U,
	 4603063884010218172U,  4605740310302420207U,
	13829112347157196015U,  4603063884010218172U,
	 4605586791482848547U,  4603293641160266722U,
	13826665678015042530U,  4605586791482848547U,
	 4595654028864046335U,  4607035262954517034U,
	13830407299809292842U,  4595654028864046335U,
	 4607124449686274900U,  4592826452951465409U,
	13816198489806241217U,  4607124449686274900U,
	 4603781852316960384U,  4605224709411790590U,
	13828596746266566398U,  4603781852316960384U,
	 4606052795787882823U,  4602406247776385022U,
	13825778284631160830U,  4606052795787882823U,
	 4598635880488956483U,  4606833664420673202U,
	13830205701275449010U,  4598635880488956483U,
	 4606753451050079834U,  4599164736579548843U,
	13822536773434324651U,  4606753451050079834U,
	 4601918851211878557U,  4606183055233559255U,
	13829555092088335063U,  4601918851211878557U,
	 4605049409688478101U,  4603995455647851249U,
	13827367492502627057U,  4605049409688478101U,
	 4590626485056654602U,  4607151534426937478U,
	13830523571281713286U,  4590626485056654602U,
	 4607163731439411601U,  4589303678145802340U,
	13812675715000578148U,  4607163731439411601U,
	 4604121000955189926U,  4604941113561600762U,
	13828313150416376570U,  4604121000955189926U,
	 4606257600839867033U,  4601622657843474729U,
	13824994694698250537U,  4606257600839867033U,
	 4599479600326345459U,  4606701442584137310U,
	13830073479438913118U,  4599479600326345459U,
	 4606877885424248132U,  4598316292140394014U,
	13821688328995169822U,  4606877885424248132U,
	 4602686793990243041U,  4605971073215153165U,
	13829343110069928973U,  4602686793990243041U,
	 4605326714874986465U,  4603651144395358093U,
	13827023181250133901U,  4605326714874986465U,
	 4593907249284540294U,  4607104153983298999U,
	13830476190838074807U,  4593907249284540294U,
	 4607063608453868552U,  4595000592312171144U,
	13818372629166946952U,  4607063608453868552U,
	 4603429196809300824U,  4605491322423429598U,
	13828863359278205406U,  4603429196809300824U,
	 4605829012964735987U,  4602923807199184054U,
	13826295844053959862U,  4605829012964735987U,
	 4597385183080791534U,  4606945027305114062U,
	13830317064159889870U,  4597385183080791534U,
	 4606608350964852124U,  4599999947619525579U,
	13823371984474301387U,  4606608350964852124U,
	 4601123065313358619U,  4606375745674388705U,
	13829747782529164513U,  4601123065313358619U,
	 4604755543975806820U,  4604325745441780828U,
	13827697782296556636U,  4604755543975806820U,
	 4585023436363055487U,  4607177290141793710U,
	13830549326996569518U,  4585023436363055487U,
	 4607175255902437396U,  4585907115494236537U,
	13809279152349012345U,  4607175255902437396U,
	 4604285253548209224U,  4604793159020491611U,
	13828165195875267419U,  4604285253548209224U,
	 4606352730697093817U,  4601223560006786057U,
	13824595596861561865U,  4606352730697093817U,
	 4599896339047301634U,  4606627607157935956U,
	13829999644012711764U,  4599896339047301634U,
	 4606932257325205256U,  4597600270510262682U,
	13820972307365038490U,  4606932257325205256U,
	 4602876755014813164U,  4605858005670328613U,
	13829230042525104421U,  4602876755014813164U,
	 4605458946901419122U,  4603473988668005304U,
	13826846025522781112U,  4605458946901419122U,
	 4594782329999411347U,  4607072388129742377U,
	13830444424984518185U,  4594782329999411347U,
	 4607096716058023245U,  4594126307716900071U,
	13817498344571675879U,  4607096716058023245U,
	 4603607160562208225U,  4605360179893335444U,
	13828732216748111252U,  4603607160562208225U,
	 4605943243960030558U,  4602734543519989142U,
	13826106580374764950U,  4605943243960030558U,
	 4598209407597805010U,  4606891971185517504U,
	13830264008040293312U,  4598209407597805010U,
	 4606683463531482757U,  4599584122834874440U,
	13822956159689650248U,  4606683463531482757U,
	 4601523323048804569U,  4606281842017099424U,
	13829653878871875232U,  4601523323048804569U,
	 4604904503566677638U,  4604162403772767740U,
	13827534440627543548U,  4604904503566677638U,
	 4588556721781247689U,  4607167120476811757U,
	13830539157331587565U,  4588556721781247689U,
	 4607146792632922887U,  4591066993883984169U,
	13814439030738759977U,  4607146792632922887U,
	 4603953166845776383U,  4605084992581147553U,
	13828457029435923361U,  4603953166845776383U,
	 4606157602458368090U,  4602016966272225497U,
	13825389003127001305U,  4606157602458368090U,
	 4599059363095165615U,  4606770142132396069U,
	13830142178987171877U,  4599059363095165615U,
	 4606818271362779153U,  4598742041476147134U,
	13822114078330922942U,  4606818271362779153U,
	 4602309411551204896U,  4606079444829232727U,
	13829451481684008535U,  4602309411551204896U,
	 4605190175055178825U,  4603825001630339212U,
	13827197038485115020U,  4605190175055178825U,
	 4592387007752762956U,  4607130541380624519U,
	13830502578235400327U,  4592387007752762956U,
	 4607025146816593591U,  4595871363584150300U,
	13819243400438926108U,  4607025146816593591U,
	 4603248068256948438U,  4605618058006716661U,
	13828990094861492469U,  4603248068256948438U,
	 4605710171610479304U,  4603110210506737381U,
	13826482247361513189U,  4605710171610479304U,
	 4596521820799644122U,  4606992800820440327U,
	13830364837675216135U,  4596521820799644122U,
	 4606528158595189433U,  4600411960456200676U,
	13823783997310976484U,  4606528158595189433U,
	 4600718319105833937U,  4606464709641375231U,
	13829836746496151039U,  4600718319105833937U,
	 4604602620643553229U,  4604485382263976838U,
	13827857419118752646U,  4604602620643553229U,
	 4576459225186735875U,  4607182037296057423U,
	13830554074150833231U,  4576459225186735875U,
	 4607182037296057423U,  4576459225186735875U,
	13799831262041511683U,  4607182037296057423U,
	 4604485382263976838U,  4604602620643553229U,
	13827974657498329037U,  4604485382263976838U,
	 4606464709641375231U,  4600718319105833937U,
	13824090355960609745U,  4606464709641375231U,
	 4600411960456200676U,  4606528158595189433U,
	13829900195449965241U,  4600411960456200676U,
	 4606992800820440327U,  4596521820799644122U,
	13819893857654419930U,  4606992800820440327U,
	 4603110210506737381U,  4605710171610479304U,
	13829082208465255112U,  4603110210506737381U,
	 4605618058006716661U,  4603248068256948438U,
	13826620105111724246U,  4605618058006716661U,
	 4595871363584150300U,  4607025146816593591U,
	13830397183671369399U,  4595871363584150300U,
	 4607130541380624519U,  4592387007752762956U,
	13815759044607538764U,  4607130541380624519U,
	 4603825001630339212U,  4605190175055178825U,
	13828562211909954633U,  4603825001630339212U,
	 4606079444829232727U,  4602309411551204896U,
	13825681448405980704U,  4606079444829232727U,
	 4598742041476147134U,  4606818271362779153U,
	13830190308217554961U,  4598742041476147134U,
	 4606770142132396069U,  4599059363095165615U,
	13822431399949941423U,  4606770142132396069U,
	 4602016966272225497U,  4606157602458368090U,
	13829529639313143898U,  4602016966272225497U,
	 4605084992581147553U,  4603953166845776383U,
	13827325203700552191U,  4605084992581147553U,
	 4591066993883984169U,  4607146792632922887U,
	13830518829487698695U,  4591066993883984169U,
	 4607167120476811757U,  4588556721781247689U,
	13811928758636023497U,  4607167120476811757U,
	 4604162403772767740U,  4604904503566677638U,
	13828276540421453446U,  4604162403772767740U,
	 4606281842017099424U,  4601523323048804569U,
	13824895359903580377U,  4606281842017099424U,
	 4599584122834874440U,  4606683463531482757U,
	13830055500386258565U,  4599584122834874440U,
	 4606891971185517504U,  4598209407597805010U,
	13821581444452580818U,  4606891971185517504U,
	 4602734543519989142U,  4605943243960030558U,
	13829315280814806366U,  4602734543519989142U,
	 4605360179893335444U,  4603607160562208225U,
	13826979197416984033U,  4605360179893335444U,
	 4594126307716900071U,  4607096716058023245U,
	13830468752912799053U,  4594126307716900071U,
	 4607072388129742377U,  4594782329999411347U,
	13818154366854187155U,  4607072388129742377U,
	 4603473988668005304U,  4605458946901419122U,
	13828830983756194930U,  4603473988668005304U,
	 4605858005670328613U,  4602876755014813164U,
	13826248791869588972U,  4605858005670328613U,
	 4597600270510262682U,  4606932257325205256U,
	13830304294179981064U,  4597600270510262682U,
	 4606627607157935956U,  4599896339047301634U,
	13823268375902077442U,  4606627607157935956U,
	 4601223560006786057U,  4606352730697093817U,
	13829724767551869625U,  4601223560006786057U,
	 4604793159020491611U,  4604285253548209224U,
	13827657290402985032U,  4604793159020491611U,
	 4585907115494236537U,  4607175255902437396U,
	13830547292757213204U,  4585907115494236537U,
	 4607177290141793710U,  4585023436363055487U,
	13808395473217831295U,  4607177290141793710U,
	 4604325745441780828U,  4604755543975806820U,
	13828127580830582628U,  4604325745441780828U,
	 4606375745674388705U,  4601123065313358619U,
	13824495102168134427U,  4606375745674388705U,
	 4599999947619525579U,  4606608350964852124U,
	13829980387819627932U,  4599999947619525579U,
	 4606945027305114062U,  4597385183080791534U,
	13820757219935567342U,  4606945027305114062U,
	 4602923807199184054U,  4605829012964735987U,
	13829201049819511795U,  4602923807199184054U,
	 4605491322423429598U,  4603429196809300824U,
	13826801233664076632U,  4605491322423429598U,
	 4595000592312171144U,  4607063608453868552U,
	13830435645308644360U,  4595000592312171144U,
	 4607104153983298999U,  4593907249284540294U,
	13817279286139316102U,  4607104153983298999U,
	 4603651144395358093U,  4605326714874986465U,
	13828698751729762273U,  4603651144395358093U,
	 4605971073215153165U,  4602686793990243041U,
	13826058830845018849U,  4605971073215153165U,
	 4598316292140394014U,  4606877885424248132U,
	13830249922279023940U,  4598316292140394014U,
	 4606701442584137310U,  4599479600326345459U,
	13822851637181121267U,  4606701442584137310U,
	 4601622657843474729U,  4606257600839867033U,
	13829629637694642841U,  4601622657843474729U,
	 4604941113561600762U,  4604121000955189926U,
	13827493037809965734U,  4604941113561600762U,
	 4589303678145802340U,  4607163731439411601U,
	13830535768294187409U,  4589303678145802340U,
	 4607151534426937478U,  4590626485056654602U,
	13813998521911430410U,  4607151534426937478U,
	 4603995455647851249U,  4605049409688478101U,
	13828421446543253909U,  4603995455647851249U,
	 4606183055233559255U,  4601918851211878557U,
	13825290888066654365U,  4606183055233559255U,
	 4599164736579548843U,  4606753451050079834U,
	13830125487904855642U,  4599164736579548843U,
	 4606833664420673202U,  4598635880488956483U,
	13822007917343732291U,  4606833664420673202U,
	 4602406247776385022U,  4606052795787882823U,
	13829424832642658631U,  4602406247776385022U,
	 4605224709411790590U,  4603781852316960384U,
	13827153889171736192U,  4605224709411790590U,
	 4592826452951465409U,  4607124449686274900U,
	13830496486541050708U,  4592826452951465409U,
	 4607035262954517034U,  4595654028864046335U,
	13819026065718822143U,  4607035262954517034U,
	 4603293641160266722U,  4605586791482848547U,
	13828958828337624355U,  4603293641160266722U,
	 4605740310302420207U,  4603063884010218172U,
	13826435920864993980U,  4605740310302420207U,
	 4596738097012783531U,  4606981354314050484U,
	13830353391168826292U,  4596738097012783531U,
	 4606548680329491866U,  4600309328230211502U,
	13823681365084987310U,  4606548680329491866U,
	 4600819913163773071U,  4606442934727379583U,
	13829814971582155391U,  4600819913163773071U,
	 4604641218080103285U,  4604445825685214043U,
	13827817862539989851U,  4604641218080103285U,
	 4579996072175835083U,  4607181359080094673U,
	13830553395934870481U,  4579996072175835083U,
	 4607180341788068727U,  4581846703643734566U,
	13805218740498510374U,  4607180341788068727U,
	 4604406033021674239U,  4604679572075463103U,
	13828051608930238911U,  4604406033021674239U,
	 4606420848538580260U,  4600921238092511730U,
	13824293274947287538U,  4606420848538580260U,
	 4600206446098256018U,  4606568886807728474U,
	13829940923662504282U,  4600206446098256018U,
	 4606969576261663845U,  4596954088216812973U,
	13820326125071588781U,  4606969576261663845U,
	 4603017373458244943U,  4605770164172969910U,
	13829142201027745718U,  4603017373458244943U,
	 4605555245917486022U,  4603339021357904144U,
	13826711058212679952U,  4605555245917486022U,
	 4595436449949385485U,  4607045045516813836U,
	13830417082371589644U,  4595436449949385485U,
	 4607118021058468598U,  4593265590854265407U,
	13816637627709041215U,  4607118021058468598U,
	 4603738491917026584U,  4605258978359093269U,
	13828631015213869077U,  4603738491917026584U,
	 4606025850160239809U,  4602502755147763107U,
	13825874792002538915U,  4606025850160239809U,
	 4598529532600161144U,  4606848731493011465U,
	13830220768347787273U,  4598529532600161144U,
	 4606736437002195879U,  4599269903251194481U,
	13822641940105970289U,  4606736437002195879U,
	 4601820425647934753U,  4606208206518262803U,
	13829580243373038611U,  4601820425647934753U,
	 4605013567986435066U,  4604037525321326463U,
	13827409562176102271U,  4605013567986435066U,
	 4590185751760970393U,  4607155938267770208U,
	13830527975122546016U,  4590185751760970393U,
	 4607160003989618959U,  4589744810590291021U,
	13813116847445066829U,  4607160003989618959U,
	 4604079374282302598U,  4604977468824438271U,
	13828349505679214079U,  4604079374282302598U,
	 4606233055365547081U,  4601721693286060937U,
	13825093730140836745U,  4606233055365547081U,
	 4599374859150636784U,  4606719100629313491U,
	13830091137484089299U,  4599374859150636784U,
	 4606863472012527185U,  4598423001813699022U,
	13821795038668474830U,  4606863472012527185U,
	 4602598930031891166U,  4605998608960791335U,
	13829370645815567143U,  4602598930031891166U,
	 4605292980606880364U,  4603694922063032361U,
	13827066958917808169U,  4605292980606880364U,
	 4593688012422887515U,  4607111255739239816U,
	13830483292594015624U,  4593688012422887515U,
	 4607054494135176056U,  4595218635031890910U,
	13818590671886666718U,  4607054494135176056U,
	 4603384207141321914U,  4605523422498301790U,
	13828895459353077598U,  4603384207141321914U,
	 4605799732098147061U,  4602970680601913687U,
	13826342717456689495U,  4605799732098147061U,
	 4597169786279785693U,  4606957467106717424U,
	13830329503961493232U,  4597169786279785693U,
	 4606588777269136769U,  4600103317933788342U,
	13823475354788564150U,  4606588777269136769U,
	 4601022290077223616U,  4606398451906509788U,
	13829770488761285596U,  4601022290077223616U,
	 4604717681185626434U,  4604366005771528720U,
	13827738042626304528U,  4604717681185626434U,
	 4583614727651146525U,  4607178985458280057U,
	13830551022313055865U,  4583614727651146525U,
	 4607172882816799076U,  4586790578280679046U,
	13810162615135454854U,  4607172882816799076U,
	 4604244531615310815U,  4604830524903495634U,
	13828202561758271442U,  4604244531615310815U,
	 4606329407841126011U,  4601323770373937522U,
	13824695807228713330U,  4606329407841126011U,
	 4599792496117920694U,  4606646545123403481U,
	13830018581978179289U,  4599792496117920694U,
	 4606919157647773535U,  4597815040470278984U,
	13821187077325054792U,  4606919157647773535U,
	 4602829525820289164U,  4605886709123365959U,
	13829258745978141767U,  4602829525820289164U,
	 4605426297151190466U,  4603518581031047189U,
	13826890617885822997U,  4605426297151190466U,
	 4594563856311064231U,  4607080832832247697U,
	13830452869687023505U,  4594563856311064231U,
	 4607088942243446236U,  4594345179472540681U,
	13817717216327316489U,  4607088942243446236U,
	 4603562972219549215U,  4605393374401988274U,
	13828765411256764082U,  4603562972219549215U,
	 4605915122243179241U,  4602782121393764535U,
	13826154158248540343U,  4605915122243179241U,
	 4598029484874872834U,  4606905728766014348U,
	13830277765620790156U,  4598029484874872834U,
	 4606665164148251002U,  4599688422741010356U,
	13823060459595786164U,  4606665164148251002U,
	 4601423692641949331U,  4606305777984577632U,
	13829677814839353440U,  4601423692641949331U,
	 4604867640218014515U,  4604203581176243359U,
	13827575618031019167U,  4604867640218014515U,
	 4587673791460508439U,  4607170170974224083U,
	13830542207828999891U,  4587673791460508439U,
	 4607141713064252300U,  4591507261658050721U,
	13814879298512826529U,  4607141713064252300U,
	 4603910660507251362U,  4605120315324767624U,
	13828492352179543432U,  4603910660507251362U,
	 4606131849150971908U,  4602114767134999006U,
	13825486803989774814U,  4606131849150971908U,
	 4598953786765296928U,  4606786509620734768U,
	13830158546475510576U,  4598953786765296928U,
	 4606802552898869248U,  4598848011564831930U,
	13822220048419607738U,  4606802552898869248U,
	 4602212250118051877U,  4606105796280968177U,
	13829477833135743985U,  4602212250118051877U,
	 4605155376589456981U,  4603867938232615808U,
	13827239975087391616U,  4605155376589456981U,
	 4591947271803021404U,  4607136295912168606U,
	13830508332766944414U,  4591947271803021404U,
	 4607014697483910382U,  4596088445927168004U,
	13819460482781943812U,  4607014697483910382U,
	 4603202304363743346U,  4605649044311923410U,
	13829021081166699218U,  4603202304363743346U,
	 4605679749231851918U,  4603156351203636159U,
	13826528388058411967U,  4605679749231851918U,
	 4596305267720071930U,  4607003915349878877U,
	13830375952204654685U,  4596305267720071930U,
	 4606507322377452870U,  4600514338912178239U,
	13823886375766954047U,  4606507322377452870U,
	 4600616459743653188U,  4606486172460753999U,
	13829858209315529807U,  4600616459743653188U,
	 4604563781218984604U,  4604524701268679793U,
	13827896738123455601U,  4604563781218984604U,
	 4569220649180767418U,  4607182376410422530U,
	13830554413265198338U,  4569220649180767418U
};

const fpr fpr_p2_tab[] = {
	4611686018427387904U,
	4607182418800017408U,
	4602678819172646912U,
	4598175219545276416U,
	4593671619917905920U,
	4589168020290535424U,
	4584664420663164928U,
	4580160821035794432U,
	4575657221408423936U,
	4571153621781053440U,
	4566650022153682944U
};

#elif FALCON_FPNATIVE // yyyFPEMU+0 yyyFPNATIVE+1

const fpr fpr_gm_tab[] = {
	{0}, {0}, /* unused */
	{-0.000000000000000000000000000}, { 1.000000000000000000000000000},
	{ 0.707106781186547524400844362}, { 0.707106781186547524400844362},
	{-0.707106781186547524400844362}, { 0.707106781186547524400844362},
	{ 0.923879532511286756128183189}, { 0.382683432365089771728459984},
	{-0.382683432365089771728459984}, { 0.923879532511286756128183189},
	{ 0.382683432365089771728459984}, { 0.923879532511286756128183189},
	{-0.923879532511286756128183189}, { 0.382683432365089771728459984},
	{ 0.980785280403230449126182236}, { 0.195090322016128267848284868},
	{-0.195090322016128267848284868}, { 0.980785280403230449126182236},
	{ 0.555570233019602224742830814}, { 0.831469612302545237078788378},
	{-0.831469612302545237078788378}, { 0.555570233019602224742830814},
	{ 0.831469612302545237078788378}, { 0.555570233019602224742830814},
	{-0.555570233019602224742830814}, { 0.831469612302545237078788378},
	{ 0.195090322016128267848284868}, { 0.980785280403230449126182236},
	{-0.980785280403230449126182236}, { 0.195090322016128267848284868},
	{ 0.995184726672196886244836953}, { 0.098017140329560601994195564},
	{-0.098017140329560601994195564}, { 0.995184726672196886244836953},
	{ 0.634393284163645498215171613}, { 0.773010453362736960810906610},
	{-0.773010453362736960810906610}, { 0.634393284163645498215171613},
	{ 0.881921264348355029712756864}, { 0.471396736825997648556387626},
	{-0.471396736825997648556387626}, { 0.881921264348355029712756864},
	{ 0.290284677254462367636192376}, { 0.956940335732208864935797887},
	{-0.956940335732208864935797887}, { 0.290284677254462367636192376},
	{ 0.956940335732208864935797887}, { 0.290284677254462367636192376},
	{-0.290284677254462367636192376}, { 0.956940335732208864935797887},
	{ 0.471396736825997648556387626}, { 0.881921264348355029712756864},
	{-0.881921264348355029712756864}, { 0.471396736825997648556387626},
	{ 0.773010453362736960810906610}, { 0.634393284163645498215171613},
	{-0.634393284163645498215171613}, { 0.773010453362736960810906610},
	{ 0.098017140329560601994195564}, { 0.995184726672196886244836953},
	{-0.995184726672196886244836953}, { 0.098017140329560601994195564},
	{ 0.998795456205172392714771605}, { 0.049067674327418014254954977},
	{-0.049067674327418014254954977}, { 0.998795456205172392714771605},
	{ 0.671558954847018400625376850}, { 0.740951125354959091175616897},
	{-0.740951125354959091175616897}, { 0.671558954847018400625376850},
	{ 0.903989293123443331586200297}, { 0.427555093430282094320966857},
	{-0.427555093430282094320966857}, { 0.903989293123443331586200297},
	{ 0.336889853392220050689253213}, { 0.941544065183020778412509403},
	{-0.941544065183020778412509403}, { 0.336889853392220050689253213},
	{ 0.970031253194543992603984207}, { 0.242980179903263889948274162},
	{-0.242980179903263889948274162}, { 0.970031253194543992603984207},
	{ 0.514102744193221726593693839}, { 0.857728610000272069902269984},
	{-0.857728610000272069902269984}, { 0.514102744193221726593693839},
	{ 0.803207531480644909806676513}, { 0.595699304492433343467036529},
	{-0.595699304492433343467036529}, { 0.803207531480644909806676513},
	{ 0.146730474455361751658850130}, { 0.989176509964780973451673738},
	{-0.989176509964780973451673738}, { 0.146730474455361751658850130},
	{ 0.989176509964780973451673738}, { 0.146730474455361751658850130},
	{-0.146730474455361751658850130}, { 0.989176509964780973451673738},
	{ 0.595699304492433343467036529}, { 0.803207531480644909806676513},
	{-0.803207531480644909806676513}, { 0.595699304492433343467036529},
	{ 0.857728610000272069902269984}, { 0.514102744193221726593693839},
	{-0.514102744193221726593693839}, { 0.857728610000272069902269984},
	{ 0.242980179903263889948274162}, { 0.970031253194543992603984207},
	{-0.970031253194543992603984207}, { 0.242980179903263889948274162},
	{ 0.941544065183020778412509403}, { 0.336889853392220050689253213},
	{-0.336889853392220050689253213}, { 0.941544065183020778412509403},
	{ 0.427555093430282094320966857}, { 0.903989293123443331586200297},
	{-0.903989293123443331586200297}, { 0.427555093430282094320966857},
	{ 0.740951125354959091175616897}, { 0.671558954847018400625376850},
	{-0.671558954847018400625376850}, { 0.740951125354959091175616897},
	{ 0.049067674327418014254954977}, { 0.998795456205172392714771605},
	{-0.998795456205172392714771605}, { 0.049067674327418014254954977},
	{ 0.999698818696204220115765650}, { 0.024541228522912288031734529},
	{-0.024541228522912288031734529}, { 0.999698818696204220115765650},
	{ 0.689540544737066924616730630}, { 0.724247082951466920941069243},
	{-0.724247082951466920941069243}, { 0.689540544737066924616730630},
	{ 0.914209755703530654635014829}, { 0.405241314004989870908481306},
	{-0.405241314004989870908481306}, { 0.914209755703530654635014829},
	{ 0.359895036534988148775104572}, { 0.932992798834738887711660256},
	{-0.932992798834738887711660256}, { 0.359895036534988148775104572},
	{ 0.975702130038528544460395766}, { 0.219101240156869797227737547},
	{-0.219101240156869797227737547}, { 0.975702130038528544460395766},
	{ 0.534997619887097210663076905}, { 0.844853565249707073259571205},
	{-0.844853565249707073259571205}, { 0.534997619887097210663076905},
	{ 0.817584813151583696504920884}, { 0.575808191417845300745972454},
	{-0.575808191417845300745972454}, { 0.817584813151583696504920884},
	{ 0.170961888760301226363642357}, { 0.985277642388941244774018433},
	{-0.985277642388941244774018433}, { 0.170961888760301226363642357},
	{ 0.992479534598709998156767252}, { 0.122410675199216198498704474},
	{-0.122410675199216198498704474}, { 0.992479534598709998156767252},
	{ 0.615231590580626845484913563}, { 0.788346427626606262009164705},
	{-0.788346427626606262009164705}, { 0.615231590580626845484913563},
	{ 0.870086991108711418652292404}, { 0.492898192229784036873026689},
	{-0.492898192229784036873026689}, { 0.870086991108711418652292404},
	{ 0.266712757474898386325286515}, { 0.963776065795439866686464356},
	{-0.963776065795439866686464356}, { 0.266712757474898386325286515},
	{ 0.949528180593036667195936074}, { 0.313681740398891476656478846},
	{-0.313681740398891476656478846}, { 0.949528180593036667195936074},
	{ 0.449611329654606600046294579}, { 0.893224301195515320342416447},
	{-0.893224301195515320342416447}, { 0.449611329654606600046294579},
	{ 0.757208846506484547575464054}, { 0.653172842953776764084203014},
	{-0.653172842953776764084203014}, { 0.757208846506484547575464054},
	{ 0.073564563599667423529465622}, { 0.997290456678690216135597140},
	{-0.997290456678690216135597140}, { 0.073564563599667423529465622},
	{ 0.997290456678690216135597140}, { 0.073564563599667423529465622},
	{-0.073564563599667423529465622}, { 0.997290456678690216135597140},
	{ 0.653172842953776764084203014}, { 0.757208846506484547575464054},
	{-0.757208846506484547575464054}, { 0.653172842953776764084203014},
	{ 0.893224301195515320342416447}, { 0.449611329654606600046294579},
	{-0.449611329654606600046294579}, { 0.893224301195515320342416447},
	{ 0.313681740398891476656478846}, { 0.949528180593036667195936074},
	{-0.949528180593036667195936074}, { 0.313681740398891476656478846},
	{ 0.963776065795439866686464356}, { 0.266712757474898386325286515},
	{-0.266712757474898386325286515}, { 0.963776065795439866686464356},
	{ 0.492898192229784036873026689}, { 0.870086991108711418652292404},
	{-0.870086991108711418652292404}, { 0.492898192229784036873026689},
	{ 0.788346427626606262009164705}, { 0.615231590580626845484913563},
	{-0.615231590580626845484913563}, { 0.788346427626606262009164705},
	{ 0.122410675199216198498704474}, { 0.992479534598709998156767252},
	{-0.992479534598709998156767252}, { 0.122410675199216198498704474},
	{ 0.985277642388941244774018433}, { 0.170961888760301226363642357},
	{-0.170961888760301226363642357}, { 0.985277642388941244774018433},
	{ 0.575808191417845300745972454}, { 0.817584813151583696504920884},
	{-0.817584813151583696504920884}, { 0.575808191417845300745972454},
	{ 0.844853565249707073259571205}, { 0.534997619887097210663076905},
	{-0.534997619887097210663076905}, { 0.844853565249707073259571205},
	{ 0.219101240156869797227737547}, { 0.975702130038528544460395766},
	{-0.975702130038528544460395766}, { 0.219101240156869797227737547},
	{ 0.932992798834738887711660256}, { 0.359895036534988148775104572},
	{-0.359895036534988148775104572}, { 0.932992798834738887711660256},
	{ 0.405241314004989870908481306}, { 0.914209755703530654635014829},
	{-0.914209755703530654635014829}, { 0.405241314004989870908481306},
	{ 0.724247082951466920941069243}, { 0.689540544737066924616730630},
	{-0.689540544737066924616730630}, { 0.724247082951466920941069243},
	{ 0.024541228522912288031734529}, { 0.999698818696204220115765650},
	{-0.999698818696204220115765650}, { 0.024541228522912288031734529},
	{ 0.999924701839144540921646491}, { 0.012271538285719926079408262},
	{-0.012271538285719926079408262}, { 0.999924701839144540921646491},
	{ 0.698376249408972853554813503}, { 0.715730825283818654125532623},
	{-0.715730825283818654125532623}, { 0.698376249408972853554813503},
	{ 0.919113851690057743908477789}, { 0.393992040061048108596188661},
	{-0.393992040061048108596188661}, { 0.919113851690057743908477789},
	{ 0.371317193951837543411934967}, { 0.928506080473215565937167396},
	{-0.928506080473215565937167396}, { 0.371317193951837543411934967},
	{ 0.978317370719627633106240097}, { 0.207111376192218549708116020},
	{-0.207111376192218549708116020}, { 0.978317370719627633106240097},
	{ 0.545324988422046422313987347}, { 0.838224705554838043186996856},
	{-0.838224705554838043186996856}, { 0.545324988422046422313987347},
	{ 0.824589302785025264474803737}, { 0.565731810783613197389765011},
	{-0.565731810783613197389765011}, { 0.824589302785025264474803737},
	{ 0.183039887955140958516532578}, { 0.983105487431216327180301155},
	{-0.983105487431216327180301155}, { 0.183039887955140958516532578},
	{ 0.993906970002356041546922813}, { 0.110222207293883058807899140},
	{-0.110222207293883058807899140}, { 0.993906970002356041546922813},
	{ 0.624859488142386377084072816}, { 0.780737228572094478301588484},
	{-0.780737228572094478301588484}, { 0.624859488142386377084072816},
	{ 0.876070094195406607095844268}, { 0.482183772079122748517344481},
	{-0.482183772079122748517344481}, { 0.876070094195406607095844268},
	{ 0.278519689385053105207848526}, { 0.960430519415565811199035138},
	{-0.960430519415565811199035138}, { 0.278519689385053105207848526},
	{ 0.953306040354193836916740383}, { 0.302005949319228067003463232},
	{-0.302005949319228067003463232}, { 0.953306040354193836916740383},
	{ 0.460538710958240023633181487}, { 0.887639620402853947760181617},
	{-0.887639620402853947760181617}, { 0.460538710958240023633181487},
	{ 0.765167265622458925888815999}, { 0.643831542889791465068086063},
	{-0.643831542889791465068086063}, { 0.765167265622458925888815999},
	{ 0.085797312344439890461556332}, { 0.996312612182778012627226190},
	{-0.996312612182778012627226190}, { 0.085797312344439890461556332},
	{ 0.998118112900149207125155861}, { 0.061320736302208577782614593},
	{-0.061320736302208577782614593}, { 0.998118112900149207125155861},
	{ 0.662415777590171761113069817}, { 0.749136394523459325469203257},
	{-0.749136394523459325469203257}, { 0.662415777590171761113069817},
	{ 0.898674465693953843041976744}, { 0.438616238538527637647025738},
	{-0.438616238538527637647025738}, { 0.898674465693953843041976744},
	{ 0.325310292162262934135954708}, { 0.945607325380521325730945387},
	{-0.945607325380521325730945387}, { 0.325310292162262934135954708},
	{ 0.966976471044852109087220226}, { 0.254865659604514571553980779},
	{-0.254865659604514571553980779}, { 0.966976471044852109087220226},
	{ 0.503538383725717558691867071}, { 0.863972856121586737918147054},
	{-0.863972856121586737918147054}, { 0.503538383725717558691867071},
	{ 0.795836904608883536262791915}, { 0.605511041404325513920626941},
	{-0.605511041404325513920626941}, { 0.795836904608883536262791915},
	{ 0.134580708507126186316358409}, { 0.990902635427780025108237011},
	{-0.990902635427780025108237011}, { 0.134580708507126186316358409},
	{ 0.987301418157858382399815802}, { 0.158858143333861441684385360},
	{-0.158858143333861441684385360}, { 0.987301418157858382399815802},
	{ 0.585797857456438860328080838}, { 0.810457198252594791726703434},
	{-0.810457198252594791726703434}, { 0.585797857456438860328080838},
	{ 0.851355193105265142261290312}, { 0.524589682678468906215098464},
	{-0.524589682678468906215098464}, { 0.851355193105265142261290312},
	{ 0.231058108280671119643236018}, { 0.972939952205560145467720114},
	{-0.972939952205560145467720114}, { 0.231058108280671119643236018},
	{ 0.937339011912574923201899593}, { 0.348418680249434568419308588},
	{-0.348418680249434568419308588}, { 0.937339011912574923201899593},
	{ 0.416429560097637182562598911}, { 0.909167983090522376563884788},
	{-0.909167983090522376563884788}, { 0.416429560097637182562598911},
	{ 0.732654271672412834615546649}, { 0.680600997795453050594430464},
	{-0.680600997795453050594430464}, { 0.732654271672412834615546649},
	{ 0.036807222941358832324332691}, { 0.999322384588349500896221011},
	{-0.999322384588349500896221011}, { 0.036807222941358832324332691},
	{ 0.999322384588349500896221011}, { 0.036807222941358832324332691},
	{-0.036807222941358832324332691}, { 0.999322384588349500896221011},
	{ 0.680600997795453050594430464}, { 0.732654271672412834615546649},
	{-0.732654271672412834615546649}, { 0.680600997795453050594430464},
	{ 0.909167983090522376563884788}, { 0.416429560097637182562598911},
	{-0.416429560097637182562598911}, { 0.909167983090522376563884788},
	{ 0.348418680249434568419308588}, { 0.937339011912574923201899593},
	{-0.937339011912574923201899593}, { 0.348418680249434568419308588},
	{ 0.972939952205560145467720114}, { 0.231058108280671119643236018},
	{-0.231058108280671119643236018}, { 0.972939952205560145467720114},
	{ 0.524589682678468906215098464}, { 0.851355193105265142261290312},
	{-0.851355193105265142261290312}, { 0.524589682678468906215098464},
	{ 0.810457198252594791726703434}, { 0.585797857456438860328080838},
	{-0.585797857456438860328080838}, { 0.810457198252594791726703434},
	{ 0.158858143333861441684385360}, { 0.987301418157858382399815802},
	{-0.987301418157858382399815802}, { 0.158858143333861441684385360},
	{ 0.990902635427780025108237011}, { 0.134580708507126186316358409},
	{-0.134580708507126186316358409}, { 0.990902635427780025108237011},
	{ 0.605511041404325513920626941}, { 0.795836904608883536262791915},
	{-0.795836904608883536262791915}, { 0.605511041404325513920626941},
	{ 0.863972856121586737918147054}, { 0.503538383725717558691867071},
	{-0.503538383725717558691867071}, { 0.863972856121586737918147054},
	{ 0.254865659604514571553980779}, { 0.966976471044852109087220226},
	{-0.966976471044852109087220226}, { 0.254865659604514571553980779},
	{ 0.945607325380521325730945387}, { 0.325310292162262934135954708},
	{-0.325310292162262934135954708}, { 0.945607325380521325730945387},
	{ 0.438616238538527637647025738}, { 0.898674465693953843041976744},
	{-0.898674465693953843041976744}, { 0.438616238538527637647025738},
	{ 0.749136394523459325469203257}, { 0.662415777590171761113069817},
	{-0.662415777590171761113069817}, { 0.749136394523459325469203257},
	{ 0.061320736302208577782614593}, { 0.998118112900149207125155861},
	{-0.998118112900149207125155861}, { 0.061320736302208577782614593},
	{ 0.996312612182778012627226190}, { 0.085797312344439890461556332},
	{-0.085797312344439890461556332}, { 0.996312612182778012627226190},
	{ 0.643831542889791465068086063}, { 0.765167265622458925888815999},
	{-0.765167265622458925888815999}, { 0.643831542889791465068086063},
	{ 0.887639620402853947760181617}, { 0.460538710958240023633181487},
	{-0.460538710958240023633181487}, { 0.887639620402853947760181617},
	{ 0.302005949319228067003463232}, { 0.953306040354193836916740383},
	{-0.953306040354193836916740383}, { 0.302005949319228067003463232},
	{ 0.960430519415565811199035138}, { 0.278519689385053105207848526},
	{-0.278519689385053105207848526}, { 0.960430519415565811199035138},
	{ 0.482183772079122748517344481}, { 0.876070094195406607095844268},
	{-0.876070094195406607095844268}, { 0.482183772079122748517344481},
	{ 0.780737228572094478301588484}, { 0.624859488142386377084072816},
	{-0.624859488142386377084072816}, { 0.780737228572094478301588484},
	{ 0.110222207293883058807899140}, { 0.993906970002356041546922813},
	{-0.993906970002356041546922813}, { 0.110222207293883058807899140},
	{ 0.983105487431216327180301155}, { 0.183039887955140958516532578},
	{-0.183039887955140958516532578}, { 0.983105487431216327180301155},
	{ 0.565731810783613197389765011}, { 0.824589302785025264474803737},
	{-0.824589302785025264474803737}, { 0.565731810783613197389765011},
	{ 0.838224705554838043186996856}, { 0.545324988422046422313987347},
	{-0.545324988422046422313987347}, { 0.838224705554838043186996856},
	{ 0.207111376192218549708116020}, { 0.978317370719627633106240097},
	{-0.978317370719627633106240097}, { 0.207111376192218549708116020},
	{ 0.928506080473215565937167396}, { 0.371317193951837543411934967},
	{-0.371317193951837543411934967}, { 0.928506080473215565937167396},
	{ 0.393992040061048108596188661}, { 0.919113851690057743908477789},
	{-0.919113851690057743908477789}, { 0.393992040061048108596188661},
	{ 0.715730825283818654125532623}, { 0.698376249408972853554813503},
	{-0.698376249408972853554813503}, { 0.715730825283818654125532623},
	{ 0.012271538285719926079408262}, { 0.999924701839144540921646491},
	{-0.999924701839144540921646491}, { 0.012271538285719926079408262},
	{ 0.999981175282601142656990438}, { 0.006135884649154475359640235},
	{-0.006135884649154475359640235}, { 0.999981175282601142656990438},
	{ 0.702754744457225302452914421}, { 0.711432195745216441522130290},
	{-0.711432195745216441522130290}, { 0.702754744457225302452914421},
	{ 0.921514039342041943465396332}, { 0.388345046698826291624993541},
	{-0.388345046698826291624993541}, { 0.921514039342041943465396332},
	{ 0.377007410216418256726567823}, { 0.926210242138311341974793388},
	{-0.926210242138311341974793388}, { 0.377007410216418256726567823},
	{ 0.979569765685440534439326110}, { 0.201104634842091911558443546},
	{-0.201104634842091911558443546}, { 0.979569765685440534439326110},
	{ 0.550457972936604802977289893}, { 0.834862874986380056304401383},
	{-0.834862874986380056304401383}, { 0.550457972936604802977289893},
	{ 0.828045045257755752067527592}, { 0.560661576197336023839710223},
	{-0.560661576197336023839710223}, { 0.828045045257755752067527592},
	{ 0.189068664149806212754997837}, { 0.981963869109555264072848154},
	{-0.981963869109555264072848154}, { 0.189068664149806212754997837},
	{ 0.994564570734255452119106243}, { 0.104121633872054579120943880},
	{-0.104121633872054579120943880}, { 0.994564570734255452119106243},
	{ 0.629638238914927025372981341}, { 0.776888465673232450040827983},
	{-0.776888465673232450040827983}, { 0.629638238914927025372981341},
	{ 0.879012226428633477831323711}, { 0.476799230063322133342158117},
	{-0.476799230063322133342158117}, { 0.879012226428633477831323711},
	{ 0.284407537211271843618310615}, { 0.958703474895871555374645792},
	{-0.958703474895871555374645792}, { 0.284407537211271843618310615},
	{ 0.955141168305770721498157712}, { 0.296150888243623824121786128},
	{-0.296150888243623824121786128}, { 0.955141168305770721498157712},
	{ 0.465976495767966177902756065}, { 0.884797098430937780104007041},
	{-0.884797098430937780104007041}, { 0.465976495767966177902756065},
	{ 0.769103337645579639346626069}, { 0.639124444863775743801488193},
	{-0.639124444863775743801488193}, { 0.769103337645579639346626069},
	{ 0.091908956497132728624990979}, { 0.995767414467659793982495643},
	{-0.995767414467659793982495643}, { 0.091908956497132728624990979},
	{ 0.998475580573294752208559038}, { 0.055195244349689939809447526},
	{-0.055195244349689939809447526}, { 0.998475580573294752208559038},
	{ 0.666999922303637506650154222}, { 0.745057785441465962407907310},
	{-0.745057785441465962407907310}, { 0.666999922303637506650154222},
	{ 0.901348847046022014570746093}, { 0.433093818853151968484222638},
	{-0.433093818853151968484222638}, { 0.901348847046022014570746093},
	{ 0.331106305759876401737190737}, { 0.943593458161960361495301445},
	{-0.943593458161960361495301445}, { 0.331106305759876401737190737},
	{ 0.968522094274417316221088329}, { 0.248927605745720168110682816},
	{-0.248927605745720168110682816}, { 0.968522094274417316221088329},
	{ 0.508830142543107036931749324}, { 0.860866938637767279344583877},
	{-0.860866938637767279344583877}, { 0.508830142543107036931749324},
	{ 0.799537269107905033500246232}, { 0.600616479383868926653875896},
	{-0.600616479383868926653875896}, { 0.799537269107905033500246232},
	{ 0.140658239332849230714788846}, { 0.990058210262297105505906464},
	{-0.990058210262297105505906464}, { 0.140658239332849230714788846},
	{ 0.988257567730749491404792538}, { 0.152797185258443427720336613},
	{-0.152797185258443427720336613}, { 0.988257567730749491404792538},
	{ 0.590759701858874228423887908}, { 0.806847553543799272206514313},
	{-0.806847553543799272206514313}, { 0.590759701858874228423887908},
	{ 0.854557988365400520767862276}, { 0.519355990165589587361829932},
	{-0.519355990165589587361829932}, { 0.854557988365400520767862276},
	{ 0.237023605994367206867735915}, { 0.971503890986251775537099622},
	{-0.971503890986251775537099622}, { 0.237023605994367206867735915},
	{ 0.939459223602189911962669246}, { 0.342660717311994397592781983},
	{-0.342660717311994397592781983}, { 0.939459223602189911962669246},
	{ 0.422000270799799685941287941}, { 0.906595704514915365332960588},
	{-0.906595704514915365332960588}, { 0.422000270799799685941287941},
	{ 0.736816568877369875090132520}, { 0.676092703575315960360419228},
	{-0.676092703575315960360419228}, { 0.736816568877369875090132520},
	{ 0.042938256934940823077124540}, { 0.999077727752645382888781997},
	{-0.999077727752645382888781997}, { 0.042938256934940823077124540},
	{ 0.999529417501093163079703322}, { 0.030674803176636625934021028},
	{-0.030674803176636625934021028}, { 0.999529417501093163079703322},
	{ 0.685083667772700381362052545}, { 0.728464390448225196492035438},
	{-0.728464390448225196492035438}, { 0.685083667772700381362052545},
	{ 0.911706032005429851404397325}, { 0.410843171057903942183466675},
	{-0.410843171057903942183466675}, { 0.911706032005429851404397325},
	{ 0.354163525420490382357395796}, { 0.935183509938947577642207480},
	{-0.935183509938947577642207480}, { 0.354163525420490382357395796},
	{ 0.974339382785575860518721668}, { 0.225083911359792835991642120},
	{-0.225083911359792835991642120}, { 0.974339382785575860518721668},
	{ 0.529803624686294668216054671}, { 0.848120344803297251279133563},
	{-0.848120344803297251279133563}, { 0.529803624686294668216054671},
	{ 0.814036329705948361654516690}, { 0.580813958095764545075595272},
	{-0.580813958095764545075595272}, { 0.814036329705948361654516690},
	{ 0.164913120489969921418189113}, { 0.986308097244598647863297524},
	{-0.986308097244598647863297524}, { 0.164913120489969921418189113},
	{ 0.991709753669099522860049931}, { 0.128498110793793172624415589},
	{-0.128498110793793172624415589}, { 0.991709753669099522860049931},
	{ 0.610382806276309452716352152}, { 0.792106577300212351782342879},
	{-0.792106577300212351782342879}, { 0.610382806276309452716352152},
	{ 0.867046245515692651480195629}, { 0.498227666972781852410983869},
	{-0.498227666972781852410983869}, { 0.867046245515692651480195629},
	{ 0.260794117915275518280186509}, { 0.965394441697689374550843858},
	{-0.965394441697689374550843858}, { 0.260794117915275518280186509},
	{ 0.947585591017741134653387321}, { 0.319502030816015677901518272},
	{-0.319502030816015677901518272}, { 0.947585591017741134653387321},
	{ 0.444122144570429231642069418}, { 0.895966249756185155914560282},
	{-0.895966249756185155914560282}, { 0.444122144570429231642069418},
	{ 0.753186799043612482483430486}, { 0.657806693297078656931182264},
	{-0.657806693297078656931182264}, { 0.753186799043612482483430486},
	{ 0.067443919563664057897972422}, { 0.997723066644191609848546728},
	{-0.997723066644191609848546728}, { 0.067443919563664057897972422},
	{ 0.996820299291165714972629398}, { 0.079682437971430121147120656},
	{-0.079682437971430121147120656}, { 0.996820299291165714972629398},
	{ 0.648514401022112445084560551}, { 0.761202385484261814029709836},
	{-0.761202385484261814029709836}, { 0.648514401022112445084560551},
	{ 0.890448723244757889952150560}, { 0.455083587126343823535869268},
	{-0.455083587126343823535869268}, { 0.890448723244757889952150560},
	{ 0.307849640041534893682063646}, { 0.951435020969008369549175569},
	{-0.951435020969008369549175569}, { 0.307849640041534893682063646},
	{ 0.962121404269041595429604316}, { 0.272621355449948984493347477},
	{-0.272621355449948984493347477}, { 0.962121404269041595429604316},
	{ 0.487550160148435954641485027}, { 0.873094978418290098636085973},
	{-0.873094978418290098636085973}, { 0.487550160148435954641485027},
	{ 0.784556597155575233023892575}, { 0.620057211763289178646268191},
	{-0.620057211763289178646268191}, { 0.784556597155575233023892575},
	{ 0.116318630911904767252544319}, { 0.993211949234794533104601012},
	{-0.993211949234794533104601012}, { 0.116318630911904767252544319},
	{ 0.984210092386929073193874387}, { 0.177004220412148756196839844},
	{-0.177004220412148756196839844}, { 0.984210092386929073193874387},
	{ 0.570780745886967280232652864}, { 0.821102514991104679060430820},
	{-0.821102514991104679060430820}, { 0.570780745886967280232652864},
	{ 0.841554977436898409603499520}, { 0.540171472729892881297845480},
	{-0.540171472729892881297845480}, { 0.841554977436898409603499520},
	{ 0.213110319916091373967757518}, { 0.977028142657754351485866211},
	{-0.977028142657754351485866211}, { 0.213110319916091373967757518},
	{ 0.930766961078983731944872340}, { 0.365612997804773870011745909},
	{-0.365612997804773870011745909}, { 0.930766961078983731944872340},
	{ 0.399624199845646828544117031}, { 0.916679059921042663116457013},
	{-0.916679059921042663116457013}, { 0.399624199845646828544117031},
	{ 0.720002507961381629076682999}, { 0.693971460889654009003734389},
	{-0.693971460889654009003734389}, { 0.720002507961381629076682999},
	{ 0.018406729905804820927366313}, { 0.999830581795823422015722275},
	{-0.999830581795823422015722275}, { 0.018406729905804820927366313},
	{ 0.999830581795823422015722275}, { 0.018406729905804820927366313},
	{-0.018406729905804820927366313}, { 0.999830581795823422015722275},
	{ 0.693971460889654009003734389}, { 0.720002507961381629076682999},
	{-0.720002507961381629076682999}, { 0.693971460889654009003734389},
	{ 0.916679059921042663116457013}, { 0.399624199845646828544117031},
	{-0.399624199845646828544117031}, { 0.916679059921042663116457013},
	{ 0.365612997804773870011745909}, { 0.930766961078983731944872340},
	{-0.930766961078983731944872340}, { 0.365612997804773870011745909},
	{ 0.977028142657754351485866211}, { 0.213110319916091373967757518},
	{-0.213110319916091373967757518}, { 0.977028142657754351485866211},
	{ 0.540171472729892881297845480}, { 0.841554977436898409603499520},
	{-0.841554977436898409603499520}, { 0.540171472729892881297845480},
	{ 0.821102514991104679060430820}, { 0.570780745886967280232652864},
	{-0.570780745886967280232652864}, { 0.821102514991104679060430820},
	{ 0.177004220412148756196839844}, { 0.984210092386929073193874387},
	{-0.984210092386929073193874387}, { 0.177004220412148756196839844},
	{ 0.993211949234794533104601012}, { 0.116318630911904767252544319},
	{-0.116318630911904767252544319}, { 0.993211949234794533104601012},
	{ 0.620057211763289178646268191}, { 0.784556597155575233023892575},
	{-0.784556597155575233023892575}, { 0.620057211763289178646268191},
	{ 0.873094978418290098636085973}, { 0.487550160148435954641485027},
	{-0.487550160148435954641485027}, { 0.873094978418290098636085973},
	{ 0.272621355449948984493347477}, { 0.962121404269041595429604316},
	{-0.962121404269041595429604316}, { 0.272621355449948984493347477},
	{ 0.951435020969008369549175569}, { 0.307849640041534893682063646},
	{-0.307849640041534893682063646}, { 0.951435020969008369549175569},
	{ 0.455083587126343823535869268}, { 0.890448723244757889952150560},
	{-0.890448723244757889952150560}, { 0.455083587126343823535869268},
	{ 0.761202385484261814029709836}, { 0.648514401022112445084560551},
	{-0.648514401022112445084560551}, { 0.761202385484261814029709836},
	{ 0.079682437971430121147120656}, { 0.996820299291165714972629398},
	{-0.996820299291165714972629398}, { 0.079682437971430121147120656},
	{ 0.997723066644191609848546728}, { 0.067443919563664057897972422},
	{-0.067443919563664057897972422}, { 0.997723066644191609848546728},
	{ 0.657806693297078656931182264}, { 0.753186799043612482483430486},
	{-0.753186799043612482483430486}, { 0.657806693297078656931182264},
	{ 0.895966249756185155914560282}, { 0.444122144570429231642069418},
	{-0.444122144570429231642069418}, { 0.895966249756185155914560282},
	{ 0.319502030816015677901518272}, { 0.947585591017741134653387321},
	{-0.947585591017741134653387321}, { 0.319502030816015677901518272},
	{ 0.965394441697689374550843858}, { 0.260794117915275518280186509},
	{-0.260794117915275518280186509}, { 0.965394441697689374550843858},
	{ 0.498227666972781852410983869}, { 0.867046245515692651480195629},
	{-0.867046245515692651480195629}, { 0.498227666972781852410983869},
	{ 0.792106577300212351782342879}, { 0.610382806276309452716352152},
	{-0.610382806276309452716352152}, { 0.792106577300212351782342879},
	{ 0.128498110793793172624415589}, { 0.991709753669099522860049931},
	{-0.991709753669099522860049931}, { 0.128498110793793172624415589},
	{ 0.986308097244598647863297524}, { 0.164913120489969921418189113},
	{-0.164913120489969921418189113}, { 0.986308097244598647863297524},
	{ 0.580813958095764545075595272}, { 0.814036329705948361654516690},
	{-0.814036329705948361654516690}, { 0.580813958095764545075595272},
	{ 0.848120344803297251279133563}, { 0.529803624686294668216054671},
	{-0.529803624686294668216054671}, { 0.848120344803297251279133563},
	{ 0.225083911359792835991642120}, { 0.974339382785575860518721668},
	{-0.974339382785575860518721668}, { 0.225083911359792835991642120},
	{ 0.935183509938947577642207480}, { 0.354163525420490382357395796},
	{-0.354163525420490382357395796}, { 0.935183509938947577642207480},
	{ 0.410843171057903942183466675}, { 0.911706032005429851404397325},
	{-0.911706032005429851404397325}, { 0.410843171057903942183466675},
	{ 0.728464390448225196492035438}, { 0.685083667772700381362052545},
	{-0.685083667772700381362052545}, { 0.728464390448225196492035438},
	{ 0.030674803176636625934021028}, { 0.999529417501093163079703322},
	{-0.999529417501093163079703322}, { 0.030674803176636625934021028},
	{ 0.999077727752645382888781997}, { 0.042938256934940823077124540},
	{-0.042938256934940823077124540}, { 0.999077727752645382888781997},
	{ 0.676092703575315960360419228}, { 0.736816568877369875090132520},
	{-0.736816568877369875090132520}, { 0.676092703575315960360419228},
	{ 0.906595704514915365332960588}, { 0.422000270799799685941287941},
	{-0.422000270799799685941287941}, { 0.906595704514915365332960588},
	{ 0.342660717311994397592781983}, { 0.939459223602189911962669246},
	{-0.939459223602189911962669246}, { 0.342660717311994397592781983},
	{ 0.971503890986251775537099622}, { 0.237023605994367206867735915},
	{-0.237023605994367206867735915}, { 0.971503890986251775537099622},
	{ 0.519355990165589587361829932}, { 0.854557988365400520767862276},
	{-0.854557988365400520767862276}, { 0.519355990165589587361829932},
	{ 0.806847553543799272206514313}, { 0.590759701858874228423887908},
	{-0.590759701858874228423887908}, { 0.806847553543799272206514313},
	{ 0.152797185258443427720336613}, { 0.988257567730749491404792538},
	{-0.988257567730749491404792538}, { 0.152797185258443427720336613},
	{ 0.990058210262297105505906464}, { 0.140658239332849230714788846},
	{-0.140658239332849230714788846}, { 0.990058210262297105505906464},
	{ 0.600616479383868926653875896}, { 0.799537269107905033500246232},
	{-0.799537269107905033500246232}, { 0.600616479383868926653875896},
	{ 0.860866938637767279344583877}, { 0.508830142543107036931749324},
	{-0.508830142543107036931749324}, { 0.860866938637767279344583877},
	{ 0.248927605745720168110682816}, { 0.968522094274417316221088329},
	{-0.968522094274417316221088329}, { 0.248927605745720168110682816},
	{ 0.943593458161960361495301445}, { 0.331106305759876401737190737},
	{-0.331106305759876401737190737}, { 0.943593458161960361495301445},
	{ 0.433093818853151968484222638}, { 0.901348847046022014570746093},
	{-0.901348847046022014570746093}, { 0.433093818853151968484222638},
	{ 0.745057785441465962407907310}, { 0.666999922303637506650154222},
	{-0.666999922303637506650154222}, { 0.745057785441465962407907310},
	{ 0.055195244349689939809447526}, { 0.998475580573294752208559038},
	{-0.998475580573294752208559038}, { 0.055195244349689939809447526},
	{ 0.995767414467659793982495643}, { 0.091908956497132728624990979},
	{-0.091908956497132728624990979}, { 0.995767414467659793982495643},
	{ 0.639124444863775743801488193}, { 0.769103337645579639346626069},
	{-0.769103337645579639346626069}, { 0.639124444863775743801488193},
	{ 0.884797098430937780104007041}, { 0.465976495767966177902756065},
	{-0.465976495767966177902756065}, { 0.884797098430937780104007041},
	{ 0.296150888243623824121786128}, { 0.955141168305770721498157712},
	{-0.955141168305770721498157712}, { 0.296150888243623824121786128},
	{ 0.958703474895871555374645792}, { 0.284407537211271843618310615},
	{-0.284407537211271843618310615}, { 0.958703474895871555374645792},
	{ 0.476799230063322133342158117}, { 0.879012226428633477831323711},
	{-0.879012226428633477831323711}, { 0.476799230063322133342158117},
	{ 0.776888465673232450040827983}, { 0.629638238914927025372981341},
	{-0.629638238914927025372981341}, { 0.776888465673232450040827983},
	{ 0.104121633872054579120943880}, { 0.994564570734255452119106243},
	{-0.994564570734255452119106243}, { 0.104121633872054579120943880},
	{ 0.981963869109555264072848154}, { 0.189068664149806212754997837},
	{-0.189068664149806212754997837}, { 0.981963869109555264072848154},
	{ 0.560661576197336023839710223}, { 0.828045045257755752067527592},
	{-0.828045045257755752067527592}, { 0.560661576197336023839710223},
	{ 0.834862874986380056304401383}, { 0.550457972936604802977289893},
	{-0.550457972936604802977289893}, { 0.834862874986380056304401383},
	{ 0.201104634842091911558443546}, { 0.979569765685440534439326110},
	{-0.979569765685440534439326110}, { 0.201104634842091911558443546},
	{ 0.926210242138311341974793388}, { 0.377007410216418256726567823},
	{-0.377007410216418256726567823}, { 0.926210242138311341974793388},
	{ 0.388345046698826291624993541}, { 0.921514039342041943465396332},
	{-0.921514039342041943465396332}, { 0.388345046698826291624993541},
	{ 0.711432195745216441522130290}, { 0.702754744457225302452914421},
	{-0.702754744457225302452914421}, { 0.711432195745216441522130290},
	{ 0.006135884649154475359640235}, { 0.999981175282601142656990438},
	{-0.999981175282601142656990438}, { 0.006135884649154475359640235},
	{ 0.999995293809576171511580126}, { 0.003067956762965976270145365},
	{-0.003067956762965976270145365}, { 0.999995293809576171511580126},
	{ 0.704934080375904908852523758}, { 0.709272826438865651316533772},
	{-0.709272826438865651316533772}, { 0.704934080375904908852523758},
	{ 0.922701128333878570437264227}, { 0.385516053843918864075607949},
	{-0.385516053843918864075607949}, { 0.922701128333878570437264227},
	{ 0.379847208924051170576281147}, { 0.925049240782677590302371869},
	{-0.925049240782677590302371869}, { 0.379847208924051170576281147},
	{ 0.980182135968117392690210009}, { 0.198098410717953586179324918},
	{-0.198098410717953586179324918}, { 0.980182135968117392690210009},
	{ 0.553016705580027531764226988}, { 0.833170164701913186439915922},
	{-0.833170164701913186439915922}, { 0.553016705580027531764226988},
	{ 0.829761233794523042469023765}, { 0.558118531220556115693702964},
	{-0.558118531220556115693702964}, { 0.829761233794523042469023765},
	{ 0.192080397049892441679288205}, { 0.981379193313754574318224190},
	{-0.981379193313754574318224190}, { 0.192080397049892441679288205},
	{ 0.994879330794805620591166107}, { 0.101069862754827824987887585},
	{-0.101069862754827824987887585}, { 0.994879330794805620591166107},
	{ 0.632018735939809021909403706}, { 0.774953106594873878359129282},
	{-0.774953106594873878359129282}, { 0.632018735939809021909403706},
	{ 0.880470889052160770806542929}, { 0.474100214650550014398580015},
	{-0.474100214650550014398580015}, { 0.880470889052160770806542929},
	{ 0.287347459544729526477331841}, { 0.957826413027532890321037029},
	{-0.957826413027532890321037029}, { 0.287347459544729526477331841},
	{ 0.956045251349996443270479823}, { 0.293219162694258650606608599},
	{-0.293219162694258650606608599}, { 0.956045251349996443270479823},
	{ 0.468688822035827933697617870}, { 0.883363338665731594736308015},
	{-0.883363338665731594736308015}, { 0.468688822035827933697617870},
	{ 0.771060524261813773200605759}, { 0.636761861236284230413943435},
	{-0.636761861236284230413943435}, { 0.771060524261813773200605759},
	{ 0.094963495329638998938034312}, { 0.995480755491926941769171600},
	{-0.995480755491926941769171600}, { 0.094963495329638998938034312},
	{ 0.998640218180265222418199049}, { 0.052131704680283321236358216},
	{-0.052131704680283321236358216}, { 0.998640218180265222418199049},
	{ 0.669282588346636065720696366}, { 0.743007952135121693517362293},
	{-0.743007952135121693517362293}, { 0.669282588346636065720696366},
	{ 0.902673318237258806751502391}, { 0.430326481340082633908199031},
	{-0.430326481340082633908199031}, { 0.902673318237258806751502391},
	{ 0.333999651442009404650865481}, { 0.942573197601446879280758735},
	{-0.942573197601446879280758735}, { 0.333999651442009404650865481},
	{ 0.969281235356548486048290738}, { 0.245955050335794611599924709},
	{-0.245955050335794611599924709}, { 0.969281235356548486048290738},
	{ 0.511468850437970399504391001}, { 0.859301818357008404783582139},
	{-0.859301818357008404783582139}, { 0.511468850437970399504391001},
	{ 0.801376171723140219430247777}, { 0.598160706996342311724958652},
	{-0.598160706996342311724958652}, { 0.801376171723140219430247777},
	{ 0.143695033150294454819773349}, { 0.989622017463200834623694454},
	{-0.989622017463200834623694454}, { 0.143695033150294454819773349},
	{ 0.988721691960323767604516485}, { 0.149764534677321517229695737},
	{-0.149764534677321517229695737}, { 0.988721691960323767604516485},
	{ 0.593232295039799808047809426}, { 0.805031331142963597922659282},
	{-0.805031331142963597922659282}, { 0.593232295039799808047809426},
	{ 0.856147328375194481019630732}, { 0.516731799017649881508753876},
	{-0.516731799017649881508753876}, { 0.856147328375194481019630732},
	{ 0.240003022448741486568922365}, { 0.970772140728950302138169611},
	{-0.970772140728950302138169611}, { 0.240003022448741486568922365},
	{ 0.940506070593268323787291309}, { 0.339776884406826857828825803},
	{-0.339776884406826857828825803}, { 0.940506070593268323787291309},
	{ 0.424779681209108833357226189}, { 0.905296759318118774354048329},
	{-0.905296759318118774354048329}, { 0.424779681209108833357226189},
	{ 0.738887324460615147933116508}, { 0.673829000378756060917568372},
	{-0.673829000378756060917568372}, { 0.738887324460615147933116508},
	{ 0.046003182130914628814301788}, { 0.998941293186856850633930266},
	{-0.998941293186856850633930266}, { 0.046003182130914628814301788},
	{ 0.999618822495178597116830637}, { 0.027608145778965741612354872},
	{-0.027608145778965741612354872}, { 0.999618822495178597116830637},
	{ 0.687315340891759108199186948}, { 0.726359155084345976817494315},
	{-0.726359155084345976817494315}, { 0.687315340891759108199186948},
	{ 0.912962190428398164628018233}, { 0.408044162864978680820747499},
	{-0.408044162864978680820747499}, { 0.912962190428398164628018233},
	{ 0.357030961233430032614954036}, { 0.934092550404258914729877883},
	{-0.934092550404258914729877883}, { 0.357030961233430032614954036},
	{ 0.975025345066994146844913468}, { 0.222093620973203534094094721},
	{-0.222093620973203534094094721}, { 0.975025345066994146844913468},
	{ 0.532403127877197971442805218}, { 0.846490938774052078300544488},
	{-0.846490938774052078300544488}, { 0.532403127877197971442805218},
	{ 0.815814410806733789010772660}, { 0.578313796411655563342245019},
	{-0.578313796411655563342245019}, { 0.815814410806733789010772660},
	{ 0.167938294974731178054745536}, { 0.985797509167567424700995000},
	{-0.985797509167567424700995000}, { 0.167938294974731178054745536},
	{ 0.992099313142191757112085445}, { 0.125454983411546238542336453},
	{-0.125454983411546238542336453}, { 0.992099313142191757112085445},
	{ 0.612810082429409703935211936}, { 0.790230221437310055030217152},
	{-0.790230221437310055030217152}, { 0.612810082429409703935211936},
	{ 0.868570705971340895340449876}, { 0.495565261825772531150266670},
	{-0.495565261825772531150266670}, { 0.868570705971340895340449876},
	{ 0.263754678974831383611349322}, { 0.964589793289812723836432159},
	{-0.964589793289812723836432159}, { 0.263754678974831383611349322},
	{ 0.948561349915730288158494826}, { 0.316593375556165867243047035},
	{-0.316593375556165867243047035}, { 0.948561349915730288158494826},
	{ 0.446868840162374195353044389}, { 0.894599485631382678433072126},
	{-0.894599485631382678433072126}, { 0.446868840162374195353044389},
	{ 0.755201376896536527598710756}, { 0.655492852999615385312679701},
	{-0.655492852999615385312679701}, { 0.755201376896536527598710756},
	{ 0.070504573389613863027351471}, { 0.997511456140303459699448390},
	{-0.997511456140303459699448390}, { 0.070504573389613863027351471},
	{ 0.997060070339482978987989949}, { 0.076623861392031492278332463},
	{-0.076623861392031492278332463}, { 0.997060070339482978987989949},
	{ 0.650846684996380915068975573}, { 0.759209188978388033485525443},
	{-0.759209188978388033485525443}, { 0.650846684996380915068975573},
	{ 0.891840709392342727796478697}, { 0.452349587233770874133026703},
	{-0.452349587233770874133026703}, { 0.891840709392342727796478697},
	{ 0.310767152749611495835997250}, { 0.950486073949481721759926101},
	{-0.950486073949481721759926101}, { 0.310767152749611495835997250},
	{ 0.962953266873683886347921481}, { 0.269668325572915106525464462},
	{-0.269668325572915106525464462}, { 0.962953266873683886347921481},
	{ 0.490226483288291154229598449}, { 0.871595086655951034842481435},
	{-0.871595086655951034842481435}, { 0.490226483288291154229598449},
	{ 0.786455213599085757522319464}, { 0.617647307937803932403979402},
	{-0.617647307937803932403979402}, { 0.786455213599085757522319464},
	{ 0.119365214810991364593637790}, { 0.992850414459865090793563344},
	{-0.992850414459865090793563344}, { 0.119365214810991364593637790},
	{ 0.984748501801904218556553176}, { 0.173983873387463827950700807},
	{-0.173983873387463827950700807}, { 0.984748501801904218556553176},
	{ 0.573297166698042212820171239}, { 0.819347520076796960824689637},
	{-0.819347520076796960824689637}, { 0.573297166698042212820171239},
	{ 0.843208239641845437161743865}, { 0.537587076295645482502214932},
	{-0.537587076295645482502214932}, { 0.843208239641845437161743865},
	{ 0.216106797076219509948385131}, { 0.976369731330021149312732194},
	{-0.976369731330021149312732194}, { 0.216106797076219509948385131},
	{ 0.931884265581668106718557199}, { 0.362755724367397216204854462},
	{-0.362755724367397216204854462}, { 0.931884265581668106718557199},
	{ 0.402434650859418441082533934}, { 0.915448716088267819566431292},
	{-0.915448716088267819566431292}, { 0.402434650859418441082533934},
	{ 0.722128193929215321243607198}, { 0.691759258364157774906734132},
	{-0.691759258364157774906734132}, { 0.722128193929215321243607198},
	{ 0.021474080275469507418374898}, { 0.999769405351215321657617036},
	{-0.999769405351215321657617036}, { 0.021474080275469507418374898},
	{ 0.999882347454212525633049627}, { 0.015339206284988101044151868},
	{-0.015339206284988101044151868}, { 0.999882347454212525633049627},
	{ 0.696177131491462944788582591}, { 0.717870045055731736211325329},
	{-0.717870045055731736211325329}, { 0.696177131491462944788582591},
	{ 0.917900775621390457642276297}, { 0.396809987416710328595290911},
	{-0.396809987416710328595290911}, { 0.917900775621390457642276297},
	{ 0.368466829953372331712746222}, { 0.929640895843181265457918066},
	{-0.929640895843181265457918066}, { 0.368466829953372331712746222},
	{ 0.977677357824509979943404762}, { 0.210111836880469621717489972},
	{-0.210111836880469621717489972}, { 0.977677357824509979943404762},
	{ 0.542750784864515906586768661}, { 0.839893794195999504583383987},
	{-0.839893794195999504583383987}, { 0.542750784864515906586768661},
	{ 0.822849781375826332046780034}, { 0.568258952670131549790548489},
	{-0.568258952670131549790548489}, { 0.822849781375826332046780034},
	{ 0.180022901405699522679906590}, { 0.983662419211730274396237776},
	{-0.983662419211730274396237776}, { 0.180022901405699522679906590},
	{ 0.993564135520595333782021697}, { 0.113270952177564349018228733},
	{-0.113270952177564349018228733}, { 0.993564135520595333782021697},
	{ 0.622461279374149972519166721}, { 0.782650596166575738458949301},
	{-0.782650596166575738458949301}, { 0.622461279374149972519166721},
	{ 0.874586652278176112634431897}, { 0.484869248000791101822951699},
	{-0.484869248000791101822951699}, { 0.874586652278176112634431897},
	{ 0.275571819310958163076425168}, { 0.961280485811320641748659653},
	{-0.961280485811320641748659653}, { 0.275571819310958163076425168},
	{ 0.952375012719765858529893608}, { 0.304929229735402406490728633},
	{-0.304929229735402406490728633}, { 0.952375012719765858529893608},
	{ 0.457813303598877221904961155}, { 0.889048355854664562540777729},
	{-0.889048355854664562540777729}, { 0.457813303598877221904961155},
	{ 0.763188417263381271704838297}, { 0.646176012983316364832802220},
	{-0.646176012983316364832802220}, { 0.763188417263381271704838297},
	{ 0.082740264549375693111987083}, { 0.996571145790554847093566910},
	{-0.996571145790554847093566910}, { 0.082740264549375693111987083},
	{ 0.997925286198596012623025462}, { 0.064382630929857460819324537},
	{-0.064382630929857460819324537}, { 0.997925286198596012623025462},
	{ 0.660114342067420478559490747}, { 0.751165131909686411205819422},
	{-0.751165131909686411205819422}, { 0.660114342067420478559490747},
	{ 0.897324580705418281231391836}, { 0.441371268731716692879988968},
	{-0.441371268731716692879988968}, { 0.897324580705418281231391836},
	{ 0.322407678801069848384807478}, { 0.946600913083283570044599823},
	{-0.946600913083283570044599823}, { 0.322407678801069848384807478},
	{ 0.966190003445412555433832961}, { 0.257831102162159005614471295},
	{-0.257831102162159005614471295}, { 0.966190003445412555433832961},
	{ 0.500885382611240786241285004}, { 0.865513624090569082825488358},
	{-0.865513624090569082825488358}, { 0.500885382611240786241285004},
	{ 0.793975477554337164895083757}, { 0.607949784967773667243642671},
	{-0.607949784967773667243642671}, { 0.793975477554337164895083757},
	{ 0.131540028702883111103387493}, { 0.991310859846115418957349799},
	{-0.991310859846115418957349799}, { 0.131540028702883111103387493},
	{ 0.986809401814185476970235952}, { 0.161886393780111837641387995},
	{-0.161886393780111837641387995}, { 0.986809401814185476970235952},
	{ 0.583308652937698294392830961}, { 0.812250586585203913049744181},
	{-0.812250586585203913049744181}, { 0.583308652937698294392830961},
	{ 0.849741768000852489471268395}, { 0.527199134781901348464274575},
	{-0.527199134781901348464274575}, { 0.849741768000852489471268395},
	{ 0.228072083170885739254457379}, { 0.973644249650811925318383912},
	{-0.973644249650811925318383912}, { 0.228072083170885739254457379},
	{ 0.936265667170278246576310996}, { 0.351292756085567125601307623},
	{-0.351292756085567125601307623}, { 0.936265667170278246576310996},
	{ 0.413638312238434547471944324}, { 0.910441292258067196934095369},
	{-0.910441292258067196934095369}, { 0.413638312238434547471944324},
	{ 0.730562769227827561177758850}, { 0.682845546385248068164596123},
	{-0.682845546385248068164596123}, { 0.730562769227827561177758850},
	{ 0.033741171851377584833716112}, { 0.999430604555461772019008327},
	{-0.999430604555461772019008327}, { 0.033741171851377584833716112},
	{ 0.999204758618363895492950001}, { 0.039872927587739811128578738},
	{-0.039872927587739811128578738}, { 0.999204758618363895492950001},
	{ 0.678350043129861486873655042}, { 0.734738878095963464563223604},
	{-0.734738878095963464563223604}, { 0.678350043129861486873655042},
	{ 0.907886116487666212038681480}, { 0.419216888363223956433010020},
	{-0.419216888363223956433010020}, { 0.907886116487666212038681480},
	{ 0.345541324963989065539191723}, { 0.938403534063108112192420774},
	{-0.938403534063108112192420774}, { 0.345541324963989065539191723},
	{ 0.972226497078936305708321144}, { 0.234041958583543423191242045},
	{-0.234041958583543423191242045}, { 0.972226497078936305708321144},
	{ 0.521975292937154342694258318}, { 0.852960604930363657746588082},
	{-0.852960604930363657746588082}, { 0.521975292937154342694258318},
	{ 0.808656181588174991946968128}, { 0.588281548222645304786439813},
	{-0.588281548222645304786439813}, { 0.808656181588174991946968128},
	{ 0.155828397654265235743101486}, { 0.987784141644572154230969032},
	{-0.987784141644572154230969032}, { 0.155828397654265235743101486},
	{ 0.990485084256457037998682243}, { 0.137620121586486044948441663},
	{-0.137620121586486044948441663}, { 0.990485084256457037998682243},
	{ 0.603066598540348201693430617}, { 0.797690840943391108362662755},
	{-0.797690840943391108362662755}, { 0.603066598540348201693430617},
	{ 0.862423956111040538690933878}, { 0.506186645345155291048942344},
	{-0.506186645345155291048942344}, { 0.862423956111040538690933878},
	{ 0.251897818154216950498106628}, { 0.967753837093475465243391912},
	{-0.967753837093475465243391912}, { 0.251897818154216950498106628},
	{ 0.944604837261480265659265493}, { 0.328209843579092526107916817},
	{-0.328209843579092526107916817}, { 0.944604837261480265659265493},
	{ 0.435857079922255491032544080}, { 0.900015892016160228714535267},
	{-0.900015892016160228714535267}, { 0.435857079922255491032544080},
	{ 0.747100605980180144323078847}, { 0.664710978203344868130324985},
	{-0.664710978203344868130324985}, { 0.747100605980180144323078847},
	{ 0.058258264500435759613979782}, { 0.998301544933892840738782163},
	{-0.998301544933892840738782163}, { 0.058258264500435759613979782},
	{ 0.996044700901251989887944810}, { 0.088853552582524596561586535},
	{-0.088853552582524596561586535}, { 0.996044700901251989887944810},
	{ 0.641481012808583151988739898}, { 0.767138911935820381181694573},
	{-0.767138911935820381181694573}, { 0.641481012808583151988739898},
	{ 0.886222530148880631647990821}, { 0.463259783551860197390719637},
	{-0.463259783551860197390719637}, { 0.886222530148880631647990821},
	{ 0.299079826308040476750336973}, { 0.954228095109105629780430732},
	{-0.954228095109105629780430732}, { 0.299079826308040476750336973},
	{ 0.959571513081984528335528181}, { 0.281464937925757984095231007},
	{-0.281464937925757984095231007}, { 0.959571513081984528335528181},
	{ 0.479493757660153026679839798}, { 0.877545290207261291668470750},
	{-0.877545290207261291668470750}, { 0.479493757660153026679839798},
	{ 0.778816512381475953374724325}, { 0.627251815495144113509622565},
	{-0.627251815495144113509622565}, { 0.778816512381475953374724325},
	{ 0.107172424956808849175529148}, { 0.994240449453187946358413442},
	{-0.994240449453187946358413442}, { 0.107172424956808849175529148},
	{ 0.982539302287441255907040396}, { 0.186055151663446648105438304},
	{-0.186055151663446648105438304}, { 0.982539302287441255907040396},
	{ 0.563199344013834115007363772}, { 0.826321062845663480311195452},
	{-0.826321062845663480311195452}, { 0.563199344013834115007363772},
	{ 0.836547727223511984524285790}, { 0.547894059173100165608820571},
	{-0.547894059173100165608820571}, { 0.836547727223511984524285790},
	{ 0.204108966092816874181696950}, { 0.978948175319062194715480124},
	{-0.978948175319062194715480124}, { 0.204108966092816874181696950},
	{ 0.927362525650401087274536959}, { 0.374164062971457997104393020},
	{-0.374164062971457997104393020}, { 0.927362525650401087274536959},
	{ 0.391170384302253888687512949}, { 0.920318276709110566440076541},
	{-0.920318276709110566440076541}, { 0.391170384302253888687512949},
	{ 0.713584868780793592903125099}, { 0.700568793943248366792866380},
	{-0.700568793943248366792866380}, { 0.713584868780793592903125099},
	{ 0.009203754782059819315102378}, { 0.999957644551963866333120920},
	{-0.999957644551963866333120920}, { 0.009203754782059819315102378},
	{ 0.999957644551963866333120920}, { 0.009203754782059819315102378},
	{-0.009203754782059819315102378}, { 0.999957644551963866333120920},
	{ 0.700568793943248366792866380}, { 0.713584868780793592903125099},
	{-0.713584868780793592903125099}, { 0.700568793943248366792866380},
	{ 0.920318276709110566440076541}, { 0.391170384302253888687512949},
	{-0.391170384302253888687512949}, { 0.920318276709110566440076541},
	{ 0.374164062971457997104393020}, { 0.927362525650401087274536959},
	{-0.927362525650401087274536959}, { 0.374164062971457997104393020},
	{ 0.978948175319062194715480124}, { 0.204108966092816874181696950},
	{-0.204108966092816874181696950}, { 0.978948175319062194715480124},
	{ 0.547894059173100165608820571}, { 0.836547727223511984524285790},
	{-0.836547727223511984524285790}, { 0.547894059173100165608820571},
	{ 0.826321062845663480311195452}, { 0.563199344013834115007363772},
	{-0.563199344013834115007363772}, { 0.826321062845663480311195452},
	{ 0.186055151663446648105438304}, { 0.982539302287441255907040396},
	{-0.982539302287441255907040396}, { 0.186055151663446648105438304},
	{ 0.994240449453187946358413442}, { 0.107172424956808849175529148},
	{-0.107172424956808849175529148}, { 0.994240449453187946358413442},
	{ 0.627251815495144113509622565}, { 0.778816512381475953374724325},
	{-0.778816512381475953374724325}, { 0.627251815495144113509622565},
	{ 0.877545290207261291668470750}, { 0.479493757660153026679839798},
	{-0.479493757660153026679839798}, { 0.877545290207261291668470750},
	{ 0.281464937925757984095231007}, { 0.959571513081984528335528181},
	{-0.959571513081984528335528181}, { 0.281464937925757984095231007},
	{ 0.954228095109105629780430732}, { 0.299079826308040476750336973},
	{-0.299079826308040476750336973}, { 0.954228095109105629780430732},
	{ 0.463259783551860197390719637}, { 0.886222530148880631647990821},
	{-0.886222530148880631647990821}, { 0.463259783551860197390719637},
	{ 0.767138911935820381181694573}, { 0.641481012808583151988739898},
	{-0.641481012808583151988739898}, { 0.767138911935820381181694573},
	{ 0.088853552582524596561586535}, { 0.996044700901251989887944810},
	{-0.996044700901251989887944810}, { 0.088853552582524596561586535},
	{ 0.998301544933892840738782163}, { 0.058258264500435759613979782},
	{-0.058258264500435759613979782}, { 0.998301544933892840738782163},
	{ 0.664710978203344868130324985}, { 0.747100605980180144323078847},
	{-0.747100605980180144323078847}, { 0.664710978203344868130324985},
	{ 0.900015892016160228714535267}, { 0.435857079922255491032544080},
	{-0.435857079922255491032544080}, { 0.900015892016160228714535267},
	{ 0.328209843579092526107916817}, { 0.944604837261480265659265493},
	{-0.944604837261480265659265493}, { 0.328209843579092526107916817},
	{ 0.967753837093475465243391912}, { 0.251897818154216950498106628},
	{-0.251897818154216950498106628}, { 0.967753837093475465243391912},
	{ 0.506186645345155291048942344}, { 0.862423956111040538690933878},
	{-0.862423956111040538690933878}, { 0.506186645345155291048942344},
	{ 0.797690840943391108362662755}, { 0.603066598540348201693430617},
	{-0.603066598540348201693430617}, { 0.797690840943391108362662755},
	{ 0.137620121586486044948441663}, { 0.990485084256457037998682243},
	{-0.990485084256457037998682243}, { 0.137620121586486044948441663},
	{ 0.987784141644572154230969032}, { 0.155828397654265235743101486},
	{-0.155828397654265235743101486}, { 0.987784141644572154230969032},
	{ 0.588281548222645304786439813}, { 0.808656181588174991946968128},
	{-0.808656181588174991946968128}, { 0.588281548222645304786439813},
	{ 0.852960604930363657746588082}, { 0.521975292937154342694258318},
	{-0.521975292937154342694258318}, { 0.852960604930363657746588082},
	{ 0.234041958583543423191242045}, { 0.972226497078936305708321144},
	{-0.972226497078936305708321144}, { 0.234041958583543423191242045},
	{ 0.938403534063108112192420774}, { 0.345541324963989065539191723},
	{-0.345541324963989065539191723}, { 0.938403534063108112192420774},
	{ 0.419216888363223956433010020}, { 0.907886116487666212038681480},
	{-0.907886116487666212038681480}, { 0.419216888363223956433010020},
	{ 0.734738878095963464563223604}, { 0.678350043129861486873655042},
	{-0.678350043129861486873655042}, { 0.734738878095963464563223604},
	{ 0.039872927587739811128578738}, { 0.999204758618363895492950001},
	{-0.999204758618363895492950001}, { 0.039872927587739811128578738},
	{ 0.999430604555461772019008327}, { 0.033741171851377584833716112},
	{-0.033741171851377584833716112}, { 0.999430604555461772019008327},
	{ 0.682845546385248068164596123}, { 0.730562769227827561177758850},
	{-0.730562769227827561177758850}, { 0.682845546385248068164596123},
	{ 0.910441292258067196934095369}, { 0.413638312238434547471944324},
	{-0.413638312238434547471944324}, { 0.910441292258067196934095369},
	{ 0.351292756085567125601307623}, { 0.936265667170278246576310996},
	{-0.936265667170278246576310996}, { 0.351292756085567125601307623},
	{ 0.973644249650811925318383912}, { 0.228072083170885739254457379},
	{-0.228072083170885739254457379}, { 0.973644249650811925318383912},
	{ 0.527199134781901348464274575}, { 0.849741768000852489471268395},
	{-0.849741768000852489471268395}, { 0.527199134781901348464274575},
	{ 0.812250586585203913049744181}, { 0.583308652937698294392830961},
	{-0.583308652937698294392830961}, { 0.812250586585203913049744181},
	{ 0.161886393780111837641387995}, { 0.986809401814185476970235952},
	{-0.986809401814185476970235952}, { 0.161886393780111837641387995},
	{ 0.991310859846115418957349799}, { 0.131540028702883111103387493},
	{-0.131540028702883111103387493}, { 0.991310859846115418957349799},
	{ 0.607949784967773667243642671}, { 0.793975477554337164895083757},
	{-0.793975477554337164895083757}, { 0.607949784967773667243642671},
	{ 0.865513624090569082825488358}, { 0.500885382611240786241285004},
	{-0.500885382611240786241285004}, { 0.865513624090569082825488358},
	{ 0.257831102162159005614471295}, { 0.966190003445412555433832961},
	{-0.966190003445412555433832961}, { 0.257831102162159005614471295},
	{ 0.946600913083283570044599823}, { 0.322407678801069848384807478},
	{-0.322407678801069848384807478}, { 0.946600913083283570044599823},
	{ 0.441371268731716692879988968}, { 0.897324580705418281231391836},
	{-0.897324580705418281231391836}, { 0.441371268731716692879988968},
	{ 0.751165131909686411205819422}, { 0.660114342067420478559490747},
	{-0.660114342067420478559490747}, { 0.751165131909686411205819422},
	{ 0.064382630929857460819324537}, { 0.997925286198596012623025462},
	{-0.997925286198596012623025462}, { 0.064382630929857460819324537},
	{ 0.996571145790554847093566910}, { 0.082740264549375693111987083},
	{-0.082740264549375693111987083}, { 0.996571145790554847093566910},
	{ 0.646176012983316364832802220}, { 0.763188417263381271704838297},
	{-0.763188417263381271704838297}, { 0.646176012983316364832802220},
	{ 0.889048355854664562540777729}, { 0.457813303598877221904961155},
	{-0.457813303598877221904961155}, { 0.889048355854664562540777729},
	{ 0.304929229735402406490728633}, { 0.952375012719765858529893608},
	{-0.952375012719765858529893608}, { 0.304929229735402406490728633},
	{ 0.961280485811320641748659653}, { 0.275571819310958163076425168},
	{-0.275571819310958163076425168}, { 0.961280485811320641748659653},
	{ 0.484869248000791101822951699}, { 0.874586652278176112634431897},
	{-0.874586652278176112634431897}, { 0.484869248000791101822951699},
	{ 0.782650596166575738458949301}, { 0.622461279374149972519166721},
	{-0.622461279374149972519166721}, { 0.782650596166575738458949301},
	{ 0.113270952177564349018228733}, { 0.993564135520595333782021697},
	{-0.993564135520595333782021697}, { 0.113270952177564349018228733},
	{ 0.983662419211730274396237776}, { 0.180022901405699522679906590},
	{-0.180022901405699522679906590}, { 0.983662419211730274396237776},
	{ 0.568258952670131549790548489}, { 0.822849781375826332046780034},
	{-0.822849781375826332046780034}, { 0.568258952670131549790548489},
	{ 0.839893794195999504583383987}, { 0.542750784864515906586768661},
	{-0.542750784864515906586768661}, { 0.839893794195999504583383987},
	{ 0.210111836880469621717489972}, { 0.977677357824509979943404762},
	{-0.977677357824509979943404762}, { 0.210111836880469621717489972},
	{ 0.929640895843181265457918066}, { 0.368466829953372331712746222},
	{-0.368466829953372331712746222}, { 0.929640895843181265457918066},
	{ 0.396809987416710328595290911}, { 0.917900775621390457642276297},
	{-0.917900775621390457642276297}, { 0.396809987416710328595290911},
	{ 0.717870045055731736211325329}, { 0.696177131491462944788582591},
	{-0.696177131491462944788582591}, { 0.717870045055731736211325329},
	{ 0.015339206284988101044151868}, { 0.999882347454212525633049627},
	{-0.999882347454212525633049627}, { 0.015339206284988101044151868},
	{ 0.999769405351215321657617036}, { 0.021474080275469507418374898},
	{-0.021474080275469507418374898}, { 0.999769405351215321657617036},
	{ 0.691759258364157774906734132}, { 0.722128193929215321243607198},
	{-0.722128193929215321243607198}, { 0.691759258364157774906734132},
	{ 0.915448716088267819566431292}, { 0.402434650859418441082533934},
	{-0.402434650859418441082533934}, { 0.915448716088267819566431292},
	{ 0.362755724367397216204854462}, { 0.931884265581668106718557199},
	{-0.931884265581668106718557199}, { 0.362755724367397216204854462},
	{ 0.976369731330021149312732194}, { 0.216106797076219509948385131},
	{-0.216106797076219509948385131}, { 0.976369731330021149312732194},
	{ 0.537587076295645482502214932}, { 0.843208239641845437161743865},
	{-0.843208239641845437161743865}, { 0.537587076295645482502214932},
	{ 0.819347520076796960824689637}, { 0.573297166698042212820171239},
	{-0.573297166698042212820171239}, { 0.819347520076796960824689637},
	{ 0.173983873387463827950700807}, { 0.984748501801904218556553176},
	{-0.984748501801904218556553176}, { 0.173983873387463827950700807},
	{ 0.992850414459865090793563344}, { 0.119365214810991364593637790},
	{-0.119365214810991364593637790}, { 0.992850414459865090793563344},
	{ 0.617647307937803932403979402}, { 0.786455213599085757522319464},
	{-0.786455213599085757522319464}, { 0.617647307937803932403979402},
	{ 0.871595086655951034842481435}, { 0.490226483288291154229598449},
	{-0.490226483288291154229598449}, { 0.871595086655951034842481435},
	{ 0.269668325572915106525464462}, { 0.962953266873683886347921481},
	{-0.962953266873683886347921481}, { 0.269668325572915106525464462},
	{ 0.950486073949481721759926101}, { 0.310767152749611495835997250},
	{-0.310767152749611495835997250}, { 0.950486073949481721759926101},
	{ 0.452349587233770874133026703}, { 0.891840709392342727796478697},
	{-0.891840709392342727796478697}, { 0.452349587233770874133026703},
	{ 0.759209188978388033485525443}, { 0.650846684996380915068975573},
	{-0.650846684996380915068975573}, { 0.759209188978388033485525443},
	{ 0.076623861392031492278332463}, { 0.997060070339482978987989949},
	{-0.997060070339482978987989949}, { 0.076623861392031492278332463},
	{ 0.997511456140303459699448390}, { 0.070504573389613863027351471},
	{-0.070504573389613863027351471}, { 0.997511456140303459699448390},
	{ 0.655492852999615385312679701}, { 0.755201376896536527598710756},
	{-0.755201376896536527598710756}, { 0.655492852999615385312679701},
	{ 0.894599485631382678433072126}, { 0.446868840162374195353044389},
	{-0.446868840162374195353044389}, { 0.894599485631382678433072126},
	{ 0.316593375556165867243047035}, { 0.948561349915730288158494826},
	{-0.948561349915730288158494826}, { 0.316593375556165867243047035},
	{ 0.964589793289812723836432159}, { 0.263754678974831383611349322},
	{-0.263754678974831383611349322}, { 0.964589793289812723836432159},
	{ 0.495565261825772531150266670}, { 0.868570705971340895340449876},
	{-0.868570705971340895340449876}, { 0.495565261825772531150266670},
	{ 0.790230221437310055030217152}, { 0.612810082429409703935211936},
	{-0.612810082429409703935211936}, { 0.790230221437310055030217152},
	{ 0.125454983411546238542336453}, { 0.992099313142191757112085445},
	{-0.992099313142191757112085445}, { 0.125454983411546238542336453},
	{ 0.985797509167567424700995000}, { 0.167938294974731178054745536},
	{-0.167938294974731178054745536}, { 0.985797509167567424700995000},
	{ 0.578313796411655563342245019}, { 0.815814410806733789010772660},
	{-0.815814410806733789010772660}, { 0.578313796411655563342245019},
	{ 0.846490938774052078300544488}, { 0.532403127877197971442805218},
	{-0.532403127877197971442805218}, { 0.846490938774052078300544488},
	{ 0.222093620973203534094094721}, { 0.975025345066994146844913468},
	{-0.975025345066994146844913468}, { 0.222093620973203534094094721},
	{ 0.934092550404258914729877883}, { 0.357030961233430032614954036},
	{-0.357030961233430032614954036}, { 0.934092550404258914729877883},
	{ 0.408044162864978680820747499}, { 0.912962190428398164628018233},
	{-0.912962190428398164628018233}, { 0.408044162864978680820747499},
	{ 0.726359155084345976817494315}, { 0.687315340891759108199186948},
	{-0.687315340891759108199186948}, { 0.726359155084345976817494315},
	{ 0.027608145778965741612354872}, { 0.999618822495178597116830637},
	{-0.999618822495178597116830637}, { 0.027608145778965741612354872},
	{ 0.998941293186856850633930266}, { 0.046003182130914628814301788},
	{-0.046003182130914628814301788}, { 0.998941293186856850633930266},
	{ 0.673829000378756060917568372}, { 0.738887324460615147933116508},
	{-0.738887324460615147933116508}, { 0.673829000378756060917568372},
	{ 0.905296759318118774354048329}, { 0.424779681209108833357226189},
	{-0.424779681209108833357226189}, { 0.905296759318118774354048329},
	{ 0.339776884406826857828825803}, { 0.940506070593268323787291309},
	{-0.940506070593268323787291309}, { 0.339776884406826857828825803},
	{ 0.970772140728950302138169611}, { 0.240003022448741486568922365},
	{-0.240003022448741486568922365}, { 0.970772140728950302138169611},
	{ 0.516731799017649881508753876}, { 0.856147328375194481019630732},
	{-0.856147328375194481019630732}, { 0.516731799017649881508753876},
	{ 0.805031331142963597922659282}, { 0.593232295039799808047809426},
	{-0.593232295039799808047809426}, { 0.805031331142963597922659282},
	{ 0.149764534677321517229695737}, { 0.988721691960323767604516485},
	{-0.988721691960323767604516485}, { 0.149764534677321517229695737},
	{ 0.989622017463200834623694454}, { 0.143695033150294454819773349},
	{-0.143695033150294454819773349}, { 0.989622017463200834623694454},
	{ 0.598160706996342311724958652}, { 0.801376171723140219430247777},
	{-0.801376171723140219430247777}, { 0.598160706996342311724958652},
	{ 0.859301818357008404783582139}, { 0.511468850437970399504391001},
	{-0.511468850437970399504391001}, { 0.859301818357008404783582139},
	{ 0.245955050335794611599924709}, { 0.969281235356548486048290738},
	{-0.969281235356548486048290738}, { 0.245955050335794611599924709},
	{ 0.942573197601446879280758735}, { 0.333999651442009404650865481},
	{-0.333999651442009404650865481}, { 0.942573197601446879280758735},
	{ 0.430326481340082633908199031}, { 0.902673318237258806751502391},
	{-0.902673318237258806751502391}, { 0.430326481340082633908199031},
	{ 0.743007952135121693517362293}, { 0.669282588346636065720696366},
	{-0.669282588346636065720696366}, { 0.743007952135121693517362293},
	{ 0.052131704680283321236358216}, { 0.998640218180265222418199049},
	{-0.998640218180265222418199049}, { 0.052131704680283321236358216},
	{ 0.995480755491926941769171600}, { 0.094963495329638998938034312},
	{-0.094963495329638998938034312}, { 0.995480755491926941769171600},
	{ 0.636761861236284230413943435}, { 0.771060524261813773200605759},
	{-0.771060524261813773200605759}, { 0.636761861236284230413943435},
	{ 0.883363338665731594736308015}, { 0.468688822035827933697617870},
	{-0.468688822035827933697617870}, { 0.883363338665731594736308015},
	{ 0.293219162694258650606608599}, { 0.956045251349996443270479823},
	{-0.956045251349996443270479823}, { 0.293219162694258650606608599},
	{ 0.957826413027532890321037029}, { 0.287347459544729526477331841},
	{-0.287347459544729526477331841}, { 0.957826413027532890321037029},
	{ 0.474100214650550014398580015}, { 0.880470889052160770806542929},
	{-0.880470889052160770806542929}, { 0.474100214650550014398580015},
	{ 0.774953106594873878359129282}, { 0.632018735939809021909403706},
	{-0.632018735939809021909403706}, { 0.774953106594873878359129282},
	{ 0.101069862754827824987887585}, { 0.994879330794805620591166107},
	{-0.994879330794805620591166107}, { 0.101069862754827824987887585},
	{ 0.981379193313754574318224190}, { 0.192080397049892441679288205},
	{-0.192080397049892441679288205}, { 0.981379193313754574318224190},
	{ 0.558118531220556115693702964}, { 0.829761233794523042469023765},
	{-0.829761233794523042469023765}, { 0.558118531220556115693702964},
	{ 0.833170164701913186439915922}, { 0.553016705580027531764226988},
	{-0.553016705580027531764226988}, { 0.833170164701913186439915922},
	{ 0.198098410717953586179324918}, { 0.980182135968117392690210009},
	{-0.980182135968117392690210009}, { 0.198098410717953586179324918},
	{ 0.925049240782677590302371869}, { 0.379847208924051170576281147},
	{-0.379847208924051170576281147}, { 0.925049240782677590302371869},
	{ 0.385516053843918864075607949}, { 0.922701128333878570437264227},
	{-0.922701128333878570437264227}, { 0.385516053843918864075607949},
	{ 0.709272826438865651316533772}, { 0.704934080375904908852523758},
	{-0.704934080375904908852523758}, { 0.709272826438865651316533772},
	{ 0.003067956762965976270145365}, { 0.999995293809576171511580126},
	{-0.999995293809576171511580126}, { 0.003067956762965976270145365}
};

const fpr fpr_p2_tab[] = {
	{ 2.00000000000 },
	{ 1.00000000000 },
	{ 0.50000000000 },
	{ 0.25000000000 },
	{ 0.12500000000 },
	{ 0.06250000000 },
	{ 0.03125000000 },
	{ 0.01562500000 },
	{ 0.00781250000 },
	{ 0.00390625000 },
	{ 0.00195312500 }
};

#elif FALCON_FXP // yyyFPNATIVE+0 yyyFPEMU+0 yyyFXP+1

/* 
 * This is efficient, but leaks a lot of information on b (not on a)
 * and would probably not work on masked implementation since we would 
 * loose to much information when normalizing b...
*/
uint128_t inner_fpr_div_not_ct(fpr a, fpr b) {
	// Shift FP to have the most significant bit set
	int shl = 0; // normalization shift
	/* NOT CT: this leaks the integer part of b */
	while ( b.v >= UINT128(1, 0) ) { b.v >>= 1; shl++; }	
	/* NOT CT: this leaks the number of 0 MSb in the deciamal part of b */
	while ( ((b.v & (1ull << 63)) == 0) ) { b.v <<= 1; shl--; }
	a.v >>= shl;
	
	fpr x = fpr_sub(c48_17 , fpr_mul(b, c32_17));
	int i;
	for (i=0;i<5;i++) // iterate 5 time is enough to get 64bits precision
		x = fpr_add(x , fpr_mul(x, fpr_sub(fpr_one, fpr_mul(b, x))));;
	a = fpr_mul(a, x);
	return a.v;
}

/* 
 * This is efficient, but would probably not work on masked implementation
 *  since we would loose to much information when normalizing b...
*/
uint128_t inner_fpr_div(fpr a, fpr b) {
	int64_t shift = 0; // normalization shift
	uint64_t tmp;
	
	// normalize b in constant time with intrinsic to count leading bits of 
	// both integer and fractional part
	tmp = fpr_trunc(b);
    cmov_int(&shift, (clz(tmp | 1U) ^ (sizeof(tmp)*8-1))+1, tmp != 0);
    b = fpr_div2e(b, shift);
	a = fpr_div2e(a, shift);

	shift = 0;
    tmp = fpr_frac(b);
    cmov_int(&shift, clz(tmp | 1U), tmp != 0);
    b = fpr_mul2e(b, shift);
	a = fpr_mul2e(a, shift);
	
	
	
	fpr x = fpr_sub(c48_17 , fpr_mul(b, c32_17));
	int i;
	for (i=0;i<5;i++) // iterate 5 time is enough to get 64-bit precision
	{
		x = fpr_add(x , fpr_mul(x, fpr_sub(fpr_one, fpr_mul(b, x))));;
	}
	a = fpr_mul(a, x);
	return a.v;
}

uint128_t
inner_fpr_div_ct(fpr x_fpr, fpr y_fpr)
{
	uint128_t x = x_fpr.v, y = y_fpr.v;
	/*
	 * Get absolute values and signs. From now on, we can suppose
	 * that x and y fit on 63 bits (we ignore edge conditions).
	 */
	uint128_t sx = x >> 127;
	x = (x ^ -sx) + sx;
	uint128_t sy = y >> 127;
	y = (y ^ -sy) + sy;

	/*
	 * Do a bit by bit division, assuming that the quotient fits.
	 * The numerator starts at x*2^63, and is shifted one bit a time.
	 */
	uint128_t q = 0;
	uint128_t num = x >> 63;
	for (int i = 127; i >= 65; i --) {
		uint128_t b = 1 - ((num - y) >> 127);
		q |= b << i;
		num -= y & -b;
		num <<= 1;
		num |= (x >> (i - 65)) & 1;
	}
	for (int i = 64; i >= 0; i --) {
		uint128_t b = 1 - ((num - y) >> 127);
		q |= b << i;
		num -= y & -b;
		num <<= 1;
	}

	/*
	 * Rounding: if the remainder is at least y/2 (scaled), we add
	 * 2^(-64) to the quotient.
	 */
	uint128_t b = 1 - ((num - y) >> 127);
	q += b;

	/*
	 * Sign management: if the original x and y had different signs,
	 * then we must negate the quotient.
	 */
	sx ^= sy;
	q = (q ^ -sx) + sx;

	return q;
}

/*
 * Compute an estimation of the square root using Newton approximation.
 * This method allows to iteration on the approximation until we reach the
 * desired precision. 
 * In Falcon, sqrt is used during hte tree computation, which is in turn used 
 * during the sampling of the signature. Therefore, we need a high level of 
 * precision to keep our distribution indistinguishable from the ideal one, an 
 * preserve security.
 * We do 10 iterations of the algorithm to reach 64 bits of accuracy, we 
 * usually reach the max accuracy within 4-5 iterations. 
*/
fpr fpr_sqrt_newton(fpr x) {
	// start our guess at 100, since input values are around 10 0000
	uint128_t ni;	
	fpr n_prev = {UINT128(100, 0)}; 
	// Converge quadratically
	for (unsigned i = 0; i < 10; i++) {
		ni = (n_prev.v + inner_fpr_div(x, n_prev)) >> 1;
		n_prev.v = ni;
	}
	x.v = ni; 
	return x;
}

/* Compute 1/sqrt(x) using Newton approximation. 
 * Given a good approximation, it converges quadratically, so only few 
 * iterations are needed. We only handle values around x = 10 000, so we can
 * give start with a very accurate approximation.
 * We compute 1/sqrt(x), then multiply the reuslt by x to avoid, using 
 * divisions during the iterative approximation.
 */
fpr fpr_sqrt(fpr x) {
	fpr ni;	
#if defined(FMLDSA_Q64) && FMLDSA_Q64
	/* RN(1/sqrt(8380417)) in Q64.64. */
	fpr n_prev = {UINT128(0, 0x0016a372eb82f3eb)};
#else
	fpr n_prev = {UINT128(0, 0x028f5c28f5c28f5d)}; /* 0.01 */
#endif
	fpr cst_3 = {UINT128(3, 0)}; 
	// Converge quadratically
	for (unsigned i = 0; i < 8; i++) {
		ni = fpr_mul(n_prev, 
				fpr_div2e(fpr_sub(cst_3, fpr_mul(x, fpr_sqr(n_prev))), 1));
		n_prev.v = ni.v;
	}
	x = fpr_mul(ni, x);
	return x;
}

/*
 * Compute an estimation of the square root with a dicothomy search.
 * In Falcon, sqrt is used during hte tree computation, which is in turn used 
 * during the sampling of the signature. Therefore, we need a high level of 
 * precision to keep our distribution indistinguishable from the ideal one, an 
 * preserve security.
 * WARNING: This implementation is not constant-time, but needs to be... 
*/
fpr fpr_sqrt_not_ct(fpr x) {
	// start our guess at 100, since input values are around 10 0000
	uint128_t ni = UINT128(100, 0), n_prev = ni>>1;	 

	// look for an interavll [2^k, 2^(k+1)) that fits
	while(((ni>>64)*(ni>>64)) < (x.v>>64)){
		n_prev = ni;
		ni <<= 1;
	}
	
	// Find the square root by dichotomy
	uint128_t diff = ni - n_prev;
	fpr tmp = {0};
	while(diff > (uint128_t)1) {
		tmp.v = (ni + n_prev) >> 1;	
		if (fpr_lt(x, fpr_sqr(tmp))) {
			ni = tmp.v;
		}
		else {
			n_prev = tmp.v;
		}
		diff = ni - n_prev;
	}

	x.v = tmp.v;
    return x;
}



uint64_t fpr_expm_p63(fpr x, fpr ccs)
{
	/*
	 * Polynomial approximation of exp(-x) is taken from FACCT:
	 *   https://eprint.iacr.org/2018/1234
	 * Specifically, values are extracted from the implementation
	 * referenced from the FACCT article, and available at:
	 *   https://github.com/raykzhao/gaussian
	 * Here, the coefficients have been scaled up by 2^63 and
	 * converted to integers.
	 *
	 * Tests over more than 24 billions of random inputs in the
	 * 0..log(2) range have never shown a deviation larger than
	 * 2^(-50) from the true mathematical value.
	 */
	static const uint64_t C[] = {
		0x00000004741183A3u,
		0x00000036548CFC06u,
		0x0000024FDCBF140Au,
		0x0000171D939DE045u,
		0x0000D00CF58F6F84u,
		0x000680681CF796E3u,
		0x002D82D8305B0FEAu,
		0x011111110E066FD0u,
		0x0555555555070F00u,
		0x155555555581FF00u,
		0x400000000002B400u,
		0x7FFFFFFFFFFF4800u,
		0x8000000000000000u
	};

	uint64_t z, y;
	unsigned u;

	y = C[0];
	z = fpr_scale_p63(x);
	print_debug(("\t\t\t\tz = %016lx\n", z));
	for (u = 1; u < (sizeof C) / sizeof(C[0]); u ++) {
		/*
		 * Compute product z * y over 128 bits, but keep only
		 * the top 64 bits.
		 */
		uint64_t c = (uint64_t) (((uint128_t)z * (uint128_t)y) >> 64);
		y = C[u] - c;
	}
	/* lower nibble doesn't match for some reasons... but it still works. The reference impelmentaiton prabably drops some value since it scale down the result by 2^64 anyway */
	print_debug(("\t\t\t\ty = %016lx\n", y));	

	/*
	 * The scaling factor must be applied at the end. Since y is now
	 * in fixed-point notation, we have to convert the factor to the
	 * same format, and do an extra integer multiplication.
	 * 
	 * Beware, when ccs >= 1, the scaling will overflow, and we loose 
	 * the MSB. To prevent errors, the scaling then returns 0, and it 
	 * is handled by the calling function.
	 * I am not a big fan of this, but th'ats how FPEMU works...
	 */
	z = fpr_scale_p63(ccs);
	print_debug(("\t\t\t\tz = %016lx\n", z));
	y = (uint64_t) (((uint128_t)z * (uint128_t)y) >> 64);

	return y;
}


const fpr fpr_gm_tab[] = {
	{UINT128(0x0000000000000000, 0x0000000000000000)}, {UINT128(0x0000000000000000, 0x0000000000000000)},
	{UINT128(0x0000000000000000, 0x0000000000000000)}, {UINT128(0x0000000000000001, 0x0000000000000000)},
	{UINT128(0x0000000000000000, 0xb504f333f9de6485)}, {UINT128(0x0000000000000000, 0xb504f333f9de6485)},
	{UINT128(0xffffffffffffffff, 0x4afb0ccc06219b7b)}, {UINT128(0x0000000000000000, 0xb504f333f9de6485)},
	{UINT128(0x0000000000000000, 0xec835e79946a3146)}, {UINT128(0x0000000000000000, 0x61f78a9abaa58b47)},
	{UINT128(0xffffffffffffffff, 0x9e087565455a74b9)}, {UINT128(0x0000000000000000, 0xec835e79946a3146)},
	{UINT128(0x0000000000000000, 0x61f78a9abaa58b47)}, {UINT128(0x0000000000000000, 0xec835e79946a3146)},
	{UINT128(0xffffffffffffffff, 0x137ca1866b95ceba)}, {UINT128(0x0000000000000000, 0x61f78a9abaa58b47)},
	{UINT128(0x0000000000000000, 0xfb14be7fbae58157)}, {UINT128(0x0000000000000000, 0x31f17078d34c156d)},
	{UINT128(0xffffffffffffffff, 0xce0e8f872cb3ea93)}, {UINT128(0x0000000000000000, 0xfb14be7fbae58157)},
	{UINT128(0x0000000000000000, 0x8e39d9cd73464365)}, {UINT128(0x0000000000000000, 0xd4db3148750d181a)},
	{UINT128(0xffffffffffffffff, 0x2b24ceb78af2e7e6)}, {UINT128(0x0000000000000000, 0x8e39d9cd73464365)},
	{UINT128(0x0000000000000000, 0xd4db3148750d181a)}, {UINT128(0x0000000000000000, 0x8e39d9cd73464365)},
	{UINT128(0xffffffffffffffff, 0x71c626328cb9bc9b)}, {UINT128(0x0000000000000000, 0xd4db3148750d181a)},
	{UINT128(0x0000000000000000, 0x31f17078d34c156d)}, {UINT128(0x0000000000000000, 0xfb14be7fbae58157)},
	{UINT128(0xffffffffffffffff, 0x04eb4180451a7ea9)}, {UINT128(0x0000000000000000, 0x31f17078d34c156d)},
	{UINT128(0x0000000000000000, 0xfec46d1e89292cf1)}, {UINT128(0x0000000000000000, 0x1917a6bc29b42be2)},
	{UINT128(0xffffffffffffffff, 0xe6e85943d64bd41e)}, {UINT128(0x0000000000000000, 0xfec46d1e89292cf1)},
	{UINT128(0x0000000000000000, 0xa267992848eeb0c1)}, {UINT128(0x0000000000000000, 0xc5e40358a8ba05a8)},
	{UINT128(0xffffffffffffffff, 0x3a1bfca75745fa58)}, {UINT128(0x0000000000000000, 0xa267992848eeb0c1)},
	{UINT128(0x0000000000000000, 0xe1c5978c05ed8692)}, {UINT128(0x0000000000000000, 0x78ad74e01bd8ec79)},
	{UINT128(0xffffffffffffffff, 0x87528b1fe4271387)}, {UINT128(0x0000000000000000, 0xe1c5978c05ed8692)},
	{UINT128(0x0000000000000000, 0x4a5018bb567c16a3)}, {UINT128(0x0000000000000000, 0xf4fa0ab6316ed2ed)},
	{UINT128(0xffffffffffffffff, 0x0b05f549ce912d13)}, {UINT128(0x0000000000000000, 0x4a5018bb567c16a3)},
	{UINT128(0x0000000000000000, 0xf4fa0ab6316ed2ed)}, {UINT128(0x0000000000000000, 0x4a5018bb567c16a3)},
	{UINT128(0xffffffffffffffff, 0xb5afe744a983e95d)}, {UINT128(0x0000000000000000, 0xf4fa0ab6316ed2ed)},
	{UINT128(0x0000000000000000, 0x78ad74e01bd8ec79)}, {UINT128(0x0000000000000000, 0xe1c5978c05ed8692)},
	{UINT128(0xffffffffffffffff, 0x1e3a6873fa12796e)}, {UINT128(0x0000000000000000, 0x78ad74e01bd8ec79)},
	{UINT128(0x0000000000000000, 0xc5e40358a8ba05a8)}, {UINT128(0x0000000000000000, 0xa267992848eeb0c1)},
	{UINT128(0xffffffffffffffff, 0x5d9866d7b7114f3f)}, {UINT128(0x0000000000000000, 0xc5e40358a8ba05a8)},
	{UINT128(0x0000000000000000, 0x1917a6bc29b42be2)}, {UINT128(0x0000000000000000, 0xfec46d1e89292cf1)},
	{UINT128(0xffffffffffffffff, 0x013b92e176d6d30f)}, {UINT128(0x0000000000000000, 0x1917a6bc29b42be2)},
	{UINT128(0x0000000000000000, 0xffb10f1bcb6bef1e)}, {UINT128(0x0000000000000000, 0x0c8fb2f886ec09f4)},
	{UINT128(0xffffffffffffffff, 0xf3704d077913f60c)}, {UINT128(0x0000000000000000, 0xffb10f1bcb6bef1e)},
	{UINT128(0x0000000000000000, 0xabeb49a46764fd16)}, {UINT128(0x0000000000000000, 0xbdaef913557d76f1)},
	{UINT128(0xffffffffffffffff, 0x425106ecaa82890f)}, {UINT128(0x0000000000000000, 0xabeb49a46764fd16)},
	{UINT128(0x0000000000000000, 0xe76bd7a1e63b9787)}, {UINT128(0x0000000000000000, 0x6d744027857300ae)},
	{UINT128(0xffffffffffffffff, 0x928bbfd87a8cff52)}, {UINT128(0x0000000000000000, 0xe76bd7a1e63b9787)},
	{UINT128(0x0000000000000000, 0x563e69d6ac7f73f9)}, {UINT128(0x0000000000000000, 0xf1090827b43725fe)},
	{UINT128(0xffffffffffffffff, 0x0ef6f7d84bc8da02)}, {UINT128(0x0000000000000000, 0x563e69d6ac7f73f9)},
	{UINT128(0x0000000000000000, 0xf853f7dc9186b953)}, {UINT128(0x0000000000000000, 0x3e33f2f642be355f)},
	{UINT128(0xffffffffffffffff, 0xc1cc0d09bd41caa1)}, {UINT128(0x0000000000000000, 0xf853f7dc9186b953)},
	{UINT128(0x0000000000000000, 0x839c3cc917ff6cb5)}, {UINT128(0x0000000000000000, 0xdb941a28cb71ec88)},
	{UINT128(0xffffffffffffffff, 0x246be5d7348e1378)}, {UINT128(0x0000000000000000, 0x839c3cc917ff6cb5)},
	{UINT128(0x0000000000000000, 0xcd9f023f9c3a059f)}, {UINT128(0x0000000000000000, 0x987fbfe70b81a709)},
	{UINT128(0xffffffffffffffff, 0x67804018f47e58f7)}, {UINT128(0x0000000000000000, 0xcd9f023f9c3a059f)},
	{UINT128(0x0000000000000000, 0x259020dd1cc27445)}, {UINT128(0x0000000000000000, 0xfd3aabf84528b50c)},
	{UINT128(0xffffffffffffffff, 0x02c55407bad74af4)}, {UINT128(0x0000000000000000, 0x259020dd1cc27445)},
	{UINT128(0x0000000000000000, 0xfd3aabf84528b50c)}, {UINT128(0x0000000000000000, 0x259020dd1cc27445)},
	{UINT128(0xffffffffffffffff, 0xda6fdf22e33d8bbb)}, {UINT128(0x0000000000000000, 0xfd3aabf84528b50c)},
	{UINT128(0x0000000000000000, 0x987fbfe70b81a709)}, {UINT128(0x0000000000000000, 0xcd9f023f9c3a059f)},
	{UINT128(0xffffffffffffffff, 0x3260fdc063c5fa61)}, {UINT128(0x0000000000000000, 0x987fbfe70b81a709)},
	{UINT128(0x0000000000000000, 0xdb941a28cb71ec88)}, {UINT128(0x0000000000000000, 0x839c3cc917ff6cb5)},
	{UINT128(0xffffffffffffffff, 0x7c63c336e800934b)}, {UINT128(0x0000000000000000, 0xdb941a28cb71ec88)},
	{UINT128(0x0000000000000000, 0x3e33f2f642be355f)}, {UINT128(0x0000000000000000, 0xf853f7dc9186b953)},
	{UINT128(0xffffffffffffffff, 0x07ac08236e7946ad)}, {UINT128(0x0000000000000000, 0x3e33f2f642be355f)},
	{UINT128(0x0000000000000000, 0xf1090827b43725fe)}, {UINT128(0x0000000000000000, 0x563e69d6ac7f73f9)},
	{UINT128(0xffffffffffffffff, 0xa9c1962953808c07)}, {UINT128(0x0000000000000000, 0xf1090827b43725fe)},
	{UINT128(0x0000000000000000, 0x6d744027857300ae)}, {UINT128(0x0000000000000000, 0xe76bd7a1e63b9787)},
	{UINT128(0xffffffffffffffff, 0x1894285e19c46879)}, {UINT128(0x0000000000000000, 0x6d744027857300ae)},
	{UINT128(0x0000000000000000, 0xbdaef913557d76f1)}, {UINT128(0x0000000000000000, 0xabeb49a46764fd16)},
	{UINT128(0xffffffffffffffff, 0x5414b65b989b02ea)}, {UINT128(0x0000000000000000, 0xbdaef913557d76f1)},
	{UINT128(0x0000000000000000, 0x0c8fb2f886ec09f4)}, {UINT128(0x0000000000000000, 0xffb10f1bcb6bef1e)},
	{UINT128(0xffffffffffffffff, 0x004ef0e4349410e2)}, {UINT128(0x0000000000000000, 0x0c8fb2f886ec09f4)},
	{UINT128(0x0000000000000000, 0xffec4304266865da)}, {UINT128(0x0000000000000000, 0x0648557de8d99f7f)},
	{UINT128(0xffffffffffffffff, 0xf9b7aa8217266081)}, {UINT128(0x0000000000000000, 0xffec4304266865da)},
	{UINT128(0x0000000000000000, 0xb085baa8e966f6db)}, {UINT128(0x0000000000000000, 0xb96841bf7ffcb21b)},
	{UINT128(0xffffffffffffffff, 0x4697be4080034de5)}, {UINT128(0x0000000000000000, 0xb085baa8e966f6db)},
	{UINT128(0x0000000000000000, 0xea09a68a6e49cd63)}, {UINT128(0x0000000000000000, 0x67bde50ea3b628b7)},
	{UINT128(0xffffffffffffffff, 0x98421af15c49d749)}, {UINT128(0x0000000000000000, 0xea09a68a6e49cd63)},
	{UINT128(0x0000000000000000, 0x5c2214c3e9167abc)}, {UINT128(0x0000000000000000, 0xeed89db66611e308)},
	{UINT128(0xffffffffffffffff, 0x1127624999ee1cf8)}, {UINT128(0x0000000000000000, 0x5c2214c3e9167abc)},
	{UINT128(0x0000000000000000, 0xf9c79d63272c4629)}, {UINT128(0x0000000000000000, 0x381704d4fc9ec5fa)},
	{UINT128(0xffffffffffffffff, 0xc7e8fb2b03613a06)}, {UINT128(0x0000000000000000, 0xf9c79d63272c4629)},
	{UINT128(0x0000000000000000, 0x88f59aa0da591422)}, {UINT128(0x0000000000000000, 0xd84852c0a80ffcdc)},
	{UINT128(0xffffffffffffffff, 0x27b7ad3f57f00324)}, {UINT128(0x0000000000000000, 0x88f59aa0da591422)},
	{UINT128(0x0000000000000000, 0xd14d3d02313c0eee)}, {UINT128(0x0000000000000000, 0x93682a66e896f545)},
	{UINT128(0xffffffffffffffff, 0x6c97d59917690abb)}, {UINT128(0x0000000000000000, 0xd14d3d02313c0eee)},
	{UINT128(0x0000000000000000, 0x2bc42889167f8caa)}, {UINT128(0x0000000000000000, 0xfc3b27d38a5d49ac)},
	{UINT128(0xffffffffffffffff, 0x03c4d82c75a2b654)}, {UINT128(0x0000000000000000, 0x2bc42889167f8caa)},
	{UINT128(0x0000000000000000, 0xfe1323870cfe9a3e)}, {UINT128(0x0000000000000000, 0x1f564e56a9730e35)},
	{UINT128(0xffffffffffffffff, 0xe0a9b1a9568cf1cb)}, {UINT128(0x0000000000000000, 0xfe1323870cfe9a3e)},
	{UINT128(0x0000000000000000, 0x9d7fd1490285c9e4)}, {UINT128(0x0000000000000000, 0xc9d1124c931fda7b)},
	{UINT128(0xffffffffffffffff, 0x362eedb36ce02585)}, {UINT128(0x0000000000000000, 0x9d7fd1490285c9e4)},
	{UINT128(0x0000000000000000, 0xdebe05637ca94cfc)}, {UINT128(0x0000000000000000, 0x7e2e936fe26ae7ee)},
	{UINT128(0xffffffffffffffff, 0x81d16c901d951812)}, {UINT128(0x0000000000000000, 0xdebe05637ca94cfc)},
	{UINT128(0x0000000000000000, 0x4447498ac7d9dd83)}, {UINT128(0x0000000000000000, 0xf6ba073b424b19e9)},
	{UINT128(0xffffffffffffffff, 0x0945f8c4bdb4e617)}, {UINT128(0x0000000000000000, 0x4447498ac7d9dd83)},
	{UINT128(0x0000000000000000, 0xf314476247088f75)}, {UINT128(0x0000000000000000, 0x504d72505d98050d)},
	{UINT128(0xffffffffffffffff, 0xafb28dafa267faf3)}, {UINT128(0x0000000000000000, 0xf314476247088f75)},
	{UINT128(0x0000000000000000, 0x7319ba64c711785b)}, {UINT128(0x0000000000000000, 0xe4aa5909a08fa7b5)},
	{UINT128(0xffffffffffffffff, 0x1b55a6f65f70584b)}, {UINT128(0x0000000000000000, 0x7319ba64c711785b)},
	{UINT128(0x0000000000000000, 0xc1d8705ffcbb6e91)}, {UINT128(0x0000000000000000, 0xa73655df1f2f489f)},
	{UINT128(0xffffffffffffffff, 0x58c9aa20e0d0b761)}, {UINT128(0x0000000000000000, 0xc1d8705ffcbb6e91)},
	{UINT128(0x0000000000000000, 0x12d52092ce19f5cd)}, {UINT128(0x0000000000000000, 0xff4e6d680c41d0aa)},
	{UINT128(0xffffffffffffffff, 0x00b19297f3be2f56)}, {UINT128(0x0000000000000000, 0x12d52092ce19f5cd)},
	{UINT128(0x0000000000000000, 0xff4e6d680c41d0aa)}, {UINT128(0x0000000000000000, 0x12d52092ce19f5cd)},
	{UINT128(0xffffffffffffffff, 0xed2adf6d31e60a33)}, {UINT128(0x0000000000000000, 0xff4e6d680c41d0aa)},
	{UINT128(0x0000000000000000, 0xa73655df1f2f489f)}, {UINT128(0x0000000000000000, 0xc1d8705ffcbb6e91)},
	{UINT128(0xffffffffffffffff, 0x3e278fa00344916f)}, {UINT128(0x0000000000000000, 0xa73655df1f2f489f)},
	{UINT128(0x0000000000000000, 0xe4aa5909a08fa7b5)}, {UINT128(0x0000000000000000, 0x7319ba64c711785b)},
	{UINT128(0xffffffffffffffff, 0x8ce6459b38ee87a5)}, {UINT128(0x0000000000000000, 0xe4aa5909a08fa7b5)},
	{UINT128(0x0000000000000000, 0x504d72505d98050d)}, {UINT128(0x0000000000000000, 0xf314476247088f75)},
	{UINT128(0xffffffffffffffff, 0x0cebb89db8f7708b)}, {UINT128(0x0000000000000000, 0x504d72505d98050d)},
	{UINT128(0x0000000000000000, 0xf6ba073b424b19e9)}, {UINT128(0x0000000000000000, 0x4447498ac7d9dd83)},
	{UINT128(0xffffffffffffffff, 0xbbb8b6753826227d)}, {UINT128(0x0000000000000000, 0xf6ba073b424b19e9)},
	{UINT128(0x0000000000000000, 0x7e2e936fe26ae7ee)}, {UINT128(0x0000000000000000, 0xdebe05637ca94cfc)},
	{UINT128(0xffffffffffffffff, 0x2141fa9c8356b304)}, {UINT128(0x0000000000000000, 0x7e2e936fe26ae7ee)},
	{UINT128(0x0000000000000000, 0xc9d1124c931fda7b)}, {UINT128(0x0000000000000000, 0x9d7fd1490285c9e4)},
	{UINT128(0xffffffffffffffff, 0x62802eb6fd7a361c)}, {UINT128(0x0000000000000000, 0xc9d1124c931fda7b)},
	{UINT128(0x0000000000000000, 0x1f564e56a9730e35)}, {UINT128(0x0000000000000000, 0xfe1323870cfe9a3e)},
	{UINT128(0xffffffffffffffff, 0x01ecdc78f30165c2)}, {UINT128(0x0000000000000000, 0x1f564e56a9730e35)},
	{UINT128(0x0000000000000000, 0xfc3b27d38a5d49ac)}, {UINT128(0x0000000000000000, 0x2bc42889167f8caa)},
	{UINT128(0xffffffffffffffff, 0xd43bd776e9807356)}, {UINT128(0x0000000000000000, 0xfc3b27d38a5d49ac)},
	{UINT128(0x0000000000000000, 0x93682a66e896f545)}, {UINT128(0x0000000000000000, 0xd14d3d02313c0eee)},
	{UINT128(0xffffffffffffffff, 0x2eb2c2fdcec3f112)}, {UINT128(0x0000000000000000, 0x93682a66e896f545)},
	{UINT128(0x0000000000000000, 0xd84852c0a80ffcdc)}, {UINT128(0x0000000000000000, 0x88f59aa0da591422)},
	{UINT128(0xffffffffffffffff, 0x770a655f25a6ebde)}, {UINT128(0x0000000000000000, 0xd84852c0a80ffcdc)},
	{UINT128(0x0000000000000000, 0x381704d4fc9ec5fa)}, {UINT128(0x0000000000000000, 0xf9c79d63272c4629)},
	{UINT128(0xffffffffffffffff, 0x0638629cd8d3b9d7)}, {UINT128(0x0000000000000000, 0x381704d4fc9ec5fa)},
	{UINT128(0x0000000000000000, 0xeed89db66611e308)}, {UINT128(0x0000000000000000, 0x5c2214c3e9167abc)},
	{UINT128(0xffffffffffffffff, 0xa3ddeb3c16e98544)}, {UINT128(0x0000000000000000, 0xeed89db66611e308)},
	{UINT128(0x0000000000000000, 0x67bde50ea3b628b7)}, {UINT128(0x0000000000000000, 0xea09a68a6e49cd63)},
	{UINT128(0xffffffffffffffff, 0x15f6597591b6329d)}, {UINT128(0x0000000000000000, 0x67bde50ea3b628b7)},
	{UINT128(0x0000000000000000, 0xb96841bf7ffcb21b)}, {UINT128(0x0000000000000000, 0xb085baa8e966f6db)},
	{UINT128(0xffffffffffffffff, 0x4f7a455716990925)}, {UINT128(0x0000000000000000, 0xb96841bf7ffcb21b)},
	{UINT128(0x0000000000000000, 0x0648557de8d99f7f)}, {UINT128(0x0000000000000000, 0xffec4304266865da)},
	{UINT128(0xffffffffffffffff, 0x0013bcfbd9979a26)}, {UINT128(0x0000000000000000, 0x0648557de8d99f7f)},
	{UINT128(0x0000000000000000, 0xfffb10b4dc96dabc)}, {UINT128(0x0000000000000000, 0x03243a3f9bd8f08d)},
	{UINT128(0xffffffffffffffff, 0xfcdbc5c064270f73)}, {UINT128(0x0000000000000000, 0xfffb10b4dc96dabc)},
	{UINT128(0x0000000000000000, 0xb2c8c92f83c1eb88)}, {UINT128(0x0000000000000000, 0xb73a22a755457449)},
	{UINT128(0xffffffffffffffff, 0x48c5dd58aaba8bb7)}, {UINT128(0x0000000000000000, 0xb2c8c92f83c1eb88)},
	{UINT128(0x0000000000000000, 0xeb4b0b9e4f345618)}, {UINT128(0x0000000000000000, 0x64dca98ef24f5cb5)},
	{UINT128(0xffffffffffffffff, 0x9b2356710db0a34b)}, {UINT128(0x0000000000000000, 0xeb4b0b9e4f345618)},
	{UINT128(0x0000000000000000, 0x5f0ea4c477339c07)}, {UINT128(0x0000000000000000, 0xedb29311c504d653)},
	{UINT128(0xffffffffffffffff, 0x124d6cee3afb29ad)}, {UINT128(0x0000000000000000, 0x5f0ea4c477339c07)},
	{UINT128(0x0000000000000000, 0xfa7301d859796672)}, {UINT128(0x0000000000000000, 0x3505404b6008a13d)},
	{UINT128(0xffffffffffffffff, 0xcafabfb49ff75ec3)}, {UINT128(0x0000000000000000, 0xfa7301d859796672)},
	{UINT128(0x0000000000000000, 0x8b9a6b1ef6da4503)}, {UINT128(0x0000000000000000, 0xd695e4f10ea88571)},
	{UINT128(0xffffffffffffffff, 0x296a1b0ef1577a8f)}, {UINT128(0x0000000000000000, 0x8b9a6b1ef6da4503)},
	{UINT128(0x0000000000000000, 0xd31848d817d70e17)}, {UINT128(0x0000000000000000, 0x90d3ccc99f5ac58c)},
	{UINT128(0xffffffffffffffff, 0x6f2c333660a53a74)}, {UINT128(0x0000000000000000, 0xd31848d817d70e17)},
	{UINT128(0x0000000000000000, 0x2edbb3bca17e628f)}, {UINT128(0x0000000000000000, 0xfbaccd1d0903bb0a)},
	{UINT128(0xffffffffffffffff, 0x045332e2f6fc44f6)}, {UINT128(0x0000000000000000, 0x2edbb3bca17e628f)},
	{UINT128(0x0000000000000000, 0xfe70afeb6d33d6a3)}, {UINT128(0x0000000000000000, 0x1c3785c79ec2d4f6)},
	{UINT128(0xffffffffffffffff, 0xe3c87a38613d2b0a)}, {UINT128(0x0000000000000000, 0xfe70afeb6d33d6a3)},
	{UINT128(0x0000000000000000, 0x9ff6ca9a2ab6a26e)}, {UINT128(0x0000000000000000, 0xc7de651f7ca0674a)},
	{UINT128(0xffffffffffffffff, 0x38219ae0835f98b6)}, {UINT128(0x0000000000000000, 0x9ff6ca9a2ab6a26e)},
	{UINT128(0x0000000000000000, 0xe046213392aa486d)}, {UINT128(0x0000000000000000, 0x7b70654bbde35623)},
	{UINT128(0xffffffffffffffff, 0x848f9ab4421ca9dd)}, {UINT128(0x0000000000000000, 0xe046213392aa486d)},
	{UINT128(0x0000000000000000, 0x474d10fd336cf747)}, {UINT128(0x0000000000000000, 0xf5dec646f85ba1c7)},
	{UINT128(0xffffffffffffffff, 0x0a2139b907a45e39)}, {UINT128(0x0000000000000000, 0x474d10fd336cf747)},
	{UINT128(0x0000000000000000, 0xf40bdd5a66886630)}, {UINT128(0x0000000000000000, 0x4d50430b860546c4)},
	{UINT128(0xffffffffffffffff, 0xb2afbcf479fab93c)}, {UINT128(0x0000000000000000, 0xf40bdd5a66886630)},
	{UINT128(0x0000000000000000, 0x75e5dd6e1b8e2556)}, {UINT128(0x0000000000000000, 0xe33c59a4439cd8ed)},
	{UINT128(0xffffffffffffffff, 0x1cc3a65bbc632713)}, {UINT128(0x0000000000000000, 0x75e5dd6e1b8e2556)},
	{UINT128(0x0000000000000000, 0xc3e2007dd175f5a5)}, {UINT128(0x0000000000000000, 0xa4d224dcd849c5b1)},
	{UINT128(0xffffffffffffffff, 0x5b2ddb2327b63a4f)}, {UINT128(0x0000000000000000, 0xc3e2007dd175f5a5)},
	{UINT128(0x0000000000000000, 0x15f6d00a9aa418c2)}, {UINT128(0x0000000000000000, 0xff0e57e5ead848d2)},
	{UINT128(0xffffffffffffffff, 0x00f1a81a1527b72e)}, {UINT128(0x0000000000000000, 0x15f6d00a9aa418c2)},
	{UINT128(0x0000000000000000, 0xff84ab2c738d6a04)}, {UINT128(0x0000000000000000, 0x0fb2b73cfc106ff7)},
	{UINT128(0xffffffffffffffff, 0xf04d48c303ef9009)}, {UINT128(0x0000000000000000, 0xff84ab2c738d6a04)},
	{UINT128(0x0000000000000000, 0xa99414951aacae5f)}, {UINT128(0x0000000000000000, 0xbfc7671ab8bb84c7)},
	{UINT128(0xffffffffffffffff, 0x403898e547447b39)}, {UINT128(0x0000000000000000, 0xa99414951aacae5f)},
	{UINT128(0x0000000000000000, 0xe60f879fe7e2e1e6)}, {UINT128(0x0000000000000000, 0x70492760047b9a7f)},
	{UINT128(0xffffffffffffffff, 0x8fb6d89ffb846581)}, {UINT128(0x0000000000000000, 0xe60f879fe7e2e1e6)},
	{UINT128(0x0000000000000000, 0x53478909e39da893)}, {UINT128(0x0000000000000000, 0xf21352595e0bf351)},
	{UINT128(0xffffffffffffffff, 0x0decada6a1f40caf)}, {UINT128(0x0000000000000000, 0x53478909e39da893)},
	{UINT128(0x0000000000000000, 0xf78bc51f239e12c7)}, {UINT128(0x0000000000000000, 0x413ee038dff6b7fe)},
	{UINT128(0xffffffffffffffff, 0xbec11fc720094802)}, {UINT128(0x0000000000000000, 0xf78bc51f239e12c7)},
	{UINT128(0x0000000000000000, 0x80e7e43a61f5b6cc)}, {UINT128(0x0000000000000000, 0xdd2d5339ac8692fe)},
	{UINT128(0xffffffffffffffff, 0x22d2acc653796d02)}, {UINT128(0x0000000000000000, 0x80e7e43a61f5b6cc)},
	{UINT128(0x0000000000000000, 0xcbbbf7a63eba0dd6)}, {UINT128(0x0000000000000000, 0x9b02c58832cf95c1)},
	{UINT128(0xffffffffffffffff, 0x64fd3a77cd306a3f)}, {UINT128(0x0000000000000000, 0xcbbbf7a63eba0dd6)},
	{UINT128(0x0000000000000000, 0x2273e19db5eaed57)}, {UINT128(0x0000000000000000, 0xfdabcb8caeba091c)},
	{UINT128(0xffffffffffffffff, 0x025434735145f6e4)}, {UINT128(0x0000000000000000, 0x2273e19db5eaed57)},
	{UINT128(0x0000000000000000, 0xfcbfc926484cd43b)}, {UINT128(0x0000000000000000, 0x28aaed62527cb3b6)},
	{UINT128(0xffffffffffffffff, 0xd755129dad834c4a)}, {UINT128(0x0000000000000000, 0xfcbfc926484cd43b)},
	{UINT128(0x0000000000000000, 0x95f6d92fd79f4fbb)}, {UINT128(0x0000000000000000, 0xcf7a1f794d7ca1b2)},
	{UINT128(0xffffffffffffffff, 0x3085e086b2835e4e)}, {UINT128(0x0000000000000000, 0x95f6d92fd79f4fbb)},
	{UINT128(0x0000000000000000, 0xd9f269f7aab88c2a)}, {UINT128(0x0000000000000000, 0x864b826aec4c74e6)},
	{UINT128(0xffffffffffffffff, 0x79b47d9513b38b1a)}, {UINT128(0x0000000000000000, 0xd9f269f7aab88c2a)},
	{UINT128(0x0000000000000000, 0x3b269fca8a8622bb)}, {UINT128(0x0000000000000000, 0xf91297bbb1d6cdbf)},
	{UINT128(0xffffffffffffffff, 0x06ed68444e293241)}, {UINT128(0x0000000000000000, 0x3b269fca8a8622bb)},
	{UINT128(0x0000000000000000, 0xeff573116df1555e)}, {UINT128(0x0000000000000000, 0x5931f774fc9f1845)},
	{UINT128(0xffffffffffffffff, 0xa6ce088b0360e7bb)}, {UINT128(0x0000000000000000, 0xeff573116df1555e)},
	{UINT128(0x0000000000000000, 0x6a9b20adb4ff262b)}, {UINT128(0x0000000000000000, 0xe8bf3ba1f1aedfbc)},
	{UINT128(0xffffffffffffffff, 0x1740c45e0e512044)}, {UINT128(0x0000000000000000, 0x6a9b20adb4ff262b)},
	{UINT128(0x0000000000000000, 0xbb8f3af81b93095d)}, {UINT128(0x0000000000000000, 0xae3bddf3280c620e)},
	{UINT128(0xffffffffffffffff, 0x51c4220cd7f39df2)}, {UINT128(0x0000000000000000, 0xbb8f3af81b93095d)},
	{UINT128(0x0000000000000000, 0x096c32baca2ae68c)}, {UINT128(0x0000000000000000, 0xffd3977ff7bae4ea)},
	{UINT128(0xffffffffffffffff, 0x002c688008451b16)}, {UINT128(0x0000000000000000, 0x096c32baca2ae68c)},
	{UINT128(0x0000000000000000, 0xffd3977ff7bae4ea)}, {UINT128(0x0000000000000000, 0x096c32baca2ae68c)},
	{UINT128(0xffffffffffffffff, 0xf693cd4535d51974)}, {UINT128(0x0000000000000000, 0xffd3977ff7bae4ea)},
	{UINT128(0x0000000000000000, 0xae3bddf3280c620e)}, {UINT128(0x0000000000000000, 0xbb8f3af81b93095d)},
	{UINT128(0xffffffffffffffff, 0x4470c507e46cf6a3)}, {UINT128(0x0000000000000000, 0xae3bddf3280c620e)},
	{UINT128(0x0000000000000000, 0xe8bf3ba1f1aedfbc)}, {UINT128(0x0000000000000000, 0x6a9b20adb4ff262b)},
	{UINT128(0xffffffffffffffff, 0x9564df524b00d9d5)}, {UINT128(0x0000000000000000, 0xe8bf3ba1f1aedfbc)},
	{UINT128(0x0000000000000000, 0x5931f774fc9f1845)}, {UINT128(0x0000000000000000, 0xeff573116df1555e)},
	{UINT128(0xffffffffffffffff, 0x100a8cee920eaaa2)}, {UINT128(0x0000000000000000, 0x5931f774fc9f1845)},
	{UINT128(0x0000000000000000, 0xf91297bbb1d6cdbf)}, {UINT128(0x0000000000000000, 0x3b269fca8a8622bb)},
	{UINT128(0xffffffffffffffff, 0xc4d960357579dd45)}, {UINT128(0x0000000000000000, 0xf91297bbb1d6cdbf)},
	{UINT128(0x0000000000000000, 0x864b826aec4c74e6)}, {UINT128(0x0000000000000000, 0xd9f269f7aab88c2a)},
	{UINT128(0xffffffffffffffff, 0x260d9608554773d6)}, {UINT128(0x0000000000000000, 0x864b826aec4c74e6)},
	{UINT128(0x0000000000000000, 0xcf7a1f794d7ca1b2)}, {UINT128(0x0000000000000000, 0x95f6d92fd79f4fbb)},
	{UINT128(0xffffffffffffffff, 0x6a0926d02860b045)}, {UINT128(0x0000000000000000, 0xcf7a1f794d7ca1b2)},
	{UINT128(0x0000000000000000, 0x28aaed62527cb3b6)}, {UINT128(0x0000000000000000, 0xfcbfc926484cd43b)},
	{UINT128(0xffffffffffffffff, 0x034036d9b7b32bc5)}, {UINT128(0x0000000000000000, 0x28aaed62527cb3b6)},
	{UINT128(0x0000000000000000, 0xfdabcb8caeba091c)}, {UINT128(0x0000000000000000, 0x2273e19db5eaed57)},
	{UINT128(0xffffffffffffffff, 0xdd8c1e624a1512a9)}, {UINT128(0x0000000000000000, 0xfdabcb8caeba091c)},
	{UINT128(0x0000000000000000, 0x9b02c58832cf95c1)}, {UINT128(0x0000000000000000, 0xcbbbf7a63eba0dd6)},
	{UINT128(0xffffffffffffffff, 0x34440859c145f22a)}, {UINT128(0x0000000000000000, 0x9b02c58832cf95c1)},
	{UINT128(0x0000000000000000, 0xdd2d5339ac8692fe)}, {UINT128(0x0000000000000000, 0x80e7e43a61f5b6cc)},
	{UINT128(0xffffffffffffffff, 0x7f181bc59e0a4934)}, {UINT128(0x0000000000000000, 0xdd2d5339ac8692fe)},
	{UINT128(0x0000000000000000, 0x413ee038dff6b7fe)}, {UINT128(0x0000000000000000, 0xf78bc51f239e12c7)},
	{UINT128(0xffffffffffffffff, 0x08743ae0dc61ed39)}, {UINT128(0x0000000000000000, 0x413ee038dff6b7fe)},
	{UINT128(0x0000000000000000, 0xf21352595e0bf351)}, {UINT128(0x0000000000000000, 0x53478909e39da893)},
	{UINT128(0xffffffffffffffff, 0xacb876f61c62576d)}, {UINT128(0x0000000000000000, 0xf21352595e0bf351)},
	{UINT128(0x0000000000000000, 0x70492760047b9a7f)}, {UINT128(0x0000000000000000, 0xe60f879fe7e2e1e6)},
	{UINT128(0xffffffffffffffff, 0x19f07860181d1e1a)}, {UINT128(0x0000000000000000, 0x70492760047b9a7f)},
	{UINT128(0x0000000000000000, 0xbfc7671ab8bb84c7)}, {UINT128(0x0000000000000000, 0xa99414951aacae5f)},
	{UINT128(0xffffffffffffffff, 0x566beb6ae55351a1)}, {UINT128(0x0000000000000000, 0xbfc7671ab8bb84c7)},
	{UINT128(0x0000000000000000, 0x0fb2b73cfc106ff7)}, {UINT128(0x0000000000000000, 0xff84ab2c738d6a04)},
	{UINT128(0xffffffffffffffff, 0x007b54d38c7295fc)}, {UINT128(0x0000000000000000, 0x0fb2b73cfc106ff7)},
	{UINT128(0x0000000000000000, 0xff0e57e5ead848d2)}, {UINT128(0x0000000000000000, 0x15f6d00a9aa418c2)},
	{UINT128(0xffffffffffffffff, 0xea092ff5655be73e)}, {UINT128(0x0000000000000000, 0xff0e57e5ead848d2)},
	{UINT128(0x0000000000000000, 0xa4d224dcd849c5b1)}, {UINT128(0x0000000000000000, 0xc3e2007dd175f5a5)},
	{UINT128(0xffffffffffffffff, 0x3c1dff822e8a0a5b)}, {UINT128(0x0000000000000000, 0xa4d224dcd849c5b1)},
	{UINT128(0x0000000000000000, 0xe33c59a4439cd8ed)}, {UINT128(0x0000000000000000, 0x75e5dd6e1b8e2556)},
	{UINT128(0xffffffffffffffff, 0x8a1a2291e471daaa)}, {UINT128(0x0000000000000000, 0xe33c59a4439cd8ed)},
	{UINT128(0x0000000000000000, 0x4d50430b860546c4)}, {UINT128(0x0000000000000000, 0xf40bdd5a66886630)},
	{UINT128(0xffffffffffffffff, 0x0bf422a5997799d0)}, {UINT128(0x0000000000000000, 0x4d50430b860546c4)},
	{UINT128(0x0000000000000000, 0xf5dec646f85ba1c7)}, {UINT128(0x0000000000000000, 0x474d10fd336cf747)},
	{UINT128(0xffffffffffffffff, 0xb8b2ef02cc9308b9)}, {UINT128(0x0000000000000000, 0xf5dec646f85ba1c7)},
	{UINT128(0x0000000000000000, 0x7b70654bbde35623)}, {UINT128(0x0000000000000000, 0xe046213392aa486d)},
	{UINT128(0xffffffffffffffff, 0x1fb9decc6d55b793)}, {UINT128(0x0000000000000000, 0x7b70654bbde35623)},
	{UINT128(0x0000000000000000, 0xc7de651f7ca0674a)}, {UINT128(0x0000000000000000, 0x9ff6ca9a2ab6a26e)},
	{UINT128(0xffffffffffffffff, 0x60093565d5495d92)}, {UINT128(0x0000000000000000, 0xc7de651f7ca0674a)},
	{UINT128(0x0000000000000000, 0x1c3785c79ec2d4f6)}, {UINT128(0x0000000000000000, 0xfe70afeb6d33d6a3)},
	{UINT128(0xffffffffffffffff, 0x018f501492cc295d)}, {UINT128(0x0000000000000000, 0x1c3785c79ec2d4f6)},
	{UINT128(0x0000000000000000, 0xfbaccd1d0903bb0a)}, {UINT128(0x0000000000000000, 0x2edbb3bca17e628f)},
	{UINT128(0xffffffffffffffff, 0xd1244c435e819d71)}, {UINT128(0x0000000000000000, 0xfbaccd1d0903bb0a)},
	{UINT128(0x0000000000000000, 0x90d3ccc99f5ac58c)}, {UINT128(0x0000000000000000, 0xd31848d817d70e17)},
	{UINT128(0xffffffffffffffff, 0x2ce7b727e828f1e9)}, {UINT128(0x0000000000000000, 0x90d3ccc99f5ac58c)},
	{UINT128(0x0000000000000000, 0xd695e4f10ea88571)}, {UINT128(0x0000000000000000, 0x8b9a6b1ef6da4503)},
	{UINT128(0xffffffffffffffff, 0x746594e10925bafd)}, {UINT128(0x0000000000000000, 0xd695e4f10ea88571)},
	{UINT128(0x0000000000000000, 0x3505404b6008a13d)}, {UINT128(0x0000000000000000, 0xfa7301d859796672)},
	{UINT128(0xffffffffffffffff, 0x058cfe27a686998e)}, {UINT128(0x0000000000000000, 0x3505404b6008a13d)},
	{UINT128(0x0000000000000000, 0xedb29311c504d653)}, {UINT128(0x0000000000000000, 0x5f0ea4c477339c07)},
	{UINT128(0xffffffffffffffff, 0xa0f15b3b88cc63f9)}, {UINT128(0x0000000000000000, 0xedb29311c504d653)},
	{UINT128(0x0000000000000000, 0x64dca98ef24f5cb5)}, {UINT128(0x0000000000000000, 0xeb4b0b9e4f345618)},
	{UINT128(0xffffffffffffffff, 0x14b4f461b0cba9e8)}, {UINT128(0x0000000000000000, 0x64dca98ef24f5cb5)},
	{UINT128(0x0000000000000000, 0xb73a22a755457449)}, {UINT128(0x0000000000000000, 0xb2c8c92f83c1eb88)},
	{UINT128(0xffffffffffffffff, 0x4d3736d07c3e1478)}, {UINT128(0x0000000000000000, 0xb73a22a755457449)},
	{UINT128(0x0000000000000000, 0x03243a3f9bd8f08d)}, {UINT128(0x0000000000000000, 0xfffb10b4dc96dabc)},
	{UINT128(0xffffffffffffffff, 0x0004ef4b23692544)}, {UINT128(0x0000000000000000, 0x03243a3f9bd8f08d)},
	{UINT128(0x0000000000000000, 0xfffec42c7454926c)}, {UINT128(0x0000000000000000, 0x01921f0fe6700712)},
	{UINT128(0xffffffffffffffff, 0xfe6de0f0198ff8ee)}, {UINT128(0x0000000000000000, 0xfffec42c7454926c)},
	{UINT128(0x0000000000000000, 0xb3e7bc248d78802f)}, {UINT128(0x0000000000000000, 0xb6206b9e0c13a893)},
	{UINT128(0xffffffffffffffff, 0x49df9461f3ec576d)}, {UINT128(0x0000000000000000, 0xb3e7bc248d78802f)},
	{UINT128(0x0000000000000000, 0xebe85815c767cb01)}, {UINT128(0x0000000000000000, 0x636a94bb2292bf47)},
	{UINT128(0xffffffffffffffff, 0x9c956b44dd6d40b9)}, {UINT128(0x0000000000000000, 0xebe85815c767cb01)},
	{UINT128(0x0000000000000000, 0x60838ec13aab0fce)}, {UINT128(0x0000000000000000, 0xed1c1d4b344c3d50)},
	{UINT128(0xffffffffffffffff, 0x12e3e2b4cbb3c2b0)}, {UINT128(0x0000000000000000, 0x60838ec13aab0fce)},
	{UINT128(0x0000000000000000, 0xfac5158bc4f42120)}, {UINT128(0x0000000000000000, 0x337b97e5b886cccc)},
	{UINT128(0xffffffffffffffff, 0xcc84681a47793334)}, {UINT128(0x0000000000000000, 0xfac5158bc4f42120)},
	{UINT128(0x0000000000000000, 0x8cead04f95cdbf67)}, {UINT128(0x0000000000000000, 0xd5b992c8b606a352)},
	{UINT128(0xffffffffffffffff, 0x2a466d3749f95cae)}, {UINT128(0x0000000000000000, 0x8cead04f95cdbf67)},
	{UINT128(0x0000000000000000, 0xd3fac294ff34e4d1)}, {UINT128(0x0000000000000000, 0x8f87845de430d778)},
	{UINT128(0xffffffffffffffff, 0x70787ba21bcf2888)}, {UINT128(0x0000000000000000, 0xd3fac294ff34e4d1)},
	{UINT128(0x0000000000000000, 0x3066cdd138c98b74)}, {UINT128(0x0000000000000000, 0xfb61fbefadddb986)},
	{UINT128(0xffffffffffffffff, 0x049e04105222467a)}, {UINT128(0x0000000000000000, 0x3066cdd138c98b74)},
	{UINT128(0x0000000000000000, 0xfe9bc8a1105c22a6)}, {UINT128(0x0000000000000000, 0x1aa7b724495c037a)},
	{UINT128(0xffffffffffffffff, 0xe55848dbb6a3fc86)}, {UINT128(0x0000000000000000, 0xfe9bc8a1105c22a6)},
	{UINT128(0x0000000000000000, 0xa12ff8bc735d8af7)}, {UINT128(0x0000000000000000, 0xc6e22998b4c6608f)},
	{UINT128(0xffffffffffffffff, 0x391dd6674b399f71)}, {UINT128(0x0000000000000000, 0xa12ff8bc735d8af7)},
	{UINT128(0x0000000000000000, 0xe106f1fd4b8d7c97)}, {UINT128(0x0000000000000000, 0x7a0f83abe1444f43)},
	{UINT128(0xffffffffffffffff, 0x85f07c541ebbb0bd)}, {UINT128(0x0000000000000000, 0xe106f1fd4b8d7c97)},
	{UINT128(0x0000000000000000, 0x48ceeeaf0eedc59a)}, {UINT128(0x0000000000000000, 0xf56d97473d446cdb)},
	{UINT128(0xffffffffffffffff, 0x0a9268b8c2bb9325)}, {UINT128(0x0000000000000000, 0x48ceeeaf0eedc59a)},
	{UINT128(0x0000000000000000, 0xf48421b0efbf939c)}, {UINT128(0x0000000000000000, 0x4bd08b6bb00e1a8b)},
	{UINT128(0xffffffffffffffff, 0xb42f74944ff1e575)}, {UINT128(0x0000000000000000, 0xf48421b0efbf939c)},
	{UINT128(0x0000000000000000, 0x774a3c5207315fc4)}, {UINT128(0x0000000000000000, 0xe28210095b483752)},
	{UINT128(0xffffffffffffffff, 0x1d7deff6a4b7c8ae)}, {UINT128(0x0000000000000000, 0x774a3c5207315fc4)},
	{UINT128(0x0000000000000000, 0xc4e3f4d26ea553b7)}, {UINT128(0x0000000000000000, 0xa39da8dcc39a38e6)},
	{UINT128(0xffffffffffffffff, 0x5c6257233c65c71a)}, {UINT128(0x0000000000000000, 0xc4e3f4d26ea553b7)},
	{UINT128(0x0000000000000000, 0x1787586a5d5b2158)}, {UINT128(0x0000000000000000, 0xfeea9cff8fa2ae55)},
	{UINT128(0xffffffffffffffff, 0x01156300705d51ab)}, {UINT128(0x0000000000000000, 0x1787586a5d5b2158)},
	{UINT128(0x0000000000000000, 0xff9c187c6abade6b)}, {UINT128(0x0000000000000000, 0x0e214689606bf168)},
	{UINT128(0xffffffffffffffff, 0xf1deb9769f940e98)}, {UINT128(0x0000000000000000, 0xff9c187c6abade6b)},
	{UINT128(0x0000000000000000, 0xaac081c4ba89ba8b)}, {UINT128(0x0000000000000000, 0xbebc1b6619ed9117)},
	{UINT128(0xffffffffffffffff, 0x4143e499e6126ee9)}, {UINT128(0x0000000000000000, 0xaac081c4ba89ba8b)},
	{UINT128(0x0000000000000000, 0xe6becc4c5997af07)}, {UINT128(0x0000000000000000, 0x6edf3c8c12f404d0)},
	{UINT128(0xffffffffffffffff, 0x9120c373ed0bfb30)}, {UINT128(0x0000000000000000, 0xe6becc4c5997af07)},
	{UINT128(0x0000000000000000, 0x54c36202bcf08e06)}, {UINT128(0x0000000000000000, 0xf18f574386712644)},
	{UINT128(0xffffffffffffffff, 0x0e70a8bc798ed9bc)}, {UINT128(0x0000000000000000, 0x54c36202bcf08e06)},
	{UINT128(0x0000000000000000, 0xf7f110605caf6391)}, {UINT128(0x0000000000000000, 0x3fb9b835bfdbf16a)},
	{UINT128(0xffffffffffffffff, 0xc04647ca40240e96)}, {UINT128(0x0000000000000000, 0xf7f110605caf6391)},
	{UINT128(0x0000000000000000, 0x8242b1357110d373)}, {UINT128(0x0000000000000000, 0xdc61c693a82745d6)},
	{UINT128(0xffffffffffffffff, 0x239e396c57d8ba2a)}, {UINT128(0x0000000000000000, 0x8242b1357110d373)},
	{UINT128(0x0000000000000000, 0xccae7976c0691178)}, {UINT128(0x0000000000000000, 0x99c200686472b4a9)},
	{UINT128(0xffffffffffffffff, 0x663dff979b8d4b57)}, {UINT128(0x0000000000000000, 0xccae7976c0691178)},
	{UINT128(0x0000000000000000, 0x24022da9d8f79d6e)}, {UINT128(0x0000000000000000, 0xfd747472367dd6c6)},
	{UINT128(0xffffffffffffffff, 0x028b8b8dc982293a)}, {UINT128(0x0000000000000000, 0x24022da9d8f79d6e)},
	{UINT128(0x0000000000000000, 0xfcfe72ad6d9641f3)}, {UINT128(0x0000000000000000, 0x271db7619b1a2774)},
	{UINT128(0xffffffffffffffff, 0xd8e2489e64e5d88c)}, {UINT128(0x0000000000000000, 0xfcfe72ad6d9641f3)},
	{UINT128(0x0000000000000000, 0x973c071f4750b49d)}, {UINT128(0x0000000000000000, 0xce8d8faf5406ab8c)},
	{UINT128(0xffffffffffffffff, 0x31727050abf95474)}, {UINT128(0x0000000000000000, 0x973c071f4750b49d)},
	{UINT128(0x0000000000000000, 0xdac44ff490a02711)}, {UINT128(0x0000000000000000, 0x84f483a0be2f0404)},
	{UINT128(0xffffffffffffffff, 0x7b0b7c5f41d0fbfc)}, {UINT128(0x0000000000000000, 0xdac44ff490a02711)},
	{UINT128(0x0000000000000000, 0x3cad943c203436a4)}, {UINT128(0x0000000000000000, 0xf8b47a9fb902e76d)},
	{UINT128(0xffffffffffffffff, 0x074b856046fd1893)}, {UINT128(0x0000000000000000, 0x3cad943c203436a4)},
	{UINT128(0x0000000000000000, 0xf08066514c055f7f)}, {UINT128(0x0000000000000000, 0x57b89cde7a9a4d64)},
	{UINT128(0xffffffffffffffff, 0xa84763218565b29c)}, {UINT128(0x0000000000000000, 0xf08066514c055f7f)},
	{UINT128(0x0000000000000000, 0x6c0835b1fd002451)}, {UINT128(0x0000000000000000, 0xe816a7f595ec9233)},
	{UINT128(0xffffffffffffffff, 0x17e9580a6a136dcd)}, {UINT128(0x0000000000000000, 0x6c0835b1fd002451)},
	{UINT128(0x0000000000000000, 0xbca002ba7aaf25eb)}, {UINT128(0x0000000000000000, 0xad146952eb9282b0)},
	{UINT128(0xffffffffffffffff, 0x52eb96ad146d7d50)}, {UINT128(0x0000000000000000, 0xbca002ba7aaf25eb)},
	{UINT128(0x0000000000000000, 0x0afe00694866a1b5)}, {UINT128(0x0000000000000000, 0xffc38ed6dc0ef98c)},
	{UINT128(0xffffffffffffffff, 0x003c712923f10674)}, {UINT128(0x0000000000000000, 0x0afe00694866a1b5)},
	{UINT128(0x0000000000000000, 0xffe128ef8e9fc17b)}, {UINT128(0x0000000000000000, 0x07da4dcc7473c040)},
	{UINT128(0xffffffffffffffff, 0xf825b2338b8c3fc0)}, {UINT128(0x0000000000000000, 0xffe128ef8e9fc17b)},
	{UINT128(0x0000000000000000, 0xaf61a4ac1b83a1df)}, {UINT128(0x0000000000000000, 0xba7ca46d46946803)},
	{UINT128(0xffffffffffffffff, 0x45835b92b96b97fd)}, {UINT128(0x0000000000000000, 0xaf61a4ac1b83a1df)},
	{UINT128(0x0000000000000000, 0xe9659107077cf610)}, {UINT128(0x0000000000000000, 0x692d049f7a879242)},
	{UINT128(0xffffffffffffffff, 0x96d2fb6085786dbe)}, {UINT128(0x0000000000000000, 0xe9659107077cf610)},
	{UINT128(0x0000000000000000, 0x5aaa75f71df85ac8)}, {UINT128(0x0000000000000000, 0xef682fbef23ecda7)},
	{UINT128(0xffffffffffffffff, 0x1097d0410dc13259)}, {UINT128(0x0000000000000000, 0x5aaa75f71df85ac8)},
	{UINT128(0x0000000000000000, 0xf96e4e4844d4e82b)}, {UINT128(0x0000000000000000, 0x399f196625650c49)},
	{UINT128(0xffffffffffffffff, 0xc660e699da9af3b7)}, {UINT128(0x0000000000000000, 0xf96e4e4844d4e82b)},
	{UINT128(0x0000000000000000, 0x87a135d95473ec8a)}, {UINT128(0x0000000000000000, 0xd91e6a38009da15b)},
	{UINT128(0xffffffffffffffff, 0x26e195c7ff625ea5)}, {UINT128(0x0000000000000000, 0x87a135d95473ec8a)},
	{UINT128(0x0000000000000000, 0xd064af55d7c9b43f)}, {UINT128(0x0000000000000000, 0x94b0393b14e54157)},
	{UINT128(0xffffffffffffffff, 0x6b4fc6c4eb1abea9)}, {UINT128(0x0000000000000000, 0xd064af55d7c9b43f)},
	{UINT128(0x0000000000000000, 0x2a37bf0b2f8be3f4)}, {UINT128(0x0000000000000000, 0xfc7eaffd720ed674)},
	{UINT128(0xffffffffffffffff, 0x038150028df1298c)}, {UINT128(0x0000000000000000, 0x2a37bf0b2f8be3f4)},
	{UINT128(0x0000000000000000, 0xfde0b0bf220c2fd5)}, {UINT128(0x0000000000000000, 0x20e5408f75063a47)},
	{UINT128(0xffffffffffffffff, 0xdf1abf708af9c5b9)}, {UINT128(0x0000000000000000, 0xfde0b0bf220c2fd5)},
	{UINT128(0x0000000000000000, 0x9c420c2eff590e60)}, {UINT128(0x0000000000000000, 0xcac77f24736eb554)},
	{UINT128(0xffffffffffffffff, 0x353880db8c914aac)}, {UINT128(0x0000000000000000, 0x9c420c2eff590e60)},
	{UINT128(0x0000000000000000, 0xddf6be249c075038)}, {UINT128(0x0000000000000000, 0x7f8bd92f9c483ed7)},
	{UINT128(0xffffffffffffffff, 0x807426d063b7c129)}, {UINT128(0x0000000000000000, 0xddf6be249c075038)},
	{UINT128(0x0000000000000000, 0x42c3673f6f6e404f)}, {UINT128(0x0000000000000000, 0xf7241712d4edde4a)},
	{UINT128(0xffffffffffffffff, 0x08dbe8ed2b1221b6)}, {UINT128(0x0000000000000000, 0x42c3673f6f6e404f)},
	{UINT128(0x0000000000000000, 0xf294f82394ffe321)}, {UINT128(0x0000000000000000, 0x51cae2955c414eff)},
	{UINT128(0xffffffffffffffff, 0xae351d6aa3beb101)}, {UINT128(0x0000000000000000, 0xf294f82394ffe321)},
	{UINT128(0x0000000000000000, 0x71b1fd265c002a42)}, {UINT128(0x0000000000000000, 0xe55e0b4d05c80389)},
	{UINT128(0xffffffffffffffff, 0x1aa1f4b2fa37fc77)}, {UINT128(0x0000000000000000, 0x71b1fd265c002a42)},
	{UINT128(0x0000000000000000, 0xc0d0d99dabd65d45)}, {UINT128(0x0000000000000000, 0xa86604facd04d96a)},
	{UINT128(0xffffffffffffffff, 0x5799fb0532fb2696)}, {UINT128(0x0000000000000000, 0xc0d0d99dabd65d45)},
	{UINT128(0x0000000000000000, 0x11440134d709b281)}, {UINT128(0x0000000000000000, 0xff6ac765b39e1e1a)},
	{UINT128(0xffffffffffffffff, 0x0095389a4c61e1e6)}, {UINT128(0x0000000000000000, 0x11440134d709b281)},
	{UINT128(0x0000000000000000, 0xff2f9d7971ca0365)}, {UINT128(0x0000000000000000, 0x14661179272095ae)},
	{UINT128(0xffffffffffffffff, 0xeb99ee86d8df6a52)}, {UINT128(0x0000000000000000, 0xff2f9d7971ca0365)},
	{UINT128(0x0000000000000000, 0xa6050a2f6000204a)}, {UINT128(0x0000000000000000, 0xc2de28d74ac6628c)},
	{UINT128(0xffffffffffffffff, 0x3d21d728b5399d74)}, {UINT128(0x0000000000000000, 0xa6050a2f6000204a)},
	{UINT128(0x0000000000000000, 0xe3f4729119e798da)}, {UINT128(0x0000000000000000, 0x74805ba3a76d6b1c)},
	{UINT128(0xffffffffffffffff, 0x8b7fa45c589294e4)}, {UINT128(0x0000000000000000, 0xe3f4729119e798da)},
	{UINT128(0x0000000000000000, 0x4ecf3be81052ddf4)}, {UINT128(0x0000000000000000, 0xf3913edb54ba2243)},
	{UINT128(0xffffffffffffffff, 0x0c6ec124ab45ddbd)}, {UINT128(0x0000000000000000, 0x4ecf3be81052ddf4)},
	{UINT128(0x0000000000000000, 0xf64d969e1dfc211a)}, {UINT128(0x0000000000000000, 0x45ca835dd945dc92)},
	{UINT128(0xffffffffffffffff, 0xba357ca226ba236e)}, {UINT128(0x0000000000000000, 0xf64d969e1dfc211a)},
	{UINT128(0x0000000000000000, 0x7cd01658ff419d07)}, {UINT128(0x0000000000000000, 0xdf83270a9bbee891)},
	{UINT128(0xffffffffffffffff, 0x207cd8f56441176f)}, {UINT128(0x0000000000000000, 0x7cd01658ff419d07)},
	{UINT128(0x0000000000000000, 0xc8d8b37ea4ed0f63)}, {UINT128(0x0000000000000000, 0x9ebc11c62c1a1dfc)},
	{UINT128(0xffffffffffffffff, 0x6143ee39d3e5e204)}, {UINT128(0x0000000000000000, 0xc8d8b37ea4ed0f63)},
	{UINT128(0x0000000000000000, 0x1dc70ecbae9fc913)}, {UINT128(0x0000000000000000, 0xfe432367f5b90a63)},
	{UINT128(0xffffffffffffffff, 0x01bcdc980a46f59d)}, {UINT128(0x0000000000000000, 0x1dc70ecbae9fc913)},
	{UINT128(0x0000000000000000, 0xfbf5314f31eb7376)}, {UINT128(0x0000000000000000, 0x2d502609ec95647d)},
	{UINT128(0xffffffffffffffff, 0xd2afd9f6136a9b83)}, {UINT128(0x0000000000000000, 0xfbf5314f31eb7376)},
	{UINT128(0x0000000000000000, 0x921eafdcc560f9c6)}, {UINT128(0x0000000000000000, 0xd233c6408cd64237)},
	{UINT128(0xffffffffffffffff, 0x2dcc39bf7329bdc9)}, {UINT128(0x0000000000000000, 0x921eafdcc560f9c6)},
	{UINT128(0x0000000000000000, 0xd77025a1e0a39d8c)}, {UINT128(0x0000000000000000, 0x8a48ad799b6759f4)},
	{UINT128(0xffffffffffffffff, 0x75b752866498a60c)}, {UINT128(0x0000000000000000, 0xd77025a1e0a39d8c)},
	{UINT128(0x0000000000000000, 0x368e65de7ace44de)}, {UINT128(0x0000000000000000, 0xfa1e842ffc96e4e1)},
	{UINT128(0xffffffffffffffff, 0x05e17bd003691b1f)}, {UINT128(0x0000000000000000, 0x368e65de7ace44de)},
	{UINT128(0x0000000000000000, 0xee46be5a0813016c)}, {UINT128(0x0000000000000000, 0x5d98d03c9063d92c)},
	{UINT128(0xffffffffffffffff, 0xa2672fc36f9c26d4)}, {UINT128(0x0000000000000000, 0xee46be5a0813016c)},
	{UINT128(0x0000000000000000, 0x664dc58506f7faeb)}, {UINT128(0x0000000000000000, 0xeaab7a9749f584ff)},
	{UINT128(0xffffffffffffffff, 0x15548568b60a7b01)}, {UINT128(0x0000000000000000, 0x664dc58506f7faeb)},
	{UINT128(0x0000000000000000, 0xb8521598bb6bce27)}, {UINT128(0x0000000000000000, 0xb1a81d18e0df488a)},
	{UINT128(0xffffffffffffffff, 0x4e57e2e71f20b776)}, {UINT128(0x0000000000000000, 0xb8521598bb6bce27)},
	{UINT128(0x0000000000000000, 0x04b64daef8c3bf4e)}, {UINT128(0x0000000000000000, 0xfff4e5a25a8d095c)},
	{UINT128(0xffffffffffffffff, 0x000b1a5da572f6a4)}, {UINT128(0x0000000000000000, 0x04b64daef8c3bf4e)},
	{UINT128(0x0000000000000000, 0xfff4e5a25a8d095c)}, {UINT128(0x0000000000000000, 0x04b64daef8c3bf4e)},
	{UINT128(0xffffffffffffffff, 0xfb49b251073c40b2)}, {UINT128(0x0000000000000000, 0xfff4e5a25a8d095c)},
	{UINT128(0x0000000000000000, 0xb1a81d18e0df488a)}, {UINT128(0x0000000000000000, 0xb8521598bb6bce27)},
	{UINT128(0xffffffffffffffff, 0x47adea67449431d9)}, {UINT128(0x0000000000000000, 0xb1a81d18e0df488a)},
	{UINT128(0x0000000000000000, 0xeaab7a9749f584ff)}, {UINT128(0x0000000000000000, 0x664dc58506f7faeb)},
	{UINT128(0xffffffffffffffff, 0x99b23a7af9080515)}, {UINT128(0x0000000000000000, 0xeaab7a9749f584ff)},
	{UINT128(0x0000000000000000, 0x5d98d03c9063d92c)}, {UINT128(0x0000000000000000, 0xee46be5a0813016c)},
	{UINT128(0xffffffffffffffff, 0x11b941a5f7ecfe94)}, {UINT128(0x0000000000000000, 0x5d98d03c9063d92c)},
	{UINT128(0x0000000000000000, 0xfa1e842ffc96e4e1)}, {UINT128(0x0000000000000000, 0x368e65de7ace44de)},
	{UINT128(0xffffffffffffffff, 0xc9719a218531bb22)}, {UINT128(0x0000000000000000, 0xfa1e842ffc96e4e1)},
	{UINT128(0x0000000000000000, 0x8a48ad799b6759f4)}, {UINT128(0x0000000000000000, 0xd77025a1e0a39d8c)},
	{UINT128(0xffffffffffffffff, 0x288fda5e1f5c6274)}, {UINT128(0x0000000000000000, 0x8a48ad799b6759f4)},
	{UINT128(0x0000000000000000, 0xd233c6408cd64237)}, {UINT128(0x0000000000000000, 0x921eafdcc560f9c6)},
	{UINT128(0xffffffffffffffff, 0x6de150233a9f063a)}, {UINT128(0x0000000000000000, 0xd233c6408cd64237)},
	{UINT128(0x0000000000000000, 0x2d502609ec95647d)}, {UINT128(0x0000000000000000, 0xfbf5314f31eb7376)},
	{UINT128(0xffffffffffffffff, 0x040aceb0ce148c8a)}, {UINT128(0x0000000000000000, 0x2d502609ec95647d)},
	{UINT128(0x0000000000000000, 0xfe432367f5b90a63)}, {UINT128(0x0000000000000000, 0x1dc70ecbae9fc913)},
	{UINT128(0xffffffffffffffff, 0xe238f134516036ed)}, {UINT128(0x0000000000000000, 0xfe432367f5b90a63)},
	{UINT128(0x0000000000000000, 0x9ebc11c62c1a1dfc)}, {UINT128(0x0000000000000000, 0xc8d8b37ea4ed0f63)},
	{UINT128(0xffffffffffffffff, 0x37274c815b12f09d)}, {UINT128(0x0000000000000000, 0x9ebc11c62c1a1dfc)},
	{UINT128(0x0000000000000000, 0xdf83270a9bbee891)}, {UINT128(0x0000000000000000, 0x7cd01658ff419d07)},
	{UINT128(0xffffffffffffffff, 0x832fe9a700be62f9)}, {UINT128(0x0000000000000000, 0xdf83270a9bbee891)},
	{UINT128(0x0000000000000000, 0x45ca835dd945dc92)}, {UINT128(0x0000000000000000, 0xf64d969e1dfc211a)},
	{UINT128(0xffffffffffffffff, 0x09b26961e203dee6)}, {UINT128(0x0000000000000000, 0x45ca835dd945dc92)},
	{UINT128(0x0000000000000000, 0xf3913edb54ba2243)}, {UINT128(0x0000000000000000, 0x4ecf3be81052ddf4)},
	{UINT128(0xffffffffffffffff, 0xb130c417efad220c)}, {UINT128(0x0000000000000000, 0xf3913edb54ba2243)},
	{UINT128(0x0000000000000000, 0x74805ba3a76d6b1c)}, {UINT128(0x0000000000000000, 0xe3f4729119e798da)},
	{UINT128(0xffffffffffffffff, 0x1c0b8d6ee6186726)}, {UINT128(0x0000000000000000, 0x74805ba3a76d6b1c)},
	{UINT128(0x0000000000000000, 0xc2de28d74ac6628c)}, {UINT128(0x0000000000000000, 0xa6050a2f6000204a)},
	{UINT128(0xffffffffffffffff, 0x59faf5d09fffdfb6)}, {UINT128(0x0000000000000000, 0xc2de28d74ac6628c)},
	{UINT128(0x0000000000000000, 0x14661179272095ae)}, {UINT128(0x0000000000000000, 0xff2f9d7971ca0365)},
	{UINT128(0xffffffffffffffff, 0x00d062868e35fc9b)}, {UINT128(0x0000000000000000, 0x14661179272095ae)},
	{UINT128(0x0000000000000000, 0xff6ac765b39e1e1a)}, {UINT128(0x0000000000000000, 0x11440134d709b281)},
	{UINT128(0xffffffffffffffff, 0xeebbfecb28f64d7f)}, {UINT128(0x0000000000000000, 0xff6ac765b39e1e1a)},
	{UINT128(0x0000000000000000, 0xa86604facd04d96a)}, {UINT128(0x0000000000000000, 0xc0d0d99dabd65d45)},
	{UINT128(0xffffffffffffffff, 0x3f2f26625429a2bb)}, {UINT128(0x0000000000000000, 0xa86604facd04d96a)},
	{UINT128(0x0000000000000000, 0xe55e0b4d05c80389)}, {UINT128(0x0000000000000000, 0x71b1fd265c002a42)},
	{UINT128(0xffffffffffffffff, 0x8e4e02d9a3ffd5be)}, {UINT128(0x0000000000000000, 0xe55e0b4d05c80389)},
	{UINT128(0x0000000000000000, 0x51cae2955c414eff)}, {UINT128(0x0000000000000000, 0xf294f82394ffe321)},
	{UINT128(0xffffffffffffffff, 0x0d6b07dc6b001cdf)}, {UINT128(0x0000000000000000, 0x51cae2955c414eff)},
	{UINT128(0x0000000000000000, 0xf7241712d4edde4a)}, {UINT128(0x0000000000000000, 0x42c3673f6f6e404f)},
	{UINT128(0xffffffffffffffff, 0xbd3c98c09091bfb1)}, {UINT128(0x0000000000000000, 0xf7241712d4edde4a)},
	{UINT128(0x0000000000000000, 0x7f8bd92f9c483ed7)}, {UINT128(0x0000000000000000, 0xddf6be249c075038)},
	{UINT128(0xffffffffffffffff, 0x220941db63f8afc8)}, {UINT128(0x0000000000000000, 0x7f8bd92f9c483ed7)},
	{UINT128(0x0000000000000000, 0xcac77f24736eb554)}, {UINT128(0x0000000000000000, 0x9c420c2eff590e60)},
	{UINT128(0xffffffffffffffff, 0x63bdf3d100a6f1a0)}, {UINT128(0x0000000000000000, 0xcac77f24736eb554)},
	{UINT128(0x0000000000000000, 0x20e5408f75063a47)}, {UINT128(0x0000000000000000, 0xfde0b0bf220c2fd5)},
	{UINT128(0xffffffffffffffff, 0x021f4f40ddf3d02b)}, {UINT128(0x0000000000000000, 0x20e5408f75063a47)},
	{UINT128(0x0000000000000000, 0xfc7eaffd720ed674)}, {UINT128(0x0000000000000000, 0x2a37bf0b2f8be3f4)},
	{UINT128(0xffffffffffffffff, 0xd5c840f4d0741c0c)}, {UINT128(0x0000000000000000, 0xfc7eaffd720ed674)},
	{UINT128(0x0000000000000000, 0x94b0393b14e54157)}, {UINT128(0x0000000000000000, 0xd064af55d7c9b43f)},
	{UINT128(0xffffffffffffffff, 0x2f9b50aa28364bc1)}, {UINT128(0x0000000000000000, 0x94b0393b14e54157)},
	{UINT128(0x0000000000000000, 0xd91e6a38009da15b)}, {UINT128(0x0000000000000000, 0x87a135d95473ec8a)},
	{UINT128(0xffffffffffffffff, 0x785eca26ab8c1376)}, {UINT128(0x0000000000000000, 0xd91e6a38009da15b)},
	{UINT128(0x0000000000000000, 0x399f196625650c49)}, {UINT128(0x0000000000000000, 0xf96e4e4844d4e82b)},
	{UINT128(0xffffffffffffffff, 0x0691b1b7bb2b17d5)}, {UINT128(0x0000000000000000, 0x399f196625650c49)},
	{UINT128(0x0000000000000000, 0xef682fbef23ecda7)}, {UINT128(0x0000000000000000, 0x5aaa75f71df85ac8)},
	{UINT128(0xffffffffffffffff, 0xa5558a08e207a538)}, {UINT128(0x0000000000000000, 0xef682fbef23ecda7)},
	{UINT128(0x0000000000000000, 0x692d049f7a879242)}, {UINT128(0x0000000000000000, 0xe9659107077cf610)},
	{UINT128(0xffffffffffffffff, 0x169a6ef8f88309f0)}, {UINT128(0x0000000000000000, 0x692d049f7a879242)},
	{UINT128(0x0000000000000000, 0xba7ca46d46946803)}, {UINT128(0x0000000000000000, 0xaf61a4ac1b83a1df)},
	{UINT128(0xffffffffffffffff, 0x509e5b53e47c5e21)}, {UINT128(0x0000000000000000, 0xba7ca46d46946803)},
	{UINT128(0x0000000000000000, 0x07da4dcc7473c040)}, {UINT128(0x0000000000000000, 0xffe128ef8e9fc17b)},
	{UINT128(0xffffffffffffffff, 0x001ed71071603e85)}, {UINT128(0x0000000000000000, 0x07da4dcc7473c040)},
	{UINT128(0x0000000000000000, 0xffc38ed6dc0ef98c)}, {UINT128(0x0000000000000000, 0x0afe00694866a1b5)},
	{UINT128(0xffffffffffffffff, 0xf501ff96b7995e4b)}, {UINT128(0x0000000000000000, 0xffc38ed6dc0ef98c)},
	{UINT128(0x0000000000000000, 0xad146952eb9282b0)}, {UINT128(0x0000000000000000, 0xbca002ba7aaf25eb)},
	{UINT128(0xffffffffffffffff, 0x435ffd458550da15)}, {UINT128(0x0000000000000000, 0xad146952eb9282b0)},
	{UINT128(0x0000000000000000, 0xe816a7f595ec9233)}, {UINT128(0x0000000000000000, 0x6c0835b1fd002451)},
	{UINT128(0xffffffffffffffff, 0x93f7ca4e02ffdbaf)}, {UINT128(0x0000000000000000, 0xe816a7f595ec9233)},
	{UINT128(0x0000000000000000, 0x57b89cde7a9a4d64)}, {UINT128(0x0000000000000000, 0xf08066514c055f7f)},
	{UINT128(0xffffffffffffffff, 0x0f7f99aeb3faa081)}, {UINT128(0x0000000000000000, 0x57b89cde7a9a4d64)},
	{UINT128(0x0000000000000000, 0xf8b47a9fb902e76d)}, {UINT128(0x0000000000000000, 0x3cad943c203436a4)},
	{UINT128(0xffffffffffffffff, 0xc3526bc3dfcbc95c)}, {UINT128(0x0000000000000000, 0xf8b47a9fb902e76d)},
	{UINT128(0x0000000000000000, 0x84f483a0be2f0404)}, {UINT128(0x0000000000000000, 0xdac44ff490a02711)},
	{UINT128(0xffffffffffffffff, 0x253bb00b6f5fd8ef)}, {UINT128(0x0000000000000000, 0x84f483a0be2f0404)},
	{UINT128(0x0000000000000000, 0xce8d8faf5406ab8c)}, {UINT128(0x0000000000000000, 0x973c071f4750b49d)},
	{UINT128(0xffffffffffffffff, 0x68c3f8e0b8af4b63)}, {UINT128(0x0000000000000000, 0xce8d8faf5406ab8c)},
	{UINT128(0x0000000000000000, 0x271db7619b1a2774)}, {UINT128(0x0000000000000000, 0xfcfe72ad6d9641f3)},
	{UINT128(0xffffffffffffffff, 0x03018d529269be0d)}, {UINT128(0x0000000000000000, 0x271db7619b1a2774)},
	{UINT128(0x0000000000000000, 0xfd747472367dd6c6)}, {UINT128(0x0000000000000000, 0x24022da9d8f79d6e)},
	{UINT128(0xffffffffffffffff, 0xdbfdd25627086292)}, {UINT128(0x0000000000000000, 0xfd747472367dd6c6)},
	{UINT128(0x0000000000000000, 0x99c200686472b4a9)}, {UINT128(0x0000000000000000, 0xccae7976c0691178)},
	{UINT128(0xffffffffffffffff, 0x335186893f96ee88)}, {UINT128(0x0000000000000000, 0x99c200686472b4a9)},
	{UINT128(0x0000000000000000, 0xdc61c693a82745d6)}, {UINT128(0x0000000000000000, 0x8242b1357110d373)},
	{UINT128(0xffffffffffffffff, 0x7dbd4eca8eef2c8d)}, {UINT128(0x0000000000000000, 0xdc61c693a82745d6)},
	{UINT128(0x0000000000000000, 0x3fb9b835bfdbf16a)}, {UINT128(0x0000000000000000, 0xf7f110605caf6391)},
	{UINT128(0xffffffffffffffff, 0x080eef9fa3509c6f)}, {UINT128(0x0000000000000000, 0x3fb9b835bfdbf16a)},
	{UINT128(0x0000000000000000, 0xf18f574386712644)}, {UINT128(0x0000000000000000, 0x54c36202bcf08e06)},
	{UINT128(0xffffffffffffffff, 0xab3c9dfd430f71fa)}, {UINT128(0x0000000000000000, 0xf18f574386712644)},
	{UINT128(0x0000000000000000, 0x6edf3c8c12f404d0)}, {UINT128(0x0000000000000000, 0xe6becc4c5997af07)},
	{UINT128(0xffffffffffffffff, 0x194133b3a66850f9)}, {UINT128(0x0000000000000000, 0x6edf3c8c12f404d0)},
	{UINT128(0x0000000000000000, 0xbebc1b6619ed9117)}, {UINT128(0x0000000000000000, 0xaac081c4ba89ba8b)},
	{UINT128(0xffffffffffffffff, 0x553f7e3b45764575)}, {UINT128(0x0000000000000000, 0xbebc1b6619ed9117)},
	{UINT128(0x0000000000000000, 0x0e214689606bf168)}, {UINT128(0x0000000000000000, 0xff9c187c6abade6b)},
	{UINT128(0xffffffffffffffff, 0x0063e78395452195)}, {UINT128(0x0000000000000000, 0x0e214689606bf168)},
	{UINT128(0x0000000000000000, 0xfeea9cff8fa2ae55)}, {UINT128(0x0000000000000000, 0x1787586a5d5b2158)},
	{UINT128(0xffffffffffffffff, 0xe878a795a2a4dea8)}, {UINT128(0x0000000000000000, 0xfeea9cff8fa2ae55)},
	{UINT128(0x0000000000000000, 0xa39da8dcc39a38e6)}, {UINT128(0x0000000000000000, 0xc4e3f4d26ea553b7)},
	{UINT128(0xffffffffffffffff, 0x3b1c0b2d915aac49)}, {UINT128(0x0000000000000000, 0xa39da8dcc39a38e6)},
	{UINT128(0x0000000000000000, 0xe28210095b483752)}, {UINT128(0x0000000000000000, 0x774a3c5207315fc4)},
	{UINT128(0xffffffffffffffff, 0x88b5c3adf8cea03c)}, {UINT128(0x0000000000000000, 0xe28210095b483752)},
	{UINT128(0x0000000000000000, 0x4bd08b6bb00e1a8b)}, {UINT128(0x0000000000000000, 0xf48421b0efbf939c)},
	{UINT128(0xffffffffffffffff, 0x0b7bde4f10406c64)}, {UINT128(0x0000000000000000, 0x4bd08b6bb00e1a8b)},
	{UINT128(0x0000000000000000, 0xf56d97473d446cdb)}, {UINT128(0x0000000000000000, 0x48ceeeaf0eedc59a)},
	{UINT128(0xffffffffffffffff, 0xb7311150f1123a66)}, {UINT128(0x0000000000000000, 0xf56d97473d446cdb)},
	{UINT128(0x0000000000000000, 0x7a0f83abe1444f43)}, {UINT128(0x0000000000000000, 0xe106f1fd4b8d7c97)},
	{UINT128(0xffffffffffffffff, 0x1ef90e02b4728369)}, {UINT128(0x0000000000000000, 0x7a0f83abe1444f43)},
	{UINT128(0x0000000000000000, 0xc6e22998b4c6608f)}, {UINT128(0x0000000000000000, 0xa12ff8bc735d8af7)},
	{UINT128(0xffffffffffffffff, 0x5ed007438ca27509)}, {UINT128(0x0000000000000000, 0xc6e22998b4c6608f)},
	{UINT128(0x0000000000000000, 0x1aa7b724495c037a)}, {UINT128(0x0000000000000000, 0xfe9bc8a1105c22a6)},
	{UINT128(0xffffffffffffffff, 0x0164375eefa3dd5a)}, {UINT128(0x0000000000000000, 0x1aa7b724495c037a)},
	{UINT128(0x0000000000000000, 0xfb61fbefadddb986)}, {UINT128(0x0000000000000000, 0x3066cdd138c98b74)},
	{UINT128(0xffffffffffffffff, 0xcf99322ec736748c)}, {UINT128(0x0000000000000000, 0xfb61fbefadddb986)},
	{UINT128(0x0000000000000000, 0x8f87845de430d778)}, {UINT128(0x0000000000000000, 0xd3fac294ff34e4d1)},
	{UINT128(0xffffffffffffffff, 0x2c053d6b00cb1b2f)}, {UINT128(0x0000000000000000, 0x8f87845de430d778)},
	{UINT128(0x0000000000000000, 0xd5b992c8b606a352)}, {UINT128(0x0000000000000000, 0x8cead04f95cdbf67)},
	{UINT128(0xffffffffffffffff, 0x73152fb06a324099)}, {UINT128(0x0000000000000000, 0xd5b992c8b606a352)},
	{UINT128(0x0000000000000000, 0x337b97e5b886cccc)}, {UINT128(0x0000000000000000, 0xfac5158bc4f42120)},
	{UINT128(0xffffffffffffffff, 0x053aea743b0bdee0)}, {UINT128(0x0000000000000000, 0x337b97e5b886cccc)},
	{UINT128(0x0000000000000000, 0xed1c1d4b344c3d50)}, {UINT128(0x0000000000000000, 0x60838ec13aab0fce)},
	{UINT128(0xffffffffffffffff, 0x9f7c713ec554f032)}, {UINT128(0x0000000000000000, 0xed1c1d4b344c3d50)},
	{UINT128(0x0000000000000000, 0x636a94bb2292bf47)}, {UINT128(0x0000000000000000, 0xebe85815c767cb01)},
	{UINT128(0xffffffffffffffff, 0x1417a7ea389834ff)}, {UINT128(0x0000000000000000, 0x636a94bb2292bf47)},
	{UINT128(0x0000000000000000, 0xb6206b9e0c13a893)}, {UINT128(0x0000000000000000, 0xb3e7bc248d78802f)},
	{UINT128(0xffffffffffffffff, 0x4c1843db72877fd1)}, {UINT128(0x0000000000000000, 0xb6206b9e0c13a893)},
	{UINT128(0x0000000000000000, 0x01921f0fe6700712)}, {UINT128(0x0000000000000000, 0xfffec42c7454926c)},
	{UINT128(0xffffffffffffffff, 0x00013bd38bab6d94)}, {UINT128(0x0000000000000000, 0x01921f0fe6700712)},
	{UINT128(0x0000000000000000, 0xffffb10b10e80e96)}, {UINT128(0x0000000000000000, 0x00c90fc5f66525d3)},
	{UINT128(0xffffffffffffffff, 0xff36f03a099ada2d)}, {UINT128(0x0000000000000000, 0xffffb10b10e80e96)},
	{UINT128(0x0000000000000000, 0xb4768f550ce389fe)}, {UINT128(0x0000000000000000, 0xb592e7697f14dd4b)},
	{UINT128(0xffffffffffffffff, 0x4a6d189680eb22b5)}, {UINT128(0x0000000000000000, 0xb4768f550ce389fe)},
	{UINT128(0x0000000000000000, 0xec3624222d227bd2)}, {UINT128(0x0000000000000000, 0x62b12e1b57b51118)},
	{UINT128(0xffffffffffffffff, 0x9d4ed1e4a84aeee8)}, {UINT128(0x0000000000000000, 0xec3624222d227bd2)},
	{UINT128(0x0000000000000000, 0x613daaabce40fcb7)}, {UINT128(0x0000000000000000, 0xecd006ec59ea3070)},
	{UINT128(0xffffffffffffffff, 0x132ff913a615cf90)}, {UINT128(0x0000000000000000, 0x613daaabce40fcb7)},
	{UINT128(0x0000000000000000, 0xfaed376a1b42e55a)}, {UINT128(0x0000000000000000, 0x32b693d36c55f3de)},
	{UINT128(0xffffffffffffffff, 0xcd496c2c93aa0c22)}, {UINT128(0x0000000000000000, 0xfaed376a1b42e55a)},
	{UINT128(0x0000000000000000, 0x8d9280b89b9df49c)}, {UINT128(0x0000000000000000, 0xd54aa3d165cc7019)},
	{UINT128(0xffffffffffffffff, 0x2ab55c2e9a338fe7)}, {UINT128(0x0000000000000000, 0x8d9280b89b9df49c)},
	{UINT128(0x0000000000000000, 0xd46b3b72a2d68fca)}, {UINT128(0x0000000000000000, 0x8ee0db26e24390f9)},
	{UINT128(0xffffffffffffffff, 0x711f24d91dbc6f07)}, {UINT128(0x0000000000000000, 0xd46b3b72a2d68fca)},
	{UINT128(0x0000000000000000, 0x312c2e4f88300850)}, {UINT128(0x0000000000000000, 0xfb3baab441e770f8)},
	{UINT128(0xffffffffffffffff, 0x04c4554bbe188f08)}, {UINT128(0x0000000000000000, 0x312c2e4f88300850)},
	{UINT128(0x0000000000000000, 0xfeb0696d3ae4f04e)}, {UINT128(0x0000000000000000, 0x19dfb6eb24a85c3e)},
	{UINT128(0xffffffffffffffff, 0xe6204914db57a3c2)}, {UINT128(0x0000000000000000, 0xfeb0696d3ae4f04e)},
	{UINT128(0x0000000000000000, 0xa1cbfad9521bfd1c)}, {UINT128(0x0000000000000000, 0xc66353a8c232a43d)},
	{UINT128(0xffffffffffffffff, 0x399cac573dcd5bc3)}, {UINT128(0x0000000000000000, 0xa1cbfad9521bfd1c)},
	{UINT128(0x0000000000000000, 0xe1668a498f1f892d)}, {UINT128(0x0000000000000000, 0x795ea1b4f360936a)},
	{UINT128(0xffffffffffffffff, 0x86a15e4b0c9f6c96)}, {UINT128(0x0000000000000000, 0xe1668a498f1f892d)},
	{UINT128(0x0000000000000000, 0x498f9a655552e91e)}, {UINT128(0x0000000000000000, 0xf5341c9f32bffeba)},
	{UINT128(0xffffffffffffffff, 0x0acbe360cd400146)}, {UINT128(0x0000000000000000, 0x498f9a655552e91e)},
	{UINT128(0x0000000000000000, 0xf4bf61b00b5982b8)}, {UINT128(0x0000000000000000, 0x4b10693a5514819d)},
	{UINT128(0xffffffffffffffff, 0xb4ef96c5aaeb7e63)}, {UINT128(0x0000000000000000, 0xf4bf61b00b5982b8)},
	{UINT128(0x0000000000000000, 0x77fbfd9aa5077659)}, {UINT128(0x0000000000000000, 0xe224198a0e002124)},
	{UINT128(0xffffffffffffffff, 0x1ddbe675f1ffdedc)}, {UINT128(0x0000000000000000, 0x77fbfd9aa5077659)},
	{UINT128(0x0000000000000000, 0xc56438f6f0ec3ccb)}, {UINT128(0x0000000000000000, 0xa302d34959951244)},
	{UINT128(0xffffffffffffffff, 0x5cfd2cb6a66aedbc)}, {UINT128(0x0000000000000000, 0xc56438f6f0ec3ccb)},
	{UINT128(0x0000000000000000, 0x184f8712c130a087)}, {UINT128(0x0000000000000000, 0xfed7d3a8a29c603c)},
	{UINT128(0xffffffffffffffff, 0x01282c575d639fc4)}, {UINT128(0x0000000000000000, 0x184f8712c130a087)},
	{UINT128(0x0000000000000000, 0xffa6e2a58df6947e)}, {UINT128(0x0000000000000000, 0x0d5880deafc18b54)},
	{UINT128(0xffffffffffffffff, 0xf2a77f21503e74ac)}, {UINT128(0x0000000000000000, 0xffa6e2a58df6947e)},
	{UINT128(0x0000000000000000, 0xab561a8cbb24f411)}, {UINT128(0x0000000000000000, 0xbe35c4e716999631)},
	{UINT128(0xffffffffffffffff, 0x41ca3b18e96669cf)}, {UINT128(0x0000000000000000, 0xab561a8cbb24f411)},
	{UINT128(0x0000000000000000, 0xe715993ccd02fe9d)}, {UINT128(0x0000000000000000, 0x6e29e053f55a4d9b)},
	{UINT128(0xffffffffffffffff, 0x91d61fac0aa5b265)}, {UINT128(0x0000000000000000, 0xe715993ccd02fe9d)},
	{UINT128(0x0000000000000000, 0x5581004bd19ed1a1)}, {UINT128(0x0000000000000000, 0xf14c7a21c8cbd0f5)},
	{UINT128(0xffffffffffffffff, 0x0eb385de37342f0b)}, {UINT128(0x0000000000000000, 0x5581004bd19ed1a1)},
	{UINT128(0x0000000000000000, 0xf822d0a67b9c5cb6)}, {UINT128(0x0000000000000000, 0x3ef6e9017a700333)},
	{UINT128(0xffffffffffffffff, 0xc10916fe858ffccd)}, {UINT128(0x0000000000000000, 0xf822d0a67b9c5cb6)},
	{UINT128(0x0000000000000000, 0x82ef9f618dc5b70f)}, {UINT128(0x0000000000000000, 0xdbfb34373c974b0f)},
	{UINT128(0xffffffffffffffff, 0x2404cbc8c368b4f1)}, {UINT128(0x0000000000000000, 0x82ef9f618dc5b70f)},
	{UINT128(0x0000000000000000, 0xcd26fd2158358e7e)}, {UINT128(0x0000000000000000, 0x99210f624d30facc)},
	{UINT128(0xffffffffffffffff, 0x66def09db2cf0534)}, {UINT128(0x0000000000000000, 0xcd26fd2158358e7e)},
	{UINT128(0x0000000000000000, 0x24c9329bfa6812d4)}, {UINT128(0x0000000000000000, 0xfd57de5867eedc3a)},
	{UINT128(0xffffffffffffffff, 0x02a821a7981123c6)}, {UINT128(0x0000000000000000, 0x24c9329bfa6812d4)},
	{UINT128(0x0000000000000000, 0xfd1cdd63d0bc8736)}, {UINT128(0x0000000000000000, 0x2656f7f28a2d4e96)},
	{UINT128(0xffffffffffffffff, 0xd9a9080d75d2b16a)}, {UINT128(0x0000000000000000, 0xfd1cdd63d0bc8736)},
	{UINT128(0x0000000000000000, 0x97de125a2080a8ee)}, {UINT128(0x0000000000000000, 0xce16888783ae13b4)},
	{UINT128(0xffffffffffffffff, 0x31e977787c51ec4c)}, {UINT128(0x0000000000000000, 0x97de125a2080a8ee)},
	{UINT128(0x0000000000000000, 0xdb2c78a7ede238aa)}, {UINT128(0x0000000000000000, 0x844889019584609a)},
	{UINT128(0xffffffffffffffff, 0x7bb776fe6a7b9f66)}, {UINT128(0x0000000000000000, 0xdb2c78a7ede238aa)},
	{UINT128(0x0000000000000000, 0x3d70d68c5bc68f21)}, {UINT128(0x0000000000000000, 0xf88485e44c7af48b)},
	{UINT128(0xffffffffffffffff, 0x077b7a1bb3850b75)}, {UINT128(0x0000000000000000, 0x3d70d68c5bc68f21)},
	{UINT128(0x0000000000000000, 0xf0c5017ee336ca10)}, {UINT128(0x0000000000000000, 0x56fb9e2e76ced87a)},
	{UINT128(0xffffffffffffffff, 0xa90461d189312786)}, {UINT128(0x0000000000000000, 0xf0c5017ee336ca10)},
	{UINT128(0x0000000000000000, 0x6cbe5c76cc65c9fb)}, {UINT128(0x0000000000000000, 0xe7c187467233d509)},
	{UINT128(0xffffffffffffffff, 0x183e78b98dcc2af7)}, {UINT128(0x0000000000000000, 0x6cbe5c76cc65c9fb)},
	{UINT128(0x0000000000000000, 0xbd27b83dfcbe927a)}, {UINT128(0x0000000000000000, 0xac800eafb91ef9aa)},
	{UINT128(0xffffffffffffffff, 0x537ff15046e10656)}, {UINT128(0x0000000000000000, 0xbd27b83dfcbe927a)},
	{UINT128(0x0000000000000000, 0x0bc6dd52c3a342ec)}, {UINT128(0x0000000000000000, 0xffba9dd8dc8b1e84)},
	{UINT128(0xffffffffffffffff, 0x004562272374e17c)}, {UINT128(0x0000000000000000, 0x0bc6dd52c3a342ec)},
	{UINT128(0x0000000000000000, 0xffe704e71533c509)}, {UINT128(0x0000000000000000, 0x071153d3394ec723)},
	{UINT128(0xffffffffffffffff, 0xf8eeac2cc6b138dd)}, {UINT128(0x0000000000000000, 0xffe704e71533c509)},
	{UINT128(0x0000000000000000, 0xaff3e5ef2b507c07)}, {UINT128(0x0000000000000000, 0xb9f2ac703cca0db4)},
	{UINT128(0xffffffffffffffff, 0x460d538fc335f24c)}, {UINT128(0x0000000000000000, 0xaff3e5ef2b507c07)},
	{UINT128(0x0000000000000000, 0xe9b7e3de5fdedc8c)}, {UINT128(0x0000000000000000, 0x6875950ed42aff7f)},
	{UINT128(0xffffffffffffffff, 0x978a6af12bd50081)}, {UINT128(0x0000000000000000, 0xe9b7e3de5fdedc8c)},
	{UINT128(0x0000000000000000, 0x5b66618e2832d7f8)}, {UINT128(0x0000000000000000, 0xef20b07b6c6c0b38)},
	{UINT128(0xffffffffffffffff, 0x10df4f849393f4c8)}, {UINT128(0x0000000000000000, 0x5b66618e2832d7f8)},
	{UINT128(0x0000000000000000, 0xf99b42d1d57781ec)}, {UINT128(0x0000000000000000, 0x38db20a6bae9b9c9)},
	{UINT128(0xffffffffffffffff, 0xc724df5945164637)}, {UINT128(0x0000000000000000, 0xf99b42d1d57781ec)},
	{UINT128(0x0000000000000000, 0x884b9246854ab50c)}, {UINT128(0x0000000000000000, 0xd8b3a1526517a48c)},
	{UINT128(0xffffffffffffffff, 0x274c5ead9ae85b74)}, {UINT128(0x0000000000000000, 0x884b9246854ab50c)},
	{UINT128(0x0000000000000000, 0xd0d93696053af099)}, {UINT128(0x0000000000000000, 0x940c5f7a69e5ce1d)},
	{UINT128(0xffffffffffffffff, 0x6bf3a085961a31e3)}, {UINT128(0x0000000000000000, 0xd0d93696053af099)},
	{UINT128(0x0000000000000000, 0x2afe010ca997bc65)}, {UINT128(0x0000000000000000, 0xfc5d39be5a5bbc4c)},
	{UINT128(0xffffffffffffffff, 0x03a2c641a5a443b4)}, {UINT128(0x0000000000000000, 0x2afe010ca997bc65)},
	{UINT128(0x0000000000000000, 0xfdfa3878546c3d29)}, {UINT128(0x0000000000000000, 0x201dd15adf70b65a)},
	{UINT128(0xffffffffffffffff, 0xdfe22ea5208f49a6)}, {UINT128(0x0000000000000000, 0xfdfa3878546c3d29)},
	{UINT128(0x0000000000000000, 0x9ce11f1eb18147b2)}, {UINT128(0x0000000000000000, 0xca4c871d625361aa)},
	{UINT128(0xffffffffffffffff, 0x35b378e29dac9e56)}, {UINT128(0x0000000000000000, 0x9ce11f1eb18147b2)},
	{UINT128(0x0000000000000000, 0xde5aa65869193806)}, {UINT128(0x0000000000000000, 0x7edd5d70934b7753)},
	{UINT128(0xffffffffffffffff, 0x8122a28f6cb488ad)}, {UINT128(0x0000000000000000, 0xde5aa65869193806)},
	{UINT128(0x0000000000000000, 0x43856d385d2736a5)}, {UINT128(0x0000000000000000, 0xf6ef5b503c32858a)},
	{UINT128(0xffffffffffffffff, 0x0910a4afc3cd7a76)}, {UINT128(0x0000000000000000, 0x43856d385d2736a5)},
	{UINT128(0x0000000000000000, 0xf2d4eaa8233e997e)}, {UINT128(0x0000000000000000, 0x510c437224dbc18c)},
	{UINT128(0xffffffffffffffff, 0xaef3bc8ddb243e74)}, {UINT128(0x0000000000000000, 0xf2d4eaa8233e997e)},
	{UINT128(0x0000000000000000, 0x7265ff0e194e229e)}, {UINT128(0x0000000000000000, 0xe50478cdce2246bd)},
	{UINT128(0xffffffffffffffff, 0x1afb873231ddb943)}, {UINT128(0x0000000000000000, 0x7265ff0e194e229e)},
	{UINT128(0x0000000000000000, 0xc154e09faa2ff69b)}, {UINT128(0x0000000000000000, 0xa7ce612e65243292)},
	{UINT128(0xffffffffffffffff, 0x58319ed19adbcd6e)}, {UINT128(0x0000000000000000, 0xc154e09faa2ff69b)},
	{UINT128(0x0000000000000000, 0x120c9674ed444c82)}, {UINT128(0x0000000000000000, 0xff5ce92982087868)},
	{UINT128(0xffffffffffffffff, 0x00a316d67df78798)}, {UINT128(0x0000000000000000, 0x120c9674ed444c82)},
	{UINT128(0x0000000000000000, 0xff3f542a416b0135)}, {UINT128(0x0000000000000000, 0x139d9f12c5a29905)},
	{UINT128(0xffffffffffffffff, 0xec6260ed3a5d66fb)}, {UINT128(0x0000000000000000, 0xff3f542a416b0135)},
	{UINT128(0x0000000000000000, 0xa69de36ac4fbfadd)}, {UINT128(0x0000000000000000, 0xc25b888d7c1fcd39)},
	{UINT128(0xffffffffffffffff, 0x3da4777283e032c7)}, {UINT128(0x0000000000000000, 0xa69de36ac4fbfadd)},
	{UINT128(0x0000000000000000, 0xe44fac3814e09857)}, {UINT128(0x0000000000000000, 0x73cd2ebb873482ef)},
	{UINT128(0xffffffffffffffff, 0x8c32d14478cb7d11)}, {UINT128(0x0000000000000000, 0xe44fac3814e09857)},
	{UINT128(0x0000000000000000, 0x4f8e6fa5bb09c582)}, {UINT128(0x0000000000000000, 0xf3530e2aea9d3967)},
	{UINT128(0xffffffffffffffff, 0x0cacf1d51562c699)}, {UINT128(0x0000000000000000, 0x4f8e6fa5bb09c582)},
	{UINT128(0x0000000000000000, 0xf6841af4cc8048a5)}, {UINT128(0x0000000000000000, 0x4508fbbf1a4de123)},
	{UINT128(0xffffffffffffffff, 0xbaf70440e5b21edd)}, {UINT128(0x0000000000000000, 0xf6841af4cc8048a5)},
	{UINT128(0x0000000000000000, 0x7d7f7b995b368ba2)}, {UINT128(0x0000000000000000, 0xdf20db088aa60405)},
	{UINT128(0xffffffffffffffff, 0x20df24f77559fbfb)}, {UINT128(0x0000000000000000, 0x7d7f7b995b368ba2)},
	{UINT128(0x0000000000000000, 0xc95520fe2d40a74c)}, {UINT128(0x0000000000000000, 0x9e1e224c0e28bc95)},
	{UINT128(0xffffffffffffffff, 0x61e1ddb3f1d7436b)}, {UINT128(0x0000000000000000, 0xc95520fe2d40a74c)},
	{UINT128(0x0000000000000000, 0x1e8eb7fde4aa3eb8)}, {UINT128(0x0000000000000000, 0xfe2b71dbecd7aefd)},
	{UINT128(0xffffffffffffffff, 0x01d48e2413285103)}, {UINT128(0x0000000000000000, 0x1e8eb7fde4aa3eb8)},
	{UINT128(0x0000000000000000, 0xfc187a52063060c3)}, {UINT128(0x0000000000000000, 0x2c8a35063b061a6c)},
	{UINT128(0xffffffffffffffff, 0xd375caf9c4f9e594)}, {UINT128(0x0000000000000000, 0xfc187a52063060c3)},
	{UINT128(0x0000000000000000, 0x92c39a65db88809e)}, {UINT128(0x0000000000000000, 0xd1c0c252c9de2c87)},
	{UINT128(0xffffffffffffffff, 0x2e3f3dad3621d379)}, {UINT128(0x0000000000000000, 0x92c39a65db88809e)},
	{UINT128(0x0000000000000000, 0xd7dc7ec4fabdb012)}, {UINT128(0x0000000000000000, 0x899f4e7f712a765f)},
	{UINT128(0xffffffffffffffff, 0x7660b1808ed589a1)}, {UINT128(0x0000000000000000, 0xd7dc7ec4fabdb012)},
	{UINT128(0x0000000000000000, 0x3752c669e2bb5946)}, {UINT128(0x0000000000000000, 0xf9f35de0dde328ac)},
	{UINT128(0xffffffffffffffff, 0x060ca21f221cd754)}, {UINT128(0x0000000000000000, 0x3752c669e2bb5946)},
	{UINT128(0x0000000000000000, 0xee8ff79c548acd10)}, {UINT128(0x0000000000000000, 0x5cdd8f2498424757)},
	{UINT128(0xffffffffffffffff, 0xa32270db67bdb8a9)}, {UINT128(0x0000000000000000, 0xee8ff79c548acd10)},
	{UINT128(0x0000000000000000, 0x6705f51037e7c8ca)}, {UINT128(0x0000000000000000, 0xea5ad8d8c3a91f06)},
	{UINT128(0xffffffffffffffff, 0x15a527273c56e0fa)}, {UINT128(0x0000000000000000, 0x6705f51037e7c8ca)},
	{UINT128(0x0000000000000000, 0xb8dd64b0720df648)}, {UINT128(0x0000000000000000, 0xb117227f6117f9fa)},
	{UINT128(0xffffffffffffffff, 0x4ee8dd809ee80606)}, {UINT128(0x0000000000000000, 0xb8dd64b0720df648)},
	{UINT128(0x0000000000000000, 0x057f53487eac8978)}, {UINT128(0x0000000000000000, 0xfff0e343865bbb14)},
	{UINT128(0xffffffffffffffff, 0x000f1cbc79a444ec)}, {UINT128(0x0000000000000000, 0x057f53487eac8978)},
	{UINT128(0x0000000000000000, 0xfff84a1e29de8572)}, {UINT128(0x0000000000000000, 0x03ed452d5732f951)},
	{UINT128(0xffffffffffffffff, 0xfc12bad2a8cd06af)}, {UINT128(0x0000000000000000, 0xfff84a1e29de8572)},
	{UINT128(0x0000000000000000, 0xb238aa1bfa9ad508)}, {UINT128(0x0000000000000000, 0xb7c654ce4adba9f3)},
	{UINT128(0xffffffffffffffff, 0x4839ab31b524560d)}, {UINT128(0x0000000000000000, 0xb238aa1bfa9ad508)},
	{UINT128(0x0000000000000000, 0xeafb8b944453f530)}, {UINT128(0x0000000000000000, 0x659556deae523d7b)},
	{UINT128(0xffffffffffffffff, 0x9a6aa92151adc285)}, {UINT128(0x0000000000000000, 0xeafb8b944453f530)},
	{UINT128(0x0000000000000000, 0x5e53d7984f7eb8c2)}, {UINT128(0x0000000000000000, 0xedfcf21cabacd3b2)},
	{UINT128(0xffffffffffffffff, 0x12030de354532c4e)}, {UINT128(0x0000000000000000, 0x5e53d7984f7eb8c2)},
	{UINT128(0x0000000000000000, 0xfa491035e55da3a4)}, {UINT128(0x0000000000000000, 0x35c9e3abe77356d6)},
	{UINT128(0xffffffffffffffff, 0xca361c54188ca92a)}, {UINT128(0x0000000000000000, 0xfa491035e55da3a4)},
	{UINT128(0x0000000000000000, 0x8af1b726df15e13d)}, {UINT128(0x0000000000000000, 0xd703479a2f6776cd)},
	{UINT128(0xffffffffffffffff, 0x28fcb865d0988933)}, {UINT128(0x0000000000000000, 0x8af1b726df15e13d)},
	{UINT128(0x0000000000000000, 0xd2a6488487a91919)}, {UINT128(0x0000000000000000, 0x91796b31609f0c55)},
	{UINT128(0xffffffffffffffff, 0x6e8694ce9f60f3ab)}, {UINT128(0x0000000000000000, 0xd2a6488487a91919)},
	{UINT128(0x0000000000000000, 0x2e15fb1a1189fe94)}, {UINT128(0x0000000000000000, 0xfbd14ce0d191516f)},
	{UINT128(0xffffffffffffffff, 0x042eb31f2e6eae91)}, {UINT128(0x0000000000000000, 0x2e15fb1a1189fe94)},
	{UINT128(0x0000000000000000, 0xfe5a381c8a1aa225)}, {UINT128(0x0000000000000000, 0x1cff533b307dc132)},
	{UINT128(0xffffffffffffffff, 0xe300acc4cf823ece)}, {UINT128(0x0000000000000000, 0xfe5a381c8a1aa225)},
	{UINT128(0x0000000000000000, 0x9f599f55f0340062)}, {UINT128(0x0000000000000000, 0xc85bca1abaf7f0a8)},
	{UINT128(0xffffffffffffffff, 0x37a435e545080f58)}, {UINT128(0x0000000000000000, 0x9f599f55f0340062)},
	{UINT128(0x0000000000000000, 0xdfe4e92d0d8a37f6)}, {UINT128(0x0000000000000000, 0x7c20641affdff6b4)},
	{UINT128(0xffffffffffffffff, 0x83df9be50020094c)}, {UINT128(0x0000000000000000, 0xdfe4e92d0d8a37f6)},
	{UINT128(0x0000000000000000, 0x468bdfefa3c90d65)}, {UINT128(0x0000000000000000, 0xf6167a58d7b59027)},
	{UINT128(0xffffffffffffffff, 0x09e985a7284a6fd9)}, {UINT128(0x0000000000000000, 0x468bdfefa3c90d65)},
	{UINT128(0x0000000000000000, 0xf3ced94d28b2ce8b)}, {UINT128(0x0000000000000000, 0x4e0fd78d4edaaa58)},
	{UINT128(0xffffffffffffffff, 0xb1f02872b12555a8)}, {UINT128(0x0000000000000000, 0xf3ced94d28b2ce8b)},
	{UINT128(0x0000000000000000, 0x753340aea1827365)}, {UINT128(0x0000000000000000, 0xe398ac4cf556b733)},
	{UINT128(0xffffffffffffffff, 0x1c6753b30aa948cd)}, {UINT128(0x0000000000000000, 0x753340aea1827365)},
	{UINT128(0x0000000000000000, 0xc36050ecd50ca831)}, {UINT128(0x0000000000000000, 0xa56bca8b391785dc)},
	{UINT128(0xffffffffffffffff, 0x5a943574c6e87a24)}, {UINT128(0x0000000000000000, 0xc36050ecd50ca831)},
	{UINT128(0x0000000000000000, 0x152e774a4d4d0a1c)}, {UINT128(0x0000000000000000, 0xff1f495f4ec430d8)},
	{UINT128(0xffffffffffffffff, 0x00e0b6a0b13bcf28)}, {UINT128(0x0000000000000000, 0x152e774a4d4d0a1c)},
	{UINT128(0x0000000000000000, 0xff780814130c893d)}, {UINT128(0x0000000000000000, 0x107b614e463063b2)},
	{UINT128(0xffffffffffffffff, 0xef849eb1b9cf9c4e)}, {UINT128(0x0000000000000000, 0xff780814130c893d)},
	{UINT128(0x0000000000000000, 0xa8fd40e6ccd52ffe)}, {UINT128(0x0000000000000000, 0xc04c5bab7297d323)},
	{UINT128(0xffffffffffffffff, 0x3fb3a4548d682cdd)}, {UINT128(0x0000000000000000, 0xa8fd40e6ccd52ffe)},
	{UINT128(0x0000000000000000, 0xe5b7105006d4c561)}, {UINT128(0x0000000000000000, 0x70fdb51c98c4a59d)},
	{UINT128(0xffffffffffffffff, 0x8f024ae3673b5a63)}, {UINT128(0x0000000000000000, 0xe5b7105006d4c561)},
	{UINT128(0x0000000000000000, 0x52894f446e0bcbb6)}, {UINT128(0x0000000000000000, 0xf2546ffc0e72f287)},
	{UINT128(0xffffffffffffffff, 0x0dab9003f18d0d79)}, {UINT128(0x0000000000000000, 0x52894f446e0bcbb6)},
	{UINT128(0x0000000000000000, 0xf7583a62852a23b3)}, {UINT128(0x0000000000000000, 0x42013817ad987955)},
	{UINT128(0xffffffffffffffff, 0xbdfec7e8526786ab)}, {UINT128(0x0000000000000000, 0xf7583a62852a23b3)},
	{UINT128(0x0000000000000000, 0x803a06415c170526)}, {UINT128(0x0000000000000000, 0xdd924d05b620678b)},
	{UINT128(0xffffffffffffffff, 0x226db2fa49df9875)}, {UINT128(0x0000000000000000, 0x803a06415c170526)},
	{UINT128(0x0000000000000000, 0xcb41fa15ebff0778)}, {UINT128(0x0000000000000000, 0x9ba298dc0bfc6a89)},
	{UINT128(0xffffffffffffffff, 0x645d6723f4039577)}, {UINT128(0x0000000000000000, 0xcb41fa15ebff0778)},
	{UINT128(0x0000000000000000, 0x21ac9b7964cf0ba4)}, {UINT128(0x0000000000000000, 0xfdc68c6b356db630)},
	{UINT128(0xffffffffffffffff, 0x02397394ca9249d0)}, {UINT128(0x0000000000000000, 0x21ac9b7964cf0ba4)},
	{UINT128(0x0000000000000000, 0xfc9f8a7c2d603c61)}, {UINT128(0x0000000000000000, 0x297162fef3f510db)},
	{UINT128(0xffffffffffffffff, 0xd68e9d010c0aef25)}, {UINT128(0x0000000000000000, 0xfc9f8a7c2d603c61)},
	{UINT128(0x0000000000000000, 0x9553b743d75ac040)}, {UINT128(0x0000000000000000, 0xcfefa7898a4ef23d)},
	{UINT128(0xffffffffffffffff, 0x3010587675b10dc3)}, {UINT128(0x0000000000000000, 0x9553b743d75ac040)},
	{UINT128(0x0000000000000000, 0xd988ad2f9bdf9bbc)}, {UINT128(0x0000000000000000, 0x86f685c25e25acf6)},
	{UINT128(0xffffffffffffffff, 0x79097a3da1da530a)}, {UINT128(0x0000000000000000, 0xd988ad2f9bdf9bbc)},
	{UINT128(0x0000000000000000, 0x3a62ee9a597bdc98)}, {UINT128(0x0000000000000000, 0xf940bfe2304e6c46)},
	{UINT128(0xffffffffffffffff, 0x06bf401dcfb193ba)}, {UINT128(0x0000000000000000, 0x3a62ee9a597bdc98)},
	{UINT128(0x0000000000000000, 0xefaf1b54dd2cdf10)}, {UINT128(0x0000000000000000, 0x59ee5272b58f2a5e)},
	{UINT128(0xffffffffffffffff, 0xa611ad8d4a70d5a2)}, {UINT128(0x0000000000000000, 0xefaf1b54dd2cdf10)},
	{UINT128(0x0000000000000000, 0x69e4334f6fcc6b41)}, {UINT128(0x0000000000000000, 0xe912ae372d27045e)},
	{UINT128(0xffffffffffffffff, 0x16ed51c8d2d8fba2)}, {UINT128(0x0000000000000000, 0x69e4334f6fcc6b41)},
	{UINT128(0x0000000000000000, 0xbb062961823b1ddd)}, {UINT128(0x0000000000000000, 0xaecef739f1a2df11)},
	{UINT128(0xffffffffffffffff, 0x513108c60e5d20ef)}, {UINT128(0x0000000000000000, 0xbb062961823b1ddd)},
	{UINT128(0x0000000000000000, 0x08a342eda160bf5b)}, {UINT128(0x0000000000000000, 0xffdaaf212fed72dc)},
	{UINT128(0xffffffffffffffff, 0x002550ded0128d24)}, {UINT128(0x0000000000000000, 0x08a342eda160bf5b)},
	{UINT128(0x0000000000000000, 0xffcbe2104600a0aa)}, {UINT128(0x0000000000000000, 0x0a351cb7fc30bc89)},
	{UINT128(0xffffffffffffffff, 0xf5cae34803cf4377)}, {UINT128(0x0000000000000000, 0xffcbe2104600a0aa)},
	{UINT128(0x0000000000000000, 0xada859327ba24152)}, {UINT128(0x0000000000000000, 0xbc17d8dc859ad584)},
	{UINT128(0xffffffffffffffff, 0x43e827237a652a7c)}, {UINT128(0x0000000000000000, 0xada859327ba24152)},
	{UINT128(0x0000000000000000, 0xe86b397ace95c470)}, {UINT128(0x0000000000000000, 0x6b51cc49737023b3)},
	{UINT128(0xffffffffffffffff, 0x94ae33b68c8fdc4d)}, {UINT128(0x0000000000000000, 0xe86b397ace95c470)},
	{UINT128(0x0000000000000000, 0x58756572230809f7)}, {UINT128(0x0000000000000000, 0xf03b36c9407aa3e9)},
	{UINT128(0xffffffffffffffff, 0x0fc4c936bf855c17)}, {UINT128(0x0000000000000000, 0x58756572230809f7)},
	{UINT128(0x0000000000000000, 0xf8e3d5f1423842a1)}, {UINT128(0x0000000000000000, 0x3bea2c7e02133380)},
	{UINT128(0xffffffffffffffff, 0xc415d381fdeccc80)}, {UINT128(0x0000000000000000, 0xf8e3d5f1423842a1)},
	{UINT128(0x0000000000000000, 0x85a02c3c7c2f5ca6)}, {UINT128(0x0000000000000000, 0xda5ba04ef3c929f5)},
	{UINT128(0xffffffffffffffff, 0x25a45fb10c36d60b)}, {UINT128(0x0000000000000000, 0x85a02c3c7c2f5ca6)},
	{UINT128(0x0000000000000000, 0xcf04176da12390ad)}, {UINT128(0x0000000000000000, 0x96999e9a74ddbde4)},
	{UINT128(0xffffffffffffffff, 0x696661658b22421c)}, {UINT128(0x0000000000000000, 0xcf04176da12390ad)},
	{UINT128(0x0000000000000000, 0x27e45eafb6912638)}, {UINT128(0x0000000000000000, 0xfcdf6be7def1464d)},
	{UINT128(0xffffffffffffffff, 0x03209418210eb9b3)}, {UINT128(0x0000000000000000, 0x27e45eafb6912638)},
	{UINT128(0x0000000000000000, 0xfd906e340eaa6402)}, {UINT128(0x0000000000000000, 0x233b12817c49ce7b)},
	{UINT128(0xffffffffffffffff, 0xdcc4ed7e83b63185)}, {UINT128(0x0000000000000000, 0xfd906e340eaa6402)},
	{UINT128(0x0000000000000000, 0x9a6292960a6f0ab1)}, {UINT128(0x0000000000000000, 0xcc35778a2bac9ca2)},
	{UINT128(0xffffffffffffffff, 0x33ca8875d453635e)}, {UINT128(0x0000000000000000, 0x9a6292960a6f0ab1)},
	{UINT128(0x0000000000000000, 0xdcc7d0fec8aaf2ab)}, {UINT128(0x0000000000000000, 0x819572af6decac85)},
	{UINT128(0xffffffffffffffff, 0x7e6a8d509213537b)}, {UINT128(0x0000000000000000, 0xdcc7d0fec8aaf2ab)},
	{UINT128(0x0000000000000000, 0x407c601ae7f746bc)}, {UINT128(0x0000000000000000, 0xf7beb728e51dfcb9)},
	{UINT128(0xffffffffffffffff, 0x084148d71ae20347)}, {UINT128(0x0000000000000000, 0x407c601ae7f746bc)},
	{UINT128(0x0000000000000000, 0xf1d19f63ae7428a3)}, {UINT128(0x0000000000000000, 0x54058f7065c11e14)},
	{UINT128(0xffffffffffffffff, 0xabfa708f9a3ee1ec)}, {UINT128(0x0000000000000000, 0xf1d19f63ae7428a3)},
	{UINT128(0x0000000000000000, 0x6f94545fff03652c)}, {UINT128(0x0000000000000000, 0xe667710616f4fc5a)},
	{UINT128(0xffffffffffffffff, 0x19988ef9e90b03a6)}, {UINT128(0x0000000000000000, 0x6f94545fff03652c)},
	{UINT128(0x0000000000000000, 0xbf41fc3d81b430dc)}, {UINT128(0x0000000000000000, 0xaa2a7fa8acefdd64)},
	{UINT128(0xffffffffffffffff, 0x55d580575310229c)}, {UINT128(0x0000000000000000, 0xbf41fc3d81b430dc)},
	{UINT128(0x0000000000000000, 0x0eea037cc0476485)}, {UINT128(0x0000000000000000, 0xff90b0a7098f6444)},
	{UINT128(0xffffffffffffffff, 0x006f4f58f6709bbc)}, {UINT128(0x0000000000000000, 0x0eea037cc0476485)},
	{UINT128(0x0000000000000000, 0xfefcc917b99839a6)}, {UINT128(0x0000000000000000, 0x16bf1b3e79b128bb)},
	{UINT128(0xffffffffffffffff, 0xe940e4c1864ed745)}, {UINT128(0x0000000000000000, 0xfefcc917b99839a6)},
	{UINT128(0x0000000000000000, 0xa4381983048ff748)}, {UINT128(0x0000000000000000, 0xc463373a40dd06a4)},
	{UINT128(0xffffffffffffffff, 0x3b9cc8c5bf22f95c)}, {UINT128(0x0000000000000000, 0xa4381983048ff748)},
	{UINT128(0x0000000000000000, 0xe2df7acff7c2cf84)}, {UINT128(0x0000000000000000, 0x76983173e8436379)},
	{UINT128(0xffffffffffffffff, 0x8967ce8c17bc9c87)}, {UINT128(0x0000000000000000, 0xe2df7acff7c2cf84)},
	{UINT128(0x0000000000000000, 0x4c907ed8e2eabc1f)}, {UINT128(0x0000000000000000, 0xf4484add6b01254c)},
	{UINT128(0xffffffffffffffff, 0x0bb7b52294fedab4)}, {UINT128(0x0000000000000000, 0x4c907ed8e2eabc1f)},
	{UINT128(0x0000000000000000, 0xf5a67a8adc4088cb)}, {UINT128(0x0000000000000000, 0x480e160f5c9ef71d)},
	{UINT128(0xffffffffffffffff, 0xb7f1e9f0a36108e3)}, {UINT128(0x0000000000000000, 0xf5a67a8adc4088cb)},
	{UINT128(0x0000000000000000, 0x7ac01a57c9588154)}, {UINT128(0x0000000000000000, 0xe0a6cee232f2bb9d)},
	{UINT128(0xffffffffffffffff, 0x1f59311dcd0d4463)}, {UINT128(0x0000000000000000, 0x7ac01a57c9588154)},
	{UINT128(0x0000000000000000, 0xc76084da43624635)}, {UINT128(0x0000000000000000, 0xa0939331e8846238)},
	{UINT128(0xffffffffffffffff, 0x5f6c6cce177b9dc8)}, {UINT128(0x0000000000000000, 0xc76084da43624635)},
	{UINT128(0x0000000000000000, 0x1b6fa6ec38f64c77)}, {UINT128(0x0000000000000000, 0xfe868ac6c3043b2f)},
	{UINT128(0xffffffffffffffff, 0x017975393cfbc4d1)}, {UINT128(0x0000000000000000, 0x1b6fa6ec38f64c77)},
	{UINT128(0x0000000000000000, 0xfb87b21a5bf5b918)}, {UINT128(0x0000000000000000, 0x2fa14f77a5963719)},
	{UINT128(0xffffffffffffffff, 0xd05eb0885a69c8e7)}, {UINT128(0x0000000000000000, 0xfb87b21a5bf5b918)},
	{UINT128(0x0000000000000000, 0x902dd50bab06b1b8)}, {UINT128(0x0000000000000000, 0xd389c6f4eb07a41d)},
	{UINT128(0xffffffffffffffff, 0x2c76390b14f85be3)}, {UINT128(0x0000000000000000, 0x902dd50bab06b1b8)},
	{UINT128(0x0000000000000000, 0xd627fde9f7d63e7f)}, {UINT128(0x0000000000000000, 0x8c42c8f9d2372645)},
	{UINT128(0xffffffffffffffff, 0x73bd37062dc8d9bb)}, {UINT128(0x0000000000000000, 0xd627fde9f7d63e7f)},
	{UINT128(0x0000000000000000, 0x34407c363063b471)}, {UINT128(0x0000000000000000, 0xfa9c58fd796837d5)},
	{UINT128(0xffffffffffffffff, 0x0563a7028697c82b)}, {UINT128(0x0000000000000000, 0x34407c363063b471)},
	{UINT128(0x0000000000000000, 0xed67a1673455c602)}, {UINT128(0x0000000000000000, 0x5fc9374dcd079094)},
	{UINT128(0xffffffffffffffff, 0xa036c8b232f86f6c)}, {UINT128(0x0000000000000000, 0xed67a1673455c602)},
	{UINT128(0x0000000000000000, 0x6423be07bdef454d)}, {UINT128(0x0000000000000000, 0xeb99fa84606ff600)},
	{UINT128(0xffffffffffffffff, 0x1466057b9f900a00)}, {UINT128(0x0000000000000000, 0x6423be07bdef454d)},
	{UINT128(0x0000000000000000, 0xb6ad7f7a557e64f3)}, {UINT128(0x0000000000000000, 0xb35879fa959c323d)},
	{UINT128(0xffffffffffffffff, 0x4ca786056a63cdc3)}, {UINT128(0x0000000000000000, 0xb6ad7f7a557e64f3)},
	{UINT128(0x0000000000000000, 0x025b2d61ca12e05e)}, {UINT128(0x0000000000000000, 0xfffd3964bc6275bb)},
	{UINT128(0xffffffffffffffff, 0x0002c69b439d8a45)}, {UINT128(0x0000000000000000, 0x025b2d61ca12e05e)},
	{UINT128(0x0000000000000000, 0xfffd3964bc6275bb)}, {UINT128(0x0000000000000000, 0x025b2d61ca12e05e)},
	{UINT128(0xffffffffffffffff, 0xfda4d29e35ed1fa2)}, {UINT128(0x0000000000000000, 0xfffd3964bc6275bb)},
	{UINT128(0x0000000000000000, 0xb35879fa959c323d)}, {UINT128(0x0000000000000000, 0xb6ad7f7a557e64f3)},
	{UINT128(0xffffffffffffffff, 0x49528085aa819b0d)}, {UINT128(0x0000000000000000, 0xb35879fa959c323d)},
	{UINT128(0x0000000000000000, 0xeb99fa84606ff600)}, {UINT128(0x0000000000000000, 0x6423be07bdef454d)},
	{UINT128(0xffffffffffffffff, 0x9bdc41f84210bab3)}, {UINT128(0x0000000000000000, 0xeb99fa84606ff600)},
	{UINT128(0x0000000000000000, 0x5fc9374dcd079094)}, {UINT128(0x0000000000000000, 0xed67a1673455c602)},
	{UINT128(0xffffffffffffffff, 0x12985e98cbaa39fe)}, {UINT128(0x0000000000000000, 0x5fc9374dcd079094)},
	{UINT128(0x0000000000000000, 0xfa9c58fd796837d5)}, {UINT128(0x0000000000000000, 0x34407c363063b471)},
	{UINT128(0xffffffffffffffff, 0xcbbf83c9cf9c4b8f)}, {UINT128(0x0000000000000000, 0xfa9c58fd796837d5)},
	{UINT128(0x0000000000000000, 0x8c42c8f9d2372645)}, {UINT128(0x0000000000000000, 0xd627fde9f7d63e7f)},
	{UINT128(0xffffffffffffffff, 0x29d802160829c181)}, {UINT128(0x0000000000000000, 0x8c42c8f9d2372645)},
	{UINT128(0x0000000000000000, 0xd389c6f4eb07a41d)}, {UINT128(0x0000000000000000, 0x902dd50bab06b1b8)},
	{UINT128(0xffffffffffffffff, 0x6fd22af454f94e48)}, {UINT128(0x0000000000000000, 0xd389c6f4eb07a41d)},
	{UINT128(0x0000000000000000, 0x2fa14f77a5963719)}, {UINT128(0x0000000000000000, 0xfb87b21a5bf5b918)},
	{UINT128(0xffffffffffffffff, 0x04784de5a40a46e8)}, {UINT128(0x0000000000000000, 0x2fa14f77a5963719)},
	{UINT128(0x0000000000000000, 0xfe868ac6c3043b2f)}, {UINT128(0x0000000000000000, 0x1b6fa6ec38f64c77)},
	{UINT128(0xffffffffffffffff, 0xe4905913c709b389)}, {UINT128(0x0000000000000000, 0xfe868ac6c3043b2f)},
	{UINT128(0x0000000000000000, 0xa0939331e8846238)}, {UINT128(0x0000000000000000, 0xc76084da43624635)},
	{UINT128(0xffffffffffffffff, 0x389f7b25bc9db9cb)}, {UINT128(0x0000000000000000, 0xa0939331e8846238)},
	{UINT128(0x0000000000000000, 0xe0a6cee232f2bb9d)}, {UINT128(0x0000000000000000, 0x7ac01a57c9588154)},
	{UINT128(0xffffffffffffffff, 0x853fe5a836a77eac)}, {UINT128(0x0000000000000000, 0xe0a6cee232f2bb9d)},
	{UINT128(0x0000000000000000, 0x480e160f5c9ef71d)}, {UINT128(0x0000000000000000, 0xf5a67a8adc4088cb)},
	{UINT128(0xffffffffffffffff, 0x0a59857523bf7735)}, {UINT128(0x0000000000000000, 0x480e160f5c9ef71d)},
	{UINT128(0x0000000000000000, 0xf4484add6b01254c)}, {UINT128(0x0000000000000000, 0x4c907ed8e2eabc1f)},
	{UINT128(0xffffffffffffffff, 0xb36f81271d1543e1)}, {UINT128(0x0000000000000000, 0xf4484add6b01254c)},
	{UINT128(0x0000000000000000, 0x76983173e8436379)}, {UINT128(0x0000000000000000, 0xe2df7acff7c2cf84)},
	{UINT128(0xffffffffffffffff, 0x1d208530083d307c)}, {UINT128(0x0000000000000000, 0x76983173e8436379)},
	{UINT128(0x0000000000000000, 0xc463373a40dd06a4)}, {UINT128(0x0000000000000000, 0xa4381983048ff748)},
	{UINT128(0xffffffffffffffff, 0x5bc7e67cfb7008b8)}, {UINT128(0x0000000000000000, 0xc463373a40dd06a4)},
	{UINT128(0x0000000000000000, 0x16bf1b3e79b128bb)}, {UINT128(0x0000000000000000, 0xfefcc917b99839a6)},
	{UINT128(0xffffffffffffffff, 0x010336e84667c65a)}, {UINT128(0x0000000000000000, 0x16bf1b3e79b128bb)},
	{UINT128(0x0000000000000000, 0xff90b0a7098f6444)}, {UINT128(0x0000000000000000, 0x0eea037cc0476485)},
	{UINT128(0xffffffffffffffff, 0xf115fc833fb89b7b)}, {UINT128(0x0000000000000000, 0xff90b0a7098f6444)},
	{UINT128(0x0000000000000000, 0xaa2a7fa8acefdd64)}, {UINT128(0x0000000000000000, 0xbf41fc3d81b430dc)},
	{UINT128(0xffffffffffffffff, 0x40be03c27e4bcf24)}, {UINT128(0x0000000000000000, 0xaa2a7fa8acefdd64)},
	{UINT128(0x0000000000000000, 0xe667710616f4fc5a)}, {UINT128(0x0000000000000000, 0x6f94545fff03652c)},
	{UINT128(0xffffffffffffffff, 0x906baba000fc9ad4)}, {UINT128(0x0000000000000000, 0xe667710616f4fc5a)},
	{UINT128(0x0000000000000000, 0x54058f7065c11e14)}, {UINT128(0x0000000000000000, 0xf1d19f63ae7428a3)},
	{UINT128(0xffffffffffffffff, 0x0e2e609c518bd75d)}, {UINT128(0x0000000000000000, 0x54058f7065c11e14)},
	{UINT128(0x0000000000000000, 0xf7beb728e51dfcb9)}, {UINT128(0x0000000000000000, 0x407c601ae7f746bc)},
	{UINT128(0xffffffffffffffff, 0xbf839fe51808b944)}, {UINT128(0x0000000000000000, 0xf7beb728e51dfcb9)},
	{UINT128(0x0000000000000000, 0x819572af6decac85)}, {UINT128(0x0000000000000000, 0xdcc7d0fec8aaf2ab)},
	{UINT128(0xffffffffffffffff, 0x23382f0137550d55)}, {UINT128(0x0000000000000000, 0x819572af6decac85)},
	{UINT128(0x0000000000000000, 0xcc35778a2bac9ca2)}, {UINT128(0x0000000000000000, 0x9a6292960a6f0ab1)},
	{UINT128(0xffffffffffffffff, 0x659d6d69f590f54f)}, {UINT128(0x0000000000000000, 0xcc35778a2bac9ca2)},
	{UINT128(0x0000000000000000, 0x233b12817c49ce7b)}, {UINT128(0x0000000000000000, 0xfd906e340eaa6402)},
	{UINT128(0xffffffffffffffff, 0x026f91cbf1559bfe)}, {UINT128(0x0000000000000000, 0x233b12817c49ce7b)},
	{UINT128(0x0000000000000000, 0xfcdf6be7def1464d)}, {UINT128(0x0000000000000000, 0x27e45eafb6912638)},
	{UINT128(0xffffffffffffffff, 0xd81ba150496ed9c8)}, {UINT128(0x0000000000000000, 0xfcdf6be7def1464d)},
	{UINT128(0x0000000000000000, 0x96999e9a74ddbde4)}, {UINT128(0x0000000000000000, 0xcf04176da12390ad)},
	{UINT128(0xffffffffffffffff, 0x30fbe8925edc6f53)}, {UINT128(0x0000000000000000, 0x96999e9a74ddbde4)},
	{UINT128(0x0000000000000000, 0xda5ba04ef3c929f5)}, {UINT128(0x0000000000000000, 0x85a02c3c7c2f5ca6)},
	{UINT128(0xffffffffffffffff, 0x7a5fd3c383d0a35a)}, {UINT128(0x0000000000000000, 0xda5ba04ef3c929f5)},
	{UINT128(0x0000000000000000, 0x3bea2c7e02133380)}, {UINT128(0x0000000000000000, 0xf8e3d5f1423842a1)},
	{UINT128(0xffffffffffffffff, 0x071c2a0ebdc7bd5f)}, {UINT128(0x0000000000000000, 0x3bea2c7e02133380)},
	{UINT128(0x0000000000000000, 0xf03b36c9407aa3e9)}, {UINT128(0x0000000000000000, 0x58756572230809f7)},
	{UINT128(0xffffffffffffffff, 0xa78a9a8ddcf7f609)}, {UINT128(0x0000000000000000, 0xf03b36c9407aa3e9)},
	{UINT128(0x0000000000000000, 0x6b51cc49737023b3)}, {UINT128(0x0000000000000000, 0xe86b397ace95c470)},
	{UINT128(0xffffffffffffffff, 0x1794c685316a3b90)}, {UINT128(0x0000000000000000, 0x6b51cc49737023b3)},
	{UINT128(0x0000000000000000, 0xbc17d8dc859ad584)}, {UINT128(0x0000000000000000, 0xada859327ba24152)},
	{UINT128(0xffffffffffffffff, 0x5257a6cd845dbeae)}, {UINT128(0x0000000000000000, 0xbc17d8dc859ad584)},
	{UINT128(0x0000000000000000, 0x0a351cb7fc30bc89)}, {UINT128(0x0000000000000000, 0xffcbe2104600a0aa)},
	{UINT128(0xffffffffffffffff, 0x00341defb9ff5f56)}, {UINT128(0x0000000000000000, 0x0a351cb7fc30bc89)},
	{UINT128(0x0000000000000000, 0xffdaaf212fed72dc)}, {UINT128(0x0000000000000000, 0x08a342eda160bf5b)},
	{UINT128(0xffffffffffffffff, 0xf75cbd125e9f40a5)}, {UINT128(0x0000000000000000, 0xffdaaf212fed72dc)},
	{UINT128(0x0000000000000000, 0xaecef739f1a2df11)}, {UINT128(0x0000000000000000, 0xbb062961823b1ddd)},
	{UINT128(0xffffffffffffffff, 0x44f9d69e7dc4e223)}, {UINT128(0x0000000000000000, 0xaecef739f1a2df11)},
	{UINT128(0x0000000000000000, 0xe912ae372d27045e)}, {UINT128(0x0000000000000000, 0x69e4334f6fcc6b41)},
	{UINT128(0xffffffffffffffff, 0x961bccb0903394bf)}, {UINT128(0x0000000000000000, 0xe912ae372d27045e)},
	{UINT128(0x0000000000000000, 0x59ee5272b58f2a5e)}, {UINT128(0x0000000000000000, 0xefaf1b54dd2cdf10)},
	{UINT128(0xffffffffffffffff, 0x1050e4ab22d320f0)}, {UINT128(0x0000000000000000, 0x59ee5272b58f2a5e)},
	{UINT128(0x0000000000000000, 0xf940bfe2304e6c46)}, {UINT128(0x0000000000000000, 0x3a62ee9a597bdc98)},
	{UINT128(0xffffffffffffffff, 0xc59d1165a6842368)}, {UINT128(0x0000000000000000, 0xf940bfe2304e6c46)},
	{UINT128(0x0000000000000000, 0x86f685c25e25acf6)}, {UINT128(0x0000000000000000, 0xd988ad2f9bdf9bbc)},
	{UINT128(0xffffffffffffffff, 0x267752d064206444)}, {UINT128(0x0000000000000000, 0x86f685c25e25acf6)},
	{UINT128(0x0000000000000000, 0xcfefa7898a4ef23d)}, {UINT128(0x0000000000000000, 0x9553b743d75ac040)},
	{UINT128(0xffffffffffffffff, 0x6aac48bc28a53fc0)}, {UINT128(0x0000000000000000, 0xcfefa7898a4ef23d)},
	{UINT128(0x0000000000000000, 0x297162fef3f510db)}, {UINT128(0x0000000000000000, 0xfc9f8a7c2d603c61)},
	{UINT128(0xffffffffffffffff, 0x03607583d29fc39f)}, {UINT128(0x0000000000000000, 0x297162fef3f510db)},
	{UINT128(0x0000000000000000, 0xfdc68c6b356db630)}, {UINT128(0x0000000000000000, 0x21ac9b7964cf0ba4)},
	{UINT128(0xffffffffffffffff, 0xde5364869b30f45c)}, {UINT128(0x0000000000000000, 0xfdc68c6b356db630)},
	{UINT128(0x0000000000000000, 0x9ba298dc0bfc6a89)}, {UINT128(0x0000000000000000, 0xcb41fa15ebff0778)},
	{UINT128(0xffffffffffffffff, 0x34be05ea1400f888)}, {UINT128(0x0000000000000000, 0x9ba298dc0bfc6a89)},
	{UINT128(0x0000000000000000, 0xdd924d05b620678b)}, {UINT128(0x0000000000000000, 0x803a06415c170526)},
	{UINT128(0xffffffffffffffff, 0x7fc5f9bea3e8fada)}, {UINT128(0x0000000000000000, 0xdd924d05b620678b)},
	{UINT128(0x0000000000000000, 0x42013817ad987955)}, {UINT128(0x0000000000000000, 0xf7583a62852a23b3)},
	{UINT128(0xffffffffffffffff, 0x08a7c59d7ad5dc4d)}, {UINT128(0x0000000000000000, 0x42013817ad987955)},
	{UINT128(0x0000000000000000, 0xf2546ffc0e72f287)}, {UINT128(0x0000000000000000, 0x52894f446e0bcbb6)},
	{UINT128(0xffffffffffffffff, 0xad76b0bb91f4344a)}, {UINT128(0x0000000000000000, 0xf2546ffc0e72f287)},
	{UINT128(0x0000000000000000, 0x70fdb51c98c4a59d)}, {UINT128(0x0000000000000000, 0xe5b7105006d4c561)},
	{UINT128(0xffffffffffffffff, 0x1a48efaff92b3a9f)}, {UINT128(0x0000000000000000, 0x70fdb51c98c4a59d)},
	{UINT128(0x0000000000000000, 0xc04c5bab7297d323)}, {UINT128(0x0000000000000000, 0xa8fd40e6ccd52ffe)},
	{UINT128(0xffffffffffffffff, 0x5702bf19332ad002)}, {UINT128(0x0000000000000000, 0xc04c5bab7297d323)},
	{UINT128(0x0000000000000000, 0x107b614e463063b2)}, {UINT128(0x0000000000000000, 0xff780814130c893d)},
	{UINT128(0xffffffffffffffff, 0x0087f7ebecf376c3)}, {UINT128(0x0000000000000000, 0x107b614e463063b2)},
	{UINT128(0x0000000000000000, 0xff1f495f4ec430d8)}, {UINT128(0x0000000000000000, 0x152e774a4d4d0a1c)},
	{UINT128(0xffffffffffffffff, 0xead188b5b2b2f5e4)}, {UINT128(0x0000000000000000, 0xff1f495f4ec430d8)},
	{UINT128(0x0000000000000000, 0xa56bca8b391785dc)}, {UINT128(0x0000000000000000, 0xc36050ecd50ca831)},
	{UINT128(0xffffffffffffffff, 0x3c9faf132af357cf)}, {UINT128(0x0000000000000000, 0xa56bca8b391785dc)},
	{UINT128(0x0000000000000000, 0xe398ac4cf556b733)}, {UINT128(0x0000000000000000, 0x753340aea1827365)},
	{UINT128(0xffffffffffffffff, 0x8accbf515e7d8c9b)}, {UINT128(0x0000000000000000, 0xe398ac4cf556b733)},
	{UINT128(0x0000000000000000, 0x4e0fd78d4edaaa58)}, {UINT128(0x0000000000000000, 0xf3ced94d28b2ce8b)},
	{UINT128(0xffffffffffffffff, 0x0c3126b2d74d3175)}, {UINT128(0x0000000000000000, 0x4e0fd78d4edaaa58)},
	{UINT128(0x0000000000000000, 0xf6167a58d7b59027)}, {UINT128(0x0000000000000000, 0x468bdfefa3c90d65)},
	{UINT128(0xffffffffffffffff, 0xb97420105c36f29b)}, {UINT128(0x0000000000000000, 0xf6167a58d7b59027)},
	{UINT128(0x0000000000000000, 0x7c20641affdff6b4)}, {UINT128(0x0000000000000000, 0xdfe4e92d0d8a37f6)},
	{UINT128(0xffffffffffffffff, 0x201b16d2f275c80a)}, {UINT128(0x0000000000000000, 0x7c20641affdff6b4)},
	{UINT128(0x0000000000000000, 0xc85bca1abaf7f0a8)}, {UINT128(0x0000000000000000, 0x9f599f55f0340062)},
	{UINT128(0xffffffffffffffff, 0x60a660aa0fcbff9e)}, {UINT128(0x0000000000000000, 0xc85bca1abaf7f0a8)},
	{UINT128(0x0000000000000000, 0x1cff533b307dc132)}, {UINT128(0x0000000000000000, 0xfe5a381c8a1aa225)},
	{UINT128(0xffffffffffffffff, 0x01a5c7e375e55ddb)}, {UINT128(0x0000000000000000, 0x1cff533b307dc132)},
	{UINT128(0x0000000000000000, 0xfbd14ce0d191516f)}, {UINT128(0x0000000000000000, 0x2e15fb1a1189fe94)},
	{UINT128(0xffffffffffffffff, 0xd1ea04e5ee76016c)}, {UINT128(0x0000000000000000, 0xfbd14ce0d191516f)},
	{UINT128(0x0000000000000000, 0x91796b31609f0c55)}, {UINT128(0x0000000000000000, 0xd2a6488487a91919)},
	{UINT128(0xffffffffffffffff, 0x2d59b77b7856e6e7)}, {UINT128(0x0000000000000000, 0x91796b31609f0c55)},
	{UINT128(0x0000000000000000, 0xd703479a2f6776cd)}, {UINT128(0x0000000000000000, 0x8af1b726df15e13d)},
	{UINT128(0xffffffffffffffff, 0x750e48d920ea1ec3)}, {UINT128(0x0000000000000000, 0xd703479a2f6776cd)},
	{UINT128(0x0000000000000000, 0x35c9e3abe77356d6)}, {UINT128(0x0000000000000000, 0xfa491035e55da3a4)},
	{UINT128(0xffffffffffffffff, 0x05b6efca1aa25c5c)}, {UINT128(0x0000000000000000, 0x35c9e3abe77356d6)},
	{UINT128(0x0000000000000000, 0xedfcf21cabacd3b2)}, {UINT128(0x0000000000000000, 0x5e53d7984f7eb8c2)},
	{UINT128(0xffffffffffffffff, 0xa1ac2867b081473e)}, {UINT128(0x0000000000000000, 0xedfcf21cabacd3b2)},
	{UINT128(0x0000000000000000, 0x659556deae523d7b)}, {UINT128(0x0000000000000000, 0xeafb8b944453f530)},
	{UINT128(0xffffffffffffffff, 0x1504746bbbac0ad0)}, {UINT128(0x0000000000000000, 0x659556deae523d7b)},
	{UINT128(0x0000000000000000, 0xb7c654ce4adba9f3)}, {UINT128(0x0000000000000000, 0xb238aa1bfa9ad508)},
	{UINT128(0xffffffffffffffff, 0x4dc755e405652af8)}, {UINT128(0x0000000000000000, 0xb7c654ce4adba9f3)},
	{UINT128(0x0000000000000000, 0x03ed452d5732f951)}, {UINT128(0x0000000000000000, 0xfff84a1e29de8572)},
	{UINT128(0xffffffffffffffff, 0x0007b5e1d6217a8e)}, {UINT128(0x0000000000000000, 0x03ed452d5732f951)},
	{UINT128(0x0000000000000000, 0xfff0e343865bbb14)}, {UINT128(0x0000000000000000, 0x057f53487eac8978)},
	{UINT128(0xffffffffffffffff, 0xfa80acb781537688)}, {UINT128(0x0000000000000000, 0xfff0e343865bbb14)},
	{UINT128(0x0000000000000000, 0xb117227f6117f9fa)}, {UINT128(0x0000000000000000, 0xb8dd64b0720df648)},
	{UINT128(0xffffffffffffffff, 0x47229b4f8df209b8)}, {UINT128(0x0000000000000000, 0xb117227f6117f9fa)},
	{UINT128(0x0000000000000000, 0xea5ad8d8c3a91f06)}, {UINT128(0x0000000000000000, 0x6705f51037e7c8ca)},
	{UINT128(0xffffffffffffffff, 0x98fa0aefc8183736)}, {UINT128(0x0000000000000000, 0xea5ad8d8c3a91f06)},
	{UINT128(0x0000000000000000, 0x5cdd8f2498424757)}, {UINT128(0x0000000000000000, 0xee8ff79c548acd10)},
	{UINT128(0xffffffffffffffff, 0x11700863ab7532f0)}, {UINT128(0x0000000000000000, 0x5cdd8f2498424757)},
	{UINT128(0x0000000000000000, 0xf9f35de0dde328ac)}, {UINT128(0x0000000000000000, 0x3752c669e2bb5946)},
	{UINT128(0xffffffffffffffff, 0xc8ad39961d44a6ba)}, {UINT128(0x0000000000000000, 0xf9f35de0dde328ac)},
	{UINT128(0x0000000000000000, 0x899f4e7f712a765f)}, {UINT128(0x0000000000000000, 0xd7dc7ec4fabdb012)},
	{UINT128(0xffffffffffffffff, 0x2823813b05424fee)}, {UINT128(0x0000000000000000, 0x899f4e7f712a765f)},
	{UINT128(0x0000000000000000, 0xd1c0c252c9de2c87)}, {UINT128(0x0000000000000000, 0x92c39a65db88809e)},
	{UINT128(0xffffffffffffffff, 0x6d3c659a24777f62)}, {UINT128(0x0000000000000000, 0xd1c0c252c9de2c87)},
	{UINT128(0x0000000000000000, 0x2c8a35063b061a6c)}, {UINT128(0x0000000000000000, 0xfc187a52063060c3)},
	{UINT128(0xffffffffffffffff, 0x03e785adf9cf9f3d)}, {UINT128(0x0000000000000000, 0x2c8a35063b061a6c)},
	{UINT128(0x0000000000000000, 0xfe2b71dbecd7aefd)}, {UINT128(0x0000000000000000, 0x1e8eb7fde4aa3eb8)},
	{UINT128(0xffffffffffffffff, 0xe17148021b55c148)}, {UINT128(0x0000000000000000, 0xfe2b71dbecd7aefd)},
	{UINT128(0x0000000000000000, 0x9e1e224c0e28bc95)}, {UINT128(0x0000000000000000, 0xc95520fe2d40a74c)},
	{UINT128(0xffffffffffffffff, 0x36aadf01d2bf58b4)}, {UINT128(0x0000000000000000, 0x9e1e224c0e28bc95)},
	{UINT128(0x0000000000000000, 0xdf20db088aa60405)}, {UINT128(0x0000000000000000, 0x7d7f7b995b368ba2)},
	{UINT128(0xffffffffffffffff, 0x82808466a4c9745e)}, {UINT128(0x0000000000000000, 0xdf20db088aa60405)},
	{UINT128(0x0000000000000000, 0x4508fbbf1a4de123)}, {UINT128(0x0000000000000000, 0xf6841af4cc8048a5)},
	{UINT128(0xffffffffffffffff, 0x097be50b337fb75b)}, {UINT128(0x0000000000000000, 0x4508fbbf1a4de123)},
	{UINT128(0x0000000000000000, 0xf3530e2aea9d3967)}, {UINT128(0x0000000000000000, 0x4f8e6fa5bb09c582)},
	{UINT128(0xffffffffffffffff, 0xb071905a44f63a7e)}, {UINT128(0x0000000000000000, 0xf3530e2aea9d3967)},
	{UINT128(0x0000000000000000, 0x73cd2ebb873482ef)}, {UINT128(0x0000000000000000, 0xe44fac3814e09857)},
	{UINT128(0xffffffffffffffff, 0x1bb053c7eb1f67a9)}, {UINT128(0x0000000000000000, 0x73cd2ebb873482ef)},
	{UINT128(0x0000000000000000, 0xc25b888d7c1fcd39)}, {UINT128(0x0000000000000000, 0xa69de36ac4fbfadd)},
	{UINT128(0xffffffffffffffff, 0x59621c953b040523)}, {UINT128(0x0000000000000000, 0xc25b888d7c1fcd39)},
	{UINT128(0x0000000000000000, 0x139d9f12c5a29905)}, {UINT128(0x0000000000000000, 0xff3f542a416b0135)},
	{UINT128(0xffffffffffffffff, 0x00c0abd5be94fecb)}, {UINT128(0x0000000000000000, 0x139d9f12c5a29905)},
	{UINT128(0x0000000000000000, 0xff5ce92982087868)}, {UINT128(0x0000000000000000, 0x120c9674ed444c82)},
	{UINT128(0xffffffffffffffff, 0xedf3698b12bbb37e)}, {UINT128(0x0000000000000000, 0xff5ce92982087868)},
	{UINT128(0x0000000000000000, 0xa7ce612e65243292)}, {UINT128(0x0000000000000000, 0xc154e09faa2ff69b)},
	{UINT128(0xffffffffffffffff, 0x3eab1f6055d00965)}, {UINT128(0x0000000000000000, 0xa7ce612e65243292)},
	{UINT128(0x0000000000000000, 0xe50478cdce2246bd)}, {UINT128(0x0000000000000000, 0x7265ff0e194e229e)},
	{UINT128(0xffffffffffffffff, 0x8d9a00f1e6b1dd62)}, {UINT128(0x0000000000000000, 0xe50478cdce2246bd)},
	{UINT128(0x0000000000000000, 0x510c437224dbc18c)}, {UINT128(0x0000000000000000, 0xf2d4eaa8233e997e)},
	{UINT128(0xffffffffffffffff, 0x0d2b1557dcc16682)}, {UINT128(0x0000000000000000, 0x510c437224dbc18c)},
	{UINT128(0x0000000000000000, 0xf6ef5b503c32858a)}, {UINT128(0x0000000000000000, 0x43856d385d2736a5)},
	{UINT128(0xffffffffffffffff, 0xbc7a92c7a2d8c95b)}, {UINT128(0x0000000000000000, 0xf6ef5b503c32858a)},
	{UINT128(0x0000000000000000, 0x7edd5d70934b7753)}, {UINT128(0x0000000000000000, 0xde5aa65869193806)},
	{UINT128(0xffffffffffffffff, 0x21a559a796e6c7fa)}, {UINT128(0x0000000000000000, 0x7edd5d70934b7753)},
	{UINT128(0x0000000000000000, 0xca4c871d625361aa)}, {UINT128(0x0000000000000000, 0x9ce11f1eb18147b2)},
	{UINT128(0xffffffffffffffff, 0x631ee0e14e7eb84e)}, {UINT128(0x0000000000000000, 0xca4c871d625361aa)},
	{UINT128(0x0000000000000000, 0x201dd15adf70b65a)}, {UINT128(0x0000000000000000, 0xfdfa3878546c3d29)},
	{UINT128(0xffffffffffffffff, 0x0205c787ab93c2d7)}, {UINT128(0x0000000000000000, 0x201dd15adf70b65a)},
	{UINT128(0x0000000000000000, 0xfc5d39be5a5bbc4c)}, {UINT128(0x0000000000000000, 0x2afe010ca997bc65)},
	{UINT128(0xffffffffffffffff, 0xd501fef35668439b)}, {UINT128(0x0000000000000000, 0xfc5d39be5a5bbc4c)},
	{UINT128(0x0000000000000000, 0x940c5f7a69e5ce1d)}, {UINT128(0x0000000000000000, 0xd0d93696053af099)},
	{UINT128(0xffffffffffffffff, 0x2f26c969fac50f67)}, {UINT128(0x0000000000000000, 0x940c5f7a69e5ce1d)},
	{UINT128(0x0000000000000000, 0xd8b3a1526517a48c)}, {UINT128(0x0000000000000000, 0x884b9246854ab50c)},
	{UINT128(0xffffffffffffffff, 0x77b46db97ab54af4)}, {UINT128(0x0000000000000000, 0xd8b3a1526517a48c)},
	{UINT128(0x0000000000000000, 0x38db20a6bae9b9c9)}, {UINT128(0x0000000000000000, 0xf99b42d1d57781ec)},
	{UINT128(0xffffffffffffffff, 0x0664bd2e2a887e14)}, {UINT128(0x0000000000000000, 0x38db20a6bae9b9c9)},
	{UINT128(0x0000000000000000, 0xef20b07b6c6c0b38)}, {UINT128(0x0000000000000000, 0x5b66618e2832d7f8)},
	{UINT128(0xffffffffffffffff, 0xa4999e71d7cd2808)}, {UINT128(0x0000000000000000, 0xef20b07b6c6c0b38)},
	{UINT128(0x0000000000000000, 0x6875950ed42aff7f)}, {UINT128(0x0000000000000000, 0xe9b7e3de5fdedc8c)},
	{UINT128(0xffffffffffffffff, 0x16481c21a0212374)}, {UINT128(0x0000000000000000, 0x6875950ed42aff7f)},
	{UINT128(0x0000000000000000, 0xb9f2ac703cca0db4)}, {UINT128(0x0000000000000000, 0xaff3e5ef2b507c07)},
	{UINT128(0xffffffffffffffff, 0x500c1a10d4af83f9)}, {UINT128(0x0000000000000000, 0xb9f2ac703cca0db4)},
	{UINT128(0x0000000000000000, 0x071153d3394ec723)}, {UINT128(0x0000000000000000, 0xffe704e71533c509)},
	{UINT128(0xffffffffffffffff, 0x0018fb18eacc3af7)}, {UINT128(0x0000000000000000, 0x071153d3394ec723)},
	{UINT128(0x0000000000000000, 0xffba9dd8dc8b1e84)}, {UINT128(0x0000000000000000, 0x0bc6dd52c3a342ec)},
	{UINT128(0xffffffffffffffff, 0xf43922ad3c5cbd14)}, {UINT128(0x0000000000000000, 0xffba9dd8dc8b1e84)},
	{UINT128(0x0000000000000000, 0xac800eafb91ef9aa)}, {UINT128(0x0000000000000000, 0xbd27b83dfcbe927a)},
	{UINT128(0xffffffffffffffff, 0x42d847c203416d86)}, {UINT128(0x0000000000000000, 0xac800eafb91ef9aa)},
	{UINT128(0x0000000000000000, 0xe7c187467233d509)}, {UINT128(0x0000000000000000, 0x6cbe5c76cc65c9fb)},
	{UINT128(0xffffffffffffffff, 0x9341a389339a3605)}, {UINT128(0x0000000000000000, 0xe7c187467233d509)},
	{UINT128(0x0000000000000000, 0x56fb9e2e76ced87a)}, {UINT128(0x0000000000000000, 0xf0c5017ee336ca10)},
	{UINT128(0xffffffffffffffff, 0x0f3afe811cc935f0)}, {UINT128(0x0000000000000000, 0x56fb9e2e76ced87a)},
	{UINT128(0x0000000000000000, 0xf88485e44c7af48b)}, {UINT128(0x0000000000000000, 0x3d70d68c5bc68f21)},
	{UINT128(0xffffffffffffffff, 0xc28f2973a43970df)}, {UINT128(0x0000000000000000, 0xf88485e44c7af48b)},
	{UINT128(0x0000000000000000, 0x844889019584609a)}, {UINT128(0x0000000000000000, 0xdb2c78a7ede238aa)},
	{UINT128(0xffffffffffffffff, 0x24d38758121dc756)}, {UINT128(0x0000000000000000, 0x844889019584609a)},
	{UINT128(0x0000000000000000, 0xce16888783ae13b4)}, {UINT128(0x0000000000000000, 0x97de125a2080a8ee)},
	{UINT128(0xffffffffffffffff, 0x6821eda5df7f5712)}, {UINT128(0x0000000000000000, 0xce16888783ae13b4)},
	{UINT128(0x0000000000000000, 0x2656f7f28a2d4e96)}, {UINT128(0x0000000000000000, 0xfd1cdd63d0bc8736)},
	{UINT128(0xffffffffffffffff, 0x02e3229c2f4378ca)}, {UINT128(0x0000000000000000, 0x2656f7f28a2d4e96)},
	{UINT128(0x0000000000000000, 0xfd57de5867eedc3a)}, {UINT128(0x0000000000000000, 0x24c9329bfa6812d4)},
	{UINT128(0xffffffffffffffff, 0xdb36cd640597ed2c)}, {UINT128(0x0000000000000000, 0xfd57de5867eedc3a)},
	{UINT128(0x0000000000000000, 0x99210f624d30facc)}, {UINT128(0x0000000000000000, 0xcd26fd2158358e7e)},
	{UINT128(0xffffffffffffffff, 0x32d902dea7ca7182)}, {UINT128(0x0000000000000000, 0x99210f624d30facc)},
	{UINT128(0x0000000000000000, 0xdbfb34373c974b0f)}, {UINT128(0x0000000000000000, 0x82ef9f618dc5b70f)},
	{UINT128(0xffffffffffffffff, 0x7d10609e723a48f1)}, {UINT128(0x0000000000000000, 0xdbfb34373c974b0f)},
	{UINT128(0x0000000000000000, 0x3ef6e9017a700333)}, {UINT128(0x0000000000000000, 0xf822d0a67b9c5cb6)},
	{UINT128(0xffffffffffffffff, 0x07dd2f598463a34a)}, {UINT128(0x0000000000000000, 0x3ef6e9017a700333)},
	{UINT128(0x0000000000000000, 0xf14c7a21c8cbd0f5)}, {UINT128(0x0000000000000000, 0x5581004bd19ed1a1)},
	{UINT128(0xffffffffffffffff, 0xaa7effb42e612e5f)}, {UINT128(0x0000000000000000, 0xf14c7a21c8cbd0f5)},
	{UINT128(0x0000000000000000, 0x6e29e053f55a4d9b)}, {UINT128(0x0000000000000000, 0xe715993ccd02fe9d)},
	{UINT128(0xffffffffffffffff, 0x18ea66c332fd0163)}, {UINT128(0x0000000000000000, 0x6e29e053f55a4d9b)},
	{UINT128(0x0000000000000000, 0xbe35c4e716999631)}, {UINT128(0x0000000000000000, 0xab561a8cbb24f411)},
	{UINT128(0xffffffffffffffff, 0x54a9e57344db0bef)}, {UINT128(0x0000000000000000, 0xbe35c4e716999631)},
	{UINT128(0x0000000000000000, 0x0d5880deafc18b54)}, {UINT128(0x0000000000000000, 0xffa6e2a58df6947e)},
	{UINT128(0xffffffffffffffff, 0x00591d5a72096b82)}, {UINT128(0x0000000000000000, 0x0d5880deafc18b54)},
	{UINT128(0x0000000000000000, 0xfed7d3a8a29c603c)}, {UINT128(0x0000000000000000, 0x184f8712c130a087)},
	{UINT128(0xffffffffffffffff, 0xe7b078ed3ecf5f79)}, {UINT128(0x0000000000000000, 0xfed7d3a8a29c603c)},
	{UINT128(0x0000000000000000, 0xa302d34959951244)}, {UINT128(0x0000000000000000, 0xc56438f6f0ec3ccb)},
	{UINT128(0xffffffffffffffff, 0x3a9bc7090f13c335)}, {UINT128(0x0000000000000000, 0xa302d34959951244)},
	{UINT128(0x0000000000000000, 0xe224198a0e002124)}, {UINT128(0x0000000000000000, 0x77fbfd9aa5077659)},
	{UINT128(0xffffffffffffffff, 0x880402655af889a7)}, {UINT128(0x0000000000000000, 0xe224198a0e002124)},
	{UINT128(0x0000000000000000, 0x4b10693a5514819d)}, {UINT128(0x0000000000000000, 0xf4bf61b00b5982b8)},
	{UINT128(0xffffffffffffffff, 0x0b409e4ff4a67d48)}, {UINT128(0x0000000000000000, 0x4b10693a5514819d)},
	{UINT128(0x0000000000000000, 0xf5341c9f32bffeba)}, {UINT128(0x0000000000000000, 0x498f9a655552e91e)},
	{UINT128(0xffffffffffffffff, 0xb670659aaaad16e2)}, {UINT128(0x0000000000000000, 0xf5341c9f32bffeba)},
	{UINT128(0x0000000000000000, 0x795ea1b4f360936a)}, {UINT128(0x0000000000000000, 0xe1668a498f1f892d)},
	{UINT128(0xffffffffffffffff, 0x1e9975b670e076d3)}, {UINT128(0x0000000000000000, 0x795ea1b4f360936a)},
	{UINT128(0x0000000000000000, 0xc66353a8c232a43d)}, {UINT128(0x0000000000000000, 0xa1cbfad9521bfd1c)},
	{UINT128(0xffffffffffffffff, 0x5e340526ade402e4)}, {UINT128(0x0000000000000000, 0xc66353a8c232a43d)},
	{UINT128(0x0000000000000000, 0x19dfb6eb24a85c3e)}, {UINT128(0x0000000000000000, 0xfeb0696d3ae4f04e)},
	{UINT128(0xffffffffffffffff, 0x014f9692c51b0fb2)}, {UINT128(0x0000000000000000, 0x19dfb6eb24a85c3e)},
	{UINT128(0x0000000000000000, 0xfb3baab441e770f8)}, {UINT128(0x0000000000000000, 0x312c2e4f88300850)},
	{UINT128(0xffffffffffffffff, 0xced3d1b077cff7b0)}, {UINT128(0x0000000000000000, 0xfb3baab441e770f8)},
	{UINT128(0x0000000000000000, 0x8ee0db26e24390f9)}, {UINT128(0x0000000000000000, 0xd46b3b72a2d68fca)},
	{UINT128(0xffffffffffffffff, 0x2b94c48d5d297036)}, {UINT128(0x0000000000000000, 0x8ee0db26e24390f9)},
	{UINT128(0x0000000000000000, 0xd54aa3d165cc7019)}, {UINT128(0x0000000000000000, 0x8d9280b89b9df49c)},
	{UINT128(0xffffffffffffffff, 0x726d7f4764620b64)}, {UINT128(0x0000000000000000, 0xd54aa3d165cc7019)},
	{UINT128(0x0000000000000000, 0x32b693d36c55f3de)}, {UINT128(0x0000000000000000, 0xfaed376a1b42e55a)},
	{UINT128(0xffffffffffffffff, 0x0512c895e4bd1aa6)}, {UINT128(0x0000000000000000, 0x32b693d36c55f3de)},
	{UINT128(0x0000000000000000, 0xecd006ec59ea3070)}, {UINT128(0x0000000000000000, 0x613daaabce40fcb7)},
	{UINT128(0xffffffffffffffff, 0x9ec2555431bf0349)}, {UINT128(0x0000000000000000, 0xecd006ec59ea3070)},
	{UINT128(0x0000000000000000, 0x62b12e1b57b51118)}, {UINT128(0x0000000000000000, 0xec3624222d227bd2)},
	{UINT128(0xffffffffffffffff, 0x13c9dbddd2dd842e)}, {UINT128(0x0000000000000000, 0x62b12e1b57b51118)},
	{UINT128(0x0000000000000000, 0xb592e7697f14dd4b)}, {UINT128(0x0000000000000000, 0xb4768f550ce389fe)},
	{UINT128(0xffffffffffffffff, 0x4b8970aaf31c7602)}, {UINT128(0x0000000000000000, 0xb592e7697f14dd4b)},
	{UINT128(0x0000000000000000, 0x00c90fc5f66525d3)}, {UINT128(0x0000000000000000, 0xffffb10b10e80e96)},
	{UINT128(0xffffffffffffffff, 0x00004ef4ef17f16a)}, {UINT128(0x0000000000000000, 0x00c90fc5f66525d3)},
};


const fpr fpr_p2_tab[] = {
    { UINT128(2, 0) },
    { UINT128(1, 0) },
    { UINT128(0, 0x8000000000000000) },
    { UINT128(0, 0x4000000000000000) },
    { UINT128(0, 0x2000000000000000) },
    { UINT128(0, 0x1000000000000000) },
    { UINT128(0, 0x0800000000000000) },
    { UINT128(0, 0x0400000000000000) },
    { UINT128(0, 0x0200000000000000) },
    { UINT128(0, 0x0100000000000000) },
};

#else // yyyFPNATIVE+0 yyyFPEMU+0 yyyFXP+0

#error No FP implementation selected

#endif // yyyFPNATIVE- yyyFPEMU- yyyFXP-
