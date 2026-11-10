#include "ng_inner.h"

#include <stdio.h>

/* 
 * We use the following constant to reject some values during the keygen.
 * They were carefully chosen to meet the precision requirement for a secure
 * fixed point implementation as described in the paper 
 * "Fixed-Point Implementation of the Falcon Signature Scheme" from CRYPTO 2026.
 * Values for n < 512 are conservative: they match the value for n = 512
*/
static const fxr fxr_q_gamma_hybrid2[] = {
	FXR(0), /* unused */ 
	FXR((12289*4*4ULL) << 32), 
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32),
	FXR((12289*4*4ULL) << 32), 
	FXR((12289*4*4ULL) << 32), /* for n = 1024, we need to figure out the best value */
};
/* the lower bound is precomputed as q/gamma² */
static const fxr fxr_q_div_gamma_hybrid2[] = {
	FXR(0), /* unused */ 
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625, 
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625,
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625, 
	FXR(0x0000030010000000ULL), // 12289/4*4 = 768.0625, /* for n = 1024, we need to figure out the best value */
};
static const fxr fxr_gamma_root_sqr[] = {
	FXR(0), /* unused */
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32),
	FXR((32*32ULL) << 32), /* for n = 1024, we need to figure out the best value */
};
static const fxr fxr_gamma_fg_sqr[] = {
	FXR(0), /* unused */
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32),
	FXR((256*256ULL) << 32), /* for n = 1024, we need to figure out the best value */
};
static const fxr fxr_gamma_FG_sqr[] = {
	FXR(0) ,/* unused */
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32),
	FXR((4096*4096ULL) << 32), /* for n = 1024, we need to figure out the best value */
};


const ntru_profile SOLVE_Falcon_256 = {
	12289,
	2, 8,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	10,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3 }
};

const ntru_profile SOLVE_Falcon_512 = {
	12289,
	9, 9,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	9,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3 }
};

const ntru_profile SOLVE_Falcon_1024 = {
	12289,
	10, 10,
	{ 1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207 },
	{ 1, 2, 3, 6, 11, 21, 40, 78, 155, 308 },
	{ 1, 1, 2, 2, 2, 3, 3, 4, 5, 7 },
	8,
	{ 0, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127 },
	{ 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 4 }
};

/* Falcon, q = 12289, n = 256 -> kmax = 24 */
const uint16_t gauss_Falcon_256[] = {
	24,
	    1,     3,     6,    11,    22,    40,    73,   129,
	  222,   371,   602,   950,  1460,  2183,  3179,  4509,
	 6231,  8395, 11032, 14150, 17726, 21703, 25995, 30487,
	35048, 39540, 43832, 47809, 51385, 54503, 57140, 59304,
	61026, 62356, 63352, 64075, 64585, 64933, 65164, 65313,
	65406, 65462, 65495, 65513, 65524, 65529, 65532, 65534
};

/* Falcon, q = 12289, n = 512 -> kmax = 17 */
const uint16_t gauss_Falcon_512[] = {
	17,
	    1,     4,    11,    28,    65,   146,   308,   615,
	 1164,  2083,  3535,  5692,  8706, 12669, 17574, 23285,
	29542, 35993, 42250, 47961, 52866, 56829, 59843, 62000,
	63452, 64371, 64920, 65227, 65389, 65470, 65507, 65524,
	65531, 65534
};

/* Falcon, q = 12289, n = 1024 -> kmax = 12 */
const uint16_t gauss_Falcon_1024[] = {
	12,
	    2,     8,    28,    94,   280,   742,  1761,  3753,
	 7197, 12472, 19623, 28206, 37329, 45912, 53063, 58338,
	61782, 63774, 64793, 65255, 65441, 65507, 65527, 65533
};

/* 
 * Compute the square of norm of the i-th element of the FFT embedding
 * of a polynomial f.
 * This function accounts for the internal layout of the FFT 
 * representation, where Re is stored at i, and Im is stored at i+(n/2).
 * Returns f[i]² + f[i+hn]² = Re(f_i)² + Im(f_i)²
*/
static inline fxr
norm_C_sqr(const fxr *f, unsigned n, unsigned i)
{
	unsigned hn = n >> 1;
	return fxr_add(fxr_sqr(f[i]), fxr_sqr(f[i + hn]));
}

/* 
 * Given the FFT embedding of a polynomial f, check that it's norm infinity is 
 * within a given range.
 * We only iterate over the n/2 first values since other are just conjugate.
 * This the bound is a predefined constant, we compute the square of the norm 
 * and compare it to the square of our bound to spare a square root computation. 
*/
static int
vect_check_norm_inf_fft(unsigned logn, const fxr *f, const fxr *bound)
{
	size_t n = (size_t)1 << (logn);
	size_t hn = (size_t)1 << (logn - 1);
	for (size_t i = 0; i < hn; i++) {
		if (fxr_lt(bound[logn], norm_C_sqr(f, n, i))) {
			return 1;
		}
	}
	return 0;
}

/* 
 * alpha_hybrid is defined in Antrag and Mitaka as the square root of the maximum, 
 * for i = 0 to n, of the max((|phi_i(f)|² + |phi_i(g)|²)/q , q/(|phi_i(f)|² + |phi_i(g)|²)
 * with phi_i being the FFT embedding.
 * We did a few optimization to avoid inversions and square root. 
 * Let x = (|phi_i(f)|² + |phi_i(g)|²)_i
 * Since alpha_hybrid > 0, we can check alpha_hybrid² < gamma_hybrid² instead,
 * avoiding a square root. This is equivalent to verifying that for every element
 *     max⁡(x[i]/q,  q/x[i]) < gamma²
 * which is equivalent to checking that 
 *     q/gamma² < x[i] < q*gamma², for all i
 * Since both q and gamma are constant, we precomput the bounds.
*/
static int
check_alpha_hybrid(const fxr *f_fft, const fxr *g_fft, unsigned logn)
{
	const unsigned n = 1u << logn;
	const unsigned hn = n >> 1;

	for (unsigned u = 0; u < hn; u++) {
		/* x = |FFT(f[i])|² + |FFT(g[i])|² */
		fxr x = fxr_add(
			norm_C_sqr(f_fft, n, u), 
			norm_C_sqr(g_fft, n, u)
		);
		/* As per alpha_hybrid definition, we should check if either a=X[u]/q 
		 * or b=q/X[u] cross bound for any index u.
		 * First, we can avoid a division in a by comparing X[u] to q*bound. 
		 * This gives us the check x > upper as a stop condition.
		 * Second, we have b = 1/a. This gives us the check 
		 * x < q/bound as the second stop condition.
		 */
		if (!fxr_lt(x, fxr_q_gamma_hybrid2[logn]) || 
			!fxr_lt(fxr_q_div_gamma_hybrid2[logn], x)) {
			return 1;
		}
	}

	return 0;
}

/* see ntrugen.h */
int
Falcon_keygen(unsigned logn,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	ntrugen_rng rng, void *restrict rng_context,
	void *restrict tmp, size_t tmp_len)
{
#ifdef DEBUG
	uint64_t rejections_alpha_hybrid = 0;
	uint64_t rejections_alpha_gpv = 0;
	uint64_t rejections_norm_fg = 0;
	uint64_t rejections_norm_FG = 0;
	uint64_t rejections_norm_k = 0;
#endif
	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 7) {
		return -1;
	}
	if (logn < 2 || logn > 10) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (void *)utmp2;
	/* Keep 2*n raw completion bytes and six n-element fxr vectors. */
	if (tmp_len < ((size_t)50 << logn)) {
		return -1;
	}

	const ntru_profile *prof;
	for (;;) {
		/*
		 * Generate f and g.
		 */
		switch (logn) {
		case 8:
			prof = &SOLVE_Falcon_256;
			gauss_sample_poly(logn, f, gauss_Falcon_256,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_Falcon_256,
				rng, rng_context);
			break;
		case 9:
			prof = &SOLVE_Falcon_512;
			gauss_sample_poly(logn, f, gauss_Falcon_512,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_Falcon_512,
				rng, rng_context);
			break;
		case 10:
			prof = &SOLVE_Falcon_1024;
			gauss_sample_poly(logn, f, gauss_Falcon_1024,
				rng, rng_context);
			gauss_sample_poly(logn, g, gauss_Falcon_1024,
				rng, rng_context);
			break;
		default:
			prof = &SOLVE_Falcon_256;
			if (logn < 2 || logn > 7) {
				return -1;
			}
			gauss_sample_poly_reduced(logn, f, gauss_Falcon_256,
				rng, rng_context);
			gauss_sample_poly_reduced(logn, g, gauss_Falcon_256,
				rng, rng_context);
			break;
		}

		/*
		 * Check that (g,-f) has an acceptable norm (norm must be
		 * less than 1.17*sqrt(q)).
		 */
		if ((poly_sqnorm(logn, f) + poly_sqnorm(logn, g)) >= 16823) {
#ifdef DEBUG
			rejections_alpha_gpv++;
			if (rejections_alpha_gpv%1000 == 0)
				printf("\trejected alpha_gpv %lu times\n", rejections_alpha_gpv);
#endif
				continue;
		}

		/*
		 * Check that FFT(f) and FFT(g) coefficient are within range
		 * of acceptable values for our fixed point implementation.
		*/
		size_t n = (size_t)1 << logn;
		fxr *rt1 = (fxr *)tt32;
		fxr *rt2 = rt1 + n;
		vect_set(logn, rt1, f);
		vect_set(logn, rt2, g);
		vect_FFT(logn, rt1);
		vect_FFT(logn, rt2);
		if (vect_check_norm_inf_fft(logn, rt1, fxr_gamma_fg_sqr) ||
				vect_check_norm_inf_fft(logn, rt2, fxr_gamma_fg_sqr)) {
#ifdef DEBUG
			rejections_norm_fg++;
			if (rejections_norm_fg%1000 == 0)
				printf("\trejected ||FFT(f), FFT(g)|| bound %lu times\n", rejections_norm_fg);
#endif
				continue;
		}

		if (check_alpha_hybrid(rt1, rt2, logn)) {
#ifdef DEBUG
			rejections_alpha_hybrid++;
			if (rejections_alpha_hybrid%1000 == 0)
				printf("\trejected rejections_alpha_hybrid %lu times\n", rejections_alpha_hybrid);
#endif
				continue;
		}

		/*
		 * For Falcon, we need to check that the orthogonalized
		 * vector also has an acceptable norm.
		 * At this point we already have the following variables layout:
		 * rt1 = FFT(f)
		 * rt2 = FFT(g)
		 */
		fxr *rt3 = rt2 + n;
		vect_invnorm_fft(logn, rt3, rt1, rt2, 0);
		vect_adj_fft(logn, rt1);
		vect_adj_fft(logn, rt2);
		vect_mul_realconst(logn, rt1, fxr_of(12289));
		vect_mul_realconst(logn, rt2, fxr_of(12289));
		vect_mul_autoadj_fft(logn, rt1, rt3);
		vect_mul_autoadj_fft(logn, rt2, rt3);
		vect_iFFT(logn, rt1);
		vect_iFFT(logn, rt2);
		fxr sn = fxr_zero;
		for (size_t u = 0; u < n; u ++) {
			sn = fxr_add(sn,
				fxr_add(fxr_sqr(rt1[u]), fxr_sqr(rt2[u])));
		}
		if (!fxr_lt(sn, fxr_of_scaled32(72251709809335))) {
#ifdef DEBUG
			rejections_alpha_gpv++;
			if (rejections_alpha_gpv%1000 == 0)
				printf("\trejected alpha_gpv %lu times\n", rejections_alpha_gpv);
#endif
				continue;
		}

		/*
		 * Check that f is invertible modulo X^n+1 and modulo q.
		 *
		 * We compute modulo the non-prime 2147465883, which is
		 * a multiple of 12289.
		 */
		if (!poly_is_invertible(logn, f,
			2147465883, 2763744365, 248710,
			12289, 2863078533, 45, tt32))
		{
			continue;
		}

		/*
		 * Solve the NTRU equation.
		 */
#if NTRUGEN_STATS
		stats_solve_attempt ++;
#endif

		int err;
        if (logn == 10) {
            ntru_profile wide = *prof;
            wide.reduce_bits = 20;
            err = solve_NTRU_wide(&wide, logn, f, g, tt32);
        } else {
            err = solve_NTRU(prof, logn, f, g, tt32);
        }
		switch (err) {
		case SOLVE_OK:
#if NTRUGEN_STATS
			stats_solve_success ++;
#endif
			break;
#if NTRUGEN_STATS
		case SOLVE_ERR_GCD:
			stats_solve_err_gcd ++;
			continue;
		case SOLVE_ERR_REDUCE:
			stats_solve_err_reduce ++;
			continue;
		case SOLVE_ERR_LIMIT:
			stats_solve_err_limit ++;
			continue;
#endif
		default:
			continue;
		}

		/*
		 * Return the computed F and G.
		 */
		int8_t *tF = (int8_t *)tt32;
		int8_t *tG = tF + n;
		const int8_t *F_chk = tF;
		const int8_t *G_chk = tG;

		/* 
		* Check FFT range for fixed point implementation.
		*/
		fxr *F_fft, *G_fft;
		F_fft = (fxr *)(void *)(tG + n);
		G_fft = F_fft + n;
		vect_set(logn, F_fft, F_chk);
		vect_set(logn, G_fft, G_chk);
		vect_FFT(logn, F_fft);
		vect_FFT(logn, G_fft);
		if (vect_check_norm_inf_fft(logn, F_fft, fxr_gamma_FG_sqr) ||
				vect_check_norm_inf_fft(logn, G_fft, fxr_gamma_FG_sqr)) {
#ifdef DEBUG
			rejections_norm_FG++;
			if (rejections_norm_FG%1000 == 0)
				printf("\ttrejected ||FFT(F), FFT(G)|| %lu times\n", rejections_norm_FG);
#endif
				continue;
		}

		// if (G != NULL) {
			/*
			* Compute k = (F*F_adj + G*g_adj) / (f*f_adj + g*g_adj) and
			* check that poly_check_inf_norm(FFT(k), logn, gamma_root).
			* Note that we still have F_fft and G_fft from earlier. 
			*/
			fxr *f_fft, *g_fft, *invnorm, *k;
			f_fft = G_fft + n;
			g_fft = f_fft + n;
			invnorm = g_fft + n;
			k = invnorm + n;
			vect_set(logn, f_fft, f);
			vect_set(logn, g_fft, g);
			vect_FFT(logn, f_fft);
			vect_FFT(logn, g_fft);
			vect_add_muladj_fft(logn, k, F_fft, G_fft, f_fft, g_fft);
			vect_invnorm_fft(logn, invnorm, f_fft, g_fft, 0);
			vect_mul_autoadj_fft(logn, k, invnorm);
			if (vect_check_norm_inf_fft(logn, k, fxr_gamma_root_sqr)) {
#ifdef DEBUG
				rejections_norm_k++;
				if (rejections_norm_k%1000 == 0)
					printf("\trejected ||FFT(k)|| > gamma_root %lu times\n", rejections_norm_k);
#endif
					continue;
			}
		// }

		if (F != NULL) {
			memmove(F, tF, n);
		}
		if (G != NULL) {
			memmove(G, tG, n);
		}
		if (tt32 != tmp) {
			memmove(tmp, tt32, 2 * n);
		}
#ifdef DEBUG
		printf("Signature success:\n\trejection_norm_fg: %lu\n\trejection_alpha_hybrid: %lu\n\trejection_alpha_GPV: %lu\n\trejection__norm_FG: %lu\n\trejection_norm_k: %lu\n\n", rejections_norm_fg, rejections_alpha_hybrid, rejections_alpha_gpv, rejections_norm_FG, rejections_norm_k);
#endif
		return 0;
	}
}

#if 0
/* remove -- this should be in a Falcon implementation, with mod q code,
   not in ntrugen */

/* see ntrugen.h */
int
Falcon_recover_G(unsigned logn,
	int8_t *restrict G,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, void *tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer is the right size; adjust
	 * alignment.
	 */
	if (tmp_len < 1) {
		return -1;
	}
	if (logn < 2 || logn > 10) {
		return -1;
	}
	const ntru_profile *prof;
	switch (logn) {
	case 9:  prof = &SOLVE_Falcon_512;  break;
	case 10: prof = &SOLVE_Falcon_1024; break;
	default: prof = &SOLVE_Falcon_256;  break;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 1) & ~(uintptr_t)1;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint16_t *tt16 = (void *)utmp2;
	if (tmp_len < ((size_t)4 << logn)) {
		return -1;
	}

	/*
	 * f*G - g*F = q, hence G = (q + g*F)/f
	 * We work modulo q, hence G = (g*F)/f mod q.
	 * Keygen should have made sure that f is invertible modulo q.
	 */
	size_t n = (size_t)1 << logn;
	uint32_t r = 0;
	uint16_t *t1 = tt16;
	uint16_t *t2 = t1 + n;
	mq12289_poly_set_small(logn, t1, g);
	mq12289_poly_set_small(logn, t2, F);
	mq12289_NTT(logn, t1);
	mq12289_NTT(logn, t2);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mq12289_tomonty(mq12289_montymul(t1[u], t2[u]));
	}
	mq12289_poly_set_small(logn, t2, f);
	mq12289_NTT(logn, t2);
	for (size_t u = 0; u < n; u ++) {
		/* top bit of r is set if t2[u] is zero (i.e. 12289) */
		r |= 12288 - (uint32_t)t2[u];
		t2[u] = mq12289_div(t1[u], t2[u]);
	}
	mq12289_iNTT(logn, t2);
	int32_t lim = prof->coeff_FG_limit[logn];
	uint8_t *Gt = (uint8_t *)tmp;
	for (size_t u = 0; u < n; u ++) {
		int32_t z = mq12289_snorm(t2[u]);
		/* top bit of r is set if z > +lim or z < -lim */
		r |= (uint32_t)(lim - z);
		r |= (uint32_t)(lim + z);
		Gt[u] = (uint8_t)z;
	}
	if (G != NULL) {
		memmove(G, Gt, n);
	}
	if ((void *)Gt != tmp) {
		memmove(tmp, Gt, n);
	}
	return -(int)(r >> 31);
}

/* see ntrugen.h */
int
Falcon_compute_public(unsigned logn,
	uint16_t *h, const int8_t *restrict f, const int8_t *restrict g,
	void *tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 1) {
		return -1;
	}
	if (logn < 2 || logn > 10) {
		return -1;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 1) & ~(uintptr_t)1;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint16_t *tt16 = (void *)utmp2;
	if (tmp_len < ((size_t)4 << logn)) {
		return -1;
	}

	size_t n = (size_t)1 << logn;
	uint16_t *t1 = tt16;
	uint16_t *t2 = t1 + n;
	mq12289_poly_set_small(logn, t1, f);
	mq12289_poly_set_small(logn, t2, g);
	mq12289_NTT(logn, t1);
	mq12289_NTT(logn, t2);
	uint32_t r = 0;
	for (size_t u = 0; u < n; u ++) {
		r |= 12288 - (uint32_t)t1[u];
		t1[u] = mq12289_div(t2[u], t1[u]);
	}
	mq12289_iNTT(logn, t1);
	for (size_t u = 0; u < n; u ++) {
		t1[u] = mq12289_unorm(t1[u]);
	}
	if (h != NULL) {
		memmove(h, t1, n * sizeof *t1);
	}
	if ((void *)t1 != tmp) {
		memmove(tmp, t1, n * sizeof *t1);
	}
	return -(int)(r >> 31);
}
#endif
