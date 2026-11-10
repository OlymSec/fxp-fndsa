#ifndef ANTRAG_HALF_FFT_H
#define ANTRAG_HALF_FFT_H
#include "wide.h"
#include "roots.h"
typedef struct {
    hw inv[512], work[512];
} half_workspace;
static void half_fft(unsigned logn, hw *f, unsigned inverse, unsigned *bad)
{
    unsigned hn = 1U << (logn - 1);
    if (!inverse) {
        unsigned t = hn;
        for (unsigned lm = 1; lm < logn; lm++) {
            unsigned m = 1U << lm, ht = t >> 1;
            for (unsigned i = 0; i < m / 2; i++) {
                hw sr = half_roots[m + i][0], si = half_roots[m + i][1];
                for (unsigned j = i * t; j < i * t + ht; j++) {
                    hw xr = f[j], xi = f[j + hn], yr, yi;
                    hw_cmul(&yr, &yi, f[j + ht], f[j + ht + hn], sr, si, bad);
                    f[j] = hw_add(xr, yr, bad);
                    f[j + hn] = hw_add(xi, yi, bad);
                    f[j + ht] = hw_sub(xr, yr, bad);
                    f[j + ht + hn] = hw_sub(xi, yi, bad);
                }
            }
            t = ht;
        }
    } else {
        unsigned ht = 1;
        for (unsigned lm = logn - 1; lm > 0; lm--) {
            unsigned m = 1U << lm, t = ht << 1;
            for (unsigned i = 0; i < m / 2; i++) {
                hw sr = half_roots[m + i][0], si = -half_roots[m + i][1];
                for (unsigned j = i * t; j < i * t + ht; j++) {
                    hw xr = f[j], xi = f[j + hn], yr = f[j + ht], yi = f[j + ht + hn];
                    f[j] = hw_half(hw_add(xr, yr, bad));
                    f[j + hn] = hw_half(hw_add(xi, yi, bad));
                    hw_cmul(&f[j + ht], &f[j + ht + hn], hw_half(hw_sub(xr, yr, bad)),
                            hw_half(hw_sub(xi, yi, bad)), sr, si, bad);
                }
            }
            ht = t;
        }
    }
}
/* Extract floor(f*2^(80-scale)); scan every limb because scale is secret.
 * Six adjacent base-2^31 digits cover a 128-bit result and a 30-bit offset.
 * The caller checks the scaled signed range before using the result. */
static void half_convert(unsigned logn, hw *d, const uint32_t *f, unsigned len, unsigned scale)
{
    unsigned n = 1U << logn, sch = (scale + 13) / 31, bits = (scale + 13) % 31;
    for (unsigned i = 0; i < n; i++) {
        uint32_t w[6] = {0}, sign = -(f[(len - 1) * n + i] >> 30) >> 1;
        for (unsigned j = 0; j < len; j++) {
            uint32_t x = f[j * n + i];
            for (unsigned k = 0; k < 6; k++)
                w[k] |= x & (uint32_t)hw_mask(j + 3 == sch + k);
        }
        for (unsigned k = 0; k < 6; k++)
            w[k] |= sign & (uint32_t)hw_mask(sch + k >= len + 3);
        hw lo = (hw)w[0] | (hw)w[1] << 31 | (hw)w[2] << 62 | (hw)w[3] << 93 | (hw)w[4] << 124;
        uint64_t hi = ((uint64_t)w[4] >> 4) | (uint64_t)w[5] << 27;
        d[i] = (lo >> bits) | (((hw)hi << 1) << (127 - bits));
    }
}
static int half_prepare(unsigned logn, const uint32_t *f, unsigned len, unsigned scale,
                        half_workspace *ws)
{
    unsigned bad = 0, hn = 1U << (logn - 1);
    half_convert(logn, ws->inv, f, len, scale);
    half_fft(logn, ws->inv, 0, &bad);
    for (unsigned i = 0; i < hn; i++)
        hw_cinv(&ws->inv[i], &ws->inv[i + hn], &bad);
    return !bad;
}
static int half_quotient(unsigned logn, int32_t *k, const uint32_t *F, unsigned len, unsigned scale,
                         half_workspace *ws)
{
    unsigned bad = 0, n = 1U << logn, hn = n / 2;
    /* A specified Q80 input-domain cutoff, before extracting the low word. */
    if (poly_max_bitlength(logn, F, len) > scale + 46)
        return 0;
    half_convert(logn, ws->work, F, len, scale);
    half_fft(logn, ws->work, 0, &bad);
    for (unsigned i = 0; i < hn; i++)
        hw_cmul(&ws->work[i], &ws->work[i + hn], ws->work[i], ws->work[i + hn], ws->inv[i],
                ws->inv[i + hn], &bad);
    half_fft(logn, ws->work, 1, &bad);
    for (unsigned i = 0; i < n; i++) {
        hw x = hw_add(ws->work[i], (hw)1 << 31, &bad);
        hw v = x >> 32;
        bad |= (v >> 31) != (hw_mask128(x >> 127) >> 63);
        k[i] = (int32_t)(uint32_t)v;
    }
    return !bad;
}
#endif
