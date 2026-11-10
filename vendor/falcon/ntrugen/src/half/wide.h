/* Fixed Q80 arithmetic for the optional F-only completion policy.
 * Values are two's-complement 128-bit words; products have 256-bit carriers.
 * A range failure is a declared rejection of this completion policy. */
#ifndef ANTRAG_HALF_WIDE_H
#define ANTRAG_HALF_WIDE_H
#include <stdint.h>
#include "product.h"
#include "division.h"
typedef unsigned __int128 hw;
#define HW_FRAC 80
static inline uint64_t hw_mask(uint64_t x)
{
    return fx_kg_mask(x);
}
static inline hw hw_mask128(unsigned x)
{
    uint64_t m = hw_mask(x);
    return ((hw)m << 64) | m;
}
static inline hw hw_abs(hw x)
{
    hw m = hw_mask128(x >> 127);
    return (x ^ m) - m;
}
static inline hw hw_add(hw x, hw y, unsigned *bad)
{
    hw z = x + y;
    *bad |= (unsigned)((~(x ^ y) & (x ^ z)) >> 127);
    return z;
}
static inline hw hw_sub(hw x, hw y, unsigned *bad)
{
    hw z = x - y;
    *bad |= (unsigned)(((x ^ y) & (x ^ z)) >> 127);
    return z;
}
static inline hw hw_half(hw x)
{
    /* ceil(x/2), without overflowing x+1. */
    return (x >> 1) + (x & 1) + (x & ((hw)1 << 127));
}
static inline hw hw_mul(hw a, hw b, unsigned *bad)
{
    uint64_t p[4];
    fx_q64_product(p, hw_abs(a), hw_abs(b));
    uint64_t m = hw_mask((a ^ b) >> 127), carry = m & 1;
    for (unsigned i = 0; i < 4; i++) {
        hw t = (hw)(p[i] ^ m) + carry;
        p[i] = (uint64_t)t;
        carry = t >> 64;
    }
    /* floor((a*b + 2^79)/2^80), ties towards positive infinity. */
    hw t = (hw)p[1] + (UINT64_C(1) << 15);
    p[1] = t;
    carry = t >> 64;
    t = (hw)p[2] + carry;
    p[2] = t;
    carry = t >> 64;
    p[3] += carry;
    hw z = ((hw)(p[2] >> 16 | p[3] << 48) << 64) | (p[1] >> 16 | p[2] << 48);
    *bad |= (unsigned)((p[3] >> 16) != (hw_mask(z >> 127) >> 16));
    return z;
}
static inline void hw_cmul(hw *zr, hw *zi, hw ar, hw ai, hw br, hw bi, unsigned *bad)
{
    hw rr = hw_mul(ar, br, bad), ii = hw_mul(ai, bi, bad);
    hw ri = hw_mul(ar, bi, bad), ir = hw_mul(ai, br, bad);
    *zr = hw_sub(rr, ii, bad);
    *zi = hw_add(ri, ir, bad);
}
/* Shift a fixed-size array left by 0..63 bits, without shift-by-64. */
static inline void hw_lshift(uint64_t *x, unsigned n, unsigned bits)
{
    for (unsigned j = n; j-- > 1;)
        x[j] = (x[j] << bits) | ((x[j - 1] >> 1) >> (63 - bits));
    x[0] <<= bits;
}
/* floor(|a|*2^112/d), sign restored. d is an exact squared magnitude. The Q80 inverse is scaled by 2^-48.
 * Fixed 6-by-4-limb division, with masked normalization. */
static inline hw hw_inv_component(hw a, const uint64_t d_in[4], unsigned *bad)
{
    hw x = hw_abs(a);
    uint64_t d[4] = {d_in[0], d_in[1], d_in[2], d_in[3]};
    uint64_t num[6] = {0, (uint64_t)x << 48, (uint64_t)(x >> 16), (uint64_t)(x >> 80), 0, 0};
    unsigned zero = (d[0] | d[1] | d[2] | d[3]) == 0;
    *bad |= zero;
    d[0] |= zero;
    for (unsigned i = 0; i < 3; i++) {
        uint64_t m = hw_mask(d[3] == 0);
        for (unsigned j = 4; j-- > 1;)
            d[j] = (d[j] & ~m) | (d[j - 1] & m);
        d[0] &= ~m;
        for (unsigned j = 6; j-- > 1;)
            num[j] = (num[j] & ~m) | (num[j - 1] & m);
        num[0] &= ~m;
    }
    unsigned bits = __builtin_clzll(d[3]);
    hw_lshift(d, 4, bits);
    hw_lshift(num, 6, bits);
    uint64_t q[2];
    uint64_t high = Zn(half_div_qr)(q, num, 6, d, 4);
    *bad |= (high | (q[1] >> 63)) != 0;
    hw z = ((hw)q[1] << 64) | q[0], m = hw_mask128(a >> 127);
    return (z ^ m) - m;
}
static inline void hw_cinv(hw *ar, hw *ai, unsigned *bad)
{
    uint64_t rr[4], ii[4], d[4], carry = 0;
    fx_q64_product(rr, hw_abs(*ar), hw_abs(*ar));
    fx_q64_product(ii, hw_abs(*ai), hw_abs(*ai));
    for (unsigned j = 0; j < 4; j++) {
        hw t = (hw)rr[j] + ii[j] + carry;
        d[j] = t;
        carry = t >> 64;
    }
    hw re = hw_inv_component(*ar, d, bad), im = hw_inv_component(*ai, d, bad);
    *ar = re;
    *ai = -im;
}
#endif
