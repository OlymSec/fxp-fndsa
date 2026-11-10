#include "ng_inner.h"
#undef mp_NTT
#undef mp_iNTT
#define mp_NTT  Zn(mp_NTT_reference)
#define mp_iNTT Zn(mp_iNTT_reference)
#include "ng_mp31.c"
#undef mp_NTT
#undef mp_iNTT

#if defined(__aarch64__)
#include <arm_neon.h>
static inline uint32x4_t add4(uint32x4_t a, uint32x4_t b, uint32x4_t p)
{
    uint32x4_t d = vsubq_u32(vaddq_u32(a, b), p);
    return vaddq_u32(
        d, vandq_u32(p, vreinterpretq_u32_s32(vshrq_n_s32(vreinterpretq_s32_u32(d), 31))));
}
static inline uint32x4_t sub4(uint32x4_t a, uint32x4_t b, uint32x4_t p)
{
    uint32x4_t d = vsubq_u32(a, b);
    return vaddq_u32(
        d, vandq_u32(p, vreinterpretq_u32_s32(vshrq_n_s32(vreinterpretq_s32_u32(d), 31))));
}
static inline uint32x4_t half4(uint32x4_t a, uint32x4_t p)
{
    uint32x4_t odd = vsubq_u32(vdupq_n_u32(0), vandq_u32(a, vdupq_n_u32(1)));
    return vshrq_n_u32(vaddq_u32(a, vandq_u32(odd, p)), 1);
}
static inline uint32x4_t mont4(uint32x4_t a, uint32x4_t b, uint32x4_t p, uint32x4_t pi)
{
    uint64x2_t lo = vmull_u32(vget_low_u32(a), vget_low_u32(b));
    uint64x2_t hi = vmull_high_u32(a, b);
    uint32x4_t w = vmulq_u32(vmulq_u32(a, b), pi);
    lo = vmlal_u32(lo, vget_low_u32(w), vget_low_u32(p));
    hi = vmlal_high_u32(hi, w, p);
    uint32x4_t q = vcombine_u32(vshrn_n_u64(lo, 32), vshrn_n_u64(hi, 32));
    return sub4(q, p, p);
}
static inline void butterfly4(uint32x4_t *a, uint32x4_t *b, uint32x4_t w, uint32x4_t p,
                              uint32x4_t pi, int inverse)
{
    uint32x4_t x = *a, y = *b;
    if (inverse) {
        *a = half4(add4(x, y, p), p);
        *b = mont4(sub4(x, y, p), w, p, pi);
    } else {
        y = mont4(y, w, p, pi);
        *a = add4(x, y, p);
        *b = sub4(x, y, p);
    }
}
/* IRS's four-lane butterflies and interleaved final stages, with NTRUGEN's
 * own prime, root ordering, and per-stage inverse normalization. */
static void transform(unsigned logn, uint32_t *a, const uint32_t *roots, uint32_t modulus,
                      uint32_t inverse_modulus, int inverse)
{
    unsigned n = 1U << logn;
    uint32x4_t p = vdupq_n_u32(modulus), pi = vdupq_n_u32(inverse_modulus);
    for (unsigned stage = 0; stage < logn; stage++) {
        unsigned half = inverse ? 1U << stage : n >> (stage + 1);
        unsigned m = inverse ? n >> (stage + 1) : 1U << stage;
        if (half == 1) {
            for (unsigned i = 0; i < m; i += 4) {
                uint32x4x2_t v = vld2q_u32(a + 2 * i);
                butterfly4(&v.val[0], &v.val[1], vld1q_u32(roots + m + i), p, pi, inverse);
                vst2q_u32(a + 2 * i, v);
            }
        } else if (half == 2) {
            for (unsigned i = 0; i < m; i += 2) {
                uint32x4_t first = vld1q_u32(a + 4 * i), second = vld1q_u32(a + 4 * i + 4);
                uint32x4_t x = vcombine_u32(vget_low_u32(first), vget_low_u32(second));
                uint32x4_t y = vcombine_u32(vget_high_u32(first), vget_high_u32(second));
                uint32x4_t w = vcombine_u32(vdup_n_u32(roots[m + i]), vdup_n_u32(roots[m + i + 1]));
                butterfly4(&x, &y, w, p, pi, inverse);
                vst1q_u32(a + 4 * i, vcombine_u32(vget_low_u32(x), vget_low_u32(y)));
                vst1q_u32(a + 4 * i + 4, vcombine_u32(vget_high_u32(x), vget_high_u32(y)));
            }
        } else {
            for (unsigned i = 0; i < m; i++) {
                uint32x4_t w = vdupq_n_u32(roots[m + i]);
                for (unsigned j = 2 * half * i; j < 2 * half * i + half; j += 4) {
                    uint32x4_t x = vld1q_u32(a + j), y = vld1q_u32(a + j + half);
                    butterfly4(&x, &y, w, p, pi, inverse);
                    vst1q_u32(a + j, x);
                    vst1q_u32(a + j + half, y);
                }
            }
        }
    }
}
#endif

void Zn(mp_NTT)(unsigned logn, uint32_t *restrict a, const uint32_t *restrict gm, uint32_t p,
                uint32_t pi)
{
#if defined(__aarch64__)
    if (logn >= 3) {
        transform(logn, a, gm, p, pi, 0);
        return;
    }
#endif
    Zn(mp_NTT_reference)(logn, a, gm, p, pi);
}
void Zn(mp_iNTT)(unsigned logn, uint32_t *restrict a, const uint32_t *restrict gm, uint32_t p,
                 uint32_t pi)
{
#if defined(__aarch64__)
    if (logn >= 3) {
        transform(logn, a, gm, p, pi, 1);
        return;
    }
#endif
    Zn(mp_iNTT_reference)(logn, a, gm, p, pi);
}
