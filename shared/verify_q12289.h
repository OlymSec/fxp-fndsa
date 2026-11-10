/* Shared PUBLIC verification arithmetic. No key generation or secret sampler
 * uses this path. Callers supply canonical input residues and heap/TLS scratch.
 * Both libraries use identical butterflies and cubic products. */
#ifndef ANSA_VERIFY_Q12289_H
#define ANSA_VERIFY_Q12289_H
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include "verify_q12289_tables.h"
#if defined(__aarch64__)
#include <arm_neon.h>
#endif
typedef union {
    uint32_t w32[1024];
    uint16_t w16[1024];
} ansa_vfy_scratch;
static inline uint32_t ansa_vfy_add(uint32_t a, uint32_t b)
{
    uint32_t x = a + b;
    return x - (12289U & (0U - (uint32_t)(x >= 12289U)));
}
static inline uint32_t ansa_vfy_sub(uint32_t a, uint32_t b)
{
    return a + (12289U & (0U - (uint32_t)(a < b))) - b;
}
static inline uint32_t ansa_vfy_mul(uint32_t a, uint32_t b)
{
    return a * b % 12289U;
}
static inline uint32_t ansa_vfy_half(uint32_t x)
{
    return (x + (12289U & (0U - (x & 1U)))) >> 1;
}
static inline uint32_t ansa_vfy_mont(uint32_t a, uint32_t b)
{
    uint32_t x = a * b, m = (x * 12287U) & 65535U;
    return ansa_vfy_sub((x + m * 12289U) >> 16, 12289U);
}
static inline int ansa_vfy_decode14(uint16_t *out, const uint8_t *in, unsigned n)
{
    unsigned invalid = 0;
    for (unsigned i = 0; i < n; i += 4) {
        uint64_t word = 0;
        for (unsigned j = 0; j < 7; j++)
            word |= (uint64_t)in[7 * (i / 4) + j] << (8 * j);
        for (unsigned j = 0; j < 4; j++) {
            out[i + j] = (uint16_t)((word >> (14 * j)) & 16383U);
            invalid |= out[i + j] >= 12289;
        }
    }
    return invalid != 0;
}
static inline int32_t ansa_vfy_signed(unsigned magnitude, unsigned negative)
{
    uint32_t mask = 0U - negative;
#if defined(__GNUC__) || defined(__clang__)
    /* Keep the 16-bit output variant from turning this into a random-sign
     * branch. The arithmetic and its 32-bit result are identical. */
    __asm__("" : "+r"(mask));
#endif
    return (int32_t)magnitude - 2 * (int32_t)(magnitude & mask);
}

/* The same buffered LSB-first response decoder, with either output carrier.
 * A nonzero accumulator has its terminator within the available low bits;
 * the magnitude bound rejects a long zero run before any wide shift. */
#define ANSA_VFY_DECODER(name, type)                                                               \
    static inline int name(type *out, const uint8_t *in, unsigned n, size_t len, unsigned bound)   \
    {                                                                                              \
        size_t next = 0;                                                                           \
        uint64_t accumulator = 0;                                                                  \
        unsigned bits = 0;                                                                         \
        for (unsigned i = 0; i < n; i++) {                                                         \
            if (bits < 24) {                                                                       \
                if (len - next >= 4) {                                                             \
                    uint32_t word;                                                                 \
                    memcpy(&word, in + next, sizeof word);                                         \
                    next += 4;                                                                     \
                    word = ansa_vfy_le32(word);                                                    \
                    accumulator |= (uint64_t)word << bits;                                         \
                    bits += 32;                                                                    \
                } else {                                                                           \
                    while (next < len) {                                                           \
                        accumulator |= (uint64_t)in[next++] << bits;                               \
                        bits += 8;                                                                 \
                    }                                                                              \
                }                                                                                  \
            }                                                                                      \
            if (bits < 8)                                                                          \
                return 1;                                                                          \
            unsigned negative = (unsigned)accumulator & 1U;                                        \
            unsigned magnitude = (unsigned)(accumulator >> 1) & 127U;                              \
            accumulator >>= 8;                                                                     \
            bits -= 8;                                                                             \
            if (!accumulator)                                                                      \
                return 1;                                                                          \
            unsigned zeros = (unsigned)__builtin_ctzll(accumulator);                               \
            magnitude += zeros << 7;                                                               \
            if (magnitude > bound || (negative && !magnitude))                                     \
                return 1;                                                                          \
            accumulator >>= zeros + 1;                                                             \
            bits -= zeros + 1;                                                                     \
            out[i] = (type)ansa_vfy_signed(magnitude, negative);                                   \
        }                                                                                          \
        if (accumulator)                                                                           \
            return 1;                                                                              \
        while (next < len)                                                                         \
            if (in[next++])                                                                        \
                return 1;                                                                          \
        return 0;                                                                                  \
    }
static inline uint32_t ansa_vfy_le32(uint32_t x)
{
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return __builtin_bswap32(x);
#else
    return x;
#endif
}
ANSA_VFY_DECODER(ansa_vfy_decode16, int16_t)
ANSA_VFY_DECODER(ansa_vfy_decode32, int32_t)
#undef ANSA_VFY_DECODER

/* The decoders and reconstructed-component check establish |a|,|b|<=2047.
 * Every pair of squares/cross products fits int32; accumulate in int64. */
#define ANSA_VFY_NORM(name, type)                                                                  \
    static inline int name(const type *a, const type *b, unsigned n, uint64_t bound)               \
    {                                                                                              \
        int64_t norm = 0;                                                                          \
        for (unsigned i = 0; i < n; i++)                                                           \
            norm += (int32_t)a[i] * a[i] + (int32_t)b[i] * b[i];                                   \
        return norm >= 0 && (uint64_t)norm <= bound;                                               \
    }
ANSA_VFY_NORM(ansa_vfy_norm16, int16_t)
ANSA_VFY_NORM(ansa_vfy_norm32, int32_t)
#undef ANSA_VFY_NORM

/* All power-of-two inputs are below q; intermediate butterfly values stay
 * canonical. R=2^16, so x+m*q <2^30 and one reduction suffices. */
#if defined(__aarch64__)
typedef uint16_t ansa_vfy_word;
static inline uint16x8_t ansa_vfy_fix8(uint16x8_t x)
{
    return vaddq_u16(x, vandq_u16(vreinterpretq_u16_s16(vshrq_n_s16(vreinterpretq_s16_u16(x), 15)),
                                  vdupq_n_u16(12289)));
}
static inline uint16x8_t ansa_vfy_add8(uint16x8_t a, uint16x8_t b)
{
    return ansa_vfy_fix8(vsubq_u16(vaddq_u16(a, b), vdupq_n_u16(12289)));
}
static inline uint16x8_t ansa_vfy_sub8(uint16x8_t a, uint16x8_t b)
{
    return ansa_vfy_fix8(vsubq_u16(a, b));
}
static inline uint16x8_t ansa_vfy_mont8(uint16x8_t a, uint16x8_t b)
{
    uint16x8_t m = vmulq_u16(vmulq_u16(a, b), vdupq_n_u16(12287));
    uint32x4_t lo =
        vmlal_n_u16(vmull_u16(vget_low_u16(a), vget_low_u16(b)), vget_low_u16(m), 12289);
    uint32x4_t hi = vmlal_n_u16(vmull_high_u16(a, b), vget_high_u16(m), 12289);
    uint16x8_t x = vcombine_u16(vshrn_n_u32(lo, 16), vshrn_n_u32(hi, 16));
    return ansa_vfy_fix8(vsubq_u16(x, vdupq_n_u16(12289)));
}
#else
typedef uint32_t ansa_vfy_word;
#endif
static inline ansa_vfy_word *ansa_vfy_pow2_buffer(ansa_vfy_scratch *s)
{
#if defined(__aarch64__)
    return s->w16;
#else
    return s->w32;
#endif
}
static inline void ansa_vfy_pow2_product(uint16_t *out, ansa_vfy_scratch *s, const uint16_t *hat,
                                         unsigned n, int antrag_order)
{
    ansa_vfy_word *a = ansa_vfy_pow2_buffer(s);
    for (size_t m = 1, t = n; m < n; m <<= 1, t >>= 1) {
        size_t half = t >> 1;
        for (size_t i = 0; i < m; i++) {
            uint32_t w = ansa_vfy_gm[m + i];
            size_t start = i * t;
#if defined(__aarch64__)
            if (half >= 8) {
                uint16x8_t root = vdupq_n_u16(w);
                for (size_t j = start; j < start + half; j += 8) {
                    uint16x8_t x = vld1q_u16(a + j),
                               y = ansa_vfy_mont8(vld1q_u16(a + j + half), root);
                    vst1q_u16(a + j, ansa_vfy_add8(x, y));
                    vst1q_u16(a + j + half, ansa_vfy_sub8(x, y));
                }
                continue;
            }
#endif
            for (size_t j = start; j < start + half; j++) {
                uint32_t x = a[j], y = ansa_vfy_mont(a[j + half], w);
                a[j] = ansa_vfy_add(x, y);
                a[j + half] = ansa_vfy_sub(x, y);
            }
        }
    }
    if (antrag_order) {
        const uint16_t *map = n == 512 ? ansa_vfy_map512 : ansa_vfy_map1024;
        for (unsigned i = 0; i < n; i++)
            a[i] = ansa_vfy_mont(a[i], hat[map[i]]);
    } else
        for (unsigned i = 0; i < n; i++)
            a[i] = ansa_vfy_mont(a[i], hat[i]);
    for (size_t m = n, t = 1; m > 1; m >>= 1, t <<= 1) {
        size_t half = m >> 1;
        for (size_t i = 0; i < half; i++) {
            uint32_t w = ansa_vfy_igm[half + i];
            size_t start = 2 * i * t;
#if defined(__aarch64__)
            if (t >= 8) {
                uint16x8_t root = vdupq_n_u16(w);
                for (size_t j = start; j < start + t; j += 8) {
                    uint16x8_t x = vld1q_u16(a + j), y = vld1q_u16(a + j + t);
                    vst1q_u16(a + j, ansa_vfy_add8(x, y));
                    vst1q_u16(a + j + t, ansa_vfy_mont8(ansa_vfy_sub8(x, y), root));
                }
                continue;
            }
#endif
            for (size_t j = start; j < start + t; j++) {
                uint32_t x = a[j], y = a[j + t];
                a[j] = ansa_vfy_add(x, y);
                a[j + t] = ansa_vfy_mont(ansa_vfy_sub(x, y), w);
            }
        }
    }
    uint32_t scale = n == 512 ? ANSA_VFY_SCALE512 : ANSA_VFY_SCALE1024;
#if defined(__aarch64__)
    uint16x8_t factor = vdupq_n_u16(scale);
    for (unsigned i = 0; i < n; i += 8)
        vst1q_u16(out + i, ansa_vfy_mont8(vld1q_u16(a + i), factor));
#else
    for (unsigned i = 0; i < n; i++)
        out[i] = (uint16_t)ansa_vfy_mont(a[i], scale);
#endif
}
#endif
