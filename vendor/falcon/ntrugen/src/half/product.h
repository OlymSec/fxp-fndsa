#ifndef FX_KEYGEN_WIDE_H
#define FX_KEYGEN_WIDE_H
#include <stdint.h>
typedef unsigned __int128 fx_kg_u128;

/* Complete unsigned 128-by-128 product, least-significant word first. */
static inline void fx_q64_product(uint64_t out[4], fx_kg_u128 a, fx_kg_u128 b)
{
    uint64_t al = (uint64_t)a, ah = (uint64_t)(a >> 64), bl = (uint64_t)b, bh = (uint64_t)(b >> 64);
    fx_kg_u128 ll = (fx_kg_u128)al * bl;
    fx_kg_u128 cross1 = (fx_kg_u128)ah * bl + (uint64_t)(ll >> 64);
    fx_kg_u128 cross2 = (fx_kg_u128)al * bh + (uint64_t)cross1;
    fx_kg_u128 high = (fx_kg_u128)ah * bh + (cross1 >> 64) + (cross2 >> 64);
    out[0] = (uint64_t)ll;
    out[1] = (uint64_t)cross2;
    out[2] = (uint64_t)high;
    out[3] = (uint64_t)(high >> 64);
}

static inline uint64_t fx_kg_mask(uint64_t bit) {
    uint64_t m=0-bit; __asm__("" : "+r"(m)); return m;
}
#endif
