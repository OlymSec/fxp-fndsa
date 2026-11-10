#ifndef XP_POW2_WINDOW_H
#define XP_POW2_WINDOW_H
#include "pow2_clean.h"
/* Exact masked-window method from Antrag. Profile cutoffs are unchanged. */
#define HALF_BARREL 316U /* cap + 5, cap <= max_bl_large + 3 = 311 */
static inline uint32_t half_mux32(uint32_t m, uint32_t a, uint32_t b) /* m ? b : a */
{
    return a ^ ((a ^ b) & m);
}
static void half_window(uint32_t w[6], const uint32_t *f, unsigned n, unsigned i, unsigned len,
                        uint32_t sch, uint32_t sign)
{
    uint32_t buf[HALF_BARREL], cap = len + 3;
    if (cap + 5 > HALF_BARREL) /* public: beyond every profile's limb count */
        __builtin_trap();
    unsigned nb = 0;
    while ((1U << nb) <= cap) /* public */
        nb++;
    uint32_t ge = (uint32_t)hw_mask((uint64_t)(sch >= cap));
    sch = half_mux32(ge, sch, cap); /* min(sch, cap) < 2^nb */
    buf[0] = buf[1] = buf[2] = 0;
    const uint32_t *x = f + i;
    for (unsigned j = 0; j < len; j++, x += n)
        buf[3 + j] = *x;
    unsigned b = nb - 1, cnt = (1U << b) + 5, step = 1U << b, t = 0;
    unsigned used = cap > cnt ? cap : cnt;
    uint32_t m = (uint32_t)hw_mask((sch >> b) & 1U);
    for (; t + step < cap; t++)
        buf[t] = half_mux32(m, buf[t], buf[t + step]);
    for (; t < cnt && t < cap; t++)
        buf[t] = half_mux32(m, buf[t], sign);
    for (; t < cnt; t++)
        buf[t] = sign;
    while (b-- > 0) {
        cnt = (1U << b) + 5;
        step = 1U << b;
        m = (uint32_t)hw_mask((sch >> b) & 1U);
        for (t = 0; t < cnt; t++)
            buf[t] = half_mux32(m, buf[t], buf[t + step]);
    }
    for (unsigned k = 0; k < 6; k++)
        w[k] = buf[k];
    xp_pow2_clean(buf, (size_t)used * sizeof *buf);
}
#endif
