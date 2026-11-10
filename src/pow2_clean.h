#ifndef XP_POW2_CLEAN_H
#define XP_POW2_CLEAN_H
#include <stdint.h>
#include <stddef.h>
/* Private word wipe, also usable by standalone imported-backend tests. */
static inline void xp_pow2_clean(void *ptr, size_t len)
{
    volatile unsigned char *p = ptr;
    if (!p)
        return;
#if (defined(__GNUC__) || defined(__clang__)) && UINTPTR_MAX > UINT32_MAX
    typedef uint64_t word __attribute__((__may_alias__));
    while (len && ((uintptr_t)p & 7U)) {
        *p++ = 0;
        len--;
    }
    volatile word *w = (volatile word *)(void *)p;
    while (len >= sizeof *w) {
        *w++ = 0;
        len -= sizeof *w;
    }
    p = (volatile unsigned char *)(void *)w;
#endif
    while (len--)
        *p++ = 0;
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : : "r"(ptr) : "memory");
#endif
}
#endif
