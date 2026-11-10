#include "ng_inner.h"
#include "half_fft_backend.h"
#include "reference_pow2_half_convert.inc"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef TEST_CT
#include "timecop/memcheck.h"
#endif
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "FAIL %d: %s\n", __LINE__, #x);                                        \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
static uint64_t state = UINT64_C(0x1234567812345678);
static uint32_t rnd(void)
{
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return (uint32_t)state & UINT32_C(0x7fffffff);
}
int main(void)
{
    uint32_t *f = malloc(4 * 308 * sizeof *f);
    hw a[4], b[4];
    CHECK(f);
    unsigned tests = 0;
    for (unsigned len = 1; len <= 308; len++)
        for (unsigned neg = 0; neg < 2; neg++) {
#ifdef TEST_CT
            if (len != 1 && len != 4 && len != 5 && len != 8 && len != 64 && len != 128 &&
                len != 207 && len != 308)
                continue;
#endif
            for (unsigned j = 0; j < 4 * len; j++)
                f[j] = rnd();
            for (unsigned j = 0; j < 4; j++)
                f[4 * (len - 1) + j] = (f[4 * (len - 1) + j] & UINT32_C(0x3fffffff)) | (neg << 30);
            const unsigned shifts[] = {0, 1, len / 2, len + 2, len + 3, len + 4, len + 17};
            for (unsigned k = 0; k < sizeof shifts / sizeof *shifts; k++)
                for (unsigned bits = 0; bits < 31; bits++) {
#ifdef TEST_CT
                    if (bits != 0 && bits != 15 && bits != 30)
                        continue;
#endif
                    unsigned value = 31 * shifts[k] + bits;
                    if (value < 13)
                        continue;
                    unsigned scale = value - 13;
                    half_convert_reference(2, b, f, len, scale);
#ifdef TEST_CT
                    CHECK(RUNNING_ON_VALGRIND);
                    VALGRIND_MAKE_MEM_UNDEFINED(f, 4 * len * sizeof *f);
                    VALGRIND_MAKE_MEM_UNDEFINED(&scale, sizeof scale);
#endif
                    half_convert(2, a, f, len, scale);
#ifdef TEST_CT
                    VALGRIND_MAKE_MEM_DEFINED(f, 4 * len * sizeof *f);
                    VALGRIND_MAKE_MEM_DEFINED(&scale, sizeof scale);
                    VALGRIND_MAKE_MEM_DEFINED(a, sizeof a);
#endif
                    CHECK(!memcmp(a, b, sizeof a));
                    tests++;
                }
        }
    free(f);
    printf("PASS %u signed/shift/limb-boundary Q80 conversions\n", tests);
    return 0;
}
