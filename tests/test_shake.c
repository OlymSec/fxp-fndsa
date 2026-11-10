#include "falcon.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Separately compiled, unchanged original SHAKE under a private prefix. */
void shake_reference_i_shake256_init(inner_shake256_context *);
void shake_reference_i_shake256_inject(inner_shake256_context *, const uint8_t *, size_t);
void shake_reference_i_shake256_flip(inner_shake256_context *);
void shake_reference_i_shake256_extract(inner_shake256_context *, uint8_t *, size_t);
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "FAIL %d: %s\n", __LINE__, #x);                                        \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
static uint64_t word(void)
{
    static uint64_t s = UINT64_C(0x24315127681024);
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    return s;
}
int main(void)
{
    static const uint8_t empty[32] = {0x46, 0xb9, 0xdd, 0x2b, 0x0b, 0xa8, 0x8d, 0x13,
                                      0x23, 0x3b, 0x3f, 0xeb, 0x74, 0x3e, 0xeb, 0x24,
                                      0x3f, 0xcd, 0x52, 0xea, 0x62, 0xb8, 0x1b, 0x82,
                                      0xb5, 0x0c, 0x27, 0x64, 0x6e, 0xd5, 0x76, 0x2f};
    uint8_t input[4096], a[8193], b[8193];
    inner_shake256_context x, y;
    inner_shake256_init(&x);
    inner_shake256_flip(&x);
    inner_shake256_extract(&x, a, sizeof empty);
    CHECK(!memcmp(a, empty, sizeof empty));
    for (unsigned trial = 0; trial < 200; trial++) {
        size_t len = trial < 10 ? (size_t[]){0, 1, 135, 136, 137, 271, 272, 273, 512, 4096}[trial]
                                : word() % 4096;
        for (size_t i = 0; i < len; i++)
            input[i] = (uint8_t)word();
        inner_shake256_init(&x);
        shake_reference_i_shake256_init(&y);
        size_t p = 0;
        while (p < len) {
            size_t n = 1 + word() % 273;
            if (n > len - p)
                n = len - p;
            inner_shake256_inject(&x, input + p, n);
            shake_reference_i_shake256_inject(&y, input + p, n);
            p += n;
        }
        inner_shake256_flip(&x);
        shake_reference_i_shake256_flip(&y);
        p = 0;
        while (p < sizeof a) {
            size_t n = trial % 2 ? 20 : 1 + word() % 733;
            if (n > sizeof a - p)
                n = sizeof a - p;
            inner_shake256_extract(&x, a + p, n);
            shake_reference_i_shake256_extract(&y, b + p, n);
            p += n;
        }
        CHECK(!memcmp(a, b, sizeof a));
        CHECK(!memcmp(&x, &y, sizeof x));
    }
    puts("PASS SHAKE256 known answer and 200 original-stream/state comparisons, including "
         "fragmented reads and rate boundaries");
    return 0;
}
