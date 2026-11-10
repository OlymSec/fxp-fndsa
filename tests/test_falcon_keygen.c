#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntrugen.h"

static unsigned rng_calls;
static void fill(void *ctx, void *out, size_t length)
{
    uint64_t *state = ctx;
    uint8_t *p = out;
    rng_calls++;
    for (size_t i = 0; i < length; i++) {
        uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        p[i] = (uint8_t)(z ^ (z >> 31));
    }
}
static void equation(unsigned n, const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G)
{
    int64_t a[2048] = {0};
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < n; j++)
            a[i + j] += (int64_t)f[i] * G[j] - (int64_t)g[i] * F[j];
    for (unsigned i = 0; i < n; i++)
        assert(a[i] - a[i + n] == (i ? 0 : 12289));
}
int main(void)
{
    for (unsigned logn = 9; logn <= 10; logn++) {
        unsigned n = 1U << logn;
        size_t length = NTRUGEN_FALCON_TMP_SIZE(logn);
        uint8_t *memory = malloc(length + 32);
        assert(memory);
        int8_t f[1024], g[1024], F[1024], G[1024], expected[4096];
        uint64_t expected_state = 0;
        for (unsigned pointers = 0; pointers < 4; pointers++)
            for (unsigned offset = 0; offset < 8; offset++) {
                uint64_t state = UINT64_C(123456789) + logn;
                rng_calls = 0;
                memset(memory, 0xA5, length + 32);
                assert(!ntrugen_Falcon_keygen(logn, f, g, (pointers & 1) ? F : NULL,
                                              (pointers & 2) ? G : NULL, fill, &state,
                                              memory + offset, length));
                assert(rng_calls > 0);
                for (unsigned i = 0; i < offset; i++)
                    assert(memory[i] == 0xA5);
                for (size_t i = offset + length; i < length + 32; i++)
                    assert(memory[i] == 0xA5);
                if (pointers & 1)
                    assert(!memcmp(F, memory + offset, n));
                if (pointers & 2)
                    assert(!memcmp(G, memory + offset + n, n));
                equation(n, f, g, (int8_t *)(memory + offset), (int8_t *)(memory + offset + n));
                if (!pointers && !offset) {
                    memcpy(expected, f, n);
                    memcpy(expected + n, g, n);
                    memcpy(expected + 2 * n, memory + offset, 2 * n);
                    expected_state = state;
                } else {
                    assert(state == expected_state);
                    assert(!memcmp(expected, f, n));
                    assert(!memcmp(expected + n, g, n));
                    assert(!memcmp(expected + 2 * n, memory + offset, 2 * n));
                }
            }
        for (unsigned offset = 0; offset < 8; offset++) {
            uint64_t state = 1;
            rng_calls = 0;
            uint8_t *tmp = memory + offset;
            size_t padding = (8 - ((uintptr_t)tmp & 7)) & 7;
            assert(ntrugen_Falcon_keygen(logn, f, g, F, G, fill, &state, tmp,
                                         50 * n + padding - 1) == -1);
            assert(!rng_calls);
        }
        free(memory);
        printf("N%u: all F/G pointer combinations, alignments, exact equations, raw outputs and "
               "short-buffer rejection PASS\n",
               n);
    }
}
