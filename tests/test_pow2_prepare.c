#define _POSIX_C_SOURCE 200809L
#include "pow2_prepare.h"
#include "falcon.h"
#include "reconstruction.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "FAIL %d: %s\n", __LINE__, #x);                                        \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
static uint64_t now(void)
{
    struct timespec t;
    CHECK(!clock_gettime(CLOCK_MONOTONIC, &t));
    return (uint64_t)t.tv_sec * UINT64_C(1000000000) + t.tv_nsec;
}
int main(int argc, char **argv)
{
    int bench = argc == 2 && !strcmp(argv[1], "--bench");
    for (unsigned logn = 9; logn <= 10; logn++) {
        size_t n = 1U << logn, sl = FALCON_PRIVKEY_SIZE(logn), el = FALCON_EXPANDEDKEY_SIZE(logn),
               cl = study_prepared_size();
        size_t wl = FALCON_TMPSIZE_KEYGEN(logn);
        if (wl < FALCON_TMPSIZE_EXPANDPRIV(logn))
            wl = FALCON_TMPSIZE_EXPANDPRIV(logn);
        uint8_t *sk = malloc(sl), *expanded = calloc(1, el), *expected = malloc(el),
                *work = malloc(wl);
        int8_t *basis = malloc(4 * n);
        uint16_t *hat = malloc(2 * n);
        void *cache = calloc(1, cl), *reference = calloc(1, cl);
        CHECK(sk && expanded && expected && work && basis && hat && cache && reference);
        for (unsigned ki = 0; ki < 4; ki++) {
            uint8_t seed[64] = {0};
            seed[0] = (uint8_t)logn;
            seed[1] = (uint8_t)ki;
            shake256_context rng;
            shake256_init_prng_from_seed(&rng, seed, sizeof seed);
            CHECK(!falcon_keygen_make(&rng, logn, sk, sl, NULL, 0, work, wl));
            CHECK(!falcon_expand_privkey(expanded, el, sk, sl, work, wl));
            CHECK(!study_prepare_expanded(reference, expanded));
            memcpy(expected, expanded, el);
            CHECK(!xp_pow2_decode(basis, hat, logn, sk, sl, work));
            CHECK(!xp_pow2_prepare_core(expanded, cache, basis, logn, work));
            CHECK(!memcmp(expected, expanded, el));
            CHECK(!memcmp(reference, cache, cl));
            /* The shared NTRU check rejects a mismatched, otherwise small basis. */
            basis[2 * n] ^= 1;
            CHECK(xp_pow2_prepare_core(expanded, cache, basis, logn, work) != 0);
            basis[2 * n] ^= 1;
            CHECK(!xp_pow2_prepare_core(expanded, cache, basis, logn, work));
            CHECK(!memcmp(reference, cache, cl));
            if (bench && ki < 2) {
                for (unsigned rep = 0; rep < 6; rep++) {
                    uint64_t start = now();
                    for (unsigned i = 0; i < 128; i++)
                        CHECK(!xp_pow2_prepare_core(expanded, cache, basis, logn, work));
                    uint64_t ns = now() - start;
                    CHECK(!memcmp(reference, cache, cl));
                    printf("PREP,%zu,%u,%u,128,%llu\n", n, ki, rep, (unsigned long long)ns);
                }
            }
            xp_clean(&rng, sizeof rng);
        }
        xp_clean(basis, 4 * n);
        xp_clean(work, wl);
        xp_clean(sk, sl);
        xp_clean(expanded, el);
        xp_clean(expected, el);
        xp_clean(cache, cl);
        xp_clean(reference, cl);
        free(sk);
        free(expanded);
        free(expected);
        free(work);
        free(basis);
        free(hat);
        free(cache);
        free(reference);
        printf("PASS %zu: exact expanded/cache handoff and mismatched-basis rejection\n", n);
    }
    return 0;
}
