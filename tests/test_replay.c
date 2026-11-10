/* The new wrapper continues an attempt stream. The pinned sampler seeded
 * from SHAKE(seed56) must equal the continuing-stream entry point when given
 * that same SHAKE(seed56) state, including a rejected candidate. */
#include "falcon.h"
#include "reconstruction.h"
#include "../src/protocol.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "FAIL %d: %s\n", __LINE__, #x);                                        \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
int main(void)
{
    for (unsigned logn = 9; logn <= 10; logn++) {
        size_t n = 1u << logn, sl = FALCON_PRIVKEY_SIZE(logn), size = FALCON_TMPSIZE_SIGNTREE(logn);
        if (size < FALCON_TMPSIZE_KEYGEN(logn))
            size = FALCON_TMPSIZE_KEYGEN(logn);
        if (size < FALCON_TMPSIZE_EXPANDPRIV(logn))
            size = FALCON_TMPSIZE_EXPANDPRIV(logn);
        uint8_t *sk = malloc(sl), *expanded = malloc(FALCON_EXPANDEDKEY_SIZE(logn)),
                *work = malloc(size);
        void *cache = calloc(1, study_prepared_size());
        CHECK(sk && expanded && work && cache);
        uint8_t seed[64] = {0};
        seed[0] = (uint8_t)logn;
        shake256_context rng;
        shake256_init_prng_from_seed(&rng, seed, 64);
        CHECK(!falcon_keygen_make(&rng, logn, sk, sl, NULL, 0, work, size));
        CHECK(!falcon_expand_privkey(expanded, FALCON_EXPANDEDKEY_SIZE(logn), sk, sl, work, size));
        CHECK(!study_prepare_expanded(cache, expanded));
        const fpr *tree =
            (const fpr *)(((uintptr_t)expanded + FPR_LEN) & ~(uintptr_t)(FPR_LEN - 1));
        uint16_t target[1024];
        int16_t a[2][1024], b[2][1024];
        for (unsigned trial = 0; trial < 32; trial++) {
            for (unsigned i = 0; i < n; i++)
                target[i] = (uint16_t)((17 * i + trial) % 12289);
            seed[1] = (uint8_t)trial;
            inner_shake256_context old, new;
            uint8_t material[56];
            inner_shake256_init(&old);
            inner_shake256_inject(&old, seed, 64);
            inner_shake256_flip(&old);
            new = old;
            inner_shake256_extract(&new, material, 56);
            inner_shake256_init(&new);
            inner_shake256_inject(&new, material, 56);
            inner_shake256_flip(&new);
            study_select_prepared(cache);
            int x = Zf(mrm_sign_tree_attempt)(a[0], a[1], &old, tree, target, logn, work);
            int y = xp_tree_attempt(b[0], b[1], &new, tree, target, logn, work);
            study_select_prepared(NULL);
            CHECK(x == y);
            if (!x)
                CHECK(!memcmp(a[0], b[0], n * 2) && !memcmp(a[1], b[1], n * 2));
        }
        study_clear_temporary(xp_clean);
        xp_clean(work, size);
        xp_clean(sk, sl);
        xp_clean(expanded, FALCON_EXPANDEDKEY_SIZE(logn));
        xp_clean(cache, study_prepared_size());
        free(work);
        free(sk);
        free(expanded);
        free(cache);
        printf("PASS %zu: 32 single-attempt stream/arithmetic replay comparisons\n", n);
    }
    return 0;
}
