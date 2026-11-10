#include "pow2_prepare.h"
#include "falcon.h"
#include "reconstruction.h"
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
void study_trace_reset(void);
uint64_t study_trace_get(void);
unsigned study_trace_count(void);
int main(void)
{
    for (unsigned logn = 9; logn <= 10; logn++) {
        size_t n = 1U << logn, sl = FALCON_PRIVKEY_SIZE(logn), el = FALCON_EXPANDEDKEY_SIZE(logn),
               cl = study_prepared_size();
        size_t wl = FALCON_TMPSIZE_KEYGEN(logn);
        if (wl < FALCON_TMPSIZE_SIGNTREE(logn))
            wl = FALCON_TMPSIZE_SIGNTREE(logn);
        uint8_t *sk = malloc(sl), *expanded = calloc(1, el), *work = malloc(wl);
        uint16_t *target = malloc(n * 2);
        int16_t *out = malloc(n * 4);
        void *cache = malloc(cl);
        CHECK(sk && expanded && work && target && out && cache);
        for (unsigned ki = 0; ki < 3; ki++) {
            uint8_t seed[64] = {0};
            seed[0] = (uint8_t)logn;
            seed[1] = (uint8_t)ki;
            shake256_context rng;
            shake256_init_prng_from_seed(&rng, seed, sizeof seed);
            CHECK(!falcon_keygen_make(&rng, logn, sk, sl, NULL, 0, work, wl));
            CHECK(!falcon_expand_privkey(expanded, el, sk, sl, work, wl));
            CHECK(!study_prepare_expanded(cache, expanded));
            const fpr *tree =
                (const fpr *)(((uintptr_t)expanded + sizeof(fpr)) & ~(uintptr_t)(sizeof(fpr) - 1));
            for (unsigned trial = 0; trial < 64; trial++) {
                for (unsigned j = 0; j < n; j++)
                    target[j] = (uint16_t)((trial * 1777 + ki * 11 + j * 19) % 12289);
                seed[2] = (uint8_t)trial;
                inner_shake256_context stream;
                inner_shake256_init(&stream);
                inner_shake256_inject(&stream, seed, sizeof seed);
                inner_shake256_flip(&stream);
                study_select_prepared(cache);
                study_trace_reset();
                int status = xp_tree_attempt(out, out + n, &stream, tree, target, logn, work);
                CHECK(study_trace_count() == 2 * n);
                uint64_t hash = study_trace_get();
                if (!status)
                    for (unsigned j = 0; j < 2 * n; j++)
                        hash = (hash ^ (uint16_t)out[j]) * UINT64_C(1099511628211);
                printf("TRACE,%zu,%u,%u,%d,%u,%016llx\n", n, ki, trial, status, study_trace_count(),
                       (unsigned long long)hash);
                study_select_prepared(NULL);
                study_clear_temporary(xp_clean);
                xp_clean(&stream, sizeof stream);
            }
            xp_clean(&rng, sizeof rng);
        }
        xp_clean(sk, sl);
        xp_clean(expanded, el);
        xp_clean(work, wl);
        xp_clean(cache, cl);
        xp_clean(out, n * 4);
        free(sk);
        free(expanded);
        free(work);
        free(target);
        free(out);
        free(cache);
    }
    return 0;
}
