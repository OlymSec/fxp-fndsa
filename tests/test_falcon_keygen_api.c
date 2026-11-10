#include "falcon.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void)
{
    for (unsigned logn = 9; logn <= 10; logn++) {
        size_t sklen = FALCON_PRIVKEY_SIZE(logn), pklen = FALCON_PUBKEY_SIZE(logn);
        size_t scratch = FALCON_TMPSIZE_SIGNDYN(logn);
        if (scratch < FALCON_TMPSIZE_KEYGEN(logn))
            scratch = FALCON_TMPSIZE_KEYGEN(logn);
        uint8_t *sk = malloc(sklen), *pk = malloc(pklen), *pub = malloc(pklen),
                *tmp = malloc(scratch), sig[2496];
        assert(sk && pk && pub && tmp);
        for (unsigned k = 0; k < 8; k++) {
            uint8_t seed[64] = {0};
            seed[0] = (uint8_t)logn;
            seed[1] = (uint8_t)k;
            shake256_context rng;
            shake256_init_prng_from_seed(&rng, seed, sizeof seed);
            assert(!falcon_keygen_make(&rng, logn, sk, sklen, pk, pklen, tmp, scratch));
            assert(!falcon_make_public(pub, pklen, sk, sklen, tmp, scratch));
            assert(!memcmp(pk, pub, pklen));
            size_t len = sizeof sig;
            assert(!falcon_sign_dyn(&rng, sig, &len, FALCON_SIG_PADDED, sk, sklen, seed,
                                    sizeof seed, tmp, scratch));
            assert(!falcon_verify(sig, len, FALCON_SIG_PADDED, pk, pklen, seed, sizeof seed, tmp,
                                  scratch));
            seed[2] ^= 1;
            assert(falcon_verify(sig, len, FALCON_SIG_PADDED, pk, pklen, seed, sizeof seed, tmp,
                                 scratch));
        }
        free(sk);
        free(pk);
        free(pub);
        free(tmp);
        printf("Falcon-%u public keygen/derivation/signing/verification PASS\n", 1U << logn);
    }
}
