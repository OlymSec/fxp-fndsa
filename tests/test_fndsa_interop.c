#include "ansa/fxp_fndsa.h"
#include "fndsa.h"
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
#ifndef CORRECTED_MRM
#define CORRECTED_MRM 0
#endif
static void check_key(unsigned n, const uint8_t *sk, const uint8_t *pk)
{
    unsigned logn = n == 512 ? 9 : 10;
    size_t sl = FNDSA_SIGN_KEY_SIZE(logn), pl = FNDSA_VRFY_KEY_SIZE(logn);
    CHECK(sl == fxp_fndsa_secret_key_bytes(n) && pl == fxp_fndsa_public_key_bytes(n));
    uint8_t derived[1793], a[2479], b[2479], seed[64] = {7}, digest[64] = {1}, ctx[255] = {3},
                                             rec[1633], out[1633];
    uint8_t *tmp = malloc(FNDSA_SIGN_TMP_SIZE(logn));
    CHECK(tmp);
    fxp_fndsa_prepared *key = NULL;
    CHECK(!fxp_fndsa_public_from_secret(sk, sl, derived, sizeof derived));
    CHECK(!memcmp(pk, derived, pl));
    CHECK(!fxp_fndsa_prepare(&key, sk, sl));
    const char *ids[] = {FXP_FNDSA_HASH_ID_RAW, FXP_FNDSA_HASH_ID_SHA256, FXP_FNDSA_HASH_ID_SHA512,
                         FXP_FNDSA_HASH_ID_EXTMU};
    for (unsigned m = 0; m < 4; m++) {
        fxp_fndsa_message msg = {digest, m == 1 ? 32 : 64, ctx, sizeof ctx, ids[m]};
        size_t al = 0;
        CHECK(!fxp_fndsa_sign_message(key, &msg, seed, a, sizeof a, &al));
        CHECK(fndsa_verify(a, al, pk, pl, ctx, sizeof ctx, ids[m], digest, msg.data_len) == 1);
        size_t bl =
            fndsa_sign_seeded_temp(sk, sl, ctx, sizeof ctx, ids[m], digest, msg.data_len, seed,
                                   sizeof seed, b, sizeof b, tmp, FNDSA_SIGN_TMP_SIZE(logn));
        CHECK(bl == al);
        CHECK(!fxp_fndsa_verify_message(pk, pl, &msg, b, bl));
        CHECK(memcmp(a, b, al) != 0);
        if (m == 0) {
            void *expanded = malloc(FNDSA_EXP_KEY_SIZE(logn));
            void *etmp = malloc(FNDSA_SIGNEXP_TMP_SIZE(logn));
            CHECK(expanded && etmp);
            void *ek = fndsa_expand_signkey(sk, sl, expanded, FNDSA_EXP_KEY_SIZE(logn));
            CHECK(ek);
            bl = fndsa_signexp_seeded_temp(ek, ctx, sizeof ctx, ids[m], digest, msg.data_len, seed,
                                           sizeof seed, b, sizeof b, etmp,
                                           FNDSA_SIGNEXP_TMP_SIZE(logn));
            CHECK(bl == al && !fxp_fndsa_verify_message(pk, pl, &msg, b, bl));
            free(expanded);
            free(etmp);
        }
    }
    for (unsigned i = 0; i < sizeof rec; i++)
        rec[i] = (uint8_t)(i * 17 + 3);
    size_t sizes[] = {0, 1, 17, fxp_fndsa_mrm_max_bytes(n)};
    CHECK(sizes[3] == FNDSA_MAX_MRM_SIZE(logn));
    for (unsigned i = 0; i < 4; i++) {
        size_t count = sizes[i], al = 0, recovered = 99;
        fxp_fndsa_message msg = {digest, 64, ctx, sizeof ctx, FXP_FNDSA_HASH_ID_SHA512};
        CHECK(!fxp_fndsa_sign_mrm_message(key, rec, count, &msg, seed, a, sizeof a, &al));
        int length = fndsa_verify_mrm_temp(a, al, pk, pl, ctx, sizeof ctx, msg.hash_id, digest, 64,
                                           out, sizeof out, tmp, FNDSA_SIGN_TMP_SIZE(logn));
        if (CORRECTED_MRM)
            CHECK(length == (int)count && !memcmp(rec, out, count));
        else
            CHECK(length == -1 || length == 0); /* Upstream failure can alias empty success. */
        size_t bl = fndsa_sign_mrm_seeded_temp(sk, sl, rec, count, ctx, sizeof ctx, msg.hash_id,
                                               digest, 64, seed, sizeof seed, b, sizeof b, tmp,
                                               FNDSA_SIGN_TMP_SIZE(logn));
        CHECK(bl == al);
        int status = fxp_fndsa_verify_mrm_message(pk, pl, &msg, b, bl, out, sizeof out, &recovered);
        if (CORRECTED_MRM)
            CHECK(!status && recovered == count && !memcmp(rec, out, count));
        else
            CHECK(status < 0 && !recovered);
    }
    fxp_fndsa_free(key);
    free(tmp);
}
int main(void)
{
    for (unsigned n = 512; n <= 1024; n *= 2) {
        unsigned logn = n == 512 ? 9 : 10;
        uint8_t sk[2369], pk[1793], seed[64] = {0};
        seed[0] = (uint8_t)logn;
        CHECK(!fxp_fndsa_keygen(n, seed, pk, sizeof pk, sk, sizeof sk));
        check_key(n, sk, pk);
        fndsa_keygen_seeded(logn, seed, sizeof seed, sk, pk);
        check_key(n, sk, pk);
        printf("PASS %u: both key generators; bidirectional ordinary, context, prehash, "
               "external-mu, expanded signing; MRM %s\n",
               n,
               CORRECTED_MRM ? "specified full-tag interoperability"
                             : "uncorrected tag mismatch detected");
    }
    return 0;
}
