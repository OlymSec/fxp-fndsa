#include "ansa/fxp_fndsa.h"
#include <pthread.h>
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
typedef struct {
    const fxp_fndsa_prepared *key;
    const uint8_t *pk, *expected;
    size_t pl, sl;
    uint8_t seed[64];
} thread_case;
static void *thread_sign(void *arg)
{
    thread_case *c = arg;
    uint8_t sig[1280];
    size_t len;
    for (unsigned i = 0; i < 8; i++) {
        CHECK(!fxp_fndsa_sign(c->key, (const uint8_t *)"abc", 3, c->seed, sig, sizeof sig, &len));
        CHECK(len == c->sl && !memcmp(sig, c->expected, len));
        CHECK(!fxp_fndsa_verify(c->pk, c->pl, (const uint8_t *)"abc", 3, sig, len));
    }
    return NULL;
}
static void test(unsigned n)
{
    uint8_t pk[1793], sk[2369], derived[1793], sig[2479], copy[2479], seed[64] = {0}, rec[1633],
                                                                      out[1633];
    size_t pl = fxp_fndsa_public_key_bytes(n), sl = fxp_fndsa_secret_key_bytes(n), len = 0, got = 0;
    size_t max = fxp_fndsa_mrm_max_bytes(n);
    seed[0] = (uint8_t)(n >> 8);
    CHECK(!fxp_fndsa_keygen(n, seed, pk, sizeof pk, sk, sizeof sk));
    CHECK(pk[0] == (n == 512 ? 9 : 10));
    CHECK(!fxp_fndsa_public_from_secret(sk, sl, derived, sizeof derived));
    CHECK(!memcmp(pk, derived, pl));
    CHECK(fxp_fndsa_keygen(n, seed, pk, pl - 1, sk, sizeof sk) == FXP_FNDSA_SIZE);
    fxp_fndsa_prepared *key = NULL;
    CHECK(!fxp_fndsa_prepare(&key, sk, sl));
    sk[sl - 1] ^= 1;
    fxp_fndsa_prepared *bad = NULL;
    CHECK(fxp_fndsa_prepare(&bad, sk, sl) == FXP_FNDSA_FORMAT && !bad);
    sk[sl - 1] ^= 1;
    CHECK(!fxp_fndsa_sign(key, (const uint8_t *)"abc", 3, seed, sig, sizeof sig, &len));
    CHECK(len == fxp_fndsa_signature_bytes(n));
    memcpy(copy, sig, len);
    CHECK(!fxp_fndsa_verify(pk, pl, (const uint8_t *)"abc", 3, sig, len));
    CHECK(fxp_fndsa_verify(pk, pl, (const uint8_t *)"abc", 2, sig, len) < 0);
    memcpy(derived, pk, pl);
    derived[0] |= 0x10;
    CHECK(fxp_fndsa_verify(derived, pl, (const uint8_t *)"abc", 3, sig, len) == FXP_FNDSA_FORMAT);
    memcpy(derived, pk, pl);
    unsigned h = 1;
    derived[h] = 255;
    derived[h + 1] |= 63;
    CHECK(fxp_fndsa_verify(derived, pl, (const uint8_t *)"abc", 3, sig, len) == FXP_FNDSA_FORMAT);
    thread_case c = {key, pk, copy, pl, len, {0}};
    memcpy(c.seed, seed, 64);
    pthread_t threads[2];
    for (unsigned i = 0; i < 2; i++)
        CHECK(!pthread_create(threads + i, NULL, thread_sign, &c));
    for (unsigned i = 0; i < 2; i++)
        CHECK(!pthread_join(threads[i], NULL));
    CHECK(fxp_fndsa_sign(key, NULL, 1, seed, sig, sizeof sig, &len) == FXP_FNDSA_BAD_ARGUMENT &&
          !len);
    CHECK(!fxp_fndsa_sign(key, NULL, 0, NULL, sig, sizeof sig, &len));
    CHECK(!fxp_fndsa_verify(pk, pl, NULL, 0, sig, len));
    uint8_t digest[64] = {9}, ctx[255] = {1};
    const char *ids[] = {FXP_FNDSA_HASH_ID_RAW, FXP_FNDSA_HASH_ID_SHA512, FXP_FNDSA_HASH_ID_EXTMU};
    for (unsigned i = 0; i < 3; i++) {
        fxp_fndsa_message m = {digest, 64, ctx, sizeof ctx, ids[i]};
        CHECK(!fxp_fndsa_sign_message(key, &m, seed, sig, sizeof sig, &len));
        CHECK(!fxp_fndsa_verify_message(pk, pl, &m, sig, len));
        ctx[0] ^= 1;
        CHECK((fxp_fndsa_verify_message(pk, pl, &m, sig, len) == 0) == (i == 2));
        ctx[0] ^= 1;
    }
    fxp_fndsa_message m = {digest, 64, ctx, 256, NULL};
    CHECK(fxp_fndsa_sign_message(key, &m, seed, sig, sizeof sig, &len) == FXP_FNDSA_BAD_ARGUMENT &&
          !len);
    for (size_t i = 0; i < max; i++)
        rec[i] = (uint8_t)(i * 97 + n);
    rec[max - 1] = 0;
    const size_t lengths[] = {0, 1, 17, max};
    for (unsigned i = 0; i < 4; i++) {
        size_t count = lengths[i];
        CHECK(!fxp_fndsa_sign_mrm(key, rec, count, (const uint8_t *)"abc", 3, seed, sig, sizeof sig,
                                  &len));
        CHECK(len == fxp_fndsa_mrm_signature_bytes(n));
        memset(out, 0xa5, sizeof out);
        got = 42;
        CHECK(!fxp_fndsa_verify_mrm(pk, pl, (const uint8_t *)"abc", 3, sig, len, out, sizeof out,
                                    &got));
        CHECK(got == count && !memcmp(rec, out, count));
        memcpy(copy, sig, len);
        CHECK(!fxp_fndsa_sign_mrm(key, rec, count, (const uint8_t *)"abc", 3, seed, sig, sizeof sig,
                                  &len));
        CHECK(!memcmp(copy, sig, len));
        unsigned old = sig[0];
        sig[0] = 0;
        got = 42;
        memset(out, 0xa5, sizeof out);
        CHECK(fxp_fndsa_verify_mrm(pk, pl, (const uint8_t *)"abc", 3, sig, len, out, sizeof out,
                                   &got) < 0 &&
              !got);
        for (unsigned j = 0; j < sizeof out; j++)
            CHECK(out[j] == 0xa5);
        sig[0] = (uint8_t)old;
        CHECK(fxp_fndsa_verify_mrm(pk, pl, (const uint8_t *)"abc", 2, sig, len, out, sizeof out,
                                   &got) < 0 &&
              !got);
        CHECK(fxp_fndsa_verify_mrm(pk, pl, (const uint8_t *)"abc", 3, sig, len - 1, out, sizeof out,
                                   &got) < 0 &&
              !got);
        if (count)
            CHECK(fxp_fndsa_verify_mrm(pk, pl, (const uint8_t *)"abc", 3, sig, len, out, count - 1,
                                       &got) == FXP_FNDSA_SIZE &&
                  !got);
    }
    CHECK(fxp_fndsa_sign_mrm(key, rec, max + 1, NULL, 0, seed, sig, sizeof sig, &len) ==
              FXP_FNDSA_SIZE &&
          !len);
    CHECK(!fxp_fndsa_sign_mrm(key, NULL, 0, NULL, 0, NULL, sig, sizeof sig, &len));
    CHECK(!fxp_fndsa_verify_mrm(pk, pl, NULL, 0, sig, len, NULL, 0, &got) && !got);
    fxp_fndsa_free(key);
    printf("PASS %u: key binding, ordinary/context/prehash/external-mu, MRM bytes/empty/errors, "
           "concurrent replay\n",
           n);
}
int main(void)
{
    CHECK(!fxp_fndsa_public_key_bytes(769));
    CHECK(!fxp_fndsa_mrm_max_bytes(0));
    fxp_fndsa_free(NULL);
    test(512);
    CHECK(!fxp_fndsa_public_key_bytes(768));
    CHECK(!fxp_fndsa_secret_key_bytes(768));
    test(1024);
    return 0;
}
