#define _GNU_SOURCE
#define _DARWIN_C_SOURCE
#include "ansa/fxp_fndsa.h"
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#if defined(__linux__)
#include <sched.h>
#endif
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
    CHECK(!clock_gettime(CLOCK_MONOTONIC_RAW, &t));
    return (uint64_t)t.tv_sec * 1000000000U + (uint64_t)t.tv_nsec;
}
static void report(unsigned degree, const char *operation, unsigned count, uint64_t elapsed)
{
    printf("%u,%s,%u,%llu,%.3f\n", degree, operation, count, (unsigned long long)elapsed,
           (double)elapsed / count);
}
int main(int argc, char **argv)
{
    unsigned selected = argc > 1 ? (unsigned)strtoul(argv[1], NULL, 10) : 0;
#if defined(__APPLE__)
    CHECK(!pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0));
#elif defined(__linux__)
    const char *cpu = getenv("FXP_BENCH_CPU");
    if (cpu) {
        cpu_set_t mask;
        CPU_ZERO(&mask);
        CPU_SET(atoi(cpu), &mask);
        CHECK(!sched_setaffinity(0, sizeof mask, &mask));
    }
#endif
    setvbuf(stdout, NULL, _IOLBF, 0);
    puts("degree,operation,count,elapsed_ns,ns_per_operation");
    for (unsigned d = 0; d < 2; d++) {
        unsigned degree = (unsigned[]){512, 1024}[d];
        if (selected && degree != selected)
            continue;
        size_t sklen = fxp_fndsa_secret_key_bytes(degree),
               pklen = fxp_fndsa_public_key_bytes(degree), count = fxp_fndsa_mrm_max_bytes(degree);
        uint8_t sk[3][2369], pk[3][1793], derived[1793], sig[2479], seed[64] = {0}, msg[32] = {0};
        uint8_t rec[1633], recovered[1633];
        size_t written, recovered_len;
        uint64_t total = 0;
        for (unsigned i = 0; i < 105; i++) {
            uint64_t start = now();
            CHECK(!fxp_fndsa_keygen(degree, NULL, pk[0], pklen, sk[0], sklen));
            uint64_t elapsed = now() - start;
            if (i >= 5)
                total += elapsed;
            CHECK(!fxp_fndsa_public_from_secret(sk[0], sklen, derived, pklen));
            CHECK(!memcmp(derived, pk[0], pklen));
        }
        report(degree, "keygen", 100, total);
        fxp_fndsa_prepared *keys[3] = {0};
        for (unsigned k = 0; k < 3; k++) {
            seed[0] = (uint8_t)degree;
            seed[1] = (uint8_t)(degree >> 8);
            seed[2] = (uint8_t)k;
            seed[3] = 73;
            CHECK(!fxp_fndsa_keygen(degree, seed, pk[k], pklen, sk[k], sklen));
            CHECK(!fxp_fndsa_prepare(keys + k, sk[k], sklen));
        }
        total = 0;
        for (unsigned k = 0; k < 3; k++)
            for (unsigned i = 0; i < 110; i++) {
                fxp_fndsa_prepared *key = NULL;
                uint64_t start = now();
                CHECK(!fxp_fndsa_prepare(&key, sk[k], sklen));
                uint64_t elapsed = now() - start;
                if (i >= 10)
                    total += elapsed;
                fxp_fndsa_free(key);
            }
        report(degree, "prepare", 300, total);
        for (unsigned mrm = 0; mrm < 2; mrm++) {
            uint64_t signing = 0, verification = 0;
            for (unsigned k = 0; k < 3; k++) {
                for (size_t i = 0; i < count; i++)
                    rec[i] = (uint8_t)(37 * i + k);
                for (unsigned i = 0; i < 1010; i++) {
                    memset(seed, 0, sizeof seed);
                    seed[0] = (uint8_t)d;
                    seed[1] = (uint8_t)k;
                    seed[2] = (uint8_t)mrm;
                    for (unsigned j = 0; j < 4; j++) {
                        seed[4 + j] = (uint8_t)(i >> (8 * j));
                        msg[j] = seed[4 + j];
                    }
                    uint64_t start = now();
                    int status = mrm ? fxp_fndsa_sign_mrm(keys[k], rec, count, msg, sizeof msg,
                                                          seed, sig, sizeof sig, &written)
                                     : fxp_fndsa_sign(keys[k], msg, sizeof msg, seed, sig,
                                                      sizeof sig, &written);
                    uint64_t elapsed = now() - start;
                    CHECK(!status);
                    if (i >= 10)
                        signing += elapsed;
                    start = now();
                    status = mrm ? fxp_fndsa_verify_mrm(pk[k], pklen, msg, sizeof msg, sig, written,
                                                        recovered, sizeof recovered, &recovered_len)
                                 : fxp_fndsa_verify(pk[k], pklen, msg, sizeof msg, sig, written);
                    elapsed = now() - start;
                    CHECK(!status);
                    if (i >= 10)
                        verification += elapsed;
                    if (mrm)
                        CHECK(recovered_len == count && !memcmp(rec, recovered, count));
                }
            }
            report(degree, mrm ? "mrm_sign" : "sign", 3000, signing);
            report(degree, mrm ? "mrm_verify" : "verify", 3000, verification);
        }
        for (unsigned k = 0; k < 3; k++)
            fxp_fndsa_free(keys[k]);
        fprintf(stderr, "degree %u: every key and ordinary/MRM signature verified\n", degree);
    }
    return 0;
}
