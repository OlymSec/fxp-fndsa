#define _GNU_SOURCE
#include "falcon_mrm.h"
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
enum { KEYS = 3, MESSAGE_BYTES = 32 };
static uint64_t now_ns(void) { struct timespec t; if (clock_gettime(CLOCK_MONOTONIC_RAW, &t)) abort(); return (uint64_t)t.tv_sec * UINT64_C(1000000000) + t.tv_nsec; }
static int run_profile(unsigned logn, unsigned count) {
	size_t pl = FALCON_PRIVKEY_SIZE(logn), ql = FALCON_PUBKEY_SIZE(logn), el = FALCON_EXPANDEDKEY_SIZE(logn), mu = falcon_mrm_message_coeffs(logn), ss = falcon_mrm_signature_size(logn); uint8_t *priv[KEYS], *pub[KEYS], *exp[KEYS], *kg = malloc(FALCON_TMPSIZE_KEYGEN(logn)), *ek = malloc(FALCON_TMPSIZE_EXPANDPRIV(logn)), *st = malloc(FALCON_TMPSIZE_SIGNTREE(logn)), *vt = malloc(FALCON_TMPSIZE_VERIFY(logn)); uint16_t *m0[KEYS]; uint8_t (*msg)[MESSAGE_BYTES] = calloc(count, sizeof *msg), (*seed)[64] = calloc(count, sizeof *seed), (*sig)[2479] = calloc(count, sizeof *sig); uint16_t *out = calloc((size_t)count * mu, sizeof *out); uint64_t ts = 0, tv = 0;
	if (!kg || !ek || !st || !vt || !msg || !seed || !sig || !out) return 1;
	for (unsigned k = 0; k < KEYS; k ++) { shake256_context r; uint8_t ks[64] = {0}; priv[k] = malloc(pl); pub[k] = malloc(ql); exp[k] = malloc(el); m0[k] = calloc(mu, sizeof **m0); if (!priv[k] || !pub[k] || !exp[k] || !m0[k]) return 1; ks[0] = (uint8_t)(logn + 11 * k); shake256_init_prng_from_seed(&r, ks, sizeof ks); if (falcon_keygen_make(&r, logn, priv[k], pl, pub[k], ql, kg, FALCON_TMPSIZE_KEYGEN(logn)) || falcon_expand_privkey(exp[k], el, priv[k], pl, ek, FALCON_TMPSIZE_EXPANDPRIV(logn))) return 1; for (size_t u = 0; u < mu; u ++) m0[k][u] = (uint16_t)((19 * u + k) % 12289); }
	for (unsigned k = 0; k < KEYS; k ++) { for (unsigned i = 0; i < count; i ++) { for (size_t j = 0; j < MESSAGE_BYTES; j ++) msg[i][j] = (uint8_t)(j + 31 * k + i); for (size_t j = 0; j < 64; j ++) seed[i][j] = (uint8_t)(j * 7 + 13 * k + i); } for (unsigned i = 0; i < 5; i ++) { size_t w; if (falcon_mrm_sign_seeded(sig[0], ss, &w, exp[k], pub[k], m0[k], mu, msg[0], MESSAGE_BYTES, seed[0], st, FALCON_TMPSIZE_SIGNTREE(logn))) return 1; } uint64_t begin = now_ns(); for (unsigned i = 0; i < count; i ++) { size_t w; if (falcon_mrm_sign_seeded(sig[i], ss, &w, exp[k], pub[k], m0[k], mu, msg[i], MESSAGE_BYTES, seed[i], st, FALCON_TMPSIZE_SIGNTREE(logn)) || w != ss) return 1; } ts += now_ns() - begin; begin = now_ns(); for (unsigned i = 0; i < count; i ++) if (falcon_mrm_verify(sig[i], ss, pub[k], msg[i], MESSAGE_BYTES, out + (size_t)i * mu, mu, vt, FALCON_TMPSIZE_VERIFY(logn))) return 1; tv += now_ns() - begin; }
	printf("Falcon-%u MRM count=%u sign_ns=%.3f sign_ops_s=%.3f verify_ns=%.3f verify_ops_s=%.3f bytes=%zu recovered=%zu\n", (unsigned)(1U << logn), count, (double)ts / (count * KEYS), (double)(count * KEYS) * 1e9 / ts, (double)tv / (count * KEYS), (double)(count * KEYS) * 1e9 / tv, ss, mu); for (unsigned k = 0; k < KEYS; k ++) { free(priv[k]); free(pub[k]); free(exp[k]); free(m0[k]); } free(kg); free(ek); free(st); free(vt); free(msg); free(seed); free(sig); free(out); return 0;
}
int main(int argc, char **argv) { unsigned count = argc > 1 ? (unsigned)strtoul(argv[1], NULL, 10) : 200; cpu_set_t set; if (!count || count > 100000) return 2; CPU_ZERO(&set); CPU_SET(7, &set); if (sched_setaffinity(0, sizeof set, &set)) return 2; puts("Falcon MRM benchmark: CPU=7; 3 expanded keys/profile; seeded fixed-point extraction"); return run_profile(9, count) || run_profile(10, count); }
