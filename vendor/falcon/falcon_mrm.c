/* Message recovery for the fixed-point Falcon sampler; see MRM.md. */
#include "falcon_mrm.h"
#include <string.h>

#if !FALCON_MRM_EXTRACT_Q64 || FALCON_FPR_MODE != 2
#error "Falcon MRM requires the fixed-point MRM build"
#endif
#if defined(FMLDSA_Q64) && FMLDSA_Q64
#error "Falcon MRM uses q=12289"
#endif

static size_t mrm_lambda(unsigned logn) { return logn == 9 ? 19 : logn == 10 ? 38 : 0; }
static size_t mrm_gamma(unsigned logn) { return logn == 9 ? 15 : logn == 10 ? 24 : 0; }
static size_t pkh_bytes(unsigned logn) { return logn == 9 ? 32 : 64; }
static size_t packed_bytes(size_t count) { return (count * 14 + 7) >> 3; }

size_t falcon_mrm_message_coeffs(unsigned logn)
{
    return mrm_lambda(logn) == 0 ? 0
        : ((size_t)1 << logn) - mrm_lambda(logn) - mrm_gamma(logn);
}
size_t falcon_mrm_signature_size(unsigned logn)
{
    if (logn != 9 && logn != 10) return 0;
    return 1 + 2 * (FALCON_SIG_PADDED_SIZE(logn) - 41);
}
static void clear_bytes(void *data, size_t length)
{
    volatile uint8_t *p = data;
    while (length --) *p ++ = 0;
}
static int public_logn(const void *key, size_t length)
{
    const uint8_t *pk = key;
    if (!pk || !length || (pk[0] != 9 && pk[0] != 10)) return -1;
    if (length != FALCON_PUBKEY_SIZE(pk[0])) return -1;
    return pk[0];
}
/* Inputs have already been validated. Big-endian fourteen-bit fields. */
static void pack_zq(uint8_t *out, const uint16_t *x, size_t count)
{
#if defined(FALCON_INTEGER_RECONSTRUCTION) && FALCON_INTEGER_RECONSTRUCTION
    while(count>=4) {
        uint64_t word=((uint64_t)x[0]<<42)|((uint64_t)x[1]<<28)|((uint64_t)x[2]<<14)|x[3];
        for(unsigned j=0;j<7;j++)out[j]=(uint8_t)(word>>(48-8*j));
        x+=4;out+=7;count-=4;
    }
#endif
    uint32_t acc = 0;
    unsigned bits = 0;
    size_t pos = 0;
    for (size_t u = 0; u < count; u ++) {
        acc = (acc << 14) | x[u];
        bits += 14;
        while (bits >= 8) {
            bits -= 8;
            out[pos ++] = (uint8_t)(acc >> bits);
        }
    }
    if (bits) out[pos] = (uint8_t)(acc << (8 - bits));
}
static void put_u64(uint8_t out[8], uint64_t value)
{
    for (unsigned u = 0; u < 8; u ++, value >>= 8) out[u] = (uint8_t)value;
}
static void hash_init(inner_shake256_context *sc, const char *tag,
    unsigned logn, const uint8_t *pkh)
{
    uint8_t params[3] = {(uint8_t)logn, (uint8_t)mrm_lambda(logn),
        (uint8_t)mrm_gamma(logn)};
    inner_shake256_init(sc);
    inner_shake256_inject(sc, (const uint8_t *)tag, strlen(tag));
    inner_shake256_inject(sc, params, sizeof params);
    inner_shake256_inject(sc, pkh, pkh_bytes(logn));
}
static void hash_point_count(inner_shake256_context *sc, uint16_t *out, size_t count)
{
    while (count) {
        uint8_t block[136];
        size_t take = count < 68 ? count : 68;
        inner_shake256_extract(sc, block, 2 * take);
        for (size_t u = 0; u < take; u ++) {
            unsigned v = ((unsigned)block[2 * u] << 8) | block[2 * u + 1];
            if (v < 61445) {
                *out ++ = (uint16_t)(v % 12289);
                count --;
            }
        }
    }
}
static void public_hash(uint8_t *out, const uint8_t *pk, size_t pk_len, unsigned logn)
{
    inner_shake256_context sc;
    uint8_t p = (uint8_t)logn;
    static const uint8_t tag[] = "FALCON-MRM-PKH-v1";
    inner_shake256_init(&sc);
    inner_shake256_inject(&sc, tag, sizeof tag - 1);
    inner_shake256_inject(&sc, &p, 1);
    inner_shake256_inject(&sc, pk, pk_len);
    inner_shake256_flip(&sc);
    inner_shake256_extract(&sc, out, pkh_bytes(logn));
}
static void hash_h0(uint16_t *out, unsigned logn, const uint8_t *pkh,
    const uint16_t *rho, const uint16_t *m0, const uint8_t *m1, size_t m1_len)
{
    inner_shake256_context sc;
    size_t gamma = mrm_gamma(logn), mu = falcon_mrm_message_coeffs(logn);
    uint8_t rb[42], mb[1684], lb[8];
    pack_zq(rb, rho, gamma);
    pack_zq(mb, m0, mu);
    hash_init(&sc, "FALCON-MRM-H0-v1", logn, pkh);
    inner_shake256_inject(&sc, rb, packed_bytes(gamma));
    inner_shake256_inject(&sc, mb, packed_bytes(mu));
    put_u64(lb, m1_len);
    inner_shake256_inject(&sc, lb, sizeof lb);
    inner_shake256_inject(&sc, m1, m1_len);
    inner_shake256_flip(&sc);
    hash_point_count(&sc, out, mrm_lambda(logn));
}
static void hash_h1(uint16_t *out, unsigned logn, const uint8_t *pkh, const uint16_t *c0)
{
    inner_shake256_context sc;
    size_t lambda = mrm_lambda(logn);
    uint8_t cb[67];
    pack_zq(cb, c0, lambda);
    hash_init(&sc, "FALCON-MRM-H1-v1", logn, pkh);
    inner_shake256_inject(&sc, cb, packed_bytes(lambda));
    inner_shake256_flip(&sc);
    hash_point_count(&sc, out, ((size_t)1 << logn) - lambda);
}
static int derive_attempt(uint16_t *rho, uint8_t root[64], unsigned logn,
    const uint8_t *pkh, const uint8_t seed[64], uint64_t attempt,
    const uint16_t *m0, const uint8_t *m1, size_t m1_len)
{
    inner_shake256_context sc;
    size_t gamma = mrm_gamma(logn), mu = falcon_mrm_message_coeffs(logn);
    uint8_t ctr[8], mb[1684], rb[42], length[8], block[2];
    size_t got = 0;
    pack_zq(mb, m0, mu);
    put_u64(ctr, attempt);
    put_u64(length, m1_len);
    hash_init(&sc, "FALCON-MRM-NONCE-v1", logn, pkh);
    inner_shake256_inject(&sc, seed, 64);
    inner_shake256_inject(&sc, ctr, sizeof ctr);
    inner_shake256_inject(&sc, mb, packed_bytes(mu));
    inner_shake256_inject(&sc, length, sizeof length);
    inner_shake256_inject(&sc, m1, m1_len);
    inner_shake256_flip(&sc);
    for (size_t i = 0; i < 256 && got < gamma; i ++) {
        inner_shake256_extract(&sc, block, 2);
        unsigned v = ((unsigned)block[0] << 8) | block[1];
        if (v < 61445) rho[got ++] = (uint16_t)(v % 12289);
    }
    clear_bytes(&sc, sizeof sc);
    if (got != gamma) return 1;
    pack_zq(rb, rho, gamma);
    hash_init(&sc, "FALCON-MRM-ROOT-v1", logn, pkh);
    inner_shake256_inject(&sc, seed, 64);
    inner_shake256_inject(&sc, ctr, sizeof ctr);
    inner_shake256_inject(&sc, rb, packed_bytes(gamma));
    inner_shake256_inject(&sc, mb, packed_bytes(mu));
    inner_shake256_inject(&sc, length, sizeof length);
    inner_shake256_inject(&sc, m1, m1_len);
    inner_shake256_flip(&sc);
    inner_shake256_extract(&sc, root, 64);
    clear_bytes(&sc, sizeof sc);
    return 0;
}
int falcon_mrm_sign_seeded(void *signature, size_t capacity, size_t *written,
    const void *expanded, size_t expanded_len, const void *public_key, size_t pk_len,
    const uint16_t *m0, size_t m0_count, const uint8_t *m1, size_t m1_len,
    const uint8_t seed[64], void *tmp, size_t tmp_len)
{
    if (written) *written = 0;
    int degree = public_logn(public_key, pk_len);
    if (!signature || !written || !expanded || !seed || !tmp || !m0
        || (m1_len && !m1) || degree < 0) return FALCON_ERR_BADARG;
    unsigned logn = (unsigned)degree;
    size_t lambda = mrm_lambda(logn), gamma = mrm_gamma(logn);
    size_t mu = falcon_mrm_message_coeffs(logn), size = falcon_mrm_signature_size(logn);
    size_t response = (size - 1) / 2;
    if (capacity < size || expanded_len < FALCON_EXPANDEDKEY_SIZE(logn)
        || tmp_len < FALCON_TMPSIZE_SIGNTREE(logn) || m0_count != mu
        || *(const uint8_t *)expanded != logn) return FALCON_ERR_BADARG;
    for (size_t u = 0; u < mu; u ++) if (m0[u] >= 12289) return FALCON_ERR_BADARG;
    uint16_t h[1024];
    if (Zf(modq_decode)(h, logn, (const uint8_t *)public_key + 1, pk_len - 1)
        != pk_len - 1) return FALCON_ERR_FORMAT;
    const fpr *key = (const fpr *)(((uintptr_t)expanded + FPR_LEN) & ~(uintptr_t)(FPR_LEN - 1));
    uint8_t *work = (uint8_t *)(((uintptr_t)tmp + FPR_LEN - 1) & ~(uintptr_t)(FPR_LEN - 1));
    uint8_t pkh[64], root[64];
    uint16_t rho[24], target[1024], y[1024];
    int16_t s0[1024], s1[1024];
    uint8_t *out = signature;
    int status = FALCON_ERR_RANDOM;
    memset(out, 0, size);
    public_hash(pkh, public_key, pk_len, logn);
    for (uint64_t attempt = 0; attempt < 65536; attempt ++) {
        if (derive_attempt(rho, root, logn, pkh, seed, attempt, m0, m1, m1_len)) continue;
        hash_h0(target, logn, pkh, rho, m0, m1, m1_len);
        hash_h1(y, logn, pkh, target);
        for (size_t u = 0; u < gamma; u ++) target[lambda + u] =
            (uint16_t)(((unsigned)y[u] + rho[u]) % 12289U);
        for (size_t u = 0; u < mu; u ++) target[lambda + gamma + u] =
            (uint16_t)(((unsigned)y[gamma + u] + m0[u]) % 12289U);
        inner_shake256_context rng;
        inner_shake256_init(&rng);
        inner_shake256_inject(&rng, root, sizeof root);
        inner_shake256_flip(&rng);
        int rejected = Zf(mrm_sign_tree_attempt)(s0, s1, &rng, key, target, logn, work);
        clear_bytes(&rng, sizeof rng);
        if (rejected) continue;
        size_t v0 = Zf(comp_encode)(out + 1, response, s0, logn);
        size_t v1 = Zf(comp_encode)(out + 1 + response, response, s1, logn);
        if (!v0 || !v1) continue;
        out[0] = (uint8_t)(0x70 + logn);
        memset(out + 1 + v0, 0, response - v0);
        memset(out + 1 + response + v1, 0, response - v1);
        *written = size;
        status = 0;
        break;
    }
    clear_bytes(root, sizeof root);
    clear_bytes(s0, sizeof s0);
    clear_bytes(s1, sizeof s1);
    if (status) memset(out, 0, size);
    return status;
}
int falcon_mrm_sign(void *signature, size_t capacity, size_t *written,
    const void *expanded, size_t expanded_len, const void *pk, size_t pk_len,
    const uint16_t *m0, size_t m0_count, const uint8_t *m1, size_t m1_len,
    void *tmp, size_t tmp_len)
{
    if (written) *written = 0;
    uint8_t seed[64];
    if (!Zf(get_seed)(seed, sizeof seed)) return FALCON_ERR_RANDOM;
    int status = falcon_mrm_sign_seeded(signature, capacity, written,
        expanded, expanded_len, pk, pk_len, m0, m0_count, m1, m1_len, seed, tmp, tmp_len);
    clear_bytes(seed, sizeof seed);
    return status;
}
int falcon_mrm_verify(const void *signature, size_t signature_len,
    const void *public_key, size_t pk_len, const uint8_t *m1, size_t m1_len,
    uint16_t *recoverable, size_t capacity, void *tmp, size_t tmp_len)
{
    int degree = public_logn(public_key, pk_len);
    if (!signature || !recoverable || !tmp || (m1_len && !m1) || degree < 0)
        return FALCON_ERR_BADARG;
    unsigned logn = (unsigned)degree;
    size_t lambda = mrm_lambda(logn), gamma = mrm_gamma(logn);
    size_t mu = falcon_mrm_message_coeffs(logn), size = falcon_mrm_signature_size(logn);
    size_t response = (size - 1) / 2;
    if (signature_len != size || capacity < mu || tmp_len < FALCON_TMPSIZE_VERIFY(logn)
        || *(const uint8_t *)signature != 0x70 + logn) return FALCON_ERR_FORMAT;
    int16_t s0[1024], s1[1024];
    const uint8_t *in = signature;
    size_t v0 = Zf(comp_decode)(s0, logn, in + 1, response);
    size_t v1 = Zf(comp_decode)(s1, logn, in + 1 + response, response);
    if (!v0 || !v1) return FALCON_ERR_FORMAT;
    for (size_t u = v0; u < response; u ++) if (in[1 + u]) return FALCON_ERR_FORMAT;
    for (size_t u = v1; u < response; u ++) if (in[1 + response + u]) return FALCON_ERR_FORMAT;
    if (!Zf(is_short)(s0, s1, logn)) return FALCON_ERR_BADSIG;
    uint16_t h[1024], target[1024], y[1024], rho[24], m0[962], expected[38];
    uint8_t pkh[64];
    if (Zf(modq_decode)(h, logn, (const uint8_t *)public_key + 1, pk_len - 1)
        != pk_len - 1) return FALCON_ERR_FORMAT;
    Zf(to_ntt_monty)(h, logn);
    Zf(mrm_recover_target)(target, s0, s1, h, logn);
    public_hash(pkh, public_key, pk_len, logn);
    hash_h1(y, logn, pkh, target);
    for (size_t u = 0; u < gamma; u ++) rho[u] =
        (uint16_t)((target[lambda + u] + 12289U - y[u]) % 12289U);
    for (size_t u = 0; u < mu; u ++) m0[u] =
        (uint16_t)((target[lambda + gamma + u] + 12289U - y[gamma + u]) % 12289U);
    hash_h0(expected, logn, pkh, rho, m0, m1, m1_len);
    if (memcmp(expected, target, lambda * sizeof *target)) return FALCON_ERR_BADSIG;
    memcpy(recoverable, m0, mu * sizeof *m0);
    return 0;
}
