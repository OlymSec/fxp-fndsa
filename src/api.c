#include "protocol.h"
#include "pow2_prepare.h"
#include "../shared/verify_q12289.h"
#include "antrag_port.h"
#include "falcon.h"
#include "reconstruction.h"
#include <stdlib.h>
#include <string.h>

struct fxp_fndsa_prepared {
    unsigned n, logn;
    size_t expanded_len, cache_len;
    uint8_t pkh[64], kh[40];
    uint8_t *expanded;
    void *cache;
};
static void clear_free(void *p, size_t n)
{
    if (p) {
        xp_clean(p, n);
        free(p);
    }
}
static unsigned log_degree(unsigned n)
{
    return n == 512 ? 9 : n == 1024 ? 10 : 0;
}
static int supported(unsigned n)
{
    return n == 512 || n == 1024;
}
static size_t raw_secret_bytes(unsigned n)
{
    return FALCON_PRIVKEY_SIZE(log_degree(n));
}

size_t fxp_fndsa_public_key_bytes(unsigned n)
{
    return supported(n) ? xp_header(n) + 14U * n / 8U : 0;
}
size_t fxp_fndsa_secret_key_bytes(unsigned n)
{
    return supported(n) ? raw_secret_bytes(n) + 64 : 0;
}
size_t fxp_fndsa_signature_bytes(unsigned n)
{
    return supported(n) ? xp_header(n) + 40 + xp_response_bytes(n) : 0;
}
size_t fxp_fndsa_mrm_signature_bytes(unsigned n)
{
    return supported(n) ? xp_header(n) + 2 * xp_response_bytes(n) : 0;
}
size_t fxp_fndsa_mrm_max_bytes(unsigned n)
{
    return supported(n) ? xp_recovery_bytes(n) : 0;
}

static unsigned be_read(const uint8_t *p, size_t *pos, unsigned bits)
{
    unsigned x = 0;
    while (bits--) {
        x = (x << 1) | ((p[*pos / 8] >> (7 - *pos % 8)) & 1U);
        (*pos)++;
    }
    return x;
}
static void be_write(uint8_t *p, size_t *pos, unsigned x, unsigned bits)
{
    while (bits--) {
        p[*pos / 8] |= (uint8_t)(((x >> bits) & 1U) << (7 - *pos % 8));
        (*pos)++;
    }
}
/* The pinned sampler retains its original private encoding internally. */
static int secret_convert(uint8_t *out, const uint8_t *in, unsigned n, int to_wire)
{
    size_t hdr = xp_header(n), pos = hdr * 8, offset = hdr;
    unsigned fg = n == 1024 ? 5 : 6;
    uint16_t *v = malloc(n * sizeof *v);
    if (!v)
        return FXP_FNDSA_MEMORY;
    memset(out, 0, raw_secret_bytes(n));
    if (to_wire)
        xp_set_header(out, n, 2);
    else
        out[0] = (uint8_t)(0x50 + log_degree(n));
    int status = 0;
    for (unsigned role = 0; role < 3; role++) {
        unsigned bits = role == 2 ? 8 : fg;
        if (to_wire) {
            for (unsigned i = 0; i < n; i++)
                v[i] = (uint16_t)be_read(in, &pos, bits);
            xp_pack(out + offset, v, n, bits);
        } else {
            status = xp_unpack(v, in + offset, n, bits, 1U << bits);
            if (status)
                break;
            for (unsigned i = 0; i < n; i++) {
                if (v[i] == (1U << (bits - 1))) {
                    status = FXP_FNDSA_FORMAT;
                    break;
                }
                be_write(out, &pos, v[i], bits);
            }
            if (status)
                break;
        }
        offset += (size_t)n * bits / 8;
    }
    clear_free(v, n * sizeof *v);
    return status;
}
static int public_convert(uint8_t *out, const uint8_t *raw, unsigned n)
{
    uint16_t *h = malloc(n * sizeof *h);
    if (!h)
        return FXP_FNDSA_MEMORY;
    size_t pos = xp_header(n) * 8;
    for (unsigned i = 0; i < n; i++)
        h[i] = (uint16_t)be_read(raw, &pos, 14);

    {
#if !FALCON_NTT_PUBLIC_KEY
        Zf(to_ntt_monty)(h, log_degree(n));
#endif
        /* Remove the backend's R=2^16 Montgomery factor. */
        for (unsigned i = 0; i < n; i++)
            h[i] = (uint16_t)((uint32_t)h[i] * 2304U % XP_Q);
    }
    xp_set_header(out, n, 1);
    xp_pack(out + xp_header(n), h, n, 14);
    free(h);
    return 0;
}
static int derive_raw(uint8_t *pk, const uint8_t *sk, unsigned n, void *tmp, size_t size)
{
    return falcon_make_public(pk, fxp_fndsa_public_key_bytes(n), sk, raw_secret_bytes(n),
                                         tmp, size);
}
int fxp_fndsa_keygen(unsigned n, const uint8_t *seed, uint8_t *pk, size_t pkcap, uint8_t *sk,
                     size_t skcap)
{
    if (!supported(n) || !pk || !sk)
        return FXP_FNDSA_BAD_ARGUMENT;
    size_t pl = fxp_fndsa_public_key_bytes(n), sl = fxp_fndsa_secret_key_bytes(n),
           rawlen = raw_secret_bytes(n);
    if (pkcap < pl || skcap < sl)
        return FXP_FNDSA_SIZE;
    memset(pk, 0, pl);
    memset(sk, 0, sl);
    size_t size = FALCON_TMPSIZE_KEYGEN(log_degree(n));
    uint8_t *rawsk = malloc(rawlen), *rawpk = malloc(pl);
    void *tmp = malloc(size);
    int status = FXP_FNDSA_MEMORY;
    if (!rawsk || !rawpk || !tmp)
        goto done;
    {
        uint8_t entropy[64];
        shake256_context rng;
        if (!seed) {
            if (!Zf(get_seed)(entropy, sizeof entropy)) {
                status = FXP_FNDSA_RANDOM;
                xp_clean(entropy, sizeof entropy);
                goto done;
            }
            seed = entropy;
        }
        shake256_init_prng_from_seed(&rng, seed, 64);
        status = falcon_keygen_make(&rng, log_degree(n), rawsk, rawlen, rawpk, pl, tmp, size);
        xp_clean(&rng, sizeof rng);
        xp_clean(entropy, sizeof entropy);
    }
    if (!status)
        status = secret_convert(sk, rawsk, n, 1);
    if (!status)
        status = public_convert(pk, rawpk, n);
    if (!status)
        xp_hash(sk + sl - 64, 64, pk, pl);
done:
    clear_free(rawsk, rawlen);
    free(rawpk);
    clear_free(tmp, size);
    if (status) {
        xp_clean(sk, sl);
        memset(pk, 0, pl);
    }
    return status;
}
static int decode_and_derive(uint8_t *rawsk, uint8_t *pk, const uint8_t *sk, size_t sl, unsigned n)
{
    int status = secret_convert(rawsk, sk, n, 0);
    if (status)
        return status;
    size_t pl = fxp_fndsa_public_key_bytes(n),
           size = FALCON_TMPSIZE_MAKEPUB(log_degree(n));
    uint8_t *rawpk = malloc(pl);
    void *tmp = malloc(size);
    uint8_t hash[64];
    if (!rawpk || !tmp) {
        status = FXP_FNDSA_MEMORY;
        goto done;
    }
    status = derive_raw(rawpk, rawsk, n, tmp, size);
    if (!status)
        status = public_convert(pk, rawpk, n);
    if (!status) {
        xp_hash(hash, 64, pk, pl);
        unsigned different = 0;
        for (unsigned i = 0; i < 64; i++)
            different |= hash[i] ^ sk[sl - 64 + i];
        if (different)
            status = FXP_FNDSA_FORMAT;
    }
done:
    xp_clean(hash, sizeof hash);
    free(rawpk);
    clear_free(tmp, size);
    return status;
}
int fxp_fndsa_public_from_secret(const uint8_t *sk, size_t sl, uint8_t *pk, size_t cap)
{
    if (!sk || !pk)
        return FXP_FNDSA_BAD_ARGUMENT;
    unsigned n = xp_degree(sk, sl, 2);
    if (!n || sl != fxp_fndsa_secret_key_bytes(n))
        return FXP_FNDSA_FORMAT;
    size_t pl = fxp_fndsa_public_key_bytes(n);
    if (cap < pl)
        return FXP_FNDSA_SIZE;
    memset(pk, 0, pl);
    size_t rawlen = raw_secret_bytes(n);
    uint8_t *raw = malloc(rawlen);
    if (!raw)
        return FXP_FNDSA_MEMORY;
    int status = decode_and_derive(raw, pk, sk, sl, n);
    clear_free(raw, rawlen);
    if (status)
        memset(pk, 0, pl);
    return status;
}
void fxp_fndsa_free(fxp_fndsa_prepared *key)
{
    if (!key)
        return;
    clear_free(key->expanded, key->expanded_len);
    clear_free(key->cache, key->cache_len);
    clear_free(key, sizeof *key);
}
int fxp_fndsa_prepare(fxp_fndsa_prepared **out, const uint8_t *sk, size_t sl)
{
    if (out)
        *out = NULL;
    if (!out || !sk)
        return FXP_FNDSA_BAD_ARGUMENT;
    unsigned n = xp_degree(sk, sl, 2);
    if (!n || sl != fxp_fndsa_secret_key_bytes(n))
        return FXP_FNDSA_FORMAT;
    size_t rawlen = raw_secret_bytes(n), pl = fxp_fndsa_public_key_bytes(n);
    fxp_fndsa_prepared *key = calloc(1, sizeof *key);
    uint8_t *raw = malloc(rawlen), *pk = malloc(pl);
    int status = FXP_FNDSA_MEMORY;
    void *tmp = NULL;
    size_t size = 0;
    if (!key || !raw || !pk)
        goto done;
    key->n = n;
    key->logn = log_degree(n);
    {
        status = secret_convert(raw, sk, n, 0);
        if (status)
            goto done;
        key->expanded_len = FALCON_EXPANDEDKEY_SIZE(key->logn);
        key->cache_len = study_prepared_size();
        key->expanded = malloc(key->expanded_len);
        key->cache = calloc(1, key->cache_len);
        size = 6U * n + 16 + FALCON_TMPSIZE_EXPANDPRIV(key->logn);
        tmp = malloc(size);
        if (!key->expanded || !key->cache || !tmp) {
            status = FXP_FNDSA_MEMORY;
            goto done;
        }
        int8_t *basis = tmp;
        uint16_t *h = (uint16_t *)(basis + 4U * n);
        void *work = (void *)(((uintptr_t)(h + n) + 15U) & ~(uintptr_t)15U);
        status = xp_pow2_decode(basis, h, key->logn, raw, rawlen, work);
        if (status)
            goto done;
        xp_set_header(pk, n, 1);
        xp_pack(pk + xp_header(n), h, n, 14);
        xp_hash(key->pkh, 64, pk, pl);
        unsigned different = 0;
        for (unsigned i = 0; i < 64; i++)
            different |= key->pkh[i] ^ sk[sl - 64 + i];
        if (different) {
            status = FXP_FNDSA_FORMAT;
            goto done;
        }
        status = xp_pow2_prepare_core(key->expanded, key->cache, basis, key->logn, work);
    }
    if (status)
        goto done;
    unsigned fg = n == 1024 ? 5 : 6;
    xp_hash(key->kh, 40, sk + xp_header(n), (size_t)2 * n * fg / 8);
done:
    clear_free(raw, rawlen);
    free(pk);
    clear_free(tmp, size);
    if (status)
        fxp_fndsa_free(key);
    else
        *out = key;
    return status;
}

typedef struct {
    uint16_t target[1024], hat[1024], product[1024], rho[24], rec[962], mask[1024], expected[38];
    int16_t response[2][1024];
    ansa_vfy_scratch verification;
    uint8_t pkh[64], mu[64], seed[40], nonce[40], bytes[XP_MAX_REC + 8];
    inner_shake256_context stream;
} workspace;
static int signing_seed(uint8_t out[40], const fxp_fndsa_prepared *key, const uint8_t mu[64],
                        const uint8_t *seed, int mrm, const uint8_t *embedded, size_t len)
{
    /* Distinct from Pornin's arithmetic and from all previous SDK samplers. */
    const uint8_t id_std[8] = {'F', 'X', '6', '4', 'S', '0', '0', '1'},
                  id_mrm[8] = {'F', 'X', '6', '4', 'M', '0', '0', '1'};
    uint8_t entropy[40];
    size_t seedlen = 64;
    inner_shake256_context s;
    if (!seed) {
        if (!Zf(get_seed)(entropy, 40)) {
            xp_clean(entropy, sizeof entropy);
            return FXP_FNDSA_RANDOM;
        }
        seed = entropy;
        seedlen = 40;
    }
    inner_shake256_init(&s);
    inner_shake256_inject(&s, mrm ? id_mrm : id_std, 8);
    inner_shake256_inject(&s, key->kh, 40);
    inner_shake256_inject(&s, mu, 64);
    inner_shake256_inject(&s, seed, seedlen);
    inner_shake256_flip(&s);
    inner_shake256_extract(&s, out, 40);
    if (mrm) {
        inner_shake256_init(&s);
        inner_shake256_inject(&s, out, 40);
        inner_shake256_inject(&s, embedded, len);
        inner_shake256_flip(&s);
        inner_shake256_extract(&s, out, 40);
    }
    xp_clean(&s, sizeof s);
    xp_clean(entropy, sizeof entropy);
    return 0;
}
static int sign_common(const fxp_fndsa_prepared *key, const fxp_fndsa_message *msg,
                       const uint8_t *seed, const uint8_t *embedded, size_t embedded_len, int mrm,
                       uint8_t *sig, size_t cap, size_t *written)
{
    if (written)
        *written = 0;
    if (!key || !msg || !sig || !written || (!embedded && embedded_len))
        return FXP_FNDSA_BAD_ARGUMENT;
    unsigned n = key->n;
    size_t size = mrm ? fxp_fndsa_mrm_signature_bytes(n) : fxp_fndsa_signature_bytes(n);
    if (cap < size)
        return FXP_FNDSA_SIZE;
    memset(sig, 0, size);
    if (mrm && embedded_len > xp_recovery_bytes(n))
        return FXP_FNDSA_SIZE;
    workspace *w = antrag_port_scratch_push(sizeof *w, 0);
    size_t sample_bytes = FALCON_TMPSIZE_SIGNTREE(key->logn);
    void *sample = antrag_port_scratch_push(sample_bytes, 1);
    int status = FXP_FNDSA_MEMORY;
    if (!w || !sample)
        goto done;
    status = xp_message_mu(w->mu, key->pkh, msg);
    if (status)
        goto done;
    status = signing_seed(w->seed, key, w->mu, seed, mrm, embedded, embedded_len);
    if (status)
        goto done;
    if (mrm) {
        status = xp_mrm_encode(w->rec, n, embedded, embedded_len, w->bytes);
        if (status)
            goto done;
    }
    status = FXP_FNDSA_EXHAUSTED;
    for (unsigned attempt = 0; attempt < 27; attempt++) {
        uint8_t ctr = (uint8_t)attempt;
        inner_shake256_init(&w->stream);
        inner_shake256_inject(&w->stream, w->seed, 40);
        inner_shake256_inject(&w->stream, &ctr, 1);
        inner_shake256_flip(&w->stream);
        if (mrm) {
            unsigned lambda = xp_lambda(n), gamma = xp_gamma(n), mu = xp_mu(n);
            xp_coefficients(&w->stream, w->rho, gamma);
            xp_mrm_h1(w->target, n, w->rho, w->mu, embedded, embedded_len);
            xp_mrm_h2(w->mask, n, w->target);
            for (unsigned i = 0; i < gamma; i++)
                w->target[lambda + i] = (uint16_t)((w->mask[i] + w->rho[i]) % XP_Q);
            for (unsigned i = 0; i < mu; i++)
                w->target[lambda + gamma + i] = (uint16_t)((w->mask[gamma + i] + w->rec[i]) % XP_Q);
        } else {
            inner_shake256_extract(&w->stream, w->nonce, 40);
            xp_target(w->target, n, w->nonce, w->mu);
        }
        {
            const fpr *expanded =
                (const fpr *)(((uintptr_t)key->expanded + FPR_LEN) & ~(uintptr_t)(FPR_LEN - 1));
            study_select_prepared(key->cache);
            status = xp_tree_attempt(w->response[0], w->response[1], &w->stream, expanded,
                                     w->target, key->logn, sample)
                         ? FXP_FNDSA_EXHAUSTED
                         : 0;
            study_select_prepared(NULL);
            study_clear_temporary(xp_clean);
        }
        if (status == FXP_FNDSA_EXHAUSTED)
            continue;
        if (status)
            break;
        if (!xp_norm(n, w->response[0], w->response[1])) {
            status = FXP_FNDSA_EXHAUSTED;
            continue;
        }
        size_t off = xp_header(n);
        xp_set_header(sig, n, mrm ? 4 : 3);
        if (mrm) {
            status = xp_response_encode(sig + off, w->response[0], n);
            off += xp_response_bytes(n);
        } else {
            memcpy(sig + off, w->nonce, 40);
            off += 40;
        }
        if (!status)
            status = xp_response_encode(sig + off, w->response[1], n);
        if (status) {
            status = FXP_FNDSA_EXHAUSTED;
            continue;
        }
        *written = size;
        break;
    }
done:
    if (status)
        xp_clean(sig, size);
    antrag_port_scratch_pop(sample, sample_bytes, 1);
    antrag_port_scratch_pop(w, sizeof *w, 1);
    return status;
}
static int verify_common(const uint8_t *pk, size_t pl, const fxp_fndsa_message *msg,
                         const uint8_t *sig, size_t sl, int mrm, uint8_t *embedded, size_t capacity,
                         size_t *recovered)
{
    if (recovered)
        *recovered = 0;
    if (!pk || !msg || !sig || (mrm && (!recovered || (!embedded && capacity))))
        return FXP_FNDSA_BAD_ARGUMENT;
    unsigned n = xp_degree(pk, pl, 1);
    if (!n || pl != fxp_fndsa_public_key_bytes(n) || xp_degree(sig, sl, mrm ? 4 : 3) != n ||
        sl != (mrm ? fxp_fndsa_mrm_signature_bytes(n) : fxp_fndsa_signature_bytes(n)))
        return FXP_FNDSA_FORMAT;
    workspace *w = antrag_port_scratch_push(sizeof *w, 0);
    if (!w)
        return FXP_FNDSA_MEMORY;
    int status = xp_unpack(w->hat, pk + xp_header(n), n, 14, XP_Q);
    if (status)
        goto done;
    xp_hash(w->pkh, 64, pk, pl);
    status = xp_message_mu(w->mu, w->pkh, msg);
    if (status)
        goto done;
    size_t off = xp_header(n);
    if (mrm) {
        status = xp_response_decode(w->response[0], sig + off, n);
        off += xp_response_bytes(n);
        if (status)
            goto done;
    } else {
        memcpy(w->nonce, sig + off, 40);
        off += 40;
    }
    status = xp_response_decode(w->response[1], sig + off, n);
    if (status)
        goto done;
    {
        ansa_vfy_word *a = ansa_vfy_pow2_buffer(&w->verification);
        for (unsigned i = 0; i < n; i++) {
            int32_t v = w->response[1][i];
            a[i] = (ansa_vfy_word)(v + (12289 & -(v < 0)));
        }
        ansa_vfy_pow2_product(w->product, &w->verification, w->hat, n, 0);
    }
    if (mrm)
        for (unsigned i = 0; i < n; i++)
            w->product[i] = (uint16_t)((w->product[i] + XP_Q + w->response[0][i]) % XP_Q);
    if (!mrm) {
        xp_target(w->target, n, w->nonce, w->mu);
        for (unsigned i = 0; i < n; i++) {
            int x = (int)w->target[i] - (int)w->product[i];
            x += XP_Q & -(x < -6144);
            x -= XP_Q & -(x > 6144);
            if (x < -(int)xp_bound(n) || x > (int)xp_bound(n)) {
                status = FXP_FNDSA_BAD_SIGNATURE;
                goto done;
            }
            w->response[0][i] = (int16_t)x;
        }
    }
    uint64_t beta = n == 512 ? 34034726 : 70265242;
    if (!ansa_vfy_norm16(w->response[0], w->response[1], n, beta)) {
        status = FXP_FNDSA_BAD_SIGNATURE;
        goto done;
    }
    if (mrm) {
        unsigned lambda = xp_lambda(n), gamma = xp_gamma(n), mu = xp_mu(n);
        xp_mrm_h2(w->mask, n, w->product);
        for (unsigned i = 0; i < gamma; i++)
            w->rho[i] = (uint16_t)((w->product[lambda + i] + XP_Q - w->mask[i]) % XP_Q);
        for (unsigned i = 0; i < mu; i++)
            w->rec[i] =
                (uint16_t)((w->product[lambda + gamma + i] + XP_Q - w->mask[gamma + i]) % XP_Q);
        size_t len = 0;
        status = xp_mrm_decode(w->bytes, &len, n, w->rec);
        if (status)
            goto done;
        xp_mrm_h1(w->expected, n, w->rho, w->mu, w->bytes, len);
        unsigned diff = 0;
        for (unsigned i = 0; i < lambda; i++)
            diff |= w->expected[i] ^ w->product[i];
        if (diff) {
            status = FXP_FNDSA_BAD_SIGNATURE;
            goto done;
        }
        if (len > capacity) {
            status = FXP_FNDSA_SIZE;
            goto done;
        }
        if (len)
            memcpy(embedded, w->bytes, len);
        *recovered = len;
    }
    status = 0;
done:
    antrag_port_scratch_pop(w, sizeof *w, 0); /* public verification data */
    return status;
}
int fxp_fndsa_sign_message(const fxp_fndsa_prepared *k, const fxp_fndsa_message *m,
                           const uint8_t *seed, uint8_t *sig, size_t cap, size_t *len)
{
    return sign_common(k, m, seed, NULL, 0, 0, sig, cap, len);
}
int fxp_fndsa_verify_message(const uint8_t *pk, size_t pl, const fxp_fndsa_message *m,
                             const uint8_t *sig, size_t sl)
{
    return verify_common(pk, pl, m, sig, sl, 0, NULL, 0, NULL);
}
int fxp_fndsa_sign_mrm_message(const fxp_fndsa_prepared *k, const uint8_t *rec, size_t rl,
                               const fxp_fndsa_message *m, const uint8_t *seed, uint8_t *sig,
                               size_t cap, size_t *len)
{
    return sign_common(k, m, seed, rec, rl, 1, sig, cap, len);
}
int fxp_fndsa_verify_mrm_message(const uint8_t *pk, size_t pl, const fxp_fndsa_message *m,
                                 const uint8_t *sig, size_t sl, uint8_t *rec, size_t cap,
                                 size_t *len)
{
    return verify_common(pk, pl, m, sig, sl, 1, rec, cap, len);
}
int fxp_fndsa_sign(const fxp_fndsa_prepared *k, const uint8_t *data, size_t dl, const uint8_t *seed,
                   uint8_t *sig, size_t cap, size_t *len)
{
    fxp_fndsa_message m = {data, dl, NULL, 0, NULL};
    return fxp_fndsa_sign_message(k, &m, seed, sig, cap, len);
}
int fxp_fndsa_verify(const uint8_t *pk, size_t pl, const uint8_t *data, size_t dl,
                     const uint8_t *sig, size_t sl)
{
    fxp_fndsa_message m = {data, dl, NULL, 0, NULL};
    return fxp_fndsa_verify_message(pk, pl, &m, sig, sl);
}
int fxp_fndsa_sign_mrm(const fxp_fndsa_prepared *k, const uint8_t *rec, size_t rl,
                       const uint8_t *data, size_t dl, const uint8_t *seed, uint8_t *sig,
                       size_t cap, size_t *len)
{
    fxp_fndsa_message m = {data, dl, NULL, 0, NULL};
    return fxp_fndsa_sign_mrm_message(k, rec, rl, &m, seed, sig, cap, len);
}
int fxp_fndsa_verify_mrm(const uint8_t *pk, size_t pl, const uint8_t *data, size_t dl,
                         const uint8_t *sig, size_t sl, uint8_t *rec, size_t cap, size_t *len)
{
    fxp_fndsa_message m = {data, dl, NULL, 0, NULL};
    return fxp_fndsa_verify_mrm_message(pk, pl, &m, sig, sl, rec, cap, len);
}
