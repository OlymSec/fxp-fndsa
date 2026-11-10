#include "protocol.h"
#include "../shared/verify_q12289.h"
#include <string.h>
#ifndef FXP_RADIX_Q2
#define FXP_RADIX_Q2 1
#endif
#ifndef FXP_PACK14_GROUPED
#define FXP_PACK14_GROUPED 1
#endif

void xp_clean(void *ptr, size_t len)
{
    if (!ptr)
        return;
    volatile uint8_t *p = ptr;
#if (defined(__GNUC__) || defined(__clang__)) && UINTPTR_MAX > UINT32_MAX
    typedef uint64_t word __attribute__((__may_alias__));
    while (len && ((uintptr_t)p & 7U)) {
        *p++ = 0;
        len--;
    }
    volatile word *w = (volatile word *)(void *)p;
    while (len >= sizeof *w) {
        *w++ = 0;
        len -= sizeof *w;
    }
    p = (volatile uint8_t *)(void *)w;
#endif
    while (len--)
        *p++ = 0;
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : : "r"(ptr) : "memory");
#endif
}
void xp_hash(uint8_t *out, size_t outlen, const void *in, size_t len)
{
    inner_shake256_context s;
    inner_shake256_init(&s);
    inner_shake256_inject(&s, in, len);
    inner_shake256_flip(&s);
    inner_shake256_extract(&s, out, outlen);
    xp_clean(&s, sizeof s);
}
size_t xp_header(unsigned n)
{
    (void)n;
    return 1;
}
void xp_set_header(uint8_t *p, unsigned n, unsigned kind)
{
    {
        p[0] = (uint8_t)((kind == 1   ? 0
                          : kind == 2 ? 0x50
                          : kind == 3 ? 0x30
                                      : 0x20) +
                         (n == 512 ? 9 : 10));
    }
}
unsigned xp_degree(const uint8_t *p, size_t len, unsigned kind)
{
    if (!p || !len)
        return 0;
    unsigned n = (p[0] & 15U) == 9 ? 512 : (p[0] & 15U) == 10 ? 1024 : 0;
    uint8_t h;
    if (!n)
        return 0;
    xp_set_header(&h, n, kind);
    return p[0] == h ? n : 0;
}
size_t xp_response_bytes(unsigned n)
{
    return n == 512 ? 625 : 1239;
}
unsigned xp_bound(unsigned n)
{
    (void)n;
    return 840;
}
unsigned xp_lambda(unsigned n)
{
    return n == 512 ? 19 : 38;
}
unsigned xp_gamma(unsigned n)
{
    return n == 512 ? 15 : 24;
}
unsigned xp_mu(unsigned n)
{
    return n - xp_lambda(n) - xp_gamma(n);
}
size_t xp_recovery_bytes(unsigned n)
{
    unsigned k = xp_mu(n), bits = 163 * (k / 12) + 27 * ((k % 12) / 2);
    return (bits - 1) / 8;
}
void xp_pack(uint8_t *out, const uint16_t *v, unsigned n, unsigned width)
{
#if FXP_PACK14_GROUPED
    if (width == 14 && !(n & 3U)) {
        for (unsigned i = 0; i < n; i += 4) {
            uint64_t word = (uint64_t)v[i] | ((uint64_t)v[i + 1] << 14) |
                            ((uint64_t)v[i + 2] << 28) | ((uint64_t)v[i + 3] << 42);
            for (unsigned j = 0; j < 7; j++)
                out[7 * (i / 4) + j] = (uint8_t)(word >> (8 * j));
        }
        return;
    }
#endif
    uint64_t acc = 0;
    unsigned bits = 0;
    size_t pos = 0;
    for (unsigned i = 0; i < n; i++) {
        acc |= (uint64_t)v[i] << bits;
        bits += width;
        while (bits >= 8) {
            out[pos++] = (uint8_t)acc;
            acc >>= 8;
            bits -= 8;
        }
    }
    if (bits)
        out[pos] = (uint8_t)acc;
}
int xp_unpack(uint16_t *v, const uint8_t *in, unsigned n, unsigned width, unsigned limit)
{
#if FXP_PACK14_GROUPED
    if (width == 14 && limit == XP_Q && !(n & 3U))
        return ansa_vfy_decode14(v, in, n) ? FXP_FNDSA_FORMAT : 0;
#endif
    uint32_t acc = 0;
    unsigned bits = 0, invalid = 0;
    size_t pos = 0;
    for (unsigned i = 0; i < n; i++) {
        while (bits < width) {
            acc |= (uint32_t)in[pos++] << bits;
            bits += 8;
        }
        v[i] = (uint16_t)(acc & ((1U << width) - 1));
        invalid |= v[i] >= limit;
        acc >>= width;
        bits -= width;
    }
    return (invalid | acc) ? FXP_FNDSA_FORMAT : 0;
}
int xp_response_encode(uint8_t *out, const int16_t *v, unsigned n)
{
    size_t len = xp_response_bytes(n), pos = 0;
    unsigned bits = 0, bound = xp_bound(n);
    uint64_t acc = 0;
    memset(out, 0, len);
    for (unsigned i = 0; i < n; i++) {
        unsigned neg = v[i] < 0, mag = neg ? -(int)v[i] : v[i];
        if (mag > bound)
            return FXP_FNDSA_FORMAT;
        acc |= (uint64_t)(neg | ((mag & 127U) << 1)) << bits;
        bits += 8;
        acc |= UINT64_C(1) << (bits + (mag >> 7));
        bits += (mag >> 7) + 1;
        while (bits >= 8) {
            if (pos == len)
                return FXP_FNDSA_SIZE;
            out[pos++] = (uint8_t)acc;
            acc >>= 8;
            bits -= 8;
        }
    }
    if (bits) {
        if (pos == len)
            return FXP_FNDSA_SIZE;
        out[pos] = (uint8_t)acc;
    }
    return 0;
}
int xp_response_decode(int16_t *v, const uint8_t *in, unsigned n)
{
    return ansa_vfy_decode16(v, in, n, xp_response_bytes(n), xp_bound(n)) ? FXP_FNDSA_FORMAT : 0;
}

int xp_message_mu(uint8_t out[64], const uint8_t pkh[64], const fxp_fndsa_message *m)
{
    if (!m || (!m->data && m->data_len))
        return FXP_FNDSA_BAD_ARGUMENT;
    const uint8_t *id = (const uint8_t *)(m->hash_id ? m->hash_id : "");
    if (id[0] == 0xfe) {
        if (m->data_len != 64)
            return FXP_FNDSA_BAD_ARGUMENT;
        memcpy(out, m->data, 64);
        return 0;
    }
    if ((!m->context && m->context_len) || m->context_len > 255)
        return FXP_FNDSA_BAD_ARGUMENT;
    size_t idlen = 0;
    uint8_t mode = 0, ctxlen = (uint8_t)m->context_len;
    if (id[0] == 6) {
        if (id[1] > 127)
            return FXP_FNDSA_BAD_ARGUMENT;
        idlen = (size_t)id[1] + 2;
        mode = 1;
    } else if (id[0])
        return FXP_FNDSA_BAD_ARGUMENT;
    inner_shake256_context s;
    inner_shake256_init(&s);
    inner_shake256_inject(&s, pkh, 64);
    inner_shake256_inject(&s, &mode, 1);
    inner_shake256_inject(&s, &ctxlen, 1);
    inner_shake256_inject(&s, m->context, m->context_len);
    inner_shake256_inject(&s, id, idlen);
    inner_shake256_inject(&s, m->data, m->data_len);
    inner_shake256_flip(&s);
    inner_shake256_extract(&s, out, 64);
    xp_clean(&s, sizeof s);
    return 0;
}
void xp_coefficients(inner_shake256_context *s, uint16_t *out, unsigned n)
{
    uint8_t b[136];
    while (n) {
        unsigned take = n < 68 ? n : 68;
        inner_shake256_extract(s, b, 2 * take);
        for (unsigned i = 0; i < take; i++) {
            unsigned x = b[2 * i] | ((unsigned)b[2 * i + 1] << 8);
            if (x < 61445) {
                *out++ = (uint16_t)(x % XP_Q);
                n--;
            }
        }
    }
    xp_clean(b, sizeof b);
}
void xp_inject_u16(inner_shake256_context *s, const uint16_t *v, unsigned n)
{
    uint8_t b[136];
    while (n) {
        unsigned take = n < 68 ? n : 68;
        for (unsigned i = 0; i < take; i++) {
            b[2 * i] = (uint8_t)v[i];
            b[2 * i + 1] = (uint8_t)(v[i] >> 8);
        }
        inner_shake256_inject(s, b, 2 * take);
        v += take;
        n -= take;
    }
    xp_clean(b, sizeof b);
}
void xp_target(uint16_t *out, unsigned n, const uint8_t nonce[40], const uint8_t mu[64])
{
    inner_shake256_context s;
    inner_shake256_init(&s);
    inner_shake256_inject(&s, nonce, 40);
    inner_shake256_inject(&s, mu, 64);
    inner_shake256_flip(&s);
    xp_coefficients(&s, out, n);
    xp_clean(&s, sizeof s);
}
void xp_mrm_h1(uint16_t *out, unsigned n, const uint16_t *rho, const uint8_t mu[64],
               const uint8_t *msg, size_t len)
{
    static const char tag[] = "Falcon MRM: H1";
    inner_shake256_context s;
    inner_shake256_init(&s);
    inner_shake256_inject(&s, (const uint8_t *)tag, sizeof tag - 1);
    xp_inject_u16(&s, rho, xp_gamma(n));
    inner_shake256_inject(&s, mu, 64);
    inner_shake256_inject(&s, msg, len);
    inner_shake256_flip(&s);
    xp_coefficients(&s, out, xp_lambda(n));
    xp_clean(&s, sizeof s);
}
void xp_mrm_h2(uint16_t *out, unsigned n, const uint16_t *c0)
{
    static const char tag[] = "Falcon MRM: H2";
    inner_shake256_context s;
    inner_shake256_init(&s);
    inner_shake256_inject(&s, (const uint8_t *)tag, sizeof tag - 1);
    xp_inject_u16(&s, c0, xp_lambda(n));
    inner_shake256_flip(&s);
    xp_coefficients(&s, out, n - xp_lambda(n));
    xp_clean(&s, sizeof s);
}
static unsigned block_bits(unsigned n)
{
    unsigned k = xp_mu(n);
    return 163 * (k / 12) + 27 * ((k % 12) / 2);
}
/* Little-endian bit extraction, bounded by the padded scratch's extra bytes. */
static uint32_t read32(const uint8_t *p, unsigned bit, unsigned count)
{
    p += bit / 8;
    uint64_t x = 0;
    for (unsigned j = 0; j < 5; j++)
        x |= (uint64_t)p[j] << (8 * j);
    return (uint32_t)(x >> (bit % 8)) & (UINT32_MAX >> (32 - count));
}
static void write32(uint8_t *p, unsigned bit, unsigned count, uint32_t x)
{
    p += bit / 8;
    uint64_t v = (uint64_t)(x & (UINT32_MAX >> (32 - count))) << (bit % 8);
    for (unsigned j = 0; j < 5; j++)
        p[j] |= (uint8_t)(v >> (8 * j));
}
int xp_mrm_encode(uint16_t *out, unsigned n, const uint8_t *msg, size_t len, uint8_t *scratch)
{
    if ((!msg && len) || len > xp_recovery_bytes(n))
        return FXP_FNDSA_SIZE;
    memset(scratch, 0, XP_MAX_REC + 8);
    if (len)
        memcpy(scratch, msg, len);
    scratch[len] = 1;
    unsigned k = xp_mu(n), bit = 0, pos = 0;
    for (unsigned block = 0; block < k / 12 + (k % 12) / 2; block++) {
        unsigned width = block < k / 12 ? 163 : 27, digits = width == 163 ? 12 : 2;
        uint32_t limbs[6] = {0};
        for (unsigned j = 0; j < (width + 31) / 32; j++) {
            unsigned take = width - 32 * j;
            if (take > 32)
                take = 32;
            limbs[j] = read32(scratch, bit + 32 * j, take);
        }
        bit += width;
#if FXP_RADIX_Q2
        /* Two base-q digits per pass, as in Antrag's byte codec. After j
         * passes, fewer than width-27*j bits remain because q^2 > 2^27. */
        for (unsigned j = 0; j < digits / 2; j++) {
            uint64_t rem = 0;
            unsigned active = (width - 27 * j + 31) / 32;
            for (unsigned x = active; x-- > 0;) {
                uint64_t value = (rem << 32) | limbs[x];
                uint32_t quotient = (uint32_t)(value / UINT64_C(151019521));
                rem = value - (uint64_t)quotient * UINT64_C(151019521);
                limbs[x] = quotient;
            }
            uint32_t high = (uint32_t)(rem / XP_Q);
            out[pos++] = (uint16_t)(rem - (uint64_t)high * XP_Q);
            out[pos++] = (uint16_t)high;
        }
#else
        for (unsigned j = 0; j < digits; j++) {
            uint64_t rem = 0;
            for (unsigned x = (width + 31) / 32; x-- > 0;) {
                uint64_t v = (rem << 32) | limbs[x];
                limbs[x] = (uint32_t)(v / XP_Q);
                rem = v % XP_Q;
            }
            out[pos++] = (uint16_t)rem;
        }
#endif
        xp_clean(limbs, sizeof limbs);
    }
    return 0;
}
int xp_mrm_decode(uint8_t *out, size_t *len, unsigned n, const uint16_t *v)
{
    unsigned k = xp_mu(n), bit = 0, pos = 0;
    memset(out, 0, XP_MAX_REC + 8);
    *len = 0;
    for (unsigned block = 0; block < k / 12 + (k % 12) / 2; block++) {
        unsigned width = block < k / 12 ? 163 : 27, digits = width == 163 ? 12 : 2;
        uint32_t limbs[6] = {0};
#if FXP_RADIX_Q2
        /* q^(2*d) < 2^(28*d); update only the limbs that can be nonzero. */
        for (unsigned j = digits / 2; j-- > 0;) {
            uint32_t low = v[pos + 2 * j], high = v[pos + 2 * j + 1];
            if (low >= XP_Q || high >= XP_Q)
                return FXP_FNDSA_FORMAT;
            uint64_t carry = low + (uint64_t)high * XP_Q;
            unsigned active = (28 * (digits / 2 - j) + 31) / 32;
            for (unsigned x = 0; x < active; x++) {
                uint64_t value = (uint64_t)limbs[x] * UINT64_C(151019521) + carry;
                limbs[x] = (uint32_t)value;
                carry = value >> 32;
            }
            if (carry)
                return FXP_FNDSA_FORMAT;
        }
#else
        for (unsigned j = digits; j-- > 0;) {
            if (v[pos + j] >= XP_Q)
                return FXP_FNDSA_FORMAT;
            uint64_t carry = v[pos + j];
            for (unsigned x = 0; x < 6; x++) {
                uint64_t z = (uint64_t)limbs[x] * XP_Q + carry;
                limbs[x] = (uint32_t)z;
                carry = z >> 32;
            }
            if (carry)
                return FXP_FNDSA_FORMAT;
        }
#endif
        if (width == 163 ? limbs[5] >= 8 : limbs[0] >= (1U << 27))
            return FXP_FNDSA_FORMAT;
        for (unsigned j = 0; j < (width + 31) / 32; j++) {
            unsigned take = width - 32 * j;
            if (take > 32)
                take = 32;
            write32(out, bit + 32 * j, take, limbs[j]);
        }
        bit += width;
        pos += digits;
    }
    unsigned bytes = (block_bits(n) + 7) / 8;
    while (bytes && !out[bytes - 1])
        bytes--;
    if (!bytes || out[bytes - 1] != 1 || bytes - 1 > xp_recovery_bytes(n))
        return FXP_FNDSA_FORMAT;
    out[bytes - 1] = 0;
    *len = bytes - 1;
    return 0;
}
int xp_norm(unsigned n, const int16_t *a, const int16_t *b)
{
    unsigned bound = xp_bound(n), valid = 1;
    int64_t norm = 0;
    for (unsigned i = 0; i < n; i++) {
        int64_t x = a[i], y = b[i];
        valid &= (x >= -(int64_t)bound) & (x <= bound) & (y >= -(int64_t)bound) & (y <= bound);
        norm += x * x + y * y;
    }
    uint64_t beta = n == 512 ? 34034726 : 70265242;
    return valid && norm >= 0 && (uint64_t)norm <= beta;
}
