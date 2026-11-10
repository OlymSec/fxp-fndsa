#include "../src/protocol.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "FAIL %d %s\n", __LINE__, #x);                                         \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
static uint32_t state = 1;
static uint32_t word(void)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}
static int bit(const uint8_t *p, size_t at)
{
    return (p[at / 8] >> (at % 8)) & 1;
}
static int reference(int16_t *out, const uint8_t *p, unsigned n)
{
    size_t at = 0, bits = xp_response_bytes(n) * 8;
    for (unsigned i = 0; i < n; i++) {
        if (at + 9 > bits)
            return -3;
        unsigned neg = bit(p, at++), mag = 0;
        for (unsigned j = 0; j < 7; j++)
            mag |= (unsigned)bit(p, at++) << j;
        for (;;) {
            if (at >= bits)
                return -3;
            if (bit(p, at++))
                break;
            mag += 128;
            if (mag > xp_bound(n))
                return -3;
        }
        if (mag > xp_bound(n) || (neg && !mag))
            return -3;
        out[i] = (int16_t)(neg ? -(int)mag : (int)mag);
    }
    while (at < bits)
        if (bit(p, at++))
            return -3;
    return 0;
}
int main(void)
{
    for (unsigned n = 512; n <= 1024; n *= 2) {
        int16_t values[1024], a[1024], b[1024];
        uint8_t encoded[1245], changed[1245], packed[1641], message[1633], decoded[1641];
        uint16_t coeff[962];
        for (unsigned trial = 0; trial < 3000; trial++) {
            for (unsigned i = 0; i < n; i++)
                values[i] = (int16_t)((int)(word() % 301) - 150);
            CHECK(!xp_response_encode(encoded, values, n));
            CHECK(!xp_response_decode(a, encoded, n));
            CHECK(!memcmp(a, values, 2 * n));
            for (unsigned j = 0; j < 8; j++) {
                memcpy(changed, encoded, xp_response_bytes(n));
                unsigned pos = word() % (xp_response_bytes(n) * 8);
                changed[pos / 8] ^= (uint8_t)(1u << (pos % 8));
                int x = xp_response_decode(a, changed, n), y = reference(b, changed, n);
                CHECK(x == y);
                if (!x)
                    CHECK(!memcmp(a, b, 2 * n));
            }
            size_t len = word() % (xp_recovery_bytes(n) + 1), got = 0;
            for (size_t i = 0; i < len; i++)
                message[i] = (uint8_t)word();
            CHECK(!xp_mrm_encode(coeff, n, message, len, packed));
            CHECK(!xp_mrm_decode(decoded, &got, n, coeff));
            CHECK(got == len && !memcmp(message, decoded, len));
        }
        memset(coeff, 0, sizeof coeff);
        size_t len = 7;
        CHECK(xp_mrm_decode(decoded, &len, n, coeff) < 0 && !len);
        coeff[0] = 1;
        CHECK(!xp_mrm_decode(decoded, &len, n, coeff) && !len);
        coeff[0] = 2;
        CHECK(xp_mrm_decode(decoded, &len, n, coeff) < 0);
        coeff[0] = 12289;
        CHECK(xp_mrm_decode(decoded, &len, n, coeff) < 0);
        for (unsigned i = 0; i < xp_mu(n); i++)
            coeff[i] = 12288;
        CHECK(xp_mrm_decode(decoded, &len, n, coeff) < 0);
        printf("PASS %u: 24000 decoder mutations and 3000 variable byte-message round trips; "
               "invalid MRM padding/ranges\n",
               n);
    }
    return 0;
}
