#include "ng_inner.h"
#include "division.h"
typedef __uint128_t ct3_u128;

/* 1 if a < b, else 0. */
static inline uint64_t ct3l_lt(uint64_t a, uint64_t b)
{
    return (uint64_t)(((ct3_u128)a - b) >> 64) & 1U;
}

/* The barrier keeps the compiler from rebuilding selects on secrets as
 * branches or address selects. */
static inline uint64_t ct3l_mask(uint64_t bit)
{
    uint64_t m = (uint64_t)0 - (bit & 1U);
    __asm__("" : "+r"(m));
    return m;
}

/* Reciprocal of a normalized limb d (d >= 2^63):
 *   v = floor((2^128 - 1) / d) - 2^64 = floor(((~d) * 2^64 + (2^64 - 1)) / d).
 * The high word ~d is < d, so the quotient fits in 64 bits.  Restoring
 * division with exactly 64 steps; each step shifts one 1-bit of the low word
 * into the partial remainder r < d, so (top:r) < 2d and one masked
 * subtraction restores r < d.  Straight-line code (no loop) so that the
 * object audit can check it for branches. */
#define CT3L_RECIP_STEP                                                                            \
    do {                                                                                           \
        uint64_t top_ = r >> 63;                                                                   \
        r = (r << 1) | 1U;                                                                         \
        uint64_t ge_ = top_ | (ct3l_lt(r, d) ^ 1U);                                                \
        r -= d & ct3l_mask(ge_);                                                                   \
        q = (q << 1) | ge_;                                                                        \
    } while (0)
#define CT3L_RECIP_STEP4                                                                           \
    CT3L_RECIP_STEP;                                                                               \
    CT3L_RECIP_STEP;                                                                               \
    CT3L_RECIP_STEP;                                                                               \
    CT3L_RECIP_STEP
#define CT3L_RECIP_STEP16                                                                          \
    CT3L_RECIP_STEP4;                                                                              \
    CT3L_RECIP_STEP4;                                                                              \
    CT3L_RECIP_STEP4;                                                                              \
    CT3L_RECIP_STEP4
static uint64_t ct3l_reciprocal(uint64_t d)
{
    uint64_t r = ~d, q = 0;
    CT3L_RECIP_STEP16;
    CT3L_RECIP_STEP16;
    CT3L_RECIP_STEP16;
    CT3L_RECIP_STEP16;
    return q;
}
#undef CT3L_RECIP_STEP16
#undef CT3L_RECIP_STEP4
#undef CT3L_RECIP_STEP

/* floor((u1 * 2^64 + u0) / d) for normalized d, u1 < d, v the reciprocal of
 * d: Moeller-Granlund "Improved division by invariant integers" (IEEE TC
 * 2011), Algorithm 4, with both adjustment branches replaced by masks.  The
 * algorithm proves that at most one decrement and then at most one increment
 * yield the exact quotient. */
static uint64_t ct3l_div_2by1(uint64_t u1, uint64_t u0, uint64_t d, uint64_t v)
{
    ct3_u128 p = (ct3_u128)v * u1 + (((ct3_u128)u1 << 64) | u0);
    uint64_t q1 = (uint64_t)(p >> 64) + 1U;
    uint64_t q0 = (uint64_t)p;
    uint64_t r = u0 - q1 * d;
    uint64_t m = ct3l_mask(ct3l_lt(q0, r));
    q1 += m;
    r += d & m;
    uint64_t m2 = ct3l_mask(ct3l_lt(r, d) ^ 1U);
    q1 -= m2;
    return q1;
}

/* 1 if a < b for 128-bit values (borrow of a two-word subtraction). */
static inline uint64_t ct3l_lt128(ct3_u128 a, ct3_u128 b)
{
    uint64_t borrow = ct3l_lt((uint64_t)a, (uint64_t)b);
    ct3_u128 t = (ct3_u128)(uint64_t)(a >> 64) - (uint64_t)(b >> 64) - borrow;
    return (uint64_t)(t >> 64) & 1U;
}

/* Trial quotient digit of Knuth's Algorithm D (TAOCP 4.3.1, step D3) for a
 * normalized divisor with top limbs d1 >= 2^63 and d0 (d0 = 0 when the
 * divisor has one limb), and window top limbs u2 <= d1, u1, u0.
 *
 * First q2 = min(floor((u2*2^64 + u1) / d1), 2^64 - 1): when u2 == d1 the
 * value 2^64 - 1 is selected and u2 is replaced by 0 in the 2-by-1 division
 * so that its precondition holds.  Then one masked step of Knuth's test
 *   rhat < 2^64 and q2*d0 > rhat*2^64 + u0,  rhat = u2*2^64 + u1 - q2*d1,
 * i.e. q2*(d1*2^64 + d0) > u2*2^128 + u1*2^64 + u0, decrements q2 by one.
 * With q the true digit and q3 = floor(U3/V2) the 3-by-2 quotient of the
 * top limbs: q <= q3 <= q + 1 and q <= q2 <= q + 2 (Theorem 4.3.1B), so the
 * result is in [q, q + 1] (proof/NTRU-COMPLETION-HARDENING-3F.md). */
static uint64_t ct3l_trial_quotient(uint64_t u2, uint64_t u1, uint64_t u0, uint64_t d1, uint64_t d0,
                                    uint64_t v)
{
    uint64_t x = u2 ^ d1;
    uint64_t equal = ct3l_mask(((x | ((uint64_t)0 - x)) >> 63) ^ 1U);
    uint64_t qhat = ct3l_div_2by1(u2 & ~equal, u1, d1, v) | equal;
    ct3_u128 rhat = (((ct3_u128)u2 << 64) | u1) - (ct3_u128)qhat * d1;
    uint64_t small = ct3l_lt((uint64_t)(rhat >> 64), 1U);
    ct3_u128 rhs = ((ct3_u128)(uint64_t)rhat << 64) | u0;
    return qhat - (small & ct3l_lt128(rhs, (ct3_u128)qhat * d0));
}

/* row[0..an+1] = row[0..an-1] + a * (b0 + b1 * 2^64), where row[an] and
 * row[an+1] are written, not read.  Pending carries are a two-limb value
 * (c0 at the current position, c1 one above); every partial sum is at most
 * (2^64-1)^2 + 2(2^64-1) < 2^128. */
static void ct3l_addmul_2(uint64_t *row, const uint64_t *ap, unsigned an, uint64_t b0, uint64_t b1)
{
    uint64_t c0 = 0, c1 = 0;
    for (unsigned i = 0; i < an; i++) {
        uint64_t ai = ap[i];
        ct3_u128 t = (ct3_u128)ai * b0 + row[i] + c0;
        row[i] = (uint64_t)t;
        ct3_u128 u = (ct3_u128)ai * b1 + (uint64_t)(t >> 64) + c1;
        c0 = (uint64_t)u;
        c1 = (uint64_t)(u >> 64);
    }
    row[an] = c0;
    row[an + 1U] = c1;
}

/* row[0..an] = row[0..an-1] + a * b0, where row[an] is written, not read. */
static void ct3l_addmul_1(uint64_t *row, const uint64_t *ap, unsigned an, uint64_t b0)
{
    uint64_t carry = 0;
    for (unsigned i = 0; i < an; i++) {
        ct3_u128 t = (ct3_u128)ap[i] * b0 + row[i] + carry;
        row[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
    row[an] = carry;
}

static void antrag_fx_ct3_limbs_mul(uint64_t *rp, const uint64_t *ap, unsigned an, const uint64_t *bp,
                             unsigned bn)
{
    for (unsigned i = 0; i < an; i++)
        rp[i] = 0;
    unsigned j = 0;
    for (; j + 1U < bn; j += 2U)
        ct3l_addmul_2(rp + j, ap, an, bp[j], bp[j + 1U]);
    if (j < bn)
        ct3l_addmul_1(rp + j, ap, an, bp[j]);
}

static void antrag_fx_ct3_limbs_mul_low(uint64_t *rp, const uint64_t *ap, const uint64_t *bp, unsigned n)
{
    for (unsigned i = 0; i < n; i++)
        rp[i] = 0;
    for (unsigned j = 0; j < n; j++) {
        uint64_t bj = bp[j], carry = 0;
        uint64_t *row = rp + j;
        for (unsigned i = 0; i < n - j; i++) {
            ct3_u128 t = (ct3_u128)ap[i] * bj + row[i] + carry;
            row[i] = (uint64_t)t;
            carry = (uint64_t)(t >> 64);
        }
    }
}

static void antrag_fx_ct3_limbs_sqr(uint64_t *rp, const uint64_t *ap, unsigned n)
{
    /* Off-diagonal sum S = sum_{i<j} a_i a_j 2^(64(i+j)), then 2S + diagonal.
     * 2S <= a^2 < 2^(128n), so the doubling shift loses no bit. */
    for (unsigned i = 0; i < 2U * n; i++)
        rp[i] = 0;
    for (unsigned i = 0; i + 1U < n; i++) {
        uint64_t ai = ap[i], carry = 0;
        for (unsigned j = i + 1U; j < n; j++) {
            ct3_u128 t = (ct3_u128)ai * ap[j] + rp[i + j] + carry;
            rp[i + j] = (uint64_t)t;
            carry = (uint64_t)(t >> 64);
        }
        rp[i + n] = carry;
    }
    uint64_t high = 0;
    for (unsigned i = 0; i < 2U * n; i++) {
        uint64_t word = rp[i];
        rp[i] = (word << 1) | high;
        high = word >> 63;
    }
    uint64_t carry = 0;
    for (unsigned i = 0; i < n; i++) {
        ct3_u128 square = (ct3_u128)ap[i] * ap[i];
        ct3_u128 low = (ct3_u128)rp[2U * i] + (uint64_t)square + carry;
        rp[2U * i] = (uint64_t)low;
        ct3_u128 top = (ct3_u128)rp[2U * i + 1U] + (uint64_t)(square >> 64) + (uint64_t)(low >> 64);
        rp[2U * i + 1U] = (uint64_t)top;
        carry = (uint64_t)(top >> 64);
    }
}

/* w[0..dn] -= q * d[0..dn-1]; returns the borrow out of w[dn]. */
static uint64_t ct3l_submul(uint64_t *w, const uint64_t *dp, unsigned dn, uint64_t q)
{
    uint64_t carry = 0, borrow = 0;
    for (unsigned i = 0; i < dn; i++) {
        ct3_u128 p = (ct3_u128)q * dp[i] + carry;
        carry = (uint64_t)(p >> 64);
        ct3_u128 t = (ct3_u128)w[i] - (uint64_t)p - borrow;
        w[i] = (uint64_t)t;
        borrow = (uint64_t)(t >> 64) & 1U;
    }
    ct3_u128 t = (ct3_u128)w[dn] - carry - borrow;
    w[dn] = (uint64_t)t;
    return (uint64_t)(t >> 64) & 1U;
}

/* w[0..dn] += d & mask; returns the carry out of w[dn]. */
static uint64_t ct3l_add_masked(uint64_t *w, const uint64_t *dp, unsigned dn, uint64_t mask)
{
    uint64_t carry = 0;
    for (unsigned i = 0; i < dn; i++) {
        ct3_u128 t = (ct3_u128)w[i] + (dp[i] & mask) + carry;
        w[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
    ct3_u128 t = (ct3_u128)w[dn] + carry;
    w[dn] = (uint64_t)t;
    return (uint64_t)(t >> 64);
}

uint64_t Zn(half_div_qr)(uint64_t *qp, uint64_t *np, unsigned nn, const uint64_t *dp,
                                    unsigned dn)
{
    const uint64_t d = dp[dn - 1U];
    const uint64_t v = ct3l_reciprocal(d);
    uint64_t *top = np + (nn - dn);

    /* Top quotient limb.  The top dn limbs T < 2^(64dn) <= 2D because D is
     * normalized, so floor(T/D) is 0 or 1 and one masked subtraction leaves
     * T mod D. */
    uint64_t borrow = 0;
    for (unsigned i = 0; i < dn; i++) {
        ct3_u128 t = (ct3_u128)top[i] - dp[i] - borrow;
        borrow = (uint64_t)(t >> 64) & 1U;
    }
    uint64_t qh = borrow ^ 1U;
    uint64_t mask = ct3l_mask(qh);
    borrow = 0;
    for (unsigned i = 0; i < dn; i++) {
        ct3_u128 t = (ct3_u128)top[i] - (dp[i] & mask) - borrow;
        top[i] = (uint64_t)t;
        borrow = (uint64_t)(t >> 64) & 1U;
    }

    /* Schoolbook digits.  Invariant: the dn limbs above position j hold a
     * partial remainder R < D, so the window W = np[j..j+dn] = R*2^64 + np[j]
     * satisfies W < 2^64 D and its top limb is <= d.  The trial quotient
     * satisfies q <= qhat <= q + 1 for the true digit q = floor(W / D), so
     * W - qhat*D lies in [-D, D).  It is formed modulo 2^(64(dn+1)) with the
     * borrow as its sign, and one masked add-back of D (with qhat decremented
     * by the sign) restores 0 <= W - q*D < D. */
    const uint64_t d0 = dn >= 2U ? dp[dn - 2U] : 0;
    for (unsigned j = nn - dn; j-- > 0;) {
        uint64_t *w = np + j;
        uint64_t u0 = dn >= 2U ? w[dn - 2U] : 0;
        uint64_t qhat = ct3l_trial_quotient(w[dn], w[dn - 1U], u0, d, d0, v);
        uint64_t negative = ct3l_submul(w, dp, dn, qhat);
        qhat -= negative;
        (void)ct3l_add_masked(w, dp, dn, ct3l_mask(negative));
        qp[j] = qhat;
    }
    return qh;
}

