#include "pow2_prepare.h"
#include "falcon.h"
#include <string.h>
#ifndef XP_POW2_HANDOFF
#define XP_POW2_HANDOFF 1
#endif

int xp_pow2_decode(int8_t *basis, uint16_t *public_ntt, unsigned logn, const uint8_t *raw,
                   size_t rawlen, void *work)
{
    if ((logn != 9 && logn != 10) || rawlen != FALCON_PRIVKEY_SIZE(logn) ||
        raw[0] != (uint8_t)(0x50 + logn))
        return FXP_FNDSA_FORMAT;
    size_t n = (size_t)1 << logn, pos = 1;
    for (unsigned role = 0; role < 3; role++) {
        unsigned bits = role == 2 ? Zf(max_FG_bits)[logn] : Zf(max_fg_bits)[logn];
        size_t len = Zf(trim_i8_decode)(basis + role * n, logn, bits, raw + pos, rawlen - pos);
        if (!len)
            return FXP_FNDSA_FORMAT;
        pos += len;
    }
    if (pos != rawlen ||
        !Zf(complete_private)(basis + 3 * n, basis, basis + n, basis + 2 * n, logn, work))
        return FXP_FNDSA_FORMAT;
    for (size_t i = 0; i < 4 * n; i++)
        if (basis[i] == -128)
            return FXP_FNDSA_FORMAT;
    if (!Zf(compute_public)(public_ntt, basis, basis + n, logn, work))
        return FXP_FNDSA_FORMAT;
    Zf(to_ntt_monty)(public_ntt, logn);
    for (size_t i = 0; i < n; i++)
        public_ntt[i] = (uint16_t)((uint32_t)public_ntt[i] * 2304U % 12289U);
    return FXP_FNDSA_OK;
}

int xp_pow2_prepare_core(void *expanded, void *cache, const int8_t *basis, unsigned logn,
                         void *work)
{
    if (logn != 9 && logn != 10)
        return FXP_FNDSA_FORMAT;
    size_t n = (size_t)1 << logn;
    fpr *tree = (fpr *)(((uintptr_t)expanded + sizeof(fpr)) & ~(uintptr_t)(sizeof(fpr) - 1));
    *(uint8_t *)expanded = (uint8_t)logn;
    Zf(expand_privkey)(tree, basis, basis + n, basis + 2 * n, basis + 3 * n, logn, work);
#if XP_POW2_HANDOFF
    return study_prepare_basis_work(cache, expanded, basis, work) ? FXP_FNDSA_FORMAT : FXP_FNDSA_OK;
#else
    return study_prepare_expanded_work(cache, expanded, work) ? FXP_FNDSA_FORMAT : FXP_FNDSA_OK;
#endif
}
