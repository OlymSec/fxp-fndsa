#ifndef XP_POW2_PREPARE_H
#define XP_POW2_PREPARE_H
#include "protocol.h"
/* basis is four consecutive n-byte vectors f,g,F,G. Decode/recovery and
 * public derivation are outside the preallocated core operation. */
int xp_pow2_decode(int8_t *basis, uint16_t *public_ntt, unsigned logn, const uint8_t *raw,
                   size_t rawlen, void *work);
int xp_pow2_prepare_core(void *expanded, void *cache, const int8_t *basis, unsigned logn,
                         void *work);
int study_prepare_basis_work(void *cache, const void *expanded, const int8_t *basis, void *work);
int study_prepare_expanded_work(void *cache, const void *expanded, void *work);
#endif
