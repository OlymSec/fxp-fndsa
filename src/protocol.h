#ifndef FXP_PROTOCOL_H
#define FXP_PROTOCOL_H
#include "ansa/fxp_fndsa.h"
#include "inner.h"
#define XP_Q       12289U
#define XP_MAX_N   1024U
#define XP_MAX_REC 1633U
/* Private protocol helpers; not exported by the library. */
void xp_clean(void *, size_t);
void xp_hash(uint8_t *, size_t, const void *, size_t);
unsigned xp_degree(const uint8_t *, size_t, unsigned kind);
size_t xp_header(unsigned n);
void xp_set_header(uint8_t *, unsigned n, unsigned kind);
size_t xp_response_bytes(unsigned n);
unsigned xp_bound(unsigned n);
unsigned xp_lambda(unsigned n);
unsigned xp_gamma(unsigned n);
unsigned xp_mu(unsigned n);
size_t xp_recovery_bytes(unsigned n);
void xp_pack(uint8_t *, const uint16_t *, unsigned, unsigned);
int xp_unpack(uint16_t *, const uint8_t *, unsigned, unsigned, unsigned);
int xp_response_encode(uint8_t *, const int16_t *, unsigned);
int xp_response_decode(int16_t *, const uint8_t *, unsigned);
int xp_message_mu(uint8_t[64], const uint8_t[64], const fxp_fndsa_message *);
void xp_coefficients(inner_shake256_context *, uint16_t *, unsigned);
void xp_inject_u16(inner_shake256_context *, const uint16_t *, unsigned);
void xp_target(uint16_t *, unsigned, const uint8_t[40], const uint8_t[64]);
void xp_mrm_h1(uint16_t *, unsigned, const uint16_t *, const uint8_t[64], const uint8_t *, size_t);
void xp_mrm_h2(uint16_t *, unsigned, const uint16_t *);
int xp_mrm_encode(uint16_t *, unsigned, const uint8_t *, size_t, uint8_t *);
int xp_mrm_decode(uint8_t *, size_t *, unsigned, const uint16_t *);
int xp_norm(unsigned, const int16_t *, const int16_t *);
int xp_tree_attempt(int16_t *, int16_t *, inner_shake256_context *, const fpr *, const uint16_t *,
                    unsigned, uint8_t *);
#endif
