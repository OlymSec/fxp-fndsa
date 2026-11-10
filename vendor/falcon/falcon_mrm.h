#ifndef FALCON_MRM_H
#define FALCON_MRM_H

#include "falcon.h"

#ifdef __cplusplus
extern "C" {
#endif

size_t falcon_mrm_message_coeffs(unsigned logn);
size_t falcon_mrm_signature_size(unsigned logn);

/* These interfaces require a build with -DFALCON_MRM_EXTRACT_Q64=1. */
int falcon_mrm_sign_seeded(
	void *signature, size_t signature_capacity, size_t *signature_length,
	const void *expanded_key, size_t expanded_key_length,
	const void *public_key, size_t public_key_length,
	const uint16_t *recoverable, size_t recoverable_count,
	const uint8_t *explicit_message, size_t explicit_length,
	const uint8_t seed[64], void *tmp, size_t tmp_len);
int falcon_mrm_sign(
	void *signature, size_t signature_capacity, size_t *signature_length,
	const void *expanded_key, size_t expanded_key_length,
	const void *public_key, size_t public_key_length,
	const uint16_t *recoverable, size_t recoverable_count,
	const uint8_t *explicit_message, size_t explicit_length,
	void *tmp, size_t tmp_len);
int falcon_mrm_verify(
	const void *signature, size_t signature_length,
	const void *public_key, size_t public_key_length, const uint8_t *explicit_message,
	size_t explicit_length, uint16_t *recoverable,
	size_t recoverable_capacity, void *tmp, size_t tmp_len);

#ifdef __cplusplus
}
#endif

#endif
