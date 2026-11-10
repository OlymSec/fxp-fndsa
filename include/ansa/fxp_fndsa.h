/* ANSA SDK: 64-bit fixed-point Falcon with the specified FN-DSA protocol.
 * 512/1024 use fndsa's NTT keys, LSB-first encodings and message representative.
 * Only degrees 512 and 1024 are supported.
 * The fixed-point sampler differs from fndsa; seeded signatures deliberately
 * differ through a distinct implementation identifier. This is not a claim
 * of final FIPS 206 conformance. Buffers must not overlap. */
#ifndef ANSA_FXP_FNDSA_H
#define ANSA_FXP_FNDSA_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
#if defined(__GNUC__) || defined(__clang__)
#define FXP_FNDSA_API __attribute__((visibility("default")))
#else
#define FXP_FNDSA_API
#endif
/* Key-generation and explicit signing seeds are exactly 64 bytes. NULL uses
 * the OS RNG (64 bytes for keygen, 40 bytes for hedged signing). */
#define FXP_FNDSA_SEED_BYTES         64
#define FXP_FNDSA_HASH_ID_RAW        ""
#define FXP_FNDSA_HASH_ID_SHA256     "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x01"
#define FXP_FNDSA_HASH_ID_SHA384     "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x02"
#define FXP_FNDSA_HASH_ID_SHA512     "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x03"
#define FXP_FNDSA_HASH_ID_SHA512_256 "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x06"
#define FXP_FNDSA_HASH_ID_SHA3_256   "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x08"
#define FXP_FNDSA_HASH_ID_SHA3_384   "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x09"
#define FXP_FNDSA_HASH_ID_SHA3_512   "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x0A"
#define FXP_FNDSA_HASH_ID_SHAKE128   "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x0B"
#define FXP_FNDSA_HASH_ID_SHAKE256   "\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x0C"
#define FXP_FNDSA_HASH_ID_EXTMU      "\xFE"
enum {
    FXP_FNDSA_OK = 0,
    FXP_FNDSA_RANDOM = -1,
    FXP_FNDSA_SIZE = -2,
    FXP_FNDSA_FORMAT = -3,
    FXP_FNDSA_BAD_SIGNATURE = -4,
    FXP_FNDSA_BAD_ARGUMENT = -5,
    FXP_FNDSA_INTERNAL = -6,
    FXP_FNDSA_MEMORY = -7,
    FXP_FNDSA_EXHAUSTED = -8
};
typedef struct fxp_fndsa_prepared fxp_fndsa_prepared;
/* hash_id is RAW (NULL also means RAW), a DER OID as in fndsa.h, or EXTMU.
 * For a prehash, the caller supplies the digest. Context length is <=255.
 * EXTMU takes exactly 64 bytes and ignores context, as fndsa does. */
typedef struct {
    const uint8_t *data;
    size_t data_len;
    const uint8_t *context;
    size_t context_len;
    const char *hash_id;
} fxp_fndsa_message;
FXP_FNDSA_API size_t fxp_fndsa_public_key_bytes(unsigned degree);
FXP_FNDSA_API size_t fxp_fndsa_secret_key_bytes(unsigned degree);
FXP_FNDSA_API size_t fxp_fndsa_signature_bytes(unsigned degree);
FXP_FNDSA_API size_t fxp_fndsa_mrm_signature_bytes(unsigned degree);
FXP_FNDSA_API size_t fxp_fndsa_mrm_max_bytes(unsigned degree);
FXP_FNDSA_API int fxp_fndsa_keygen(unsigned degree, const uint8_t *seed, uint8_t *public_key,
                                   size_t public_capacity, uint8_t *secret_key,
                                   size_t secret_capacity);
FXP_FNDSA_API int fxp_fndsa_public_from_secret(const uint8_t *secret_key, size_t secret_len,
                                               uint8_t *public_key, size_t public_capacity);
/* Preparation owns immutable sampler caches. Temporary secret storage is
 * on the heap and wiped. Prepared keys may be shared across signing threads. */
FXP_FNDSA_API int fxp_fndsa_prepare(fxp_fndsa_prepared **out, const uint8_t *secret_key,
                                    size_t secret_len);
FXP_FNDSA_API void fxp_fndsa_free(fxp_fndsa_prepared *key);
/* Raw-message convenience APIs use an empty context. All APIs return 0 on
 * success and a negative status on error. written is zero on failure. */
FXP_FNDSA_API int fxp_fndsa_sign(const fxp_fndsa_prepared *, const uint8_t *message,
                                 size_t message_len, const uint8_t *seed, uint8_t *signature,
                                 size_t capacity, size_t *written);
FXP_FNDSA_API int fxp_fndsa_verify(const uint8_t *public_key, size_t public_len,
                                   const uint8_t *message, size_t message_len,
                                   const uint8_t *signature, size_t signature_len);
/* MRM uses variable-length BYTES and the trailing-one/base-q encoding in
 * fndsa's ds-spec.txt. Full H1/H2 tags are used. Thus MRM does not interoperate
 * with the uncorrected two-byte-tag bug in the pinned fndsa implementation.
 * Recovery output is written only after authentication; empty recovery is
 * success with recovered_len=0, distinct from every negative error status. */
FXP_FNDSA_API int fxp_fndsa_sign_mrm(const fxp_fndsa_prepared *, const uint8_t *embedded,
                                     size_t embedded_len, const uint8_t *message,
                                     size_t message_len, const uint8_t *seed, uint8_t *signature,
                                     size_t capacity, size_t *written);
FXP_FNDSA_API int fxp_fndsa_verify_mrm(const uint8_t *public_key, size_t public_len,
                                       const uint8_t *message, size_t message_len,
                                       const uint8_t *signature, size_t signature_len,
                                       uint8_t *embedded, size_t embedded_capacity,
                                       size_t *recovered_len);
FXP_FNDSA_API int fxp_fndsa_sign_message(const fxp_fndsa_prepared *, const fxp_fndsa_message *,
                                         const uint8_t *seed, uint8_t *signature, size_t capacity,
                                         size_t *written);
FXP_FNDSA_API int fxp_fndsa_verify_message(const uint8_t *public_key, size_t public_len,
                                           const fxp_fndsa_message *, const uint8_t *signature,
                                           size_t signature_len);
FXP_FNDSA_API int fxp_fndsa_sign_mrm_message(const fxp_fndsa_prepared *, const uint8_t *embedded,
                                             size_t embedded_len, const fxp_fndsa_message *,
                                             const uint8_t *seed, uint8_t *signature,
                                             size_t capacity, size_t *written);
FXP_FNDSA_API int fxp_fndsa_verify_mrm_message(const uint8_t *public_key, size_t public_len,
                                               const fxp_fndsa_message *, const uint8_t *signature,
                                               size_t signature_len, uint8_t *embedded,
                                               size_t embedded_capacity, size_t *recovered_len);
#ifdef __cplusplus
}
#endif
#endif
