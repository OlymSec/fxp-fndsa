#ifndef FALCON_EXACT_NTT_H
#define FALCON_EXACT_NTT_H
#define FALCON_REC_UNSUPPORTED 1
#define FALCON_REC_INTERNAL 2
#include <stdint.h>
__extension__ typedef __int128 falcon_rec_i128;
__extension__ typedef unsigned __int128 falcon_rec_u128;
#define FALCON_REC_MAX_N 1024
typedef struct {
    unsigned n;
    uint32_t p,inv_n,p0i,rho,irho,rho_inverse_difference;
    uint64_t reciprocal;
    uint16_t reverse[FALCON_REC_MAX_N];
    uint32_t twist[FALCON_REC_MAX_N],untwist[FALCON_REC_MAX_N];
    uint32_t roots[FALCON_REC_MAX_N],iroots[FALCON_REC_MAX_N];
    uint32_t stages[FALCON_REC_MAX_N],istages[FALCON_REC_MAX_N];
    uint32_t shoup[FALCON_REC_MAX_N],ishoup[FALCON_REC_MAX_N];
    /* Bit-reversed-spectrum schedule (n=512/1024, AVX2): psi^brv(k), inverses,
     * Shoup factors, and per-lane twiddles for the three in-register layers. */
    uint32_t ct_w[FALCON_REC_MAX_N],ct_f[FALCON_REC_MAX_N],gs_w[FALCON_REC_MAX_N],gs_f[FALCON_REC_MAX_N];
    uint32_t ct_lw[3][FALCON_REC_MAX_N/2],ct_lf[3][FALCON_REC_MAX_N/2],gs_lw[3][FALCON_REC_MAX_N/2],gs_lf[3][FALCON_REC_MAX_N/2];
    /* Mixed degree uses ct_w/gs_w for Montgomery roots and ct_f/gs_f for
     * its public cubic roots, rho and interpolation constants. */
    uint32_t gs_ninv,gs_ninv_f,gs_last,gs_last_f;
} falcon_rec_ntt_plan;
int falcon_rec_ntt_init(falcon_rec_ntt_plan *c,unsigned n,uint32_t p);
void falcon_rec_ntt(const falcon_rec_ntt_plan *p,uint32_t *a,int inverse);
void falcon_rec_ntt_pair(const falcon_rec_ntt_plan *,uint32_t *,uint32_t *);
void falcon_rec_linear(const falcon_rec_ntt_plan *,uint32_t *,const uint32_t *,const uint32_t *,
               const uint32_t *,const uint32_t *,const uint32_t *,const uint32_t *,
               const uint32_t *,const uint32_t *);
void falcon_rec_response(const falcon_rec_ntt_plan *,uint32_t *,const uint32_t *,const uint32_t *,
                  const uint32_t *,const uint32_t *,const uint32_t *);
uint32_t falcon_rec_mod(falcon_rec_i128 a,uint32_t p);
uint32_t falcon_rec_mul(const falcon_rec_ntt_plan *p,uint32_t a,uint32_t b);
uint32_t falcon_rec_add(const falcon_rec_ntt_plan *p,uint32_t a,uint32_t b);
uint32_t falcon_rec_sub(const falcon_rec_ntt_plan *p,uint32_t a,uint32_t b);
#endif
