#ifndef STUDY_REC_NTT_H
#define STUDY_REC_NTT_H
#include <stdint.h>
#include <stddef.h>
#include "opt_config.h"
#if FALCON_IMPL_VECTOR_NTT
#include "exact_ntt.h"
#endif
typedef struct {
    unsigned logn;
#if FALCON_IMPL_VECTOR_NTT
    const falcon_rec_ntt_plan *plan;
    uint32_t key[4][1024];
#else
    uint32_t gm[1024],igm[1024],key[4][1024];
#endif
} study_rec_key;
int study_rec_prepare(study_rec_key *,const int32_t [4][1024],unsigned);
/* Admitted signing coordinates only; fixed path, no runtime range guard. */
int study_rec_pair(const study_rec_key *,int32_t *,int32_t *,
                   const int32_t *,const int32_t *,const uint16_t *);
void study_rec_clear_temporary(void (*clear)(void *,size_t));
#endif
