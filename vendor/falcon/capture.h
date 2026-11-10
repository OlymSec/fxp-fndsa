#ifndef STUDY_CAPTURE_H
#define STUDY_CAPTURE_H
#include "rec_ntt.h"
#include <assert.h>
#include "certify.h"
typedef struct {
    study_rec_key key;
    const fpr *expanded;
    unsigned order[2048];
#if FALCON_IMPL_MIRROR_CACHE
    fpr mirrors[512]; /* DFS order of the stored left-subtree leaves. */
#endif
} study_prepared;
static _Thread_local const study_prepared *study_active;
#if FALCON_IMPL_MIRROR_CACHE
static const fpr *study_mirrored_widths(void) {
    return study_active ? study_active->mirrors : NULL;
}
#endif
size_t study_prepared_size(void) { return sizeof(study_prepared); }
void study_select_prepared(const void *cache) { study_active=cache; }

int study_prepare_expanded(void *memory,const void *encoded_expanded)
{
    study_prepared *prepared=memory;
    if(!prepared) return -1;
    memset(prepared,0,sizeof *prepared);
    if(!encoded_expanded)return -1;
    unsigned logn=*(const uint8_t *)encoded_expanded;
    if(logn!=9 && logn!=10)return -1;
    uintptr_t ptr=((uintptr_t)encoded_expanded+sizeof(fpr))&~(uintptr_t)(sizeof(fpr)-1);
    const fpr *candidate=(const fpr *)ptr;
    unsigned n=1U<<logn;
    int32_t key[4][1024];fpr a[1024];
    unsigned slotmap[4]={0,2,1,3};
    for(unsigned slot=0;slot<4;slot++) {
        memcpy(a,candidate+slotmap[slot]*n,n*sizeof *a);
        Zf(iFFT)(a,logn);
        for(unsigned j=0;j<n;j++) {
            int64_t x=fpr_rint(a[j]);if(slot>=2)x=-x;
            if(x< -127 || x>127)return -2;key[slot][j]=(int32_t)x;
        }
        for(unsigned j=0;j<n;j++)a[j]=fpr_of(slot>=2?-(int64_t)key[slot][j]:key[slot][j]);
        Zf(FFT)(a,logn);
        if(memcmp(a,candidate+slotmap[slot]*n,n*sizeof *a))return -3;
    }
    fpr *mirrors=NULL;
#if FALCON_IMPL_MIRROR_CACHE
    mirrors=prepared->mirrors;
#endif
    if(!cert_basis_tree(candidate,logn,mirrors))return -4;
    if(!study_rec_prepare(&prepared->key,key,logn))return -5;
    for(unsigned k=0;k<2*n;k++) {
        unsigned local=k%n,rev=0,v=n/2-1-local/2;
        for(unsigned bit=0;bit<logn-1;bit++){rev=(rev<<1)|(v&1);v>>=1;}
        prepared->order[k]=(k<n?n:0)+rev+(local&1)*(n/2);
    }
    prepared->expanded=candidate;return 0;
}
typedef struct {
    samplerZ sampler;void *ctx;unsigned count;int32_t z[2048];
} study_capture;
typedef struct {
    study_capture capture;
    int64_t shift0[1024],shift1[1024];
    int32_t full0[1024],full1[1024],exact0[1024],exact1[1024];
} study_workspace;
static _Thread_local study_workspace study_work;
void study_clear_temporary(void (*clear)(void *,size_t)) {
    clear(&study_work,sizeof study_work);
    study_rec_clear_temporary(clear);
}
static int study_capture_sample(void *ctx,fpr mu,fpr sigma)
{
    study_capture *c=ctx;
    int z=c->sampler(c->ctx,mu,sigma);
    c->z[study_active->order[c->count++]]=z;
    return z;
}
#endif
