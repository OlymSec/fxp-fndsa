#include "rec_ntt.h"
#include <string.h>
#if FALCON_IMPL_VECTOR_NTT
#include <pthread.h>
#define REC_PRIME UINT32_C(2147473409)
static falcon_rec_ntt_plan rec_plans[2];
static pthread_once_t rec_once[2]={PTHREAD_ONCE_INIT,PTHREAD_ONCE_INIT};
static int rec_status[2];
static void rec_init512(void){rec_status[0]=falcon_rec_ntt_init(rec_plans,512,REC_PRIME);}
static void rec_init1024(void){rec_status[1]=falcon_rec_ntt_init(rec_plans+1,1024,REC_PRIME);}
/* Full I32 input conversion, without signed negation or data-dependent division. */
static uint32_t rec_reduce(int32_t x) {
    uint32_t sign=(uint32_t)x>>31,mask=0U-sign;
    uint32_t mag=((uint32_t)x^mask)+sign,over=0U-(uint32_t)(mag>=REC_PRIME);
#if defined(__GNUC__) || defined(__clang__)
    __asm__("" : "+r"(over));
#endif
    uint32_t r=mag-(REC_PRIME&over),negative=(REC_PRIME-r)&(0U-(uint32_t)(r!=0));
    return (negative&mask)|(r&~mask);
}
int study_rec_prepare(study_rec_key *c,const int32_t key[4][1024],unsigned logn) {
    if(!c||!key||(logn!=9&&logn!=10))return 0;
    memset(c,0,sizeof *c);c->logn=logn;
    unsigned slot=logn-9,n=1U<<logn;
    if(pthread_once(rec_once+slot,slot?rec_init1024:rec_init512)||rec_status[slot])return 0;
    c->plan=rec_plans+slot;const falcon_rec_ntt_plan *p=c->plan;
    for(unsigned k=0;k<4;k++) {
        for(unsigned j=0;j<n;j++)c->key[k][j]=rec_reduce(key[k][j]);
        falcon_rec_ntt(p,c->key[k],0);
    }
    /* The same exact NTRU check on the already admitted small basis. */
    uint32_t v[1024];
    for(unsigned j=0;j<n;j++)v[j]=falcon_rec_sub(p,
        falcon_rec_mul(p,c->key[2][j],c->key[1][j]),
        falcon_rec_mul(p,c->key[0][j],c->key[3][j]));
    falcon_rec_ntt(p,v,1);unsigned valid=1;
    for(unsigned j=0;j<n;j++)valid&=v[j]==(j?0U:12289U);
    const uint32_t mont=(uint32_t)((UINT64_C(1)<<32)%REC_PRIME);
    for(unsigned k=0;k<4;k++)for(unsigned j=0;j<n;j++)c->key[k][j]=falcon_rec_mul(p,c->key[k][j],mont);
    return (int)valid;
}
static _Thread_local uint32_t study_ntt_workspace[3][1024];
void study_rec_clear_temporary(void (*clear)(void *,size_t)) {
    clear(study_ntt_workspace,sizeof study_ntt_workspace);
}
int study_rec_pair(const study_rec_key *c,int32_t *s0,int32_t *s1,
                   const int32_t *z0,const int32_t *z1,const uint16_t *hm) {
    const falcon_rec_ntt_plan *p=c->plan;unsigned n=1U<<c->logn;
    uint32_t *a=study_ntt_workspace[0],*b=study_ntt_workspace[1],*v=study_ntt_workspace[2];
    for(unsigned j=0;j<n;j++){a[j]=rec_reduce(z0[j]);b[j]=rec_reduce(z1[j]);}
    falcon_rec_ntt_pair(p,a,b);
    for(unsigned k=0;k<2;k++) {
        falcon_rec_response(p,v,a,c->key[2*k],b,c->key[2*k+1],NULL);
        falcon_rec_ntt(p,v,1);
        for(unsigned j=0;j<n;j++) {
            uint32_t mask=0U-(uint32_t)(v[j]>REC_PRIME/2);
#if defined(__GNUC__) || defined(__clang__)
            __asm__("" : "+r"(mask));
#endif
            int32_t x=(int32_t)((int64_t)v[j]-(int64_t)(mask&REC_PRIME));
            if(k)s1[j]=x;else s0[j]=(int32_t)hm[j]-x;
        }
    }
    return 1;
}
#else
#define NTRUGEN_PREFIX ntrugen
#define NTRUGEN_AVX2 0
#define NTRUGEN_ASM_CORTEXM4 0
#define NTRUGEN_STATS 0
#include "ng_inner.h"
int study_rec_prepare(study_rec_key *c,const int32_t key[4][1024],unsigned logn)
{
    memset(c,0,sizeof *c);c->logn=logn;
    unsigned n=1U<<logn;const small_prime *p=PRIMES;
    mp_mkgmigm(logn,c->gm,c->igm,p->g,p->ig,p->p,p->p0i);
    for(unsigned slot=0;slot<4;slot++) {
        for(unsigned j=0;j<n;j++) {
            int32_t x=key[slot][j];
            c->key[slot][j]=mp_set(x,p->p);
        }
        mp_NTT(logn,c->key[slot],c->gm,p->p,p->p0i);
        for(unsigned j=0;j<n;j++)
            c->key[slot][j]=mp_montymul(c->key[slot][j],p->R2,p->p,p->p0i);
    }
    /* Exact NTRU check: products of small key coefficients are bounded
     * by 2*n*127^2 < p/2, so modular equality is integer equality. */
    uint32_t v[1024];
    for(unsigned j=0;j<n;j++) {
        uint32_t fG=mp_montymul(c->key[2][j],c->key[1][j],p->p,p->p0i);
        uint32_t gF=mp_montymul(c->key[0][j],c->key[3][j],p->p,p->p0i);
        v[j]=mp_sub(fG,gF,p->p);
    }
    mp_iNTT(logn,v,c->igm,p->p,p->p0i);
    unsigned valid=1;
    for(unsigned j=0;j<n;j++)valid&=mp_montymul(v[j],1,p->p,p->p0i)==(j?0:12289);
    return (int)valid;
}
static _Thread_local uint32_t study_ntt_workspace[3][1024];
void study_rec_clear_temporary(void (*clear)(void *,size_t)) {
    clear(study_ntt_workspace,sizeof study_ntt_workspace);
}
int study_rec_pair(const study_rec_key *c,int32_t *s0,int32_t *s1,
                   const int32_t *z0,const int32_t *z1,const uint16_t *hm)
{
    unsigned n=1U<<c->logn;const small_prime *p=PRIMES;
    /* Fixed single-prime path; inherited bounds are listed in ../EXACTNESS.md. */
    uint32_t *a=study_ntt_workspace[0],*b=study_ntt_workspace[1],*v=study_ntt_workspace[2];
    for(unsigned j=0;j<n;j++) {a[j]=mp_set(z0[j],p->p);b[j]=mp_set(z1[j],p->p);}
    mp_NTT(c->logn,a,c->gm,p->p,p->p0i);
    mp_NTT(c->logn,b,c->gm,p->p,p->p0i);
    for(unsigned i=0;i<2;i++) {
        for(unsigned j=0;j<n;j++) {
            uint32_t x=mp_montymul(a[j],c->key[2*i][j],p->p,p->p0i);
            uint32_t y=mp_montymul(b[j],c->key[2*i+1][j],p->p,p->p0i);
            v[j]=mp_add(x,y,p->p);
        }
        mp_iNTT(c->logn,v,c->igm,p->p,p->p0i);
        for(unsigned j=0;j<n;j++) {
            int32_t x=mp_norm(v[j],p->p);
            if(i) s1[j]=x;else s0[j]=(int32_t)hm[j]-x;
        }
    }
    return 1;
}

#endif
