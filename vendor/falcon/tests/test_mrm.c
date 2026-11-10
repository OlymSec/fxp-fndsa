#include "falcon_mrm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(x) do { if (!(x)) { fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#x); return 1; } } while (0)

static int run_profile(unsigned logn)
{
    size_t pl=FALCON_PRIVKEY_SIZE(logn), ql=FALCON_PUBKEY_SIZE(logn);
    size_t el=FALCON_EXPANDEDKEY_SIZE(logn), size=falcon_mrm_signature_size(logn);
    size_t mu=falcon_mrm_message_coeffs(logn), tl=FALCON_TMPSIZE_KEYGEN(logn);
    if(tl<FALCON_TMPSIZE_EXPANDPRIV(logn)) tl=FALCON_TMPSIZE_EXPANDPRIV(logn);
    if(tl<FALCON_TMPSIZE_SIGNTREE(logn)) tl=FALCON_TMPSIZE_SIGNTREE(logn);
    uint8_t *priv=malloc(pl), *pk=malloc(ql), *expanded=malloc(el);
    uint8_t *scratch=malloc(tl+1), *sig=malloc(size+1), *repeat=malloc(size);
    uint16_t m0[962], decoded[962], sentinel[962];
    uint8_t seed[64]={0}, msg[32]={0};
    shake256_context rng;
    CHECK(priv && pk && expanded && scratch && sig && repeat);
    seed[0]=(uint8_t)logn;
    shake256_init_prng_from_seed(&rng,seed,sizeof seed);
    CHECK(!falcon_keygen_make(&rng,logn,priv,pl,pk,ql,scratch,tl));
    CHECK(!falcon_expand_privkey(expanded,el,priv,pl,scratch,tl));
    for(size_t u=0;u<mu;u++) m0[u]=(uint16_t)((u*17+logn)%12289);
    memset(sentinel,0xa5,sizeof sentinel);
    for(unsigned trial=0;trial<20;trial++) {
        size_t written=0, written2=0;
        seed[1]=(uint8_t)trial; msg[0]=(uint8_t)trial;
        m0[trial]=(trial&1)?12288:0;
        CHECK(!falcon_mrm_sign_seeded(sig,size,&written,expanded,el,pk,ql,
            m0,mu,msg,sizeof msg,seed,scratch+1,tl));
        CHECK(written==size);
        CHECK(!falcon_mrm_verify(sig,size,pk,ql,msg,sizeof msg,decoded,mu,scratch+1,tl));
        CHECK(!memcmp(m0,decoded,mu*sizeof *m0));
        CHECK(!falcon_mrm_sign_seeded(repeat,size,&written2,expanded,el,pk,ql,
            m0,mu,msg,sizeof msg,seed,scratch,tl));
        CHECK(written2==size && !memcmp(sig,repeat,size));
        msg[0]^=128;
        memcpy(decoded,sentinel,sizeof decoded);
        CHECK(falcon_mrm_verify(sig,size,pk,ql,msg,sizeof msg,decoded,mu,scratch,tl));
        CHECK(!memcmp(decoded,sentinel,sizeof decoded));
        msg[0]^=128;
        CHECK(falcon_mrm_verify(sig,size-1,pk,ql,msg,sizeof msg,decoded,mu,scratch,tl));
        CHECK(falcon_mrm_verify(sig,size+1,pk,ql,msg,sizeof msg,decoded,mu,scratch,tl));
        CHECK(falcon_mrm_verify(sig,size,pk,ql-1,msg,sizeof msg,decoded,mu,scratch,tl));
        CHECK(falcon_mrm_verify(sig,size,pk,ql,msg,sizeof msg,decoded,mu-1,scratch,tl));
        sig[trial%size]^=1;
        CHECK(falcon_mrm_verify(sig,size,pk,ql,msg,sizeof msg,decoded,mu,scratch,tl));
        m0[0]=12289;
        CHECK(falcon_mrm_sign_seeded(sig,size,&written,expanded,el,pk,ql,
            m0,mu,msg,sizeof msg,seed,scratch,tl)==FALCON_ERR_BADARG && !written);
        m0[0]=0;
    }
    /* Every possible one-byte public header, including degrees >10. */
    uint8_t saved=pk[0];
    for(unsigned i=0;i<256;i++) {
        if(i==saved) continue;
        pk[0]=(uint8_t)i;
        CHECK(falcon_mrm_verify(sig,size,pk,ql,msg,sizeof msg,decoded,mu,scratch,tl));
    }
    pk[0]=saved;
    /* Compare the exported modular operation with independent negacyclic multiplication. */
    uint16_t h[1024], hntt[1024], c[1024];
    int16_t s0[1024], s1[1024];
    size_t n=(size_t)1<<logn;
    for(size_t i=0;i<n;i++) { h[i]=(uint16_t)((i*31)%12289); s0[i]=(int16_t)((int)(i%13)-6); s1[i]=(int16_t)((int)(i%11)-5); }
    memcpy(hntt,h,n*sizeof *h);
    Zf(to_ntt_monty)(hntt,logn);
    Zf(mrm_recover_target)(c,s0,s1,hntt,logn);
    for(size_t i=0;i<n;i++) {
        int64_t v=s0[i];
        for(size_t j=0;j<n;j++) v+=(j<=i?1:-1)*(int64_t)s1[j]*h[(i+n-j)%n];
        v%=12289; if(v<0) v+=12289;
        CHECK(c[i]==v);
    }
    printf("Falcon-%u MRM: replay/recovery/malformed-input/unaligned-scratch/independent-ring PASS (%zu bytes)\n",1U<<logn,size);
    free(priv); free(pk); free(expanded); free(scratch); free(sig); free(repeat);
    return 0;
}
int main(void) { return run_profile(9)||run_profile(10); }
