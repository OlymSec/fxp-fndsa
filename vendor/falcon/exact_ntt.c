#define _POSIX_C_SOURCE 200809L
/* Exact fixed-schedule kernels; Shoup AVX2 and Montgomery NEON butterflies. */
#ifndef FALCON_REC_SHOUP_NTT
#if defined(__AVX2__)
#define FALCON_REC_SHOUP_NTT 1
#else
#define FALCON_REC_SHOUP_NTT 0
#endif
#endif
#include "exact_ntt.h"
#include <string.h>
#if defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
#include <arm_neon.h>
#endif
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
#include <immintrin.h>
#endif
#ifndef FALCON_REC_NTT_CT
#if (defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD) && FALCON_REC_SHOUP_NTT) || (defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD))
#define FALCON_REC_NTT_CT 1
#else
#define FALCON_REC_NTT_CT 0
#endif
#endif
uint32_t falcon_rec_mod(falcon_rec_i128 a,uint32_t p) {
    /* Fixed four-word reduction; p is a public odd 31-bit prime. */
    falcon_rec_u128 u=(falcon_rec_u128)a,negative=u>>127,mask=(falcon_rec_u128)0-negative;
    falcon_rec_u128 mag=(u^mask)+negative;
    uint64_t reciprocal=UINT64_MAX/p,r=0;
    for(unsigned j=4;j-->0;) {
        uint64_t x=(r<<32)|(uint32_t)(mag>>(32*j));
        uint64_t q=(uint64_t)(((falcon_rec_u128)x*reciprocal)>>64);
        r=x-q*p;r-=p & (0U-(uint32_t)(r>=p));
    }
    uint32_t neg=(p-(uint32_t)r)&(0U-(uint32_t)(r!=0));
    return (neg&(uint32_t)mask)|((uint32_t)r&~(uint32_t)mask);
}
uint32_t falcon_rec_add(const falcon_rec_ntt_plan *c,uint32_t a,uint32_t b) {
    uint32_t r=a+b; return r-(c->p & (0U-(uint32_t)(r>=c->p)));
}
uint32_t falcon_rec_sub(const falcon_rec_ntt_plan *c,uint32_t a,uint32_t b) {
    return a-b+(c->p & (0U-(uint32_t)(a<b)));
}
uint32_t falcon_rec_mul(const falcon_rec_ntt_plan *c,uint32_t a,uint32_t b) {
    uint64_t z=(uint64_t)a*b;
    uint64_t q=(uint64_t)(((falcon_rec_u128)z*c->reciprocal)>>64);
    uint64_t r=z-q*c->p;
    return (uint32_t)(r-(c->p & (0U-(uint32_t)(r>=c->p))));
}
static uint32_t power(const falcon_rec_ntt_plan *c,uint32_t a,uint32_t e) {
    uint32_t r=1;
    while(e) {if(e&1) r=falcon_rec_mul(c,r,a);a=falcon_rec_mul(c,a,a);e>>=1;}
    return r;
}
int falcon_rec_ntt_init(falcon_rec_ntt_plan *c,unsigned n,uint32_t p) {
    if((n!=512 && n!=768 && n!=1024) || (p-1)%(n==768?2304:2*n)) return FALCON_REC_UNSUPPORTED;
    memset(c,0,sizeof *c);c->n=n;c->p=p;
    c->reciprocal=(uint64_t)(((falcon_rec_u128)1<<64)/p);
    if(n==768) {
        uint32_t psi=0;
        for(uint32_t g=2;g<1000;g++) {
            psi=power(c,g,(p-1)/2304);
            if(power(c,psi,1152)!=1 && power(c,psi,768)!=1)break;
        }
        if(power(c,psi,1152)==1 || power(c,psi,768)==1)return FALCON_REC_INTERNAL;
        uint32_t ipsi=power(c,psi,p-2),omega=power(c,psi,6),iomega=power(c,omega,p-2);
        c->rho=power(c,psi,384);c->irho=power(c,c->rho,p-2);
        c->rho_inverse_difference=power(c,falcon_rec_sub(c,c->rho,c->irho),p-2);
        c->inv_n=power(c,384,p-2);uint32_t inv=p;
        for(unsigned i=0;i<5;i++)inv*=2U-p*inv;
        c->p0i=0U-inv;uint32_t mont=(uint32_t)((UINT64_C(1)<<32)%p);
        c->roots[0]=c->iroots[0]=1;
        for(unsigned i=1;i<384;i++){c->roots[i]=falcon_rec_mul(c,c->roots[i-1],omega);c->iroots[i]=falcon_rec_mul(c,c->iroots[i-1],iomega);}
        for(unsigned i=0;i<128;i++){unsigned r=0;for(unsigned j=0;j<7;j++)r=(r<<1)|((i>>j)&1);c->reverse[i]=(uint16_t)r;}
        for(unsigned len=2;len<=128;len*=2)for(unsigned k=0;k<len/2;k++){
            c->stages[len/2+k]=c->roots[k*(384/len)];
            c->istages[len/2+k]=c->iroots[k*(384/len)];
            c->shoup[len/2+k]=(uint32_t)(((uint64_t)c->stages[len/2+k]<<32)/p);
            c->ishoup[len/2+k]=(uint32_t)(((uint64_t)c->istages[len/2+k]<<32)/p);
#if !FALCON_REC_SHOUP_NTT
            c->stages[len/2+k]=falcon_rec_mul(c,c->stages[len/2+k],mont);
            c->istages[len/2+k]=falcon_rec_mul(c,c->istages[len/2+k],mont);
#endif
        }
        uint32_t pos=1,neg=1;
        for(unsigned i=0;i<384;i++){
            c->twist[i]=falcon_rec_mul(c,pos,mont);c->twist[i+384]=falcon_rec_mul(c,neg,mont);
            c->untwist[i]=falcon_rec_mul(c,falcon_rec_mul(c,neg,c->inv_n),mont);
            c->untwist[i+384]=falcon_rec_mul(c,falcon_rec_mul(c,pos,c->inv_n),mont);
            pos=falcon_rec_mul(c,pos,psi);neg=falcon_rec_mul(c,neg,ipsi);
            c->ct_w[i]=falcon_rec_mul(c,c->roots[i],mont);
            c->gs_w[i]=falcon_rec_mul(c,c->iroots[i],mont);
        }
        c->ct_f[0]=c->ct_w[128];c->gs_f[0]=c->gs_w[128];
        c->ct_f[1]=falcon_rec_mul(c,c->rho,mont);
        c->ct_f[2]=falcon_rec_mul(c,c->irho,mont);
        c->ct_f[3]=falcon_rec_mul(c,c->rho_inverse_difference,mont);
        return 0;
    }
    uint32_t psi=0;
    for(uint32_t g=2;g<1000;g++) {
        psi=power(c,g,(p-1)/(2*n));
        if(power(c,psi,n)==p-1) break;
    }
    if(power(c,psi,n)!=p-1) return FALCON_REC_INTERNAL;
    uint32_t ipsi=power(c,psi,p-2),omega=falcon_rec_mul(c,psi,psi),iomega=falcon_rec_mul(c,ipsi,ipsi);
    c->inv_n=power(c,n,p-2);
    uint32_t inv=p;for(unsigned i=0;i<5;i++)inv*=2U-p*inv;c->p0i=0U-inv;
    c->twist[0]=c->untwist[0]=c->roots[0]=c->iroots[0]=1;
    unsigned bits=n==512?9:10;
    for(unsigned i=0;i<n;i++) {
        unsigned r=0;for(unsigned j=0;j<bits;j++)r=(r<<1)|((i>>j)&1);
        c->reverse[i]=(uint16_t)r;
        if(i) {
            c->twist[i]=falcon_rec_mul(c,c->twist[i-1],psi);
            c->untwist[i]=falcon_rec_mul(c,c->untwist[i-1],ipsi);
            c->roots[i]=falcon_rec_mul(c,c->roots[i-1],omega);
            c->iroots[i]=falcon_rec_mul(c,c->iroots[i-1],iomega);
        }
    }
    uint32_t mont=(uint32_t)((UINT64_C(1)<<32)%p);
    for(unsigned len=2;len<=n;len*=2)for(unsigned k=0;k<len/2;k++) {
        c->stages[len/2+k]=c->roots[k*(n/len)];
        c->istages[len/2+k]=c->iroots[k*(n/len)];
        c->shoup[len/2+k]=(uint32_t)(((uint64_t)c->stages[len/2+k]<<32)/p);
        c->ishoup[len/2+k]=(uint32_t)(((uint64_t)c->istages[len/2+k]<<32)/p);
#if !FALCON_REC_SHOUP_NTT
        c->stages[len/2+k]=falcon_rec_mul(c,c->stages[len/2+k],mont);
        c->istages[len/2+k]=falcon_rec_mul(c,c->istages[len/2+k],mont);
#endif
    }
    for(unsigned i=0;i<n;i++) {
        c->twist[i]=falcon_rec_mul(c,c->twist[i],mont);
        c->untwist[i]=falcon_rec_mul(c,falcon_rec_mul(c,c->untwist[i],c->inv_n),mont);
    }
    /* Negacyclic Cooley-Tukey/Gentleman-Sande twiddles: block k of every layer
     * uses psi^brv(k); spectra are produced in bit-reversed evaluation order. */
    for(unsigned k=1;k<n;k++) {
        unsigned r=0;for(unsigned j=0;j<bits;j++)r=(r<<1)|((k>>j)&1);
        c->ct_w[k]=power(c,psi,r);c->gs_w[k]=power(c,ipsi,r);
        c->ct_f[k]=(uint32_t)(((uint64_t)c->ct_w[k]<<32)/p);
        c->gs_f[k]=(uint32_t)(((uint64_t)c->gs_w[k]<<32)/p);
    }
#if defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
    /* NEON: chunks of 8, u-lane positions for the in-register len=2,1 layers. */
    enum{FALCON_REC_LAYERS=2,FALCON_REC_CHUNK=8,FALCON_REC_LANES=4,FALCON_REC_TOP=2};
    static const uint8_t lanes[2][8]={{0,1,4,5},{0,4,2,6}};
#else
    /* AVX2: chunks of 16, u-lane positions for the in-register len=4,2,1 layers. */
    enum{FALCON_REC_LAYERS=3,FALCON_REC_CHUNK=16,FALCON_REC_LANES=8,FALCON_REC_TOP=4};
    static const uint8_t lanes[3][8]={{0,1,2,3,8,9,10,11},{0,1,4,5,8,9,12,13},{0,4,2,6,8,12,10,14}};
#endif
    for(unsigned layer=0;layer<FALCON_REC_LAYERS;layer++) {
        unsigned len=FALCON_REC_TOP>>layer;
        for(unsigned chunk=0;chunk<n/FALCON_REC_CHUNK;chunk++)for(unsigned l=0;l<FALCON_REC_LANES;l++) {
            unsigned k=(n+FALCON_REC_CHUNK*chunk+lanes[layer][l])/(2*len),at=FALCON_REC_LANES*chunk+l;
            c->ct_lw[layer][at]=c->ct_w[k];c->ct_lf[layer][at]=c->ct_f[k];
            c->gs_lw[layer][at]=c->gs_w[k];c->gs_lf[layer][at]=c->gs_f[k];
        }
    }
    c->gs_ninv=c->inv_n;c->gs_ninv_f=(uint32_t)(((uint64_t)c->inv_n<<32)/p);
    c->gs_last=falcon_rec_mul(c,c->gs_w[1],c->inv_n);c->gs_last_f=(uint32_t)(((uint64_t)c->gs_last<<32)/p);
    return 0;
}
static inline uint32_t montmul(const falcon_rec_ntt_plan *c,uint32_t a,uint32_t b) {
    uint64_t z=(uint64_t)a*b;uint32_t m=(uint32_t)z*c->p0i;
    uint32_t r=(uint32_t)((z+(uint64_t)m*c->p)>>32);
    return r-(c->p & (0U-(uint32_t)(r>=c->p)));
}
#if defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
static inline uint32x4_t mont4(const falcon_rec_ntt_plan *c,uint32x4_t a,uint32x4_t b) {
    uint32x4_t mm=vmulq_n_u32(vmulq_u32(a,b),c->p0i);
    uint64x2_t lo=vmull_u32(vget_low_u32(a),vget_low_u32(b));
    uint64x2_t hi=vmull_u32(vget_high_u32(a),vget_high_u32(b));
    lo=vmlal_n_u32(lo,vget_low_u32(mm),c->p);hi=vmlal_n_u32(hi,vget_high_u32(mm),c->p);
    uint32x4_t r=vcombine_u32(vshrn_n_u64(lo,32),vshrn_n_u64(hi,32)),p=vdupq_n_u32(c->p);
    return vsubq_u32(r,vandq_u32(vcgeq_u32(r,p),p));
}
#endif
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
static inline __m256i lt8(__m256i a,__m256i b) {
    __m256i sign=_mm256_set1_epi32((int)0x80000000U);
    return _mm256_cmpgt_epi32(_mm256_xor_si256(b,sign),_mm256_xor_si256(a,sign));
}
static inline __m256i mont8(const falcon_rec_ntt_plan *c,__m256i a,__m256i b) {
    __m256i p=_mm256_set1_epi32((int)c->p);
    __m256i m=_mm256_mullo_epi32(_mm256_mullo_epi32(a,b),_mm256_set1_epi32((int)c->p0i));
    __m256i lo=_mm256_add_epi64(_mm256_mul_epu32(a,b),_mm256_mul_epu32(m,p));
    __m256i hi=_mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(a,32),_mm256_srli_epi64(b,32)),_mm256_mul_epu32(_mm256_srli_epi64(m,32),p));
    __m256i r=_mm256_blend_epi32(_mm256_srli_epi64(lo,32),hi,0xaa);
    return _mm256_sub_epi32(r,_mm256_andnot_si256(lt8(r,p),p));
}
#endif
/* Mixed-radix operations on canonical residues. The public root factors
 * are Montgomery-scaled in the otherwise unused mixed-profile CT cache. */
#if defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
typedef uint32x4_t falcon_rec_vec;
#define FALCON_REC_VL 4U
static inline falcon_rec_vec mv_load(const uint32_t *p){return vld1q_u32(p);}
static inline void mv_store(uint32_t *p,falcon_rec_vec x){vst1q_u32(p,x);}
static inline falcon_rec_vec mv_set(uint32_t x){return vdupq_n_u32(x);}
static inline falcon_rec_vec mv_add(const falcon_rec_ntt_plan *c,falcon_rec_vec a,falcon_rec_vec b){falcon_rec_vec s=vaddq_u32(a,b),p=mv_set(c->p);return vsubq_u32(s,vandq_u32(vcgeq_u32(s,p),p));}
static inline falcon_rec_vec mv_sub(const falcon_rec_ntt_plan *c,falcon_rec_vec a,falcon_rec_vec b){return vaddq_u32(vsubq_u32(a,b),vandq_u32(vcltq_u32(a,b),mv_set(c->p)));}
static inline falcon_rec_vec mv_mul(const falcon_rec_ntt_plan *c,falcon_rec_vec a,falcon_rec_vec b){return mont4(c,a,b);}
static inline falcon_rec_vec mv_roots2(const uint32_t *p){return (falcon_rec_vec){p[0],p[2],p[4],p[6]};}
#elif defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
typedef __m256i falcon_rec_vec;
#define FALCON_REC_VL 8U
static inline falcon_rec_vec mv_load(const uint32_t *p){return _mm256_loadu_si256((const __m256i *)p);}
static inline void mv_store(uint32_t *p,falcon_rec_vec x){_mm256_storeu_si256((__m256i *)p,x);}
static inline falcon_rec_vec mv_set(uint32_t x){return _mm256_set1_epi32((int)x);}
static inline falcon_rec_vec mv_add(const falcon_rec_ntt_plan *c,falcon_rec_vec a,falcon_rec_vec b){falcon_rec_vec s=_mm256_add_epi32(a,b),p=mv_set(c->p);return _mm256_sub_epi32(s,_mm256_andnot_si256(lt8(s,p),p));}
static inline falcon_rec_vec mv_sub(const falcon_rec_ntt_plan *c,falcon_rec_vec a,falcon_rec_vec b){return _mm256_add_epi32(_mm256_sub_epi32(a,b),_mm256_and_si256(lt8(a,b),mv_set(c->p)));}
static inline falcon_rec_vec mv_mul(const falcon_rec_ntt_plan *c,falcon_rec_vec a,falcon_rec_vec b){return mont8(c,a,b);}
static inline falcon_rec_vec mv_roots2(const uint32_t *p){return _mm256_setr_epi32((int)p[0],(int)p[2],(int)p[4],(int)p[6],(int)p[8],(int)p[10],(int)p[12],(int)p[14]);}
#endif
/* All cached keys are Montgomery-scaled; live spectra remain canonical.
 * Fused public expressions avoid per-coefficient external function calls. */
void falcon_rec_linear(const falcon_rec_ntt_plan *restrict p,uint32_t *restrict out,
    const uint32_t *a,const uint32_t *A,const uint32_t *b,const uint32_t *B,
    const uint32_t *c,const uint32_t *C,const uint32_t *d,const uint32_t *D) {
    unsigned i=0;
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
    __m256i pp=_mm256_set1_epi32((int)p->p);
#define FALCON_REC_LOAD(x) _mm256_loadu_si256((const __m256i *)((x)+i))
#define FALCON_REC_SUB(v,x) do{__m256i x_=(x);(v)=_mm256_add_epi32(_mm256_sub_epi32((v),x_),_mm256_and_si256(lt8((v),x_),pp));}while(0)
    for(;i+8<=p->n;i+=8) {
        __m256i v=a?mont8(p,FALCON_REC_LOAD(a),FALCON_REC_LOAD(A)):_mm256_setzero_si256();
        FALCON_REC_SUB(v,mont8(p,FALCON_REC_LOAD(b),FALCON_REC_LOAD(B)));FALCON_REC_SUB(v,mont8(p,FALCON_REC_LOAD(c),FALCON_REC_LOAD(C)));
        if(d)FALCON_REC_SUB(v,mont8(p,FALCON_REC_LOAD(d),FALCON_REC_LOAD(D)));
        _mm256_storeu_si256((__m256i *)(out+i),v);
    }
#undef FALCON_REC_LOAD
#undef FALCON_REC_SUB
#elif defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
    uint32x4_t pp=vdupq_n_u32(p->p);
#define FALCON_REC_SUB(v,x) do{uint32x4_t x_=(x);(v)=vaddq_u32(vsubq_u32((v),x_),vandq_u32(vcltq_u32((v),x_),pp));}while(0)
    for(;i+4<=p->n;i+=4) {
        uint32x4_t v=a?mont4(p,vld1q_u32(a+i),vld1q_u32(A+i)):vdupq_n_u32(0);
        FALCON_REC_SUB(v,mont4(p,vld1q_u32(b+i),vld1q_u32(B+i)));FALCON_REC_SUB(v,mont4(p,vld1q_u32(c+i),vld1q_u32(C+i)));
        if(d)FALCON_REC_SUB(v,mont4(p,vld1q_u32(d+i),vld1q_u32(D+i)));
        vst1q_u32(out+i,v);
    }
#undef FALCON_REC_SUB
#endif
    for(;i<p->n;i++) {
        uint32_t v=a?montmul(p,a[i],A[i]):0;
        v=falcon_rec_sub(p,v,montmul(p,b[i],B[i]));v=falcon_rec_sub(p,v,montmul(p,c[i],C[i]));
        if(d)v=falcon_rec_sub(p,v,montmul(p,d[i],D[i]));
        out[i]=v;
    }
}
void falcon_rec_response(const falcon_rec_ntt_plan *restrict p,uint32_t *restrict out,const uint32_t *a,const uint32_t *A,
                  const uint32_t *b,const uint32_t *B,const uint32_t *target) {
    unsigned i=0;
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
    __m256i pp=_mm256_set1_epi32((int)p->p);
#define FALCON_REC_LOAD(x) _mm256_loadu_si256((const __m256i *)((x)+i))
    for(;i+8<=p->n;i+=8) {
        __m256i v=_mm256_add_epi32(mont8(p,FALCON_REC_LOAD(a),FALCON_REC_LOAD(A)),mont8(p,FALCON_REC_LOAD(b),FALCON_REC_LOAD(B)));
        v=_mm256_sub_epi32(v,_mm256_andnot_si256(lt8(v,pp),pp));
        if(target){__m256i t=FALCON_REC_LOAD(target);v=_mm256_add_epi32(_mm256_sub_epi32(t,v),_mm256_and_si256(lt8(t,v),pp));}
        _mm256_storeu_si256((__m256i *)(out+i),v);
    }
#undef FALCON_REC_LOAD
#elif defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
    uint32x4_t pp=vdupq_n_u32(p->p);
    for(;i+4<=p->n;i+=4) {
        uint32x4_t v=vaddq_u32(mont4(p,vld1q_u32(a+i),vld1q_u32(A+i)),mont4(p,vld1q_u32(b+i),vld1q_u32(B+i)));
        v=vsubq_u32(v,vandq_u32(vcgeq_u32(v,pp),pp));
        if(target){uint32x4_t t=vld1q_u32(target+i);v=vaddq_u32(vsubq_u32(t,v),vandq_u32(vcltq_u32(t,v),pp));}
        vst1q_u32(out+i,v);
    }
#endif
    for(;i<p->n;i++) {uint32_t v=falcon_rec_add(p,montmul(p,a[i],A[i]),montmul(p,b[i],B[i]));out[i]=target?falcon_rec_sub(p,target[i],v):v;}
}
static inline uint32_t shoup1(const falcon_rec_ntt_plan *p,uint32_t a,uint32_t b,uint32_t factor) {
    uint32_t q=(uint32_t)(((uint64_t)a*factor)>>32),r=a*b-q*p->p;
    return r-(p->p & (0U-(uint32_t)(r>=p->p)));
}
#if defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
static inline uint32x4_t shoup4(const falcon_rec_ntt_plan *p,uint32x4_t a,uint32x4_t b,uint32x4_t factor) {
    uint64x2_t lo=vmull_u32(vget_low_u32(a),vget_low_u32(factor)),hi=vmull_u32(vget_high_u32(a),vget_high_u32(factor));
    uint32x4_t q=vcombine_u32(vshrn_n_u64(lo,32),vshrn_n_u64(hi,32)),pp=vdupq_n_u32(p->p);
    uint32x4_t r=vsubq_u32(vmulq_u32(a,b),vmulq_u32(q,pp));
    return vsubq_u32(r,vandq_u32(vcgeq_u32(r,pp),pp));
}
#endif
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
static inline __m256i shoup8(const falcon_rec_ntt_plan *p,__m256i a,__m256i b,__m256i factor) {
    __m256i lo=_mm256_mul_epu32(a,factor),hi=_mm256_mul_epu32(_mm256_srli_epi64(a,32),_mm256_srli_epi64(factor,32));
    __m256i q=_mm256_blend_epi32(_mm256_srli_epi64(lo,32),hi,0xaa),pp=_mm256_set1_epi32((int)p->p);
    __m256i r=_mm256_sub_epi32(_mm256_mullo_epi32(a,b),_mm256_mullo_epi32(q,pp));
    return _mm256_sub_epi32(r,_mm256_andnot_si256(lt8(r,pp),pp));
}
#endif
static void radix2(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a,unsigned n,int inverse) {
    for(unsigned i=0;i<n;i++)if(i<c->reverse[i]) {
        uint32_t t=a[i];a[i]=a[c->reverse[i]];a[c->reverse[i]]=t;
    }
    const uint32_t *roots=inverse?c->istages:c->stages;
#if FALCON_REC_SHOUP_NTT
    const uint32_t *factors=inverse?c->ishoup:c->shoup;
#endif
    unsigned first_length=2;
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD) && FALCON_REC_SHOUP_NTT
    __m256i modulus=_mm256_set1_epi32((int)c->p);
    __m256i w4=_mm256_setr_epi32((int)roots[2],(int)roots[3],(int)roots[2],(int)roots[3],(int)roots[2],(int)roots[3],(int)roots[2],(int)roots[3]);
    __m256i f4=_mm256_setr_epi32((int)factors[2],(int)factors[3],(int)factors[2],(int)factors[3],(int)factors[2],(int)factors[3],(int)factors[2],(int)factors[3]);
    __m256i w8=_mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(roots+4)));
    __m256i f8=_mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(factors+4)));
    __m256i evens=_mm256_setr_epi32(0,2,4,6,0,2,4,6),odds=_mm256_setr_epi32(1,3,5,7,1,3,5,7);
    __m256i lowpairs=_mm256_setr_epi32(0,1,4,5,0,1,4,5),hipairs=_mm256_setr_epi32(2,3,6,7,2,3,6,7);
#define FALCON_REC_SUMDIFF() do { \
        sum=_mm256_add_epi32(u,v);dif=_mm256_sub_epi32(u,v); \
        sum=_mm256_sub_epi32(sum,_mm256_andnot_si256(lt8(sum,modulus),modulus)); \
        dif=_mm256_add_epi32(dif,_mm256_and_si256(lt8(u,v),modulus)); \
    }while(0)
    for(unsigned start=0;start<n;start+=16) {
        __m256i x=_mm256_loadu_si256((const __m256i *)(a+start)),y=_mm256_loadu_si256((const __m256i *)(a+start+8));
        __m256i u=_mm256_permute2x128_si256(_mm256_permutevar8x32_epi32(x,evens),_mm256_permutevar8x32_epi32(y,evens),0x20);
        __m256i v=_mm256_permute2x128_si256(_mm256_permutevar8x32_epi32(x,odds),_mm256_permutevar8x32_epi32(y,odds),0x20),sum,dif;
        FALCON_REC_SUMDIFF();
        __m256i lo=_mm256_unpacklo_epi32(sum,dif),hi=_mm256_unpackhi_epi32(sum,dif);
        x=_mm256_permute2x128_si256(lo,hi,0x20);y=_mm256_permute2x128_si256(lo,hi,0x31);
        u=_mm256_permute2x128_si256(_mm256_permutevar8x32_epi32(x,lowpairs),_mm256_permutevar8x32_epi32(y,lowpairs),0x20);
        v=_mm256_permute2x128_si256(_mm256_permutevar8x32_epi32(x,hipairs),_mm256_permutevar8x32_epi32(y,hipairs),0x20);
        v=shoup8(c,v,w4,f4);FALCON_REC_SUMDIFF();
        lo=_mm256_unpacklo_epi64(sum,dif);hi=_mm256_unpackhi_epi64(sum,dif);
        x=_mm256_permute2x128_si256(lo,hi,0x20);y=_mm256_permute2x128_si256(lo,hi,0x31);
        u=_mm256_permute2x128_si256(x,y,0x20);v=shoup8(c,_mm256_permute2x128_si256(x,y,0x31),w8,f8);
        FALCON_REC_SUMDIFF();
        _mm256_storeu_si256((__m256i *)(a+start),_mm256_permute2x128_si256(sum,dif,0x20));
        _mm256_storeu_si256((__m256i *)(a+start+8),_mm256_permute2x128_si256(sum,dif,0x31));
    }
#undef FALCON_REC_SUMDIFF
    first_length=16;
#endif
    for(unsigned len=first_length;len<=n;len*=2)for(unsigned start=0;start<n;start+=len) {
        unsigned k=0;
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
        __m256i p=_mm256_set1_epi32((int)c->p);
        for(;k+8<=len/2;k+=8) {
            __m256i u=_mm256_loadu_si256((const __m256i *)(a+start+k));
#if FALCON_REC_SHOUP_NTT
            __m256i v=shoup8(c,_mm256_loadu_si256((const __m256i *)(a+start+k+len/2)),_mm256_loadu_si256((const __m256i *)(roots+len/2+k)),_mm256_loadu_si256((const __m256i *)(factors+len/2+k)));
#else
            __m256i v=mont8(c,_mm256_loadu_si256((const __m256i *)(a+start+k+len/2)),_mm256_loadu_si256((const __m256i *)(roots+len/2+k)));
#endif
            __m256i sum=_mm256_add_epi32(u,v),dif=_mm256_sub_epi32(u,v);
            sum=_mm256_sub_epi32(sum,_mm256_andnot_si256(lt8(sum,p),p));
            dif=_mm256_add_epi32(dif,_mm256_and_si256(lt8(u,v),p));
            _mm256_storeu_si256((__m256i *)(a+start+k),sum);_mm256_storeu_si256((__m256i *)(a+start+k+len/2),dif);
        }
#elif defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
        uint32x4_t p=vdupq_n_u32(c->p);
        for(;k+4<=len/2;k+=4) {
            uint32x4_t u=vld1q_u32(a+start+k);
#if FALCON_REC_SHOUP_NTT
            uint32x4_t v=shoup4(c,vld1q_u32(a+start+k+len/2),vld1q_u32(roots+len/2+k),vld1q_u32(factors+len/2+k));
#else
            uint32x4_t v=mont4(c,vld1q_u32(a+start+k+len/2),vld1q_u32(roots+len/2+k));
#endif
            uint32x4_t sum=vaddq_u32(u,v),dif=vsubq_u32(u,v);
            sum=vsubq_u32(sum,vandq_u32(vcgeq_u32(sum,p),p));
            dif=vaddq_u32(dif,vandq_u32(vcltq_u32(u,v),p));
            vst1q_u32(a+start+k,sum);vst1q_u32(a+start+k+len/2,dif);
        }
#endif
        for(;k<len/2;k++) {
            uint32_t u=a[start+k];
#if FALCON_REC_SHOUP_NTT
            uint32_t v=shoup1(c,a[start+k+len/2],roots[len/2+k],factors[len/2+k]);
#else
            uint32_t v=montmul(c,a[start+k+len/2],roots[len/2+k]);
#endif
            a[start+k]=falcon_rec_add(c,u,v);a[start+k+len/2]=falcon_rec_sub(c,u,v);
        }
    }

}
#include "exact_ntt_pair.inc"
#if FALCON_REC_NTT_CT && defined(__AVX2__)
/* Butterflies on canonical residues in [0,p), p<2^31; each output is canonical. */
static inline __m256i ct_add(__m256i a,__m256i b,__m256i p){__m256i s=_mm256_add_epi32(a,b);return _mm256_sub_epi32(s,_mm256_andnot_si256(lt8(s,p),p));}
static inline __m256i ct_sub(__m256i a,__m256i b,__m256i p){return _mm256_add_epi32(_mm256_sub_epi32(a,b),_mm256_and_si256(lt8(a,b),p));}
#define FALCON_REC_LD(x) _mm256_loadu_si256((const __m256i *)(x))
#define FALCON_REC_ST(x,v) _mm256_storeu_si256((__m256i *)(x),(v))
#define FALCON_REC_SH(a,b) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a),_mm256_castsi256_ps(b),(int)0x88))
#define FALCON_REC_SHO(a,b) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a),_mm256_castsi256_ps(b),(int)0xdd))
/* Forward: natural coefficients in, bit-reversed spectrum out; m arrays share twiddle loads. */
static inline __attribute__((always_inline)) void ct_forward(const falcon_rec_ntt_plan *restrict c,uint32_t *const *arr,unsigned m) {
    const unsigned n=c->n;const __m256i p=_mm256_set1_epi32((int)c->p);
    unsigned k=1;   /* block k of layer len is n/(2len)+start/(2len) */
    for(unsigned len=n/2;len>=8;len>>=1)for(unsigned start=0;start<n;start+=2*len,k++) {
        __m256i w=_mm256_set1_epi32((int)c->ct_w[k]),f=_mm256_set1_epi32((int)c->ct_f[k]);
        for(unsigned j=start;j<start+len;j+=8)
#pragma GCC unroll 2
            for(unsigned t=0;t<m;t++) {
                uint32_t *a=arr[t];
                __m256i u=FALCON_REC_LD(a+j),v=shoup8(c,FALCON_REC_LD(a+j+len),w,f);
                FALCON_REC_ST(a+j,ct_add(u,v,p));FALCON_REC_ST(a+j+len,ct_sub(u,v,p));
            }
    }
    for(unsigned chunk=0;chunk<n/16;chunk++) {
        __m256i w0=FALCON_REC_LD(c->ct_lw[0]+8*chunk),f0=FALCON_REC_LD(c->ct_lf[0]+8*chunk);
        __m256i w1=FALCON_REC_LD(c->ct_lw[1]+8*chunk),f1=FALCON_REC_LD(c->ct_lf[1]+8*chunk);
        __m256i w2=FALCON_REC_LD(c->ct_lw[2]+8*chunk),f2=FALCON_REC_LD(c->ct_lf[2]+8*chunk);
#pragma GCC unroll 2
        for(unsigned t=0;t<m;t++) {
            uint32_t *a=arr[t]+16*chunk;
            __m256i x=FALCON_REC_LD(a),y=FALCON_REC_LD(a+8);
            __m256i u=_mm256_permute2x128_si256(x,y,0x20),v=_mm256_permute2x128_si256(x,y,0x31),r;
            r=shoup8(c,v,w0,f0);v=ct_sub(u,r,p);u=ct_add(u,r,p);
            x=_mm256_unpacklo_epi64(u,v);y=_mm256_unpackhi_epi64(u,v);
            r=shoup8(c,y,w1,f1);v=ct_sub(x,r,p);u=ct_add(x,r,p);
            x=FALCON_REC_SH(u,v);y=FALCON_REC_SHO(u,v);
            r=shoup8(c,y,w2,f2);v=ct_sub(x,r,p);u=ct_add(x,r,p);
            x=_mm256_unpacklo_epi32(u,v);y=_mm256_unpackhi_epi32(u,v);
            u=_mm256_unpacklo_epi64(x,y);v=_mm256_unpackhi_epi64(x,y);
            FALCON_REC_ST(a,_mm256_permute2x128_si256(u,v,0x20));FALCON_REC_ST(a+8,_mm256_permute2x128_si256(u,v,0x31));
        }
    }
}
/* Inverse: bit-reversed spectrum in, natural coefficients out, including 1/n. */
static void gs_inverse(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a) {
    const unsigned n=c->n;const __m256i p=_mm256_set1_epi32((int)c->p);
    for(unsigned chunk=0;chunk<n/16;chunk++) {
        uint32_t *b=a+16*chunk;
        __m256i x=FALCON_REC_LD(b),y=FALCON_REC_LD(b+8);
        __m256i u=_mm256_permute2x128_si256(x,y,0x20),v=_mm256_permute2x128_si256(x,y,0x31);
        x=_mm256_unpacklo_epi64(u,v);y=_mm256_unpackhi_epi64(u,v);
        u=FALCON_REC_SH(x,y);v=FALCON_REC_SHO(x,y);
        x=ct_add(u,v,p);y=shoup8(c,ct_sub(u,v,p),FALCON_REC_LD(c->gs_lw[2]+8*chunk),FALCON_REC_LD(c->gs_lf[2]+8*chunk));
        u=_mm256_unpacklo_epi32(x,y);v=_mm256_unpackhi_epi32(x,y);
        x=ct_add(u,v,p);y=shoup8(c,ct_sub(u,v,p),FALCON_REC_LD(c->gs_lw[1]+8*chunk),FALCON_REC_LD(c->gs_lf[1]+8*chunk));
        u=_mm256_unpacklo_epi64(x,y);v=_mm256_unpackhi_epi64(x,y);
        x=ct_add(u,v,p);y=shoup8(c,ct_sub(u,v,p),FALCON_REC_LD(c->gs_lw[0]+8*chunk),FALCON_REC_LD(c->gs_lf[0]+8*chunk));
        FALCON_REC_ST(b,_mm256_permute2x128_si256(x,y,0x20));FALCON_REC_ST(b+8,_mm256_permute2x128_si256(x,y,0x31));
    }
    for(unsigned len=8,base=n/16;len<n/2;len<<=1,base>>=1)for(unsigned start=0,k=base;start<n;start+=2*len,k++) {
        __m256i w=_mm256_set1_epi32((int)c->gs_w[k]),f=_mm256_set1_epi32((int)c->gs_f[k]);
        for(unsigned j=start;j<start+len;j+=8) {
            __m256i u=FALCON_REC_LD(a+j),v=FALCON_REC_LD(a+j+len);
            FALCON_REC_ST(a+j,ct_add(u,v,p));FALCON_REC_ST(a+j+len,shoup8(c,ct_sub(u,v,p),w,f));
        }
    }
    const unsigned half=n/2;
    __m256i wn=_mm256_set1_epi32((int)c->gs_ninv),fn=_mm256_set1_epi32((int)c->gs_ninv_f);
    __m256i wl=_mm256_set1_epi32((int)c->gs_last),fl=_mm256_set1_epi32((int)c->gs_last_f);
    for(unsigned j=0;j<half;j+=8) {
        __m256i u=FALCON_REC_LD(a+j),v=FALCON_REC_LD(a+j+half);
        FALCON_REC_ST(a+j,shoup8(c,ct_add(u,v,p),wn,fn));FALCON_REC_ST(a+j+half,shoup8(c,ct_sub(u,v,p),wl,fl));
    }
}
#undef FALCON_REC_LD
#undef FALCON_REC_ST
#undef FALCON_REC_SH
#undef FALCON_REC_SHO
#elif FALCON_REC_NTT_CT
/* NEON: same transform and bit-reversed spectrum order, four lanes, chunks of 8. */
static inline uint32x4_t ct_add4(uint32x4_t a,uint32x4_t b,uint32x4_t p){uint32x4_t s=vaddq_u32(a,b);return vsubq_u32(s,vandq_u32(vcgeq_u32(s,p),p));}
static inline uint32x4_t ct_sub4(uint32x4_t a,uint32x4_t b,uint32x4_t p){return vaddq_u32(vsubq_u32(a,b),vandq_u32(vcltq_u32(a,b),p));}
static inline __attribute__((always_inline)) void ct_forward(const falcon_rec_ntt_plan *restrict c,uint32_t *const *arr,unsigned m) {
    const unsigned n=c->n;const uint32x4_t p=vdupq_n_u32(c->p);
    unsigned k=1;   /* block k of layer len is n/(2len)+start/(2len) */
    for(unsigned len=n/2;len>=4;len>>=1)for(unsigned start=0;start<n;start+=2*len,k++) {
        uint32x4_t w=vdupq_n_u32(c->ct_w[k]),f=vdupq_n_u32(c->ct_f[k]);
        for(unsigned j=start;j<start+len;j+=4)
#pragma GCC unroll 2
            for(unsigned t=0;t<m;t++) {
                uint32_t *a=arr[t];
                uint32x4_t u=vld1q_u32(a+j),v=shoup4(c,vld1q_u32(a+j+len),w,f);
                vst1q_u32(a+j,ct_add4(u,v,p));vst1q_u32(a+j+len,ct_sub4(u,v,p));
            }
    }
    for(unsigned chunk=0;chunk<n/8;chunk++) {
        uint32x4_t w0=vld1q_u32(c->ct_lw[0]+4*chunk),f0=vld1q_u32(c->ct_lf[0]+4*chunk);
        uint32x4_t w1=vld1q_u32(c->ct_lw[1]+4*chunk),f1=vld1q_u32(c->ct_lf[1]+4*chunk);
#pragma GCC unroll 2
        for(unsigned t=0;t<m;t++) {
            uint32_t *a=arr[t]+8*chunk;
            uint32x4_t x=vld1q_u32(a),y=vld1q_u32(a+4);
            uint32x4_t u=vcombine_u32(vget_low_u32(x),vget_low_u32(y)),v=vcombine_u32(vget_high_u32(x),vget_high_u32(y)),r;
            r=shoup4(c,v,w0,f0);v=ct_sub4(u,r,p);u=ct_add4(u,r,p);
            x=vuzp1q_u32(u,v);y=vuzp2q_u32(u,v);
            r=shoup4(c,y,w1,f1);v=ct_sub4(x,r,p);u=ct_add4(x,r,p);
            x=vzip1q_u32(u,v);y=vzip2q_u32(u,v);
            vst1q_u32(a,vcombine_u32(vget_low_u32(x),vget_low_u32(y)));
            vst1q_u32(a+4,vcombine_u32(vget_high_u32(x),vget_high_u32(y)));
        }
    }
}
static void gs_inverse(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a) {
    const unsigned n=c->n;const uint32x4_t p=vdupq_n_u32(c->p);
    for(unsigned chunk=0;chunk<n/8;chunk++) {
        uint32_t *b=a+8*chunk;
        uint32x4_t x=vld1q_u32(b),y=vld1q_u32(b+4);
        uint32x4_t u=vcombine_u32(vget_low_u32(x),vget_low_u32(y)),v=vcombine_u32(vget_high_u32(x),vget_high_u32(y));
        x=vuzp1q_u32(u,v);y=vuzp2q_u32(u,v);
        u=ct_add4(x,y,p);v=shoup4(c,ct_sub4(x,y,p),vld1q_u32(c->gs_lw[1]+4*chunk),vld1q_u32(c->gs_lf[1]+4*chunk));
        x=vzip1q_u32(u,v);y=vzip2q_u32(u,v);
        u=ct_add4(x,y,p);v=shoup4(c,ct_sub4(x,y,p),vld1q_u32(c->gs_lw[0]+4*chunk),vld1q_u32(c->gs_lf[0]+4*chunk));
        vst1q_u32(b,vcombine_u32(vget_low_u32(u),vget_low_u32(v)));
        vst1q_u32(b+4,vcombine_u32(vget_high_u32(u),vget_high_u32(v)));
    }
    for(unsigned len=4,base=n/8;len<n/2;len<<=1,base>>=1)for(unsigned start=0,k=base;start<n;start+=2*len,k++) {
        uint32x4_t w=vdupq_n_u32(c->gs_w[k]),f=vdupq_n_u32(c->gs_f[k]);
        for(unsigned j=start;j<start+len;j+=4) {
            uint32x4_t u=vld1q_u32(a+j),v=vld1q_u32(a+j+len);
            vst1q_u32(a+j,ct_add4(u,v,p));vst1q_u32(a+j+len,shoup4(c,ct_sub4(u,v,p),w,f));
        }
    }
    const unsigned half=n/2;
    uint32x4_t wn=vdupq_n_u32(c->gs_ninv),fn=vdupq_n_u32(c->gs_ninv_f);
    uint32x4_t wl=vdupq_n_u32(c->gs_last),fl=vdupq_n_u32(c->gs_last_f);
    for(unsigned j=0;j<half;j+=4) {
        uint32x4_t u=vld1q_u32(a+j),v=vld1q_u32(a+j+half);
        vst1q_u32(a+j,shoup4(c,ct_add4(u,v,p),wn,fn));vst1q_u32(a+j+half,shoup4(c,ct_sub4(u,v,p),wl,fl));
    }
}
#endif
static void cyclic384(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a,int inverse) {
    uint32_t lanes[3][128];
#if !defined(FALCON_REC_VL)
    const uint32_t *roots=inverse?c->iroots:c->roots;
#endif
    for(unsigned r=0;r<3;r++) {
        for(unsigned i=0;i<128;i++)lanes[r][i]=a[3*i+r];
        radix2(c,lanes[r],128,inverse);
    }
#if defined(FALCON_REC_VL)
    const uint32_t *mr=inverse?c->gs_w:c->ct_w;
    falcon_rec_vec w=mv_set(inverse?c->gs_f[0]:c->ct_f[0]);
    for(unsigned k=0;k<128;k+=FALCON_REC_VL) {
        falcon_rec_vec x=mv_load(lanes[0]+k);
        falcon_rec_vec y=mv_mul(c,mv_load(lanes[1]+k),mv_load(mr+k));
        falcon_rec_vec z=mv_mul(c,mv_load(lanes[2]+k),mv_roots2(mr+2*k));
        /* w^2=-1-w: one product gives both nontrivial outputs. */
        falcon_rec_vec yz=mv_mul(c,mv_sub(c,y,z),w);
        mv_store(a+k,mv_add(c,x,mv_add(c,y,z)));
        mv_store(a+128+k,mv_add(c,mv_sub(c,x,z),yz));
        mv_store(a+256+k,mv_sub(c,mv_sub(c,x,y),yz));
    }
#else
    for(unsigned r=1;r<3;r++)for(unsigned k=0;k<128;k++)lanes[r][k]=falcon_rec_mul(c,lanes[r][k],roots[r*k]);
    for(unsigned ell=0;ell<3;ell++)for(unsigned k=0;k<128;k++) {
        uint32_t v=lanes[0][k];
        for(unsigned r=1;r<3;r++)v=falcon_rec_add(c,v,falcon_rec_mul(c,lanes[r][k],roots[(128*r*ell)%384]));
        a[128*ell+k]=v;
    }
#endif
}
static void twist_words(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a,const uint32_t *restrict tw) {
    unsigned i=0;
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
    for(;i+8<=c->n;i+=8)_mm256_storeu_si256((__m256i *)(a+i),mont8(c,_mm256_loadu_si256((const __m256i *)(a+i)),_mm256_loadu_si256((const __m256i *)(tw+i))));
#elif defined(__aarch64__) && !defined(FALCON_REC_NO_SIMD)
    for(;i+4<=c->n;i+=4)vst1q_u32(a+i,mont4(c,vld1q_u32(a+i),vld1q_u32(tw+i)));
#endif
    for(;i<c->n;i++)a[i]=montmul(c,a[i],tw[i]);
}
#ifdef FALCON_REC_COMPONENT_PROFILE
static void falcon_rec_ntt_body(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a,int inverse)
#else
void falcon_rec_ntt(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a,int inverse)
#endif
{
    unsigned n=c->n;
    if(n==768) {
#if defined(FALCON_REC_VL)
        if(!inverse)for(unsigned i=0;i<384;i+=FALCON_REC_VL) {
            falcon_rec_vec lo=mv_load(a+i),hi=mv_load(a+384+i);
            falcon_rec_vec v0=mv_add(c,lo,mv_mul(c,hi,mv_set(c->ct_f[1])));
            falcon_rec_vec v1=mv_add(c,lo,mv_mul(c,hi,mv_set(c->ct_f[2])));
            mv_store(a+i,mv_mul(c,v0,mv_load(c->twist+i)));
            mv_store(a+384+i,mv_mul(c,v1,mv_load(c->twist+384+i)));
        }
#else
        if(!inverse)for(unsigned i=0;i<384;i++) {
            uint32_t lo=a[i],hi=a[i+384];
            a[i]=montmul(c,falcon_rec_add(c,lo,falcon_rec_mul(c,c->rho,hi)),c->twist[i]);
            a[i+384]=montmul(c,falcon_rec_add(c,lo,falcon_rec_mul(c,c->irho,hi)),c->twist[i+384]);
        }
#endif
        cyclic384(c,a,inverse);cyclic384(c,a+384,inverse);
#if defined(FALCON_REC_VL)
        if(inverse)for(unsigned i=0;i<384;i+=FALCON_REC_VL) {
            falcon_rec_vec v0=mv_mul(c,mv_load(a+i),mv_load(c->untwist+i));
            falcon_rec_vec v1=mv_mul(c,mv_load(a+384+i),mv_load(c->untwist+384+i));
            falcon_rec_vec hi=mv_mul(c,mv_sub(c,v0,v1),mv_set(c->ct_f[3]));
            mv_store(a+i,mv_sub(c,v0,mv_mul(c,hi,mv_set(c->ct_f[1]))));
            mv_store(a+384+i,hi);
        }
#else
        if(inverse)for(unsigned i=0;i<384;i++) {
            uint32_t v0=montmul(c,a[i],c->untwist[i]),v1=montmul(c,a[i+384],c->untwist[i+384]);
            uint32_t hi=falcon_rec_mul(c,falcon_rec_sub(c,v0,v1),c->rho_inverse_difference);
            a[i]=falcon_rec_sub(c,v0,falcon_rec_mul(c,c->rho,hi));a[i+384]=hi;
        }
#endif
        return;
    }
#if FALCON_REC_NTT_CT
    if(inverse)gs_inverse(c,a);
    else {uint32_t *arr[1]={a};ct_forward(c,arr,1);}
    return;
#endif
    if(!inverse)twist_words(c,a,c->twist);
    radix2(c,a,n,inverse);
    if(inverse)twist_words(c,a,c->untwist);
}

#ifdef FALCON_REC_COMPONENT_PROFILE
#include <time.h>
_Thread_local uint64_t falcon_rec_ntt_ns[2],falcon_rec_ntt_calls[2];
static uint64_t ntt_clock(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return (uint64_t)t.tv_sec*1000000000U+t.tv_nsec;}
void falcon_rec_ntt(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a,int inverse) {
#ifdef FALCON_REC_COMPONENT_PROFILE
    uint64_t before=ntt_clock();
#endif
    falcon_rec_ntt_body(c,a,inverse);
#ifdef FALCON_REC_COMPONENT_PROFILE
    falcon_rec_ntt_ns[inverse]+=ntt_clock()-before;falcon_rec_ntt_calls[inverse]++;
#endif
}
#endif

void falcon_rec_ntt_pair(const falcon_rec_ntt_plan *c,uint32_t *a,uint32_t *b) {
#if FALCON_REC_NTT_CT
    if(c->n!=768) {
#ifdef FALCON_REC_COMPONENT_PROFILE
        uint64_t before=ntt_clock();
#endif
        uint32_t *arr[2]={a,b};
#if FALCON_REC_SCHEDULE
        ct_forward(c,arr,2);
#else
        ct_forward(c,arr,1);ct_forward(c,arr+1,1);
#endif
#ifdef FALCON_REC_COMPONENT_PROFILE
        falcon_rec_ntt_ns[0]+=ntt_clock()-before;falcon_rec_ntt_calls[0]+=2;
#endif
        return;
    }
#endif
#if FALCON_REC_SCHEDULE && defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD) && FALCON_REC_SHOUP_NTT
    /* The paired radix-two schedule requires a power-of-two degree.
     * Mixed-ring calls use the native cubic decomposition in falcon_rec_ntt. */
    if(c->n==768){falcon_rec_ntt(c,a,0);falcon_rec_ntt(c,b,0);return;}
#ifdef FALCON_REC_COMPONENT_PROFILE
    uint64_t before=ntt_clock();
#endif
    for(unsigned i=0;i<c->n;i+=8) {
        __m256i tw=_mm256_loadu_si256((const __m256i *)(c->twist+i));
        _mm256_storeu_si256((__m256i *)(a+i),mont8(c,_mm256_loadu_si256((const __m256i *)(a+i)),tw));
        _mm256_storeu_si256((__m256i *)(b+i),mont8(c,_mm256_loadu_si256((const __m256i *)(b+i)),tw));
    }
    radix2_pair(c,a,b,c->n);
#ifdef FALCON_REC_COMPONENT_PROFILE
    falcon_rec_ntt_ns[0]+=ntt_clock()-before;falcon_rec_ntt_calls[0]+=2;
#endif
#else
    falcon_rec_ntt(c,a,0);falcon_rec_ntt(c,b,0);
#endif
}
