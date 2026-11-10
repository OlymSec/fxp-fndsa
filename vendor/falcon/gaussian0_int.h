#ifndef FALCON_GAUSSIAN0_INT_H
#define FALCON_GAUSSIAN0_INT_H
#include <stdint.h>
#if defined(__aarch64__)
#include <arm_neon.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#endif
/* The original 72-bit thresholds, low/middle/high 24-bit limbs. */
static const uint32_t falcon_gaussian0_thresholds[3][24] = {
    {3741698U,8248194U,2736639U,10046180U,4136815U,7650655U,7826148U,11363290U,8086568U,265321U,13644283U,9111839U,6138264U,12545723U,3104126U,28824U,198U,1U,0U,0U,0U,0U,0U,0U},
    {3068844U,1580863U,13669192U,4421575U,7122675U,13063405U,14505003U,16768101U,8444042U,12844466U,1232676U,38047U,870U,14U,0U,0U,0U,0U,0U,0U,0U,0U,0U,0U},
    {10745844U,5559083U,2260429U,708981U,169348U,30538U,4132U,417U,31U,1U,0U,0U,0U,0U,0U,0U,0U,0U,0U,0U,0U,0U,0U,0U},
};
static inline int falcon_gaussian0_rank(uint64_t lo,unsigned hi) {
    uint32_t a=(uint32_t)lo&0xffffffU,b=(uint32_t)(lo>>24)&0xffffffU,c=(uint32_t)(lo>>48)|(hi<<16);
#if defined(__aarch64__)
    uint32x4_t x=vdupq_n_u32(a),y=vdupq_n_u32(b),z=vdupq_n_u32(c),sum=vdupq_n_u32(0);
    for(unsigned i=0;i<24;i+=4) {
        uint32x4_t u=vld1q_u32(falcon_gaussian0_thresholds[0]+i),v=vld1q_u32(falcon_gaussian0_thresholds[1]+i),w=vld1q_u32(falcon_gaussian0_thresholds[2]+i);
        uint32x4_t lower=vorrq_u32(vcgtq_u32(v,y),vandq_u32(vceqq_u32(v,y),vcgtq_u32(u,x)));
        sum=vsubq_u32(sum,vorrq_u32(vcgtq_u32(w,z),vandq_u32(vceqq_u32(w,z),lower)));
    }
    return (int)vaddvq_u32(sum);
#elif defined(__AVX2__)
    __m256i x=_mm256_set1_epi32((int)a),y=_mm256_set1_epi32((int)b),z=_mm256_set1_epi32((int)c),sum=_mm256_setzero_si256();
    /* All compared limbs are <2^24, so signed comparisons are exact. */
    for(unsigned i=0;i<24;i+=8) {
        __m256i u=_mm256_loadu_si256((const __m256i *)(falcon_gaussian0_thresholds[0]+i)),v=_mm256_loadu_si256((const __m256i *)(falcon_gaussian0_thresholds[1]+i)),w=_mm256_loadu_si256((const __m256i *)(falcon_gaussian0_thresholds[2]+i));
        __m256i lower=_mm256_or_si256(_mm256_cmpgt_epi32(v,y),_mm256_and_si256(_mm256_cmpeq_epi32(v,y),_mm256_cmpgt_epi32(u,x)));
        sum=_mm256_sub_epi32(sum,_mm256_or_si256(_mm256_cmpgt_epi32(w,z),_mm256_and_si256(_mm256_cmpeq_epi32(w,z),lower)));
    }
    __m128i t=_mm_add_epi32(_mm256_castsi256_si128(sum),_mm256_extracti128_si256(sum,1));
    t=_mm_add_epi32(t,_mm_srli_si128(t,8));t=_mm_add_epi32(t,_mm_srli_si128(t,4));return _mm_cvtsi128_si32(t);
#else
    unsigned sum=0;
    for(unsigned i=0;i<18;i++) {
        uint32_t borrow=(a-falcon_gaussian0_thresholds[0][i])>>31;
        borrow=(b-falcon_gaussian0_thresholds[1][i]-borrow)>>31;
        sum+=(c-falcon_gaussian0_thresholds[2][i]-borrow)>>31;
    }
    return (int)sum;
#endif
}
#endif
