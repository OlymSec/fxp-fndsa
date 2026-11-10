/* Vector Montgomery butterflies modulo 12289, same roots and scaling. */
#if defined(__aarch64__)
#include <arm_neon.h>
#define MQ_LANES 8
#define MQV uint16x8_t
#define MQ_LOAD(p) vld1q_u16(p)
#define MQ_STORE(p,x) vst1q_u16(p,x)
#define MQ_SET(x) vdupq_n_u16(x)
#define MQ_ADD(x,y) vaddq_u16(x,y)
#define MQ_SUB(x,y) vsubq_u16(x,y)
static MQV mqv_fix(MQV x){return vaddq_u16(x,vandq_u16(vreinterpretq_u16_s16(vshrq_n_s16(vreinterpretq_s16_u16(x),15)),MQ_SET(Q)));}
static MQV mqv_mul(MQV a,MQV b){
 uint32x4_t lo=vmull_u16(vget_low_u16(a),vget_low_u16(b));
 uint32x4_t hi=vmull_high_u16(a,b);
 MQV m=vmulq_u16(vmulq_u16(a,b),MQ_SET(Q0I));
 lo=vmlal_n_u16(lo,vget_low_u16(m),Q);hi=vmlal_n_u16(hi,vget_high_u16(m),Q);
 MQV x=vcombine_u16(vshrn_n_u32(lo,16),vshrn_n_u32(hi,16));
 return mqv_fix(MQ_SUB(x,MQ_SET(Q)));
}
#elif defined(__AVX2__)
#include <immintrin.h>
#define MQ_LANES 16
#define MQV __m256i
#define MQ_LOAD(p) _mm256_loadu_si256((const __m256i *)(p))
#define MQ_STORE(p,x) _mm256_storeu_si256((__m256i *)(p),x)
#define MQ_SET(x) _mm256_set1_epi16(x)
#define MQ_ADD(x,y) _mm256_add_epi16(x,y)
#define MQ_SUB(x,y) _mm256_sub_epi16(x,y)
static MQV mqv_fix(MQV x){return MQ_ADD(x,_mm256_and_si256(_mm256_srai_epi16(x,15),MQ_SET(Q)));}
static MQV mqv_mul(MQV a,MQV b){
 MQV low=_mm256_mullo_epi16(a,b),m=_mm256_mullo_epi16(low,MQ_SET(Q0I));
 MQV carry=_mm256_andnot_si256(_mm256_cmpeq_epi16(low,MQ_SET(0)),MQ_SET(1));
 MQV x=MQ_ADD(MQ_ADD(_mm256_mulhi_epu16(a,b),_mm256_mulhi_epu16(m,MQ_SET(Q))),carry);
 return mqv_fix(MQ_SUB(x,MQ_SET(Q)));
}
#endif
#ifdef MQ_LANES
static MQV mqv_add(MQV a,MQV b){return mqv_fix(MQ_SUB(MQ_ADD(a,b),MQ_SET(Q)));}
static MQV mqv_sub(MQV a,MQV b){return mqv_fix(MQ_SUB(a,b));}
#endif
