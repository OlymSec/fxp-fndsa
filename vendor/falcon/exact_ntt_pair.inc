/* Paired radix-two schedule: identical modular butterflies with shared public twiddle loads. */
#if FALCON_REC_SCHEDULE && defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD) && FALCON_REC_SHOUP_NTT
static void radix2_pair(const falcon_rec_ntt_plan *restrict c,uint32_t *restrict a0,uint32_t *restrict a1,unsigned n) {
    for(unsigned i=0;i<n;i++)if(i<c->reverse[i]) {
        uint32_t t=a0[i];a0[i]=a0[c->reverse[i]];a0[c->reverse[i]]=t;
        t=a1[i];a1[i]=a1[c->reverse[i]];a1[c->reverse[i]]=t;
    }
    const uint32_t *roots=c->stages;
#if FALCON_REC_SHOUP_NTT
    const uint32_t *factors=c->shoup;
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
#pragma GCC unroll 2
      for(unsigned pass=0;pass<2;pass++) {
        uint32_t *a=pass?a1:a0;
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
    }
#undef FALCON_REC_SUMDIFF
    first_length=16;
#endif
    for(unsigned len=first_length;len<=n;len*=2)for(unsigned start=0;start<n;start+=len) {
        unsigned k=0;
#if defined(__AVX2__) && !defined(FALCON_REC_NO_SIMD)
        __m256i p=_mm256_set1_epi32((int)c->p);
        for(;k+8<=len/2;k+=8) {
            __m256i root=_mm256_loadu_si256((const __m256i *)(roots+len/2+k));
            __m256i factor=_mm256_loadu_si256((const __m256i *)(factors+len/2+k));
#pragma GCC unroll 2
          for(unsigned pass=0;pass<2;pass++) {
            uint32_t *a=pass?a1:a0;
            __m256i u=_mm256_loadu_si256((const __m256i *)(a+start+k));
#if FALCON_REC_SHOUP_NTT
            __m256i v=shoup8(c,_mm256_loadu_si256((const __m256i *)(a+start+k+len/2)),root,factor);
#else
            __m256i v=mont8(c,_mm256_loadu_si256((const __m256i *)(a+start+k+len/2)),_mm256_loadu_si256((const __m256i *)(roots+len/2+k)));
#endif
            __m256i sum=_mm256_add_epi32(u,v),dif=_mm256_sub_epi32(u,v);
            sum=_mm256_sub_epi32(sum,_mm256_andnot_si256(lt8(sum,p),p));
            dif=_mm256_add_epi32(dif,_mm256_and_si256(lt8(u,v),p));
            _mm256_storeu_si256((__m256i *)(a+start+k),sum);_mm256_storeu_si256((__m256i *)(a+start+k+len/2),dif);
          }
        }
#endif
        for(;k<len/2;k++) {
          for(unsigned pass=0;pass<2;pass++) {
            uint32_t *a=pass?a1:a0;
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

}
#endif
