#include "inner.h"
#include <stdio.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #c); exit(1); } } while (0)
void falcon_support_reference_i_shake256_init(inner_shake256_context *);
void falcon_support_reference_i_shake256_inject(inner_shake256_context *,const uint8_t *,size_t);
void falcon_support_reference_i_shake256_flip(inner_shake256_context *);
void falcon_support_reference_i_shake256_extract(inner_shake256_context *,uint8_t *,size_t);
#if defined(FALCON_SUPPORT_WRAP_AUXV)
#include <sys/auxv.h>
#include <asm/hwcap.h>
static int force_portable;
unsigned long __real_getauxval(unsigned long);
unsigned long __wrap_getauxval(unsigned long type) {
    unsigned long value=__real_getauxval(type);
    return force_portable && type==AT_HWCAP ? value & ~(unsigned long)HWCAP_SHA3 : value;
}
#endif
static uint64_t random_state=UINT64_C(0x202609071a4c923b);
static uint32_t word(void) {
    random_state^=random_state>>12; random_state^=random_state<<25; random_state^=random_state>>27;
    return (uint32_t)((random_state*UINT64_C(0x2545f4914f6cdd1d))>>32);
}
static void state_equal(const inner_shake256_context *a,const inner_shake256_context *b) {
    CHECK(a->dptr==b->dptr && !memcmp(a->st.A,b->st.A,sizeof a->st.A));
}
static void shake_tests(void) {
    static const char kat[]="46b9dd2b0ba88d13233b3feb743eeb243fcd52ea62b81b82b50c27646ed5762fd75dc4ddd8c0f200cb05019d67b592f6fc821c49479ab48640292eacb3b7c4be";
    uint8_t input[4096],a[4096],b[4096]; inner_shake256_context x,y,copy;
    inner_shake256_init(&x); inner_shake256_flip(&x); inner_shake256_extract(&x,a,64);
    for(unsigned i=0;i<64;i++) { unsigned expected; CHECK(sscanf(kat+2*i,"%2x",&expected)==1); CHECK(a[i]==expected); }
    for(unsigned trial=0;trial<2000;trial++) {
        size_t len=trial<276?trial:word()%sizeof input;
        for(size_t i=0;i<len;i++) input[i]=(uint8_t)word();
        inner_shake256_init(&x); falcon_support_reference_i_shake256_init(&y);
        for(size_t pos=0;pos<len;) {
            size_t take=1+word()%173; if(take>len-pos) take=len-pos;
            inner_shake256_inject(&x,input+pos,take);
            falcon_support_reference_i_shake256_inject(&y,input+pos,take);
            pos+=take; state_equal(&x,&y);
        }
        inner_shake256_flip(&x); falcon_support_reference_i_shake256_flip(&y); state_equal(&x,&y);
        copy=x;
        for(size_t pos=0;pos<sizeof a;) {
            size_t take=1+word()%173; if(take>sizeof a-pos) take=sizeof a-pos;
            inner_shake256_extract(&x,a+pos,take);
            falcon_support_reference_i_shake256_extract(&y,b+pos,take);
            CHECK(!memcmp(a+pos,b+pos,take)); state_equal(&x,&y); pos+=take;
        }
        inner_shake256_extract(&copy,b,sizeof b); CHECK(!memcmp(a,b,sizeof a)); state_equal(&x,&copy);
    }
    for(unsigned trial=0;trial<5000;trial++) {
        for(unsigned i=0;i<25;i++) x.st.A[i]=((uint64_t)word()<<32)|word();
        x.dptr=136; y=x;
        inner_shake256_extract(&x,a,1); falcon_support_reference_i_shake256_extract(&y,b,1);
        CHECK(a[0]==b[0]); state_equal(&x,&y);
    }
    puts("SHAKE: standard KAT, 2000 fragmented streams/copies, 5000 full random states PASS");
}
static void reference_hash(inner_shake256_context *sc,uint16_t *out,unsigned logn) {
    for(size_t n=(size_t)1<<logn;n;) {
        uint8_t bytes[2]; falcon_support_reference_i_shake256_extract(sc,bytes,2);
        unsigned value=((unsigned)bytes[0]<<8)|bytes[1];
        if(value<61445) { while(value>=12289) value-=12289; *out++=(uint16_t)value; n--; }
    }
}
static void hash_tests(void) {
    inner_shake256_context a,b; uint16_t x[1024],y[1024]; uint8_t seed[65],tail_a[200],tail_b[200];
    for(unsigned trial=0;trial<3000;trial++) {
        unsigned logn=trial%11;
        for(unsigned i=0;i<sizeof seed;i++) seed[i]=(uint8_t)word();
        inner_shake256_init(&a); inner_shake256_inject(&a,seed,sizeof seed); inner_shake256_flip(&a);
        if(trial<100) {
            memset(a.st.dbuf,255,136); a.dptr=0;
            for(unsigned i=0;i<68;i++) if(i%7==0) {
                unsigned v=(unsigned[]){0,12288,12289,61444,61445,65535}[i%6];
                a.st.dbuf[2*i]=(uint8_t)(v>>8); a.st.dbuf[2*i+1]=(uint8_t)v;
            }
        }
        b=a; Zf(hash_to_point_vartime)(&a,x,logn); reference_hash(&b,y,logn);
        CHECK(!memcmp(x,y,((size_t)1<<logn)*sizeof *x)); state_equal(&a,&b);
        inner_shake256_extract(&a,tail_a,sizeof tail_a);
        falcon_support_reference_i_shake256_extract(&b,tail_b,sizeof tail_b);
        CHECK(!memcmp(tail_a,tail_b,sizeof tail_a));
    }
    puts("Hash-to-point: 3000 cases, rejection boundaries, identical full state and suffix PASS");
}
static unsigned bit(const uint8_t *in,size_t pos) { return (in[pos/8]>>(7-pos%8))&1; }
static size_t reference_decode(int16_t *out,unsigned logn,const uint8_t *in,size_t len) {
    size_t pos=0;
    for(size_t i=0;i<((size_t)1<<logn);i++) {
        if(pos/8>=len || len-pos/8<1+(pos%8!=0)) return 0;
        unsigned negative=bit(in,pos++),m=0;
        for(unsigned j=0;j<7;j++) m=(m<<1)|bit(in,pos++);
        for(;;) {
            if(pos/8>=len) return 0;
            if(bit(in,pos++)) break;
            m+=128; if(m>2047) return 0;
        }
        if(negative && !m) return 0;
        out[i]=(int16_t)(negative?-(int)m:(int)m);
    }
    size_t consumed=(pos+7)/8;
    while(pos<8*consumed) if(bit(in,pos++)) return 0;
    return consumed;
}
static void check_decode(unsigned logn,const uint8_t *bytes,size_t len) {
    int16_t a[1024],b[1024];
    size_t got=Zf(comp_decode)(a,logn,bytes,len),expected=reference_decode(b,logn,bytes,len);
    CHECK(got==expected); if(got) CHECK(!memcmp(a,b,((size_t)1<<logn)*sizeof *a));
}
static void codec_tests(void) {
    uint8_t bytes[4096]; int16_t coefficients[1024],decoded[1024];
    for(unsigned logn=0;logn<=10;logn++) for(unsigned trial=0;trial<2000;trial++) {
        size_t n=(size_t)1<<logn;
        for(size_t i=0;i<n;i++) coefficients[i]=(int16_t)((int)(word()%4095)-2047);
        if(trial<16) for(size_t i=0;i<n;i++) coefficients[i]=(int16_t)((int[]){0,1,-1,127,-127,128,-128,255,-255,256,-256,1023,-1023,2047,-2047,0}[trial]);
        size_t len=Zf(comp_encode)(bytes,sizeof bytes,coefficients,logn); CHECK(len);
        CHECK(Zf(comp_decode)(decoded,logn,bytes,len)==len && !memcmp(decoded,coefficients,n*sizeof *decoded));
        check_decode(logn,bytes,len);
        if(trial<16) for(size_t limit=0;limit<len;limit++) check_decode(logn,bytes,limit);
        else check_decode(logn,bytes,word()%(len+1));
        for(unsigned mutation=0;mutation<4;mutation++) {
            size_t index=word()%len; uint8_t mask=(uint8_t)(1U<<(word()%8));
            bytes[index]^=mask; check_decode(logn,bytes,len); bytes[index]^=mask;
        }
        for(size_t i=0;i<len+1;i++) bytes[i]=(uint8_t)word();
        check_decode(logn,bytes,len+1);
    }
    memset(bytes,0,sizeof bytes); bytes[0]=0x80; bytes[1]=0x80; CHECK(!Zf(comp_decode)(decoded,0,bytes,2));
    puts("Decoder: 22000 valid vectors, truncations, mutations, arbitrary inputs and negative zero PASS");
}
int main(int argc,char **argv) {
#if defined(FALCON_SUPPORT_WRAP_AUXV)
    force_portable=argc==2 && !strcmp(argv[1],"--portable");
    CHECK(argc==1 || force_portable);
    printf("ARM SHA3 supported=%u; forced portable=%d\n",(unsigned)((__real_getauxval(AT_HWCAP)&HWCAP_SHA3)!=0),force_portable);
#else
    (void)argc; (void)argv;
#endif
    shake_tests(); hash_tests(); codec_tests(); return 0;
}
