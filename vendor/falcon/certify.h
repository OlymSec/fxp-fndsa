#ifndef STUDY_CERTIFY_H
#define STUDY_CERTIFY_H
/* Preparation-only checks. All tests are exact integer inequalities with
 * outward Q16 bounds; no sampling-centre error assumption is made. */
static uint64_t cert_abs_q16(fpr x,uint64_t cap)
{
    uint128_t a=(x.v>>127)?-x.v:x.v;
    uint128_t v=(a+(((uint128_t)1<<48)-1))>>48;
    return v>cap?cap+1:(uint64_t)v;
}
static int cert_norm(const fpr *a,unsigned count,unsigned limit,unsigned margin)
{
    uint64_t cap=(uint64_t)limit<<16,sum=0;
    for(unsigned i=0;i<count;i++) {
        uint64_t x=cert_abs_q16(a[i],cap)+margin;
        if(x>cap)return 0;
        sum+=x*x;
    }
    return sum<=cap*cap;
}
static int cert_children(const fpr *tree,unsigned logn,unsigned orig_logn,fpr *mirrors)
{
    if(!logn) {
        /* Also keep the existing scalar sampler/mirrored-leaf arithmetic
         * within a plainly representable range. */
        fpr lo=fpr_div2e(fpr_of(1),2),hi=fpr_of(2);
        if((int128_t)tree[0].v<(int128_t)lo.v || tree[0].v>hi.v)return 0;
        fpr mirror=compute_mirror_leaf(tree[0],orig_logn);
        if(mirrors)*mirrors=mirror;
        return (int128_t)mirror.v>=(int128_t)lo.v && mirror.v<=hi.v;
    }
    size_t n=(size_t)1<<logn;
    for(size_t j=0;j<n/2;j++) {
        uint64_t re=cert_abs_q16(tree[j],UINT64_C(1)<<17);
        uint64_t im=cert_abs_q16(tree[j+n/2],UINT64_C(1)<<17);
        /* |L| <= 11/10, with outward rounding. */
        if(100*(re*re+im*im)>(UINT64_C(121)<<32))return 0;
    }
    return cert_children(tree+n,logn-1,orig_logn,mirrors)
        && cert_children(tree+n+((size_t)logn<<(logn-1)),logn-1,orig_logn,mirrors ? mirrors+n/2 : NULL);
}
static void cert_cmul(fpr *re,fpr *im,fpr ar,fpr ai,fpr br,fpr bi)
{
    *re=fpr_sub(fpr_mul(ar,br),fpr_mul(ai,bi));
    *im=fpr_add(fpr_mul(ar,bi),fpr_mul(ai,br));
}
static int cert_basis_tree(const fpr *b,unsigned logn,fpr *mirrors)
{
    size_t n=(size_t)1<<logn,hn=n/2;const fpr *L=b+4*n;
    if(!cert_children(L+n,logn-1,logn,mirrors))return 0;
    for(size_t j=0;j<hn;j++) {
        fpr l[2]={L[j],L[j+hn]};
        if(!cert_norm(l,2,64,0))return 0;
        fpr first[4]={b[j],b[j+hn],b[n+j],b[n+j+hn]};
        /* Four Q16 units cover the proven FFT table/rounding error. */
        if(!cert_norm(first,4,512,4))return 0;
        fpr orth[4];
        for(unsigned k=0;k<2;k++) {
            fpr br=b[(2+k)*n+j],bi=b[(2+k)*n+j+hn],r,i;
            uint64_t raw_cap=UINT64_C(262144)<<16;
            if(cert_abs_q16(br,raw_cap)>raw_cap || cert_abs_q16(bi,raw_cap)>raw_cap)return 0;
            cert_cmul(&r,&i,l[0],l[1],first[2*k],first[2*k+1]);
            orth[2*k]=fpr_sub(br,r);orth[2*k+1]=fpr_sub(bi,i);
        }
        if(!cert_norm(orth,4,512,4))return 0;
    }
    return 1;
}
#endif
