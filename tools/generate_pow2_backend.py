#!/usr/bin/env python3
'''Asserted adapters for preallocated power-of-two preparation and Q80 conversion.'''
import argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def replace(s,old,new,count=1):
    if s.count(old)!=count: raise RuntimeError(f'Expected {count} occurrences of {old!r}')
    return s.replace(old,new)
def capture():
    s=(ROOT/'vendor/falcon/capture.h').read_text()
    s=replace(s,'#include "certify.h"','#include "certify.h"\n#include "pow2_clean.h"\n#include <stdlib.h>\nint study_rec_prepare_work(study_rec_key *,const int32_t *,unsigned,unsigned,uint32_t *);')
    start=s.index('int study_prepare_expanded(');end=s.index('\ntypedef struct {\n    samplerZ',start)
    block=s[start:end]
    # The common finish preserves the existing certificate and NTRU check.
    tail=block[block.index('    fpr *mirrors=NULL;'):]
    tail=replace(tail,'study_rec_prepare(&prepared->key,key,logn)','study_rec_prepare_work(&prepared->key,key,n,logn,ntt_work)')
    finish='''static int study_prepare_finish(study_prepared *prepared,const fpr *candidate,
    const int32_t *key,unsigned logn,uint32_t *ntt_work)
{
    unsigned n=1U<<logn;
'''+tail
    general='''int study_prepare_expanded_work(void *memory,const void *encoded_expanded,void *work)
{
    study_prepared *prepared=memory;
    if(!prepared)return -1;
    memset(prepared,0,sizeof *prepared);
    if(!encoded_expanded)return -1;
    unsigned logn=*(const uint8_t *)encoded_expanded;
    if(logn!=9 && logn!=10)return -1;
    unsigned n=1U<<logn;
    uintptr_t ptr=((uintptr_t)encoded_expanded+sizeof(fpr))&~(uintptr_t)(sizeof(fpr)-1);
    const fpr *candidate=(const fpr *)ptr;
    fpr *a=work;
    int32_t *key=(int32_t *)(a+n);
    uint32_t *ntt_work=(uint32_t *)(key+4*n);
    unsigned slotmap[4]={0,2,1,3};
    int status=0;
    for(unsigned slot=0;slot<4;slot++) {
        memcpy(a,candidate+slotmap[slot]*n,n*sizeof *a);
        Zf(iFFT)(a,logn);
        for(unsigned j=0;j<n;j++) {
            int64_t x=fpr_rint(a[j]);if(slot>=2)x=-x;
            if(x< -127 || x>127){status=-2;goto done;}
            key[slot*n+j]=(int32_t)x;
        }
        for(unsigned j=0;j<n;j++)a[j]=fpr_of(slot>=2?-(int64_t)key[slot*n+j]:key[slot*n+j]);
        Zf(FFT)(a,logn);
        if(memcmp(a,candidate+slotmap[slot]*n,n*sizeof *a)){status=-3;goto done;}
    }
    status=study_prepare_finish(prepared,candidate,key,logn,ntt_work);
done:
    return status;
}
int study_prepare_expanded(void *memory,const void *encoded_expanded)
{
    if(!memory || !encoded_expanded)return -1;
    unsigned logn=*(const uint8_t *)encoded_expanded;
    if(logn!=9 && logn!=10)return -1;
    size_t bytes=((size_t)1<<logn)*(sizeof(fpr)+5*sizeof(uint32_t));
    void *work=malloc(bytes);
    if(!work)return -6;
    int status=study_prepare_expanded_work(memory,encoded_expanded,work);
    xp_pow2_clean(work,bytes);
    free(work);
    return status;
}
/* Internal trusted handoff: expanded was just produced by expand_privkey
 * from this exact bounded basis. The general importer above retains all
 * representation checks. Certificate and NTRU checks are common to both. */
int study_prepare_basis_work(void *memory,const void *encoded_expanded,
    const int8_t *basis,void *work)
{
    study_prepared *prepared=memory;
    unsigned logn=*(const uint8_t *)encoded_expanded;
    if(logn!=9 && logn!=10)return -1;
    unsigned n=1U<<logn;
    int32_t *key=work;
    const unsigned slots[4]={1,3,0,2};
    for(unsigned k=0;k<4;k++)for(unsigned j=0;j<n;j++) {
        int32_t x=basis[slots[k]*n+j];
        if(x< -127 || x>127)return -2;
        key[k*n+j]=x;
    }
    memset(prepared,0,sizeof *prepared);
    uintptr_t ptr=((uintptr_t)encoded_expanded+sizeof(fpr))&~(uintptr_t)(sizeof(fpr)-1);
    return study_prepare_finish(prepared,(const fpr *)ptr,key,logn,(uint32_t *)(key+4*n));
}
'''
    s=s[:start]+finish+general+s[end:]
    trace='''#if defined(XP_POW2_TRACE) && XP_POW2_TRACE
static _Thread_local uint64_t study_trace_value;
static _Thread_local unsigned study_trace_calls;
void study_trace_reset(void) { study_trace_value=UINT64_C(1469598103934665603); study_trace_calls=0; }
uint64_t study_trace_get(void) { return study_trace_value; }
unsigned study_trace_count(void) { return study_trace_calls; }
static void study_trace_bytes(const void *memory,size_t len)
{
    const uint8_t *p=memory;
    for(size_t i=0;i<len;i++)study_trace_value=(study_trace_value^p[i])*UINT64_C(1099511628211);
}
#endif
'''
    s=replace(s,'static int study_capture_sample(',trace+'static int study_capture_sample(')
    s=replace(s,'    int z=c->sampler(c->ctx,mu,sigma);','''    int z=c->sampler(c->ctx,mu,sigma);
#if defined(XP_POW2_TRACE) && XP_POW2_TRACE
    const sampler_context *spc=c->ctx;
    study_trace_bytes(&mu,sizeof mu);study_trace_bytes(&sigma,sizeof sigma);
    study_trace_bytes(&z,sizeof z);study_trace_bytes(&spc->p.ptr,sizeof spc->p.ptr);
    study_trace_bytes(spc->p.state.d,sizeof spc->p.state.d);
    study_trace_bytes(spc->p.buf.d,sizeof spc->p.buf.d);
    study_trace_calls++;
#endif''')
    return s
def rec():
    s=(ROOT/'vendor/falcon/rec_ntt.c').read_text()
    s=replace(s,'#include <string.h>','#include <string.h>\n#include <stdlib.h>\n#include "pow2_clean.h"')
    s=replace(s,'int study_rec_prepare(study_rec_key *c,const int32_t key[4][1024],unsigned logn)',
              'int study_rec_prepare_work(study_rec_key *c,const int32_t *key,unsigned stride,unsigned logn,uint32_t *v)',2)
    s=replace(s,'rec_reduce(key[k][j])','rec_reduce(key[k*stride+j])')
    s=replace(s,'int32_t x=key[slot][j];','int32_t x=key[slot*stride+j];')
    s=replace(s,'uint32_t v[1024];','/* n-word validation workspace supplied by caller */',2)
    s+='''
int study_rec_prepare(study_rec_key *c,const int32_t key[4][1024],unsigned logn)
{
    if(logn!=9 && logn!=10)return 0;
    size_t bytes=((size_t)1<<logn)*sizeof(uint32_t);
    uint32_t *v=malloc(bytes);
    if(!v)return 0;
    int status=study_rec_prepare_work(c,&key[0][0],1024,logn,v);
    xp_pow2_clean(v,bytes);
    free(v);
    return status;
}
'''
    return s
def half_fft():
    s=(ROOT/'vendor/falcon/ntrugen/src/half/fft.h').read_text()
    s=replace(s,'#include "wide.h"','#include "half/wide.h"')
    s=replace(s,'#include "roots.h"','#include "half/roots.h"\n#include "pow2_window.h"\n#ifndef XP_POW2_WINDOW\n#define XP_POW2_WINDOW 1\n#endif')
    start=s.index('        for (unsigned j = 0; j < len; j++) {',s.index('static void half_convert'))
    end=s.index('        hw lo =',start)
    old=s[start:end]
    new='''#if XP_POW2_WINDOW
        if (len > 4) { /* public input length */
            half_window(w,f,n,i,len,sch,sign);
        } else
#endif
        {
'''+old+'''        }
'''
    return s[:start]+new+s[end:]
def half_solver():
    s=(ROOT/'vendor/falcon/ntrugen/src/half/solver.c').read_text()
    return replace(s,'#include "fft.h"','#include "half_fft_backend.h"')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['capture','rec','half-fft','half-solver']);p.add_argument('--output',type=Path);p.add_argument('--check',action='store_true');a=p.parse_args()
    if a.check:capture();rec();half_fft();half_solver();print('PASS uniquely matched preparation, reconstruction and Q80-window adapters')
    if a.output:a.output.write_text({'capture':capture,'rec':rec,'half-fft':half_fft,'half-solver':half_solver}[a.kind]())
