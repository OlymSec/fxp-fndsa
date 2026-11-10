#!/usr/bin/env python3
"""Directed root intervals and a rational bound for trusted integer-basis handoff."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,re
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args()
source=ROOT/'vendor/falcon/fpr.c'
s=source.read_text().rsplit('const fpr fpr_gm_tab[] = {',1)[1].split('};',1)[0]
words=[int(x,16)*2**64+int(y,16) for x,y in re.findall(r'UINT128\((0x[0-9a-f]+),\s*(0x[0-9a-f]+)\)',s)]
words=[x-2**128 if x>=2**127 else x for x in words]
assert len(words)==2048
mp.iv.dps=90
for k in range(1,1024):
    rev=int(f'{k:010b}'[::-1],2);angle=mp.iv.pi*rev/1024
    for j,fn in enumerate([mp.iv.cos,mp.iv.sin]):
        x=words[2*k+j];e=fn(angle)*2**64-x
        assert bool(e>-2) and bool(e<2),(k,j,e)
        assert abs(x)<=2**64
u=F(1,2**64);rows=[]
for logn in [9,10]:
    # Forward/inverse each have logn-1 nontrivial butterfly stages.
    # Represented magnitude <= 127*4^j + (4*u/3)*(4^j-1) < M_j=128*4^j.
    # Component error: E_{j+1} <= 4 E_j + (8 M_j+8)*u.
    steps=2*(logn-1);m=F(128);e=F(0);represented=F(127)
    for j in range(steps):
        e=4*e+(8*m+8)*u;m*=4
        represented=4*represented+4*u
        assert represented < m
        assert e <= (j+1)*4**(j+1)*258*u
    final=e/F(2**(logn-1))+u/2
    assert final<F(1,2**24)
    assert final<F(1,2)
    assert 4*m < 2**63 # sums and real factors, before Q64 scaling
    rows.append({'degree':2**logn,'stages':steps,'roundtrip_absolute_error_upper':str(final),
                 'roundtrip_error_less_than': '2^-24',
                 'integer_recovery_unique':True,'largest_component_carrier_bound':str(m)})
result={'status':'PASS','scope':'Q64 table error and integer FFT roundtrip, coefficients in [-127,127]; preserves old import recovery test for internally expanded keys','interval_decimal_precision':90,'checked_nonzero_root_components':2046,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':rows}
if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
