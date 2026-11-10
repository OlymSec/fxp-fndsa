#!/usr/bin/env python3
"""Generate the narrow XKCP hook without editing the imported Falcon code."""
from pathlib import Path
import argparse
import re

ROOT=Path(__file__).resolve().parents[1]


def adapter():
    source=(ROOT/'vendor/xkcp/KeccakP-1600-opt64.c').read_text()
    constants=re.findall(r'static const uint64_t KeccakF1600RoundConstants\[24\] = \{.*?\};',source,re.S)
    if len(constants)!=1:
        raise RuntimeError('XKCP round-constant definition changed')
    header=source[:source.index('---')].rstrip()+'\n*/\n'
    return header+'''#include <stdint.h>
#define ROL64(a, offset) ((((uint64_t)(a)) << (offset)) ^ (((uint64_t)(a)) >> (64-(offset))))
#define FullUnrolling
#include "KeccakP-1600-64.macros"
#include "KeccakP-1600-unrolling.macros"
'''+constants[0]+'''

/* ANSA adapter: ordinary-layout, 24-round permutation only. The upstream
 * copyFromState/rounds24/copyToState body uses Falcon's state pointer. */
void fxp_keccak_permute(uint64_t stateAsLanes[25])
{
    declareABCDE
    copyFromState(A, stateAsLanes)
    rounds24
    copyToState(stateAsLanes, A)
}
'''


def shake():
    source=(ROOT/'vendor/falcon/shake.c').read_text()
    include='#include "inner.h"'
    definition='static void\nprocess_block(uint64_t *A)\n{\n'
    if source.count(include)!=1 or source.count(definition)!=1:
        raise RuntimeError('Falcon SHAKE hook no longer matches exactly once')
    source=source.replace(include,include+'\n#if FALCON_XKCP_SHAKE\nvoid fxp_keccak_permute(uint64_t state[25]);\n#endif',1)
    return source.replace(definition,definition+'''#if FALCON_XKCP_SHAKE
    fxp_keccak_permute(A);
    return;
#endif
''',1)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path);p.add_argument('--write-adapter',action='store_true');p.add_argument('--check',action='store_true');a=p.parse_args()
    if a.write_adapter:(ROOT/'vendor/xkcp/keccak_permute.c').write_text(adapter())
    if a.output:a.output.write_text(shake())
    if a.check:
        if (ROOT/'vendor/xkcp/keccak_permute.c').read_text()!=adapter():raise SystemExit('stale XKCP permutation adapter')
        shake()
        print('PASS exact XKCP adapter and uniquely matched Falcon permutation hook')
    if not (a.write_adapter or a.output or a.check):p.error('select --output, --write-adapter or --check')
