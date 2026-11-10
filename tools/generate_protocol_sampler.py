#!/usr/bin/env python3
"""Expose one Q64 tree attempt with a continuing sequential SHAKE stream."""
from pathlib import Path
import argparse
ROOT=Path(__file__).resolve().parents[1]

def generate():
    source=(ROOT/'vendor/falcon/sign.c').read_text()
    source='#ifndef FALCON_IMPL_SKIP_DEAD_OUTPUT\n#define FALCON_IMPL_SKIP_DEAD_OUTPUT 1\n#endif\n'+source
    assert source.count('#include "capture.h"') == 1
    source=source.replace('#include "capture.h"','#include "capture_backend.h"')
    wipe = "\tvolatile uint8_t *clear = (volatile uint8_t *)&spc;\n\tfor (size_t u = 0; u < sizeof spc; u ++) clear[u] = 0;"
    assert source.count(wipe) == 1
    source=source.replace(wipe, "\txp_pow2_clean(&spc, sizeof spc);")
    start=source.index('int\nZf(mrm_sign_tree_attempt)(')
    end=source.index('\n#endif',start)
    function=source[start:end]
    assert function.count('Zf(prng_init)(&spc.p, rng);')==1
    function=function.replace('Zf(mrm_sign_tree_attempt)','xp_tree_attempt',1)
    function=function.replace('Zf(prng_init)(&spc.p, rng);', '''/* Resume the attempt stream immediately after its nonce. */
    _Static_assert(sizeof *rng <= sizeof spc.p.state.d, "SHAKE state capacity");
    memset(spc.p.state.d, 0, sizeof spc.p.state.d);
    memcpy(spc.p.state.d, rng, sizeof *rng);
    spc.p.type = 1;
    Zf(prng_refill)(&spc.p);''')
    return source+'\n/* FN-DSA protocol adapter; pinned sampler arithmetic above is unchanged. */\n'+function+'\n'

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--check',action='store_true');a=p.parse_args()
    s=generate()
    if a.output:a.output.write_text(s)
    if a.check:print('PASS uniquely matched single-attempt sampler adapter')
