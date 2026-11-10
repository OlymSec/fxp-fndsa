#!/usr/bin/env python3
"""Generate the NTT-key adapters with asserted matches against the pinned import."""
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace(source, old, new, count=1):
    if source.count(old) != count:
        raise RuntimeError(f'Public-key adapter expected {count} matches: {old!r}')
    return source.replace(old, new)


def generate(kind, ntt):
    source = (ROOT / 'vendor/falcon' / (kind + '.c')).read_text()
    if not ntt:
        return source
    if kind == 'falcon':
        source = replace(source, '\tZf(to_ntt_monty)(h, logn);\n',
                         '\t/* The public-key payload is already NTT/Montgomery. */\n')
        source = replace(source, 'pk[0] = 0x00 + logn;',
                         'Zf(to_ntt_monty)(h, logn);\n\tpk[0] = 0x10 + logn;', 2)
        source = replace(source, 'if ((pk[0] & 0xF0) != 0x00)',
                         'if ((pk[0] & 0xF0) != 0x10)')
    elif kind == 'falcon_mrm':
        source = replace(source, '''    if (!pk || !length || (pk[0] != 9 && pk[0] != 10)) return -1;
    if (length != FALCON_PUBKEY_SIZE(pk[0])) return -1;
    return pk[0];''', '''    if (!pk || !length || (pk[0] != 0x19 && pk[0] != 0x1a)) return -1;
    unsigned logn = pk[0] & 15U;
    if (length != FALCON_PUBKEY_SIZE(logn)) return -1;
    return (int)logn;''')
        source = replace(source, '    Zf(to_ntt_monty)(h, logn);\n',
                         '    /* The public-key payload is already NTT/Montgomery. */\n')
    else:
        raise ValueError(kind)
    return '/* Generated NTT public-key adapter; do not edit. */\n' + source


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kind', choices=['falcon', 'falcon_mrm'])
    p.add_argument('--ntt', type=int, choices=[0, 1], default=1)
    p.add_argument('--output', type=Path)
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    if a.check:
        for kind in ['falcon', 'falcon_mrm']:
            generate(kind, True)
            assert generate(kind, False) == (ROOT / 'vendor/falcon' / (kind + '.c')).read_text()
        print('PASS uniquely matched NTT-key adapters and unchanged coefficient backend')
    if a.output:
        if not a.kind:
            p.error('--output requires --kind')
        a.output.write_text(generate(a.kind, a.ntt))
    if not a.check and not a.output:
        p.error('select --output or --check')


if __name__ == '__main__':
    main()
