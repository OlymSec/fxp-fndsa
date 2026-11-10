#!/usr/bin/env python3
"""Export only the SDK public interface, including in the static archive."""
import re
import sys
from pathlib import Path


def names(header):
    return sorted(set(re.findall(r'\b(fxp_fndsa_\w+)\s*\(', Path(header).read_text())))


if __name__ == '__main__':
    mode, header = sys.argv[1:]
    api = names(header)
    if mode == 'names':
        print(' '.join(api))
    elif mode == 'map':
        print('{ global:\n'+'\n'.join('    '+name+';' for name in api)+'\nlocal: *; };')
    elif mode == 'darwin':
        print('\n'.join('_'+name for name in api))
    else:
        raise SystemExit('unknown export format')
