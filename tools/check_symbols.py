#!/usr/bin/env python3
"""Check static and shared exports and absence of OpenSSL/GMP dependencies."""
import subprocess
import sys
from exports import names

header, static, shared = sys.argv[1:]
expected = set(names(header))
for library, dynamic in [(static, False), (shared, True)]:
    cmd = ['nm', '-gU', library] if sys.platform == 'darwin' else ['nm', *(['-D'] if dynamic else []), '--defined-only', '-g', library]
    output = subprocess.check_output(cmd, text=True)
    actual = {line.split()[-1].lstrip('_') for line in output.splitlines() if len(line.split()) == 3}
    if actual != expected:
        raise SystemExit(f'{library}: missing={expected-actual}, extra={actual-expected}')
    undefined = subprocess.check_output(['nm', '-u', library], text=True)
    if any(s in undefined for s in ['OPENSSL_', 'EVP_', '__gmp']):
        raise SystemExit(f'{library}: unexpected arithmetic or crypto-library dependency')
print(f'PASS static and shared exports: exactly {len(expected)} fxp_fndsa functions')
