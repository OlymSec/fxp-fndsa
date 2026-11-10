#!/usr/bin/env python3
"""Verify the unmodified import of the measured 512/1024 implementation."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root/'provenance/import-sha256.json').read_text())
actual = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in (root/'vendor/falcon').rglob('*') if p.is_file()}
if actual != manifest:
    changed = sorted(k for k in actual.keys() | manifest.keys() if actual.get(k) != manifest.get(k))
    raise SystemExit('Imported source differs: '+', '.join(changed))
print(f'PASS {len(actual)} imported files match the port snapshot')

xkcp=root/'vendor/xkcp'
manifest=json.loads((xkcp/'provenance.json').read_text())
for name,digest in manifest['sha256'].items():
    if hashlib.sha256((xkcp/name).read_bytes()).hexdigest()!=digest:
        raise SystemExit('XKCP upstream source differs: '+name)
print('PASS pinned XKCP single-state upstream sources')
