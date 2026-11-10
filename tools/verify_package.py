#!/usr/bin/env python3
"""Verify all source-package files against SHA256SUMS."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[1]
count=0
for line in (root/'SHA256SUMS').read_text().splitlines():
    digest,name=line.split('  ',1)
    path=Path(name)
    assert not path.is_absolute() and '..' not in path.parts
    assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest,name
    count+=1
print(f'PASS: {count} source-package files')
