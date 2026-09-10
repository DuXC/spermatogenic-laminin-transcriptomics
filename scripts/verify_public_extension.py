"""Verify frozen v0.3 numerical data before running reproduction."""
from pathlib import Path
import hashlib, json
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'public_extension_manifest.json').read_text())
for row in m['files']:
    p=root/row['path']
    assert p.is_file(), f"Missing file: {row['path']}"
    assert p.stat().st_size==row['bytes'], f"Size mismatch: {row['path']}"
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'], f"Checksum mismatch: {row['path']}"
print(f"PASS: {len(m['files'])} v0.3 numerical source files match the frozen manifest")
