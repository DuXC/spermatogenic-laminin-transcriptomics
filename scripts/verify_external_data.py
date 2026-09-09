from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'external_data_manifest.json').read_text())
for f in manifest['files']:
 p=root/f['path'];assert p.exists(),f['path'];assert p.stat().st_size==f['bytes'],f['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256'],f['path']
print('PASS',len(manifest['files']),'external release files')
