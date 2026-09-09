"""Verify the released data bytes and the complete-member scientific structure."""
from pathlib import Path
import csv,hashlib,json
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'data_manifest.json').read_text())
bad=[]
for x in m['files']:
    p=r/x['file']
    if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=x['sha256']:bad.append(x['file'])
assert not bad, bad
t=r/'08_manuscript_refinement_20260909/03_tables'
with (t/'Table_S3_all_30_laminin_gene_evidence.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==30 and len({x['gene'] for x in rows})==30
print('PASS',len(m['files']),'data checksums and 30 complete laminin members')
