from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import gzip, hashlib, json, urllib.request

S=Path(__file__).resolve().parents[1]
d=json.loads((S/'01_metadata/GSE153947_metadata.json').read_text())
jobs=[]
for x in d['samples']:
    for k,vs in x.items():
        if k.startswith('supplementary_file'):
            for url in vs:
                if url.endswith('_counts.tsv.gz'): jobs.append((url.replace('ftp://','https://'),S/'02_raw'/url.rsplit('/',1)[1]))
for url in d['series']['supplementary_file']:
    if 'Cell_metadata' in url:jobs.append((url.replace('ftp://','https://'),S/'01_metadata'/url.rsplit('/',1)[1]))
assert len(jobs)==7
def run(job):
    url,path=job
    try:
        if not path.exists():
            req=urllib.request.Request(url,headers={'User-Agent':'Scientific-public-data-reanalysis/1.0'})
            with urllib.request.urlopen(req,timeout=60) as r:
                assert r.status==200
                payload=r.read()
            assert payload[:2]==b'\x1f\x8b'
            path.write_bytes(payload)
        with gzip.open(path,'rb') as f:
            while f.read(1024*1024):pass
        return {'status':'VERIFIED_GZIP','url':url,'path':str(path.relative_to(S)),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    except Exception as e:return {'status':'ERROR','url':url,'error':str(e)}
rows=[]
with ThreadPoolExecutor(max_workers=3) as pool:
    for r in pool.map(run,jobs):rows.append(r);print(json.dumps(r),flush=True)
(S/'00_admin/sc_download_receipt.json').write_text(json.dumps({'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'files':rows},indent=2)+'\n')
assert all(r['status']=='VERIFIED_GZIP' for r in rows)
