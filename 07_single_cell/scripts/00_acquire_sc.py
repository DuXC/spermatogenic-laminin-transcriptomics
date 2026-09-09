"""Acquire individual public count matrices with byte-level provenance."""
from pathlib import Path
import csv,hashlib,json,urllib.request,gzip,concurrent.futures,datetime
SC=Path(__file__).resolve().parents[1];ROOT=SC.parent
meta=list(csv.DictReader((ROOT/'02_annotation/GSE149512_donor_metadata_audit.tsv').open(),delimiter='\t'))
selected=[r for r in meta if r['paper_group'] in ['adult_normal_spermatogenesis_OA','iNOA']]
receipt_path=SC/'00_admin/download_receipt.json'
old={r['file']:r for r in json.loads(receipt_path.read_text())} if receipt_path.exists() else {}
def download(r):
    url=r['matrix_url'];name=url.rsplit('/',1)[-1];p=SC/'01_raw'/name
    if not p.exists():
        q=p.with_suffix(p.suffix+'.partial');n=0
        with urllib.request.urlopen(url,timeout=60) as src,q.open('wb') as dst:
            expected=src.headers.get('Content-Length')
            while True:
                b=src.read(1<<20)
                if not b:break
                n+=len(b)
                if n>500_000_000:raise ValueError('Individual matrix exceeds 500 MB cap')
                dst.write(b)
        if expected:assert n==int(expected),(name,n,expected)
        assert q.read_bytes()[:2]==b'\x1f\x8b'
        q.rename(p)
    sha=hashlib.sha256(p.read_bytes()).hexdigest()
    if name in old:assert sha==old[name]['sha256']
    with gzip.open(p,'rt') as f:
        rd=csv.reader(f);header=next(rd);row=next(rd)
    return {'file':name,'donor_id':r['donor_id'],'group':r['paper_group'],'technology':r['capture_technology'],'url':url,'bytes':p.stat().st_size,'sha256':sha,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'header_fields':len(header),'first_header_fields':header[:4],'first_data_fields':row[:4]}
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    fs={pool.submit(download,r):r for r in selected}
    for fut in concurrent.futures.as_completed(fs):
        item=fut.result();rows.append(item);print(json.dumps(item),flush=True)
        receipt_path.write_text(json.dumps(sorted(rows,key=lambda r:r['donor_id']),indent=2))
assert len(rows)==8
print('COMPLETE: 8 donor matrices')
