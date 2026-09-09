"""Download the exact data asset and extract it beneath the repository root."""
from pathlib import Path
import hashlib,json,urllib.request,zipfile
root=Path(__file__).resolve().parents[1]
url='https://github.com/DuXC/spermatogenic-laminin-transcriptomics/releases/download/v0.2.0/spermatogenic-laminin-external-data-v0.2.0.zip'
expected='5a639c8f34ffb78c74192ea578b80777e9995a2d3dc84598bbd89e097805f096'
target=root/'spermatogenic-laminin-external-data-v0.2.0.zip'
if not target.exists():
    req=urllib.request.Request(url,headers={'User-Agent':'spermatogenic-laminin-reproduction'})
    with urllib.request.urlopen(req,timeout=60) as src,target.with_suffix('.partial').open('wb') as dst:
        while True:
            b=src.read(1048576)
            if not b:break
            dst.write(b)
    target.with_suffix('.partial').rename(target)
assert hashlib.sha256(target.read_bytes()).hexdigest()==expected,'Asset checksum mismatch'
with zipfile.ZipFile(target) as z:
    for name in z.namelist():
        p=(root/name).resolve()
        assert p.is_relative_to(root),'Unsafe archive path'
    z.extractall(root)
for rel in ['00_admin','logs','env/R_library','07_single_cell/00_admin','07_single_cell/01_raw','07_single_cell/03_processed','07_single_cell/04_results','07_single_cell/05_figures','08_manuscript_refinement_20260909/00_admin']:
    (root/rel).mkdir(parents=True,exist_ok=True)
for rel in ['00_admin','01_metadata','02_raw','03_processed/single_cell','04_results/bulk','04_results/single_cell','05_figures','logs']:
    (root/'10_external_validation_20260910'/rel).mkdir(parents=True,exist_ok=True)
print('Downloaded and extracted verified external data asset')
