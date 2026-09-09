from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import gzip, hashlib, json, urllib.parse, urllib.request, sys

S = Path(__file__).resolve().parents[1]
out = S / '01_metadata'
receipts = []

def download(url, path):
    if path.exists():
        payload = path.read_bytes()
    else:
        req = urllib.request.Request(url, headers={'User-Agent': 'Scientific-public-data-reanalysis/1.0'})
        with urllib.request.urlopen(req, timeout=60) as r:
            assert r.status == 200
            payload = r.read()
        path.write_bytes(payload)
    return payload

def fetch(acc):
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz'
    path = out / (acc + '_family.soft.gz')
    try:
        payload = download(url, path)
        raw = gzip.decompress(payload).decode('utf-8', errors='replace')
        assert '^SERIES = ' + acc in raw
        # Extract record metadata only; expression tables are not inspected at this stage.
        metadata = '\n'.join(line for line in raw.splitlines() if line.startswith(('^', '!')) and not line.startswith(('!sample_table_', '!platform_table_')))
        (out / (acc + '_metadata.soft')).write_text(metadata + '\n')
        samples = []
        current = None
        for line in metadata.splitlines():
            if line.startswith('^SAMPLE = '):
                current = {'gsm': line.split(' = ',1)[1]}; samples.append(current)
            elif line.startswith('^'):
                current = None
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                key, value = line.split(' = ', 1)
                current.setdefault(key[8:], []).append(value)
        series = {}
        for line in metadata.splitlines():
            if line.startswith('!Series_') and ' = ' in line:
                key, value = line.split(' = ', 1); series.setdefault(key[8:], []).append(value)
        detail = {'accession': acc, 'series': series, 'samples': samples}
        (out / (acc + '_metadata.json')).write_text(json.dumps(detail, ensure_ascii=False, indent=2)+'\n')
        return {'accession': acc, 'status': 'VERIFIED_METADATA', 'url': url, 'bytes': len(payload),
                'sha256': hashlib.sha256(payload).hexdigest(), 'samples': len(samples),
                'title': series.get('title'), 'design': series.get('overall_design'),
                'pmid': series.get('pubmed_id')}
    except Exception as exc:
        return {'accession': acc, 'status': 'ERROR', 'url': url, 'error': str(exc)}

extra = len(sys.argv) > 1
accessions = sys.argv[1:] if extra else ['GSE45885', 'GSE9210', 'GSE45887', 'GSE108886', 'GSE154535', 'GSE202647', 'GSE106487', 'GSE157421']
with ThreadPoolExecutor(max_workers=3) as pool:
    for receipt in pool.map(fetch, accessions):
        receipts.append(receipt); print(json.dumps(receipt, ensure_ascii=False), flush=True)

query = '(testis OR testicular) AND ("single cell" OR "single-cell" OR scRNA-seq) AND (azoospermia OR cryptozoospermia OR "spermatogenic failure") AND "Homo sapiens"[Organism]'
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urllib.parse.urlencode({'db':'gds','term':query,'retmode':'json','retmax':100})
try:
    payload = download(url, out / 'GEO_native_search.json')
    record = json.loads(payload); ids = record['esearchresult']['idlist']
    summary_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?' + urllib.parse.urlencode({'db':'gds','id':','.join(ids),'retmode':'json'})
    download(summary_url, out / 'GEO_native_search_summaries.json')
    receipts.append({'status':'NATIVE_SEARCH_SAVED','query':query,'count':record['esearchresult']['count'],'ids':len(ids),'url':url})
except Exception as exc:
    receipts.append({'status':'SEARCH_ERROR','query':query,'error':str(exc)})
(S/('00_admin/metadata_extension_receipt.json' if extra else '00_admin/metadata_acquisition_receipt.json')).write_text(json.dumps({'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'receipts':receipts}, ensure_ascii=False, indent=2)+'\n')
print(json.dumps(receipts[-1], ensure_ascii=False), flush=True)
