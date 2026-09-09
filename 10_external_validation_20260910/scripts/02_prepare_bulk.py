"""Prepare the prespecified external bulk cohorts without fitting outcomes."""
from pathlib import Path
from collections import defaultdict
import gzip, io, json, re, hashlib
import numpy as np
import pandas as pd

S=Path(__file__).resolve().parents[1];P=S.parent;M=S/'01_metadata';O=S/'03_processed'
frozen=json.loads((S/'00_admin/analysis_freeze_receipt.json').read_text())
for r in frozen['files']:assert hashlib.sha256((S/Path(r['file'])).read_bytes()).hexdigest()==r['sha256']
h=pd.read_csv(P/'01_raw/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False)
h=h[h.status=='Approved']; approved=set(h.symbol)
by_id=defaultdict(set);by_symbol=defaultdict(set);by_ref=defaultdict(set)
def tokens(s):return [v for v in re.split(r'\s*///\s*|\||;|,',str(s)) if v and v not in ['---','NA','nan']]
for _,r in h.iterrows():
    if r.entrez_id:by_id[r.entrez_id].add(r.symbol)
    by_symbol[r.symbol].add(r.symbol)
    for c in ['prev_symbol','alias_symbol']:
        for v in tokens(r[c]):by_symbol[v].add(r.symbol)
    for v in tokens(r.refseq_accession):by_ref[v.split('.')[0]].add(r.symbol)
def resolve(ids,syms,refs):
    a=set().union(*(by_id[x] for x in tokens(ids)))
    b=set().union(*({x} if x in approved else by_symbol[x] for x in tokens(syms)))
    c=set().union(*(by_ref[x.split('.')[0]] for x in tokens(refs)))
    gs=a or b or c
    method='Entrez' if a else 'symbol' if b else 'RefSeq' if c else 'unmapped'
    if len(gs)!=1:return '',method,'unmapped' if not gs else 'ambiguous'
    gene=next(iter(gs))
    if a and b and gene not in b:return '',method,'Entrez_symbol_conflict'
    return gene,method,'mapped_unique'

audit={};coverage=[];sets=json.loads((M/'fixed_gene_sets.json').read_text())
allsets={**sets['primary'],**sets['exploratory_markers']}
for acc in ['GSE9210','GSE108886']:
    meta=json.loads((M/(acc+'_metadata.json')).read_text())
    text=gzip.open(M/(acc+'_family.soft.gz'),'rt').read()
    platform=pd.read_csv(io.StringIO(text.split('!platform_table_begin\n')[1].split('!platform_table_end')[0]),sep='\t',dtype=str,keep_default_na=False).set_index('ID')
    blocks=text.split('!sample_table_begin\n')[1:]
    assert len(blocks)==len(meta['samples'])
    matrices=[];samples=[];det=[]
    for source,block in zip(meta['samples'],blocks):
        tab=pd.read_csv(io.StringIO(block.split('!sample_table_end')[0]),sep='\t',dtype={'ID_REF':str}).set_index('ID_REF')
        gsm=source['gsm'];title=source['title'][0]
        chars={v.split(':',1)[0].strip():v.split(':',1)[1].strip() for v in source.get('characteristics_ch1',[]) if ':' in v}
        group=('NOA' if '_NOA_' in title else 'OA') if acc=='GSE9210' else ('NOA' if chars['disease state'].startswith('Non obstructive') else 'OA' if chars['disease state']=='Obstructive azoospermia' else 'pooled_control')
        included=group!='pooled_control'
        samples.append({'cohort':acc,'gsm':gsm,'title':title,'group':group,'included':included,
                        'source':source.get('source_name_ch1',[''])[0], 'original_diagnosis':chars.get('disease state','; '.join(source.get('characteristics_ch1',[]))),
                        'unit':'individual patient' if included else 'pooled reference excluded', 'platform':source['platform_id'][0]})
        matrices.append(tab.VALUE.rename(gsm))
        if 'DETECTION P-VALUE' in tab:det.append(tab['DETECTION P-VALUE'].rename(gsm))
    sm=pd.DataFrame(samples);sm.to_csv(M/(acc+'_sample_audit.tsv'),sep='\t',index=False)
    sm=sm[sm.included].reset_index(drop=True)
    mat=pd.concat(matrices,axis=1).loc[:,sm.gsm]
    assert mat.index.is_unique and mat.columns.is_unique
    platform=platform.reindex(mat.index)
    assert not platform.isna().all(axis=1).any()
    platform.to_csv(M/(acc+'_original_platform.tsv.gz'),sep='\t',compression='gzip')
    if acc=='GSE9210':
        ids=platform.GENE;syms=platform.GENE_SYMBOL;refs=platform.REFSEQ;biological=platform.CONTROL_TYPE.eq('FALSE')
        # VALUE is ln(Cy5/Cy3), verified against the native normalized channel intensities.
        mat=mat/np.log(2)
        detection=pd.Series(True,index=mat.index)
        unit='log2 common-reference ratio'
    else:
        ids=platform.Entrez_Gene_ID;syms=platform.Symbol;refs=platform.RefSeq_ID
        biological=platform.Species.str.lower().eq('homo sapiens') & (platform.Entrez_Gene_ID.ne('')|platform.Symbol.ne(''))
        pdet=pd.concat(det,axis=1).loc[:,sm.gsm]
        detection=pd.Series(False,index=mat.index)
        for group in ['OA','NOA']:
            columns=sm.loc[sm.group==group,'gsm'];detection|=(pdet[columns]<.01).sum(axis=1)>=np.ceil(len(columns)/2)
        pdet.to_csv(O/(acc+'_detection_p.tsv.gz'),sep='\t',compression='gzip')
        unit='deposited normalized signal'
    mp=pd.DataFrame({'probe_id':mat.index,'original_entrez':ids.values,'original_symbol':syms.values,'original_refseq':refs.values})
    mp[['gene','mapping_method','mapping_status']]=pd.DataFrame([resolve(i,j,k) for i,j,k in zip(ids,syms,refs)])
    mp['biological_probe']=biological.values;mp['complete']=mat.notna().all(axis=1).values;mp['detected']=detection.values
    present80=pd.Series(True,index=mat.index)
    for group in ['OA','NOA']:
        cols=sm.loc[sm.group==group,'gsm'];present80&=mat[cols].notna().sum(axis=1)>=np.ceil(.8*len(cols))
    mp['present_80pct_each_group']=present80.values
    mapped=mp.mapping_status.eq('mapped_unique')&mp.biological_probe
    mp['primary_eligible']=mapped&mp.complete&mp.detected
    sensitivity=mapped&mp.present_80pct_each_group if acc=='GSE9210' else mapped&mp.complete
    mp['sensitivity_eligible']=sensitivity
    mp.to_csv(M/(acc+'_probe_mapping.tsv.gz'),sep='\t',index=False,compression='gzip')
    variants={'primary':mp.primary_eligible,'available_case' if acc=='GSE9210' else 'all_mapped':sensitivity}
    for variant,keep in variants.items():
        selected=mat.loc[keep.values].copy();selected.index=mp.loc[keep,'gene']
        # Average observed probes equally; no missing expression values are filled.
        gene=selected.groupby(level=0,sort=True).mean()
        gene=gene.loc[gene.std(axis=1,skipna=True)>0]
        if variant=='primary':assert not gene.isna().any().any()
        gene.to_csv(O/(acc+'_'+variant+'_gene_expression.tsv.gz'),sep='\t',index_label='gene',compression='gzip')
        assert len(gene)>1000
        for name,members in allsets.items():
            measured=sorted(set(members)&set(gene.index))
            coverage.append({'cohort':acc,'variant':variant,'program':name,'n_total':len(members),'n_measured':len(measured),'genes':'|'.join(measured),'missing':'|'.join(sorted(set(members)-set(measured)))})
    sm['gene_effect_unit']=unit;sm.to_csv(M/(acc+'_analysis_samples.tsv'),sep='\t',index=False)
    audit[acc]={'patients':len(sm),'groups':sm.group.value_counts().to_dict(),'source_probes':len(mp),'primary_probes':int(mp.primary_eligible.sum()),
                'sensitivity_probes':int(mp.sensitivity_eligible.sum()),'mapping_status':mp.mapping_status.value_counts().to_dict(),'gene_effect_unit':unit}
    print(acc,json.dumps(audit[acc]),flush=True)
pd.DataFrame(coverage).to_csv(M/'external_program_coverage.tsv',sep='\t',index=False)
(S/'00_admin/bulk_preparation_QA.json').write_text(json.dumps(audit,indent=2)+'\n')

# Record explicit native reanalysis links rather than assuming independence from new accessions.
a=json.loads((M/'GSE45885_metadata.json').read_text())['samples'];b=json.loads((M/'GSE45887_metadata.json').read_text())['samples']
links=[]
for source in a:
    for rel in source.get('relation',[]):
        if rel.startswith('Reanalyzed by: '):
            target=rel.split(': ',1)[1]
            if target in {x['gsm'] for x in b}:links.append({'original_gsm':source['gsm'],'reanalysis_gsm':target,'evidence':'native GEO Sample_relation'})
pd.DataFrame(links).to_csv(M/'GSE45885_GSE45887_reuse.tsv',sep='\t',index=False)
print('EXPLICIT_REANALYSIS_LINKS',len(links),flush=True)
