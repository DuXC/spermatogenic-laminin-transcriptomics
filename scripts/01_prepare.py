"""Deterministic GEO parsing, current-symbol resolution, and frozen gene sets."""
from pathlib import Path
import csv,gzip,io,json,re,hashlib,datetime
from collections import defaultdict
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'01_raw';ANN=ROOT/'02_annotation';PRO=ROOT/'03_processed'
def soft_table(p):
    opener=gzip.open if p.suffix=='.gz' else open
    lines=[];inside=False
    with opener(p,'rt') as f:
        for l in f:
            if l.startswith('!platform_table_begin'):inside=True;continue
            if l.startswith('!platform_table_end'):break
            if inside:lines.append(l)
    return pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str,keep_default_na=False).set_index('ID')
def metadata(p):
    out=defaultdict(list)
    with gzip.open(p,'rt') as f:
        for row in csv.reader(f,delimiter='\t'):
            if row and row[0]=='!series_matrix_table_begin':break
            if row and row[0].startswith('!Sample_'):out[row[0][8:]].append(row[1:])
    return out
def tokens(s):return [x for x in re.split(r'\s*///\s*|\||;|,',str(s)) if x and x not in ['---','NA','nan']]
h=pd.read_csv(RAW/'hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False)
h=h[h.status=='Approved']; by_id=defaultdict(set);by_symbol=defaultdict(set);by_ref=defaultdict(set)
approved=set(h.symbol)
for _,r in h.iterrows():
    if r.entrez_id:by_id[r.entrez_id].add(r.symbol)
    by_symbol[r.symbol].add(r.symbol)
    for col in ['prev_symbol','alias_symbol']:
        for v in tokens(r[col]):by_symbol[v].add(r.symbol)
    for v in tokens(r.refseq_accession):by_ref[v.split('.')[0]].add(r.symbol)
def resolve(ids,syms,refs):
    idset=set().union(*(by_id[x] for x in tokens(ids)))
    # Approved current symbols take precedence over aliases of the same string.
    symset=set().union(*({x} if x in approved else by_symbol[x] for x in tokens(syms)))
    refset=set().union(*(by_ref[x.split('.')[0]] for x in tokens(refs)))
    chosen=idset or symset or refset
    method='Entrez_HGNC' if idset else 'symbol_HGNC' if symset else 'RefSeq_HGNC' if refset else 'unmapped'
    if len(chosen)>1:return '',method,'ambiguous_multiple_current_genes','|'.join(sorted(chosen))
    if not chosen:return '',method,'unmapped',''
    g=next(iter(chosen))
    if idset and symset and g not in symset:return '',method,'Entrez_symbol_conflict','|'.join(sorted(idset|symset))
    return g,method,'mapped_unique',g

audit={};all_samples=[]
for acc,plat in [('GSE4797','GPL2891'),('GSE145467','GPL4133')]:
    meta=metadata(RAW/(acc+'_series_matrix.txt.gz'))
    samp=pd.DataFrame({'gsm':meta['geo_accession'][0],'title':meta['title'][0]})
    for field in ['source_name_ch1','description','platform_id','label_ch1','label_ch2']:
        if field in meta:samp[field]=meta[field][0]
    for j,v in enumerate(meta['characteristics_ch1']):samp['characteristics_ch1_'+str(j+1)]=v
    samp['cohort']=acc
    if acc=='GSE4797':
        samp['johnsen']=samp.title.str.split('.').str[0].astype(int)
        samp['group']='JS'+samp.johnsen.astype(str)
        samp['series_pathology']=samp.johnsen.map({10:'full spermatogenesis',8:'arrest at the spermatid stage',5:'arrest at spermatocyte stage',2:'Sertoli-cell-only syndrome'})
        samp['clinical_NOA_status']='not assigned from histology'
        samp['sample_unit']='one biopsy per man; GEO series states 28 men'
    else:
        samp['group']=np.where(samp.title.str.contains('impaired'),'impaired','normal')
        samp['johnsen']=np.nan
        samp['paper_patient_id']=''
        samp['clinical_link_status']='UNLINKED: paper patient numbers differ from GEO Sample numbers'
        samp['sample_unit']='one expression column per study sample; n=20'
    samp.to_csv(ANN/(acc+'_samples.tsv'),sep='\t',index=False)
    all_samples.append(samp)
    mat=pd.read_csv(RAW/(acc+'_series_matrix.txt.gz'),sep='\t',comment='!',index_col=0)
    mat.index=mat.index.astype(str); assert list(mat.columns)==list(samp.gsm)
    assert mat.index.is_unique and samp.gsm.is_unique
    p=soft_table(RAW/(plat+'_table.txt'))
    assert p.index.is_unique
    p=p.reindex(mat.index)
    if acc=='GSE4797':
        auto=soft_table(RAW/'GPL2891.annot.gz').reindex(mat.index).fillna('')
        mp=pd.DataFrame({'probe_id':mat.index,'original_symbol':auto['Gene symbol'].values,'original_entrez':auto['Gene ID'].values,'original_refseq_genbank':p.GB_ACC.values,'probe_type':p.PROBE_TYPE.values,'target_scope':p.PUB_PROBE_TARGETS.values})
        # Retain every original platform column and every historical automatic annotation.
        p.to_csv(ANN/(plat+'_original.tsv.gz'),sep='\t',compression='gzip')
        auto.to_csv(ANN/(plat+'_GEO_20071123_annotation.tsv.gz'),sep='\t',compression='gzip')
    else:
        mp=pd.DataFrame({'probe_id':mat.index,'original_symbol':p.GENE_SYMBOL.values,'original_entrez':p.GENE.values,'original_refseq_genbank':p.REFSEQ.values,'probe_type':p.CONTROL_TYPE.values})
        p.to_csv(ANN/(plat+'_original.tsv.gz'),sep='\t',compression='gzip')
    resolved=[resolve(r.original_entrez,r.original_symbol,r.original_refseq_genbank) for _,r in mp.iterrows()]
    mp[['gene','mapping_method','mapping_status','mapping_candidates']]=pd.DataFrame(resolved,index=mp.index)
    mp['complete_all_samples']=mat.notna().all(axis=1).values
    if acc=='GSE4797':
        above=[]
        for group in ['JS10','JS8','JS5','JS2']:
            sub=mat.loc[:,samp.loc[samp.group==group,'gsm']]
            above.append((sub>=.2).sum(axis=1)>=np.ceil(sub.shape[1]/2))
        mp['passes_low_signal_filter']=np.logical_or.reduce(above)
    else:mp['passes_low_signal_filter']=True
    mp['analysis_eligible']=(mp.mapping_status=='mapped_unique') & mp.complete_all_samples & mp.passes_low_signal_filter
    if acc=='GSE145467':mp['analysis_eligible'] &= mp.probe_type.str.upper().eq('FALSE')
    else:mp['analysis_eligible'] &= mp.probe_type.eq('DISCOVERY')
    mp.to_csv(ANN/(acc+'_probe_to_current_gene.tsv.gz'),sep='\t',index=False,compression='gzip')
    mat.to_csv(PRO/(acc+'_uploaded_probe_values.tsv.gz'),sep='\t',index_label='probe_id',compression='gzip')
    audit[acc]={'probes':len(mp),'samples':len(samp),'groups':samp.group.value_counts().to_dict(),'mapping_status':mp.mapping_status.value_counts().to_dict(),'eligible_probes':int(mp.analysis_eligible.sum()),'eligible_genes':int(mp.loc[mp.analysis_eligible,'gene'].nunique()),'incomplete_probes':int((~mp.complete_all_samples).sum()),'low_signal_probes':int((~mp.passes_low_signal_filter).sum()),'missing_entries':int(mat.isna().sum().sum()),'nonpositive_entries':int((mat<=0).sum().sum()),'quantiles':np.nanquantile(mat,[0,.01,.25,.5,.75,.99,1]).tolist()}
    print(acc,json.dumps(audit[acc]),flush=True)
pd.concat(all_samples,ignore_index=True).to_csv(ANN/'all_samples.tsv',sep='\t',index=False)

fixed=['GOCC_BASEMENT_MEMBRANE','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_LAMININ_INTERACTIONS','REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
sets={};rows=[]
for p in [RAW/'c2.cp.reactome.v2025.1.Hs.symbols.gmt',RAW/'c5.go.cc.v2025.1.Hs.symbols.gmt']:
    for line in p.read_text().splitlines():
        name,url,*members=line.split('\t')
        if name not in fixed:continue
        mapped=set()
        for m in members:
            gene,method,status,candidates=resolve('',m,'')
            rows.append(dict(set=name,source_gene=m,current_gene=gene,status=status,source_url=url,collection=p.name))
            if gene:mapped.add(gene)
        sets[name]=sorted(mapped)
assert set(sets)==set(fixed)
markers={'Sertoli_markers':['SOX9','WT1','FSHR','CLDN11','GATA4'],'Leydig_markers':['STAR','CYP11A1','CYP17A1','HSD3B2','INSL3'],'Peritubular_markers':['ACTA2','MYH11','TAGLN','CNN1','DES'],'Meiotic_markers':['SYCP1','SYCP3','DMC1','SPO11','MEIOB'],'Postmeiotic_markers':['PRM1','PRM2','TNP1','TNP2','ACRV1']}
pd.DataFrame(rows).to_csv(ANN/'fixed_pathway_membership.tsv',sep='\t',index=False)
(ANN/'fixed_gene_sets.json').write_text(json.dumps({'primary':sets,'exploratory_markers':markers},indent=2))
(ROOT/'00_admin/preparation_audit.json').write_text(json.dumps(audit,indent=2))
for p in [ANN/'fixed_gene_sets.json', ROOT/'data_manifest.json']:
    print('FROZEN',p.name,hashlib.sha256(p.read_bytes()).hexdigest())
