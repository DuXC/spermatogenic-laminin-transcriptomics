#!/usr/bin/env python3
"""Chunked original-count QC, donor-aware adult annotation embedding.

Run from any directory using the isolated SC environment. No bulk file is edited.
"""
from pathlib import Path
from collections import defaultdict
import os, json, re, time, hashlib, importlib.metadata
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import scanpy as sc
import harmonypy

S = Path(__file__).resolve().parents[1]
P = S.parent
CACHE = (P / '.cache' / 'single_cell')
CACHE.mkdir(parents=True, exist_ok=True)
SEED = 20260909
sc.settings.n_jobs = 4
sc.settings.verbosity = 2

CURATED = {
 'Sertoli': 'SOX9 WT1 AMH FSHR CLU INHA GATA4 KRT18 KRT8 SOX8 DMRT1 FATE1',
 'Leydig': 'INSL3 STAR CYP17A1 HSD3B2 CYP11A1 LHCGR NR5A1',
 'Peritubular_myoid': 'DPEP1 TSHZ2 OSR2 FHL2 PTGDS CRISPLD2 C7 MEG3',
 'Vascular_mural': 'RGS5 RERGL MCAM NOTCH3 CSPG4 PLN ADIRF MYL9',
 'Shared_contractile': 'ACTA2 TAGLN MYH11 CNN1 DES TPM2',
 'Endothelial': 'PECAM1 VWF KDR EMCN CLDN5 RAMP2 GNG11 SELE',
 'Macrophage': 'LST1 TYROBP FCER1G AIF1 CD68 CSF1R C1QA C1QB C1QC CD74',
 'T_NK': 'CD3D CD3E TRAC NKG7 GNLY KLRD1',
 'B_cell': 'MS4A1 CD79A CD79B CD37',
 'Spermatogonia': 'UTF1 GFRA1 ZBTB16 FGFR3 PIWIL4 UCHL1 MAGEA4 KIT SALL4',
 'Spermatocyte': 'SYCP1 SYCP2 SYCP3 DMC1 HORMAD1 SPO11 MEIOB REC8',
 'Spermatid': 'PRM1 PRM2 TNP1 TNP2 ACRV1 SPACA1 IZUMO1 ODF1 ODF2',
}
PRIMARY = json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary']
ECM = set().union(*map(set, PRIMARY.values()))
CURATED = {k: [x for x in v.split() if x not in ECM] for k,v in CURATED.items()}

def say(s):
    print(time.strftime('%H:%M:%S'), s, flush=True)

def symbol_mapper():
    h = pd.read_csv(P/'01_raw/hgnc_complete_set.txt', sep='\t', dtype=str).fillna('')
    h = h[h.status.eq('Approved')]
    official = set(h.symbol)
    transformed, aliases = defaultdict(set), defaultdict(set)
    for _, r in h.iterrows():
        sym = r.symbol
        transformed[re.sub('[^A-Za-z0-9_.]', '.', sym)].add(sym)
        for fld in ['prev_symbol','alias_symbol']:
            for term in r[fld].split('|'):
                if term:
                    aliases[term].add(sym)
                    aliases[re.sub('[^A-Za-z0-9_.]', '.', term)].add(sym)
    def lookup(x):
        if x in official: return x, 'approved_exact'
        t = transformed.get(x, set())
        if len(t) == 1: return next(iter(t)), 'approved_R_name'
        if len(t) > 1: return '', 'ambiguous_approved_R_name'
        t = aliases.get(x, set())
        if len(t) == 1: return next(iter(t)), 'unique_alias_or_previous'
        return '', 'ambiguous_alias' if t else 'unmapped'
    return lookup

def reference_panels(lookup):
    df = pd.read_excel(P/'references/GSE149512_Supplementary_Data_1.xlsx')
    labels = {'endotheliocyte':'Endothelial','Leydig_cells':'Leydig',
              'macrophages':'Macrophage','PTM_cells':'Peritubular_myoid',
              'VSM_cells':'Vascular_mural','Sertoli_cells':'Sertoli',
              'spermatid':'Spermatid','spermatocyte':'Spermatocyte','spermatogonia':'Spermatogonia'}
    ref = {}
    for raw, g in df.groupby('cluster'):
        out = []
        for x in g.sort_values('avg_logFC', ascending=False).gene:
            sym, _ = lookup(x)
            if not sym or sym in ECM or re.match(r'^(MT-|RPL|RPS)', sym) or sym in out: continue
            out.append(sym)
            if len(out) == 30: break
        ref[labels[raw]] = out
    return ref

def prepare_one(row, lookup):
    donor = row.donor_id
    cached = CACHE/f'{donor}_counts_qc.h5ad'
    if cached.exists():
        say('reuse verified cache '+donor)
        return ad.read_h5ad(cached)
    path = S/'01_raw'/f'{row.gsm}_{donor}matrix.csv.gz'
    say('parse '+path.name)
    chunks, names, cells = [], [], None
    max_value = 0
    for ch in pd.read_csv(path, index_col=0, chunksize=128):
        if cells is None: cells = ch.columns.astype(str).to_numpy()
        assert np.array_equal(cells, ch.columns.astype(str).to_numpy())
        v = ch.to_numpy()
        assert np.isfinite(v).all() and (v >= 0).all() and np.equal(v, np.floor(v)).all()
        max_value = max(max_value, int(v.max()))
        assert max_value < np.iinfo(np.int32).max
        chunks.append(sparse.csr_matrix(v, dtype=np.int32))
        names.extend(ch.index.astype(str).tolist())
    X = sparse.vstack(chunks, format='csc').T.tocsr()
    del chunks
    assert len(set(cells)) == len(cells)
    nfeat = np.diff(X.indptr)
    counts = np.asarray(X.sum(axis=1)).ravel()
    mito = np.array([bool(re.match(r'^MT[-.]',x)) for x in names])
    pct = np.asarray(X[:,mito].sum(axis=1)).ravel()/np.maximum(counts,1)*100
    obs = pd.DataFrame(index=pd.Index(cells, name='cell_id'))
    obs['donor_id'] = donor
    obs['group'] = 'OA' if row.paper_group == 'adult_normal_spermatogenesis_OA' else 'iNOA'
    obs['technology'] = row.capture_technology
    obs['age_years'] = row.age_years
    obs['n_features_original'] = nfeat
    obs['total_counts_original'] = counts
    obs['pct_mito_original'] = pct
    obs['fail_low_features'] = nfeat <= 500
    obs['fail_high_features'] = nfeat >= 9000
    obs['fail_high_counts'] = counts >= 80000
    obs['fail_mito_40'] = pct >= 40
    obs['primary_qc_pass'] = ~(obs[['fail_low_features','fail_high_features','fail_high_counts','fail_mito_40']].any(axis=1))
    obs['mt20_pass'] = obs.primary_qc_pass & (pct < 20)
    obs.to_csv(S/'02_annotation'/f'{donor}_all_cell_QC.tsv.gz', sep='\t')
    mapping = pd.DataFrame({'original_feature': names})
    mapping[['symbol','mapping_status']] = [lookup(x) for x in names]
    mapping.to_csv(S/'02_annotation'/f'{donor}_feature_mapping.tsv.gz', sep='\t', index=False)
    approved = sorted(set(mapping.symbol)-{''})
    idx = {x:i for i,x in enumerate(approved)}
    valid = mapping.symbol.ne('').to_numpy()
    M = sparse.coo_matrix((np.ones(valid.sum(),dtype=np.int32),
          (np.flatnonzero(valid), [idx[x] for x in mapping.symbol[valid]])),
          shape=(len(names),len(approved))).tocsr()
    X = (X[obs.primary_qc_pass.to_numpy()] @ M).tocsr()
    obs = obs[obs.primary_qc_pass].copy()
    obs['total_counts_mapped'] = np.asarray(X.sum(axis=1)).ravel()
    a = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=pd.Index(approved,name='symbol')))
    try:
        sc.pp.scrublet(a, expected_doublet_rate=.06, random_state=SEED, verbose=True)
        a.obs['scrublet_status'] = 'complete'
        status = {'donor_id':donor,'status':'complete',
                  'threshold':float(a.uns['scrublet']['threshold']),
                  'predicted_doublets':int(a.obs.predicted_doublet.sum())}
    except Exception as e:
        a.obs['doublet_score'] = np.nan
        a.obs['predicted_doublet'] = False
        a.obs['scrublet_status'] = 'failed'
        status = {'donor_id':donor,'status':'failed','error':repr(e)}
    status.update({'raw_cells':len(cells),'primary_qc_cells':a.n_obs,'raw_features':len(names),
                   'mapped_genes':a.n_vars,'unmapped_features':int((~valid).sum()),
                   'mitochondrial_features':int(mito.sum()),'max_original_count':max_value})
    (S/'00_admin'/f'{donor}_prepare_receipt.json').write_text(json.dumps(status,indent=2))
    a.write_h5ad(cached, compression='gzip')
    say(str(status))
    return a

def main():
    lookup = symbol_mapper()
    ref = reference_panels(lookup)
    panels = {'curated_'+k:v for k,v in CURATED.items()} | {'paper_'+k:v for k,v in ref.items()}
    assert not ECM.intersection(set().union(*map(set,panels.values())))
    (S/'02_annotation/annotation_marker_panels.json').write_text(json.dumps(panels,indent=2))
    d = pd.read_csv(P/'02_annotation/GSE149512_donor_metadata_audit.tsv', sep='\t')
    # SI p.2 clinical table re-inspected in this extension: LZ011 total 7896.
    d['prior_release_total_cells'] = d.published_total_cells
    d.loc[d.donor_id.eq('LZ011'),'published_total_cells'] = 7896
    d['sc_extension_use'] = 'separate biological stratum'
    keep = d.paper_group.isin(['adult_normal_spermatogenesis_OA','iNOA'])
    d.loc[keep,'sc_extension_use'] = 'adult source extension'
    d.loc[keep,'matrix_analysis_status'] = 'downloaded; processed in 07_single_cell'
    d.to_csv(S/'02_annotation/GSE149512_donor_metadata_v2.tsv',sep='\t',index=False)
    arrays = [prepare_one(r,lookup) for _,r in d[keep].iterrows()]
    a = ad.concat(arrays,join='outer',merge='first',fill_value=0)
    del arrays
    assert a.obs_names.is_unique
    a.X = a.X.astype(np.float32)
    a.layers['counts'] = a.X.astype(np.int32)
    a.var['n_cells'] = np.asarray((a.X > 0).sum(axis=0)).ravel()
    sc.pp.normalize_total(a,target_sum=1e4)
    sc.pp.log1p(a)
    all_markers = set().union(*map(set,panels.values()))
    eligible = np.array([x not in ECM and x not in all_markers and not re.match(r'^(MT-|RPL|RPS)',x) for x in a.var_names]) & (a.var.n_cells.to_numpy() >= 3)
    temp = a[:,eligible].copy()
    sc.pp.highly_variable_genes(temp,n_top_genes=3000,flavor='seurat',batch_key='donor_id')
    selected = temp.var_names[temp.var.highly_variable]
    a.var['highly_variable'] = a.var_names.isin(selected)
    a.var.to_csv(S/'02_annotation/gene_coverage_and_hvg.tsv.gz',sep='\t')
    b = temp[:,temp.var.highly_variable].copy()
    del temp
    b.layers.clear()
    sc.pp.scale(b,max_value=10)
    sc.tl.pca(b,n_comps=30,svd_solver='arpack',random_state=SEED)
    a.obsm['X_pca_uncorrected'] = b.obsm['X_pca'].copy()
    a.uns['pca_variance_ratio'] = b.uns['pca']['variance_ratio']
    del b
    say('Harmony annotation embedding')
    ho = harmonypy.run_harmony(a.obsm['X_pca_uncorrected'], a.obs, 'donor_id',
                             max_iter_harmony=20,random_state=SEED)
    z = np.asarray(ho.Z_corr)
    a.obsm['X_pca_harmony'] = z if z.shape[0] == a.n_obs else z.T
    sc.pp.neighbors(a,n_neighbors=30,use_rep='X_pca_harmony',random_state=SEED)
    sc.tl.leiden(a,resolution=1.0,random_state=SEED,key_added='leiden',flavor='igraph',n_iterations=2,directed=False)
    sc.tl.umap(a,random_state=SEED)
    score_rows = []
    for label, genes in panels.items():
        measured = [x for x in genes if x in a.var_names]
        if len(measured) >= 2:
            sc.tl.score_genes(a,measured,score_name=label,random_state=SEED,ctrl_size=50,
                              gene_pool=[x for x in a.var_names if x not in ECM])
        else:
            a.obs[label] = np.nan
        score_rows.append({'panel':label,'expected':len(genes),'measured':len(measured),'genes':';'.join(measured)})
    pd.DataFrame(score_rows).to_csv(S/'02_annotation/annotation_marker_coverage.tsv',sep='\t',index=False)
    a.obs.to_csv(S/'02_annotation/cells_preannotation.tsv.gz',sep='\t')
    cols = list(panels)
    a.obs.groupby('leiden',observed=True)[cols].mean().to_csv(S/'04_results/cluster_marker_scores.tsv',sep='\t')
    pd.crosstab(a.obs.leiden,a.obs.donor_id).to_csv(S/'04_results/cluster_donor_counts.tsv',sep='\t')
    # Cluster markers support annotation only; cell-level p-values are not patient evidence.
    sc.tl.rank_genes_groups(a,groupby='leiden',method='wilcoxon',n_genes=40,use_raw=False)
    sc.get.rank_genes_groups_df(a,group=None).to_csv(S/'04_results/cluster_annotation_markers.tsv.gz',sep='\t',index=False)
    pd.DataFrame(a.obsm['X_umap'],index=a.obs_names,columns=['UMAP1','UMAP2']).to_csv(S/'03_processed/umap.tsv.gz',sep='\t')
    a.write_h5ad(CACHE/'adult_preannotation.h5ad',compression='gzip')
    a.write_h5ad(S/'03_processed/adult_preannotation.h5ad',compression='gzip')
    versions = {p:importlib.metadata.version(p) for p in ['scanpy','anndata','numpy','scipy','pandas','scikit-learn','igraph','leidenalg','harmonypy']}
    (S/'00_admin/software_versions.json').write_text(json.dumps(versions,indent=2))
    say(f'COMPLETE: {a.n_obs} cells, {a.n_vars} mapped genes, {a.obs.leiden.nunique()} clusters')

if __name__ == '__main__': main()
