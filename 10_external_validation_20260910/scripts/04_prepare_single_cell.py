"""Aggregate original Alevin counts by donor, preserving published cell identities."""
from pathlib import Path
from collections import defaultdict
import gzip, json, hashlib, re, sys, gc
import numpy as np
import pandas as pd
from scipy import sparse

S=Path(__file__).resolve().parents[1];P=S.parent;O=S/'03_processed/single_cell';O.mkdir(exist_ok=True)
h=pd.read_csv(P/'01_raw/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False);h=h[h.status=='Approved']
approved=set(h.symbol);aliases=defaultdict(set)
for _,r in h.iterrows():
    for col in ['prev_symbol','alias_symbol']:
        for token in r[col].split('|'):
            if token:aliases[token].add(r.symbol)
def resolve(g):
    if g in approved:return g,'approved'
    vals=aliases[g]
    return (next(iter(vals)),'unique_alias') if len(vals)==1 else ('','ambiguous' if vals else 'unmapped')
meta=pd.read_csv(S/'01_metadata/GSE153947_Cell_metadata.tsv.gz',sep='\t',index_col=0)
assert len(meta)==28666 and meta.index.is_unique
sets=json.loads((S/'01_metadata/fixed_gene_sets.json').read_text());ecm=set().union(*map(set,sets['primary'].values()))
markers={
 'contractile_mural':['ACTA2','MYH11','TAGLN','CNN1','DES'],
 'perivascular':['RGS5','CSPG4','MCAM','NOTCH3'],
 'Sertoli':['SOX9','WT1','FSHR','CLDN11','GATA4'],
 'Leydig':['STAR','CYP11A1','CYP17A1','HSD3B2','INSL3'],
 'endothelial':['PECAM1','VWF','EMCN','KDR'],
 'immune':['PTPRC','LST1','TYROBP'],
 'germline':['DDX4','DAZL','MAGEA4','PRM1','SYCP3']}
marker_genes=set().union(*map(set,markers.values()))-ecm
source_genes=ecm|marker_genes
compartments={
 'Peritubular_compartment':['PMCs','Fibrotic peritubular myoid cells'],
 'Sertoli':['Sertoli cells'],'Leydig':['Leydig cells'],'Endothelial':['Endothelial cells'],
 'Macrophage':['Macrophages'],'Perivascular':['Perivascular cells'],
 'Author_PMC_state':['PMCs'],'Author_fibrotic_PMC_state':['Fibrotic peritubular myoid cells']}
donors=sys.argv[1:] or ['N1','N2','N3','Cr1','Cr2','Cr3']
for donor in donors:
    ds=meta[meta.Sample_ID==donor].copy();ds['source_cell_id']=ds.index;ds['raw_barcode']=ds.index.str.split('_').str[0]
    assert ds.raw_barcode.is_unique
    paths=list((S/'02_raw').glob('*_'+donor+'_counts.tsv.gz'));assert len(paths)==1;path=paths[0]
    with gzip.open(path,'rt') as f:barcodes=f.readline().rstrip('\n').split('\t')
    assert len(set(barcodes))==len(barcodes)
    positions=pd.Index(barcodes).get_indexer(ds.raw_barcode);assert min(positions)>=0
    chunks=[];original_genes=[]
    iterator=pd.read_csv(path,sep='\t',header=None,skiprows=1,names=['gene']+barcodes,index_col=0,chunksize=400,
                         dtype={**{b:'float64' for b in barcodes},'gene':str})
    for k,block in enumerate(iterator):
        values=block.iloc[:,positions].to_numpy(dtype=np.float64)
        assert np.isfinite(values).all() and values.min()>=0
        original_genes.extend(block.index.astype(str));chunks.append(sparse.csr_matrix(values))
        if k%25==0:print(donor,'parsed_genes',len(original_genes),flush=True)
    counts=sparse.vstack(chunks,format='csr');del chunks;gc.collect()
    assert counts.shape==(len(original_genes),len(ds))
    total=np.asarray(counts.sum(axis=0)).ravel();detected=np.asarray((counts>0).sum(axis=0)).ravel()
    mt=np.asarray([g.startswith('MT-') for g in original_genes]);mt_total=np.asarray(counts[mt].sum(axis=0)).ravel();mt_fraction=mt_total/total
    ds['original_feature_total']=total;ds['original_detected_features']=detected;ds['computed_mito_fraction']=mt_fraction;ds['mt20']=mt_fraction<.2
    ds['provided_mito_value']=ds.Percent_mitochondrial
    ds.to_csv(O/(donor+'_cell_QC.tsv.gz'),sep='\t',index=False,compression='gzip')
    mapping=pd.DataFrame({'original_gene':original_genes})
    mapping[['gene','status']]=pd.DataFrame([resolve(g) for g in original_genes])
    mapping.to_csv(O/(donor+'_feature_mapping.tsv.gz'),sep='\t',index=False,compression='gzip')
    kept=mapping.gene.ne('');genes=sorted(mapping.loc[kept,'gene'].unique());gidx={g:i for i,g in enumerate(genes)}
    transform=sparse.csr_matrix((np.ones(kept.sum()),([gidx[g] for g in mapping.loc[kept,'gene']],np.flatnonzero(kept))),shape=(len(genes),len(original_genes)))
    mapped=(transform@counts).tocsr();del counts;gc.collect()
    pb=[];pbs=[];summaries=[];coverage=[]
    for variant in ['primary','mt20']:
        base=np.ones(len(ds),dtype=bool) if variant=='primary' else ds.mt20.to_numpy()
        for name,labels in compartments.items():
            mask=base&ds.Cluster_identity.isin(labels).to_numpy();n=int(mask.sum())
            sid=donor+'__'+variant+'__'+name
            pb.append(np.asarray(mapped[:,mask].sum(axis=1)).ravel())
            pbs.append({'sample_id':sid,'donor':donor,'group':'OA' if donor.startswith('N')else'Crypto','variant':variant,'compartment':name,'n_cells':n,'source_labels':'|'.join(labels)})
        for name in sorted(meta.Cluster_identity.unique()):
            mask=base&ds.Cluster_identity.eq(name).to_numpy();n=int(mask.sum())
            coverage.append({'donor':donor,'group':'OA' if donor.startswith('N')else'Crypto','variant':variant,'author_cell_type':name,'n_cells':n})
            if not n:continue
            selected=[gidx[g] for g in sorted(source_genes) if g in gidx]
            sub=mapped[selected][:,mask].multiply(10000/total[mask]).tocsr()
            positive=np.asarray((sub>0).sum(axis=1)).ravel()/n
            mean_normalized=np.asarray(sub.sum(axis=1)).ravel()/n
            sub.data=np.log1p(sub.data);mean_log=np.asarray(sub.sum(axis=1)).ravel()/n
            for j,idx in enumerate(selected):
                summaries.append({'donor':donor,'group':'OA' if donor.startswith('N')else'Crypto','variant':variant,'author_cell_type':name,'gene':genes[idx],'n_cells':n,
                                  'mean_log1p_CPTT':mean_log[j],'mean_CPTT':mean_normalized[j],'estimated_count_positive_fraction':positive[j],
                                  'independent_annotation_marker':genes[idx]in marker_genes})
    pd.DataFrame(np.column_stack(pb),index=genes,columns=[x['sample_id']for x in pbs]).to_csv(O/(donor+'_pseudobulk_counts.tsv.gz'),sep='\t',index_label='gene',compression='gzip')
    pd.DataFrame(pbs).to_csv(O/(donor+'_pseudobulk_samples.tsv'),sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(O/(donor+'_author_label_coverage.tsv'),sep='\t',index=False)
    summary=pd.DataFrame(summaries);summary.to_csv(O/(donor+'_source_summaries.tsv.gz'),sep='\t',index=False,compression='gzip')
    summary[summary.independent_annotation_marker].to_csv(O/(donor+'_independent_marker_review.tsv'),sep='\t',index=False)
    receipt={'donor':donor,'status':'COUNTS_AGGREGATED_BEFORE_ECM_INFERENCE','cells':len(ds),'mt20_cells':int(ds.mt20.sum()),'original_features':len(original_genes),'mapped_genes':len(genes),
             'counts_are_fractional_estimates':True,'rounding_applied':False,'min_detected_features':int(detected.min()),'mito_fraction_range':[float(mt_fraction.min()),float(mt_fraction.max())],
             'provided_mito_fraction_max_difference':float(np.max(np.abs(mt_fraction-ds.Percent_mitochondrial.to_numpy()))),
             'raw_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'published_cell_mapping_complete':True}
    (S/'00_admin'/(donor+'_external_prepare_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    print('COMPLETE',json.dumps(receipt),flush=True)
    del mapped;gc.collect()
