#!/usr/bin/env python3
"""Recount detection as integer positives / cells, avoiding sparse-mean drift."""
from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
S=Path(__file__).resolve().parents[1]
a=ad.read_h5ad(S/'03_processed/adult_annotated.h5ad')
src=pd.read_csv(S/'04_results/ECM_all_genes_donor_sources.tsv.gz',sep='\t')
before=src.fraction_detected.to_numpy().copy()
genes=sorted(src.gene.unique());ix=a.var_names.get_indexer(genes);x=(a.layers['counts'][:,ix]>0).tocsr()
variants={'primary':np.ones(a.n_obs,dtype=bool),'mt20':a.obs.mt20_pass.to_numpy(dtype=bool),'no_predicted_doublets':~a.obs.predicted_doublet.to_numpy(dtype=bool)}
for (variant,donor,typ),loc in src.groupby(['variant','donor_id','cell_type']).indices.items():
 mask=variants[variant]&a.obs.donor_id.eq(donor).to_numpy()&a.obs.cell_type.eq(typ).to_numpy()
 frac=pd.Series(np.asarray(x[mask].sum(axis=0)).ravel()/int(mask.sum()),index=genes)
 src.loc[loc,'fraction_detected']=src.loc[loc,'gene'].map(frac).to_numpy()
src.to_csv(S/'04_results/ECM_all_genes_donor_sources.tsv.gz',sep='\t',index=False)
ref=src[src.variant.eq('primary')&src.group.eq('OA')&src.source_eligible]
ref.groupby(['cell_type','gene']).agg(mean_CP10K=('mean_CP10K','mean'),fraction_detected=('fraction_detected','mean'),n_donors=('donor_id','nunique')).reset_index().to_csv(S/'04_results/reference_ECM_gene_source_summary.tsv',sep='\t',index=False)
for fn,keys in [('annotation_marker_donor_means.tsv.gz',['donor_id','cell_type']),('cluster_independent_marker_expression.tsv',['cluster','group'])]:
 p=S/'04_results'/fn;t=pd.read_csv(p,sep='\t');gg=sorted(t.gene.unique());xx=(a.layers['counts'][:,a.var_names.get_indexer(gg)]>0).tocsr()
 for unit,loc in t.groupby(keys).indices.items():
  mask=np.ones(a.n_obs,dtype=bool)
  for key,val in zip(keys,unit):
   k='leiden' if key=='cluster' else key
   mask &= a.obs[k].astype(str).eq(str(val)).to_numpy()
  f=pd.Series(np.asarray(xx[mask].sum(axis=0)).ravel()/int(mask.sum()),index=gg)
  t.loc[loc,'fraction_detected']=t.loc[loc,'gene'].map(f).to_numpy()
 t.to_csv(p,sep='\t',index=False)
receipt={'status':'EXACT_INTEGER_RECOUNT_COMPLETE','original_bound_overshoots':int((before>1).sum()),'maximum_absolute_fraction_change':float(np.max(np.abs(before-src.fraction_detected.to_numpy()))),'all_fractions_in_unit_interval':bool(src.fraction_detected.between(0,1).all()),'counts_labels_scores_and_pseudobulk_effects_changed':False}
(S/'00_admin/detection_fraction_QA.json').write_text(json.dumps(receipt,indent=2))
print(receipt,flush=True)
