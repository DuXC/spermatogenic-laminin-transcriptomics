#!/usr/bin/env python3
"""Apply reviewed cluster identities; export donor-level sources and counts."""
from pathlib import Path
import json, itertools
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

S=Path(__file__).resolve().parents[1];P=S.parent
CACHE=(P / '.cache' / 'single_cell')

def main():
 a=ad.read_h5ad(CACHE/'adult_preannotation.h5ad')
 labels=pd.read_csv(S/'02_annotation/cluster_labels.tsv',sep='\t',dtype={'cluster':str})
 assert labels.cluster.is_unique and set(labels.cluster)==set(a.obs.leiden.astype(str))
 labels=labels.set_index('cluster')
 a.obs['cell_type']=a.obs.leiden.astype(str).map(labels.cell_type).to_numpy()
 a.obs['annotation_confidence']=a.obs.leiden.astype(str).map(labels.confidence).to_numpy()
 # Optional reviewed cell overrides from targeted mural subclustering.
 overrides=S/'02_annotation/cell_annotation_overrides.tsv.gz'
 if overrides.exists():
  o=pd.read_csv(overrides,sep='\t',index_col=0)
  assert o.index.is_unique and o.index.isin(a.obs_names).all()
  a.obs.loc[o.index,['cell_type','annotation_confidence']]=o[['cell_type','annotation_confidence']]
 a.obs.to_csv(S/'02_annotation/cells_final.tsv.gz',sep='\t')
 a.write_h5ad(CACHE/'adult_annotated.h5ad',compression='gzip')
 # Final counts, normalized expression, embeddings and cell labels share one object.
 a.write_h5ad(S/'03_processed/adult_annotated.h5ad',compression='gzip')
 sets=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary']
 ecm=sorted(set().union(*map(set,sets.values())))
 present=[x for x in ecm if x in a.var_names]
 ix=a.var_names.get_indexer(present)
 logX=a.X[:,ix].tocsr();norm=logX.copy();norm.data=np.expm1(norm.data)
 raw=a.layers['counts']
 samples=[];vectors=[];sources=[];countrows=[]
 types=sorted(a.obs.cell_type.unique());donors=sorted(a.obs.donor_id.unique())
 variants={'primary':np.ones(a.n_obs,dtype=bool),
           'mt20':a.obs.mt20_pass.to_numpy(dtype=bool),
           'no_predicted_doublets':~a.obs.predicted_doublet.to_numpy(dtype=bool)}
 for variant,keep in variants.items():
  for donor,typ in itertools.product(donors,types):
   m=keep & a.obs.donor_id.eq(donor).to_numpy() & a.obs.cell_type.eq(typ).to_numpy()
   n=int(m.sum());r=a.obs[a.obs.donor_id.eq(donor)].iloc[0]
   ann='high' if n and a.obs.loc[m,'annotation_confidence'].eq('high').all() else 'provisional'
   unit=dict(variant=variant,donor_id=donor,cell_type=typ,group=str(r['group']),technology=str(r.technology),n_cells=n,annotation_confidence=ann)
   countrows.append(unit)
   if not n:continue
   sample_id=f'{variant}__{donor}__{typ}'
   samples.append(dict(sample_id=sample_id,**unit))
   vectors.append(np.asarray(raw[m].sum(axis=0)).ravel())
   mean=np.asarray(norm[m].mean(axis=0)).ravel()
   detect=np.asarray((raw[m][:,ix]>0).sum(axis=0)).ravel()/n
   d=pd.DataFrame({'gene':present,'mean_CP10K':mean,'fraction_detected':detect})
   for k,v in unit.items():d[k]=v
   d['source_eligible']=(n>=30 and ann=='high')
   sources.append(d)
 counts=pd.DataFrame(countrows)
 counts['capture_fraction']=counts.n_cells/counts.groupby(['variant','donor_id']).n_cells.transform('sum')
 counts.to_csv(S/'04_results/donor_celltype_coverage.tsv',sep='\t',index=False)
 meta=pd.DataFrame(samples)
 meta.to_csv(S/'03_processed/pseudobulk_samples.tsv',sep='\t',index=False)
 pd.DataFrame(np.column_stack(vectors),index=pd.Index(a.var_names,name='gene'),columns=meta.sample_id).to_csv(S/'03_processed/pseudobulk_counts.tsv.gz',sep='\t')
 src=pd.concat(sources,ignore_index=True)
 src.to_csv(S/'04_results/ECM_all_genes_donor_sources.tsv.gz',sep='\t',index=False)
 gene_cov=[]
 all_detection=np.asarray((raw>0).sum(axis=0)).ravel()
 reference_detection=np.asarray((raw[a.obs.group.eq('OA').to_numpy()]>0).sum(axis=0)).ravel()
 for name,genes in sets.items():
  for gene in genes:
   j=a.var_names.get_indexer([gene])[0]
   gene_cov.append(dict(pathway=name,gene=gene,measured=j>=0,
                        reference_cells_detected=int(reference_detection[j]) if j>=0 else 0,
                        all_cells_detected=int(all_detection[j]) if j>=0 else 0))
 pd.DataFrame(gene_cov).to_csv(S/'04_results/ECM_program_gene_coverage.tsv',sep='\t',index=False)
 # One gene standardization learned from primary eligible OA donor/type units.
 ref=src[(src.variant=='primary') & (src.group=='OA') & src.source_eligible].copy()
 ref['log1p_mean']=np.log1p(ref.mean_CP10K)
 param=ref.groupby('gene').log1p_mean.agg(['mean','std','count'])
 param.to_csv(S/'04_results/source_score_standardization.tsv',sep='\t')
 src['z']=(np.log1p(src.mean_CP10K)-src.gene.map(param['mean']))/src.gene.map(param['std']).replace(0,np.nan)
 scores=[]
 keys=['variant','donor_id','cell_type','group','technology','n_cells','annotation_confidence','source_eligible']
 for name,genes in sets.items():
  sub=src[src.gene.isin(genes)]
  d=sub.groupby(keys,observed=True).z.agg(['mean','count']).reset_index().rename(columns={'mean':'score','count':'scored_genes'})
  d['pathway']=name;d['fixed_set_size']=len(genes)
  scores.append(d)
 scores=pd.concat(scores,ignore_index=True)
 scores.to_csv(S/'04_results/ECM_program_donor_source_scores.tsv',sep='\t',index=False)
 ref_scores=scores[(scores.variant=='primary') & (scores.group=='OA') & scores.source_eligible]
 summary=ref_scores.groupby(['cell_type','pathway']).score.agg(['mean','std','count']).reset_index()
 summary['SE']=summary['std']/np.sqrt(summary['count'])
 summary['CI_low']=summary['mean']-stats.t.ppf(.975,summary['count']-1)*summary.SE
 summary['CI_high']=summary['mean']+stats.t.ppf(.975,summary['count']-1)*summary.SE
 summary.to_csv(S/'04_results/reference_program_source_summary.tsv',sep='\t',index=False)
 genes=ref.groupby(['cell_type','gene']).agg(mean_CP10K=('mean_CP10K','mean'),fraction_detected=('fraction_detected','mean'),n_donors=('donor_id','nunique')).reset_index()
 genes.to_csv(S/'04_results/reference_ECM_gene_source_summary.tsv',sep='\t',index=False)
 # All-donor independent marker summaries for visual label checking.
 panels=json.loads((S/'02_annotation/annotation_marker_panels.json').read_text())
 mark=sorted(set().union(*map(set,panels.values())).intersection(a.var_names))
 midx=a.var_names.get_indexer(mark);markerrows=[]
 for donor,typ in itertools.product(donors,types):
  m=a.obs.donor_id.eq(donor).to_numpy() & a.obs.cell_type.eq(typ).to_numpy()
  if m.sum()<30:continue
  val=np.asarray(a.X[m][:,midx].mean(axis=0)).ravel()
  det=np.asarray((raw[m][:,midx]>0).sum(axis=0)).ravel()/int(m.sum())
  t=pd.DataFrame(dict(gene=mark,mean_log1p_CP10K=val,fraction_detected=det,donor_id=donor,cell_type=typ,n_cells=int(m.sum())))
  markerrows.append(t)
 pd.concat(markerrows).to_csv(S/'04_results/annotation_marker_donor_means.tsv.gz',sep='\t',index=False)
 print('COMPLETE',a.shape,'pseudobulk units',len(meta),'source rows',len(src),flush=True)

if __name__=='__main__':main()
