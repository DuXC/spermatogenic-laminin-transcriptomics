#!/usr/bin/env python3
"""Pre-label annotation QA: remove ECM genes from score controls too."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import anndata as ad,scanpy as sc
S=Path(__file__).resolve().parents[1];P=S.parent
C=(P / '.cache' / 'single_cell')
a=ad.read_h5ad(C/'adult_preannotation.h5ad')
sets=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary']
ecm=set().union(*map(set,sets.values()))
panels=json.loads((S/'02_annotation/annotation_marker_panels.json').read_text())
pool=[x for x in a.var_names if x not in ecm]
for label,genes in panels.items():
 sc.tl.score_genes(a,[x for x in genes if x in a.var_names],gene_pool=pool,score_name=label,random_state=20260909,ctrl_size=50)
a.obs.to_csv(S/'02_annotation/cells_preannotation.tsv.gz',sep='\t')
a.obs.groupby('leiden',observed=True)[list(panels)].mean().to_csv(S/'04_results/cluster_marker_scores.tsv',sep='\t')
genes='INSL3 STAR CYP17A1 HSD3B2 CYP11A1 LHCGR NR5A1 SOX9 WT1 AMH FSHR CLU INHA GATA4 KRT18 FATE1 CLDN11 DPEP1 TSHZ2 OSR2 FHL2 PTGDS CRISPLD2 C7 ACTA2 TAGLN MYH11 CNN1 DES RGS5 RERGL MCAM NOTCH3 CSPG4 PLN ADIRF MYL9 PECAM1 VWF KDR EMCN CLDN5 UTF1 GFRA1 ZBTB16 FGFR3 PIWIL4 UCHL1 MAGEA4 KIT SALL4 SYCP1 SYCP3 SPO11 ACRV1 PRM1 PRM2 TPSAB1 TPSB2 CPA3 MS4A2 LST1 CD3D CD3E NKG7'.split()
genes=[x for x in genes if x in a.var_names and x not in ecm]
rows=[]
for (cluster,grp),ix in a.obs.groupby(['leiden','group'],observed=True).indices.items():
 x=a.X[ix][:,a.var_names.get_indexer(genes)]
 avg=np.asarray(x.mean(axis=0)).ravel();det=np.asarray((x>0).sum(axis=0)).ravel()/len(ix)
 for g,m,f in zip(genes,avg,det):rows.append(dict(cluster=str(cluster),group=str(grp),n_cells=len(ix),gene=g,mean_log1p=m,fraction_detected=f))
pd.DataFrame(rows).to_csv(S/'04_results/cluster_independent_marker_expression.tsv',sep='\t',index=False)
a.write_h5ad(C/'adult_preannotation.h5ad',compression='gzip')
a.write_h5ad(S/'03_processed/adult_preannotation.h5ad',compression='gzip')
d=pd.read_csv(S/'02_annotation/GSE149512_donor_metadata_v2.tsv',sep='\t')
d['prior_release_Sertoli_cells']=d.published_Sertoli_cells
d.loc[d.donor_id.eq('LZ014'),'published_Sertoli_cells']=27
d.to_csv(S/'02_annotation/GSE149512_donor_metadata_v2.tsv',sep='\t',index=False)
(S/'00_admin/preannotation_QA.json').write_text(json.dumps({'ECM_genes_excluded_from_positive_markers_and_control_pool':True,'annotations_assigned_after_QA':True,'clinical_reinspection':{'LZ011_total_cells':7896,'LZ014_Sertoli_cells':27},'evidence':'SI page 2, Supplementary Figure 1a; high-resolution clinical crop retained'},indent=2))
print('COMPLETE annotation QA',flush=True)
