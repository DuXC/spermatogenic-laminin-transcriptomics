#!/usr/bin/env python3
"""Pooled-cell identity display only; gene selection excludes the ECM union."""
from pathlib import Path
import json
import anndata as ad
import numpy as np,pandas as pd
S=Path(__file__).resolve().parents[1];P=S.parent
a=ad.read_h5ad(S/'03_processed/adult_annotated.h5ad')
ecm=set().union(*map(set,json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary'].values()))
genes='UTF1 MAGEA4 SYCP3 HORMAD1 PRM1 ACRV1 SOX9 FATE1 CLDN11 INSL3 STAR CYP17A1 DPEP1 MYH11 RGS5 MCAM NOTCH3 CLDN5 KDR GNG11 TYROBP FCER1G LST1 CPA3 MS4A2 CD3D NKG7'.split()
genes=[x for x in genes if x in a.var_names and x not in ecm]
x=a.X[:,a.var_names.get_indexer(genes)].tocsr();rows=[]
for typ,ix in a.obs.groupby('cell_type',observed=True).indices.items():
 avg=np.asarray(x[ix].mean(axis=0)).ravel();frac=np.asarray((x[ix]>0).sum(axis=0)).ravel()/len(ix)
 for gene,m,f in zip(genes,avg,frac):rows.append(dict(cell_type=typ,gene=gene,mean=m,fraction=f,n_cells=len(ix)))
pd.DataFrame(rows).to_csv(S/'04_results/annotation_display_means.tsv',sep='\t',index=False)
print('COMPLETE',len(genes),'independent displayed markers',flush=True)
