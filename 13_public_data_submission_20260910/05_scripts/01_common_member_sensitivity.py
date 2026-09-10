"""Frozen common-coverage sensitivity; no new gene or cohort selection."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy.stats import t

Q=Path(__file__).resolve().parents[1];P=Q.parent;S=P/'10_external_validation_20260910'
O=Q/'02_sources/common_member';O.mkdir(exist_ok=True)
cohorts=['GSE4797','GSE145467','GSE9210','GSE108886']
inputs={}
mat={};meta={}
for c in cohorts:
 f=P/'03_processed'/f'{c}_gene_log_expression.tsv.gz' if c in cohorts[:2] else S/'03_processed'/f'{c}_primary_gene_expression.tsv.gz'
 m=P/'02_annotation'/f'{c}_samples.tsv' if c in cohorts[:2] else S/'01_metadata'/f'{c}_analysis_samples.tsv'
 mat[c]=pd.read_csv(f,sep='\t',index_col=0)
 meta[c]=pd.read_csv(m,sep='\t').set_index('gsm').loc[mat[c].columns]
 assert not mat[c].index.duplicated().any()
 assert np.isfinite(mat[c].values).all()
 for x in [f,m]:inputs[str(x.relative_to(P))]=hashlib.sha256(x.read_bytes()).hexdigest()
membership=P/'02_annotation/fixed_pathway_membership.tsv'
sets=pd.read_csv(membership,sep='\t')
inputs[str(membership.relative_to(P))]=hashlib.sha256(membership.read_bytes()).hexdigest()
shared=set.intersection(*(set(x.index) for x in mat.values()))
def estimate(c,y):
 g=meta[c]['group'].astype(str)
 levels=['JS10','JS8','JS5','JS2'] if c=='GSE4797' else ['normal','impaired'] if c=='GSE145467' else ['OA','NOA']
 assert set(g)==set(levels),(c,set(g))
 X=np.column_stack([(g==x).astype(float) for x in levels]);v=np.zeros(len(levels));v[0]=-1;v[-1]=1
 xtx=np.linalg.inv(X.T@X);b=xtx@X.T@y;res=y-X@b;df=len(y)-len(levels)
 est=float(v@b);se=float(np.sqrt((res@res/df)*(v@xtx@v)))
 return {'effect':est,'CI_low':est-t.ppf(.975,df)*se,'CI_high':est+t.ppf(.975,df)*se,'SE':se,'df':df,'p':2*t.sf(abs(est/se),df),'n':len(y)}
rows=[];scores=[];members=[];loo=[]
for name,g in sets.groupby('set',sort=True):
 allgenes=set(g.current_gene.dropna());common=sorted(allgenes&shared)
 assert len(common)>=2,(name,len(common))
 for gene in sorted(allgenes):members.append({'program':name,'gene':gene,'common_to_four':gene in common})
 for c in cohorts:
  e=mat[c].loc[common].T
  assert (e.std(ddof=1)>0).all()
  z=(e-e.mean())/e.std(ddof=1);y=z.mean(axis=1)
  rows.append({'cohort':c,'program':name,'n_common':len(common),'n_fixed':len(allgenes),**estimate(c,y.values)})
  for gsm,val in y.items():scores.append({'cohort':c,'program':name,'gsm':gsm,'group':meta[c].loc[gsm,'group'],'score':val})
  if name=='REACTOME_LAMININ_INTERACTIONS':
   for gene in common:loo.append({'cohort':c,'omitted_gene':gene,'n_retained':len(common)-1,**estimate(c,z.drop(columns=gene).mean(axis=1).values)})
res=pd.DataFrame(rows);pvals=res.p.values;order=np.argsort(pvals);q=np.minimum.accumulate((pvals[order]*len(pvals)/np.arange(1,len(pvals)+1))[::-1])[::-1];out=np.empty(len(q));out[order]=np.minimum(q,1);res['BH_FDR_24_tests']=out
for name,frame in [('all_program_contrasts',res),('common_member_audit',pd.DataFrame(members)),('sample_scores',pd.DataFrame(scores)),('laminin_leave_one_gene_out',pd.DataFrame(loo))]:frame.to_csv(O/f'{name}.tsv',sep='\t',index=False)
receipt={'status':'COMPLETE','plan_sha256':hashlib.sha256((Q/'00_admin/COMMON_MEMBER_SENSITIVITY_PLAN.md').read_bytes()).hexdigest(),'inputs':inputs,'common_background_genes':len(shared),'post_hoc':True,'tests':len(res),'cohort_effects_not_pooled':True,'laminin':res[res.program=='REACTOME_LAMININ_INTERACTIONS'].to_dict('records')}
(O/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt['laminin'],indent=2))
