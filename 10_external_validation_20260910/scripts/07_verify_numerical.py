from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
from scipy import stats

S=Path(__file__).resolve().parents[1];checks=[]
def check(name,value):
    assert bool(value),name
    checks.append({'check':name,'status':'PASS'})
def bh(values,n=None):
    a=np.asarray(values);order=np.argsort(a);ranked=a[order];m=len(a) if n is None else n
    q=np.minimum.accumulate((ranked*m/np.arange(1,len(a)+1))[::-1])[::-1]
    result=np.empty(len(a));result[order]=np.minimum(q,1);return result
freeze=json.loads((S/'00_admin/analysis_freeze_receipt.json').read_text())
for r in freeze['files']:check('Frozen '+Path(r['file']).name,hashlib.sha256((S/Path(r['file'])).read_bytes()).hexdigest()==r['sha256'])
names=list(json.loads((S/'01_metadata/fixed_gene_sets.json').read_text())['primary'])
ef=pd.read_csv(S/'04_results/bulk/external_score_effects.tsv',sep='\t')
for acc,noa,oa in [('GSE9210',47,11),('GSE108886',8,3)]:
    sm=pd.read_csv(S/'01_metadata'/(acc+'_analysis_samples.tsv'),sep='\t')
    check(acc+' unique patient groups',sm.gsm.is_unique and (sm.group=='NOA').sum()==noa and (sm.group=='OA').sum()==oa)
    check(acc+' excludes pooled control','GSM2915444'not in set(sm.gsm))
    scores=pd.read_csv(S/'04_results/bulk'/(acc+'_scores.tsv'),sep='\t',index_col=0)
    tab=ef[ef.cohort==acc]
    for r in tab.itertuples():
        a=scores.loc[r.program,sm.loc[sm.group=='NOA','gsm']].to_numpy();b=scores.loc[r.program,sm.loc[sm.group=='OA','gsm']].to_numpy()
        est=a.mean()-b.mean();df=len(a)+len(b)-2
        variance=((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/df
        se=np.sqrt(variance*(1/len(a)+1/len(b)));p=2*stats.t.sf(abs(est/se),df)
        check(acc+' independent score calculation '+r.program,np.allclose([est,se,p,est-stats.t.ppf(.975,df)*se,est+stats.t.ppf(.975,df)*se],[r.effect,r.SE,r.P,r.CI_low,r.CI_high],rtol=1e-9,atol=1e-11))
    primary=tab[tab.program.isin(names)];check(acc+' BH six',np.allclose(bh(primary.P,6),primary.FDR_six))
    genes=pd.read_csv(S/'04_results/bulk'/(acc+'_all_genes.tsv.gz'),sep='\t')
    check(acc+' all gene intervals and FDR',genes.gene.is_unique and (genes['CI.L']<=genes.expression_difference).all() and (genes['CI.R']>=genes.expression_difference).all() and genes['adj.P.Val'].between(0,1).all())
primary=ef[ef.program.isin(names)];check('External 12-test score family',np.allclose(bh(primary.P,12),primary.FDR_external_12))
audit=pd.read_csv(S/'01_metadata/GSE45885_GSE45887_reuse.tsv',sep='\t');check('Twenty native sample reuse links',len(audit)==20 and audit.reanalysis_gsm.is_unique)
pc=pd.read_csv(S/'03_processed/single_cell/all_pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
sm=pd.read_csv(S/'03_processed/single_cell/all_pseudobulk_samples.tsv',sep='\t')
check('Pseudobulk column identity',pc.columns.tolist()==sm.sample_id.tolist())
check('Fractional count preservation',(pc.to_numpy()>=0).all() and ((pc.to_numpy()%1)>0).any())
for d in ['N1','N2','N3','Cr1','Cr2','Cr3']:
    for v in ['primary','mt20']:
        a=pc[d+'__'+v+'__Peritubular_compartment'];b=pc[d+'__'+v+'__Author_PMC_state']+pc[d+'__'+v+'__Author_fibrotic_PMC_state']
        check(d+' '+v+' count conservation',np.allclose(a,b,rtol=1e-12,atol=1e-7))
    for c in sm.compartment.unique():check(d+' '+c+' mitochondrial subset count bound',(pc[d+'__mt20__'+c]<=pc[d+'__primary__'+c]+1e-7).all())
    r=json.loads((S/'00_admin'/(d+'_external_prepare_receipt.json')).read_text())
    source=next((S/'02_raw').glob('*_'+d+'_counts.tsv.gz'))
    check(d+' downloaded source unchanged',r['raw_sha256']==hashlib.sha256(source.read_bytes()).hexdigest())
elig=pd.read_csv(S/'04_results/single_cell/donor_compartment_eligibility.tsv',sep='\t')
check('Eligible compartments are prespecified',set(elig.loc[elig.eligible,'compartment'])=={'Peritubular_compartment','Author_PMC_state'})
ca=pd.read_csv(S/'04_results/single_cell/program_camera_and_effects.tsv',sep='\t')
for (fam,mode),g in ca.groupby(['family','mode']):check('SC BH '+fam+' '+mode,np.allclose(bh(g.PValue,len(g.compartment.unique())*6),g.FDR_family))
genes=pd.read_csv(S/'04_results/single_cell/all_30_laminin_gene_effects.tsv',sep='\t')
check('Every SC scenario reports all 30 members',genes.groupby('scenario').gene.nunique().eq(30).all())
check('SC intervals contain all measured estimates',(genes.loc[genes.measured,'CI.L']<=genes.loc[genes.measured,'logFC']).all() and (genes.loc[genes.measured,'CI.R']>=genes.loc[genes.measured,'logFC']).all())
coverage=pd.read_csv(S/'04_results/external_cell_coverage.tsv',sep='\t')
check('Published cells counted once',coverage[coverage.variant=='primary'].n_cells.sum()==28666)
(S/'00_admin/external_numerical_QA.json').write_text(json.dumps({'status':'PASS','checks_passed':len(checks),'checks':checks},indent=2)+'\n')
print('PASS',len(checks),'scientific checks')
