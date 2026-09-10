from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def bh(p):
    p=np.array(p,dtype=float);order=np.argsort(p);q=np.empty_like(p)
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q

N=Path(__file__).resolve().parents[1]
O=N/'03_analysis/PXD011817'
F=N.parent/'13_public_data_submission_20260910/03_figures';F.mkdir(parents=True,exist_ok=True)
score=pd.read_csv(O/'program_score_group_effects.tsv',sep='\t').set_index('program')
camera=pd.read_csv(O/'camera_estimated_correlation.tsv',sep='\t').set_index('program')
people=pd.read_csv(O/'participant_program_scores.tsv',sep='\t')
gene=pd.read_csv(O/'all_gene_results.tsv',sep='\t').set_index('gene')
x=pd.read_csv(O/'gene_log2_LFQ_min3_each_group.tsv',sep='\t',index_col=0)
fixed=json.loads((N.parent/'10_external_validation_20260910/01_metadata/fixed_gene_sets.json').read_text())['primary']
lam='REACTOME_LAMININ_INTERACTIONS'
audit=[]
for name in score.index:
 a=people.loc[people.group.eq('SRpos'),name].to_numpy();b=people.loc[people.group.eq('SRneg'),name].to_numpy()
 delta=b.mean()-a.mean()
 p=stats.ttest_ind(b,a,equal_var=True).pvalue
 se=np.sqrt(((len(a)-1)*np.var(a,ddof=1)+(len(b)-1)*np.var(b,ddof=1))/(len(a)+len(b)-2)*(1/len(a)+1/len(b)))
 low,high=delta+np.array([-1,1])*stats.t.ppf(.975,8)*se
 assert np.allclose([delta,p,low,high],score.loc[name,['effect','PValue','CI_low','CI_high']].to_numpy(float),rtol=1e-8,atol=1e-10)
 audit.append(dict(check='program_effect_CI_p_independent_scipy',program=name,pass_check=True))
delta=x.iloc[:,5:].mean(axis=1)-x.iloc[:,:5].mean(axis=1)
assert np.allclose(delta.loc[gene.index],gene.logFC)
assert np.allclose(bh(gene['P.Value']),gene['adj.P.Val'])
assert np.allclose(bh(camera.PValue),camera.FDR_six)
complete=pd.read_csv(O/'gene_log2_LFQ_complete_10.tsv',sep='\t',index_col=0)
measured=[g for g in sorted(fixed[lam]) if g in gene.index]
all30=pd.DataFrame(index=sorted(fixed[lam]));all30.index.name='gene'
all30=all30.join(gene[['logFC','CI.L','CI.R','P.Value','adj.P.Val','SRpos_n','SRneg_n']])
all30['status']=np.where(all30.index.isin(gene.index),'EVALUATED','NOT_EVALUABLE_AFTER_QC')
all30.to_csv(O/'all_30_laminin_members_results.tsv',sep='\t',na_rep='NA')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
fig=plt.figure(figsize=(183/25.4,150/25.4))
gs=fig.add_gridspec(2,2,width_ratios=[1,1.12],height_ratios=[.8,1.15],hspace=.53,wspace=.72,left=.22,right=.98,bottom=.14,top=.91)
ax=fig.add_subplot(gs[0,:])
names=list(score.index)
labels=['Collagen formation','ECM degradation','ECM organization','Integrin interactions','Laminin interactions','Basement membrane']
y=np.arange(len(names))[::-1]
v=score.loc[names]
ax.errorbar(v.effect,y,xerr=[v.effect-v.CI_low,v.CI_high-v.effect],fmt='o',color='#385F77',ecolor='#8C9BA5',capsize=2,markersize=3,linewidth=.8)
ax.axvline(0,color='#555555',linewidth=.6,linestyle='--')
ax.set_yticks(y,labels);ax.set_xlim(-1.2,1.2);ax.set_xlabel('Program-score difference (95% CI)')
ax.set_title('a   Six fixed ECM programs',loc='left',pad=8)
ax.text(1.0,1.08,'5 retrieval-negative vs 5 retrieval-positive participants',ha='right',transform=ax.transAxes,fontsize=7)
ax.text(.02,-.35,'Estimated-correlation CAMERA: all six q = 0.930',transform=ax.transAxes,fontsize=7)

bx=fig.add_subplot(gs[1,0])
colors=['#718DA2','#BE9C7D']
for i,g in enumerate(['SRpos','SRneg']):
 vv=people.loc[people.group.eq(g),lam].to_numpy()
 bx.scatter(i+np.linspace(-.08,.08,len(vv)),vv,s=17,color=colors[i],edgecolor='white',linewidth=.3,zorder=3)
 bx.plot([i-.2,i+.2],[vv.mean(),vv.mean()],color='#333333',linewidth=1.1)
bx.set_xticks([0,1],['Retrieval +\nn = 5','Retrieval −\nn = 5']);bx.set_xlim(-.55,1.55)
bx.set_ylabel('Mean standardized laminin score')
bx.set_title('b   Participant scores',loc='left',pad=8)
bx.text(.02,1.00,'15 measured members',va='top',transform=bx.transAxes,fontsize=7)

cx=fig.add_subplot(gs[1,1])
gg=gene.loc[measured];yy=np.arange(len(gg))[::-1]
cx.errorbar(gg.logFC,yy,xerr=[gg.logFC-gg['CI.L'],gg['CI.R']-gg.logFC],fmt='o',color='#385F77',ecolor='#8C9BA5',capsize=1.5,markersize=2.8,linewidth=.7)
cx.set_yticks(yy,gg.index);cx.axvline(0,color='#555555',linestyle='--',linewidth=.6);cx.set_xlim(-2.6,2.6)
cx.set_xlabel('log2 LFQ difference (95% CI)')
cx.set_title('c   Measured laminin members',loc='left',pad=8)
fig.text(.015,.975,'Independent ECM proteomics in idiopathic NOA',fontsize=10,va='top')
fig.text(.015,.045,'Positive differences: retrieval-negative minus retrieval-positive. Technical triplicates aggregated per participant.\nAll 15 gene q values > 0.79 after correction across 412 genes; intervals are unadjusted. Full 30-member table supplied.',fontsize=7,va='bottom')
base=F/'Figure_6'
for ext in ['svg','pdf','png','tiff']:
 fig.savefig(base.with_suffix('.'+ext),dpi=600 if ext=='tiff' else 200,facecolor='white')
plt.close(fig)
receipt=dict(status='NUMERICAL_QA_PASS',independent_checks=audit,gene_mean_difference_check=True,gene_BH_check=True,camera_BH_check=True,n_participants=10,n_gene_tests=len(gene),all_30_members_exported=len(all30),program='REACTOME_LAMININ_INTERACTIONS',camera_q=float(camera.loc[lam,'FDR_six']),program_effect=float(score.loc[lam,'effect']),program_CI=[float(score.loc[lam,'CI_low']),float(score.loc[lam,'CI_high'])],figure_visual_QA='PENDING',source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.glob('*.tsv')})
(N.parent/'13_public_data_submission_20260910/00_admin/protein_numerical_QA.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','independent_checks']},indent=2))
