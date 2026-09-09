"""Authoritative publication figure generator; reads frozen bulk-analysis TSVs."""
from pathlib import Path
import json, textwrap, importlib.metadata
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import t

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'05_figures';RES=ROOT/'04_results'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
rng=np.random.default_rng(20260909)
COL={'JS10':'#617E98','JS8':'#A6B6C5','JS5':'#D3AE7A','JS2':'#A35B55','normal':'#617E98','impaired':'#A35B55'}
COHORT={'GSE4797':'#52738B','GSE145467':'#AE6B55'}
P=['GOCC_BASEMENT_MEMBRANE','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_LAMININ_INTERACTIONS','REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
LABEL=dict(zip(P,['Basement membrane','ECM organization','Collagen formation','Laminin interactions','Integrin interactions','ECM degradation']))
ACCS=['GSE4797','GSE145467'];sources={};figures=[]
def read(rel):return pd.read_csv(ROOT/rel,sep='\t')
S={a:read(f'02_annotation/{a}_samples.tsv') for a in ACCS}
SC={a:read(f'04_results/primary/{a}_module_scores.tsv').set_index('set') for a in ACCS}
E={a:read(f'04_results/primary/{a}_module_effects.tsv') for a in ACCS}
CA={a:read(f'04_results/primary/{a}_camera.tsv') for a in ACCS}
def letter(ax,x):ax.text(-.14,1.04,x,transform=ax.transAxes,weight='bold',fontsize=10,va='bottom')
def save(fig,name,src):
    for ext in ['svg','pdf','png']:fig.savefig(OUT/f'{name}.{ext}',dpi=300)
    figures.append((name,fig));sources[name]=src
def labels(ax,groups):
    ax.set_xticks(range(len(groups)),[g.replace('JS','JS ') for g in groups]);ax.set_xlim(-.5,len(groups)-.5)

# 1. Analyzed units and unsupervised sample QC.
fig,axs=plt.subplots(2,2,figsize=(183/25.4,135/25.4),layout='constrained')
ax=axs[0,0];counts=[12,6,5,5,10,10];ys=[5,4,3,2,.5,-.5];names=['JS 10','JS 8','JS 5','JS 2','Normal','Impaired'];colors=[COL[k] for k in ['JS10','JS8','JS5','JS2','normal','impaired']]
ax.barh(ys,counts,color=colors,height=.65);ax.set_yticks(ys,names);ax.set_xlim(0,15);ax.set_xlabel('Men / study samples');ax.set_title('GSE4797: 28 men; GSE145467: 20 samples',loc='left')
for x,y in zip(counts,ys):ax.text(x+.2,y,str(x),va='center')
letter(ax,'a')
for a,ax,let in [('GSE4797',axs[0,1],'b'),('GSE145467',axs[1,0],'c')]:
    d=read(f'04_results/qc/{a}_PCA.tsv');v=read(f'04_results/qc/{a}_PCA_variance.tsv').variance_fraction
    for g,sub in d.groupby('group',sort=False):ax.scatter(sub.PC1,sub.PC2,s=23,color=COL[g],label=f'{g} (n={len(sub)})',edgecolor='white',linewidth=.35)
    ax.set_xlabel(f'PC1 ({100*v.iloc[0]:.1f}%)');ax.set_ylabel(f'PC2 ({100*v.iloc[1]:.1f}%)');ax.set_title(a,loc='left');ax.legend(loc='best',frameon=False,ncol=2);letter(ax,let)
ax=axs[1,1]
for i,a in enumerate(ACCS):
    d=read(f'04_results/qc/{a}_sample_QC.tsv');ax.scatter(i+rng.uniform(-.16,.16,len(d)),d.within_group_median_correlation,s=22,c=d.group.map(COL),edgecolor='white',linewidth=.35)
    ax.plot([i-.21,i+.21],[d.within_group_median_correlation.median()]*2,color='black',lw=1)
ax.set_xticks([0,1],ACCS);ax.set_ylim(.4,1);ax.set_ylabel('Median correlation with same-group samples');ax.set_title('Within-group expression similarity',loc='left');letter(ax,'d')
save(fig,'Figure_1_Cohorts_QC',['02_annotation/*_samples.tsv','04_results/qc/*'])

# 2. Six fixed pathway scores, every patient shown.
fig,axs=plt.subplots(2,3,figsize=(183/25.4,125/25.4),layout='constrained')
s=S['GSE4797'];groups=['JS10','JS8','JS5','JS2']
for ax,nm,let in zip(axs.flat,P,'abcdef'):
    means=[];cis=[]
    for i,g in enumerate(groups):
        vals=SC['GSE4797'].loc[nm,s.loc[s.group==g,'gsm']].to_numpy(float)
        ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=16,color=COL[g],alpha=.85,edgecolor='white',linewidth=.3)
        means.append(vals.mean());cis.append(t.ppf(.975,len(vals)-1)*vals.std(ddof=1)/np.sqrt(len(vals)))
    ax.errorbar(range(4),means,yerr=cis,fmt='o-',color='#303B43',markersize=3,lw=.8,capsize=2)
    labels(ax,groups);ax.set_title(LABEL[nm],loc='left');ax.axhline(0,color='#D4D4D4',lw=.6,zorder=0);ax.set_ylabel('Mean gene z score');letter(ax,let)
fig.supxlabel('Cross-sectional histological categories: n = 12 / 6 / 5 / 5',fontsize=7)
save(fig,'Figure_2_Pathology_ECM',['04_results/primary/GSE4797_module_scores.tsv','02_annotation/GSE4797_samples.tsv'])

# 3. Score effect intervals and competitive enrichment.
fig=plt.figure(figsize=(183/25.4,105/25.4),layout='constrained');gs=fig.add_gridspec(1,2,width_ratios=[1.35,1])
ax=fig.add_subplot(gs[0,0]);ay=fig.add_subplot(gs[0,1]);y=np.arange(6)
for a,off in [('GSE4797',-.13),('GSE145467',.13)]:
    cn='JS2_vs_JS10' if a=='GSE4797' else 'impaired_vs_normal';d=E[a].query('contrast==@cn').set_index('set').loc[P]
    ax.errorbar(d.effect,y+off,xerr=[d.effect-d.CI_low,d.CI_high-d.effect],fmt='o',ms=4,color=COHORT[a],lw=.9,capsize=2,label='JS 2 vs JS 10' if a=='GSE4797' else 'Impaired vs normal')
ax.set_yticks(y,[LABEL[n] for n in P]);ax.invert_yaxis();ax.axvline(0,color='#BEBEBE',lw=.7);ax.set_xlabel('Difference in mean gene z score (95% CI)');ax.set_title('Pathway score differences',loc='left');ax.legend(frameon=False,loc='lower right');letter(ax,'a')
qs=np.stack([CA['GSE4797'].query('contrast==@cn').set_index('set').loc[P].FDR_family.to_numpy() for cn in ['JS8_vs_JS10','JS5_vs_JS10','JS2_vs_JS10']]+[CA['GSE145467'].set_index('set').loc[P].FDR_family.to_numpy()],axis=1)
im=ay.imshow(-np.log10(qs),cmap='YlOrBr',vmin=0,vmax=5,aspect='auto')
for i in range(6):
    for j in range(4):ay.text(j,i,f'{qs[i,j]:.3g}',ha='center',va='center',fontsize=6,color='white' if -np.log10(qs[i,j])>3 else '#292929',weight='bold' if qs[i,j]<.05 else 'normal')
ay.set_yticks(y,[LABEL[n] for n in P]);ay.yaxis.tick_right();ay.set_xticks(range(4),['JS 8','JS 5','JS 2','Impaired'],rotation=35,ha='right');ay.set_title('Competitive enrichment: FDR',loc='left');ay.set_xlabel('vs JS 10                  vs normal',fontsize=6);letter(ay,'b')
cb=fig.colorbar(im,ax=ay,orientation='horizontal',shrink=.8,pad=.13);cb.set_label('−log10(FDR); fixed gene correlation = 0.01',fontsize=6)
save(fig,'Figure_3_Cross_Cohort_Pathways',['04_results/primary/*_module_effects.tsv','04_results/primary/*_camera.tsv'])

# 4. Full shared ECM universe and every measured laminin member.
m=read('04_results/primary/cross_cohort_gene_effects.tsv.gz');mm=m.query('contrast_GSE4797=="JS2_vs_JS10" and in_fixed_ECM_union').copy()
fixed=json.loads((ROOT/'02_annotation/fixed_gene_sets.json').read_text());lam=fixed['primary']['REACTOME_LAMININ_INTERACTIONS']
fig=plt.figure(figsize=(183/25.4,175/25.4),layout='constrained');gs=fig.add_gridspec(1,2,width_ratios=[1.35,1])
ax=fig.add_subplot(gs[0,0]);ay=fig.add_subplot(gs[0,1]);both=(mm['adj.P.Val_GSE4797']<.05)&(mm['adj.P.Val_GSE145467']<.05)&(np.sign(mm.logFC_GSE4797)==np.sign(mm.logFC_GSE145467))
ax.scatter(mm.logFC_GSE4797,mm.logFC_GSE145467,s=13,c=np.where(both,'#53758B','#C6CED3'),alpha=.8,edgecolors='none')
highlight=mm[mm.gene.isin(lam)];ax.scatter(highlight.logFC_GSE4797,highlight.logFC_GSE145467,s=30,facecolor='none',edgecolor='#A56149',lw=.8,label='Laminin pathway members')
ax.axhline(0,color='#B7B7B7',lw=.7);ax.axvline(0,color='#B7B7B7',lw=.7);ax.set_xlabel('GSE4797: JS 2 − JS 10 (log2 difference)');ax.set_ylabel('GSE145467: impaired − normal (log2 ratio difference)');ax.set_title('Common measurable ECM genes',loc='left')
summary=read('04_results/primary/cross_cohort_direction_summary.tsv').query('contrast=="JS2_vs_JS10" and scope=="fixed_ECM_union"').iloc[0]
ax.text(.04,.97,f'n = {summary.genes}\nSame direction: {summary.same_sign_n}/{summary.genes} ({summary.same_sign_fraction:.1%})\nSpearman ρ = {summary.spearman_logFC:.2f}',transform=ax.transAxes,va='top',fontsize=7)
ax.legend(frameon=False,loc='lower right',fontsize=6);letter(ax,'a');ax.set_box_aspect(1.1)
dat=[]
for cn in ['JS8_vs_JS10','JS5_vs_JS10','JS2_vs_JS10']:
    z=m[m.contrast_GSE4797==cn].set_index('gene');dat.append(z.logFC_GSE4797)
dat.append(m[m.contrast_GSE4797=='JS2_vs_JS10'].set_index('gene').logFC_GSE145467)
h=pd.concat(dat,axis=1);h=h.loc[h.index.isin(lam)];h.columns=['JS 8','JS 5','JS 2','Impaired'];h=h.sort_values('JS 2',ascending=False)
im=ay.imshow(h,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto');ay.set_yticks(range(len(h)),h.index,fontsize=6);ay.set_xticks(range(4),h.columns,rotation=35,ha='right');ay.set_title(f'Laminin pathway: {len(h)} common genes',loc='left');letter(ay,'b')
cb=fig.colorbar(im,ax=ay,orientation='horizontal',pad=.1,shrink=.8);cb.set_label('log2 expression difference (color capped at ±3)',fontsize=6)
h.to_csv(OUT/'Figure_4_laminin_source_data.tsv',sep='\t',index_label='gene')
save(fig,'Figure_4_Gene_Direction',['04_results/primary/cross_cohort_gene_effects.tsv.gz','02_annotation/fixed_gene_sets.json','04_results/primary/cross_cohort_direction_summary.tsv'])

# 5. Composition clues and statistical sensitivity.
fig,axs=plt.subplots(2,2,figsize=(183/25.4,145/25.4),layout='constrained')
markers=['Sertoli_markers','Leydig_markers','Peritubular_markers','Meiotic_markers','Postmeiotic_markers']
ax=axs[0,0];pal=['#4F6F89','#A0785E','#9B929C','#749486','#485555']
for nm,c in zip(markers,pal):
    means=[SC['GSE4797'].loc[nm,S['GSE4797'].loc[S['GSE4797'].group==g,'gsm']].mean() for g in groups]
    ax.plot(range(4),means,'o-',ms=3,lw=1,color=c,label=nm.replace('_markers',''))
labels(ax,groups);ax.set_ylabel('Mean gene z score');ax.set_title('Exploratory cell-marker programs',loc='left');ax.legend(frameon=False,fontsize=6,ncol=2,loc='lower left');letter(ax,'a')
for ax,a,let in [(axs[0,1],'GSE4797','b'),(axs[1,0],'GSE145467','c')]:
    ss=S[a];xx=SC[a].loc['Postmeiotic_markers',ss.gsm].to_numpy(float);yy=SC[a].loc['REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION',ss.gsm].to_numpy(float)
    for g in ss.group.unique():
        mask=ss.group.to_numpy()==g;ax.scatter(xx[mask],yy[mask],s=22,color=COL[g],label=g,edgecolor='white',linewidth=.3)
    ax.set_xlabel('Postmeiotic-marker score');ax.set_ylabel('ECM organization score');ax.set_title(f'{a}: Pearson r = {np.corrcoef(xx,yy)[0,1]:.2f}',loc='left');ax.legend(frameon=False,fontsize=6,ncol=2);letter(ax,let)
ax=axs[1,1]
for i,a in enumerate(ACCS):
    cn='JS2_vs_JS10' if a=='GSE4797' else 'impaired_vs_normal'
    q1=CA[a].query('contrast==@cn').set_index('set').loc[P].FDR_family
    q2=read(f'04_results/sensitivity/{a}_camera_estimated_correlation.tsv').query('contrast==@cn').set_index('set').loc[P].FDR_family
    x=-np.log10(q1.to_numpy());y2=-np.log10(q2.to_numpy());ax.scatter(x,y2,color=COHORT[a],s=25,label=a)
    j=P.index('REACTOME_LAMININ_INTERACTIONS');ax.annotate('Laminin',(x[j],y2[j]),xytext=(3,5),textcoords='offset points',fontsize=6)
ax.axhline(-np.log10(.05),ls='--',lw=.7,color='#969696');ax.axvline(-np.log10(.05),ls='--',lw=.7,color='#969696');ax.set_xlabel('−log10(FDR), fixed correlation 0.01');ax.set_ylabel('−log10(FDR), estimated correlation');ax.set_title('Sensitivity to gene correlation',loc='left');ax.legend(frameon=False,fontsize=6);letter(ax,'d')
save(fig,'Figure_5_Composition_Robustness',['04_results/primary/*_module_scores.tsv','04_results/primary/*_camera.tsv','04_results/sensitivity/*_camera_estimated_correlation.tsv','02_annotation/*_samples.tsv'])

with PdfPages(OUT/'Bulk_Figures_v1_20260909.pdf') as pdf:
    for name,fig in figures:pdf.savefig(fig)
(OUT/'figure_source_manifest.json').write_text(json.dumps(sources,indent=2))
(ROOT/'env/python_versions.json').write_text(json.dumps({m:importlib.metadata.version(m) for m in ['numpy','pandas','scipy','matplotlib']},indent=2))
plt.close('all')
print('Generated',len(figures),'figures and combined PDF')
