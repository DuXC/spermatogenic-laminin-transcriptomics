from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
Q=Path(__file__).resolve().parents[1]
C=Q/'02_sources/common_member'
F=Q/'03_figures';F.mkdir(parents=True,exist_ok=True)
effects=pd.read_csv(C/'all_program_contrasts.tsv',sep='\t');effects=effects[effects.program.eq('REACTOME_LAMININ_INTERACTIONS')].set_index('cohort')
loo=pd.read_csv(C/'laminin_leave_one_gene_out.tsv',sep='\t');order=list(effects.index);lab=['JS2 − JS10','Impaired − normal','NOA − OA','NOA − OA']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(183/25.4,108/25.4),gridspec_kw={'left':.19,'right':.97,'top':.84,'bottom':.31,'wspace':.30})
y=np.arange(4)[::-1];e=effects.loc[order]
axs[0].errorbar(e.effect,y,xerr=[e.effect-e.CI_low,e.CI_high-e.effect],fmt='o',ms=4,color='#385F77',capsize=2,lw=.9)
axs[0].set_yticks(y,[c+'\n'+l for c,l in zip(order,lab)]);axs[0].set_title('a   Shared eight-member score',loc='left',fontsize=9);axs[0].set_xlabel('Difference with 95% CI')
for i,c in enumerate(order):
 vals=loo.loc[loo.cohort.eq(c),'effect'];axs[1].plot([vals.min(),vals.max()],[y[i],y[i]],color='#385F77',lw=2);axs[1].scatter(vals,[y[i]]*len(vals),s=13,color='#718DA2',edgecolor='white',linewidth=.3,zorder=3);axs[1].plot(e.loc[c,'effect'],y[i],'|',color='black',ms=10,mew=1.2,zorder=4)
axs[1].set_yticks(y,[]);axs[1].set_title('b   Gene omission influence',loc='left',fontsize=9);axs[1].set_xlabel('Range of point estimates')
for ax in axs:ax.axvline(0,color='#777777',lw=.6,ls='--');ax.set_xlim(-.08,2.5);ax.set_ylim(-.6,3.6)
genes=sorted(loo.omitted_gene.unique())
fig.text(.03,.965,'Laminin tissue direction with identical coverage across four cohorts',fontsize=10,va='top')
fig.text(.03,.19,'Common members: '+', '.join(genes),fontsize=7)
fig.text(.03,.13,'Panel b: each point omits one gene; blue segment spans all eight point estimates. Black mark: all eight genes.\nThese spans are influence ranges, not confidence intervals. Full omission intervals are in Table S24.\nPost hoc coverage sensitivity; cohort models and biological comparisons remain separate.',fontsize=7,va='top',linespacing=1.5)
for ext in ['svg','pdf','png']:fig.savefig(F/('Figure_S10.'+ext),dpi=200,facecolor='white')
plt.close(fig)
