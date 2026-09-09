from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import t

S=Path(__file__).resolve().parents[1];O=S/'05_figures'
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':8,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.5,'lines.linewidth':.8})
blue='#4E7189';ochre='#B27654';colors=[blue,ochre]
names=['GOCC_BASEMENT_MEMBRANE','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_LAMININ_INTERACTIONS','REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
labels=['Basement membrane','ECM organization','Collagen formation','Laminin interactions','Integrin interactions','ECM degradation']
lam='REACTOME_LAMININ_INTERACTIONS';manifest={}
def panel(ax,label,title):ax.text(-.13,1.06,label,transform=ax.transAxes,fontweight='bold',fontsize=10);ax.set_title(title,loc='left',pad=8)
def qfmt(q):return 'NA' if pd.isna(q)else '<0.001' if q<.001 else f'{q:.3f}'
def save(fig,name,sources):
    for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=300,facecolor='white')
    manifest[name]={'width_mm':round(fig.get_size_inches()[0]*25.4,2),'height_mm':round(fig.get_size_inches()[1]*25.4,2),'sources':sources,'files':{ext:hashlib.sha256((O/(name+'.'+ext)).read_bytes()).hexdigest()for ext in ['pdf','svg','png']}}
    plt.close(fig)
def qheat(ax,values,rowlabels,collabels):
    values=np.asarray(values,dtype=float);plot=-np.log10(np.maximum(values,1e-5));masked=np.ma.masked_invalid(np.minimum(plot,5));cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#eeeeee')
    ax.imshow(masked,cmap=cmap,vmin=0,vmax=5,aspect='auto')
    ax.set_xticks(range(len(collabels)),collabels);ax.set_yticks(range(len(rowlabels)),rowlabels)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):ax.text(j,i,qfmt(values[i,j]),ha='center',va='center',fontsize=6.1,color='white'if plot[i,j]>2.9 else 'black',fontweight='bold'if values[i,j]<.05 else 'normal')
    ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)

def figure7():
    ef=pd.read_csv(S/'04_results/bulk/external_score_effects.tsv',sep='\t');ca=pd.read_csv(S/'04_results/bulk/external_camera.tsv',sep='\t')
    fig=plt.figure(figsize=(183/25.4,178/25.4));gs=fig.add_gridspec(2,2,height_ratios=[1.05,1.05],left=.25,right=.98,top=.89,bottom=.17,hspace=.58,wspace=.55)
    ax=fig.add_subplot(gs[0,:]);panel(ax,'a','Fixed-program effects in two external biopsy cohorts')
    for j,acc in enumerate(['GSE9210','GSE108886']):
        a=ef[ef.cohort==acc].set_index('program').loc[names];y=np.arange(6)+(j-.5)*.2
        ax.errorbar(a.effect,y,xerr=np.vstack([a.effect-a.CI_low,a.CI_high-a.effect]),fmt='o',ms=3,color=colors[j],capsize=1.8,label=acc+(' · 47 NOA / 11 OA'if j==0 else' · 8 NOA / 3 OA'))
    ax.set_yticks(range(6),labels);ax.invert_yaxis();ax.axvline(0,color='#aaaaaa',lw=.6);ax.set_xlabel('Program-score difference (NOA − OA), 95% CI');ax.texts[0].set_position((-.13,1.20));ax.legend(loc='lower center',bbox_to_anchor=(.5,1.00),ncol=2,fontsize=6,frameon=False);ax.set_title('Fixed-program effects in two external biopsy cohorts',loc='left',pad=25)
    ax=fig.add_subplot(gs[1,0]);panel(ax,'b','Individual laminin scores')
    rng=np.random.default_rng(20260910)
    for j,acc in enumerate(['GSE9210','GSE108886']):
        sm=pd.read_csv(S/'01_metadata'/(acc+'_analysis_samples.tsv'),sep='\t');scores=pd.read_csv(S/'04_results/bulk'/(acc+'_scores.tsv'),sep='\t',index_col=0).loc[lam]
        for k,grp in enumerate(['OA','NOA']):
            values=scores[sm.loc[sm.group==grp,'gsm']].to_numpy();pos=j*3+k
            ax.scatter(pos+rng.uniform(-.16,.16,len(values)),values,s=7,alpha=.75,color='#a3afb8'if k==0 else colors[j],edgecolors='none')
            mean=values.mean();half=t.ppf(.975,len(values)-1)*values.std(ddof=1)/np.sqrt(len(values));ax.errorbar(pos,mean,yerr=half,fmt='_',ms=9,color='black',capsize=3,lw=1)
    ax.set_xticks([0,1,3,4],['OA\n11','NOA\n47','OA\n3','NOA\n8']);ax.set_ylabel('Mean member-gene z score')
    ax.text(.5,-.29,'GSE9210\n9 / 30 members',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=6.5);ax.text(3.5,-.29,'GSE108886\n28 / 30 members',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=6.5)
    ax=fig.add_subplot(gs[1,1]);panel(ax,'c','Competitive enrichment FDR')
    vals=[]
    for nm in names:
        row=[]
        for acc in ['GSE9210','GSE108886']:
            for mode in ['estimated','fixed_0.01']:
                a=ca[(ca.cohort==acc)&(ca.program==nm)&(ca['mode']==mode)];row.append(a.FDR_six.iloc[0]if len(a)else np.nan)
        vals.append(row)
    qheat(ax,vals,['BM','ECM','Collagen','Laminin','Integrin','Degradation'],['9210\nEst.','9210\nFixed','108886\nEst.','108886\nFixed'])
    ax.text(.5,-.23,'NA: fewer than 10 measured members\nBold: FDR < 0.05; six-test families',ha='center',va='top',transform=ax.transAxes,fontsize=6.2)
    save(fig,'Figure_7_External_Bulk_Validation',['04_results/bulk/external_score_effects.tsv','04_results/bulk/external_camera.tsv','04_results/bulk/*_scores.tsv','01_metadata/external_program_coverage.tsv'])

def figure8():
    donors=['N1','N2','N3','Cr1','Cr2','Cr3'];sm=pd.read_csv(S/'03_processed/single_cell/all_pseudobulk_samples.tsv',sep='\t');frac=pd.read_csv(S/'04_results/peritubular_captured_state_fractions.tsv',sep='\t').set_index('donor').loc[donors]
    ca=pd.read_csv(S/'04_results/single_cell/program_camera_and_effects.tsv',sep='\t');g=pd.read_csv(S/'04_results/single_cell/all_30_laminin_gene_effects.tsv',sep='\t')
    fig=plt.figure(figsize=(183/25.4,244/25.4));gs=fig.add_gridspec(3,2,height_ratios=[.72,.74,2.0],left=.22,right=.98,top=.95,bottom=.09,hspace=.68,wspace=.35)
    ax=fig.add_subplot(gs[0,0]);panel(ax,'a','Somatic donor coverage')
    cs=['Peritubular_compartment','Sertoli','Leydig','Endothelial','Macrophage','Perivascular'];a=sm[(sm.variant=='primary')&sm.compartment.isin(cs)].pivot(index='compartment',columns='donor',values='n_cells').reindex(index=cs,columns=donors)
    ax.imshow((a>=30).astype(int),cmap=ListedColormap(['#eeeeee','#b9ccdb']),vmin=0,vmax=1,aspect='auto')
    for i in range(len(cs)):
        for j in range(6):ax.text(j,i,str(a.iloc[i,j]),ha='center',va='center',fontsize=6)
    ax.set_xticks(range(6),donors,rotation=45);ax.set_yticks(range(6),['Peritubular¹','Sertoli','Leydig','Endothelial','Macrophage','Perivascular']);ax.tick_params(length=0);ax.spines[['left','bottom']].set_visible(False)
    ax=fig.add_subplot(gs[0,1]);panel(ax,'b','Captured peritubular states')
    y=frac.author_fibrotic_state_fraction_within_peritubular.to_numpy();ax.bar(range(6),1-y,color=blue,label='Author PMC state');ax.bar(range(6),y,bottom=1-y,color=ochre,label='Author fibrotic PMC state');ax.set_xticks(range(6),donors,rotation=45);ax.set_ylim(0,1);ax.set_ylabel('Captured-cell fraction');ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),fontsize=5.8,frameon=False,ncol=1)
    ax=fig.add_subplot(gs[1,:]);panel(ax,'c','Correlation-estimated competitive FDR across prespecified analyses')
    scenarios=['primary__Peritubular_compartment__unadjusted','primary__Author_PMC_state__unadjusted','primary__Peritubular_compartment__age_adjusted','primary__Author_PMC_state__age_adjusted','mt20__Peritubular_compartment__unadjusted','mt20__Author_PMC_state__unadjusted']
    q=ca[ca['mode']=='estimated'].pivot(index='program',columns='scenario',values='FDR_family').reindex(index=names,columns=scenarios)
    qheat(ax,q.to_numpy(),labels,['Peritubular¹','PMC state','Peritubular¹\n+ age','PMC state\n+ age','Peritubular¹\nMT < 20%','PMC state\nMT < 20%'])
    genes=sorted(json.loads((S/'01_metadata/fixed_gene_sets.json').read_text())['primary'][lam])
    for j,scenario in enumerate(scenarios[:2]):
        ax=fig.add_subplot(gs[2,j]);panel(ax,'d'if j==0 else'e','All laminin members · '+('peritubular¹'if j==0 else'PMC state'))
        a=g[g.scenario==scenario].set_index('gene').reindex(genes);yy=np.arange(30);ok=a.measured.fillna(False).to_numpy(dtype=bool)
        ax.errorbar(a.loc[ok,'logFC'],yy[ok],xerr=np.vstack([a.loc[ok,'logFC']-a.loc[ok,'CI.L'],a.loc[ok,'CI.R']-a.loc[ok,'logFC']]),fmt='o',ms=2.4,color=colors[j],capsize=1,lw=.6)
        ax.axvline(0,color='#aaaaaa',lw=.6);ax.set_yticks(yy,genes if j==0 else ['']*30);ax.invert_yaxis();ax.set_xlim(-4.8,3.6);ax.set_xlabel('Crypto − OA log2 fold change\nModerated 95% CI; 3 donors per group')
        for i in yy[~ok]:ax.text(.03,i,'NA',transform=ax.get_yaxis_transform(),color='#888888',fontsize=6)
    fig.text(.02,.012,'¹ Author-defined mixture of PMC states. Blue coverage cells meet n ≥ 30; state fractions are captured-cell fractions.',fontsize=6)
    save(fig,'Figure_8_External_Peritubular_Validation',['03_processed/single_cell/all_pseudobulk_samples.tsv','04_results/peritubular_captured_state_fractions.tsv','04_results/single_cell/program_camera_and_effects.tsv','04_results/single_cell/all_30_laminin_gene_effects.tsv'])

def figures6():
    fig=plt.figure(figsize=(183/25.4,200/25.4));gs=fig.add_gridspec(3,2,left=.20,right=.97,top=.95,bottom=.14,hspace=.85,wspace=.5,height_ratios=[1,1,1.15])
    for j,acc in enumerate(['GSE9210','GSE108886']):
        ax=fig.add_subplot(gs[0,j]);panel(ax,'ab'[j],acc+' expression PCA')
        a=pd.read_csv(S/('04_results/bulk/'+acc+'_PCA.tsv'),sep='\t')
        for grp,col in [('OA',blue),('NOA',ochre)]:
            x=a[a.group==grp];ax.scatter(x.PC1,x.PC2,s=8,color=col,label=grp+' (n = '+str(len(x))+')',alpha=.8)
        ax.set_xlabel('PC1');ax.set_ylabel('PC2');ax.legend(fontsize=6,frameon=False,loc='best')
    ax=fig.add_subplot(gs[1,0]);panel(ax,'c','Published-cell QC sensitivity')
    a=pd.read_csv(S/'03_processed/single_cell/all_cell_QC.tsv.gz',sep='\t');donors=['N1','N2','N3','Cr1','Cr2','Cr3']
    aa=a.groupby('Sample_ID').agg(primary=('mt20','size'),mt20=('mt20','sum')).reindex(donors)
    ax.bar(np.arange(6)-.18,aa.primary/1000,width=.36,color=blue,label='Published curated');ax.bar(np.arange(6)+.18,aa.mt20/1000,width=.36,color=ochre,label='MT < 20%')
    ax.set_xticks(range(6),donors,rotation=45);ax.set_ylabel('Cells (thousands)');ax.set_ylim(0,7.5);ax.legend(fontsize=5.7,frameon=False,loc='upper right')
    ax=fig.add_subplot(gs[1,1]);panel(ax,'d','Measured program coverage')
    b=pd.read_csv(S/'01_metadata/external_program_coverage.tsv',sep='\t');c=pd.read_csv(S/'04_results/single_cell/program_member_coverage.tsv',sep='\t');v=[];ns=[]
    for nm in names:
        rs=[b[(b.cohort==acc)&(b.variant=='primary')&(b.program==nm)].iloc[0]for acc in ['GSE9210','GSE108886']]+[c[(c.scenario==sc)&(c.program==nm)].iloc[0]for sc in ['primary__Peritubular_compartment__unadjusted','primary__Author_PMC_state__unadjusted']]
        v.append([r.n_measured/r.n_total for r in rs]);ns.append([r.n_measured for r in rs])
    ax.imshow(v,cmap='Blues',vmin=0,vmax=1,aspect='auto')
    for i in range(6):
        for j in range(4):ax.text(j,i,str(ns[i][j]),ha='center',va='center',fontsize=6,color='white'if v[i][j]>.7 else 'black')
    ax.set_xticks(range(4),['9210','108886','Peritub.','PMC'],rotation=35);ax.set_yticks(range(6),['BM','ECM','Collagen','Laminin','Integrin','Degradation']);ax.tick_params(length=0)
    ax=fig.add_subplot(gs[2,:]);panel(ax,'e','Independent markers across author-defined somatic identities')
    m=pd.read_csv(S/'04_results/external_independent_marker_review.tsv',sep='\t');types=['PMCs','Fibrotic peritubular myoid cells','Perivascular cells','Sertoli cells','Leydig cells','Endothelial cells','Macrophages'];gsym=['MYH11','TAGLN','CNN1','DES','RGS5','CSPG4','NOTCH3','SOX9','CLDN11','CYP11A1','STAR','EMCN','PTPRC','LST1','TYROBP']
    m=m[(m.variant=='primary')&(m.n_cells>=30)&m.author_cell_type.isin(types)&m.gene.isin(gsym)]
    mat=m.groupby(['author_cell_type','gene']).mean_log1p_CPTT.mean().unstack().reindex(index=types,columns=gsym);n=m.groupby('author_cell_type').donor.nunique().reindex(types)
    scale=mat.div(mat.max(axis=0).replace(0,np.nan),axis=1).fillna(0);ax.imshow(scale,cmap='Blues',vmin=0,vmax=1,aspect='auto');ax.set_xticks(range(len(gsym)),gsym,rotation=55,ha='right');ax.set_yticks(range(7),[l+' ('+str(n[t])+')'for l,t in zip(['PMC','Fibrotic PMC','Perivascular','Sertoli','Leydig','Endothelial','Macrophage'],types)]);ax.tick_params(length=0)
    fig.text(.03,.015,'Marker color: equal-donor mean log1p CP10K, divided by each gene’s maximum. Parentheses: eligible donors, both groups.\nOnly donor/identity units with ≥30 cells are displayed; all marker genes are outside the fixed ECM union. Descriptive QC.',fontsize=6)
    save(fig,'Figure_S6_External_Validation_QC',['04_results/bulk/*_PCA.tsv','03_processed/single_cell/all_cell_QC.tsv.gz','01_metadata/external_program_coverage.tsv','04_results/single_cell/program_member_coverage.tsv','04_results/external_independent_marker_review.tsv'])

if __name__=='__main__':
    figure7();figure8();figures6();(S/'00_admin/external_figure_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('GENERATED',list(manifest))
