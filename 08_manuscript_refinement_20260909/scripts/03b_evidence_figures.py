#!/usr/bin/env python3
"""Integrative figures from complete frozen gene and program estimates."""
from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pypdf import PdfReader,PdfWriter
S=Path(__file__).resolve().parents[1];P=S.parent;SC=P/'07_single_cell';O=S/'04_figures';T=S/'03_tables'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
sets=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary']
order=['GOCC_BASEMENT_MEMBRANE','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_LAMININ_INTERACTIONS','REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
labels=['Basement membrane','ECM organization','Collagen formation','Laminin interactions','Integrin interactions','ECM degradation']
palette=['#56758E','#A76450'];manifest={}
def save(fig,name,sources):
    for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=300)
    manifest[name]={'sources':sources,'width_mm':183,'height_mm':float(fig.get_size_inches()[1]*25.4)}
    plt.close(fig)
def letter(ax,s):ax.text(-.10,1.035,s,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')

# Effect magnitude and competitive inference answer separate questions.
d=pd.read_csv(T/'Table_2_program_effects_and_enrichment.tsv',sep='\t')
fig=plt.figure(figsize=(183/25.4,120/25.4))
gs=fig.add_gridspec(1,2,width_ratios=[1.45,1],left=.24,right=.98,bottom=.22,top=.88,wspace=.26)
ax=fig.add_subplot(gs[0,0]);ay=fig.add_subplot(gs[0,1]);y=np.arange(6)
for k,c in enumerate(['GSE4797','GSE145467']):
    q=d[d.cohort.eq(c)].set_index('pathway').loc[order]
    ax.errorbar(q.effect,y+(k-.5)*.23,xerr=np.vstack([q.effect-q.CI_low,q.CI_high-q.effect]),fmt='o',ms=3.5,lw=.9,capsize=2,color=palette[k],label=['GSE4797: JS2 − JS10 (5 vs 12)','GSE145467: impaired − normal (10 vs 10)'][k])
ax.axvline(0,color='#B8B8B8',lw=.6);ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Program score difference (95% CI)');ax.set_title('Severe-group effect sizes',loc='left');letter(ax,'a')
ax.legend(loc='upper left',bbox_to_anchor=(-.46,-.16),frameon=False,fontsize=6.3)
cols=[]
for c in ['GSE4797','GSE145467']:
    q=d[d.cohort.eq(c)].set_index('pathway').loc[order]
    cols.extend([q.camera_fixed_FDR.values,q.camera_estimated_FDR.values])
qv=np.array(cols).T
im=ay.imshow(-np.log10(qv),cmap='Blues',vmin=0,vmax=5,aspect='auto')
for i in range(6):
    for j in range(4):
        q=qv[i,j];tx=f'{q:.3f}' if q>=.001 else f'{q:.1e}'
        ay.text(j,i,tx,ha='center',va='center',fontsize=6.6,color='white' if -np.log10(q)>2.7 else '#242424',weight='bold' if q<.05 else 'normal')
ay.set_yticks(y,[]);ay.set_xticks(range(4),['Fixed\nGSE4797','Estimated\nGSE4797','Fixed\nGSE145467','Estimated\nGSE145467'],rotation=90)
ay.axvline(1.5,color='white',lw=2);ay.set_title('Competitive camera FDR',loc='left');letter(ay,'b')
fig.text(.98,.02,'Bold: FDR <0.05. Fixed ρ = 0.01; estimated ρ is set-specific.\nFDR families: 18 primary-cohort tests; 6 independent-cohort tests.',ha='right',fontsize=6.2)
save(fig,'Figure_2_Cross_Cohort_Evidence',['03_tables/Table_2_program_effects_and_enrichment.tsv'])

# Every fixed member, ordered by molecular group, keeps distinct data layers visible.
z=pd.read_csv(T/'Table_S3_all_30_laminin_gene_evidence.tsv',sep='\t')
src=pd.read_csv(T/'Table_S4_all_30_laminin_gene_sources.tsv',sep='\t')
types=['Spermatogonia','Spermatocyte','Spermatid','Sertoli','Interstitial_stromal','Peritubular_myoid','Endothelial']
short=['SPG (5)','SPC (5)','SPT (5)','Sertoli (4)','Interstitial (5)','Myoid (4)','Endothelial (3)']
fig=plt.figure(figsize=(183/25.4,225/25.4))
gs=fig.add_gridspec(1,4,width_ratios=[1.05,1.05,1.62,1.1],left=.12,right=.975,bottom=.15,top=.925,wspace=.24)
axs=[fig.add_subplot(gs[0,j]) for j in range(4)];yy=np.arange(30)
for j,c in [(0,'GSE4797'),(1,'GSE145467'),(3,'interstitial')]:
    ax=axs[j];est=z[c+'_logFC'];lo=z[c+'_CI.L'];hi=z[c+'_CI.R'];ok=est.notna()
    ax.errorbar(est[ok],yy[ok],xerr=np.vstack([est[ok]-lo[ok],hi[ok]-est[ok]]),fmt='o',ms=2.4,lw=.65,capsize=1.2,color=palette[j] if j<2 else '#515D64')
    for i in yy[~ok]:ax.text(0,i,'NA',ha='center',va='center',fontsize=5.8,color='#888888')
    ax.axvline(0,color='#AAA',lw=.6);ax.set_xlabel('log2 fold change\n(95% CI)',fontsize=6.5)
    ax.set_ylim(29.7,-.7);ax.set_yticks(yy,z.gene if j==0 else []);ax.xaxis.set_major_locator(plt.MaxNLocator(4));letter(ax,'abcd'[j])
axs[0].set_title('Primary bulk\nJS2 − JS10\n5 vs 12 men',loc='left',fontsize=7)
axs[1].set_title('Independent bulk\nImpaired − normal\n10 vs 10 samples',loc='left',fontsize=7)
axs[3].set_title('Interstitial pseudobulk\niNOA − OA\n3 vs 3 BD donors',loc='left',fontsize=7)
ax=axs[2];m=src.pivot(index='gene',columns='cell_type',values='mean_CP10K').loc[z.gene,types]
raw=np.log1p(m.to_numpy());maxi=raw.max(axis=1,keepdims=True);scaled=np.divide(raw,maxi,out=np.zeros_like(raw),where=maxi>0)
im=ax.imshow(scaled,cmap='YlOrBr',vmin=0,vmax=1,aspect='auto')
ax.set_yticks(yy,[]);ax.set_xticks(range(7),short,rotation=90,fontsize=6);ax.set_title('Adult reference sources\nEqual-donor expression\nRelative within each gene',loc='left',fontsize=7);letter(ax,'c')
for bound in [10.5,18.5,21.5]:
    for ax in axs:ax.axhline(bound,color='#D1D1D1',lw=.55)
cbax=fig.add_axes([.72,.058,.20,.008]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal');cb.set_ticks([0,.5,1]);cb.set_label('Relative log1p mean CP10K',fontsize=6)
fig.text(.12,.058,'Rows: laminin chains, collagens, matrix partners,\nintegrin receptors. NA: unmeasured or not tested.',fontsize=6.2)
fig.text(.12,.024,'SPG / SPC / SPT: germ-cell stages. Source parentheses: eligible reference donors. Full FDR values: Table S3.',fontsize=6)
save(fig,'Figure_6_Laminin_Gene_Evidence',['03_tables/Table_S3_all_30_laminin_gene_evidence.tsv','03_tables/Table_S4_all_30_laminin_gene_sources.tsv'])

# Program overlap and measured coverage are descriptive audit quantities.
over=pd.read_csv(T/'Table_S6_program_Jaccard_overlap.tsv',sep='\t',index_col=0).loc[order,order]
fig,(ax,ay)=plt.subplots(1,2,figsize=(183/25.4,108/25.4),gridspec_kw={'width_ratios':[1,1.1]},layout='constrained')
im=ax.imshow(over,cmap='Blues',vmin=0,vmax=1)
for i in range(6):
    for j in range(6):ax.text(j,i,f'{over.iloc[i,j]:.2f}',ha='center',va='center',fontsize=6.5,color='white' if over.iloc[i,j]>.6 else '#333')
lab=['Basement membrane','ECM organization','Collagen formation','Laminin interactions','Integrin interactions','ECM degradation']
ax.set_xticks(range(6),lab,rotation=90,fontsize=6);ax.set_yticks(range(6),lab,fontsize=6);ax.set_title('Fixed-set Jaccard overlap',loc='left');letter(ax,'a')
cov=pd.read_csv(T/'Table_S7_single_cell_program_coverage.tsv',sep='\t')
grid=[]
for k in order:
    ds=d[d.pathway.eq(k)].set_index('cohort');sc=cov[cov.pathway.eq(k)]
    grid.append([len(sets[k]),int(ds.loc['GSE4797','measured_members']),int(ds.loc['GSE145467','measured_members']),int(sc.measured.sum()),int((sc.reference_cells_detected>0).sum())])
grid=np.array(grid);den=grid[:,[0]];ay.imshow(grid/den,cmap='Greys',vmin=0,vmax=1,aspect='auto')
for i in range(6):
    for j in range(5):ay.text(j,i,str(grid[i,j]),ha='center',va='center',fontsize=7,color='white' if grid[i,j]/grid[i,0]>.6 else '#222')
ay.set_yticks(range(6),[]);ay.set_xticks(range(5),['Fixed set','GSE4797','GSE145467','SC measured','SC detected\nreference'],rotation=90,fontsize=6);ay.set_title('Complete membership and coverage',loc='left');letter(ay,'b')
save(fig,'Figure_S5_Program_Overlap_Coverage',['03_tables/Table_S6_program_Jaccard_overlap.tsv','03_tables/Table_S7_single_cell_program_coverage.tsv'])

# Reuse validated panels with a clear manuscript numbering map.
mapping=[(P/'05_figures','Figure_2_Pathology_ECM','Figure_1_Histological_Programs'),(P/'05_figures','Figure_4_Gene_Direction','Figure_3_Gene_Direction'),(P/'05_figures','Figure_5_Composition_Robustness','Figure_4_Composition_Robustness'),(P/'05_figures','Figure_1_Cohorts_QC','Figure_S1_Bulk_Cohorts_QC'),(SC/'05_figures','Figure_S1_Single_Cell_QC','Figure_S2_Single_Cell_QC'),(SC/'05_figures','Figure_S3_Source_Sensitivity','Figure_S4_Source_Sensitivity')]
for root,old,new in mapping:
    for ext in ['pdf','svg','png']:
        path=root/(old+'.'+ext)
        if path.exists():shutil.copy2(path,O/(new+'.'+ext))
    manifest[new]={'reused_from':str((root/old).relative_to(P)),'evidence_and_geometry':'unchanged'}
for prefix,title in [('Figure_','Main_Figures_1_to_6_v3_20260909'),('Figure_S','Supplementary_Figures_S1_to_S5_v3_20260909')]:
    files=sorted([p for p in O.glob(prefix+'*.pdf') if (p.stem.split('_')[1].isdigit() if prefix=='Figure_' else True)])
    w=PdfWriter()
    for path in files:
        for pg in PdfReader(path).pages:w.add_page(pg)
    with (O/(title+'.pdf')).open('wb') as h:w.write(h)
(S/'00_admin/evidence_figure_manifest_v3.json').write_text(json.dumps(manifest,indent=2))
print('COMPLETE evidence figures and main / supplementary compilations')
