#!/usr/bin/env python3
"""Authoritative adult SC figures from frozen tables, with all members shown."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr
from pypdf import PdfReader,PdfWriter
S=Path(__file__).resolve().parents[1];P=S.parent;O=S/'05_figures';R=S/'04_results'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.2,'ytick.labelsize':6.5,'legend.fontsize':6,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
SET=['GOCC_BASEMENT_MEMBRANE','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_LAMININ_INTERACTIONS','REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
SL=['Basement\nmembrane','ECM\norganization','Collagen\nformation','Laminin\ninteractions','Integrin\ninteractions','ECM\ndegradation']
SHORT={'Interstitial_stromal':'Interstitial','Peritubular_myoid':'Peritubular myoid','Vascular_mural':'Vascular mural','Sertoli':'Sertoli','Endothelial':'Endothelial','Macrophage':'Macrophage','Spermatogonia':'Spermatogonia','Spermatocyte':'Spermatocyte','Spermatid':'Spermatid','Mast_cell':'Mast cell','T_NK':'T / NK','Steroidogenic_unresolved':'Steroidogenic unresolved','Sertoli_mural_mixed':'Sertoli–mural mixed'}
TYPES=['Spermatogonia','Spermatocyte','Spermatid','Sertoli','Interstitial_stromal','Peritubular_myoid','Vascular_mural','Endothelial','Macrophage']
ALLT=TYPES+['Mast_cell','T_NK','Steroidogenic_unresolved','Sertoli_mural_mixed']
COL=dict(zip(ALLT,['#B4BE8A','#85AFA4','#A6BED0','#AE84AA','#C9976D','#B35F57','#887BA7','#638FA7','#9C9983','#9C785C','#6D957D','#D0D0D0','#9C9C9C']))
DONORS=['LZ003','LZ007','LZ013','LZ014','LZ015','LZ017','LZ018','LZ019']
DC=['#617E98']*5+['#A35B55']*3
def read(name):return pd.read_csv(R/name,sep='\t')
def letter(ax,s):ax.text(-.10,1.04,s,transform=ax.transAxes,weight='bold',fontsize=10,va='bottom')
figs=[];manifest={}
def save(fig,name,sources):
 for ext in ['svg','pdf','png']:fig.savefig(O/f'{name}.{ext}',dpi=300)
 figs.append((name,fig));manifest[name]={'sources':sources,'width_mm':183,'height_mm':round(fig.get_size_inches()[1]*25.4,2)}

cells=pd.read_csv(S/'02_annotation/cells_final.tsv.gz',sep='\t',index_col=0)
um=pd.read_csv(S/'03_processed/umap.tsv.gz',sep='\t',index_col=0)
cells=cells.join(um)
cov=read('donor_celltype_coverage.tsv');primary=cov[cov.variant.eq('primary')]
refsum=read('reference_program_source_summary.tsv');genes=read('reference_ECM_gene_source_summary.tsv')
scores=read('ECM_program_donor_source_scores.tsv')
fixed=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary'];lam=fixed['REACTOME_LAMININ_INTERACTIONS']
nref=refsum.groupby('cell_type')['count'].first().to_dict()

# Figure 6: reference embedding, six fixed programs and all 30 laminin members.
fig=plt.figure(figsize=(183/25.4,175/25.4),layout='constrained')
gs=fig.add_gridspec(2,2,height_ratios=[1,1.08],width_ratios=[1,1.12],hspace=.14,wspace=.08)
ax=fig.add_subplot(gs[0,0]);ay=fig.add_subplot(gs[0,1]);az=fig.add_subplot(gs[1,:])
ref=cells[cells.group.eq('OA')]
for typ in ALLT:
 d=ref[ref.cell_type.eq(typ)]
 ax.scatter(d.UMAP1,d.UMAP2,s=1.8,color=COL[typ],rasterized=True,lw=0,alpha=.7)
abbr=dict(zip(TYPES,['SPG','SPC','SPT','SER','INT','PTM','VSM','END','MAC']))
for typ in TYPES:
 d=ref[ref.cell_type.eq(typ)]
 ax.text(d.UMAP1.median(),d.UMAP2.median(),abbr[typ],fontsize=6.4,weight='bold',ha='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.75,pad=.5))
ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.set_title('Adult reference: 26,005 cells / 5 donors',loc='left');letter(ax,'a')
ax.text(.02,-.10,'SPG / SPC / SPT: germ-cell stages\nSER: Sertoli; INT: interstitial; PTM: myoid\nVSM: vascular mural; END: endothelial; MAC: macrophage',transform=ax.transAxes,fontsize=5.8,va='top')
h=refsum.pivot(index='cell_type',columns='pathway',values='mean').loc[TYPES,SET]
im=ay.imshow(h,cmap='RdBu_r',vmin=-.85,vmax=.85,aspect='auto')
ay.set_yticks(range(len(TYPES)),[f'{abbr[t]} ({nref[t]})' for t in TYPES]);ay.set_xticks(range(6),SL,rotation=90,fontsize=6)
ay.set_title('Program sources (eligible reference donors)',loc='left');letter(ay,'b')
cb=fig.colorbar(im,ax=ay,orientation='horizontal',shrink=.75,pad=.02,fraction=.045);cb.set_label('Mean member-gene source z score',fontsize=6)
g=genes[genes.gene.isin(lam)].set_index(['cell_type','gene'])
X=[];Y=[];color=[];sizes=[]
for i,typ in enumerate(TYPES):
 for j,gene in enumerate(lam):
  r=g.loc[(typ,gene)];X.append(j);Y.append(i);color.append(np.log1p(r.mean_CP10K));sizes.append(40*r.fraction_detected)
norm=Normalize(0,max(color))
dp=az.scatter(X,Y,s=sizes,c=color,cmap='YlOrBr',norm=norm,edgecolors='#6A5741',linewidth=.15)
az.set_xlim(-.7,len(lam)-.3);az.set_ylim(len(TYPES)-.4,-.7)
az.set_xticks(range(len(lam)),lam,rotation=90,fontsize=6.5)
az.set_yticks(range(len(TYPES)),[f'{SHORT[t]} ({nref[t]})' for t in TYPES]);az.set_title('Every laminin-interaction member',loc='left');letter(az,'c')
az.grid(axis='y',color='#E7E7E7',lw=.4);az.set_axisbelow(True)
cb=fig.colorbar(dp,ax=az,orientation='horizontal',shrink=.30,pad=.01,fraction=.045);cb.set_label('log1p mean CP10K (equal donor weight)',fontsize=6)
handles=[az.scatter([],[],s=40*f,facecolor='#A89070',edgecolor='#6A5741',lw=.15,label=f'{100*f:.0f}%') for f in [.1,.5,.9]]
az.legend(handles=handles,title='Detected cells',ncol=3,loc='upper right',bbox_to_anchor=(1.0,-.31),frameon=False,title_fontsize=6,handletextpad=.5,columnspacing=1)
h.to_csv(O/'Figure_6_program_source_data.tsv',sep='\t')
genes[genes.gene.isin(lam)].to_csv(O/'Figure_6_all_laminin_source_data.tsv',sep='\t',index=False)
save(fig,'Figure_6_Adult_ECM_Sources',['cells_final.tsv.gz','umap.tsv.gz','reference_program_source_summary.tsv','reference_ECM_gene_source_summary.tsv'])

# S1: original-cell metrics and explicit filtering/doublet flags.
raw=[];receipts=[]
for donor in DONORS:
 raw.append(pd.read_csv(S/'02_annotation'/f'{donor}_all_cell_QC.tsv.gz',sep='\t'))
 receipts.append(json.loads((S/'00_admin'/f'{donor}_prepare_receipt.json').read_text()))
q=pd.concat(raw,ignore_index=True);qr=pd.DataFrame(receipts).set_index('donor_id').loc[DONORS]
fig,axs=plt.subplots(3,2,figsize=(183/25.4,190/25.4),layout='constrained')
ax=axs[0,0];xx=np.arange(8)
ax.bar(xx,qr.raw_cells,color='#E0E3E6',label='Uploaded cells');ax.bar(xx,qr.primary_qc_cells,color=DC,width=.55,label='Pass primary QC')
for i,v in enumerate(qr.primary_qc_cells):ax.text(i,v+120,str(v),ha='center',fontsize=5.8)
ax.set_ylabel('Cells');ax.set_title('Raw-to-QC retention',loc='left');ax.set_ylim(0,12000);ax.legend(frameon=False,fontsize=5.8);letter(ax,'a')
def bx(ax,col,title,ylabel,log=False):
 bp=ax.boxplot([q.loc[q.donor_id.eq(d),col] for d in DONORS],positions=xx,widths=.55,showfliers=False,patch_artist=True,medianprops={'color':'#2F3438','linewidth':.7})
 for patch,c in zip(bp['boxes'],DC):patch.set_facecolor(c);patch.set_alpha(.55)
 ax.set_title(title,loc='left');ax.set_ylabel(ylabel)
 if log:ax.set_yscale('log')
bx(axs[0,1],'total_counts_original','Original library counts','Counts / cell',True);letter(axs[0,1],'b')
bx(axs[1,0],'n_features_original','Detected original features','Features / cell');axs[1,0].axhline(500,color='#777',ls='--',lw=.7);axs[1,0].axhline(9000,color='#777',ls='--',lw=.7);letter(axs[1,0],'c')
bx(axs[1,1],'pct_mito_original','Mitochondrial counts','Mitochondrial counts (%)');axs[1,1].axhline(40,color='#777',ls='--',lw=.7);axs[1,1].axhline(20,color='#A35B55',ls=':',lw=.8);letter(axs[1,1],'d')
ax=axs[2,0];rate=[cells.loc[cells.donor_id.eq(d),'predicted_doublet'].mean()*100 for d in DONORS]
ax.scatter(xx,rate,c=DC,s=24);ax.set_ylim(0,max(rate)+1.3);ax.set_ylabel('Predicted doublets (%)');ax.set_title('Scrublet flags after primary QC',loc='left');letter(ax,'e')
ax=axs[2,1];bottom=np.zeros(8)
for typ in ALLT:
 d=primary[primary.cell_type.eq(typ)].set_index('donor_id').loc[DONORS]
 ax.bar(xx,d.capture_fraction,bottom=bottom,color=COL[typ],label=SHORT[typ]);bottom+=d.capture_fraction.to_numpy()
ax.set_ylim(0,1);ax.set_ylabel('Capture fraction');ax.set_title('Captured cell composition',loc='left');letter(ax,'f')
for ax in axs.flat:ax.set_xticks(xx,DONORS,rotation=45,ha='right');ax.tick_params(axis='x',labelsize=6)
fig.legend(*axs[2,1].get_legend_handles_labels(),loc='outside lower center',ncol=4,frameon=False,fontsize=5.5)
save(fig,'Figure_S1_Single_Cell_QC',['*_all_cell_QC.tsv.gz','*_prepare_receipt.json','donor_celltype_coverage.tsv','cells_final.tsv.gz'])

# S2: independently chosen annotation markers and all donor/type counts.
fig,(ax,ay)=plt.subplots(2,1,figsize=(183/25.4,210/25.4),layout='constrained',gridspec_kw={'height_ratios':[1,1.1],'hspace':.10})
m=read('annotation_display_means.tsv')
b=m.set_index(['cell_type','gene'])
marks='UTF1 MAGEA4 SYCP3 HORMAD1 PRM1 ACRV1 SOX9 FATE1 CLDN11 INSL3 STAR CYP17A1 DPEP1 MYH11 RGS5 MCAM NOTCH3 CLDN5 KDR GNG11 TYROBP FCER1G LST1 CPA3 MS4A2 CD3D NKG7'.split()
marks=[x for x in marks if x in b.index.get_level_values('gene')]
mx=b['mean'].groupby('gene').max()
for i,typ in enumerate(ALLT):
 for j,gene in enumerate(marks):
  v=b.loc[(typ,gene)];ax.scatter(j,i,s=38*v['fraction'],c=[v['mean']/max(mx[gene],1e-9)],cmap='Blues',vmin=0,vmax=1,edgecolors='#657789',lw=.1)
ax.set_xlim(-.7,len(marks)-.3);ax.set_ylim(len(ALLT)-.4,-.7);ax.set_yticks(range(len(ALLT)),[SHORT[t] for t in ALLT]);ax.set_xticks(range(len(marks)),marks,rotation=90);ax.set_title('Independent markers for annotation (all eight donors)',loc='left');letter(ax,'a')
ax.text(0,-.29,'Color: expression relative to each gene maximum; size: detected-cell fraction. Cell-level annotation display.',transform=ax.transAxes,fontsize=6)
co=primary.pivot(index='cell_type',columns='donor_id',values='n_cells').loc[ALLT,DONORS]
im=ay.imshow(np.log1p(co),cmap='Greys',aspect='auto',vmin=0,vmax=np.log1p(co.to_numpy().max()))
for i in range(len(ALLT)):
 for j in range(8):
  n=int(co.iloc[i,j]);ay.text(j,i,str(n),ha='center',va='center',fontsize=6.7,color='white' if n>100 else '#333333',weight='bold' if n>=30 else 'normal')
ay.set_xticks(range(8),[d+'\n'+('10X' if j<2 else 'BD') for j,d in enumerate(DONORS)],fontsize=6.5)
ay.set_yticks(range(len(ALLT)),[SHORT[t] for t in ALLT]);ay.axvline(4.5,color='#A35B55',lw=1.2)
ay.set_title('Donor coverage (bold: ≥30 cells; rightmost three donors: iNOA)',loc='left');letter(ay,'b')
save(fig,'Figure_S2_Annotation_Coverage',['annotation_display_means.tsv','cluster_labels.tsv','donor_celltype_coverage.tsv'])

# S3: reference donor variation and sensitivity; within-interstitial pseudobulk.
fig=plt.figure(figsize=(183/25.4,150/25.4),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1,1.2],hspace=.13,wspace=.16)
ax=fig.add_subplot(gs[:,0]);ay=fig.add_subplot(gs[0,1]);az=fig.add_subplot(gs[1,1])
lamid='REACTOME_LAMININ_INTERACTIONS'
d=scores[scores.variant.eq('primary') & scores.group.eq('OA') & scores.source_eligible & scores.pathway.eq(lamid)]
for i,typ in enumerate(TYPES):
 sub=d[d.cell_type.eq(typ)]
 for k,r in enumerate(sub.itertuples()):
  ax.scatter(r.score,i+(k-(len(sub)-1)/2)*.10,s=22,marker='o' if r.technology=='10X' else '^',color=COL[typ],edgecolor='white',lw=.3)
 ax.plot([sub.score.mean()]*2,[i-.23,i+.23],color='black',lw=.8)
ax.axvline(0,color='#CCC',lw=.6);ax.set_yticks(range(9),[f'{SHORT[t]} ({nref[t]})' for t in TYPES]);ax.invert_yaxis();ax.set_xlabel('Laminin source z score');ax.set_title('Every eligible reference donor',loc='left');letter(ax,'a')
for marker,label in [('o','10X'),('^','BD Rhapsody')]:ax.scatter([],[],marker=marker,color='#666',label=label,s=22)
ax.legend(frameon=False,loc='lower right')
base=scores[scores.variant.eq('primary')&scores.group.eq('OA')&scores.source_eligible].groupby(['cell_type','pathway']).score.mean()
sens=[]
for variant,color,mark,label in [('mt20','#AE6B55','o','Mito <20%'),('no_predicted_doublets','#617E98','^','Doublet flags excluded')]:
 z=scores[scores.variant.eq(variant)&scores.group.eq('OA')&scores.source_eligible].groupby(['cell_type','pathway']).score.mean()
 pair=pd.concat([base.rename('primary'),z.rename('sensitivity')],axis=1).dropna()
 rho=spearmanr(pair.primary,pair.sensitivity).statistic
 ay.scatter(pair.primary,pair.sensitivity,s=11,color=color,marker=mark,alpha=.7,label=f'{label}: ρ={rho:.3f}')
 sens.append(dict(variant=variant,unit='equal-donor cell-type/program source summary',paired_units=len(pair),spearman=float(rho)))
 pair.to_csv(O/f'Figure_S3_{variant}_source_data.tsv',sep='\t')
ay.plot([-.65,.95],[-.65,.95],color='#AAA',lw=.7);ay.set_xlabel('Primary source summary');ay.set_ylabel('Sensitivity source summary');ay.set_title('Source patterns across QC choices',loc='left');ay.legend(frameon=False,fontsize=5.8);letter(ay,'b')
cam=read('pseudobulk_camera.tsv')
for variant,color,off in [('primary','#303B43',-.16),('mt20','#AE6B55',0),('no_predicted_doublets','#617E98',.16)]:
 z=cam[cam.variant.eq(variant)].set_index('pathway').loc[SET]
 az.scatter(z.mean_member_logFC,np.arange(6)+off,s=16,color=color,label={'primary':'Primary','mt20':'Mito <20%','no_predicted_doublets':'Doublet flags excluded'}[variant])
az.axvline(0,color='#AAA',lw=.7);az.set_yticks(range(6),[x.replace('\n',' ') for x in SL],fontsize=6);az.invert_yaxis();az.set_xlabel('Mean measured-member log2 fold change');az.set_title('Interstitial pseudobulk: iNOA − OA (3 vs 3)',loc='left');letter(az,'c')
az.text(0,-.32,'Primary camera FDR = 0.929 for all six sets.\nDots summarize member-gene effects; they are not pathway CIs.',transform=az.transAxes,fontsize=5.9)
az.legend(frameon=False,loc='upper left',fontsize=5.4,handletextpad=.3)
save(fig,'Figure_S3_Source_Sensitivity',['ECM_program_donor_source_scores.tsv','pseudobulk_camera.tsv'])

with PdfPages(O/'Single_Cell_Figures_v2_20260909.pdf') as pdf:
 for _,fig in figs:pdf.savefig(fig)
writer=PdfWriter()
for p in [P/'05_figures/Bulk_Figures_v1_20260909.pdf',O/'Figure_6_Adult_ECM_Sources.pdf']:
 for page in PdfReader(p).pages:writer.add_page(page)
with (O/'Main_Figures_1_to_6_v2_20260909.pdf').open('wb') as h:writer.write(h)
supp=PdfWriter()
for name in ['Figure_S1_Single_Cell_QC','Figure_S2_Annotation_Coverage','Figure_S3_Source_Sensitivity']:
 for page in PdfReader(O/(name+'.pdf')).pages:supp.add_page(page)
with (O/'Supplementary_Figures_S1_to_S3_v2_20260909.pdf').open('wb') as h:supp.write(h)
(S/'00_admin/figure_manifest.json').write_text(json.dumps(manifest,indent=2))
(S/'04_results/source_sensitivity_summary.json').write_text(json.dumps(sens,indent=2))
for _,fig in figs:plt.close(fig)
print('COMPLETE four SC figures and combined main figures',flush=True)
