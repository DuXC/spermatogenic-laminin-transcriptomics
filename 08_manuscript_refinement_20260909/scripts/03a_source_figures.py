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
V=Path(__file__).resolve().parents[1];P=V.parent;S=P/'07_single_cell';O=V/'04_figures';R=S/'04_results'
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
fixed=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary'];lam=[g for v in json.loads((V/'00_admin/laminin_order.json').read_text()).values() for g in v]
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
h.to_csv(O/'Figure_5_program_source_data.tsv',sep='\t')
genes[genes.gene.isin(lam)].to_csv(O/'Figure_5_all_laminin_source_data.tsv',sep='\t',index=False)
save(fig,'Figure_5_Adult_ECM_Sources',['cells_final.tsv.gz','umap.tsv.gz','reference_program_source_summary.tsv','reference_ECM_gene_source_summary.tsv'])


# S2: independently chosen annotation markers and all donor/type counts.
fig,(ax,ay)=plt.subplots(2,1,figsize=(183/25.4,210/25.4),layout='constrained',gridspec_kw={'height_ratios':[1,1.1],'hspace':.10})
m=read('annotation_display_means.tsv')
b=m.set_index(['cell_type','gene'])
marks='UTF1 MAGEA4 SYCP3 HORMAD1 PRM1 ACRV1 SOX9 FATE1 CLDN11 INSL3 STAR CYP17A1 DPEP1 MYH11 RGS5 MCAM NOTCH3 CLDN5 KDR GNG11 TYROBP FCER1G LST1 CPA3 MS4A2 CD3D NKG7'.split()
marks=[x for x in marks if x in b.index.get_level_values('gene')]
mx=b['mean'].groupby('gene').max()
for i,typ in enumerate(ALLT):
 for j,gene in enumerate(marks):
  v=b.loc[(typ,gene)];dp=ax.scatter(j,i,s=38*v['fraction'],c=[v['mean']/max(mx[gene],1e-9)],cmap='Blues',vmin=0,vmax=1,edgecolors='#657789',lw=.1)
ax.set_xlim(-.7,len(marks)-.3);ax.set_ylim(len(ALLT)-.4,-.7);ax.set_yticks(range(len(ALLT)),[SHORT[t] for t in ALLT]);ax.set_xticks(range(len(marks)),marks,rotation=90);ax.set_title('Independent markers for annotation (all eight donors)',loc='left');letter(ax,'a')
cb=fig.colorbar(dp,ax=ax,orientation='horizontal',shrink=.30,pad=.015,fraction=.035);cb.set_label('Expression relative to gene maximum',fontsize=6);cb.set_ticks([0,.5,1])
handles=[ax.scatter([],[],s=38*f,color='#6B91B1',label=f'{100*f:.0f}%') for f in [.1,.5,.9]]
ax.legend(handles=handles,title='Detected cells',ncol=3,loc='upper right',bbox_to_anchor=(1.0,-.25),frameon=False,title_fontsize=6)
co=primary.pivot(index='cell_type',columns='donor_id',values='n_cells').loc[ALLT,DONORS]
im=ay.imshow(np.log1p(co),cmap='Greys',aspect='auto',vmin=0,vmax=np.log1p(co.to_numpy().max()))
for i in range(len(ALLT)):
 for j in range(8):
  n=int(co.iloc[i,j]);ay.text(j,i,str(n),ha='center',va='center',fontsize=6.7,color='white' if n>100 else '#333333',weight='bold' if n>=30 else 'normal')
ay.set_xticks(range(8),[d+'\n'+('10X' if j<2 else 'BD') for j,d in enumerate(DONORS)],fontsize=6.5)
ay.set_yticks(range(len(ALLT)),[SHORT[t] for t in ALLT]);ay.axvline(4.5,color='#A35B55',lw=1.2)
ay.set_title('Donor coverage (bold: ≥30 cells; rightmost three donors: iNOA)',loc='left');letter(ay,'b')
save(fig,'Figure_S3_Annotation_Coverage',['annotation_display_means.tsv','cluster_labels.tsv','donor_celltype_coverage.tsv'])


(V/'00_admin/SC_figure_manifest_v3.json').write_text(json.dumps(manifest,indent=2))
print('COMPLETE Figure 5 and Figure S3')
