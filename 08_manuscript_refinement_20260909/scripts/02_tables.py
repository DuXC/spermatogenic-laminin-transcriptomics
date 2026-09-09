#!/usr/bin/env python3
"""Assemble review tables from frozen numerical evidence without refitting models."""
from pathlib import Path
import json,shutil
import numpy as np,pandas as pd
S=Path(__file__).resolve().parents[1];P=S.parent;SC=P/'07_single_cell';O=S/'03_tables'
sets=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary']
order=['GOCC_BASEMENT_MEMBRANE','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_COLLAGEN_FORMATION','REACTOME_LAMININ_INTERACTIONS','REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
labels=['Basement membrane','ECM organization','Collagen formation','Laminin interactions','Integrin interactions','ECM degradation']
groups={'Laminin chains':'LAMA1 LAMA2 LAMA3 LAMA4 LAMA5 LAMB1 LAMB2 LAMB3 LAMC1 LAMC2 LAMC3'.split(),'Collagens':'COL4A1 COL4A2 COL4A3 COL4A4 COL4A5 COL4A6 COL7A1 COL18A1'.split(),'Matrix partners':'HSPG2 NID1 NID2'.split(),'Integrin receptors':'ITGA1 ITGA2 ITGA3 ITGA6 ITGA7 ITGAV ITGB1 ITGB4'.split()}
lam=[g for genes in groups.values() for g in genes];assert set(lam)==set(sets['REACTOME_LAMININ_INTERACTIONS']) and len(lam)==30
rows=[];bulk={}
for cohort,contrast in [('GSE4797','JS2_vs_JS10'),('GSE145467','impaired_vs_normal')]:
    raw=pd.read_csv(P/f'04_results/primary/{cohort}_limma_all_contrasts.tsv.gz',sep='\t')
    bulk[cohort]=raw[raw.contrast.eq(contrast)].set_index('gene')
    mod=pd.read_csv(P/f'04_results/primary/{cohort}_module_effects.tsv',sep='\t').query('contrast==@contrast').set_index('set')
    fixed=pd.read_csv(P/f'04_results/primary/{cohort}_camera.tsv',sep='\t').query('contrast==@contrast').set_index('set')
    est=pd.read_csv(P/f'04_results/sensitivity/{cohort}_camera_estimated_correlation.tsv',sep='\t').query('contrast==@contrast').set_index('set')
    for key,label in zip(order,labels):
        a,b,c=mod.loc[key],fixed.loc[key],est.loc[key]
        rows.append(dict(cohort=cohort,contrast=contrast,pathway=key,label=label,fixed_members=len(sets[key]),measured_members=int(b.NGenes),effect=a.effect,CI_low=a.CI_low,CI_high=a.CI_high,module_P=a.P,module_FDR=a.FDR_family,camera_fixed_rho=.01,camera_fixed_FDR=b.FDR_family,camera_estimated_rho=c.Correlation,camera_estimated_FDR=c.FDR_family,camera_test_family=18 if cohort=='GSE4797' else 6))
pd.DataFrame(rows).to_csv(O/'Table_2_program_effects_and_enrichment.tsv',sep='\t',index=False)
sources=pd.read_csv(SC/'04_results/reference_ECM_gene_source_summary.tsv',sep='\t')
pb=pd.read_csv(SC/'04_results/pseudobulk_all_genes.tsv.gz',sep='\t');pb=pb[pb.variant.eq('primary')].set_index('gene')
rows=[]
for gene in lam:
    r=dict(gene=gene,molecular_group=next(k for k,v in groups.items() if gene in v))
    for cohort,b in bulk.items():
        r[cohort+'_measured']=gene in b.index
        for col in ['logFC','CI.L','CI.R','P.Value','FDR_primary_family']:
            r[cohort+'_'+col]=b.loc[gene,col] if gene in b.index else np.nan
    d=sources[sources.gene.eq(gene)&sources.n_donors.ge(3)]
    t=d.loc[d.mean_CP10K.idxmax()]
    r.update(max_mean_source_at_least_3_donors=t.cell_type,max_source_n_donors=int(t.n_donors),max_source_CP10K=t.mean_CP10K,max_source_detection_fraction=t.fraction_detected)
    r['interstitial_tested']=gene in pb.index
    for col in ['logFC','CI.L','CI.R','P.Value','adj.P.Val']:
        r['interstitial_'+col]=pb.loc[gene,col] if gene in pb.index else np.nan
    r['same_bulk_direction']=np.sign(r['GSE4797_logFC'])==np.sign(r['GSE145467_logFC']) if all(r[c+'_measured'] for c in bulk) else np.nan
    rows.append(r)
pd.DataFrame(rows).to_csv(O/'Table_S3_all_30_laminin_gene_evidence.tsv',sep='\t',index=False,na_rep='NA')
sources[sources.gene.isin(lam)].to_csv(O/'Table_S4_all_30_laminin_gene_sources.tsv',sep='\t',index=False)
universe=set().union(*(set(v) for v in sets.values()))
members=pd.DataFrame([dict(gene=g,**{k:g in sets[k] for k in order},measured_GSE4797=g in bulk['GSE4797'].index,measured_GSE145467=g in bulk['GSE145467'].index) for g in sorted(universe)])
members.to_csv(O/'Table_S5_complete_program_memberships.tsv',sep='\t',index=False)
overlap=pd.DataFrame([[len(set(sets[a])&set(sets[b]))/len(set(sets[a])|set(sets[b])) for b in order] for a in order],index=order,columns=order)
overlap.to_csv(O/'Table_S6_program_Jaccard_overlap.tsv',sep='\t')
for src,dst in [('GSE149512_donor_metadata_v2.tsv','Table_S1_original_donor_metadata.tsv')]:shutil.copy2(SC/'02_annotation'/src,O/dst)
shutil.copy2(SC/'04_results/donor_celltype_coverage.tsv',O/'Table_S2_donor_cell_coverage.tsv')
shutil.copy2(SC/'04_results/ECM_program_gene_coverage.tsv',O/'Table_S7_single_cell_program_coverage.tsv')
shutil.copy2(SC/'04_results/pseudobulk_camera.tsv',O/'Table_S8_interstitial_pathway_sensitivities.tsv')
shutil.copy2(SC/'04_results/ECM_program_donor_source_scores.tsv',O/'Table_S9_donor_program_source_scores.tsv')
cohorts=[dict(cohort='GSE4797',role='Primary histological contrasts',unit='28 men',groups='JS10 12; JS8 6; JS5 5; JS2 5',platform='CodeLink GPL2891',measured_genes=12931,scope='Cross-sectional biopsy categories'),dict(cohort='GSE145467',role='Independent bulk directional support',unit='20 samples',groups='Normal 10; impaired 10',platform='Agilent GPL4133; two-channel common reference',measured_genes=18884,scope='Broad groups; patient-level clinical linkage unresolved'),dict(cohort='GSE149512 adult subset',role='Adult source reference and eligible pseudobulk',unit='8 donors; 47003 cells',groups='OA 5; iNOA 3',platform='10X 2 OA; BD Rhapsody 3 OA and 3 iNOA',measured_genes=28529,scope='Only broad interstitial cells eligible for BD 3 vs 3')]
pd.DataFrame(cohorts).to_csv(O/'Table_1_analyzed_cohorts.tsv',sep='\t',index=False)
qa={'fixed_programs':len(sets),'union_genes':len(universe),'laminin_genes':len(lam),'laminin_both_bulk_measured':int(sum(all(r[c+'_measured'] for c in bulk) for r in rows)),'laminin_pseudobulk_tested':int(sum(r['interstitial_tested'] for r in rows)),'integration_is_descriptive':True,'model_refits':0}
(S/'00_admin/table_build_QA.json').write_text(json.dumps(qa,indent=2))
(S/'00_admin/laminin_order.json').write_text(json.dumps(groups,indent=2))
print(json.dumps(qa,indent=2))
