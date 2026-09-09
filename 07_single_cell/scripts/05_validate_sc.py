#!/usr/bin/env python3
"""Scientific data-integrity and numerical checks for the adult SC release."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import stats
S=Path(__file__).resolve().parents[1];P=S.parent
checks=[]
def check(name,passed,details=None):
 checks.append(dict(check=name,passed=bool(passed),details=details))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 return h.hexdigest()
a=ad.read_h5ad(S/'03_processed/adult_annotated.h5ad')
cells=pd.read_csv(S/'02_annotation/cells_final.tsv.gz',sep='\t',index_col=0)
check('Final cell labels match archived matrix',a.obs_names.equals(cells.index) and (a.obs.cell_type.astype(str).to_numpy()==cells.cell_type.to_numpy()).all())
check('Unique adult cells and approved genes',a.obs_names.is_unique and a.var_names.is_unique and a.shape==(47003,28529),list(a.shape))
check('Only eight authorized adult donors',set(a.obs.donor_id)=={'LZ003','LZ007','LZ013','LZ014','LZ015','LZ017','LZ018','LZ019'})
check('Count layer remains nonnegative integers',np.issubdtype(a.layers['counts'].dtype,np.integer) and (a.layers['counts'].data>=0).all())
check('Primary original-feature filters hold',((cells.n_features_original>500)&(cells.n_features_original<9000)&(cells.total_counts_original<80000)&(cells.pct_mito_original<40)).all())
check('Count mapping does not add molecules',(cells.total_counts_mapped<=cells.total_counts_original).all())
receipts=json.loads((S/'00_admin/download_receipt.json').read_text())
check('Eight original matrix hashes verified',len(receipts)==8 and all(sha(S/'01_raw'/r['file'])==r['sha256'] for r in receipts))
original=pd.concat([pd.read_csv(S/'02_annotation'/f'{r["donor_id"]}_all_cell_QC.tsv.gz',sep='\t') for r in receipts])
check('Raw and QC cell denominators reconcile',len(original)==48815 and int(original.primary_qc_pass.sum())==a.n_obs)
check('Doublet detection completed for every donor',cells.scrublet_status.eq('complete').all(),int(cells.predicted_doublet.sum()))
sets=json.loads((P/'02_annotation/fixed_gene_sets.json').read_text())['primary'];ecm=set().union(*map(set,sets.values()))
panels=json.loads((S/'02_annotation/annotation_marker_panels.json').read_text())
check('Annotation positive markers independent of ECM',not ecm.intersection(set().union(*map(set,panels.values()))))
check('Clustering HVGs independent of ECM',not ecm.intersection(a.var_names[a.var.highly_variable]))
check('Annotation control pools corrected before assignment',json.loads((S/'00_admin/preannotation_QA.json').read_text())['ECM_genes_excluded_from_positive_markers_and_control_pool'])
cov=pd.read_csv(S/'04_results/donor_celltype_coverage.tsv',sep='\t')
tot=cov[cov.variant.eq('primary')].groupby('donor_id').n_cells.sum()
check('Cell-type census reconciles for every donor',tot.equals(cells.groupby('donor_id').size().rename('n_cells')))
check('Capture fractions sum to one',np.allclose(cov.groupby(['variant','donor_id']).capture_fraction.sum(),1))
src=pd.read_csv(S/'04_results/ECM_all_genes_donor_sources.tsv.gz',sep='\t')
check('All measured ECM union genes exported',set(src.gene)==ecm and len(ecm)==364)
check('Expression and detection bounds valid',(src.mean_CP10K>=0).all() and src.fraction_detected.between(0,1).all())
check('Source eligibility follows fixed coverage',(src.source_eligible==((src.n_cells>=30)&src.annotation_confidence.eq('high'))).all())
summ=pd.read_csv(S/'04_results/reference_ECM_gene_source_summary.tsv',sep='\t')
ref=src[src.variant.eq('primary')&src.group.eq('OA')&src.source_eligible]
recalc=ref.groupby(['cell_type','gene']).mean_CP10K.mean()
got=summ.set_index(['cell_type','gene']).mean_CP10K.reindex(recalc.index)
check('Reference means give donors equal weight',np.allclose(got,recalc))
coverage=pd.read_csv(S/'04_results/ECM_program_gene_coverage.tsv',sep='\t')
lam=coverage[coverage.pathway.eq('REACTOME_LAMININ_INTERACTIONS')]
check('All 30 laminin genes measured and detected in references',len(lam)==30 and lam.measured.all() and lam.reference_cells_detected.gt(0).all())
# Rebuild one entire donor/type pseudobulk directly from the archived counts.
pb=pd.read_csv(S/'03_processed/pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
ps=pd.read_csv(S/'03_processed/pseudobulk_samples.tsv',sep='\t')
check('Pseudobulk sample order explicit',list(pb.columns)==ps.sample_id.tolist())
unit='primary__LZ015__Interstitial_stromal'
mask=a.obs.donor_id.eq('LZ015')&a.obs.cell_type.eq('Interstitial_stromal')
check('Independent recount of smallest reference interstitial unit',np.array_equal(np.asarray(a.layers['counts'][mask.to_numpy()].sum(axis=0)).ravel(),pb[unit].reindex(a.var_names).to_numpy()),int(mask.sum()))
elig=pd.read_csv(S/'04_results/pseudobulk_eligibility.tsv',sep='\t')
check('Only eligible broad interstitial comparisons performed',set(elig[elig.eligible].cell_type)=={'Interstitial_stromal'} and len(elig[elig.eligible])==3)
e=ps[ps.variant.eq('primary')&ps.cell_type.eq('Interstitial_stromal')&ps.technology.eq('BD_Rhapsody')]
check('Same-platform pseudobulk uses 3 versus 3 donors',len(e)==6 and e.n_cells.ge(30).all() and e.groupby('group').size().to_dict()=={'OA':3,'iNOA':3})
cam=pd.read_csv(S/'04_results/pseudobulk_camera.tsv',sep='\t')
check('All six fixed sets retained in every pseudobulk sensitivity',cam.groupby('variant').pathway.nunique().eq(6).all() and len(cam)==18)
for variant,d in cam.groupby('variant'):
 check('BH family '+variant,np.allclose(stats.false_discovery_control(d.PValue),d.FDR_all_types_six_sets))
g=pd.read_csv(S/'04_results/pseudobulk_all_genes.tsv.gz',sep='\t')
check('Pseudobulk confidence intervals contain point effects',((g['CI.L']<=g.logFC)&(g.logFC<=g['CI.R'])).all())
base=json.loads((P/'00_admin/release_manifest.json').read_text())['files']
bad=[r['file'] for r in base if not (P/r['file']).exists() or sha(P/r['file'])!=r['sha256']]
check('Frozen first bulk release unchanged',not bad,{'files_checked':len(base),'differences':bad})
clinical=pd.read_csv(S/'02_annotation/GSE149512_donor_metadata_v2.tsv',sep='\t').set_index('donor_id')
check('High-resolution clinical values propagated',clinical.loc['LZ011','published_total_cells']==7896 and clinical.loc['LZ014','published_Sertoli_cells']==27)
manifest=json.loads((S/'00_admin/fulltext_manifest.json').read_text())
receipts=list((S/'00_admin/fulltext_receipts').glob('*Commit*.json'))
check('Central fulltext commit receipt retained',len(receipts)==1 and len(manifest['downloaded'])==2)
result={'status':'PASS' if all(x['passed'] for x in checks) else 'FAIL','checks':checks,'summary':{'raw_cells':len(original),'retained_cells':a.n_obs,'reference_cells':int(cells.group.eq('OA').sum()),'iNOA_cells':int(cells.group.eq('iNOA').sum()),'mapped_genes':a.n_vars,'ECM_union_genes':len(ecm),'predicted_doublets':int(cells.predicted_doublet.sum()),'mt20_cells':int(cells.mt20_pass.sum()),'median_count_retention':float((cells.total_counts_mapped/cells.total_counts_original).median()),'confident_population_categories':int(cells[cells.annotation_confidence.eq('high')].cell_type.nunique()),'provisional_population_categories':int(cells[~cells.annotation_confidence.eq('high')].cell_type.nunique()),'pseudobulk_eligible_population':'Interstitial_stromal','primary_camera_FDR':float(cam[cam.variant.eq('primary')].FDR_all_types_six_sets.iloc[0])}}
(S/'00_admin/numerical_QA.json').write_text(json.dumps(result,indent=2))
print(result['status'],len(checks),'checks',flush=True)
for c in checks:
 if not c['passed']:print('FAILED',c,flush=True)
assert result['status']=='PASS'
