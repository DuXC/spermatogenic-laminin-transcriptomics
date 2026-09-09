from pathlib import Path
from datetime import datetime,timezone
import pandas as pd,json,hashlib
S=Path(__file__).resolve().parents[1];O=S/'03_processed/single_cell';D=['N1','N2','N3','Cr1','Cr2','Cr3']
counts=[];samples=[];coverage=[];markers=[];qc=[]
age={'N1':31,'N2':33,'N3':55,'Cr1':39,'Cr2':25,'Cr3':36}
for d in D:
    receipt=json.loads((S/'00_admin'/(d+'_external_prepare_receipt.json')).read_text());assert receipt['published_cell_mapping_complete']
    counts.append(pd.read_csv(O/(d+'_pseudobulk_counts.tsv.gz'),sep='\t',index_col=0))
    x=pd.read_csv(O/(d+'_pseudobulk_samples.tsv'),sep='\t');x['age']=age[d];samples.append(x)
    coverage.append(pd.read_csv(O/(d+'_author_label_coverage.tsv'),sep='\t'))
    markers.append(pd.read_csv(O/(d+'_independent_marker_review.tsv'),sep='\t'))
    qc.append(pd.read_csv(O/(d+'_cell_QC.tsv.gz'),sep='\t'))
x=pd.concat(counts,axis=1);s=pd.concat(samples,ignore_index=True)
assert x.columns.tolist()==s.sample_id.tolist() and not x.isna().any().any()
x.to_csv(O/'all_pseudobulk_counts.tsv.gz',sep='\t',index_label='gene',compression='gzip')
s.to_csv(O/'all_pseudobulk_samples.tsv',sep='\t',index=False)
pd.concat(coverage,ignore_index=True).to_csv(S/'04_results/external_cell_coverage.tsv',sep='\t',index=False)
pd.concat(markers,ignore_index=True).to_csv(S/'04_results/external_independent_marker_review.tsv',sep='\t',index=False)
cells=pd.concat(qc,ignore_index=True);assert len(cells)==28666
cells.to_csv(O/'all_cell_QC.tsv.gz',sep='\t',index=False,compression='gzip')
label_coverage=pd.concat(coverage,ignore_index=True)
fractions=[]
for d in D:
    a=label_coverage[(label_coverage.donor==d)&(label_coverage.variant=='primary')].set_index('author_cell_type').n_cells
    total=int(a.sum());pmc=int(a['PMCs']);fib=int(a['Fibrotic peritubular myoid cells'])
    fractions.append({'donor':d,'group':'OA' if d.startswith('N') else 'Crypto','total_captured_cells':total,'author_PMC_cells':pmc,'author_fibrotic_PMC_cells':fib,
                      'peritubular_captured_cell_fraction':(pmc+fib)/total,'author_fibrotic_state_fraction_within_peritubular':fib/(pmc+fib)})
pd.DataFrame(fractions).to_csv(S/'04_results/peritubular_captured_state_fractions.tsv',sep='\t',index=False)
review=S/'00_admin/EXTERNAL_ANNOTATION_REVIEW_v1.md'
(S/'00_admin/external_annotation_freeze_receipt.json').write_text(json.dumps({'status':'REVIEWED_BEFORE_SC_ECM_INFERENCE','recorded_at_utc':datetime.now(timezone.utc).isoformat(),
 'review_sha256':hashlib.sha256(review.read_bytes()).hexdigest(),'cells':len(cells),'pseudobulk_units':len(s),'mapped_genes':len(x),
 'all_source_rows_retained':True,'broad_compartment_interpretation':'author-defined mixture of peritubular stromal states'},indent=2)+'\n')
print('JOIN_COMPLETE',x.shape,'cells',len(cells))
