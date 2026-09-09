# Independent public data validation scope

Frozen on 2026-09-10 before examining any new expression contrasts. This is a local prospective analysis record, not a public preregistration.

The existing v3.1 manuscript and its original data, six ECM programs, 30 laminin members, donor assignments and negative results remain the reference analysis. This extension asks whether its tissue-level laminin direction is reproducible in additional independent studies and whether additional adult single-cell studies can support comparable within-lineage tests.

## Cohort selection

Screen primary GEO records and original publications for adult human testicular tissue, original sample-to-donor mapping, pathology/etiology, assay and processing, study reuse and available measurements. Candidate bulk studies identified before new expression analysis are GSE45885, GSE9210, GSE45887 and GSE108886. GSE45885 and GSE45887 require explicit overlap assessment because they originate from the same research group; they cannot automatically be counted as two independent cohorts. Candidate single-cell studies include GSE154535, GSE202647, GSE106487, GSE157421 and additional original studies found by the documented native database search. GSE149512 is already analyzed and is not new external validation. Adult reference atlases are evaluated separately from disease cohorts.

The eligibility decision is based on sample identity and metadata, never on significance or direction. All screened candidates and reasons are retained. Bulk independent tests require at least three independent donors in both comparison groups and no complete disease-by-assay or disease-by-study confounding. Single-cell tests require at least three independent adult donors per group, at least 30 confidently assigned cells per donor in that lineage, comparable assays, and no complete disease-by-assay or disease-by-study confounding. Pooled libraries require explicit donor resolution. This minimum permits an exploratory donor-level test; it does not establish adequate power. If a study fails a gate, retain a descriptive eligibility result without loosening the threshold.

## Fixed questions and analysis family

The laminin program is the focal external endpoint, with all six original ECM programs reported as one multiplicity family for each prespecified cohort contrast. Severe histology versus normal spermatogenesis is primary where the original individual labels support it; heterogeneous NOA versus OA is a distinct clinical comparison. Histological categories are not converted into unsupported clinical diagnoses. The same original laminin members and exploratory cell-stage marker sets are retained. No hub-gene selection or new classifier is part of this extension.

Before model fitting, write and hash a cohort-specific plan documenting exact donor lists, contrasts, measurement scale, probe mapping, missing-value rules, multiplicity families and sensitivity analyses. Processing is guided by original documentation and data structure; it must not be selected according to pathway results. Studies are modeled separately. Report score effects and 95% confidence intervals, all member-gene effects and coverage, correlation-estimated competitive enrichment, fixed-correlation sensitivity, and prespecified sample influence checks. No data-dependent sample removal is permitted.

## Composition and interpretation

Any reference-based composition analysis requires an adult independent reference with defensible cell identities and donor-aware held-out or synthetic-mixture validation before patient deconvolution. If this requirement cannot be met, do not describe a marker score as a measured cell proportion. Existing marker adjustment remains a sensitivity analysis with its original limitations. Public spatial or protein data require independent participants and appropriate disease/tissue context before being described as orthogonal validation.

A successful bulk result strengthens tissue-level replication. It does not establish within-cell-type activation or causal ECM remodeling. A failed or unstable external result is reported with the same priority as support. New hospital data or tissue experiments are not initiated in this computational extension.
