# Protein and common-member data dictionary

Version 0.3.0 accompanies manuscript v5. Paths retain their analysis-stage directory names. No new participant enrollment or clinical prediction is represented.

## PXD011817 isolated extracellular matrix

Primary source: Alfano et al., Fertility and Sterility 2019, [10.1016/j.fertnstert.2018.12.002](https://doi.org/10.1016/j.fertnstert.2018.12.002); [PRIDE PXD011817](https://www.ebi.ac.uk/pride/archive/projects/PXD011817).

The folder `12_higher_tier_evidence_plan_20260910/03_analysis/PXD011817` contains ten participant columns. ECM1–5 are retrieval-positive iNOA; ECM6–10 are retrieval-negative iNOA. Each originated from three technical measurements. Biological replication is ten men, not thirty measurements. Effects are retrieval-negative minus retrieval-positive. They do not estimate a NOA–OA contrast.

- `gene_log2_LFQ_min3_each_group.tsv`: 412 tested genes, positive LFQ values log2 transformed and aggregated as the participant median when at least two technical measurements were observed. Each retained gene has at least three observed participants per group. NA means unmeasured.
- `gene_log2_LFQ_complete_10.tsv`: 300 genes measured in every participant; background for competitive tests and scores.
- `participants.tsv`, `technical_run_to_participant.tsv`: original study labels and the verified replicate/group mapping.
- `all_protein_group_annotation_audit.tsv`: original and elaborated protein identifiers, quality flags, HGNC gene mapping and the outcome-independent duplicate selection rule.
- `supplement_to_PRIDE_value_check.tsv`: 105 supplementary groups and 3150 LFQ values reconciled exactly to the deposit.
- `fixed_program_coverage.tsv`: all six fixed programs and their complete-case measured members. Laminin interactions have 15 of 30 members.
- `all_gene_results.tsv`: robust trend-adjusted limma effects, moderated 95% intervals and BH FDR across 412 genes.
- `all_30_laminin_members_results.tsv`: all fixed members, including explicit unevaluable status.
- `camera_estimated_correlation.tsv`, `camera_fixed_correlation_001.tsv`: primary estimated-correlation and fixed-0.01 sensitivity results, BH family of six programs.
- `participant_program_scores.tsv`, `program_score_group_effects.tsv`: within-cohort mean gene z scores, OLS group differences, unadjusted 95% t intervals and BH family of six score tests.
- `laminin_leave_one_participant_out.tsv`: every participant omission. The full-cohort standardization is retained for score influence analysis; CAMERA is refitted per omission.

The original workbooks can be obtained lawfully from their deposit or publisher and placed under `12_higher_tier_evidence_plan_20260910/02_acquired`:

| Local filename | Source | Expected SHA-256 |
|---|---|---|
| PXD011817_QUANT.xlsx | [PRIDE QUANT.xlsx](https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/01/PXD011817/QUANT.xlsx) | 7d2e64dc6cbf8a231d49fdb1bd249c6da849d242efd5b4c5f15852a4785a98b1 |
| Alfano2019_mmc5.xlsx | [Publisher Table S5](https://ars.els-cdn.com/content/image/1-s2.0-S0015028218322519-mmc5.xlsx) | 68a28be6584f875929d728560e8c8af0095832c88c56b1925eec49bc2bd189f2 |

After the v0.1.0 and v0.2.0 data downloads, run `python3 12_higher_tier_evidence_plan_20260910/05_scripts/04_prepare_ecm_proteomics.py` to rebuild the participant inputs. The script reads the frozen HGNC snapshot and fixed gene sets from the earlier assets. Python requires pandas, numpy and openpyxl. Inference requires R, limma and statmod; figures require scipy and matplotlib. Publisher access controls may require an ordinary authorized browser download.

## Four-cohort common-member sensitivity

`13_public_data_submission_20260910/02_sources/common_member` contains all 24 tests (four cohorts by six fixed programs), every original set member with shared-coverage status, all participant scores and all 32 laminin gene omissions. The plan explicitly records that this is post hoc; it was frozen before these new calculations, after the original full-member effects were known.

The common laminin-interaction subset has eight members: five collagen genes (COL4A1, COL4A2, COL4A5, COL4A6, COL7A1), ITGA3, ITGA6 and LAMC3. It is not a set of eight laminin chains. Group designs remain four-level Johnsen categories for GSE4797, impaired/normal for GSE145467 and NOA/OA for each external cohort. The displayed GSE4797 contrast is JS2 minus JS10. Genes are z standardized within each full cohort, using the sample standard deviation.

The primary sensitivity family comprises 24 score tests. Gene omissions are influence analyses with all results retained; their range of point estimates is not a confidence interval. Contrasts are not pooled. S23–S25 contain score tests, gene omissions and complete coverage; sample scores are available in the repository. These findings do not replace full-program reporting or competitive CAMERA tests.

## Source tables and candidate-resource audit

The extension's `04_supplement/source_tables` folder contains updated Table 1 and S18–S25. Earlier S1–S17 remain in the prior data assets. Table S22 records unresolved eligibility of PXD023979 (TMT channel mapping) and PXD032722 (group-separated date strata and unresolved acquisition comparability). No new group effects are reported for these two datasets.

Di Persio et al. 2021 Table S1, [10.1016/j.xcrm.2021.100395](https://doi.org/10.1016/j.xcrm.2021.100395), identifies PXD023979 participants N1–N5 and Cr2–Cr5. Five (N1, N2, N3, Cr2, Cr3) are shared with GSE153947. Assay overlap is retained in the lineage audit rather than counted as independent donor replication.
