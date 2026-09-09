# Spermatogenic laminin transcriptomics

Code and numerical source data for **Laminin related transcription across human spermatogenic histologies and adult peritubular cell states**.

Version v0.2.0 includes manuscript v4: four separately analyzed bulk cohorts (117 biopsies or study samples) and two adult single-cell cohorts (14 donors). New analyses add GSE9210, GSE108886 and GSE153947 to the original GSE4797, GSE145467 and GSE149512 analysis. Donors supply independent replication in cellular comparisons. Complete fixed programs, measured and unmeasured members, and negative results are retained.

## Authors

Xiancheng Du, Shuchun Tao, Yongkun Zhu, Chunhui Liu and Chao Sun. Department of Urology, Zhongda Hospital, School of Medicine, Southeast University, Nanjing, China. Chunhui Liu and Chao Sun are corresponding authors; Chao Sun is the last author.

Funding: China Postdoctoral Science Foundation, 2024M750457 (Chunhui Liu); Jiangsu Provincial Research Project on Traditional Chinese Medicine and Integrated Chinese-Western Medicine, ZXFZ2026021 (Chao Sun). The authors declare no competing interests.

## Obtain and verify the data

Clone this repository at tag `v0.2.0`, then run from its root:

```bash
python3 scripts/download_release_data.py
python3 scripts/verify_release_data.py
python3 scripts/download_external_data.py
python3 scripts/verify_external_data.py
```

Code v0.2.0 uses the original v0.1.0 data snapshot and a separate v0.2.0 external-validation asset. The original download script still retrieves the v0.1.0 snapshot. The release asset contains public bulk inputs and frozen annotation snapshots, normalized bulk expression, original-publication donor/feature mappings, per-cell annotations, donor pseudobulk counts, complete effects and sensitivities, and numerical figure source tables. File-level SHA-256 checksums are in `data_manifest.json`. GEO retains the original sequencing/count deposits. Recreate full cell matrices from the public eight-donor downloads when reprocessing them.

## Environment and reproduction

Python package versions for the cellular analysis are recorded in `07_single_cell/00_admin/python_requirements.lock.txt`. The original execution used R 4.6.0, limma 3.68.5 and edgeR 4.10.5; complete R session information is retained. Install R dependencies `limma`, `edgeR`, `statmod` and `jsonlite` using the compatible Bioconductor release. Figure generation additionally uses `pypdf`; Word authoring is not part of this code release. Arial or DejaVu Sans is used for plots, so appearance can vary with installed fonts.

```bash
bash scripts/run_reproduction.sh verify
bash scripts/run_reproduction.sh bulk
bash scripts/run_reproduction.sh pseudobulk
bash scripts/run_reproduction.sh figures
```

Set `SPERM_PYTHON` and `SPERM_RSCRIPT` to the desired interpreters. Reproduction commands overwrite generated results in the cloned working copy; use a fresh checkout/data extraction to retain the exact release bytes.

For full single-cell reprocessing, run `07_single_cell/scripts/00_acquire_sc.py`, `01_prepare_cluster.py`, and `01b_annotation_review.py` in that order. The cache lives beneath `.cache/single_cell`. Inspect donor/cluster identities before applying frozen labels with `02_annotate_aggregate.py`; cluster numbers alone are not transferable across a new clustering run. The released cell annotations provide the reference mapping. Whole-cell reprocessing is computationally larger than the donor-count comparison and may be sensitive to numerical libraries. It was not rerun during repository packaging.

## External validation reproduction

After both data downloads, use:

```bash
bash scripts/run_external_reproduction.sh prepare
bash scripts/run_external_reproduction.sh bulk
bash scripts/run_external_reproduction.sh pseudobulk
bash scripts/run_external_reproduction.sh verify
bash scripts/run_external_reproduction.sh figures
```

`prepare` reconstructs external bulk expression from the native GEO SOFT deposits and the frozen HGNC reference. `pseudobulk` starts from released donor-compartment counts. `aggregate` additionally reconstructs those counts from all six original Alevin matrices and the published cell annotation, preserving fractional estimates; this optional step requires more time and memory. `all` runs preparation, models, scientific checks and figures from the released pseudobulks. Python requires numpy, pandas, scipy and matplotlib; original full-cell annotation additionally uses the packages listed above. Use compatible R/Bioconductor dependencies and, where needed, set `R_LIBS_USER` to their installation directory.

The external entrypoints are rooted at `10_external_validation_20260910`. The native eligibility audit, timestamped local plans, coverage criteria and annotation review were frozen before new ECM inference. They are local prospective records, not claims of public preregistration. Source labels and fixed pathways are unchanged across sensitivities. A fresh data extraction is required before checksum verification if reproduction has overwritten generated outputs.

## Release verification

The portable R code was rerun from the frozen uploaded bulk expression and donor pseudobulk counts. All 44 generated numerical result tables matched the original frozen results within relative tolerance 1e-10 and absolute tolerance 1e-12. The bulk preparation entrypoint reproduced all 12 annotation tables and fixed memberships; the public single-cell validator passed 30 scientific checks on the frozen cell matrix. The complete figure entrypoint reproduced all 11 main/supplementary PNGs byte-for-byte and all 11 integrative source tables. These checks reproduce the original v3 results. For the v0.2.0 extension, public bulk preparation and both model entrypoints reproduced all 40 checked numerical tables at relative tolerance 1e-10 and absolute tolerance 1e-12. The portable scientific validator passed 117 checks. All three new PNGs reproduced byte-for-byte. The six raw single-cell aggregation runs were completed during analysis; they were not repeated solely for repository packaging.

## Scientific interpretation

Tissue laminin score directions recurred across four bulk cohorts, with important measurement differences. GSE9210 covers only 9 of 30 laminin members and does not meet the competitive-test minimum. GSE108886 covers 28 members; its estimated-correlation camera result is nonsignificant, and postmeiotic-marker adjustment attenuates the score association. Gene effects there are deposited normalized-signal differences because the exact logarithm base was not documented.

In GSE153947, estimated-correlation donor pseudobulk tests support an upward peritubular ECM shift in cryptozoospermia, including within the author PMC-state sensitivity. Age and mitochondrial sensitivities use the same six donors and are not independent replications. All six overlapping ECM programs moved upward. The author fibrotic-state label is retained provenance, not a new histological diagnosis, and captured-cell fractions are not tissue proportions.

The original GSE149512 iNOA interstitial comparison remains negative for all six competitive tests. It involves a different clinical context and population from the new peritubular analysis. Transcript patterns do not establish deposited matrix protein, causality, treatment response or a diagnostic model.

GSE145467 patient numbers in the paper could not be linked safely to expression-column numbering. Individual clinical covariates have not been assigned to those columns. Developmental samples and distinct etiologies are retained as separate metadata strata for GSE149512. Additional public atlases are not assumed to provide new independent donors.

## Navigation

- `scripts/`: bulk preparation, models and plotting.
- `07_single_cell/scripts/`: acquisition, annotation, donor aggregation and pseudobulk models.
- `08_manuscript_refinement_20260909/scripts/`: complete source tables and integrative figures.
- `08_manuscript_refinement_20260909/03_tables/`: the two main tables and Tables S1–S9 after data extraction.
- `08_manuscript_refinement_20260909/04_figures/`: six main and five supplementary figures in SVG/PDF/PNG after extraction.
- `10_external_validation_20260910/`: new cohort metadata, complete results, eight external supplemental source tables and three new figures after data extraction.
- `EXTERNAL_DATA_DICTIONARY.md`: new contrasts, units and multiplicity families.
- `DATA_DICTIONARY.md`: units, missingness, contrasts and multiplicity.
- `THIRD_PARTY_SOURCES.md`: original datasets, reference versions and rights.

## Citation and licensing

Cite this version using `CITATION.cff` and cite the original data-generating studies. Project code is MIT licensed; original project tables, figures and documentation are CC BY 4.0. Third-party data retain their original terms, listed in `THIRD_PARTY_SOURCES.md`. This repository release is a data/code deposit associated with a manuscript under author review.
