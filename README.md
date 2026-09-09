# Spermatogenic laminin transcriptomics

Code and numerical source data for **Histological patterns and adult cellular sources of laminin related transcription in human spermatogenic dysfunction**.

Code version v0.1.1 preserves the scientific results used in manuscript v3. The work analyzes tissue-level patterns in GSE4797 and GSE145467 and adult cell sources in GSE149512. Donors supply the independent replication in cellular comparisons. Results include complete gene families and nonsignificant sensitivity analyses.

## Authors

Xiancheng Du, Shuchun Tao, Yongkun Zhu, Chunhui Liu and Chao Sun. Department of Urology, Zhongda Hospital, School of Medicine, Southeast University, Nanjing, China. Chunhui Liu and Chao Sun are corresponding authors; Chao Sun is the last author.

Funding: China Postdoctoral Science Foundation, 2024M750457 (Chunhui Liu); Jiangsu Provincial Research Project on Traditional Chinese Medicine and Integrated Chinese-Western Medicine, ZXFZ2026021 (Chao Sun). The authors declare no competing interests.

## Obtain and verify the data

Clone this repository at tag `v0.1.1`, then run from its root:

```bash
python3 scripts/download_release_data.py
python3 scripts/verify_release_data.py
```

Code v0.1.1 uses the versioned v0.1.0 data snapshot. The release asset contains public bulk inputs and frozen annotation snapshots, normalized bulk expression, original-publication donor/feature mappings, per-cell annotations, donor pseudobulk counts, complete effects and sensitivities, and numerical figure source tables. File-level SHA-256 checksums are in `data_manifest.json`. GEO retains the original sequencing/count deposits. Recreate full cell matrices from the public eight-donor downloads when reprocessing them.

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

## Release verification

The portable R code was rerun from the frozen uploaded bulk expression and donor pseudobulk counts. All 44 generated numerical result tables matched the original frozen results within relative tolerance 1e-10 and absolute tolerance 1e-12. The bulk preparation entrypoint reproduced all 12 annotation tables and fixed memberships; the public single-cell validator passed 30 scientific checks on the frozen cell matrix. The complete figure entrypoint reproduced all 11 main/supplementary PNGs byte-for-byte and all 11 integrative source tables. These checks reproduce existing results.

## Scientific interpretation

The severe-histology direction of the tissue laminin score was reproducible between bulk cohorts. Competitive enrichment was sensitive to residual inter-gene correlation. Only broad interstitial stromal cells met the same-platform three-versus-three coverage criterion, and none of six ECM programs was competitively enriched there. Cellular sources are equal-donor descriptive expression summaries, not matrix protein measurements or validated cell proportions. These findings do not establish within-lineage ECM activation, causality, treatment response or a diagnostic model.

GSE145467 patient numbers in the paper could not be linked safely to expression-column numbering. Individual clinical covariates have not been assigned to those columns. Developmental samples and distinct etiologies are retained as separate metadata strata for GSE149512. Additional public atlases are not assumed to provide new independent donors.

## Navigation

- `scripts/`: bulk preparation, models and plotting.
- `07_single_cell/scripts/`: acquisition, annotation, donor aggregation and pseudobulk models.
- `08_manuscript_refinement_20260909/scripts/`: complete source tables and integrative figures.
- `08_manuscript_refinement_20260909/03_tables/`: the two main tables and Tables S1–S9 after data extraction.
- `08_manuscript_refinement_20260909/04_figures/`: six main and five supplementary figures in SVG/PDF/PNG after extraction.
- `DATA_DICTIONARY.md`: units, missingness, contrasts and multiplicity.
- `THIRD_PARTY_SOURCES.md`: original datasets, reference versions and rights.

## Citation and licensing

Cite this version using `CITATION.cff` and cite the original data-generating studies. Project code is MIT licensed; original project tables, figures and documentation are CC BY 4.0. Third-party data retain their original terms, listed in `THIRD_PARTY_SOURCES.md`. This repository release is a data/code deposit associated with a manuscript under author review.
