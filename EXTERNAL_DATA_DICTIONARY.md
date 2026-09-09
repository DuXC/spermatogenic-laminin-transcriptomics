# External validation data dictionary

The v0.2.0 asset adds `10_external_validation_20260910`; the v0.1.0 base asset supplies the original frozen HGNC reference in `01_raw`. All source TSVs are tab-separated, UTF-8. NA is missing or untested, never zero. Donor and GSM identifiers are local to their original study.

| Layer | Unit and interpretation |
|---|---|
| GSE9210 expression | Log2 common-reference ratio, converted from deposited natural-log ratios |
| GSE108886 expression | Deposited quantile-normalized signal; exact logarithm base not documented |
| Bulk program scores | Mean within-cohort member-gene z score; member coverage differs by platform |
| Bulk contrast | NOA minus OA, at individual-biopsy level |
| Single-cell original counts | Fractional Alevin estimates; retained without rounding |
| Single-cell pseudobulk | Sum within each donor and prespecified compartment or state |
| Single-cell contrast | Cryptozoospermia minus OA, three donors per group |
| Cell-source `mean_CPTT` | Mean counts per ten thousand; source files retain equal-donor inputs |
| Cell-source `mean_log1p_CPTT` | Mean of cell-level log1p CP10K, used only for descriptive annotation review |
| `n_cells` and captured fractions | Sequenced-cell recovery, not unbiased tissue proportions |
| Gene confidence intervals | Moderated limma 95% intervals, with the effect unit stated in each table |
| `mean_member_log2FC` | Descriptive arithmetic mean of member effects; no pathway CI |

Bulk `FDR_six` is BH across six fixed programs in a cohort; `FDR_external_12` combines both external cohorts. Competitive and score tests are separate families. Estimated nonnegative camera correlation is primary in this extension; fixed 0.01 is sensitivity. Competitive tests require ten measured members; program scores require three. GSE9210 laminin competitive tests are not evaluable at 9/30 coverage.

Single-cell `FDR_family` covers all six programs in each primary or explicitly named sensitivity family. The eligible primary somatic comparison is the author-defined peritubular mixture. Author PMC state, donor age and mitochondrial <20% analyses are prespecified sensitivities sharing donors. `adj.P.Val` is BH across all tested genes within a model. All 30 laminin members remain in every corresponding table, even when filtered from inference.

`08_tables` retains main Tables 1–3 and Tables S1–S17; S1–S9 preserve the original analysis. `09_references` supplies the updated 26-record RIS and BibTeX bibliography. `external_data_manifest.json` at repository root records SHA-256 for every external data asset file. `00_admin` records frozen analysis and annotation decisions plus execution environments. Publication-source full texts and private author manuscripts are not distributed in the data asset.
