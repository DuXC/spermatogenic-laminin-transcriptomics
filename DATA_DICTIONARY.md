# Data dictionary

All scientific text tables use UTF-8 and tab separators; `.gz` files are gzip compressed. Blank/NA estimates denote unavailable measurement or failed predefined expression/coverage eligibility, not zero expression.

| Quantity | Meaning |
|---|---|
| Bulk `logFC` | Separate-cohort log2-scale gene effects; JS2 minus JS10 or impaired minus normal |
| `CI.L`, `CI.R` | Moderated 95% gene-effect limits; program-score ordinary-model intervals have their own columns |
| `FDR_primary_family` | GSE4797: all genes over three primary contrasts; GSE145467: its tested gene family |
| `adj.P.Val` in pseudobulk | BH adjustment over all 13,130 tested interstitial genes |
| Competitive camera FDR | Six programs across three GSE4797 primary contrasts (18 tests), six GSE145467 tests; six eligible interstitial tests within each QC analysis |
| `mean_CP10K` | Cell-normalized expression aggregated within donor/type, then equally weighted across eligible reference donors |
| `fraction_detected` | Nonzero original-count cells divided by the applicable cell denominator |
| `n_donors` | Biological donor coverage, not independent cells |
| `source_eligible` | Confident population label and at least 30 cells in that donor/population |
| Capture fraction | Observed recovered cells; not an unbiased histological cell fraction |
| Program source score | Reference-standardized relative transcriptional source score |

Source populations are not mutually exclusive gene-expression compartments. Figure 6 uses gene-wise scaled log1p equal-donor mean expression among populations with at least three eligible reference donors. Complete unscaled source values and coverage are in Table S4. All 30 laminin members are retained; 29 were shared by the bulk platforms and 25 passed the interstitial expression filter.

Primary SC QC uses >500 and <9000 detected features, <80000 counts and <40% mitochondrial counts. Sensitivities restrict mitochondrial counts to <20% or remove predicted doublets. None of the primary thresholds was relaxed for the repository release.
