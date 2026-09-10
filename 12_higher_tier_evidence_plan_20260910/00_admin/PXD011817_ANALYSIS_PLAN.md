# PXD011817 secondary proteomic analysis

Frozen 2026-09-10 after inspecting the original supplementary spreadsheet headers and methods, before any new group-effect calculation. This is a prospectively specified extension of the current reanalysis, not a registered original protocol.

## Population and provenance

The 2019 study included ten independently sampled men with idiopathic NOA in its proteomics experiment, five with positive sperm retrieval and five with negative sperm retrieval. Three mass-spectrometry runs per ECM preparation are technical repeats. Original mmc5.xlsx explicitly labels ECM 1–5 positive and 6–10 negative. Verify all 105 proteins × 30 reported LFQ measurements against QUANT.xlsx before using that mapping. The source describes negative retrieval as SCO and positive retrieval as heterogeneous residual spermatogenesis. The contrast is retrieval-negative minus retrieval-positive within iNOA.

## Preparation

- Use the full deposited proteinGroups_elaborated sheet, not just the 105 selected matrisome proteins. Reconcile each row with the original proteinGroups sheet by its numeric id and exact protein-group identifiers.
- Respect MaxQuant LFQ normalization already applied by the depositors. Treat zero as an unmeasured abundance, log2 positive LFQ, and take the median of available technical-replicate log2 abundances per participant, requiring at least two of three measured repeats for that participant.
- Require at least two peptides, at least one unique peptide, and no reverse, contaminant or only-by-site flag in the source MaxQuant group. Map all supplied gene names through the project's frozen HGNC approved-symbol/previous-symbol/alias annotation. Keep a group for gene-level analysis only if all supplied names resolve to the same approved gene. This is a single-gene-annotation rule; it does not claim independent verification of every accession. Export all annotation exclusions.
- Where several eligible groups map to the same gene, choose the group with the highest median LFQ abundance across all ten participants, then lexicographic protein ids for an exact tie. This rule is fixed without reference to direction or statistical significance.
- Gene-level descriptive/inferential tables require at least three observed participants per clinical group. Primary competitive program testing and program scores use the same complete-case gene universe (a qualified measurement for all ten participants), so no protein imputation is required.

### Source-format clarification before effect calculation

The elaborated sheet modifies contaminant accession strings in 20 rows. Two of these rows retain the original contaminant strings in the paper's supplementary spreadsheet. Reconcile these rows by numeric source id, canonical accession-set containment and exact LFQ values for every technical measurement. Preserve both accession strings in the audit. Apply the original Contaminant/Reverse flags and exclude any original group containing a CON__ or REV__ accession even if its flag is blank. The complete 105-by-30 supplementary LFQ grid must still match exactly. This clarification addresses source formatting and is recorded before any inferential analysis.

## Statistical specification

- Keep all six frozen ECM programs and all thirty laminin-interaction members from 10_external_validation_20260910/01_metadata/fixed_gene_sets.json.
- A competitive program is evaluable when at least ten of its members occur in the complete-case gene universe. Retain smaller programs as NOT_EVALUABLE_COVERAGE. Count all six programs in the multiple-testing family, assigning p=1 to non-evaluable entries solely when computing six-program BH adjustment.
- Primary competitive analysis: limma camera with an intercept and retrieval-negative indicator; estimate within-set correlation (inter.gene.cor=NA). Retain the fixed-correlation 0.01 analysis as a prespecified sensitivity. No data-dependent batch correction or covariate assignment.
- Gene-level analysis: limma with robust empirical Bayes and trend, with BH adjustment across every tested gene. Export effects, confidence intervals, raw p values, adjusted p values and observed participant counts.
- Descriptive program scores: center and scale each complete-case gene across all ten participants, average the available fixed-program genes, and estimate the group difference with an ordinary two-group linear model and 95% CI. These scores describe coordinated abundance rather than competitive specificity. Adjust six score tests as one family.
- Laminin sensitivity: leave out each participant in turn; preserve the complete-case member list and original score scaling. Report every leave-one-out group estimate and estimated-correlation camera result.

## Interpretation

This is ECM extracted from tissue and associated with sperm-retrieval/histology strata. It does not measure transcription within peritubular cells, establish a causal role, or provide a diagnostic validation cohort. Individual clinical covariates and acquisition randomization are not supplied in the downloaded files. Retain these limitations and any discordance with RNA. PXD032722 and PXD023979 remain held for their separate batch/channel questions.
