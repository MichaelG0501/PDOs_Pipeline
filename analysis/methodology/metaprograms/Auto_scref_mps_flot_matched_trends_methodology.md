# Methodology: scRef Finalised Metaprograms Matched FLOT PDO Response

## Overview
This methodology evaluates the expression and activity dynamics of the 17 finalized, centred-refined single-cell reference (scRef) metaprograms (MPs) across 4 matched pairs of FLOT chemotherapy-treated versus untreated patient-derived organoids (PDOs).

## Datasets and Inputs
1. **scRef Metaprograms (17 MPs)**:
   - Gene lists: `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`
   - Canonical grouping & annotations: `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`
   - Functional states represented: `Cell cycle` (MP1, MP5, MP13+), `Classic proliferation` (MP2+), `Squamous-to-intestinal` (MP14, MP3+, MP6+, MP11+, MP9+, MP10+), `Glandular-to-intestinal` (MP18b, MP16, MP17, MP8b, MP8+), `Stress-adaptive` (MP12), and `Cancer-cell immune mimicry` (MP15).
2. **Matched FLOT PDO Samples (8 samples across 4 patients)**:
   - Pre-treatment responder: `SUR1070` (Untreated vs Treated)
   - Post-treatment non-responder: `SUR1072` (Untreated vs Treated)
   - Post-treatment responder: `SUR1090` (Untreated vs Treated)
   - Pre-treatment non-responder: `SUR1181` (Untreated vs Treated)
   - Source data: `/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/PDOs_list_PDOs.rds`

## Scoring Methodology
- **Signature scoring**: Per-cell signature scores are computed using `UCell::ScoreSignatures_UCell(maxRank = 1500)`.
- **Sample pseudobulk summaries**: For each MP and sample, cell-level score distributions are summarized via sample mean, median, Q1 (25th percentile), Q3 (75th percentile), and cell counts.

## Paired Trend Analysis & Statistical Unit
- **Statistical unit**: Patient-level paired response across the 4 independent matched pairs.
- **Direction concordance**:
  - `increase`: Both mean and median delta (`Treated - Untreated`) > 0 in $\ge 3$ of 4 pairs.
  - `decrease`: Both mean and median delta (`Treated - Untreated`) < 0 in $\ge 3$ of 4 pairs.
  - `mixed`: Inconsistent response direction across pairs.
- **Hypothesis testing**:
  - Paired Wilcoxon signed-rank test comparing Treated vs Untreated sample means and medians.
  - Benjamini-Hochberg (BH) false discovery rate adjustment.

## Visualizations
- **Mean & Median Pair Trend Plots**: Individual panel per MP showing sample-level mean (solid line) and median (dashed line) trajectories across each patient pair, with paired Wilcoxon p-value annotations and visual spacing between patient pairs.
- **Sample Activity Boxplots**: Per-sample boxplots showing full single-cell score distributions.
- **Activity Heatmap**: Matrix of sample mean scores annotated by MP State, Direction, and Pair Support.
