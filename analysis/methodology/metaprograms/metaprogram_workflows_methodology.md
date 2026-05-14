# PDO Metaprogram Workflow Methodology

This document describes scripts under `analysis/metaprograms/`.

## 1. Aim

The metaprogram workflow identifies robust malignant PDO transcriptional
programs with GeneNMF, selects an optimal number of metaprograms, scores cells
with UCell, and compares PDO metaprograms with scATLAS and 3CA references.

The current chosen optimal PDO result is:

- `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`

## 2. Core Inputs

- `PDOs_outs/geneNMF_outs.rds`
- `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_{k}.rds`
- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/UCell_scores_filtered.rds`
- scRef metaprogram objects under `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/`
- 3CA MP reference `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`

## 3. MP Filtering

Downstream PDO scripts must apply both filters before interpreting MPs:

1. Silhouette filter: remove MPs with silhouette `< 0`.
2. Sample-coverage filter: remove MPs with sample coverage `< 0.25`.

This removes sparse MPs from downstream state annotation and correlation
analyses.

## 4. Script Roles

### `find_optimal_nmf.R`

Reads available `geneNMF_metaprograms_nMP_{k}.rds` files, computes silhouette
and WSS diagnostics, applies kneedle-style inflection detection, and reports the
optimal PDO nMP. Current selected nMP is 13.

### `extend_nMP_range.R`

Extends `getMetaPrograms()` over additional nMP values from the raw
`geneNMF_outs.rds` object. Run in the `gnmf` environment.

### `update_optimal_mp.R`

Updates the default optimal metaprogram object after nMP selection. Run in the
`gnmf` environment.

### `mp_ucell_scoring.R`

Scores PDO cells for retained MPs. These scores are consumed by state definition
and downstream MP analyses.

### Correlation and comparison scripts

`PDO_mp_correlation_crossdata.R`, `mp_correlation_pdo.R`,
`Auto_3CA_pseudobulk_correlation_crossdata.R`, and related scripts are terminal
comparison analyses unless explicitly documented otherwise. Their outputs should
not replace the canonical nMP=13 object.

## 5. Output Tiers And Replotting

Long-running metaprogram analyses should cache:

- raw GeneNMF or `getMetaPrograms()` objects in `intermediate/`
- score matrices in `intermediate/`
- correlation or enrichment tables in `tables/`
- heatmaps/scatters in `figures/`
- run summaries in `logs/`

Replot-only changes should not rerun GeneNMF or UCell scoring when cached score
matrices are available.
