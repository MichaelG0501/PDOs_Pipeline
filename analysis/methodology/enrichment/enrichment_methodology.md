# PDO Enrichment Methodology

This document describes scripts under `analysis/enrichment/` and enrichment
stages embedded in cell-state workflows.

## 1. Aim

Enrichment scripts annotate PDO metaprograms or treatment-response signatures
against curated references:

- MSigDB Hallmark
- GO Biological Process
- 3CA metaprograms
- developmental stage references
- focused pathway references such as WNT

## 2. Core Inputs

- `PDOs_outs/MP_outs_default.rds` or the selected
  `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`
- `PDOs_outs/cluster_enrich.rds` where enrichment has already been computed
- 3CA MPs: `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`
- developmental references:
  `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/*.rds`
- MSigDB Hallmark via `msigdbr`
- GO mappings via `clusterProfiler`, `org.Hs.eg.db`, and related packages

Any script that downloads or reads external reference files must document the
exact path or URL in its header.

## 3. Script Roles

### `enrichment_annotation.R`

Terminal annotation workflow for PDO metaprograms. It consumes the optimal MP
result and enrichment tables, then writes heatmaps and annotation summaries.

### `enrichment_extract.R`

Utility script to extract or reshape enrichment results for downstream tables.

### `enrichment_plotting.R` and `enrich_plot.R`

Plotting/helper scripts for enrichment figures. New code should avoid copying
plot helpers into unrelated scripts; reusable helpers should move to
`analysis/shared/` or a dedicated enrichment helper.

### `create_mp_excel.R`

Exports metaprogram annotation tables for manual review or Excel-based
summaries.

### `wnt_enrich.R` and `scGSEA.R`

Focused pathway or scoring analyses. Treat as terminal unless another script is
documented as consuming their outputs.

## 4. Output Standards

New enrichment workflows should write:

- `intermediate/` for enrichment result RDS objects
- `tables/` for long enrichment result CSVs
- `figures/` for heatmaps and dot plots
- `logs/` for run summaries
- `reports/` for multi-page enrichment annotation PDFs

Large enrichment computations should support cached result reuse so heatmap
font sizes, label rotation, or color scales can be changed without recomputing
GO or Hallmark statistics.
