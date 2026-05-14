# PDO CNV Workflow Methodology

This document describes CNV-related workflows under `analysis/cnv/`.

## 1. Aim

CNV scripts infer expression-derived CNA profiles, define malignant CNA
subclones, compare subclones with PDO states/metaprograms, and optionally run
Numbat for allele-aware copy-number validation.

## 2. Core Inputs

InferCNA workflow inputs:

- `PDOs_outs/by_samples/<sample>/<sample>.rds`
- Carroll 2023 non-malignant reference:
  `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds`
- gene order:
  `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`

Subclone workflow inputs:

- `PDOs_outs/cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds`
- `PDOs_outs/cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds`
- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/Auto_PDO_final_states.rds`
- `PDOs_outs/Auto_PDO_mp_adj_noreg.rds`
- `PDOs_outs/UCell_scores_filtered.rds`
- `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`

Numbat workflow inputs:

- velocity/demultiplex CellRanger BAMs
- per-sample QC barcode files
- per-sample count RDS/cell maps
- phased allele counts from `pileup_and_phase.R`
- official Numbat Singularity image

## 3. Script Roles

### `Auto_PDO_infercna.R`

Active InferCNA workflow. It writes target/full InferCNA matrices, sample-level
heatmap caches, scatter plots, and per-sample CNV matrices. It supports quick
replotting by reusing intermediate files when present.

### `Auto_PDO_cnv_subclone_mp_heatmap.R`

Active CNA subclone/state/MP workflow. It intersects InferCNA cells with final
PDO state and MP matrices, infers per-sample CNA subclones, and writes cohort
and sample-level figures and tables.

### `Auto_PDO_cna_diagnostics_SUR1121_SUR1141.R`

Diagnostic-only audit for the SUR1121/SUR1141 expression-derived CNA similarity.
It must not overwrite canonical InferCNA or subclone outputs.

### Numbat scripts

The Numbat export, sample-run, pileup, and concordance scripts are currently
untracked. They should remain unstaged unless the user explicitly asks to stage
the Numbat workflow. Their methodology should stay under
`analysis/methodology/cnv/`.

## 4. Output Standards

New CNV scripts should write:

- `intermediate/` for InferCNA/Numbat matrices and plotting caches
- `tables/` for subclone, concordance, and diagnostic summaries
- `figures/` for heatmaps and scatter plots
- `logs/` for run summaries
- `reports/` for multi-page sample reports

CNV heatmaps and slide figures should use readable labels and legends at the
final PDF dimensions. If a heatmap is generated from a long-running matrix
calculation, the matrix must be cached before plotting.
