# PDO Clinical Association Methodology

This document describes clinical association and survival scripts under
`analysis/clinical/` and the newer final clinical plotting scripts currently
kept under `analysis/plotting/`.

## 1. Aim

The clinical scripts test whether finalized PDO states or PDO metaprogram
scores are associated with clinical variables such as gender, age group,
clinical response, tumour type, histology, T stage, collection timepoint, and
batch.

For presentation work, final clinical figures should use the finalized five
state vector:

- `PDOs_outs/Auto_PDO_final_states.rds`

and should avoid older pre-final noreg states except when explicitly labelled
as a comparison.

## 2. Main Inputs

Core inputs:

- `PDOs_outs/PDOs_merged.rds` or `PDOs_outs/PDOs_final.rds`
- `PDOs_outs/Auto_PDO_final_states.rds`
- `PDOs_outs/UCell_scores_filtered.rds`
- optional `PDOs_outs/Auto_PDO_mp_adj_noreg.rds`

Clinical workbook:

- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx`

Bulk survival inputs, where used:

- TCGA TPM and metadata files under `/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/`

## 3. Current Script Roles

### `clinical_variable_plots.R`

Historical clinical state-composition plotting script adapted from scRef. It
builds stacked state-composition summaries by clinical group. It should be
treated as legacy once the final stacked/boxplot clinical scripts are merged
into one canonical script.

### `clinical_mp_ucell_plots.R`

Historical MP UCell clinical association plotting script. It remains useful as
an MP-focused diagnostic, but final slide figures should use the merged final
clinical association workflow once staged.

### `survival_clinical_mps.R`

Survival/clinical association script using MP and state summaries. It is a
terminal analysis and should not create state vectors for downstream use.

### `analysis/plotting/Auto_PDO_clinical_assoc_stacked_final.R`

Untracked final stacked-bar clinical association workflow. It uses final PDO
states and the cleaned scRef-style stacked layout. It should be moved or merged
into `analysis/clinical/` before being committed.

### `analysis/plotting/Auto_PDO_clinical_assoc_boxplots_scref_style.R`

Untracked final boxplot workflow for MP/state clinical associations. It should
be merged with the stacked final workflow into one canonical clinical figure
script before being committed.

Recommended future canonical script name:

- `analysis/clinical/clinical_association_final_figures.R`

## 4. Plotting Standards

Clinical plots are intended for PowerPoint slides. They should therefore use:

- readable axis and legend fonts
- explicit sample/cell counts where possible
- fixed canonical state order
- fixed canonical state colors
- PDF output as the primary format
- high-resolution PNG only when needed for slide software compatibility

Legends, row names, and point sizes must be sized relative to final figure
dimensions so labels remain readable after insertion into slides.

## 5. Output Tiers And Replotting

Final clinical scripts should write:

- `tables/` for association statistics and plotting data
- `figures/` for PDFs/PNGs
- `logs/` for run summaries
- `intermediate/` for cached joined metadata or score matrices

When only aesthetics change, the script should support replotting from cached
plotting data rather than rebuilding Seurat or clinical joins.
