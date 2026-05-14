# Shared PDO Analysis Configuration And Logging Methodology

This document describes the shared analysis infrastructure in:

- `analysis/shared/Auto_pdo_analysis_config.R`
- `analysis/shared/Auto_pdo_analysis_helpers.R`

These files do not run a biological analysis by themselves. They define the
defaults that new downstream scripts should reuse instead of copying constants
from older scripts.

## 1. Current Defaults

The preferred state definition is:

- `Approach B, noreg`

The preferred finalized state vector is:

- `PDOs_outs/Auto_PDO_final_states.rds`

The pre-final four-state noreg vector remains available for methods that
explicitly need it, such as pseudotime workflows that should avoid the later
3CA unresolved relabeling:

- `PDOs_outs/Auto_PDO_states_noreg.rds`

The preferred noreg MP activity matrix is:

- `PDOs_outs/Auto_PDO_mp_adj_noreg.rds`

`SUR843T3_PDO` is excluded from NMF and downstream analyses.

## 2. Constants Centralized In `Auto_pdo_analysis_config.R`

The config file centralizes:

- project, analysis, output, and temp directories
- canonical PDO state order and colors
- optional scATLAS-style `Immune Infiltrating` placement
- MP descriptions
- MP-to-state mappings
- cell-cycle MP IDs
- common thresholds:
  - MP silhouette cutoff
  - sample-coverage cutoff
  - state assignment threshold
  - hybrid gap
  - minimum cells for state-sample pseudobulk analyses
- metadata column names
- plot dimensions and font defaults for presentation slides
- external reference paths
- cache/replot environment-variable names

New scripts should source this file before defining workflow-specific
parameters.

## 3. Shared Helper Functions

`Auto_pdo_analysis_helpers.R` provides reusable functions for:

- parsing boolean, numeric, and integer environment variables
- resolving cache policy with `PDO_FORCE_REBUILD` and `PDO_REPLOT_ONLY`
- creating output tier folders
- checking required input files
- preserving names when coercing final-state vectors to character
- retrieving Seurat assay matrices across Seurat v4/v5 APIs
- applying a slide-readable ggplot theme
- writing lightweight run-summary logs

The state-vector helper exists because downstream lookups such as
`state_vec[Cells(seurat_obj)]` silently fail if barcode names are dropped.

## 4. Output Tiers

Long-running workflows should create these subfolders under their own output
directory:

- `intermediate/` for heavy RDS/cache/model objects
- `tables/` for final CSV/TSV/XLSX tables
- `figures/` for PDF/PNG plot exports
- `logs/` for run summaries and session information
- `reports/` for multi-page narrative PDFs or markdown reports

Existing older scripts may still write flat outputs. When they are updated,
move new outputs into tiers without deleting the previous output names unless
the user explicitly requests a breaking change.

## 5. Cache And Replot Policy

Heavy scripts should implement:

- `PDO_FORCE_REBUILD=1` to ignore cached intermediates and recompute.
- `PDO_REPLOT_ONLY=1` to reuse cached intermediates and regenerate only plots
  and reports.

The run summary must record whether cached objects were reused.

## 6. Run Summary Logging

Long-running scripts should write a lightweight text log to `logs/` containing:

- script name
- start/end or write time
- input files
- output files
- parameters used
- cache/replot settings
- session/package versions when relevant

The helper `pdo_write_run_summary()` provides the default implementation.
