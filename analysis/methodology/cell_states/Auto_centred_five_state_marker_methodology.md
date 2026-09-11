# Centred Five-State Marker Methodology

## Scope

This methodology covers `Auto_five_state_markers.R`. It derives recurrent
markers for the five canonical states assigned by
`Auto_06_centred_refined_state_definition_noreg.R`. It does not use the legacy
`Auto_PDO_final_states.rds` vector or uncentred nMP=13 scores.

## Inputs and eligible cells

The script reads the live `PDOs_merged.rds` count matrix and the named live
`centred_mp_refinement/centred_refined_noreg_states.rds` vector. Cell barcodes
are intersected exactly. Only cells in the five canonical biological states are
retained; `Unresolved` and `Hybrid` cells are excluded from marker discovery.
`SUR843T3_PDO` is absent from the canonical upstream analysis.

The state order is sourced from `Auto_pdo_analysis_config.R`: Classic
proliferation, Basal to intestinal metaplasia, SMG to intestinal metaplasia,
Stress adaptive, then PDO medium induced.

## Expression preprocessing and embedding

Genes detected in fewer than 10 retained cells are discarded. A lean Seurat
object is created from raw RNA counts, normalized, and processed with 3,000
variable features, 30 principal components, and clustering resolution 0.5.
The embedding is descriptive QC; cluster identity does not replace the supplied
centred state vector.

## Marker screens and threshold choices

The global screen tests each state against all other eligible cells. Candidate
genes then enter a sample-aware recurrence analysis. A state/sample comparison
is eligible only when both the target state and pooled other states contain at
least 10 cells in that sample. This minimum avoids tests driven by tiny groups
while retaining rarer PDO states. Per-sample tests preserve effect size,
adjusted P value, target/off-target prevalence, and cell counts.

Candidate pools are limited to the top 1,000 genes per state after the global
screen so recurrence calculations remain tractable without restricting the
final result to a small presentation list. The recurrence summary records how
many eligible samples support each gene/state association and retains the
prevalence difference (`pct_state - pct_other`). Final ranking prioritizes
recurrent, state-specific positive markers. The full CSV preserves all ranking
columns; top-five tables are presentation subsets only.

## Persistent outputs and cache behavior

State-by-gene normalized-expression and detection-prevalence matrices support
the marker heatmap, surface-marker workflow, and downstream review. They are
stored under live `PDOs_outs/Auto_five_state_markers/`, with reusable RDS in
`intermediate/`, source/result tables in `tables/`, figures in `figures/`, and
the generated compact method note in `reports/`. The intermediate RDS remain
in live storage because they support replotting and downstream review.

Set `PDO_FORCE_REBUILD=1` to ignore cached embedding, global-marker,
per-sample-DGE, and specificity objects. A forced rebuild is mandatory after a
state-vector or threshold change. Run the workflow through PBS, not on a login
node.

## Interpretation

Cells are measurement units in within-sample marker tests, whereas recurrence
across independently cultured PDO samples is the robustness criterion. The
workflow prioritizes descriptive markers and does not treat cell-level P values
as independent patient-level evidence. Surface suitability is assessed
separately using UniProt topology and surfaceome references.
