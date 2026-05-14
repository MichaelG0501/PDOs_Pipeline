# PDO Cell-State Workflow Methodology

This document describes the current cell-state workflow under
`analysis/cell_states/`. It is written against the current code and the
current downstream convention that `Approach B, noreg` is the chosen state
definition.

## 1. Current State-Definition Choice

The active state-definition route is:

1. Build PDO metaprogram activity scores from the optimal GeneNMF nMP=13 result.
2. Apply the PDO two-step MP filter:
   - remove MPs with silhouette `< 0`
   - remove MPs with sample coverage `< 0.25`
3. Assign the four PDO-intrinsic states using Approach B with no cell-cycle
   regression:
   - `Classic Proliferative`
   - `Basal to Intest. Meta`
   - `SMG-like Metaplasia`
   - `Stress-adaptive`
4. Mark low-confidence cells as `Unresolved`.
5. Mark near-tied cells as `Hybrid`.
6. Relabel unresolved cells with retained 3CA support where appropriate.
7. Collapse finalized labels into the five-state vector:
   - `Classic Proliferative`
   - `Basal to Intest. Meta`
   - `SMG-like Metaplasia`
   - `Stress-adaptive`
   - `3CA_EMT_and_Protein_maturation`

The preferred final object for downstream analysis is:

- `PDOs_outs/Auto_PDO_final_states.rds`

The pre-final noreg state vector remains available for direct method comparison
and selected trajectory workflows:

- `PDOs_outs/Auto_PDO_states_noreg.rds`

## 2. Active Upstream State Scripts

### `PDO_states_analysis.R`

Status: active upstream, but name is historical.

Purpose:

- runs the PDO-adapted Approach B noreg state definition
- writes the noreg state vector, noreg MP matrix, and top-MP assignments
- generates noreg QC heatmaps and proportion summaries

Key inputs:

- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/UCell_scores_filtered.rds`
- `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`

Key outputs:

- `PDOs_outs/Auto_PDO_states_noreg.rds`
- `PDOs_outs/Auto_PDO_mp_adj_noreg.rds`
- `PDOs_outs/Auto_PDO_top_mp.rds`
- noreg heatmap and state-proportion PDFs

Recommended future rename:

- `state_definition_approach_b_noreg.R`

### `PDO_unresolved_relabel.R`

Status: active upstream, but name is historical.

Purpose:

- examines `Unresolved` cells from the noreg Approach B state vector
- uses PDO MPs plus 3CA MP support to rescue unresolved cells into the finalized
  3CA EMT/protein-maturation compartment where supported
- writes an intermediate relabeled state vector and coverage table

Key inputs:

- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/UCell_scores_filtered.rds`
- `PDOs_outs/UCell_3CA_MPs.rds`
- `PDOs_outs/Auto_PDO_states_noreg.rds`
- `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`
- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`

Key outputs:

- `PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_states.rds`
- `PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_mp_coverage.csv`
- relabel heatmaps, proportion plots, and optional survival summaries

Recommended future rename:

- `final_state_unresolved_relabel_3ca.R`

### `PDO_finalize_states.R`

Status: active upstream, but name is historical.

Purpose:

- converts the unresolved-relabel state vector into the final five-state PDO
  vector used by downstream analysis
- preserves original barcode names during normalization
- redraws final-state heatmaps and proportion plots

Key inputs:

- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/UCell_scores_filtered.rds`
- `PDOs_outs/UCell_3CA_MPs.rds`
- `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`
- `PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_states.rds`

Key outputs:

- `PDOs_outs/Auto_PDO_final_states.rds`
- `PDOs_outs/Auto_PDO_final_states_heatmap.pdf`
- `PDOs_outs/Auto_PDO_final_states_proportion.pdf`

Recommended future rename:

- `final_state_vector_finalize.R`

## 3. Terminal State-Figure Scripts

Scripts such as `pdo_overall_state_proportions.R`,
`sample_abundance_pdo.R`, `Auto_pdo_sn_matched_pair_comparison.R`,
`Auto_PDO_final_mp_scenic.R`, marker workflows, and SCENIC comparison scripts
consume `Auto_PDO_final_states.rds` and should not write state vectors that are
used downstream unless explicitly documented in `analysis/ANALYSIS_MAP.md`.

## 4. Legacy And Comparison Scripts

Scripts that compare alternative state definitions or older scoring methods
should be clearly marked as `legacy` in their header and in the analysis map.
They may produce comparison plots, but their state vectors should not be used
by downstream analysis.

Known legacy/comparison categories:

- scRef-derived state calls on PDOs
- pairwise hybrid plots based only on pre-final noreg states
- MP scoring comparisons that write alternative state calls
- old matched-FLOT DGE/survival analyses superseded by
  `Auto_pdo_flot_matched_response.R`

## 5. Output And Cache Standards

New or substantially revised cell-state scripts should write output tiers under
their workflow output directory:

- `intermediate/` for heavy state, score, or model objects
- `tables/` for cell-level or summary CSVs
- `figures/` for plot files
- `logs/` for run summaries
- `reports/` for multi-page presentation PDFs

Heavy scripts should support `PDO_FORCE_REBUILD=1` and `PDO_REPLOT_ONLY=1`
where a cached intermediate object makes replotting possible.
