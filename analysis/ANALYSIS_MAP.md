# PDO Analysis Map

This document is the canonical map for `analysis/`. Update it whenever a script
is added, renamed, superseded, moved, or given a new downstream dependency.

## Current Defaults

- Preferred state definition: `Approach B, noreg`
- Preferred final state object: `PDOs_outs/Auto_PDO_final_states.rds`
- Pre-final noreg state object: `PDOs_outs/Auto_PDO_states_noreg.rds`
- Pre-final noreg MP activity object: `PDOs_outs/Auto_PDO_mp_adj_noreg.rds`
- Optimal PDO GeneNMF result: `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`
- Shared constants: `analysis/shared/Auto_pdo_analysis_config.R`
- Shared helper functions: `analysis/shared/Auto_pdo_analysis_helpers.R`
- Output tiers for new long-running analyses: `intermediate/`, `tables/`, `figures/`, `logs/`, `reports/`

## Run Order

1. Core QC and GeneNMF pipeline
   - `QC_Pipeline.R`
   - `geneNMF.R`
   - optional classical per-sample `NMF.R`

2. Metaprogram selection and scoring
   - `analysis/metaprograms/extend_nMP_range.R` if more nMP values are needed
   - `analysis/metaprograms/find_optimal_nmf.R`
   - `analysis/metaprograms/update_optimal_mp.R`
   - `analysis/metaprograms/mp_ucell_scoring.R`

3. Current PDO state definition
   - `analysis/cell_states/PDO_states_analysis.R`
   - Use the noreg outputs for final-state relabeling and selected trajectory workflows.

4. Final-state relabeling and final-state vector
   - `analysis/cell_states/PDO_unresolved_relabel.R`
   - `analysis/cell_states/PDO_finalize_states.R`
   - Preferred downstream state vector: `PDOs_outs/Auto_PDO_final_states.rds`

5. Final-state terminal figures and tables
   - `analysis/cell_states/pdo_overall_state_proportions.R`
   - `analysis/cell_states/sample_abundance_pdo.R`
   - `analysis/cell_states/Auto_pdo_sn_matched_pair_comparison.R`
   - `analysis/cell_states/Auto_compare_untreated_proportions.R`
   - `analysis/cell_states/Auto_five_state_markers.R`
   - `analysis/cell_states/Auto_five_state_surface_markers.R`
   - `analysis/cell_states/Auto_marker_comparison_excel.R`
   - `analysis/cell_states/Auto_PDO_scAtlas_scenic_comparison.R`
   - `analysis/cell_states/Auto_PDO_final_mp_scenic.R`

6. Matched-FLOT response analyses
   - Canonical current script: `analysis/cell_states/Auto_pdo_flot_matched_response.R`
   - Optional high-resolution GeneNMF route:
     `Auto_pdo_flot_matched_geneNMF.R` ->
     `Auto_pdo_flot_matched_highres_mp_trend_filter.R` ->
     `Auto_pdo_flot_highres_enrichment_annotation.R`

7. Clinical association and survival
   - Historical scripts: `analysis/clinical/clinical_variable_plots.R`,
     `analysis/clinical/clinical_mp_ucell_plots.R`, and
     `analysis/clinical/survival_clinical_mps.R`
   - Final clinical plotting uses the merged canonical script:
     `analysis/clinical/clinical_association_final_figures.R`

8. CNV, demultiplex, trajectory, drug reversal, and enrichment
   - See the dependency table and folder methodology files below.

## Dependency Map

| Script | Status | Key inputs | Key outputs | Downstream use |
| :--- | :--- | :--- | :--- | :--- |
| `shared/Auto_pdo_analysis_config.R` | active shared | none | constants | all new scripts |
| `shared/Auto_pdo_analysis_helpers.R` | active shared | shared config | helper functions | all new scripts |
| `metaprograms/find_optimal_nmf.R` | active upstream | `Metaprogrammes_Results/geneNMF_metaprograms_nMP_{k}.rds` | nMP diagnostic plots/tables | selects nMP=13 |
| `metaprograms/extend_nMP_range.R` | active upstream | `geneNMF_outs.rds` | additional `geneNMF_metaprograms_nMP_{k}.rds` | nMP selection |
| `metaprograms/update_optimal_mp.R` | active upstream | selected nMP object, `PDOs_merged.rds` | `MP_outs_default.rds`, updated MP tables/plots | enrichment and state workflows |
| `metaprograms/mp_ucell_scoring.R` | active upstream | `PDOs_merged.rds`, optimal MP object | `UCell_scores_filtered.rds` | state definition, MP figures |
| `metaprograms/PDO_mp_correlation_crossdata.R` | terminal comparison | PDO/scRef MP objects and UCell scores | cross-data Jaccard/scatter/bar plots | terminal |
| `metaprograms/mp_correlation_pdo.R` | terminal comparison | `PDOs_final.rds`, MP object | PDO MP correlation plots/RDS | terminal |
| `metaprograms/Auto_3CA_pseudobulk_correlation_crossdata.R` | terminal comparison | `PDOs_merged.rds`, scRef object, 3CA MPs, OSCC GEO files | pseudobulk/bulk 3CA correlation plots/tables | terminal |
| `metaprograms/Auto_mp_chromosomal_mapping_pdo.R` | untracked terminal | MP gene lists and gene coordinates | chromosomal mapping figures | terminal; do not stage unless requested |
| `metaprograms/Find_NMF.R`, `robust_NMF.R`, `robust_nmf_scref.R`, `MP_analysis_pdos.R`, `MP_dist.R`, `nmf_plot.R`, `compute_pdo_mp_in_scref.R` | historical utilities | varies | NMF diagnostics/utilities | no canonical downstream dependency documented |
| `cell_states/PDO_states_analysis.R` | active upstream, historical name | `PDOs_merged.rds`, `UCell_scores_filtered.rds`, nMP=13 object, cell-cycle genes | `Auto_PDO_states_noreg.rds`, `Auto_PDO_mp_adj_noreg.rds`, `Auto_PDO_top_mp.rds` | unresolved relabel, final state, selected trajectory |
| `cell_states/PDO_unresolved_relabel.R` | active upstream, historical name | noreg states, PDO/3CA UCell, nMP=13 object, 3CA MPs | `unresolved_states/Auto_PDO_unresolved_relabel_states.rds`, coverage table, figures | final state |
| `cell_states/PDO_finalize_states.R` | active upstream, historical name | unresolved relabel states, PDO/3CA UCell, nMP=13 object | `Auto_PDO_final_states.rds`, final-state figures | preferred state vector for downstream |
| `cell_states/pdo_overall_state_proportions.R` | terminal figure | `Auto_PDO_final_states.rds` | overall final-state proportion plot | terminal |
| `cell_states/sample_abundance_pdo.R` | terminal figure/table | `PDOs_merged.rds`, final states, UCell, clinical workbook | sample/state abundance figures and tables | terminal |
| `cell_states/Auto_pdo_sn_matched_pair_comparison.R` | terminal comparison | PDO final states/MPs, snRNA-seq malignant object/states/MPs | two-page matched PDO/snRNA-seq comparison PDF and tables | terminal |
| `cell_states/Auto_compare_untreated_proportions.R` | terminal figure/table | `PDOs_all_meta.rds`, final states | SUR1090/SUR1072 untreated comparison | terminal |
| `cell_states/Auto_five_state_markers.R` | terminal marker workflow | `PDOs_merged.rds`, final states | five-state marker tables, heatmaps, caches | marker comparison/surface-marker workflows |
| `cell_states/Auto_five_state_surface_markers.R` | terminal marker workflow | marker tables, final states, UniProt/surfaceome/GO refs | FACS-oriented surface-marker workbook/tables | terminal |
| `cell_states/Auto_marker_comparison_excel.R` | terminal comparison | scRef and PDO marker caches/tables | cross-dataset marker workbook and heatmaps | terminal |
| `cell_states/Auto_marker_sample_expression_report.R` | terminal report | final states and marker sets | marker expression report | terminal |
| `cell_states/Auto_marker_selection_simulation.R` | terminal simulation | five-state marker cache or final states | marker-gate and qPCR simulation outputs | terminal |
| `cell_states/Auto_scATLAS_four_marker_specificity.R` | terminal figure/table | scATLAS marker specificity cache, PDO final states | four-marker scATLAS specificity heatmap/table | terminal |
| `cell_states/Auto_PDO_scAtlas_scenic_comparison.R` | terminal comparison | PDO/scATLAS SCENIC or specificity outputs | RSS comparison heatmaps/workbook | terminal |
| `cell_states/Auto_PDO_final_mp_scenic.R` | heavy terminal workflow | final states, PDO/3CA UCell, nMP=13 object, cistarget DBs | SCENIC selected cells, regulons, networks | terminal |
| `cell_states/Auto_pdo_flot_matched_response.R` | canonical matched-FLOT workflow | final states, `PDOs_merged.rds`, UCell, noreg MP matrix, Hallmark/cell-cycle refs | cached response object, tables, final presentation PDFs | terminal; replaces older matched scripts |
| `cell_states/Auto_pdo_flot_matched_geneNMF.R` | optional upstream | `PDOs_list_PDOs.rds` | matched-sample GeneNMF object | high-resolution FLOT MP trend workflow |
| `cell_states/Auto_pdo_flot_matched_highres_mp_trend_filter.R` | optional upstream/terminal | matched GeneNMF object, per-sample RDS, optional 3CA/cell-cycle refs | retained high-res MPs, UCell scores, trend plots | high-resolution enrichment |
| `cell_states/Auto_pdo_flot_highres_enrichment_annotation.R` | optional terminal | retained high-res MP genes/trends, GO/Hallmark/3CA/developmental refs | high-res MP enrichment tables/plots | terminal |
| `cell_states/legacy_compare_mp_scoring_state_definition.R` | legacy comparison | `PDOs_merged.rds`, optimal MP object, optional UCell | alternative activity/state call files and comparison plots | no downstream use; should be legacy-prefixed |
| `cell_states/legacy_states_scref_pairwise_nodeplot.R` | legacy comparison | `PDOs_final.rds`, scRef MPs | scRef-derived PDO node plots | no downstream use; should be legacy-prefixed |
| `cell_states/legacy_state_hybrid_subtyping_noreg.R` | legacy comparison | pre-final noreg states/MP matrix | hybrid subtype plots/tables | no downstream use; should be legacy-prefixed |
| `cell_states/legacy_state_hybrid_pairwise_nodeplot_noreg.R` | legacy comparison | pre-final noreg states/MP matrix | noreg hybrid pairwise nodeplot | no downstream use; should be legacy-prefixed |
| `cell_states/legacy_pdo_flot_matched_dge_findmarkers.R` | legacy matched-FLOT DGE | final states, matched samples | old Seurat DGE/enrichment outputs | superseded by canonical matched-FLOT response |
| `cell_states/legacy_pdo_flot_matched_survival_and_state_plots.R` | legacy matched-FLOT/survival | final states, UCell, TCGA inputs | old survival and matched-sample plots | superseded by canonical matched-FLOT response |
| `clinical/clinical_variable_plots.R` | legacy/terminal | final states, clinical workbook | older stacked clinical plots | superseded once final clinical merge is staged |
| `clinical/clinical_mp_ucell_plots.R` | legacy/terminal | MP UCell, clinical workbook | older MP clinical plots | superseded once final clinical merge is staged |
| `clinical/survival_clinical_mps.R` | terminal survival | final states, UCell, clinical/survival fields | survival/association CSV/PDF | terminal |
| `clinical/clinical_association_final_figures.R` | active final clinical script | final states, MP scores, clinical workbook | final stacked and boxplot PDFs/tables | terminal; fully merged and canonical |
| `plotting/heatmap.R` | historical plotting utility | varies | heatmaps | no canonical downstream dependency documented |
| `enrichment/enrichment_annotation.R` | terminal enrichment | optimal MP object, enrichment refs/results | PDO enrichment annotation PDFs | terminal |
| `enrichment/enrichment_extract.R`, `enrichment_plotting.R`, `enrich_plot.R`, `create_mp_excel.R`, `wnt_enrich.R`, `scGSEA.R` | terminal/utilities | enrichment RDS/reference sets | extracted tables and figures | terminal |
| `cnv/Auto_PDO_infercna.R` | active upstream/terminal | per-sample PDO RDS, Carroll ref, gene order | InferCNA matrices, heatmaps, scatter plots, caches | CNA subclone workflow |
| `cnv/Auto_PDO_cnv_subclone_mp_heatmap.R` | active terminal | InferCNA target matrices, final states, MP scores | CNA subclone/state/MP PDFs and tables | terminal; Numbat concordance |
| `cnv/Auto_PDO_cna_diagnostics_SUR1121_SUR1141.R` | untracked diagnostic | InferCNA caches, untreated RDS, gene order | SUR1121/SUR1141 diagnostic tables | diagnostic only; do not stage unless requested |
| `cnv/Auto_PDO_numbat_export_inputs.R`, `Auto_PDO_numbat_run_sample.R`, `Auto_PDO_numbat_concordance_heatmaps.R` | untracked optional CNV validation | velocity/demux BAMs, allele counts, Numbat outputs | Numbat manifests, clone calls, concordance heatmaps | optional; do not stage unless requested |
| `cnv/CNV_filter.R`, `cnv_profile.R`, `plot_CNV.R` | historical CNV utilities | varies | older CNV outputs | no new downstream use documented |
| `demultiplex/Auto_*` | organized demultiplex workflow | FASTQs, CellRanger/Souporcell, WES VCFs | demultiplex rerun outputs and verification | external staging; see methodology |
| `trajectory/Auto_*` | trajectory/velocity workflow | noreg/final states, BAMs, velocyto/scVelo refs | pseudotime/velocity tables and figures | mostly untracked; see methodology |
| `cell_states/Auto_drug_reversal/*` | organized drug-reversal workflow | final states, DEG signatures, ASGARD/scDrugPrio/CLUE refs | drug reversal tables and figures | terminal; see drug methodology |

## Superseded Or No-Downstream Scripts

These scripts are retained for comparison or file-safety history. They should
not be used by new downstream analysis. If renamed, use the `legacy_` prefix:

- `analysis/cell_states/legacy_compare_mp_scoring_state_definition.R`
- `analysis/cell_states/legacy_states_scref_pairwise_nodeplot.R`
- `analysis/cell_states/legacy_state_hybrid_subtyping_noreg.R`
- `analysis/cell_states/legacy_state_hybrid_pairwise_nodeplot_noreg.R`
- `analysis/cell_states/legacy_pdo_flot_matched_dge_findmarkers.R`
- `analysis/cell_states/legacy_pdo_flot_matched_survival_and_state_plots.R`
- `analysis/clinical/clinical_variable_plots.R` after final clinical merge
- `analysis/clinical/clinical_mp_ucell_plots.R` after final clinical merge

Delete-candidates are not deleted by agents. If a script is confirmed to have
no unique purpose after the map review, rename it with a `delete_` prefix and
leave final removal to the user.

## Outdated Downstream Pointers To Avoid

- Do not use alternative state vectors written by
  `legacy_compare_mp_scoring_state_definition.R` for downstream analysis.
- Do not use scRef-derived PDO states from
  `legacy_states_scref_pairwise_nodeplot.R` for downstream analysis.
- Prefer `PDOs_merged.rds` plus `Auto_PDO_final_states.rds` over
  `PDOs_final.rds` in new scripts, unless the existing script explicitly needs
  the old object.
- Use `Auto_PDO_final_states.rds` for terminal state abundance, marker,
  clinical, SCENIC, drug-reversal, and matched-FLOT workflows.
- Use `Auto_PDO_states_noreg.rds` only for method comparison and workflows that
  explicitly require pre-final four-state calls.

## Cache And Replot Policy

Long-running scripts should write heavy intermediates before plotting. Plotting
changes should be reproducible by reading cached `intermediate/` or `tables/`
outputs and regenerating only `figures/` or `reports/`.

Recommended environment toggles:

- `PDO_FORCE_REBUILD=1`: ignore cached intermediates and recompute.
- `PDO_REPLOT_ONLY=1`: read cached intermediates and regenerate plots/reports.

Run summaries should be written to `logs/` and record start/end time, inputs,
outputs, parameters, cached-object reuse, and session/package versions when
relevant.

## External Data And Download Requirements

- Surface marker workflow downloads/caches UniProt reviewed human topology data
  and ETH Zurich human surfaceome Table S3 if local copies are absent.
- Drug-reversal workflows require ASGARD/LINCS, scDrugPrio PPI/drug-target
  resources, and optional CLUE API access as documented in the drug methodology.
- FLOT high-resolution enrichment uses Hallmark via `msigdbr`, GO via
  `org.Hs.eg.db`/`clusterProfiler`, 3CA MPs, and developmental references.
- Numbat requires the official `pkharchenkolab/numbat-rbase:latest` container
  and allele counts generated by `pileup_and_phase.R`.

## Untracked Files Not To Stage

The following paths were untracked at the start of this cleanup and should not
be staged unless the user explicitly asks:

- `Auto_pdo_cnv_subclone_mp.sh`
- `Auto_pdo_flot_highres_enrichment_annotation.sh`
- `Auto_pdo_flot_highres_metaprogram_trends.sh`
- `Auto_pdo_infercna.sh`
- `Auto_qsub_cnv_heatmap.sh`

- `Auto_run_compare_mp_scoring_state_definition.sh`
- `Auto_run_compare_mp_scoring_state_definition_4core.sh`
- `analysis/cnv/Auto_00_submit_pdo_numbat.sh`
- `analysis/cnv/Auto_PDO_cna_diagnostics_SUR1121_SUR1141.R`
- `analysis/cnv/Auto_PDO_numbat_concordance_heatmaps.R`
- `analysis/cnv/Auto_PDO_numbat_export_inputs.R`
- `analysis/cnv/Auto_PDO_numbat_run_sample.R`
- `analysis/cnv/Auto_prepare_pdo_numbat_container.sh`
- `analysis/cnv/Auto_run_pdo_numbat_concordance.sh`
- `analysis/cnv/Auto_run_pdo_numbat_pileup.sh`
- `analysis/cnv/Auto_run_pdo_numbat_sample.sh`
- `analysis/metaprograms/Auto_mp_chromosomal_mapping_pdo.R`
- `analysis/methodology/Auto_PDO_cnv_subclone_methodology.md`

- `analysis/trajectory/`
