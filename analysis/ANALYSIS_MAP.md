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
| `cnv/Auto_PDO_numbat_export_inputs.R`, `Auto_PDO_numbat_run_sample.R`, `Auto_PDO_numbat_concordance_heatmaps.R`, `Auto_PDO_numbat_concordance_summary_plots.R`, `Auto_PDO_numbat_subclone_mp_heatmap.R` | untracked optional CNV validation | velocity/demux BAMs, allele counts, Numbat outputs; optional conservative clone layer via `PDO_NUMBAT_CLONE_MODE=conservative` | Numbat manifests, clone calls, raw/conservative concordance heatmaps, summary plots, MP/state clone reports | optional; do not stage unless requested |
| `cnv/Auto_PDO_numbat_phylogeny_visualisation.R` | active terminal Numbat tree visualization | Numbat manifest, per-sample `tree_final_<iter>.rds`, `clones_<iter>.rds` | per-sample phylogeny PDF and tree summary CSV | terminal audit figure |
| `cnv/Auto_PDO_numbat_conservative_recut.R` | active optional conservative Numbat clone layer | cached Numbat `treeML`, genotype matrix, expression posterior, allele posterior | conservative clone posterior CSVs, tree-cut sweep, conservative phylogeny PDF | optional downstream clone layer for concordance and MP/state analyses |
| `cnv/wes_subclone/Auto_*` | active WES FACETS/CNVkit input workflow; PyClone outputs legacy | Sarek tumour-normal CRAMs/Mutect2 VCFs, original Sarek launch script, Broad/GATK GRCh38 FASTA, UCSC hg38 `snp151Common` FACETS VCF | FACETS allele-specific segments/purity/ploidy and high-resolution CNVkit-grid `.cns` inputs under `PDOs_outs/Auto_wes_subclone/` | upstream input tier for absolute WES visualization and HATCHet/THetA2 clone-CNA audits; PyClone-VI SNV clusters are not clone-specific CNA |
| `cnv/wes_subclone/legacy_Auto_wes_scrna_subclone_highres_audit.R` | legacy diagnostic/replot | live WES subclone FACETS/PyClone outputs, Sarek CNVkit `.cns`, ephemeral Numbat and inferCNA outputs | previously wrote `figures_highres/`, `tables/visualisation_highres/`, and `PDOs_outs/cnv/cnv_compare_highres/` | superseded by conditional-shift absolute WES plots and HATCHet clone-CNA comparison; retained as method history only |
| `cnv/wes_subclone/Auto_wes_absolute_highres_subclone_compare.R`, `Auto_prepare_hatchet_clone_cna_env.sh`, `Auto_run_hatchet_clone_cna_sample.sh`, `Auto_hatchet_clone_cna_compare.R` | active WES absolute/clone-CNA workflow | FACETS allele-specific segments/purity/ploidy, CNVkit-resolution `.cns`, HATCHet `.bb` inputs, native Numbat outputs from ephemeral | conditional-shift absolute `.cns` and figures under `PDOs_outs/Auto_wes_absolute_cna/`; HATCHet `best/chosen` UCN tables and comparison figures under `PDOs_outs/Auto_wes_clone_cna/` | preferred current WES/scRNA CNA visualization and model-based WES clone-CNA audit; HATCHet selected one tumour clone per WES sample, so output is not evidence for multiple WES CNA subclones |
| `cnv/wes_subclone/Auto_prepare_theta2_clone_cna_env.sh`, `Auto_run_theta2_clone_cna_sample.sh`, `Auto_make_theta2_snp_counts.R`, `Auto_theta2_n3_clone_cna_compare.R` | active diagnostic WES clone-CNA sensitivity | conditional-shift WES `.cns` when present, `reference.cnn`, FACETS `snp-pileup`, THetA2 outputs | forced THetA2 n3 imported `.cns`, summary, and comparison figures under `PDOs_outs/Auto_wes_clone_cna/` | secondary audit only; THetA2 solutions did not improve the main WES/scRNA match over conditional-shift CNVkit bulk |
| `cnv/wes_subclone/Auto_make_phylowgs_inputs.R`, `Auto_prepare_phylowgs_env.sh`, `Auto_run_phylowgs_sample.sh`, `Auto_summarise_phylowgs_results.R`, `Auto_plot_phylowgs_numbat_compare.R` | active integrated WES SNV+CNA phylogeny and terminal visualization | FACETS purity/ploidy/allele-specific segments plus FACETS-aware PyClone-VI SNV count/results tables; native Numbat `bulk_clones` outputs from ephemeral | PhyloWGS SSM/CNV inputs, result bundles, top-tree population/CNV/SSM/inherited-clone-CNA tables under `PDOs_outs/Auto_wes_subclone/tables/phylowgs/` and `reports/phylowgs/`; PhyloWGS-vs-Numbat figures under `PDOs_outs/Auto_wes_subclone/figures/phylowgs/` and visualisation tables under `tables/phylowgs_visualisation/` | current best integrated attempt to assign FACETS CNA events onto SNV-defined phylogenetic populations; visualization shows sparse assigned FACETS CNA events against genome-wide Numbat profiles and is not proof of unique full clone CNA genomes from single WES when events are shared or weakly resolved |
| `cnv/CNV_filter.R`, `cnv_profile.R`, `plot_CNV.R` | historical CNV utilities | varies | older CNV outputs | no new downstream use documented |
| `demultiplex/Auto_*` | organized demultiplex workflow | FASTQs, CellRanger/Souporcell, WES VCFs | demultiplex rerun outputs and verification | external staging; see methodology |

####################
| `demultiplex/Auto_00_submit_demultiplex_pool.sh`, `Auto_01_cellranger_pdo_pool.sh`, `Auto_02_souporcell_pdo_pool.sh`, `Auto_03_reference_and_assign.sh`, `Auto_04_write_demultiplexed_counts.R` | active generic multiplexed-PDO workflow | paired FASTQ folder, matching FASTA, optional WES/VCF reference calls | ephemeral Cell Ranger/Souporcell intermediates; live assignment audits, submission manifests, and donor or temporary-cluster count CSVs | `reference` mode produces donor matrices; no-VCF `temporary` mode emits clearly provisional `TEMP_<pool>_SouporcellCluster<id>_PDO.csv` singlet matrices pending later genotype replacement |
####################
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

####################

####################
## 2026-07-08 Conditional WES Shift Update

`cnv/wes_subclone/Auto_wes_absolute_highres_subclone_compare.R` now writes
`PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_<sample>_*_conditional_shift_absolute.cns`.
The displayed WES absolute row applies the FACETS ploidy offset only when native
Numbat `bulk_clones` median log2 indicates a global amplified baseline. Current
policy: `PDO_1090_vs_NT_1090` applied shift 0; `PDO_1181_vs_NT_1181` applied
shift 0.681. Numbat pseudo-bulk is plotted from gene-level weighted
`bulk_clones` values when available.
####################
## 2026-07-08 CNV Map Update

| Script | Status | Inputs | Outputs | Downstream Use |
| :--- | :--- | :--- | :--- | :--- |
| `cnv/wes_subclone/Auto_wes_absolute_highres_subclone_compare.R` | active terminal diagnostic/replot | live `PDOs_outs/Auto_wes_subclone/tables/cns_highres/Auto_<sample>_*.cns`, live FACETS purity/ploidy tables, ephemeral Numbat by-sample outputs, ephemeral inferCNA matrix/metadata when available | corrected high-resolution conditional-shift `.cns` files under `PDOs_outs/Auto_wes_absolute_cna/tables/cns/`, corrected `Auto_wes_absolute_cna_compare_<sample>.pdf/.png`, and summary/correlation CSVs under `PDOs_outs/Auto_wes_absolute_cna/tables/` | preferred current absolute WES/scRNA CNA visualization; preserves CNVkit resolution, applies `log2(FACETS ploidy / 2)` only when native Numbat indicates a globally amplified baseline, and shows projected WES SNV-cluster CNA rows next to native Numbat/inferCNA tracks |
| `cnv/wes_subclone/Auto_hatchet_clone_cna_compare.R` | active terminal diagnostic/replot | live HATCHet `best.bbc.ucn`, conditional-shift WES bulk `.cns`, native Numbat `bulk_clones` outputs | HATCHet-vs-Numbat figures, clone summary, segment table, and correlations under `PDOs_outs/Auto_wes_clone_cna/` | current model-based clone-specific WES CNA audit; HATCHet uses shifted read depth and FACETS purity but selected one tumour clone per pair, with SUR1181 discordant from the all-amplified Numbat/CNVkit bulk baseline |

The `analysis/cnv/wes_subclone/Auto_run_wes_absolute_cna_compare.sh` wrapper
now runs `Auto_wes_absolute_highres_subclone_compare.R`. The older
`legacy_Auto_wes_absolute_cna_compare.R` FACETS-only output is retained for method
history but should not be used as the presentation comparison because its WES
track has too few segments.

Cleanup note: obsolete PyClone-projection figures/tables, low-resolution `.cns`,
old globally ploidy-adjusted `.cns`, copied HATCHet intermediates, and root PBS
debug logs were removed on 2026-07-08. The removal manifest is
`PDOs_outs/Auto_wes_clone_cna/tables/Auto_wes_cleanup_manifest.tsv`.
####################

####################
## 2026-07-10 WES Subclone CNA Correction

`cnv/wes_subclone/Auto_prepare_facets_snp_genome_order.sh` creates the current
FACETS common-SNP resource
`PDOs_outs/Auto_wes_subclone/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.genome_order.vcf.gz`.
Use this genome-order VCF for `snp-pileup`; the older lexicographic VCF caused
chr1/chr10-limited pileups after CRAM traversal.

`cnv/wes_subclone/Auto_run_phylowgs_sample.sh` supports
`PHYLOWGS_RUN_SUFFIX` for fresh MCMC run directories without deleting prior
outputs. The current accepted PhyloWGS run used
`PHYLOWGS_RUN_SUFFIX=_genome_order_20260709`.

`cnv/wes_subclone/Auto_plot_phylowgs_numbat_compare.R` now plots PhyloWGS
final clone CNA profiles directly: inherited FACETS/PhyloWGS CNA events are
shown on the conditional WES baseline, and non-event intervals are filled as
neutral rather than copied from the WES bulk backbone. Native high-resolution
Numbat `bulk_clones_final.tsv.gz` from the ephemeral sample output is the
preferred scRNA source; conservative live Numbat files are only a fallback and
are treated as already log2-scaled when used.
####################
