# PDO Analysis Map

This is the authoritative status, dependency, and run-order map for
`analysis/`. Script headers contain exact file-level I/O. Methodology files are
required only for workflows with non-trivial scientific choices, thresholds, or
external references.

## Status vocabulary

- **active upstream**: creates a current downstream input.
- **active terminal**: creates current figures/tables but no canonical input.
- **active support**: wrapper, helper, configuration, or setup component.
- **legacy**: retained for provenance; must not feed current analysis.
- **delete-candidate**: retained for manual review; agents do not delete it.

## Canonical current objects

| Purpose | Live path |
| :--- | :--- |
| Merged PDO Seurat object | `PDOs_outs/PDOs_merged.rds` |
| Final centred-refined MP genes | `PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds` |
| Final centred-refined MP weights | `PDOs_outs/centred_mp_refinement/merged_refined_mp_gene_weights.rds` |
| Final centred-refined UCell scores | `PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds` |
| Final centred state vector | `PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds` |
| Adjusted centred MP matrix | `PDOs_outs/centred_mp_refinement/centred_refined_noreg_mp_adj.rds` |
| State-group maximum matrix | `PDOs_outs/centred_mp_refinement/centred_refined_noreg_group_max.rds` |
| MP display grouping/order | `PDOs_outs/centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv`; `centred_refined_mp_strict_order.rds` |
| Centred enrichment | `PDOs_outs/centred_mp_refinement/cluster_enrich_centred.rds` |
| External 3CA UCell scores | `PDOs_outs/UCell_3CA_MPs.rds` |

The uncentred nMP=13 objects, `UCell_scores_filtered.rds`,
`Auto_PDO_states_noreg.rds`, and `Auto_PDO_final_states.rds` are legacy and
must not be used by new downstream work.

## Canonical run order

| Order | Script | Status | Inputs | Outputs / downstream use |
| :--- | :--- | :--- | :--- | :--- |
| 1 | `metaprograms/centred/Auto_01_centred_geneNMF.R` | active upstream | `PDOs_list_PDOs.rds`; excludes SUR843T3_PDO | raw centred multiNMF cache; nMP 4:25 objects |
| 2 | `metaprograms/centred/Auto_02_nmf_rank_selection_diagnostics.R` | active upstream/QC | nMP 4:25 objects | `optimal_nMP.rds`; rank/enrichment/batch QC |
| 3 | `metaprograms/centred/Auto_03_mp_refinement_submp.R` | active upstream | optimal object, raw NMF programs, merged PDO | refined genes/weights/assignments/UCell and split audit |
| 4 | `metaprograms/centred/Auto_04_mp_refinement_merge_correlated_submps.R` | active canonical upstream | step-03 objects | final filtered genes, weights, assignments, UCell, enrichment |
| 4a | `metaprograms/centred/Auto_07a_ucell_scoring.R` | active score utility | final genes, 3CA genes, merged PDO | refreshes current PDO and 3CA UCell matrices before state calling |
| 5 | `metaprograms/centred/Auto_05_centred_refined_mp_ordered_heatmaps.R` | active upstream/terminal | final step-04 objects | persistent biological grouping/order; QC heatmaps |
| 6 | `metaprograms/centred/Auto_06_centred_refined_state_definition_noreg.R` | active canonical upstream | merged PDO, final UCell | current named state vector and normalized score matrices |
| 7b | `metaprograms/centred/Auto_07_3ca_vs_refined_mp_correlation.R` | active terminal | current UCell/state matrices | terminal all-cell/unresolved-cell comparison |
| 8 | `metaprograms/centred/Auto_08_tcga_mp_survival_volcano_centred.R` | active terminal | final genes, TCGA ESCA TPM/metadata | survival CSV and PDF |
| table | `metaprograms/centred/Auto_06_centred_refined_mp_annotation_excel_export.R` | active terminal | final genes and optimal object | current annotated MP workbook |

Step 04 is the sole final MP filter: at least three observed PDO samples and at
least five genes after refinement/merging. Step 05 must not impose another
gene-count or coverage filter. State definition uses fixed normalized-score
thresholds 0.5 (minimum best group) and 0.3 (best-minus-second hybrid gap).
See the two centred methodology files for full rationale.

## Current centred-state downstream analyses

| Script | Status | Current inputs | Outputs / downstream use |
| :--- | :--- | :--- | :--- |
| `cell_states/Auto_five_state_markers.R` | active upstream/terminal | merged PDO + centred state vector | live recurrent marker tables/caches; feeds surface markers |
| `cell_states/Auto_selected_state_marker_expression.R` | active terminal | merged PDO + centred state vector; four matched untreated/FLOT-treated pairs | prespecified four-state marker heatmaps by state/sample/treatment, paired patient-state deltas, and paired per-cell boxplots/source tables |
| `cell_states/Auto_five_state_surface_markers.R` | active terminal | current marker tables + centred states + UniProt/surfaceome | FACS-ranked CSV/workbook/figure |
| `cell_states/sample_abundance_pdo.R` | active terminal | centred state/UCell + clinical workbook | `Auto_sample_abundance_pdo/` PDF |
| `cell_states/centred_refined_pdo_flot_matched_response.R` | active terminal | centred states, UCell, adjusted/group-max matrices, strict order | live matched-FLOT cache, source tables, figures, logs |
| `cell_states/Auto_PDO_scRef_centred_state_concordance.R` | active terminal comparison | independently defined PDO/scRef centred states and MPs | concordance cache, tables, report |
| `cell_states/Auto_PDO_treatment_scenic.R` | active terminal | eight matched PDO samples + cisTarget databases | live `treatment_scenic/` regulon matrices/tables/figures; ephemeral SCENIC work cache |
| `cell_states/Auto_PDO_personalised_state_regulon_rss.R` | active terminal | merged PDO + centred states + four fixed scATLAS regulons | per-specimen projected AUCell activity, binary/multiclass RSS gaps, source tables, and explanatory figures |
| `cell_states/Auto_crispr_selective_dependency_analysis.R` | active terminal | final centred PDO MP genes/weights, merged PDO + centred states, oesophageal differential-dependency workbook, current scATLAS epithelial atlas + centred MP genes/weights/states + ranked state markers | independent PDO and scATLAS MP/state overlap and expression-matched depletion-burden scores; fresh PDO six-state sample-blocked edgeR signatures; scATLAS ranked top-150 state signatures; four Classic-proliferation scATLAS/CRISPR/PDO three-way overlaps; gene-level source tables, workbook, and figures |
| `clinical/Auto_centred_clinical_auc_response_associations.R` | active terminal | merged PDO + current centred states/UCell/MP genes/order + `PDO_ClinicalMetadata_V3.xlsx` | categorical clinical state/MP associations, 12-PDO AUCrel linear patterns, and sensitive-versus-resistant state/MP/pathway outputs |

Current marker and surface-marker runs must follow step 06 and be forced after a
state-vector change. The matched-FLOT and SCENIC workflows are terminal and do
not redefine states.

## Centred high-resolution matched-FLOT MP route

This terminal route deliberately preserves smaller, diverse treatment-trend
MPs and does not replace the canonical centred-refined MPs or states.

| Order | Script | Status | Inputs | Outputs / downstream use |
| :--- | :--- | :--- | :--- | :--- |
| F1 | `cell_states/Auto_pdo_flot_matched_geneNMF.R` | active upstream | eight matched samples from `PDOs_list_PDOs.rds` | centred (`center=TRUE`) NMF programmes |
| F2 | `cell_states/Auto_pdo_flot_matched_highres_mp_trend_filter.R` | active upstream | F1 programmes; matched post-QC PDOs; 3CA/cell-cycle references | nMP=`round(total programmes/2)` MPs, UCell scores, paired trend audit and retained MP genes |
| F3 | `cell_states/Auto_pdo_flot_highres_enrichment_annotation.R` | active terminal/support | F2 retained genes and trend table | GO/Hallmark/3CA/developmental enrichment tables and figures |
| F4 | `cell_states/Auto_pdo_flot_highres_cluster_heatmap.R` | active terminal | F2 scores/metadata/labels; current centred state vector | sample-balanced MP correlation clustering plus individual-MP state-resolved response and absolute-score heatmaps; no manual grouping |
| F5 | `cell_states/Auto_pdo_flot_matched_survival_and_state_plots.R` | active terminal | F2 retained genes/labels; TCGA ESCA TPM/metadata | selected-MP Cox sensitivity table and volcano PDF |

F2 retains an MP when its sample mean and median UCell scores change in the
same direction in at least three of four treated-versus-untreated pairs. MP
display names use the best non-cell-cycle 3CA enrichment match. No manual
functional grouping is part of this route. See the dedicated methodology file.

## Other active metaprogram comparisons

| Script | Status | Purpose |
| :--- | :--- | :--- |
| `metaprograms/Auto_3CA_state_correlation_crossdata.R` | active terminal | state-resolved 3CA comparison between current PDO and scRef states |
| `metaprograms/Auto_3CA_pseudobulk_correlation_crossdata.R` | active terminal, input-blocked | scATLAS/OAC-PDO/OSCC-PDO 3CA score comparison; do not submit until the documented live GSE269447 organoid TPM files are staged |
| `metaprograms/Auto_parse_external_centred_mp_timecourse.R` | active terminal | projects current PDO/scRef centred MPs into Parse time-course data |
| `metaprograms/Auto_mp_chromosomal_mapping_pdo.R` | active terminal | maps final centred MP genes to hg38 positions |
| `metaprograms/Auto_scref_mps_flot_matched_trends.R` | active terminal | evaluates and visualizes mean and median expression trends of the 17 finalized scRef MPs across matched FLOT PDO pairs |

## CNV and demultiplex workflows

CNV analyses are optional diagnostics and do not define PDO cell types or
states.

| Workflow | Status | Scripts / notes |
| :--- | :--- | :--- |
| InferCNA/CNA diagnostics | active optional | `cnv/Auto_PDO_infercna.R`, `Auto_PDO_cnv_compare.R`, `Auto_PDO_cna_diagnostics_SUR1121_SUR1141.R` |
| Numbat core | active optional | `Auto_PDO_numbat_export_inputs.R`, `Auto_PDO_numbat_run_sample.R`, `Auto_PDO_numbat_conservative_recut.R`, `Auto_PDO_numbat_phylogeny_visualisation.R` and non-legacy wrappers |
| Indel visualization | active terminal | `cnv/Auto_indel_visualisation.R` and `Auto_run_indel_visualisation.sh` |
| WES subclone/CNA | active optional/terminal | all non-`legacy_` scripts under `cnv/wes_subclone/`; exact order and limitations are in `Auto_wes_subclone_methodology.md` |
| Demultiplex | active upstream/QC | non-legacy scripts under `demultiplex/`; Cell Ranger → Souporcell → reference assignment → count export → verification/publish |
| New-four-sample Strelka | active upstream | `demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_strelka.sh` before reference assignment |

Heavy Cell Ranger, Souporcell, Numbat, Sarek, FACETS, HATCHet, THetA2, and
PhyloWGS objects belong in ephemeral storage. Assignment keys, manifests needed
by another script, final segment tables, summaries, reports, and figures belong
in live storage. Active wrappers execute from the live repository even when an
explicit heavy work directory is under ephemeral storage.

## Shared support

| Script | Status | Purpose |
| :--- | :--- | :--- |
| `shared/Auto_pdo_analysis_config.R` | active support | current live/ephemeral paths, state/MP order, colors, thresholds |
| `shared/Auto_pdo_analysis_helpers.R` | active support | output tiers, state normalization, logging/cache helpers |

## Legacy inventory

Every `legacy_` file is provenance-only even if an older internal comment says
“active.” The authoritative registry block at the top and this map take
precedence.

- `cell_states/legacy_*.R`: uncentred state definition/finalization, old
  scRef-derived labeling, old marker/simulation/SCENIC, superseded matched-FLOT
  DGE/general response analyses, and old state-proportion comparisons. The
  centred high-resolution matched-FLOT MP scripts listed above are active and
  no longer carry the `legacy_` prefix.
- `clinical/legacy_*.R`: clinical/state/MP analyses tied to
  `Auto_PDO_final_states.rds` or `UCell_scores_filtered.rds`.
- `metaprograms/legacy_*.R`: uncentred nMP selection/scoring/correlation and
  exploratory NMF utilities.
- `enrichment/legacy_*.R`: uncentred or external scRef enrichment utilities;
  current centred enrichment is produced in step 04.
- `plotting/legacy_heatmap.R`: historical standalone NMF heatmap fragment.
- `cnv/legacy_*.R/.sh` and `cnv/wes_subclone/legacy_*`: superseded
  state/MP-coupled CNA summaries or earlier WES methods.
- `cell_states/Auto_drug_reversal/*`: legacy as a workflow because its inputs
  and manual panels use the superseded state/marker route; filenames are
  retained to avoid breaking the multi-script bundle.
- `trajectory/*`: legacy as a workflow because exported metadata, root-state
  definitions, pseudotime, and node plots use the superseded four/five-state
  route. Raw velocity concepts may be reusable only after a centred-state
  redesign.
- `cell_states/legacy_pdo_flot_matched_response.R` is superseded specifically
  by `centred_refined_pdo_flot_matched_response.R`.
- `cell_states/legacy_PDO_state_concordance.R` is superseded by
  `Auto_PDO_scRef_centred_state_concordance.R`.

Legacy wrappers are retained with their target names repaired where practical.
They must not be submitted as part of the current run order.

## Methodology index

Active complex workflows and their methodology:

- centred NMF/refinement:
  `methodology/metaprograms/centred/Auto_centred_metaprogram_refinement_methodology.md`;
- centred ordering/state/3CA/survival:
  `methodology/metaprograms/centred/Auto_centred_ordering_and_state_methodology.md`;
- current marker discovery:
  `methodology/cell_states/Auto_centred_five_state_marker_methodology.md`;
- current surface markers:
  `methodology/cell_states/Auto_five_state_surface_marker_methodology.md`;
- matched-FLOT:
  `methodology/cell_states/centred_refined_pdo_flot_matched_response_methodology.md`;
- current-centred clinical/AUCrel and response extremes:
  `methodology/clinical/Auto_centred_clinical_auc_response_associations_methodology.md`;
- centred high-resolution matched-FLOT MPs:
  `methodology/cell_states/Auto_pdo_flot_centred_highres_metaprogram_methodology.md`;
- treatment SCENIC:
  `methodology/cell_states/Auto_PDO_treatment_scenic_methodology.md`;
- personalised state-regulon RSS:
  `methodology/cell_states/Auto_PDO_personalised_state_regulon_rss_methodology.md`;
- oesophageal CRISPR selective dependency:
  `methodology/cell_states/Auto_crispr_selective_dependency_methodology.md`;
- PDO/scRef concordance:
  `methodology/cell_states/Auto_PDO_scRef_centred_state_concordance_methodology.md`;
- Parse and 3CA comparisons: the corresponding files under
  `methodology/metaprograms/`;
- CNV/WES, demultiplex/Strelka: the corresponding mirrored methodology folders.

Simple direct exports and plot-only scripts use `Methodology: none` in their
headers; no separate file is required.

## PBS and replot rules

All non-trivial runs use PBS and `#PBS -koed`. Current centred wrappers are at
the repository root. Long-running scripts should honor `PDO_FORCE_REBUILD=1`
and/or `PDO_REPLOT_ONLY=1` when their headers document those modes. Never
submit legacy scripts merely to refresh an old figure.
