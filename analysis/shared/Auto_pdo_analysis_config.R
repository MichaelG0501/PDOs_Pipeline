####################
# Analysis registry:
#   Status: active shared configuration
#   Script: analysis/shared/Auto_pdo_analysis_config.R
#   Methodology: analysis/methodology/shared/shared_config_and_logging_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs: none
#   Outputs: shared constants for PDO downstream analysis scripts
####################

####################
# Central PDO analysis configuration
####################

PDO_PROJECT_DIR <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
PDO_EPHEMERAL_PROJECT_DIR <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
PDO_ANALYSIS_DIR <- file.path(PDO_PROJECT_DIR, "analysis")
PDO_OUTPUT_DIR <- file.path(PDO_PROJECT_DIR, "PDOs_outs")
PDO_TEMP_DIR <- file.path(PDO_PROJECT_DIR, "temp")
PDO_LIVE_OUTS <- PDO_OUTPUT_DIR
PDO_EPHEMERAL_OUTS <- file.path(PDO_EPHEMERAL_PROJECT_DIR, "PDOs_outs")

PDO_PREFERRED_STATE_DEFINITION <- "centred refined noreg Approach-B analogue"
PDO_PREFERRED_STATE_VECTOR <- file.path(
  "centred_mp_refinement", "centred_refined_noreg_states.rds"
)
PDO_PREFERRED_MP_MATRIX <- file.path(
  "centred_mp_refinement", "centred_refined_noreg_mp_adj.rds"
)
PDO_PREFERRED_GROUP_MAX <- file.path(
  "centred_mp_refinement", "centred_refined_noreg_group_max.rds"
)
PDO_PREFERRED_UCELL <- file.path(
  "centred_mp_refinement", "merged_refined_ucell_scores.rds"
)
PDO_PREFERRED_MP_GENES <- file.path(
  "centred_mp_refinement", "merged_refined_mp_genes.rds"
)
PDO_LEGACY_STATE_VECTOR <- "Auto_PDO_final_states.rds"
PDO_LEGACY_PRE_FINAL_STATE_VECTOR <- "Auto_PDO_states_noreg.rds"
PDO_LEGACY_MP_MATRIX <- "Auto_PDO_mp_adj_noreg.rds"
PDO_EXCLUDED_SAMPLE <- "SUR843T3_PDO"

PDO_OUTPUT_TIERS <- c("intermediate", "tables", "figures", "logs", "reports")

PDO_STATE_ORDER <- c(
  "Classic proliferation",
  "Columnar-to-intestinal",
  "Glandular differentiation",
  "Stress-adaptive",
  "ECM-remodelling",
  "Motile-cilia differentiation"
)

PDO_STATE_ORDER_WITH_OPTIONAL <- c(
  PDO_STATE_ORDER,
  "Unresolved",
  "Hybrid"
)

PDO_STATE_COLORS <- c(
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Motile-cilia differentiation" = "#F781BF",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

PDO_MP_DESCRIPTIONS <- c(
  "MP1" = "G2/M cell cycle",
  "MP2" = "G1/S cell cycle",
  "MP11" = "Single-nucleus-associated cell cycle",
  "MP3" = "Replication-dependent histones",
  "MP19+" = "MYC-associated proliferation",
  "MP14b" = "Proliferative epithelial plasticity",
  "MP13b" = "Metabolic-detox columnar epithelium",
  "MP5+" = "Inflammatory-reactive columnar epithelium",
  "MP12" = "KRAS-active columnar epithelium",
  "MP15" = "Intestinal metaplasia",
  "MP17+" = "Ciliated progenitor epithelium",
  "MP8+" = "Secretory-transport glandular epithelium",
  "MP16b" = "EMT/KRAS adaptive plasticity",
  "MP9" = "ECM-remodelling epithelium",
  "MP18" = "Motile-cilia differentiation"
)

PDO_MP_STATE_GROUPS <- list(
  "Classic proliferation" = c("MP19+"),
  "Columnar-to-intestinal" = c("MP14b", "MP13b", "MP5+", "MP12", "MP15"),
  "Glandular differentiation" = c("MP17+", "MP8+"),
  "Stress-adaptive" = c("MP16b"),
  "ECM-remodelling" = c("MP9"),
  "Motile-cilia differentiation" = c("MP18")
)

PDO_CELL_CYCLE_MPS <- c("MP11", "MP1", "MP2", "MP3")

PDO_THRESHOLDS <- list(
  parent_mp_min_silhouette = 0,
  mp_min_sample_n = 3L,
  mp_min_genes = 5L,
  state_assignment_threshold = 0.5,
  hybrid_gap = 0.3,
  marker_min_cells_state = 10L,
  marker_min_cells_rest = 10L,
  marker_specificity_gap = 0,
  marker_min_hit_sample_n = 1
)

PDO_METADATA_COLUMNS <- list(
  sample = "orig.ident",
  batch = "Batch",
  patient = "SUR",
  state = "state",
  final_state = "final_state",
  top_mp = "top_mp",
  treatment = "Treatment"
)

PDO_PLOT_DEFAULTS <- list(
  base_size = 13,
  axis_text_size = 11,
  legend_text_size = 10,
  legend_title_size = 11,
  strip_text_size = 12,
  heatmap_row_font_size = 9,
  heatmap_column_font_size = 9,
  pdf_width = 12,
  pdf_height = 8,
  slide_pdf_width = 13.333,
  slide_pdf_height = 7.5,
  dpi = 300
)

PDO_EXTERNAL_PATHS <- list(
  cell_cycle_genes = "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  three_ca_mps = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv",
  clinical_workbook = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx",
  developmental_reference_dir = "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage",
  sc_ref_pipeline = "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
)

PDO_CACHE_ENV <- list(
  force_rebuild = "PDO_FORCE_REBUILD",
  replot_only = "PDO_REPLOT_ONLY"
)
