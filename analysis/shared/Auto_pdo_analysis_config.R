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

PDO_PROJECT_DIR <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
PDO_ANALYSIS_DIR <- file.path(PDO_PROJECT_DIR, "analysis")
PDO_OUTPUT_DIR <- file.path(PDO_PROJECT_DIR, "PDOs_outs")
PDO_TEMP_DIR <- file.path(PDO_PROJECT_DIR, "temp")

PDO_PREFERRED_STATE_DEFINITION <- "Approach B, noreg"
PDO_PREFERRED_STATE_VECTOR <- "Auto_PDO_final_states.rds"
PDO_PREFERRED_PRE_FINAL_STATE_VECTOR <- "Auto_PDO_states_noreg.rds"
PDO_PREFERRED_MP_MATRIX <- "Auto_PDO_mp_adj_noreg.rds"
PDO_EXCLUDED_SAMPLE <- "SUR843T3_PDO"

PDO_OUTPUT_TIERS <- c("intermediate", "tables", "figures", "logs", "reports")

PDO_STATE_ORDER <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

PDO_STATE_ORDER_WITH_OPTIONAL <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation",
  "Unresolved",
  "Hybrid"
)

PDO_STATE_COLORS <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Immune Infiltrating" = "#A65628",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

PDO_MP_DESCRIPTIONS <- c(
  "MP6" = "G2M Cell Cycle",
  "MP7" = "DNA repair",
  "MP5" = "MYC-related Proliferation",
  "MP1" = "G2M checkpoint",
  "MP3" = "G1S Cell Cycle",
  "MP8" = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9" = "ECM Remodeling Epi.",
  "MP4" = "Intestinal Metaplasia"
)

PDO_MP_STATE_GROUPS <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "SMG-like Metaplasia" = c("MP8"),
  "Stress-adaptive" = c("MP10", "MP9")
)

PDO_CELL_CYCLE_MPS <- c("MP6", "MP7", "MP1", "MP3")

PDO_THRESHOLDS <- list(
  mp_min_silhouette = 0,
  mp_min_sample_coverage = 0.25,
  state_assignment_threshold = 0.5,
  hybrid_gap = 0.3,
  min_cells_state_sample = 20,
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
  sc_ref_pipeline = "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline",
  snseq_pipeline = "/rds/general/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline"
)

PDO_CACHE_ENV <- list(
  force_rebuild = "PDO_FORCE_REBUILD",
  replot_only = "PDO_REPLOT_ONLY"
)
