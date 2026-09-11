####################
# Analysis registry:
#   Status: active terminal centred-state marker visualization
#   Script: analysis/cell_states/Auto_selected_state_marker_expression.R
#   Methodology: none (direct deterministic expression summaries and plots)
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Visualizes 12 prespecified markers for four canonical centred PDO states
#     across four matched untreated/FLOT-treated patient pairs. Provides state,
#     sample, sample-by-state, treatment-by-state, paired-delta, and per-cell views.
#   Inputs:
#     - PDOs_outs/PDOs_merged.rds
#     - PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     - analysis/shared/Auto_pdo_analysis_config.R
#     - analysis/shared/Auto_pdo_analysis_helpers.R
#   Outputs:
#     - PDOs_outs/Auto_matched_state_marker_expression/intermediate/Auto_selected_state_marker_per_cell_expression.rds
#     - PDOs_outs/Auto_matched_state_marker_expression/tables/Auto_selected_state_marker_*.csv*
#     - PDOs_outs/Auto_matched_state_marker_expression/figures/Auto_selected_state_marker_*.pdf
#     - PDOs_outs/Auto_matched_state_marker_expression/logs/Auto_selected_state_marker_expression_run_summary.txt
#   Downstream use: none; source tables and per-cell expression reproduce all figures.
#   Cache/replot behavior:
#     Reuses the live per-cell cache when it is newer than both primary inputs.
#     Set PDO_FORCE_REBUILD=1 to rebuild it or PDO_REPLOT_ONLY=1 to require it.
#   Run command: qsub Auto_run_selected_state_marker_expression.sh
#   Conda env: dmtcp
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(Cairo)
  library(grid)
})

####################
# setup and constants
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
script_path <- file.path(
  project_dir,
  "analysis",
  "cell_states",
  "Auto_selected_state_marker_expression.R"
)

source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_helpers.R"))

input_paths <- c(
  merged_pdo = file.path(PDO_OUTPUT_DIR, "PDOs_merged.rds"),
  state_vector = file.path(PDO_OUTPUT_DIR, PDO_PREFERRED_STATE_VECTOR)
)
pdo_require_files(input_paths, names(input_paths))

out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_matched_state_marker_expression")
output_tiers <- pdo_ensure_output_tiers(out_dir)
figures_dir <- output_tiers[["figures"]]
tables_dir <- output_tiers[["tables"]]
intermediate_dir <- output_tiers[["intermediate"]]

cache_policy <- pdo_cache_policy()
min_cells_sample_state <- 1L

target_states <- PDO_STATE_ORDER[seq_len(4)]
state_colors <- PDO_STATE_COLORS[target_states]
state_numbers <- setNames(seq_along(target_states), target_states)

patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
treatment_order <- c("Untreated", "Treated")
matched_samples <- as.vector(rbind(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
))
sample_patient <- setNames(
  sub("_(Treated|Untreated)_PDO$", "", matched_samples),
  matched_samples
)
sample_treatment <- setNames(
  ifelse(grepl("_Untreated_", matched_samples), "Untreated", "Treated"),
  matched_samples
)
patient_colors <- c(
  SUR1070 = "#4C78A8",
  SUR1072 = "#59A14F",
  SUR1090 = "#B07AA1",
  SUR1181 = "#F28E2B"
)
treatment_colors <- c(Untreated = "#D58B2D", Treated = "#374151")

marker_definition <- data.table(
  marker_state = rep(target_states, each = 3),
  gene = c(
    "MKI67", "TOP2A", "DLGAP5",
    "CEACAM5", "LGALS4", "TFF1",
    "PLCB4", "ROR1", "SOBP",
    "ID2", "CREB5", "LGALS1"
  )
)
marker_definition[, state_number := state_numbers[marker_state]]
marker_definition[, marker_rank := seq_len(.N)]
marker_order <- marker_definition$gene
marker_state_by_gene <- setNames(marker_definition$marker_state, marker_definition$gene)

cache_path <- file.path(
  intermediate_dir,
  "Auto_selected_state_marker_per_cell_expression.rds"
)
per_cell_csv <- file.path(
  tables_dir,
  "Auto_selected_state_marker_per_cell_expression.csv.gz"
)

message("=== Auto_selected_state_marker_expression.R ===")
message("Output directory: ", out_dir)
message("Matched samples: ", paste(matched_samples, collapse = ", "))

####################
# helpers
####################
row_zscore <- function(mat) {
  result <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )

  for (i in seq_len(nrow(mat))) {
    values <- mat[i, ]
    finite <- is.finite(values)
    if (!any(finite)) next
    value_sd <- stats::sd(values[finite])
    if (!is.finite(value_sd) || value_sd == 0) {
      result[i, finite] <- 0
    } else {
      result[i, finite] <- (values[finite] - mean(values[finite])) / value_sd
    }
  }
  result
}

natural_sample_order <- function(samples) {
  order_table <- data.table(sample = unique(as.character(samples)))
  order_table[, patient_number := suppressWarnings(
    as.integer(sub("^SUR([0-9]+).*$", "\\1", sample))
  )]
  order_table[is.na(patient_number), patient_number := .Machine$integer.max]
  order_table[, treatment_rank := fifelse(
    grepl("_Untreated_", sample),
    1L,
    fifelse(grepl("_Treated_", sample), 2L, 0L)
  )]
  setorder(order_table, patient_number, treatment_rank, sample)
  order_table$sample
}

sample_class <- function(samples) {
  factor(
    fifelse(
      grepl("_Untreated_", samples),
      "Untreated",
      fifelse(grepl("_Treated_", samples), "Treated", "Other")
    ),
    levels = c("Other", "Untreated", "Treated")
  )
}

build_matrix <- function(summary_table, column_name, value_name, column_order) {
  result <- matrix(
    NA_real_,
    nrow = length(marker_order),
    ncol = length(column_order),
    dimnames = list(marker_order, column_order)
  )
  row_index <- match(as.character(summary_table$gene), marker_order)
  column_index <- match(as.character(summary_table[[column_name]]), column_order)
  valid <- !is.na(row_index) & !is.na(column_index)
  if (anyDuplicated(paste(row_index[valid], column_index[valid], sep = "|"))) {
    stop("Summary table has duplicated gene-column combinations: ", column_name)
  }
  result[cbind(row_index[valid], column_index[valid])] <- summary_table[[value_name]][valid]
  result
}

marker_row_annotation <- function() {
  marker_state <- factor(
    marker_state_by_gene[marker_order],
    levels = target_states
  )
  names(marker_state) <- marker_order
  rowAnnotation(
    Marker_state = marker_state,
    col = list(Marker_state = state_colors),
    show_annotation_name = FALSE,
    show_legend = TRUE,
    simple_anno_size = unit(4, "mm")
  )
}

save_marker_heatmap <- function(
    matrix_unscaled,
    filename,
    title,
    subtitle,
    width,
    height,
    top_annotation = NULL,
    column_split = NULL,
    column_labels = colnames(matrix_unscaled),
    column_names_rot = 45,
    scale_rows = TRUE,
    heatmap_name = "Row Z-score",
    heatmap_colors = colorRamp2(
      c(-2, 0, 2),
      c("#1D4E89", "#F8F4EC", "#B22222")
    ),
    legend_title = "Relative\nexpression") {
  matrix_scaled <- if (scale_rows) row_zscore(matrix_unscaled) else matrix_unscaled
  marker_split <- factor(
    marker_state_by_gene[rownames(matrix_scaled)],
    levels = target_states
  )

  heatmap <- Heatmap(
    matrix_scaled,
    name = heatmap_name,
    col = heatmap_colors,
    na_col = "#D9D9D9",
    left_annotation = marker_row_annotation(),
    top_annotation = top_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = marker_split,
    column_split = column_split,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_gap = unit(2.5, "mm"),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    column_labels = column_labels,
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_rot = column_names_rot,
    border = TRUE,
    heatmap_legend_param = list(
      title = legend_title,
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )

  pdf(filename, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  pushViewport(viewport(
    layout = grid.layout(
      nrow = 2,
      ncol = 1,
      heights = unit(c(1.7, 1), c("cm", "null"))
    )
  ))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(
    title,
    x = unit(0.5, "npc"),
    y = unit(0.72, "npc"),
    gp = gpar(fontsize = 15, fontface = "bold")
  )
  grid.text(
    subtitle,
    x = unit(0.5, "npc"),
    y = unit(0.25, "npc"),
    gp = gpar(fontsize = 9)
  )
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(
    heatmap,
    newpage = FALSE,
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  popViewport(2)
  invisible(filename)
}

####################
# extract or reuse per-cell expression
####################
cache_exists <- file.exists(cache_path)
cache_is_current <- cache_exists &&
  file.info(cache_path)$mtime >= max(file.info(input_paths)$mtime)

if (cache_policy$replot_only && !cache_exists) {
  stop("PDO_REPLOT_ONLY=1 but the required per-cell cache is missing: ", cache_path)
}

reuse_cache <- cache_exists && (
  cache_policy$replot_only ||
    (!cache_policy$force_rebuild && cache_is_current)
)

if (reuse_cache) {
  message("Reusing live per-cell expression cache.")
  cache_object <- readRDS(cache_path)
  per_cell <- as.data.table(cache_object$per_cell)
  sample_order <- cache_object$sample_order
  if (!setequal(unique(as.character(per_cell$sample)), matched_samples)) {
    stop("Matched-pair cache does not contain exactly the eight required samples.")
  }
  rm(cache_object)
} else {
  message("Loading canonical merged PDO object and centred state vector.")
  pdos <- readRDS(input_paths[["merged_pdo"]])
  state_vector <- pdo_normalise_final_state_vector(
    readRDS(input_paths[["state_vector"]])
  )

  if (is.null(names(state_vector)) || anyDuplicated(names(state_vector))) {
    stop("Canonical state vector must have unique cell-barcode names.")
  }
  if (!"orig.ident" %in% colnames(pdos@meta.data)) {
    stop("PDOs_merged.rds is missing required metadata column: orig.ident")
  }

  sample_by_cell <- setNames(as.character(pdos$orig.ident), colnames(pdos))
  common_cells <- intersect(colnames(pdos), names(state_vector))
  missing_samples <- setdiff(matched_samples, unique(sample_by_cell[common_cells]))
  if (length(missing_samples) > 0) {
    stop("Matched sample(s) missing from canonical inputs: ", paste(missing_samples, collapse = ", "))
  }
  keep_cells <- common_cells[
    state_vector[common_cells] %in% target_states &
      sample_by_cell[common_cells] %in% matched_samples
  ]

  if (length(keep_cells) == 0) {
    stop("No cells remain after canonical-state filtering and sample exclusion.")
  }
  if (any(sample_by_cell[keep_cells] == PDO_EXCLUDED_SAMPLE)) {
    stop("Excluded sample remains in selected cells: ", PDO_EXCLUDED_SAMPLE)
  }

  message("Extracting log-normalized RNA expression for ", length(marker_order), " markers.")
  expression_all <- pdo_get_assay_matrix(pdos, assay = "RNA", layer = "data")
  missing_markers <- setdiff(marker_order, rownames(expression_all))
  if (length(missing_markers) > 0) {
    stop("Markers missing from the RNA data layer: ", paste(missing_markers, collapse = ", "))
  }

  expression_matrix <- as.matrix(
    expression_all[marker_order, keep_cells, drop = FALSE]
  )
  if (any(!is.finite(expression_matrix))) {
    stop("Non-finite values found in selected marker expression matrix.")
  }

  cell_expression_wide <- as.data.table(
    t(expression_matrix),
    keep.rownames = "cell"
  )
  cell_metadata <- data.table(
    cell = keep_cells,
    sample = sample_by_cell[keep_cells],
    state = state_vector[keep_cells],
    patient = sample_patient[sample_by_cell[keep_cells]],
    treatment = sample_treatment[sample_by_cell[keep_cells]]
  )
  per_cell <- melt(
    cell_expression_wide,
    id.vars = "cell",
    measure.vars = marker_order,
    variable.name = "gene",
    value.name = "expression",
    variable.factor = FALSE
  )
  per_cell <- cell_metadata[per_cell, on = "cell"]
  per_cell[, marker_state := marker_state_by_gene[gene]]
  per_cell[, detected := expression > 0]
  per_cell[, state := factor(state, levels = target_states)]
  per_cell[, marker_state := factor(marker_state, levels = target_states)]
  per_cell[, sample := factor(sample, levels = matched_samples)]
  per_cell[, patient := factor(patient, levels = patient_order)]
  per_cell[, treatment := factor(treatment, levels = treatment_order)]

  if (anyNA(per_cell$sample) || anyNA(per_cell$state) ||
      anyNA(per_cell$marker_state) || anyNA(per_cell$patient) ||
      anyNA(per_cell$treatment)) {
    stop("Missing sample/patient/treatment/state annotations after per-cell expression assembly.")
  }

  sample_order <- matched_samples
  setorder(per_cell, marker_state, gene, state, patient, treatment, cell)
  saveRDS(
    list(
      per_cell = per_cell,
      marker_definition = marker_definition,
      sample_order = sample_order,
      expression_scale = "RNA data layer (log-normalized expression)",
      matched_patients = patient_order,
      matched_samples = matched_samples,
      excluded_sample = PDO_EXCLUDED_SAMPLE,
      created = Sys.time()
    ),
    cache_path
  )

  rm(
    pdos,
    state_vector,
    expression_all,
    expression_matrix,
    cell_expression_wide,
    cell_metadata
  )
  invisible(gc())
}

if (!all(marker_order %in% unique(per_cell$gene))) {
  stop("Per-cell data do not contain all requested markers.")
}
if (any(as.character(per_cell$sample) == PDO_EXCLUDED_SAMPLE)) {
  stop("Excluded sample found in per-cell cache: ", PDO_EXCLUDED_SAMPLE)
}

fwrite(per_cell, per_cell_csv)
fwrite(
  marker_definition[, .(state_number, marker_state, gene, marker_rank)],
  file.path(tables_dir, "Auto_selected_state_marker_definition.csv")
)

####################
# auditable expression summaries
####################
message("Computing per-sample, per-state, and sample-by-state summaries.")

cell_counts <- unique(per_cell[, .(cell, patient, treatment, sample, state)])[ 
  , .(n_cells = .N),
  by = .(patient, treatment, sample, state)
]
cell_counts[, eligible_for_heatmap := n_cells >= min_cells_sample_state]
setorder(cell_counts, state, sample)

sample_state_summary <- per_cell[
  , .(
    mean_expression = mean(expression),
    median_expression = median(expression),
    pct_detected = 100 * mean(detected),
    n_cells = .N
  ),
  by = .(patient, treatment, sample, state, marker_state, gene)
]
sample_state_summary[, eligible_for_heatmap := n_cells >= min_cells_sample_state]
sample_state_summary[, sample_state_id := paste(state, sample, sep = " || ")]

sample_summary <- per_cell[
  , .(
    mean_expression = mean(expression),
    median_expression = median(expression),
    pct_detected = 100 * mean(detected),
    n_cells = .N,
    n_states_present = uniqueN(state)
  ),
  by = .(patient, treatment, sample, marker_state, gene)
]

state_summary <- sample_state_summary[
  eligible_for_heatmap == TRUE,
  .(
    mean_expression = mean(mean_expression),
    median_sample_mean = median(mean_expression),
    mean_pct_detected = mean(pct_detected),
    total_n_cells = sum(n_cells),
    n_samples = .N
  ),
  by = .(state, marker_state, gene)
]

treatment_state_summary <- sample_state_summary[
  eligible_for_heatmap == TRUE,
  .(
    mean_expression = mean(mean_expression),
    median_sample_mean = median(mean_expression),
    mean_pct_detected = mean(pct_detected),
    total_n_cells = sum(n_cells),
    n_patients = uniqueN(patient)
  ),
  by = .(state, treatment, marker_state, gene)
]
treatment_state_summary[, treatment_state_id := paste(state, treatment, sep = " || ")]

paired_delta_summary <- dcast(
  sample_state_summary[eligible_for_heatmap == TRUE],
  patient + state + marker_state + gene ~ treatment,
  value.var = "mean_expression"
)
for (treatment_name in treatment_order) {
  if (!treatment_name %in% colnames(paired_delta_summary)) {
    paired_delta_summary[, (treatment_name) := NA_real_]
  }
}
paired_delta_summary[, treated_minus_untreated := Treated - Untreated]
paired_delta_summary[, patient_state_id := paste(state, patient, sep = " || ")]

paired_delta_state_summary <- paired_delta_summary[
  is.finite(treated_minus_untreated),
  .(
    mean_delta = mean(treated_minus_untreated),
    median_delta = median(treated_minus_untreated),
    min_delta = min(treated_minus_untreated),
    max_delta = max(treated_minus_untreated),
    n_pairs = .N
  ),
  by = .(state, marker_state, gene)
]

setorder(sample_summary, marker_state, gene, sample)
setorder(state_summary, marker_state, gene, state)
setorder(sample_state_summary, marker_state, gene, state, sample)
setorder(treatment_state_summary, marker_state, gene, state, treatment)
setorder(paired_delta_summary, marker_state, gene, state, patient)
setorder(paired_delta_state_summary, marker_state, gene, state)

fwrite(
  cell_counts,
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_cell_counts.csv")
)
fwrite(
  sample_summary,
  file.path(tables_dir, "Auto_selected_state_marker_sample_summary.csv")
)
fwrite(
  state_summary,
  file.path(tables_dir, "Auto_selected_state_marker_state_summary.csv")
)
fwrite(
  sample_state_summary,
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_summary.csv")
)
fwrite(
  treatment_state_summary,
  file.path(tables_dir, "Auto_selected_state_marker_treatment_state_summary.csv")
)
fwrite(
  paired_delta_summary,
  file.path(tables_dir, "Auto_selected_state_marker_paired_delta_summary.csv")
)
fwrite(
  paired_delta_state_summary,
  file.path(tables_dir, "Auto_selected_state_marker_paired_delta_state_summary.csv")
)

####################
# heatmap matrices and source tables
####################
state_matrix <- build_matrix(
  state_summary,
  column_name = "state",
  value_name = "mean_expression",
  column_order = target_states
)

sample_matrix <- build_matrix(
  sample_summary,
  column_name = "sample",
  value_name = "mean_expression",
  column_order = sample_order
)

sample_state_plot <- copy(sample_state_summary)
sample_state_plot[
  eligible_for_heatmap == FALSE,
  heatmap_expression := NA_real_
]
sample_state_plot[
  eligible_for_heatmap == TRUE,
  heatmap_expression := mean_expression
]

combined_keys <- unlist(
  lapply(target_states, function(state_name) {
    paste(state_name, sample_order, sep = " || ")
  }),
  use.names = FALSE
)
combined_matrix <- build_matrix(
  sample_state_plot,
  column_name = "sample_state_id",
  value_name = "heatmap_expression",
  column_order = combined_keys
)

treatment_state_keys <- unlist(
  lapply(target_states, function(state_name) {
    paste(state_name, treatment_order, sep = " || ")
  }),
  use.names = FALSE
)
treatment_state_matrix <- build_matrix(
  treatment_state_summary,
  column_name = "treatment_state_id",
  value_name = "mean_expression",
  column_order = treatment_state_keys
)

patient_state_keys <- unlist(
  lapply(target_states, function(state_name) {
    paste(state_name, patient_order, sep = " || ")
  }),
  use.names = FALSE
)
paired_delta_matrix <- build_matrix(
  paired_delta_summary,
  column_name = "patient_state_id",
  value_name = "treated_minus_untreated",
  column_order = patient_state_keys
)

fwrite(
  as.data.table(state_matrix, keep.rownames = "gene"),
  file.path(tables_dir, "Auto_selected_state_marker_state_mean_expression_matrix.csv")
)
fwrite(
  as.data.table(sample_matrix, keep.rownames = "gene"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_mean_expression_matrix.csv")
)
fwrite(
  as.data.table(combined_matrix, keep.rownames = "gene"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_mean_expression_matrix.csv")
)
fwrite(
  as.data.table(treatment_state_matrix, keep.rownames = "gene"),
  file.path(tables_dir, "Auto_selected_state_marker_treatment_state_mean_expression_matrix.csv")
)
fwrite(
  as.data.table(paired_delta_matrix, keep.rownames = "gene"),
  file.path(tables_dir, "Auto_selected_state_marker_paired_delta_matrix.csv")
)

####################
# separated marker heatmaps
####################
message("Writing marker heatmaps.")

state_annotation <- HeatmapAnnotation(
  State = factor(target_states, levels = target_states),
  col = list(State = state_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)

state_heatmap_path <- file.path(
  figures_dir,
  "Auto_selected_state_marker_per_state_heatmap.pdf"
)
save_marker_heatmap(
  state_matrix,
  filename = state_heatmap_path,
  title = "Selected marker expression by canonical PDO state",
  subtitle = paste0(
    "Sample-balanced mean log-normalized RNA expression; rows scaled independently. ",
    "Marker blocks denote their prespecified state."
  ),
  width = 12,
  height = 8.5,
  top_annotation = state_annotation,
  column_split = factor(target_states, levels = target_states),
  column_labels = target_states,
  column_names_rot = 35
)

sample_annotation <- HeatmapAnnotation(
  Treatment = factor(sample_treatment[sample_order], levels = treatment_order),
  Patient = factor(sample_patient[sample_order], levels = patient_order),
  col = list(Treatment = treatment_colors, Patient = patient_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm")
)

sample_heatmap_path <- file.path(
  figures_dir,
  "Auto_selected_state_marker_per_sample_heatmap.pdf"
)
save_marker_heatmap(
  sample_matrix,
  filename = sample_heatmap_path,
  title = "Selected marker expression across four matched PDO pairs",
  subtitle = paste0(
    "Untreated and treated samples are adjacent within patient; means pool the four displayed states. ",
    "Rows are scaled independently."
  ),
  width = max(14, min(24, 7 + 0.42 * ncol(sample_matrix))),
  height = 9,
  top_annotation = sample_annotation,
  column_labels = sample_order,
  column_names_rot = 55
)

combined_state <- factor(
  rep(target_states, each = length(sample_order)),
  levels = target_states
)
combined_sample <- rep(sample_order, times = length(target_states))
combined_annotation <- HeatmapAnnotation(
  State = combined_state,
  Treatment = factor(sample_treatment[combined_sample], levels = treatment_order),
  Patient = factor(sample_patient[combined_sample], levels = patient_order),
  col = list(
    State = state_colors,
    Treatment = treatment_colors,
    Patient = patient_colors
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm")
)

combined_heatmap_path <- file.path(
  figures_dir,
  "Auto_selected_state_marker_per_sample_per_state_heatmap.pdf"
)
save_marker_heatmap(
  combined_matrix,
  filename = combined_heatmap_path,
  title = "Selected marker expression by sample and canonical PDO state",
  subtitle = paste0(
    "Untreated and treated samples are adjacent within each patient and state; ",
    "grey indicates an absent group. Rows are scaled independently."
  ),
  width = max(18, min(32, 8 + 0.24 * ncol(combined_matrix))),
  height = 9.5,
  top_annotation = combined_annotation,
  column_split = combined_state,
  column_labels = combined_sample,
  column_names_rot = 60
)

treatment_state_state <- factor(
  rep(target_states, each = length(treatment_order)),
  levels = target_states
)
treatment_state_treatment <- factor(
  rep(treatment_order, times = length(target_states)),
  levels = treatment_order
)
treatment_state_annotation <- HeatmapAnnotation(
  State = treatment_state_state,
  Treatment = treatment_state_treatment,
  col = list(State = state_colors, Treatment = treatment_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm")
)

treatment_state_heatmap_path <- file.path(
  figures_dir,
  "Auto_selected_state_marker_treated_vs_untreated_by_state_heatmap.pdf"
)
save_marker_heatmap(
  treatment_state_matrix,
  filename = treatment_state_heatmap_path,
  title = "Untreated versus treated marker expression within each PDO state",
  subtitle = paste0(
    "Each value is the equally weighted mean of four patient-specific sample-state means; ",
    "rows are scaled independently."
  ),
  width = 14,
  height = 9,
  top_annotation = treatment_state_annotation,
  column_split = treatment_state_state,
  column_labels = rep(treatment_order, times = length(target_states)),
  column_names_rot = 35
)

paired_delta_state <- factor(
  rep(target_states, each = length(patient_order)),
  levels = target_states
)
paired_delta_patient <- factor(
  rep(patient_order, times = length(target_states)),
  levels = patient_order
)
paired_delta_annotation <- HeatmapAnnotation(
  State = paired_delta_state,
  Patient = paired_delta_patient,
  col = list(State = state_colors, Patient = patient_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(4, "mm")
)
delta_limit <- max(abs(paired_delta_matrix), na.rm = TRUE)
if (!is.finite(delta_limit) || delta_limit == 0) delta_limit <- 1

paired_delta_heatmap_path <- file.path(
  figures_dir,
  "Auto_selected_state_marker_paired_treated_minus_untreated_heatmap.pdf"
)
save_marker_heatmap(
  paired_delta_matrix,
  filename = paired_delta_heatmap_path,
  title = "Paired treatment change in marker expression by patient and PDO state",
  subtitle = paste0(
    "Unscaled treated-minus-untreated log-normalized mean expression; ",
    "red indicates higher expression after FLOT and blue indicates lower expression."
  ),
  width = 17,
  height = 9,
  top_annotation = paired_delta_annotation,
  column_split = paired_delta_state,
  column_labels = rep(patient_order, times = length(target_states)),
  column_names_rot = 45,
  scale_rows = FALSE,
  heatmap_name = "Treated -\nUntreated",
  heatmap_colors = colorRamp2(
    c(-delta_limit, 0, delta_limit),
    c("#2166AC", "white", "#B2182B")
  ),
  legend_title = "Mean log-expression\ndifference"
)

####################
# per-cell distribution boxplots
####################
message("Writing per-cell expression boxplots.")

sample_medians <- per_cell[
  , .(sample_median = median(expression)),
  by = .(patient, treatment, sample, state, marker_state, gene)
]
sample_median_pairs <- dcast(
  sample_medians,
  patient + state + marker_state + gene ~ treatment,
  value.var = "sample_median"
)
sample_median_pairs[, treated_minus_untreated := Treated - Untreated]
fwrite(
  sample_medians,
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_medians.csv")
)
fwrite(
  sample_median_pairs,
  file.path(tables_dir, "Auto_selected_state_marker_paired_median_delta.csv")
)

boxplot_paths <- setNames(
  file.path(
    figures_dir,
    paste0(
      "Auto_selected_state_marker_per_cell_boxplot_marker_state_",
      seq_along(target_states),
      ".pdf"
    )
  ),
  target_states
)
boxplot_multipage_path <- file.path(
  figures_dir,
  "Auto_selected_state_marker_per_cell_boxplots.pdf"
)
boxplot_pages <- vector("list", length(target_states))
names(boxplot_pages) <- target_states

for (marker_state_name in target_states) {
  page_data <- per_cell[
    marker_state == marker_state_name & state == marker_state_name
  ]
  page_medians <- sample_medians[
    marker_state == marker_state_name & state == marker_state_name
  ]
  gene_levels <- marker_definition[
    marker_state == marker_state_name,
    gene
  ]
  page_data[, gene := factor(gene, levels = gene_levels)]
  page_medians[, gene := factor(gene, levels = gene_levels)]

  page_plot <- ggplot(
    page_data,
    aes(x = treatment, y = expression, fill = treatment)
  ) +
    geom_boxplot(
      width = 0.72,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.75
    ) +
    geom_point(
      data = page_medians,
      aes(
        x = treatment,
        y = sample_median,
        group = patient,
        colour = patient
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.1,
      stroke = 0.5,
      fill = "white",
      alpha = 0.85
    ) +
    geom_line(
      data = page_medians,
      aes(
        x = treatment,
        y = sample_median,
        group = patient,
        colour = patient
      ),
      inherit.aes = FALSE,
      linewidth = 0.45,
      alpha = 0.75
    ) +
    facet_wrap(~gene, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = treatment_colors, drop = FALSE) +
    scale_colour_manual(values = patient_colors, drop = FALSE) +
    labs(
      title = paste0(
        "State ",
        state_numbers[[marker_state_name]],
        ": ",
        marker_state_name,
        " markers"
      ),
      subtitle = paste0(
        "Within-state cell distributions by treatment; ",
        "points and lines connect matched patient sample medians"
      ),
      x = NULL,
      y = "Log-normalized RNA expression",
      fill = "Treatment",
      colour = "Patient"
    ) +
    pdo_theme_slide(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "italic", size = 13),
      panel.spacing.x = unit(8, "mm"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  boxplot_pages[[marker_state_name]] <- page_plot

  CairoPDF(
    file = boxplot_paths[[marker_state_name]],
    width = 13.333,
    height = 7.5,
    onefile = TRUE
  )
  print(page_plot)
  dev.off()
}

CairoPDF(
  file = boxplot_multipage_path,
  width = 13.333,
  height = 7.5,
  onefile = TRUE
)
for (marker_state_name in target_states) {
  print(boxplot_pages[[marker_state_name]])
}
dev.off()

####################
# validation and run record
####################
required_outputs <- c(
  cache_path,
  per_cell_csv,
  file.path(tables_dir, "Auto_selected_state_marker_definition.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_cell_counts.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_summary.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_state_summary.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_summary.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_treatment_state_summary.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_paired_delta_summary.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_paired_delta_state_summary.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_medians.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_paired_median_delta.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_state_mean_expression_matrix.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_mean_expression_matrix.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_sample_state_mean_expression_matrix.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_treatment_state_mean_expression_matrix.csv"),
  file.path(tables_dir, "Auto_selected_state_marker_paired_delta_matrix.csv"),
  state_heatmap_path,
  sample_heatmap_path,
  combined_heatmap_path,
  treatment_state_heatmap_path,
  paired_delta_heatmap_path,
  boxplot_multipage_path,
  unname(boxplot_paths)
)

missing_outputs <- required_outputs[
  !file.exists(required_outputs) | file.info(required_outputs)$size <= 0
]
if (length(missing_outputs) > 0) {
  stop("Missing or empty required output(s): ", paste(missing_outputs, collapse = ", "))
}

run_log <- pdo_write_run_summary(
  script = script_path,
  out_dir = out_dir,
  inputs = unname(input_paths),
  outputs = required_outputs,
  parameters = list(
    expression_scale = "RNA data layer (log-normalized expression)",
    matched_patients = paste(patient_order, collapse = "; "),
    matched_samples = paste(matched_samples, collapse = "; "),
    target_states = paste(target_states, collapse = "; "),
    markers = paste(marker_order, collapse = "; "),
    min_cells_sample_state = min_cells_sample_state,
    excluded_sample = PDO_EXCLUDED_SAMPLE,
    n_cells = uniqueN(per_cell$cell),
    n_patients = uniqueN(per_cell$patient),
    n_samples = length(sample_order),
    per_cell_cache_reused = reuse_cache
  ),
  cache = cache_policy,
  status = "completed"
)

message("Completed selected marker visualization.")
message("Cells: ", uniqueN(per_cell$cell), "; samples: ", length(sample_order))
message("Run log: ", run_log)
