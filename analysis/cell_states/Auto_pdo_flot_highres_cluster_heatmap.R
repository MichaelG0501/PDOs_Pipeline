####################
# Analysis registry:
#   Status: active terminal; centred high-resolution matched-FLOT heatmaps
#   Script: analysis/cell_states/Auto_pdo_flot_highres_cluster_heatmap.R
#   Methodology: analysis/methodology/cell_states/Auto_pdo_flot_centred_highres_metaprogram_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Displays retained centred high-resolution MPs directly across the current
#     five PDO states. MPs are ordered by selected treatment direction and
#     data-derived similarity; labels use each MP's best non-cell-cycle 3CA
#     enrichment match. Also computes a sample-balanced MP-by-MP UCell
#     correlation heatmap using Fisher-Z-averaged within-sample Spearman
#     correlations. No manual MP grouping or functional cluster is used.
#   Inputs:
#     - selected-MP UCell matrix, cell metadata, trend table and enrichment
#       labels from Auto_pdo_flot_matched_highres_mp_trend_filter.R
#     - live: centred_mp_refinement/centred_refined_noreg_states.rds
#   Outputs:
#     - live tables: state/patient/treatment MP summaries and paired deltas
#     - live tables: MP-by-MP mean correlation and correlation p-value matrices
#     - live figures: MP correlation, MP-by-state delta and absolute-score
#       heatmaps (PDF/PNG)
#   Downstream use: none; terminal source tables and figures.
#   Cache/replot behavior: deterministic replot from persistent live inputs.
#   Run command: use PBS in dmtcp after the trend filter.
#   Conda env: dmtcp
####################

library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

####################
# Persistent paths and inputs
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

out_dir <- file.path(PDO_LIVE_OUTS, "Auto_pdo_flot_centred_highres_metaprogram_trends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

config_path <- file.path(out_dir, "Auto_pdo_flot_highres_current_config.csv")
pdo_require_files(config_path)
config <- read.csv(config_path, check.names = FALSE, stringsAsFactors = FALSE)
n_mp <- as.integer(config$nMP[[1]])

input_paths <- c(
  ucell = file.path(out_dir, paste0("Auto_pdo_flot_highres_UCell_scores_nMP", n_mp, ".rds")),
  cell_meta = file.path(out_dir, paste0("Auto_pdo_flot_highres_cell_metadata_nMP", n_mp, ".rds")),
  trend = file.path(out_dir, paste0("Auto_pdo_flot_highres_trend_summary_nMP", n_mp, ".csv")),
  selected_genes = file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mp_genes_nMP", n_mp, ".rds")),
  labels = file.path(out_dir, paste0("Auto_pdo_flot_highres_top_3CA_noncellcycle_nMP", n_mp, ".csv")),
  states = file.path(PDO_LIVE_OUTS, PDO_PREFERRED_STATE_VECTOR)
)
pdo_require_files(input_paths, names(input_paths))

ucell <- readRDS(input_paths[["ucell"]])
cell_meta <- readRDS(input_paths[["cell_meta"]])
trend_summary <- read.csv(input_paths[["trend"]], check.names = FALSE, stringsAsFactors = FALSE)
selected_genes <- readRDS(input_paths[["selected_genes"]])
label_table <- read.csv(input_paths[["labels"]], check.names = FALSE, stringsAsFactors = FALSE)
state_vector <- pdo_normalise_final_state_vector(readRDS(input_paths[["states"]]))

if (is.null(names(state_vector))) stop("Current centred state vector is not cell-named.")
if (!all(c("cell", "patient", "treatment") %in% colnames(cell_meta))) {
  stop("High-resolution cell metadata lacks cell, patient, or treatment.")
}
if (!all(rownames(ucell) %in% cell_meta$cell)) {
  stop("High-resolution UCell cells are not fully represented in cell metadata.")
}

trend_summary$retained <- trend_summary$retained %in% c(TRUE, "TRUE")
retained_order <- trend_summary |>
  filter(retained, MP %in% names(selected_genes), MP %in% colnames(ucell)) |>
  mutate(direction_order = match(treatment_direction, c("increase", "decrease"))) |>
  arrange(direction_order, desc(pair_support_n), trend_p_value, MP) |>
  pull(MP)
if (length(retained_order) == 0L) stop("No retained high-resolution MPs are available.")

label_map <- setNames(label_table$top_3ca_noncc, label_table$MP)
label_map <- label_map[retained_order]
label_map[is.na(label_map) | !nzchar(label_map)] <- "no non-cell-cycle 3CA match"
display_labels <- setNames(
  paste0(retained_order, " | ", unname(label_map)),
  retained_order
)

common_cells <- intersect(intersect(rownames(ucell), cell_meta$cell), names(state_vector))
cell_df <- cell_meta[match(common_cells, cell_meta$cell), , drop = FALSE] |>
  mutate(
    state = state_vector[cell],
    state = factor(state, levels = PDO_STATE_ORDER),
    patient = factor(as.character(patient), levels = c("SUR1070", "SUR1072", "SUR1090", "SUR1181")),
    treatment = factor(as.character(treatment), levels = c("Untreated", "Treated"))
  ) |>
  filter(!is.na(state), !is.na(patient), !is.na(treatment))

if (nrow(cell_df) == 0L) {
  stop("No matched cells overlap the retained-MP scores and current five-state vector.")
}
if (anyDuplicated(cell_df$cell)) stop("Cell metadata contains duplicated cell identifiers.")

####################
# Sample-balanced MP-by-MP UCell correlation and clustering
####################
correlation_dir <- file.path(out_dir, "clustering")
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

correlation_scores <- as.matrix(ucell[cell_df$cell, retained_order, drop = FALSE])
sample_ids <- paste(as.character(cell_df$patient), as.character(cell_df$treatment), sep = " | ")
sample_levels <- unique(sample_ids)
n_retained <- length(retained_order)

sample_correlations <- array(
  NA_real_,
  dim = c(n_retained, n_retained, length(sample_levels)),
  dimnames = list(retained_order, retained_order, sample_levels)
)
for (sample_id in sample_levels) {
  sample_cells <- which(sample_ids == sample_id)
  if (length(sample_cells) < 10L) next
  sample_correlations[, , sample_id] <- suppressWarnings(
    cor(correlation_scores[sample_cells, , drop = FALSE], method = "spearman")
  )
}

fisher_z <- atanh(pmin(pmax(sample_correlations, -0.999), 0.999))
mean_rho <- matrix(
  NA_real_,
  nrow = n_retained,
  ncol = n_retained,
  dimnames = list(retained_order, retained_order)
)
correlation_p <- mean_rho
for (i in seq_len(n_retained)) {
  mean_rho[i, i] <- 1
  correlation_p[i, i] <- 0
  if (i == n_retained) next
  for (j in seq.int(i + 1L, n_retained)) {
    pair_z <- fisher_z[i, j, ]
    pair_z <- pair_z[is.finite(pair_z)]
    if (length(pair_z) < 3L) next
    pair_rho <- tanh(mean(pair_z))
    pair_test <- tryCatch(t.test(pair_z), error = function(e) NULL)
    pair_p <- if (is.null(pair_test)) NA_real_ else pair_test$p.value
    mean_rho[i, j] <- mean_rho[j, i] <- pair_rho
    correlation_p[i, j] <- correlation_p[j, i] <- pair_p
  }
}

mean_rho_csv <- file.path(
  correlation_dir,
  paste0("Auto_pdo_flot_highres_mp_mean_spearman_correlation_nMP", n_mp, ".csv")
)
correlation_p_csv <- file.path(
  correlation_dir,
  paste0("Auto_pdo_flot_highres_mp_correlation_p_values_nMP", n_mp, ".csv")
)
write.csv(mean_rho, mean_rho_csv, row.names = TRUE)
write.csv(correlation_p, correlation_p_csv, row.names = TRUE)

rho_for_clustering <- mean_rho
rho_for_clustering[!is.finite(rho_for_clustering)] <- 0
diag(rho_for_clustering) <- 1
rho_distance_matrix <- 1 - rho_for_clustering
rho_distance_matrix[rho_distance_matrix < 0] <- 0
rho_distance <- as.dist(rho_distance_matrix)
mp_clustering <- hclust(rho_distance, method = "average")

trend_map <- setNames(trend_summary$treatment_direction, trend_summary$MP)
trend_annotation <- factor(
  unname(trend_map[retained_order]),
  levels = c("increase", "decrease"),
  labels = c("Increased with treatment", "Decreased with treatment")
)
names(trend_annotation) <- retained_order
trend_colors <- c(
  "Increased with treatment" = "#B63E2F",
  "Decreased with treatment" = "#245F7B"
)

finite_off_diagonal <- abs(mean_rho[row(mean_rho) != col(mean_rho) & is.finite(mean_rho)])
rho_limit <- if (length(finite_off_diagonal) == 0L) {
  0.5
} else {
  min(0.95, max(0.4, unname(quantile(finite_off_diagonal, 0.98))))
}
rho_col_fun <- colorRamp2(c(-rho_limit, 0, rho_limit), c("#245F7B", "white", "#B63E2F"))

correlation_heatmap <- Heatmap(
  mean_rho,
  name = "Mean rho",
  col = rho_col_fun,
  cluster_rows = mp_clustering,
  cluster_columns = mp_clustering,
  row_labels = display_labels[rownames(mean_rho)],
  column_labels = display_labels[colnames(mean_rho)],
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 55,
  rect_gp = gpar(col = "white", lwd = 0.15),
  top_annotation = HeatmapAnnotation(
    FLOT_trend = trend_annotation[colnames(mean_rho)],
    col = list(FLOT_trend = trend_colors)
  ),
  left_annotation = rowAnnotation(
    FLOT_trend = trend_annotation[rownames(mean_rho)],
    col = list(FLOT_trend = trend_colors),
    show_legend = FALSE
  ),
  heatmap_legend_param = list(
    title = paste0("Fisher-Z mean\nSpearman rho\n(", length(sample_levels), " samples)")
  ),
  width = unit(11, "inch"),
  height = unit(11, "inch")
)

correlation_pdf <- file.path(
  correlation_dir,
  paste0("Auto_pdo_flot_highres_mp_correlation_heatmap_nMP", n_mp, ".pdf")
)
correlation_png <- sub("\\.pdf$", ".png", correlation_pdf)
pdf(correlation_pdf, width = 22, height = 22, useDingbats = FALSE)
draw(
  correlation_heatmap,
  column_title = "Retained centred high-resolution MP correlation across matched PDO samples",
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom"
)
dev.off()
png(correlation_png, width = 6600, height = 6600, res = 300)
draw(
  correlation_heatmap,
  column_title = "Retained centred high-resolution MP correlation across matched PDO samples",
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom"
)
dev.off()

####################
# MP-level state-resolved summaries
####################
score_long <- as.data.frame(ucell[cell_df$cell, retained_order, drop = FALSE]) |>
  rownames_to_column("cell") |>
  left_join(cell_df[, c("cell", "state", "patient", "treatment")], by = "cell") |>
  pivot_longer(all_of(retained_order), names_to = "MP", values_to = "ucell_score")

sample_state_summary <- score_long |>
  group_by(MP, state, patient, treatment) |>
  summarise(
    mean_score = mean(ucell_score, na.rm = TRUE),
    median_score = median(ucell_score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) |>
  complete(
    MP = retained_order,
    state = factor(PDO_STATE_ORDER, levels = PDO_STATE_ORDER),
    patient = factor(
      c("SUR1070", "SUR1072", "SUR1090", "SUR1181"),
      levels = c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
    ),
    treatment = factor(c("Untreated", "Treated"), levels = c("Untreated", "Treated")),
    fill = list(n_cells = 0L)
  )

paired_delta <- sample_state_summary |>
  select(MP, state, patient, treatment, mean_score, median_score, n_cells) |>
  pivot_wider(
    names_from = treatment,
    values_from = c(mean_score, median_score, n_cells)
  ) |>
  mutate(
    mean_delta = mean_score_Treated - mean_score_Untreated,
    median_delta = median_score_Treated - median_score_Untreated
  )

state_delta_summary <- paired_delta |>
  group_by(MP, state) |>
  summarise(
    n_pairs = sum(is.finite(mean_delta)),
    mean_delta = mean(mean_delta, na.rm = TRUE),
    median_delta = median(mean_delta, na.rm = TRUE),
    paired_wilcox_p = if (sum(is.finite(mean_delta)) >= 3L) {
      wilcox.test(mean_delta[is.finite(mean_delta)], mu = 0, exact = TRUE)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  group_by(state) |>
  mutate(paired_wilcox_fdr = p.adjust(paired_wilcox_p, method = "BH")) |>
  ungroup()

absolute_summary <- sample_state_summary |>
  group_by(MP, state, treatment) |>
  summarise(mean_score = mean(mean_score, na.rm = TRUE), .groups = "drop")

write.csv(
  sample_state_summary,
  file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_state_sample_summary_nMP", n_mp, ".csv")),
  row.names = FALSE
)
write.csv(
  paired_delta,
  file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_state_patient_deltas_nMP", n_mp, ".csv")),
  row.names = FALSE
)
write.csv(
  state_delta_summary,
  file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_state_delta_summary_nMP", n_mp, ".csv")),
  row.names = FALSE
)
write.csv(
  absolute_summary,
  file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_state_absolute_summary_nMP", n_mp, ".csv")),
  row.names = FALSE
)

####################
# MP-by-state treated-minus-untreated heatmap
####################
delta_mat <- state_delta_summary |>
  select(MP, state, mean_delta) |>
  pivot_wider(names_from = state, values_from = mean_delta) |>
  column_to_rownames("MP") |>
  as.matrix()
delta_mat <- delta_mat[retained_order, PDO_STATE_ORDER, drop = FALSE]

sig_mat <- state_delta_summary |>
  mutate(label = ifelse(is.finite(paired_wilcox_fdr), sprintf("q=%.2f", paired_wilcox_fdr), "")) |>
  select(MP, state, label) |>
  pivot_wider(names_from = state, values_from = label) |>
  column_to_rownames("MP") |>
  as.matrix()
sig_mat <- sig_mat[rownames(delta_mat), colnames(delta_mat), drop = FALSE]
sig_mat[is.na(sig_mat)] <- ""

direction_map <- setNames(trend_summary$treatment_direction, trend_summary$MP)
direction_split <- factor(
  unname(direction_map[rownames(delta_mat)]),
  levels = c("increase", "decrease"),
  labels = c("Increased in treated pairs", "Decreased in treated pairs")
)
finite_delta <- abs(delta_mat[is.finite(delta_mat)])
delta_limit <- if (length(finite_delta) == 0L) 0.01 else max(0.002, quantile(finite_delta, 0.95))
delta_col_fun <- colorRamp2(
  c(-delta_limit, 0, delta_limit),
  c("#245F7B", "white", "#B63E2F")
)

delta_heatmap <- Heatmap(
  delta_mat,
  name = "Mean delta",
  col = delta_col_fun,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  row_split = direction_split,
  row_labels = display_labels[rownames(delta_mat)],
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_rot = 35,
  top_annotation = HeatmapAnnotation(
    State = colnames(delta_mat),
    col = list(State = PDO_STATE_COLORS[PDO_STATE_ORDER]),
    show_annotation_name = FALSE
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.4f", delta_mat[i, j]), x, y - unit(1.2, "mm"), gp = gpar(fontsize = 6))
    if (nzchar(sig_mat[i, j])) {
      grid.text(sig_mat[i, j], x, y + unit(1.8, "mm"), gp = gpar(fontsize = 5.5))
    }
  },
  row_gap = unit(3, "mm"),
  width = unit(10, "cm"),
  height = unit(max(12, 0.42 * nrow(delta_mat)), "cm")
)

delta_pdf <- file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_state_delta_heatmap_nMP", n_mp, ".pdf"))
pdf(delta_pdf, width = 13, height = max(8, 0.22 * nrow(delta_mat) + 3), useDingbats = FALSE)
draw(delta_heatmap, column_title = "Centred high-resolution MP response within current PDO states")
dev.off()

delta_png <- sub("\\.pdf$", ".png", delta_pdf)
png(delta_png, width = 3900, height = max(2400, 70 * nrow(delta_mat) + 900), res = 300)
draw(delta_heatmap, column_title = "Centred high-resolution MP response within current PDO states")
dev.off()

####################
# Absolute untreated/treated activity heatmap
####################
absolute_wide <- absolute_summary |>
  mutate(state_treatment = paste(state, treatment, sep = " | ")) |>
  select(MP, state_treatment, mean_score) |>
  pivot_wider(names_from = state_treatment, values_from = mean_score) |>
  column_to_rownames("MP")

absolute_order <- unlist(lapply(PDO_STATE_ORDER, function(state_name) {
  paste(state_name, c("Untreated", "Treated"), sep = " | ")
}), use.names = FALSE)
absolute_mat <- as.matrix(absolute_wide[retained_order, absolute_order, drop = FALSE])
scaled_absolute <- t(scale(t(absolute_mat)))
scaled_absolute[!is.finite(scaled_absolute)] <- 0

absolute_limit <- max(1.5, quantile(abs(scaled_absolute), 0.98, na.rm = TRUE))
absolute_col_fun <- colorRamp2(
  c(-absolute_limit, 0, absolute_limit),
  c("#245F7B", "white", "#B63E2F")
)
column_states <- sub(" \\| (Untreated|Treated)$", "", colnames(absolute_mat))
column_treatments <- sub("^.* \\| ", "", colnames(absolute_mat))

absolute_heatmap <- Heatmap(
  scaled_absolute,
  name = "Row z-score",
  col = absolute_col_fun,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  row_split = direction_split,
  column_split = factor(column_states, levels = PDO_STATE_ORDER),
  cluster_column_slices = FALSE,
  row_labels = display_labels[rownames(scaled_absolute)],
  column_labels = column_treatments,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 35,
  top_annotation = HeatmapAnnotation(
    State = column_states,
    Treatment = column_treatments,
    col = list(
      State = PDO_STATE_COLORS[PDO_STATE_ORDER],
      Treatment = c(Untreated = "#B8B8B8", Treated = "#B43C3C")
    ),
    show_annotation_name = FALSE
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.4f", absolute_mat[i, j]), x, y, gp = gpar(fontsize = 5.5))
  },
  row_gap = unit(3, "mm"),
  column_gap = unit(3, "mm"),
  width = unit(18, "cm"),
  height = unit(max(12, 0.42 * nrow(scaled_absolute)), "cm")
)

absolute_pdf <- file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_state_absolute_heatmap_nMP", n_mp, ".pdf"))
pdf(absolute_pdf, width = 16, height = max(8, 0.22 * nrow(scaled_absolute) + 3), useDingbats = FALSE)
draw(absolute_heatmap, column_title = "Centred high-resolution MP activity by current PDO state")
dev.off()

absolute_png <- sub("\\.pdf$", ".png", absolute_pdf)
png(absolute_png, width = 4800, height = max(2400, 70 * nrow(scaled_absolute) + 900), res = 300)
draw(absolute_heatmap, column_title = "Centred high-resolution MP activity by current PDO state")
dev.off()

pdo_write_run_summary(
  script = "analysis/cell_states/Auto_pdo_flot_highres_cluster_heatmap.R",
  out_dir = out_dir,
  inputs = unname(input_paths),
  outputs = c(
    mean_rho_csv,
    correlation_p_csv,
    correlation_pdf,
    correlation_png,
    delta_pdf,
    delta_png,
    absolute_pdf,
    absolute_png
  ),
  parameters = list(
    n_mp = n_mp,
    retained_mp_n = length(retained_order),
    valid_cell_n = nrow(cell_df),
    correlation_samples = length(sample_levels),
    correlation_method = "within-sample Spearman; Fisher-Z mean",
    manual_grouping = FALSE,
    label_source = "best non-cell-cycle 3CA enrichment"
  )
)

message("Centred high-resolution MP state heatmaps completed.")
