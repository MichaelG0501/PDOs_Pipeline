####################
# Analysis registry:
#   Status: active canonical state-definition workflow
#   Script: analysis/metaprograms/centred/Auto_06_centred_refined_state_definition_noreg.R
#   Methodology: analysis/methodology/metaprograms/centred/Auto_centred_ordering_and_state_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
#   Description:
#     Defines per-cell states from centred refined MP UCell scores using only the
#     requested non-cell-cycle centred refined MP groups for PDOs. This is a noreg-only
#     Approach B analogue: MP scores are sample-centred and study-scaled directly,
#     with no cell-cycle regression branch and no cell-cycle MPs in the state
#     definition. Cell-cycle MPs are kept only as plotted heatmap rows.
#
#   Inputs:
#     - live: PDOs_outs/PDOs_merged.rds
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds
#     - live: /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#
#   Outputs (live: PDOs_outs/centred_mp_refinement/):
#     tables/centred_refined_noreg_state_counts.csv
#     tables/centred_refined_noreg_sample_state_abundance.csv
#     figures/centred_refined_noreg_state_heatmap.pdf
#     figures/centred_refined_noreg_state_proportion_withpie.pdf
#     figures/centred_refined_noreg_ccscore_boxplot.pdf
#     figures/centred_refined_noreg_sample_state_abundance.pdf
#   Outputs (live and ephemeral cache):
#     centred_refined_noreg_states.rds
#     centred_refined_noreg_mp_adj.rds
#     centred_refined_noreg_group_max.rds
#
#   Downstream use:
#     All current MP/state downstream analyses must use the live state vector and
#     normalized matrices above. Figure and CSV outputs are terminal summaries.
#
#   Cache/replot behavior:
#     Lightweight deterministic rebuild from the live merged UCell matrix. No
#     NMF, refinement, or UCell scoring is rerun.
#
#   Run command: qsub Auto_centred_06.sh
#   Conda env: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
})

# === Paths ===
live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
ephemeral_base <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

outdir_live <- file.path(live_base, "centred_mp_refinement")
outdir_ephemeral <- file.path(ephemeral_base, "centred_mp_refinement", "intermediate")

for (subdir in c("tables", "figures")) {
  dir.create(file.path(outdir_live, subdir), recursive = TRUE, showWarnings = FALSE)
}
dir.create(outdir_ephemeral, recursive = TRUE, showWarnings = FALSE)

seurat_file <- file.path(live_base, "PDOs_merged.rds")
ucell_file <- file.path(outdir_live, "merged_refined_ucell_scores.rds")
if (!file.exists(ucell_file)) {
  ucell_file <- file.path(outdir_ephemeral, "merged_refined_ucell_scores.rds")
}
cell_cycle_file <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"

required_inputs <- c(seurat_file, ucell_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))
}

cat("Loading centred refined MP state inputs...\n")
pdos_merged <- readRDS(seurat_file)
ucell_scores <- readRDS(ucell_file)
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")

####################
# Requested centred refined MP grouping. Cell-cycle MPs remain
# visible in heatmap rows, but only non-CC state_groups define state labels.
####################
cc_mps <- PDO_CELL_CYCLE_MPS
state_groups <- PDO_MP_STATE_GROUPS
excluded_mps <- character(0)

state_cols <- PDO_STATE_COLORS

mp_group_cols <- c(
  "Cell cycle" = "#6B7280",
  state_cols[names(state_groups)],
  "Excluded" = "grey80",
  "Other" = "grey70"
)

state_level_order <- c(names(state_groups), "Unresolved", "Hybrid")
plot_mp_order <- c(cc_mps, unlist(state_groups, use.names = FALSE), excluded_mps)

# Apply threshold filter logic: only use MPs that are actually retained (if any filtering happened in upstream)
# But since ucell_scores might contain all, we just check what's requested
missing_mps <- setdiff(plot_mp_order, colnames(ucell_scores))
if (length(missing_mps) > 0) {
  cat("WARNING: Merged centred refined UCell scores are missing requested MP(s): ", paste(missing_mps, collapse = ", "), "\n")
  # Filter them out
  cc_mps <- intersect(cc_mps, colnames(ucell_scores))
  state_groups <- lapply(state_groups, function(x) intersect(x, colnames(ucell_scores)))
  state_groups <- state_groups[lengths(state_groups) > 0]
  plot_mp_order <- c(cc_mps, unlist(state_groups, use.names = FALSE), excluded_mps)
}

common_cells <- intersect(rownames(ucell_scores), Cells(pdos_merged))
if (length(common_cells) == 0) {
  stop("No overlapping cells between centred refined UCell scores and PDOs_merged.rds")
}
pdos_merged <- pdos_merged[, common_cells]
ucell_scores <- ucell_scores[common_cells, plot_mp_order, drop = FALSE]

sample_var <- pdos_merged$orig.ident
names(sample_var) <- Cells(pdos_merged)
# For PDOs, study_var is typically 'Batch'
study_var <- pdos_merged$Batch
names(study_var) <- Cells(pdos_merged)

####################
# Helpers copied from the noreg branch of the Approach B workflow.
####################
z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

mp_desc_map <- setNames(
  paste(names(PDO_MP_DESCRIPTIONS), PDO_MP_DESCRIPTIONS, sep = "_"),
  names(PDO_MP_DESCRIPTIONS)
)

display_mp_label <- function(mp_vec, sep = ": ") {
  desc <- unname(mp_desc_map[mp_vec])
  desc[is.na(desc)] <- mp_vec[is.na(desc)]
  ifelse(desc == mp_vec, mp_vec, desc)
}

make_prop_data <- function(label_vec, sample_vec, all_labels) {
  df <- data.frame(
    orig.ident = as.character(sample_vec),
    label = as.character(label_vec),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$label), , drop = FALSE]
  df <- df[df$label %in% all_labels, , drop = FALSE]
  all_samples <- sort(unique(as.character(sample_vec)))
  counts <- df %>% count(orig.ident, label, name = "n")

  tidyr::expand_grid(
    orig.ident = all_samples,
    label = as.character(all_labels)
  ) %>%
    left_join(counts, by = c("orig.ident", "label")) %>%
    mutate(n = ifelse(is.na(n), 0L, as.integer(n))) %>%
    group_by(orig.ident) %>%
    mutate(pct = {
      total_n <- sum(n, na.rm = TRUE)
      if (total_n > 0) 100 * n / total_n else rep(0, dplyr::n())
    }) %>%
    ungroup()
}

plot_abundance <- function(prop_data, sample_order, col_map, title_text, totals_df) {
  plot_df <- prop_data %>%
    filter(orig.ident %in% sample_order) %>%
    mutate(orig.ident = factor(orig.ident, levels = sample_order))

  totals_plot <- totals_df %>%
    filter(orig.ident %in% sample_order) %>%
    mutate(orig.ident = factor(orig.ident, levels = sample_order))

  scale_factor <- max(totals_plot$total_n, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

  plot_df$label <- factor(plot_df$label, levels = rev(names(col_map)))

  ggplot(plot_df, aes(x = orig.ident, y = pct, fill = label)) +
    geom_col(width = 0.75) +
    geom_point(
      data = totals_plot,
      aes(x = orig.ident, y = total_n / scale_factor),
      color = "black",
      size = 1.5,
      shape = 18,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = totals_plot,
      aes(x = orig.ident, y = total_n / scale_factor, group = 1),
      color = "black",
      alpha = 0.4,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = col_map, breaks = names(col_map), drop = FALSE) +
    scale_y_continuous(
      name = "Proportion (%)",
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * scale_factor, name = "Total Cell Count (N)", labels = comma)
    ) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(title = title_text, x = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    guides(fill = guide_legend(ncol = 4, reverse = FALSE))
}
####################

cat("Defining centred refined noreg states...\n")
mp_adj_all <- z_normalise(as.matrix(ucell_scores[, plot_mp_order, drop = FALSE]), sample_var, study_var)

group_max <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_all))
  if (length(mps_avail) == 0) return(rep(NA_real_, nrow(mp_adj_all)))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_all[, mps_avail]))
  apply(mp_adj_all[, mps_avail, drop = FALSE], 1, max)
})
group_max <- as.matrix(group_max)
rownames(group_max) <- rownames(mp_adj_all)
group_max[!is.finite(group_max)] <- NA_real_

####################
# State-call thresholds inherited from the noreg Approach-B analogue. A cell is
# unresolved when its best normalized state-group score is < 0.5. Otherwise it
# is hybrid when the best-minus-second-best group score is < 0.3. These are
# fixed biological decision thresholds, not sample-specific fitted cutoffs.
####################
threshold <- PDO_THRESHOLDS$state_assignment_threshold
hybrid_gap <- PDO_THRESHOLDS$hybrid_gap

best_group_idx <- max.col(group_max, ties.method = "first")
best_group_val <- apply(group_max, 1, max, na.rm = TRUE)
state_vec <- names(state_groups)[best_group_idx]
state_vec[!is.finite(best_group_val) | best_group_val < threshold] <- "Unresolved"

sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
gap <- sorted_groups[, 1] - sorted_groups[, 2]
state_vec[(gap < hybrid_gap) & (state_vec != "Unresolved")] <- "Hybrid"
names(state_vec) <- rownames(group_max)

saveRDS(state_vec, file.path(outdir_ephemeral, "centred_refined_noreg_states.rds"))
saveRDS(mp_adj_all, file.path(outdir_ephemeral, "centred_refined_noreg_mp_adj.rds"))
saveRDS(group_max, file.path(outdir_ephemeral, "centred_refined_noreg_group_max.rds"))
####################
# Also save to live for downstream replotting (AGENTS.md storage policy)
saveRDS(state_vec, file.path(outdir_live, "centred_refined_noreg_states.rds"))
saveRDS(mp_adj_all, file.path(outdir_live, "centred_refined_noreg_mp_adj.rds"))
saveRDS(group_max, file.path(outdir_live, "centred_refined_noreg_group_max.rds"))
####################

state_counts <- data.frame(state = as.character(state_vec), stringsAsFactors = FALSE) %>%
  count(state, name = "cells") %>%
  mutate(pct = 100 * cells / sum(cells)) %>%
  arrange(match(state, state_level_order))
write.csv(
  state_counts,
  file.path(outdir_live, "tables", "centred_refined_noreg_state_counts.csv"),
  row.names = FALSE
)

####################
# Cell-cycle score for the per-cell heatmap.
####################
cc_score <- rep(NA_real_, length(state_vec))
names(cc_score) <- names(state_vec)
if (file.exists(cell_cycle_file)) {
  cell_cycle_genes <- read.csv(cell_cycle_file, header = TRUE, stringsAsFactors = FALSE)[, 1:3]
  cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
  expr_mat <- tryCatch(
    GetAssayData(pdos_merged, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(pdos_merged, assay = "RNA", slot = "data")
  )
  cc_consensus <- intersect(cc_consensus, rownames(expr_mat))
  if (length(cc_consensus) > 0) {
    cc_top <- names(sort(Matrix::rowMeans(expr_mat[cc_consensus, , drop = FALSE], na.rm = TRUE),
                         decreasing = TRUE))[seq_len(min(50, length(cc_consensus)))]
    cc_score <- Matrix::colMeans(expr_mat[cc_top, names(state_vec), drop = FALSE])
    cc_score <- as.numeric(cc_score)
    names(cc_score) <- names(state_vec)
  }
}

####################
# Per-cell heatmap, matching Approach B styling and sampling behavior.
####################
make_state_heatmap <- function(state_vec, mp_adj_all) {
  set.seed(42)
  max_cells_total <- 8000
  state_counts_local <- table(state_vec)
  state_fracs <- state_counts_local / sum(state_counts_local)
  cells_per_state <- pmax(round(state_fracs * max_cells_total), 20)
  state_cells <- split(names(state_vec), state_vec)
  cells_to_plot <- unlist(
    mapply(function(cells, n) sample(cells, min(length(cells), n)),
           state_cells,
           cells_per_state[names(state_cells)],
           SIMPLIFY = FALSE),
    use.names = FALSE
  )
  cells_to_plot <- unique(cells_to_plot)

  sub_scores <- t(mp_adj_all[cells_to_plot, plot_mp_order, drop = FALSE])
  rownames(sub_scores) <- display_mp_label(rownames(sub_scores))

  present_states <- intersect(state_level_order, unique(as.character(state_vec[cells_to_plot])))
  if (length(present_states) == 0) present_states <- unique(as.character(state_vec[cells_to_plot]))
  split_vec <- factor(as.character(state_vec[cells_to_plot]), levels = present_states)

  local_state_cols <- state_cols[present_states]
  local_state_cols[is.na(local_state_cols)] <- "grey80"

  study_vals <- pdos_merged@meta.data[cells_to_plot, "Batch"]
  study_levels <- unique(as.character(pdos_merged$Batch))
  study_cols <- setNames(hue_pal()(length(study_levels)), study_levels)

  cc_max <- max(cc_score[cells_to_plot], na.rm = TRUE)
  if (!is.finite(cc_max) || cc_max <= 0) cc_max <- 1

  col_ann <- HeatmapAnnotation(
    State = split_vec,
    CC_score = cc_score[cells_to_plot],
    Study = study_vals,
    col = list(
      State = local_state_cols,
      CC_score = colorRamp2(c(0, cc_max), c("white", "darkgreen")),
      Study = study_cols
    ),
    annotation_name_side = "left",
    show_legend = TRUE,
    na_col = "white"
  )

  mp_to_group <- rep("Other", length(plot_mp_order))
  names(mp_to_group) <- plot_mp_order
  mp_to_group[intersect(cc_mps, names(mp_to_group))] <- "Cell cycle"
  for (grp in names(state_groups)) {
    mp_to_group[intersect(state_groups[[grp]], names(mp_to_group))] <- grp
  }
  if (length(excluded_mps) > 0) mp_to_group[intersect(excluded_mps, names(mp_to_group))] <- "Excluded"
  mp_group_label <- mp_to_group
  names(mp_group_label) <- rownames(sub_scores)

  row_ann <- rowAnnotation(
    MP_group = factor(mp_group_label, levels = c("Cell cycle", names(state_groups), "Excluded", "Other")),
    col = list(MP_group = mp_group_cols),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      MP_group = list(title_gp = gpar(fontsize = 16, fontface = "bold"),
                      labels_gp = gpar(fontsize = 14))
    )
  )

  lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

  Heatmap(
    sub_scores,
    name = "Adj score",
    col = col_fun_sc,
    top_annotation = col_ann,
    left_annotation = row_ann,
    column_split = split_vec,
    column_order = (function() {
      col_order_list <- lapply(levels(split_vec), function(lvl) {
        idx <- which(as.character(split_vec) == lvl)
        if (length(idx) <= 1) return(idx)
        mat_lvl <- sub_scores[, idx, drop = FALSE]
        dcols <- dist(t(mat_lvl))
        hc <- hclust(dcols, method = "ward.D2")
        idx[hc$order]
      })
      full_ord <- unlist(col_order_list, use.names = FALSE)
      if (length(full_ord) != ncol(sub_scores) || !setequal(full_ord, seq_len(ncol(sub_scores)))) {
        return(seq_len(ncol(sub_scores)))
      }
      full_ord
    })(),
    column_gap = unit(1.5, "mm"),
    row_split = factor(mp_group_label, levels = c("Cell cycle", names(state_groups), "Excluded", "Other")),
    row_title = NULL,
    row_gap = unit(1.8, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
    row_names_max_width = unit(130, "mm"),
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 5,
    border = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  )
}

heatmap_pdf <- file.path(outdir_live, "figures", "centred_refined_noreg_state_heatmap.pdf")
cairo_pdf(heatmap_pdf, width = 24, height = 13)
draw(
  make_state_heatmap(state_vec, mp_adj_all),
  merge_legend = TRUE,
  padding = unit(c(12, 24, 12, 24), "mm")
)
grid.text(
  "PDO Centred refined states - noreg",
  x = unit(2, "mm"),
  y = unit(1, "npc") - unit(2, "mm"),
  just = c("left", "top"),
  gp = gpar(fontsize = 20, fontface = "bold")
)
dev.off()

####################
# Overall and per-study state proportions, plus CC-score boxplot.
####################
make_prop_and_pie <- function(state_vec) {
  prop_df <- data.frame(
    state = as.character(state_vec),
    study = as.character(pdos_merged@meta.data[names(state_vec), "Batch"]),
    stringsAsFactors = FALSE
  )
  overall <- prop_df %>% count(state) %>% mutate(study = "Total", pct = 100 * n / sum(n))
  per_study <- prop_df %>%
    count(study, state) %>%
    group_by(study) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  plot_df <- bind_rows(per_study, overall)
  plot_df$state <- factor(plot_df$state, levels = state_level_order)
  study_levels <- c(sort(unique(per_study$study)), "Total")
  plot_df$study <- factor(plot_df$study, levels = study_levels)
  plot_df$is_total <- factor(ifelse(plot_df$study == "Total", "Total", "Studies"),
                             levels = c("Studies", "Total"))

  p_bar <- ggplot(plot_df, aes(study, pct, fill = state)) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(pct > 3, sprintf("%.1f%%", pct), "")),
      position = position_stack(vjust = 0.5),
      size = 4.5
    ) +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    facet_grid(~ is_total, scales = "free_x", space = "free_x") +
    labs(title = "PDO Centred refined noreg: state proportions", x = NULL, y = "% of cells", fill = "State") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.title.y = element_text(size = 16, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      panel.spacing = unit(1, "lines")
    )

  pie_df <- overall %>% mutate(label = paste0(state, "\n", sprintf("%.1f%%", pct)))
  p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = ifelse(pct > 3, label, "")),
      position = position_stack(vjust = 0.5),
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    labs(title = "PDO Centred refined noreg: overall pie", fill = "State") +
    theme_void(base_size = 16) +
    theme(legend.position = "none",
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

  list(bar = p_bar, pie = p_pie, overall = overall)
}

make_cc_box <- function(state_vec) {
  cc_df <- data.frame(
    state = factor(as.character(state_vec), levels = state_level_order),
    cc_score = as.numeric(cc_score[names(state_vec)]),
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(state), is.finite(cc_score))

  ggplot(cc_df, aes(state, cc_score, fill = state)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.15, alpha = 0.2) +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    labs(title = "PDO Centred refined noreg: cell-cycle score by state", x = NULL, y = "Cell-cycle score") +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold")
    )
}

prop_result <- make_prop_and_pie(state_vec)
ggsave(
  file.path(outdir_live, "figures", "centred_refined_noreg_state_proportion_withpie.pdf"),
  prop_result$bar / prop_result$pie,
  width = 10,
  height = 18
)
ggsave(
  file.path(outdir_live, "figures", "centred_refined_noreg_ccscore_boxplot.pdf"),
  make_cc_box(state_vec),
  width = 11,
  height = 6
)

####################
# Per-sample abundance, state-only analogue of state_mp_sample_abundance.R.
####################
sample_by_cell <- pdos_merged$orig.ident
names(sample_by_cell) <- Cells(pdos_merged)

totals_df <- data.frame(
  orig.ident = as.character(sample_by_cell),
  stringsAsFactors = FALSE
) %>%
  count(orig.ident, name = "total_n")

state_df <- data.frame(
  cell = names(state_vec),
  state = as.character(state_vec),
  orig.ident = as.character(sample_by_cell[names(state_vec)]),
  stringsAsFactors = FALSE
)

target_states <- names(state_groups)
counts_long <- state_df %>%
  filter(state %in% target_states) %>%
  count(orig.ident, state, .drop = FALSE) %>%
  complete(orig.ident, state = target_states, fill = list(n = 0))

rank_df <- counts_long %>%
  group_by(orig.ident) %>%
  summarise(
    target_n = sum(n),
    geo_mean_score = exp(mean(log(n + 1))),
    .groups = "drop"
  ) %>%
  left_join(totals_df, by = "orig.ident") %>%
  arrange(desc(geo_mean_score), orig.ident) %>%
  mutate(rank = row_number())

diversity_order <- rank_df$orig.ident
diversity_order <- c(diversity_order, setdiff(sort(unique(sample_by_cell)), diversity_order))

study_map <- pdos_merged@meta.data %>%
  group_by(orig.ident) %>%
  summarise(study = first(Batch), .groups = "drop")
study_order <- study_map %>% arrange(study, orig.ident) %>% pull(orig.ident)

state_levels <- state_level_order
prop_states <- make_prop_data(state_vec, sample_by_cell[names(state_vec)], state_levels)
col_states <- state_cols[state_levels]

write.csv(
  prop_states,
  file.path(outdir_live, "tables", "centred_refined_noreg_sample_state_abundance.csv"),
  row.names = FALSE
)

n_samples <- length(unique(sample_by_cell))
pdf_w <- max(16, 0.15 * n_samples)
pdf_h <- 8

p_abund_div <- plot_abundance(
  prop_data = prop_states,
  sample_order = diversity_order,
  col_map = col_states,
  title_text = "PDO Centred refined states (noreg) | Sort: Diversity",
  totals_df = totals_df
)
p_abund_study <- plot_abundance(
  prop_data = prop_states,
  sample_order = study_order,
  col_map = col_states,
  title_text = "PDO Centred refined states (noreg) | Sort: Batch",
  totals_df = totals_df
)

sample_abundance_pdf <- file.path(outdir_live, "figures", "centred_refined_noreg_sample_state_abundance.pdf")
pdf(sample_abundance_pdf, width = pdf_w, height = pdf_h, onefile = TRUE)
print(p_abund_div)
print(p_abund_study)
dev.off()

cat("Saved:", heatmap_pdf, "\n")
cat("Saved:", sample_abundance_pdf, "\n")
cat("Saved: intermediate/centred_refined_noreg_states.rds\n")

cat("=== Auto_06 PDO Centred Refined State Definition Complete ===\n")
