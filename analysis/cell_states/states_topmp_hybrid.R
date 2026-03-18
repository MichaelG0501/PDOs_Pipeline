####################
# Auto_states_topmp_hybrid.R
# Adapted from scRef_Pipeline/analysis/cell_states/states_topmp_hybrid.R
# PDO-specific top-MP state assignment with grouped MPs, hybrid detection,
# CC annotations, and hybrid subclassification.
#
# Input:  PDOs_outs/PDOs_merged.rds
#         PDOs_outs/MP_outs_default.rds
#         PDOs_outs/UCell_scores_filtered.rds
#         /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
# Output: PDOs_outs/Auto_topmp_v2_states_A.rds
#         PDOs_outs/Auto_topmp_v2_states_B.rds
#         PDOs_outs/Auto_topmp_v2_mp_adj.rds
#         PDOs_outs/Auto_topmp_v2_heatmap_A.pdf
#         PDOs_outs/Auto_topmp_v2_heatmap_B.pdf
#         PDOs_outs/Auto_topmp_v2_proportion_A.pdf
#         PDOs_outs/Auto_topmp_v2_proportion_B.pdf
#         PDOs_outs/Auto_topmp_v2_mean_heatmap.pdf
#         PDOs_outs/Auto_topmp_v2_summary.csv
####################
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

cat("=== Auto states top-MP v2 ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

# ============================================================================
# 1. Load data
# ============================================================================
cat("Loading data...\n")
tmdata_all <- readRDS("PDOs_merged.rds")
geneNMF.metaprograms <- readRDS("MP_outs_default.rds") # UPDATE: Change to optimal nMP result after running Auto_find_optimal_nmf.R
ucell_scores <- readRDS("UCell_scores_filtered.rds") # UPDATE: Verify UCell output filename after geneNMF pipeline

cell_cycle_genes <- read.csv("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
                             header = TRUE, stringsAsFactors = FALSE)[, 1:3]
####################
# NOTE: Cancer signature scoring removed — not applicable to PDO-only analysis
####################

# ============================================================================
# 2. Silhouette filtering
# ============================================================================
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  cat("Removing low-quality MPs (silhouette < 0):", bad_mp_names, "\n")
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
retained_mps <- names(mp.genes)
cat("Retained MPs:", retained_mps, "\n")

# Derive canonical MP ordering from geneNMF program tree
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)
cat("Canonical MP tree order:", mp_tree_order_names, "\n")

# ============================================================================
# 3. MP annotations, CC identification, group definitions
# ============================================================================
####################
# UPDATE AFTER geneNMF: Fill in MP descriptions based on enrichment results
# Run Auto_enrichment_annotation.R first, then annotate each MP
####################
mp_descriptions <- c(
  "MP6"  = "G2M Cell Cycle",
  "MP7"  = "DNA repair",
  "MP5"  = "MYC-related Proliferation",
  "MP1"  = "G2M checkpoint",
  "MP3"  = "G1S Cell Cycle",
  "MP8"  = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9"  = "ECM Remodeling Epi.",
  "MP4"  = "Intestinal Metaplasia"
)
# TEMPORARY: auto-generate names if not yet annotated
if (length(mp_descriptions) == 0) {
  mp_descriptions <- setNames(retained_mps, retained_mps)
}

####################
# UPDATE AFTER geneNMF: Identify which MPs correspond to cell cycle
# Look for MPs enriched in G2M, S-phase, G1/S genes
####################
cc_mps <- c("MP6", "MP7", "MP1", "MP3")
non_cc_mps <- setdiff(retained_mps, cc_mps)

####################
# UPDATE AFTER geneNMF: Define state groupings based on MP biology
# Example structure (adapt to actual PDO MPs):
####################
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "Basal to Intest. Meta" = c("MP4")
)

stopifnot(all(unlist(state_groups) %in% non_cc_mps))

# ============================================================================
# 4. Align UCell scores to Seurat object
# ============================================================================
common_cells <- intersect(rownames(ucell_scores), Cells(tmdata_all))
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
tmdata_all <- tmdata_all[, common_cells]
cat(sprintf("Aligned %d cells\n", length(common_cells)))

retained_in_ucell <- intersect(retained_mps, colnames(ucell_scores))
cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
non_cc_in_ucell <- intersect(non_cc_mps, colnames(ucell_scores))
cat("Non-CC MPs available in UCell:", non_cc_in_ucell, "\n")

# ============================================================================
# 5. Compute CC scores (top 50 genes)
# ============================================================================
cat("Computing CC scores...\n")
cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(tmdata_all))
gene_means_cc <- rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE)
cc_top50 <- names(sort(gene_means_cc, decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))
names(cc_score) <- colnames(tmdata_all)

# ============================================================================
# 6. Regress cell-cycle signal from non-CC MPs
# ============================================================================
cat("Regressing cell-cycle signal...\n")
ucell_mat <- as.matrix(ucell_scores[, retained_in_ucell, drop = FALSE])
X_cc <- ucell_mat[, cc_in_ucell, drop = FALSE]
Y_other <- ucell_mat[, non_cc_in_ucell, drop = FALSE]

X <- cbind(Intercept = 1, X_cc)
XtX_inv <- solve(crossprod(X))
B <- XtX_inv %*% crossprod(X, Y_other)
Y_hat <- X %*% B
Y_resid <- Y_other - Y_hat
cat("Residualization complete.\n")

# ============================================================================
# 7. Z-score normalisation (centre per sample, scale by within-batch SD)
# ============================================================================
z_normalise <- function(mat, sample_var, batch_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.batch <- batch_var[rownames(mat)]

  batch_sd <- clust_df %>%
    group_by(.batch) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".batch") %>%
    as.matrix()

  batch_sd[is.na(batch_sd) | batch_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell

  for (mp in colnames(mp_adj)) {
    cell_batches <- clust_centered$.batch
    mp_adj[, mp] <- mp_adj[, mp] / batch_sd[cell_batches, mp]
  }

  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

cat("Z-normalising non-CC MPs...\n")
sample_var <- tmdata_all$orig.ident
batch_var <- tmdata_all$Batch
mp_adj_noncc <- z_normalise(Y_resid, sample_var, batch_var)

cat("Z-normalising CC MPs (no regression)...\n")
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, batch_var)

mp_adj_all <- cbind(mp_adj_noncc, mp_adj_cc)

# ============================================================================
# 8. Group-level state assignment
# ============================================================================
cat("Assigning group-level states...\n")
group_max <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_noncc))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc[, mps_avail]))
  apply(mp_adj_noncc[, mps_avail, drop = FALSE], 1, max)
})
group_max <- as.matrix(group_max)

THRESHOLD <- 0.5
best_group_idx <- max.col(group_max, ties.method = "first")
best_group_val <- apply(group_max, 1, max)
group_names <- names(state_groups)
base_state <- group_names[best_group_idx]
base_state[best_group_val < THRESHOLD] <- "Unresolved"

####################
# Hybrid Approach A (updated threshold)
####################
HYBRID_THR_A <- 1
sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
top1_val <- sorted_groups[, 1]
top2_val <- sorted_groups[, 2]

state_A <- base_state
hybrid_mask_A <- (top1_val > HYBRID_THR_A) & (top2_val > HYBRID_THR_A)
state_A[hybrid_mask_A & base_state != "Unresolved"] <- "Hybrid"

####################
# Hybrid Approach B (updated gap threshold)
####################
HYBRID_GAP_B <- 0.3
gap <- top1_val - top2_val
hybrid_mask_B <- (gap < HYBRID_GAP_B) & (base_state != "Unresolved")
state_B <- base_state
state_B[hybrid_mask_B] <- "Hybrid"

names(state_A) <- rownames(mp_adj_noncc)
names(state_B) <- rownames(mp_adj_noncc)

# ============================================================================
# 9. State order from mp_tree_order
# ============================================================================
group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))
####################
# State order: tree-ordered groups, Unresolved second last, Hybrid last
####################
state_level_order <- c(ordered_group_names, "Unresolved", "Hybrid")

# ============================================================================
# 10. Colour palettes
# ============================================================================
group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "SMG-like Metaplasia"   = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "Basal to Intest. Meta" = "#FF7F00",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

  batch_cols <- setNames(
    DiscretePalette(length(unique(tmdata_all$Batch)), palette = "polychrome"),
    unique(tmdata_all$Batch)
  )

# ============================================================================
# 11. Plotting helpers
# ============================================================================
plot_for_state <- function(state_vec, suffix) {
  cat(sprintf("Generating plots for approach %s...\n", suffix))
  meta <- tmdata_all@meta.data

  ####################
  # Build comprehensive state occurrence/diversity score
  ####################
  state_df_full <- data.frame(
    cell = names(state_vec),
    state = as.character(state_vec),
    sample = as.character(meta[names(state_vec), "orig.ident"]),
    batch = as.character(meta[names(state_vec), "Batch"]),
    stringsAsFactors = FALSE
  )
  total_samples <- length(unique(state_df_full$sample))
  total_batches <- length(unique(state_df_full$batch))
  state_div_df <- state_df_full %>%
    dplyr::group_by(state) %>%
    dplyr::summarise(
      sample_cov = dplyr::n_distinct(sample) / max(total_samples, 1),
      batch_cov = dplyr::n_distinct(batch) / max(total_batches, 1),
      diversity_score = 0.5 * sample_cov + 0.5 * batch_cov,
      .groups = "drop"
    )
  state_div_map <- setNames(state_div_df$diversity_score, state_div_df$state)

  # Proportional downsampling
  set.seed(42)
  MAX_CELLS_TOTAL <- 8000
  state_counts <- table(state_vec)
  state_fracs <- state_counts / sum(state_counts)
  cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)
  state_cells <- split(names(state_vec), state_vec)
  cells_to_plot <- unlist(
    mapply(function(cells, n) sample(cells, min(length(cells), n)),
           state_cells,
           cells_per_state[names(state_cells)],
           SIMPLIFY = FALSE),
    use.names = FALSE
  )

  # Heatmap matrix: rows = MPs, columns = cells
  sub_scores_orig <- t(mp_adj_all[cells_to_plot, , drop = FALSE])
  ####################
  # MP row order: three cell-cycle MPs first, then a gap, then non-CC MPs in tree order
  ####################
  cc_block_order <- cc_mps[cc_mps %in% rownames(sub_scores_orig)]
  non_cc_block_order <- mp_tree_order_names[
    mp_tree_order_names %in% rownames(sub_scores_orig) &
      !(mp_tree_order_names %in% cc_mps)
  ]
  mp_row_order <- c(cc_block_order, non_cc_block_order)
  sub_scores <- sub_scores_orig[mp_row_order, , drop = FALSE]

  mp_label_map <- mp_descriptions
  mp_label_map[setdiff(rownames(sub_scores), names(mp_label_map))] <- setdiff(rownames(sub_scores), names(mp_label_map))
  rownames(sub_scores) <- mp_label_map[rownames(sub_scores)]

  # Create a column split that only includes states present in the sampled cells
  present_states <- intersect(state_level_order, unique(as.character(state_vec[cells_to_plot])))
  if (length(present_states) == 0) present_states <- unique(as.character(state_vec[cells_to_plot]))
  split_vec <- factor(as.character(state_vec[cells_to_plot]), levels = present_states)
  batch_vals <- meta[cells_to_plot, "Batch"]

  max_cc <- max(cc_score[cells_to_plot], na.rm = TRUE)

  # map diversity values onto cells; missing -> 0
  state_div_vals <- state_div_map[as.character(split_vec)]
  state_div_vals[is.na(state_div_vals)] <- 0
  names(state_div_vals) <- cells_to_plot

  # Build a local palette for only the present states to avoid mismatch length errors
  local_group_cols <- group_cols[names(group_cols) %in% present_states]
  # ensure every present state has a colour (fallback grey)
  for (st in present_states) if (!st %in% names(local_group_cols)) local_group_cols[st] <- "grey80"
  # preserve order
  local_group_cols <- local_group_cols[present_states]

  col_list_ann <- list(
    State = local_group_cols,
    CC_score = colorRamp2(c(0, max_cc), c("white", "darkgreen")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Batch = batch_cols
  )

  col_ann <- HeatmapAnnotation(
    State = split_vec,
    CC_score = cc_score[cells_to_plot],
    Diversity = state_div_vals,
    Batch = batch_vals,
    col = col_list_ann,
    annotation_name_side = "left",
    show_legend = TRUE,
    na_col = "white"
  )

  ####################
  # Row grouping annotation: related MP groups (replacing CC/non-CC MP_type)
  ####################
  mp_to_group <- rep("Other", length(mp_row_order))
  names(mp_to_group) <- mp_row_order
  mp_to_group[cc_mps[cc_mps %in% names(mp_to_group)]] <- "Cell_cycle"
  for (grp in names(state_groups)) {
    grp_mps <- intersect(state_groups[[grp]], names(mp_to_group))
    mp_to_group[grp_mps] <- grp
  }
  group_colors_row <- c(group_cols[ordered_group_names], Cell_cycle = "gold", Other = "grey70")
  mp_group_label <- mp_to_group
  names(mp_group_label) <- rownames(sub_scores)

  row_ann <- rowAnnotation(
    MP_group = factor(mp_group_label, levels = c("Cell_cycle", ordered_group_names, "Other")),
    col = list(MP_group = group_colors_row),
    show_annotation_name = FALSE
  )

  lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
  col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

  ht_sc <- Heatmap(
    sub_scores,
    name = "Adj score",
    col = col_fun_sc,
    top_annotation = col_ann,
    left_annotation = row_ann,
    column_split = split_vec,
    # enforce strict column-split order: compute an ordered index vector per split
    column_order = (function() {
      cols_all <- colnames(sub_scores)
      col_order_list <- lapply(levels(split_vec), function(lvl) {
        idx <- which(as.character(split_vec) == lvl)
        if (length(idx) <= 1) return(idx)
        mat_lvl <- sub_scores[, idx, drop = FALSE]
        dcols <- dist(t(mat_lvl))
        hc <- hclust(dcols, method = "ward.D2")
        idx[hc$order]
      })
      # flatten to a single ordering vector (ComplexHeatmap accepts a vector)
      full_ord <- unlist(col_order_list, use.names = FALSE)
      # ensure we have all columns exactly once; if not, fall back to default order
      if (length(full_ord) != ncol(sub_scores) || !setequal(full_ord, seq_len(ncol(sub_scores)))) {
        warning('Computed column_order incomplete or invalid; using default column order')
        return(seq_len(ncol(sub_scores)))
      }
      full_ord
    })(),
    column_gap = unit(1.5, "mm"),
    row_split = factor(ifelse(mp_row_order %in% cc_mps, "Cell_cycle_MPs", "Other_MPs"),
                       levels = c("Cell_cycle_MPs", "Other_MPs")),
    row_gap = unit(2.5, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    use_raster = TRUE,
    raster_quality = 5,
    border = FALSE,
    rect_gp = gpar(col = NA)
  )

  pdf(sprintf("Auto_topmp_v2_heatmap_%s.pdf", suffix), width = 18,
      height = max(7, length(mp_row_order) * 0.5), useDingbats = FALSE)
  draw(ht_sc, merge_legend = TRUE)
  dev.off()

  # Proportion barplot (overall + per-batch)
  overall_df <- data.frame(state = state_vec) %>%
    count(state) %>%
    mutate(pct = 100 * n / sum(n), batch = "Overall")

  per_batch_df <- data.frame(state = state_vec, batch = meta$Batch) %>%
    count(state, batch) %>%
    group_by(batch) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  plot_df <- bind_rows(overall_df, per_batch_df)
  plot_df$state <- factor(plot_df$state, levels = state_level_order)
  plot_df$batch <- factor(plot_df$batch, levels = c("Overall", sort(unique(meta$Batch))))

  p_prop <- ggplot(plot_df, aes(x = batch, y = pct, fill = state)) +
    geom_col(colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = group_cols) +
    labs(title = sprintf("Cell state proportions (Top-MP v2, %s)", suffix),
         x = NULL, y = "% of cells", fill = "State") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(sprintf("Auto_topmp_v2_proportion_%s.pdf", suffix), p_prop, width = 14, height = 7)
}

mean_heatmap <- function(state_vec, title_suffix) {
  mean_scores_df <- as.data.frame(mp_adj_noncc)
  mean_scores_df$state <- state_vec
  mean_mat <- mean_scores_df %>%
    group_by(state) %>%
    summarise(across(all_of(colnames(mp_adj_noncc)), mean, na.rm = TRUE), .groups = "drop") %>%
    tibble::column_to_rownames("state") %>%
    as.matrix()

  ####################
  # Enforce state order from mp_tree_order with Unresolved then Hybrid at end
  ####################
  ordered_states_present <- state_level_order[state_level_order %in% rownames(mean_mat)]
  mean_mat <- mean_mat[ordered_states_present, , drop = FALSE]

  non_cc_tree_order <- mp_tree_order_names[mp_tree_order_names %in% colnames(mp_adj_noncc)]
  mean_mat <- mean_mat[, non_cc_tree_order, drop = FALSE]
  colnames(mean_mat) <- mp_descriptions[colnames(mean_mat)]

  cell_counts <- table(state_vec)
  row_anno <- rowAnnotation(
    Cells = anno_barplot(as.numeric(cell_counts[rownames(mean_mat)]),
                         gp = gpar(fill = "steelblue"),
                         width = unit(2, "cm"))
  )

  col_fun_z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

  Heatmap(
    mean_mat,
    name = "Mean Z-score",
    col = col_fun_z,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    clustering_method_rows = "ward.D2",
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10, fontface = "italic"),
    column_names_rot = 45,
    column_title = sprintf("Mean MP activity per state (%s)", title_suffix),
    rect_gp = gpar(col = "grey80", lwd = 0.3),
    right_annotation = row_anno,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", mean_mat[i, j]), x, y, gp = gpar(fontsize = 7))
    }
  )
}

# ============================================================================
# 12. Plotting for approaches A and B
# ============================================================================
plot_for_state(state_A, "A")
plot_for_state(state_B, "B")

pdf("Auto_topmp_v2_mean_heatmap.pdf", width = 14, height = 8, useDingbats = FALSE)
draw(mean_heatmap(state_A, "Approach A"), merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right")
draw(mean_heatmap(state_B, "Approach B"), merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# ============================================================================
# 13. Save outputs
# ============================================================================
saveRDS(state_A, "Auto_topmp_v2_states_A.rds")
saveRDS(state_B, "Auto_topmp_v2_states_B.rds")
saveRDS(mp_adj_noncc, "Auto_topmp_v2_mp_adj.rds")

summary_dir <- file.path("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "Auto_topmp_v2_summary.csv")

summary_df <- data.frame(
  approach = c(rep("A", length(table(state_A))), rep("B", length(table(state_B)))),
  state = c(names(table(state_A)), names(table(state_B))),
  cells = c(as.integer(table(state_A)), as.integer(table(state_B))),
  stringsAsFactors = FALSE
)
summary_df$pct <- ave(summary_df$cells, summary_df$approach, FUN = function(x) 100 * x / sum(x))
write.csv(summary_df, summary_path, row.names = FALSE)

cat("Saved: Auto_topmp_v2_states_A.rds\n")
cat("Saved: Auto_topmp_v2_states_B.rds\n")
cat("Saved: Auto_topmp_v2_mp_adj.rds\n")
cat("Saved: Auto_topmp_v2_mean_heatmap.pdf\n")
cat(sprintf("Saved: %s\n", summary_path))

cat("=== Auto states top-MP v2 complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")


###################################Hybrid subclassification############################

####################
# Auto_states_topmp_v2_hybridB.R
# Hybrid-only subdivision and visualisation based on Approach B from
# analysis/cell_states/Auto_states_topmp_hybrid.R
#
# Input:  PDOs_outs/PDOs_merged.rds
#         PDOs_outs/MP_outs_default.rds
#         PDOs_outs/UCell_scores_filtered.rds
#         PDOs_outs/Auto_topmp_v2_states_B.rds
#         PDOs_outs/Auto_topmp_v2_mp_adj.rds
#         /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
# Output: PDOs_outs/Auto_topmp_v2_hybridB_subtypes.rds
#         PDOs_outs/Auto_topmp_v2_hybridB_states_expanded.rds
#         PDOs_outs/Auto_topmp_v2_hybridB_heatmap.pdf
#         PDOs_outs/Auto_topmp_v2_hybridB_mean_heatmap.pdf
#         PDOs_outs/Auto_topmp_v2_hybridB_proportion.pdf
#         PDOs_outs/Auto_topmp_v2_hybridB_umap_top12.pdf
#         PDOs_outs/Auto_topmp_v2_realstates_plushybrid_umap_top12.pdf
#         PDOs_outs/Auto_topmp_v2_hybridB_batch_umap_top12.pdf
#         PDOs_outs/Auto_topmp_v2_hybridB_summary.csv
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

cat("=== Auto topMP v2 hybrid-B subtyping ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load inputs
####################
tmdata_all <- readRDS("PDOs_merged.rds")
geneNMF.metaprograms <- readRDS("MP_outs_default.rds") # UPDATE: Change to optimal nMP result after running Auto_find_optimal_nmf.R
ucell_scores <- readRDS("UCell_scores_filtered.rds") # UPDATE: Verify UCell output filename after geneNMF pipeline
state_B <- readRDS("Auto_topmp_v2_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_mp_adj.rds")

cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]
####################
# NOTE: Cancer signature scoring removed — not applicable to PDO-only analysis
####################

####################
# 2) Canonical MP setup and ordering
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
retained_mps <- names(mp.genes)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

####################
# UPDATE AFTER geneNMF: Identify which MPs correspond to cell cycle
# Look for MPs enriched in G2M, S-phase, G1/S genes
####################
cc_mps <- c("MP6", "MP7", "MP1", "MP3")

####################
# UPDATE AFTER geneNMF: Define state groupings based on MP biology
# Example structure (adapt to actual PDO MPs):
####################
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "Basal to Intest. Meta" = c("MP4")
)
# TEMPORARY: if empty, create one group per non-CC MP
if (length(state_groups) == 0) {
  non_cc <- setdiff(retained_mps, cc_mps)
  state_groups <- setNames(as.list(non_cc), paste0("Group_", seq_along(non_cc)))
}

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))

pair_levels <- combn(ordered_group_names, 2, simplify = FALSE)
pair_labels <- vapply(pair_levels, function(x) paste(x, collapse = "__"), character(1))

####################
# UPDATE AFTER geneNMF: Fill in MP descriptions based on enrichment results
# Run Auto_enrichment_annotation.R first, then annotate each MP
####################
mp_descriptions <- c(
  "MP6"  = "G2M Cell Cycle",
  "MP7"  = "DNA repair",
  "MP5"  = "MYC-related Proliferation",
  "MP1"  = "G2M checkpoint",
  "MP3"  = "G1S Cell Cycle",
  "MP8"  = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9"  = "ECM Remodeling Epi.",
  "MP4"  = "Intestinal Metaplasia"
)
# TEMPORARY: auto-generate names if not yet annotated
if (length(mp_descriptions) == 0) {
  mp_descriptions <- setNames(retained_mps, retained_mps)
}

####################
# 3) Align matrices and compute helper scores
####################
common_cells <- intersect(intersect(Cells(tmdata_all), rownames(ucell_scores)), names(state_B))
tmdata_all <- tmdata_all[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
state_B <- state_B[common_cells]

cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))

z_normalise <- function(mat, sample_var, batch_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.batch <- batch_var[rownames(mat)]

  batch_sd <- clust_df %>%
    group_by(.batch) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".batch") %>%
    as.matrix()

  batch_sd[is.na(batch_sd) | batch_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell

  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / batch_sd[clust_centered$.batch, mp]
  }

  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

sample_var <- tmdata_all$orig.ident
batch_var <- tmdata_all$Batch
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, batch_var)

mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc[common_cells, , drop = FALSE])

cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))


####################
# 4) Hybrid-B subdivision
####################
HYBRID_GAP_B <- 0.3

group_max <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_noncc))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc[common_cells, mps_avail]))
  apply(mp_adj_noncc[common_cells, mps_avail, drop = FALSE], 1, max)
})
group_max <- as.matrix(group_max)
rownames(group_max) <- common_cells

hybrid_cells <- names(state_B)[state_B == "Hybrid"]
cat(sprintf("Hybrid cells found in Approach B: %d\n", length(hybrid_cells)))
if (length(hybrid_cells) == 0) stop("No Hybrid cells found in state_B")

assign_hybrid_subtype <- function(score_vec, order_vec, gap_thr = 0.3) {
  score_vec <- score_vec[order_vec]
  top1 <- max(score_vec, na.rm = TRUE)
  active <- names(score_vec)[score_vec >= (top1 - gap_thr)]

  if (length(active) > 2) return("MultiHybrid_3plus")

  if (length(active) == 2) {
    active <- order_vec[order_vec %in% active]
    return(paste(active, collapse = "__"))
  }

  ord <- names(sort(score_vec, decreasing = TRUE))[1:2]
  ord <- order_vec[order_vec %in% ord]
  paste(ord, collapse = "__")
}

hybrid_subtype <- vapply(
  hybrid_cells,
  function(cl) assign_hybrid_subtype(group_max[cl, ], ordered_group_names, HYBRID_GAP_B),
  character(1)
)

hybrid_levels <- c(pair_labels, "MultiHybrid_3plus")
hybrid_subtype <- factor(hybrid_subtype, levels = hybrid_levels)
names(hybrid_subtype) <- hybrid_cells

expanded_state <- as.character(state_B)
names(expanded_state) <- names(state_B)
expanded_state[names(hybrid_subtype)] <- as.character(hybrid_subtype)
expanded_state <- factor(expanded_state)

saveRDS(hybrid_subtype, "Auto_topmp_v2_hybridB_subtypes.rds")
saveRDS(expanded_state, "Auto_topmp_v2_hybridB_states_expanded.rds")

####################
# 4b) 5 real states + Hybrid view (for requested UMAP display)
####################
state5_plus_hybrid <- as.character(state_B)
names(state5_plus_hybrid) <- names(state_B)
state5_plus_hybrid[state5_plus_hybrid %in% c("Unresolved")] <- NA_character_
state5_plus_hybrid <- factor(
  state5_plus_hybrid,
  levels = c(ordered_group_names, "Hybrid")
)

####################
# 5) Hybrid-only per-cell heatmap
####################
set.seed(42)
MAX_HYBRID_CELLS <- 6000
subtype_cells <- split(names(hybrid_subtype), hybrid_subtype)
subtype_counts <- table(hybrid_subtype)
subtype_fracs <- subtype_counts / sum(subtype_counts)
cells_per_subtype <- pmax(round(subtype_fracs * MAX_HYBRID_CELLS), 30)

cells_to_plot <- unlist(
  mapply(
    function(cells, n) sample(cells, min(length(cells), n)),
    subtype_cells,
    cells_per_subtype[names(subtype_cells)],
    SIMPLIFY = FALSE
  ),
  use.names = FALSE
)

sub_scores_orig <- t(mp_adj_all[cells_to_plot, , drop = FALSE])
cc_block_order <- cc_mps[cc_mps %in% rownames(sub_scores_orig)]
non_cc_block_order <- mp_tree_order_names[
  mp_tree_order_names %in% rownames(sub_scores_orig) & !(mp_tree_order_names %in% cc_mps)
]
mp_row_order <- c(cc_block_order, non_cc_block_order)
sub_scores <- sub_scores_orig[mp_row_order, , drop = FALSE]

mp_label_map <- mp_descriptions
missing_mps <- setdiff(rownames(sub_scores), names(mp_label_map))
if (length(missing_mps) > 0) mp_label_map[missing_mps] <- missing_mps
rownames(sub_scores) <- mp_label_map[rownames(sub_scores)]

split_vec <- factor(as.character(hybrid_subtype[cells_to_plot]), levels = hybrid_levels)
split_vec <- droplevels(split_vec)

meta <- tmdata_all@meta.data
batch_vals <- as.character(meta[cells_to_plot, "Batch"])

state_df <- data.frame(
  cell = names(hybrid_subtype),
  subtype = as.character(hybrid_subtype),
  sample = as.character(meta[names(hybrid_subtype), "orig.ident"]),
  batch = as.character(meta[names(hybrid_subtype), "Batch"]),
  stringsAsFactors = FALSE
)

total_samples <- length(unique(state_df$sample))
total_batches <- length(unique(state_df$batch))
subtype_div_df <- state_df %>%
  group_by(subtype) %>%
  summarise(
    sample_cov = n_distinct(sample) / max(total_samples, 1),
    batch_cov = n_distinct(batch) / max(total_batches, 1),
    diversity_score = 0.5 * sample_cov + 0.5 * batch_cov,
    .groups = "drop"
  )
subtype_div_map <- setNames(subtype_div_df$diversity_score, subtype_div_df$subtype)
div_vals <- subtype_div_map[as.character(split_vec)]
div_vals[is.na(div_vals)] <- 0
names(div_vals) <- cells_to_plot

hybrid_cols <- setNames(
  c(hue_pal()(length(pair_labels)), "black"),
  c(pair_labels, "MultiHybrid_3plus")
)
local_hybrid_cols <- hybrid_cols[levels(split_vec)]

batch_cols <- setNames(
  DiscretePalette(length(unique(meta$Batch)), palette = "polychrome"),
  unique(meta$Batch)
)

col_ann <- HeatmapAnnotation(
  HybridSubtype = split_vec,
  CC_score = cc_score[cells_to_plot],
  Diversity = div_vals,
  Batch = batch_vals,
  col = list(
    HybridSubtype = local_hybrid_cols,
    CC_score = colorRamp2(c(0, max(cc_score[cells_to_plot], na.rm = TRUE)), c("white", "darkgreen")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Batch = batch_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE,
  na_col = "white"
)

mp_to_group <- rep("Other", length(mp_row_order))
names(mp_to_group) <- mp_row_order
mp_to_group[cc_mps[cc_mps %in% names(mp_to_group)]] <- "Cell_cycle"
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], names(mp_to_group))
  mp_to_group[grp_mps] <- grp
}

group_cols <- setNames(
  c(hue_pal()(length(ordered_group_names)), "gold", "grey70"),
  c(ordered_group_names, "Cell_cycle", "Other")
)

mp_group_label <- mp_to_group
names(mp_group_label) <- rownames(sub_scores)

row_ann <- rowAnnotation(
  MP_group = factor(mp_group_label, levels = c("Cell_cycle", ordered_group_names, "Other")),
  col = list(MP_group = group_cols),
  show_annotation_name = FALSE
)

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

ht_sc <- Heatmap(
  sub_scores,
  name = "Adj score",
  col = col_fun_sc,
  top_annotation = col_ann,
  left_annotation = row_ann,
  column_split = split_vec,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  cluster_rows = FALSE,
  row_split = factor(ifelse(mp_row_order %in% cc_mps, "Cell_cycle_MPs", "Other_MPs"),
                     levels = c("Cell_cycle_MPs", "Other_MPs")),
  row_gap = unit(2.5, "mm"),
  column_gap = unit(1.5, "mm"),
  show_row_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8.5, fontface = "italic"),
  show_column_names = FALSE,
  column_title_rot = 30,
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = gpar(col = NA)
)

pdf("Auto_topmp_v2_hybridB_heatmap.pdf", width = 19, height = max(7, length(mp_row_order) * 0.5), useDingbats = FALSE)
draw(ht_sc, merge_legend = TRUE)
dev.off()

####################
# 6) Hybrid-only proportion plot
####################
overall_df <- data.frame(subtype = hybrid_subtype) %>%
  count(subtype) %>%
  mutate(pct = 100 * n / sum(n), batch = "Overall")

per_batch_df <- data.frame(
  subtype = hybrid_subtype,
  batch = meta[names(hybrid_subtype), "Batch"],
  stringsAsFactors = FALSE
) %>%
  count(subtype, batch) %>%
  group_by(batch) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

plot_df <- bind_rows(overall_df, per_batch_df)
plot_df$subtype <- factor(plot_df$subtype, levels = hybrid_levels)
plot_df$batch <- factor(plot_df$batch, levels = c("Overall", sort(unique(as.character(meta$Batch)))))

p_prop <- ggplot(plot_df, aes(x = batch, y = pct, fill = subtype)) +
  geom_col(colour = "black", linewidth = 0.15) +
  scale_fill_manual(values = hybrid_cols) +
  labs(
    title = "Hybrid-B subtype proportions",
    subtitle = "10 double-hybrid classes + MultiHybrid_3plus",
    x = NULL,
    y = "% of hybrid cells",
    fill = "Hybrid subtype"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Auto_topmp_v2_hybridB_proportion.pdf", p_prop, width = 14, height = 7)

####################
# 7) Mean heatmap per hybrid subtype (non-CC MPs)
####################
mean_scores_df <- as.data.frame(mp_adj_noncc[names(hybrid_subtype), , drop = FALSE])
mean_scores_df$subtype <- hybrid_subtype
mean_mat <- mean_scores_df %>%
  group_by(subtype) %>%
  summarise(across(all_of(colnames(mp_adj_noncc)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames("subtype") %>%
  as.matrix()

ordered_subtypes_present <- hybrid_levels[hybrid_levels %in% rownames(mean_mat)]
mean_mat <- mean_mat[ordered_subtypes_present, , drop = FALSE]

non_cc_tree_order <- mp_tree_order_names[
  mp_tree_order_names %in% colnames(mp_adj_noncc) & !(mp_tree_order_names %in% cc_mps)
]
mean_mat <- mean_mat[, non_cc_tree_order, drop = FALSE]
colnames(mean_mat) <- mp_descriptions[colnames(mean_mat)]

cell_counts <- table(hybrid_subtype)
row_anno <- rowAnnotation(
  Cells = anno_barplot(as.numeric(cell_counts[rownames(mean_mat)]), gp = gpar(fill = "steelblue"), width = unit(2, "cm"))
)

col_fun_z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mean <- Heatmap(
  mean_mat,
  name = "Mean Z-score",
  col = col_fun_z,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9.5, fontface = "italic"),
  column_names_rot = 45,
  column_title = "Mean non-CC MP activity per hybrid-B subtype",
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  right_annotation = row_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mean_mat[i, j]), x, y, gp = gpar(fontsize = 6.5))
  }
)

pdf("Auto_topmp_v2_hybridB_mean_heatmap.pdf", width = 15, height = 8, useDingbats = FALSE)
draw(ht_mean, merge_legend = TRUE)
dev.off()

####################
# 8) UMAPs in per-sample style (separate real-state and hybrid views)
####################
if ("umap" %in% names(tmdata_all@reductions)) {
  # Real states + Hybrid (single Hybrid class), per-sample facets
  state_real_hybrid <- as.character(state_B)
  names(state_real_hybrid) <- names(state_B)
  state_real_hybrid[state_real_hybrid == "Unresolved"] <- NA_character_
  state_real_hybrid <- factor(state_real_hybrid, levels = c(ordered_group_names, "Hybrid"))
  names(state_real_hybrid) <- names(state_B)
  tmdata_all$Auto_state_real_hybrid <- state_real_hybrid[Cells(tmdata_all)]

  state5_cols <- setNames(
    c(hue_pal()(length(ordered_group_names)), "orchid"),
    c(ordered_group_names, "Hybrid")
  )

  # Rank samples by abundance of defined states (excluding unresolved)
  md_rank <- tmdata_all@meta.data
  rank_df <- md_rank %>%
    dplyr::mutate(state_plot = Auto_state_real_hybrid) %>%
    dplyr::filter(!is.na(state_plot)) %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::summarise(total_defined = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_defined))

  # safe slicing: compute n_top outside of dplyr mask
  n_top <- min(12, nrow(rank_df))
  if (n_top == 0) n_top <- 0
  rank_df <- rank_df %>% dplyr::slice_head(n = n_top)

  top_ids <- as.character(rank_df$orig.ident)

  # Build one-page-per-sample list for real-state+hybrid view
  real_plots <- vector("list", length(top_ids))
  for (i in seq_along(top_ids)) {
    sid <- top_ids[i]
    sub_obj <- subset(tmdata_all, subset = orig.ident == sid)
    # Recompute a fast per-sample UMAP (Normalize -> HVG -> Scale -> PCA -> UMAP)
    sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sub_obj <- ScaleData(sub_obj, verbose = FALSE)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
    sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

    p <- DimPlot(
      sub_obj,
      reduction = "umap",
      group.by = "Auto_state_real_hybrid",
      cols = state5_cols,
      pt.size = 0.7
    ) +
      labs(
        title = sid,
        subtitle = sprintf("Defined cells: %d", rank_df$total_defined[i])
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8)
      )
    real_plots[[i]] <- p
  }

  pdf("Auto_topmp_v2_realstates_plushybrid_umap_top12.pdf", width = 12, height = 9, useDingbats = FALSE)
  for (p in real_plots) print(p)
  dev.off()

  # Batch-colored UMAP view per sample
  batch_levels <- unique(as.character(tmdata_all$Batch))
  batch_cols_umap <- setNames(
    DiscretePalette(length(batch_levels), palette = "polychrome"),
    batch_levels
  )
  batch_plots <- vector("list", length(top_ids))
  for (i in seq_along(top_ids)) {
    sid <- top_ids[i]
    sub_obj <- subset(tmdata_all, subset = orig.ident == sid)
    sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sub_obj <- ScaleData(sub_obj, verbose = FALSE)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
    sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

    p <- DimPlot(
      sub_obj,
      reduction = "umap",
      group.by = "Batch",
      cols = batch_cols_umap,
      pt.size = 0.7
    ) +
      labs(
        title = sid,
        subtitle = "Batch-colored UMAP",
        colour = "Batch"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8)
      )
    batch_plots[[i]] <- p
  }

  pdf("Auto_topmp_v2_hybridB_batch_umap_top12.pdf", width = 12, height = 9, useDingbats = FALSE)
  for (p in batch_plots) print(p)
  dev.off()

  # Hybrid-only subtype view, per-sample facets
  hybrid_plot_label <- rep(NA_character_, length(common_cells))
  names(hybrid_plot_label) <- common_cells
  hybrid_plot_label[names(hybrid_subtype)] <- as.character(hybrid_subtype)
  hybrid_plot_label <- factor(hybrid_plot_label, levels = hybrid_levels)
  tmdata_all$Auto_hybrid_subtype <- hybrid_plot_label[Cells(tmdata_all)]

  # Build one-page-per-sample list for hybrid-only subtype view
  hybrid_plots <- vector("list", length(top_ids))
  for (i in seq_along(top_ids)) {
    sid <- top_ids[i]
    sub_obj <- subset(tmdata_all, subset = orig.ident == sid)
    # Recompute UMAP per sample to ensure locality and comparable layout
    sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sub_obj <- ScaleData(sub_obj, verbose = FALSE)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), verbose = FALSE)
    sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)

    p <- DimPlot(
      sub_obj,
      reduction = "umap",
      group.by = "Auto_hybrid_subtype",
      cols = hybrid_cols,
      pt.size = 0.7
    ) +
      labs(
        title = sid,
        subtitle = sprintf("Hybrid-labelled cells in sample"),
        colour = "Hybrid subtype"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8)
      )
    hybrid_plots[[i]] <- p
  }

  pdf("Auto_topmp_v2_hybridB_umap_top12.pdf", width = 12, height = 9, useDingbats = FALSE)
  for (p in hybrid_plots) print(p)
  dev.off()
}

####################
# 9) Summary output
####################
summary_df <- data.frame(
  subtype = names(table(hybrid_subtype)),
  cells = as.integer(table(hybrid_subtype)),
  stringsAsFactors = FALSE
)
summary_df$pct <- 100 * summary_df$cells / sum(summary_df$cells)

summary_dir <- file.path("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "Auto_topmp_v2_hybridB_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

cat("Saved: Auto_topmp_v2_hybridB_subtypes.rds\n")
cat("Saved: Auto_topmp_v2_hybridB_states_expanded.rds\n")
cat("Saved: Auto_topmp_v2_hybridB_heatmap.pdf\n")
cat("Saved: Auto_topmp_v2_hybridB_mean_heatmap.pdf\n")
cat("Saved: Auto_topmp_v2_hybridB_proportion.pdf\n")
cat("Saved: Auto_topmp_v2_realstates_plushybrid_umap_top12.pdf\n")
cat("Saved: Auto_topmp_v2_hybridB_umap_top12.pdf\n")
cat("Saved: Auto_topmp_v2_hybridB_batch_umap_top12.pdf\n")
cat(sprintf("Saved: %s\n", summary_path))
cat("=== Auto topMP v2 hybrid-B subtyping complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
