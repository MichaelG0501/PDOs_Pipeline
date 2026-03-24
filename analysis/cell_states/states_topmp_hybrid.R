####################
# Auto_PDO_hybrid_subtyping_noreg.R
# Hybrid-only subdivision and visualisation
# Dependent on outputs from Auto_PDO_states_analysis.R (noreg mode)
#
# Input:  PDOs_outs/PDOs_final.rds
#         PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#         PDOs_outs/UCell_scores_filtered.rds
#         PDOs_outs/Auto_PDO_states_noreg.rds
#         PDOs_outs/Auto_PDO_mp_adj_noreg.rds
#         /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
# Output: PDOs_outs/Auto_PDO_hybrid_subtypes_noreg.rds
#         PDOs_outs/Auto_PDO_hybrid_states_expanded_noreg.rds
#         PDOs_outs/Auto_PDO_hybrid_heatmap_noreg.pdf
#         PDOs_outs/Auto_PDO_hybrid_proportion_noreg.pdf
#         PDOs_outs/Auto_PDO_hybrid_mean_heatmap_noreg.pdf
#         PDOs_outs/Auto_PDO_hybrid_summary_noreg.csv
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

cat("=== Auto PDO hybrid subtyping (noreg) ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load inputs
####################
tmdata_all <- readRDS("PDOs_final.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds") 
ucell_scores <- readRDS("UCell_scores_filtered.rds") 
state_B <- readRDS("Auto_PDO_states_noreg.rds")
mp_adj_noncc <- readRDS("Auto_PDO_mp_adj_noreg.rds")

cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]

# Fix batch: only two batches to match main noreg script
tmdata_all$Batch <- ifelse(tmdata_all$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")

####################
# 2) Canonical MP setup and ordering (Two-Step Filter)
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes

# Step 1: Silhouette filtering
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)

# Step 2: Sample coverage filtering
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]

# Apply both filters
mp.genes <- mp.genes[!names(mp.genes) %in% c(bad_mp_names, low_coverage_mps)]
retained_mps <- names(mp.genes)

# Get MP tree order (reversed for PDO)
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- rev(unique(ordered_clusters))
mp_tree_order_names <- paste0("MP", mp_tree_order)
mp_tree_order_names <- mp_tree_order_names[mp_tree_order_names %in% retained_mps]

# CC vs Non-CC MPs
cc_mps <- c("MP6", "MP7", "MP1", "MP3")
cc_mps <- cc_mps[cc_mps %in% retained_mps]
non_cc_mps <- setdiff(retained_mps, cc_mps)

# State groups (scRef nomenclature)
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "SMG-like Metaplasia"   = c("MP8")
)

# Filter state groups to only include available MPs
state_groups <- lapply(state_groups, function(mps) mps[mps %in% retained_mps])
state_groups <- state_groups[sapply(state_groups, length) > 0]

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))

pair_levels <- combn(ordered_group_names, 2, simplify = FALSE)
pair_labels <- vapply(pair_levels, function(x) paste(x, collapse = "__"), character(1))

# MP descriptions for plotting
mp_descriptions <- c(
  "MP6"  = "MP6_G2M Cell Cycle",
  "MP7"  = "MP7_DNA repair",
  "MP5"  = "MP5_MYC-related Proliferation",
  "MP1"  = "MP1_G2M checkpoint",
  "MP3"  = "MP3_G1S Cell Cycle",
  "MP8"  = "MP8_Columnar Progenitor",
  "MP10" = "MP10_Inflammatory Stress Epi.",
  "MP9"  = "MP9_ECM Remodeling Epi.",
  "MP4"  = "MP4_Intestinal Metaplasia"
)

####################
# 3) Align matrices and compute helper scores
####################
common_cells <- intersect(intersect(Cells(tmdata_all), rownames(ucell_scores)), names(state_B))
tmdata_all <- tmdata_all[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]

cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))

# Need CC scores to attach to mp_adj_all for heatmap
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

mp_adj_all <- cbind(mp_adj_noncc, mp_adj_cc[common_cells, , drop = FALSE])

cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(tmdata_all))
cc_top50 <- names(sort(rowMeans(tmdata_all@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(tmdata_all@assays$RNA$data[cc_top50, , drop = FALSE]))


####################
# 4) Hybrid subdivision (noreg logic)
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
cat(sprintf("Hybrid cells found: %d\n", length(hybrid_cells)))
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

saveRDS(hybrid_subtype, "Auto_PDO_hybrid_subtypes_noreg.rds")
saveRDS(expanded_state, "Auto_PDO_hybrid_states_expanded_noreg.rds")

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

batch_cols_cmap <- setNames(c("brown", "darkgreen"), c("Cynthia_batch", "New_batch"))

col_ann <- HeatmapAnnotation(
  HybridSubtype = split_vec,
  CC_score = cc_score[cells_to_plot],
  Diversity = div_vals,
  Batch = batch_vals,
  col = list(
    HybridSubtype = local_hybrid_cols,
    CC_score = colorRamp2(c(0, max(cc_score[cells_to_plot], na.rm = TRUE)), c("white", "darkgreen")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Batch = batch_cols_cmap
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

group_cols_ann <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Cell_cycle"            = "gold",
  "Other"                 = "grey70"
)

mp_group_label <- mp_to_group
names(mp_group_label) <- rownames(sub_scores)

row_ann <- rowAnnotation(
  MP_group = factor(mp_group_label, levels = c("Cell_cycle", ordered_group_names, "Other")),
  col = list(MP_group = group_cols_ann),
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

pdf("Auto_PDO_hybrid_heatmap_noreg.pdf", width = 19, height = max(7, length(mp_row_order) * 0.5), useDingbats = FALSE)
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
    title = "Hybrid subtype proportions (noreg)",
    x = NULL,
    y = "% of hybrid cells",
    fill = "Hybrid subtype"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Auto_PDO_hybrid_proportion_noreg.pdf", p_prop, width = 14, height = 7)

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
  column_title = "Mean non-CC MP activity per hybrid subtype (noreg)",
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  right_annotation = row_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mean_mat[i, j]), x, y, gp = gpar(fontsize = 6.5))
  }
)

pdf("Auto_PDO_hybrid_mean_heatmap_noreg.pdf", width = 15, height = 8, useDingbats = FALSE)
draw(ht_mean, merge_legend = TRUE)
dev.off()

####################
# 8) Summary output
####################
summary_df <- data.frame(
  subtype = names(table(hybrid_subtype)),
  cells = as.integer(table(hybrid_subtype)),
  stringsAsFactors = FALSE
)
summary_df$pct <- 100 * summary_df$cells / sum(summary_df$cells)

summary_path <- file.path(getwd(), "Auto_PDO_hybrid_summary_noreg.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

cat("Saved outputs to PDOs_outs/ :\n")
cat(" - Auto_PDO_hybrid_subtypes_noreg.rds\n")
cat(" - Auto_PDO_hybrid_states_expanded_noreg.rds\n")
cat(" - Auto_PDO_hybrid_heatmap_noreg.pdf\n")
cat(" - Auto_PDO_hybrid_proportion_noreg.pdf\n")
cat(" - Auto_PDO_hybrid_mean_heatmap_noreg.pdf\n")
cat(" - Auto_PDO_hybrid_summary_noreg.csv\n")
cat("=== Auto PDO hybrid subtyping complete ===\n")