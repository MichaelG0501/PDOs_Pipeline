####################
# Analysis registry:
#   Status: active
#   Script: analysis/cell_states/Auto_PDO_unresolved_pan_cancer_noreg_comparison.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     PDOs_outs/PDOs_merged.rds
#     PDOs_outs/Auto_PDO_states_noreg.rds
#     /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
#     /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#   Outputs:
#     PDOs_outs/UCell_3CA_MPs.rds
#     PDOs_outs/Auto_PDO_noreg_unresolved_subclass.csv
#     PDOs_outs/Auto_PDO_noreg_unresolved_subclass.rds
#     PDOs_outs/Auto_PDO_noreg_unresolved_mp_adj.rds
#     PDOs_outs/Auto_PDO_noreg_unresolved_heatmap.pdf
#     PDOs_outs/Auto_PDO_noreg_unresolved_barplot.pdf
####################

library(Seurat)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

args <- commandArgs(trailingOnly = TRUE)
ncores <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 8L

message("Loading PDOs_merged.rds...")
pdos <- readRDS("PDOs_merged.rds")
state_noreg <- readRDS("Auto_PDO_states_noreg.rds")

# Fix batch for standardisation
pdos$Batch <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")

message("Loading 3CA MP signatures...")
mp_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
if (!file.exists(mp_path)) stop("Missing 3CA MP gene list: ", mp_path)

mp_df <- read.csv(mp_path, check.names = FALSE)
mp_list <- as.list(mp_df)
mp_list <- lapply(mp_list, function(x) unique(x[x != "" & !is.na(x)]))
mp_list <- mp_list[lengths(mp_list) > 0]
names(mp_list) <- make.names(sub("^MP", "3CA_mp_", names(mp_list)))

message("Running AddModuleScore_UCell for 3CA MPs...")
pdos <- AddModuleScore_UCell(pdos, features = mp_list, ncores = ncores, name = "")
score_cols <- grep("^X3CA_mp|^3CA_mp", colnames(pdos@meta.data), value = TRUE)
if (length(score_cols) == 0) stop("No 3CA UCell score columns were generated.")

ucell_3ca <- pdos@meta.data[, score_cols, drop = FALSE]
saveRDS(ucell_3ca, file = "UCell_3CA_MPs.rds")

message("Saved UCell_3CA_MPs.rds with ", nrow(ucell_3ca), " cells and ", ncol(ucell_3ca), " MPs.")

cc_genes_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"
cell_cycle_genes <- read.csv(cc_genes_path, header = TRUE, stringsAsFactors = FALSE)[, 1:3]

common_cells <- intersect(Cells(pdos), rownames(ucell_3ca))
pdos <- pdos[, common_cells]
ucell_3ca <- ucell_3ca[common_cells, , drop = FALSE]
state_noreg <- state_noreg[common_cells]

cc_consensus <- intersect(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1], rownames(pdos))
cc_top50 <- names(sort(rowMeans(pdos@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(pdos@assays$RNA$data[cc_top50, , drop = FALSE]))

# Pre-defined colors from scRef pan-cancer analysis to ensure exact matching
scref_colors <- c(
  "X3CA_mp_30.Respiration.1" = "#D078FF", "X3CA_mp_17.EMT.III" = "#A9A400", 
  "X3CA_mp_12.Protein.maturation" = "#F8766D", "X3CA_mp_34.Secreted.II" = "#00C094", 
  "X3CA_mp_43.PDAC.classical" = "#C49A00", "X3CA_mp_6.Stress.1" = "#EB8335", 
  "X3CA_mp_4.Chromatin" = "#EC69EF", "X3CA_mp_19.EMT.V" = "#FB61D7", 
  "X3CA_mp_9.Stress.2" = "#00BA38", "X3CA_mp_3.Cell.Cylce.HMG.rich" = "#53B400", 
  "X3CA_mp_33.Secreted.I" = "#FF6B96", "X3CA_mp_66.Unassigned.1" = "#00BE6D", 
  "X3CA_mp_63.PDAC.related.5" = "#FF63B9", "X3CA_mp_5.Cell.cycle.single.nucleus" = "#86AC00", 
  "X3CA_mp_22.Interferon.MHC.II..I." = "#00C0B5", "X3CA_mp_20.EMT.VI" = "#A58AFF", 
  "X3CA_mp_27.MYC" = "#00B6EB", "X3CA_mp_10.Proteasomal.degradation" = "#DA8F00", 
  "X3CA_mp_59.PDAC.related.1" = "#00BDD2", "X3CA_mp_60.PDAC.related.2" = "#00ABFD", 
  "X3CA_mp_1.Cell.Cycle...G2.M" = "#619CFF"
)

THRESHOLD <- 0.10

classify_mode <- function(state_vec, mode_name) {
  unresolved_cells <- names(state_vec)[state_vec == "Unresolved"]
  unresolved_cells <- intersect(unresolved_cells, rownames(ucell_3ca))
  if (length(unresolved_cells) == 0) return(NULL)

  sub_scores <- ucell_3ca[unresolved_cells, , drop = FALSE]
  
  sample_var   <- pdos$orig.ident[unresolved_cells]
  study_var    <- pdos$Batch[unresolved_cells]
  
  clust_df <- as.data.frame(sub_scores)
  clust_df$.cell   <- rownames(sub_scores)
  clust_df$.sample <- sample_var
  clust_df$.study  <- study_var

  mps <- colnames(sub_scores)

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(mps), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()

  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(mps), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, mps])
  rownames(mp_adj) <- clust_centered$.cell

  for (mp in mps) {
    cell_studies <- clust_centered$.study
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[cell_studies, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0

  # Use raw UCell scores for classification (no normalisation) as requested
  sub_scores_raw <- ucell_3ca[unresolved_cells, mps, drop = FALSE]
  top_mp_idx <- max.col(sub_scores_raw, ties.method = "first")
  subclass <- mps[top_mp_idx]

  out <- data.frame(
    cell = unresolved_cells,
    mode = mode_name,
    unresolved_subclass = subclass,
    stringsAsFactors = FALSE
  )

  write.csv(out, paste0("Auto_PDO_", mode_name, "_unresolved_subclass.csv"), row.names = FALSE)
  saveRDS(out, paste0("Auto_PDO_", mode_name, "_unresolved_subclass.rds"))
  saveRDS(mp_adj, paste0("Auto_PDO_", mode_name, "_unresolved_mp_adj.rds"))
  list(df = out, mp_adj = mp_adj)
}

make_unresolved_heatmap <- function(res_list, mode_name) {
  if (is.null(res_list) || is.null(res_list$df)) return(NULL)

  df_mode <- res_list$df
  mp_adj <- res_list$mp_adj
  cells <- df_mode$cell
  scores <- mp_adj[cells, , drop = FALSE]
  split_vec <- factor(df_mode$unresolved_subclass)

  set.seed(42)
  max_cells <- 6000
  subtype_cells <- split(cells, split_vec)
  subtype_counts <- table(split_vec)
  subtype_fracs <- subtype_counts / sum(subtype_counts)
  cells_per_subtype <- pmax(round(subtype_fracs * max_cells), 20)
  cells_to_plot <- unlist(mapply(
    function(cset, n) sample(cset, min(length(cset), n)),
    subtype_cells,
    cells_per_subtype[names(subtype_cells)],
    SIMPLIFY = FALSE
  ), use.names = FALSE)

  sub_scores <- t(scores[cells_to_plot, , drop = FALSE])
  rownames(sub_scores) <- gsub("\\.", " ", gsub("^X3CA_", "3CA_", rownames(sub_scores)))
  split_plot <- factor(df_mode$unresolved_subclass[match(cells_to_plot, df_mode$cell)])
  split_plot <- droplevels(split_plot)

  subtype_cols <- scref_colors[levels(split_plot)]
  missing_ht <- is.na(subtype_cols)
  if(any(missing_ht)) subtype_cols[missing_ht] <- scales::hue_pal()(sum(missing_ht))
  names(subtype_cols) <- levels(split_plot)

  study_vals <- pdos@meta.data[cells_to_plot, "Batch"]
  study_cols <- setNames(
    DiscretePalette(length(unique(pdos$Batch)), palette = "polychrome"),
    unique(pdos$Batch)
  )

  col_ann <- HeatmapAnnotation(
    Subclass = split_plot,
    CC_score = cc_score[cells_to_plot],
    Study = study_vals,
    col = list(
      Subclass = subtype_cols,
      CC_score = colorRamp2(c(0, max(cc_score[cells_to_plot], na.rm = TRUE)), c("white", "darkgreen")),
      Study = study_cols
    ),
    annotation_name_side = "left",
    na_col = "white"
  )

  lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
  Heatmap(
    sub_scores,
    name = "Adj score",
    col = colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3")),
    top_annotation = col_ann,
    column_split = split_plot,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    show_row_dend = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    use_raster = TRUE,
    raster_quality = 5,
    column_title = paste0("Unresolved subclassification (", mode_name, ")")
  )
}

res_noreg <- classify_mode(state_noreg, "noreg")

pdf("Auto_PDO_noreg_unresolved_heatmap.pdf", width = 18, height = 10, useDingbats = FALSE)
if (!is.null(res_noreg)) {
  ht_noreg <- make_unresolved_heatmap(res_noreg, "noreg")
  if (!is.null(ht_noreg)) draw(ht_noreg, merge_legend = TRUE)
}
dev.off()

summary_df <- res_noreg$df %>%
  count(mode, unresolved_subclass, name = "cells") %>%
  group_by(mode) %>%
  mutate(pct = 100 * cells / sum(cells)) %>%
  ungroup() %>%
  arrange(desc(pct))

write.csv(summary_df, "Auto_PDO_noreg_unresolved_summary.csv", row.names = FALSE)

# Subset for the plot to show only the top half of MPs
summary_df_plot <- head(summary_df, n = ceiling(nrow(summary_df) / 2))
summary_df_plot$unresolved_subclass <- factor(summary_df_plot$unresolved_subclass, levels = summary_df_plot$unresolved_subclass)

# Clean subclass labels for plot
summary_df_plot$unresolved_subclass_label <- gsub("\\.", " ", gsub("^X3CA_mp_|^3CA_mp_", "", as.character(summary_df_plot$unresolved_subclass)))
summary_df_plot$unresolved_subclass_label <- factor(summary_df_plot$unresolved_subclass_label, levels = unique(summary_df_plot$unresolved_subclass_label))

# Map colors to the cleaned labels
my_colors <- setNames(scref_colors[as.character(summary_df_plot$unresolved_subclass)], as.character(summary_df_plot$unresolved_subclass_label))

# Fallback for classes not in the top 21 of scRef
missing_classes <- is.na(my_colors)
if (any(missing_classes)) {
  fallback_colors <- scales::hue_pal()(sum(missing_classes))
  my_colors[missing_classes] <- fallback_colors
}

p_bar <- ggplot(summary_df_plot, aes(x = unresolved_subclass_label, y = pct, fill = unresolved_subclass_label)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 5.3, fontface = "bold") +
  theme_classic(base_size = 20) +
  labs(
    title = "Unresolved pan-cancer subclass proportion",
    subtitle = "PDO noreg",
    x = "",
    y = "% unresolved cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 15)),
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 20)),
    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = my_colors)

pdf("Auto_PDO_noreg_unresolved_barplot.pdf", width = 12, height = 12)
print(p_bar)
dev.off()

message("Saved unified unresolved outputs (noreg only).")
