####################
# Auto_PDO_mp_correlation_crossdata.R
#
# Compare PDO MPs (nMP=13) against scRef MPs
# Generates: gene overlap heatmap (Jaccard), dataset-level 3CA correlation
#
# Input:
#   PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#   PDOs_outs/UCell_3CA_MPs.rds
#   scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   scRef_Pipeline/ref_outs/UCell_3CA_MPs.rds
#
# Output:
#   PDOs_outs/Auto_PDO_mp_correlation_jaccard.pdf
#   PDOs_outs/Auto_PDO_mp_correlation_crossdata_scatter.pdf
#   PDOs_outs/Auto_PDO_mp_correlation_crossdata_bar.pdf
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# Load PDO data
####################
message("Loading PDO metaprograms...")
MP_pdo <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")

message("Loading PDO 3CA UCell scores...")
ucell_3ca_pdo <- readRDS("UCell_3CA_MPs.rds")

####################
# Load scRef data
####################
message("Loading scRef metaprograms...")
MP_sc <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

message("Loading scRef 3CA UCell scores...")
ucell_3ca_sc <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_3CA_MPs.rds")

####################
# MP descriptions - STRICT
####################
pdo_mp_descriptions <- c(
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

sc_mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

pdo_mp_descriptions <- setNames(
  paste("PDOs", names(pdo_mp_descriptions), pdo_mp_descriptions, sep = "_"),
  names(pdo_mp_descriptions)
)

sc_mp_descriptions <- setNames(
  paste("scATLAS", names(sc_mp_descriptions), sc_mp_descriptions, sep = "_"),
  names(sc_mp_descriptions)
)

####################
# Filter and get tree order for PDO
# PDO: silhouette >= 0 AND sample coverage >= 0.25
####################
pdo_list <- MP_pdo$metaprograms.genes

pdo_sil <- MP_pdo$metaprograms.metrics$silhouette
names(pdo_sil) <- paste0("MP", seq_along(pdo_sil))
pdo_cov <- MP_pdo$metaprograms.metrics$sampleCoverage
names(pdo_cov) <- paste0("MP", seq_along(pdo_cov))

bad_pdo <- names(pdo_sil)[pdo_sil < 0]
low_cov_pdo <- names(pdo_cov)[pdo_cov < 0.25]
drop_mps_pdo <- unique(c(bad_pdo, low_cov_pdo))

# Get tree order from PDO
pdo_tree_order <- MP_pdo$programs.tree$order
pdo_ordered_clusters <- MP_pdo$programs.clusters[pdo_tree_order]
pdo_mp_tree_order <- rev(unique(pdo_ordered_clusters))
pdo_mp_tree_order <- paste0("MP", pdo_mp_tree_order)

# FILTER: only keep MPs that have descriptions AND pass filters
pdo_mp_tree_order <- pdo_mp_tree_order[pdo_mp_tree_order %in% names(pdo_mp_descriptions)]
pdo_mp_tree_order <- pdo_mp_tree_order[!pdo_mp_tree_order %in% drop_mps_pdo]

message("PDO MPs (sil>=0 & cov>=0.25, with desc): ", paste(pdo_mp_tree_order, collapse = ", "))

####################
# Filter and get tree order for scRef
# scRef: silhouette >= 0 only
####################
sc_list <- MP_sc$metaprograms.genes

sc_sil <- MP_sc$metaprograms.metrics$silhouette
names(sc_sil) <- paste0("MP", seq_along(sc_sil))

bad_sc <- names(sc_sil)[sc_sil < 0]
drop_mps_sc <- unique(bad_sc)

# Get tree order from scRef
sc_tree_order <- MP_sc$programs.tree$order
sc_ordered_clusters <- MP_sc$programs.clusters[sc_tree_order]
sc_mp_tree_order <- rev(unique(sc_ordered_clusters))
sc_mp_tree_order <- paste0("MP", sc_mp_tree_order)

# FILTER: only keep MPs that have descriptions AND pass silhouette filter
sc_mp_tree_order <- sc_mp_tree_order[sc_mp_tree_order %in% names(sc_mp_descriptions)]
sc_mp_tree_order <- rev(sc_mp_tree_order[!sc_mp_tree_order %in% drop_mps_sc])

message("scRef MPs (sil>=0 only, with desc): ", paste(sc_mp_tree_order, collapse = ", "))

####################
# JACCARD: Filter gene lists and rename
####################
pdo_list <- pdo_list[pdo_mp_tree_order]
sc_list <- sc_list[sc_mp_tree_order]

# Use descriptions directly
names(pdo_list) <- pdo_mp_descriptions[names(pdo_list)]
names(sc_list) <- sc_mp_descriptions[names(sc_list)]

####################
# Gene overlap (Jaccard)
####################
message("Computing gene overlap...")
universe <- unique(c(unlist(pdo_list), unlist(sc_list)))
n_pdo <- length(pdo_list)
n_sc <- length(sc_list)

jaccard_mat <- matrix(NA_real_, n_pdo, n_sc, dimnames = list(names(pdo_list), names(sc_list)))
overlap_n_mat <- jaccard_mat
pval_mat <- jaccard_mat

for (i in seq_len(n_pdo)) {
  A <- pdo_list[[i]]
  for (j in seq_len(n_sc)) {
    B <- sc_list[[j]]
    inter <- length(intersect(A, B))
    uni <- length(union(A, B))
    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j] <- if (uni == 0) NA_real_ else inter / uni
    
    a <- inter
    b <- length(setdiff(A, B))
    c <- length(setdiff(B, A))
    d <- length(setdiff(universe, union(A, B)))
    pval_mat[i, j] <- if (any(c(a,b,c,d) < 0)) NA_real_
    else fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
  }
}

padj_mat <- matrix(p.adjust(as.vector(pval_mat), method = "BH"), nrow = n_pdo, ncol = n_sc, dimnames = dimnames(pval_mat))

stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat), dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05] <- "*"
stars_mat[padj_mat < 0.01] <- "**"
stars_mat[padj_mat < 0.001] <- "***"

display_mat <- paste0(overlap_n_mat, "\n", stars_mat)
dim(display_mat) <- dim(overlap_n_mat)
dimnames(display_mat) <- dimnames(overlap_n_mat)

pdf("Auto_PDO_mp_correlation_jaccard.pdf", width = 12, height = 10, useDingbats = FALSE)
pheatmap(t(jaccard_mat),
         cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey85",
         main = "Gene sets overlap (Jaccard)",
         angle_col = "90",
         display_numbers = t(display_mat),
         fontsize_number = 8,
         number_color = "black",
         fontsize_row = 10,
         fontsize_col = 10,
         color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100))
dev.off()

####################
# 3CA correlation (dataset-level)
####################
message("Computing 3CA UCell score correlation...")
pdo_cols <- colnames(ucell_3ca_pdo)
sc_cols <- colnames(ucell_3ca_sc)
common_cols <- intersect(pdo_cols, sc_cols)

pdo_mean_scores <- colMeans(ucell_3ca_pdo[, common_cols, drop = FALSE], na.rm = TRUE)
sc_mean_scores <- colMeans(ucell_3ca_sc[, common_cols, drop = FALSE], na.rm = TRUE)

common_mps <- intersect(names(pdo_mean_scores), names(sc_mean_scores))
comp_df <- data.frame(MP = common_mps, PDO_score = pdo_mean_scores[common_mps], scRef_score = sc_mean_scores[common_mps])

cor_val <- cor(comp_df$PDO_score, comp_df$scRef_score, method = "spearman")

# Add status and labeling threshold
comp_df$Label <- ifelse(comp_df$PDO_score >= 0.1 | comp_df$scRef_score >= 0.1, sub("^X3CA_mp_", "", comp_df$MP), NA)
comp_df$Status <- ifelse(comp_df$PDO_score < 0.1 & comp_df$scRef_score < 0.1, "Low", "Significant")

# Determine max limit for synced axes
max_limit <- max(c(comp_df$scRef_score, comp_df$PDO_score), na.rm = TRUE) * 1.05

p_scatter <- ggplot(comp_df, aes(x = scRef_score, y = PDO_score)) +
  # Threshold lines
  geom_vline(xintercept = 0.1, linetype = "dotted", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "black", linewidth = 0.4, alpha = 0.5) +
  # Points
  geom_point(aes(color = Status), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Low" = "grey60", "Significant" = "black")) +
  # Repel labels - Reverted to simpler original style
  geom_text_repel(aes(label = Label), size = 2.5, max.overlaps = 20, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", fill = "red", alpha = 0.1) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  # Sync axes
  xlim(0, max_limit) + 
  ylim(0, max_limit) +
  coord_fixed() +
  labs(title = "PDO vs scAtlas 3CA MP mean scores",
       subtitle = "Threshold: score >= 0.1 in at least one dataset",
       x = "scAtlas mean Score",
       y = "PDO mean Score") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
ggsave("Auto_PDO_mp_correlation_crossdata_scatter.pdf", p_scatter, width = 9, height = 9, useDingbats = FALSE)

# Sample ID extraction function
extract_sample_id <- function(x) {
  ifelse(
    grepl("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", x),
    sub("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", "", x),
    x
  )
}

# Compute Sample-level stats for each MP
message("Computing sample-level distributions for 3CA MPs...")
pdo_samples <- extract_sample_id(rownames(ucell_3ca_pdo))
sc_samples <- extract_sample_id(rownames(ucell_3ca_sc))

# For each MP, calculate mean score per sample
pdo_sample_dist <- list()
sc_sample_dist <- list()
p_vals_dist <- c()

for (m in common_mps) {
  # PDO sample means
  p_df <- data.frame(score = ucell_3ca_pdo[, m], sample = pdo_samples)
  p_means <- p_df %>% group_by(sample) %>% summarize(mean_score = mean(score, na.rm=TRUE)) %>% pull(mean_score)
  pdo_sample_dist[[m]] <- p_means
  
  # scAtlas sample means
  s_df <- data.frame(score = ucell_3ca_sc[, m], sample = sc_samples)
  s_means <- s_df %>% group_by(sample) %>% summarize(mean_score = mean(score, na.rm=TRUE)) %>% pull(mean_score)
  sc_sample_dist[[m]] <- s_means
  
  # Wilcoxon test between sample means
  p_vals_dist[m] <- if(length(p_means) >= 3 && length(s_means) >= 3) wilcox.test(p_means, s_means)$p.value else NA
}

adj_p_dist <- p.adjust(p_vals_dist, method = "BH")

# Create plot data for boxplots (sample-level means)
all_sample_means <- list()
for (m in common_mps) {
  df_m <- data.frame(
    MP = m,
    Score = c(pdo_sample_dist[[m]], sc_sample_dist[[m]]),
    Dataset = c(rep("PDO", length(pdo_sample_dist[[m]])), rep("scAtlas", length(sc_sample_dist[[m]])))
  )
  all_sample_means[[m]] <- df_m
}
plot_df_dist <- do.call(rbind, all_sample_means)

library(ggplot2)

# 1. Prepare significance labels dataframe
sig_df <- data.frame(
  MP = common_mps,
  adj_p = adj_p_dist,
  PDO_mean = sapply(pdo_sample_dist, mean),
  scRef_mean = sapply(sc_sample_dist, mean)
)

# Calculate stars
sig_df$stars <- cut(sig_df$adj_p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), labels=c("***", "**", "*", ""))
sig_df$stars <- as.character(sig_df$stars)
sig_df$stars[is.na(sig_df$stars)] <- ""

# Clean MP names
sig_df$MP_label <- sub("^X3CA_mp_", "", sig_df$MP)

# 🔥 FORMATTING FIX: Append stars in brackets directly to the MP name (if significant)
sig_df$MP_annot <- paste0(sig_df$MP_label, ifelse(sig_df$stars == "", "", paste0(" (", sig_df$stars, ")")))

# Sort by PDO_mean to establish factor levels for plotting
sig_df <- sig_df[order(sig_df$PDO_mean, decreasing = TRUE), ]
sig_df$MP_annot <- factor(sig_df$MP_annot, levels = rev(sig_df$MP_annot))

# 2. Merge annotated labels into your main plot dataframe
plot_df_dist$MP_label <- sub("^X3CA_mp_", "", plot_df_dist$MP)
plot_df_dist$MP_annot <- sig_df$MP_annot[match(plot_df_dist$MP_label, sig_df$MP_label)]

# 3. Build the Plot
p_box <- ggplot(plot_df_dist, aes(x = MP_annot, y = Score, fill = Dataset)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    position = position_dodge(0.75),     # Tighter dodging for a compact look
    color = "black",
    linewidth = 0.4,                     # Replaced 'size' with modern 'linewidth'
    alpha = 0.85,
    coef = 0                             # Kept as requested to remove whiskers
  ) +
  # Note: stat_summary (white dot) and geom_text (floating stars) are completely removed
  scale_fill_manual(values = c(PDO = "#E41A1C", scAtlas = "#377EB8")) +
  coord_flip() +
  # 🔥 AESTHETICS FIX: Cleaner, modern, and more beautiful theme
  theme_classic(base_size = 12) +
  labs(
    title = "PDO vs scAtlas 3CA MP Distributions",
    subtitle = "Significance: *** p<0.001, ** p<0.01, * p<0.05",
    x = NULL,                            # Removed redundant y-axis title
    y = "Mean Score per Sample"
  ) +
  theme(
    legend.position = "top",             # Moves legend to top to save lateral space
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "lines"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"), # Soft vertical guides
    axis.line.y = element_blank(),       # Removes harsh vertical axis line
    axis.ticks.y = element_blank()       # Removes cluttered y-axis ticks
  )

# 4. Save
# Note: 12x15 inches is massive and causes awkward spacing. 
# Reduced dimensions to 8x10 to force the plot to be physically compact.
ggsave("Auto_PDO_mp_correlation_crossdata_bar.pdf", p_box, width = 8, height = 10, useDingbats = FALSE)
####################
# Cross-correlation: PDO MPs vs scRef MPs in scRef cells (SAMPLE-AVERAGED)
####################
message("=== DONE ===")
