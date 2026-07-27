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

setwd("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs")

####################
# Load PDO data
####################
message("Loading PDO metaprograms...")
MP_pdo <- readRDS("centred_mp_refinement/geneNMF_metaprograms_nMP_13.rds")

message("Loading 3CA MP definitions from CSV...")
library(UCell)
csv_path <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
mp_df <- read.csv(csv_path, stringsAsFactors = FALSE)
mp_3ca_list <- as.list(mp_df)
mp_3ca_list <- lapply(mp_3ca_list, function(x) x[x != "" & !is.na(x)])
# Prefix to ensure valid column names
names(mp_3ca_list) <- paste0("X3CA_", names(mp_3ca_list))

message("Loading PDO 3CA UCell scores...")
pdo_3ca_path <- "UCell_3CA_MPs.rds"
if (file.exists(pdo_3ca_path)) {
  ucell_3ca_pdo <- readRDS(pdo_3ca_path)
} else {
  message("PDO 3CA scores not found. Computing on the fly...")
  pdo_obj_3ca <- if (file.exists("PDOs_final.rds")) readRDS("PDOs_final.rds") else readRDS("PDOs_merged.rds")
  pdo_obj_3ca <- AddModuleScore_UCell(pdo_obj_3ca, features = mp_3ca_list, ncores = 1, name = "")
  ucell_3ca_pdo <- pdo_obj_3ca@meta.data[, names(mp_3ca_list), drop = FALSE]
  saveRDS(ucell_3ca_pdo, pdo_3ca_path)
  rm(pdo_obj_3ca)
  gc()
}

####################
# Load scRef data
####################
message("Loading scRef metaprograms...")
sc_mp_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/centred_mp_refinement/merged_refined_mp_genes.rds"
if(file.exists(sc_mp_path)){
  merged_mp_genes <- readRDS(sc_mp_path)
} else {
  # Fallback to older path just in case
  merged_mp_genes <- readRDS("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds")
}

message("Loading scRef 3CA UCell scores...")
sc_3ca_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/UCell_3CA_MPs.rds"
if (file.exists(sc_3ca_path)) {
  ucell_3ca_sc <- readRDS(sc_3ca_path)
} else {
  message("scRef 3CA scores not found. Computing on the fly...")
  sc_obj_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds"
  if (file.exists(sc_obj_path)) {
    sc_obj_3ca <- readRDS(sc_obj_path)
    sc_obj_3ca <- AddModuleScore_UCell(sc_obj_3ca, features = mp_3ca_list, ncores = 1, name = "")
    ucell_3ca_sc <- sc_obj_3ca@meta.data[, names(mp_3ca_list), drop = FALSE]
    saveRDS(ucell_3ca_sc, sc_3ca_path)
    rm(sc_obj_3ca)
    gc()
  } else {
    message("Could not find scRef Seurat object at ", sc_obj_path)
    ucell_3ca_sc <- NULL
  }
}

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

state_groups_sc <- list(
  "Cell cycle" = c("MP1", "MP5", "MP13+"),
  "Classic proliferation" = c("MP2+"),
  "Basal to intestinal metaplasia" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "SMG to intestinal metaplasia" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17", "MP2x"),
  "Stress adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15"),
  "Excluded" = c("MP11c", "MP18a")
)

sc_mp_descriptions_raw <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "replication-stress-associated cell cycling",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Stress-reactive columnar epithelium",
  "MP11+" = "Epithelial antiviral interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP2x" = "Wnt-active glandular stem/progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry",
  "MP11c" = "Excluded",
  "MP18a" = "Excluded"
)

pdo_mp_descriptions <- setNames(
  paste("PDOs", names(pdo_mp_descriptions), pdo_mp_descriptions, sep = "_"),
  names(pdo_mp_descriptions)
)

sc_mp_descriptions <- setNames(
  paste("scATLAS", names(sc_mp_descriptions_raw), sc_mp_descriptions_raw, sep = "_"),
  names(sc_mp_descriptions_raw)
)

####################
# Filter and get tree order for PDO
# Use refined merged MPs and unsupervised clustering order
####################
####################
# Load refined PDO MPs and apply threshold filter (coverage >= 3/24 samples, ngenes >= 10)
pdo_merged_mp_genes <- readRDS("centred_mp_refinement/merged_refined_mp_genes.rds")
pdo_list <- pdo_merged_mp_genes

# Load correlation matrix to get exact unsupervised clustering order
cached_cor <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/centred_mp_refinement/intermediate/merged_refined_mp_correlation_matrices.rds"
cor_matrices <- readRDS(cached_cor)
mean_rho <- cor_matrices$mean_rho

metrics_path <- "centred_mp_refinement/tables/merged_refined_mp_metrics.rds"
if (!file.exists(metrics_path)) {
  metrics_path <- file.path(dirname(cached_cor), "merged_refined_mp_metrics.rds")
}
if (file.exists(metrics_path)) {
  metrics_df <- readRDS(metrics_path)
  cov_threshold <- 3 / 24 - 1e-5
  min_genes <- 10
  retained <- rownames(metrics_df)[!is.na(metrics_df$sampleCoverage) & metrics_df$sampleCoverage >= cov_threshold &
                                    !is.na(metrics_df$numberGenes) & metrics_df$numberGenes >= min_genes]
  pdo_list <- pdo_list[names(pdo_list) %in% retained]
  valid_cor_mps <- intersect(colnames(mean_rho), retained)
  mean_rho <- mean_rho[valid_cor_mps, valid_cor_mps, drop = FALSE]
}

# Extract exact unsupervised clustering order using ComplexHeatmap
ht_cor_unsup <- ComplexHeatmap::Heatmap(mean_rho, cluster_rows = TRUE, cluster_columns = TRUE)
pdf(NULL)
hm_drawn <- ComplexHeatmap::draw(ht_cor_unsup)
dev.off()
final_col_order <- colnames(mean_rho)[ComplexHeatmap::column_order(hm_drawn)]

# Ensure we only use available MPs
pdo_mp_tree_order <- final_col_order[final_col_order %in% names(pdo_list)]

# Note: Descriptions for PDO MPs are just their raw names currently
pdo_mp_descriptions <- setNames(pdo_mp_tree_order, pdo_mp_tree_order)

message("PDO Refined MPs (ordered): ", paste(pdo_mp_tree_order, collapse = ", "))
####################

####################
# Filter and get tree order for scRef
# scRef: use predefined ordered groups (already filtered)
####################
sc_list <- merged_mp_genes

strict_refined_mp_order <- unlist(state_groups_sc, use.names = FALSE)
sc_mp_tree_order <- strict_refined_mp_order[strict_refined_mp_order %in% names(sc_list)]
sc_mp_tree_order <- sc_mp_tree_order[sc_mp_tree_order %in% names(sc_mp_descriptions_raw)]

message("scRef MPs (ordered): ", paste(sc_mp_tree_order, collapse = ", "))

####################
# JACCARD: Filter gene lists and rename
####################
pdo_list <- pdo_list[pdo_mp_tree_order]
sc_list <- sc_list[sc_mp_tree_order]

pdo_list_for_ucell <- pdo_list
sc_list_for_ucell <- sc_list

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
# Helper functions
####################
# Sample ID extraction function
extract_sample_id <- function(x) {
  ifelse(
    grepl("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", x),
    sub("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", "", x),
    x
  )
}

####################
# 3CA correlation (dataset-level)
####################
if (!is.null(ucell_3ca_pdo) && !is.null(ucell_3ca_sc)) {
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
} else {
  message("Skipping 3CA UCell correlation because input files are missing.")
}
####################
# Cross-correlation: PDO MPs vs scRef MPs in PDO cells (SAMPLE-AVERAGED)
####################
message("=== Scoring scATLAS MPs in PDOs ===")
library(UCell)
pdo_obj <- if (file.exists("PDOs_final.rds")) readRDS("PDOs_final.rds") else readRDS("PDOs_merged.rds")

# AddModuleScore_UCell for scATLAS MPs
sc_mps_to_score <- sc_list_for_ucell
# Safely prefix names to prevent metadata column collision
names(sc_mps_to_score) <- paste0("scATLAS_", names(sc_mps_to_score))

existing_cols <- intersect(colnames(pdo_obj@meta.data), names(sc_mps_to_score))
if (length(existing_cols) > 0) {
  pdo_obj@meta.data <- pdo_obj@meta.data[, !colnames(pdo_obj@meta.data) %in% existing_cols, drop = FALSE]
}
pdo_obj <- AddModuleScore_UCell(pdo_obj, features = sc_mps_to_score, ncores = 1, name = "")

# AddModuleScore_UCell for PDO MPs
message("=== Scoring PDO MPs in PDOs ===")
pdo_mps_to_score <- pdo_list_for_ucell
names(pdo_mps_to_score) <- paste0("PDO_", names(pdo_mps_to_score))

existing_pdo_cols <- intersect(colnames(pdo_obj@meta.data), names(pdo_mps_to_score))
if (length(existing_pdo_cols) > 0) {
  pdo_obj@meta.data <- pdo_obj@meta.data[, !colnames(pdo_obj@meta.data) %in% existing_pdo_cols, drop = FALSE]
}
pdo_obj <- AddModuleScore_UCell(pdo_obj, features = pdo_mps_to_score, ncores = 1, name = "")

# Save UCell scores locally (in live as requested)
sc_ucell_scores <- t(as.matrix(pdo_obj@meta.data[, names(sc_mps_to_score), drop = FALSE]))
saveRDS(sc_ucell_scores, "Auto_PDO_ucell_scATLAS_MPs_in_PDOs.rds")
message("Saved UCell scores for scATLAS MPs in PDOs to Auto_PDO_ucell_scATLAS_MPs_in_PDOs.rds")

message("=== Computing sample-averaged cross-correlation ===")
sc_cols <- names(sc_mps_to_score)
pdo_cols <- names(pdo_mps_to_score)

pdo_mat <- as.matrix(pdo_obj@meta.data[, pdo_cols, drop = FALSE])
sc_mat <- as.matrix(pdo_obj@meta.data[, sc_cols, drop = FALSE])

mod_mat <- t(pdo_mat)
ref_mat <- t(sc_mat)

# Map back to descriptions (we need to strip the prefixes we added for UCell)
rownames(mod_mat) <- pdo_mp_descriptions[sub("^PDO_", "", names(pdo_mps_to_score))]
rownames(ref_mat) <- sc_mp_descriptions[sub("^scATLAS_", "", names(sc_mps_to_score))]

sample_ids <- extract_sample_id(colnames(mod_mat))
unique_samples <- unique(sample_ids)
unique_samples <- unique_samples[!is.na(unique_samples) & unique_samples != ""]

n_pdo <- nrow(mod_mat)
n_sc <- nrow(ref_mat)

cor_list <- list()
for (smp in unique_samples) {
  idx <- which(sample_ids == smp)
  if (length(idx) >= 10) {
    pdo_scores <- mod_mat[, idx, drop = FALSE]
    ref_scores <- ref_mat[, idx, drop = FALSE]
    
    cors <- matrix(NA, n_sc, n_pdo) # Rows = scRef, Cols = PDO
    for (i in seq_len(n_sc)) {
      for (j in seq_len(n_pdo)) {
        if (sd(ref_scores[i,]) > 0 && sd(pdo_scores[j,]) > 0) {
          cors[i, j] <- cor(ref_scores[i,], pdo_scores[j,], method = "spearman")
        }
      }
    }
    rownames(cors) <- rownames(ref_mat)
    colnames(cors) <- rownames(mod_mat)
    cor_list[[smp]] <- cors
  }
}

if (length(cor_list) == 0) {
  stop("No samples with >= 10 cells found for cross-correlation.")
}

sum_z <- matrix(0, n_sc, n_pdo)
count_z <- matrix(0, n_sc, n_pdo)
z_list <- list()

for (smp in names(cor_list)) {
  z <- atanh(pmin(pmax(cor_list[[smp]], -0.999), 0.999))
  z_list[[smp]] <- z
  valid <- !is.na(z)
  sum_z[valid] <- sum_z[valid] + z[valid]
  count_z[valid] <- count_z[valid] + 1
}

mean_z <- sum_z / count_z
mean_rho <- tanh(mean_z)

p_vals <- matrix(1, n_sc, n_pdo)
for (i in seq_len(n_sc)) {
  for (j in seq_len(n_pdo)) {
    vals <- sapply(z_list, function(z) z[i, j])
    vals <- vals[!is.na(vals)]
    if (length(vals) >= 3) {
      if (sd(vals, na.rm = TRUE) == 0) {
        p_vals[i, j] <- if (mean(vals, na.rm = TRUE) == 0) 1 else 1e-16
      } else {
        p_vals[i, j] <- t.test(vals, mu = 0)$p.value
      }
    }
  }
}

dimnames(mean_rho) <- list(rownames(ref_mat), rownames(mod_mat))
dimnames(p_vals) <- dimnames(mean_rho)

col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

pdf("Auto_PDO_mp_correlation_crossdata_expression_heatmap.pdf", width = 14, height = 10, useDingbats = FALSE)
ht <- Heatmap(mean_rho, name = "Mean Spearman\n(Meta-analysis)", col = col_cor,
  cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_vals[i, j]
    r_val <- mean_rho[i, j]
    lvl <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    grid.text(sprintf("%.2f\n%s", r_val, lvl), x, y, gp = gpar(fontsize = 9, fontface = "bold"))
  }, 
  row_names_gp = gpar(fontsize = 10), 
  column_names_gp = gpar(fontsize = 10),
  column_title = "PDO Metaprograms",
  row_title = "scATLAS Metaprograms")
draw(ht)
dev.off()

message("=== DONE ===")
