####################
# Analysis registry (authoritative override):
#   Status: legacy; retained for provenance, no current downstream use
#   Script: analysis/metaprograms/legacy_PDO_mp_correlation_crossdata.R
#   Methodology: historical method only; see analysis/ANALYSIS_MAP.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Preserves a superseded implementation or analysis tied to superseded
#     inputs. Do not use its outputs as current centred-MP/state inputs. The
#     original historical inputs, outputs, and method notes remain below.
####################

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
pdo_mp_descriptions_raw <- c(
  "MP11"  = "Single-nucleus-associated cell cycle",
  "MP1"   = "G2/M cell cycle",
  "MP2"   = "G1/S cell cycle",
  "MP3"   = "Replication-dependent histones",
  "MP19+" = "MYC-associated proliferation",
  "MP15"  = "Intestinal metaplasia",
  "MP5+"  = "Inflammatory-reactive columnar epithelium",
  "MP12"  = "KRAS-active columnar epithelium",
  "MP13b" = "Metabolic-detox columnar epithelium",
  "MP14b" = "Proliferative epithelial plasticity",
  "MP16b" = "EMT/KRAS adaptive plasticity",
  "MP17+" = "Ciliated progenitor epithelium",
  "MP8+"  = "Secretory-transport glandular epithelium",
  "MP9"   = "ECM-remodelling epithelium",
  "MP18"  = "Motile-cilia differentiation"
)

sc_mp_descriptions_raw <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "Single-nucleus cell cycle",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Inflammatory-reactive columnar epithelium",
  "MP11+" = "Epithelial type I interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP18b" = "Mucous-secretory differentiation",
  "MP17" = "Immune-interactive glandular progenitor",
  "MP12" = "Hypoxic inflammatory adaptive plasticity",
  "MP15" = "T/NK-like cancer-cell immune mimicry"
)

state_groups_pdo <- list(
  "Cell cycle" = c("MP11", "MP1", "MP2", "MP3"),
  "Classic proliferation" = c("MP19+"),
  "Columnar-to-intestinal" = c("MP14b", "MP13b", "MP5+", "MP12", "MP15"),
  "Glandular differentiation" = c("MP17+", "MP8+"),
  "Stress-adaptive" = c("MP16b"),
  "ECM-remodelling" = c("MP9"),
  "Motile-cilia differentiation" = c("MP18")
)

state_cols_pdo <- c(
  "Cell cycle" = "#6B7280",
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Motile-cilia differentiation" = "#F781BF"
)

state_groups_sc <- list(
  "Cell cycle" = c("MP1", "MP5", "MP13+"),
  "Classic proliferation" = c("MP2+"),
  "Squamous-to-intestinal" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "Glandular-to-intestinal" = c("MP8+", "MP8b", "MP16", "MP18b", "MP17"),
  "Stress-adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)

state_cols_sc <- c(
  "Cell cycle" = "#6B7280",
  "Classic proliferation" = "#E41A1C",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)

pdo_mp_descriptions <- setNames(
  paste(names(pdo_mp_descriptions_raw), pdo_mp_descriptions_raw, sep = ": "),
  names(pdo_mp_descriptions_raw)
)

sc_mp_descriptions <- setNames(
  paste(names(sc_mp_descriptions_raw), sc_mp_descriptions_raw, sep = ": "),
  names(sc_mp_descriptions_raw)
)

####################
# Filter and get tree order for PDO
# Use strict state-based order (matching Auto_05)
####################
# Load refined PDO MPs and apply threshold filter (coverage >= 3/24 samples, ngenes >= 10)
pdo_merged_mp_genes <- readRDS("centred_mp_refinement/merged_refined_mp_genes.rds")
pdo_list <- pdo_merged_mp_genes

metrics_path <- "centred_mp_refinement/tables/merged_refined_mp_metrics.rds"
if (file.exists(metrics_path)) {
  metrics_df <- readRDS(metrics_path)
  cov_threshold <- 3 / 24 - 1e-5
  min_genes <- 10
  retained <- rownames(metrics_df)[!is.na(metrics_df$sampleCoverage) & metrics_df$sampleCoverage >= cov_threshold &
                                    !is.na(metrics_df$numberGenes) & metrics_df$numberGenes >= min_genes]
  pdo_list <- pdo_list[names(pdo_list) %in% retained]
}

strict_pdo_order <- unlist(state_groups_pdo, use.names = FALSE)
pdo_mp_tree_order <- strict_pdo_order[strict_pdo_order %in% names(pdo_list)]

message("PDO Refined MPs (strict order): ", paste(pdo_mp_tree_order, collapse = ", "))

####################
# Filter and get tree order for scRef
# scRef: use predefined ordered groups (already filtered)
####################
sc_list <- merged_mp_genes

strict_sc_order <- unlist(state_groups_sc, use.names = FALSE)
sc_mp_tree_order <- strict_sc_order[strict_sc_order %in% names(sc_list)]
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

clean_3ca_label <- function(x) {
  x <- gsub("^X3CA_mp_", "3CA_", x)
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("^XMP", "3CA_MP", x)
  x <- gsub("^MP", "3CA_MP", x)
  x <- gsub("^3CA_3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x <- gsub(" +", " ", x)
  x
}

pdo_cols_clean <- clean_3ca_label(pdo_cols)
sc_cols_clean <- clean_3ca_label(sc_cols)

pdo_base <- sub("^3CA_", "", pdo_cols_clean)
sc_base <- sub("^3CA_", "", sc_cols_clean)
common_bases <- intersect(pdo_base, sc_base)

if (length(common_bases) > 0) {
  pdo_match_cols <- pdo_cols[match(common_bases, pdo_base)]
  sc_match_cols <- sc_cols[match(common_bases, sc_base)]
  
  pdo_mean_scores <- colMeans(ucell_3ca_pdo[, pdo_match_cols, drop = FALSE], na.rm = TRUE)
  sc_mean_scores <- colMeans(ucell_3ca_sc[, sc_match_cols, drop = FALSE], na.rm = TRUE)
  
  # Align names to bases for dataframe
  names(pdo_mean_scores) <- common_bases
  names(sc_mean_scores) <- common_bases
  
  comp_df <- data.frame(MP = common_bases, PDO_score = pdo_mean_scores, scRef_score = sc_mean_scores)


cor_val <- cor(comp_df$PDO_score, comp_df$scRef_score, method = "spearman")

# Add status and labeling threshold
comp_df$Label <- ifelse(comp_df$PDO_score >= 0.1 | comp_df$scRef_score >= 0.1, sub("^mp_", "", comp_df$MP), NA)
comp_df$Status <- ifelse(comp_df$PDO_score < 0.1 & comp_df$scRef_score < 0.1, "Low", "Significant")

# Determine max limit for synced axes safely
if (nrow(comp_df) > 0 && sum(!is.na(comp_df$scRef_score)) > 0) {
  max_limit <- max(c(comp_df$scRef_score, comp_df$PDO_score), na.rm = TRUE) * 1.05
} else {
  max_limit <- 1.0
}

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

for (i in seq_along(common_bases)) {
  base_m <- common_bases[i]
  p_col <- pdo_match_cols[i]
  s_col <- sc_match_cols[i]
  
  # PDO sample means
  p_df <- data.frame(score = ucell_3ca_pdo[, p_col], sample = pdo_samples)
  p_means <- p_df %>% group_by(sample) %>% summarize(mean_score = mean(score, na.rm=TRUE)) %>% pull(mean_score)
  pdo_sample_dist[[base_m]] <- p_means
  
  # scAtlas sample means
  s_df <- data.frame(score = ucell_3ca_sc[, s_col], sample = sc_samples)
  s_means <- s_df %>% group_by(sample) %>% summarize(mean_score = mean(score, na.rm=TRUE)) %>% pull(mean_score)
  sc_sample_dist[[base_m]] <- s_means
  
  # Wilcoxon test between sample means
  p_vals_dist[base_m] <- if(length(p_means) >= 3 && length(s_means) >= 3) wilcox.test(p_means, s_means)$p.value else NA
}

adj_p_dist <- p.adjust(p_vals_dist, method = "BH")

# Create plot data for boxplots (sample-level means)
all_sample_means <- list()
for (m in common_bases) {
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
  MP = common_bases,
  adj_p = adj_p_dist,
  PDO_mean = sapply(pdo_sample_dist, mean),
  scRef_mean = sapply(sc_sample_dist, mean)
)

# Calculate stars
sig_df$stars <- cut(sig_df$adj_p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), labels=c("***", "**", "*", ""))
sig_df$stars <- as.character(sig_df$stars)
sig_df$stars[is.na(sig_df$stars)] <- ""

# Clean MP names
sig_df$MP_label <- sub("^mp_", "", sig_df$MP)

# 🔥 FORMATTING FIX: Append stars in brackets directly to the MP name (if significant)
sig_df$MP_annot <- paste0(sig_df$MP_label, ifelse(sig_df$stars == "", "", paste0(" (", sig_df$stars, ")")))

# Sort by PDO_mean to establish factor levels for plotting
sig_df <- sig_df[order(sig_df$PDO_mean, decreasing = TRUE), ]
sig_df$MP_annot <- factor(sig_df$MP_annot, levels = rev(sig_df$MP_annot))

# 2. Merge annotated labels into your main plot dataframe
plot_df_dist$MP_label <- sub("^mp_", "", plot_df_dist$MP)
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
  message("No common 3CA MPs found between datasets. Skipping scatter/bar plots.")
}
} else {
  message("Skipping 3CA UCell correlation because input files are missing.")
}
####################
# Cross-correlation: PDO MPs vs scRef MPs in PDO cells (SAMPLE-AVERAGED)
####################
sc_ucell_path <- "Auto_PDO_ucell_scATLAS_MPs_in_PDOs.rds"

sc_mps_to_score <- sc_list_for_ucell
names(sc_mps_to_score) <- paste0("scATLAS_", names(sc_mps_to_score))

if (file.exists(sc_ucell_path)) {
  message("=== Loading saved scATLAS UCell scores ===")
  sc_ucell_scores <- readRDS(sc_ucell_path) # MP x cells
  sc_mat <- t(sc_ucell_scores) # cells x MP
} else {
  message("=== Scoring scATLAS MPs in PDOs ===")
  library(UCell)
  pdo_obj <- if (file.exists("PDOs_final.rds")) readRDS("PDOs_final.rds") else readRDS("PDOs_merged.rds")
  
  existing_cols <- intersect(colnames(pdo_obj@meta.data), names(sc_mps_to_score))
  if (length(existing_cols) > 0) {
    pdo_obj@meta.data <- pdo_obj@meta.data[, !colnames(pdo_obj@meta.data) %in% existing_cols, drop = FALSE]
  }
  pdo_obj <- AddModuleScore_UCell(pdo_obj, features = sc_mps_to_score, ncores = 1, name = "")
  
  sc_ucell_scores <- t(as.matrix(pdo_obj@meta.data[, names(sc_mps_to_score), drop = FALSE]))
  saveRDS(sc_ucell_scores, sc_ucell_path)
  message("Saved UCell scores for scATLAS MPs in PDOs to ", sc_ucell_path)
  sc_mat <- t(sc_ucell_scores)
}

pdo_ucell_path <- "centred_mp_refinement/merged_refined_ucell_scores.rds"
pdo_mps_to_score <- pdo_list_for_ucell
names(pdo_mps_to_score) <- paste0("PDO_", names(pdo_mps_to_score))

if (file.exists(pdo_ucell_path)) {
  message("=== Loading saved PDO UCell scores ===")
  pdo_ucell_scores <- readRDS(pdo_ucell_path) # cells x MP
  # match original MP names
  intersect_cols <- intersect(colnames(pdo_ucell_scores), names(pdo_list_for_ucell))
  pdo_mat <- pdo_ucell_scores[, intersect_cols, drop = FALSE]
  # prefix with PDO_ to match downstream code
  colnames(pdo_mat) <- paste0("PDO_", colnames(pdo_mat))
} else {
  message("=== Scoring PDO MPs in PDOs ===")
  if (!exists("pdo_obj")) pdo_obj <- if (file.exists("PDOs_final.rds")) readRDS("PDOs_final.rds") else readRDS("PDOs_merged.rds")
  
  existing_pdo_cols <- intersect(colnames(pdo_obj@meta.data), names(pdo_mps_to_score))
  if (length(existing_pdo_cols) > 0) {
    pdo_obj@meta.data <- pdo_obj@meta.data[, !colnames(pdo_obj@meta.data) %in% existing_pdo_cols, drop = FALSE]
  }
  pdo_obj <- AddModuleScore_UCell(pdo_obj, features = pdo_mps_to_score, ncores = 1, name = "")
  pdo_mat <- as.matrix(pdo_obj@meta.data[, names(pdo_mps_to_score), drop = FALSE])
}

message("=== Computing sample-averaged cross-correlation ===")
sc_cols <- names(sc_mps_to_score)
pdo_cols <- names(pdo_mps_to_score)

# ensure we only take the columns we need
pdo_mat <- pdo_mat[, pdo_cols, drop = FALSE]
sc_mat <- sc_mat[, sc_cols, drop = FALSE]

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

mp_to_state_pdo <- unlist(lapply(names(state_groups_pdo), function(state) {
  setNames(rep(state, length(state_groups_pdo[[state]])), state_groups_pdo[[state]])
}))
mp_to_state_sc <- unlist(lapply(names(state_groups_sc), function(state) {
  setNames(rep(state, length(state_groups_sc[[state]])), state_groups_sc[[state]])
}))

combined_state_order <- c(
  "Cell cycle",
  "Classic proliferation",
  "Columnar-to-intestinal",
  "Squamous-to-intestinal",
  "Glandular differentiation",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "ECM-remodelling",
  "Cancer-cell immune mimicry",
  "Motile-cilia differentiation"
)

state_vec_for_mps_pdo <- factor(mp_to_state_pdo[pdo_mp_tree_order], levels = combined_state_order)
state_vec_for_mps_sc <- factor(mp_to_state_sc[sc_mp_tree_order], levels = combined_state_order)

combined_colors <- c(
  "Cell cycle" = "#6B7280",
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Motile-cilia differentiation" = "#F781BF"
)

ha_left <- rowAnnotation(
  State = state_vec_for_mps_sc,
  col = list(State = combined_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(title = "State", at = names(combined_colors))
)
ha_top <- HeatmapAnnotation(
  State = state_vec_for_mps_pdo,
  col = list(State = combined_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))
hm_width <- unit(10.5, "inch")
hm_height <- unit(10.5, "inch")

ht_cor <- Heatmap(
  mean_rho,
  name = "Mean Spearman\n(Meta-analysis)",
  col = col_cor,
  rect_gp = gpar(col = "white", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = ha_left,
  top_annotation = ha_top,
  row_split = state_vec_for_mps_sc,
  column_split = state_vec_for_mps_pdo,
  column_title_rot = 20,
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  row_title = NULL,
  row_names_side = "right",
  column_names_side = "bottom",
  column_names_rot = 30,
  row_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
  row_names_max_width = unit(128, "mm"),
  column_names_max_height = unit(128, "mm"),
  width = hm_width,
  height = hm_height,
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_vals[i, j]
    rho <- mean_rho[i, j]
    if (is.na(p) || is.na(rho)) {
      grid.text("NA", x, y, gp = gpar(fontsize = 8.5, col = "grey50"))
    } else if (p < 0.001) {
      grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 8.5))
    } else if (p < 0.01) {
      grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 8.5))
    } else if (p < 0.05) {
      grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 8.5))
    } else {
      grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 8.5))
    }
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
)

pdf("Auto_PDO_mp_correlation_crossdata_expression_heatmap.pdf", width = 22, height = 18, useDingbats = FALSE)
draw(
  ht_cor,
  heatmap_legend_side = "left",
  padding = unit(c(20, 20, 20, 20), "mm")
)
dev.off()

message("=== DONE ===")
