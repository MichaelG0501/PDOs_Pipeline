####################
# Auto_PDO_mp_correlation_crossdata.R
#
# Compare PDO MPs (nMP=13) against scRef MPs
# Generates: gene overlap heatmap (Jaccard), expression correlation
#
# Input:
#   PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#   PDOs_outs/UCell_3CA_MPs.rds
#   scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds
#   scRef_Pipeline/ref_outs/UCell_3CA_MPs.rds
#   scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds
#   scRef_Pipeline/ref_outs/UCell_ref_terms_v2.rds
#
# Output:
#   PDOs_outs/Auto_PDO_mp_correlation_jaccard.pdf
#   PDOs_outs/Auto_PDO_mp_correlation_crossdata_scatter.pdf
#   PDOs_outs/Auto_PDO_mp_correlation_crossdata_bar.pdf
#   PDOs_outs/Auto_PDO_mp_crossref_correlation_meta.pdf
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
  "MP1" = "G2M_checkpoint",
  "MP3" = "G1S_Cycle",
  "MP4" = "Intest_diff",
  "MP5" = "MYC_Biosynth",
  "MP6" = "G2M_mitotic",
  "MP7" = "DNA",
  "MP8" = "Columnar_progenitor",
  "MP9" = "EMT_related",
  "MP10" = "Stress_induced_plasticity"
)

sc_mp_descriptions <- c(
  "MP1"  = "G2M_cycle",
  "MP2"  = "MYC_prolif",
  "MP5"  = "IFN_response",
  "MP7"  = "S_cycle",
  "MP8"  = "Intestinal_diff",
  "MP9"  = "G1S_cycle",
  "MP10" = "Columnar_diff",
  "MP12" = "Neuro_epithelial",
  "MP13" = "Partial_EMT",
  "MP14" = "Hypoxia_epithelial",
  "MP15" = "T_NK_infiltration",
  "MP16" = "Secretory_diff",
  "MP17" = "Squamous_transition",
  "MP18" = "Adaptive_secretory"
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

p_scatter <- ggplot(comp_df, aes(x = scRef_score, y = PDO_score, label = MP)) +
  geom_point(size = 3, alpha = 0.7) + geom_text_repel(size = 3, max.overlaps = 15) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(title = paste0("PDO vs scRef 3CA MP mean scores\nSpearman r = ", round(cor_val, 3)), x = "scRef mean score", y = "PDO mean score") +
  theme_minimal(base_size = 12)
ggsave("Auto_PDO_mp_correlation_crossdata_scatter.pdf", p_scatter, width = 8, height = 6, useDingbats = FALSE)

comp_df_long <- reshape2::melt(comp_df, id.vars = "MP", variable.name = "Dataset", value.name = "Score")
p_bar <- ggplot(comp_df_long, aes(x = MP, y = Score, fill = Dataset)) + geom_bar(position = "dodge", stat = "identity") +
  coord_flip() + theme_minimal(base_size = 11) + labs(title = "PDO vs scRef 3CA MP mean scores") +
  scale_fill_manual(values = c(PDO_score = "#E41A1C", scRef_score = "#377AB8"))
ggsave("Auto_PDO_mp_correlation_crossdata_bar.pdf", p_bar, width = 10, height = 8, useDingbats = FALSE)

####################
# Cross-correlation: PDO MPs vs scRef MPs in scRef cells (SAMPLE-AVERAGED)
####################
message("Computing sample-averaged cross-correlation...")

ucell_pdo_in_scref <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds")
ucell_scref_nmp19 <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds")

# Get MP rows from PDO - USE SAME FILTERED LIST
pdo_mp_rows <- pdo_mp_tree_order

# Get MP cols from scRef - USE SAME FILTERED LIST (MP1,2,5,7,8,9,10,12,13,14,15,16,17,18)
sc_mp_cols <- colnames(ucell_scref_nmp19)
sc_mp_cols <- sc_mp_cols[sc_mp_cols %in% sc_mp_tree_order]

message("PDO MPs for correlation: ", paste(pdo_mp_rows, collapse = ", "))
message("scRef MPs for correlation: ", paste(sc_mp_cols, collapse = ", "))

if (length(pdo_mp_rows) >= 3 && length(sc_mp_cols) >= 3) {
  # ucell_pdo_in_scref: rows=PDO MPs, cols=cells (already transposed)
  # ucell_scref_nmp19: rows=cells, cols=scRef MPs (need to transpose)
  mod_mat <- ucell_pdo_in_scref[pdo_mp_rows, , drop = FALSE]  # Already: MPs x cells
  ref_mat <- t(ucell_scref_nmp19[, sc_mp_cols, drop = FALSE])  # Transpose: MPs x cells
  
  # Rename with descriptions
  rownames(mod_mat) <- pdo_mp_descriptions[rownames(mod_mat)]
  rownames(ref_mat) <- sc_mp_descriptions[rownames(ref_mat)]
  
  # Align columns (cell barcodes)
  common_cells <- intersect(colnames(mod_mat), colnames(ref_mat))
  mod_mat <- mod_mat[, common_cells, drop = FALSE]
  ref_mat <- ref_mat[, common_cells, drop = FALSE]
  
  message("Common cells: ", ncol(mod_mat))
  
  extract_sample_id <- function(x) {
    ifelse(
      grepl("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", x),
      sub("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", "", x),
      x
    )
  }
  
  sample_ids <- extract_sample_id(common_cells)
  sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]
  
  sample_counts <- table(sample_ids)
  samples <- names(sample_counts)
  
  message("Number of samples: ", length(samples))
  
  # Per-sample correlation: compute correlation for each sample, then average
  # After transpose: rows = MPs, cols = cells
  n_pdo <- nrow(mod_mat)  # 9 PDO MPs
  n_sc <- nrow(ref_mat)   # 14 scRef MPs
  
  # Store per-sample correlations
  cor_list <- list()
  for (smp in samples) {
    idx <- which(sample_ids == smp)
    if (length(idx) >= 10) {
      # Correlations between each PDO MP and each scRef MP within this sample
      # mod_mat: rows=PDO MPs, cols=cells; ref_mat: rows=scRef MPs, cols=cells
      pdo_scores <- mod_mat[, idx, drop = FALSE]   # 9 x n_cells
      ref_scores <- ref_mat[, idx, drop = FALSE]   # 14 x n_cells
      
      # Compute Spearman correlation for each pair
      cors <- matrix(NA, n_pdo, n_sc)
      for (i in seq_len(n_pdo)) {
        for (j in seq_len(n_sc)) {
          if (sd(pdo_scores[i,]) > 0 && sd(ref_scores[j,]) > 0) {
            cors[i, j] <- cor(pdo_scores[i,], ref_scores[j,], method = "spearman")
          }
        }
      }
      rownames(cors) <- rownames(mod_mat)
      colnames(cors) <- rownames(ref_mat)
      cor_list[[smp]] <- cors
    }
  }
  
  # Stack all sample correlations and compute meta-analysis
  allCors <- abind::abind(cor_list, along = 3)
  
  # Fisher z-transformation, then mean, then back-transform
  z_scores <- atanh(pmin(pmax(allCors, -0.999), 0.999))
  
  # Mean z across samples, then convert back to r
  mean_z <- apply(z_scores, c(1, 2), mean, na.rm = TRUE)
  mean_rho <- tanh(mean_z)
  
  # Compute p-value: one-sample t-test if multiple samples, else use individual p-values
  if (length(samples) > 1) {
    # t-test that mean correlation is different from 0
    p_vals <- matrix(1, nrow = n_pdo, ncol = n_sc)
    for (i in seq_len(n_pdo)) {
      for (j in seq_len(n_sc)) {
        z_vals <- z_scores[i, j, ]
        z_vals <- z_vals[!is.na(z_vals)]
        if (length(z_vals) >= 3) {
          # One-sample t-test against 0
          test_res <- t.test(z_vals, mu = 0)
          p_vals[i, j] <- test_res$p.value
        }
      }
    }
  } else {
    # Use approximation for single sample
    p_vals <- 2 * pnorm(-abs(mean_rho * sqrt(ncol(mod_mat))))
  }
  
  dimnames(mean_rho) <- list(rownames(mod_mat), rownames(ref_mat))
  dimnames(p_vals) <- dimnames(mean_rho)
  
  # Get the correct order from the descriptions (matching Jaccard plot)
  # Jaccard uses: pdo_mp_descriptions[pdo_mp_tree_order] and sc_mp_descriptions[sc_mp_tree_order]
  pdo_order_for_plot <- pdo_mp_descriptions[pdo_mp_tree_order[pdo_mp_tree_order %in% names(pdo_mp_descriptions) & !pdo_mp_tree_order %in% drop_mps_pdo]]
  sc_order_for_plot <- sc_mp_descriptions[sc_mp_tree_order[sc_mp_tree_order %in% names(sc_mp_descriptions) & !sc_mp_tree_order %in% drop_mps_sc]]
  
  # Filter to only include MPs that exist in mean_rho
  pdo_order_for_plot <- pdo_order_for_plot[pdo_order_for_plot %in% rownames(mean_rho)]
  sc_order_for_plot <- sc_order_for_plot[sc_order_for_plot %in% colnames(mean_rho)]
  
  # Reorder mean_rho: rows=PDO, cols=scRef
  mean_rho_ordered <- mean_rho[pdo_order_for_plot, sc_order_for_plot, drop = FALSE]
  p_vals_ordered <- p_vals[pdo_order_for_plot, sc_order_for_plot, drop = FALSE]
  
  # Plot meta-analysis heatmap
  # mean_rho_ordered: rows = PDO MPs, cols = scRef MPs
  # After t(): rows = scRef MPs, cols = PDO MPs
  # In cell_fun: i = scRef (row after transpose), j = PDO (col after transpose)
  # So need to access p_vals_ordered[j, i] not p_vals_ordered[i, j]
  col_cor <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  pdf("Auto_PDO_mp_crossref_correlation_meta.pdf", width = 14, height = 10, useDingbats = FALSE)
  ht <- Heatmap(t(mean_rho_ordered), name = "Mean Spearman\n(Meta-analysis)", col = col_cor,
    cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 1),
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- p_vals_ordered[j, i]; r_val <- mean_rho_ordered[j, i]
      lvl <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      grid.text(sprintf("%.2f\n%s", r_val, lvl), x, y, gp = gpar(fontsize = 9, fontface = "bold"))
    }, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
  draw(ht)
  dev.off()
  
  message("Meta-analysis complete with significance levels")
}

message("=== DONE ===")
