####################
# MP correlation analysis within PDOs
# Adapted from scRef_Pipeline/analysis/metaprograms/mp_correlation_sc.R
# and mp_correlation_crossdata.R for PDO-vs-PDO use
####################
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# Load objects
####################
tmdata_pdos <- readRDS("PDOs_final.rds")
MP_pdo <- readRDS("MP_outs_default.rds")

####################
# Filters: silhouette >= 0 and sample coverage >= 0.25 (>= 5/20 samples)
####################
metrics <- MP_pdo$metaprograms.metrics
rownames(metrics) <- paste0("MP", seq_len(nrow(metrics)))

sil_bad <- rownames(metrics)[metrics$silhouette < 0]
cov_bad <- rownames(metrics)[metrics$sampleCoverage < 0.25]
drop_mps <- unique(c(sil_bad, cov_bad))

keep_mps <- setdiff(paste0("MP", seq_len(nrow(metrics))), drop_mps)
keep_mps <- keep_mps[keep_mps %in% colnames(tmdata_pdos@meta.data)]

message("Silhouette-filtered MPs: ", paste(sil_bad, collapse = ", "))
message("Coverage-filtered MPs (<25%): ", paste(cov_bad, collapse = ", "))
message("Retained MPs: ", paste(keep_mps, collapse = ", "))

####################
# PDO order requirement: rev(mp_tree_order)
####################
tree_order <- MP_pdo$programs.tree$order
ordered_clusters <- MP_pdo$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order)]
mp_tree_order <- rev(mp_tree_order)
mp_tree_order <- paste0("MP", mp_tree_order)
mp_tree_order <- mp_tree_order[mp_tree_order %in% keep_mps]

####################
# MP description mapping provided by user
####################
mp_descriptions <- c(
  "MP6"  = "MP6_G2M_mitotic",
  "MP7"  = "MP7_DNA",
  "MP5"  = "MP5_MYC Biosynth",
  "MP1"  = "MP1_G2M_checkpoint",
  "MP3"  = "MP3_G1S_Cycle",
  "MP8"  = "MP8_Columnar progenitor",
  "MP10" = "MP10_Stress-induced plasticity",
  "MP9"  = "MP9_EMT_related",
  "MP4"  = "MP4_Intest diff"
)

####################
# Build MP matrix (MP x cell)
####################
mod_mat <- t(as.matrix(tmdata_pdos@meta.data[, mp_tree_order, drop = FALSE]))
rownames(mod_mat) <- ifelse(
  rownames(mod_mat) %in% names(mp_descriptions),
  mp_descriptions[rownames(mod_mat)],
  rownames(mod_mat)
)

####################
# 1) Global MP-MP Spearman correlation (all cells)
####################
cor_global <- cor(t(mod_mat), method = "spearman", use = "pairwise.complete.obs")

max_abs <- max(abs(cor_global), na.rm = TRUE)
if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
breaks <- seq(-max_abs, max_abs, length.out = 101)

pdf("mp_correlation_pdo_global_pheatmap.pdf", width = 8, height = 7)
pheatmap(
  cor_global,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = breaks,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 11,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "PDO MP-MP Spearman correlation (all cells)"
)
dev.off()

####################
# 2) Per-sample meta-correlation (Fisher Z + t-test)
####################
sample_ids <- tmdata_pdos$orig.ident
samples <- unique(sample_ids)

n_mps <- nrow(mod_mat)
cor_array <- array(
  NA_real_,
  dim = c(n_mps, n_mps, length(samples)),
  dimnames = list(rownames(mod_mat), rownames(mod_mat), samples)
)

for (smp in samples) {
  idx <- which(sample_ids == smp)
  if (length(idx) < 10) next
  sub_mat <- mod_mat[, idx, drop = FALSE]
  cor_array[, , smp] <- cor(t(sub_mat), method = "spearman", use = "pairwise.complete.obs")
}

z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))

mean_rho <- matrix(NA_real_, n_mps, n_mps, dimnames = list(rownames(mod_mat), rownames(mod_mat)))
p_vals <- matrix(NA_real_, n_mps, n_mps, dimnames = list(rownames(mod_mat), rownames(mod_mat)))

for (i in seq_len(n_mps)) {
  for (j in seq_len(n_mps)) {
    z_scores <- z_array[i, j, ]
    z_scores <- z_scores[is.finite(z_scores)]
    if (length(z_scores) >= 3) {
      mean_rho[i, j] <- tanh(mean(z_scores))
      ####################
      # Guard against constant vectors in t.test
      ####################
      if (sd(z_scores) == 0) {
        p_vals[i, j] <- 1
      } else {
        tt <- tryCatch(t.test(z_scores), error = function(e) NULL)
        p_vals[i, j] <- if (is.null(tt)) NA_real_ else tt$p.value
      }
    }
  }
}

if (all(is.na(mean_rho))) {
  stop("No valid per-sample correlations computed. Check sample sizes or MP matrix.")
}

col_cor <- colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))

pdf("mp_correlation_pdo_meta_complexheatmap.pdf", width = 9, height = 8)
ht <- Heatmap(
  mean_rho,
  name = paste0("Mean Rho\n(", length(samples), " samples)"),
  col = col_cor,
  rect_gp = gpar(col = "white", lwd = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_vals[i, j]
    r <- mean_rho[i, j]
    if (is.na(p) || is.na(r)) {
      grid.text("NA", x, y, gp = gpar(fontsize = 8, col = "grey50"))
    } else if (p < 0.001) {
      grid.text(paste0(sprintf("%.2f", r), "\n***"), x, y, gp = gpar(fontsize = 9))
    } else if (p < 0.01) {
      grid.text(paste0(sprintf("%.2f", r), "\n**"), x, y, gp = gpar(fontsize = 9))
    } else if (p < 0.05) {
      grid.text(paste0(sprintf("%.2f", r), "\n*"), x, y, gp = gpar(fontsize = 9))
    } else {
      grid.text(sprintf("%.2f", r), x, y, gp = gpar(fontsize = 9))
    }
  },
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold")
)
draw(ht, heatmap_legend_side = "right")
dev.off()

####################
# 3) Gene-set overlap Jaccard self-heatmap
####################
mp_genes <- MP_pdo$metaprograms.genes[keep_mps]
mp_genes <- mp_genes[mp_tree_order]
names(mp_genes) <- ifelse(names(mp_genes) %in% names(mp_descriptions), mp_descriptions[names(mp_genes)], names(mp_genes))
mp_genes <- lapply(mp_genes, unique)

universe <- unique(unlist(mp_genes))
n <- length(mp_genes)

jaccard_mat <- matrix(NA_real_, n, n, dimnames = list(names(mp_genes), names(mp_genes)))
overlap_n_mat <- jaccard_mat
pval_mat <- jaccard_mat

for (i in seq_len(n)) {
  A <- mp_genes[[i]]
  for (j in seq_len(n)) {
    B <- mp_genes[[j]]
    inter <- length(intersect(A, B))
    uni <- length(union(A, B))
    overlap_n_mat[i, j] <- inter
    jaccard_mat[i, j] <- if (uni == 0) NA_real_ else inter / uni

    a <- inter
    b <- length(setdiff(A, B))
    cc <- length(setdiff(B, A))
    d <- length(setdiff(universe, union(A, B)))
    pval_mat[i, j] <- if (any(c(a, b, cc, d) < 0)) NA_real_ else {
      fisher.test(matrix(c(a, b, cc, d), nrow = 2), alternative = "greater")$p.value
    }
  }
}

padj_mat <- matrix(
  p.adjust(as.vector(pval_mat), method = "BH"),
  nrow = nrow(pval_mat),
  ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)

stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat), dimnames = dimnames(padj_mat))
stars_mat[padj_mat < 0.05] <- "*"
stars_mat[padj_mat < 0.01] <- "**"
stars_mat[padj_mat < 0.001] <- "***"

display_mat <- matrix(
  paste0(overlap_n_mat, "\n", stars_mat),
  nrow = nrow(overlap_n_mat),
  ncol = ncol(overlap_n_mat),
  dimnames = dimnames(overlap_n_mat)
)

pdf("mp_jaccard_self_pdo.pdf", width = 8.5, height = 7.5)
pheatmap(
  jaccard_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "PDO MP gene-set overlap (Jaccard)",
  angle_col = "90",
  display_numbers = display_mat,
  fontsize_number = 8,
  number_color = "black",
  fontsize_row = 10,
  fontsize_col = 10,
  color = colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
)
dev.off()

####################
# Save matrices for downstream use
####################
saveRDS(
  list(
    retained_mps = mp_tree_order,
    retained_mp_labels = rownames(mod_mat),
    global_correlation = cor_global,
    meta_correlation = mean_rho,
    meta_p_values = p_vals,
    jaccard = jaccard_mat,
    overlap_n = overlap_n_mat,
    jaccard_padj = padj_mat
  ),
  file = "mp_correlation_pdo_outputs.rds"
)

message("Saved: mp_correlation_pdo_global_pheatmap.pdf")
message("Saved: mp_correlation_pdo_meta_complexheatmap.pdf")
message("Saved: mp_jaccard_self_pdo.pdf")
message("Saved: mp_correlation_pdo_outputs.rds")
