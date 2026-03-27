####################
# Auto_mp_correlation_in_sc.R
#
# Compare PDO MPs (nMP=13) against scRef MPs (nMP=19) 
# using scAtlas samples (scRef cells).
# Generates: Sample-averaged Spearman correlation heatmap (meta-analysis)
#
# Input:
#   /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds
#   /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds (PDO MP genes/metrics)
#   /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds (scRef MP genes/metrics)
#
# Output:
#   PDOs_outs/Auto_PDO_mp_crossref_in_sc_correlation_meta.pdf
####################

cat("=== Script Started ===\n")
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(data.table)
library(ggplot2)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

# Sample ID extraction function
extract_sample_id <- function(x) {
  ifelse(
    grepl("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", x),
    sub("_[ACGTN]+(?:[-._][A-Za-z0-9]+)*$", "", x),
    x
  )
}

message("=== Loading Metaprograms ===")
# PDO MPs (nMP=13)
MP_pdo <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")
# scRef MPs (nMP=19)
MP_sc <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

message("=== Loading pre-computed UCell scores (scRef cells) ===")
# PDO MPs in scRef cells (transposed: MPs x cells)
ucell_pdo_in_sc <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds")
# scRef MPs in scRef cells (cells x MPs)
ucell_sc_in_sc <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds")

####################
# MP descriptions - STRICT AS REQUESTED
####################
# State Groupings for ordering
state_groups_sc <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intest. Meta" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
)

state_groups_pdo <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "SMG-like Metaplasia"   = c("MP8")
)

pdo_mp_descriptions_raw <- c(
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

sc_mp_descriptions_raw <- c(
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

# Full labels for plotting
pdo_mp_descriptions <- setNames(
  paste("PDOs", names(pdo_mp_descriptions_raw), pdo_mp_descriptions_raw, sep = "_"),
  names(pdo_mp_descriptions_raw)
)

sc_mp_descriptions <- setNames(
  paste("scATLAS", names(sc_mp_descriptions_raw), sc_mp_descriptions_raw, sep = "_"),
  names(sc_mp_descriptions_raw)
)

####################
# Filter and Order PDO MPs
####################
pdo_sil <- MP_pdo$metaprograms.metrics$silhouette
names(pdo_sil) <- paste0("MP", seq_along(pdo_sil))
pdo_cov <- MP_pdo$metaprograms.metrics$sampleCoverage
names(pdo_cov) <- paste0("MP", seq_along(pdo_cov))

bad_pdo <- names(pdo_sil)[pdo_sil < 0]
low_cov_pdo <- names(pdo_cov)[pdo_cov < 0.25]
drop_mps_pdo <- unique(c(bad_pdo, low_cov_pdo))

pdo_order_ordered <- unlist(state_groups_pdo)
pdo_remaining <- setdiff(names(pdo_mp_descriptions), pdo_order_ordered)
pdo_mp_order <- c(pdo_order_ordered, pdo_remaining)

pdo_mps_to_use <- pdo_mp_order[pdo_mp_order %in% names(pdo_mp_descriptions)]
pdo_mps_to_use <- pdo_mps_to_use[!pdo_mps_to_use %in% drop_mps_pdo]

message("Selected PDO MPs (ordered): ", paste(pdo_mps_to_use, collapse = ", "))

####################
# Filter and Order scRef MPs
####################
sc_sil <- MP_sc$metaprograms.metrics$silhouette
names(sc_sil) <- paste0("MP", seq_along(sc_sil))
bad_sc <- names(sc_sil)[sc_sil < 0]

sc_order_ordered <- unlist(state_groups_sc)
sc_remaining <- setdiff(names(sc_mp_descriptions), sc_order_ordered)
sc_mp_order <- c(sc_order_ordered, sc_remaining)

sc_mps_to_use <- sc_mp_order[sc_mp_order %in% names(sc_mp_descriptions)]
sc_mps_to_use <- sc_mps_to_use[!sc_mps_to_use %in% bad_sc]

message("Selected scRef MPs (ordered): ", paste(sc_mps_to_use, collapse = ", "))

####################
# Align Data
####################
# ucell_pdo_in_sc: rows = MPs, cols = cells
# ucell_sc_in_sc: rows = cells, cols = MPs

common_cells <- intersect(colnames(ucell_pdo_in_sc), rownames(ucell_sc_in_sc))
message("Common cells (scRef): ", length(common_cells))

mod_mat <- ucell_pdo_in_sc[pdo_mps_to_use, common_cells, drop = FALSE]
ref_mat <- t(ucell_sc_in_sc[common_cells, sc_mps_to_use, drop = FALSE])

# Rename rows to full descriptions
rownames(mod_mat) <- pdo_mp_descriptions[rownames(mod_mat)]
rownames(ref_mat) <- sc_mp_descriptions[rownames(ref_mat)]

####################
# Meta-analysis Spearman Correlation
####################
message("=== Computing sample-averaged cross-correlation ===")

sample_ids <- extract_sample_id(colnames(mod_mat))
unique_samples <- unique(sample_ids)
# Filter out empty or NA samples
unique_samples <- unique_samples[!is.na(unique_samples) & unique_samples != ""]
message("Number of samples found: ", length(unique_samples))

n_pdo <- nrow(mod_mat)
n_sc <- nrow(ref_mat)

# Store per-sample correlations
cor_list <- list()
for (smp in unique_samples) {
  idx <- which(sample_ids == smp)
  if (length(idx) >= 10) {
    pdo_scores <- mod_mat[, idx, drop = FALSE]
    ref_scores <- ref_mat[, idx, drop = FALSE]
    
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

if (length(cor_list) == 0) {
  stop("No samples with >= 10 cells found for correlation.")
}

# Average correlations across samples
sum_z <- matrix(0, n_pdo, n_sc)
count_z <- matrix(0, n_pdo, n_sc)
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

# P-values: one-sample t-test across samples
p_vals <- matrix(1, n_pdo, n_sc)
for (i in seq_len(n_pdo)) {
  for (j in seq_len(n_sc)) {
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

dimnames(mean_rho) <- list(rownames(mod_mat), rownames(ref_mat))
dimnames(p_vals) <- dimnames(mean_rho)

####################
# Heatmap
####################
col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

pdf("Auto_PDO_mp_crossref_in_sc_correlation_meta.pdf", width = 14, height = 10, useDingbats = FALSE)
# Rows = scRef, Cols = PDO (to match layout of the PDO-cell plot)
ht <- Heatmap(t(mean_rho), name = "Mean Spearman\n(Meta-analysis)", col = col_cor,
  cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_vals[j, i]
    r_val <- mean_rho[j, i]
    lvl <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    grid.text(sprintf("%.2f\n%s", r_val, lvl), x, y, gp = gpar(fontsize = 9, fontface = "bold"))
  }, 
  row_names_gp = gpar(fontsize = 10), 
  column_names_gp = gpar(fontsize = 10),
  column_title = "PDO Metaprograms",
  row_title = "scRef Metaprogrammes (nMP=19)")
draw(ht)
dev.off()

message("=== Meta-analysis complete for scRef samples ===")
