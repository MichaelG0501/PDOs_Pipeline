####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/centred/Auto_01_centred_geneNMF.R
#   Description:
#     Replicates geneNMF.R but uses center=TRUE in multiNMF. This natively
#     transforms the log normalised matrix by making it centered per gene
#     (subtracting gene mean) and sets negative values to zero before running
#     NMF factorization. Extracts metaprograms for nMP=4:25.
#   Inputs:
#     - PDOs_outs/PDOs_list_PDOs.rds  (named list of post-QC Seurat objects)
#   Outputs:
#     - ephemeral: PDOs_outs/centred_mp_refinement/intermediate/geneNMF_outs.rds
#     - live: PDOs_outs/centred_mp_refinement/geneNMF_metaprograms_nMP_{k}.rds
#     - live: PDOs_outs/centred_mp_refinement/figures/metaprograms_heatmap_nMP_{k}.png
#   Conda env: gnmf
####################

library(GeneNMF)
library(RColorBrewer)
library(Seurat)

# === Paths ===
live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
ephemeral_base <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

outdir_live <- file.path(live_base, "centred_mp_refinement")
outdir_ephemeral <- file.path(ephemeral_base, "centred_mp_refinement", "intermediate")
fig_dir <- file.path(outdir_live, "figures")

dir.create(outdir_live, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_ephemeral, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# === Load input data ===
cat("Loading PDOs_list_PDOs.rds...\n")
pdos_list_path <- file.path(live_base, "PDOs_list_PDOs.rds")
if (!file.exists(pdos_list_path)) stop("Input file not found: ", pdos_list_path)
pdos.list <- readRDS(pdos_list_path)

# Exclude SUR843T3_PDO
pdos.list$SUR843T3_PDO <- NULL
cat("Loaded", length(pdos.list), "samples (after excluding SUR843T3_PDO)\n")

# === Step 1: multiNMF with center=TRUE ===
# Internally, GeneNMF:::getDataMatrix applies:
#   mat <- t(scale(Matrix::t(mat), center=TRUE, scale=FALSE))
#   mat[mat < 0] <- 0
rds_outs <- file.path(outdir_ephemeral, "geneNMF_outs.rds")
if (file.exists(rds_outs)) {
  cat("Loading existing geneNMF.programs...\n")
  geneNMF.programs <- readRDS(rds_outs)
} else {
  cat("Running multiNMF with center=TRUE, k=4:9...\n")
  geneNMF.programs <- multiNMF(pdos.list, assay = "RNA", k = 4:9,
                                min.exp = 0.05, center = TRUE)
  saveRDS(geneNMF.programs, file = rds_outs)
  cat("Saved:", rds_outs, "\n")
}

# === Step 2: Extract Metaprograms for nMP = 4 to 25 ===
k_vals <- 4:25
for (k in k_vals) {
  rds_path <- file.path(outdir_live, paste0("geneNMF_metaprograms_nMP_", k, ".rds"))
  png_path <- file.path(fig_dir, paste0("metaprograms_heatmap_nMP_", k, ".png"))

  if (file.exists(rds_path)) {
    cat(paste("nMP =", k, "already exists, skipping.\n"))
    next
  }

  cat(paste("Extracting Metaprograms for nMP =", k, "\n"))
  geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                           metric = "cosine",
                                           specificity.weight = 5,
                                           weight.explained = 0.5,
                                           nMP = k,
                                           min.confidence = 0.5)
  saveRDS(geneNMF.metaprograms, file = rds_path)

  if (!file.exists(png_path)) {
    cat(paste("Plotting heatmap for nMP =", k, "\n"))
    png(png_path, width = 3000, height = 2500, res = 300)
    plotMetaPrograms(geneNMF.metaprograms, similarity.cutoff = c(0, 1))
    dev.off()
  }

  cat(paste("Saved nMP =", k, "\n"))
}

cat("Auto_01_centred_geneNMF.R completed successfully.\n")
