#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
demux_root <- "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex"
out_dir <- file.path(out_root, "Auto_demultiplex_verification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(
  "SUR1070_Untreated_PDO",
  "SUR1090_Untreated_PDO",
  "SUR1072_Untreated_PDO",
  "SUR1121_Untreated_PDO",
  "SUR1141_Untreated_PDO",
  "SUR1181_Untreated_PDO"
)

sample_paths <- file.path(out_root, "by_samples", samples, paste0(samples, ".rds"))
names(sample_paths) <- samples
missing_samples <- names(sample_paths)[!file.exists(sample_paths)]
if (length(missing_samples) > 0) stop("Missing sample RDS files: ", paste(missing_samples, collapse = ", "))

objs <- lapply(sample_paths, readRDS)
names(objs) <- samples

sample_dims <- rbindlist(lapply(samples, function(sample) {
  obj <- objs[[sample]]
  data.table(
    sample = sample,
    cells = ncol(obj),
    genes = nrow(obj),
    orig_ident = paste(unique(obj$orig.ident), collapse = ";")
  )
}))
fwrite(sample_dims, file.path(out_dir, "Auto_existing_untreated_sample_dimensions.csv"))

common_genes <- Reduce(intersect, lapply(objs, rownames))
pseudobulk <- sapply(samples, function(sample) {
  counts <- get_counts(objs[[sample]])[common_genes, , drop = FALSE]
  lib_size <- Matrix::colSums(counts)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
  cpm <- Matrix::t(Matrix::t(counts) * (1e6 / lib_size))
  Matrix::rowMeans(log1p(cpm))
})

expr_cor <- as.data.table(round(cor(pseudobulk, method = "pearson"), 5), keep.rownames = "sample")
fwrite(expr_cor, file.path(out_dir, "Auto_existing_untreated_pseudobulk_pearson.csv"))

merged_path <- file.path(out_root, "PDOs_merged.rds")
if (file.exists(merged_path)) {
  merged <- readRDS(merged_path)
  emb <- Embeddings(merged, "umap")
  meta <- merged@meta.data
  centroids <- rbindlist(lapply(samples, function(sample) {
    cells <- rownames(meta)[meta$orig.ident == sample]
    vals <- colMeans(emb[cells, , drop = FALSE])
    data.table(sample = sample, umap_1 = vals[1], umap_2 = vals[2], cells = length(cells))
  }))
  fwrite(centroids, file.path(out_dir, "Auto_existing_untreated_umap_centroids.csv"))
  dist_mat <- as.matrix(dist(as.matrix(centroids[, .(umap_1, umap_2)])))
  rownames(dist_mat) <- centroids$sample
  colnames(dist_mat) <- centroids$sample
  fwrite(as.data.table(round(dist_mat, 5), keep.rownames = "sample"), file.path(out_dir, "Auto_existing_untreated_umap_centroid_distances.csv"))
}

cnv_paths <- file.path(out_root, "cnv", "by_samples", samples, paste0("Auto_", samples, "_infercna_outs_Carroll_2023.rds"))
names(cnv_paths) <- samples
if (all(file.exists(cnv_paths))) {
  cnv <- lapply(cnv_paths, readRDS)
  names(cnv) <- samples
  common_cnv_genes <- Reduce(intersect, lapply(cnv, rownames))
  cnv_profiles <- sapply(samples, function(sample) {
    Matrix::rowMeans(cnv[[sample]][common_cnv_genes, , drop = FALSE], na.rm = TRUE)
  })
  cnv_cor <- as.data.table(round(cor(cnv_profiles, method = "pearson"), 5), keep.rownames = "sample")
  fwrite(cnv_cor, file.path(out_dir, "Auto_existing_untreated_infercna_profile_pearson.csv"))
}

assignment_path <- file.path(demux_root, "assignment_audit", "PDOs_Untreated", "Auto_PDOs_Untreated_barcode_assignment.csv")
if (file.exists(assignment_path)) {
  assignment <- fread(assignment_path)
  assignment <- assignment[keep_for_counts == TRUE]
  assignment[, barcode_core := barcode]

  remap <- rbindlist(lapply(samples, function(sample) {
    old_barcodes <- sub(paste0("^", sample, "_"), "", colnames(objs[[sample]]))
    data.table(current_sample = sample, barcode_core = old_barcodes)
  }))
  remap <- merge(remap, assignment[, .(barcode_core, new_donor = donor, new_status = status)], by = "barcode_core", all.x = TRUE)
  remap_summary <- remap[, .N, by = .(current_sample, new_donor)]
  setorder(remap_summary, current_sample, -N)
  fwrite(remap, file.path(out_dir, "Auto_existing_vs_rerun_barcode_remap.csv"))
  fwrite(remap_summary, file.path(out_dir, "Auto_existing_vs_rerun_barcode_remap_summary.csv"))
} else {
  writeLines(
    paste0("No rerun barcode assignment found at ", assignment_path),
    file.path(out_dir, "Auto_existing_vs_rerun_barcode_remap_NOT_RUN.txt")
  )
}

focus <- data.table(
  metric = c("expression_pseudobulk_pearson", "infercna_profile_pearson"),
  SUR1121_vs_SUR1141 = c(
    cor(pseudobulk[, "SUR1121_Untreated_PDO"], pseudobulk[, "SUR1141_Untreated_PDO"]),
    NA_real_
  )
)

cnv_cor_file <- file.path(out_dir, "Auto_existing_untreated_infercna_profile_pearson.csv")
if (file.exists(cnv_cor_file)) {
  cnv_cor_dt <- fread(cnv_cor_file)
  focus[metric == "infercna_profile_pearson", SUR1121_vs_SUR1141 := cnv_cor_dt[sample == "SUR1121_Untreated_PDO", SUR1141_Untreated_PDO]]
}
fwrite(focus, file.path(out_dir, "Auto_existing_SUR1121_SUR1141_focus_summary.csv"))
