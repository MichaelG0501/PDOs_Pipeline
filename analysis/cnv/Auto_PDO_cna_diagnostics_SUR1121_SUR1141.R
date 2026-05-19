####################
# Quantify why SUR1121_Untreated_PDO and SUR1141_Untreated_PDO look nearly
# identical in expression-derived CNA profiles.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

out_dir <- file.path("cnv", "Auto_PDO_cna_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
focus_samples <- c("SUR1121_Untreated_PDO", "SUR1141_Untreated_PDO")
untreated_samples <- c(
  "SUR1070_Untreated_PDO",
  "SUR1090_Untreated_PDO",
  "SUR1072_Untreated_PDO",
  "SUR1121_Untreated_PDO",
  "SUR1141_Untreated_PDO",
  "SUR1181_Untreated_PDO"
)
treated_samples <- c(
  "SUR1070_Treated_PDO",
  "SUR1090_Treated_PDO",
  "SUR1072_Treated_PDO",
  "SUR1181_Treated_PDO"
)

pair_metrics <- function(mat) {
  samples <- colnames(mat)
  pairs <- combn(samples, 2, simplify = FALSE)
  bind_rows(lapply(pairs, function(pair) {
    x <- mat[, pair[1]]
    y <- mat[, pair[2]]
    delta <- x - y
    data.frame(
      sample_a = pair[1],
      sample_b = pair[2],
      n_features = sum(is.finite(x) & is.finite(y)),
      pearson = suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "pearson")),
      spearman = suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "spearman")),
      mean_abs_diff = mean(abs(delta), na.rm = TRUE),
      median_abs_diff = median(abs(delta), na.rm = TRUE),
      rmse = sqrt(mean(delta^2, na.rm = TRUE)),
      max_abs_diff = max(abs(delta), na.rm = TRUE),
      exactly_identical = identical(unname(x), unname(y)),
      stringsAsFactors = FALSE
    )
  }))
}

read_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

mean_log_cpm <- function(counts) {
  lib_size <- Matrix::colSums(counts)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
  cpm <- Matrix::t(Matrix::t(counts) * (1e6 / lib_size))
  Matrix::rowMeans(log2((cpm / 10) + 1))
}

make_bins <- function(feature_names, bin_size = 50L) {
  gene_order <- read.table(
    gene_order_path,
    header = FALSE,
    col.names = c("gene_id", "chromosome", "start", "end")
  )
  chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
  go <- gene_order %>%
    filter(.data$gene_id %in% feature_names, .data$chromosome %in% chrom_levels) %>%
    mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
    arrange(.data$chromosome, .data$start) %>%
    group_by(.data$chromosome) %>%
    mutate(
      g_rank = row_number(),
      bin_in_chr = ((.data$g_rank - 1L) %/% bin_size) + 1L,
      bin_key = paste(.data$chromosome, .data$bin_in_chr, sep = "_")
    ) %>%
    ungroup()
  split(go$gene_id, factor(go$bin_key, levels = unique(go$bin_key)))
}

bin_vector_matrix <- function(mat, bin_size = 50L) {
  bins <- make_bins(rownames(mat), bin_size = bin_size)
  out <- do.call(rbind, lapply(bins, function(genes) {
    colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(out) <- names(bins)
  out
}

sample_mean_path <- "cnv/Auto_PDO_cnv_sample_mean_profiles_Carroll_2023.rds"
heatmap_input_path <- "cnv/Auto_PDO_cnv_heatmap_all_samples_Carroll_2023_input.rds"
target_meta_path <- "cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds"

sample_mean <- readRDS(sample_mean_path)
target_meta <- readRDS(target_meta_path)

all_profile_metrics <- pair_metrics(sample_mean)
write.csv(all_profile_metrics, file.path(out_dir, "Auto_infercna_sample_mean_pair_metrics.csv"), row.names = FALSE)
write.csv(
  all_profile_metrics %>% filter(.data$sample_a %in% untreated_samples, .data$sample_b %in% untreated_samples),
  file.path(out_dir, "Auto_infercna_sample_mean_pair_metrics_untreated.csv"),
  row.names = FALSE
)
write.csv(
  all_profile_metrics %>% filter(.data$sample_a %in% treated_samples, .data$sample_b %in% treated_samples),
  file.path(out_dir, "Auto_infercna_sample_mean_pair_metrics_treated.csv"),
  row.names = FALSE
)

heatmap_input <- readRDS(heatmap_input_path)
binned_mat <- heatmap_input$binned_matrix
binned_meta <- heatmap_input$metadata
binned_sample_mean <- sapply(unique(as.character(binned_meta$sample)), function(sample) {
  cells <- rownames(binned_meta)[binned_meta$sample == sample]
  rowMeans(binned_mat[, cells, drop = FALSE], na.rm = TRUE)
})
write.csv(pair_metrics(binned_sample_mean), file.path(out_dir, "Auto_infercna_binned_pair_metrics.csv"), row.names = FALSE)

per_sample_dims <- bind_rows(lapply(focus_samples, function(sample) {
  path <- file.path("cnv", "by_samples", sample, paste0("Auto_", sample, "_infercna_outs_Carroll_2023.rds"))
  x <- readRDS(path)
  data.frame(
    sample = sample,
    n_genes = nrow(x),
    n_cells = ncol(x),
    md5_path = path,
    stringsAsFactors = FALSE
  )
}))
write.csv(per_sample_dims, file.path(out_dir, "Auto_focus_per_sample_infercna_dimensions.csv"), row.names = FALSE)

focus_delta <- sample_mean[, focus_samples[1]] - sample_mean[, focus_samples[2]]
top_focus_delta <- data.frame(
  gene_id = rownames(sample_mean),
  SUR1121_mean = sample_mean[, focus_samples[1]],
  SUR1141_mean = sample_mean[, focus_samples[2]],
  delta_1121_minus_1141 = focus_delta,
  abs_delta = abs(focus_delta),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(.data$abs_delta)) %>%
  head(200)
write.csv(top_focus_delta, file.path(out_dir, "Auto_SUR1121_SUR1141_top_infercna_gene_differences.csv"), row.names = FALSE)

sample_expr <- list()
sample_cells <- list()
for (sample in untreated_samples) {
  obj <- readRDS(file.path("by_samples", sample, paste0(sample, ".rds")))
  counts <- read_counts(obj)
  sample_expr[[sample]] <- mean_log_cpm(counts)
  sample_cells[[sample]] <- ncol(counts)
  rm(obj, counts)
  gc()
}

common_expr_genes <- Reduce(intersect, lapply(sample_expr, names))
expr_mat <- sapply(sample_expr, function(x) x[common_expr_genes])
rownames(expr_mat) <- common_expr_genes
write.csv(pair_metrics(expr_mat), file.path(out_dir, "Auto_raw_logcpm_pseudobulk_pair_metrics_untreated.csv"), row.names = FALSE)

expr_binned <- bin_vector_matrix(expr_mat)
write.csv(pair_metrics(expr_binned), file.path(out_dir, "Auto_raw_logcpm_genomic_binned_pair_metrics_untreated.csv"), row.names = FALSE)

sample_counts <- data.frame(
  sample = names(sample_cells),
  post_qc_cells = as.integer(unlist(sample_cells)),
  infercna_cells = as.integer(table(target_meta$sample)[names(sample_cells)]),
  stringsAsFactors = FALSE
)
write.csv(sample_counts, file.path(out_dir, "Auto_untreated_postqc_vs_infercna_cell_counts.csv"), row.names = FALSE)

focus_summary <- data.frame(
  metric = c(
    "infercna_sample_mean_pearson",
    "infercna_binned_heatmap_pearson",
    "raw_logcpm_pseudobulk_pearson",
    "raw_logcpm_genomic_binned_pearson",
    "infercna_sample_mean_exactly_identical",
    "infercna_sample_mean_mean_abs_diff",
    "infercna_sample_mean_max_abs_diff",
    "SUR1121_post_qc_cells",
    "SUR1141_post_qc_cells"
  ),
  value = c(
    suppressWarnings(cor(sample_mean[, focus_samples[1]], sample_mean[, focus_samples[2]], use = "pairwise.complete.obs")),
    suppressWarnings(cor(binned_sample_mean[, focus_samples[1]], binned_sample_mean[, focus_samples[2]], use = "pairwise.complete.obs")),
    suppressWarnings(cor(expr_mat[, focus_samples[1]], expr_mat[, focus_samples[2]], use = "pairwise.complete.obs")),
    suppressWarnings(cor(expr_binned[, focus_samples[1]], expr_binned[, focus_samples[2]], use = "pairwise.complete.obs")),
    identical(unname(sample_mean[, focus_samples[1]]), unname(sample_mean[, focus_samples[2]])),
    mean(abs(focus_delta), na.rm = TRUE),
    max(abs(focus_delta), na.rm = TRUE),
    sample_counts$post_qc_cells[sample_counts$sample == focus_samples[1]],
    sample_counts$post_qc_cells[sample_counts$sample == focus_samples[2]]
  ),
  stringsAsFactors = FALSE
)
write.csv(focus_summary, file.path(out_dir, "Auto_SUR1121_SUR1141_cna_diagnostic_summary.csv"), row.names = FALSE)

print(focus_summary)
message("Wrote diagnostics to: ", file.path(out_root, out_dir))
