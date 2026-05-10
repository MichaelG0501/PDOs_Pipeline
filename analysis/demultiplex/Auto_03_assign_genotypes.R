#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0 || hit[length(hit)] == length(args)) return(default)
  args[hit[length(hit)] + 1]
}

pool <- get_arg("--pool", "PDOs")
ref_gt_path <- get_arg("--ref_gt")
cluster_gt_path <- get_arg("--cluster_gt")
ref_samples_path <- get_arg("--ref_samples")
cluster_samples_path <- get_arg("--cluster_samples")
outdir <- get_arg("--outdir", ".")

if (is.null(ref_gt_path) || is.null(cluster_gt_path) ||
    is.null(ref_samples_path) || is.null(cluster_samples_path)) {
  stop("Required arguments: --ref_gt --cluster_gt --ref_samples --cluster_samples")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_samples <- function(path) {
  x <- readLines(path, warn = FALSE)
  x[nzchar(x)]
}

gt_to_dosage <- function(dt, sample_cols) {
  mat <- as.matrix(dt[, ..sample_cols])
  storage.mode(mat) <- "character"
  mat <- sub(":.*$", "", mat)
  mat <- gsub("\\|", "/", mat)
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  out[mat == "0/0"] <- 0
  out[mat == "0/1"] <- 1
  out[mat == "1/0"] <- 1
  out[mat == "1/1"] <- 2
  out
}

make_site_id <- function(dt) {
  paste(dt$chrom, dt$pos, dt$ref, dt$alt, sep = ":")
}

ref_samples <- read_samples(ref_samples_path)
cluster_samples <- read_samples(cluster_samples_path)

ref_dt <- fread(ref_gt_path, header = FALSE, sep = "\t", data.table = TRUE)
cluster_dt <- fread(cluster_gt_path, header = FALSE, sep = "\t", data.table = TRUE)

expected_ref_cols <- 4 + length(ref_samples)
expected_cluster_cols <- 4 + length(cluster_samples)
if (ncol(ref_dt) != expected_ref_cols) {
  stop("Reference GT table has ", ncol(ref_dt), " columns; expected ", expected_ref_cols)
}
if (ncol(cluster_dt) != expected_cluster_cols) {
  stop("Cluster GT table has ", ncol(cluster_dt), " columns; expected ", expected_cluster_cols)
}

setnames(ref_dt, c("chrom", "pos", "ref", "alt", ref_samples))
setnames(cluster_dt, c("chrom", "pos", "ref", "alt", cluster_samples))

ref_dt[, site_id := make_site_id(.SD), .SDcols = c("chrom", "pos", "ref", "alt")]
cluster_dt[, site_id := make_site_id(.SD), .SDcols = c("chrom", "pos", "ref", "alt")]
ref_dt <- ref_dt[!duplicated(site_id)]
cluster_dt <- cluster_dt[!duplicated(site_id)]

common_sites <- intersect(ref_dt$site_id, cluster_dt$site_id)
if (length(common_sites) < 1000) {
  stop("Too few shared SNP sites for genotype assignment: ", length(common_sites))
}

setkey(ref_dt, site_id)
setkey(cluster_dt, site_id)
ref_dt <- ref_dt[common_sites]
cluster_dt <- cluster_dt[common_sites]

ref_mat <- gt_to_dosage(ref_dt, ref_samples)
cluster_mat <- gt_to_dosage(cluster_dt, cluster_samples)

cor_mat <- matrix(
  NA_real_,
  nrow = length(cluster_samples),
  ncol = length(ref_samples),
  dimnames = list(cluster_samples, ref_samples)
)
n_mat <- matrix(
  0L,
  nrow = length(cluster_samples),
  ncol = length(ref_samples),
  dimnames = list(cluster_samples, ref_samples)
)

for (cluster_id in cluster_samples) {
  x <- cluster_mat[, cluster_id]
  for (donor_id in ref_samples) {
    y <- ref_mat[, donor_id]
    keep <- is.finite(x) & is.finite(y)
    n_mat[cluster_id, donor_id] <- sum(keep)
    if (sum(keep) >= 100) {
      cor_mat[cluster_id, donor_id] <- suppressWarnings(cor(x[keep], y[keep], method = "pearson"))
    }
  }
}

cor_out <- data.table(Cluster = rownames(cor_mat), as.data.table(cor_mat))
fwrite(cor_out, file.path(outdir, paste0("Auto_", pool, "_ref_clust_pearson_correlations.tsv")), sep = "\t")

n_out <- data.table(Cluster = rownames(n_mat), as.data.table(n_mat))
fwrite(n_out, file.path(outdir, paste0("Auto_", pool, "_ref_clust_overlap_sites.tsv")), sep = "\t")

cluster_key <- rbindlist(lapply(rownames(cor_mat), function(cluster_id) {
  vals <- cor_mat[cluster_id, ]
  ord <- order(vals, decreasing = TRUE, na.last = NA)
  best <- if (length(ord) >= 1) colnames(cor_mat)[ord[1]] else NA_character_
  second <- if (length(ord) >= 2) colnames(cor_mat)[ord[2]] else NA_character_
  data.table(
    Cluster_ID = cluster_id,
    Genotype_ID = best,
    Correlation = if (!is.na(best)) vals[best] else NA_real_,
    Second_Genotype_ID = second,
    Second_Correlation = if (!is.na(second)) vals[second] else NA_real_,
    Correlation_Margin = if (!is.na(best) && !is.na(second)) vals[best] - vals[second] else NA_real_,
    Sites_Compared = if (!is.na(best)) n_mat[cluster_id, best] else NA_integer_
  )
}))

donor_key <- rbindlist(lapply(colnames(cor_mat), function(donor_id) {
  vals <- cor_mat[, donor_id]
  ord <- order(vals, decreasing = TRUE, na.last = NA)
  best <- if (length(ord) >= 1) rownames(cor_mat)[ord[1]] else NA_character_
  second <- if (length(ord) >= 2) rownames(cor_mat)[ord[2]] else NA_character_
  data.table(
    Genotype_ID = donor_id,
    Cluster_ID = best,
    Correlation = if (!is.na(best)) vals[best] else NA_real_,
    Second_Cluster_ID = second,
    Second_Correlation = if (!is.na(second)) vals[second] else NA_real_,
    Correlation_Margin = if (!is.na(best) && !is.na(second)) vals[best] - vals[second] else NA_real_,
    Sites_Compared = if (!is.na(best)) n_mat[best, donor_id] else NA_integer_
  )
}))

cluster_best <- setNames(cluster_key$Genotype_ID, cluster_key$Cluster_ID)
donor_key[, Reciprocal_Best := !is.na(Cluster_ID) & cluster_best[Cluster_ID] == Genotype_ID]
donor_key[, Cluster_ID_Reciprocal := fifelse(Reciprocal_Best, Cluster_ID, "unassigned")]

fwrite(cluster_key, file.path(outdir, paste0("Auto_", pool, "_cluster_to_donor_key.tsv")), sep = "\t")
fwrite(donor_key, file.path(outdir, paste0("Auto_", pool, "_donor_to_cluster_key.tsv")), sep = "\t")

official_key <- donor_key[, .(
  Genotype_ID,
  Cluster_ID = Cluster_ID_Reciprocal,
  Correlation = fifelse(Reciprocal_Best, Correlation, NA_real_),
  Correlation_Margin,
  Sites_Compared,
  Reciprocal_Best
)]
fwrite(official_key, file.path(outdir, paste0("Auto_", pool, "_Genotype_ID_key.txt")), sep = "\t")

summary_lines <- c(
  paste0("pool\t", pool),
  paste0("reference_gt\t", ref_gt_path),
  paste0("cluster_gt\t", cluster_gt_path),
  paste0("shared_sites\t", length(common_sites)),
  paste0("reference_samples\t", paste(ref_samples, collapse = ",")),
  paste0("cluster_samples\t", paste(cluster_samples, collapse = ","))
)
writeLines(summary_lines, file.path(outdir, paste0("Auto_", pool, "_assignment_summary.tsv")))

png(file.path(outdir, paste0("Auto_", pool, "_ref_clust_pearson_correlation.png")),
    width = 1200, height = 900, res = 150)
heatmap(
  cor_mat,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  col = colorRampPalette(c("white", "red"))(101),
  margins = c(10, 10),
  main = paste(pool, "cluster-reference genotype correlation")
)
dev.off()
