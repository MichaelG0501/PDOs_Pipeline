#!/usr/bin/env Rscript

####################
# Parameterized genotyping_save.R-style assignment.
#
# This mirrors the official Souporcell/Demuxafy Assign_Indiv_by_Geno.R logic
# used in the old genotyping_save.R workflow: match reference and cluster VCF
# sites by CHROM:POS:REF:ALT, convert GT calls to dosages, calculate Pearson
# correlations between every donor and Souporcell cluster, then write both the
# full correlation matrix and reciprocal-best Genotype_ID_key output.
####################

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

calculate_gt_dosage <- function(dt, sample_cols) {
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

ref_gt <- fread(ref_gt_path, header = FALSE, sep = "\t", data.table = TRUE)
cluster_gt <- fread(cluster_gt_path, header = FALSE, sep = "\t", data.table = TRUE)

expected_ref_cols <- 4 + length(ref_samples)
expected_cluster_cols <- 4 + length(cluster_samples)
if (ncol(ref_gt) != expected_ref_cols) {
  stop("Reference genotype table has ", ncol(ref_gt), " columns; expected ", expected_ref_cols)
}
if (ncol(cluster_gt) != expected_cluster_cols) {
  stop("Cluster genotype table has ", ncol(cluster_gt), " columns; expected ", expected_cluster_cols)
}

setnames(ref_gt, c("chrom", "pos", "ref", "alt", ref_samples))
setnames(cluster_gt, c("chrom", "pos", "ref", "alt", cluster_samples))

ref_gt[, site_id := make_site_id(.SD), .SDcols = c("chrom", "pos", "ref", "alt")]
cluster_gt[, site_id := make_site_id(.SD), .SDcols = c("chrom", "pos", "ref", "alt")]
ref_gt <- ref_gt[!duplicated(site_id)]
cluster_gt <- cluster_gt[!duplicated(site_id)]

shared_sites <- intersect(ref_gt$site_id, cluster_gt$site_id)
if (length(shared_sites) < 1000) {
  stop("Too few shared SNP sites for genotype assignment: ", length(shared_sites))
}

setkey(ref_gt, site_id)
setkey(cluster_gt, site_id)
ref_gt <- ref_gt[shared_sites]
cluster_gt <- cluster_gt[shared_sites]

ref_dosage <- calculate_gt_dosage(ref_gt, ref_samples)
cluster_dosage <- calculate_gt_dosage(cluster_gt, cluster_samples)

pearson_correlations <- matrix(
  NA_real_,
  nrow = length(cluster_samples),
  ncol = length(ref_samples),
  dimnames = list(cluster_samples, ref_samples)
)
overlap_sites <- matrix(
  0L,
  nrow = length(cluster_samples),
  ncol = length(ref_samples),
  dimnames = list(cluster_samples, ref_samples)
)

for (cluster_id in cluster_samples) {
  cluster_vec <- cluster_dosage[, cluster_id]
  for (donor_id in ref_samples) {
    donor_vec <- ref_dosage[, donor_id]
    keep <- is.finite(cluster_vec) & is.finite(donor_vec)
    overlap_sites[cluster_id, donor_id] <- sum(keep)
    if (sum(keep) >= 100) {
      pearson_correlations[cluster_id, donor_id] <- suppressWarnings(
        cor(donor_vec[keep], cluster_vec[keep], method = "pearson", use = "complete.obs")
      )
    }
  }
}

correlation_out <- data.table(Cluster = rownames(pearson_correlations), as.data.table(pearson_correlations))
fwrite(correlation_out, file.path(outdir, paste0("Auto_", pool, "_ref_clust_pearson_correlations.tsv")), sep = "\t")

overlap_out <- data.table(Cluster = rownames(overlap_sites), as.data.table(overlap_sites))
fwrite(overlap_out, file.path(outdir, paste0("Auto_", pool, "_ref_clust_overlap_sites.tsv")), sep = "\t")

cluster_key <- rbindlist(lapply(rownames(pearson_correlations), function(cluster_id) {
  vals <- pearson_correlations[cluster_id, ]
  ord <- order(vals, decreasing = TRUE, na.last = NA)
  best <- if (length(ord) >= 1) colnames(pearson_correlations)[ord[1]] else NA_character_
  second <- if (length(ord) >= 2) colnames(pearson_correlations)[ord[2]] else NA_character_
  data.table(
    Cluster_ID = cluster_id,
    Genotype_ID = best,
    Correlation = if (!is.na(best)) vals[best] else NA_real_,
    Second_Genotype_ID = second,
    Second_Correlation = if (!is.na(second)) vals[second] else NA_real_,
    Correlation_Margin = if (!is.na(best) && !is.na(second)) vals[best] - vals[second] else NA_real_,
    Sites_Compared = if (!is.na(best)) overlap_sites[cluster_id, best] else NA_integer_
  )
}))

donor_key <- rbindlist(lapply(colnames(pearson_correlations), function(donor_id) {
  vals <- pearson_correlations[, donor_id]
  ord <- order(vals, decreasing = TRUE, na.last = NA)
  best <- if (length(ord) >= 1) rownames(pearson_correlations)[ord[1]] else NA_character_
  second <- if (length(ord) >= 2) rownames(pearson_correlations)[ord[2]] else NA_character_
  data.table(
    Genotype_ID = donor_id,
    Cluster_ID = best,
    Correlation = if (!is.na(best)) vals[best] else NA_real_,
    Second_Cluster_ID = second,
    Second_Correlation = if (!is.na(second)) vals[second] else NA_real_,
    Correlation_Margin = if (!is.na(best) && !is.na(second)) vals[best] - vals[second] else NA_real_,
    Sites_Compared = if (!is.na(best)) overlap_sites[best, donor_id] else NA_integer_
  )
}))

cluster_best <- setNames(cluster_key$Genotype_ID, cluster_key$Cluster_ID)
donor_key[, Reciprocal_Best := !is.na(Cluster_ID) & cluster_best[Cluster_ID] == Genotype_ID]
donor_key[, Cluster_ID_Reciprocal := fifelse(Reciprocal_Best, Cluster_ID, "unassigned")]

official_key <- donor_key[, .(
  Genotype_ID,
  Cluster_ID = Cluster_ID_Reciprocal,
  Correlation = fifelse(Reciprocal_Best, Correlation, NA_real_),
  Correlation_Margin,
  Sites_Compared,
  Reciprocal_Best
)]

fwrite(cluster_key, file.path(outdir, paste0("Auto_", pool, "_cluster_to_donor_key.tsv")), sep = "\t")
fwrite(donor_key, file.path(outdir, paste0("Auto_", pool, "_donor_to_cluster_key.tsv")), sep = "\t")
fwrite(official_key, file.path(outdir, paste0("Auto_", pool, "_Genotype_ID_key.txt")), sep = "\t")

writeLines(
  c(
    paste0("pool\t", pool),
    paste0("method\tgenotyping_save_style_gt_pearson"),
    paste0("reference_gt\t", ref_gt_path),
    paste0("cluster_gt\t", cluster_gt_path),
    paste0("shared_sites\t", length(shared_sites)),
    paste0("reference_samples\t", paste(ref_samples, collapse = ",")),
    paste0("cluster_samples\t", paste(cluster_samples, collapse = ","))
  ),
  file.path(outdir, paste0("Auto_", pool, "_assignment_summary.tsv"))
)

png(
  file.path(outdir, paste0("Auto_", pool, "_ref_clust_pearson_correlation.png")),
  width = 1200,
  height = 900,
  res = 150
)
heatmap(
  pearson_correlations,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  col = colorRampPalette(c("white", "red"))(101),
  margins = c(10, 10),
  main = paste(pool, "genotyping_save-style correlation")
)
dev.off()
