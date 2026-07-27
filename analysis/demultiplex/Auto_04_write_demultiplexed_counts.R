#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0 || hit[length(hit)] == length(args)) return(default)
  args[hit[length(hit)] + 1]
}

pool <- get_arg("--pool")
ephemeral_root <- get_arg("--ephemeral_root", "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate")
live_root <- get_arg("--live_root", "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/demultiplex")
min_margin <- as.numeric(get_arg("--min_margin", "0"))

####################
# A pool without WES/VCF references can still be exported safely as four
# temporary, cluster-labelled matrices.  These labels deliberately record
# that they are Souporcell clusters, rather than biological donor identities.
temporary_cluster_prefix <- get_arg("--temporary_cluster_prefix", "")
counts_out_dir_override <- get_arg("--counts_out_dir", "")
temporary_cluster_mode <- nzchar(temporary_cluster_prefix)
####################

if (is.null(pool) || !nzchar(pool)) stop("Usage: Rscript Auto_04_write_demultiplexed_counts.R --pool PDOs_Untreated")

condition <- sub("^PDOs_", "", pool)
matrix_dir <- file.path(ephemeral_root, "cellranger", pool, "outs", "filtered_feature_bc_matrix")
clusters_path <- file.path(ephemeral_root, "souporcell", pool, "clusters.tsv")
key_path <- file.path(live_root, "genotype_assignment", pool, paste0("Auto_", pool, "_cluster_to_donor_key.tsv"))
####################
counts_out_dir <- if (nzchar(counts_out_dir_override)) {
  counts_out_dir_override
} else {
  file.path(live_root, "counts_csv", pool)
}
####################
audit_out_dir <- file.path(live_root, "assignment_audit", pool)

if (!dir.exists(matrix_dir)) stop("Missing CellRanger matrix directory: ", matrix_dir)
if (!file.exists(clusters_path)) stop("Missing Souporcell clusters.tsv: ", clusters_path)
####################
if (!temporary_cluster_mode && !file.exists(key_path)) {
  stop("Missing genotype assignment key: ", key_path)
}
####################

dir.create(counts_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_out_dir, recursive = TRUE, showWarnings = FALSE)

clusters <- fread(clusters_path)
####################
if (temporary_cluster_mode) {
  cluster_ids <- sort(unique(unlist(strsplit(as.character(clusters$assignment), "/", fixed = TRUE))))
  cluster_ids <- cluster_ids[nzchar(cluster_ids)]
  cluster_to_donor <- setNames(
    paste0(temporary_cluster_prefix, cluster_ids),
    cluster_ids
  )
  cluster_to_margin <- setNames(rep(NA_real_, length(cluster_ids)), cluster_ids)
} else {
  key <- fread(key_path)
  required_key <- c("Cluster_ID", "Genotype_ID", "Correlation_Margin")
  if (!all(required_key %in% colnames(key))) stop("Assignment key missing columns: ", paste(setdiff(required_key, colnames(key)), collapse = ", "))

  key[, Cluster_ID := as.character(Cluster_ID)]
  key[, Genotype_ID := as.character(Genotype_ID)]
  key[!is.finite(Correlation_Margin), Correlation_Margin := NA_real_]
  cluster_to_donor <- setNames(key$Genotype_ID, key$Cluster_ID)
  cluster_to_margin <- setNames(key$Correlation_Margin, key$Cluster_ID)
}
####################

map_assignment <- function(x, lookup) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  mapped <- lookup[parts]
  if (any(is.na(mapped))) return(NA_character_)
  paste(mapped, collapse = "/")
}

clusters[, donor := vapply(assignment, map_assignment, character(1), lookup = cluster_to_donor)]
clusters[, assignment_margin := vapply(as.character(assignment), function(x) {
  parts <- strsplit(x, "/", fixed = TRUE)[[1]]
  vals <- cluster_to_margin[parts]
  if (any(is.na(vals))) return(NA_real_)
  min(vals)
}, numeric(1))]

clusters[, keep_for_counts := status == "singlet" &
           !is.na(donor) &
####################
           (temporary_cluster_mode | grepl("^SUR[0-9]+$", donor)) &
####################
           (is.na(assignment_margin) | assignment_margin >= min_margin)]

counts <- Read10X(data.dir = matrix_dir)
if (is.list(counts)) counts <- counts[[1]]

barcode_assignment <- clusters[, .(
  barcode,
  status,
  assignment,
  donor,
  assignment_margin,
  keep_for_counts
)]
fwrite(barcode_assignment, file.path(audit_out_dir, paste0("Auto_", pool, "_barcode_assignment.csv")))

sample_summary <- clusters[, .N, by = .(status, donor, keep_for_counts)]
setorder(sample_summary, status, donor)
fwrite(sample_summary, file.path(audit_out_dir, paste0("Auto_", pool, "_sample_assignment_summary.csv")))

kept <- clusters[keep_for_counts == TRUE]
if (nrow(kept) == 0) stop("No singlet donor-assigned cells passed filters for ", pool)

for (donor_id in sort(unique(kept$donor))) {
  barcodes <- kept[donor == donor_id, barcode]
  barcodes <- intersect(barcodes, colnames(counts))
####################
  sample_name <- if (temporary_cluster_mode) {
    paste0(donor_id, "_PDO")
  } else {
    paste0(donor_id, "_", condition, "_PDO")
  }
####################
  output_file <- file.path(counts_out_dir, paste0(sample_name, ".csv"))

  if (length(barcodes) == 0) {
    warning("No matrix barcodes found for ", sample_name)
    next
  }

  message("Writing ", sample_name, ": ", length(barcodes), " cells")
  sample_counts <- counts[, barcodes, drop = FALSE]
  fwrite(as.matrix(sample_counts), output_file, row.names = TRUE)
}

writeLines(
  c(
    paste0("pool\t", pool),
    paste0("matrix_dir\t", matrix_dir),
    paste0("clusters_path\t", clusters_path),
    paste0("key_path\t", key_path),
    paste0("counts_out_dir\t", counts_out_dir),
    paste0("min_margin\t", min_margin),
####################
    paste0("temporary_cluster_prefix\t", temporary_cluster_prefix),
    paste0("temporary_cluster_mode\t", temporary_cluster_mode)
####################
  ),
  file.path(audit_out_dir, paste0("Auto_", pool, "_counts_export_summary.tsv"))
)
