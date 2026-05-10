#!/usr/bin/env Rscript

####################
# Compare the currently assigned SUR1121/SUR1141 barcode groups against the
# rerun Souporcell cluster probabilities and genotyping_save-style donor key.
####################

suppressPackageStartupMessages({
  library(data.table)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
pdos_out <- file.path(root_dir, "PDOs_outs")
demux_root <- "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex"
pool <- "PDOs_Untreated"
out_dir <- file.path(pdos_out, "Auto_demultiplex_verification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

current_samples <- c("SUR1121_Untreated_PDO", "SUR1141_Untreated_PDO")
current_csvs <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all",
  paste0(current_samples, ".csv")
)
names(current_csvs) <- current_samples

clusters_path <- file.path(demux_root, "souporcell", pool, "clusters.tsv")
cluster_key_path <- file.path(demux_root, "genotype_assignment", pool, paste0("Auto_", pool, "_cluster_to_donor_key.tsv"))
cor_path <- file.path(demux_root, "genotype_assignment", pool, paste0("Auto_", pool, "_ref_clust_pearson_correlations.tsv"))

if (!all(file.exists(current_csvs))) {
  stop("Missing current input CSVs: ", paste(current_csvs[!file.exists(current_csvs)], collapse = ", "))
}
if (!file.exists(clusters_path)) stop("Missing rerun clusters.tsv: ", clusters_path)
if (!file.exists(cluster_key_path)) stop("Missing rerun cluster key: ", cluster_key_path)
if (!file.exists(cor_path)) stop("Missing rerun genotype correlation matrix: ", cor_path)

read_csv_barcodes <- function(path) {
  header <- readLines(path, n = 1, warn = FALSE)
  fields <- strsplit(header, ",", fixed = TRUE)[[1]]
  fields <- gsub('^"|"$', "", fields)
  fields <- fields[nzchar(fields)]
  fields
}

current_barcodes <- rbindlist(lapply(names(current_csvs), function(sample) {
  data.table(
    current_sample = sample,
    barcode = read_csv_barcodes(current_csvs[[sample]])
  )
}))

clusters <- fread(clusters_path)
cluster_key <- fread(cluster_key_path)
cor_mat <- fread(cor_path)
setnames(cor_mat, "Cluster", "assignment")

cluster_cols <- grep("^cluster[0-9]+$", colnames(clusters), value = TRUE)
if (length(cluster_cols) == 0) stop("No cluster probability columns found in clusters.tsv")

clusters[, assignment := as.character(assignment)]
cluster_key[, Cluster_ID := as.character(Cluster_ID)]
cluster_key[, Genotype_ID := as.character(Genotype_ID)]

setkey(clusters, barcode)
matched <- clusters[current_barcodes, on = "barcode"]
####################
if ("i.current_sample" %in% colnames(matched)) {
  matched[, current_sample := i.current_sample]
  matched[, i.current_sample := NULL]
}
####################

matched[, rerun_cluster := assignment]
matched <- merge(
  matched,
  cluster_key[, .(rerun_cluster = Cluster_ID, rerun_genotype = Genotype_ID, rerun_cluster_correlation = Correlation, rerun_correlation_margin = Correlation_Margin)],
  by = "rerun_cluster",
  all.x = TRUE
)

prob_long <- melt(
  matched[, c("current_sample", "barcode", cluster_cols), with = FALSE],
  id.vars = c("current_sample", "barcode"),
  variable.name = "cluster_col",
  value.name = "log_probability"
)
prob_long[, rerun_cluster := sub("^cluster", "", as.character(cluster_col))]
prob_long <- merge(
  prob_long,
  cluster_key[, .(rerun_cluster = Cluster_ID, donor = Genotype_ID)],
  by = "rerun_cluster",
  all.x = TRUE
)

donor_probability_summary <- prob_long[
  !is.na(donor),
  .(
    mean_log_probability = mean(log_probability, na.rm = TRUE),
    median_log_probability = median(log_probability, na.rm = TRUE),
    cells_with_probability = sum(is.finite(log_probability))
  ),
  by = .(current_sample, donor)
]
setorder(donor_probability_summary, current_sample, -mean_log_probability)

assignment_summary <- matched[, .N, by = .(current_sample, status, rerun_cluster, rerun_genotype)]
setorder(assignment_summary, current_sample, -N)

sample_match_summary <- matched[, .(
  current_barcodes = .N,
  matched_to_rerun_clusters = sum(!is.na(rerun_cluster)),
  singlets = sum(status == "singlet", na.rm = TRUE),
  donor_SUR1121 = sum(rerun_genotype == "SUR1121", na.rm = TRUE),
  donor_SUR1141 = sum(rerun_genotype == "SUR1141", na.rm = TRUE),
  donor_other = sum(!is.na(rerun_genotype) & !rerun_genotype %in% c("SUR1121", "SUR1141")),
  unmatched = sum(is.na(rerun_cluster))
), by = current_sample]

fwrite(matched, file.path(out_dir, "Auto_current_SUR1121_SUR1141_barcodes_vs_rerun_clusters.csv"))
fwrite(assignment_summary, file.path(out_dir, "Auto_current_SUR1121_SUR1141_assignment_summary.csv"))
fwrite(sample_match_summary, file.path(out_dir, "Auto_current_SUR1121_SUR1141_sample_match_summary.csv"))
fwrite(donor_probability_summary, file.path(out_dir, "Auto_current_SUR1121_SUR1141_donor_probability_summary.csv"))

print(sample_match_summary)
print(donor_probability_summary[current_sample %in% current_samples & donor %in% c("SUR1121", "SUR1141")])
