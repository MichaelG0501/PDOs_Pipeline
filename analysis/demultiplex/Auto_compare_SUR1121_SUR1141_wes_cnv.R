#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_dir <- file.path(root_dir, "PDOs_outs", "Auto_demultiplex_verification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  SUR1121 = "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit/PDO_1121_vs_NT_1121/PDO_1121.cns",
  SUR1141 = "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit/PDO_1141_vs_NT_1141/PDO_1141.cns"
)

if (!all(file.exists(paths))) {
  stop("Missing CNVkit segment files: ", paste(paths[!file.exists(paths)], collapse = ", "))
}

chrom_levels <- paste0("chr", c(1:22, "X", "Y"))
segments <- lapply(paths, function(path) fread(path)[chromosome %in% chrom_levels])
max_end <- rbindlist(segments)[, .(max_end = max(end)), by = chromosome]

bins <- max_end[, .(start = seq(0, max_end, by = 1e6)), by = chromosome]
bins[, end := start + 1e6]
bins[, mid := start + 5e5]

assign_segment_profile <- function(seg) {
  out <- copy(bins)
  out[, log2 := NA_real_]
  for (chr in unique(out$chromosome)) {
    idx <- which(out$chromosome == chr)
    chr_seg <- seg[chromosome == chr]
    hits <- findInterval(out$mid[idx], chr_seg$start)
    ok <- hits > 0 & out$mid[idx] <= chr_seg$end[pmax(hits, 1)]
    out$log2[idx[ok]] <- chr_seg$log2[hits[ok]]
  }
  out$log2
}

profiles <- data.table(
  bin = paste(bins$chromosome, bins$start, bins$end, sep = ":"),
  chromosome = bins$chromosome,
  start = bins$start,
  end = bins$end,
  SUR1121 = assign_segment_profile(segments$SUR1121),
  SUR1141 = assign_segment_profile(segments$SUR1141)
)
profiles <- profiles[is.finite(SUR1121) & is.finite(SUR1141)]
profiles[, abs_diff := abs(SUR1121 - SUR1141)]
profiles[, SUR1121_major := abs(SUR1121) > 0.25]
profiles[, SUR1141_major := abs(SUR1141) > 0.25]

summary <- data.table(
  metric = c(
    "bins_compared",
    "pearson",
    "spearman",
    "mean_abs_diff",
    "median_abs_diff",
    "max_abs_diff",
    "major_cna_overlap_bins",
    "major_cna_SUR1121_only_bins",
    "major_cna_SUR1141_only_bins"
  ),
  value = c(
    nrow(profiles),
    cor(profiles$SUR1121, profiles$SUR1141, method = "pearson"),
    cor(profiles$SUR1121, profiles$SUR1141, method = "spearman"),
    mean(profiles$abs_diff),
    median(profiles$abs_diff),
    max(profiles$abs_diff),
    sum(profiles$SUR1121_major & profiles$SUR1141_major),
    sum(profiles$SUR1121_major & !profiles$SUR1141_major),
    sum(!profiles$SUR1121_major & profiles$SUR1141_major)
  )
)

fwrite(profiles, file.path(out_dir, "Auto_SUR1121_SUR1141_wes_cnv_1Mb_profiles.csv"))
fwrite(summary, file.path(out_dir, "Auto_SUR1121_SUR1141_wes_cnv_summary.csv"))

print(summary)
