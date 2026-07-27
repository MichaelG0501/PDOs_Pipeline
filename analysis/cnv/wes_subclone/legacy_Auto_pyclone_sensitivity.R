####################
# Auto_pyclone_sensitivity.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_pyclone_sensitivity.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_input.tsv
#   - PDOs_outs/Auto_wes_subclone/conda_env/bin/pyclone-vi
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_pyclone_sensitivity_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_pyclone_sensitivity_cluster_summary.csv
#   - PDOs_outs/Auto_wes_subclone/intermediate/pyclone_sensitivity/<sample>/k<max_clusters>/
# Downstream use: terminal sensitivity check for PyClone-VI clone-count stability.
####################

suppressPackageStartupMessages({
  library(data.table)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
if (dir.exists(wd)) setwd(wd)

samples <- c("PDO_1090_vs_NT_1090", "PDO_1181_vs_NT_1181")
max_cluster_grid <- as.integer(strsplit(Sys.getenv("PYCLONE_SENSITIVITY_K", "5,10,20,40"), ",", fixed = TRUE)[[1]])
restarts <- as.integer(Sys.getenv("PYCLONE_SENSITIVITY_RESTARTS", "3"))
density <- Sys.getenv("PYCLONE_DENSITY", "beta-binomial")

out_root <- "PDOs_outs/Auto_wes_subclone"
pyclone_bin <- file.path(out_root, "conda_env/bin/pyclone-vi")
if (!file.exists(pyclone_bin)) stop("Missing PyClone-VI binary: ", pyclone_bin)

input_dir <- file.path(out_root, "tables/pyclone")
intermediate_dir <- file.path(out_root, "intermediate/pyclone_sensitivity")
table_dir <- file.path(out_root, "tables/visualisation")
log_dir <- file.path(out_root, "logs")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

min_ccf_sep <- function(x) {
  x <- sort(unique(round(x[is.finite(x)], 4)))
  if (length(x) < 2) return(NA_real_)
  min(diff(x))
}

parse_final_elbo <- function(log_file) {
  if (!file.exists(log_file)) return(NA_real_)
  lines <- readLines(log_file, warn = FALSE)
  elbo_lines <- grep("^Final ELBO:", lines, value = TRUE)
  if (length(elbo_lines) == 0) return(NA_real_)
  suppressWarnings(as.numeric(sub("^Final ELBO:[[:space:]]*", "", tail(elbo_lines, 1))))
}

run_one <- function(sample, max_clusters) {
  sample_input <- file.path(input_dir, paste0("Auto_", sample, "_pyclone_vi_input.tsv"))
  if (!file.exists(sample_input)) stop("Missing PyClone-VI input: ", sample_input)

  run_dir <- file.path(intermediate_dir, sample, paste0("k", max_clusters))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  h5 <- file.path(run_dir, paste0("Auto_", sample, "_k", max_clusters, "_r", restarts, ".h5"))
  result_file <- file.path(run_dir, paste0("Auto_", sample, "_k", max_clusters, "_r", restarts, "_results.tsv"))
  fit_log <- file.path(run_dir, paste0("Auto_", sample, "_k", max_clusters, "_r", restarts, "_fit.log"))
  write_log <- file.path(run_dir, paste0("Auto_", sample, "_k", max_clusters, "_r", restarts, "_write_results.log"))

  if (Sys.getenv("PDO_FORCE_REBUILD", "0") == "1" || !file.exists(result_file) || file.size(result_file) == 0) {
    fit_args <- c(
      "fit",
      "-i", sample_input,
      "-o", h5,
      "-c", as.character(max_clusters),
      "-d", density,
      "-r", as.character(restarts)
    )
    status <- system2(pyclone_bin, fit_args, stdout = fit_log, stderr = fit_log)
    if (!identical(status, 0L)) stop("PyClone-VI sensitivity fit failed for ", sample, " k=", max_clusters, ". See ", fit_log)

    write_args <- c("write-results-file", "-i", h5, "-o", result_file)
    status <- system2(pyclone_bin, write_args, stdout = write_log, stderr = write_log)
    if (!identical(status, 0L)) stop("PyClone-VI sensitivity write-results failed for ", sample, " k=", max_clusters, ". See ", write_log)
  }

  result <- fread(result_file)
  result[, cluster_id := as.character(cluster_id)]
  cluster_summary <- result[, .(
    n_mutations = .N,
    median_ccf = as.numeric(median(cellular_prevalence, na.rm = TRUE)),
    mean_ccf = as.numeric(mean(cellular_prevalence, na.rm = TRUE)),
    median_assignment_prob = as.numeric(median(cluster_assignment_prob, na.rm = TRUE)),
    mean_assignment_prob = as.numeric(mean(cluster_assignment_prob, na.rm = TRUE)),
    frac_assignment_prob_ge_0_8 = as.numeric(mean(cluster_assignment_prob >= 0.8, na.rm = TRUE))
  ), by = cluster_id]
  cluster_summary[, `:=`(
    sample = sample,
    max_clusters = max_clusters,
    restarts = restarts
  )]

  summary <- data.table(
    sample = sample,
    max_clusters = max_clusters,
    restarts = restarts,
    n_clusters_result = uniqueN(result$cluster_id),
    n_pyclone_variants = nrow(result),
    final_elbo = parse_final_elbo(fit_log),
    median_assignment_prob = as.numeric(median(result$cluster_assignment_prob, na.rm = TRUE)),
    mean_assignment_prob = as.numeric(mean(result$cluster_assignment_prob, na.rm = TRUE)),
    frac_assignment_prob_ge_0_8 = as.numeric(mean(result$cluster_assignment_prob >= 0.8, na.rm = TRUE)),
    min_cluster_n = min(cluster_summary$n_mutations),
    min_cluster_fraction = min(cluster_summary$n_mutations) / sum(cluster_summary$n_mutations),
    min_median_ccf_separation = min_ccf_sep(cluster_summary$median_ccf),
    result_file = result_file,
    fit_log = fit_log,
    write_log = write_log
  )

  list(summary = summary, clusters = cluster_summary)
}

runs <- list()
i <- 1L
for (sample in samples) {
  for (max_clusters in max_cluster_grid) {
    message("Running PyClone-VI sensitivity: ", sample, " max_clusters=", max_clusters)
    runs[[i]] <- run_one(sample, max_clusters)
    i <- i + 1L
  }
}

sensitivity_summary <- rbindlist(lapply(runs, `[[`, "summary"), use.names = TRUE)
sensitivity_clusters <- rbindlist(lapply(runs, `[[`, "clusters"), use.names = TRUE)
setorder(sensitivity_summary, sample, max_clusters)
setorder(sensitivity_clusters, sample, max_clusters, -median_ccf)

summary_path <- file.path(table_dir, "Auto_wes_subclone_pyclone_sensitivity_summary.csv")
cluster_path <- file.path(table_dir, "Auto_wes_subclone_pyclone_sensitivity_cluster_summary.csv")
fwrite(sensitivity_summary, summary_path)
fwrite(sensitivity_clusters, cluster_path)

fwrite(data.table(
  finished = as.character(Sys.time()),
  samples = paste(samples, collapse = ","),
  max_cluster_grid = paste(max_cluster_grid, collapse = ","),
  restarts = restarts,
  summary = summary_path,
  cluster_summary = cluster_path
), file.path(log_dir, "Auto_wes_subclone_pyclone_sensitivity_summary.tsv"), sep = "\t")

message("Wrote PyClone-VI sensitivity summary: ", summary_path)
