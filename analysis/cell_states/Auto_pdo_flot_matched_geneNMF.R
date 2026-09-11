####################
# Analysis registry:
#   Status: active upstream; centred high-resolution matched-FLOT GeneNMF
#   Script: analysis/cell_states/Auto_pdo_flot_matched_geneNMF.R
#   Methodology: analysis/methodology/cell_states/Auto_pdo_flot_centred_highres_metaprogram_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Runs GeneNMF on the eight matched untreated/FLOT-treated PDO samples.
#     As in the canonical centred route, multiNMF uses center=TRUE so each
#     gene is centred before negative values are set to zero. This script
#     preserves the individual NMF programmes for the deliberately
#     high-resolution metaprogram extraction and paired trend filter.
#   Inputs:
#     - live: PDOs_outs/PDOs_list_PDOs.rds
#   Outputs:
#     - live intermediate: PDOs_outs/Auto_pdo_flot_centred_highres_metaprogram_trends/intermediate/Auto_pdo_flot_matched_centred_geneNMF_outs.rds
#     - live tables: matched-sample and NMF-programme summaries
#     - live logs: run summary and session information
#   Downstream use:
#     - Auto_pdo_flot_matched_highres_mp_trend_filter.R reads the programme
#       object and extracts nMP=round(total NMF programmes / 2).
#   Cache/replot behavior:
#     - Reuses the programme object unless --force or PDO_FORCE_REBUILD=1.
#     - PDO_REPLOT_ONLY is not applicable because this script has no plots.
#   Run command:
#     qsub -l select=1:ncpus=8:mem=120gb -l walltime=24:00:00 -N flot_hr_nmf -koed -- /bin/bash -lc 'module purge; module load tools/dev; eval "$(~/miniforge3/bin/conda shell.bash hook)"; conda activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf; cd /rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline; Rscript analysis/cell_states/Auto_pdo_flot_matched_geneNMF.R'
#   Conda env: gnmf
####################

library(GeneNMF)
library(Seurat)

####################
# Configuration and persistent paths
####################
args <- commandArgs(trailingOnly = TRUE)
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

out_dir <- file.path(PDO_LIVE_OUTS, "Auto_pdo_flot_centred_highres_metaprogram_trends")
tier_dirs <- pdo_ensure_output_tiers(out_dir)
input_path <- file.path(PDO_LIVE_OUTS, "PDOs_list_PDOs.rds")
program_path <- file.path(
  tier_dirs[["intermediate"]],
  "Auto_pdo_flot_matched_centred_geneNMF_outs.rds"
)
sample_summary_path <- file.path(
  tier_dirs[["tables"]],
  "Auto_pdo_flot_matched_centred_geneNMF_sample_summary.csv"
)
program_summary_path <- file.path(
  tier_dirs[["tables"]],
  "Auto_pdo_flot_matched_centred_geneNMF_program_summary.csv"
)

patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
matched_samples <- as.vector(rbind(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
))
force_rebuild <- "--force" %in% args || pdo_get_env_flag("PDO_FORCE_REBUILD", FALSE)
start_time <- Sys.time()

program_count <- function(program_object) {
  if (is.list(program_object) && !is.null(program_object$w)) {
    return(ncol(program_object$w))
  }
  if (is.matrix(program_object)) {
    return(ncol(program_object))
  }
  NA_integer_
}

####################
# Load and validate the matched PDO input
####################
pdo_require_files(input_path)
message("Loading matched PDO input: ", input_path)
pdos_list <- readRDS(input_path)
pdos_list[[PDO_EXCLUDED_SAMPLE]] <- NULL

missing_samples <- setdiff(matched_samples, names(pdos_list))
if (length(missing_samples) > 0L) {
  stop(
    "Matched sample(s) missing from PDOs_list_PDOs.rds: ",
    paste(missing_samples, collapse = ", ")
  )
}
pdos_list <- pdos_list[matched_samples]
if (length(pdos_list) != length(matched_samples)) {
  stop("Matched PDO list did not retain exactly eight samples.")
}

sample_summary <- data.frame(
  sample = names(pdos_list),
  patient = sub("_(Untreated|Treated)_PDO$", "", names(pdos_list)),
  treatment = ifelse(grepl("_Treated_", names(pdos_list)), "Treated", "Untreated"),
  n_cells = vapply(pdos_list, ncol, numeric(1)),
  n_genes = vapply(pdos_list, nrow, numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(sample_summary, sample_summary_path, row.names = FALSE)

####################
# Centred GeneNMF
####################
if (file.exists(program_path) && !force_rebuild) {
  message("Reusing centred GeneNMF programme cache: ", program_path)
  gene_nmf_programs <- readRDS(program_path)
} else {
  message("Running centred multiNMF on ", length(pdos_list), " matched PDO samples.")
  gene_nmf_programs <- GeneNMF::multiNMF(
    pdos_list,
    assay = "RNA",
    k = 4:9,
    min.exp = 0.05,
    center = TRUE
  )
  saveRDS(gene_nmf_programs, program_path, compress = FALSE)
}

expected_fit_n <- length(pdos_list) * length(4:9)
if (length(gene_nmf_programs) != expected_fit_n) {
  stop(
    "multiNMF returned ", length(gene_nmf_programs),
    " sample-rank fits; expected ", expected_fit_n, "."
  )
}
if (is.null(names(gene_nmf_programs)) || any(!grepl("\\.k[4-9]$", names(gene_nmf_programs)))) {
  stop("multiNMF sample-rank fits do not have the expected '<sample>.k<rank>' names.")
}
program_summary <- data.frame(
  nmf_fit = names(gene_nmf_programs),
  sample = sub("\\.k[4-9]$", "", names(gene_nmf_programs)),
  k = suppressWarnings(as.integer(sub("^.*\\.k", "", names(gene_nmf_programs)))),
  n_programmes = vapply(gene_nmf_programs, program_count, numeric(1)),
  stringsAsFactors = FALSE
)
if (anyNA(program_summary$n_programmes) || any(program_summary$n_programmes < 1L)) {
  stop("At least one matched sample has an invalid NMF programme count.")
}
program_summary$total_nmf_programmes <- sum(program_summary$n_programmes)
program_summary$center <- TRUE
program_summary$k_range <- "4:9"
write.csv(program_summary, program_summary_path, row.names = FALSE)

pdo_write_run_summary(
  script = "analysis/cell_states/Auto_pdo_flot_matched_geneNMF.R",
  out_dir = out_dir,
  inputs = input_path,
  outputs = c(program_path, sample_summary_path, program_summary_path),
  parameters = list(
    matched_sample_n = length(matched_samples),
    sample_rank_fit_n = length(gene_nmf_programs),
    k = "4:9",
    min_exp = 0.05,
    center = TRUE,
    total_nmf_programmes = unique(program_summary$total_nmf_programmes),
    elapsed_minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)
  ),
  cache = list(force_rebuild = force_rebuild, cache_reused = file.exists(program_path) && !force_rebuild)
)

message("Centred matched-FLOT GeneNMF completed.")
message("Programme object: ", program_path)
message("Total NMF programmes: ", unique(program_summary$total_nmf_programmes))
