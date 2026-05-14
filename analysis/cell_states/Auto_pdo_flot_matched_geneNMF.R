####################
# Auto_pdo_flot_matched_geneNMF.R
#
# Runs GeneNMF multiNMF on the eight matched untreated/FLOT-treated PDO
# samples used by legacy_pdo_flot_matched_survival_and_state_plots.R. This script only generates the
# sample-level NMF programme object; downstream high-resolution MP selection,
# UCell scoring, and paired trend filtering are handled by
# Auto_pdo_flot_matched_highres_mp_trend_filter.R.
#
# Env: gnmf
####################

suppressPackageStartupMessages({
  library(GeneNMF)
  library(Seurat)
})

####################
# setup
####################
args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_pdo_flot_highres_metaprogram_trends"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
matched_samples <- as.vector(rbind(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
))

geneNMF_program_path <- file.path(out_dir, "Auto_pdo_flot_matched_geneNMF_outs.rds")
sample_summary_path <- file.path(out_dir, "Auto_pdo_flot_matched_geneNMF_sample_summary.csv")
program_summary_path <- file.path(out_dir, "Auto_pdo_flot_matched_geneNMF_program_summary.csv")

program_count <- function(x) {
  if (is.list(x) && !is.null(x$w)) {
    return(ncol(x$w))
  }
  if (is.matrix(x)) {
    return(ncol(x))
  }
  NA_integer_
}

####################
# load matched samples
####################
if (!file.exists("PDOs_list_PDOs.rds")) {
  stop("Input file PDOs_list_PDOs.rds not found in PDOs_outs.")
}

message("Loading PDOs_list_PDOs.rds")
pdos.list <- readRDS("PDOs_list_PDOs.rds")

missing_samples <- setdiff(matched_samples, names(pdos.list))
if (length(missing_samples) > 0) {
  stop("Matched sample(s) missing from PDOs_list_PDOs.rds: ", paste(missing_samples, collapse = ", "))
}

pdos.list <- pdos.list[matched_samples]

sample_summary <- data.frame(
  sample = names(pdos.list),
  patient = sub("_(Untreated|Treated)_PDO$", "", names(pdos.list)),
  treatment = ifelse(grepl("_Treated_", names(pdos.list)), "Treated", "Untreated"),
  n_cells = vapply(pdos.list, ncol, numeric(1)),
  n_genes = vapply(pdos.list, nrow, numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(sample_summary, sample_summary_path, row.names = FALSE)

####################
# run multiNMF
####################
if (file.exists(geneNMF_program_path) && !force) {
  message("Existing GeneNMF programme object found; loading: ", geneNMF_program_path)
  geneNMF.programs <- readRDS(geneNMF_program_path)
} else {
  message("Running multiNMF on ", length(pdos.list), " matched PDO samples.")
  geneNMF.programs <- GeneNMF::multiNMF(
    pdos.list,
    assay = "RNA",
    k = 4:9,
    min.exp = 0.05
  )
  saveRDS(geneNMF.programs, geneNMF_program_path, compress = FALSE)
}

program_summary <- data.frame(
  sample = names(geneNMF.programs),
  n_programmes = vapply(geneNMF.programs, program_count, numeric(1)),
  stringsAsFactors = FALSE
)
program_summary$total_nmf_programmes <- sum(program_summary$n_programmes, na.rm = TRUE)
write.csv(program_summary, program_summary_path, row.names = FALSE)

message("GeneNMF programme object written to: ", geneNMF_program_path)
message("Total NMF programmes: ", unique(program_summary$total_nmf_programmes))
message("Auto_pdo_flot_matched_geneNMF.R completed successfully.")
