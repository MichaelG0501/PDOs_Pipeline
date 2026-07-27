####################
# Auto_run_facets_sample.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_run_facets_sample.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - FACETS snp-pileup output: PDOs_outs/Auto_wes_subclone/intermediate/facets/<sample>/<sample>.snp_pileup.gz
#   - Validated Sarek launch/reference metadata written by Auto_validate_sarek_reference.sh
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/intermediate/facets/<sample>/Auto_<sample>_facets_fit.rds
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_segments.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_purity_ploidy.tsv
# Downstream use: allele-specific major/minor copy-number input for PyClone-VI.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(facets)
})

args <- commandArgs(trailingOnly = TRUE)
sample <- if (length(args) >= 1 && nzchar(args[1])) args[1] else Sys.getenv("sample")
out_root <- if (length(args) >= 2 && nzchar(args[2])) args[2] else Sys.getenv("OUT_ROOT")
cval <- as.integer(if (length(args) >= 3 && nzchar(args[3])) args[3] else Sys.getenv("FACETS_CVAL", "150"))
min_nhet <- as.integer(Sys.getenv("FACETS_MIN_NHET", "15"))

if (!nzchar(sample)) stop("Supply sample as first argument or environment variable sample.")
if (!nzchar(out_root)) {
  out_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
}

facets_dir <- file.path(out_root, "intermediate/facets", sample)
tables_dir <- file.path(out_root, "tables/facets")
logs_dir <- file.path(out_root, "logs")
dir.create(facets_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

pileup_file <- file.path(facets_dir, paste0(sample, ".snp_pileup.gz"))
if (!file.exists(pileup_file)) stop("Missing FACETS snp-pileup file: ", pileup_file)

message("Reading FACETS pileup: ", pileup_file)
rcmat <- readSnpMatrix(pileup_file)
if (nrow(rcmat) < min_nhet) {
  stop("Too few heterozygous SNP rows for FACETS in ", sample, ": ", nrow(rcmat), " < ", min_nhet)
}

message("Running FACETS preProcSample/procSample/emcncf with cval=", cval)
preprocessed <- preProcSample(rcmat, cval = cval)
processed <- procSample(preprocessed, cval = cval)
fit <- emcncf(processed)

purity <- suppressWarnings(as.numeric(fit$purity))
ploidy <- suppressWarnings(as.numeric(fit$ploidy))
if (!is.finite(purity) || !is.finite(ploidy)) {
  stop("FACETS did not return finite purity/ploidy for ", sample)
}

cncf <- as.data.table(fit$cncf)
if (nrow(cncf) == 0) stop("FACETS returned no copy-number segments for ", sample)

choose_col <- function(dt, candidates, label) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) stop("FACETS output lacks required ", label, " column. Tried: ", paste(candidates, collapse = ", "))
  hit[1]
}

chrom_col <- choose_col(cncf, c("chrom", "Chromosome", "chr"), "chromosome")
start_col <- choose_col(cncf, c("loc.start", "start", "Start"), "segment start")
end_col <- choose_col(cncf, c("loc.end", "end", "End"), "segment end")
tcn_col <- choose_col(cncf, c("tcn.em", "tcn", "tcn.fit"), "total copy number")
lcn_col <- choose_col(cncf, c("lcn.em", "lcn", "lcn.fit"), "minor copy number")
cf_col <- choose_col(cncf, c("cf.em", "cf", "cf.fit"), "cellular fraction")

segments <- data.table(
  sample = sample,
  chrom = as.character(cncf[[chrom_col]]),
  start = as.integer(cncf[[start_col]]),
  end = as.integer(cncf[[end_col]]),
  total_cn = pmax(0L, as.integer(round(as.numeric(cncf[[tcn_col]])))),
  minor_cn = pmax(0L, as.integer(round(as.numeric(cncf[[lcn_col]])))),
  cf = as.numeric(cncf[[cf_col]]),
  purity = purity,
  ploidy = ploidy
)
segments[, major_cn := pmax(0L, total_cn - minor_cn)]
segments[, normal_cn := 2L]
segments <- segments[is.finite(start) & is.finite(end) & end >= start]
if (nrow(segments) == 0) stop("No usable FACETS segments after filtering for ", sample)

fit_path <- file.path(facets_dir, paste0("Auto_", sample, "_facets_fit.rds"))
segments_path <- file.path(tables_dir, paste0("Auto_", sample, "_facets_segments.tsv"))
purity_path <- file.path(tables_dir, paste0("Auto_", sample, "_facets_purity_ploidy.tsv"))
raw_cncf_path <- file.path(tables_dir, paste0("Auto_", sample, "_facets_cncf_raw.tsv"))

saveRDS(
  list(sample = sample, rcmat = rcmat, preprocessed = preprocessed, processed = processed, fit = fit),
  fit_path
)
fwrite(cncf, raw_cncf_path, sep = "\t")
fwrite(segments, segments_path, sep = "\t")
fwrite(data.table(sample = sample, purity = purity, ploidy = ploidy, cval = cval, n_snp_rows = nrow(rcmat)), purity_path, sep = "\t")

writeLines(
  c(
    paste0("sample\t", sample),
    paste0("pileup_file\t", pileup_file),
    paste0("fit_rds\t", fit_path),
    paste0("segments_tsv\t", segments_path),
    paste0("purity_ploidy_tsv\t", purity_path),
    paste0("purity\t", purity),
    paste0("ploidy\t", ploidy),
    paste0("finished\t", as.character(Sys.time()))
  ),
  file.path(logs_dir, paste0("Auto_", sample, "_facets_run_summary.tsv"))
)

message("FACETS complete for ", sample)
