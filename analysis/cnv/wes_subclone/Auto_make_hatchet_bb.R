####################
# Auto_make_hatchet_bb.R
#
# Analysis registry
# Status: active helper
# Script: analysis/cnv/wes_subclone/Auto_make_hatchet_bb.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Sarek CNVkit tumour .cnr target/antitarget bin table
#   - FACETS snp-pileup from matched normal/tumour WES pair
# Outputs:
#   - live/ephemeral HATCHet .bb input table with RDR/BAF fields
#   - live HATCHet .bb preparation summary TSV
# Downstream use: HATCHet cluster-bins and compute-cn.
####################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript Auto_make_hatchet_bb.R <wgs_id> <out_bb> <summary_tsv> <bin_size_bp>")
}

wgs_id <- args[[1]]
out_bb <- args[[2]]
summary_tsv <- args[[3]]
bin_size <- as.integer(args[[4]])
if (!is.finite(bin_size) || bin_size <= 0) stop("Invalid bin size: ", args[[4]])

sarek_root <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect"
live_project <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
ephemeral_project <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"

pair_for_sample <- function(sample_id) {
  switch(
    sample_id,
    "PDO_1090_vs_NT_1090" = list(tumor = "PDO_1090", normal = "NT_1090", hatchet_sample = "PDO_1090"),
    "PDO_1181_vs_NT_1181" = list(tumor = "PDO_1181", normal = "NT_1181", hatchet_sample = "PDO_1181"),
    stop("Unsupported WES pair for HATCHet BB preparation: ", sample_id)
  )
}

pair <- pair_for_sample(wgs_id)
cnr_path <- file.path(sarek_root, "variant_calling/cnvkit", wgs_id, paste0(pair$tumor, ".cnr"))
shift_summary_path <- file.path(
  live_project,
  "PDOs_outs/Auto_wes_absolute_cna/tables/Auto_wes_absolute_cna_profile_summary.csv"
)
pileup_path <- file.path(
  ephemeral_project,
  "PDOs_outs/Auto_wes_subclone/intermediate/facets",
  wgs_id,
  paste0(wgs_id, ".snp_pileup.gz")
)

if (!file.exists(cnr_path)) stop("Missing CNVkit CNR: ", cnr_path)
if (!file.exists(pileup_path)) stop("Missing FACETS snp-pileup: ", pileup_path)
dir.create(dirname(out_bb), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)

applied_shift_log2 <- 0
shift_source <- "none"
if (file.exists(shift_summary_path)) {
  shift_dt <- fread(shift_summary_path)
  target_wgs_id <- wgs_id
  shift_row <- shift_dt[wgs_id == target_wgs_id & profile == "bulk_highres"]
  if (nrow(shift_row) > 0 && is.finite(as.numeric(shift_row$applied_shift_log2[1]))) {
    applied_shift_log2 <- as.numeric(shift_row$applied_shift_log2[1])
    shift_source <- shift_summary_path
  }
}

autosomes_chr <- paste0("chr", 1:22)
autosomes <- as.character(1:22)

cnr <- fread(cnr_path)
required_cnr <- c("chromosome", "start", "end", "depth", "log2", "weight")
missing_cnr <- setdiff(required_cnr, names(cnr))
if (length(missing_cnr) > 0) stop("CNVkit CNR missing columns: ", paste(missing_cnr, collapse = ", "))
cnr <- cnr[chromosome %chin% autosomes_chr]
cnr[, `:=`(
  CHR = sub("^chr", "", chromosome),
  window_start = floor(as.numeric(start) / bin_size) * bin_size,
  window_end = floor(as.numeric(start) / bin_size) * bin_size + bin_size,
  rd_raw = 2 ^ (as.numeric(log2) + applied_shift_log2),
  depth_num = as.numeric(depth),
  weight_num = as.numeric(weight)
)]
cnr <- cnr[CHR %chin% autosomes & is.finite(rd_raw) & rd_raw > 0 & is.finite(depth_num)]
cnr[!is.finite(weight_num) | weight_num <= 0, weight_num := 1]

rd_bins <- cnr[, .(
  START = min(window_start),
  END = max(window_end),
  RD = weighted.mean(rd_raw, weight_num, na.rm = TRUE),
  COV = weighted.mean(depth_num, weight_num, na.rm = TRUE),
  N_CNVKIT_BINS = .N
), by = .(CHR, window_start)]

pileup <- fread(pileup_path)
required_pileup <- c("Chromosome", "Position", "Ref", "Alt", "File1R", "File1A", "File2R", "File2A")
missing_pileup <- setdiff(required_pileup, names(pileup))
if (length(missing_pileup) > 0) stop("FACETS pileup missing columns: ", paste(missing_pileup, collapse = ", "))

pileup <- pileup[Chromosome %chin% autosomes_chr]
pileup <- pileup[Ref %chin% c("A", "C", "G", "T") & Alt %chin% c("A", "C", "G", "T")]
pileup[, `:=`(
  normal_ref = as.integer(File1R),
  normal_alt = as.integer(File1A),
  tumor_ref = as.integer(File2R),
  tumor_alt = as.integer(File2A)
)]
pileup[, `:=`(
  normal_depth = normal_ref + normal_alt,
  tumor_depth = tumor_ref + tumor_alt
)]
pileup <- pileup[normal_depth >= 20 & tumor_depth >= 20]
pileup[, normal_alt_frac := normal_alt / pmax(1, normal_depth)]
het <- pileup[normal_alt_frac >= 0.25 & normal_alt_frac <= 0.75]
het[, `:=`(
  CHR = sub("^chr", "", Chromosome),
  window_start = floor(as.numeric(Position) / bin_size) * bin_size,
  allele_a = pmax(tumor_ref, tumor_alt),
  allele_b = pmin(tumor_ref, tumor_alt)
)]

baf_bins <- het[, .(
  `#SNPS` = .N,
  ALPHA = sum(allele_a, na.rm = TRUE),
  BETA = sum(allele_b, na.rm = TRUE)
), by = .(CHR, window_start)]
baf_bins[, BAF := BETA / pmax(1, ALPHA + BETA)]

bb <- merge(rd_bins, baf_bins, by = c("CHR", "window_start"), all = FALSE)
bb <- bb[`#SNPS` > 0 & is.finite(RD) & RD > 0 & is.finite(BAF)]
bb[, SAMPLE := pair$hatchet_sample]
bb <- bb[order(as.integer(CHR), START)]
bb_out <- bb[, .(
  `#CHR` = CHR,
  START = as.integer(START),
  END = as.integer(END),
  SAMPLE,
  RD,
  `#SNPS`,
  COV,
  ALPHA,
  BETA,
  BAF
)]

fwrite(bb_out, out_bb, sep = "\t", quote = FALSE)

summary <- data.table(
  wgs_id = wgs_id,
  tumor = pair$tumor,
  normal = pair$normal,
  hatchet_sample = pair$hatchet_sample,
  bin_size_bp = bin_size,
  cnr_path = cnr_path,
  pileup_path = pileup_path,
  n_cnr_rows = nrow(cnr),
  n_het_snps = nrow(het),
  n_hatchet_bins = nrow(bb_out),
  median_rd = median(bb_out$RD, na.rm = TRUE),
  median_baf = median(bb_out$BAF, na.rm = TRUE),
  applied_shift_log2 = applied_shift_log2,
  rd_scale_factor = 2 ^ applied_shift_log2,
  shift_source = shift_source,
  output_bb = out_bb,
  finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
fwrite(summary, summary_tsv, sep = "\t")

if (nrow(bb_out) < 100) {
  stop("Too few HATCHet bins with both RDR and BAF support: ", nrow(bb_out))
}
