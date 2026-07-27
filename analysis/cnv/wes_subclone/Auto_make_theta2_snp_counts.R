####################
# Auto_make_theta2_snp_counts.R
#
# Analysis registry
# Status: active helper
# Script: analysis/cnv/wes_subclone/Auto_make_theta2_snp_counts.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - FACETS snp-pileup CSV.gz for a matched normal/tumour WES pair
# Outputs:
#   - THetA2 normal.snp_formatted.txt and tumor.snp_formatted.txt files
#   - lightweight SNP-count summary TSV
# Downstream use: helper for THetA2 clone-specific WES CNA deconvolution.
####################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript Auto_make_theta2_snp_counts.R <sample> <pileup.csv.gz> <out_dir> <summary_tsv>")
}

sample_id <- args[[1]]
pileup_path <- args[[2]]
out_dir <- args[[3]]
summary_path <- args[[4]]

if (!file.exists(pileup_path)) stop("Missing FACETS snp-pileup: ", pileup_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)

dt <- fread(pileup_path)
required <- c("Chromosome", "Position", "Ref", "Alt", "File1R", "File1A", "File2R", "File2A")
missing <- setdiff(required, names(dt))
if (length(missing) > 0) stop("FACETS pileup missing columns: ", paste(missing, collapse = ", "))

autosomes <- paste0("chr", 1:22)
dt <- dt[Chromosome %in% autosomes]
dt <- dt[Ref %chin% c("A", "C", "G", "T") & Alt %chin% c("A", "C", "G", "T")]
dt[, `:=`(
  normal_ref = as.integer(File1R),
  normal_alt = as.integer(File1A),
  tumor_ref = as.integer(File2R),
  tumor_alt = as.integer(File2A)
)]
dt[, `:=`(
  normal_depth = normal_ref + normal_alt,
  tumor_depth = tumor_ref + tumor_alt
)]
dt <- dt[normal_depth >= 20 & tumor_depth >= 20]
dt[, normal_alt_frac := normal_alt / pmax(1, normal_depth)]

# THetA2 BAF input should represent likely germline heterozygous SNPs. The
# matched normal is used for this filter; somatic Mutect2 SNVs are not used.
het <- dt[normal_alt_frac >= 0.25 & normal_alt_frac <= 0.75]
het[, chrom_theta := sub("^chr", "", Chromosome)]

normal_out <- het[, .(
  `#Chrm` = chrom_theta,
  Pos = as.integer(Position),
  Ref_Allele = normal_ref,
  Mut_Allele = normal_alt
)]
tumor_out <- het[, .(
  `#Chrm` = chrom_theta,
  Pos = as.integer(Position),
  Ref_Allele = tumor_ref,
  Mut_Allele = tumor_alt
)]

normal_path <- file.path(out_dir, paste0("Auto_", sample_id, ".normal.snp_formatted.txt"))
tumor_path <- file.path(out_dir, paste0("Auto_", sample_id, ".tumor.snp_formatted.txt"))
fwrite(normal_out, normal_path, sep = "\t")
fwrite(tumor_out, tumor_path, sep = "\t")

summary <- data.table(
  sample = sample_id,
  pileup_path = pileup_path,
  n_raw_rows = nrow(dt),
  n_likely_germline_het_snps = nrow(het),
  normal_snp_file = normal_path,
  tumor_snp_file = tumor_path,
  normal_depth_min = if (nrow(het) > 0) min(het$normal_depth) else NA_integer_,
  tumor_depth_min = if (nrow(het) > 0) min(het$tumor_depth) else NA_integer_,
  normal_alt_frac_median = if (nrow(het) > 0) median(het$normal_alt_frac) else NA_real_,
  finished = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
fwrite(summary, summary_path, sep = "\t")

if (nrow(het) < 100) {
  stop("Too few likely germline heterozygous SNPs for THetA2 BAF input: ", nrow(het))
}
