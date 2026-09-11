####################
# Analysis registry (authoritative override):
#   Status: legacy; retained for provenance, no current downstream use
#   Script: analysis/cnv/wes_subclone/legacy_Auto_make_pyclone_input.R
#   Methodology: historical method only; see analysis/ANALYSIS_MAP.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Preserves a superseded implementation or analysis tied to superseded
#     inputs. Do not use its outputs as current centred-MP/state inputs. The
#     original historical inputs, outputs, and method notes remain below.
####################

####################
# Auto_make_pyclone_input.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_make_pyclone_input.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Mutect2 filtered VCF from Sarek, PASS somatic SNVs only
#   - FACETS allele-specific segments from Auto_run_facets_sample.R
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_input.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_mutation_map.tsv
# Downstream use: direct PyClone-VI input. CNVkit total CN is not used.
####################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
sample <- if (length(args) >= 1 && nzchar(args[1])) args[1] else Sys.getenv("sample")
out_root <- if (length(args) >= 2 && nzchar(args[2])) args[2] else Sys.getenv("OUT_ROOT")
mutect_vcf <- if (length(args) >= 3 && nzchar(args[3])) args[3] else Sys.getenv("MUTECT_VCF")

if (!nzchar(sample)) stop("Supply sample as first argument or environment variable sample.")
if (!nzchar(out_root)) {
  out_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
}
if (!nzchar(mutect_vcf) || !file.exists(mutect_vcf)) stop("Missing Mutect2 VCF: ", mutect_vcf)
if (!file.exists(paste0(mutect_vcf, ".tbi")) && !file.exists(paste0(mutect_vcf, ".csi"))) {
  stop("Mutect2 VCF is not indexed: ", mutect_vcf)
}

exclude_sex <- identical(Sys.getenv("EXCLUDE_SEX_CHROMS_DEFAULT", "1"), "1")
min_tumour_depth <- as.integer(Sys.getenv("PYCLONE_MIN_TUMOUR_DEPTH", "20"))
min_alt_count <- as.integer(Sys.getenv("PYCLONE_MIN_ALT_COUNT", "3"))

segments_path <- file.path(out_root, "tables/facets", paste0("Auto_", sample, "_facets_segments.tsv"))
purity_path <- file.path(out_root, "tables/facets", paste0("Auto_", sample, "_facets_purity_ploidy.tsv"))
if (!file.exists(segments_path)) stop("Missing FACETS segments: ", segments_path)
if (!file.exists(purity_path)) stop("Missing FACETS purity/ploidy table: ", purity_path)

segments <- fread(segments_path)
purity_tbl <- fread(purity_path)
if (!all(c("chrom", "start", "end", "major_cn", "minor_cn", "normal_cn", "purity") %in% names(segments))) {
  stop("FACETS segments table lacks required PyClone columns: ", segments_path)
}
tumour_content <- as.numeric(purity_tbl$purity[1])
if (!is.finite(tumour_content) || tumour_content <= 0 || tumour_content > 1) {
  stop("Invalid FACETS tumour_content/purity for ", sample, ": ", tumour_content)
}

message("Extracting Mutect2 AD fields with bcftools: ", mutect_vcf)
sample_names <- system2("bcftools", c("query", "-l", mutect_vcf), stdout = TRUE)
if (length(sample_names) != 2) stop("Expected exactly two VCF samples in ", mutect_vcf, "; found: ", paste(sample_names, collapse = ", "))

query_cmd <- paste(
  "bcftools query -f",
  shQuote("%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n"),
  shQuote(mutect_vcf)
)
vcf_dt <- fread(cmd = query_cmd, header = FALSE, sep = "\t", quote = "")
if (ncol(vcf_dt) != 7) stop("Unexpected bcftools query column count for ", mutect_vcf, ": ", ncol(vcf_dt))
setnames(vcf_dt, c("chrom", "pos", "ref", "alt", "filter", "normal_ad", "tumour_ad"))

parse_ad <- function(x, which) {
  parts <- tstrsplit(x, ",", fixed = TRUE)
  suppressWarnings(as.integer(parts[[which]]))
}

vcf_dt[, tumour_ref_count := parse_ad(tumour_ad, 1)]
vcf_dt[, tumour_alt_count := parse_ad(tumour_ad, 2)]
vcf_dt[, normal_ref_count := parse_ad(normal_ad, 1)]
vcf_dt[, normal_alt_count := parse_ad(normal_ad, 2)]
vcf_dt[, chrom_clean := sub("^chr", "", as.character(chrom))]

autosomes <- as.character(1:22)
snvs <- vcf_dt[
  filter == "PASS" &
    nchar(ref) == 1 &
    nchar(alt) == 1 &
    !grepl(",", alt, fixed = TRUE) &
    is.finite(tumour_ref_count) &
    is.finite(tumour_alt_count)
]
if (exclude_sex) snvs <- snvs[chrom_clean %in% autosomes]
snvs <- snvs[(tumour_ref_count + tumour_alt_count) >= min_tumour_depth & tumour_alt_count >= min_alt_count]
if (nrow(snvs) == 0) stop("No PASS somatic SNVs remain after PyClone filters for ", sample)

segments[, chrom_clean := sub("^chr", "", as.character(chrom))]
if (exclude_sex) segments <- segments[chrom_clean %in% autosomes]
segments <- segments[is.finite(start) & is.finite(end)]
segments[, `:=`(start = as.integer(start), end = as.integer(end))]

seg_ov <- segments[, .(
  chrom_clean,
  seg_start = start,
  seg_end = end,
  major_cn = as.integer(major_cn),
  minor_cn = as.integer(minor_cn),
  normal_cn = as.integer(normal_cn),
  tumour_content = tumour_content
)]
setnames(seg_ov, c("chrom_clean", "seg_start", "seg_end"), c("chrom_clean", "start", "end"))

var_ov <- snvs[, .(
  chrom_clean,
  start = as.integer(pos),
  end = as.integer(pos),
  chrom,
  pos = as.integer(pos),
  ref,
  alt,
  tumour_ref_count,
  tumour_alt_count,
  normal_ref_count,
  normal_alt_count
)]

setkey(seg_ov, chrom_clean, start, end)
setkey(var_ov, chrom_clean, start, end)
mapped <- foverlaps(var_ov, seg_ov, nomatch = 0L)
if (nrow(mapped) == 0) stop("No PASS somatic SNVs overlap FACETS segments for ", sample)

mapped[, mutation_id := paste(chrom, pos, ref, alt, sep = ":")]
pyclone_input <- mapped[, .(
  mutation_id,
  sample_id = sample,
  ref_counts = as.integer(tumour_ref_count),
  alt_counts = as.integer(tumour_alt_count),
  major_cn = pmax(0L, as.integer(major_cn)),
  minor_cn = pmax(0L, as.integer(minor_cn)),
  normal_cn = pmax(0L, as.integer(normal_cn)),
  tumour_content = tumour_content
)]
pyclone_input <- unique(pyclone_input, by = "mutation_id")
pyclone_input <- pyclone_input[major_cn >= minor_cn & (major_cn + minor_cn) > 0]
if (nrow(pyclone_input) < 50) {
  warning("Only ", nrow(pyclone_input), " variants remain for PyClone-VI in ", sample, "; inference may be unstable.")
}

tables_dir <- file.path(out_root, "tables/pyclone")
logs_dir <- file.path(out_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(tables_dir, paste0("Auto_", sample, "_pyclone_vi_input.tsv"))
map_path <- file.path(tables_dir, paste0("Auto_", sample, "_pyclone_vi_mutation_map.tsv"))
summary_path <- file.path(logs_dir, paste0("Auto_", sample, "_pyclone_input_summary.tsv"))

fwrite(pyclone_input, input_path, sep = "\t")
fwrite(mapped[, .(
  mutation_id, chrom, pos, ref, alt,
  tumour_ref_count, tumour_alt_count, normal_ref_count, normal_alt_count,
  segment_start = i.start, segment_end = i.end,
  major_cn, minor_cn, normal_cn, tumour_content
)], map_path, sep = "\t")
fwrite(data.table(
  sample = sample,
  mutect_vcf = mutect_vcf,
  n_pass_snv_after_depth_filter = nrow(snvs),
  n_pyclone_variants = nrow(pyclone_input),
  exclude_sex_chromosomes = exclude_sex,
  min_tumour_depth = min_tumour_depth,
  min_alt_count = min_alt_count,
  tumour_content = tumour_content,
  finished = as.character(Sys.time())
), summary_path, sep = "\t")

message("Wrote PyClone-VI input: ", input_path)
