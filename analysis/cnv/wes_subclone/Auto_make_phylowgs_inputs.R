####################
# Auto_make_phylowgs_inputs.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_make_phylowgs_inputs.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - FACETS allele-specific segments and purity/ploidy tables.
#   - FACETS-aware PyClone-VI input/results tables.
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/ssm_data.txt
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/cnv_data.txt
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_*_audit.tsv
# Downstream use: direct PhyloWGS input and provenance tables.
####################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
sample <- if (length(args) >= 1 && nzchar(args[1])) args[1] else Sys.getenv("sample")
live_root <- if (length(args) >= 2 && nzchar(args[2])) args[2] else Sys.getenv("OUT_ROOT")
ephemeral_root <- if (length(args) >= 3 && nzchar(args[3])) args[3] else Sys.getenv("EPHEMERAL_OUT_ROOT")

if (!nzchar(sample)) stop("Supply sample as first argument or environment variable sample.")
if (!nzchar(live_root)) {
  live_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
}
if (!nzchar(ephemeral_root)) {
  ephemeral_root <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
}

find_existing <- function(paths, label) {
  hits <- paths[file.exists(paths)]
  if (!length(hits)) {
    stop("Missing ", label, ". Checked: ", paste(paths, collapse = " | "))
  }
  hits[[1]]
}

copy_to_live_if_missing <- function(src, dest) {
  if (!file.exists(src)) return(FALSE)
  if (file.exists(dest)) return(FALSE)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dest, overwrite = FALSE)
  if (!ok) stop("Failed to copy critical input from ", src, " to ", dest)
  TRUE
}

live_pyclone_dir <- file.path(live_root, "tables/pyclone")
eph_pyclone_dir <- file.path(ephemeral_root, "tables/pyclone")
critical_pyclone_files <- c(
  paste0("Auto_", sample, "_pyclone_vi_input.tsv"),
  paste0("Auto_", sample, "_pyclone_vi_results.tsv"),
  paste0("Auto_", sample, "_pyclone_vi_mutation_map.tsv")
)
copied <- vapply(
  critical_pyclone_files,
  function(fname) copy_to_live_if_missing(file.path(eph_pyclone_dir, fname), file.path(live_pyclone_dir, fname)),
  logical(1)
)

segments_path <- find_existing(c(
  file.path(live_root, "tables/facets", paste0("Auto_", sample, "_facets_segments.tsv")),
  file.path(ephemeral_root, "tables/facets", paste0("Auto_", sample, "_facets_segments.tsv"))
), "FACETS segments")
purity_path <- find_existing(c(
  file.path(live_root, "tables/facets", paste0("Auto_", sample, "_facets_purity_ploidy.tsv")),
  file.path(ephemeral_root, "tables/facets", paste0("Auto_", sample, "_facets_purity_ploidy.tsv"))
), "FACETS purity/ploidy")
pyclone_input_path <- find_existing(c(
  file.path(live_pyclone_dir, paste0("Auto_", sample, "_pyclone_vi_input.tsv")),
  file.path(eph_pyclone_dir, paste0("Auto_", sample, "_pyclone_vi_input.tsv"))
), "PyClone-VI input")
pyclone_results_path <- find_existing(c(
  file.path(live_pyclone_dir, paste0("Auto_", sample, "_pyclone_vi_results.tsv")),
  file.path(eph_pyclone_dir, paste0("Auto_", sample, "_pyclone_vi_results.tsv"))
), "PyClone-VI results")

segments <- fread(segments_path, na.strings = c("", "NA", "NaN"))
purity_tbl <- fread(purity_path)
pyclone_input <- fread(pyclone_input_path)
pyclone_results <- fread(pyclone_results_path)

required_segments <- c("chrom", "start", "end", "total_cn", "minor_cn", "cf", "purity", "ploidy", "major_cn", "normal_cn")
missing_segments <- setdiff(required_segments, names(segments))
if (length(missing_segments)) stop("FACETS segments missing columns: ", paste(missing_segments, collapse = ", "))
required_pyclone <- c("mutation_id", "ref_counts", "alt_counts")
missing_pyclone <- setdiff(required_pyclone, names(pyclone_input))
if (length(missing_pyclone)) stop("PyClone input missing columns: ", paste(missing_pyclone, collapse = ", "))

purity <- as.numeric(purity_tbl$purity[[1]])
ploidy <- as.numeric(purity_tbl$ploidy[[1]])
if (!is.finite(purity) || purity <= 0 || purity > 1) stop("Invalid FACETS purity for ", sample)
if (!is.finite(ploidy) || ploidy <= 0) stop("Invalid FACETS ploidy for ", sample)

clean_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
num <- function(x) suppressWarnings(as.numeric(x))

segments[, `:=`(
  chrom_clean = clean_chr(chrom),
  start = as.integer(start),
  end = as.integer(end),
  total_cn_num = num(total_cn),
  major_cn_num = num(major_cn),
  minor_cn_num = num(minor_cn),
  normal_cn_num = num(normal_cn),
  cf_num = num(cf)
)]
segments[!is.finite(normal_cn_num), normal_cn_num := 2]
segments[!is.finite(cf_num), cf_num := purity]
segments[, cf_num := pmin(pmax(cf_num, 0.001), purity)]

segments[, allele_cn_imputed := !is.finite(major_cn_num) | !is.finite(minor_cn_num)]
segments[!is.finite(major_cn_num) & is.finite(minor_cn_num) & is.finite(total_cn_num),
         major_cn_num := total_cn_num - minor_cn_num]
segments[!is.finite(minor_cn_num) & is.finite(major_cn_num) & is.finite(total_cn_num),
         minor_cn_num := total_cn_num - major_cn_num]
segments[!is.finite(major_cn_num) & !is.finite(minor_cn_num) & is.finite(total_cn_num),
         `:=`(major_cn_num = ceiling(total_cn_num / 2), minor_cn_num = floor(total_cn_num / 2))]
segments <- segments[
  chrom_clean %in% as.character(1:22) &
    is.finite(start) & is.finite(end) & end > start &
    is.finite(major_cn_num) & is.finite(minor_cn_num)
]
segments[, `:=`(
  major_cn_int = pmax(0L, as.integer(round(major_cn_num))),
  minor_cn_int = pmax(0L, as.integer(round(minor_cn_num))),
  total_cn_int = pmax(0L, as.integer(round(major_cn_num + minor_cn_num))),
  normal_cn_int = pmax(1L, as.integer(round(normal_cn_num)))
)]
swap_idx <- which(segments$major_cn_int < segments$minor_cn_int)
if (length(swap_idx)) {
  old_major <- segments$major_cn_int[swap_idx]
  segments$major_cn_int[swap_idx] <- segments$minor_cn_int[swap_idx]
  segments$minor_cn_int[swap_idx] <- old_major
}

split_id <- tstrsplit(pyclone_input$mutation_id, ":", fixed = TRUE)
if (length(split_id) < 4) stop("Mutation IDs are not chr:pos:ref:alt formatted.")
pyclone_input[, `:=`(
  chrom = split_id[[1]],
  pos = as.integer(split_id[[2]]),
  ref = split_id[[3]],
  alt = split_id[[4]],
  chrom_clean = clean_chr(split_id[[1]]),
  ref_counts = as.integer(ref_counts),
  alt_counts = as.integer(alt_counts)
)]
ssm <- unique(pyclone_input[
  chrom_clean %in% as.character(1:22) &
    is.finite(pos) &
    is.finite(ref_counts) &
    is.finite(alt_counts) &
    (ref_counts + alt_counts) > 0
], by = "mutation_id")
ssm[, chrom_order := suppressWarnings(as.integer(chrom_clean))]
setorder(ssm, chrom_order, pos, mutation_id)
ssm[, `:=`(
  ssm_id = paste0("s", seq_len(.N) - 1L),
  gene = paste(chrom_clean, pos, sep = "_"),
  total_counts = ref_counts + alt_counts,
  mu_r = 0.999,
  mu_v = 0.499
)]

if (nrow(ssm) < 50) {
  warning("Only ", nrow(ssm), " SSMs available for PhyloWGS in ", sample, "; tree inference may be unstable.")
}

ssm_data <- ssm[, .(
  id = ssm_id,
  gene,
  a = ref_counts,
  d = total_counts,
  mu_r,
  mu_v
)]

is_neutral <- segments$total_cn_int == segments$normal_cn_int &
  segments$major_cn_int == 1L &
  segments$minor_cn_int == 1L
cnv_segments <- copy(segments[!is_neutral])

median_depth <- median(ssm$total_counts, na.rm = TRUE)
if (!is.finite(median_depth) || median_depth <= 0) median_depth <- 60
het_snp_rate <- as.numeric(Sys.getenv("PHYLOWGS_HETSNP_RATE", "0.0007"))
max_equiv_ssms <- as.numeric(Sys.getenv("PHYLOWGS_MAX_EQUIV_SSMS", "3000"))
if (!is.finite(het_snp_rate) || het_snp_rate <= 0) het_snp_rate <- 0.0007
if (!is.finite(max_equiv_ssms) || max_equiv_ssms <= 0) max_equiv_ssms <- 3000

format_overlapping_ssms <- function(chrom_value, seg_start, seg_end, minor_cn_value, major_cn_value) {
  hits <- ssm[chrom_clean == chrom_value & pos >= seg_start & pos <= seg_end, ssm_id]
  if (!length(hits)) return("")
  paste(paste(hits, minor_cn_value, major_cn_value, sep = ","), collapse = ";")
}

if (nrow(cnv_segments)) {
  cnv_segments[, `:=`(
    cnv = paste0("c", seq_len(.N) - 1L),
    pseudo_total_reads = as.integer(round(pmax(
      1,
      pmin((end - start + 1) * het_snp_rate * median_depth, max_equiv_ssms * median_depth)
    )))
  )]
  cnv_segments[, pseudo_ref_reads := as.integer(round((1 - (cf_num / 2)) * pseudo_total_reads))]
  cnv_segments[, ssms := mapply(
    format_overlapping_ssms,
    chrom_clean,
    start,
    end,
    minor_cn_int,
    major_cn_int,
    USE.NAMES = FALSE
  )]
  cnv_segments[, physical_cnvs := paste0(
    "chrom=", chrom_clean,
    ",start=", start,
    ",end=", end,
    ",major_cn=", major_cn_int,
    ",minor_cn=", minor_cn_int,
    ",cell_prev=", signif(cf_num, 6)
  )]
  cnv_data <- cnv_segments[, .(
    cnv,
    a = pseudo_ref_reads,
    d = pseudo_total_reads,
    ssms,
    physical_cnvs
  )]
} else {
  cnv_data <- data.table(cnv = character(), a = integer(), d = integer(), ssms = character(), physical_cnvs = character())
}

pyclone_map <- merge(
  ssm[, .(mutation_id, ssm_id, chrom, pos, ref, alt, ref_counts, alt_counts, total_counts)],
  pyclone_results,
  by = "mutation_id",
  all.x = TRUE
)
setorder(pyclone_map, ssm_id)

out_dir <- file.path(live_root, "tables/phylowgs", sample)
log_dir <- file.path(live_root, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

ssm_path <- file.path(out_dir, "ssm_data.txt")
cnv_path <- file.path(out_dir, "cnv_data.txt")
ssm_audit_path <- file.path(out_dir, paste0("Auto_", sample, "_phylowgs_ssm_pyclone_map.tsv"))
cnv_audit_path <- file.path(out_dir, paste0("Auto_", sample, "_phylowgs_cnv_segment_audit.tsv"))
summary_path <- file.path(out_dir, paste0("Auto_", sample, "_phylowgs_input_summary.tsv"))
manifest_path <- file.path(log_dir, paste0("Auto_", sample, "_phylowgs_input_manifest.tsv"))

fwrite(ssm_data, ssm_path, sep = "\t", quote = FALSE)
fwrite(cnv_data, cnv_path, sep = "\t", quote = FALSE)
fwrite(pyclone_map, ssm_audit_path, sep = "\t", quote = FALSE)
fwrite(cnv_segments[, .(
  cnv,
  chrom = chrom_clean,
  start,
  end,
  major_cn = major_cn_int,
  minor_cn = minor_cn_int,
  total_cn = total_cn_int,
  normal_cn = normal_cn_int,
  cell_prev = cf_num,
  purity,
  ploidy,
  allele_cn_imputed,
  pseudo_ref_reads,
  pseudo_total_reads,
  n_overlapping_ssms = fifelse(nchar(ssms) == 0, 0L, lengths(strsplit(ssms, ";", fixed = TRUE)))
)], cnv_audit_path, sep = "\t", quote = FALSE)

cluster_count <- if ("cluster_id" %in% names(pyclone_results)) uniqueN(pyclone_results$cluster_id) else NA_integer_
summary <- data.table(
  sample = sample,
  facets_segments = segments_path,
  purity_ploidy = purity_path,
  pyclone_input = pyclone_input_path,
  pyclone_results = pyclone_results_path,
  copied_critical_pyclone_files_to_live = paste(names(copied)[copied], collapse = ";"),
  purity = purity,
  ploidy = ploidy,
  n_ssms = nrow(ssm_data),
  n_pyclone_clusters = cluster_count,
  n_facets_segments = nrow(segments),
  n_cnv_events = nrow(cnv_data),
  n_cnv_events_with_ssms = if (nrow(cnv_data)) sum(nchar(cnv_data$ssms) > 0) else 0L,
  n_cnv_events_with_imputed_allele_cn = if (nrow(cnv_segments)) sum(cnv_segments$allele_cn_imputed) else 0L,
  median_ssm_depth = median_depth,
  het_snp_rate = het_snp_rate,
  phylowgs_ssm_data = ssm_path,
  phylowgs_cnv_data = cnv_path,
  finished = as.character(Sys.time())
)
fwrite(summary, summary_path, sep = "\t", quote = FALSE)
manifest <- data.table(
  field = names(summary),
  value = vapply(summary, function(x) paste(as.character(x), collapse = ";"), character(1))
)
fwrite(manifest, manifest_path, sep = "\t", quote = FALSE)

message("Wrote PhyloWGS inputs for ", sample, " to ", out_dir)
