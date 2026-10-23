####################
# Current analysis registry (takes precedence over the historical block below):
#   Status: active terminal visualization and MUC16 comparison
#   Script: analysis/cnv/Auto_indel_visualisation.R
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Restore live SNP/indel metrics, export exact-allele MUC16 status
#     in both PDOs and PASS-SNP pileup support, then regenerate historical plots.
#   Inputs: Sarek Mutect2 annotated VCFs and PDO CRAMs for 1072/1090; genome.fa;
#     VCFs/CRAMs under /rds/general/project/spatialtranscriptomics/live/sarek_mutect/;
#     exact filenames are declared below.
#   Outputs: PDOs_outs/indel_analysis/ (live): metrics TSVs, Auto_PDO_* VCFs,
#     Auto_variant_prioritization_details.csv, Auto_indel_chromatogram_style.pdf,
#     Auto_indel_sequence_alignments.png, Auto_variant_clusters_plot.pdf;
#     tables/Auto_indel_depth_profiles.csv; intermediate/Auto_plot_sources.rds;
#     tables/Auto_MUC16_all_variant_status.csv and
#     tables/Auto_MUC16_PASS_SNP_read_support.csv; preparation logs/report.
#     tables/Auto_variant_plot_selection.csv; intermediate/Auto_variant_plot_data.rds.
#   Plot scope: first 18 ranked variants, grouped into gene panels; includes
#     MUC16 at rank 17 and TP53 at rank 18 in the current data.
#   Downstream: terminal figures and review tables; no state-analysis dependency.
#   Cache/replot: complete deterministic rebuild, no cache/fallback.
#   Run: Rscript analysis/cnv/Auto_indel_visualisation.R
#   PBS: qsub job_scripts/cnv/Auto_indel_visualisation.sh from project root;
#     4 CPUs, 32 GB, 02:00:00; dmtcp, SAMtools 1.22.1, BCFtools 1.22.
####################

####################
# Analysis registry:
#   Status: active visualization
#   Script: analysis/cnv/Auto_indel_visualisation.R
#   Methodology: None (Visualization and Depth verification script)
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     PDOs_outs/indel_analysis/PDO_1072_variants_metrics.tsv
#     PDOs_outs/indel_analysis/PDO_1090_variants_metrics.tsv
#     Sarek recalibrated CRAMs
#   Outputs:
#     PDOs_outs/indel_analysis/Auto_indel_chromatogram_style.pdf
#     PDOs_outs/indel_analysis/Auto_indel_sequence_alignments.png
#     PDOs_outs/indel_analysis/Auto_variant_prioritization_table.csv
#     PDOs_outs/indel_analysis/Auto_variant_clusters_plot.png
#   Downstream: Terminal figures
####################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

WD <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(WD)

OUT_DIR <- file.path(WD, "PDOs_outs/indel_analysis")
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive=TRUE)

# 1. Load metrics TSVs
cols <- c("chrom", "pos", "ref", "alt", "filter", "tlod", "dp", "germq", "nalod", "nlod", "mbq", "mpos", "csq", "ad_normal", "af_normal", "ad_tumor", "af_tumor")

read_variant_metrics <- function(filepath, sample_name) {
  if (!file.exists(filepath)) return(data.frame())
  
  df <- read.csv(filepath, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  if(nrow(df) == 0) return(data.frame())
  if(ncol(df) > 17) df <- df[, 1:17]
  if(ncol(df) < 17) return(data.frame()) # Incomplete data
  
  colnames(df) <- cols
  
  df$sample <- sample_name
  df$af_tumor <- as.numeric(gsub("]", "", df$af_tumor))
  df$af_normal <- as.numeric(gsub("\\[", "", df$af_normal))
  df$dp <- as.numeric(df$dp)
  df$tlod <- as.numeric(df$tlod)
  
  df$indel_len <- abs(nchar(df$ref) - nchar(df$alt))
  df$is_snp <- df$indel_len == 0 & nchar(df$ref) == 1 & nchar(df$alt) == 1
  
  # Parse CSQ
  first_csq <- sapply(strsplit(as.character(df$csq), ","), `[`, 1)
  csq_parts <- strsplit(first_csq, "\\|")
  
  df$consequence <- sapply(csq_parts, function(x) if(length(x) >= 2) x[2] else NA)
  df$gene <- sapply(csq_parts, function(x) if(length(x) >= 4) x[4] else NA)
  
  df$is_exon <- grepl("missense|synonymous|stop|splice|frameshift|inframe|exon|coding", df$consequence, ignore.case=TRUE)
  
  return(df)
}

####################
source("analysis/shared/Auto_pdo_analysis_config.R")
source("analysis/shared/Auto_pdo_analysis_helpers.R")
out_dir <- file.path(PDO_OUTPUT_DIR, "indel_analysis")
for (tier in c("", "intermediate", "tables", "logs", "reports")) {
  dir.create(file.path(out_dir, tier), recursive = TRUE, showWarnings = FALSE)
}
command_log <- file.path(out_dir, "logs/Auto_metrics_commands.log")
writeLines(paste("Started", Sys.time()), command_log)
run_checked <- function(tool, args, stdout = "") {
  if (!nzchar(Sys.which(tool))) stop("Required executable missing: ", tool)
  cat(paste(Sys.time(), tool, paste(args, collapse = " "), "\n"),
      file = command_log, append = TRUE)
  stderr_path <- file.path(out_dir, "logs/Auto_last_command_stderr.log")
  status <- system2(tool, args, stdout = stdout, stderr = stderr_path)
  cat(readLines(stderr_path, warn = FALSE), sep = "\n", file = command_log, append = TRUE)
  if (!identical(as.integer(status), 0L)) stop(tool, " failed: ", status)
}
sarek <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect"
reference <- "/rds/general/project/tumourheterogeneity1/live/demultiplex/genome.fa"
pdo_require_files(c(reference, paste0(reference, ".fai")))
samples <- c("1072", "1090")
vcf_paths <- setNames(file.path(sarek, "annotation/mutect2",
  paste0("PDO_", samples, "_vs_NT_", samples),
  paste0("PDO_", samples, "_vs_NT_", samples, ".mutect2.filtered_VEP.ann.vcf.gz")), samples)
cram_paths <- setNames(file.path(sarek, "preprocessing/recalibrated",
  paste0("PDO_", samples), paste0("PDO_", samples, ".recal.cram")), samples)
pdo_require_files(c(vcf_paths, cram_paths, paste0(cram_paths, ".crai")))
records <- list()
metric_format <- paste0("%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TLOD\t%INFO/DP\t",
  "%INFO/GERMQ\t%INFO/NALOD\t%INFO/NLOD\t%INFO/MBQ\t%INFO/MPOS\t%INFO/CSQ\t[%AD\t%AF\t]\n")
for (sample_id in samples) {
  sample_file <- file.path(out_dir, "intermediate", paste0("Auto_samples_", sample_id, ".tsv"))
  run_checked("bcftools", c("query", "-l", shQuote(vcf_paths[[sample_id]])), sample_file)
  expected <- c(paste0("SUR", sample_id, "_NT_", sample_id),
                paste0("SUR", sample_id, "_PDO_", sample_id))
  if (!identical(readLines(sample_file), expected)) stop("Unexpected VCF sample order")
  output_vcf <- file.path(out_dir, paste0("Auto_PDO_", sample_id, "_variants.vcf.gz"))
  run_checked("bcftools", c("view", "-v", "snps,indels", shQuote(vcf_paths[[sample_id]]),
                           "-Oz", "-o", shQuote(output_vcf)))
  run_checked("bcftools", c("index", "-f", shQuote(output_vcf)))
  metrics <- file.path(out_dir, paste0("PDO_", sample_id, "_variants_metrics.tsv"))
  run_checked("bcftools", c("query", "-f", shQuote(metric_format), shQuote(output_vcf)), metrics)
  if (file.info(metrics)$size == 0) stop("Empty metrics: ", metrics)
  # Indexed region query is only an efficient extraction. SYMBOL determines membership.
  region_file <- file.path(out_dir, "intermediate", paste0("Auto_MUC16_", sample_id, ".vcf"))
  run_checked("bcftools", c("view", "-r", "chr19:8800000-9100000",
                           shQuote(vcf_paths[[sample_id]])), region_file)
  lines <- readLines(region_file)
  csq_header <- lines[grepl("^##INFO=<ID=CSQ,", lines)]
  if (length(csq_header) != 1L) stop("Missing CSQ header")
  csq_fields <- strsplit(sub('".*$', '', sub('.*Format: ', '', csq_header)), "|", fixed = TRUE)[[1]]
  body <- lines[!startsWith(lines, "#")]
  for (line in body) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    info <- strsplit(fields[8], ";", fixed = TRUE)[[1]]
    csq <- sub("^CSQ=", "", info[startsWith(info, "CSQ=")])
    if (length(csq) != 1L) stop("Missing variant CSQ")
    annotations <- lapply(strsplit(csq, ",", fixed = TRUE)[[1]], function(x) {
      values <- strsplit(x, "|", fixed = TRUE)[[1]]
      setNames(c(values, rep("", length(csq_fields) - length(values))), csq_fields)
    })
    muc <- Filter(function(x) x[["SYMBOL"]] == "MUC16", annotations)
    if (!length(muc)) next
    format_fields <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
    normal <- setNames(strsplit(fields[10], ":", fixed = TRUE)[[1]], format_fields)
    tumor <- setNames(strsplit(fields[11], ":", fixed = TRUE)[[1]], format_fields)
    alts <- strsplit(fields[5], ",", fixed = TRUE)[[1]]
    for (j in seq_along(alts)) {
      get_allele <- function(x, key, index) {
        if (!key %in% names(x)) stop("Missing FORMAT field: ", key)
        values <- strsplit(x[[key]], ",", fixed = TRUE)[[1]]
        if (length(values) < index) stop("Incomplete FORMAT field: ", key)
        values[index]
      }
      records[[length(records) + 1L]] <- data.frame(
        chrom = fields[1], pos = as.integer(fields[2]), ref = fields[4], alt = alts[j],
        sample = sample_id, filter = fields[7],
        consequence = paste(unique(vapply(muc, function(x) x[["Consequence"]], "")), collapse = ";"),
        normal_ref_reads = get_allele(normal, "AD", 1), normal_alt_reads = get_allele(normal, "AD", j + 1L),
        normal_af = get_allele(normal, "AF", j),
        pdo_ref_reads = get_allele(tumor, "AD", 1), pdo_alt_reads = get_allele(tumor, "AD", j + 1L),
        pdo_af = get_allele(tumor, "AF", j), pdo_depth = tumor[["DP"]], stringsAsFactors = FALSE)
    }
  }
}
if (!length(records)) stop("No MUC16 calls")
all_records <- do.call(rbind, records)
all_records$key <- with(all_records, paste(chrom, pos, ref, alt, sep = ":"))
if (anyDuplicated(all_records[c("key", "sample")])) stop("Duplicate MUC16 allele/sample")
keys <- unique(all_records$key)
comparison <- do.call(rbind, lapply(keys, function(key) {
  rows <- all_records[all_records$key == key, ]
  result <- rows[1, c("chrom", "pos", "ref", "alt")]
  result$consequence <- paste(unique(rows$consequence), collapse = ";")
  for (sample_id in samples) {
    row <- rows[rows$sample == sample_id, ]
    result[[paste0("PDO_", sample_id, "_status")]] <- if (!nrow(row)) "NOT_CALLED (absence unproven)" else if (row$filter == "PASS") "PASS" else "FILTERED"
    for (field in c("filter", "normal_ref_reads", "normal_alt_reads", "normal_af", "pdo_ref_reads", "pdo_alt_reads", "pdo_af", "pdo_depth")) {
      result[[paste0("PDO_", sample_id, "_", field)]] <- if (nrow(row)) row[[field]] else NA_character_
    }
  }
  result
}))
comparison <- comparison[order(comparison$chrom, comparison$pos, comparison$ref, comparison$alt), ]
write.csv(comparison, file.path(out_dir, "tables/Auto_MUC16_all_variant_status.csv"), row.names = FALSE)
pass_snps <- unique(all_records[all_records$filter == "PASS" & nchar(all_records$ref) == 1 &
                              nchar(all_records$alt) == 1, c("chrom", "pos", "ref", "alt")])
if (!nrow(pass_snps)) stop("No PASS MUC16 SNPs for read verification")
bed <- file.path(out_dir, "intermediate/Auto_MUC16_PASS_sites.bed")
write.table(unique(data.frame(pass_snps$chrom, pass_snps$pos - 1L, pass_snps$pos)),
            bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
support <- list()
for (sample_id in samples) {
  pileup <- file.path(out_dir, "intermediate", paste0("Auto_MUC16_PASS_", sample_id, ".pileup"))
  # Region limits CRAM I/O; BED restricts observations to the requested SNPs.
  run_checked("samtools", c("mpileup", "-B", "-q", "20", "-Q", "20", "-d", "100000",
    "-f", shQuote(reference), "-r", "chr19:8800000-9100000", "-l", shQuote(bed),
    shQuote(cram_paths[[sample_id]])), pileup)
  pile_lines <- readLines(pileup)
  for (i in seq_len(nrow(pass_snps))) {
    variant <- pass_snps[i, ]
    match <- Filter(function(x) x[1] == variant$chrom && as.integer(x[2]) == variant$pos,
                    strsplit(pile_lines, "\t", fixed = TRUE))
    counts <- setNames(integer(4), c("A", "C", "G", "T"))
    depth <- 0L
    if (length(match) > 1L) stop("Duplicate pileup position")
    if (length(match)) {
      x <- match[[1]]
      if (toupper(x[3]) != variant$ref) stop("Reference mismatch at ", variant$pos)
      depth <- as.integer(x[4])
      bases <- x[5]
      k <- 1L
      while (k <= nchar(bases)) {
        b <- substr(bases, k, k)
        if (b == "^") { k <- k + 2L; next }
        if (b %in% c("+", "-")) {
          k <- k + 1L; digits <- ""
          while (k <= nchar(bases) && grepl("[0-9]", substr(bases, k, k))) {
            digits <- paste0(digits, substr(bases, k, k)); k <- k + 1L
          }
          if (!nzchar(digits)) stop("Invalid pileup indel")
          k <- k + as.integer(digits); next
        }
        if (b %in% c(".", ",")) b <- variant$ref
        b <- toupper(b)
        if (b %in% names(counts)) counts[b] <- counts[b] + 1L
        k <- k + 1L
      }
    }
    called <- all_records[all_records$sample == sample_id & all_records$chrom == variant$chrom &
      all_records$pos == variant$pos & all_records$ref == variant$ref & all_records$alt == variant$alt, ]
    support[[length(support) + 1L]] <- data.frame(variant, sample = sample_id,
      vcf_status = if (nrow(called)) called$filter else "NOT_CALLED",
      pileup_depth = depth, ref_reads = counts[variant$ref], alt_reads = counts[variant$alt],
      other_acgt_reads = sum(counts) - counts[variant$ref] - counts[variant$alt],
      observed_alt_fraction = if (sum(counts) > 0) unname(counts[variant$alt] / sum(counts)) else NA_real_,
      read_status = if (!sum(counts)) "NO_USABLE_ACGT_COVERAGE" else if (!counts[variant$alt])
        "NO_ALT_READS_OBSERVED" else "ALT_READS_OBSERVED", row.names = NULL)
  }
}
support <- do.call(rbind, support)
write.csv(support, file.path(out_dir, "tables/Auto_MUC16_PASS_SNP_read_support.csv"), row.names = FALSE)
writeLines(c("# MUC16 comparison", "", paste("Generated:", Sys.time()), "",
  "All-variant status joins exact CHROM/POS/REF/ALT alleles from both annotated VCFs, including failed calls.",
  "NOT_CALLED means no identical VCF allele record; it does not establish a wild-type genotype.",
  "Multiallelic calls retain their record-level FILTER; transcript consequences are retained from MUC16 annotations.",
  "PASS SNPs are checked in both PDO CRAMs with mapping/base quality >=20, BAQ disabled, maximum depth 100000.",
  "Pileup uses default samtools exclusions and overlap handling; its counts can differ from Mutect2 AD/AF.",
  "No alternate reads at finite depth does not exclude a lower-frequency variant. This is not a new variant call.",
  "The MUC16 target plot uses measured pileup; other noncall bars show reference context only."),
  file.path(out_dir, "reports/Auto_MUC16_interpretation.md"))
writeLines(c(capture.output(sessionInfo()), "", capture.output(system2("bcftools", "--version", stdout = TRUE)),
             capture.output(system2("samtools", "--version", stdout = TRUE))),
           file.path(out_dir, "logs/Auto_metrics_session.txt"))
print(support)
cat("MUC16 allele rows:", nrow(comparison), "\n")
for (sample_id in c("1072", "1090")) {
  metrics_path <- file.path(OUT_DIR, paste0("PDO_", sample_id, "_variants_metrics.tsv"))
  pdo_require_files(metrics_path)
  field_counts <- count.fields(metrics_path, sep = "\t", quote = "", comment.char = "")
  # bcftools emits a trailing tab after the final sample field.
  if (!length(field_counts) || any(!field_counts %in% c(17L, 18L))) stop("Invalid metrics columns")
}
####################
df_1072 <- read_variant_metrics(file.path(OUT_DIR, "PDO_1072_variants_metrics.tsv"), "PDO_1072")
df_1090 <- read_variant_metrics(file.path(OUT_DIR, "PDO_1090_variants_metrics.tsv"), "PDO_1090")

# Apply Quality Filters
if (nrow(df_1072) > 0) df_1072 <- df_1072 %>% filter(dp >= 10 | is.na(dp))
if (nrow(df_1090) > 0) df_1090 <- df_1090 %>% filter(dp >= 10 | is.na(dp))

# ---------------------------------------------------------
# INDEL ANALYSIS (Original Logic for Sequence Alignments)
# ---------------------------------------------------------

top_1072 <- data.frame()
if (nrow(df_1072) > 0) {
  top_1072 <- df_1072 %>%
    filter(!is_snp, indel_len >= 5, grepl("PASS", filter)) %>%
    arrange(desc(af_tumor), desc(indel_len)) %>%
    head(3) %>%
    mutate(label = paste0("Intergenic (", chrom, ") | Deletion -", indel_len, "bp (Exclusive to 1072)\nAllele Freq: ", round(af_tumor*100, 1), "%"))
}

top_1090 <- data.frame()
if (nrow(df_1090) > 0) {
  top_1090 <- df_1090 %>%
    filter(!is_snp, indel_len >= 5, grepl("PASS", filter)) %>%
    arrange(desc(af_tumor), desc(indel_len)) %>%
    head(3) %>%
    mutate(label = paste0("Intergenic (", chrom, ") | Deletion -", indel_len, "bp (Exclusive to 1090)\nAllele Freq: ", round(af_tumor*100, 1), "%"))
}

# Define controls (hardcoded known clean regions)
controls <- data.frame(
  sample = c("Control_WT", "Control_Shared"),
  chrom = c("chr22", "chr2"),
  pos = c(29489100, 6861116),
  ref = c("N", "ACACACACACACC"),
  alt = c("N", "A"),
  indel_len = c(0, 12),
  af_tumor = c(NA, NA),
  label = c("Wild-Type Control (chr22) | Absent in both samples\nContinuous, normal depth in both samples",
            "Shared Control (chr2) | Deletion -12bp present in BOTH samples\nIdentical drop in sequencing depth is visible in both samples")
)

# Combine for sequence alignment plotting
all_indels <- bind_rows(top_1072, top_1090, controls[controls$sample == "Control_Shared", ])

# Depth Extraction using samtools
CRAM1072 <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/preprocessing/recalibrated/PDO_1072/PDO_1072.recal.cram"
CRAM1090 <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/preprocessing/recalibrated/PDO_1090/PDO_1090.recal.cram"

depth_list <- list()
extract_depth <- function(chr, pos, len, region_name, label_text) {
  start_coord <- max(1, pos - 100)
  end_coord <- pos + max(100, len + 50)
  reg_coord <- paste0(chr, ":", start_coord, "-", end_coord)
  
  cmd_1072 <- sprintf("samtools depth -r %s %s", reg_coord, CRAM1072)
  cmd_1090 <- sprintf("samtools depth -r %s %s", reg_coord, CRAM1090)
  
####################
  depth_1072 <- system(cmd_1072, intern = TRUE)
  depth_1090 <- system(cmd_1090, intern = TRUE)
  if (!is.null(attr(depth_1072, "status")) || !is.null(attr(depth_1090, "status"))) {
    stop("samtools depth failed for ", reg_coord)
  }
  if (!length(depth_1072) || !length(depth_1090)) stop("Missing depth for ", reg_coord)
  d1 <- read.table(text = depth_1072, stringsAsFactors = FALSE)
  d2 <- read.table(text = depth_1090, stringsAsFactors = FALSE)
####################
  
  if (nrow(d1) > 0) { colnames(d1) <- c("chrom", "pos", "depth"); d1$sample <- "PDO_1072"; d1$region <- region_name; d1$label <- label_text }
  if (nrow(d2) > 0) { colnames(d2) <- c("chrom", "pos", "depth"); d2$sample <- "PDO_1090"; d2$region <- region_name; d2$label <- label_text }
  
  res <- bind_rows(if(nrow(d1)>0) d1 else NULL, if(nrow(d2)>0) d2 else NULL)
  
  if (region_name != "Control_WT") {
    if(region_name == "Control_Shared") {
      res$indel_start <- pos - 23
      res$indel_end <- pos - 11
    } else {
      res$indel_start <- pos
      res$indel_end <- pos + len
    }
  } else {
    res$indel_start <- NA
    res$indel_end <- NA
  }
  return(res)
}

if (nrow(top_1072) > 0) {
  for (i in 1:nrow(top_1072)) {
    depth_list[[length(depth_list)+1]] <- extract_depth(top_1072$chrom[i], top_1072$pos[i], top_1072$indel_len[i], paste0("1072_Excl_", i), top_1072$label[i])
  }
}
if (nrow(top_1090) > 0) {
  for (i in 1:nrow(top_1090)) {
    depth_list[[length(depth_list)+1]] <- extract_depth(top_1090$chrom[i], top_1090$pos[i], top_1090$indel_len[i], paste0("1090_Excl_", i), top_1090$label[i])
  }
}
depth_list[[length(depth_list)+1]] <- extract_depth(controls$chrom[1], controls$pos[1], controls$indel_len[1], controls$sample[1], controls$label[1])
depth_list[[length(depth_list)+1]] <- extract_depth(controls$chrom[2], controls$pos[2], controls$indel_len[2], controls$sample[2], controls$label[2])

df_depth <- bind_rows(depth_list)
####################
if (!nrow(df_depth)) stop("No depth profiles produced")
write.csv(df_depth, file.path(OUT_DIR, "tables/Auto_indel_depth_profiles.csv"), row.names = FALSE)
saveRDS(list(metrics_1072 = df_1072, metrics_1090 = df_1090, selected_indels = all_indels,
             controls = controls, depth_profiles = df_depth),
        file.path(OUT_DIR, "intermediate/Auto_plot_sources.rds"))
####################

if (nrow(df_depth) > 0) {
  df_depth$sample <- factor(df_depth$sample, levels = c("PDO_1090", "PDO_1072"))
  pdf_out <- file.path(OUT_DIR, "Auto_indel_chromatogram_style.pdf")
  pdf(pdf_out, width = 11, height = 7, onefile = TRUE)
  for (reg in unique(df_depth$region)) {
    df_sub <- df_depth %>% filter(region == reg)
    if(nrow(df_sub) == 0) next
    p <- ggplot(df_sub, aes(x = pos, y = depth)) +
      geom_area(fill = "#5D8AA8", alpha = 0.5) +
      geom_line(color = "#1D428A", linewidth = 0.8) +
      facet_wrap(~sample, ncol = 1, scales = "free_y") +
      labs(title = df_sub$label[1],
           subtitle = "Top panel: Sample 1090 | Bottom panel: Sample 1072\nA structural deletion appears as a sudden drop in sequencing depth.",
           x = paste("Genomic Position (", df_sub$chrom[1], ")", sep=""),
           y = "Sequencing Depth (Reads)") +
      theme_minimal(base_size = 14) +
      theme(
        strip.text = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "grey80"),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40", size = 11),
        plot.margin = margin(20, 20, 20, 20)
      )
    if(!is.na(df_sub$indel_start[1])) {
      start_v <- df_sub$indel_start[1]
      end_v <- df_sub$indel_end[1]
      p <- p + geom_segment(aes(x = start_v, xend = end_v, y = 0, yend = 0), color = "red", linewidth = 4) +
               geom_vline(xintercept = c(start_v, end_v), color="red", linetype="dashed", alpha=0.5)
    }
    print(p)
  }
  dev.off()
}

plot_indel_alignment <- function(row_data) {
  ref <- row_data$ref; alt <- row_data$alt
  len <- abs(nchar(ref) - nchar(alt))
  type <- ifelse(nchar(ref) > nchar(alt), "deletion", "insertion")
  max_viz_len <- 40
  if (type == "deletion") {
    viz_len <- min(nchar(ref), max_viz_len)
    ref_viz <- substr(ref, 1, viz_len)
    alt_viz <- paste0(substr(alt, 1, 1), strrep("-", viz_len - 1))
    if (nchar(ref) > max_viz_len) { ref_viz <- paste0(ref_viz, "..."); alt_viz <- paste0(alt_viz, "...") }
  } else {
    viz_len <- min(nchar(alt), max_viz_len)
    ref_viz <- paste0(substr(ref, 1, 1), strrep("-", viz_len - 1))
    alt_viz <- substr(alt, 1, viz_len)
    if (nchar(alt) > max_viz_len) { ref_viz <- paste0(ref_viz, "..."); alt_viz <- paste0(alt_viz, "...") }
  }
  ref_chars <- strsplit(ref_viz, "")[[1]]; alt_chars <- strsplit(alt_viz, "")[[1]]
  n_pos <- max(length(ref_chars), length(alt_chars))
  df <- data.frame(pos = rep(1:n_pos, 2), allele = rep(c("REF", "ALT"), each = n_pos), base = c(ref_chars, alt_chars))
  base_colors <- c("A" = "#4daf4a", "C" = "#377eb8", "G" = "#ff7f00", "T" = "#e41a1c", "-" = "#cccccc", "."="white")
  p <- ggplot(df, aes(x = pos, y = allele, fill = base)) +
    geom_tile(color = "white", linewidth = 1, width = 0.9, height = 0.8) +
    geom_text(aes(label = base), size = 5, fontface = "bold", color = ifelse(df$base == "-", "transparent", "white")) +
    scale_fill_manual(values = base_colors) + scale_y_discrete(limits = c("ALT", "REF")) +
    labs(title = paste0("Sample: ", row_data$sample, " | ", row_data$chrom, ":", row_data$pos), subtitle = paste0(toupper(type), " of ", len, "bp"), x = NULL, y = NULL) +
    theme_minimal(base_size = 14) + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "none")
  return(p)
}

if (nrow(all_indels) > 0) {
  plots <- lapply(1:nrow(all_indels), function(i) plot_indel_alignment(all_indels[i, ]))
  combined_plot <- wrap_plots(plots, ncol = 1) +
    plot_annotation(title = "Somatic Indels: Sequence-Level Visualizations", subtitle = "Reference vs Alternate allele alignments for top exclusive indels")
  ggsave(file.path(OUT_DIR, "Auto_indel_sequence_alignments.png"), combined_plot, width = 12, height = 2.5 * nrow(all_indels), dpi = 300)
}

# ---------------------------------------------------------
# NEW SNP CLUSTER & PASS PRIORITIZATION
# ---------------------------------------------------------
all_vars <- bind_rows(df_1072, df_1090)

if (nrow(all_vars) > 0) {
  # 1. Filter for PASS variants and strict Exon coding
  all_vars <- all_vars %>% filter(grepl("PASS", filter))
  
  all_vars$var_id <- paste(all_vars$chrom, all_vars$pos, all_vars$ref, all_vars$alt, sep="_")
  all_vars$is_coding <- grepl("missense_variant|synonymous_variant|stop_gained|stop_lost|frameshift_variant|inframe_insertion|inframe_deletion|splice_acceptor_variant|splice_donor_variant", all_vars$consequence)
  
  all_vars <- all_vars %>% filter(gene != "" & !is.na(gene) & is_coding)
  
  # Identify variants unique to one sample
  counts <- all_vars %>% group_by(var_id) %>% summarise(n_samples = n_distinct(sample))
  all_vars <- left_join(all_vars, counts, by="var_id")
  unique_vars <- all_vars %>% filter(n_samples == 1)
  
  # 2. Identify 500bp SNP clusters (Coding)
  snp_clusters <- unique_vars %>%
    filter(is_snp) %>%
    group_by(gene, sample, chrom) %>%
    arrange(pos) %>%
    mutate(
      prev_pos = lag(pos),
      dist_to_prev = pos - prev_pos,
      next_pos = lead(pos),
      dist_to_next = next_pos - pos
    ) %>%
    filter((!is.na(dist_to_prev) & dist_to_prev <= 500) | (!is.na(dist_to_next) & dist_to_next <= 500))
  
  cluster_genes <- unique(snp_clusters$gene)
  
  # All unique PASS variants (Cluster + Single), sorted by priority
  prioritized_vars <- unique_vars %>%
    mutate(
      priority_reason = ifelse(gene %in% cluster_genes & is_snp, "500bp SNP Cluster", "High-Confidence PASS")
    ) %>%
    distinct(var_id, sample, .keep_all=TRUE) %>%
    arrange(desc(priority_reason == "500bp SNP Cluster"), desc(af_tumor), gene, pos) %>%
    select(gene, chrom, pos, ref, alt, sample, priority_reason, consequence, af_tumor, af_normal, dp, tlod, filter, everything())
  
  # Output single detailed table sorted by priority reason -> gene -> position
  write.csv(prioritized_vars, file.path(OUT_DIR, "Auto_variant_prioritization_details.csv"), row.names=FALSE)
  
  if(nrow(prioritized_vars) > 0) {
####################
    selected_variants <- prioritized_vars %>% mutate(variant_rank = row_number()) %>% head(18)
    top_genes <- unique(selected_variants$gene)
    write.csv(selected_variants, file.path(OUT_DIR, "tables/Auto_variant_plot_selection.csv"), row.names = FALSE)
    muc16_support <- read.csv(file.path(OUT_DIR, "tables/Auto_MUC16_PASS_SNP_read_support.csv"),
                             stringsAsFactors = FALSE)
    variant_plot_data <- list()
####################
    
    fasta_path <- "/rds/general/project/tumourheterogeneity1/live/demultiplex/genome.fa"
    base_colors <- c("A" = "#4daf4a", "C" = "#377eb8", "G" = "#ff7f00", "T" = "#e41a1c", "-" = "#cccccc", "N" = "#999999")
    
    pdf_path <- file.path(OUT_DIR, "Auto_variant_clusters_plot.pdf")
    pdf(pdf_path, width=20, height=8)
    
    for (g in top_genes) {
####################
      g_vars <- selected_variants %>% filter(gene == g) %>% arrange(pos)
####################
      chrom <- g_vars$chrom[1]
      min_pos <- min(g_vars$pos)
      max_pos <- max(g_vars$pos)
      
      # 20bp padding
      pad <- 20
      start_pos <- min_pos - pad
      end_pos <- max_pos + pad
      
      # Fetch sequence
      cmd <- sprintf("samtools faidx %s %s:%d-%d | grep -v '>' | tr -d '\\n'", fasta_path, chrom, start_pos, end_pos)
####################
      ref_seq <- system(cmd, intern = TRUE)
      
      if (length(ref_seq) == 0 || nchar(ref_seq) == 0) {
        stop(paste("Could not fetch sequence for", g))
      }
      
      ref_chars <- strsplit(toupper(ref_seq), "")[[1]]
      pos_seq <- start_pos:end_pos
      
      if(length(ref_chars) != length(pos_seq)) {
        stop("Reference length mismatch for ", g)
      }
####################
      
      plot_rows <- list()
      mut_labels <- c()
      
      for (s in c("PDO_1072", "PDO_1090")) {
        for (idx in 1:length(pos_seq)) {
          p <- pos_seq[idx]
          ref_b <- ref_chars[idx]
          
          # Check if this position is mutated in this sample
          mut_match <- g_vars[g_vars$pos == p & g_vars$sample == s, ]
          
          if (nrow(mut_match) > 0) {
            row <- mut_match[1,]
            alt_b <- row$alt
            af <- row$af_tumor
            
            # Alt bar
            plot_rows[[length(plot_rows)+1]] <- data.frame(
              Sample = s, Position = p, Nucleotide = alt_b, Proportion = af
            )
            # Ref bar
            plot_rows[[length(plot_rows)+1]] <- data.frame(
              Sample = s, Position = p, Nucleotide = ref_b, Proportion = 1 - af
            )
            mut_labels <- c(mut_labels, paste0(p, " (", row$ref, ">", alt_b, " AF:", round(af,2), ")"))
          } else {
            # Not mutated
            plot_rows[[length(plot_rows)+1]] <- data.frame(
              Sample = s, Position = p, Nucleotide = ref_b, Proportion = 1.0
            )
          }
        }
      }
      
      plot_data <- bind_rows(plot_rows)
####################
      evidence_label <- "Variant bars use Mutect2 AF; other bases show reference context, not measured wild-type calls."
      if (g == "MUC16") {
        muc16_target <- g_vars[1, ]
        target_support <- muc16_support %>%
          filter(chrom == muc16_target$chrom, pos == muc16_target$pos,
                 ref == muc16_target$ref, alt == muc16_target$alt)
        if (nrow(target_support) != 2L || any(target_support$pileup_depth <= 0L)) {
          stop("MUC16 target requires usable pileup evidence in both PDOs")
        }
        for (sample_id in c("PDO_1072", "PDO_1090")) {
          evidence <- target_support[target_support$sample == sub("PDO_", "", sample_id), ]
          if (nrow(evidence) != 1L) stop("Missing MUC16 sample evidence")
          base_total <- evidence$ref_reads + evidence$alt_reads + evidence$other_acgt_reads
          if (base_total <= 0L || evidence$other_acgt_reads != 0L) {
            stop("MUC16 two-base plot requires positive coverage and no third-base observations")
          }
          plot_data <- plot_data[!(plot_data$Sample == sample_id & plot_data$Position == muc16_target$pos), ]
          plot_data <- bind_rows(plot_data, data.frame(Sample = sample_id,
            Position = muc16_target$pos, Nucleotide = c(muc16_target$ref, muc16_target$alt),
            Proportion = c(evidence$ref_reads, evidence$alt_reads) / base_total))
        }
        evidence_label <- paste0("MUC16 target: Q20/MQ20 pileup; ",
          paste0("PDO_", target_support$sample, " ALT ", target_support$alt_reads, "/",
                 target_support$ref_reads + target_support$alt_reads, collapse = "; "),
          ". Flanks show reference context.")
      }
####################
      plot_data$Sample <- factor(plot_data$Sample, levels=c("PDO_1072", "PDO_1090"))
      
      plot_data$DisplayBase <- sapply(plot_data$Nucleotide, function(x) ifelse(nchar(x) > 1, substr(x, 1, 2), x))
      plot_data$ColorBase <- sapply(plot_data$Nucleotide, function(x) {
         ch <- substr(x, 1, 1)
         if (ch %in% names(base_colors)) ch else "-"
      })
      
      mut_labels <- unique(mut_labels)
      priority_label <- g_vars$priority_reason[1]
      title_str <- paste0("Gene: ", g, " | ", priority_label, " | ", chrom, ":", min_pos, "-", max_pos)
      sub_str <- paste("Variants:", paste(mut_labels, collapse="; "))
####################
      sub_str <- paste0("Priority rank(s): ", paste(g_vars$variant_rank, collapse = ", "),
                        " | ", sub_str, "\n", evidence_label)
      variant_plot_data[[g]] <- list(variants = g_vars, data = plot_data,
                                     title = title_str, subtitle = sub_str)
####################
      
      p_plot <- ggplot(plot_data, aes(x = Position, y = Proportion, fill = ColorBase)) +
        geom_col(width = 0.9, color = "black", linewidth = 0.2) +
        geom_text(aes(label = DisplayBase), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "white") +
        scale_fill_manual(values = base_colors) +
        facet_wrap(~Sample, ncol=1) +
        scale_x_continuous(breaks = seq(start_pos, end_pos, by=5)) +
        theme_minimal(base_size = 14) +
        labs(title = title_str, subtitle = sub_str, x = "Genomic Position", y = "Allele Proportion") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          strip.text = element_text(face = "bold", size=12),
          panel.grid.major.x = element_blank(),
          legend.position = "none"
        )
        
      # Draw boxes around the mutated positions
      for (pos in unique(g_vars$pos)) {
        p_plot <- p_plot + annotate("rect", xmin = pos - 0.5, xmax = pos + 0.5, ymin = 0, ymax = 1,
                          color = "black", fill = NA, linewidth = 1)
      }
        
      print(p_plot)
    }
    
####################
    saveRDS(variant_plot_data, file.path(OUT_DIR, "intermediate/Auto_variant_plot_data.rds"))
####################
    dev.off()
  }
}

cat("Successfully generated sequence visualizations, priority tables, and depth profiles in", OUT_DIR, "\n")
