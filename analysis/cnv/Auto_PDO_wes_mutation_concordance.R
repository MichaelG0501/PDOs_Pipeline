####################
# Analysis registry
# Status: active terminal; analysis/cnv/Auto_PDO_wes_mutation_concordance.R
# Methodology: analysis/methodology/cnv/Auto_wes_mutation_concordance.md
# Map: analysis/ANALYSIS_MAP.md
# Purpose: Li et al. 2018 Figure 2a/b-style plots of existing Mutect2 PASS calls.
# Inputs: /rds/general/project/spatialtranscriptomics/live/sarek_mutect/
#   variant_calling/mutect2/{sample}_vs_NT_{patient}/*.mutect2.filtered.vcf.gz;
#   annotation/mutect2/{sample}_vs_NT_{patient}/*.mutect2.filtered_VEP.ann.vcf.gz;
# Reference: published gene list and digitised approximate frequencies are hard-coded below.
# Samples: PDO 1070,1072,1090,1121,1141,1181; PT 1072,1090,1121.
# Outputs: PDOs_outs/Auto_wes_mutation_concordance/
#   tables/Auto_{all_pass_variants,coding_variants,driver_variants,sample_qc,
#   concordance,concordance_variants,gene_frequencies,reference_frequencies}.csv;
#   figures/Auto_final/Auto_01_driver_mutations_and_concordance.{pdf,png};
#   reports/Auto_results.md; logs/Auto_sessionInfo.txt.
# Downstream: terminal figure and reusable source tables; no state-analysis dependency.
# Cache: none; reruns regenerate the documented current tables, figures and run log.
# Environment: dmtcp; base R only. Resources: 1 CPU, 8 GB, 00:30:00.
# Run: Rscript analysis/cnv/Auto_PDO_wes_mutation_concordance.R
# PBS: qsub job_scripts/cnv/Auto_PDO_wes_mutation_concordance.sh
####################

source("analysis/shared/Auto_pdo_analysis_config.R")
source("analysis/shared/Auto_pdo_analysis_helpers.R")
input_root <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect"
####################
out <- file.path(PDO_OUTPUT_DIR, "Auto_wes_mutation_concordance")
dirs <- pdo_ensure_output_tiers(out, tiers = c("tables", "figures", "logs", "reports"))
####################
final_dir <- file.path(dirs[["figures"]], "Auto_final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
####################
# Approximate frequencies digitised from Li et al. 2018 Figure 2, not exact cohort counts.
# Source: https://www.nature.com/articles/s41467-018-05190-9/figures/2
# Published 2053x887 PNG: black fraction against grey, rounded to percentage points;
# approximate uncertainty 2 percentage points. Pixel coordinates are recorded below.
reference <- data.frame(
  gene = c("TP53", "PIK3CA", "LRP1B", "CDH11", "ERBB4", "PTPRT", "KCNQ3",
           "ARID1A", "PREX2", "CDKN2A", "SMARCA4", "SMAD4", "GNAS", "APC",
           "FAT4", "SACS", "NRG1"),
  mutsig_significant = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
                        FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  secrier_pct_approx = c(76, 7, 21, 6, 7, 11, 10, 15, 6, 15, 8, 16, 7, 9, 21, 8, 12),
  tcga_pct_approx = c(80, 6, 14, 14, 8, 16, 9, 12, 28, 8, 12, 20, 18, 20, 13, 20, 14),
  source = "https://www.nature.com/articles/s41467-018-05190-9/figures/2",
  method = sprintf(paste0("Digitised published 2053x887 PNG; y=%d; Secrier x=918:999; ",
                          "TCGA x=1069:1150; black fraction against local grey; ",
                          "nearest percentage point; approximate uncertainty 2 percentage points"),
                   164 + 37 * (0:16)))
####################
genes <- reference$gene
samples <- c("PT_1072", "PDO_1072", "PT_1090", "PDO_1090", "PT_1121", "PDO_1121",
             "PDO_1070", "PDO_1141", "PDO_1181")
coding_terms <- c("frameshift_variant", "stop_gained", "stop_lost", "start_lost",
                  "inframe_insertion", "inframe_deletion", "missense_variant",
                  "protein_altering_variant")
colours <- setNames(c("#377EB8", "#7B3294", "#A65628", "#E78AC3",
                      "#FFD966", "#E6AB02", "#E41A1C", "#66A61E"), coding_terms)

read_vcf <- function(path) {
  pdo_require_files(path)
  con <- gzfile(path, "rt")
  lines <- readLines(con)
  close(con)
  header <- lines[startsWith(lines, "#CHROM")]
  if (length(header) != 1L) stop("Missing VCF header: ", path)
  columns <- strsplit(sub("^#", "", header), "\t", fixed = TRUE)[[1]]
  records <- lines[!startsWith(lines, "#")]
  if (!length(records)) stop("Empty VCF: ", path)
  tab <- read.delim(text = paste(records, collapse = "\n"), header = FALSE,
                    sep = "\t", quote = "", comment.char = "", colClasses = "character")
  names(tab) <- columns
  tab[tab$FILTER == "PASS", , drop = FALSE]
}

# VEP removes common reference/alternate prefixes and suffixes from its Allele field.
vep_allele <- function(ref, alt) {
  while (nchar(ref) > 0 && nchar(alt) > 0 && substr(ref, 1, 1) == substr(alt, 1, 1)) {
    ref <- substring(ref, 2); alt <- substring(alt, 2)
  }
  while (nchar(ref) > 0 && nchar(alt) > 0 &&
         substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))) {
    ref <- substr(ref, 1, nchar(ref) - 1); alt <- substr(alt, 1, nchar(alt) - 1)
  }
  if (alt == "") "-" else alt
}

all_rows <- coding_rows <- qc_rows <- list()
for (sample in samples) {
  patient <- sub(".*_", "", sample)
  pair <- paste0(sample, "_vs_NT_", patient)
  raw_path <- file.path(input_root, "variant_calling/mutect2", pair,
                        paste0(pair, ".mutect2.filtered.vcf.gz"))
  ann_path <- file.path(input_root, "annotation/mutect2", pair,
                        paste0(pair, ".mutect2.filtered_VEP.ann.vcf.gz"))
  raw <- read_vcf(raw_path)
  ann <- read_vcf(ann_path)
  sample_column <- paste0("SUR", patient, "_", sample)
  if (!sample_column %in% names(raw)) stop("Missing sample column: ", sample_column)
  sample_rows <- list()
  for (i in seq_len(nrow(raw))) {
    r <- raw[i, ]
    alts <- strsplit(r$ALT, ",", fixed = TRUE)[[1]]
    fmt <- setNames(strsplit(r[[sample_column]], ":", fixed = TRUE)[[1]],
                    strsplit(r$FORMAT, ":", fixed = TRUE)[[1]])
    if (!all(c("AF", "AD", "DP") %in% names(fmt))) stop("Missing AF/AD/DP")
    af <- as.numeric(strsplit(fmt[["AF"]], ",", fixed = TRUE)[[1]])
    ad <- as.numeric(strsplit(fmt[["AD"]], ",", fixed = TRUE)[[1]])
    if (length(af) != length(alts) || length(ad) != length(alts) + 1L) stop("Invalid allele fields")
    for (a in seq_along(alts)) {
      if (!grepl("^[ACGT]+$", r$REF) || !grepl("^[ACGT]+$", alts[a])) stop("Unsupported allele")
      sample_rows[[length(sample_rows) + 1L]] <- data.frame(
        sample, patient, type = sub("_.*", "", sample), chrom = r$CHROM,
        pos = as.integer(r$POS), ref = r$REF, alt = alts[a],
        key = paste(r$CHROM, r$POS, r$REF, alts[a], sep = ":"),
        vaf = af[a], ref_reads = ad[1], alt_reads = ad[a + 1], depth = as.integer(fmt[["DP"]]))
    }
  }
  variants <- do.call(rbind, sample_rows)
  if (is.null(variants) || anyDuplicated(variants$key) ||
      anyNA(variants$vaf) || any(variants$vaf < 0 | variants$vaf > 1)) stop("Invalid PASS variants")
  all_rows[[sample]] <- variants
  for (i in seq_len(nrow(ann))) {
    r <- ann[i, ]
    csq <- strsplit(r$INFO, ";", fixed = TRUE)[[1]]
    csq <- csq[startsWith(csq, "CSQ=")]
    if (length(csq) != 1L) stop("Missing VEP annotation")
    alts <- strsplit(r$ALT, ",", fixed = TRUE)[[1]]
    va <- vapply(alts, function(a) vep_allele(r$REF, a), character(1))
    if (anyDuplicated(va)) stop("Ambiguous VEP allele mapping")
    for (entry in strsplit(sub("^CSQ=", "", csq), ",", fixed = TRUE)[[1]]) {
      fields <- strsplit(entry, "|", fixed = TRUE)[[1]]
      terms <- intersect(coding_terms, strsplit(fields[2], "&", fixed = TRUE)[[1]])
      if (!length(terms) || fields[8] != "protein_coding" || !nzchar(fields[4])) next
      a <- match(fields[1], va)
      if (is.na(a)) stop("Unmatched VEP allele")
      key <- paste(r$CHROM, r$POS, r$REF, alts[a], sep = ":")
      v <- variants[match(key, variants$key), , drop = FALSE]
      if (anyNA(v$key)) stop("Annotated PASS variant absent from original PASS VCF")
      coding_rows[[length(coding_rows) + 1L]] <- cbind(v, gene = fields[4],
                                                     consequence = terms[1])
    }
  }
  qc_rows[[sample]] <- data.frame(sample, raw_pass_records = nrow(raw),
                                  raw_pass_alleles = nrow(variants), vep_pass_records = nrow(ann))
  message(sample, ": ", nrow(variants), " PASS alleles")
}
variants <- do.call(rbind, all_rows)
coding <- unique(do.call(rbind, coding_rows))
if (is.null(coding) || !nrow(coding)) stop("No nonsynonymous coding variants")
coding <- coding[order(match(coding$consequence, coding_terms)), ]
coding <- coding[!duplicated(paste(coding$sample, coding$key, coding$gene)), ]
drivers <- coding[coding$gene %in% genes, ]
if (!nrow(drivers)) stop("No mutations in the published driver list")
concordance <- membership <- list()
for (patient in c("1072", "1090", "1121")) {
  for (subset in c("Total", "Nonsynonymous")) {
    x <- if (subset == "Total") variants else coding
    pdo <- unique(x$key[x$sample == paste0("PDO_", patient)])
    tumour <- unique(x$key[x$sample == paste0("PT_", patient)])
    keys <- union(pdo, tumour)
    if (!length(keys)) stop("Empty comparison: ", patient, " ", subset)
    category <- ifelse(keys %in% pdo & keys %in% tumour, "Shared",
                       ifelse(keys %in% pdo, "PDO only", "Tumour only"))
    membership[[length(membership) + 1L]] <- data.frame(patient, subset, key = keys, category)
    counts <- table(factor(category, levels = c("PDO only", "Shared", "Tumour only")))
    concordance[[length(concordance) + 1L]] <- data.frame(
      patient, subset, category = names(counts), count = as.integer(counts),
      union_count = length(keys), proportion = as.integer(counts) / length(keys))
  }
}
concordance <- do.call(rbind, concordance)
frequency <- data.frame(gene = genes)
for (type in c("PDO", "PT")) {
  denominator <- sum(startsWith(samples, paste0(type, "_")))
  frequency[[paste0(type, "_n")]] <- vapply(genes, function(g)
    length(unique(drivers$sample[drivers$gene == g & drivers$type == type])), integer(1))
  frequency[[paste0(type, "_pct")]] <- 100 * frequency[[paste0(type, "_n")]] / denominator
}
tables <- list(all_pass_variants = variants, coding_variants = coding, driver_variants = drivers,
               sample_qc = do.call(rbind, qc_rows), concordance = concordance,
               concordance_variants = do.call(rbind, membership), gene_frequencies = frequency,
               reference_frequencies = reference)
for (name in names(tables))
  write.csv(tables[[name]], file.path(dirs[["tables"]], paste0("Auto_", name, ".csv")), row.names = FALSE)
stopifnot(all(abs(aggregate(proportion ~ patient + subset, concordance, sum)$proportion - 1) < 1e-10))

draw_figure <- function() {
  x <- c(1, 1.95, 3.35, 4.3, 5.7, 6.65, 8.3, 9.6, 10.9)
  x_limits <- c(0.5, 11.4)
  layout(matrix(c(1, 7, 7, 7, 0, 2, 3, 4, 5, 6, 8, 8, 8, 8, 8),
                nrow = 3, byrow = TRUE), widths = c(5.6, 1.2, 1.35, 1.35, 3.5),
         heights = c(1.2, 6.8, 0.65))
  par(mar = c(0.3, 6, 2.4, 0.6), xaxs = "i")
  counts <- table(factor(drivers$consequence, levels = coding_terms),
                  factor(drivers$sample, levels = samples))
  plot(NA, xlim = x_limits, ylim = c(0, max(colSums(counts)) + 1),
       xaxt = "n", xlab = "", ylab = "Variants", bty = "n",
       main = "a  Driver SNVs/indels", cex.main = 1, adj = 0)
  for (j in seq_along(samples)) {
    h <- c(0, cumsum(counts[, j]))
    rect(x[j] - 0.42, head(h, -1), x[j] + 0.42, tail(h, -1), col = colours, border = NA)
  }
  par(mar = c(4.4, 6, 1.7, 0.6))
  plot(NA, xlim = x_limits, ylim = c(length(genes) + 0.5, 0.5),
       axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  for (i in seq_along(genes)) for (j in seq_along(samples)) {
    rect(x[j] - 0.42, i - 0.44, x[j] + 0.42, i + 0.44, col = "#D9D9D9", border = NA)
    d <- drivers[drivers$gene == genes[i] & drivers$sample == samples[j], ]
    if (nrow(d)) {
      edges <- seq(i - 0.44, i + 0.44, length.out = nrow(d) + 1)
      rect(x[j] - 0.42, head(edges, -1), x[j] + 0.42, tail(edges, -1),
           col = colours[d$consequence], border = "white")
      text(x[j], (head(edges, -1) + tail(edges, -1))/2, sprintf("%.2f", d$vaf),
           col = ifelse(d$consequence %in% c("inframe_insertion", "inframe_deletion"), "black", "white"),
           cex = min(0.7, 1.2/nrow(d)))
    }
  }
  axis(2, seq_along(genes), paste0(genes, ifelse(reference$mutsig_significant, "*", "")),
       las = 2, tick = FALSE, cex.axis = 0.85)
  axis(1, x, ifelse(startsWith(samples, "PT"), "T", "O"), tick = FALSE,
       cex.axis = 0.85, line = -0.3)
  text(c(mean(x[1:2]), mean(x[3:4]), mean(x[5:6]), x[7:9]), length(genes) + 1.7,
       c("SUR1072", "SUR1090", "SUR1121", "SUR1070", "SUR1141", "SUR1181"),
       srt = 45, adj = 1, xpd = NA, cex = 0.8)
  bar_panel <- function(values, title) {
    par(mar = c(4.4, 0.4, 1.7, 0.4))
    plot(NA, xlim = c(0, 105), ylim = c(length(genes) + 0.5, 0.5),
         axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    rect(0, seq_along(genes) - 0.44, 100, seq_along(genes) + 0.44, col = "#D9D9D9", border = NA)
    rect(0, seq_along(genes) - 0.44, values, seq_along(genes) + 0.44, col = "#202020", border = NA)
    axis(3, c(0, 50, 100), c("0", "50", "100%"), cex.axis = 0.65)
    title_lines <- strsplit(title, "\n", fixed = TRUE)[[1]]
    for (i in seq_along(title_lines))
      mtext(title_lines[i], side = 1, line = i, cex = 0.72)
  }
  bar_panel(frequency$PDO_pct, "PDOs\nn = 6")
  bar_panel(reference$secrier_pct_approx, "Secrier & Li\nn = 129\napproximate")
  bar_panel(reference$tcga_pct_approx, "TCGA OAC\nn = 89\napproximate")
  par(mar = c(4.4, 5, 1.7, 1))
  plot(NA, xlim = c(0, 1), ylim = c(8.2, 0), axes = FALSE, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i")
  mtext("b  SNVs/indels", side = 3, line = 2.6,
        font = 2, cex = 0.95, xpd = NA)
  axis(3, seq(0, 1, 0.25), paste0(seq(0, 100, 25), "%"), cex.axis = 0.75)
  overlap_colours <- c("PDO only" = "#FC8D62", "Shared" = "#808080", "Tumour only" = "#67A9CF")
  for (k in seq_along(c("1072", "1090", "1121"))) {
    patient <- c("1072", "1090", "1121")[k]
    y <- (k - 1) * 2.6 + c(0.8, 1.5)
    text(0, y[1] - 0.5, paste0("SUR", patient), adj = 0, font = 2, cex = 0.85)
    for (s in seq_along(c("Total", "Nonsynonymous"))) {
      d <- concordance[concordance$patient == patient & concordance$subset == c("Total", "Nonsynonymous")[s], ]
      edges <- c(0, cumsum(d$proportion))
      rect(head(edges, -1), y[s] - 0.24, tail(edges, -1), y[s] + 0.24,
           col = overlap_colours[d$category], border = "white")
      text(-0.04, y[s], c("Total", "Nonsyn.")[s], adj = 1, xpd = NA, cex = 0.75)
      text(1, y[s] + 0.38, paste0("union n = ", d$union_count[1]), adj = 1, cex = 0.6)
    }
  }
  legend("bottom", names(overlap_colours), fill = overlap_colours, bty = "n", cex = 0.73)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  used <- coding_terms[coding_terms %in% drivers$consequence]
  legend("center", c(gsub("_", " ", used), "No qualifying call"),
         fill = c(colours[used], "#D9D9D9"), bty = "n", cex = 0.85, ncol = 2)
  plot.new()
  text(0.5, 0.68, "T  Parent tumour     O  PDO     Cell labels: VAF     * Published MutSigCV significance", cex = 0.8)
  text(0.5, 0.25, "Existing PASS calls; overlap is relative to the pair's union. External frequencies approximate (Li et al., 2018).", cex = 0.75)
}
stem <- file.path(final_dir, "Auto_01_driver_mutations_and_concordance")
pdf(paste0(stem, ".pdf"), width = 13.5, height = 8.7, useDingbats = FALSE)
draw_figure()
dev.off()
png(paste0(stem, ".png"), width = 2700, height = 1740, res = 200, type = "cairo")
draw_figure()
dev.off()
####################
writeLines(c("# OAC PDO–tumour mutation results",
             "",
             "- [Final figure PDF](../figures/Auto_final/Auto_01_driver_mutations_and_concordance.pdf)",
             "- [PNG preview](../figures/Auto_final/Auto_01_driver_mutations_and_concordance.png)",
             "- [Variant and frequency source tables](../tables/)",
             "",
             "Panel a shows 17 published genes in six PDOs and three matched tumours. Complete pairs come first; standalone PDOs follow. Cell labels are caller VAFs. PDO mutation frequencies use six patients; external OAC reference frequencies are approximate.",
             "",
             "Panel b shows PDO-only, shared and tumour-only variants as fractions of each pair’s union, for all Mutect2 PASS alleles and the VEP nonsynonymous coding subset. This is callset overlap, not proof of biological absence.",
             "",
             "Inputs: existing matched-normal Mutect2 calls and VEP annotations from sarek_mutect. No signature analysis or variant recall is performed.",
             "",
             "Reference: Li et al., Nature Communications 9, 2983 (2018), https://doi.org/10.1038/s41467-018-05190-9."),
           file.path(dirs[["reports"]], "Auto_results.md"))
####################
writeLines(c(paste("Completed:", Sys.time()), capture.output(sessionInfo())),
           file.path(dirs[["logs"]], "Auto_sessionInfo.txt"))
stopifnot(all(file.info(paste0(stem, c(".pdf", ".png")))$size > 1000))
message("SUCCESS: ", out)
