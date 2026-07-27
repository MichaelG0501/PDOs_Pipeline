####################
# Auto_plot_wes_subclone_results.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_plot_wes_subclone_results.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_input.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_results.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_purity_ploidy.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_segments.tsv
#   - Optional scRNA Numbat conservative summary for context
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/figures/Auto_wes_subclone_visual_summary.pdf
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_cluster_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_reliability_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_scRNA_clone_context.csv
#   - PDOs_outs/Auto_wes_subclone/logs/Auto_wes_subclone_visualisation_summary.tsv
# Downstream use: terminal visual QC and reliability assessment for WES subclone calls.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

wd <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
if (dir.exists(wd)) setwd(wd)

samples <- c("PDO_1090_vs_NT_1090", "PDO_1181_vs_NT_1181")
out_root <- "PDOs_outs/Auto_wes_subclone"
fig_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables/visualisation")
log_dir <- file.path(out_root, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pyclone_dir <- file.path(out_root, "tables/pyclone")
facets_dir <- file.path(out_root, "tables/facets")

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  fread(path)
}

read_sample <- function(sample) {
  input <- read_required(file.path(pyclone_dir, paste0("Auto_", sample, "_pyclone_vi_input.tsv")))
  result <- read_required(file.path(pyclone_dir, paste0("Auto_", sample, "_pyclone_vi_results.tsv")))
  purity <- read_required(file.path(facets_dir, paste0("Auto_", sample, "_facets_purity_ploidy.tsv")))
  segments <- read_required(file.path(facets_dir, paste0("Auto_", sample, "_facets_segments.tsv")))
  dt <- merge(result, input, by = c("mutation_id", "sample_id"), all.x = TRUE)
  dt[, tumour_vaf := alt_counts / pmax(1, ref_counts + alt_counts)]
  dt[, tumour_depth := ref_counts + alt_counts]
  dt[, sample := sample]
  dt[, cluster_id := as.character(cluster_id)]
  list(input = input, result = dt, purity = purity, segments = segments)
}

sample_data <- setNames(lapply(samples, read_sample), samples)
pyclone_results <- rbindlist(lapply(sample_data, `[[`, "result"), use.names = TRUE, fill = TRUE)
purity_tbl <- rbindlist(lapply(sample_data, `[[`, "purity"), use.names = TRUE, fill = TRUE)
segments_tbl <- rbindlist(lapply(sample_data, `[[`, "segments"), use.names = TRUE, fill = TRUE)

cluster_summary <- pyclone_results[, .(
  n_mutations = .N,
  median_ccf = as.numeric(median(cellular_prevalence, na.rm = TRUE)),
  mean_ccf = as.numeric(mean(cellular_prevalence, na.rm = TRUE)),
  median_assignment_prob = as.numeric(median(cluster_assignment_prob, na.rm = TRUE)),
  mean_assignment_prob = as.numeric(mean(cluster_assignment_prob, na.rm = TRUE)),
  frac_assignment_prob_ge_0_8 = as.numeric(mean(cluster_assignment_prob >= 0.8, na.rm = TRUE)),
  median_tumour_vaf = as.numeric(median(tumour_vaf, na.rm = TRUE)),
  median_tumour_depth = as.numeric(median(tumour_depth, na.rm = TRUE))
), by = .(sample, cluster_id)]
setorder(cluster_summary, sample, -median_ccf)

cluster_summary[, cluster_label := paste0(
  "cluster ", cluster_id,
  "\nn=", n_mutations,
  "\nmedian CCF=", round(median_ccf, 2),
  "\nmedian P=", round(median_assignment_prob, 2)
)]

min_ccf_sep <- function(x) {
  x <- sort(unique(round(x[is.finite(x)], 4)))
  if (length(x) < 2) return(NA_real_)
  min(diff(x))
}

reliability_summary <- pyclone_results[, .(
  n_pyclone_variants = .N,
  n_pyclone_clusters = uniqueN(cluster_id),
  median_assignment_prob = median(cluster_assignment_prob, na.rm = TRUE),
  mean_assignment_prob = mean(cluster_assignment_prob, na.rm = TRUE),
  frac_assignment_prob_ge_0_8 = mean(cluster_assignment_prob >= 0.8, na.rm = TRUE),
  frac_assignment_prob_ge_0_7 = mean(cluster_assignment_prob >= 0.7, na.rm = TRUE),
  min_cluster_n = min(.N), # replaced below after cluster join
  min_cluster_fraction = NA_real_,
  min_median_ccf_separation = NA_real_
), by = sample]

cluster_metrics <- cluster_summary[, .(
  min_cluster_n = min(n_mutations),
  min_cluster_fraction = min(n_mutations) / sum(n_mutations),
  min_median_ccf_separation = min_ccf_sep(median_ccf)
), by = sample]
reliability_summary[, c("min_cluster_n", "min_cluster_fraction", "min_median_ccf_separation") := NULL]
reliability_summary <- merge(reliability_summary, cluster_metrics, by = "sample")

if (all(c("sample", "purity", "ploidy") %in% names(purity_tbl))) {
  reliability_summary <- merge(
    reliability_summary,
    unique(purity_tbl[, .(sample, facets_purity = purity, facets_ploidy = ploidy)]),
    by = "sample",
    all.x = TRUE
  )
}

reliability_summary[, reliability_call := fifelse(
  n_pyclone_variants >= 100 &
    min_cluster_n >= 20 &
    median_assignment_prob >= 0.80 &
    frac_assignment_prob_ge_0_8 >= 0.65 &
    min_median_ccf_separation >= 0.15,
  "strong",
  fifelse(
    n_pyclone_variants >= 75 &
      min_cluster_n >= 10 &
      median_assignment_prob >= 0.60 &
      frac_assignment_prob_ge_0_8 >= 0.40 &
      min_median_ccf_separation >= 0.10,
    "moderate",
    "weak/provisional"
  )
)]

reliability_summary[, reliability_notes := fifelse(
  reliability_call == "strong",
  "Supported by variant count, cluster size, assignment probability, and CCF separation.",
  fifelse(
    reliability_call == "moderate",
    "Usable but should be interpreted with copy-number and scRNA concordance checks.",
    "Do not treat the exact cluster count as definitive without sensitivity/re-run support."
  )
)]

numbat_context_path <- "PDOs_outs/Auto_PDO_numbat_subclone_mp_conservative/Auto_PDO_numbat_subclone_summary.csv"
numbat_context <- data.table()
if (file.exists(numbat_context_path)) {
  numbat <- fread(numbat_context_path)
  keep_samples <- c(
    "SUR1090_Untreated_PDO", "SUR1090_Treated_PDO",
    "SUR1181_Untreated_PDO", "SUR1181_Treated_PDO"
  )
  numbat_context <- numbat[sample %in% keep_samples, .(
    scRNA_sample = sample,
    scRNA_status = status,
    n_cells,
    n_raw_numbat_clones,
    n_display_clones,
    clone_mode,
    median_p_cnv,
    state_p_value,
    state_cramers_v
  )]
  numbat_context[, wes_pair := fifelse(grepl("SUR1090", scRNA_sample), "PDO_1090_vs_NT_1090", "PDO_1181_vs_NT_1181")]
}

autosomes <- as.character(1:22)
chr_lengths <- data.table(
  chrom_clean = autosomes,
  chr_len = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468
  )
)
chr_lengths[, chr_start := shift(cumsum(chr_len), fill = 0)]
chr_lengths[, chr_mid := chr_start + chr_len / 2]

segments_plot <- copy(segments_tbl)
segments_plot[, chrom_clean := sub("^chr", "", as.character(chrom))]
segments_plot <- segments_plot[chrom_clean %in% autosomes]
segments_plot <- merge(segments_plot, chr_lengths, by = "chrom_clean", all.x = TRUE)
segments_plot[, genome_start := chr_start + pmax(0, start)]
segments_plot[, genome_end := chr_start + pmax(start, end)]
segments_plot[, total_cn_plot := pmin(total_cn, 8)]

cluster_palette <- c(
  "0" = "#2B6CB0",
  "1" = "#D95F02",
  "2" = "#1B9E77",
  "3" = "#7570B3",
  "4" = "#E7298A",
  "5" = "#666666"
)

make_ccf_plot <- function(sample_id) {
  dt <- pyclone_results[sample == sample_id]
  labels <- cluster_summary[sample == sample_id]
  dt <- merge(dt, labels[, .(cluster_id, median_ccf, cluster_label)], by = "cluster_id")
  dt[, cluster_label := factor(cluster_label, levels = labels[order(median_ccf), cluster_label])]
  ggplot(dt, aes(cellular_prevalence, cluster_label, colour = cluster_id, alpha = cluster_assignment_prob)) +
    geom_jitter(height = 0.14, width = 0, size = 1.8) +
    geom_vline(data = labels, aes(xintercept = median_ccf, colour = cluster_id), linewidth = 0.6, show.legend = FALSE) +
    scale_x_continuous(limits = c(0, 1.02), labels = percent_format(accuracy = 1)) +
    scale_alpha_continuous(range = c(0.2, 1), limits = c(0, 1)) +
    scale_colour_manual(values = cluster_palette, drop = FALSE) +
    labs(title = paste0(sample_id, ": PyClone-VI CCF clusters"), x = "Cancer cell fraction", y = NULL, alpha = "Assignment P", colour = "Cluster") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

make_vaf_ccf_plot <- function(sample_id) {
  dt <- pyclone_results[sample == sample_id]
  ggplot(dt, aes(tumour_vaf, cellular_prevalence, colour = cluster_id, alpha = cluster_assignment_prob)) +
    geom_point(size = 1.8) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.02)) +
    scale_alpha_continuous(range = c(0.2, 1), limits = c(0, 1)) +
    scale_colour_manual(values = cluster_palette, drop = FALSE) +
    labs(title = "Raw VAF versus copy-number-corrected CCF", x = "Tumour VAF", y = "Cancer cell fraction", alpha = "Assignment P", colour = "Cluster") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

make_probability_plot <- function(sample_id) {
  dt <- pyclone_results[sample == sample_id]
  ggplot(dt, aes(cluster_id, cluster_assignment_prob, fill = cluster_id)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey35") +
    geom_boxplot(width = 0.7, outlier.alpha = 0.35) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = cluster_palette, drop = FALSE) +
    labs(title = "Mutation assignment confidence", x = "Cluster", y = "PyClone-VI assignment probability") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
}

make_cluster_size_plot <- function(sample_id) {
  dt <- cluster_summary[sample == sample_id]
  ggplot(dt, aes(reorder(cluster_id, -median_ccf), n_mutations, fill = cluster_id)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0("median CCF ", round(median_ccf, 2))), vjust = -0.35, size = 3) +
    scale_fill_manual(values = cluster_palette, drop = FALSE) +
    labs(title = "Cluster mutation support", x = "Cluster", y = "Mutations") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")
}

make_cna_plot <- function(sample_id) {
  dt <- segments_plot[sample == sample_id]
  ggplot(dt, aes(x = genome_start, xend = genome_end, y = total_cn_plot, yend = total_cn_plot, colour = minor_cn == 0)) +
    geom_segment(linewidth = 1.1, lineend = "butt") +
    scale_x_continuous(breaks = chr_lengths$chr_mid, labels = chr_lengths$chrom_clean, expand = expansion(mult = c(0.005, 0.005))) +
    scale_y_continuous(breaks = 0:8, limits = c(0, 8.5)) +
    scale_colour_manual(values = c("FALSE" = "#4C566A", "TRUE" = "#BF616A"), labels = c("retained minor", "LOH/minor=0")) +
    labs(title = "FACETS autosomal total copy number", x = "Chromosome", y = "Total CN (capped at 8)", colour = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(size = 7), legend.position = "bottom")
}

make_reliability_text <- function(sample_id) {
  rel <- reliability_summary[sample == sample_id]
  paste0(
    "Reliability call: ", rel$reliability_call,
    "\nPyClone-VI clusters: ", rel$n_pyclone_clusters,
    "\nInput variants: ", rel$n_pyclone_variants,
    "\nMedian assignment P: ", round(rel$median_assignment_prob, 3),
    "\nFraction assignment P >= 0.8: ", round(rel$frac_assignment_prob_ge_0_8, 3),
    "\nSmallest cluster: ", rel$min_cluster_n, " variants (", percent(rel$min_cluster_fraction, accuracy = 0.1), ")",
    "\nMinimum median CCF separation: ", round(rel$min_median_ccf_separation, 3),
    "\nFACETS purity/ploidy: ", round(rel$facets_purity, 3), " / ", round(rel$facets_ploidy, 3),
    "\n", rel$reliability_notes
  )
}

summary_pdf <- file.path(fig_dir, "Auto_wes_subclone_visual_summary.pdf")
pdf(summary_pdf, width = 13, height = 9)
for (sample in samples) {
  rel_text <- make_reliability_text(sample)
  text_plot <- ggplot() +
    annotate("text", x = 0, y = 1, label = rel_text, hjust = 0, vjust = 1, size = 4, lineheight = 1.05) +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = "Reliability summary") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0))

  grid.arrange(
    make_ccf_plot(sample),
    make_vaf_ccf_plot(sample),
    make_probability_plot(sample),
    make_cluster_size_plot(sample),
    make_cna_plot(sample),
    text_plot,
    ncol = 2,
    layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6)),
    top = sample
  )
}
dev.off()

fwrite(cluster_summary, file.path(table_dir, "Auto_wes_subclone_cluster_summary.csv"))
fwrite(reliability_summary, file.path(table_dir, "Auto_wes_subclone_reliability_summary.csv"))
if (nrow(numbat_context) > 0) {
  fwrite(numbat_context, file.path(table_dir, "Auto_wes_subclone_scRNA_clone_context.csv"))
}

fwrite(data.table(
  finished = as.character(Sys.time()),
  samples = paste(samples, collapse = ","),
  summary_pdf = summary_pdf,
  cluster_summary = file.path(table_dir, "Auto_wes_subclone_cluster_summary.csv"),
  reliability_summary = file.path(table_dir, "Auto_wes_subclone_reliability_summary.csv"),
  scRNA_context = ifelse(nrow(numbat_context) > 0, file.path(table_dir, "Auto_wes_subclone_scRNA_clone_context.csv"), NA_character_)
), file.path(log_dir, "Auto_wes_subclone_visualisation_summary.tsv"), sep = "\t")

message("Wrote WES subclone visual summary: ", summary_pdf)
