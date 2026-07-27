####################
# Auto_wes_absolute_cna_compare.R
#
# Analysis registry
# Status: active terminal diagnostic
# Script: analysis/cnv/wes_subclone/Auto_wes_absolute_cna_compare.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_segments.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_purity_ploidy.tsv
#   - Sarek CNVkit .somatic.call.cns/.cns under spatialtranscriptomics live
#   - ephemeral PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat outputs
# Outputs:
#   - PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_<sample>_facets_absolute*.cns
#   - PDOs_outs/Auto_wes_absolute_cna/tables/Auto_wes_absolute_cna_*.csv
#   - PDOs_outs/Auto_wes_absolute_cna/figures/Auto_wes_absolute_cna_compare_<sample>.pdf/.png
# Downstream use: terminal WES/scRNA CNA scale validation. The FACETS absolute
#   WES track is diploid-scaled log2(total_cn / 2), not CNVkit centered log2.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

live_project <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
ephemeral_project <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
wes_root <- file.path(live_project, "PDOs_outs/Auto_wes_subclone")
out_root <- file.path(live_project, "PDOs_outs/Auto_wes_absolute_cna")
cnvkit_root <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit"
numbat_root <- file.path(ephemeral_project, "PDOs_outs/Auto_PDO_numbat/by_samples")

fig_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
cns_dir <- file.path(table_dir, "cns")
log_dir <- file.path(out_root, "logs")
for (d in c(fig_dir, table_dir, cns_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

wgs_sample_map <- list(
  "PDO_1090_vs_NT_1090" = c("SUR1090_Untreated_PDO", "SUR1090_Treated_PDO"),
  "PDO_1181_vs_NT_1181" = c("SUR1181_Untreated_PDO", "SUR1181_Treated_PDO")
)

chr_order <- c(as.character(1:22), "X")
chr_sizes <- c(
  `1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
  `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
  `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
  `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
  `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
  `21` = 46709983, `22` = 50818468, X = 156040895
)
chr_cumstart <- cumsum(c(0, chr_sizes[chr_order][-length(chr_order)]))
names(chr_cumstart) <- chr_order
chr_cumend <- chr_cumstart + chr_sizes[chr_order]
chr_mids <- (chr_cumstart + chr_cumend) / 2
genome_len <- max(chr_cumend)

chr_bands <- data.table(
  genome_start = chr_cumstart[chr_order],
  genome_end = chr_cumend[chr_order],
  chr = chr_order,
  fill = ifelse(seq_along(chr_order) %% 2 == 1, "grey96", "white")
)

base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.margin = margin(2, 6, 2, 6)
  )

chr_x_scale <- scale_x_continuous(
  limits = c(0, genome_len),
  expand = c(0, 0),
  breaks = chr_mids,
  labels = chr_order
)

as_segment_profile <- function(dt, value_col = "log2_absolute_diploid") {
  out <- copy(dt)
  if (!"chr" %in% names(out)) out[, chr := sub("^chr", "", as.character(chromosome))]
  out <- out[chr %in% chr_order]
  out[, genome_start := as.numeric(start) + chr_cumstart[chr]]
  out[, genome_end := as.numeric(end) + chr_cumstart[chr]]
  out[, value := as.numeric(get(value_col))]
  out[is.finite(genome_start) & is.finite(genome_end) & genome_end >= genome_start]
}

load_cnvkit <- function(wgs_id) {
  tumor_id <- sub("_vs_.*", "", wgs_id)
  candidates <- file.path(cnvkit_root, wgs_id, paste0(tumor_id, c(".somatic.call.cns", ".cns")))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) return(NULL)
  cns <- fread(candidates[1])
  cns[, chr := sub("^chr", "", as.character(chromosome))]
  cns <- cns[chr %in% chr_order]
  cns[, source_path := candidates[1]]
  cns
}

load_facets_absolute <- function(wgs_id) {
  seg_path <- file.path(wes_root, "tables/facets", paste0("Auto_", wgs_id, "_facets_segments.tsv"))
  pp_path <- file.path(wes_root, "tables/facets", paste0("Auto_", wgs_id, "_facets_purity_ploidy.tsv"))
  if (!file.exists(seg_path) || !file.exists(pp_path)) stop("Missing FACETS output for ", wgs_id)
  seg <- fread(seg_path)
  pp <- fread(pp_path)
  purity <- as.numeric(pp$purity[1])
  ploidy <- as.numeric(pp$ploidy[1])
  seg[, chr := sub("^chr", "", as.character(chrom))]
  seg <- seg[chr %in% chr_order]
  seg[, chromosome := paste0("chr", chr)]
  seg[, `:=`(
    start = as.integer(start),
    end = as.integer(end),
    cna_ccf = pmin(1, pmax(0, as.numeric(cf) / purity)),
    log2_absolute_diploid = log2(pmax(total_cn, 0.001) / 2),
    log2_relative_ploidy = log2(pmax(total_cn, 0.001) / ploidy),
    log2_expected_mixed_diploid = log2((purity * total_cn + (1 - purity) * 2) / 2),
    purity = purity,
    ploidy = ploidy
  )]
  seg[]
}

write_facets_cns <- function(wgs_id, facets, cnvkit) {
  facets_cns <- facets[, .(
    chromosome,
    start,
    end,
    gene = "-",
    log2 = round(log2_absolute_diploid, 6),
    depth = 0,
    probes = 1L,
    weight = 1,
    total_cn,
    major_cn,
    minor_cn,
    cna_ccf,
    purity,
    ploidy
  )]
  facets_path <- file.path(cns_dir, paste0("Auto_", wgs_id, "_facets_absolute_diploid.cns"))
  fwrite(facets_cns, facets_path, sep = "\t")

  if (is.null(cnvkit)) {
    return(list(facets_path = facets_path, cnvkit_grid_path = NA_character_, cnvkit_grid = NULL))
  }

  cnv_mid <- cnvkit[, .(
    chromosome,
    chr,
    start = as.integer((start + end) / 2),
    end = as.integer((start + end) / 2),
    cnvkit_row = .I
  )]
  fac_ov <- facets[, .(
    chr,
    start,
    end,
    total_cn,
    major_cn,
    minor_cn,
    cna_ccf,
    log2_absolute_diploid,
    purity,
    ploidy
  )]
  setkey(cnv_mid, chr, start, end)
  setkey(fac_ov, chr, start, end)
  mapped <- foverlaps(cnv_mid, fac_ov, nomatch = NA)[order(cnvkit_row)]

  grid <- copy(cnvkit)
  grid[, `:=`(
    total_cn = mapped$total_cn,
    major_cn = mapped$major_cn,
    minor_cn = mapped$minor_cn,
    cna_ccf = mapped$cna_ccf,
    log2_absolute_diploid = mapped$log2_absolute_diploid,
    purity = mapped$purity,
    ploidy = mapped$ploidy
  )]
  grid_cns <- grid[, .(
    chromosome,
    start = as.integer(start),
    end = as.integer(end),
    gene = if ("gene" %in% names(grid)) gene else "-",
    log2 = round(log2_absolute_diploid, 6),
    depth = if ("depth" %in% names(grid)) depth else 0,
    probes = if ("probes" %in% names(grid)) probes else 1L,
    weight = if ("weight" %in% names(grid)) weight else 1,
    total_cn,
    major_cn,
    minor_cn,
    cna_ccf,
    purity,
    ploidy
  )]
  grid_path <- file.path(cns_dir, paste0("Auto_", wgs_id, "_facets_absolute_diploid_on_cnvkit_grid.cns"))
  fwrite(grid_cns, grid_path, sep = "\t")
  list(facets_path = facets_path, cnvkit_grid_path = grid_path, cnvkit_grid = grid)
}

final_iter_from <- function(numbat_dir, prefix = "treeML", ext = "rds") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.", ext)))
  if (length(files) == 0) return(NA_integer_)
  regex <- paste0("^", prefix, "_([0-9]+)\\.", gsub("\\.", "\\\\.", ext), "$")
  iter <- suppressWarnings(as.integer(sub(regex, "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

load_numbat_native <- function(scrna_sample) {
  numbat_dir <- file.path(numbat_root, scrna_sample, "numbat")
  if (!dir.exists(numbat_dir)) return(NULL)

  segs_file <- file.path(numbat_dir, paste0("Auto_", scrna_sample, "_numbat_segs_consensus.csv"))
  if (!file.exists(segs_file)) {
    iter_seg <- final_iter_from(numbat_dir, "segs_consensus", "tsv")
    if (is.finite(iter_seg)) segs_file <- file.path(numbat_dir, paste0("segs_consensus_", iter_seg, ".tsv"))
  }
  bulk_segment <- NULL
  if (file.exists(segs_file)) {
    segs <- fread(segs_file)
    segs[, chr := as.character(CHROM)]
    segs <- segs[chr %in% chr_order]
    segs[, phi_val := suppressWarnings(as.numeric(phi_mle))]
    segs[, value := ifelse(is.finite(phi_val) & phi_val > 0, log2(phi_val), 0)]
    segs[, genome_start := as.numeric(seg_start) + chr_cumstart[chr]]
    segs[, genome_end := as.numeric(seg_end) + chr_cumstart[chr]]
    bulk_segment <- segs[is.finite(genome_start) & is.finite(genome_end), .(genome_start, genome_end, value, cnv_state_post)]
  }

  bulk_file <- file.path(numbat_dir, "bulk_clones_final.tsv.gz")
  if (!file.exists(bulk_file)) {
    iter <- final_iter_from(numbat_dir, "bulk_clones", "tsv.gz")
    if (is.finite(iter)) bulk_file <- file.path(numbat_dir, paste0("bulk_clones_", iter, ".tsv.gz"))
  }
  if (!file.exists(bulk_file)) return(list(bulk_segment = bulk_segment, clone_profiles = list(), clone_info = data.table()))

  nb <- fread(bulk_file, select = c("CHROM", "gene_start", "gene_end", "n_cells", "members", "phi_mle_roll", "gene_index"))
  nb[, chr := as.character(CHROM)]
  nb <- nb[chr %in% chr_order]
  nb[, genome_pos := (gene_start + gene_end) / 2 + chr_cumstart[chr]]
  # Native Numbat scale: phi_mle_roll is copy-number ratio relative to diploid.
  nb[, value := log2(as.numeric(phi_mle_roll))]
  nb <- nb[is.finite(genome_pos) & is.finite(value)]
  clone_info <- nb[nchar(members) > 0 & members != '""', .(n_cells = n_cells[1]), by = members][order(-n_cells)]
  if (nrow(clone_info) > 0) {
    clone_info <- clone_info[seq_len(min(5L, nrow(clone_info)))]
    clone_info[, clone_label := paste0("Numbat Clone ", seq_len(.N), " n=", n_cells)]
  } else {
    clone_info[, clone_label := character()]
  }
  clone_profiles <- lapply(seq_len(nrow(clone_info)), function(i) {
    d <- copy(nb[members == clone_info$members[i]])
    d[order(genome_pos), .(genome_pos, value)]
  })
  names(clone_profiles) <- clone_info$clone_label
  list(bulk_segment = bulk_segment, clone_profiles = clone_profiles, clone_info = clone_info)
}

make_bins <- function(bin_size) {
  bins <- data.table(genome_start = seq(0, genome_len - 1, by = bin_size))
  bins[, genome_end := pmin(genome_start + bin_size, genome_len)]
  bins[, genome_mid := (genome_start + genome_end) / 2]
  bins[, chr := {
    out <- rep(NA_character_, .N)
    for (ch in chr_order) out[genome_mid >= chr_cumstart[ch] & genome_mid < chr_cumend[ch]] <- ch
    out
  }]
  bins[!is.na(chr)]
}

bin_segment <- function(dt, bins) {
  sapply(seq_len(nrow(bins)), function(i) {
    hit <- dt[genome_start <= bins$genome_mid[i] & genome_end >= bins$genome_mid[i]]
    if (nrow(hit) == 0) NA_real_ else hit$value[1]
  })
}

bin_point <- function(dt, bins) {
  sapply(seq_len(nrow(bins)), function(i) {
    hit <- dt[genome_pos >= bins$genome_start[i] & genome_pos < bins$genome_end[i]]
    if (nrow(hit) == 0) NA_real_ else mean(hit$value, na.rm = TRUE)
  })
}

correlate_to_numbat <- function(wes_abs, numbat_data, bin_size = 5e6) {
  bins <- make_bins(bin_size)
  wes_vec <- bin_segment(wes_abs, bins)
  profiles <- list()
  if (!is.null(numbat_data$bulk_segment)) profiles[["Numbat pseudo-bulk native"]] <- numbat_data$bulk_segment
  profiles <- c(profiles, numbat_data$clone_profiles)
  if (length(profiles) == 0) return(data.table())
  rbindlist(lapply(names(profiles), function(nm) {
    profile <- profiles[[nm]]
    num_vec <- if (all(c("genome_start", "genome_end") %in% names(profile))) bin_segment(profile, bins) else bin_point(profile, bins)
    valid <- is.finite(wes_vec) & is.finite(num_vec)
    data.table(
      scrna_profile = nm,
      bin_size_mb = bin_size / 1e6,
      n_bins = sum(valid),
      correlation = if (sum(valid) >= 20) cor(wes_vec[valid], num_vec[valid]) else NA_real_
    )
  }))
}

make_ribbon_panel <- function(dt, title, is_point = FALSE, y_limits = c(-1.5, 2), show_chr_labels = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity()
  if (is_point) {
    p <- p +
      geom_ribbon(data = dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
      geom_ribbon(data = dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
      geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.25, color = "grey25")
  } else {
    p <- p +
      geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
      geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
      geom_segment(data = dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.28, color = "grey25")
  }
  p <- p +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = "log2 vs diploid") +
    base_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_heatmap <- function(cor_dt) {
  if (nrow(cor_dt) == 0) return(NULL)
  plot_dt <- copy(cor_dt[bin_size_mb == 5])
  if (nrow(plot_dt) == 0) return(NULL)
  plot_dt[, scrna_profile := factor(scrna_profile, levels = unique(scrna_profile))]
  ggplot(plot_dt, aes(x = scrna_profile, y = "WES FACETS absolute", fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey85") +
    labs(title = "5 Mb binned correlation", x = "scRNA profile", y = NULL, fill = "r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

make_plot <- function(wgs_id, scrna_sample, wes_abs, numbat_data, cor_dt) {
  panels <- list(
    wes = make_ribbon_panel(
      wes_abs,
      paste0("WES FACETS absolute CNA - ", wgs_id, " (log2 total CN / 2)")
    )
  )
  if (!is.null(numbat_data$bulk_segment)) {
    panels$numbat_bulk <- make_ribbon_panel(
      numbat_data$bulk_segment,
      paste0("Numbat pseudo-bulk native - ", scrna_sample, " (log2 phi)"),
      y_limits = c(-1.5, 2)
    )
  }
  if (length(numbat_data$clone_profiles) > 0) {
    clone_dt <- rbindlist(lapply(names(numbat_data$clone_profiles), function(nm) {
      d <- copy(numbat_data$clone_profiles[[nm]])
      d[, profile := nm]
      d
    }))
    clone_dt[, profile := factor(profile, levels = names(numbat_data$clone_profiles))]
    panels$numbat_clones <- ggplot() +
      geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
      scale_fill_identity() +
      geom_ribbon(data = clone_dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
      geom_ribbon(data = clone_dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
      geom_line(data = clone_dt, aes(x = genome_pos, y = value), linewidth = 0.2, color = "grey25") +
      geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
      facet_wrap(~profile, ncol = 1, strip.position = "right") +
      chr_x_scale +
      coord_cartesian(ylim = c(-1.5, 2)) +
      labs(title = "Numbat clone CNA profiles native scale (not centered)", x = NULL, y = "log2 phi") +
      base_theme +
      theme(strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"))
  }
  heatmap <- make_heatmap(cor_dt)
  if (!is.null(heatmap)) panels$heatmap <- heatmap

  heights <- vapply(names(panels), function(nm) {
    if (nm == "numbat_clones") max(1.5, length(numbat_data$clone_profiles) * 0.65)
    else if (nm == "heatmap") 2.2
    else 1
  }, numeric(1))
  Reduce(`/`, panels) + plot_layout(heights = heights) +
    plot_annotation(
      title = paste0(scrna_sample, " - Absolute WES CNA vs native Numbat CNA"),
      subtitle = "WES uses FACETS total copy number scaled to diploid; Numbat uses native log2(phi)."
    )
}

all_summary <- list()
all_cor <- list()

method_assessment <- data.table(
  conclusion = c(
    "PyClone-VI does not infer clone-specific CNA profiles from WES; it clusters SNVs using copy number as input.",
    "FACETS provides allele-specific integer CN and CNA cellular fraction per segment, sufficient for absolute bulk/event CNA visualization.",
    "True clone-specific WES CNA requires a clone-CNA deconvolution method such as THetA2/HATCHet/HATCHet2/CloneHD, not PyClone output alone.",
    "This script writes FACETS absolute diploid-scaled WES CNA profiles and native-scale Numbat comparisons."
  )
)
fwrite(method_assessment, file.path(table_dir, "Auto_wes_absolute_cna_method_assessment.csv"))

message("Starting FACETS absolute CNA comparison")
for (wgs_id in names(wgs_sample_map)) {
  message("Processing ", wgs_id)
  cnvkit <- load_cnvkit(wgs_id)
  facets <- load_facets_absolute(wgs_id)
  cns_paths <- write_facets_cns(wgs_id, facets, cnvkit)

  wes_abs <- if (!is.null(cns_paths$cnvkit_grid)) {
    as_segment_profile(cns_paths$cnvkit_grid, "log2_absolute_diploid")
  } else {
    as_segment_profile(facets, "log2_absolute_diploid")
  }

  summary_row <- facets[, .(
    wgs_id = wgs_id,
    purity = purity[1],
    ploidy = ploidy[1],
    facets_segments = .N,
    median_total_cn = median(total_cn, na.rm = TRUE),
    q10_total_cn = quantile(total_cn, 0.1, na.rm = TRUE),
    q90_total_cn = quantile(total_cn, 0.9, na.rm = TRUE),
    median_log2_absolute_diploid = median(log2_absolute_diploid, na.rm = TRUE),
    q10_log2_absolute_diploid = quantile(log2_absolute_diploid, 0.1, na.rm = TRUE),
    q90_log2_absolute_diploid = quantile(log2_absolute_diploid, 0.9, na.rm = TRUE),
    facets_cns = cns_paths$facets_path,
    cnvkit_grid_cns = cns_paths$cnvkit_grid_path
  )]
  all_summary[[wgs_id]] <- summary_row

  for (scrna_sample in wgs_sample_map[[wgs_id]]) {
    message("  Comparing ", scrna_sample)
    numbat_data <- load_numbat_native(scrna_sample)
    if (is.null(numbat_data)) next
    cor_5 <- correlate_to_numbat(wes_abs, numbat_data, 5e6)
    cor_1 <- correlate_to_numbat(wes_abs, numbat_data, 1e6)
    cor_dt <- rbindlist(list(cor_5, cor_1), use.names = TRUE, fill = TRUE)
    cor_dt[, `:=`(wgs_id = wgs_id, sample = scrna_sample)]
    all_cor[[scrna_sample]] <- cor_dt
    fwrite(cor_dt, file.path(table_dir, paste0("Auto_wes_absolute_cna_correlations_", scrna_sample, ".csv")))

    plot_obj <- make_plot(wgs_id, scrna_sample, wes_abs, numbat_data, cor_dt)
    height <- max(10, 4 + if (length(numbat_data$clone_profiles) > 0) length(numbat_data$clone_profiles) * 1.25 else 0)
    pdf_path <- file.path(fig_dir, paste0("Auto_wes_absolute_cna_compare_", scrna_sample, ".pdf"))
    png_path <- file.path(fig_dir, paste0("Auto_wes_absolute_cna_compare_", scrna_sample, ".png"))
    ggsave(pdf_path, plot_obj, width = 16, height = height, limitsize = FALSE)
    ggsave(png_path, plot_obj, width = 16, height = height, dpi = 200, limitsize = FALSE)
  }
}

summary_dt <- rbindlist(all_summary, use.names = TRUE, fill = TRUE)
cor_dt <- rbindlist(all_cor, use.names = TRUE, fill = TRUE)
fwrite(summary_dt, file.path(table_dir, "Auto_wes_absolute_cna_summary.csv"))
fwrite(cor_dt, file.path(table_dir, "Auto_wes_absolute_cna_correlations.csv"))
fwrite(
  data.table(
    finished = as.character(Sys.time()),
    n_wgs = length(all_summary),
    n_correlation_rows = nrow(cor_dt)
  ),
  file.path(log_dir, "Auto_wes_absolute_cna_compare_summary.tsv"),
  sep = "\t"
)

message("Done.")
