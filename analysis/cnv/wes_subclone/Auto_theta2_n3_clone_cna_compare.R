####################
# Auto_theta2_n3_clone_cna_compare.R
#
# Analysis registry
# Status: active terminal diagnostic
# Script: analysis/cnv/wes_subclone/Auto_theta2_n3_clone_cna_compare.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_clone_cna/tables/Auto_theta2_n3_clone_cna_manifest.csv
#   - PDOs_outs/Auto_wes_clone_cna/tables/cns_theta2_n3/Auto_<sample>_theta2_n3_clone*.cns
#   - PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_<sample>_bulk_highres_conditional_shift_absolute.cns
#   - ephemeral PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat bulk_clones outputs
# Outputs:
#   - PDOs_outs/Auto_wes_clone_cna/figures/Auto_theta2_n3_clone_cna_compare_<sample>.pdf/.png
#   - PDOs_outs/Auto_wes_clone_cna/tables/Auto_theta2_n3_clone_cna_summary.csv
#   - PDOs_outs/Auto_wes_clone_cna/tables/Auto_theta2_n3_clone_cna_correlations.csv
#   - PDOs_outs/Auto_wes_clone_cna/logs/Auto_theta2_n3_clone_cna_compare_summary.tsv
# Downstream use: terminal assessment of forced THetA2 n=3 WES clone-specific CNA profiles.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

live_project <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
ephemeral_project <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
clone_root <- file.path(live_project, "PDOs_outs/Auto_wes_clone_cna")
absolute_root <- file.path(live_project, "PDOs_outs/Auto_wes_absolute_cna")
numbat_root <- file.path(ephemeral_project, "PDOs_outs/Auto_PDO_numbat/by_samples")

fig_dir <- file.path(clone_root, "figures")
table_dir <- file.path(clone_root, "tables")
log_dir <- file.path(clone_root, "logs")
for (d in c(fig_dir, table_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

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

final_iter_from <- function(numbat_dir, prefix = "treeML", ext = "rds") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.", ext)))
  if (length(files) == 0) return(NA_integer_)
  regex <- paste0("^", prefix, "_([0-9]+)\\.", gsub("\\.", "\\\\.", ext), "$")
  iter <- suppressWarnings(as.integer(sub(regex, "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

bulk_clones_path <- function(scrna_sample) {
  numbat_dir <- file.path(numbat_root, scrna_sample, "numbat")
  if (!dir.exists(numbat_dir)) return(NA_character_)
  bulk_file <- file.path(numbat_dir, "bulk_clones_final.tsv.gz")
  if (!file.exists(bulk_file)) {
    iter <- final_iter_from(numbat_dir, "bulk_clones", "tsv.gz")
    if (is.finite(iter)) bulk_file <- file.path(numbat_dir, paste0("bulk_clones_", iter, ".tsv.gz"))
  }
  if (file.exists(bulk_file)) bulk_file else NA_character_
}

read_numbat_bulk_clones <- function(scrna_sample) {
  bulk_file <- bulk_clones_path(scrna_sample)
  if (is.na(bulk_file) || !file.exists(bulk_file)) return(NULL)
  fread(bulk_file, select = c("CHROM", "gene_start", "gene_end", "n_cells", "members", "phi_mle_roll", "gene_index"))
}

load_numbat_native <- function(scrna_sample) {
  nb <- read_numbat_bulk_clones(scrna_sample)
  if (is.null(nb) || nrow(nb) == 0) return(NULL)

  nb[, chr := as.character(CHROM)]
  nb <- nb[chr %in% chr_order]
  nb[, `:=`(
    genome_pos = (gene_start + gene_end) / 2 + chr_cumstart[chr],
    phi_val = suppressWarnings(as.numeric(phi_mle_roll)),
    n_cells_num = suppressWarnings(as.numeric(n_cells))
  )]
  nb <- nb[is.finite(genome_pos) & is.finite(phi_val) & phi_val > 0]
  nb_clone <- nb[nchar(members) > 0 & members != '""']
  if (nrow(nb_clone) == 0) return(NULL)

  bulk_point <- nb_clone[, .(
    phi_weighted = if (sum(n_cells_num, na.rm = TRUE) > 0) {
      sum(phi_val * n_cells_num, na.rm = TRUE) / sum(n_cells_num, na.rm = TRUE)
    } else {
      mean(phi_val, na.rm = TRUE)
    }
  ), by = .(chr, gene_index, genome_pos)]
  bulk_point[, value := log2(phi_weighted)]
  bulk_point <- bulk_point[is.finite(value)][order(genome_pos), .(genome_pos, value)]

  nb_clone[, value := log2(phi_val)]
  clone_info <- nb_clone[, .(n_cells = n_cells_num[1]), by = members][order(-n_cells)]
  clone_info <- clone_info[seq_len(min(6L, nrow(clone_info)))]
  clone_info[, clone_label := paste0("Numbat Clone ", seq_len(.N), " n=", n_cells)]
  clone_profiles <- lapply(seq_len(nrow(clone_info)), function(i) {
    d <- copy(nb_clone[members == clone_info$members[i]])
    d[order(genome_pos), .(genome_pos, value)]
  })
  names(clone_profiles) <- clone_info$clone_label

  list(bulk_point = bulk_point, clone_profiles = clone_profiles, clone_info = clone_info)
}

as_segment_profile <- function(dt, value_col = "log2") {
  out <- copy(dt)
  out[, chr := sub("^chr", "", as.character(chromosome))]
  out <- out[chr %in% chr_order]
  out[, genome_start := as.numeric(start) + chr_cumstart[chr]]
  out[, genome_end := as.numeric(end) + chr_cumstart[chr]]
  out[, value := as.numeric(get(value_col))]
  out[is.finite(genome_start) & is.finite(genome_end) & is.finite(value)]
}

parse_theta_meta <- function(wgs_id) {
  theta_file <- file.path(table_dir, "theta2", paste0("Auto_", wgs_id, ".n3.results"))
  if (!file.exists(theta_file)) return(data.table())
  theta <- fread(theta_file)
  setnames(theta, sub("^#", "", names(theta)))
  mu <- as.numeric(strsplit(theta$mu[1], ",", fixed = TRUE)[[1]])
  tumor_mu <- mu[-1]
  data.table(
    wgs_id = wgs_id,
    theta_nll = as.numeric(theta$NLL[1]),
    normal_fraction = mu[1],
    clone_index = seq_along(tumor_mu),
    theta_clone_fraction = tumor_mu
  )
}

load_theta2_n3_clones <- function(wgs_id) {
  manifest_path <- file.path(table_dir, "Auto_theta2_n3_clone_cna_manifest.csv")
  if (!file.exists(manifest_path)) stop("Missing THetA2 n3 manifest: ", manifest_path)
  manifest <- fread(manifest_path)[sample == wgs_id]
  if (nrow(manifest) == 0) stop("No THetA2 n3 imported clones for ", wgs_id)
  theta_meta <- parse_theta_meta(wgs_id)
  manifest[, clone_index := as.integer(clone_index)]
  manifest <- merge(manifest, theta_meta, by = "clone_index", all.x = TRUE, sort = FALSE)
  setorder(manifest, clone_index)

  profiles <- list()
  summary_rows <- list()
  for (i in seq_len(nrow(manifest))) {
    path <- manifest$output_cns[i]
    if (!file.exists(path)) stop("Missing imported THetA2 n3 CNS: ", path)
    cns <- fread(path)
    cns[, cn := suppressWarnings(as.numeric(cn))]
    cns[, log2 := suppressWarnings(as.numeric(log2))]
    label <- paste0(
      "THetA2 n3 Clone ", manifest$clone_index[i],
      " frac=", ifelse(is.finite(manifest$theta_clone_fraction[i]), sprintf("%.2f", manifest$theta_clone_fraction[i]), "NA")
    )
    profiles[[label]] <- as_segment_profile(cns, "log2")
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      wgs_id = wgs_id,
      clone_index = manifest$clone_index[i],
      display_label = label,
      theta_nll = manifest$theta_nll[i],
      normal_fraction = manifest$normal_fraction[i],
      theta_clone_fraction = manifest$theta_clone_fraction[i],
      n_segments = nrow(cns),
      median_log2 = median(cns$log2, na.rm = TRUE),
      mean_log2 = mean(cns$log2, na.rm = TRUE),
      n_positive = sum(cns$log2 > 0.01, na.rm = TRUE),
      n_neutral_cn2 = sum(cns$cn == 2, na.rm = TRUE),
      n_gain_cn_ge3 = sum(cns$cn >= 3, na.rm = TRUE),
      n_loss_cn_le1 = sum(cns$cn <= 1, na.rm = TRUE),
      output_cns = path,
      source_theta_results = manifest$source_theta_results[i]
    )
  }
  list(profiles = profiles, summary = rbindlist(summary_rows, use.names = TRUE, fill = TRUE))
}

load_wes_conditional_bulk <- function(wgs_id) {
  path <- file.path(
    absolute_root,
    "tables/cns",
    paste0("Auto_", wgs_id, "_bulk_highres_conditional_shift_absolute.cns")
  )
  if (!file.exists(path)) stop("Missing conditional-shift bulk CNS: ", path)
  cns <- fread(path)
  list(
    profile = as_segment_profile(cns, "log2"),
    summary = data.table(
      wgs_id = wgs_id,
      n_segments = nrow(cns),
      median_log2 = median(as.numeric(cns$log2), na.rm = TRUE),
      applied_shift_log2 = if ("applied_shift_log2" %in% names(cns)) unique(cns$applied_shift_log2)[1] else NA_real_,
      source_cns = path
    )
  )
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

correlate_profiles <- function(wes_profiles, scrna_profiles, bin_size) {
  bins <- make_bins(bin_size)
  wes_binned <- lapply(wes_profiles, bin_segment, bins = bins)
  scrna_binned <- lapply(scrna_profiles, function(x) {
    if (all(c("genome_start", "genome_end") %in% names(x))) {
      bin_segment(x, bins)
    } else {
      bin_point(x, bins)
    }
  })
  rows <- list()
  for (wn in names(wes_binned)) {
    for (sn in names(scrna_binned)) {
      w <- wes_binned[[wn]]
      s <- scrna_binned[[sn]]
      valid <- is.finite(w) & is.finite(s)
      rows[[length(rows) + 1]] <- data.table(
        wes_profile = wn,
        scrna_profile = sn,
        bin_size_mb = bin_size / 1e6,
        n_bins = sum(valid),
        correlation = if (sum(valid) >= 20) cor(w[valid], s[valid]) else NA_real_
      )
    }
  }
  rbindlist(rows)
}

make_segment_panel <- function(dt, title, y_lab = "log2", y_limits = c(-2.2, 2.8), show_chr_labels = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_segment(data = dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.22, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_point_panel <- function(dt, title, y_lab = "log2", y_limits = c(-2.2, 2.8), show_chr_labels = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.18, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_facet_segments <- function(profile_list, title, y_lab = "log2", y_limits = c(-2.2, 2.8)) {
  if (length(profile_list) == 0) return(NULL)
  dt <- rbindlist(lapply(names(profile_list), function(nm) {
    d <- copy(profile_list[[nm]])
    d[, profile := nm]
    d
  }))
  dt[, profile := factor(profile, levels = names(profile_list))]
  ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_segment(data = dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.16, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
    facet_wrap(~profile, ncol = 1, strip.position = "right") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme +
    theme(strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"))
}

make_facet_points <- function(profile_list, title, y_lab = "log2", y_limits = c(-2.2, 2.8)) {
  if (length(profile_list) == 0) return(NULL)
  dt <- rbindlist(lapply(names(profile_list), function(nm) {
    d <- copy(profile_list[[nm]])
    d[, profile := nm]
    d
  }))
  dt[, profile := factor(profile, levels = names(profile_list))]
  ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.15, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
    facet_wrap(~profile, ncol = 1, strip.position = "right") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme +
    theme(strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"))
}

make_heatmap <- function(cor_dt) {
  plot_dt <- copy(cor_dt[bin_size_mb == 5])
  if (nrow(plot_dt) == 0) return(NULL)
  plot_dt[, wes_profile := factor(wes_profile, levels = unique(wes_profile))]
  plot_dt[, scrna_profile := factor(scrna_profile, levels = unique(scrna_profile))]
  ggplot(plot_dt, aes(x = scrna_profile, y = wes_profile, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = ifelse(is.finite(correlation), sprintf("%.2f", correlation), "NA")), size = 2.6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey85") +
    labs(title = "5 Mb binned Pearson correlations", x = "scRNA profile", y = "WES profile", fill = "r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

make_plot <- function(wgs_id, scrna_sample, wes_bulk, theta_profiles, numbat_data, cor_dt) {
  panels <- list(
    wes_bulk = make_segment_panel(
      wes_bulk,
      paste0("WES CNVkit high-res conditional-shift bulk - ", wgs_id),
      y_lab = "log2 vs diploid",
      y_limits = c(-2.2, 2.8)
    ),
    theta_clones = make_facet_segments(
      theta_profiles,
      paste0("THetA2 forced n=3 WES clone-specific CNA - ", wgs_id),
      "log2(CN/2)",
      c(-2.2, 2.8)
    ),
    numbat_bulk = make_point_panel(
      numbat_data$bulk_point,
      paste0("Numbat pseudo-bulk native high-res - ", scrna_sample),
      "log2 phi",
      c(-2.2, 2.8)
    ),
    numbat_clones = make_facet_points(
      numbat_data$clone_profiles,
      paste0("Numbat clone CNA profiles native scale - ", scrna_sample),
      "log2 phi",
      c(-2.2, 2.8)
    )
  )
  heatmap <- make_heatmap(cor_dt)
  if (!is.null(heatmap)) panels$heatmap <- heatmap

  heights <- c(
    wes_bulk = 1.0,
    theta_clones = max(1.4, length(theta_profiles) * 0.8),
    numbat_bulk = 1.0,
    numbat_clones = max(1.8, length(numbat_data$clone_profiles) * 0.55),
    heatmap = 2.5
  )[names(panels)]

  Reduce(`/`, panels) + plot_layout(heights = heights) +
    plot_annotation(
      title = paste0(scrna_sample, " - THetA2 n3 WES clone-CNA vs scRNA CNA"),
      subtitle = "THetA2 n=3 is a forced higher-clone solution; compare against model-selected BEST and native Numbat before biological interpretation."
    )
}

message("Starting THetA2 n3 clone-CNA comparison")
all_theta_summary <- list()
all_bulk_summary <- list()
all_cor <- list()

for (wgs_id in names(wgs_sample_map)) {
  message("Processing ", wgs_id)
  theta <- load_theta2_n3_clones(wgs_id)
  wes_bulk <- load_wes_conditional_bulk(wgs_id)
  all_theta_summary[[wgs_id]] <- theta$summary
  all_bulk_summary[[wgs_id]] <- wes_bulk$summary

  for (scrna_sample in wgs_sample_map[[wgs_id]]) {
    message("  Comparing ", scrna_sample)
    numbat_data <- load_numbat_native(scrna_sample)
    if (is.null(numbat_data)) {
      warning("Skipping ", scrna_sample, ": missing usable native Numbat bulk_clones output")
      next
    }

    wes_profiles <- c(list("WES CNVkit conditional-shift bulk" = wes_bulk$profile), theta$profiles)
    scrna_profiles <- c(list("Numbat pseudo-bulk native high-res" = numbat_data$bulk_point), numbat_data$clone_profiles)
    corr_5 <- correlate_profiles(wes_profiles, scrna_profiles, 5e6)
    corr_1 <- correlate_profiles(wes_profiles, scrna_profiles, 1e6)
    cor_dt <- rbindlist(list(corr_5, corr_1), use.names = TRUE, fill = TRUE)
    cor_dt[, `:=`(wgs_id = wgs_id, sample = scrna_sample)]
    all_cor[[scrna_sample]] <- cor_dt

    plot_obj <- make_plot(wgs_id, scrna_sample, wes_bulk$profile, theta$profiles, numbat_data, cor_dt)
    height <- max(12, 5.0 + 0.9 * length(theta$profiles) + 0.7 * length(numbat_data$clone_profiles))
    pdf_path <- file.path(fig_dir, paste0("Auto_theta2_n3_clone_cna_compare_", scrna_sample, ".pdf"))
    png_path <- file.path(fig_dir, paste0("Auto_theta2_n3_clone_cna_compare_", scrna_sample, ".png"))
    ggsave(pdf_path, plot_obj, width = 16, height = height, limitsize = FALSE)
    ggsave(png_path, plot_obj, width = 16, height = height, dpi = 200, limitsize = FALSE)
  }
}

theta_summary <- rbindlist(all_theta_summary, use.names = TRUE, fill = TRUE)
bulk_summary <- rbindlist(all_bulk_summary, use.names = TRUE, fill = TRUE)
cor_dt <- rbindlist(all_cor, use.names = TRUE, fill = TRUE)

fwrite(theta_summary, file.path(table_dir, "Auto_theta2_n3_clone_cna_summary.csv"))
fwrite(bulk_summary, file.path(table_dir, "Auto_theta2_n3_bulk_reference_summary.csv"))
fwrite(cor_dt, file.path(table_dir, "Auto_theta2_n3_clone_cna_correlations.csv"))

fwrite(
  data.table(
    finished = as.character(Sys.time()),
    n_theta_clone_rows = nrow(theta_summary),
    n_correlation_rows = nrow(cor_dt),
    output_fig_dir = fig_dir,
    output_table_dir = table_dir
  ),
  file.path(log_dir, "Auto_theta2_n3_clone_cna_compare_summary.tsv"),
  sep = "\t"
)

message("Done.")
