####################
# Auto_plot_phylowgs_numbat_compare.R
#
# Analysis registry
# Status: active terminal visualization
# Script: analysis/cnv/wes_subclone/Auto_plot_phylowgs_numbat_compare.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/Auto_phylowgs_top_tree_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_population_top_tree.csv
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_clone_cna_inherited_top_tree.csv
#   - PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_<sample>_bulk_highres_conditional_shift_absolute.cns when present
#   - ephemeral PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/bulk_clones_final.tsv.gz
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/figures/phylowgs/Auto_phylowgs_clone_cna_compare_<scrna_sample>.pdf/.png
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_visualisation_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_correlations.csv
#   - PDOs_outs/Auto_wes_subclone/logs/Auto_phylowgs_clone_cna_visualisation.tsv
# Downstream use: terminal visual audit of PhyloWGS final clone CNA profiles
#   against WES bulk and native high-resolution Numbat profiles.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

live_project <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
ephemeral_project <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
wes_root <- file.path(live_project, "PDOs_outs/Auto_wes_subclone")
absolute_root <- file.path(live_project, "PDOs_outs/Auto_wes_absolute_cna")
numbat_root <- file.path(ephemeral_project, "PDOs_outs/Auto_PDO_numbat/by_samples")

fig_dir <- file.path(wes_root, "figures/phylowgs")
table_dir <- file.path(wes_root, "tables/phylowgs_visualisation")
log_dir <- file.path(wes_root, "logs")
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

final_iter_from <- function(numbat_dir, prefix = "bulk_clones", ext = "tsv.gz") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.", ext)))
  if (length(files) == 0) return(NA_integer_)
  regex <- paste0("^", prefix, "_([0-9]+)\\.", gsub("\\.", "\\\\.", ext), "$")
  iter <- suppressWarnings(as.integer(sub(regex, "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

bulk_clones_path <- function(scrna_sample) {
  native_file <- file.path(numbat_root, scrna_sample, "numbat", "bulk_clones_final.tsv.gz")
  if (file.exists(native_file)) {
    return(list(path = native_file, scale = "native_phi"))
  }
  PDO_OUTPUT_DIR <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
  conservative_file <- file.path(
    PDO_OUTPUT_DIR,
    "Auto_PDO_numbat/conservative_clones/by_samples",
    scrna_sample,
    paste0("Auto_", scrna_sample, "_numbat_conservative_bulk_clones.csv.gz")
  )
  if (file.exists(conservative_file)) {
    return(list(path = conservative_file, scale = "log2_phi"))
  }
  NULL
}

load_numbat_native <- function(scrna_sample, max_clones = 6L) {
  bulk_source <- bulk_clones_path(scrna_sample)
  if (is.null(bulk_source) || !file.exists(bulk_source$path)) return(NULL)
  nb <- fread(bulk_source$path)
  keep_cols <- c("CHROM", "gene_start", "gene_end", "gene_index", "n_cells", "members", "phi_mle_roll")
  missing_cols <- setdiff(keep_cols, names(nb))
  if (length(missing_cols) > 0) stop("Numbat bulk clone file missing columns: ", paste(missing_cols, collapse = ", "))
  nb <- nb[, ..keep_cols]
  nb[, chr := sub("^chr", "", as.character(CHROM))]
  nb <- nb[chr %in% chr_order]
  nb[, `:=`(
    genome_pos = (as.numeric(gene_start) + as.numeric(gene_end)) / 2 + chr_cumstart[chr],
    phi_val = suppressWarnings(as.numeric(phi_mle_roll)),
    n_cells_num = suppressWarnings(as.numeric(n_cells)),
    members_chr = as.character(members)
  )]
  nb <- nb[is.finite(genome_pos) & is.finite(phi_val)]
  if (bulk_source$scale == "native_phi") nb <- nb[phi_val > 0]
  nb_clone <- nb[!is.na(members_chr) & nchar(members_chr) > 0 & members_chr != '""']
  if (nrow(nb_clone) == 0) return(NULL)

  bulk_point <- nb_clone[, .(
    phi_weighted = if (sum(n_cells_num, na.rm = TRUE) > 0) {
      sum(phi_val * n_cells_num, na.rm = TRUE) / sum(n_cells_num, na.rm = TRUE)
    } else {
      mean(phi_val, na.rm = TRUE)
    }
  ), by = .(chr, gene_index, genome_pos)]
  if (bulk_source$scale == "native_phi") {
    bulk_point[, value := log2(phi_weighted)]
  } else {
    bulk_point[, value := phi_weighted]
  }
  bulk_point <- bulk_point[is.finite(value)][order(genome_pos), .(genome_pos, value)]

  if (bulk_source$scale == "native_phi") {
    nb_clone[, value := log2(phi_val)]
  } else {
    nb_clone[, value := phi_val]
  }
  clone_info <- nb_clone[, .(n_cells = n_cells_num[1]), by = members_chr][order(-n_cells)]
  clone_info <- clone_info[seq_len(min(max_clones, nrow(clone_info)))]
  clone_info[, clone_label := paste0("Numbat Clone ", seq_len(.N), " n=", n_cells)]
  clone_profiles <- lapply(seq_len(nrow(clone_info)), function(i) {
    d <- copy(nb_clone[members_chr == clone_info$members_chr[i]])
    d[order(genome_pos), .(genome_pos, value)]
  })
  names(clone_profiles) <- clone_info$clone_label

  list(
    bulk_point = bulk_point,
    clone_profiles = clone_profiles,
    clone_info = clone_info,
    source = bulk_source$path,
    source_scale = bulk_source$scale
  )
}

as_segment_profile <- function(dt, chr_col = "chrom", start_col = "start", end_col = "end", value_col = "value") {
  out <- copy(dt)
  out[, chr := sub("^chr", "", as.character(get(chr_col)))]
  out <- out[chr %in% chr_order]
  out[, `:=`(
    genome_start = as.numeric(get(start_col)) + chr_cumstart[chr],
    genome_end = as.numeric(get(end_col)) + chr_cumstart[chr],
    value = as.numeric(get(value_col))
  )]
  out[is.finite(genome_start) & is.finite(genome_end) & is.finite(value) & genome_end >= genome_start]
}

load_wes_bulk <- function(wgs_id) {
  path <- file.path(
    absolute_root,
    "tables/cns",
    paste0("Auto_", wgs_id, "_bulk_highres_conditional_shift_absolute.cns")
  )
  if (!file.exists(path)) return(NULL)
  cns <- fread(path)
  if (!all(c("chromosome", "start", "end", "log2") %in% names(cns))) return(NULL)
  prof <- as_segment_profile(cns, "chromosome", "start", "end", "log2")
  shift_policy <- data.table(
    facets_ploidy_shift_log2 = if ("facets_ploidy_shift_log2" %in% names(cns)) {
      median(as.numeric(cns$facets_ploidy_shift_log2), na.rm = TRUE)
    } else {
      0
    },
    applied_shift_log2 = if ("applied_shift_log2" %in% names(cns)) {
      median(as.numeric(cns$applied_shift_log2), na.rm = TRUE)
    } else {
      0
    },
    shift_applied = if ("shift_applied" %in% names(cns)) {
      any(as.logical(cns$shift_applied))
    } else {
      FALSE
    }
  )
  list(
    profile = prof,
    shift_policy = shift_policy,
    summary = data.table(
      wgs_id = wgs_id,
      source_type = "conditional_shift_wes_bulk",
      display_label = "WES CNVkit high-res conditional-shift bulk",
      n_segments = nrow(prof),
      median_log2 = median(prof$value, na.rm = TRUE),
      mean_log2 = mean(prof$value, na.rm = TRUE),
      source = path
    )
  )
}

fill_neutral_segments <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(data.table(
      chr = character(),
      genome_start = numeric(),
      genome_end = numeric(),
      value = numeric()
    ))
  }
  event_dt <- copy(dt)[order(chr, genome_start, genome_end)]
  out <- list()
  for (ch in chr_order) {
    ch_start <- chr_cumstart[ch]
    ch_end <- chr_cumend[ch]
    d <- event_dt[chr == ch][order(genome_start, genome_end)]
    cursor <- ch_start
    if (nrow(d) > 0) {
      for (i in seq_len(nrow(d))) {
        seg_start <- max(ch_start, as.numeric(d$genome_start[i]))
        seg_end <- min(ch_end, as.numeric(d$genome_end[i]))
        if (!is.finite(seg_start) || !is.finite(seg_end) || seg_end < ch_start || seg_start > ch_end) next
        if (seg_start > cursor) {
          out[[length(out) + 1]] <- data.table(
            chr = ch,
            genome_start = cursor,
            genome_end = seg_start,
            value = 0,
            source_type = "phylowgs_neutral_fill"
          )
        }
        out[[length(out) + 1]] <- copy(d[i])[, source_type := "phylowgs_inherited_cna_event"]
        cursor <- max(cursor, seg_end)
      }
    }
    if (cursor < ch_end) {
      out[[length(out) + 1]] <- data.table(
        chr = ch,
        genome_start = cursor,
        genome_end = ch_end,
        value = 0,
        source_type = "phylowgs_neutral_fill"
      )
    }
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

load_phylowgs_clones <- function(wgs_id, shift_policy = NULL) {
  sample_dir <- file.path(wes_root, "tables/phylowgs", wgs_id)
  pop_path <- file.path(sample_dir, paste0("Auto_", wgs_id, "_phylowgs_population_top_tree.csv"))
  cna_path <- file.path(sample_dir, paste0("Auto_", wgs_id, "_phylowgs_clone_cna_inherited_top_tree.csv"))
  if (!file.exists(pop_path)) stop("Missing PhyloWGS population table: ", pop_path)
  if (!file.exists(cna_path)) stop("Missing PhyloWGS inherited clone-CNA table: ", cna_path)

  pop <- fread(pop_path)
  cna <- fread(cna_path)
  cna[, `:=`(
    total_cn_num = suppressWarnings(as.numeric(total_cn)),
    clone_population = as.integer(clone_population)
  )]
  cna <- cna[is.finite(total_cn_num) & total_cn_num > 0]
  if (is.null(shift_policy)) {
    facets_ploidy_shift <- median(log2(suppressWarnings(as.numeric(cna$ploidy)) / 2), na.rm = TRUE)
    if (!is.finite(facets_ploidy_shift)) facets_ploidy_shift <- 0
    applied_shift <- 0
  } else {
    facets_ploidy_shift <- as.numeric(shift_policy$facets_ploidy_shift_log2[1])
    applied_shift <- as.numeric(shift_policy$applied_shift_log2[1])
    if (!is.finite(facets_ploidy_shift)) facets_ploidy_shift <- 0
    if (!is.finite(applied_shift)) applied_shift <- 0
  }
  cna[, `:=`(
    absolute_log2_cn2 = log2(total_cn_num / 2),
    baseline_aligned_log2 = log2(total_cn_num / 2) - facets_ploidy_shift + applied_shift
  )]
  cna[, value := baseline_aligned_log2]

  event_pops <- pop[population != 0][order(population)]
  event_profiles <- list()
  final_profiles <- list()
  summary_rows <- list()
  segment_rows <- list()

  for (i in seq_len(nrow(event_pops))) {
    clone_pop <- as.integer(event_pops$population[i])
    d <- cna[clone_population == clone_pop]
    label <- paste0(
      "PhyloWGS Pop ", clone_pop,
      " CP=", sprintf("%.2f", as.numeric(event_pops$cellular_prevalence[i])),
      " CNV=", as.integer(event_pops$num_cnvs[i]),
      " SSM=", as.integer(event_pops$num_ssms[i])
    )
    event_profiles[[label]] <- as_segment_profile(d, "chrom", "start", "end", "value")
    final_profiles[[label]] <- fill_neutral_segments(event_profiles[[label]])
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      wgs_id = wgs_id,
      source_type = "phylowgs_final_clone_cna_neutral_filled",
      clone_population = clone_pop,
      display_label = label,
      cellular_prevalence = as.numeric(event_pops$cellular_prevalence[i]),
      num_ssms = as.integer(event_pops$num_ssms[i]),
      num_cnvs_direct = as.integer(event_pops$num_cnvs[i]),
      n_inherited_cna_segments = nrow(d),
      n_final_profile_segments = nrow(final_profiles[[label]]),
      n_unique_inherited_cna_events = uniqueN(d$cnv),
      median_total_cn = median(as.numeric(d$total_cn_num), na.rm = TRUE),
      median_absolute_log2_cn2 = median(as.numeric(d$absolute_log2_cn2), na.rm = TRUE),
      median_log2_cn2 = median(as.numeric(final_profiles[[label]]$value), na.rm = TRUE),
      mean_log2_cn2 = mean(as.numeric(final_profiles[[label]]$value), na.rm = TRUE),
      facets_ploidy_shift_log2 = facets_ploidy_shift,
      applied_shift_log2 = applied_shift,
      source = cna_path
    )
    if (nrow(final_profiles[[label]]) > 0) {
      seg_dt <- cbind(
        data.table(wgs_id = wgs_id, clone_population = clone_pop, display_label = label),
        final_profiles[[label]]
      )
      seg_dt[, source_type := ifelse(is.na(source_type), "phylowgs_final_clone_cna_neutral_filled", source_type)]
      segment_rows[[length(segment_rows) + 1]] <- seg_dt
    }
  }

  list(
    event_profiles = event_profiles,
    final_profiles = final_profiles,
    summary = rbindlist(summary_rows, use.names = TRUE, fill = TRUE),
    segments = rbindlist(segment_rows, use.names = TRUE, fill = TRUE),
    populations = pop
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

bin_segment <- function(dt, bins, missing_value = 0) {
  sapply(seq_len(nrow(bins)), function(i) {
    hit <- dt[genome_start <= bins$genome_mid[i] & genome_end >= bins$genome_mid[i]]
    if (nrow(hit) == 0) missing_value else hit$value[1]
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
  wes_binned <- lapply(wes_profiles, bin_segment, bins = bins, missing_value = 0)
  scrna_binned <- lapply(scrna_profiles, function(x) {
    if (all(c("genome_start", "genome_end") %in% names(x))) {
      bin_segment(x, bins, missing_value = NA_real_)
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
    geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.16, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_facet_segments <- function(profile_list, title, y_lab = "log2(CN/2)", y_limits = c(-2.2, 2.8), overlay_list = NULL) {
  if (length(profile_list) == 0) return(NULL)
  dt <- rbindlist(lapply(names(profile_list), function(nm) {
    d <- copy(profile_list[[nm]])
    d[, profile := nm]
    d
  }), use.names = TRUE, fill = TRUE)
  dt[, profile := factor(profile, levels = names(profile_list))]
  overlay_dt <- NULL
  if (!is.null(overlay_list) && length(overlay_list) > 0) {
    overlay_dt <- rbindlist(lapply(names(overlay_list), function(nm) {
      d <- copy(overlay_list[[nm]])
      d[, profile := nm]
      d
    }), use.names = TRUE, fill = TRUE)
    overlay_dt <- overlay_dt[profile %in% names(profile_list)]
    if (nrow(overlay_dt) > 0) overlay_dt[, profile := factor(profile, levels = names(profile_list))]
  }
  p <- ggplot() +
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
  if (!is.null(overlay_dt) && nrow(overlay_dt) > 0) {
    p <- p +
      geom_segment(
        data = overlay_dt,
        aes(x = genome_start, xend = genome_end, y = value, yend = value),
        linewidth = 0.34,
        color = "black"
      ) +
      geom_point(
        data = overlay_dt[, .(genome_start = genome_start[1], value = value[1]), by = profile],
        aes(x = genome_start, y = value),
        size = 0.8,
        color = "black"
      )
  }
  p
}

make_facet_points <- function(profile_list, title, y_lab = "log2 phi", y_limits = c(-2.2, 2.8)) {
  if (length(profile_list) == 0) return(NULL)
  dt <- rbindlist(lapply(names(profile_list), function(nm) {
    d <- copy(profile_list[[nm]])
    d[, profile := nm]
    d
  }), use.names = TRUE, fill = TRUE)
  dt[, profile := factor(profile, levels = names(profile_list))]
  ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.14, color = "grey25") +
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
    labs(title = "5 Mb binned Pearson correlations", x = "scRNA profile", y = "WES PhyloWGS profile", fill = "r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

make_plot <- function(wgs_id, scrna_sample, phylowgs_profiles, phylowgs_event_profiles, numbat_data, cor_dt, wes_bulk = NULL) {
  panels <- list()
  if (!is.null(wes_bulk)) {
    panels$wes_bulk <- make_segment_panel(
      wes_bulk,
      paste0("WES CNVkit high-res conditional-shift bulk - ", wgs_id),
      y_lab = "log2 vs diploid"
    )
  }
  panels$phylowgs_clones <- make_facet_segments(
    phylowgs_profiles,
    paste0("PhyloWGS final clone CNA profiles - ", wgs_id),
    "baseline-aligned log2 CNA",
    overlay_list = NULL
  )
  panels$numbat_bulk <- make_point_panel(
    numbat_data$bulk_point,
    paste0("Numbat pseudo-bulk native high-res - ", scrna_sample),
    "log2 phi"
  )
  panels$numbat_clones <- make_facet_points(
    numbat_data$clone_profiles,
    paste0("Numbat clone CNA profiles native scale - ", scrna_sample),
    "log2 phi"
  )
  heatmap <- make_heatmap(cor_dt)
  if (!is.null(heatmap)) panels$heatmap <- heatmap

  heights <- c(
    wes_bulk = 1.0,
    phylowgs_clones = max(1.8, length(phylowgs_profiles) * 0.72),
    numbat_bulk = 1.0,
    numbat_clones = max(1.8, length(numbat_data$clone_profiles) * 0.55),
    heatmap = 2.5
  )[names(panels)]

  Reduce(`/`, panels) + plot_layout(heights = heights) +
    plot_annotation(
      title = paste0(scrna_sample, " - PhyloWGS WES clone-CNA vs Numbat"),
      subtitle = "PhyloWGS rows show inherited CNA events with neutral-filled non-event intervals on the conditional WES baseline; Numbat is native high-resolution phi."
    )
}

start_time <- Sys.time()
message("Starting PhyloWGS/Numbat CNA visualisation")
all_summary <- list()
all_cor <- list()
all_segments <- list()
written_figures <- character()

top_summary_path <- file.path(wes_root, "tables/phylowgs/Auto_phylowgs_top_tree_summary.csv")
if (!file.exists(top_summary_path)) stop("Missing PhyloWGS top-tree summary: ", top_summary_path)
top_summary <- fread(top_summary_path)

for (wgs_id in names(wgs_sample_map)) {
  message("Processing ", wgs_id)
  wes_bulk_obj <- load_wes_bulk(wgs_id)
  if (is.null(wes_bulk_obj)) stop("Missing conditional-shift WES bulk reference for ", wgs_id)
  phylowgs <- load_phylowgs_clones(wgs_id, wes_bulk_obj$shift_policy)
  all_summary[[paste0(wgs_id, "_phylowgs")]] <- phylowgs$summary
  all_segments[[paste0(wgs_id, "_events")]] <- phylowgs$segments

  if (!is.null(wes_bulk_obj)) {
    all_summary[[paste0(wgs_id, "_bulk")]] <- wes_bulk_obj$summary
  }

  for (scrna_sample in wgs_sample_map[[wgs_id]]) {
    message("  Comparing ", scrna_sample)
    numbat_data <- load_numbat_native(scrna_sample)
    if (is.null(numbat_data)) {
      warning("Skipping ", scrna_sample, ": missing usable Numbat native bulk_clones output")
      next
    }

    scrna_profiles <- c(
      setNames(list(numbat_data$bulk_point), "Numbat pseudo-bulk"),
      numbat_data$clone_profiles
    )
    cor_dt <- rbindlist(lapply(c(5e6, 10e6), function(bin_size) {
      out <- correlate_profiles(phylowgs$final_profiles, scrna_profiles, bin_size)
      out[, `:=`(wgs_id = wgs_id, scrna_sample = scrna_sample)]
      out
    }))
    all_cor[[paste(wgs_id, scrna_sample, sep = "__")]] <- cor_dt

    plot_obj <- make_plot(
      wgs_id = wgs_id,
      scrna_sample = scrna_sample,
      phylowgs_profiles = phylowgs$final_profiles,
      phylowgs_event_profiles = phylowgs$event_profiles,
      numbat_data = numbat_data,
      cor_dt = cor_dt,
      wes_bulk = wes_bulk_obj$profile
    )
    safe_sample <- gsub("[^A-Za-z0-9]+", "_", scrna_sample)
    pdf_path <- file.path(fig_dir, paste0("Auto_phylowgs_clone_cna_compare_", safe_sample, ".pdf"))
    png_path <- file.path(fig_dir, paste0("Auto_phylowgs_clone_cna_compare_", safe_sample, ".png"))
    height <- 4.5 + length(phylowgs$final_profiles) * 0.72 + length(numbat_data$clone_profiles) * 0.55 + 2.5
    ggsave(pdf_path, plot_obj, width = 16, height = height, limitsize = FALSE)
    ggsave(png_path, plot_obj, width = 16, height = height, dpi = 200, limitsize = FALSE)
    written_figures <- c(written_figures, pdf_path, png_path)
  }
}

summary_dt <- rbindlist(all_summary, use.names = TRUE, fill = TRUE)
cor_dt <- rbindlist(all_cor, use.names = TRUE, fill = TRUE)
segment_dt <- rbindlist(all_segments, use.names = TRUE, fill = TRUE)

summary_out <- file.path(table_dir, "Auto_phylowgs_clone_cna_visualisation_summary.csv")
cor_out <- file.path(table_dir, "Auto_phylowgs_clone_cna_correlations.csv")
segment_out <- file.path(table_dir, "Auto_phylowgs_clone_cna_segments_for_plot.csv")
fwrite(summary_dt, summary_out)
fwrite(cor_dt, cor_out)
fwrite(segment_dt, segment_out)

end_time <- Sys.time()
log_dt <- data.table(
  started = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
  finished = format(end_time, "%Y-%m-%d %H:%M:%S %Z"),
  status = "success",
  top_summary = top_summary_path,
  n_top_summary_rows = nrow(top_summary),
  n_summary_rows = nrow(summary_dt),
  n_correlation_rows = nrow(cor_dt),
  n_segment_rows = nrow(segment_dt),
  n_figures = length(written_figures),
  figure_dir = fig_dir,
  table_dir = table_dir,
  figures = paste(written_figures, collapse = ";")
)
fwrite(log_dt, file.path(log_dir, "Auto_phylowgs_clone_cna_visualisation.tsv"), sep = "\t")

message("Finished PhyloWGS/Numbat CNA visualisation")
