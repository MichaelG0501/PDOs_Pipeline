####################
# Auto_wes_absolute_highres_subclone_compare.R
#
# Analysis registry
# Status: active terminal diagnostic
# Script: analysis/cnv/wes_subclone/Auto_wes_absolute_highres_subclone_compare.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/cns_highres/Auto_<sample>_*.cns
#   - PDOs_outs/Auto_wes_subclone/tables/cns_highres/Auto_wes_subclone_cns_highres_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_purity_ploidy.tsv
#   - ephemeral PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat outputs
#   - ephemeral InferCNA outputs when available
# Outputs:
#   - PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_<sample>_*_conditional_shift_absolute.cns
#   - PDOs_outs/Auto_wes_absolute_cna/tables/Auto_wes_absolute_cna_*.csv
#   - PDOs_outs/Auto_wes_absolute_cna/figures/Auto_wes_absolute_cna_compare_<sample>.pdf/.png
# Downstream use: terminal WES/scRNA CNA comparison. CNVkit high-resolution
#   profiles are shown both as centered shape tracks and conditional-shift
#   tracks. The FACETS ploidy offset is applied only when native Numbat indicates
#   a global amplified baseline.
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
numbat_root <- file.path(ephemeral_project, "PDOs_outs/Auto_PDO_numbat/by_samples")
infercna_outs_path <- file.path(ephemeral_project, "PDOs_outs/cnv/Auto_PDO_infercna_outs_Carroll_2023.rds")
infercna_meta_path <- file.path(ephemeral_project, "PDOs_outs/cnv/Auto_PDO_infercna_meta_Carroll_2023.csv")
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"

fig_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables")
cns_dir <- file.path(table_dir, "cns")
log_dir <- file.path(out_root, "logs")
for (d in c(fig_dir, table_dir, cns_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

wgs_sample_map <- list(
  "PDO_1090_vs_NT_1090" = c("SUR1090_Untreated_PDO", "SUR1090_Treated_PDO"),
  "PDO_1181_vs_NT_1181" = c("SUR1181_Untreated_PDO", "SUR1181_Treated_PDO")
)

####################
# Apply the global ploidy offset only when native Numbat indicates a real
# sample-wide amplification baseline. This keeps SUR1090 on the CNVkit-centered
# shape scale while retaining the SUR1181 polyploid shift requested for display.
####################
numbat_shift_trigger_log2 <- 0.25

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

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

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
  if (!is.na(bulk_file) && file.exists(bulk_file)) {
    fread(bulk_file, select = c("CHROM", "gene_start", "gene_end", "n_cells", "members", "phi_mle_roll", "gene_index"))
  } else {
    NULL
  }
}

numbat_native_median <- function(scrna_sample) {
  nb <- read_numbat_bulk_clones(scrna_sample)
  if (is.null(nb) || nrow(nb) == 0) return(NA_real_)
  nb[, phi_val := suppressWarnings(as.numeric(phi_mle_roll))]
  nb <- nb[nchar(members) > 0 & members != '""' & is.finite(phi_val) & phi_val > 0]
  if (nrow(nb) == 0) return(NA_real_)
  median(log2(nb$phi_val), na.rm = TRUE)
}

choose_shift_policy <- function(wgs_id) {
  pp <- load_ploidy(wgs_id)
  sample_medians <- data.table(
    sample = wgs_sample_map[[wgs_id]],
    numbat_native_median_log2 = vapply(wgs_sample_map[[wgs_id]], numbat_native_median, numeric(1))
  )
  median_numbat <- median(sample_medians$numbat_native_median_log2, na.rm = TRUE)
  if (!is.finite(median_numbat)) median_numbat <- NA_real_
  apply_shift <- is.finite(median_numbat) &&
    median_numbat >= numbat_shift_trigger_log2 &&
    pp$ploidy_shift_log2[1] > 0
  data.table(
    wgs_id = wgs_id,
    facets_purity = pp$purity[1],
    facets_ploidy = pp$ploidy[1],
    facets_ploidy_shift_log2 = pp$ploidy_shift_log2[1],
    numbat_native_median_log2 = median_numbat,
    applied_shift_log2 = if (apply_shift) pp$ploidy_shift_log2[1] else 0,
    shift_applied = apply_shift,
    shift_reason = if (apply_shift) {
      paste0("native Numbat median >= ", numbat_shift_trigger_log2, "; apply FACETS ploidy offset")
    } else {
      paste0("native Numbat median < ", numbat_shift_trigger_log2, " or unavailable; keep CNVkit centered baseline")
    }
  )
}

as_segment_profile <- function(dt, value_col = "log2") {
  out <- copy(dt)
  out[, chr := sub("^chr", "", as.character(chromosome))]
  out <- out[chr %in% chr_order]
  out[, genome_start := as.numeric(start) + chr_cumstart[chr]]
  out[, genome_end := as.numeric(end) + chr_cumstart[chr]]
  out[, value := as.numeric(get(value_col))]
  out[is.finite(genome_start) & is.finite(genome_end) & genome_end >= genome_start]
}

load_ploidy <- function(wgs_id) {
  pp_path <- file.path(wes_root, "tables/facets", paste0("Auto_", wgs_id, "_facets_purity_ploidy.tsv"))
  if (!file.exists(pp_path)) stop("Missing FACETS purity/ploidy: ", pp_path)
  pp <- fread(pp_path)
  purity <- as.numeric(pp$purity[1])
  ploidy <- as.numeric(pp$ploidy[1])
  if (!is.finite(ploidy) || ploidy <= 0) stop("Invalid ploidy for ", wgs_id)
  data.table(wgs_id = wgs_id, purity = purity, ploidy = ploidy, ploidy_shift_log2 = log2(ploidy / 2))
}

load_wes_profiles <- function(wgs_id, shift_policy) {
  summary_path <- file.path(wes_root, "tables/cns_highres", "Auto_wes_subclone_cns_highres_summary.csv")
  if (!file.exists(summary_path)) stop("Missing high-resolution WES CNS summary: ", summary_path)
  target_wgs_id <- wgs_id
  summary_dt <- fread(summary_path)[wgs_id == target_wgs_id]
  if (nrow(summary_dt) == 0) stop("No high-resolution WES CNS profiles for ", wgs_id)
  summary_dt[, is_bulk := profile == "bulk_highres"]
  summary_dt[, median_ccf_sort := fifelse(is.finite(median_ccf), median_ccf, Inf)]
  summary_dt <- summary_dt[order(!is_bulk, -median_ccf_sort, cluster_id)]
  pp <- load_ploidy(wgs_id)
  shift <- shift_policy$applied_shift_log2[1]
  facets_shift <- pp$ploidy_shift_log2[1]

  profiles <- list()
  profile_meta <- list()
  cns_manifest <- list()

  for (i in seq_len(nrow(summary_dt))) {
    path <- summary_dt$cns_path[i]
    if (!file.exists(path)) stop("Missing high-resolution CNS file: ", path)
    cns <- fread(path)
    cns[, centered_log2 := as.numeric(log2)]
    cns[, adjusted_log2 := centered_log2 + shift]
    cns[, absolute_log2 := adjusted_log2]
    cns[, `:=`(
      facets_purity = pp$purity[1],
      facets_ploidy = pp$ploidy[1],
      facets_ploidy_shift_log2 = facets_shift,
      applied_shift_log2 = shift,
      shift_applied = shift_policy$shift_applied[1],
      numbat_native_median_log2 = shift_policy$numbat_native_median_log2[1]
    )]

    is_bulk <- summary_dt$profile[i] == "bulk_highres"
    label <- if (is_bulk) {
      "WES CNVkit high-res"
    } else {
      paste0(
        "WES Cluster ", summary_dt$cluster_id[i],
        " CCF=", round(summary_dt$median_ccf[i], 2),
        " n=", summary_dt$n_mutations[i]
      )
    }

    centered <- as_segment_profile(cns, "centered_log2")
    absolute <- as_segment_profile(cns, "adjusted_log2")
    profiles[[paste0(label, " centered")]] <- centered
    profiles[[paste0(label, " baseline-adjusted")]] <- absolute
    profile_meta[[length(profile_meta) + 1]] <- data.table(
      wgs_id = wgs_id,
      profile = summary_dt$profile[i],
      display_label = label,
      cluster_id = as.character(summary_dt$cluster_id[i]),
      median_ccf = as.numeric(summary_dt$median_ccf[i]),
      n_mutations = as.integer(summary_dt$n_mutations[i]),
      n_segments = nrow(cns),
      n_nonzero_centered = sum(abs(cns$centered_log2) > 0.01, na.rm = TRUE),
      n_positive_absolute = sum(cns$adjusted_log2 > 0.01, na.rm = TRUE),
      median_centered_log2 = median(cns$centered_log2, na.rm = TRUE),
      median_absolute_log2 = median(cns$adjusted_log2, na.rm = TRUE),
      facets_purity = pp$purity[1],
      facets_ploidy = pp$ploidy[1],
      facets_ploidy_shift_log2 = facets_shift,
      applied_shift_log2 = shift,
      shift_applied = shift_policy$shift_applied[1],
      numbat_native_median_log2 = shift_policy$numbat_native_median_log2[1],
      shift_reason = shift_policy$shift_reason[1],
      source_cns = path
    )

    out_cns <- copy(cns)
    out_cns[, log2 := round(adjusted_log2, 6)]
    out_cns[, `:=`(
      centered_log2 = round(centered_log2, 6),
      absolute_log2 = round(adjusted_log2, 6),
      adjusted_log2 = round(adjusted_log2, 6)
    )]
    out_name <- paste0("Auto_", wgs_id, "_", sanitize_name(summary_dt$profile[i]), "_conditional_shift_absolute.cns")
    out_path <- file.path(cns_dir, out_name)
    fwrite(out_cns[, .(
      chromosome, start, end, gene, log2, depth, probes, weight,
      centered_log2, adjusted_log2, absolute_log2,
      facets_purity, facets_ploidy, facets_ploidy_shift_log2,
      applied_shift_log2, shift_applied, numbat_native_median_log2
    )], out_path, sep = "\t")
    cns_manifest[[length(cns_manifest) + 1]] <- data.table(
      wgs_id = wgs_id,
      profile = summary_dt$profile[i],
      display_label = label,
      output_cns = out_path
    )
  }

  list(
    profiles = profiles,
    meta = rbindlist(profile_meta, use.names = TRUE, fill = TRUE),
    cns_manifest = rbindlist(cns_manifest, use.names = TRUE, fill = TRUE),
    ploidy = pp,
    shift_policy = shift_policy
  )
}

load_numbat_native <- function(scrna_sample) {
  numbat_dir <- file.path(numbat_root, scrna_sample, "numbat")
  if (!dir.exists(numbat_dir)) return(NULL)

  nb <- read_numbat_bulk_clones(scrna_sample)
  bulk_point <- NULL
  clone_profiles <- list()
  clone_info <- data.table()
  if (!is.null(nb) && nrow(nb) > 0) {
    nb[, chr := as.character(CHROM)]
    nb <- nb[chr %in% chr_order]
    nb[, `:=`(
      genome_pos = (gene_start + gene_end) / 2 + chr_cumstart[chr],
      phi_val = suppressWarnings(as.numeric(phi_mle_roll)),
      n_cells_num = suppressWarnings(as.numeric(n_cells))
    )]
    nb <- nb[is.finite(genome_pos) & is.finite(phi_val) & phi_val > 0]
    nb_clone <- nb[nchar(members) > 0 & members != '""']
    if (nrow(nb_clone) > 0) {
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
      clone_info <- clone_info[seq_len(min(5L, nrow(clone_info)))]
      clone_info[, clone_label := paste0("Numbat Clone ", seq_len(.N), " n=", n_cells)]
      clone_profiles <- lapply(seq_len(nrow(clone_info)), function(i) {
        d <- copy(nb_clone[members == clone_info$members[i]])
        d[order(genome_pos), .(genome_pos, value)]
      })
      names(clone_profiles) <- clone_info$clone_label
    }
  }

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
    bulk_segment <- segs[is.finite(genome_start) & is.finite(genome_end), .(genome_start, genome_end, value)]
  }

  list(
    bulk_point = bulk_point,
    bulk_segment = bulk_segment,
    clone_profiles = clone_profiles,
    clone_info = clone_info
  )
}

load_infercna_profile <- function(scrna_sample) {
  if (!file.exists(infercna_outs_path) || !file.exists(infercna_meta_path) || !file.exists(gene_order_path)) return(NULL)
  gene_order <- as.data.table(read.table(gene_order_path, header = FALSE, col.names = c("gene_id", "chromosome", "start", "end")))
  gene_order[, chr := sub("^chr", "", chromosome)]
  gene_order <- gene_order[chr %in% chr_order]
  gene_order[, genome_pos := (start + end) / 2 + chr_cumstart[chr]]

  outs <- readRDS(infercna_outs_path)
  meta <- fread(infercna_meta_path)
  sample_cells <- meta[sample == scrna_sample, cell]
  sample_cells <- intersect(sample_cells, colnames(outs))
  if (length(sample_cells) < 10) return(NULL)

  mean_profile <- rowMeans(outs[, sample_cells, drop = FALSE], na.rm = TRUE)
  common_genes <- intersect(names(mean_profile), gene_order$gene_id)
  if (length(common_genes) < 1000) return(NULL)
  go <- gene_order[gene_id %in% common_genes]
  go[, value_raw := mean_profile[gene_id]]
  cna_sd <- sd(go$value_raw, na.rm = TRUE)
  if (is.finite(cna_sd) && cna_sd > 0) {
    go[, value := (value_raw - mean(value_raw, na.rm = TRUE)) / cna_sd * 0.3]
  } else {
    go[, value := value_raw]
  }
  setorder(go, genome_pos)
  window <- min(100L, nrow(go) %/% 5)
  go[, value := frollmean(value, n = window, align = "center", na.rm = TRUE)]
  go[!is.finite(value), value := value_raw]
  go[, .(genome_pos, value)]
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

make_segment_panel <- function(dt, title, y_lab = "log2", y_limits = c(-1.5, 2.2), show_chr_labels = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_segment(data = dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.25, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_point_panel <- function(dt, title, y_lab = "log2", y_limits = c(-1.5, 2.2), show_chr_labels = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity() +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.65) +
    geom_ribbon(data = dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.65) +
    geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.22, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_facet_segments <- function(profile_list, title, y_lab = "log2", y_limits = c(-1.5, 2.2)) {
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
    geom_segment(data = dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.18, color = "grey25") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
    facet_wrap(~profile, ncol = 1, strip.position = "right") +
    chr_x_scale +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = y_lab) +
    base_theme +
    theme(strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"))
}

make_facet_points <- function(profile_list, title, y_lab = "log2", y_limits = c(-1.5, 2.2)) {
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
    geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.18, color = "grey25") +
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
    geom_text(aes(label = ifelse(is.finite(correlation), sprintf("%.2f", correlation), "NA")), size = 2.8) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey85") +
    labs(title = "5 Mb binned Pearson correlations", x = "scRNA profile", y = "WES profile", fill = "r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

make_plot <- function(wgs_id, scrna_sample, wes_data, numbat_data, infercna_profile, cor_dt) {
  centered_bulk <- wes_data$profiles[["WES CNVkit high-res centered"]]
  absolute_bulk <- wes_data$profiles[["WES CNVkit high-res baseline-adjusted"]]
  sub_abs <- wes_data$profiles[grepl("^WES Cluster .* baseline-adjusted$", names(wes_data$profiles))]
  names(sub_abs) <- sub(" baseline-adjusted$", "", names(sub_abs))

  panels <- list(
    wes_centered = make_segment_panel(
      centered_bulk,
      paste0("WES CNVkit high-res centered shape - ", wgs_id),
      y_lab = "centered log2",
      y_limits = c(-1.5, 1.5)
    ),
    wes_absolute = make_segment_panel(
      absolute_bulk,
      paste0(
        "WES CNVkit conditional-shift absolute - ", wgs_id,
        " (ploidy=", round(wes_data$ploidy$ploidy[1], 2),
        ", FACETS shift=", round(wes_data$shift_policy$facets_ploidy_shift_log2[1], 2),
        ", applied shift=", round(wes_data$shift_policy$applied_shift_log2[1], 2), ")"
      ),
      y_lab = "log2 vs diploid",
      y_limits = c(-1.5, 2.2)
    )
  )

  sub_panel <- make_facet_segments(sub_abs, "WES projected subclone CNA profiles, conditional-shift", "log2 vs diploid", c(-1.5, 2.2))
  if (!is.null(sub_panel)) panels$wes_subclones <- sub_panel

  if (!is.null(numbat_data$bulk_point)) {
    panels$numbat_bulk <- make_point_panel(
      numbat_data$bulk_point,
      paste0("Numbat pseudo-bulk native high-res - ", scrna_sample, " (weighted log2 phi)"),
      y_lab = "log2 phi",
      y_limits = c(-1.5, 2.2)
    )
  } else if (!is.null(numbat_data$bulk_segment)) {
    panels$numbat_bulk <- make_segment_panel(
      numbat_data$bulk_segment,
      paste0("Numbat pseudo-bulk native consensus - ", scrna_sample, " (log2 phi)"),
      y_lab = "log2 phi",
      y_limits = c(-1.5, 2.2)
    )
  }
  if (length(numbat_data$clone_profiles) > 0) {
    panels$numbat_clones <- make_facet_points(
      numbat_data$clone_profiles,
      paste0("Numbat clone CNA profiles native scale - ", scrna_sample),
      "log2 phi",
      c(-1.5, 2.2)
    )
  }
  if (!is.null(infercna_profile)) {
    panels$infercna <- make_point_panel(
      infercna_profile,
      paste0("inferCNA mean - ", scrna_sample),
      y_lab = "scaled CNA",
      y_limits = c(-1.5, 1.5),
      show_chr_labels = TRUE
    )
  }
  heatmap <- make_heatmap(cor_dt)
  if (!is.null(heatmap)) panels$heatmap <- heatmap

  heights <- vapply(names(panels), function(nm) {
    if (nm == "wes_subclones") max(1.4, length(sub_abs) * 0.75)
    else if (nm == "numbat_clones") max(1.8, length(numbat_data$clone_profiles) * 0.65)
    else if (nm == "heatmap") 2.6
    else 1
  }, numeric(1))

  Reduce(`/`, panels) + plot_layout(heights = heights) +
    plot_annotation(
      title = paste0(scrna_sample, " - high-resolution WES/scRNA CNA comparison"),
      subtitle = "CNVkit shape is preserved for matching; WES shift is applied only when native Numbat indicates a global amplified baseline."
    )
}

message("Starting high-resolution conditional-shift WES/scRNA CNA comparison")
all_profile_meta <- list()
all_cns_manifest <- list()
all_cor <- list()
all_summary <- list()

for (wgs_id in names(wgs_sample_map)) {
  message("Processing ", wgs_id)
  shift_policy <- choose_shift_policy(wgs_id)
  message(
    "  shift policy: applied_shift=", round(shift_policy$applied_shift_log2[1], 3),
    " facets_shift=", round(shift_policy$facets_ploidy_shift_log2[1], 3),
    " numbat_median=", round(shift_policy$numbat_native_median_log2[1], 3)
  )
  wes_data <- load_wes_profiles(wgs_id, shift_policy)
  all_profile_meta[[wgs_id]] <- wes_data$meta
  all_cns_manifest[[wgs_id]] <- wes_data$cns_manifest

  bulk_meta <- wes_data$meta[profile == "bulk_highres"]
  all_summary[[wgs_id]] <- data.table(
    wgs_id = wgs_id,
    purity = wes_data$ploidy$purity[1],
    ploidy = wes_data$ploidy$ploidy[1],
    facets_ploidy_shift_log2 = shift_policy$facets_ploidy_shift_log2[1],
    applied_shift_log2 = shift_policy$applied_shift_log2[1],
    shift_applied = shift_policy$shift_applied[1],
    numbat_native_median_log2 = shift_policy$numbat_native_median_log2[1],
    shift_reason = shift_policy$shift_reason[1],
    cnvkit_segments = bulk_meta$n_segments[1],
    median_centered_log2 = bulk_meta$median_centered_log2[1],
    median_absolute_log2 = bulk_meta$median_absolute_log2[1],
    n_positive_absolute = bulk_meta$n_positive_absolute[1]
  )

  for (scrna_sample in wgs_sample_map[[wgs_id]]) {
    message("  Comparing ", scrna_sample)
    numbat_data <- load_numbat_native(scrna_sample)
    if (is.null(numbat_data)) next
    infercna_profile <- load_infercna_profile(scrna_sample)

    scrna_profiles <- list()
    if (!is.null(numbat_data$bulk_point)) {
      scrna_profiles[["Numbat pseudo-bulk native high-res"]] <- numbat_data$bulk_point
    } else if (!is.null(numbat_data$bulk_segment)) {
      scrna_profiles[["Numbat pseudo-bulk native consensus"]] <- numbat_data$bulk_segment
    }
    scrna_profiles <- c(scrna_profiles, numbat_data$clone_profiles)
    if (!is.null(infercna_profile)) scrna_profiles[["inferCNA mean"]] <- infercna_profile
    if (length(scrna_profiles) == 0) next

    plot_wes_profiles <- wes_data$profiles[
      names(wes_data$profiles) == "WES CNVkit high-res centered" |
        names(wes_data$profiles) == "WES CNVkit high-res baseline-adjusted" |
        grepl("^WES Cluster .* baseline-adjusted$", names(wes_data$profiles))
    ]
    corr_5 <- correlate_profiles(plot_wes_profiles, scrna_profiles, 5e6)
    corr_1 <- correlate_profiles(plot_wes_profiles, scrna_profiles, 1e6)
    cor_dt <- rbindlist(list(corr_5, corr_1), use.names = TRUE, fill = TRUE)
    cor_dt[, `:=`(wgs_id = wgs_id, sample = scrna_sample)]
    all_cor[[scrna_sample]] <- cor_dt
    fwrite(cor_dt, file.path(table_dir, paste0("Auto_wes_absolute_cna_correlations_", scrna_sample, ".csv")))

    plot_obj <- make_plot(wgs_id, scrna_sample, wes_data, numbat_data, infercna_profile, cor_dt)
    n_sub <- sum(grepl("^WES Cluster .* baseline-adjusted$", names(wes_data$profiles)))
    height <- max(12, 4.5 + 0.9 * n_sub + 0.9 * length(numbat_data$clone_profiles))
    pdf_path <- file.path(fig_dir, paste0("Auto_wes_absolute_cna_compare_", scrna_sample, ".pdf"))
    png_path <- file.path(fig_dir, paste0("Auto_wes_absolute_cna_compare_", scrna_sample, ".png"))
    ggsave(pdf_path, plot_obj, width = 16, height = height, limitsize = FALSE)
    ggsave(png_path, plot_obj, width = 16, height = height, dpi = 200, limitsize = FALSE)
  }
}

profile_meta <- rbindlist(all_profile_meta, use.names = TRUE, fill = TRUE)
cns_manifest <- rbindlist(all_cns_manifest, use.names = TRUE, fill = TRUE)
summary_dt <- rbindlist(all_summary, use.names = TRUE, fill = TRUE)
cor_dt <- rbindlist(all_cor, use.names = TRUE, fill = TRUE)

fwrite(summary_dt, file.path(table_dir, "Auto_wes_absolute_cna_summary.csv"))
fwrite(profile_meta, file.path(table_dir, "Auto_wes_absolute_cna_profile_summary.csv"))
fwrite(cns_manifest, file.path(table_dir, "Auto_wes_absolute_cna_cns_manifest.csv"))
fwrite(cor_dt, file.path(table_dir, "Auto_wes_absolute_cna_correlations.csv"))
fwrite(
  data.table(
    conclusion = c(
      "Corrected figures use CNVkit-resolution WES tracks, not FACETS-only coarse segments.",
      "Displayed WES absolute values are centered CNVkit log2 plus the conditional applied shift.",
      "The FACETS ploidy offset is applied only when native Numbat median log2 indicates a global amplified baseline.",
      "SUR1090 therefore keeps applied shift 0, while SUR1181 keeps the FACETS ploidy offset and remains globally amplified.",
      "Numbat pseudo-bulk and clone profiles are plotted from gene-level bulk_clones phi_mle_roll when available.",
      "WES projected subclone tracks are displayed for comparison, but true clone-specific WES CNA remains limited by single-sample bulk deconvolution."
    )
  ),
  file.path(table_dir, "Auto_wes_absolute_cna_method_assessment.csv")
)
fwrite(
  data.table(
    finished = as.character(Sys.time()),
    n_wgs = length(all_summary),
    n_correlation_rows = nrow(cor_dt),
    output_fig_dir = fig_dir,
    output_cns_dir = cns_dir
  ),
  file.path(log_dir, "Auto_wes_absolute_cna_compare_summary.tsv"),
  sep = "\t"
)

message("Done.")
