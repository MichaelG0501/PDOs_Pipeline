####################
# Auto_wes_scrna_subclone_highres_audit.R
#
# Analysis registry
# Status: active diagnostic/replot
# Script: analysis/cnv/wes_subclone/Auto_wes_scrna_subclone_highres_audit.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Sarek CNVkit .cns/.somatic.call.cns under spatialtranscriptomics live
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_segments.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_purity_ploidy.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_results.tsv
#   - ephemeral PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat outputs
#   - ephemeral PDOs_outs/cnv/Auto_PDO_infercna_* outputs
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/tables/cns_highres/Auto_<sample>_*.cns
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation_highres/*.csv
#   - PDOs_outs/Auto_wes_subclone/figures_highres/Auto_wes_scrna_subclone_match_highres_<sample>.pdf
#   - PDOs_outs/cnv/cnv_compare_highres/Auto_PDO_cnv_compare_highres_*.csv
# Downstream use: terminal audit of WES/scRNA CNA concordance and diagnosis of
#   resolution loss in per-WES-subclone CNA projections.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

live_project <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
ephemeral_project <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(live_project, "PDOs_outs/Auto_wes_subclone")
cnv_compare_dir <- file.path(live_project, "PDOs_outs/cnv/cnv_compare_highres")
cnvkit_root <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit"
numbat_root_eph <- file.path(ephemeral_project, "PDOs_outs/Auto_PDO_numbat/by_samples")
infercna_outs_path <- file.path(ephemeral_project, "PDOs_outs/cnv/Auto_PDO_infercna_outs_Carroll_2023.rds")
infercna_meta_path <- file.path(ephemeral_project, "PDOs_outs/cnv/Auto_PDO_infercna_meta_Carroll_2023.csv")
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"

cns_highres_dir <- file.path(out_root, "tables/cns_highres")
fig_highres_dir <- file.path(out_root, "figures_highres")
table_highres_dir <- file.path(out_root, "tables/visualisation_highres")
log_dir <- file.path(out_root, "logs")
for (d in c(cns_highres_dir, fig_highres_dir, table_highres_dir, cnv_compare_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

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
  `21` = 46709983, `22` = 50818468, `X` = 156040895
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

base_cnv_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin = margin(2, 6, 2, 6)
  )

chr_x_scale <- scale_x_continuous(
  limits = c(0, genome_len),
  expand = c(0, 0),
  breaks = chr_mids,
  labels = chr_order
)

as_segment_profile <- function(dt) {
  out <- copy(dt)
  out[, chr := sub("^chr", "", as.character(chromosome))]
  out <- out[chr %in% chr_order]
  out[, genome_start := as.numeric(start) + chr_cumstart[chr]]
  out[, genome_end := as.numeric(end) + chr_cumstart[chr]]
  out[, value := as.numeric(log2)]
  out[is.finite(genome_start) & is.finite(genome_end) & genome_end >= genome_start]
}

load_cnvkit_segments <- function(wgs_id) {
  tumor_id <- sub("_vs_.*", "", wgs_id)
  candidate_paths <- file.path(
    cnvkit_root,
    wgs_id,
    paste0(tumor_id, c(".somatic.call.cns", ".cns", ".call.cns"))
  )
  candidate_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(candidate_paths) == 0) stop("No CNVkit .cns found for ", wgs_id)

  cns <- fread(candidate_paths[1])
  required <- c("chromosome", "start", "end", "log2")
  if (!all(required %in% names(cns))) {
    stop("CNVkit file lacks required columns: ", candidate_paths[1])
  }
  cns <- as_segment_profile(cns)
  cns[, source_path := candidate_paths[1]]
  cns[]
}

load_facets_segments <- function(wgs_id) {
  seg_path <- file.path(out_root, "tables/facets", paste0("Auto_", wgs_id, "_facets_segments.tsv"))
  pp_path <- file.path(out_root, "tables/facets", paste0("Auto_", wgs_id, "_facets_purity_ploidy.tsv"))
  if (!file.exists(seg_path) || !file.exists(pp_path)) stop("Missing FACETS outputs for ", wgs_id)
  seg <- fread(seg_path)
  pp <- fread(pp_path)
  purity <- as.numeric(pp$purity[1])
  ploidy <- as.numeric(pp$ploidy[1])
  seg[, chr := sub("^chr", "", as.character(chrom))]
  seg <- seg[chr %in% chr_order]
  seg[, `:=`(
    start = as.integer(start),
    end = as.integer(end),
    facets_ccf = pmin(1, pmax(0, as.numeric(cf) / purity)),
    event_log2 = log2((purity * total_cn + (1 - purity) * 2) / (purity * ploidy + (1 - purity) * 2))
  )]
  list(segments = seg, purity = purity, ploidy = ploidy)
}

load_cluster_summary <- function(wgs_id) {
  res_path <- file.path(out_root, "tables/pyclone", paste0("Auto_", wgs_id, "_pyclone_vi_results.tsv"))
  if (!file.exists(res_path)) stop("Missing PyClone-VI results: ", res_path)
  res <- fread(res_path)
  res[, cluster_id := as.character(cluster_id)]
  out <- res[, .(
    median_ccf = median(cellular_prevalence, na.rm = TRUE),
    mean_assignment_prob = mean(cluster_assignment_prob, na.rm = TRUE),
    n_mutations = .N
  ), by = cluster_id]
  setorder(out, -median_ccf, cluster_id)
  out[]
}

map_facets_to_cnvkit <- function(cnvkit, facets) {
  cnv_mid <- cnvkit[, .(
    chromosome,
    chr,
    start = as.integer((start + end) / 2),
    end = as.integer((start + end) / 2),
    cnvkit_row = .I
  )]
  fac <- facets[, .(
    chr,
    start,
    end,
    total_cn,
    minor_cn,
    major_cn,
    cf,
    facets_ccf,
    event_log2
  )]
  setkey(fac, chr, start, end)
  setkey(cnv_mid, chr, start, end)
  mapped <- foverlaps(cnv_mid, fac, nomatch = NA)
  mapped <- mapped[order(cnvkit_row)]
  mapped[, .(
    cnvkit_row,
    total_cn,
    minor_cn,
    major_cn,
    cf,
    facets_ccf,
    event_log2
  )]
}

write_highres_cns <- function(wgs_id, cnvkit, facets_data, cluster_summary, ccf_tolerance = 0.15) {
  facets_map <- map_facets_to_cnvkit(cnvkit, facets_data$segments)
  base_out <- cnvkit[, .(
    chromosome,
    start = as.integer(start),
    end = as.integer(end),
    gene = if ("gene" %in% names(cnvkit)) gene else "-",
    log2 = round(value, 6),
    depth = if ("depth" %in% names(cnvkit)) depth else 0,
    probes = if ("probes" %in% names(cnvkit)) probes else 0,
    weight = if ("weight" %in% names(cnvkit)) weight else 1
  )]
  bulk_path <- file.path(cns_highres_dir, paste0("Auto_", wgs_id, "_bulk_highres.cns"))
  fwrite(base_out, bulk_path, sep = "\t")

  summary_rows <- list(data.table(
    wgs_id = wgs_id,
    profile = "bulk_highres",
    cluster_id = NA_character_,
    median_ccf = NA_real_,
    n_mutations = NA_integer_,
    n_segments = nrow(base_out),
    n_nonzero_segments = sum(abs(base_out$log2) > 0.01, na.rm = TRUE),
    cns_path = bulk_path
  ))

  profiles <- list("WES Bulk CNVkit highres" = as_segment_profile(base_out))
  for (i in seq_len(nrow(cluster_summary))) {
    cid <- cluster_summary$cluster_id[i]
    cluster_ccf <- cluster_summary$median_ccf[i]
    threshold <- max(0, cluster_ccf - ccf_tolerance)
    present <- is.na(facets_map$facets_ccf) | facets_map$facets_ccf >= threshold

    sub_out <- copy(base_out)
    sub_out[, log2 := fifelse(present, log2, 0)]
    sub_out[, probes := fifelse(abs(log2) > 0.01, probes, 0)]
    out_path <- file.path(cns_highres_dir, paste0("Auto_", wgs_id, "_cluster", cid, "_highres.cns"))
    fwrite(sub_out, out_path, sep = "\t")

    label <- paste0(
      "WES Cluster ", cid,
      " CCF=", round(cluster_ccf, 2),
      " n=", cluster_summary$n_mutations[i]
    )
    profiles[[label]] <- as_segment_profile(sub_out)

    summary_rows[[length(summary_rows) + 1]] <- data.table(
      wgs_id = wgs_id,
      profile = paste0("cluster", cid, "_highres"),
      cluster_id = cid,
      median_ccf = cluster_ccf,
      n_mutations = cluster_summary$n_mutations[i],
      n_segments = nrow(sub_out),
      n_nonzero_segments = sum(abs(sub_out$log2) > 0.01, na.rm = TRUE),
      cns_path = out_path
    )
  }

  list(profiles = profiles, summary = rbindlist(summary_rows), facets_map = facets_map)
}

final_iter_from <- function(numbat_dir, prefix = "treeML", ext = "rds") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.", ext)))
  if (length(files) == 0) return(NA_integer_)
  regex <- paste0("^", prefix, "_([0-9]+)\\.", gsub("\\.", "\\\\.", ext), "$")
  iter <- suppressWarnings(as.integer(sub(regex, "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

load_numbat_profiles <- function(scrna_sample) {
  numbat_dir <- file.path(numbat_root_eph, scrna_sample, "numbat")
  if (!dir.exists(numbat_dir)) return(NULL)

  bulk_file <- file.path(numbat_dir, "bulk_clones_final.tsv.gz")
  if (!file.exists(bulk_file)) {
    iter <- final_iter_from(numbat_dir, "bulk_clones", "tsv.gz")
    if (is.finite(iter)) bulk_file <- file.path(numbat_dir, paste0("bulk_clones_", iter, ".tsv.gz"))
  }
  if (!file.exists(bulk_file)) return(NULL)

  nb <- fread(
    bulk_file,
    select = c("CHROM", "gene_start", "gene_end", "n_cells", "members", "phi_mle_roll", "gene_index")
  )
  nb[, chr := as.character(CHROM)]
  nb <- nb[chr %in% chr_order]
  nb[, genome_pos := (gene_start + gene_end) / 2 + chr_cumstart[chr]]
  # Numbat defines phi_mle_roll as total copy-number ratio relative to diploid;
  # keep this scale and do not clone-center, otherwise global aneuploidy is lost.
  nb[, value := log2(as.numeric(phi_mle_roll))]
  nb <- nb[is.finite(genome_pos) & is.finite(value)]

  clone_info <- nb[nchar(members) > 0 & members != '""', .(n_cells = n_cells[1]), by = members]
  clone_info <- clone_info[order(-n_cells)]
  if (nrow(clone_info) == 0) return(NULL)

  cons_path <- file.path(
    ephemeral_project,
    "PDOs_outs/Auto_PDO_numbat/conservative_clones/by_samples",
    scrna_sample,
    paste0("Auto_", scrna_sample, "_clones_conservative.rds")
  )
  if (file.exists(cons_path)) {
    cons <- readRDS(cons_path)
    cons_df <- data.table(
      clone_opt = as.character(sapply(cons, `[[`, "sample")),
      cons_n = as.integer(sapply(cons, `[[`, "size"))
    )
    cons_df <- cons_df[order(-cons_n)]
    n_match <- min(nrow(clone_info), nrow(cons_df))
    clone_info <- clone_info[seq_len(n_match)]
    cons_df <- cons_df[seq_len(n_match)]
    clone_info[, n_cells := cons_df$cons_n]
    clone_info[, clone_label := paste0("Numbat ", cons_df$clone_opt, " n=", n_cells)]
  } else {
    clone_info <- clone_info[seq_len(min(5, nrow(clone_info)))]
    clone_info[, clone_label := paste0("Numbat Clone ", seq_len(.N), " n=", n_cells)]
  }

  nb_mapped <- merge(nb[members %in% clone_info$members], clone_info[, .(members, clone_label, n_cells)], by = "members")
  clone_profiles <- lapply(seq_len(nrow(clone_info)), function(i) {
    d <- copy(nb_mapped[members == clone_info$members[i]])
    d[order(genome_pos), .(genome_pos, value)]
  })
  names(clone_profiles) <- clone_info$clone_label

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

  list(clone_profiles = clone_profiles, bulk_segment = bulk_segment, clone_info = clone_info)
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

bin_segment_profile <- function(dt, bins) {
  sapply(seq_len(nrow(bins)), function(i) {
    hit <- dt[genome_start <= bins$genome_mid[i] & genome_end >= bins$genome_mid[i]]
    if (nrow(hit) == 0) NA_real_ else hit$value[1]
  })
}

bin_point_profile <- function(dt, bins) {
  sapply(seq_len(nrow(bins)), function(i) {
    hit <- dt[genome_pos >= bins$genome_start[i] & genome_pos < bins$genome_end[i]]
    if (nrow(hit) == 0) NA_real_ else mean(hit$value, na.rm = TRUE)
  })
}

correlate_profiles <- function(wes_profiles, scrna_profiles, bin_size) {
  bins <- make_bins(bin_size)
  wes_binned <- lapply(wes_profiles, bin_segment_profile, bins = bins)
  scrna_binned <- lapply(scrna_profiles, function(x) {
    if (all(c("genome_start", "genome_end") %in% names(x))) {
      bin_segment_profile(x, bins)
    } else {
      bin_point_profile(x, bins)
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

make_ribbon_panel <- function(dt, title, is_point = FALSE, show_chr_labels = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
    scale_fill_identity()
  if (is_point) {
    p <- p +
      geom_ribbon(data = dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.6) +
      geom_ribbon(data = dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.6) +
      geom_line(data = dt, aes(x = genome_pos, y = value), linewidth = 0.25, color = "grey25")
  } else {
    p <- p +
      geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.6) +
      geom_rect(data = dt, aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.6) +
      geom_segment(data = dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.25, color = "grey25")
  }
  p <- p +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    labs(title = title, x = NULL, y = "log2") +
    base_cnv_theme
  if (show_chr_labels) p <- p + theme(axis.text.x = element_text(size = 7))
  p
}

make_heatmap <- function(cor_dt, title) {
  plot_dt <- copy(cor_dt[bin_size_mb == 5])
  if (nrow(plot_dt) == 0) return(NULL)
  plot_dt[, wes_profile := factor(wes_profile, levels = unique(wes_profile))]
  plot_dt[, scrna_profile := factor(scrna_profile, levels = unique(scrna_profile))]
  ggplot(plot_dt, aes(x = scrna_profile, y = wes_profile, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = ifelse(is.finite(correlation), sprintf("%.2f", correlation), "NA")), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey85") +
    labs(title = title, x = "scRNA profile", y = "WES profile", fill = "r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
}

make_sample_plot <- function(scrna_sample, wgs_id, wes_profiles, numbat_data, infercna_profile, cor_dt) {
  panels <- list()
  panels[["wes_bulk"]] <- make_ribbon_panel(wes_profiles[[1]], paste0("WES bulk CNVkit highres - ", wgs_id))

  subclone_profiles <- wes_profiles[!grepl("Bulk", names(wes_profiles))]
  if (length(subclone_profiles) > 0) {
    sub_dt <- rbindlist(lapply(names(subclone_profiles), function(nm) {
      d <- copy(subclone_profiles[[nm]])
      d[, profile := nm]
      d
    }))
    sub_dt[, profile := factor(profile, levels = names(subclone_profiles))]
    panels[["wes_subclones"]] <- ggplot() +
      geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
      scale_fill_identity() +
      geom_rect(data = sub_dt, aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.6) +
      geom_rect(data = sub_dt, aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.6) +
      geom_segment(data = sub_dt, aes(x = genome_start, xend = genome_end, y = value, yend = value), linewidth = 0.2, color = "grey25") +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
      facet_wrap(~profile, ncol = 1, strip.position = "right") +
      chr_x_scale +
      coord_cartesian(ylim = c(-1.5, 1.5)) +
      labs(title = "WES subclone CNA projections on CNVkit grid", x = NULL, y = "log2") +
      base_cnv_theme +
      theme(strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"))
  }

  if (!is.null(numbat_data) && length(numbat_data$clone_profiles) > 0) {
    numb_dt <- rbindlist(lapply(names(numbat_data$clone_profiles), function(nm) {
      d <- copy(numbat_data$clone_profiles[[nm]])
      d[, profile := nm]
      d
    }))
    numb_dt[, profile := factor(profile, levels = names(numbat_data$clone_profiles))]
    panels[["numbat"]] <- ggplot() +
      geom_rect(data = chr_bands, aes(xmin = genome_start, xmax = genome_end, ymin = -Inf, ymax = Inf, fill = fill), color = NA) +
      scale_fill_identity() +
      geom_ribbon(data = numb_dt, aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)), fill = "#D9534F", alpha = 0.6) +
      geom_ribbon(data = numb_dt, aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0), fill = "#5B9BD5", alpha = 0.6) +
      geom_line(data = numb_dt, aes(x = genome_pos, y = value), linewidth = 0.2, color = "grey25") +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
      facet_wrap(~profile, ncol = 1, strip.position = "right") +
      chr_x_scale +
      coord_cartesian(ylim = c(-1.5, 1.5)) +
      labs(title = paste0("Numbat clone CNA profiles - ", scrna_sample), x = NULL, y = "log2") +
      base_cnv_theme +
      theme(strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"))
  }

  if (!is.null(infercna_profile)) {
    panels[["infercna"]] <- make_ribbon_panel(infercna_profile, paste0("inferCNA mean - ", scrna_sample), is_point = TRUE, show_chr_labels = TRUE)
  }
  heatmap <- make_heatmap(cor_dt, "5 Mb binned Pearson correlations")
  if (!is.null(heatmap)) panels[["heatmap"]] <- heatmap

  heights <- vapply(names(panels), function(nm) {
    if (nm == "wes_subclones") max(1.5, length(subclone_profiles) * 0.7)
    else if (nm == "numbat") max(1.5, length(numbat_data$clone_profiles) * 0.65)
    else if (nm == "heatmap") 2.5
    else 1
  }, numeric(1))
  Reduce(`/`, panels) + plot_layout(heights = heights) +
    plot_annotation(
      title = paste0(scrna_sample, " - WES/scRNA high-resolution CNA audit"),
      subtitle = "WES subclone profiles preserve CNVkit segment resolution; clone-specific CNA calls remain projections from bulk FACETS CCF."
    )
}

all_cns_summaries <- list()
all_correlations <- list()
all_bulk_compare <- list()
all_diagnostics <- list()

message("Starting high-resolution WES/scRNA CNA audit")
for (wgs_id in names(wgs_sample_map)) {
  message("Processing ", wgs_id)
  cnvkit <- load_cnvkit_segments(wgs_id)
  facets_data <- load_facets_segments(wgs_id)
  cluster_summary <- load_cluster_summary(wgs_id)
  highres <- write_highres_cns(wgs_id, cnvkit, facets_data, cluster_summary)
  wes_profiles <- highres$profiles
  all_cns_summaries[[wgs_id]] <- highres$summary

  facets_n <- nrow(facets_data$segments)
  cnvkit_n <- nrow(cnvkit)
  all_diagnostics[[wgs_id]] <- data.table(
    wgs_id = wgs_id,
    cnvkit_source = unique(cnvkit$source_path)[1],
    cnvkit_segments = cnvkit_n,
    facets_segments = facets_n,
    resolution_gain_vs_facets = round(cnvkit_n / facets_n, 2),
    pyclone_clusters = nrow(cluster_summary),
    pyclone_mutations = sum(cluster_summary$n_mutations),
    min_cluster_mutations = min(cluster_summary$n_mutations),
    min_mean_assignment_prob = min(cluster_summary$mean_assignment_prob, na.rm = TRUE)
  )

  for (scrna_sample in wgs_sample_map[[wgs_id]]) {
    message("  Matching ", scrna_sample)
    numbat_data <- load_numbat_profiles(scrna_sample)
    infercna_profile <- load_infercna_profile(scrna_sample)

    scrna_profiles <- list()
    if (!is.null(numbat_data)) {
      scrna_profiles <- c(scrna_profiles, numbat_data$clone_profiles)
      if (!is.null(numbat_data$bulk_segment)) scrna_profiles[["Numbat pseudo-bulk"]] <- numbat_data$bulk_segment
    }
    if (!is.null(infercna_profile)) scrna_profiles[["inferCNA mean"]] <- infercna_profile
    if (length(scrna_profiles) == 0) next

    corr_5 <- correlate_profiles(wes_profiles, scrna_profiles, 5e6)
    corr_1 <- correlate_profiles(wes_profiles, scrna_profiles, 1e6)
    corr <- rbindlist(list(corr_5, corr_1))
    corr[, `:=`(sample = scrna_sample, wgs_id = wgs_id)]
    all_correlations[[scrna_sample]] <- corr

    if ("Numbat pseudo-bulk" %in% names(scrna_profiles) || "inferCNA mean" %in% names(scrna_profiles)) {
      bulk_compare <- corr[wes_profile == "WES Bulk CNVkit highres" & scrna_profile %in% c("Numbat pseudo-bulk", "inferCNA mean")]
      all_bulk_compare[[scrna_sample]] <- bulk_compare
      fwrite(
        bulk_compare,
        file.path(cnv_compare_dir, paste0("Auto_PDO_cnv_compare_highres_", scrna_sample, ".csv"))
      )
    }

    plot_obj <- make_sample_plot(scrna_sample, wgs_id, wes_profiles, numbat_data, infercna_profile, corr)
    pdf_path <- file.path(fig_highres_dir, paste0("Auto_wes_scrna_subclone_match_highres_", scrna_sample, ".pdf"))
    png_path <- file.path(fig_highres_dir, paste0("Auto_wes_scrna_subclone_match_highres_", scrna_sample, ".png"))
    h <- max(12, 2.5 + 1.7 * length(wes_profiles) + if (!is.null(numbat_data)) 1.3 * length(numbat_data$clone_profiles) else 0)
    ggsave(pdf_path, plot_obj, width = 16, height = h, limitsize = FALSE)
    ggsave(png_path, plot_obj, width = 16, height = h, dpi = 200, limitsize = FALSE)
  }
}

cns_summary <- rbindlist(all_cns_summaries, use.names = TRUE, fill = TRUE)
diagnostics <- rbindlist(all_diagnostics, use.names = TRUE, fill = TRUE)
correlations <- rbindlist(all_correlations, use.names = TRUE, fill = TRUE)
bulk_compare <- rbindlist(all_bulk_compare, use.names = TRUE, fill = TRUE)

fwrite(cns_summary, file.path(cns_highres_dir, "Auto_wes_subclone_cns_highres_summary.csv"))
fwrite(diagnostics, file.path(table_highres_dir, "Auto_wes_scrna_subclone_highres_diagnostics.csv"))
fwrite(correlations, file.path(table_highres_dir, "Auto_wes_scrna_subclone_highres_correlations.csv"))
fwrite(bulk_compare, file.path(cnv_compare_dir, "Auto_PDO_cnv_compare_highres_summary.csv"))

best_matches <- correlations[!grepl("Bulk", wes_profile) & bin_size_mb == 5]
if (nrow(best_matches) > 0) {
  best_matches <- best_matches[order(sample, wes_profile, -correlation)]
  best_matches <- best_matches[, .SD[1], by = .(sample, wgs_id, wes_profile)]
  fwrite(best_matches, file.path(table_highres_dir, "Auto_wes_scrna_subclone_highres_best_matches.csv"))
}

fwrite(
  data.table(
    finished = as.character(Sys.time()),
    live_project = live_project,
    ephemeral_project = ephemeral_project,
    n_cns_profiles = nrow(cns_summary),
    n_correlation_rows = nrow(correlations)
  ),
  file.path(log_dir, "Auto_wes_scrna_subclone_highres_audit_summary.tsv"),
  sep = "\t"
)

message("Done.")
