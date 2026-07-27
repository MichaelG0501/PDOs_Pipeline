####################
# Auto_wes_scrna_subclone_match.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_wes_scrna_subclone_match.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/cns/Auto_<sample>_bulk.cns
#   - PDOs_outs/Auto_wes_subclone/tables/cns/Auto_<sample>_cluster<N>.cns
#   - PDOs_outs/Auto_wes_subclone/tables/cns/Auto_wes_subclone_cns_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_results.tsv
#   - Numbat bulk_clones_final.tsv.gz per sample (ephemeral, per-clone CNA)
#   - Optional: Numbat conservative clones .rds
#   - Optional: inferCNA outputs
#   - Gene order: hg38_gencode_v27.txt
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/figures/Auto_wes_scrna_subclone_match_<scrna_sample>.pdf
#   - PDOs_outs/Auto_wes_subclone/figures/Auto_wes_scrna_subclone_match_<scrna_sample>.png
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_scrna_subclone_match_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_scrna_subclone_correlation_matrix.csv
# Downstream use: terminal validation of which Numbat/inferCNA scRNA clones
#   correspond to which WES PyClone-VI subclones.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── Configuration ─────────────────────────────────────────────────────────────
out_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
numbat_root_eph <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_PDO_numbat/by_samples"
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"

# inferCNA paths (ephemeral)
infercna_outs_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/cnv/Auto_PDO_infercna_outs_Carroll_2023.rds"
infercna_meta_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/cnv/Auto_PDO_infercna_meta_Carroll_2023.csv"

cns_dir <- file.path(out_root, "tables/cns")
fig_dir <- file.path(out_root, "figures")
table_dir <- file.path(out_root, "tables/visualisation")
log_dir <- file.path(out_root, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# WGS sample → scRNA sample mapping
wgs_sample_map <- list(
  "PDO_1090_vs_NT_1090" = c("SUR1090_Untreated_PDO", "SUR1090_Treated_PDO"),
  "PDO_1181_vs_NT_1181" = c("SUR1181_Untreated_PDO", "SUR1181_Treated_PDO")
)

# ── Chromosome layout (hg38) ─────────────────────────────────────────────────
chr_order <- c(as.character(1:22), "X")
chr_sizes <- c(
  `1`=248956422, `2`=242193529, `3`=198295559, `4`=190214555, `5`=181538259,
  `6`=170805979, `7`=159345973, `8`=145138636, `9`=138394717, `10`=133797422,
  `11`=135086622, `12`=133275309, `13`=114364328, `14`=107043718, `15`=101991189,
  `16`=90338345, `17`=83257441, `18`=80373285, `19`=58617616, `20`=64444167,
  `21`=46709983, `22`=50818468, `X`=156040895
)

chr_cumstart <- cumsum(c(0, chr_sizes[chr_order][-length(chr_order)]))
names(chr_cumstart) <- chr_order
chr_cumend <- chr_cumstart + chr_sizes[chr_order]
chr_mids   <- (chr_cumstart + chr_cumend) / 2
genome_len <- max(chr_cumend)

chr_bands <- data.table(
  xmin = chr_cumstart[chr_order],
  xmax = chr_cumend[chr_order],
  chr  = chr_order,
  band_fill = ifelse(seq_along(chr_order) %% 2 == 1, "grey96", "white")
)

# ── Shared theme ──────────────────────────────────────────────────────────────
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
  limits = c(0, genome_len), expand = c(0, 0),
  breaks = chr_mids, labels = chr_order
)

# ── Helper: build a single ribbon panel ───────────────────────────────────────
make_ribbon_panel <- function(dt, title_str, y_lab = "log2 ratio",
                              ylim_range = c(-1.5, 1.5), show_chr_labels = FALSE,
                              line_col = "grey20", gain_col = "#D9534F",
                              loss_col = "#5B9BD5", alpha_val = 0.7,
                              is_point_data = FALSE) {
  p <- ggplot() +
    geom_rect(data = chr_bands,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_fill),
              color = NA, show.legend = FALSE) +
    scale_fill_identity()
    
  if (is_point_data) {
    p <- p +
      geom_ribbon(data = dt,
                  aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)),
                  fill = gain_col, alpha = alpha_val) +
      geom_ribbon(data = dt,
                  aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0),
                  fill = loss_col, alpha = alpha_val) +
      geom_line(data = dt,
                aes(x = genome_pos, y = value),
                linewidth = 0.35, color = line_col)
  } else {
    p <- p +
      geom_rect(data = dt,
                aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)),
                fill = gain_col, alpha = alpha_val) +
      geom_rect(data = dt,
                aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0),
                fill = loss_col, alpha = alpha_val) +
      geom_segment(data = dt,
                   aes(x = genome_start, xend = genome_end, y = value, yend = value),
                   linewidth = 0.35, color = line_col)
  }
  
  p <- p +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    chr_x_scale +
    coord_cartesian(ylim = ylim_range) +
    labs(title = title_str, y = y_lab, x = NULL) +
    base_cnv_theme

  if (show_chr_labels) {
    p <- p + theme(axis.text.x = element_text(size = 7, angle = 0))
  }
  p
}

# ── 1. Load per-subclone .cns files ──────────────────────────────────────────
load_wes_subclone_cns <- function(wgs_id) {
  # Load all .cns files for this sample
  cns_files <- list.files(cns_dir, pattern = paste0("Auto_", wgs_id, "_.*\\.cns$"),
                          full.names = TRUE)
  if (length(cns_files) == 0) {
    message("  No .cns files found for ", wgs_id)
    return(NULL)
  }

  # Load PyClone results for cluster CCF info
  res_path <- file.path(out_root, "tables/pyclone",
                        paste0("Auto_", wgs_id, "_pyclone_vi_results.tsv"))
  cluster_ccf <- NULL
  if (file.exists(res_path)) {
    res <- fread(res_path)
    res[, cluster_id := as.character(cluster_id)]
    cluster_ccf <- res[, .(median_ccf = median(cellular_prevalence, na.rm = TRUE),
                           n_mutations = .N),
                       by = cluster_id]
    setorder(cluster_ccf, -median_ccf)
  }

  profiles <- list()
  for (f in cns_files) {
    bn <- basename(f)
    cns <- fread(f)
    cns[, chr := sub("^chr", "", chromosome)]
    cns <- cns[chr %in% chr_order]
    cns[, genome_start := start + chr_cumstart[chr]]
    cns[, genome_end   := end   + chr_cumstart[chr]]
    cns[, value := log2]

    # Parse label from filename
    if (grepl("_bulk\\.cns$", bn)) {
      label <- "WES Bulk"
    } else {
      cid <- sub(".*_cluster(\\d+)\\.cns$", "\\1", bn)
      ccf_str <- ""
      n_str <- ""
      if (!is.null(cluster_ccf) && cid %in% cluster_ccf$cluster_id) {
        ccf_str <- paste0(", CCF=", round(cluster_ccf[cluster_id == cid, median_ccf], 2))
        n_str <- paste0(", n=", cluster_ccf[cluster_id == cid, n_mutations])
      }
      label <- paste0("WES Cluster ", cid, ccf_str, n_str)
    }

    profiles[[label]] <- cns
  }

  list(profiles = profiles, cluster_ccf = cluster_ccf)
}

# ── 2. Load Numbat per-clone profiles ────────────────────────────────────────
load_numbat_clones <- function(scrna_sample) {
  numbat_dir <- file.path(numbat_root_eph, scrna_sample, "numbat")
  if (!dir.exists(numbat_dir)) {
    message("  Numbat dir missing: ", numbat_dir)
    return(NULL)
  }

  bulk_file <- file.path(numbat_dir, "bulk_clones_final.tsv.gz")
  if (!file.exists(bulk_file)) {
    # Try iteration files
    files <- Sys.glob(file.path(numbat_dir, "bulk_clones_*.tsv.gz"))
    iters <- suppressWarnings(as.integer(sub(".*bulk_clones_(\\d+)\\.tsv\\.gz$", "\\1", basename(files))))
    iters <- iters[is.finite(iters)]
    if (length(iters) > 0) {
      bulk_file <- file.path(numbat_dir, paste0("bulk_clones_", max(iters), ".tsv.gz"))
    }
  }
  if (!file.exists(bulk_file)) {
    message("  Numbat bulk_clones not found")
    return(NULL)
  }

  nb <- fread(bulk_file,
              select = c("CHROM", "gene_start", "gene_end", "n_cells",
                         "members", "phi_mle_roll", "gene_index"))
  nb[, chr := as.character(CHROM)]
  nb <- nb[chr %in% chr_order]
  nb[, genome_mid := (gene_start + gene_end) / 2 + chr_cumstart[chr]]
  nb[, log2_phi := log2(phi_mle_roll)]

  # Get clone info
  clone_info <- nb[nchar(members) > 0 & members != '""',
                   .(n_cells = n_cells[1]), by = members]
  clone_info <- clone_info[order(-n_cells)]

  # Try conservative clones
  cons_path <- file.path(numbat_root_eph, "..", "conservative_clones/by_samples",
                         scrna_sample,
                         paste0("Auto_", scrna_sample, "_clones_conservative.rds"))
  if (file.exists(cons_path)) {
    message("  Using conservative Numbat clones")
    cons <- readRDS(cons_path)
    cons_sizes <- as.integer(sapply(cons, `[[`, "size"))
    cons_names <- as.character(sapply(cons, `[[`, "sample"))
    cons_df <- data.table(clone_opt = cons_names, cons_n = cons_sizes)
    cons_df <- cons_df[order(-cons_n)]
    n_match <- min(nrow(clone_info), nrow(cons_df))
    clone_info <- clone_info[1:n_match]
    cons_df <- cons_df[1:n_match]
    clone_info[, n_cells := cons_df$cons_n]
    clone_info[, clone_label := paste0("Numbat ", cons_df$clone_opt, " (n=", n_cells, ")")]
  } else {
    top_n <- min(5, nrow(clone_info))
    clone_info <- clone_info[1:top_n]
    clone_info[, clone_label := paste0("Numbat Clone ", seq_len(.N), " (n=", n_cells, ")")]
  }

  if (nrow(clone_info) == 0) {
    message("  No Numbat clones with members")
    return(NULL)
  }

  # Build per-clone profiles
  nb_mapped <- nb[members %in% clone_info$members]
  nb_mapped <- merge(nb_mapped[, -"n_cells", with = FALSE],
                     clone_info[, .(members, n_cells)], by = "members")

  clone_profiles <- list()
  for (i in seq_len(nrow(clone_info))) {
    d <- nb_mapped[members == clone_info$members[i]]
    d <- d[order(genome_mid)]
    clone_profiles[[clone_info$clone_label[i]]] <- data.table(
      genome_pos = d$genome_mid,
      value = d$log2_phi
    )
  }

  # Pseudo-bulk from all clones (weighted by n_cells)
  segs_file <- file.path(numbat_dir,
                         paste0("Auto_", scrna_sample, "_numbat_segs_consensus.csv"))
  bulk_ribbon <- NULL
  if (file.exists(segs_file)) {
    segs <- fread(segs_file)
    segs[, chr := as.character(CHROM)]
    segs <- segs[chr %in% chr_order]
    segs[, phi_val := suppressWarnings(as.numeric(phi_mle))]
    segs[, log2_seg := ifelse(is.finite(phi_val) & phi_val > 0, log2(phi_val), 0)]
    segs[, genome_start := seg_start + chr_cumstart[chr]]
    segs[, genome_end   := seg_end   + chr_cumstart[chr]]
    segs[, value := log2_seg]
    bulk_ribbon <- segs
  }

  list(
    clone_profiles = clone_profiles,
    bulk_ribbon = bulk_ribbon,
    clone_info = clone_info
  )
}

# ── 3. Load inferCNA ─────────────────────────────────────────────────────────
load_infercna <- function(scrna_sample) {
  if (!file.exists(infercna_outs_path) || !file.exists(infercna_meta_path)) {
    return(NULL)
  }

  gene_order <- as.data.table(read.table(gene_order_path, header = FALSE,
                                         col.names = c("gene_id", "chromosome", "start", "end")))
  gene_order[, chr := sub("^chr", "", chromosome)]
  gene_order <- gene_order[chr %in% chr_order]
  gene_order[, chr := factor(chr, levels = chr_order)]
  gene_order <- gene_order[order(chr, start)]
  gene_order[, genome_mid := (start + end) / 2 + chr_cumstart[as.character(chr)]]

  outs <- readRDS(infercna_outs_path)
  meta <- read.csv(infercna_meta_path)

  sample_cells <- meta$cell[meta$sample == scrna_sample]
  if (length(sample_cells) == 0) return(NULL)
  sample_cells <- intersect(sample_cells, colnames(outs))
  if (length(sample_cells) < 10) return(NULL)

  sample_mat <- outs[, sample_cells, drop = FALSE]
  mean_profile <- rowMeans(sample_mat, na.rm = TRUE)

  common_genes <- intersect(names(mean_profile), gene_order$gene_id)
  if (length(common_genes) < 1000) return(NULL)

  go <- gene_order[gene_id %in% common_genes]
  go[, cna_value := mean_profile[gene_id]]

  cna_sd <- sd(go$cna_value, na.rm = TRUE)
  if (!is.na(cna_sd) && cna_sd > 0) {
    go[, cna_value := (cna_value - mean(cna_value, na.rm = TRUE)) / cna_sd * 0.3]
  }

  go <- go[order(genome_mid)]
  window <- min(100L, nrow(go) %/% 5)
  go[, cna_smooth := frollmean(cna_value, n = window, align = "center", na.rm = TRUE)]
  go[is.na(cna_smooth), cna_smooth := cna_value]

  data.table(genome_pos = go$genome_mid, value = go$cna_smooth)
}

# ── 4. Bin profiles for correlation ──────────────────────────────────────────
bin_profile <- function(dt, bin_size = 5e6, is_point_data = FALSE) {
  bins <- data.table(bin_start = seq(0, genome_len - 1, by = bin_size))
  bins[, bin_end := bin_start + bin_size]
  bins[, bin_mid := (bin_start + bin_end) / 2]

  sapply(seq_len(nrow(bins)), function(i) {
    if (is_point_data) {
      pts <- dt[genome_pos >= bins$bin_start[i] & genome_pos < bins$bin_end[i]]
      if (nrow(pts) == 0) NA_real_ else mean(pts$value, na.rm = TRUE)
    } else {
      seg <- dt[genome_start <= bins$bin_mid[i] & genome_end >= bins$bin_mid[i]]
      if (nrow(seg) == 0) NA_real_ else seg$value[1]
    }
  })
}

# ── 5. Compute pairwise correlation matrix ───────────────────────────────────
compute_match_matrix <- function(wes_profiles, scrna_profiles) {
  wes_names <- names(wes_profiles)
  scrna_names <- names(scrna_profiles)

  wes_binned <- lapply(wes_profiles, bin_profile, is_point_data = FALSE)
  scrna_binned <- lapply(scrna_profiles, bin_profile, is_point_data = TRUE)

  mat <- matrix(NA_real_, nrow = length(wes_names), ncol = length(scrna_names),
                dimnames = list(wes_names, scrna_names))

  for (i in seq_along(wes_names)) {
    for (j in seq_along(scrna_names)) {
      w <- wes_binned[[i]]
      s <- scrna_binned[[j]]
      valid <- is.finite(w) & is.finite(s)
      if (sum(valid) >= 20) {
        mat[i, j] <- cor(w[valid], s[valid])
      }
    }
  }

  mat
}

# ── 6. Correlation heatmap ───────────────────────────────────────────────────
make_correlation_heatmap <- function(cor_mat, best_matches = NULL) {
  if (is.null(cor_mat) || all(is.na(cor_mat))) return(NULL)

  cor_dt <- as.data.table(reshape2::melt(cor_mat, na.rm = FALSE))
  setnames(cor_dt, c("wes_clone", "scrna_clone", "correlation"))
  cor_dt[, wes_clone := factor(wes_clone, levels = rownames(cor_mat))]
  cor_dt[, scrna_clone := factor(scrna_clone, levels = colnames(cor_mat))]

  # Mark best matches
  cor_dt[, is_best := FALSE]
  if (!is.null(best_matches)) {
    for (k in seq_len(nrow(best_matches))) {
      cor_dt[wes_clone == best_matches$wes_clone[k] &
               scrna_clone == best_matches$scrna_clone[k],
             is_best := TRUE]
    }
  }

  ggplot(cor_dt, aes(x = scrna_clone, y = wes_clone, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = ifelse(is.finite(correlation),
                                 sprintf("%.3f", correlation), "NA")),
              size = 3.5, fontface = ifelse(cor_dt$is_best, "bold", "plain")) +
    geom_tile(data = cor_dt[is_best == TRUE],
              aes(x = scrna_clone, y = wes_clone),
              fill = NA, color = "#E41A1C", linewidth = 1.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Pearson r", na.value = "grey80") +
    labs(title = "WES Subclone \u2194 scRNA Clone Correlation",
         x = "scRNA Clone", y = "WES Subclone") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 8, 4, 4)
    )
}

# ── 7. Scatter plot for matched pair ─────────────────────────────────────────
make_match_scatter <- function(wes_dt, scrna_dt, wes_label, scrna_label) {
  w_binned <- bin_profile(wes_dt, is_point_data = FALSE)
  s_binned <- bin_profile(scrna_dt, is_point_data = TRUE)
  valid <- is.finite(w_binned) & is.finite(s_binned)
  if (sum(valid) < 20) return(NULL)

  dt <- data.table(wes = w_binned[valid], scrna = s_binned[valid])
  r_val <- cor(dt$wes, dt$scrna)

  ggplot(dt, aes(x = wes, y = scrna)) +
    geom_point(size = 1.2, alpha = 0.4, color = "grey30") +
    geom_smooth(method = "lm", se = TRUE, color = "#D9534F", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("r = ", round(r_val, 3)),
             hjust = 1.2, vjust = 1.5, size = 4, fontface = "bold", color = "#D9534F") +
    labs(x = wes_label, y = scrna_label,
         title = paste0("Best Match: r = ", round(r_val, 3))) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ── 8. Main loop ─────────────────────────────────────────────────────────────
match_summary_rows <- list()
all_cor_matrices <- list()

message("=== WES \u2194 scRNA Subclone Matching ===")

for (wgs_id in names(wgs_sample_map)) {
  scrna_samples <- wgs_sample_map[[wgs_id]]

  # Load WES subclone profiles
  message("\nLoading WES subclone profiles: ", wgs_id)
  wes_data <- load_wes_subclone_cns(wgs_id)
  if (is.null(wes_data)) {
    message("  No WES .cns files. Run Auto_generate_wes_subclone_cns.R first.")
    next
  }

  wes_profiles <- wes_data$profiles
  wes_subclone_profiles <- wes_profiles[!grepl("Bulk", names(wes_profiles))]
  message("  WES profiles loaded: ", paste(names(wes_profiles), collapse = ", "))

  for (scrna_sample in scrna_samples) {
    message("\n  === Matching against: ", scrna_sample, " ===")

    # Load Numbat clones
    numbat_data <- load_numbat_clones(scrna_sample)
    has_numbat <- !is.null(numbat_data) && length(numbat_data$clone_profiles) > 0

    # Load inferCNA
    infercna_ribbon <- tryCatch(
      load_infercna(scrna_sample),
      error = function(e) { message("  inferCNA error: ", e$message); NULL }
    )
    has_infercna <- !is.null(infercna_ribbon)

    if (!has_numbat && !has_infercna) {
      message("  No scRNA clone data available. Skipping.")
      next
    }

    # ── Build scRNA profiles for correlation ──────────────────────────────
    scrna_profiles <- list()
    if (has_numbat) {
      scrna_profiles <- c(scrna_profiles, numbat_data$clone_profiles)
    }
    if (has_infercna) {
      scrna_profiles[["inferCNA Mean"]] <- infercna_ribbon
    }

    # ── Compute correlation matrix ────────────────────────────────────────
    if (length(wes_subclone_profiles) > 0 && length(scrna_profiles) > 0) {
      cor_mat <- compute_match_matrix(wes_subclone_profiles, scrna_profiles)
      all_cor_matrices[[scrna_sample]] <- cor_mat

      # Find best match per WES subclone
      best_matches <- data.table(
        wes_clone = rownames(cor_mat),
        scrna_clone = colnames(cor_mat)[apply(cor_mat, 1, function(x) {
          if (all(is.na(x))) NA_integer_ else which.max(x)
        })],
        correlation = apply(cor_mat, 1, function(x) {
          if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        }),
        sample = scrna_sample,
        wgs_id = wgs_id
      )
      message("  Best matches:")
      for (k in seq_len(nrow(best_matches))) {
        message("    ", best_matches$wes_clone[k], " \u2192 ",
                best_matches$scrna_clone[k],
                " (r = ", round(best_matches$correlation[k], 3), ")")
      }
    } else {
      cor_mat <- NULL
      best_matches <- data.table()
    }

    # Also correlate WES BULK against all scRNA profiles
    if ("WES Bulk" %in% names(wes_profiles) && length(scrna_profiles) > 0) {
      bulk_cor_mat <- compute_match_matrix(
        wes_profiles["WES Bulk"], scrna_profiles
      )
      bulk_best <- data.table(
        wes_clone = "WES Bulk",
        scrna_clone = colnames(bulk_cor_mat)[which.max(bulk_cor_mat[1, ])],
        correlation = max(bulk_cor_mat[1, ], na.rm = TRUE),
        sample = scrna_sample,
        wgs_id = wgs_id
      )
      if (nrow(best_matches) > 0) {
        best_matches <- rbind(best_matches, bulk_best, fill = TRUE)
      } else {
        best_matches <- bulk_best
      }
    }

    match_summary_rows[[scrna_sample]] <- best_matches

    # ── Build visualization ───────────────────────────────────────────────
    panels <- list()
    panel_heights <- c()

    # Panel 1: WES Bulk
    if ("WES Bulk" %in% names(wes_profiles)) {
      panels[["wes_bulk"]] <- make_ribbon_panel(
        wes_profiles[["WES Bulk"]],
        title_str = paste0("WES Bulk CNA \u2014 ", sub("_vs_.*", "", wgs_id)),
        show_chr_labels = FALSE,
        is_point_data = FALSE
      )
      panel_heights <- c(panel_heights, 1)
    }

    # Panel 2: Per-WES-subclone ribbons (faceted)
    if (length(wes_subclone_profiles) > 0) {
      wes_sub_dt <- rbindlist(lapply(names(wes_subclone_profiles), function(nm) {
        dt <- copy(wes_subclone_profiles[[nm]])
        dt[, clone_label := nm]
        dt
      }))
      wes_sub_dt[, clone_label := factor(clone_label, levels = names(wes_subclone_profiles))]

      panels[["wes_subclones"]] <- ggplot() +
        geom_rect(data = chr_bands,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_fill),
                  color = NA, show.legend = FALSE) +
        scale_fill_identity() +
        geom_rect(data = wes_sub_dt,
                  aes(xmin = genome_start, xmax = genome_end, ymin = 0, ymax = pmax(value, 0)),
                  fill = "#D9534F", alpha = 0.6) +
        geom_rect(data = wes_sub_dt,
                  aes(xmin = genome_start, xmax = genome_end, ymin = pmin(value, 0), ymax = 0),
                  fill = "#5B9BD5", alpha = 0.6) +
        geom_segment(data = wes_sub_dt,
                  aes(x = genome_start, xend = genome_end, y = value, yend = value),
                  linewidth = 0.25, color = "grey30") +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
        facet_wrap(~clone_label, ncol = 1, strip.position = "right") +
        chr_x_scale +
        coord_cartesian(ylim = c(-1.5, 1.5)) +
        labs(title = "WES Per-Subclone CNA Profiles (PyClone-VI)", y = "log2 ratio", x = NULL) +
        base_cnv_theme +
        theme(
          strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"),
          strip.background = element_rect(fill = "#FFF3CD", color = NA)
        )
      panel_heights <- c(panel_heights, max(1.5, length(wes_subclone_profiles) * 0.7))
    }

    # Panel 3: Per-Numbat-clone ribbons (faceted)
    if (has_numbat && length(numbat_data$clone_profiles) > 0) {
      numb_dt <- rbindlist(lapply(names(numbat_data$clone_profiles), function(nm) {
        dt <- copy(numbat_data$clone_profiles[[nm]])
        dt[, clone_label := nm]
        dt
      }))
      numb_dt[, clone_label := factor(clone_label, levels = names(numbat_data$clone_profiles))]

      panels[["numbat_clones"]] <- ggplot() +
        geom_rect(data = chr_bands,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_fill),
                  color = NA, show.legend = FALSE) +
        scale_fill_identity() +
        geom_ribbon(data = numb_dt,
                    aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)),
                    fill = "#C0392B", alpha = 0.6) +
        geom_ribbon(data = numb_dt,
                    aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0),
                    fill = "#2980B9", alpha = 0.6) +
        geom_line(data = numb_dt,
                  aes(x = genome_pos, y = value),
                  linewidth = 0.25, color = "grey30") +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
        facet_wrap(~clone_label, ncol = 1, strip.position = "right") +
        chr_x_scale +
        coord_cartesian(ylim = c(-1.5, 1.5)) +
        labs(title = paste0("Numbat Per-Clone CNA Profiles \u2014 ", scrna_sample),
             y = "log2 ratio", x = NULL) +
        base_cnv_theme +
        theme(
          axis.text.x = element_text(size = 7),
          strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"),
          strip.background = element_rect(fill = "#D5E8D4", color = NA)
        )
      n_numb <- length(numbat_data$clone_profiles)
      panel_heights <- c(panel_heights, max(1.5, n_numb * 0.6))
    }

    # Panel 4: inferCNA if available
    if (has_infercna) {
      panels[["infercna"]] <- make_ribbon_panel(
        infercna_ribbon,
        title_str = paste0("inferCNA Mean Profile \u2014 ", scrna_sample),
        show_chr_labels = TRUE,
        gain_col = "#C0392B", loss_col = "#2980B9",
        is_point_data = TRUE
      )
      panel_heights <- c(panel_heights, 1)
    }

    # Panel 5: Correlation heatmap
    heatmap_plot <- NULL
    if (!is.null(cor_mat)) {
      heatmap_plot <- make_correlation_heatmap(
        cor_mat,
        best_matches[wes_clone != "WES Bulk"]
      )
    }

    # Panel 6: Scatter plots for best matches
    scatter_panels <- list()
    if (nrow(best_matches) > 0) {
      sub_matches <- best_matches[wes_clone != "WES Bulk" & is.finite(correlation)]
      for (k in seq_len(nrow(sub_matches))) {
        wes_label <- sub_matches$wes_clone[k]
        scrna_label <- sub_matches$scrna_clone[k]
        if (wes_label %in% names(wes_subclone_profiles) &&
            scrna_label %in% names(scrna_profiles)) {
          sp <- make_match_scatter(
            wes_subclone_profiles[[wes_label]],
            scrna_profiles[[scrna_label]],
            wes_label, scrna_label
          )
          if (!is.null(sp)) scatter_panels[[paste0("sc_", k)]] <- sp
        }
      }
    }

    # ── Assemble final figure ─────────────────────────────────────────────
    ribbon_stack <- Reduce(`/`, panels)
    ribbon_stack <- ribbon_stack + plot_layout(heights = panel_heights)

    bottom_panels <- list()
    bottom_heights <- c()

    if (!is.null(heatmap_plot)) {
      bottom_panels[["heatmap"]] <- heatmap_plot
      bottom_heights <- c(bottom_heights, 2)
    }
    if (length(scatter_panels) > 0) {
      scatter_grid <- wrap_plots(scatter_panels, ncol = min(3, length(scatter_panels)))
      bottom_panels[["scatter"]] <- scatter_grid
      bottom_heights <- c(bottom_heights, 2)
    }

    if (length(bottom_panels) > 0) {
      bottom_stack <- wrap_plots(bottom_panels, ncol = length(bottom_panels))
      final_plot <- ribbon_stack / wrap_elements(full = bottom_stack) +
        plot_layout(heights = c(panel_heights, sum(bottom_heights)))
    } else {
      final_plot <- ribbon_stack
    }

    final_plot <- final_plot +
      plot_annotation(
        title = paste0(scrna_sample, " \u2014 WES Subclone \u2194 scRNA Clone Matching"),
        subtitle = paste0(
          wgs_id, " | ",
          length(wes_subclone_profiles), " WES subclones | ",
          if (has_numbat) paste0(length(numbat_data$clone_profiles), " Numbat clones") else "No Numbat",
          if (has_infercna) " | inferCNA \u2713" else ""
        ),
        theme = theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
        )
      )

    # ── Save ──────────────────────────────────────────────────────────────
    total_h <- sum(panel_heights) * 2.5 + sum(bottom_heights) * 2
    total_h <- max(total_h, 12)

    pdf_path <- file.path(fig_dir, paste0("Auto_wes_scrna_subclone_match_", scrna_sample, ".pdf"))
    png_path <- file.path(fig_dir, paste0("Auto_wes_scrna_subclone_match_", scrna_sample, ".png"))

    ggsave(pdf_path, final_plot, width = 16, height = total_h, limitsize = FALSE)
    ggsave(png_path, final_plot, width = 16, height = total_h, dpi = 200, limitsize = FALSE)
    message("  Saved: ", pdf_path)
  }
}

# ── 9. Summary tables ────────────────────────────────────────────────────────
if (length(match_summary_rows) > 0) {
  match_summary <- rbindlist(match_summary_rows, use.names = TRUE, fill = TRUE)
  match_path <- file.path(table_dir, "Auto_wes_scrna_subclone_match_summary.csv")
  fwrite(match_summary, match_path)
  message("\n=== Match Summary ===")
  print(match_summary)
  message("Wrote: ", match_path)
}

# Save correlation matrices as CSV (wide format, one per sample)
if (length(all_cor_matrices) > 0) {
  cor_rows <- list()
  for (samp in names(all_cor_matrices)) {
    m <- all_cor_matrices[[samp]]
    if (is.null(m)) next
    for (i in seq_len(nrow(m))) {
      for (j in seq_len(ncol(m))) {
        cor_rows[[length(cor_rows) + 1]] <- data.table(
          sample = samp,
          wes_clone = rownames(m)[i],
          scrna_clone = colnames(m)[j],
          correlation = m[i, j]
        )
      }
    }
  }
  if (length(cor_rows) > 0) {
    cor_dt <- rbindlist(cor_rows)
    cor_path <- file.path(table_dir, "Auto_wes_scrna_subclone_correlation_matrix.csv")
    fwrite(cor_dt, cor_path)
    message("Wrote: ", cor_path)
  }
}

# ── Log ───────────────────────────────────────────────────────────────────────
fwrite(data.table(
  finished = as.character(Sys.time()),
  wgs_samples = paste(names(wgs_sample_map), collapse = ","),
  scrna_samples = paste(unlist(wgs_sample_map), collapse = ","),
  n_match_rows = if (exists("match_summary")) nrow(match_summary) else 0L,
  n_cor_matrices = length(all_cor_matrices)
), file.path(log_dir, "Auto_wes_scrna_subclone_match_summary.tsv"), sep = "\t")

message("\nDone.")
