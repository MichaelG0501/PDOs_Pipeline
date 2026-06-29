####################
# Auto_PDO_cnv_compare.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/Auto_PDO_cnv_compare.R
# Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - CNVkit .cns segments from sarek_mutect (WGS bulk)
#   - Numbat bulk_clones_final.tsv.gz per sample (scRNA-seq, per-clone)
#   - inferCNA output: PDOs_outs/cnv/Auto_PDO_infercna_outs_Carroll_2023.rds + meta
#   - Gene order: hg38_gencode_v27.txt
# Outputs:
#   - PDOs_outs/cnv/cnv_compare/Auto_PDO_cnv_compare_<sample>.pdf  (per sample)
#   - PDOs_outs/cnv/cnv_compare/Auto_PDO_cnv_compare_summary.csv
# Downstream use: terminal figures for cross-platform CNV validation
####################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source(file.path(
  "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline",
  "analysis/shared/Auto_pdo_analysis_config.R"
))

setwd(PDO_OUTPUT_DIR)

# ── Configuration ─────────────────────────────────────────────────────────────
cnvkit_root   <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit"
numbat_root   <- file.path(PDO_OUTPUT_DIR, "Auto_PDO_numbat/by_samples")
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"

infercna_outs_path  <- file.path(PDO_OUTPUT_DIR, "cnv/Auto_PDO_infercna_outs_Carroll_2023.rds")
infercna_meta_path  <- file.path(PDO_OUTPUT_DIR, "cnv/Auto_PDO_infercna_meta_Carroll_2023.csv")

out_dir <- "cnv/cnv_compare"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

force_rebuild <- identical(Sys.getenv(PDO_CACHE_ENV$force_rebuild), "1")
replot_only   <- identical(Sys.getenv(PDO_CACHE_ENV$replot_only), "1")

# ── Sample mapping: WGS CNVkit ID → scRNA PDO sample names ───────────────────
# CNVkit uses "PDO_XXXX" while scRNA uses "SURXXXX_{Treated,Untreated}_PDO"
wgs_sample_map <- list(
  "PDO_1070_vs_NT_1070" = c("SUR1070_Untreated_PDO", "SUR1070_Treated_PDO"),
  "PDO_1072_vs_NT_1072" = c("SUR1072_Untreated_PDO", "SUR1072_Treated_PDO"),
  "PDO_1090_vs_NT_1090" = c("SUR1090_Untreated_PDO", "SUR1090_Treated_PDO"),
  "PDO_1121_vs_NT_1121" = c("SUR1121_Untreated_PDO"),
  "PDO_1141_vs_NT_1141" = c("SUR1141_Untreated_PDO"),
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

# ── Helper: expand segment table to paired start/end points for ribbon ────────
expand_seg <- function(dt, val_col, pos_start = "genome_start", pos_end = "genome_end") {
  rbindlist(lapply(seq_len(nrow(dt)), function(i) {
    data.table(
      genome_pos = c(dt[[pos_start]][i], dt[[pos_end]][i]),
      value      = rep(dt[[val_col]][i], 2)
    )
  }))
}

# ── Helper: build a single ribbon panel ───────────────────────────────────────
make_ribbon_panel <- function(ribbon_dt, title_str, y_lab = "log2 ratio",
                              ylim_range = c(-1.5, 1.5), show_chr_labels = FALSE,
                              line_col = "grey20", gain_col = "#D9534F",
                              loss_col = "#5B9BD5", alpha_val = 0.7) {
  p <- ggplot() +
    geom_rect(data = chr_bands,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_fill),
              color = NA, show.legend = FALSE) +
    scale_fill_identity() +
    geom_ribbon(data = ribbon_dt,
                aes(x = genome_pos, ymin = 0, ymax = pmax(value, 0)),
                fill = gain_col, alpha = alpha_val) +
    geom_ribbon(data = ribbon_dt,
                aes(x = genome_pos, ymin = pmin(value, 0), ymax = 0),
                fill = loss_col, alpha = alpha_val) +
    geom_line(data = ribbon_dt,
              aes(x = genome_pos, y = value),
              linewidth = 0.35, color = line_col) +
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

# ── 1. Load CNVkit segments ──────────────────────────────────────────────────
load_cnvkit <- function(wgs_id) {
  cns_path <- file.path(cnvkit_root, wgs_id,
                        paste0(sub("_vs_.*", "", wgs_id), ".cns"))
  if (!file.exists(cns_path)) {
    message("  CNVkit .cns not found: ", cns_path)
    return(NULL)
  }
  cns <- fread(cns_path)
  cns[, chr := sub("^chr", "", chromosome)]
  cns <- cns[chr %in% chr_order]
  cns[, genome_start := start + chr_cumstart[chr]]
  cns[, genome_end   := end   + chr_cumstart[chr]]
  expand_seg(cns, "log2")
}

# ── 2. Load Numbat ───────────────────────────────────────────────────────────
final_iter_from <- function(numbat_dir, prefix = "treeML", ext = "rds") {
  pattern <- paste0(prefix, "_*.", ext)
  files <- Sys.glob(file.path(numbat_dir, pattern))
  if (length(files) == 0) return(NA_integer_)
  regex <- paste0("^", prefix, "_([0-9]+)\\.", gsub("\\.", "\\\\.", ext), "$")
  iter <- suppressWarnings(as.integer(sub(regex, "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

load_numbat <- function(scrna_sample) {
  numbat_dir <- file.path(numbat_root, scrna_sample, "numbat")
  if (!dir.exists(numbat_dir)) {
    message("  Numbat dir missing: ", numbat_dir)
    return(NULL)
  }

  # Find final iteration — try bulk_clones_final first, then highest iteration
  bulk_file <- file.path(numbat_dir, "bulk_clones_final.tsv.gz")
  if (!file.exists(bulk_file)) {
    iter <- final_iter_from(numbat_dir, "bulk_clones", ext = "tsv.gz")
    if (is.finite(iter)) {
      bulk_file <- file.path(numbat_dir, paste0("bulk_clones_", iter, ".tsv.gz"))
    }
  }
  if (!file.exists(bulk_file)) {
    message("  Numbat bulk_clones not found: ", bulk_file)
    return(NULL)
  }

  nb <- fread(bulk_file,
              select = c("CHROM", "gene_start", "gene_end", "n_cells",
                         "members", "phi_mle_roll", "gene_index"))
  nb[, chr := as.character(CHROM)]
  nb <- nb[chr %in% chr_order]
  nb[, genome_mid := (gene_start + gene_end) / 2 + chr_cumstart[chr]]
  nb[, log2_phi := log2(phi_mle_roll)]

  clone_info <- nb[nchar(members) > 0, .(original_n = n_cells[1]), by = members]
  clone_info <- clone_info[order(-original_n)]

  # Check for conservative clones
  conservative_clones_path <- file.path(
    PDO_OUTPUT_DIR, "Auto_PDO_numbat/conservative_clones/by_samples", scrna_sample,
    paste0("Auto_", scrna_sample, "_clones_conservative.rds")
  )
  
  if (file.exists(conservative_clones_path)) {
    message("  Using conservative Numbat clones")
    clones_conservative <- readRDS(conservative_clones_path)
    
    cons_sizes <- as.integer(sapply(clones_conservative, `[[`, "size"))
    cons_names <- as.character(sapply(clones_conservative, `[[`, "sample"))
    cons_df <- data.table(clone_opt = cons_names, cons_n = cons_sizes)
    cons_df <- cons_df[order(-cons_n)]
    
    # Match by rank
    n_match <- min(nrow(clone_info), nrow(cons_df))
    clone_info <- clone_info[1:n_match]
    cons_df <- cons_df[1:n_match]
    
    clone_info[, n_cells := cons_df$cons_n]
    clone_info[, clone_label := paste0("Clone ", cons_df$clone_opt, " (n=", n_cells, ")")]
  } else {
    clone_info[, n_cells := original_n]
    top_n <- min(5, nrow(clone_info))
    clone_info <- clone_info[1:top_n]
    clone_info[, clone_label := paste0("Clone ", seq_len(.N), " (n=", n_cells, ")")]
  }
  
  nb_mapped <- nb[members %in% clone_info$members]
  nb_mapped <- merge(nb_mapped[, -"n_cells", with = FALSE], clone_info[, .(members, n_cells)], by = "members")

  # Pseudo-bulk: weighted mean across matched clones
  nb_bulk <- nb_mapped[, .(
    log2_phi = weighted.mean(log2_phi, w = n_cells, na.rm = TRUE)
  ), by = .(genome_mid)]
  nb_bulk <- nb_bulk[order(genome_mid)]

  # Per-clone profiles
  top_n <- nrow(clone_info)
  top_members <- clone_info$members
  top_labels <- clone_info$clone_label
  clone_profiles <- rbindlist(lapply(seq_len(top_n), function(i) {
    d <- nb_mapped[members == top_members[i]]
    d[, clone_label := top_labels[i]]
    d[order(genome_mid)]
  }))
  clone_profiles$clone_label <- factor(clone_profiles$clone_label, levels = top_labels)

  list(
    bulk_ribbon = data.table(genome_pos = nb_bulk$genome_mid, value = nb_bulk$log2_phi),
    clone_df    = clone_profiles,
    clone_info  = clone_info
  )
}

# ── 3. Load inferCNA ─────────────────────────────────────────────────────────
load_infercna <- function(scrna_sample) {
  if (!file.exists(infercna_outs_path) || !file.exists(infercna_meta_path)) {
    message("  inferCNA outputs not found. Skipping. Run Auto_PDO_infercna.R first.")
    return(NULL)
  }

  # Load gene order for genomic positions
  gene_order <- as.data.table(read.table(gene_order_path, header = FALSE,
                                         col.names = c("gene_id", "chromosome", "start", "end")))
  gene_order[, chr := sub("^chr", "", chromosome)]
  gene_order <- gene_order[chr %in% chr_order]
  gene_order[, chr := factor(chr, levels = chr_order)]
  gene_order <- gene_order[order(chr, start)]
  gene_order[, genome_mid := (start + end) / 2 + chr_cumstart[as.character(chr)]]

  # Load inferCNA matrix (genes × cells) and meta
  outs <- readRDS(infercna_outs_path)
  meta <- read.csv(infercna_meta_path)

  # Select cells for this sample
  sample_cells <- meta$cell[meta$sample == scrna_sample]
  if (length(sample_cells) == 0) {
    message("  No inferCNA cells for sample: ", scrna_sample)
    return(NULL)
  }
  sample_cells <- intersect(sample_cells, colnames(outs))
  if (length(sample_cells) < 10) {
    message("  Too few inferCNA cells (", length(sample_cells), ") for: ", scrna_sample)
    return(NULL)
  }

  # Per-sample mean CNA profile
  sample_mat <- outs[, sample_cells, drop = FALSE]
  mean_profile <- rowMeans(sample_mat, na.rm = TRUE)

  # Map genes to genomic position
  common_genes <- intersect(names(mean_profile), gene_order$gene_id)
  if (length(common_genes) < 1000) {
    message("  Too few genes with positions (", length(common_genes), ") for: ", scrna_sample)
    return(NULL)
  }

  go <- gene_order[gene_id %in% common_genes]
  go[, cna_value := mean_profile[gene_id]]
  
  # Normalize inferCNA signal strength to be comparable with Numbat/CNVkit log2 ratios
  cna_sd <- sd(go$cna_value, na.rm = TRUE)
  if (!is.na(cna_sd) && cna_sd > 0) {
    go[, cna_value := (cna_value - mean(cna_value, na.rm = TRUE)) / cna_sd * 0.3]
  }
  
  go <- go[order(genome_mid)]

  # Smooth with rolling mean for cleaner visualisation (window = 100 genes)
  window <- min(100L, nrow(go) %/% 5)
  go[, cna_smooth := frollmean(cna_value, n = window, align = "center", na.rm = TRUE)]
  go[is.na(cna_smooth), cna_smooth := cna_value]

  data.table(genome_pos = go$genome_mid, value = go$cna_smooth)
}

# ── 4. Per-chromosome segment-level concordance ─────────────────────────────
compute_concordance <- function(cnvkit_ribbon, numbat_ribbon, infercna_ribbon) {
  # Bin genome into fixed windows and compute per-window mean from each source
  bin_size <- 5e6  # 5 Mb bins
  bins <- data.table(
    bin_start = seq(0, genome_len - 1, by = bin_size)
  )
  bins[, bin_end := bin_start + bin_size]
  bins[, bin_mid := (bin_start + bin_end) / 2]

  # Assign chromosome to each bin
  bins[, chr := {
    chr_assign <- character(.N)
    for (ch in chr_order) {
      chr_assign[bin_mid >= chr_cumstart[ch] & bin_mid < chr_cumend[ch]] <- ch
    }
    chr_assign
  }]
  bins <- bins[nchar(chr) > 0]

  # Mean value per bin for each source
  bin_mean <- function(ribbon_dt, bins_dt) {
    sapply(seq_len(nrow(bins_dt)), function(i) {
      pts <- ribbon_dt[genome_pos >= bins_dt$bin_start[i] & genome_pos < bins_dt$bin_end[i]]
      if (nrow(pts) == 0) NA_real_ else mean(pts$value, na.rm = TRUE)
    })
  }

  result <- copy(bins)

  if (!is.null(cnvkit_ribbon))  result[, cnvkit  := bin_mean(cnvkit_ribbon, bins)]
  if (!is.null(numbat_ribbon))  result[, numbat  := bin_mean(numbat_ribbon, bins)]
  if (!is.null(infercna_ribbon)) result[, infercna := bin_mean(infercna_ribbon, bins)]

  result
}

# ── 5. Correlation scatter panel ─────────────────────────────────────────────
make_scatter_panel <- function(binned, x_col, y_col, x_lab, y_lab) {
  dt <- binned[is.finite(get(x_col)) & is.finite(get(y_col))]
  if (nrow(dt) < 10) return(NULL)

  r_val <- cor(dt[[x_col]], dt[[y_col]], use = "complete.obs")

  ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(size = 1, alpha = 0.4, color = "grey30") +
    geom_smooth(method = "lm", se = TRUE, color = "#D9534F", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("r = ", round(r_val, 3)),
             hjust = 1.2, vjust = 1.5, size = 4, fontface = "bold", color = "#D9534F") +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = 11) +
    theme(
      plot.margin = margin(4, 4, 4, 4),
      axis.title = element_text(face = "bold")
    )
}

# ── 6. Confusion matrix panel ────────────────────────────────────────────────
make_confusion_panel <- function(binned, x_col, y_col, x_lab, y_lab, threshold = 0.05) {
  dt <- binned[is.finite(get(x_col)) & is.finite(get(y_col))]
  if (nrow(dt) < 10) return(NULL)

  classify <- function(vals, thr) {
    ifelse(vals > thr, "gain", ifelse(vals < -thr, "loss", "neutral"))
  }

  dt[, truth := classify(get(x_col), threshold)]
  dt[, predicted := classify(get(y_col), threshold)]

  levels_use <- c("loss", "neutral", "gain")
  dt[, truth     := factor(truth,     levels = levels_use)]
  dt[, predicted := factor(predicted, levels = levels_use)]

  conf <- as.data.frame(table(Truth = dt$truth, Predicted = dt$predicted))
  conf$frac <- ave(conf$Freq, conf$Truth, FUN = function(x) x / sum(x))
  conf$label <- paste0(conf$Freq, "\n(", round(conf$frac * 100, 1), "%)")

  ggplot(conf, aes(x = Truth, y = Predicted, fill = frac)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = label), size = 3.5) +
    scale_fill_gradient(low = "white", high = "#E8C547", guide = "none") +
    labs(x = paste0("Truth (", x_lab, ")"), y = paste0("Predicted (", y_lab, ")")) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ── 7. Main loop ─────────────────────────────────────────────────────────────
summary_rows <- list()
message("=== CNV Cross-Platform Comparison ===")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  wgs_id_arg <- args[1]
  if (wgs_id_arg %in% names(wgs_sample_map)) {
    wgs_sample_map <- wgs_sample_map[wgs_id_arg]
  } else {
    # If the user passed a scRNA sample directly (like SUR1072_Untreated_PDO)
    # Find which wgs_id it belongs to and only keep that sample
    found <- FALSE
    for (wgs_id in names(wgs_sample_map)) {
      if (wgs_id_arg %in% wgs_sample_map[[wgs_id]]) {
        wgs_sample_map <- list()
        wgs_sample_map[[wgs_id]] <- wgs_id_arg
        found <- TRUE
        break
      }
    }
    if (!found) stop("Argument must be a wgs_id or scrna_sample name from wgs_sample_map.")
  }
}

for (wgs_id in names(wgs_sample_map)) {
  scrna_samples <- wgs_sample_map[[wgs_id]]

  # Load CNVkit once per WGS pair
  message("\nLoading CNVkit: ", wgs_id)
  cnvkit_ribbon <- load_cnvkit(wgs_id)

  for (scrna_sample in scrna_samples) {
    message("  Processing scRNA sample: ", scrna_sample)

    # Load Numbat
    numbat_data <- load_numbat(scrna_sample)

    # Load inferCNA
    infercna_ribbon <- tryCatch(
      load_infercna(scrna_sample),
      error = function(e) { message("  inferCNA error: ", e$message); NULL }
    )

    has_cnvkit  <- !is.null(cnvkit_ribbon)
    has_numbat  <- !is.null(numbat_data)
    has_infercna <- !is.null(infercna_ribbon)

    if (!has_cnvkit && !has_numbat && !has_infercna) {
      message("  No data available for ", scrna_sample, ". Skipping.")
      summary_rows[[scrna_sample]] <- data.frame(
        sample = scrna_sample, wgs_id = wgs_id,
        has_cnvkit = FALSE, has_numbat = FALSE, has_infercna = FALSE,
        r_cnvkit_numbat = NA, r_cnvkit_infercna = NA, r_numbat_infercna = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    # ── Build ribbon panels ───────────────────────────────────────────────────
    panels <- list()
    panel_heights <- c()

    if (has_cnvkit) {
      panels[["cnvkit"]] <- make_ribbon_panel(
        cnvkit_ribbon,
        title_str = paste0("CNVkit \u2014 WES Bulk (", sub("_vs_.*", "", wgs_id), ")"),
        show_chr_labels = TRUE
      )
      panel_heights <- c(panel_heights, 1)
    }

    if (has_infercna) {
      panels[["infercna"]] <- make_ribbon_panel(
        infercna_ribbon,
        title_str = paste0("inferCNA \u2014 scRNA Mean Profile (", scrna_sample, ")"),
        show_chr_labels = TRUE,
        gain_col = "#C0392B", loss_col = "#2980B9"
      )
      panel_heights <- c(panel_heights, 1)
    }

    if (has_numbat) {
      panels[["numbat_bulk"]] <- make_ribbon_panel(
        numbat_data$bulk_ribbon,
        title_str = paste0("Numbat \u2014 scRNA Pseudo-bulk (", scrna_sample, ")"),
        show_chr_labels = TRUE
      )
      panel_heights <- c(panel_heights, 1)

      # Per-clone faceted panel
      clone_df <- numbat_data$clone_df
      panels[["numbat_clones"]] <- ggplot() +
        geom_rect(data = chr_bands,
                  aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band_fill),
                  color = NA, show.legend = FALSE) +
        scale_fill_identity() +
        geom_ribbon(data = clone_df,
                    aes(x = genome_mid, ymin = 0, ymax = pmax(log2_phi, 0)),
                    fill = "#D9534F", alpha = 0.6) +
        geom_ribbon(data = clone_df,
                    aes(x = genome_mid, ymin = pmin(log2_phi, 0), ymax = 0),
                    fill = "#5B9BD5", alpha = 0.6) +
        geom_line(data = clone_df,
                  aes(x = genome_mid, y = log2_phi),
                  linewidth = 0.25, color = "grey30") +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
        facet_wrap(~clone_label, ncol = 1, strip.position = "right") +
        chr_x_scale +
        coord_cartesian(ylim = c(-1.5, 1.5)) +
        labs(title = "Numbat \u2014 Per-Clone CNV Profiles", y = "log2 ratio", x = NULL) +
        base_cnv_theme +
        theme(
          axis.text.x = element_text(size = 7),
          strip.text.y.right = element_text(size = 8, angle = 0, face = "bold"),
          strip.background = element_rect(fill = "grey95", color = NA)
        )
      n_clones <- nlevels(clone_df$clone_label)
      panel_heights <- c(panel_heights, max(1.5, n_clones * 0.6))
    }

    # ── Concordance metrics ───────────────────────────────────────────────────
    numbat_bulk_ribbon <- if (has_numbat) numbat_data$bulk_ribbon else NULL
    binned <- compute_concordance(cnvkit_ribbon, numbat_bulk_ribbon, infercna_ribbon)

    r_cnvkit_numbat  <- NA_real_
    r_cnvkit_infercna <- NA_real_
    r_numbat_infercna <- NA_real_

    scatter_panels <- list()

    if (has_cnvkit && has_numbat && "cnvkit" %in% names(binned) && "numbat" %in% names(binned)) {
      r_cnvkit_numbat <- cor(binned$cnvkit, binned$numbat, use = "complete.obs")
      scatter_panels[["cn"]] <- make_scatter_panel(binned, "cnvkit", "numbat",
                                                    "CNVkit (WES)", "Numbat (scRNA)")
      scatter_panels[["cn_conf"]] <- make_confusion_panel(binned, "cnvkit", "numbat",
                                                           "CNVkit", "Numbat")
    }
    if (has_cnvkit && has_infercna && "cnvkit" %in% names(binned) && "infercna" %in% names(binned)) {
      r_cnvkit_infercna <- cor(binned$cnvkit, binned$infercna, use = "complete.obs")
      scatter_panels[["ci"]] <- make_scatter_panel(binned, "cnvkit", "infercna",
                                                    "CNVkit (WES)", "inferCNA (scRNA)")
      scatter_panels[["ci_conf"]] <- make_confusion_panel(binned, "cnvkit", "infercna",
                                                           "CNVkit", "inferCNA")
    }
    if (has_numbat && has_infercna && "numbat" %in% names(binned) && "infercna" %in% names(binned)) {
      r_numbat_infercna <- cor(binned$numbat, binned$infercna, use = "complete.obs")
      scatter_panels[["ni"]] <- make_scatter_panel(binned, "numbat", "infercna",
                                                    "Numbat (scRNA)", "inferCNA (scRNA)")
      scatter_panels[["ni_conf"]] <- make_confusion_panel(binned, "numbat", "infercna",
                                                           "Numbat", "inferCNA")
    }

    # ── Assemble final figure ─────────────────────────────────────────────────
    # Stack ribbon panels vertically
    ribbon_stack <- Reduce(`/`, panels)
    ribbon_stack <- ribbon_stack + plot_layout(heights = panel_heights)

    # If scatter panels exist, combine them using wrap_elements to avoid
    # patchwork flattening nested layouts
    scatter_panels <- Filter(Negate(is.null), scatter_panels)
    scatter_height <- 0
    if (length(scatter_panels) > 0) {
      scatter_grid <- wrap_plots(scatter_panels, ncol = 2)
      scatter_height <- max(2, length(scatter_panels) * 0.6)
      final_plot <- ribbon_stack / wrap_elements(full = scatter_grid) +
        plot_layout(heights = c(panel_heights, scatter_height))
    } else {
      final_plot <- ribbon_stack
    }

    final_plot <- final_plot +
      plot_annotation(
        title = paste0(scrna_sample, " \u2014 Cross-Platform CNV Comparison"),
        subtitle = paste0(
          "CNVkit (WES bulk)", if (has_cnvkit) " \u2713" else " \u2717",
          " | Numbat (scRNA)", if (has_numbat) paste0(" \u2713 (", nrow(numbat_data$clone_info), " clones)") else " \u2717",
          " | inferCNA (scRNA)", if (has_infercna) " \u2713" else " \u2717"
        ),
        theme = theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
        )
      )

    # ── Save ──────────────────────────────────────────────────────────────────
    total_h <- sum(panel_heights) * 2.5 + if (length(scatter_panels) > 0) scatter_height * 2.5 else 0
    pdf_path <- file.path(out_dir, paste0("Auto_PDO_cnv_compare_", scrna_sample, ".pdf"))
    png_path <- file.path(out_dir, paste0("Auto_PDO_cnv_compare_", scrna_sample, ".png"))

    ggsave(pdf_path, final_plot, width = 16, height = total_h, limitsize = FALSE)
    ggsave(png_path, final_plot, width = 16, height = total_h, dpi = 200, limitsize = FALSE)
    message("  Saved: ", pdf_path)

    # ── Save binned concordance table ─────────────────────────────────────────
    fwrite(binned, file.path(out_dir, paste0("Auto_PDO_cnv_compare_binned_", scrna_sample, ".csv")))

    summary_rows[[scrna_sample]] <- data.frame(
      sample = scrna_sample,
      wgs_id = wgs_id,
      has_cnvkit = has_cnvkit,
      has_numbat = has_numbat,
      has_infercna = has_infercna,
      n_numbat_clones = if (has_numbat) nrow(numbat_data$clone_info) else NA_integer_,
      r_cnvkit_numbat = r_cnvkit_numbat,
      r_cnvkit_infercna = r_cnvkit_infercna,
      r_numbat_infercna = r_numbat_infercna,
      stringsAsFactors = FALSE
    )
  }
}

# ── 8. Summary table ─────────────────────────────────────────────────────────
summary_df <- bind_rows(summary_rows)
summary_path <- file.path(out_dir, "Auto_PDO_cnv_compare_summary.csv")
fwrite(summary_df, summary_path)
message("\n=== Summary ===")
print(summary_df)
message("Wrote: ", file.path(PDO_OUTPUT_DIR, summary_path))
message("Done.")
