####################
# Auto_PDO_infercna.R
#
# Run InferCNA on PDO single-cell RNA-seq samples with the Carroll 2023
# non-malignant scATLAS reference used by Auto_parse_infercna.R, then draw
# one all-sample CNV-profile heatmap.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(infercna)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(scales)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

args <- commandArgs(trailingOnly = TRUE)
sample_arg <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "all"
max_cells_per_sample <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 750L

out_dir <- "cnv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reference_name <- "Carroll_2023"
reference_path <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds"
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
excluded_samples <- c("SUR843T3_PDO")

if (!file.exists(reference_path)) stop("Missing reference RDS: ", reference_path)
if (!file.exists(gene_order_path)) stop("Missing gene order file: ", gene_order_path)

sample_dirs <- list.dirs("by_samples", recursive = FALSE, full.names = FALSE)
sample_dirs <- sample_dirs[grepl("_PDO$", sample_dirs)]
sample_dirs <- sample_dirs[!sample_dirs %in% excluded_samples]
sample_files <- file.path("by_samples", sample_dirs, paste0(sample_dirs, ".rds"))
names(sample_files) <- sample_dirs
sample_files <- sample_files[file.exists(sample_files)]

if (!identical(sample_arg, "all")) {
  requested <- trimws(unlist(strsplit(sample_arg, ",")))
  sample_files <- sample_files[names(sample_files) %in% requested]
}
if (length(sample_files) == 0) stop("No PDO sample RDS files found for argument: ", sample_arg)

sample_order <- names(sample_files)
message("PDO samples: ", paste(sample_order, collapse = ", "))

get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

to_cpm <- function(counts) {
  lib_size <- Matrix::colSums(counts)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- 1
  cpm <- Matrix::t(Matrix::t(counts) * (1e6 / lib_size))
  dimnames(cpm) <- dimnames(counts)
  cpm
}

message("Loading reference: ", reference_path)
reference <- readRDS(reference_path)
ref_prefix <- paste0(reference_name, "_ref")

outs_path <- file.path(out_dir, "Auto_PDO_infercna_outs_Carroll_2023.rds")
target_outs_path <- file.path(out_dir, "Auto_PDO_infercna_target_outs_Carroll_2023.rds")
target_meta_path <- file.path(out_dir, "Auto_PDO_infercna_target_meta_Carroll_2023.rds")
meta_path <- file.path(out_dir, "Auto_PDO_infercna_meta_Carroll_2023.csv")

if (file.exists(outs_path) && file.exists(target_meta_path) && file.exists(meta_path)) {
  message("Existing InferCNA outputs found. Loading to replot...")
  outs <- readRDS(outs_path)
  target_meta <- readRDS(target_meta_path)
  combined_meta <- read.csv(meta_path)
  rownames(combined_meta) <- combined_meta$cell
  target_outs <- outs[, target_meta$cell, drop = FALSE]
} else {
  message("Finding common genes across PDO samples and reference")
  sample_genes <- lapply(sample_files, function(path) rownames(readRDS(path)))
  common_genes <- Reduce(intersect, c(list(rownames(reference$matrix)), sample_genes))
  if (length(common_genes) < 5000) stop("Too few common genes for InferCNA: ", length(common_genes))
  message("Common genes: ", length(common_genes))

  target_mats <- list()
  target_meta <- list()

  for (sample in sample_order) {
    message("Preparing PDO CPM matrix for ", sample)
    obj <- readRDS(sample_files[[sample]])
    counts <- get_counts(obj)[common_genes, , drop = FALSE]
    cpm <- to_cpm(counts)
    new_cells <- paste(sample, colnames(cpm), sep = "__")
    target_mats[[sample]] <- cpm
    colnames(target_mats[[sample]]) <- new_cells
    target_meta[[sample]] <- data.frame(
      cell = new_cells,
      original_cell = colnames(obj),
      sample = sample,
      compartment = "PDO_target",
      stringsAsFactors = FALSE
    )
    rm(obj, counts, cpm)
    gc()
  }

  target_matrix <- do.call(cbind, target_mats)
  target_meta <- bind_rows(target_meta)
  rownames(target_meta) <- target_meta$cell
  rm(target_mats)
  gc()

  reference_matrix <- reference$matrix[common_genes, , drop = FALSE]
  original_ref_cells <- colnames(reference_matrix)
  colnames(reference_matrix) <- paste(ref_prefix, original_ref_cells, sep = "__")
  ref_cells <- lapply(reference$ref, function(cells) paste(ref_prefix, cells, sep = "__"))
  ref_cells <- lapply(ref_cells, intersect, colnames(reference_matrix))
  ref_cells <- ref_cells[lengths(ref_cells) > 0]
  if (length(ref_cells) < 1) stop("No valid reference cell groups after renaming reference cells")

  reference_meta <- reference$meta[match(original_ref_cells, rownames(reference$meta)), , drop = FALSE]
  reference_meta$cell <- colnames(reference_matrix)
  reference_meta$original_cell <- original_ref_cells
  reference_meta$sample <- ref_prefix
  reference_meta$compartment <- paste0("reference_", reference_meta$celltype_update)
  rownames(reference_meta) <- reference_meta$cell

  message("Combining PDO target and Carroll reference matrices")
  combined_matrix <- cbind(as.matrix(target_matrix), as.matrix(reference_matrix))
  rm(target_matrix, reference_matrix)
  gc()

  combined_meta <- bind_rows(
    target_meta,
    reference_meta[, union(colnames(target_meta), colnames(reference_meta)), drop = FALSE]
  )
  rownames(combined_meta) <- combined_meta$cell

  message("Running infercna()")
  outs <- infercna(
    combined_matrix,
    refCells = ref_cells,
    isLog = FALSE,
    verbose = TRUE
  )

  target_outs <- outs[, target_meta$cell, drop = FALSE]

  saveRDS(outs, outs_path)
  saveRDS(target_outs, target_outs_path)

  message("Calculating target CNA scatter metrics for downstream QC summaries")
  target_coord <- as.data.frame(cnaScatterPlot(target_outs))
  target_meta$cna_signal <- target_coord$cna.signal[match(target_meta$cell, rownames(target_coord))]
  target_meta$cna_cor <- target_coord$cna.cor[match(target_meta$cell, rownames(target_coord))]
  combined_meta[target_meta$cell, "cna_signal"] <- target_meta$cna_signal
  combined_meta[target_meta$cell, "cna_cor"] <- target_meta$cna_cor

  saveRDS(target_meta, target_meta_path)
  write.csv(combined_meta, meta_path, row.names = FALSE)

  for (sample in sample_order) {
    sample_out_dir <- file.path(out_dir, "by_samples", sample)
    dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
    sample_cells <- target_meta$cell[target_meta$sample == sample]
    saveRDS(
      target_outs[, sample_cells, drop = FALSE],
      file.path(sample_out_dir, paste0("Auto_", sample, "_infercna_outs_Carroll_2023.rds"))
    )
  }

  writeLines(
    c(
      paste0("reference_name\t", reference_name),
      paste0("reference_path\t", reference_path),
      paste0("reference_groups\t", paste(names(ref_cells), collapse = ",")),
      paste0("reference_cells\t", length(unlist(ref_cells, use.names = FALSE))),
      paste0("common_genes\t", if (exists("common_genes")) length(common_genes) else "NA"),
      paste0("target_cells\t", nrow(target_meta)),
      paste0("sample_count\t", length(sample_order)),
      paste0("infercna_outs\t", outs_path),
      paste0("target_outs\t", target_outs_path),
      paste0("target_meta\t", target_meta_path)
    ),
    con = file.path(out_dir, "Auto_PDO_infercna_reference_summary.tsv")
  )
}

bin_cna_matrix <- function(cna_mat, gene_order, bin_size = 50L) {
  chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
  common <- intersect(rownames(cna_mat), gene_order$gene_id)
  go <- gene_order %>%
    filter(.data$gene_id %in% common, .data$chromosome %in% chrom_levels) %>%
    mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
    arrange(.data$chromosome, .data$start)
  cna_mat <- cna_mat[go$gene_id, , drop = FALSE]
  go <- go %>%
    group_by(.data$chromosome) %>%
    mutate(
      g_rank = row_number(),
      bin_in_chr = ((.data$g_rank - 1L) %/% bin_size) + 1L,
      bin_key = paste(.data$chromosome, .data$bin_in_chr, sep = "_")
    ) %>%
    ungroup()
  bins_idx <- split(seq_len(nrow(go)), factor(go$bin_key, levels = unique(go$bin_key)))
  binned <- do.call(rbind, lapply(bins_idx, function(ix) colMeans(cna_mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins_idx)
  list(
    matrix = binned,
    row_chr = factor(sub("_.*$", "", rownames(binned)), levels = chrom_levels),
    row_chr_labels = sub("_.*$", "", rownames(binned))
  )
}

make_cna_scale <- function(mat) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  nz <- vals[vals != 0]
  scale_vals <- if (length(nz) >= 10) nz else vals
  lim <- as.numeric(quantile(abs(scale_vals), 0.985, na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- max(abs(scale_vals), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 0.1
  lim <- min(max(lim, 0.05), 1.5)
  list(
    fun = colorRamp2(c(-lim, 0, lim), c("#053061", "white", "#67001F")),
    lim = lim
  )
}

plot_cna_heatmap <- function(cna_mat, meta_df, samples_use, output_file, max_cells_per_sample) {
  set.seed(123)
  cells_use <- unlist(lapply(samples_use, function(s) {
    x <- rownames(meta_df)[meta_df$sample == s]
    if (length(x) > max_cells_per_sample) sample(x, max_cells_per_sample) else x
  }), use.names = FALSE)
  cells_use <- intersect(cells_use, colnames(cna_mat))
  cna_mat <- cna_mat[, cells_use, drop = FALSE]
  meta_plot <- meta_df[colnames(cna_mat), , drop = FALSE]

  gene_order <- read.table(
    gene_order_path,
    header = FALSE,
    col.names = c("gene_id", "chromosome", "start", "end")
  )
  binned <- bin_cna_matrix(cna_mat, gene_order)
  mat <- binned$matrix
  row_chr <- binned$row_chr
  row_chr_labels <- binned$row_chr_labels

  sample_factor <- factor(meta_plot$sample, levels = samples_use)
  sample_cols <- setNames(hcl.colors(length(samples_use), palette = "Dark 3"), samples_use)
  top_ha <- HeatmapAnnotation(
    sample = sample_factor,
    col = list(sample = sample_cols),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 14, fontface = "bold"),
    show_legend = FALSE,
    simple_anno_size = unit(5, "mm")
  )

  chr_used <- levels(droplevels(row_chr))
  base_cols <- c(
    brewer.pal(12, "Paired"),
    brewer.pal(8, "Dark2"),
    brewer.pal(9, "Set1"),
    brewer.pal(12, "Set3")
  )
  chr_cols <- setNames(base_cols[seq_along(chr_used)], chr_used)
  left_chr_bar <- rowAnnotation(
    chr = row_chr,
    col = list(chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    gp = gpar(col = NA),
    width = unit(4, "mm")
  )

  cna_scale <- make_cna_scale(mat)
  chr_bounds <- which(head(row_chr_labels, -1L) != tail(row_chr_labels, -1L))
  line_gp <- gpar(col = "black", lwd = 1.5, lineend = "square")

  ht <- Heatmap(
    mat,
    name = "InferCNA",
    col = cna_scale$fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    column_split = sample_factor,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = top_ha,
    left_annotation = left_chr_bar,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    row_title_rot = 0,
    column_title_rot = 30,
    column_gap = unit(2, "mm"),
    rect_gp = gpar(col = NA),
    border = NA,
    use_raster = TRUE,
    heatmap_legend_param = list(
      title = "CNA",
      title_gp = gpar(fontsize = 13, fontface = "bold"),
      labels_gp = gpar(fontsize = 11)
    ),
    layer_fun = function(j, i, x, y, w, h, fill) {
      hits <- intersect(i, chr_bounds)
      if (length(hits)) {
        id <- match(hits, i)
        yy <- y[id] - h[id] / 2
        grid.segments(
          x0 = unit(0, "npc"), x1 = unit(1, "npc"),
          y0 = yy, y1 = yy,
          gp = line_gp
        )
      }
    }
  )

  saveRDS(
    list(
      binned_matrix = mat,
      metadata = meta_plot,
      samples = samples_use,
      colour_clip = cna_scale$lim
    ),
    sub("\\.pdf$", "_input.rds", output_file)
  )
  pdf(output_file, width = max(12, length(samples_use) * 1.15), height = 9, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

plot_cna_heatmap(
  cna_mat = target_outs,
  meta_df = target_meta,
  samples_use = sample_order,
  output_file = file.path(out_dir, "Auto_PDO_cnv_heatmap_all_samples_Carroll_2023.pdf"),
  max_cells_per_sample = max_cells_per_sample
)

plot_cna_scatter <- function(outs, meta_df, ref_prefix, sample_order, output_file) {
  coord <- cnaScatterPlot(outs)
  scatter_df <- as.data.frame(coord)
  
  # Focus on concentrated points location by capping X-axis (signal) at 99.5th percentile
  target_vals <- scatter_df$cna.signal[meta_df[rownames(scatter_df), "sample"] != ref_prefix]
  x_lim <- quantile(target_vals, 0.995, na.rm = TRUE)
  message("Scatter plot X-axis limit (99.5th percentile signal): ", x_lim)
  scatter_df$cell <- rownames(scatter_df)
  scatter_df$sample <- meta_df[scatter_df$cell, "sample"]
  scatter_df$is_reference <- scatter_df$sample == ref_prefix
  scatter_df$plot_group <- ifelse(scatter_df$is_reference, "Reference", scatter_df$sample)
  scatter_df$plot_group <- factor(scatter_df$plot_group, levels = c("Reference", sample_order))
  
  ref_df <- scatter_df[scatter_df$is_reference, , drop = FALSE]
  thr_cor <- mean(ref_df$cna.cor, na.rm = TRUE) + 2 * sd(ref_df$cna.cor, na.rm = TRUE)
  thr_sig <- mean(ref_df$cna.signal, na.rm = TRUE) + 2 * sd(ref_df$cna.signal, na.rm = TRUE)
  if (!is.finite(thr_cor)) thr_cor <- max(ref_df$cna.cor, na.rm = TRUE)
  if (!is.finite(thr_sig)) thr_sig <- max(ref_df$cna.signal, na.rm = TRUE)
  if (!is.finite(thr_cor)) thr_cor <- 0
  if (!is.finite(thr_sig)) thr_sig <- 0
  
  scatter_summary <- scatter_df %>%
    mutate(above_reference = cna.cor > thr_cor & cna.signal > thr_sig) %>%
    group_by(sample) %>%
    summarise(
      n_cells = n(),
      median_cna_cor = median(cna.cor, na.rm = TRUE),
      median_cna_signal = median(cna.signal, na.rm = TRUE),
      q95_cna_cor = quantile(cna.cor, 0.95, na.rm = TRUE),
      q95_cna_signal = quantile(cna.signal, 0.95, na.rm = TRUE),
      pct_above_reference = 100 * mean(above_reference, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(scatter_summary, sub("\\.pdf$", "_summary.csv", output_file), row.names = FALSE)
  
  scatter_cols <- c("Reference" = "#D9D9D9", setNames(hcl.colors(length(sample_order), palette = "Dark 3"), sample_order))
  target_scatter <- scatter_df[!scatter_df$is_reference, , drop = FALSE]
  ref_scatter <- scatter_df[scatter_df$is_reference, , drop = FALSE]
  set.seed(1)
  if (nrow(ref_scatter) > 3000) ref_scatter <- ref_scatter[sample(seq_len(nrow(ref_scatter)), 3000), , drop = FALSE]
  if (nrow(target_scatter) > 30000) target_scatter <- target_scatter[sample(seq_len(nrow(target_scatter)), 30000), , drop = FALSE]
  
  p_all <- ggplot() +
    geom_point(
      data = ref_scatter,
      aes(x = cna.signal, y = cna.cor),
      color = "#D9D9D9",
      size = 0.5,
      alpha = 0.45,
      position = position_jitter(width = 0.000002, height = 0.002)
    ) +
    geom_point(
      data = target_scatter,
      aes(x = cna.signal, y = cna.cor, color = plot_group),
      size = 0.7,
      alpha = 0.65
    ) +
    geom_vline(xintercept = thr_sig, linetype = "dashed", color = "grey30", linewidth = 0.45) +
    geom_hline(yintercept = thr_cor, linetype = "dashed", color = "grey30", linewidth = 0.45) +
    scale_color_manual(values = scatter_cols, drop = FALSE, name = "Sample") +
    scale_x_continuous(labels = scientific) +
    coord_cartesian(xlim = c(0, x_lim)) +
    labs(
      title = "InferCNA Scatter",
      x = "CNA signal",
      y = "CNA correlation"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 15),
      legend.text = element_text(size = 14),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  ref_facet <- do.call(rbind, lapply(sample_order, function(s) {
    tmp <- ref_scatter
    tmp$facet_sample <- s
    tmp
  }))
  target_facet <- scatter_df[!scatter_df$is_reference & scatter_df$sample %in% sample_order, , drop = FALSE]
  target_facet$facet_sample <- target_facet$sample
  facet_df <- rbind(ref_facet, target_facet)
  facet_df$facet_sample <- factor(facet_df$facet_sample, levels = sample_order)
  
  p_facet <- ggplot() +
    geom_point(
      data = facet_df[facet_df$is_reference, , drop = FALSE],
      aes(x = cna.signal, y = cna.cor),
      color = "#D9D9D9",
      size = 0.3,
      alpha = 0.35,
      position = position_jitter(width = 0.000002, height = 0.002)
    ) +
    geom_point(
      data = facet_df[!facet_df$is_reference, , drop = FALSE],
      aes(x = cna.signal, y = cna.cor, color = sample),
      size = 0.45,
      alpha = 0.55
    ) +
    facet_wrap(~facet_sample, ncol = 3) +
    geom_vline(xintercept = thr_sig, linetype = "dashed", color = "grey30", linewidth = 0.3) +
    geom_hline(yintercept = thr_cor, linetype = "dashed", color = "grey30", linewidth = 0.3) +
    scale_color_manual(values = scatter_cols[sample_order], guide = "none") +
    scale_x_continuous(labels = scientific) +
    coord_cartesian(xlim = c(0, x_lim)) +
    labs(
      title = "InferCNA Scatter by Sample",
      x = "CNA signal",
      y = "CNA correlation"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 10)
    )
  
  pdf(output_file, width = 13, height = 8.5, useDingbats = FALSE)
  print(p_all)
  print(p_facet)
  dev.off()
  
  scatter_summary
}

sample_mean <- sapply(sample_order, function(s) {
  cells <- rownames(target_meta)[target_meta$sample == s]
  rowMeans(target_outs[, cells, drop = FALSE])
})
saveRDS(sample_mean, file.path(out_dir, "Auto_PDO_cnv_sample_mean_profiles_Carroll_2023.rds"))

scatter_summary <- plot_cna_scatter(
  outs = outs,
  meta_df = combined_meta,
  ref_prefix = ref_prefix,
  sample_order = sample_order,
  output_file = file.path(out_dir, "Auto_PDO_cnv_scatter_Carroll_2023.pdf")
)
print(scatter_summary)

message("Done. Outputs written to: ", file.path(out_root, out_dir))
