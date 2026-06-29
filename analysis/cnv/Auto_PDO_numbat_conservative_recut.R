####################
# Auto_PDO_numbat_conservative_recut.R
#
# Analysis registry
# Status: active terminal/audit workflow.
# Script: analysis/cnv/Auto_PDO_numbat_conservative_recut.R
# Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
# - PDOs_outs/Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv
# - PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/treeML_<iter>.rds
# - PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/geno_<iter>.tsv
# - PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/exp_post_<iter>.tsv
# - PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/allele_post_<iter>.tsv
# Outputs:
# - PDOs_outs/Auto_PDO_numbat/conservative_clones/Auto_PDO_numbat_tree_cut_sweep.csv
# - PDOs_outs/Auto_PDO_numbat/conservative_clones/Auto_PDO_numbat_conservative_clone_summary.csv
# - PDOs_outs/Auto_PDO_numbat/conservative_clones/Auto_PDO_numbat_conservative_phylogenetic_trees.pdf
# - PDOs_outs/Auto_PDO_numbat/conservative_clones/by_samples/<sample>/Auto_<sample>_numbat_conservative_clone_post.csv
# - PDOs_outs/Auto_PDO_numbat/conservative_clones/by_samples/<sample>/Auto_<sample>_tree_final_conservative.rds
# Downstream use: optional conservative Numbat clone layer for concordance
# plots and MP/state analyses; original Numbat outputs are left untouched.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(numbat)
  library(RColorBrewer)
  library(scales)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

preferred_n_cut <- as.integer(Sys.getenv("PDO_NUMBAT_CONSERVATIVE_N_CUT", "3"))
min_clone_frac <- as.numeric(Sys.getenv("PDO_NUMBAT_CONSERVATIVE_MIN_FRAC", "0.03"))
min_clone_cells_floor <- as.integer(Sys.getenv("PDO_NUMBAT_CONSERVATIVE_MIN_CELLS", "20"))
sweep_n_cuts <- seq_len(max(5L, preferred_n_cut))

manifest_path <- "Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv"
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)

manifest <- fread(manifest_path)
manifest <- manifest[sample != "SUR843T3_PDO"]
if (nrow(manifest) == 0) stop("No samples in Numbat manifest after exclusion.")

out_dir <- "Auto_PDO_numbat/conservative_clones"
by_sample_dir <- file.path(out_dir, "by_samples")
dir.create(by_sample_dir, recursive = TRUE, showWarnings = FALSE)
sweep_path <- file.path(out_dir, "Auto_PDO_numbat_tree_cut_sweep.csv")
reuse_sweep <- file.exists(sweep_path) && !identical(Sys.getenv("PDO_FORCE_REBUILD"), "1")

final_iter_from <- function(numbat_dir, prefix = "treeML") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.rds")))
  if (length(files) == 0) return(NA_integer_)
  iter <- suppressWarnings(as.integer(sub(paste0("^", prefix, "_([0-9]+)\\.rds$"), "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

as_vertex_df <- function(g) {
  attrs <- vertex_attr(g)
  as.data.frame(lapply(attrs, function(x) {
    if (is.list(x)) {
      vapply(x, function(y) paste(y, collapse = ","), character(1))
    } else {
      x
    }
  }), stringsAsFactors = FALSE)
}

read_genotype_matrix <- function(path) {
  geno <- fread(path)
  if (!("cell" %in% colnames(geno))) stop("Missing cell column in genotype matrix: ", path)
  cells <- geno$cell
  geno$cell <- NULL
  mat <- as.matrix(geno)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- cells
  mat
}

tree_clone_counts <- function(gtree) {
  vertex_df <- as_vertex_df(gtree)
  leaf <- as.logical(vertex_df$leaf)
  clone <- as.character(vertex_df$clone)
  clone[is.na(clone) | !nzchar(clone)] <- "unknown"
  sort(table(clone[leaf]), decreasing = TRUE)
}

summarise_counts <- function(counts) {
  if (length(counts) == 0) {
    return(data.frame(n_clones = 0L, min_clone_cells = NA_integer_, max_clone_cells = NA_integer_, stringsAsFactors = FALSE))
  }
  data.frame(
    n_clones = length(counts),
    min_clone_cells = as.integer(min(counts)),
    max_clone_cells = as.integer(max(counts)),
    stringsAsFactors = FALSE
  )
}

make_palette <- function(values) {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(character(0))
  base <- brewer.pal(max(3, min(12, length(values))), if (length(values) <= 8) "Dark2" else "Paired")
  setNames(colorRampPalette(base)(length(values)), values)
}

merge_minor_clone_post <- function(clone_post, gtree, min_clone_cells) {
  clone_post <- as.data.frame(clone_post)
  clone_post$clone_opt_raw_recut <- clone_post$clone_opt
  clone_post$GT_opt_raw_recut <- clone_post$GT_opt
  clone_post$p_opt_raw_recut <- clone_post$p_opt

  clone_raw <- as.character(clone_post$clone_opt_raw_recut)
  clone_counts <- sort(table(clone_raw), decreasing = TRUE)
  major <- names(clone_counts)[as.integer(clone_counts) >= min_clone_cells]
  if (length(major) == 0) major <- names(clone_counts)[1]

  p_cols <- paste0("p_", major)
  available_p_cols <- p_cols[p_cols %in% colnames(clone_post)]
  major_assign <- clone_raw
  minor_idx <- which(!clone_raw %in% major)

  if (length(minor_idx) > 0) {
    if (length(available_p_cols) > 0) {
      p_mat <- as.matrix(clone_post[minor_idx, available_p_cols, drop = FALSE])
      storage.mode(p_mat) <- "numeric"
      best_col <- max.col(p_mat, ties.method = "first")
      major_assign[minor_idx] <- sub("^p_", "", available_p_cols[best_col])
    } else {
      major_assign[minor_idx] <- major[1]
    }
  }

  vertex_df <- as_vertex_df(gtree)
  clone_gt <- vertex_df %>%
    filter(!is.na(.data$clone), nzchar(as.character(.data$clone))) %>%
    group_by(clone = as.character(.data$clone)) %>%
    summarise(GT = dplyr::first(as.character(.data$GT)), .groups = "drop")
  clone_gt_map <- setNames(clone_gt$GT, clone_gt$clone)

  clone_post$clone_opt <- major_assign
  clone_post$minor_clone_merged <- clone_raw != major_assign
  clone_post$GT_opt <- unname(clone_gt_map[as.character(clone_post$clone_opt)])
  clone_post$GT_opt[is.na(clone_post$GT_opt)] <- clone_post$GT_opt_raw_recut[is.na(clone_post$GT_opt)]

  merged_p_cols <- paste0("p_", clone_post$clone_opt)
  has_p <- merged_p_cols %in% colnames(clone_post)
  clone_post$p_opt <- clone_post$p_opt_raw_recut
  clone_post$p_opt[has_p] <- mapply(function(row_i, col_i) clone_post[row_i, col_i], which(has_p), merged_p_cols[has_p])
  clone_post
}

plot_tree_panel <- function(gtree, sample_id, iter, n_cut, clone_counts, cell_clone_map = NULL) {
  vertex_df <- as_vertex_df(gtree)
  root_idx <- which(as.logical(vertex_df$root))
  if (length(root_idx) == 0) root_idx <- 1L
  root_idx <- root_idx[1]
  leaf <- as.logical(vertex_df$leaf)
  clone <- as.character(vertex_df$clone)
  clone[is.na(clone) | !nzchar(clone)] <- "unknown"
  if (!is.null(cell_clone_map) && "name" %in% colnames(vertex_df)) {
    mapped <- unname(cell_clone_map[as.character(vertex_df$name)])
    clone[leaf & !is.na(mapped)] <- mapped[leaf & !is.na(mapped)]
  }
  clone_cols <- make_palette(clone)
  v_cols <- ifelse(leaf, clone_cols[clone], adjustcolor("grey55", alpha.f = 0.35))
  v_sizes <- ifelse(leaf, ifelse(sum(leaf) > 5000, 0.45, 0.75), 0.08)
  lay <- layout_as_tree(gtree, root = root_idx, circular = TRUE, mode = "out")
  plot(
    gtree,
    layout = lay,
    vertex.label = NA,
    vertex.size = v_sizes,
    vertex.color = v_cols,
    vertex.frame.color = NA,
    edge.arrow.mode = 0,
    edge.color = adjustcolor("grey45", alpha.f = 0.25),
    edge.width = 0.18,
    margin = -0.08,
    main = paste0(sample_id, "\nconservative merged Numbat tree, iter ", iter, ", n_cut=", n_cut)
  )
  legend(
    "bottomleft",
    legend = paste0("C", names(clone_counts), "  ", comma(as.integer(clone_counts))),
    col = clone_cols[names(clone_counts)],
    pch = 16,
    pt.cex = 1.1,
    bty = "n",
    cex = 0.72,
    title = "Clone cells"
  )
}

plot_clone_bar <- function(clone_counts) {
  frac <- as.numeric(clone_counts) / sum(clone_counts)
  cols <- make_palette(names(clone_counts))[names(clone_counts)]
  bp <- barplot(
    rev(frac),
    horiz = TRUE,
    col = rev(cols),
    border = NA,
    names.arg = paste0("C", rev(names(clone_counts))),
    xlim = c(0, max(frac, na.rm = TRUE) * 1.25),
    las = 1,
    cex.names = 0.75,
    xlab = "Fraction of cells",
    main = "Conservative clone sizes"
  )
  text(rev(frac), bp, labels = paste0(comma(rev(as.integer(clone_counts))), " (", percent(rev(frac), accuracy = 0.1), ")"), pos = 4, cex = 0.72, xpd = TRUE)
}

sweep_rows <- list()
summary_rows <- list()
pdf_path <- file.path(out_dir, "Auto_PDO_numbat_conservative_phylogenetic_trees.pdf")
pdf(pdf_path, width = 14, height = 8, useDingbats = FALSE)

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample[i]
  numbat_dir <- manifest$numbat_dir[i]
  iter <- final_iter_from(numbat_dir, "treeML")
  message("Conservative re-cut: ", sample_id)

  if (!is.finite(iter)) {
    plot.new()
    title(sample_id)
    text(0.5, 0.5, "No treeML_*.rds found")
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_treeML", stringsAsFactors = FALSE)
    next
  }

  tree_file <- file.path(numbat_dir, paste0("treeML_", iter, ".rds"))
  geno_file <- file.path(numbat_dir, paste0("geno_", iter, ".tsv"))
  exp_file <- file.path(numbat_dir, paste0("exp_post_", iter, ".tsv"))
  allele_file <- file.path(numbat_dir, paste0("allele_post_", iter, ".tsv"))
  original_tree_file <- file.path(numbat_dir, paste0("tree_final_", iter, ".rds"))

  missing <- c(tree_file, geno_file, exp_file, allele_file, original_tree_file)[!file.exists(c(tree_file, geno_file, exp_file, allele_file, original_tree_file))]
  if (length(missing) > 0) {
    plot.new()
    title(sample_id)
    text(0.5, 0.5, paste("Missing:", paste(basename(missing), collapse = ", ")))
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_input", missing = paste(missing, collapse = ";"), stringsAsFactors = FALSE)
    next
  }

  tree_ml <- readRDS(tree_file)
  P <- read_genotype_matrix(geno_file)
  n_cells <- nrow(P)
  min_clone_cells <- max(min_clone_cells_floor, ceiling(n_cells * min_clone_frac))

  original_counts <- tree_clone_counts(readRDS(original_tree_file))
  original_summary <- summarise_counts(original_counts)

  if (reuse_sweep) {
    sweep_cached <- fread(sweep_path)
    sample_sweep <- sweep_cached[sample == sample_id]
  } else {
    sample_sweep <- lapply(sweep_n_cuts, function(n_cut) {
      g_cut <- get_gtree(tree_ml, P, n_cut = n_cut, max_cost = 0)
      counts <- tree_clone_counts(g_cut)
      cbind(
        data.frame(sample = sample_id, iter = iter, cut_mode = "n_cut", n_cut = n_cut, n_cells = n_cells, stringsAsFactors = FALSE),
        summarise_counts(counts)
      )
    }) %>% bind_rows()
    sample_sweep$passes_min_clone <- sample_sweep$min_clone_cells >= min_clone_cells
    sample_sweep$min_clone_threshold <- min_clone_cells
    sweep_rows[[sample_id]] <- sample_sweep
  }

  candidates <- sample_sweep %>%
    filter(.data$n_cut <= preferred_n_cut, .data$n_clones <= preferred_n_cut + 1L, .data$passes_min_clone) %>%
    arrange(desc(.data$n_cut))
  selected_n_cut <- if (nrow(candidates) > 0) candidates$n_cut[1] else 1L

  g_conservative <- get_gtree(tree_ml, P, n_cut = selected_n_cut, max_cost = 0)
  recut_tree_counts <- tree_clone_counts(g_conservative)
  recut_tree_summary <- summarise_counts(recut_tree_counts)

  exp_post <- fread(exp_file)
  allele_post <- fread(allele_file)
  clone_post <- numbat:::get_clone_post(g_conservative, exp_post, allele_post)
  clone_post <- merge_minor_clone_post(clone_post, g_conservative, min_clone_cells)
  conservative_counts <- sort(table(as.character(clone_post$clone_opt)), decreasing = TRUE)
  conservative_summary <- summarise_counts(conservative_counts)
  minor_cells_merged <- sum(clone_post$minor_clone_merged, na.rm = TRUE)
  cell_clone_map <- setNames(as.character(clone_post$clone_opt), clone_post$cell)

  sample_out_dir <- file.path(by_sample_dir, sample_id)
  dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(g_conservative, file.path(sample_out_dir, paste0("Auto_", sample_id, "_tree_final_conservative.rds")))
  fwrite(clone_post, file.path(sample_out_dir, paste0("Auto_", sample_id, "_numbat_conservative_clone_post.csv")))
  clone_groups <- split(clone_post, clone_post$clone_opt)
  clones <- lapply(clone_groups, function(x) {
    list(sample = unique(x$clone_opt), members = unique(x$GT_opt), cells = x$cell, size = nrow(x))
  })
  saveRDS(clones, file.path(sample_out_dir, paste0("Auto_", sample_id, "_clones_conservative.rds")))

  layout(matrix(c(1, 2), nrow = 1), widths = c(2.1, 1))
  par(mar = c(0.4, 0.4, 2.4, 0.4))
  plot_tree_panel(g_conservative, sample_id, iter, selected_n_cut, conservative_counts, cell_clone_map)
  par(mar = c(4, 6.5, 2.4, 2))
  plot_clone_bar(conservative_counts)

  summary_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "recut",
    iter = iter,
    n_cells = n_cells,
    min_clone_threshold = min_clone_cells,
    original_n_clones = original_summary$n_clones,
    original_min_clone_cells = original_summary$min_clone_cells,
    original_max_clone_cells = original_summary$max_clone_cells,
    selected_n_cut = selected_n_cut,
    recut_tree_n_clones = recut_tree_summary$n_clones,
    recut_tree_min_clone_cells = recut_tree_summary$min_clone_cells,
    minor_cells_merged = minor_cells_merged,
    conservative_n_clones = conservative_summary$n_clones,
    conservative_min_clone_cells = conservative_summary$min_clone_cells,
    conservative_max_clone_cells = conservative_summary$max_clone_cells,
    conservative_clone_post = file.path(out_root, sample_out_dir, paste0("Auto_", sample_id, "_numbat_conservative_clone_post.csv")),
    conservative_tree = file.path(out_root, sample_out_dir, paste0("Auto_", sample_id, "_tree_final_conservative.rds")),
    stringsAsFactors = FALSE
  )
}

dev.off()

sweep_df <- if (reuse_sweep) fread(sweep_path) else bind_rows(sweep_rows)
summary_df <- bind_rows(summary_rows)
fwrite(sweep_df, file.path(out_dir, "Auto_PDO_numbat_tree_cut_sweep.csv"))
fwrite(summary_df, file.path(out_dir, "Auto_PDO_numbat_conservative_clone_summary.csv"))

message("Wrote: ", file.path(out_root, out_dir, "Auto_PDO_numbat_tree_cut_sweep.csv"))
message("Wrote: ", file.path(out_root, out_dir, "Auto_PDO_numbat_conservative_clone_summary.csv"))
message("Wrote: ", file.path(out_root, pdf_path))
