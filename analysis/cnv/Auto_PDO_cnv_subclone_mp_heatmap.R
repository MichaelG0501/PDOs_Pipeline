####################
# Auto_PDO_cnv_subclone_mp_heatmap.R
#
# CNA subclone versus PDO metaprogram/state heterogeneity.
#
# Run after analysis/cnv/Auto_PDO_infercna.R.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

args <- commandArgs(trailingOnly = TRUE)
sample_arg <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "all"
min_cells <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 40L
min_subclone_cells <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else 20L
min_subclone_frac <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 0.05
max_plot_cells <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1200L
max_subclones <- 6L
min_distinct_arm_delta <- 0.03
# min_subclone_silhouette <- 0.12
####################
louvain_k <- 15L
arm_call_threshold <- 0.10
strong_arm_delta <- 0.08
moderate_arm_delta <- 0.08
min_distinct_arms <- 2L
merge_max_arm_delta <- 0.06
same_pattern_cor <- 0.80
same_pattern_mean_delta <- 0.015
####################
mp_score_limit <- 2
cna_colour_limit <- 0.15
mp_mean_colour_limit <- 0.75

out_dir <- "Auto_PDO_cnv_subclone_mp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

infer_outs_path <- "cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds"
infer_meta_path <- "cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds"
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
excluded_samples <- c("SUR843T3_PDO")

for (path in c(infer_outs_path, infer_meta_path, gene_order_path, "PDOs_merged.rds",
               "UCell_scores_filtered.rds", "Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds",
               "Auto_PDO_final_states.rds", "Auto_PDO_mp_adj_noreg.rds")) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

gene_order <- read.table(
  gene_order_path,
  header = FALSE,
  col.names = c("gene", "chromosome", "start", "end"),
  stringsAsFactors = FALSE
)
chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- gene_order %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
  arrange(.data$chromosome, .data$start)

message("Loading PDO metadata and scores")
pdos <- readRDS("PDOs_merged.rds")
meta_full <- pdos@meta.data
if ("Batch" %in% colnames(meta_full)) {
  meta_full$Batch <- ifelse(meta_full$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", as.character(meta_full$Batch))
}

ucell_scores <- readRDS("UCell_scores_filtered.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")
mp_adj_noncc <- as.matrix(readRDS("Auto_PDO_mp_adj_noreg.rds"))
state_vec <- readRDS("Auto_PDO_final_states.rds")
state_names <- names(state_vec)
state_vec <- as.character(state_vec)
names(state_vec) <- state_names
state_vec[state_vec == "Basal to Intestinal Metaplasia"] <- "Basal to Intest. Meta"

mp_descriptions <- c(
  "MP6"  = "G2M Cell Cycle",
  "MP7"  = "DNA repair",
  "MP5"  = "MYC-related Proliferation",
  "MP1"  = "G2M checkpoint",
  "MP3"  = "G1S Cell Cycle",
  "MP8"  = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9"  = "ECM Remodeling Epi.",
  "MP4"  = "Intestinal Metaplasia"
)

state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "SMG-like Metaplasia" = c("MP8"),
  "Stress-adaptive" = c("MP10", "MP9")
)

cc_mps <- c("MP6", "MP7", "MP1", "MP3")
extra_state_order <- c("3CA_EMT_and_Protein_maturation")
state_level_order <- c(names(state_groups), extra_state_order, "Unresolved", "Hybrid")

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

####################
subclone_palette <- c(
  "Subclone A" = "#D73027",
  "Subclone B" = "#4575B4",
  "Subclone C" = "#1A9850",
  "Subclone D" = "#984EA3",
  "Subclone E" = "#FF7F00",
  "Subclone F" = "#A65628"
)
####################

mp_cols <- c(
  "MP6_G2M Cell Cycle" = "#B0B0B0",
  "MP7_DNA repair" = "#999999",
  "MP1_G2M checkpoint" = "#808080",
  "MP3_G1S Cell Cycle" = "#C0C0C0",
  "MP5_MYC-related Proliferation" = "#E41A1C",
  "MP4_Intestinal Metaplasia" = "#4DAF4A",
  "MP8_Columnar Progenitor" = "#FF7F00",
  "MP10_Inflammatory Stress Epi." = "#984EA3",
  "MP9_ECM Remodeling Epi." = "#C77CFF"
)

label_mp <- function(mps) {
  desc <- mp_descriptions[mps]
  desc[is.na(desc)] <- mps[is.na(desc)]
  paste0(mps, "_", desc)
}

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
if (!is.null(coverage_tbl)) {
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
  mp.genes <- mp.genes[!names(mp.genes) %in% low_coverage_mps]
}
retained_mps <- names(mp.genes)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- paste0("MP", rev(unique(ordered_clusters)))
mp_tree_order <- mp_tree_order[mp_tree_order %in% retained_mps]
ordered_state_mps <- unlist(lapply(state_groups, function(mps) {
  x <- mp_tree_order[mp_tree_order %in% mps]
  c(x, setdiff(mps, x))
}), use.names = FALSE)
mp_names <- unique(c(cc_mps[cc_mps %in% retained_mps], ordered_state_mps, mp_tree_order))
mp_names <- mp_names[mp_names %in% names(mp_descriptions)]
mp_labels <- setNames(label_mp(mp_names), mp_names)

infer_study <- function(sample_id) {
  sub("^([^_]+).*$", "\\1", as.character(sample_id))
}

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  study_sd <- clust_df %>%
    group_by(.data$.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop")
  study_names <- study_sd$.study
  study_sd <- as.matrix(study_sd[, colnames(mat), drop = FALSE])
  rownames(study_sd) <- study_names
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  clust_centered <- clust_df %>%
    group_by(.data$.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()
  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

score_sample_var <- as.character(meta_full$orig.ident)
names(score_sample_var) <- rownames(meta_full)
score_study_var <- if ("Batch" %in% colnames(meta_full)) as.character(meta_full$Batch) else infer_study(score_sample_var)
names(score_study_var) <- rownames(meta_full)

cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_common_cells <- intersect(rownames(ucell_scores), names(score_sample_var))
mp_adj_cc <- matrix(nrow = length(cc_common_cells), ncol = 0, dimnames = list(cc_common_cells, character(0)))
if (length(cc_in_ucell) > 0) {
  mp_adj_cc <- z_normalise(
    as.matrix(ucell_scores[cc_common_cells, cc_in_ucell, drop = FALSE]),
    score_sample_var,
    score_study_var
  )
}
score_common_cells <- intersect(rownames(mp_adj_noncc), cc_common_cells)
mp_score_mat <- cbind(
  mp_adj_noncc[score_common_cells, , drop = FALSE],
  mp_adj_cc[score_common_cells, , drop = FALSE]
)
mp_score_mat <- mp_score_mat[, intersect(mp_names, colnames(mp_score_mat)), drop = FALSE]
mp_names <- mp_names[mp_names %in% colnames(mp_score_mat)]
mp_labels <- setNames(label_mp(mp_names), mp_names)
topmp_mps <- setdiff(mp_names, cc_mps)
topmp_mps <- topmp_mps[topmp_mps %in% colnames(mp_adj_noncc)]
if (length(topmp_mps) == 0) stop("No non-cell-cycle PDO MP scores available for top-MP assignment.")

message("Loading InferCNA target matrix")
target_outs <- readRDS(infer_outs_path)
target_meta <- readRDS(infer_meta_path)
rownames(target_meta) <- target_meta$cell
target_meta <- target_meta[target_meta$sample != "SUR843T3_PDO", , drop = FALSE]

sample_dirs <- unique(as.character(target_meta$sample))
sample_dirs <- sample_dirs[!sample_dirs %in% excluded_samples]
if (!identical(sample_arg, "all")) {
  requested <- trimws(unlist(strsplit(sample_arg, ",")))
  sample_dirs <- intersect(sample_dirs, requested)
}
if (length(sample_dirs) == 0) stop("No samples found for argument: ", sample_arg)

make_palette <- function(values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values)]
  if (length(values) == 0) return(character(0))
  base <- suppressWarnings(brewer.pal(max(3, min(12, length(values))), palette))
  cols <- colorRampPalette(base)(length(values))
  stats::setNames(cols, values)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  cols <- as.character(cols)
  names(cols) <- names(cols)
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  missing_values <- setdiff(values, names(cols))
  if (length(missing_values) > 0) {
    cols <- c(cols, make_palette(missing_values, palette))
  }
  cols[values]
}

assert_named_palette <- function(cols, label) {
  if (length(cols) == 0 || is.null(names(cols)) || any(is.na(names(cols))) || any(!nzchar(names(cols)))) {
    stop("Invalid annotation palette for ", label, ": ",
         paste0("[", paste(names(cols), collapse = ","), "]"), call. = FALSE)
  }
  cols
}

####################
subclone_colours <- function(values) {
  assert_named_palette(complete_palette(subclone_palette, values, "Set2"), "Subclone")
}

finite_axis_lims <- function(values, probs = c(0.01, 0.99), pad = 0.05) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) == 0) return(c(0, 1))
  lims <- as.numeric(quantile(values, probs = probs, na.rm = TRUE))
  if (!all(is.finite(lims))) lims <- range(values, na.rm = TRUE)
  if (!is.finite(diff(lims)) || diff(lims) <= 0) {
    delta <- max(abs(lims[1]), 1) * pad
    lims <- lims + c(-delta, delta)
  } else {
    delta <- diff(lims) * pad
    lims <- lims + c(-delta, delta)
  }
  lims
}

percent_mt_cache <- new.env(parent = emptyenv())

find_percent_mt_column <- function(meta_df) {
  candidates <- c("percent.mt", "percent_mt", "percent.mito", "percent_mito", "pct_mito", "pct.mt")
  direct <- candidates[candidates %in% colnames(meta_df)]
  if (length(direct) > 0) return(direct[1])
  fuzzy <- grep("percent.*mt|mt.*percent|mito|mitochond", colnames(meta_df), ignore.case = TRUE, value = TRUE)
  if (length(fuzzy) > 0) fuzzy[1] else NA_character_
}

load_sample_percent_mt <- function(sample_id) {
  if (exists(sample_id, envir = percent_mt_cache, inherits = FALSE)) {
    return(get(sample_id, envir = percent_mt_cache, inherits = FALSE))
  }
  sample_rds <- file.path("by_samples", sample_id, paste0(sample_id, ".rds"))
  out <- numeric(0)
  if (file.exists(sample_rds)) {
    sample_obj <- readRDS(sample_rds)
    sample_meta <- sample_obj@meta.data
    mt_col <- find_percent_mt_column(sample_meta)
    if (!is.na(mt_col)) {
      vals <- as.numeric(sample_meta[[mt_col]])
      names(vals) <- rownames(sample_meta)
      pdo_names <- ifelse(grepl(paste0("^", sample_id, "_"), names(vals)), names(vals), paste(sample_id, names(vals), sep = "_"))
      vals <- c(vals, setNames(vals, pdo_names))
      out <- vals[!duplicated(names(vals))]
    }
  }
  assign(sample_id, out, envir = percent_mt_cache)
  out
}

lookup_percent_mt <- function(sample_id, target_metrics, keep_cells, meta_epi) {
  if ("percent.mt" %in% colnames(meta_epi)) {
    vals <- as.numeric(meta_epi[keep_cells, "percent.mt"])
    names(vals) <- keep_cells
    return(vals)
  }
  sample_vals <- load_sample_percent_mt(sample_id)
  vals <- sample_vals[keep_cells]
  missing <- is.na(vals)
  if (any(missing) && "original_cell" %in% colnames(target_metrics)) {
    original_cells <- as.character(target_metrics[keep_cells, "original_cell"])
    vals[missing] <- sample_vals[original_cells[missing]]
  }
  names(vals) <- keep_cells
  as.numeric(vals)
}
####################

compute_sample_mp_scores <- function(cells) {
  cells <- intersect(cells, rownames(mp_score_mat))
  z <- as.matrix(mp_score_mat[cells, mp_names, drop = FALSE])
  top_use <- as.matrix(mp_adj_noncc[cells, topmp_mps, drop = FALSE])
  top_mp <- colnames(top_use)[max.col(top_use, ties.method = "first")]
  names(top_mp) <- rownames(top_use)
  list(z = z, top_mp = top_mp)
}

centromere_pos <- c(
  chr1 = 121700000, chr2 = 91800000, chr3 = 87900000, chr4 = 50600000,
  chr5 = 48400000, chr6 = 61000000, chr7 = 59900000, chr8 = 45600000,
  chr9 = 49000000, chr10 = 40200000, chr11 = 53400000, chr12 = 35500000,
  chr13 = 17700000, chr14 = 17200000, chr15 = 19000000, chr16 = 36800000,
  chr17 = 25100000, chr18 = 18500000, chr19 = 26200000, chr20 = 28100000,
  chr21 = 12000000, chr22 = 15000000, chrX = 61000000
)

make_arm_labels <- function(go) {
  chr <- as.character(go$chromosome)
  arm <- ifelse(go$start <= centromere_pos[chr], "p", "q")
  paste0(chr, arm)
}

compute_arm_means <- function(cna_mat, cluster, arm_labels) {
  cl_levels <- sort(unique(as.character(cluster)))
  arm_levels <- unique(arm_labels)
  out <- do.call(rbind, lapply(cl_levels, function(cl) {
    cells <- names(cluster)[cluster == cl]
    vals <- tapply(seq_len(nrow(cna_mat)), factor(arm_labels, levels = arm_levels), function(ix) {
      mean(cna_mat[ix, cells, drop = FALSE], na.rm = TRUE)
    })
    vals[arm_levels]
  }))
  out <- as.matrix(out)
  rownames(out) <- cl_levels
  out[is.na(out)] <- 0
  out
}

####################
call_arm_cna <- function(arm_mean) {
  arm_call <- matrix(0L, nrow = nrow(arm_mean), ncol = ncol(arm_mean), dimnames = dimnames(arm_mean))
  arm_call[arm_mean > arm_call_threshold] <- 1L
  arm_call[arm_mean < -arm_call_threshold] <- -1L
  arm_call
}

summarise_subclone_separation <- function(cna_mat, cluster, arm_labels) {
  counts <- table(cluster)
  if (length(counts) < 2) {
    return(data.frame(
      min_max_arm_delta = NA_real_,
      max_max_arm_delta = NA_real_,
      n_pairs_same_arm_calls = NA_integer_,
      n_pairs_below_delta = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  arm_mean <- compute_arm_means(cna_mat, cluster, arm_labels)
  arm_call <- call_arm_cna(arm_mean)
  cl_levels <- rownames(arm_mean)
  pair_rows <- list()
  for (i in seq_len(length(cl_levels) - 1L)) {
    for (j in seq.int(i + 1L, length(cl_levels))) {
      delta <- abs(arm_mean[cl_levels[i], ] - arm_mean[cl_levels[j], ])
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        same_arm_calls = identical(unname(arm_call[cl_levels[i], ]), unname(arm_call[cl_levels[j], ])),
        max_arm_delta = max(delta, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  pair_df <- bind_rows(pair_rows)
  data.frame(
    min_max_arm_delta = min(pair_df$max_arm_delta, na.rm = TRUE),
    max_max_arm_delta = max(pair_df$max_arm_delta, na.rm = TRUE),
    n_pairs_same_arm_calls = sum(pair_df$same_arm_calls, na.rm = TRUE),
    n_pairs_below_delta = sum(pair_df$max_arm_delta < merge_max_arm_delta, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

louvain_cna_clusters <- function(cna_mat) {
  n_cells <- ncol(cna_mat)
  if (n_cells <= 2L) return(setNames(rep("1", n_cells), colnames(cna_mat)))
  k_use <- min(louvain_k, n_cells - 1L)
  cell_mat <- t(cna_mat)
  nn <- RANN::nn2(data = cell_mat, query = cell_mat, k = k_use + 1L)$nn.idx[, -1, drop = FALSE]
  edges <- cbind(rep(seq_len(n_cells), each = k_use), as.vector(t(nn)))
  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
  membership <- igraph::membership(igraph::cluster_louvain(graph))
  out <- as.character(membership)
  names(out) <- colnames(cna_mat)
  out
}

minimum_subclone_cells <- function(n_cells) {
  max(min_subclone_cells, ceiling(min_subclone_frac * n_cells))
}

drop_small_clusters <- function(cluster, n_cells) {
  min_keep <- minimum_subclone_cells(n_cells)
  counts <- table(cluster)
  keep_levels <- names(counts)[counts >= min_keep]
  if (length(keep_levels) == 0) keep_levels <- names(which.max(counts))
  cluster[cluster %in% keep_levels]
}

has_distinct_arm_evidence <- function(mean_a, mean_b, call_a, call_b) {
  delta <- abs(mean_a - mean_b)
  call_diff <- call_a != call_b
  any(call_diff & delta >= strong_arm_delta, na.rm = TRUE) ||
    sum(delta >= moderate_arm_delta, na.rm = TRUE) >= min_distinct_arms
}

has_same_scaled_pattern <- function(mean_a, mean_b) {
  centered_a <- mean_a - mean(mean_a, na.rm = TRUE)
  centered_b <- mean_b - mean(mean_b, na.rm = TRUE)
  profile_cor <- suppressWarnings(cor(centered_a, centered_b, use = "pairwise.complete.obs"))
  mean_delta <- mean(abs(mean_a - mean_b), na.rm = TRUE)
  is.finite(profile_cor) &&
    profile_cor >= same_pattern_cor &&
    is.finite(mean_delta) &&
    mean_delta <= same_pattern_mean_delta
}

merge_indistinct_clusters <- function(cna_mat, cluster, arm_labels) {
  cluster_names <- names(cluster)
  cluster <- as.character(cluster)
  names(cluster) <- cluster_names
  repeat {
    cl_levels <- sort(unique(cluster))
    if (length(cl_levels) <= 1) break
    arm_mean <- compute_arm_means(cna_mat, cluster, arm_labels)
    arm_call <- call_arm_cna(arm_mean)

    merge_pair <- NULL
    merge_score <- Inf
    for (i in seq_len(length(cl_levels) - 1L)) {
      for (j in seq.int(i + 1L, length(cl_levels))) {
        max_delta <- max(abs(arm_mean[cl_levels[i], ] - arm_mean[cl_levels[j], ]), na.rm = TRUE)
        distinct_evidence <- has_distinct_arm_evidence(
          arm_mean[cl_levels[i], ],
          arm_mean[cl_levels[j], ],
          arm_call[cl_levels[i], ],
          arm_call[cl_levels[j], ]
        )
        same_scaled_pattern <- has_same_scaled_pattern(
          arm_mean[cl_levels[i], ],
          arm_mean[cl_levels[j], ]
        )
        indistinct <- same_scaled_pattern ||
          !distinct_evidence ||
          (is.finite(max_delta) && max_delta < merge_max_arm_delta)
        if (indistinct && max_delta < merge_score) {
          merge_pair <- c(cl_levels[i], cl_levels[j])
          merge_score <- max_delta
        }
      }
    }
    if (is.null(merge_pair)) break
    counts <- table(cluster[cluster %in% merge_pair])
    keep <- names(which.max(counts))
    cluster[cluster %in% merge_pair] <- keep
  }
  cluster
}
####################

make_binned_cna <- function(cna_mat, go, bin_size = 100L) {
  go2 <- go %>%
    mutate(.row = seq_len(n())) %>%
    group_by(.data$chromosome) %>%
    mutate(bin = paste0(.data$chromosome, "_", ((row_number() - 1L) %/% bin_size) + 1L)) %>%
    ungroup()
  bin_levels <- unique(go2$bin)
  bins <- split(go2$.row, factor(go2$bin, levels = bin_levels))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(cna_mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  row_chr <- sub("_.*$", "", rownames(binned))
  list(mat = binned, chr = row_chr)
}

infer_cna_subclones <- function(cna_mat, go) {
  n_cells <- ncol(cna_mat)
  if (n_cells < min_cells) {
    out <- setNames(rep("Subclone A", n_cells), colnames(cna_mat))
    attr(out, "silhouette") <- NA_real_
    attr(out, "diagnostics") <- data.frame()
    attr(out, "selected_separation") <- data.frame()
    return(out)
  }
  arm_labels <- make_arm_labels(go)
  initial_cluster <- louvain_cna_clusters(cna_mat)
  initial_counts <- table(initial_cluster)
  initial_cluster <- drop_small_clusters(initial_cluster, n_cells)
  if (length(initial_cluster) < min_cells) {
    out <- setNames(rep("Subclone A", length(initial_cluster)), names(initial_cluster))
    attr(out, "silhouette") <- NA_real_
    attr(out, "diagnostics") <- data.frame(
      initial_k = length(initial_counts),
      final_k = 1L,
      initial_counts = paste(paste0(names(initial_counts), "=", as.integer(initial_counts)), collapse = ";"),
      final_counts = paste0("1=", length(initial_cluster)),
      dropped_small_cells = n_cells - length(initial_cluster),
      min_subclone_cells = minimum_subclone_cells(n_cells),
      stringsAsFactors = FALSE
    )
    attr(out, "selected_separation") <- data.frame()
    return(out)
  }
  cna_mat <- cna_mat[, names(initial_cluster), drop = FALSE]
  if (length(unique(initial_cluster)) <= 1L) {
    out <- setNames(rep("Subclone A", length(initial_cluster)), names(initial_cluster))
    attr(out, "silhouette") <- NA_real_
    attr(out, "diagnostics") <- data.frame(
      initial_k = length(initial_counts),
      final_k = 1L,
      initial_counts = paste(paste0(names(initial_counts), "=", as.integer(initial_counts)), collapse = ";"),
      final_counts = paste0("1=", length(initial_cluster)),
      dropped_small_cells = n_cells - length(initial_cluster),
      min_subclone_cells = minimum_subclone_cells(n_cells),
      stringsAsFactors = FALSE
    )
    attr(out, "selected_separation") <- data.frame()
    return(out)
  }

  merged_cluster <- merge_indistinct_clusters(cna_mat, initial_cluster, arm_labels)
  separation_df <- summarise_subclone_separation(cna_mat, merged_cluster, arm_labels)
  diagnostics <- cbind(
    data.frame(
      initial_k = length(initial_counts),
      final_k = length(unique(merged_cluster)),
      initial_counts = paste(paste0(names(initial_counts), "=", as.integer(initial_counts)), collapse = ";"),
      final_counts = paste(paste0(names(table(merged_cluster)), "=", as.integer(table(merged_cluster))), collapse = ";"),
      dropped_small_cells = n_cells - length(initial_cluster),
      min_subclone_cells = minimum_subclone_cells(n_cells),
      louvain_k = louvain_k,
      arm_call_threshold = arm_call_threshold,
      strong_arm_delta = strong_arm_delta,
      moderate_arm_delta = moderate_arm_delta,
      min_distinct_arms = min_distinct_arms,
      merge_max_arm_delta = merge_max_arm_delta,
      same_pattern_cor = same_pattern_cor,
      same_pattern_mean_delta = same_pattern_mean_delta,
      stringsAsFactors = FALSE
    ),
    separation_df
  )

  counts <- sort(table(merged_cluster), decreasing = TRUE)
  message("  inferred CNA clusters: ",
          paste(paste0(names(counts), "=", as.integer(counts)), collapse = "; "),
          "; louvain_k=", louvain_k)
  clone_letters <- LETTERS[seq_along(counts)]
  names(clone_letters) <- names(counts)
  out <- paste0("Subclone ", clone_letters[merged_cluster])
  names(out) <- names(merged_cluster)
  out <- out[colnames(cna_mat)]
  attr(out, "silhouette") <- NA_real_
  attr(out, "diagnostics") <- diagnostics
  attr(out, "selected_separation") <- separation_df
  out
}

prepare_cna_matrix <- function(outs, cells) {
  cells <- intersect(cells, colnames(outs))
  common_genes <- intersect(rownames(outs), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(outs[go$gene, cells, drop = FALSE])
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  signal <- rowMeans(abs(mat), na.rm = TRUE)
  keep_signal <- signal >= as.numeric(quantile(signal, probs = 1 / 3, na.rm = TRUE))
  list(mat = mat[keep_signal, , drop = FALSE], gene_order = go[keep_signal, , drop = FALSE])
}

order_cells_by_subclone <- function(cna_mat, subclone) {
  split_cells <- split(names(subclone), factor(subclone, levels = unique(subclone)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    d <- dist(t(cna_mat[, cells, drop = FALSE]))
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}

sample_plot_cells <- function(cells, subclone, max_cells = 1200L) {
  cells <- intersect(cells, names(subclone))
  if (length(cells) <= max_cells) return(cells)
  split_cells <- split(cells, factor(subclone[cells], levels = unique(subclone[cells])))
  target <- pmax(20L, floor(max_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  sampled <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

make_cna_heatmap <- function(binned, meta_plot, sample_id) {
  ord <- rownames(meta_plot)
  mat <- binned$mat[, ord, drop = FALSE]
  row_chr <- factor(binned$chr, levels = unique(binned$chr))
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(levels(row_chr))), levels(row_chr))
  subclone_cols <- subclone_colours(meta_plot$subclone)
  topmp_cols <- mp_cols[names(mp_cols) %in% unique(as.character(meta_plot$top_mp_label))]
  local_state_cols <- state_cols[names(state_cols) %in% unique(as.character(meta_plot$state_label))]
  missing_states <- setdiff(unique(meta_plot$state_label), names(local_state_cols))
  if (length(missing_states) > 0) {
    local_state_cols <- c(local_state_cols, setNames(hue_pal()(length(missing_states)), missing_states))
  }
  missing_topmp <- setdiff(unique(as.character(meta_plot$top_mp_label)), names(topmp_cols))
  if (length(missing_topmp) > 0) {
    topmp_cols <- c(topmp_cols, setNames(hue_pal()(length(missing_topmp)), missing_topmp))
  }
  subclone_cols <- assert_named_palette(subclone_cols, "Subclone")
  topmp_cols <- assert_named_palette(topmp_cols, "TopMP")
  local_state_cols <- assert_named_palette(local_state_cols, "State")

  top_ha <- HeatmapAnnotation(
    Subclone = meta_plot$subclone,
    TopMP = meta_plot$top_mp_label,
    State = meta_plot$state_label,
    col = list(Subclone = subclone_cols, TopMP = topmp_cols, State = local_state_cols),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(3, "mm"),
    na_col = "grey90"
  )

  left_ha <- rowAnnotation(
    Chr = row_chr,
    col = list(Chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(3, "mm")
  )

  Heatmap(
    mat,
    name = "CNA",
    col = colorRamp2(c(-cna_colour_limit, 0, cna_colour_limit), c("#2166AC", "white", "#B2182B")),
    top_annotation = top_ha,
    left_annotation = left_ha,
    row_split = row_chr,
    row_gap = unit(0, "mm"),
    column_split = factor(meta_plot$subclone, levels = unique(meta_plot$subclone)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title_rot = 30,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_title_gp = gpar(fontsize = 7),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA),
    column_title = sample_id
  )
}

make_mean_mp_heatmap <- function(mp_z, subclone) {
  mean_mat <- sapply(unique(subclone), function(cl) colMeans(mp_z[names(subclone)[subclone == cl], , drop = FALSE], na.rm = TRUE))
  if (is.null(dim(mean_mat))) mean_mat <- matrix(mean_mat, ncol = 1, dimnames = list(colnames(mp_z), unique(subclone)))
  mean_mat <- mean_mat[intersect(mp_names, rownames(mean_mat)), , drop = FALSE]
  rownames(mean_mat) <- mp_labels[rownames(mean_mat)]
  subclone_cols <- subclone_colours(colnames(mean_mat))
  top_ha <- HeatmapAnnotation(
    Subclone = factor(colnames(mean_mat), levels = colnames(mean_mat)),
    col = list(Subclone = subclone_cols),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    simple_anno_size = unit(3, "mm")
  )
  Heatmap(
    mean_mat,
    name = "Mean MP z",
    col = colorRamp2(c(-mp_mean_colour_limit, 0, mp_mean_colour_limit), c("#2166AC", "white", "#B2182B")),
    top_annotation = top_ha,
    width = unit(max(30, 18 * ncol(mean_mat)), "mm"),
    height = unit(max(50, 5.5 * nrow(mean_mat)), "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.3)
  )
}

make_corr_heatmap <- function(mp_z, subclone) {
  mean_mat <- sapply(unique(subclone), function(cl) colMeans(mp_z[names(subclone)[subclone == cl], , drop = FALSE], na.rm = TRUE))
  if (is.null(dim(mean_mat)) || ncol(mean_mat) < 2) {
    return(textGrob("Only one CNA subclone", gp = gpar(fontsize = 12)))
  }
  cm <- suppressWarnings(cor(mean_mat, method = "spearman", use = "pairwise.complete.obs"))
  cm[!is.finite(cm)] <- NA
  Heatmap(
    cm,
    name = "rho",
    col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
    cluster_rows = nrow(cm) > 2,
    cluster_columns = ncol(cm) > 2,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.3)
  )
}

test_mps_by_subclone <- function(mp_z, subclone, sample_id) {
  df <- as.data.frame(mp_z)
  df$cell <- rownames(df)
  df$subclone <- subclone[df$cell]
  long <- df %>%
    pivot_longer(all_of(colnames(mp_z)), names_to = "mp", values_to = "score_z") %>%
    mutate(mp_label = mp_labels[.data$mp])
  test_long <- long %>% filter(.data$subclone != "Unresolved")

  tests <- test_long %>%
    group_by(.data$mp, .data$mp_label) %>%
    summarise(
      p_value = if (n_distinct(.data$subclone) >= 2) {
        ms <- tapply(.data$score_z, .data$subclone, mean, na.rm = TRUE)
        hi_cl <- names(ms)[which.max(ms)]
        lo_cl <- names(ms)[which.min(ms)]
        tryCatch(wilcox.test(.data$score_z[.data$subclone == hi_cl], .data$score_z[.data$subclone == lo_cl])$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(sample = sample_id, p_adj = p.adjust(.data$p_value, method = "BH"))

  sub_rows <- list()
  row_i <- 1L
  for (mp_id in unique(test_long$mp)) {
    mp_df <- test_long[test_long$mp == mp_id, , drop = FALSE]
    mp_label <- unique(mp_df$mp_label)[1]
    for (cl in unique(mp_df$subclone)) {
      x <- mp_df$score_z[mp_df$subclone == cl]
      y <- mp_df$score_z[mp_df$subclone != cl]
      mean_score <- mean(x, na.rm = TRUE)
      rest_mean <- mean(y, na.rm = TRUE)
      delta_mean <- mean_score - rest_mean
      p_val <- if (length(x) >= 2 && length(y) >= 2) {
        tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
      sub_rows[[row_i]] <- data.frame(
        sample = sample_id,
        mp = mp_id,
        mp_label = mp_label,
        subclone = cl,
        n_cells = length(x),
        mean_score = mean_score,
        rest_mean = rest_mean,
        delta_mean = delta_mean,
        p_value = p_val,
        frac_large_effect = ifelse(delta_mean >= 0,
                                   mean(x > 1, na.rm = TRUE),
                                   mean(x < -1, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  sub_tests <- bind_rows(sub_rows)
  sub_tests <- sub_tests %>%
    group_by(.data$sample) %>%
    mutate(p_adj = p.adjust(.data$p_value, method = "BH"),
           significant = !is.na(.data$p_adj) & .data$p_adj < 0.05 & abs(.data$delta_mean) >= 0.25) %>%
    ungroup()

  list(long = long, tests = tests, sub_tests = sub_tests)
}

test_states_by_subclone <- function(meta_all, sample_id) {
  tab <- table(meta_all$subclone, meta_all$state_label)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(data.frame(sample = sample_id, test = NA_character_, p_value = NA_real_, cramers_v = NA_real_, stringsAsFactors = FALSE))
  }
  test_name <- "chisq"
  p_val <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5, na.rm = TRUE)) {
    test_name <- "chisq_simulated"
    p_val <- tryCatch(chisq.test(tab, simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA_real_)
  }
  chi <- suppressWarnings(chisq.test(tab)$statistic)
  n <- sum(tab)
  denom <- n * (min(dim(tab)) - 1)
  v <- if (is.finite(chi) && denom > 0) as.numeric(sqrt(chi / denom)) else NA_real_
  data.frame(sample = sample_id, test = test_name, p_value = p_val, cramers_v = v, stringsAsFactors = FALSE)
}

make_boxplot <- function(score_df, mp_test_df, sample_id) {
  set.seed(42)
  point_df <- score_df %>%
    group_by(.data$mp_label, .data$subclone) %>%
    group_modify(~ .x[sample(seq_len(nrow(.x)), min(nrow(.x), 200L)), , drop = FALSE]) %>%
    ungroup()

  label_df <- mp_test_df %>%
    mutate(sig_label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "NS"
    ))
  label_df$mp_label <- factor(label_df$mp_label, levels = unique(score_df$mp_label))

  y_pos <- score_df %>%
    group_by(.data$mp_label) %>%
    summarise(y = max(.data$score_z, na.rm = TRUE) + 0.25, .groups = "drop") %>%
    left_join(label_df, by = "mp_label")

  ggplot(score_df, aes(.data$subclone, .data$score_z, fill = .data$subclone)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
    geom_jitter(data = point_df, width = 0.12, size = 0.15, alpha = 0.15) +
    geom_text(data = y_pos, aes(x = 1, y = .data$y, label = .data$sig_label), inherit.aes = FALSE, size = 2.4) +
    facet_wrap(~mp_label, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = subclone_colours(score_df$subclone), drop = FALSE) +
    labs(title = paste0(sample_id, ": MP scores by CNA subclone"), x = NULL, y = "MP score z") +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6)
    )
}

make_qc_boxplot <- function(meta_plot, sample_id) {
  top_subs <- meta_plot %>% count(subclone) %>% arrange(desc(n)) %>% head(2) %>% pull(subclone)
  if (length(top_subs) < 2) return(ggplot() + theme_void() + ggtitle("Only one subclone"))

  qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "cna_signal", "cna_cor"), colnames(meta_plot))
  plot_data <- meta_plot %>%
    filter(.data$subclone %in% top_subs) %>%
    select(subclone, all_of(qc_cols)) %>%
    pivot_longer(cols = all_of(qc_cols), names_to = "QC_Metric", values_to = "Value") %>%
    mutate(QC_Metric = factor(.data$QC_Metric, levels = qc_cols))

  stats_df <- plot_data %>%
    group_by(.data$QC_Metric) %>%
    summarise(
      p_value = tryCatch({
        sub1_vals <- Value[subclone == top_subs[1]]
        sub2_vals <- Value[subclone == top_subs[2]]
        if (length(na.omit(sub1_vals)) > 2 && length(na.omit(sub2_vals)) > 2) {
          wilcox.test(sub1_vals, sub2_vals)$p.value
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_),
      max_val = ifelse(any(is.finite(Value)), max(Value, na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(label = ifelse(!is.na(p_value) & p_value < 0.05, sprintf("p=%.2e", p_value), "NS"))

  ggplot(plot_data, aes(x = subclone, y = Value, fill = subclone)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8, linewidth = 0.3) +
    facet_wrap(~QC_Metric, scales = "free_y", nrow = 2) +
    geom_text(data = stats_df, aes(x = 1.5, y = max_val * 1.05, label = label), inherit.aes = FALSE, size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = subclone_colours(plot_data$subclone), drop = FALSE) +
    theme_classic(base_size = 8) +
    labs(title = "QC and CNA Metrics (Top 2 Subclones)", x = NULL, y = "Value") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", strip.text = element_text(size = 7))
}

make_state_distribution_plot <- function(meta_plot, state_test_df = NULL) {
  df <- meta_plot %>%
    count(.data$subclone, .data$state_label, name = "n") %>%
    group_by(.data$subclone) %>%
    mutate(pct = 100 * .data$n / sum(.data$n)) %>%
    ungroup()
  local_state_cols <- state_cols[names(state_cols) %in% unique(as.character(df$state_label))]
  missing_states <- setdiff(unique(as.character(df$state_label)), names(local_state_cols))
  if (length(missing_states) > 0) local_state_cols <- c(local_state_cols, setNames(hue_pal()(length(missing_states)), missing_states))
  subtitle <- NULL
  if (!is.null(state_test_df) && nrow(state_test_df) > 0 && !is.na(state_test_df$p_value[1])) {
    subtitle <- paste0("State x subclone p=", signif(state_test_df$p_value[1], 3),
                       ", Cramer's V=", signif(state_test_df$cramers_v[1], 3))
  }
  ggplot(df, aes(.data$subclone, .data$pct, fill = .data$state_label)) +
    geom_col(color = "black", linewidth = 0.15) +
    scale_fill_manual(values = local_state_cols, breaks = state_level_order, drop = FALSE) +
    labs(title = "State distribution", subtitle = subtitle, x = NULL, y = "% PDO cells", fill = "State") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
}

blank_page <- function(sample_id, reason) {
  grid.newpage()
  grid.text(sample_id, x = 0.03, y = 0.96, just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))
  grid.text(reason, x = 0.03, y = 0.88, just = c("left", "top"), gp = gpar(fontsize = 11))
}

cell_rows <- list()
sample_rows <- list()
mp_tests_all <- list()
sub_tests_all <- list()
state_tests_all <- list()
qc_tests_all <- list()
cluster_diagnostics_all <- list()

sample_pdf <- file.path(out_dir, "Auto_PDO_cnv_subclone_mp_sample_pages.pdf")
pdf(sample_pdf, width = 18, height = 12, useDingbats = FALSE)

for (sample_id in sample_dirs) {
  message("Processing ", sample_id)
  sample_meta <- target_meta[target_meta$sample == sample_id, , drop = FALSE]
  sample_meta <- sample_meta[sample_meta$cell %in% colnames(target_outs), , drop = FALSE]
  sample_meta <- sample_meta[!duplicated(sample_meta$original_cell), , drop = FALSE]
  sample_meta$pdo_cell <- ifelse(
    sample_meta$original_cell %in% rownames(meta_full),
    sample_meta$original_cell,
    paste(sample_meta$sample, sample_meta$original_cell, sep = "_")
  )

  outs <- target_outs[, sample_meta$cell, drop = FALSE]
  colnames(outs) <- sample_meta$pdo_cell
  target_metrics <- sample_meta
  rownames(target_metrics) <- target_metrics$pdo_cell

  sample_cells <- Reduce(intersect, list(colnames(outs), rownames(mp_score_mat), names(state_vec), rownames(meta_full)))
  if (length(sample_cells) < min_cells) {
    reason <- paste0("skipped_low_cells_with_cna_state_and_mp_scores: ", length(sample_cells))
    blank_page(sample_id, reason)
    sample_rows[[sample_id]] <- data.frame(sample = sample_id, status = reason, n_pdo_cells = length(sample_cells), n_subclones = 0, stringsAsFactors = FALSE)
    next
  }

  cna_prepped <- prepare_cna_matrix(outs, sample_cells)
  if (nrow(cna_prepped$mat) < 100 || ncol(cna_prepped$mat) < min_cells) {
    reason <- paste0("skipped_low_cna_signal_matrix: ", nrow(cna_prepped$mat), " genes x ", ncol(cna_prepped$mat), " cells")
    blank_page(sample_id, reason)
    sample_rows[[sample_id]] <- data.frame(sample = sample_id, status = reason, n_pdo_cells = length(sample_cells), n_subclones = 0, stringsAsFactors = FALSE)
    next
  }

  subclone <- infer_cna_subclones(cna_prepped$mat, cna_prepped$gene_order)
  subclone_silhouette <- attr(subclone, "silhouette")
  if (is.null(subclone_silhouette)) subclone_silhouette <- NA_real_
  subclone_diagnostics <- attr(subclone, "diagnostics")
  if (!is.null(subclone_diagnostics) && nrow(subclone_diagnostics) > 0) {
    cluster_diagnostics_all[[sample_id]] <- subclone_diagnostics %>% mutate(sample = sample_id, .before = 1)
  }
  subclone_separation <- attr(subclone, "selected_separation")
  if (is.null(subclone_separation) || nrow(subclone_separation) == 0) {
    subclone_separation <- data.frame(
      min_max_arm_delta = NA_real_,
      max_max_arm_delta = NA_real_,
      n_pairs_same_arm_calls = NA_integer_,
      n_pairs_below_delta = NA_integer_
    )
  }

  mp <- compute_sample_mp_scores(names(subclone))
  keep_cells <- Reduce(intersect, list(names(subclone), rownames(mp$z), names(state_vec), rownames(meta_full)))
  if (length(keep_cells) < min_cells) {
    reason <- paste0("skipped_low_cells_after_subclone_intersection: ", length(keep_cells))
    blank_page(sample_id, reason)
    sample_rows[[sample_id]] <- data.frame(
      sample = sample_id,
      status = reason,
      n_pdo_cells = length(sample_cells),
      n_subclones = length(unique(subclone)),
      subclone_silhouette = subclone_silhouette,
      stringsAsFactors = FALSE
    )
    next
  }

  subclone <- subclone[keep_cells]
  meta_epi <- meta_full[keep_cells, , drop = FALSE]
  state_label <- as.character(state_vec[keep_cells])
  names(state_label) <- keep_cells
  top_mp_label <- mp_labels[mp$top_mp[keep_cells]]
  names(top_mp_label) <- keep_cells
  percent_mt_vals <- lookup_percent_mt(sample_id, target_metrics, keep_cells, meta_epi)
  subclone_order <- paste0("Subclone ", LETTERS[seq_along(unique(subclone))])
  subclone_order <- subclone_order[subclone_order %in% unique(subclone)]
  state_order <- c(state_level_order, sort(setdiff(unique(state_label), state_level_order)))
  state_order <- state_order[state_order %in% unique(state_label)]
  topmp_order <- mp_labels[topmp_mps]
  topmp_order <- topmp_order[topmp_order %in% unique(top_mp_label)]

  meta_all <- data.frame(
    cell = keep_cells,
    sample = sample_id,
    subclone = factor(subclone[keep_cells], levels = subclone_order),
    top_mp = mp$top_mp[keep_cells],
    top_mp_label = factor(top_mp_label[keep_cells], levels = topmp_order),
    state_label = factor(state_label[keep_cells], levels = state_order),
    cna_signal = if ("cna_signal" %in% colnames(target_metrics)) as.numeric(target_metrics[keep_cells, "cna_signal"]) else NA_real_,
    cna_cor = if ("cna_cor" %in% colnames(target_metrics)) as.numeric(target_metrics[keep_cells, "cna_cor"]) else NA_real_,
    nCount_RNA = if ("nCount_RNA" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "nCount_RNA"]) else NA_real_,
    nFeature_RNA = if ("nFeature_RNA" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "nFeature_RNA"]) else NA_real_,
    percent.mt = percent_mt_vals,
    stringsAsFactors = FALSE,
    row.names = keep_cells
  )

  set.seed(42)
  plot_cells <- sample_plot_cells(keep_cells, subclone, max_plot_cells)
  cna_order <- order_cells_by_subclone(cna_prepped$mat[, plot_cells, drop = FALSE], factor(subclone[plot_cells], levels = subclone_order))
  meta_plot <- meta_all[cna_order, , drop = FALSE]
  binned <- make_binned_cna(cna_prepped$mat[, cna_order, drop = FALSE], cna_prepped$gene_order)
  test_res <- test_mps_by_subclone(mp$z[keep_cells, , drop = FALSE], subclone, sample_id)
  state_test <- test_states_by_subclone(meta_all, sample_id)

  cell_rows[[sample_id]] <- meta_all %>% as.data.frame() %>% select(-cell) %>% mutate(cell = rownames(meta_all), .before = 1)
  sample_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "analysed",
    n_pdo_cells = length(sample_cells),
    n_cells_used = length(keep_cells),
    n_subclones = length(unique(subclone)),
    subclone_silhouette = subclone_silhouette,
    subclone_min_max_arm_delta = subclone_separation$min_max_arm_delta[1],
    subclone_max_max_arm_delta = subclone_separation$max_max_arm_delta[1],
    subclone_n_pairs_same_arm_calls = subclone_separation$n_pairs_same_arm_calls[1],
    subclone_n_pairs_below_delta = subclone_separation$n_pairs_below_delta[1],
    state_source = "Auto_PDO_final_states.rds",
    stringsAsFactors = FALSE
  )
  mp_tests_all[[sample_id]] <- test_res$tests
  sub_tests_all[[sample_id]] <- test_res$sub_tests
  state_tests_all[[sample_id]] <- state_test

  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "cna_signal", "cna_cor")
  qc_test_rows <- list()
  for (q in qc_metrics) {
    if (q %in% colnames(meta_all) && n_distinct(meta_all$subclone) >= 2) {
      ms <- tapply(meta_all[[q]], meta_all$subclone, mean, na.rm = TRUE)
      ms <- ms[is.finite(ms)]
      if (length(ms) < 2) next
      hi_cl <- names(ms)[which.max(ms)]
      lo_cl <- names(ms)[which.min(ms)]
      p <- tryCatch(wilcox.test(meta_all[[q]][meta_all$subclone == hi_cl],
                                meta_all[[q]][meta_all$subclone == lo_cl])$p.value,
                    error = function(e) NA_real_)
      qc_test_rows[[q]] <- data.frame(sample = sample_id, metric = q, p_value = p, stringsAsFactors = FALSE)
    }
  }
  qc_tests_all[[sample_id]] <- bind_rows(qc_test_rows)

  score_df <- test_res$long %>%
    mutate(subclone = factor(.data$subclone, levels = subclone_order),
           mp_label = factor(.data$mp_label, levels = mp_labels[mp_names]))

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 4,
    widths = unit(c(5.9, 2.7, 3.6, 3.6), "null"),
    heights = unit(c(5.8, 5.2), "null")
  )))

  cna_ht <- make_cna_heatmap(binned, meta_plot, sample_id)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  draw(cna_ht, newpage = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right")
  popViewport()

  mean_ht <- make_mean_mp_heatmap(mp$z[keep_cells, , drop = FALSE], subclone)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  draw(mean_ht, newpage = FALSE, heatmap_legend_side = "right")
  popViewport()

  corr_obj <- make_corr_heatmap(mp$z[keep_cells, , drop = FALSE], subclone)
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
  if (inherits(corr_obj, "Heatmap")) {
    draw(corr_obj, newpage = FALSE, heatmap_legend_side = "right")
  } else {
    grid.draw(corr_obj)
  }
  popViewport()

  print(make_boxplot(score_df, test_res$tests, sample_id), vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))
  print(make_state_distribution_plot(meta_all, state_test), vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(make_qc_boxplot(meta_plot, sample_id), vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
  popViewport()
}

dev.off()

cell_df <- bind_rows(cell_rows)
sample_df <- bind_rows(sample_rows)
mp_tests_df <- bind_rows(mp_tests_all)
sub_tests_df <- bind_rows(sub_tests_all)
state_tests_df <- bind_rows(state_tests_all)
qc_tests_df <- bind_rows(qc_tests_all)
cluster_diagnostics_df <- bind_rows(cluster_diagnostics_all)
if (nrow(qc_tests_df) > 0) {
  qc_tests_df <- qc_tests_df %>% mutate(p_adj = p.adjust(.data$p_value, method = "BH"))
}
if (nrow(state_tests_df) > 0 && "p_value" %in% colnames(state_tests_df)) {
  state_tests_df <- state_tests_df %>%
    mutate(p_adj = p.adjust(.data$p_value, method = "BH"),
           significant = !is.na(.data$p_adj) & .data$p_adj < 0.05)
}

write.csv(cell_df, file.path(out_dir, "Auto_PDO_cnv_subclone_cells.csv"), row.names = FALSE)
write.csv(sample_df, file.path(out_dir, "Auto_PDO_cnv_subclone_summary.csv"), row.names = FALSE)
write.csv(mp_tests_df, file.path(out_dir, "Auto_PDO_cnv_subclone_mp_tests.csv"), row.names = FALSE)
write.csv(sub_tests_df, file.path(out_dir, "Auto_PDO_cnv_subclone_mp_subclone_tests.csv"), row.names = FALSE)
write.csv(state_tests_df, file.path(out_dir, "Auto_PDO_cnv_subclone_state_tests.csv"), row.names = FALSE)
write.csv(qc_tests_df, file.path(out_dir, "Auto_PDO_cnv_subclone_qc_tests.csv"), row.names = FALSE)
write.csv(cluster_diagnostics_df, file.path(out_dir, "Auto_PDO_cnv_subclone_cluster_diagnostics.csv"), row.names = FALSE)

if (nrow(sub_tests_df) > 0 && nrow(mp_tests_df) > 0) {
  multi_subclone_samples <- sample_df$sample[sample_df$n_subclones >= 2]
  if (length(multi_subclone_samples) > 0) {
    sig_counts_sample <- mp_tests_df %>%
      filter(.data$sample %in% multi_subclone_samples) %>%
      mutate(significant = !is.na(.data$p_adj) & .data$p_adj < 0.05) %>%
      group_by(.data$sample) %>%
      summarise(n_sig_mps = sum(.data$significant, na.rm = TRUE), .groups = "drop") %>%
      mutate(category = case_when(
        n_sig_mps == 0 ~ "None",
        n_sig_mps == 1 ~ "One significant",
        TRUE ~ "More than one"
      ))

    mp_summary_sample <- sub_tests_df %>%
      group_by(.data$mp, .data$mp_label) %>%
      summarise(
        n_subclone_tests = n(),
        n_significant_subclone_tests = sum(.data$significant, na.rm = TRUE),
        pct_significant_subclone_tests = 100 * mean(.data$significant, na.rm = TRUE),
        median_abs_delta = median(abs(.data$delta_mean), na.rm = TRUE),
        max_abs_delta = max(abs(.data$delta_mean), na.rm = TRUE),
        .groups = "drop"
      )

    write.csv(sig_counts_sample, file.path(out_dir, "Auto_PDO_cnv_subclone_sig_count_summary.csv"), row.names = FALSE)
    write.csv(mp_summary_sample, file.path(out_dir, "Auto_PDO_cnv_subclone_mp_cohort_summary.csv"), row.names = FALSE)

    p_counts <- sig_counts_sample %>%
      count(.data$category, name = "n") %>%
      mutate(category = factor(.data$category, levels = c("None", "One significant", "More than one")),
             pct = 100 * .data$n / sum(.data$n)) %>%
      ggplot(aes(.data$category, .data$pct, fill = .data$category)) +
      geom_col(color = "black", linewidth = 0.3) +
      geom_text(aes(label = paste0(round(.data$pct, 1), "%")), vjust = -0.3, size = 4) +
      scale_fill_manual(values = c("None" = "grey70", "One significant" = "#FDB863", "More than one" = "#B2182B")) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(title = "Significant MP differences per sample", x = NULL, y = "Percentage of samples") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none")

    target_mps <- mp_names[mp_names %in% names(mp_labels)]
    mp_plot_df <- mp_tests_df %>%
      filter(.data$sample %in% multi_subclone_samples, .data$mp %in% target_mps) %>%
      mutate(mp_label = factor(.data$mp_label, levels = mp_labels[target_mps]),
             val = -log10(pmax(.data$p_adj, .Machine$double.xmin)),
             val_plot = pmax(.data$val, 1e-3))

    mp_pcts <- mp_plot_df %>%
      group_by(.data$mp_label) %>%
      summarise(pct = 100 * mean(.data$p_adj < 0.05, na.rm = TRUE), .groups = "drop")

    p_mp <- ggplot(mp_plot_df, aes(.data$mp_label, .data$val_plot, fill = .data$mp_label)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_text(data = mp_pcts, aes(x = .data$mp_label, y = max(mp_plot_df$val_plot, na.rm = TRUE) * 1.15, label = sprintf("%.1f%%", .data$pct)), inherit.aes = FALSE, size = 3) +
      scale_y_log10(expand = expansion(mult = c(0.1, 0.3))) +
      scale_fill_manual(values = mp_cols[mp_labels[target_mps]]) +
      labs(title = NULL, x = NULL, y = "-log10(p_adj)") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

    plot_list_qc <- list()
    qc_metrics <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "cna_signal", "cna_cor"), colnames(cell_df))
    qc_summary_rows <- list()
    for (samp in multi_subclone_samples) {
      samp_df <- cell_df %>% filter(.data$sample == samp)
      top_subs <- samp_df %>% count(.data$subclone) %>% arrange(desc(.data$n)) %>% head(2) %>% pull(.data$subclone)
      if (length(top_subs) < 2) next
      row_data <- data.frame(sample = samp, clone1 = top_subs[1], clone2 = top_subs[2], stringsAsFactors = FALSE)
      for (q in qc_metrics) {
        cl_means <- samp_df %>%
          filter(.data$subclone %in% top_subs) %>%
          group_by(.data$subclone) %>%
          summarise(m = mean(.data[[q]], na.rm = TRUE), .groups = "drop") %>%
          filter(is.finite(.data$m))
        if (nrow(cl_means) < 2) next
        row_data[[paste0("X_", q)]] <- max(cl_means$m, na.rm = TRUE)
        row_data[[paste0("Y_", q)]] <- min(cl_means$m, na.rm = TRUE)
      }
      qc_summary_rows[[length(qc_summary_rows) + 1L]] <- row_data
    }
    qc_summary_df <- bind_rows(qc_summary_rows)
    for (q in qc_metrics) {
      x_col <- paste0("X_", q)
      y_col <- paste0("Y_", q)
      if (!all(c(x_col, y_col) %in% colnames(qc_summary_df))) next
      wt <- tryCatch(wilcox.test(qc_summary_df[[x_col]], qc_summary_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
      diff_stat <- mean(qc_summary_df[[x_col]] - qc_summary_df[[y_col]], na.rm = TRUE)
      p_val_display <- if (is.na(wt$p.value)) "NA" else if (wt$p.value < 0.001) sprintf("%.2e", wt$p.value) else sprintf("%.3f", wt$p.value)
      qc_lims <- finite_axis_lims(unlist(qc_summary_df[, c(x_col, y_col), drop = FALSE]))
      plot_list_qc[[q]] <- ggplot(qc_summary_df, aes(.data[[x_col]], .data[[y_col]])) +
        geom_point(alpha = 0.5, size = 1.2, color = "black") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        coord_cartesian(xlim = qc_lims, ylim = qc_lims) +
        labs(title = q, subtitle = sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat), x = "Highest subclone", y = "Lowest subclone") +
        theme_classic(base_size = 9) +
        theme(plot.title = element_text(size = 8, face = "bold"),
              plot.subtitle = element_text(size = 7.5))
    }
    if (nrow(qc_tests_df) > 0) {
      qc_plot_df <- qc_tests_df %>%
        filter(.data$sample %in% multi_subclone_samples) %>%
        mutate(val = -log10(pmax(.data$p_adj, .Machine$double.xmin)),
               val_plot = pmax(.data$val, 1e-3))
      qc_pcts <- qc_plot_df %>% group_by(.data$metric) %>% summarise(pct = 100 * mean(.data$p_adj < 0.05, na.rm = TRUE), .groups = "drop")
      plot_list_qc[["QC Significance"]] <- ggplot(qc_plot_df, aes(.data$metric, .data$val_plot, fill = .data$metric)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_text(data = qc_pcts, aes(x = .data$metric, y = max(qc_plot_df$val_plot, na.rm = TRUE) * 1.15, label = sprintf("%.1f%%", .data$pct)), inherit.aes = FALSE, size = 2.5) +
        scale_y_log10(expand = expansion(mult = c(0.1, 0.3))) +
        labs(title = "QC/CNA significance", x = NULL, y = "-log10(p_adj)") +
        theme_classic(base_size = 9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.title = element_text(size = 8, face = "bold"))
    }
    if (length(plot_list_qc) == 0) {
      plot_list_qc[["No QC"]] <- ggplot() + theme_void() + ggtitle("No QC tests available")
    }

    valid_cells <- intersect(cell_df$cell, rownames(mp_score_mat))
    cell_df_valid <- cell_df[match(valid_cells, cell_df$cell), ]
    mp_scores_valid <- mp_score_mat[valid_cells, target_mps, drop = FALSE]
    subclone_means <- cell_df_valid %>%
      select(cell, sample, subclone) %>%
      bind_cols(as.data.frame(mp_scores_valid)) %>%
      filter(.data$sample %in% multi_subclone_samples) %>%
      group_by(.data$sample, .data$subclone) %>%
      summarise(across(all_of(target_mps), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

    pair_rows <- list()
    for (samp in unique(subclone_means$sample)) {
      samp_df <- subclone_means %>% filter(.data$sample == samp) %>% arrange(.data$subclone)
      if (nrow(samp_df) >= 2) {
        combos <- combn(nrow(samp_df), 2)
        for (i in seq_len(ncol(combos))) {
          idx1 <- combos[1, i]
          idx2 <- combos[2, i]
          row_data <- data.frame(sample = samp, clone1 = samp_df$subclone[idx1], clone2 = samp_df$subclone[idx2])
          for (mp in target_mps) {
            val1 <- samp_df[[mp]][idx1]
            val2 <- samp_df[[mp]][idx2]
            row_data[[paste0("X_", mp)]] <- max(val1, val2, na.rm = TRUE)
            row_data[[paste0("Y_", mp)]] <- min(val1, val2, na.rm = TRUE)
          }
          pair_rows[[length(pair_rows) + 1]] <- row_data
        }
      }
    }
    pairs_df <- bind_rows(pair_rows)

    plot_list_mp <- list()
    if (nrow(pairs_df) > 0) {
      mp_x_cols <- paste0("X_", target_mps)
      mp_y_cols <- paste0("Y_", target_mps)
      mp_all_vals <- unlist(pairs_df[, c(mp_x_cols, mp_y_cols)])
      mp_global_lims <- quantile(mp_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
      for (mp in target_mps) {
        x_col <- paste0("X_", mp)
        y_col <- paste0("Y_", mp)
        wt <- tryCatch(wilcox.test(pairs_df[[x_col]], pairs_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
        diff_stat <- mean(pairs_df[[x_col]] - pairs_df[[y_col]], na.rm = TRUE)
        p_val_display <- if (is.na(wt$p.value)) "NA" else if (wt$p.value < 0.001) sprintf("%.2e", wt$p.value) else sprintf("%.3f", wt$p.value)
        plot_list_mp[[mp]] <- ggplot(pairs_df, aes(.data[[x_col]], .data[[y_col]])) +
          geom_point(alpha = 0.5, color = mp_cols[mp_labels[mp]], size = 1.2) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
          labs(title = mp_labels[mp], subtitle = sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat), x = "Higher clone", y = "Lower clone") +
          coord_cartesian(xlim = mp_global_lims, ylim = mp_global_lims) +
          theme_classic(base_size = 9) +
          theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 7.5))
      }
    }

    target_states <- state_level_order
    state_counts <- cell_df %>%
      filter(.data$sample %in% multi_subclone_samples) %>%
      count(.data$sample, .data$subclone, .data$state_label) %>%
      group_by(.data$sample, .data$subclone) %>%
      mutate(prop = .data$n / sum(.data$n)) %>%
      ungroup()

    pair_rows_state <- list()
    for (samp in unique(state_counts$sample)) {
      samp_df_st <- state_counts %>% filter(.data$sample == samp)
      subs <- sort(unique(samp_df_st$subclone))
      if (length(subs) >= 2) {
        combos <- combn(length(subs), 2)
        for (i in seq_len(ncol(combos))) {
          sub1 <- subs[combos[1, i]]
          sub2 <- subs[combos[2, i]]
          row_data <- data.frame(sample = samp, clone1 = sub1, clone2 = sub2)
          for (st in target_states) {
            p1 <- samp_df_st$prop[samp_df_st$subclone == sub1 & samp_df_st$state_label == st]
            p2 <- samp_df_st$prop[samp_df_st$subclone == sub2 & samp_df_st$state_label == st]
            val1 <- if (length(p1) > 0) p1[1] else 0
            val2 <- if (length(p2) > 0) p2[1] else 0
            row_data[[paste0("X_", st)]] <- max(val1, val2, na.rm = TRUE)
            row_data[[paste0("Y_", st)]] <- min(val1, val2, na.rm = TRUE)
          }
          pair_rows_state[[length(pair_rows_state) + 1]] <- row_data
        }
      }
    }
    pairs_state_df <- bind_rows(pair_rows_state)

    plot_list_state <- list()
    if (nrow(pairs_state_df) > 0) {
      state_x_cols <- paste0("X_", target_states)
      state_y_cols <- paste0("Y_", target_states)
      state_all_vals <- unlist(pairs_state_df[, c(state_x_cols, state_y_cols)])
      state_global_lims <- quantile(state_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
      for (st in target_states) {
        x_col <- paste0("X_", st)
        y_col <- paste0("Y_", st)
        wt <- tryCatch(wilcox.test(pairs_state_df[[x_col]], pairs_state_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
        diff_stat <- mean(pairs_state_df[[x_col]] - pairs_state_df[[y_col]], na.rm = TRUE)
        p_val_display <- if (is.na(wt$p.value)) "NA" else if (wt$p.value < 0.001) sprintf("%.2e", wt$p.value) else sprintf("%.3f", wt$p.value)
        st_col <- if (st %in% names(state_cols)) state_cols[[st]] else "grey50"
        plot_list_state[[st]] <- ggplot(pairs_state_df, aes(.data[[x_col]], .data[[y_col]])) +
          geom_point(alpha = 0.5, color = st_col, size = 1.2) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
          labs(title = st, subtitle = sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat), x = "Higher clone", y = "Lower clone") +
          coord_cartesian(xlim = state_global_lims, ylim = state_global_lims) +
          scale_x_continuous(labels = scales::percent) +
          scale_y_continuous(labels = scales::percent) +
          theme_classic(base_size = 9) +
          theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 7.5))
      }
    }

    pdf(file.path(out_dir, "Auto_PDO_cnv_subclone_mp_cohort_summary.pdf"), width = 15, height = 9, useDingbats = FALSE)
    grid.arrange(p_counts, p_mp, ncol = 2, widths = c(1, 2))
    grid.arrange(grobs = plot_list_qc, ncol = 3)
    if (length(plot_list_mp) > 0) grid.arrange(grobs = plot_list_mp, ncol = 4)
    if (length(plot_list_state) > 0) grid.arrange(grobs = plot_list_state, ncol = 4)
    dev.off()
  }
}

message("Saved sample pages to: ", sample_pdf)
message("Saved tables and cohort summary to: ", out_dir)
