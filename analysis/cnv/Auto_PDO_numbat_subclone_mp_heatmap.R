####################
# Auto_PDO_numbat_subclone_mp_heatmap.R
#
# Numbat-derived subclone versus PDO metaprogram/state heterogeneity.
# This mirrors the InferCNA subclone MP/state workflow but uses Numbat clone
# calls and Numbat posterior CNV segments projected onto genome-wide bins.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
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
min_subclone_frac <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 0.02
max_plot_cells <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1200L
max_display_clones <- if (length(args) >= 6 && nzchar(args[6])) as.integer(args[6]) else 12L

out_dir <- "Auto_PDO_numbat_subclone_mp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- "Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv"
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
required <- c(
  manifest_path,
  gene_order_path,
  "PDOs_merged.rds",
  "UCell_scores_filtered.rds",
  "Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds",
  "Auto_PDO_final_states.rds",
  "Auto_PDO_mp_adj_noreg.rds"
)
for (path in required) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

manifest <- fread(manifest_path)
if (!identical(sample_arg, "all")) {
  requested <- trimws(unlist(strsplit(sample_arg, ",")))
  manifest <- manifest[sample %in% requested]
}
manifest <- manifest[sample != "SUR843T3_PDO"]
if (nrow(manifest) == 0) stop("No samples found for argument: ", sample_arg)

chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- fread(gene_order_path, header = FALSE, fill = TRUE)
gene_order <- gene_order[, seq_len(4), with = FALSE]
setnames(gene_order, c("gene", "chromosome", "start", "end"))
gene_order <- gene_order[!is.na(gene) & nzchar(gene)]
gene_order <- gene_order %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(start = as.numeric(.data$start), end = as.numeric(.data$end)) %>%
  filter(is.finite(.data$start), is.finite(.data$end)) %>%
  mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
  arrange(.data$chromosome, .data$start)

make_genome_bins <- function(bin_size = 100L) {
  gene_order %>%
    mutate(.row = row_number()) %>%
    group_by(.data$chromosome) %>%
    mutate(bin_index = ((row_number() - 1L) %/% bin_size) + 1L,
           bin = paste0(.data$chromosome, "_", .data$bin_index)) %>%
    ungroup() %>%
    group_by(.data$bin, .data$chromosome, .data$bin_index) %>%
    summarise(start = min(.data$start, na.rm = TRUE),
              end = max(.data$end, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(.data$chromosome, .data$bin_index)
}
genome_bins <- make_genome_bins()

message("Loading PDO metadata and MP/state scores")
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
state_level_order <- c(names(state_groups), "3CA_EMT_and_Protein_maturation", "Unresolved", "Hybrid")

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black",
  "Unassigned" = "grey85"
)

mp_cols <- c(
  "MP6_G2M Cell Cycle" = "#B0B0B0",
  "MP7_DNA repair" = "#999999",
  "MP1_G2M checkpoint" = "#808080",
  "MP3_G1S Cell Cycle" = "#C0C0C0",
  "MP5_MYC-related Proliferation" = "#E41A1C",
  "MP4_Intestinal Metaplasia" = "#4DAF4A",
  "MP8_Columnar Progenitor" = "#FF7F00",
  "MP10_Inflammatory Stress Epi." = "#984EA3",
  "MP9_ECM Remodeling Epi." = "#C77CFF",
  "Unassigned" = "grey85"
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

infer_study <- function(sample_id) sub("^([^_]+).*$", "\\1", as.character(sample_id))

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

make_palette <- function(values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(character(0))
  base <- suppressWarnings(brewer.pal(max(3, min(12, length(values))), palette))
  setNames(colorRampPalette(base)(length(values)), values)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values[is.na(values) | !nzchar(values)] <- "Unassigned"
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(c(Unassigned = "grey85"))
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  missing_values <- setdiff(values, names(cols))
  if (length(missing_values) > 0) cols <- c(cols, make_palette(missing_values, palette))
  out <- cols[values]
  names(out) <- values
  out
}

subclone_colours <- function(values) complete_palette(make_palette(values, "Dark2"), values, "Dark2")

score_numbat_joint <- function(joint_post) {
  p_cols <- colnames(joint_post)
  if (all(c("p_amp", "p_bamp", "p_del", "p_bdel") %in% p_cols)) {
    return(as.numeric(joint_post$p_amp) + as.numeric(joint_post$p_bamp) -
             as.numeric(joint_post$p_del) - as.numeric(joint_post$p_bdel))
  }
  state <- as.character(joint_post$cnv_state)
  base <- case_when(
    state %in% c("amp", "bamp") ~ 1,
    state %in% c("del", "bdel") ~ -1,
    state %in% c("loh") ~ 0.5,
    TRUE ~ 0
  )
  prob <- if ("p_cnv" %in% p_cols) as.numeric(joint_post$p_cnv) else 1
  base * prob
}

read_numbat_clone_post <- function(row) {
  sample_id <- row$sample
  clone_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_clone_post.csv"))
  map_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_cell_map.csv"))
  if (!file.exists(clone_file) || !file.exists(map_file)) return(NULL)
  clone_post <- fread(clone_file)
  cell_map <- fread(map_file)
  clone_post <- clone_post %>%
    left_join(cell_map %>% select(.data$numbat_cell, .data$cell_id, .data$raw_barcode), by = c("cell" = "numbat_cell"))
  missing_cell <- is.na(clone_post$cell_id)
  if (any(missing_cell)) {
    raw_guess <- sub(paste0("^", sample_id, "_"), "", clone_post$cell[missing_cell])
    fallback <- cell_map$cell_id[match(raw_guess, cell_map$raw_barcode)]
    clone_post$cell_id[missing_cell] <- fallback
  }
  clone_post %>%
    transmute(
      cell = .data$cell_id,
      sample = sample_id,
      numbat_cell = .data$cell,
      subclone = paste0("Numbat clone ", .data$clone_opt),
      clone_opt = as.character(.data$clone_opt),
      compartment = if ("compartment_opt" %in% colnames(clone_post)) as.character(.data$compartment_opt) else NA_character_,
      p_cnv = if ("p_cnv" %in% colnames(clone_post)) as.numeric(.data$p_cnv) else NA_real_,
      p_cnv_expr = if ("p_cnv_x" %in% colnames(clone_post)) as.numeric(.data$p_cnv_x) else NA_real_,
      p_cnv_allele = if ("p_cnv_y" %in% colnames(clone_post)) as.numeric(.data$p_cnv_y) else NA_real_
    ) %>%
    filter(!is.na(.data$cell))
}

read_numbat_binned <- function(row, cells_use) {
  sample_id <- row$sample
  joint_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_joint_post.csv.gz"))
  map_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_cell_map.csv"))
  if (!file.exists(joint_file) || !file.exists(map_file)) return(NULL)
  joint_post <- fread(joint_file)
  cell_map <- fread(map_file)
  available_cells <- intersect(cells_use, cell_map$cell_id)
  if (length(available_cells) == 0) return(NULL)
  mat <- matrix(0, nrow = nrow(genome_bins), ncol = length(available_cells),
                dimnames = list(genome_bins$bin, available_cells))
  joint_post <- joint_post %>%
    left_join(cell_map %>% select(.data$numbat_cell, .data$cell_id), by = c("cell" = "numbat_cell")) %>%
    filter(.data$cell_id %in% available_cells)
  if (nrow(joint_post) == 0) return(list(mat = mat, chr = as.character(genome_bins$chromosome)))
  joint_post$score <- score_numbat_joint(joint_post)
  joint_post$CHROM <- paste0("chr", sub("^chr", "", as.character(joint_post$CHROM)))
  joint_post$seg_start <- as.numeric(joint_post$seg_start)
  joint_post$seg_end <- as.numeric(joint_post$seg_end)
  joint_post <- joint_post %>%
    filter(.data$CHROM %in% chrom_levels,
           is.finite(.data$seg_start),
           is.finite(.data$seg_end),
           is.finite(.data$score)) %>%
    group_by(.data$CHROM, .data$seg, .data$seg_start, .data$seg_end, .data$cell_id) %>%
    summarise(score = mean(.data$score, na.rm = TRUE), .groups = "drop")
  if (nrow(joint_post) == 0) return(list(mat = mat, chr = as.character(genome_bins$chromosome)))
  segment_keys <- joint_post %>% distinct(.data$CHROM, .data$seg, .data$seg_start, .data$seg_end)
  for (i in seq_len(nrow(segment_keys))) {
    seg_row <- segment_keys[i, ]
    bin_ix <- which(
      as.character(genome_bins$chromosome) == seg_row$CHROM &
        genome_bins$end >= seg_row$seg_start &
        genome_bins$start <= seg_row$seg_end
    )
    if (length(bin_ix) == 0) next
    seg_scores <- joint_post %>%
      filter(.data$CHROM == seg_row$CHROM,
             .data$seg == seg_row$seg,
             .data$seg_start == seg_row$seg_start,
             .data$seg_end == seg_row$seg_end)
    cell_ix <- match(seg_scores$cell_id, colnames(mat))
    cell_keep <- is.finite(cell_ix)
    if (!any(cell_keep)) next
    mat[bin_ix, cell_ix[cell_keep]] <- matrix(rep(seg_scores$score[cell_keep], each = length(bin_ix)), nrow = length(bin_ix))
  }
  list(mat = mat, chr = as.character(genome_bins$chromosome))
}

compute_sample_mp_scores <- function(cells) {
  cells <- intersect(cells, rownames(mp_score_mat))
  z <- as.matrix(mp_score_mat[cells, mp_names, drop = FALSE])
  top_use <- as.matrix(mp_adj_noncc[cells, topmp_mps, drop = FALSE])
  top_mp <- colnames(top_use)[max.col(top_use, ties.method = "first")]
  names(top_mp) <- rownames(top_use)
  list(z = z, top_mp = top_mp)
}

sample_plot_cells <- function(cells, subclone, max_cells = 1200L) {
  cells <- intersect(cells, names(subclone))
  if (length(cells) <= max_cells) return(cells)
  split_cells <- split(cells, factor(subclone[cells], levels = unique(subclone[cells])))
  target <- pmax(10L, floor(max_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  sampled <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

order_cells_by_subclone <- function(cna_mat, subclone) {
  if (is.null(names(subclone)) || any(is.na(names(subclone))) || any(!nzchar(names(subclone)))) {
    names(subclone) <- colnames(cna_mat)
  }
  subclone <- subclone[colnames(cna_mat)]
  split_cells <- split(colnames(cna_mat), factor(subclone, levels = unique(subclone)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    mat_sub <- cna_mat[, cells, drop = FALSE]
    if (all(mat_sub == 0, na.rm = TRUE)) return(cells)
    d <- dist(t(mat_sub))
    if (length(d) == 0 || any(!is.finite(d))) return(cells)
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}

make_numbat_heatmap <- function(binned, meta_plot, sample_id) {
  for (anno_col in c("subclone", "compartment", "top_mp_label", "state_label")) {
    x <- as.character(meta_plot[[anno_col]])
    x[is.na(x) | !nzchar(x)] <- "Unassigned"
    meta_plot[[anno_col]] <- factor(x, levels = unique(x))
  }
  ord <- rownames(meta_plot)
  mat <- binned$mat[, ord, drop = FALSE]
  row_chr <- factor(binned$chr, levels = unique(binned$chr))
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(levels(row_chr))), levels(row_chr))
  subclone_cols <- subclone_colours(meta_plot$subclone)
  topmp_cols <- complete_palette(mp_cols, meta_plot$top_mp_label, "Paired")
  local_state_cols <- complete_palette(state_cols, meta_plot$state_label, "Set3")
  comp_cols <- complete_palette(c(tumor = "#7B3294", normal = "#008837", unknown = "grey70"), meta_plot$compartment, "Set2")

  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  lim <- as.numeric(quantile(abs(vals[vals != 0]), 0.985, na.rm = TRUE))
  if (!is.finite(lim) || lim <= 0) lim <- 1
  lim <- min(max(lim, 0.25), 1.5)

  top_ha <- HeatmapAnnotation(
    Clone = meta_plot$subclone,
    Compartment = meta_plot$compartment,
    TopMP = meta_plot$top_mp_label,
    State = meta_plot$state_label,
    col = list(Clone = subclone_cols, Compartment = comp_cols, TopMP = topmp_cols, State = local_state_cols),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(2.5, "mm"),
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
    name = "Numbat CNV",
    col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
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
    column_title_gp = gpar(fontsize = 7, fontface = "bold"),
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
  Heatmap(
    mean_mat,
    name = "Mean MP z",
    col = colorRamp2(c(-0.75, 0, 0.75), c("#2166AC", "white", "#B2182B")),
    top_annotation = HeatmapAnnotation(Clone = factor(colnames(mean_mat), levels = colnames(mean_mat)),
                                       col = list(Clone = subclone_colours(colnames(mean_mat))),
                                       show_annotation_name = TRUE,
                                       annotation_name_side = "left",
                                       simple_anno_size = unit(3, "mm")),
    width = unit(max(35, 15 * ncol(mean_mat)), "mm"),
    height = unit(max(50, 5.5 * nrow(mean_mat)), "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7, fontface = "bold"),
    border = TRUE,
    rect_gp = gpar(col = "grey85", lwd = 0.3)
  )
}

make_corr_heatmap <- function(mp_z, subclone) {
  mean_mat <- sapply(unique(subclone), function(cl) colMeans(mp_z[names(subclone)[subclone == cl], , drop = FALSE], na.rm = TRUE))
  if (is.null(dim(mean_mat)) || ncol(mean_mat) < 2) return(textGrob("Only one Numbat clone", gp = gpar(fontsize = 12)))
  cm <- suppressWarnings(cor(mean_mat, method = "spearman", use = "pairwise.complete.obs"))
  cm[!is.finite(cm)] <- NA
  Heatmap(
    cm,
    name = "rho",
    col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
    cluster_rows = nrow(cm) > 2,
    cluster_columns = ncol(cm) > 2,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
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
  tests <- long %>%
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
  for (mp_id in unique(long$mp)) {
    mp_df <- long[long$mp == mp_id, , drop = FALSE]
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
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1L
    }
  }
  sub_tests <- bind_rows(sub_rows) %>%
    group_by(.data$sample) %>%
    mutate(p_adj = p.adjust(.data$p_value, method = "BH"),
           significant = !is.na(.data$p_adj) & .data$p_adj < 0.05 & abs(.data$delta_mean) >= 0.25) %>%
    ungroup()
  list(long = long, tests = tests, sub_tests = sub_tests)
}

test_states_by_subclone <- function(meta_all, sample_id) {
  tab <- table(meta_all$subclone, meta_all$state_label)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (sum(tab) == 0 || nrow(tab) < 2 || ncol(tab) < 2) {
    return(data.frame(sample = sample_id, test = NA_character_, p_value = NA_real_, cramers_v = NA_real_, stringsAsFactors = FALSE))
  }
  test_name <- "chisq"
  p_val <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
  expected <- tryCatch(suppressWarnings(chisq.test(tab)$expected), error = function(e) matrix(NA_real_, nrow = nrow(tab), ncol = ncol(tab)))
  if (any(expected < 5, na.rm = TRUE)) {
    test_name <- "chisq_simulated"
    p_val <- tryCatch(chisq.test(tab, simulate.p.value = TRUE, B = 10000)$p.value, error = function(e) NA_real_)
  }
  chi <- tryCatch(suppressWarnings(chisq.test(tab)$statistic), error = function(e) NA_real_)
  n <- sum(tab)
  denom <- n * (min(dim(tab)) - 1)
  v <- if (is.finite(chi) && denom > 0) as.numeric(sqrt(chi / denom)) else NA_real_
  data.frame(sample = sample_id, test = test_name, p_value = p_val, cramers_v = v, stringsAsFactors = FALSE)
}

make_boxplot <- function(score_df, mp_test_df, sample_id) {
  set.seed(42)
  point_df <- score_df %>%
    group_by(.data$mp_label, .data$subclone) %>%
    group_modify(~ .x[sample(seq_len(nrow(.x)), min(nrow(.x), 160L)), , drop = FALSE]) %>%
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
    geom_jitter(data = point_df, width = 0.12, size = 0.12, alpha = 0.12) +
    geom_text(data = y_pos, aes(x = 1, y = .data$y, label = .data$sig_label), inherit.aes = FALSE, size = 2.2) +
    facet_wrap(~mp_label, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = subclone_colours(score_df$subclone), drop = FALSE) +
    labs(title = paste0(sample_id, ": MP scores by Numbat clone"), x = NULL, y = "MP score z") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none", strip.text = element_text(size = 6.2),
          axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5))
}

make_state_distribution_plot <- function(meta_plot, state_test_df = NULL) {
  df <- meta_plot %>%
    count(.data$subclone, .data$state_label, name = "n") %>%
    group_by(.data$subclone) %>%
    mutate(pct = 100 * .data$n / sum(.data$n)) %>%
    ungroup()
  local_state_cols <- complete_palette(state_cols, df$state_label, "Set3")
  subtitle <- NULL
  if (!is.null(state_test_df) && nrow(state_test_df) > 0 && !is.na(state_test_df$p_value[1])) {
    subtitle <- paste0("State x clone p=", signif(state_test_df$p_value[1], 3),
                       ", Cramer's V=", signif(state_test_df$cramers_v[1], 3))
  }
  ggplot(df, aes(.data$subclone, .data$pct, fill = .data$state_label)) +
    geom_col(color = "black", linewidth = 0.15) +
    scale_fill_manual(values = local_state_cols, breaks = state_level_order, drop = FALSE) +
    labs(title = "State abundance", subtitle = subtitle, x = NULL, y = "% PDO cells", fill = "State") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6), legend.position = "right")
}

make_clone_qc_plot <- function(meta_plot) {
  qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "p_cnv", "p_cnv_expr", "p_cnv_allele"), colnames(meta_plot))
  if (length(qc_cols) == 0) return(textGrob("No QC columns available", gp = gpar(fontsize = 10)))
  plot_data <- meta_plot %>%
    select(subclone, all_of(qc_cols)) %>%
    pivot_longer(cols = all_of(qc_cols), names_to = "QC_Metric", values_to = "Value")
  ggplot(plot_data, aes(.data$subclone, .data$Value, fill = .data$subclone)) +
    geom_boxplot(outlier.size = 0.25, alpha = 0.85, linewidth = 0.25) +
    facet_wrap(~QC_Metric, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = subclone_colours(plot_data$subclone), drop = FALSE) +
    labs(title = "QC and Numbat posterior support", x = NULL, y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom", strip.text = element_text(size = 7))
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
clone_comp_all <- list()

sample_pdf <- file.path(out_dir, "Auto_PDO_numbat_subclone_mp_sample_pages.pdf")
pdf(sample_pdf, width = 18, height = 12, useDingbats = FALSE)

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  sample_id <- row$sample
  message("Processing ", sample_id)
  clone_post <- read_numbat_clone_post(row)
  if (is.null(clone_post) || nrow(clone_post) == 0) {
    blank_page(sample_id, "missing_numbat_clone_post")
    sample_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_numbat_clone_post", stringsAsFactors = FALSE)
    next
  }

  sample_cells <- Reduce(intersect, list(clone_post$cell, rownames(mp_score_mat), names(state_vec), rownames(meta_full)))
  if (length(sample_cells) < min_cells) {
    reason <- paste0("skipped_low_cells_with_numbat_state_and_mp_scores: ", length(sample_cells))
    blank_page(sample_id, reason)
    sample_rows[[sample_id]] <- data.frame(sample = sample_id, status = reason, n_cells = length(sample_cells), stringsAsFactors = FALSE)
    next
  }

  clone_post <- clone_post[match(sample_cells, clone_post$cell), ]
  clone_counts <- sort(table(clone_post$subclone), decreasing = TRUE)
  min_keep <- max(min_subclone_cells, ceiling(min_subclone_frac * length(sample_cells)))
  keep_clones <- names(clone_counts)[clone_counts >= min_keep]
  if (length(keep_clones) == 0) keep_clones <- names(clone_counts)[1]
  if (length(keep_clones) > max_display_clones) keep_clones <- names(clone_counts)[seq_len(max_display_clones)]

  clone_post <- clone_post %>%
    mutate(subclone = ifelse(.data$subclone %in% keep_clones, .data$subclone, "Other Numbat clones"),
           compartment = ifelse(is.na(.data$compartment) | !nzchar(.data$compartment), "unknown", .data$compartment))
  subclone_counts <- sort(table(clone_post$subclone), decreasing = TRUE)
  subclone_order <- c(names(subclone_counts)[names(subclone_counts) != "Other Numbat clones"],
                      intersect("Other Numbat clones", names(subclone_counts)))

  mp <- compute_sample_mp_scores(sample_cells)
  keep_cells <- Reduce(intersect, list(clone_post$cell, rownames(mp$z), names(state_vec), rownames(meta_full)))
  if (length(keep_cells) < min_cells) {
    reason <- paste0("skipped_low_cells_after_intersection: ", length(keep_cells))
    blank_page(sample_id, reason)
    sample_rows[[sample_id]] <- data.frame(sample = sample_id, status = reason, n_cells = length(sample_cells), stringsAsFactors = FALSE)
    next
  }

  clone_post <- clone_post[match(keep_cells, clone_post$cell), ]
  subclone <- setNames(clone_post$subclone, clone_post$cell)
  subclone <- factor(subclone, levels = subclone_order)
  names(subclone) <- clone_post$cell
  meta_epi <- meta_full[keep_cells, , drop = FALSE]
  subclone_label <- as.character(subclone[keep_cells])
  names(subclone_label) <- keep_cells
  state_label <- as.character(state_vec[keep_cells])
  names(state_label) <- keep_cells
  state_order <- c(state_level_order, sort(setdiff(unique(state_label), state_level_order)))
  state_order <- state_order[state_order %in% unique(state_label)]
  top_mp <- as.character(mp$top_mp[keep_cells])
  names(top_mp) <- keep_cells
  top_mp_label <- as.character(mp_labels[top_mp])
  names(top_mp_label) <- keep_cells
  topmp_order <- mp_labels[topmp_mps]
  topmp_order <- topmp_order[topmp_order %in% unique(top_mp_label)]

  meta_all <- data.frame(
    cell = keep_cells,
    sample = sample_id,
    subclone = factor(subclone_label, levels = subclone_order),
    clone_opt = clone_post$clone_opt,
    compartment = factor(clone_post$compartment, levels = c("tumor", "normal", "unknown")),
    top_mp = top_mp,
    top_mp_label = factor(top_mp_label, levels = topmp_order),
    state_label = factor(state_label, levels = state_order),
    p_cnv = clone_post$p_cnv,
    p_cnv_expr = clone_post$p_cnv_expr,
    p_cnv_allele = clone_post$p_cnv_allele,
    nCount_RNA = if ("nCount_RNA" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "nCount_RNA"]) else NA_real_,
    nFeature_RNA = if ("nFeature_RNA" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "nFeature_RNA"]) else NA_real_,
    percent.mt = if ("percent.mt" %in% colnames(meta_epi)) as.numeric(meta_epi[keep_cells, "percent.mt"]) else NA_real_,
    stringsAsFactors = FALSE,
    row.names = keep_cells
  )

  binned <- read_numbat_binned(row, keep_cells)
  if (is.null(binned)) {
    blank_page(sample_id, "missing_numbat_joint_posterior")
    sample_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_numbat_joint_posterior", n_cells = length(keep_cells), stringsAsFactors = FALSE)
    next
  }

  set.seed(42)
  plot_cells <- sample_plot_cells(keep_cells, subclone, max_plot_cells)
  cna_order <- order_cells_by_subclone(binned$mat[, plot_cells, drop = FALSE], factor(as.character(subclone[plot_cells]), levels = subclone_order))
  meta_plot <- meta_all[cna_order, , drop = FALSE]
  binned_plot <- list(mat = binned$mat[, cna_order, drop = FALSE], chr = binned$chr)
  test_res <- test_mps_by_subclone(mp$z[keep_cells, , drop = FALSE], subclone, sample_id)
  state_test <- test_states_by_subclone(meta_all, sample_id)

  cell_rows[[sample_id]] <- meta_all %>% as.data.frame() %>% select(-cell) %>% mutate(cell = rownames(meta_all), .before = 1)
  clone_comp_all[[sample_id]] <- meta_all %>%
    count(.data$sample, .data$subclone, .data$compartment, name = "n_cells") %>%
    group_by(.data$sample, .data$subclone) %>%
    mutate(frac = .data$n_cells / sum(.data$n_cells)) %>%
    ungroup()
  sample_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "analysed",
    n_cells = length(keep_cells),
    n_raw_numbat_clones = length(clone_counts),
    n_display_clones = length(unique(meta_all$subclone)),
    min_clone_cells_display = min_keep,
    median_p_cnv = median(meta_all$p_cnv, na.rm = TRUE),
    state_p_value = state_test$p_value[1],
    state_cramers_v = state_test$cramers_v[1],
    stringsAsFactors = FALSE
  )
  mp_tests_all[[sample_id]] <- test_res$tests
  sub_tests_all[[sample_id]] <- test_res$sub_tests
  state_tests_all[[sample_id]] <- state_test

  score_df <- test_res$long %>%
    mutate(subclone = factor(.data$subclone, levels = subclone_order),
           mp_label = factor(.data$mp_label, levels = mp_labels[mp_names]))

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 4,
    widths = unit(c(5.9, 2.9, 3.4, 3.6), "null"),
    heights = unit(c(5.8, 5.2), "null")
  )))

  cna_ht <- make_numbat_heatmap(binned_plot, meta_plot, sample_id)
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
  print(make_clone_qc_plot(meta_plot), vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
  popViewport()
}

dev.off()

cell_df <- bind_rows(cell_rows)
sample_df <- bind_rows(sample_rows)
mp_tests_df <- bind_rows(mp_tests_all)
sub_tests_df <- bind_rows(sub_tests_all)
state_tests_df <- bind_rows(state_tests_all)
clone_comp_df <- bind_rows(clone_comp_all)

if (nrow(state_tests_df) > 0 && "p_value" %in% colnames(state_tests_df)) {
  state_tests_df <- state_tests_df %>%
    mutate(p_adj = p.adjust(.data$p_value, method = "BH"),
           significant = !is.na(.data$p_adj) & .data$p_adj < 0.05)
}

write.csv(cell_df, file.path(out_dir, "Auto_PDO_numbat_subclone_cells.csv"), row.names = FALSE)
write.csv(sample_df, file.path(out_dir, "Auto_PDO_numbat_subclone_summary.csv"), row.names = FALSE)
write.csv(mp_tests_df, file.path(out_dir, "Auto_PDO_numbat_subclone_mp_tests.csv"), row.names = FALSE)
write.csv(sub_tests_df, file.path(out_dir, "Auto_PDO_numbat_subclone_mp_subclone_tests.csv"), row.names = FALSE)
write.csv(state_tests_df, file.path(out_dir, "Auto_PDO_numbat_subclone_state_tests.csv"), row.names = FALSE)
write.csv(clone_comp_df, file.path(out_dir, "Auto_PDO_numbat_subclone_compartment_summary.csv"), row.names = FALSE)

if (nrow(mp_tests_df) > 0) {
  multi_clone_samples <- sample_df$sample[sample_df$n_display_clones >= 2]
  sig_counts_sample <- mp_tests_df %>%
    filter(.data$sample %in% multi_clone_samples) %>%
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
      n_clone_tests = n(),
      n_significant_clone_tests = sum(.data$significant, na.rm = TRUE),
      pct_significant_clone_tests = 100 * mean(.data$significant, na.rm = TRUE),
      median_abs_delta = median(abs(.data$delta_mean), na.rm = TRUE),
      max_abs_delta = max(abs(.data$delta_mean), na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(sig_counts_sample, file.path(out_dir, "Auto_PDO_numbat_subclone_sig_count_summary.csv"), row.names = FALSE)
  write.csv(mp_summary_sample, file.path(out_dir, "Auto_PDO_numbat_subclone_mp_cohort_summary.csv"), row.names = FALSE)

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
    filter(.data$sample %in% multi_clone_samples, .data$mp %in% target_mps) %>%
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
    labs(title = "MP association with Numbat clones", x = NULL, y = "-log10(BH p)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

  p_state <- state_tests_df %>%
    filter(.data$sample %in% multi_clone_samples) %>%
    mutate(sample = factor(.data$sample, levels = .data$sample[order(.data$cramers_v, decreasing = TRUE)])) %>%
    ggplot(aes(.data$cramers_v, .data$sample, color = .data$p_adj < 0.05)) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey50"), na.value = "grey70") +
    labs(title = "State abundance association with Numbat clones", x = "Cramer's V", y = NULL, color = "BH p < 0.05") +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_text(size = 7), legend.position = "top")

  pdf(file.path(out_dir, "Auto_PDO_numbat_subclone_mp_cohort_summary.pdf"), width = 15, height = 9, useDingbats = FALSE)
  grid.arrange(p_counts, p_mp, p_state, ncol = 3, widths = c(0.9, 1.8, 1.1))
  dev.off()
}

message("Saved sample pages to: ", sample_pdf)
message("Saved tables and cohort summary to: ", out_dir)
