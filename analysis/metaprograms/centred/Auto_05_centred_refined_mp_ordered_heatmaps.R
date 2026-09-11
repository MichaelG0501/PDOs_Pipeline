####################
# Analysis registry:
#   Status: active terminal QC/ordering workflow
#   Script: analysis/metaprograms/centred/Auto_05_centred_refined_mp_ordered_heatmaps.R
#   Methodology: analysis/methodology/metaprograms/centred/Auto_centred_ordering_and_state_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#
#   Description:
#     Replots the centred GeneNMF program similarity matrix ordered by the final
#     centred refined MPs, then plots centred refined MP UCell score correlations
#     using the requested state grouping and MP order. Adapted from scRef
#     centred pipeline 05_centred_refined_mp_ordered_heatmaps.R. The canonical
#     step-04 filtering audit is authoritative; this script does not impose an
#     additional gene-count or coverage threshold.
#
#   Inputs:
#     - live: PDOs_outs/centred_mp_refinement/optimal_nMP.rds
#     - live: PDOs_outs/centred_mp_refinement/geneNMF_metaprograms_nMP_{optimal}.rds
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_mp_assignments.rds
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds
#     - live: PDOs_outs/centred_mp_refinement/refined_ucell_scores.rds
#     - live: PDOs_outs/PDOs_merged.rds
#
#   Outputs (live: PDOs_outs/centred_mp_refinement/):
#     figures/centred_refined_mp_nmf_ordered_heatmap.pdf
#     figures/centred_refined_mp_nmf_ordered_heatmap.png
#     figures/centred_refined_mp_correlation_heatmap.pdf
#     figures/centred_refined_mp_correlation_heatmap.png
#     figures/centred_refined_mp_split_correlation_ordered_heatmap.pdf
#     figures/centred_refined_mp_split_correlation_ordered_heatmap.png
#     tables/centred_refined_mp_state_grouping.csv
#     tables/centred_refined_mp_nmf_ordered_blocks.csv
#     tables/centred_refined_mp_correlation_mean_rho.csv
#     tables/centred_refined_mp_correlation_p_values.csv
#     tables/centred_refined_mp_split_correlation_ordered_final_blocks.csv
#     tables/centred_refined_mp_split_correlation_ordered_sub_blocks.csv
#   Outputs (ephemeral):
#     intermediate/centred_refined_mp_ordered_heatmaps.rds
#
#   Downstream use:
#     - The grouping table and strict-order RDS are current inputs to centred
#       state definition, matched-FLOT, and cross-dataset projection scripts.
#     - Heatmaps are terminal QC figures.
#
#   Cache/replot behavior:
#     Rebuilds plot matrices from existing step-04 objects; it never reruns NMF,
#     refinement, or UCell scoring.
#
#   Run command: qsub Auto_replot_05.sh
#   Conda env: dmtcp
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

# === Paths ===
live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
ephemeral_base <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

outdir_live <- file.path(live_base, "centred_mp_refinement")
outdir_ephemeral <- file.path(ephemeral_base, "centred_mp_refinement", "intermediate")

for (subdir in c("tables", "figures")) {
  dir.create(file.path(outdir_live, subdir), recursive = TRUE, showWarnings = FALSE)
}
dir.create(outdir_ephemeral, recursive = TRUE, showWarnings = FALSE)
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")

optimal_file <- file.path(outdir_live, "optimal_nMP.rds")
assignment_file <- file.path(outdir_live, "merged_refined_mp_assignments.rds")
if (!file.exists(assignment_file)) {
  assignment_file <- file.path(outdir_ephemeral, "merged_refined_mp_assignments.rds")
}
gene_file <- file.path(outdir_live, "merged_refined_mp_genes.rds")
ucell_file <- file.path(outdir_live, "merged_refined_ucell_scores.rds")
if (!file.exists(ucell_file)) {
  ucell_file <- file.path(outdir_ephemeral, "merged_refined_ucell_scores.rds")
}
refined_ucell_file <- file.path(outdir_live, "refined_ucell_scores.rds")
if (!file.exists(refined_ucell_file)) {
  refined_ucell_file <- file.path(outdir_ephemeral, "refined_ucell_scores.rds")
}
seurat_file <- file.path(live_base, "PDOs_merged.rds")

required_inputs <- c(optimal_file, assignment_file, gene_file, ucell_file, refined_ucell_file, seurat_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing required input(s). Run centred GeneNMF and centred MP refinement first: ",
    paste(missing_inputs, collapse = ", ")
  )
}

optimal_nMP <- readRDS(optimal_file)
nmf_file <- file.path(outdir_live, paste0("geneNMF_metaprograms_nMP_", optimal_nMP, ".rds"))
if (!file.exists(nmf_file)) {
  stop("Missing centred GeneNMF metaprogram object: ", nmf_file)
}

cat("Loading centred refined MP inputs...\n")
geneNMF.metaprograms <- readRDS(nmf_file)
merged_assignments <- readRDS(assignment_file)
merged_mp_genes <- readRDS(gene_file)
merged_ucell <- readRDS(ucell_file)
pdos_merged <- readRDS(seurat_file)

required_assignment_cols <- c("program", "original_mp", "refined_submp", "merged_refined_mp")
missing_cols <- setdiff(required_assignment_cols, colnames(merged_assignments))
if (length(missing_cols) > 0) {
  stop("Assignment table is missing column(s): ", paste(missing_cols, collapse = ", "))
}

####################
# PDO centred refined MP state grouping and strict plotting order.
# Dotted lines before: MP19+, MP15, MP13b, MP9, MP18
####################
state_groups <- c(list("Cell cycle" = PDO_CELL_CYCLE_MPS), PDO_MP_STATE_GROUPS)

# Use provided updated MP descriptions
mp_desc_map <- setNames(
  paste(names(PDO_MP_DESCRIPTIONS), PDO_MP_DESCRIPTIONS, sep = "_"),
  names(PDO_MP_DESCRIPTIONS)
)

state_cols <- c("Cell cycle" = "#6B7280", PDO_STATE_COLORS[PDO_STATE_ORDER])

# Dotted line positions: draw dotted lines BEFORE these MPs
label_break_before <- vapply(PDO_MP_STATE_GROUPS, `[`, character(1), 1)
####################

strict_refined_mp_order <- unlist(state_groups, use.names = FALSE)
mp_to_state <- unlist(lapply(names(state_groups), function(state) {
  setNames(rep(state, length(state_groups[[state]])), state_groups[[state]])
}))

####################
# Step 04 is the single authoritative filtering point. Its persisted final gene
# list contains only MPs that passed coverage in at least three observed samples
# and the documented five-gene minimum. Do not introduce a second threshold here.
####################
finalized_mps <- intersect(
  unique(as.character(merged_assignments$merged_refined_mp)),
  names(merged_mp_genes)
)

# Validate strict order matches finalized MPs
missing_from_order <- setdiff(finalized_mps, strict_refined_mp_order)
extra_in_order <- setdiff(strict_refined_mp_order, finalized_mps)
if (length(missing_from_order) > 0) {
  cat("WARNING: Finalized MPs not in strict order (will be skipped):", paste(missing_from_order, collapse = ", "), "\n")
}
if (length(extra_in_order) > 0) {
  cat("WARNING: Strict order MPs not in finalized set (will be removed):", paste(extra_in_order, collapse = ", "), "\n")
  strict_refined_mp_order <- strict_refined_mp_order[strict_refined_mp_order %in% finalized_mps]
  # Rebuild state_groups to remove missing MPs
  state_groups <- lapply(state_groups, function(mps) mps[mps %in% finalized_mps])
  state_groups <- state_groups[lengths(state_groups) > 0]
  mp_to_state <- unlist(lapply(names(state_groups), function(state) {
    setNames(rep(state, length(state_groups[[state]])), state_groups[[state]])
  }))
}

if (!all(strict_refined_mp_order %in% colnames(merged_ucell))) {
  stop(
    "Merged UCell cache is missing requested MP(s): ",
    paste(setdiff(strict_refined_mp_order, colnames(merged_ucell)), collapse = ", ")
  )
}

mp_grouping_table <- data.frame(
  state = factor(mp_to_state[strict_refined_mp_order], levels = names(state_groups)),
  mp = strict_refined_mp_order,
  description = unname(mp_desc_map[strict_refined_mp_order]),
  plot_label = paste(strict_refined_mp_order, unname(mp_desc_map[strict_refined_mp_order])),
  stringsAsFactors = FALSE
)
write.csv(
  mp_grouping_table,
  file.path(outdir_live, "tables", "centred_refined_mp_state_grouping.csv"),
  row.names = FALSE
)
cat("Saved: tables/centred_refined_mp_state_grouping.csv\n")

# Save strict order for downstream scripts
saveRDS(
  list(
    strict_refined_mp_order = strict_refined_mp_order,
    state_groups = state_groups,
    state_cols = state_cols,
    mp_desc_map = mp_desc_map,
    mp_to_state = mp_to_state,
    label_break_before = label_break_before
  ),
  file.path(outdir_live, "centred_refined_mp_strict_order.rds")
)
cat("Saved: centred_refined_mp_strict_order.rds\n")

display_label <- function(x) {
  desc <- unname(mp_desc_map[x])
  desc[is.na(desc)] <- x[is.na(desc)]
  ifelse(desc == x, x, paste0(x, ": ", desc))
}

####################
# 1. Ordered NMF program similarity heatmap.
####################
program_similarity <- geneNMF.metaprograms$programs.similarity
original_tree_order <- rownames(program_similarity)[geneNMF.metaprograms$programs.tree$order]

merged_assignments <- merged_assignments |>
  mutate(
    program = as.character(program),
    original_mp = as.character(original_mp),
    refined_submp = as.character(refined_submp),
    merged_refined_mp = as.character(merged_refined_mp)
  ) |>
  filter(program %in% rownames(program_similarity))

# Keep only programs assigned to retained MPs
merged_assignments <- merged_assignments |>
  filter(merged_refined_mp %in% strict_refined_mp_order)

if (anyDuplicated(merged_assignments$program) > 0) {
  duplicated_programs <- unique(merged_assignments$program[duplicated(merged_assignments$program)])
  stop("Programs have duplicate finalized centred refined MP assignments: ",
       paste(head(duplicated_programs, 20), collapse = ", "))
}

program_rank <- setNames(seq_along(original_tree_order), original_tree_order)
ordered_programs <- unlist(lapply(strict_refined_mp_order, function(feature) {
  progs <- merged_assignments$program[merged_assignments$merged_refined_mp == feature]
  progs[order(program_rank[progs])]
}), use.names = FALSE)

missing_programs <- setdiff(merged_assignments$program, ordered_programs)
if (length(missing_programs) > 0) {
  stop("Assigned program(s) were not placed in the ordered heatmap: ",
       paste(head(missing_programs, 20), collapse = ", "))
}

sim_ordered <- program_similarity[ordered_programs, ordered_programs, drop = FALSE]
feature_by_program <- setNames(
  merged_assignments$merged_refined_mp[match(ordered_programs, merged_assignments$program)],
  ordered_programs
)

runs <- rle(unname(feature_by_program))
run_end <- cumsum(runs$lengths)
run_start <- run_end - runs$lengths + 1L

block_table <- data.frame(
  refined_mp = runs$values,
  state = unname(mp_to_state[runs$values]),
  strict_order_rank = match(runs$values, strict_refined_mp_order),
  start = run_start,
  end = run_end,
  n_programs = runs$lengths,
  n_genes = vapply(runs$values, function(x) length(merged_mp_genes[[x]]), integer(1)),
  label = display_label(runs$values),
  stringsAsFactors = FALSE
)

write.csv(
  block_table,
  file.path(outdir_live, "tables", "centred_refined_mp_nmf_ordered_blocks.csv"),
  row.names = FALSE
)

cat("Generating ordered centred refined MP NMF heatmap...\n")
n_prog <- nrow(sim_ordered)
col_fun_nmf <- colorRamp2(c(0, 0.16, 0.42), c("#FFF7F3", "#FB6A4A", "#67000D"))

ht_nmf <- Heatmap(
  sim_ordered,
  name = "Jaccard",
  col = col_fun_nmf,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = NULL,
  column_title = paste0("Centred NMF programmes ordered by final refined MPs (n = ", n_prog, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 4,
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    direction = "horizontal",
    legend_width = unit(8, "cm")
  ),
  width = unit(205, "mm"),
  height = unit(205, "mm")
)

adjust_label_positions <- function(raw_y, labels, base_gap = 0.026,
                                   break_gap = 0.052, low = 0.018, high = 0.982) {
  n <- length(raw_y)
  label_y <- pmax(low, pmin(high, raw_y))
  gap_before <- rep(base_gap, n)
  gap_before[labels %in% label_break_before] <- break_gap
  gap_before[1] <- 0

  for (i in seq(2L, n)) {
    label_y[i] <- min(label_y[i], label_y[i - 1L] - gap_before[i])
  }
  if (label_y[n] < low) {
    label_y <- label_y + (low - label_y[n])
  }
  for (i in seq(n - 1L, 1L)) {
    label_y[i] <- max(label_y[i], label_y[i + 1L] + gap_before[i + 1L])
  }
  if (label_y[1] > high) {
    label_y <- label_y - (label_y[1] - high)
  }
  if (label_y[n] < low) {
    scaled_gaps <- gap_before[-1L]
    scaled_gaps <- scaled_gaps * ((high - low) / sum(scaled_gaps))
    label_y <- c(high, high - cumsum(scaled_gaps))
  }

  pmax(low, pmin(high, label_y))
}

draw_nmf_heatmap <- function() {
  draw(ht_nmf, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 170), "mm"))
  decorate_heatmap_body("Jaccard", {
    raw_y <- 1 - ((block_table$start + block_table$end - 1) / (2 * n_prog))
    label_y <- adjust_label_positions(raw_y, block_table$refined_mp)

    for (rr in seq_len(nrow(block_table))) {
      s <- block_table$start[rr]
      e <- block_table$end[rr]
      centre <- (s + e - 1) / (2 * n_prog)
      extent <- (e - s + 1) / n_prog
      box_y <- 1 - centre
      ly <- label_y[rr]

      grid.rect(
        x = unit(centre, "npc"),
        y = unit(box_y, "npc"),
        width = unit(extent, "npc"),
        height = unit(extent, "npc"),
        gp = gpar(fill = NA, col = "#111111", lwd = 2.2, lty = "dotted")
      )
      grid.lines(
        x = unit(c(centre + extent / 2, 1.008), "npc"),
        y = unit(c(box_y, ly), "npc"),
        gp = gpar(col = "#111111", lwd = 0.55, lty = "dotted")
      )
      grid.text(
        block_table$label[rr],
        x = unit(1.012, "npc"),
        y = unit(ly, "npc"),
        just = "left",
        gp = gpar(fontsize = 8.8, col = "#111111", lineheight = 0.86)
      )
    }
  })
}

nmf_pdf_path <- file.path(outdir_live, "figures", "centred_refined_mp_nmf_ordered_heatmap.pdf")
nmf_png_path <- file.path(outdir_live, "figures", "centred_refined_mp_nmf_ordered_heatmap.png")

cairo_pdf(nmf_pdf_path, width = 15.2, height = 9.8)
draw_nmf_heatmap()
dev.off()

png(nmf_png_path, width = 15.2, height = 9.8, units = "in", res = 450, bg = "white")
draw_nmf_heatmap()
dev.off()

cat("Saved:", nmf_pdf_path, "\n")
cat("Saved:", nmf_png_path, "\n")

####################
# 2. Program-resolution refined MP correlation heatmap.
####################
cat("Computing centred program-resolution refined MP correlations by sample...\n")

compute_fisher_cor_p <- function(score_mat, sample_vec, min_cells = 10) {
  score_mat <- scale(as.matrix(score_mat))
  feature_names <- colnames(score_mat)
  samples <- unique(sample_vec)
  n_features <- length(feature_names)

  cor_array <- array(
    NA_real_,
    dim = c(n_features, n_features, length(samples)),
    dimnames = list(feature_names, feature_names, samples)
  )

  for (samp in samples) {
    idx <- which(sample_vec == samp)
    if (length(idx) < min_cells) next
    cor_array[, , samp] <- cor(score_mat[idx, , drop = FALSE], method = "spearman")
  }

  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
  mean_rho <- matrix(NA_real_, n_features, n_features,
                     dimnames = list(feature_names, feature_names))
  p_vals <- matrix(NA_real_, n_features, n_features,
                   dimnames = list(feature_names, feature_names))

  for (i in seq_len(n_features)) {
    for (j in seq_len(n_features)) {
      if (i == j) {
        mean_rho[i, j] <- 1
        p_vals[i, j] <- 0
        next
      }
      zs <- z_array[i, j, ]
      zs <- zs[is.finite(zs)]
      if (length(zs) < 3) next
      mean_rho[i, j] <- tanh(mean(zs))
      tt <- tryCatch(t.test(zs), error = function(e) NULL)
      p_vals[i, j] <- if (!is.null(tt)) tt$p.value else NA_real_
    }
  }

  list(mean_rho = mean_rho, p_values = p_vals, samples = samples)
}

ordered_split_programs <- unlist(lapply(strict_refined_mp_order, function(feature) {
  feature_rows <- merged_assignments[merged_assignments$merged_refined_mp == feature, , drop = FALSE]
  split_order <- unique(feature_rows$refined_submp[order(program_rank[feature_rows$program])])
  unlist(lapply(split_order, function(split_feature) {
    progs <- feature_rows$program[feature_rows$refined_submp == split_feature]
    progs[order(program_rank[progs])]
  }), use.names = FALSE)
}), use.names = FALSE)

missing_split_programs <- setdiff(merged_assignments$program, ordered_split_programs)
if (length(missing_split_programs) > 0) {
  stop("Assigned program(s) were not placed in the program-resolution correlation heatmap: ",
       paste(head(missing_split_programs, 20), collapse = ", "))
}

ordered_split_assignments <- merged_assignments[
  match(ordered_split_programs, merged_assignments$program),
  ,
  drop = FALSE
]
final_feature_by_split_program <- setNames(
  ordered_split_assignments$merged_refined_mp,
  ordered_split_programs
)
split_feature_by_program <- setNames(
  ordered_split_assignments$refined_submp,
  ordered_split_programs
)

final_split_runs <- rle(unname(final_feature_by_split_program))
final_split_end <- cumsum(final_split_runs$lengths)
final_split_start <- final_split_end - final_split_runs$lengths + 1L

sub_split_runs <- rle(paste(
  unname(final_feature_by_split_program),
  unname(split_feature_by_program),
  sep = "||"
))
sub_split_end <- cumsum(sub_split_runs$lengths)
sub_split_start <- sub_split_end - sub_split_runs$lengths + 1L
sub_split_parts <- strsplit(sub_split_runs$values, "\\|\\|")

split_final_block_table <- data.frame(
  refined_mp = final_split_runs$values,
  state = unname(mp_to_state[final_split_runs$values]),
  strict_order_rank = match(final_split_runs$values, strict_refined_mp_order),
  start = final_split_start,
  end = final_split_end,
  n_programs = final_split_runs$lengths,
  n_split_blocks = vapply(final_split_runs$values, function(x) {
    length(unique(ordered_split_assignments$refined_submp[
      ordered_split_assignments$merged_refined_mp == x
    ]))
  }, integer(1)),
  n_genes = vapply(final_split_runs$values, function(x) length(merged_mp_genes[[x]]), integer(1)),
  label = display_label(final_split_runs$values),
  stringsAsFactors = FALSE
)

split_sub_block_table <- data.frame(
  finalized_mp = vapply(sub_split_parts, `[`, character(1), 1),
  split_mp = vapply(sub_split_parts, `[`, character(1), 2),
  start = sub_split_start,
  end = sub_split_end,
  n_programs = sub_split_runs$lengths,
  stringsAsFactors = FALSE
)

write.csv(
  split_final_block_table,
  file.path(outdir_live, "tables", "centred_refined_mp_split_correlation_ordered_final_blocks.csv"),
  row.names = FALSE
)
write.csv(
  split_sub_block_table,
  file.path(outdir_live, "tables", "centred_refined_mp_split_correlation_ordered_sub_blocks.csv"),
  row.names = FALSE
)

refined_ucell <- readRDS(refined_ucell_file)
display_split_features <- unique(unname(split_feature_by_program))
missing_split_features <- setdiff(display_split_features, colnames(refined_ucell))
if (length(missing_split_features) > 0) {
  stop(
    "refined_ucell_scores.rds is missing split/full feature(s) needed for plotting: ",
    paste(missing_split_features, collapse = ", ")
  )
}

cell_meta <- pdos_merged@meta.data
sample_vec_refined <- cell_meta$orig.ident[match(rownames(refined_ucell), rownames(cell_meta))]
sample_vec_refined[is.na(sample_vec_refined)] <- rownames(refined_ucell)[is.na(sample_vec_refined)]
split_cor <- compute_fisher_cor_p(
  refined_ucell[, display_split_features, drop = FALSE],
  sample_vec_refined,
  min_cells = 10
)
split_cor$display_split_features <- display_split_features

mean_rho_features <- split_cor$mean_rho[display_split_features, display_split_features, drop = FALSE]
split_idx <- match(unname(split_feature_by_program), display_split_features)
rho_program_ordered <- mean_rho_features[split_idx, split_idx, drop = FALSE]
rownames(rho_program_ordered) <- ordered_split_programs
colnames(rho_program_ordered) <- ordered_split_programs

saveRDS(
  list(
    rho_program_resolution = rho_program_ordered,
    ordered_programs = ordered_split_programs,
    final_feature_by_program = final_feature_by_split_program,
    split_feature_by_program = split_feature_by_program,
    final_blocks = split_final_block_table,
    split_blocks = split_sub_block_table,
    feature_correlation = split_cor
  ),
  file.path(outdir_ephemeral, "centred_refined_mp_split_correlation_ordered_matrix.rds")
)

cat("Generating centred program-resolution refined MP correlation heatmap...\n")
n_prog_split <- nrow(rho_program_ordered)
col_fun_split <- colorRamp2(c(-0.4, 0, 0.4), c("#2166AC", "#FFFFFF", "#B2182B"))

ht_split_cor <- Heatmap(
  rho_program_ordered,
  name = "Mean rho",
  col = col_fun_split,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = NULL,
  column_title = paste0("Centred refined MP correlation at NMF-program resolution (n = ", n_prog_split, ")"),
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 4,
  border = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14),
    direction = "horizontal",
    legend_width = unit(8, "cm")
  ),
  width = unit(205, "mm"),
  height = unit(205, "mm")
)

draw_split_correlation_heatmap <- function() {
  draw(ht_split_cor, heatmap_legend_side = "bottom", padding = unit(c(12, 10, 14, 170), "mm"))
  decorate_heatmap_body("Mean rho", {
    for (rr in seq_len(nrow(split_sub_block_table))) {
      s <- split_sub_block_table$start[rr]
      e <- split_sub_block_table$end[rr]
      centre <- (s + e - 1) / (2 * n_prog_split)
      extent <- (e - s + 1) / n_prog_split
      grid.rect(
        x = unit(centre, "npc"),
        y = unit(1 - centre, "npc"),
        width = unit(extent, "npc"),
        height = unit(extent, "npc"),
        gp = gpar(fill = NA, col = "#444444", lwd = 1.0, lty = "solid")
      )
    }

    raw_y <- 1 - ((split_final_block_table$start + split_final_block_table$end - 1) / (2 * n_prog_split))
    label_y <- adjust_label_positions(raw_y, split_final_block_table$refined_mp)

    for (rr in seq_len(nrow(split_final_block_table))) {
      s <- split_final_block_table$start[rr]
      e <- split_final_block_table$end[rr]
      centre <- (s + e - 1) / (2 * n_prog_split)
      extent <- (e - s + 1) / n_prog_split
      box_y <- 1 - centre
      ly <- label_y[rr]

      grid.rect(
        x = unit(centre, "npc"),
        y = unit(box_y, "npc"),
        width = unit(extent, "npc"),
        height = unit(extent, "npc"),
        gp = gpar(fill = NA, col = "#111111", lwd = 2.2, lty = "dotted")
      )
      grid.lines(
        x = unit(c(centre + extent / 2, 1.008), "npc"),
        y = unit(c(box_y, ly), "npc"),
        gp = gpar(col = "#111111", lwd = 0.55, lty = "dotted")
      )
      grid.text(
        split_final_block_table$label[rr],
        x = unit(1.012, "npc"),
        y = unit(ly, "npc"),
        just = "left",
        gp = gpar(fontsize = 8.8, col = "#111111", lineheight = 0.86)
      )
    }
  })
}

split_cor_pdf_path <- file.path(outdir_live, "figures", "centred_refined_mp_split_correlation_ordered_heatmap.pdf")
split_cor_png_path <- file.path(outdir_live, "figures", "centred_refined_mp_split_correlation_ordered_heatmap.png")

cairo_pdf(split_cor_pdf_path, width = 15.2, height = 9.8)
draw_split_correlation_heatmap()
dev.off()

png(split_cor_png_path, width = 15.2, height = 9.8, units = "in", res = 450, bg = "white")
draw_split_correlation_heatmap()
dev.off()

cat("Saved:", split_cor_pdf_path, "\n")
cat("Saved:", split_cor_png_path, "\n")

####################
# 3. MP UCell Fisher-Z correlation heatmap with state annotations.
####################
cat("Computing centred refined MP correlations by sample...\n")

ucell_scores <- merged_ucell[, strict_refined_mp_order, drop = FALSE]
module_scores <- scale(as.matrix(ucell_scores))
score_mat <- module_scores

cell_meta <- pdos_merged@meta.data
cell_ids <- rownames(score_mat)
sample_vec <- cell_meta$orig.ident[match(cell_ids, rownames(cell_meta))]
sample_vec[is.na(sample_vec)] <- cell_ids[is.na(sample_vec)]
samples <- unique(sample_vec)
mps <- colnames(score_mat)
n_mps <- length(mps)

cor_array <- array(
  NA_real_,
  dim = c(n_mps, n_mps, length(samples)),
  dimnames = list(mps, mps, samples)
)

for (samp in samples) {
  cells_in_sample <- which(sample_vec == samp)
  if (length(cells_in_sample) < 10) next
  sub_mat <- score_mat[cells_in_sample, , drop = FALSE]
  cor_array[, , samp] <- cor(sub_mat, method = "spearman")
}

z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
mean_rho <- matrix(NA_real_, n_mps, n_mps, dimnames = list(mps, mps))
p_vals <- matrix(NA_real_, n_mps, n_mps, dimnames = list(mps, mps))

for (i in seq_len(n_mps)) {
  for (j in seq_len(n_mps)) {
    if (i == j) {
      mean_rho[i, j] <- 1
      p_vals[i, j] <- 0
      next
    }
    z_scores <- z_array[i, j, ]
    z_scores <- z_scores[is.finite(z_scores)]
    if (length(z_scores) < 3) next
    mean_rho[i, j] <- tanh(mean(z_scores))
    test_res <- tryCatch(t.test(z_scores), error = function(e) NULL)
    p_vals[i, j] <- if (!is.null(test_res)) test_res$p.value else NA_real_
  }
}

write.csv(
  mean_rho,
  file.path(outdir_live, "tables", "centred_refined_mp_correlation_mean_rho.csv")
)
write.csv(
  p_vals,
  file.path(outdir_live, "tables", "centred_refined_mp_correlation_p_values.csv")
)

mp_names_with_desc <- setNames(display_label(strict_refined_mp_order), strict_refined_mp_order)
mean_rho_plot <- mean_rho
rownames(mean_rho_plot) <- mp_names_with_desc[rownames(mean_rho_plot)]
colnames(mean_rho_plot) <- mp_names_with_desc[colnames(mean_rho_plot)]

state_vec_for_mps <- factor(mp_to_state[strict_refined_mp_order], levels = names(state_groups))

ha_left <- rowAnnotation(
  State = state_vec_for_mps,
  col = list(State = state_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)
ha_top <- HeatmapAnnotation(
  State = state_vec_for_mps,
  col = list(State = state_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm_width <- unit(10.5, "inch")
hm_height <- unit(10.5, "inch")
col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

draw_correlation_heatmap <- function() {
  ht_cor <- Heatmap(
    mean_rho_plot,
    name = paste0("Mean Rho\n(", sum(apply(cor_array, 3, function(x) any(is.finite(x)))), " Samples)"),
    col = col_cor,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    row_split = state_vec_for_mps,
    column_split = state_vec_for_mps,
    column_title_rot = 20,
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_title = NULL,
    row_names_side = "right",
    column_names_side = "bottom",
    column_names_rot = 30,
    row_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
    row_names_max_width = unit(128, "mm"),
    column_names_max_height = unit(128, "mm"),
    width = hm_width,
    height = hm_height,
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- p_vals[i, j]
      rho <- mean_rho[i, j]
      if (is.na(p) || is.na(rho)) {
        grid.text("NA", x, y, gp = gpar(fontsize = 8.5, col = "grey50"))
      } else if (p < 0.001) {
        grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 8.5))
      } else if (p < 0.01) {
        grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 8.5))
      } else if (p < 0.05) {
        grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 8.5))
      } else {
        grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 8.5))
      }
    },
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  )
  draw(
    ht_cor,
    heatmap_legend_side = "left",
    padding = unit(c(20, 20, 20, 20), "mm")
  )
}

cor_pdf_path <- file.path(outdir_live, "figures", "centred_refined_mp_correlation_heatmap.pdf")
cor_png_path <- file.path(outdir_live, "figures", "centred_refined_mp_correlation_heatmap.png")

cairo_pdf(cor_pdf_path, width = 22, height = 18)
draw_correlation_heatmap()
dev.off()

png(cor_png_path, width = 22, height = 18, units = "in", res = 350, bg = "white")
draw_correlation_heatmap()
dev.off()

saveRDS(
  list(
    strict_refined_mp_order = strict_refined_mp_order,
    state_groups = state_groups,
    state_cols = state_cols,
    mp_desc_map = mp_desc_map,
    ordered_program_similarity = sim_ordered,
    ordered_programs = ordered_programs,
    feature_by_program = feature_by_program,
    nmf_blocks = block_table,
    correlation_mean_rho = mean_rho,
    correlation_p_values = p_vals,
    split_correlation_program_resolution = rho_program_ordered,
    split_correlation_final_blocks = split_final_block_table,
    split_correlation_sub_blocks = split_sub_block_table,
    split_correlation_features = split_cor,
    sample_count = length(samples)
  ),
  file.path(outdir_ephemeral, "centred_refined_mp_ordered_heatmaps.rds")
)

cat("Saved:", cor_pdf_path, "\n")
cat("Saved:", cor_png_path, "\n")
cat("Saved:", split_cor_pdf_path, "\n")
cat("Saved:", split_cor_png_path, "\n")
cat("Saved: tables/centred_refined_mp_state_grouping.csv\n")
cat("Saved: intermediate/centred_refined_mp_ordered_heatmaps.rds\n")
cat("\n=== Auto_05 PDO Ordered Heatmaps Complete ===\n")
