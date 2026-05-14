####################
# Analysis registry:
#   Status: legacy comparison; no downstream use
#   Script: analysis/cell_states/legacy_compare_mp_scoring_state_definition.R
#   Recommended legacy name: analysis/cell_states/legacy_compare_mp_scoring_state_definition.R
#   Methodology: analysis/methodology/cell_states/state_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     PDOs_outs/PDOs_merged.rds
#     PDOs_outs/MP_outs_default.rds or Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#     optional PDOs_outs/UCell_scores_filtered.rds
#   Outputs:
#     PDOs_outs/Auto_compare_mp_scoring_state_definition/*
#   Downstream:
#     Comparison-only. Do not use alternative state vectors from this script
#     as current downstream state definitions.
####################

####################
# legacy_compare_mp_scoring_state_definition.R
# Compare PDO MP activity/state calls from full UCell, cumulative-weight UCell,
# and weighted-rank scoring using the optimal PDO GeneNMF metaprogram object.
####################

library(Seurat)
library(UCell)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

get_env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) default else parsed
}

get_env_integer <- function(name, default) {
  as.integer(round(get_env_numeric(name, default)))
}

cum_weight_threshold <- get_env_numeric("AUTO_MP_CUM_WEIGHT_THRESHOLD", 0.70)
max_rank <- get_env_integer("AUTO_MP_MAX_RANK", 1500)
ncores <- get_env_integer("AUTO_MP_NCORES", 4)
plot_max_cells <- get_env_integer("AUTO_MP_PLOT_MAX_CELLS", 12000)
state_threshold <- get_env_numeric("AUTO_MP_STATE_THRESHOLD", 0.5)
hybrid_gap <- get_env_numeric("AUTO_MP_HYBRID_GAP", 0.3)
random_seed <- get_env_integer("AUTO_MP_SEED", 42)

if (cum_weight_threshold <= 0 || cum_weight_threshold > 1) {
  stop("AUTO_MP_CUM_WEIGHT_THRESHOLD must be in (0, 1].")
}
if (max_rank < 10) stop("AUTO_MP_MAX_RANK is too small.")
if (ncores < 1) ncores <- 1
if (plot_max_cells < 100) plot_max_cells <- 100

set.seed(random_seed)

out_dir <- "Auto_compare_mp_scoring_state_definition"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("=== Loading PDO data and metaprograms ===")

mp_path <- Sys.getenv("AUTO_MP_OBJECT", unset = "")
if (mp_path == "") {
  if (file.exists("MP_outs_default.rds")) {
    mp_path <- "MP_outs_default.rds"
  } else {
    mp_path <- file.path("Metaprogrammes_Results", "geneNMF_metaprograms_nMP_13.rds")
  }
}
if (!file.exists(mp_path)) stop("Metaprogram object not found: ", mp_path)

geneNMF.metaprograms <- readRDS(mp_path)
pdos <- readRDS("PDOs_merged.rds")

if ("orig.ident" %in% colnames(pdos@meta.data)) {
  pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
}
pdos$Batch <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
  mat[, Cells(obj), drop = FALSE]
}

expression_matrix <- get_assay_matrix(pdos, assay = "RNA", layer = "data")

message("Cells retained after SUR843T3_PDO exclusion: ", ncol(expression_matrix))
message("Genes in expression matrix: ", nrow(expression_matrix))

####################
# MP filtering and gene-list construction
####################

if (is.null(geneNMF.metaprograms$metaprograms.genes)) {
  stop("The metaprogram object has no metaprograms.genes slot.")
}
if (is.null(geneNMF.metaprograms$metaprograms.genes.weights)) {
  stop("The metaprogram object has no metaprograms.genes.weights slot.")
}

all_mp_genes <- geneNMF.metaprograms$metaprograms.genes
all_mp_weights <- geneNMF.metaprograms$metaprograms.genes.weights

bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)

coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
exclude_mps <- unique(c(bad_mp_names, low_coverage_mps))

if (length(exclude_mps) > 0) {
  message("Excluding MPs by silhouette/sample coverage: ", paste(exclude_mps, collapse = ", "))
}

retained_mps <- setdiff(names(all_mp_genes), exclude_mps)
mp_genes_full <- all_mp_genes[retained_mps]
mp_weights <- all_mp_weights[retained_mps]
mp_weights <- mp_weights[names(mp_genes_full)]

mp_weights <- Map(function(weights, genes) {
  weights <- weights[is.finite(weights) & weights > 0]
  weights[intersect(names(weights), genes)]
}, mp_weights, mp_genes_full)

select_cumulative_weight_genes <- function(weights, threshold) {
  weights <- weights[is.finite(weights) & weights > 0]
  if (length(weights) == 0) return(character(0))
  weights <- sort(weights, decreasing = TRUE)
  cumulative_weight <- cumsum(weights / sum(weights))
  keep_n <- which(cumulative_weight >= threshold)[1]
  names(weights)[seq_len(keep_n)]
}

mp_genes_cum <- lapply(mp_weights, select_cumulative_weight_genes, threshold = cum_weight_threshold)

clean_gene_list <- function(gene_list, available_genes, min_genes = 3) {
  cleaned <- lapply(gene_list, function(genes) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    intersect(genes, available_genes)
  })
  too_short <- names(cleaned)[lengths(cleaned) < min_genes]
  if (length(too_short) > 0) {
    stop("Too few genes present for: ", paste(too_short, collapse = ", "))
  }
  cleaned
}

mp_genes_full_present <- clean_gene_list(mp_genes_full, rownames(expression_matrix))
mp_genes_cum_present <- clean_gene_list(mp_genes_cum, rownames(expression_matrix))

mp_descriptions <- c(
  "MP6"  = "MP6 G2M Cell Cycle",
  "MP7"  = "MP7 DNA repair",
  "MP5"  = "MP5 MYC-related Proliferation",
  "MP1"  = "MP1 G2M checkpoint",
  "MP3"  = "MP3 G1S Cell Cycle",
  "MP8"  = "MP8 Columnar Progenitor",
  "MP10" = "MP10 Inflammatory Stress Epi.",
  "MP9"  = "MP9 ECM Remodeling Epi.",
  "MP4"  = "MP4 Intestinal Metaplasia"
)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- rev(unique(ordered_clusters))
mp_tree_order_names <- paste0("MP", mp_tree_order)
mp_order <- mp_tree_order_names[mp_tree_order_names %in% retained_mps]
mp_order <- c(mp_order, setdiff(retained_mps, mp_order))

weight_summary <- bind_rows(lapply(retained_mps, function(mp) {
  weights <- sort(mp_weights[[mp]], decreasing = TRUE)
  selected <- mp_genes_cum[[mp]]
  selected_present <- mp_genes_cum_present[[mp]]
  data.frame(
    MP = mp,
    silhouette = geneNMF.metaprograms$metaprograms.metrics[mp, "silhouette"],
    sampleCoverage = geneNMF.metaprograms$metaprograms.metrics[mp, "sampleCoverage"],
    full_gene_n = length(mp_genes_full[[mp]]),
    full_present_gene_n = length(mp_genes_full_present[[mp]]),
    cumulative_weight_threshold = cum_weight_threshold,
    cumulative_gene_n = length(selected),
    cumulative_present_gene_n = length(selected_present),
    cumulative_weight_sum = sum(weights[selected], na.rm = TRUE) / sum(weights, na.rm = TRUE),
    min_selected_weight = min(weights[selected], na.rm = TRUE),
    max_selected_weight = max(weights[selected], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

weight_long <- bind_rows(lapply(retained_mps, function(mp) {
  weights <- sort(mp_weights[[mp]], decreasing = TRUE)
  data.frame(
    MP = mp,
    gene = names(weights),
    weight = as.numeric(weights),
    weight_rank = seq_along(weights),
    cumulative_weight = cumsum(weights / sum(weights)),
    selected_cumulative = names(weights) %in% mp_genes_cum[[mp]],
    present_in_expression_matrix = names(weights) %in% rownames(expression_matrix),
    stringsAsFactors = FALSE
  )
}))

write.csv(weight_summary, file.path(out_dir, "Auto_mp_cumulative_weight_gene_counts.csv"), row.names = FALSE)
write.csv(weight_long, file.path(out_dir, "Auto_mp_gene_weight_selection_long.csv"), row.names = FALSE)

message("Cumulative weight threshold: ", cum_weight_threshold)
message("Cumulative genes retained per MP: ",
        paste(paste0(weight_summary$MP, "=", weight_summary$cumulative_present_gene_n), collapse = ", "))

####################
# Activity scoring
####################

method_labels <- c(
  full_ucell = "Full gene-list UCell",
  cum_weight_ucell = paste0("Cumulative ", round(cum_weight_threshold * 100), "% weight UCell"),
  weighted_rank = "Weighted rank score"
)

cached_scores_path <- file.path(out_dir, "Auto_mp_activity_scores_all_methods.rds")
use_cache <- FALSE
if (file.exists(cached_scores_path)) {
  score_list_cached <- readRDS(cached_scores_path)
  cached_mps <- sort(colnames(score_list_cached[[1]]))
  needed_mps <- sort(mp_order)
  if (identical(cached_mps, needed_mps) &&
      all(c("full_ucell", "cum_weight_ucell", "weighted_rank") %in% names(score_list_cached))) {
    message("=== Using cached activity scores (MPs match) ===")
    score_list <- score_list_cached
    use_cache <- TRUE
    rm(score_list_cached)
  } else {
    message("Cached scores exist but MPs differ — recomputing.")
    rm(score_list_cached)
  }
}

if (!use_cache) {
  message("=== Precomputing UCell rank matrix ===")
  rankings <- StoreRankings_UCell(
    expression_matrix,
    maxRank = max_rank,
    ncores = ncores,
    ties.method = "average",
    force.gc = TRUE
  )

  score_ucell_from_ranks <- function(rankings, features, max_rank) {
    scores <- ScoreSignatures_UCell(
      precalc.ranks = rankings,
      features = features,
      maxRank = max_rank,
      name = ""
    )
    scores <- as.data.frame(scores)
    missing_cols <- setdiff(names(features), colnames(scores))
    if (length(missing_cols) > 0) {
      stop("UCell did not return scores for: ", paste(missing_cols, collapse = ", "))
    }
    scores <- scores[, names(features), drop = FALSE]
    as.matrix(scores)
  }

  weighted_rank_score_from_ranks <- function(rankings, weight_list, max_rank) {
    all_weighted_genes <- unique(unlist(lapply(weight_list, names), use.names = FALSE))
    genes_present <- intersect(all_weighted_genes, rownames(rankings))
    if (length(genes_present) == 0) stop("No weighted genes found in rank matrix.")

    rank_sub <- rankings[genes_present, , drop = FALSE]
    rank_sub@x <- (max_rank - rank_sub@x) / max_rank
    rank_sub@x[rank_sub@x < 0] <- 0

    weight_i <- integer(0)
    weight_j <- integer(0)
    weight_x <- numeric(0)
    mp_names <- names(weight_list)

    for (j in seq_along(weight_list)) {
      weights <- weight_list[[j]]
      weights <- weights[intersect(names(weights), genes_present)]
      weights <- weights[is.finite(weights) & weights > 0]
      if (length(weights) < 3) {
        stop("Too few weighted genes present for ", mp_names[j])
      }
      weights <- weights / sum(weights)
      weight_i <- c(weight_i, match(names(weights), genes_present))
      weight_j <- c(weight_j, rep(j, length(weights)))
      weight_x <- c(weight_x, as.numeric(weights))
    }

    weight_matrix <- sparseMatrix(
      i = weight_i,
      j = weight_j,
      x = weight_x,
      dims = c(length(genes_present), length(mp_names)),
      dimnames = list(genes_present, mp_names)
    )

    scores <- Matrix::t(rank_sub) %*% weight_matrix
    scores <- as.matrix(scores)
    scores[, mp_names, drop = FALSE]
  }

  message("=== Scoring full gene-list UCell ===")
  scores_full_ucell <- score_ucell_from_ranks(rankings, mp_genes_full_present, max_rank)

  message("=== Scoring cumulative-weight gene-list UCell ===")
  scores_cum_ucell <- score_ucell_from_ranks(rankings, mp_genes_cum_present, max_rank)

  message("=== Scoring weighted-rank activity ===")
  scores_weighted_rank <- weighted_rank_score_from_ranks(rankings, mp_weights, max_rank)

  score_list <- list(
    full_ucell = scores_full_ucell,
    cum_weight_ucell = scores_cum_ucell,
    weighted_rank = scores_weighted_rank
  )
  score_list <- lapply(score_list, function(mat) mat[, mp_order, drop = FALSE])

  saveRDS(score_list, cached_scores_path)
  rm(rankings)
  gc(verbose = FALSE)
}

if (!use_cache && file.exists("UCell_scores_filtered.rds")) {
  message("=== Auditing newly computed full UCell scores against UCell_scores_filtered.rds ===")
  scores_full_ucell <- score_list[["full_ucell"]]
  old_ucell <- readRDS("UCell_scores_filtered.rds")
  audit_cells <- intersect(rownames(old_ucell), rownames(scores_full_ucell))
  audit_mps <- intersect(colnames(old_ucell), colnames(scores_full_ucell))
  existing_audit <- bind_rows(lapply(audit_mps, function(mp) {
    x <- old_ucell[audit_cells, mp]
    y <- scores_full_ucell[audit_cells, mp]
    data.frame(
      MP = mp,
      n_cells = length(audit_cells),
      pearson_r = suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs")),
      spearman_r = suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs")),
      max_abs_difference = max(abs(x - y), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(existing_audit, file.path(out_dir, "Auto_existing_full_ucell_audit.csv"), row.names = FALSE)
}

####################
# Activity concordance plots and statistics
####################

scale_vector <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

calc_activity_stats <- function(method_x, method_y, mp) {
  x <- score_list[[method_x]][, mp]
  y <- score_list[[method_y]][, mp]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  top_n <- max(1, ceiling(0.10 * length(x)))
  x_top <- rank(-x, ties.method = "first") <= top_n
  y_top <- rank(-y, ties.method = "first") <= top_n
  top_union <- x_top | y_top
  top_intersection <- x_top & y_top

  data.frame(
    method_x = method_x,
    method_y = method_y,
    MP = mp,
    n_cells = length(x),
    pearson_r = suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs")),
    spearman_r = suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs")),
    top10_jaccard = ifelse(sum(top_union) == 0, NA_real_, sum(top_intersection) / sum(top_union)),
    top10_overlap_of_x = ifelse(sum(x_top) == 0, NA_real_, sum(top_intersection) / sum(x_top)),
    top10_overlap_of_y = ifelse(sum(y_top) == 0, NA_real_, sum(top_intersection) / sum(y_top)),
    mean_abs_z_difference = mean(abs(scale_vector(x) - scale_vector(y)), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

method_pairs <- combn(names(score_list), 2, simplify = FALSE)

activity_stats <- bind_rows(lapply(method_pairs, function(pair) {
  bind_rows(lapply(mp_order, function(mp) calc_activity_stats(pair[1], pair[2], mp)))
}))
write.csv(activity_stats, file.path(out_dir, "Auto_activity_pairwise_correlation_stats.csv"), row.names = FALSE)

make_activity_scatter <- function(method_x, method_y) {
  common_cells <- intersect(rownames(score_list[[method_x]]), rownames(score_list[[method_y]]))
  if (length(common_cells) > plot_max_cells) {
    common_cells <- sample(common_cells, plot_max_cells)
  }

  plot_df <- bind_rows(lapply(mp_order, function(mp) {
    data.frame(
      cell = common_cells,
      MP = mp,
      x = score_list[[method_x]][common_cells, mp],
      y = score_list[[method_y]][common_cells, mp],
      stringsAsFactors = FALSE
    )
  }))

  label_levels <- ifelse(is.na(mp_descriptions[mp_order]), mp_order, mp_descriptions[mp_order])
  names(label_levels) <- mp_order
  plot_df$MP_label <- factor(label_levels[plot_df$MP], levels = label_levels)

  annot_df <- activity_stats %>%
    filter(method_x == !!method_x, method_y == !!method_y) %>%
    mutate(
      MP_label = factor(label_levels[MP], levels = label_levels),
      label = paste0(
        "rho=", sprintf("%.2f", spearman_r),
        "\nJ10=", sprintf("%.2f", top10_jaccard)
      )
    )

  ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(size = 0.12, alpha = 0.18, color = "#2B2B2B") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.35, color = "#B2182B") +
    geom_text(
      data = annot_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.08,
      size = 2.7
    ) +
    facet_wrap(~ MP_label, scales = "free", ncol = 3) +
    labs(
      title = paste(method_labels[[method_x]], "vs", method_labels[[method_y]]),
      x = method_labels[[method_x]],
      y = method_labels[[method_y]]
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 8.5),
      panel.grid.minor = element_blank()
    )
}

pdf(file.path(out_dir, "Auto_activity_pairwise_scatter.pdf"), width = 13, height = 10, useDingbats = FALSE)
for (pair in method_pairs) {
  print(make_activity_scatter(pair[1], pair[2]))
}
dev.off()

####################
# State assignment from each scoring method
####################

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
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

cc_mps <- intersect(c("MP6", "MP7", "MP1", "MP3"), retained_mps)
non_cc_mps <- setdiff(retained_mps, cc_mps)

state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "SMG-like Metaplasia"   = c("MP8")
)
state_groups <- lapply(state_groups, function(mps) intersect(mps, non_cc_mps))
state_groups <- state_groups[lengths(state_groups) > 0]

state_level_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "Unresolved",
  "Hybrid"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Stress-adaptive"       = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Immune Infiltrating"   = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

sample_var <- setNames(as.character(pdos$orig.ident), Cells(pdos))
study_var <- setNames(as.character(pdos$Batch), Cells(pdos))

assign_states_from_scores <- function(scores, method_name) {
  common_cells <- intersect(rownames(scores), Cells(pdos))
  scores <- scores[common_cells, retained_mps, drop = FALSE]

  cc_present <- intersect(cc_mps, colnames(scores))
  non_cc_present <- intersect(non_cc_mps, colnames(scores))
  if (length(non_cc_present) == 0) stop("No non-cell-cycle MPs available for ", method_name)

  mp_adj_noncc <- z_normalise(scores[, non_cc_present, drop = FALSE], sample_var, study_var)
  mp_adj_cc <- NULL
  if (length(cc_present) > 0) {
    mp_adj_cc <- z_normalise(scores[, cc_present, drop = FALSE], sample_var, study_var)
  }
  mp_adj_all <- cbind(mp_adj_noncc, mp_adj_cc)

  group_max <- sapply(state_groups, function(mps) {
    mps_avail <- intersect(mps, colnames(mp_adj_noncc))
    if (length(mps_avail) == 0) return(rep(0, nrow(mp_adj_noncc)))
    if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc[, mps_avail]))
    apply(mp_adj_noncc[, mps_avail, drop = FALSE], 1, max)
  })
  group_max <- as.matrix(group_max)
  rownames(group_max) <- rownames(mp_adj_noncc)

  best_group_idx <- max.col(group_max, ties.method = "first")
  best_group_val <- apply(group_max, 1, max)
  base_state <- colnames(group_max)[best_group_idx]
  base_state[best_group_val < state_threshold] <- "Unresolved"

  state_vec <- base_state
  sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
  gap <- sorted_groups[, 1] - sorted_groups[, 2]
  state_vec[(gap < hybrid_gap) & (base_state != "Unresolved")] <- "Hybrid"
  names(state_vec) <- rownames(group_max)

  top_mp <- colnames(mp_adj_all)[max.col(mp_adj_all, ties.method = "first")]
  names(top_mp) <- rownames(mp_adj_all)

  list(
    state = state_vec,
    top_mp = top_mp,
    group_max = group_max,
    mp_adj_noncc = mp_adj_noncc,
    mp_adj_all = mp_adj_all
  )
}

message("=== Assigning states for each scoring method ===")
state_results <- lapply(names(score_list), function(method_name) {
  message("State assignment: ", method_name)
  assign_states_from_scores(score_list[[method_name]], method_name)
})
names(state_results) <- names(score_list)

for (method_name in names(state_results)) {
  saveRDS(
    state_results[[method_name]]$state,
    file.path(out_dir, paste0("Auto_PDO_states_", method_name, ".rds"))
  )
  saveRDS(
    state_results[[method_name]]$top_mp,
    file.path(out_dir, paste0("Auto_PDO_top_mp_", method_name, ".rds"))
  )
}
saveRDS(state_results, file.path(out_dir, "Auto_state_assignment_results_all_methods.rds"))

state_summary <- bind_rows(lapply(names(state_results), function(method_name) {
  state_vec <- state_results[[method_name]]$state
  as.data.frame(table(state = state_vec), stringsAsFactors = FALSE) %>%
    mutate(
      method = method_name,
      method_label = method_labels[[method_name]],
      n_cells = sum(Freq),
      pct = 100 * Freq / n_cells
    ) %>%
    rename(n = Freq)
}))
write.csv(state_summary, file.path(out_dir, "Auto_state_assignment_summary.csv"), row.names = FALSE)

composition_df <- state_summary %>%
  mutate(
    state = factor(as.character(state), levels = state_level_order),
    method_label = factor(method_label, levels = method_labels[names(score_list)])
  )

p_comp <- ggplot(composition_df, aes(x = method_label, y = pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 3, sprintf("%.1f", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.8) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  labs(x = NULL, y = "% of cells", fill = "State", title = "PDO state composition by MP scoring method") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.major.x = element_blank())

####################
# Connected-points plot: 4 main states × 3 methods, per-sample paired data
####################

message("=== Building connected-points (paired sample) state proportion plot ===")

four_main_states <- c("Classic Proliferative", "Basal to Intest. Meta",
                      "SMG-like Metaplasia", "Stress-adaptive")

paired_df <- bind_rows(lapply(names(state_results), function(method_name) {
  state_vec <- state_results[[method_name]]$state
  sample_vec <- sample_var[names(state_vec)]
  data.frame(
    method = method_name,
    method_label = method_labels[[method_name]],
    sample = sample_vec,
    state = state_vec,
    stringsAsFactors = FALSE
  )
})) %>%
  filter(state %in% four_main_states) %>%
  count(method, method_label, sample, state, name = "n_state")

sample_method_totals <- bind_rows(lapply(names(state_results), function(method_name) {
  state_vec <- state_results[[method_name]]$state
  sample_vec <- sample_var[names(state_vec)]
  data.frame(method = method_name, sample = sample_vec, stringsAsFactors = FALSE)
})) %>%
  count(method, sample, name = "n_total")

paired_df <- paired_df %>%
  left_join(sample_method_totals, by = c("method", "sample")) %>%
  mutate(pct = 100 * n_state / n_total)

all_combos <- expand.grid(
  method = unique(paired_df$method),
  sample = unique(paired_df$sample),
  state = four_main_states,
  stringsAsFactors = FALSE
)
all_combos$method_label <- method_labels[all_combos$method]
paired_df <- all_combos %>%
  left_join(paired_df %>% select(-method_label), by = c("method", "sample", "state")) %>%
  left_join(sample_method_totals, by = c("method", "sample"), suffix = c("", ".y")) %>%
  mutate(
    n_total = coalesce(n_total, n_total.y),
    n_state = ifelse(is.na(n_state), 0L, n_state),
    pct = ifelse(is.na(pct), 0, pct)
  ) %>%
  select(-n_total.y)

paired_df$state <- factor(paired_df$state, levels = four_main_states)
paired_df$method_label <- factor(paired_df$method_label, levels = method_labels[names(score_list)])

method_colors <- c(
  "Full gene-list UCell"               = "#E41A1C",
  "Cumulative 70% weight UCell"        = "#4DAF4A",
  "Weighted rank score"                = "#377EB8"
)
names(method_colors) <- method_labels[names(score_list)]

p_paired <- ggplot(paired_df, aes(x = method_label, y = pct, group = sample)) +
  geom_line(alpha = 0.35, linewidth = 0.4, color = "grey40") +
  geom_point(aes(color = method_label), size = 2.2, alpha = 0.85) +
  facet_wrap(~ state, ncol = 4, scales = "free_y") +
  scale_color_manual(values = method_colors, guide = "none") +
  labs(
    title = "Per-sample state proportion across scoring methods",
    x = NULL, y = "% of cells in sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

####################
# Per-sample state composition (20 pages)
####################

message("=== Building per-sample state composition plots ===")

sample_state_summary <- bind_rows(lapply(names(state_results), function(method_name) {
  state_vec <- state_results[[method_name]]$state
  sample_vec <- sample_var[names(state_vec)]
  data.frame(method = method_name, method_label = method_labels[[method_name]],
             sample = sample_vec, state = state_vec, stringsAsFactors = FALSE)
})) %>%
  count(method, method_label, sample, state, name = "n") %>%
  group_by(method, method_label, sample) %>%
  mutate(n_total = sum(n), pct = 100 * n / n_total) %>%
  ungroup() %>%
  mutate(
    state = factor(as.character(state), levels = state_level_order),
    method_label = factor(method_label, levels = method_labels[names(score_list)])
  )

all_samples <- sort(unique(as.character(sample_state_summary$sample)))

pdf(file.path(out_dir, "Auto_state_composition_by_method.pdf"), width = 9, height = 6, useDingbats = FALSE)
print(p_paired)
print(p_comp)
for (samp in all_samples) {
  sub_df <- sample_state_summary %>% filter(sample == samp)
  p_samp <- ggplot(sub_df, aes(x = method_label, y = pct, fill = state)) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(aes(label = ifelse(pct >= 3, sprintf("%.1f", pct), "")),
              position = position_stack(vjust = 0.5), size = 2.8) +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    labs(x = NULL, y = "% of cells", fill = "State",
         title = paste0("State composition: ", samp)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold"))
  print(p_samp)
}
dev.off()

####################
# State concordance plots and statistics
####################

comb2 <- function(x) {
  x <- as.numeric(x)
  x * (x - 1) / 2
}

adjusted_rand_from_table <- function(tab) {
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  sum_comb <- sum(comb2(tab))
  row_comb <- sum(comb2(rowSums(tab)))
  col_comb <- sum(comb2(colSums(tab)))
  total_comb <- comb2(n)
  expected <- row_comb * col_comb / total_comb
  max_index <- (row_comb + col_comb) / 2
  if (isTRUE(all.equal(max_index, expected))) return(NA_real_)
  (sum_comb - expected) / (max_index - expected)
}

cohen_kappa_from_vectors <- function(x, y, levels) {
  tab <- table(factor(x, levels = levels), factor(y, levels = levels))
  n <- as.numeric(sum(tab))
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(as.numeric(rowSums(tab)) * as.numeric(colSums(tab))) / (n * n)
  if (isTRUE(all.equal(1, pe))) return(NA_real_)
  (po - pe) / (1 - pe)
}

state_pair_stats <- bind_rows(lapply(method_pairs, function(pair) {
  state_x <- state_results[[pair[1]]]$state
  state_y <- state_results[[pair[2]]]$state
  cells <- intersect(names(state_x), names(state_y))
  state_x <- state_x[cells]
  state_y <- state_y[cells]
  levels_pair <- c(state_level_order, sort(setdiff(unique(c(state_x, state_y)), state_level_order)))
  tab <- table(factor(state_x, levels = levels_pair), factor(state_y, levels = levels_pair))
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  row_max <- apply(sweep(tab, 1, rowSums(tab), FUN = "/"), 1, max, na.rm = TRUE)

  data.frame(
    method_x = pair[1],
    method_y = pair[2],
    method_x_label = method_labels[[pair[1]]],
    method_y_label = method_labels[[pair[2]]],
    n_cells = length(cells),
    exact_match = mean(state_x == state_y),
    adjusted_rand_index = adjusted_rand_from_table(tab),
    cohen_kappa = cohen_kappa_from_vectors(state_x, state_y, levels_pair),
    weighted_row_max_concordance = sum(rowSums(tab) * row_max) / sum(tab),
    mean_row_max_concordance = mean(row_max),
    stringsAsFactors = FALSE
  )
}))
write.csv(state_pair_stats, file.path(out_dir, "Auto_state_pairwise_concordance_stats.csv"), row.names = FALSE)

state_counts_long <- bind_rows(lapply(method_pairs, function(pair) {
  state_x <- state_results[[pair[1]]]$state
  state_y <- state_results[[pair[2]]]$state
  cells <- intersect(names(state_x), names(state_y))
  df <- data.frame(
    method_x = pair[1],
    method_y = pair[2],
    state_x = state_x[cells],
    state_y = state_y[cells],
    stringsAsFactors = FALSE
  )
  df %>%
    count(method_x, method_y, state_x, state_y, name = "n") %>%
    group_by(method_x, method_y, state_x) %>%
    mutate(row_pct = 100 * n / sum(n)) %>%
    ungroup()
}))
write.csv(state_counts_long, file.path(out_dir, "Auto_state_pairwise_concordance_counts.csv"), row.names = FALSE)

make_state_heatmap <- function(method_x, method_y) {
  state_x <- state_results[[method_x]]$state
  state_y <- state_results[[method_y]]$state
  cells <- intersect(names(state_x), names(state_y))
  state_x <- state_x[cells]
  state_y <- state_y[cells]

  levels_pair <- c(state_level_order, sort(setdiff(unique(c(state_x, state_y)), state_level_order)))
  tab <- table(factor(state_x, levels = levels_pair), factor(state_y, levels = levels_pair))
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  tab_pct <- sweep(tab, 1, rowSums(tab), FUN = "/") * 100

  stat_row <- state_pair_stats %>% filter(method_x == !!method_x, method_y == !!method_y)
  title_text <- paste0(method_labels[[method_x]], " vs ", method_labels[[method_y]])

  col_fun <- colorRamp2(c(0, 50, 100), c("#FFFFFF", "#E31A1C", "#7F0000"))
  Heatmap(
    tab_pct,
    name = "% of row state",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = col_fun,
    rect_gp = gpar(col = "white", lwd = 1),
    row_names_side = "left",
    row_title = method_labels[[method_x]],
    column_title = title_text,
    heatmap_legend_param = list(title = "%"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 9),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (tab[i, j] == 0) return(NULL)
      lab <- sprintf("%.1f%%\nN=%d", tab_pct[i, j], tab[i, j])
      grid.text(
        lab,
        x,
        y,
        gp = gpar(fontsize = 8.5, col = ifelse(tab_pct[i, j] > 50, "white", "black"))
      )
    }
  )
}

pdf(file.path(out_dir, "Auto_state_pairwise_concordance_heatmaps.pdf"), width = 10, height = 7.5, useDingbats = FALSE)
for (i in seq_along(method_pairs)) {
  pair <- method_pairs[[i]]
  draw(make_state_heatmap(pair[1], pair[2]), heatmap_legend_side = "right")
}
dev.off()

method_config <- data.frame(
  parameter = c(
    "metaprogram_object",
    "cumulative_weight_threshold",
    "max_rank",
    "ncores",
    "plot_max_cells",
    "state_threshold",
    "hybrid_gap",
    "random_seed",
    "excluded_mps",
    "retained_mps"
  ),
  value = c(
    mp_path,
    as.character(cum_weight_threshold),
    as.character(max_rank),
    as.character(ncores),
    as.character(plot_max_cells),
    as.character(state_threshold),
    as.character(hybrid_gap),
    as.character(random_seed),
    paste(exclude_mps, collapse = ","),
    paste(retained_mps, collapse = ",")
  ),
  stringsAsFactors = FALSE
)
write.csv(method_config, file.path(out_dir, "Auto_run_configuration.csv"), row.names = FALSE)

message("=== SUCCESS ===")
message("Outputs written to: ", file.path(getwd(), out_dir))
