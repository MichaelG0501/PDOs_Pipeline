####################
# Auto_parse_pdo_mp_scoring_method_trend_check.R
# Scores PDO-derived nMP13 metaprograms on six Parse trajectory/recovery samples
# using three activity methods, then compares MP activity and Approach B/noreg
# state-abundance trends against the expected T1 -> T4 depletion/recovery pattern.
#
# Supports cached execution: if intermediate RDS and CSV files exist, skips the
# heavy raw data loading and scoring, proceeding directly to plotting.
# Force recompute with environment variable: PDO_FORCE_REBUILD=1
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_parse_pdo_mp_scoring_method_trend_check"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Check for cache
score_rds_path <- file.path(out_dir, "Auto_parse_pdo_mp_scores_all_methods.rds")
assignments_csv_path <- file.path(out_dir, "Auto_parse_pdo_state_assignments_all_methods.csv")
force_rebuild <- Sys.getenv("PDO_FORCE_REBUILD", unset = "0") == "1"
use_cache <- file.exists(score_rds_path) && file.exists(assignments_csv_path) && !force_rebuild

parse_samples <- c("T0", "T1", "T2", "T4", "R4", "eR4")
sample_cols <- c(
  "T0" = "#0072B2",
  "T1" = "#E69F00",
  "T2" = "#009E73",
  "T4" = "#D55E00",
  "R4" = "#CC79A7",
  "eR4" = "#56B4E9"
)

method_labels <- c(
  full_ucell = "Full gene-list UCell",
  cum_weight_ucell = "Cumulative 70% weight UCell",
  weighted_rank = "Weighted rank score"
)

mp_desc_full <- c(
  "MP6"  = "MP6_G2M Cell Cycle",
  "MP7"  = "MP7_DNA repair",
  "MP5"  = "MP5_MYC-related Proliferation",
  "MP1"  = "MP1_G2M checkpoint",
  "MP3"  = "MP3_G1S Cell Cycle",
  "MP8"  = "MP8_Columnar Progenitor",
  "MP10" = "MP10_Inflammatory Stress Epi.",
  "MP9"  = "MP9_ECM Remodeling Epi.",
  "MP4"  = "MP4_Intestinal Metaplasia"
)

mp_expected_class <- c(
  "MP5" = "lineage",
  "MP4" = "lineage",
  "MP8" = "lineage",
  "MP10" = "stress",
  "MP9" = "stress",
  "MP6" = "cell_cycle",
  "MP7" = "cell_cycle",
  "MP1" = "cell_cycle",
  "MP3" = "cell_cycle"
)

state_expected_class <- c(
  "Classic Proliferative" = "lineage",
  "Basal to Intest. Meta" = "lineage",
  "SMG-like Metaplasia" = "lineage",
  "Stress-adaptive" = "stress",
  "Unresolved" = "unresolved",
  "Hybrid" = "hybrid"
)

expected_class_cols <- c(
  lineage = "#4DAF4A",
  stress = "#984EA3",
  cell_cycle = "#999999",
  unresolved = "grey70",
  hybrid = "black"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)
state_order <- names(state_cols)

score_expected_trend <- function(sample_values, expected_class) {
  names(sample_values) <- as.character(names(sample_values))
  needed <- c("T1", "T4", "R4", "eR4")
  if (!all(needed %in% names(sample_values)) || any(!is.finite(sample_values[needed]))) {
    return(data.frame(
      expected_direction_score = NA_integer_,
      t4_minus_t1 = NA_real_,
      r4_minus_t4 = NA_real_,
      er4_minus_t4 = NA_real_,
      follows_t1_t4 = NA,
      follows_r4_recovery = NA,
      follows_er4_recovery = NA
    ))
  }

  t4_minus_t1 <- unname(sample_values["T4"] - sample_values["T1"])
  r4_minus_t4 <- unname(sample_values["R4"] - sample_values["T4"])
  er4_minus_t4 <- unname(sample_values["eR4"] - sample_values["T4"])

  if (expected_class == "lineage") {
    follows <- c(t4_minus_t1 < 0, r4_minus_t4 > 0, er4_minus_t4 > 0)
  } else if (expected_class == "stress") {
    follows <- c(t4_minus_t1 > 0, r4_minus_t4 < 0, er4_minus_t4 < 0)
  } else {
    follows <- c(NA, NA, NA)
  }

  data.frame(
    expected_direction_score = if (all(is.na(follows))) NA_integer_ else sum(follows, na.rm = TRUE),
    t4_minus_t1 = t4_minus_t1,
    r4_minus_t4 = r4_minus_t4,
    er4_minus_t4 = er4_minus_t4,
    follows_t1_t4 = follows[1],
    follows_r4_recovery = follows[2],
    follows_er4_recovery = follows[3]
  )
}

if (use_cache) {
  message("Found cached intermediate scoring and assignment files. Skipping heavy single-cell processing and UCell scoring.")
  score_list <- readRDS(score_rds_path)
  assignments <- read.csv(assignments_csv_path, check.names = FALSE)
  
  cell_meta <- data.frame(
    cell = rownames(score_list[[1]]),
    sample = sub("__.*$", "", rownames(score_list[[1]])),
    stringsAsFactors = FALSE
  )
  cell_meta$sample <- factor(cell_meta$sample, levels = parse_samples)
  
  pdo_mp_order <- colnames(score_list[[1]])
  
} else {
  message("Cached intermediates not found (or rebuild forced). Running heavy UCell scoring and state assignment...")
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(UCell)
    library(Matrix)
    library(patchwork)
  })

  parse_root <- Sys.getenv(
    "AUTO_PARSE_ROOT_DIR",
    unset = "/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"
  )
  if (!dir.exists(parse_root)) {
    parse_root <- "/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline"
  }
  if (!dir.exists(parse_root)) stop("Parse_Pipeline root not found. Set AUTO_PARSE_ROOT_DIR.")

  get_env_numeric <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || value == "") return(default)
    parsed <- suppressWarnings(as.numeric(value))
    if (is.na(parsed)) default else parsed
  }

  get_env_integer <- function(name, default) {
    as.integer(round(get_env_numeric(name, default)))
  }

  max_rank <- get_env_integer("AUTO_PARSE_PDO_MP_MAX_RANK", 1500)
  ncores <- get_env_integer("AUTO_PARSE_PDO_MP_NCORES", 4)
  cum_weight_threshold <- get_env_numeric("AUTO_PARSE_PDO_MP_CUM_WEIGHT_THRESHOLD", 0.70)
  state_threshold <- get_env_numeric("AUTO_PARSE_PDO_MP_STATE_THRESHOLD", 0.5)
  hybrid_gap <- get_env_numeric("AUTO_PARSE_PDO_MP_HYBRID_GAP", 0.3)

  if (cum_weight_threshold <= 0 || cum_weight_threshold > 1) {
    stop("AUTO_PARSE_PDO_MP_CUM_WEIGHT_THRESHOLD must be in (0, 1].")
  }
  if (ncores < 1) ncores <- 1

  message("Parse root: ", parse_root)
  
  # Helpers specific to heavy block
  get_counts <- function(obj) {
    suppressWarnings({
      tryCatch(
        GetAssayData(obj, assay = "RNA", layer = "counts"),
        error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
      )
    })
  }

  ordered_block <- function(mps, tree_order_names) {
    out <- intersect(tree_order_names, mps)
    c(out, setdiff(mps, out))
  }

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

  select_cumulative_weight_genes <- function(weights, threshold) {
    weights <- weights[is.finite(weights) & weights > 0]
    if (length(weights) == 0) return(character(0))
    weights <- sort(weights, decreasing = TRUE)
    keep_n <- which(cumsum(weights / sum(weights)) >= threshold)[1]
    names(weights)[seq_len(keep_n)]
  }

  clean_gene_list <- function(gene_list, available_genes, min_genes = 3) {
    cleaned <- lapply(gene_list, function(genes) {
      genes <- unique(as.character(genes))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      intersect(genes, available_genes)
    })
    too_short <- names(cleaned)[lengths(cleaned) < min_genes]
    if (length(too_short) > 0) {
      stop("Too few genes present in Parse matrix for: ", paste(too_short, collapse = ", "))
    }
    cleaned
  }

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
    as.matrix(scores[, names(features), drop = FALSE])
  }

  weighted_rank_score_from_ranks <- function(rankings, weight_list, max_rank) {
    all_weighted_genes <- unique(unlist(lapply(weight_list, names), use.names = FALSE))
    genes_present <- intersect(all_weighted_genes, rownames(rankings))
    if (length(genes_present) == 0) stop("No weighted PDO MP genes found in Parse rank matrix.")

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
      if (length(weights) < 3) stop("Too few weighted genes present for ", mp_names[j])
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
    as.matrix(scores[, mp_names, drop = FALSE])
  }

  mp_desc <- sub("^MP[0-9]+_", "", mp_desc_full)
  names(mp_desc) <- names(mp_desc_full)

  # Load Parse counts
  message("Loading six Parse samples")
  sample_files <- setNames(
    file.path(parse_root, "parse_outs", "by_samples", parse_samples, paste0("Auto_", parse_samples, "_final.rds")),
    parse_samples
  )
  missing_files <- sample_files[!file.exists(sample_files)]
  if (length(missing_files) > 0) {
    stop("Missing Parse sample RDS files: ", paste(missing_files, collapse = ", "))
  }

  sample_genes <- lapply(sample_files, function(path) rownames(readRDS(path)))
  common_genes <- Reduce(intersect, sample_genes)
  if (length(common_genes) == 0) stop("No common genes across Parse samples.")

  counts_list <- list()
  meta_list <- list()
  for (sample in parse_samples) {
    message("Loading counts for ", sample)
    obj <- readRDS(sample_files[[sample]])
    old_cells <- colnames(obj)
    new_cells <- paste(sample, old_cells, sep = "__")
    counts <- get_counts(obj)[common_genes, , drop = FALSE]
    colnames(counts) <- new_cells
    counts_list[[sample]] <- counts
    meta_list[[sample]] <- data.frame(
      cell = new_cells,
      original_cell = old_cells,
      sample = sample,
      study = "Parse_SUR1090",
      stringsAsFactors = FALSE
    )
    rm(obj, counts)
    gc()
  }

  cell_meta <- bind_rows(meta_list)
  rownames(cell_meta) <- cell_meta$cell
  sample_var <- setNames(cell_meta$sample, cell_meta$cell)
  study_var <- setNames(cell_meta$study, cell_meta$cell)
  counts_all <- do.call(cbind, counts_list)
  rm(counts_list)
  gc()

  write.csv(
    cell_meta %>% count(sample, name = "n_cells"),
    file.path(out_dir, "Auto_parse_sample_cell_counts.csv"),
    row.names = FALSE
  )

  # PDO MP setup
  message("Loading PDO metaprograms")
  pdo_mp_path <- Sys.getenv(
    "AUTO_PARSE_PDO_MP_OBJECT",
    unset = "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds"
  )
  if (!file.exists(pdo_mp_path)) {
    pdo_mp_path <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds"
  }
  if (!file.exists(pdo_mp_path)) stop("PDO metaprogram object not found: ", pdo_mp_path)
  pdo_geneNMF <- readRDS(pdo_mp_path)

  bad_mps <- which(pdo_geneNMF$metaprograms.metrics$silhouette < 0)
  bad_mp_names <- paste0("MP", bad_mps)
  coverage_tbl <- pdo_geneNMF$metaprograms.metrics$sampleCoverage
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
  exclude_mps <- unique(c(bad_mp_names, low_coverage_mps))

  mp_genes_full <- pdo_geneNMF$metaprograms.genes
  mp_genes_full <- mp_genes_full[!names(mp_genes_full) %in% exclude_mps]
  retained_mps <- names(mp_genes_full)

  mp_weights <- pdo_geneNMF$metaprograms.genes.weights[retained_mps]
  mp_weights <- Map(function(weights, genes) {
    weights <- weights[is.finite(weights) & weights > 0]
    weights[intersect(names(weights), genes)]
  }, mp_weights, mp_genes_full)

  mp_genes_cum <- lapply(mp_weights, select_cumulative_weight_genes, threshold = cum_weight_threshold)
  mp_genes_full_present <- clean_gene_list(mp_genes_full, rownames(counts_all))
  mp_genes_cum_present <- clean_gene_list(mp_genes_cum, rownames(counts_all))

  pdo_tree_order <- pdo_geneNMF$programs.clusters[pdo_geneNMF$programs.tree$order]
  pdo_tree_order <- paste0("MP", rev(unique(pdo_tree_order)))
  pdo_tree_order <- pdo_tree_order[pdo_tree_order %in% retained_mps]
  pdo_cc <- intersect(c("MP6", "MP7", "MP1", "MP3"), retained_mps)
  pdo_noncc <- setdiff(retained_mps, pdo_cc)

  state_groups <- list(
    "Classic Proliferative" = c("MP5"),
    "Basal to Intest. Meta" = c("MP4"),
    "SMG-like Metaplasia" = c("MP8"),
    "Stress-adaptive" = c("MP10", "MP9")
  )
  state_groups <- lapply(state_groups, function(mps) intersect(mps, pdo_noncc))
  state_groups <- state_groups[lengths(state_groups) > 0]

  pdo_mp_order <- ordered_block(pdo_cc, pdo_tree_order)
  for (state in names(state_groups)) {
    pdo_mp_order <- c(pdo_mp_order, ordered_block(intersect(state_groups[[state]], retained_mps), pdo_tree_order))
  }
  pdo_mp_order <- unique(c(pdo_mp_order, setdiff(retained_mps, pdo_mp_order)))

  weight_summary <- bind_rows(lapply(retained_mps, function(mp) {
    weights <- sort(mp_weights[[mp]], decreasing = TRUE)
    selected <- mp_genes_cum[[mp]]
    data.frame(
      MP = mp,
      full_gene_n = length(mp_genes_full[[mp]]),
      full_present_gene_n = length(mp_genes_full_present[[mp]]),
      cumulative_weight_threshold = cum_weight_threshold,
      cumulative_gene_n = length(selected),
      cumulative_present_gene_n = length(mp_genes_cum_present[[mp]]),
      cumulative_weight_sum = sum(weights[selected], na.rm = TRUE) / sum(weights, na.rm = TRUE),
      expected_class = unname(mp_expected_class[mp]),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(weight_summary, file.path(out_dir, "Auto_parse_pdo_mp_weight_gene_counts.csv"), row.names = FALSE)

  # Score three methods
  message("Precomputing UCell ranks for Parse counts")
  rankings <- StoreRankings_UCell(
    counts_all,
    maxRank = max_rank,
    ncores = ncores,
    ties.method = "average",
    force.gc = TRUE
  )

  message("Scoring full gene-list UCell")
  scores_full_ucell <- score_ucell_from_ranks(rankings, mp_genes_full_present, max_rank)

  message("Scoring cumulative-weight UCell")
  scores_cum_ucell <- score_ucell_from_ranks(rankings, mp_genes_cum_present, max_rank)

  message("Scoring weighted-rank activity")
  scores_weighted_rank <- weighted_rank_score_from_ranks(rankings, mp_weights, max_rank)

  score_list <- list(
    full_ucell = scores_full_ucell,
    cum_weight_ucell = scores_cum_ucell,
    weighted_rank = scores_weighted_rank
  )
  score_list <- lapply(score_list, function(mat) mat[, pdo_mp_order, drop = FALSE])
  saveRDS(score_list, file.path(out_dir, "Auto_parse_pdo_mp_scores_all_methods.rds"))
  rm(rankings, counts_all)
  gc()

  # State assignment
  assign_states_from_scores <- function(scores) {
    scores <- scores[cell_meta$cell, retained_mps, drop = FALSE]
    cc_present <- intersect(pdo_cc, colnames(scores))
    noncc_present <- intersect(pdo_noncc, colnames(scores))
    mp_adj_noncc <- z_normalise(scores[, noncc_present, drop = FALSE], sample_var, study_var)
    mp_adj_cc <- z_normalise(scores[, cc_present, drop = FALSE], sample_var, study_var)
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
    state_vec <- colnames(group_max)[best_group_idx]
    state_vec[best_group_val < state_threshold] <- "Unresolved"

    sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
    gap <- sorted_groups[, 1] - sorted_groups[, 2]
    state_vec[(gap < hybrid_gap) & (state_vec != "Unresolved")] <- "Hybrid"
    names(state_vec) <- rownames(group_max)

    top_noncc <- colnames(mp_adj_noncc)[max.col(mp_adj_noncc, ties.method = "first")]
    names(top_noncc) <- rownames(mp_adj_noncc)
    top_all <- colnames(mp_adj_all)[max.col(mp_adj_all, ties.method = "first")]
    names(top_all) <- rownames(mp_adj_all)

    list(
      state = state_vec,
      top_noncc = top_noncc,
      top_all = top_all,
      mp_adj_noncc = mp_adj_noncc,
      mp_adj_all = mp_adj_all,
      group_max = group_max
    )
  }

  message("Assigning Approach B/noreg states")
  state_results <- lapply(score_list, assign_states_from_scores)
  saveRDS(state_results, file.path(out_dir, "Auto_parse_pdo_state_results_all_methods.rds"))

  assignments <- bind_rows(lapply(names(state_results), function(method_name) {
    result <- state_results[[method_name]]
    data.frame(
      method = method_name,
      method_label = method_labels[[method_name]],
      cell = names(result$state),
      original_cell = cell_meta[names(result$state), "original_cell"],
      sample = cell_meta[names(result$state), "sample"],
      pdo_state = unname(result$state),
      top_noncc_mp = unname(result$top_noncc),
      top_noncc_label = unname(mp_desc_full[result$top_noncc]),
      top_all_mp = unname(result$top_all),
      top_all_label = unname(mp_desc_full[result$top_all]),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(assignments, file.path(out_dir, "Auto_parse_pdo_state_assignments_all_methods.csv"), row.names = FALSE)

  config <- data.frame(
    parameter = c(
      "parse_root",
      "pdo_metaprogram_object",
      "samples",
      "max_rank",
      "ncores",
      "cumulative_weight_threshold",
      "state_threshold",
      "hybrid_gap",
      "excluded_mps",
      "retained_mps"
    ),
    value = c(
      parse_root,
      pdo_mp_path,
      paste(parse_samples, collapse = ","),
      as.character(max_rank),
      as.character(ncores),
      as.character(cum_weight_threshold),
      as.character(state_threshold),
      as.character(hybrid_gap),
      paste(exclude_mps, collapse = ","),
      paste(retained_mps, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
  write.csv(config, file.path(out_dir, "Auto_parse_pdo_scoring_method_trend_config.csv"), row.names = FALSE)
}

# Shared Downstream Summaries & Plots
message("Computing summaries and generating figures")

mp_axis_levels <- unname(mp_desc_full[pdo_mp_order])
names(mp_axis_levels) <- pdo_mp_order

score_long <- bind_rows(lapply(names(score_list), function(method_name) {
  as.data.frame(score_list[[method_name]][, pdo_mp_order, drop = FALSE]) %>%
    tibble::rownames_to_column("cell") %>%
    left_join(cell_meta[, c("cell", "sample")], by = "cell") %>%
    pivot_longer(cols = all_of(pdo_mp_order), names_to = "MP", values_to = "score") %>%
    mutate(
      method = method_name,
      method_label = method_labels[[method_name]],
      MP_label = unname(mp_desc_full[MP]),
      expected_class = unname(mp_expected_class[MP])
    )
}))

mp_sample_summary <- score_long %>%
  group_by(method, method_label, MP, MP_label, expected_class, sample) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    q1 = quantile(score, 0.25, na.rm = TRUE),
    q3 = quantile(score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sample = factor(sample, levels = parse_samples))
write.csv(mp_sample_summary, file.path(out_dir, "Auto_parse_pdo_mp_activity_sample_summary.csv"), row.names = FALSE)

state_abundance_summary <- assignments %>%
  count(method, method_label, sample, pdo_state, name = "n") %>%
  right_join(
    expand_grid(method = names(method_labels), sample = parse_samples, pdo_state = state_order),
    by = c("method", "sample", "pdo_state")
  ) %>%
  mutate(n = replace_na(n, 0L), method_label = unname(method_labels[method])) %>%
  group_by(method, method_label, sample) %>%
  mutate(total_cells = sum(n), pct = 100 * n / pmax(total_cells, 1)) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = parse_samples),
    expected_class = unname(state_expected_class[pdo_state])
  )
write.csv(state_abundance_summary, file.path(out_dir, "Auto_parse_pdo_state_abundance_sample_summary.csv"), row.names = FALSE)

mp_trend_eval <- mp_sample_summary %>%
  filter(expected_class %in% c("lineage", "stress")) %>%
  select(method, method_label, MP, MP_label, expected_class, sample, mean_score) %>%
  pivot_wider(names_from = sample, values_from = mean_score) %>%
  rowwise() %>%
  mutate(
    trend = list(score_expected_trend(c(T1 = T1, T4 = T4, R4 = R4, eR4 = eR4), expected_class))
  ) %>%
  unnest(trend) %>%
  ungroup()
write.csv(mp_trend_eval, file.path(out_dir, "Auto_parse_pdo_mp_expected_trend_eval.csv"), row.names = FALSE)

state_trend_eval <- state_abundance_summary %>%
  filter(expected_class %in% c("lineage", "stress")) %>%
  select(method, method_label, pdo_state, expected_class, sample, pct) %>%
  pivot_wider(names_from = sample, values_from = pct) %>%
  rowwise() %>%
  mutate(
    trend = list(score_expected_trend(c(T1 = T1, T4 = T4, R4 = R4, eR4 = eR4), expected_class))
  ) %>%
  unnest(trend) %>%
  ungroup()
write.csv(state_trend_eval, file.path(out_dir, "Auto_parse_pdo_state_expected_trend_eval.csv"), row.names = FALSE)

method_trend_summary <- bind_rows(
  mp_trend_eval %>% mutate(feature_type = "MP", feature = MP_label),
  state_trend_eval %>% mutate(feature_type = "State", feature = pdo_state)
) %>%
  group_by(method, method_label, feature_type, expected_class) %>%
  summarise(
    n_features = n(),
    mean_expected_direction_score = mean(expected_direction_score, na.rm = TRUE),
    full_match_n = sum(expected_direction_score == 3, na.rm = TRUE),
    partial_or_full_match_n = sum(expected_direction_score >= 2, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(method_trend_summary, file.path(out_dir, "Auto_parse_pdo_method_expected_trend_summary.csv"), row.names = FALSE)

mp_dot_df <- mp_sample_summary %>%
  group_by(method, MP) %>%
  mutate(
    mean_delta_t1 = mean_score - mean_score[sample == "T1"][1],
    mean_z_by_mp_method = as.numeric(scale(mean_score))
  ) %>%
  ungroup() %>%
  mutate(
    method_label = factor(method_label, levels = method_labels[names(method_labels)]),
    sample = factor(sample, levels = parse_samples),
    MP_label = factor(MP_label, levels = mp_axis_levels[pdo_mp_order]),
    expected_class = factor(expected_class, levels = c("cell_cycle", "lineage", "stress"))
  )

p_mp_dot <- ggplot(mp_dot_df, aes(x = sample, y = MP_label)) +
  geom_point(aes(size = mean_score, fill = mean_delta_t1), shape = 21, color = "black", stroke = 0.25) +
  facet_grid(expected_class ~ method_label, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Mean score\nminus T1") +
  scale_size_continuous(range = c(2.2, 7.5), name = "Mean score") +
  labs(
    title = "PDO-derived MP activity on Parse samples",
    subtitle = "Dot fill is mean activity change relative to T1; T1 -> T4 depletion and R4/eR4 recovery are expected for lineage MPs, opposite for stress MPs.",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 10.5, colour = "grey35"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )
ggsave(file.path(out_dir, "Auto_parse_pdo_mp_activity_method_dotplot.pdf"), p_mp_dot, width = 16, height = 9, useDingbats = FALSE)
ggsave(file.path(out_dir, "Auto_parse_pdo_mp_activity_method_dotplot.png"), p_mp_dot, width = 16, height = 9, dpi = 300)

p_mp_line <- mp_dot_df %>%
  filter(expected_class %in% c("lineage", "stress")) %>%
  ggplot(aes(x = sample, y = mean_score, group = MP_label, color = expected_class)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(aes(fill = sample), shape = 21, color = "black", size = 2.8, stroke = 0.25) +
  facet_grid(method_label ~ MP_label, scales = "free_y") +
  scale_color_manual(values = expected_class_cols, guide = "none") +
  scale_fill_manual(values = sample_cols, name = "Sample") +
  labs(
    title = "Lineage/stress PDO MP mean activity trends on Parse samples",
    x = NULL,
    y = "Mean activity score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "top"
  )
ggsave(file.path(out_dir, "Auto_parse_pdo_mp_activity_line_trends.pdf"), p_mp_line, width = 18, height = 9, useDingbats = FALSE)

state_dot_df <- state_abundance_summary %>%
  group_by(method, pdo_state) %>%
  mutate(
    pct_delta_t1 = pct - pct[sample == "T1"][1],
    method_label = unname(method_labels[method])
  ) %>%
  ungroup() %>%
  mutate(
    method_label = factor(method_label, levels = method_labels[names(method_labels)]),
    sample = factor(sample, levels = parse_samples),
    pdo_state = factor(pdo_state, levels = state_order),
    expected_class = factor(expected_class, levels = c("lineage", "stress", "unresolved", "hybrid"))
  )

p_state_dot <- ggplot(state_dot_df, aes(x = sample, y = pdo_state)) +
  geom_point(aes(size = pct, fill = pct_delta_t1), shape = 21, color = "black", stroke = 0.25) +
  facet_grid(expected_class ~ method_label, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Percent\nminus T1") +
  scale_size_continuous(range = c(2.2, 8.5), name = "% cells") +
  labs(
    title = "PDO Approach B/noreg state abundance on Parse samples",
    subtitle = "Dot fill is abundance change relative to T1; lineage states are expected to fall at T4 then recover, stress state is expected to do the opposite.",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 10.5, colour = "grey35"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )
ggsave(file.path(out_dir, "Auto_parse_pdo_state_abundance_method_dotplot.pdf"), p_state_dot, width = 16, height = 8, useDingbats = FALSE)
ggsave(file.path(out_dir, "Auto_parse_pdo_state_abundance_method_dotplot.png"), p_state_dot, width = 16, height = 8, dpi = 300)

p_state_line <- state_dot_df %>%
  filter(expected_class %in% c("lineage", "stress")) %>%
  ggplot(aes(x = sample, y = pct, group = pdo_state, color = pdo_state)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(aes(fill = sample), shape = 21, color = "black", size = 2.8, stroke = 0.25) +
  facet_grid(method_label ~ pdo_state, scales = "free_y") +
  scale_color_manual(values = state_cols, guide = "none") +
  scale_fill_manual(values = sample_cols, name = "Sample") +
  labs(
    title = "Lineage/stress PDO state abundance trends on Parse samples",
    x = NULL,
    y = "% cells"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "top"
  )
ggsave(file.path(out_dir, "Auto_parse_pdo_state_abundance_line_trends.pdf"), p_state_line, width = 16, height = 8, useDingbats = FALSE)

p_summary <- method_trend_summary %>%
  mutate(
    method_label = factor(method_label, levels = method_labels[names(method_labels)]),
    expected_class = factor(expected_class, levels = c("lineage", "stress"))
  ) %>%
  ggplot(aes(x = method_label, y = mean_expected_direction_score, fill = expected_class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.25) +
  facet_wrap(~ feature_type, nrow = 1) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "grey35", linewidth = 0.35) +
  scale_fill_manual(values = expected_class_cols, name = "Expected class") +
  scale_y_continuous(limits = c(0, 3), breaks = 0:3) +
  labs(
    title = "Expected trend agreement by scoring method",
    subtitle = "Score counts T1->T4 direction plus R4/eR4 recovery directions; maximum = 3.",
    x = NULL,
    y = "Mean expected-direction score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10.5, colour = "grey35"),
    axis.text.x = element_text(angle = 25, hjust = 1, colour = "black"),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(out_dir, "Auto_parse_pdo_expected_trend_method_summary.pdf"), p_summary, width = 12, height = 6, useDingbats = FALSE)
ggsave(file.path(out_dir, "Auto_parse_pdo_expected_trend_method_summary.png"), p_summary, width = 12, height = 6, dpi = 300)

message("Done. Outputs written to: ", file.path(getwd(), out_dir))
