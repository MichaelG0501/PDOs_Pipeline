####################
# Analysis registry:
#   Status: active terminal comparison; supersedes the analytical approach in
#           analysis/cell_states/PDO_state_concordance.R
#   Script: analysis/cell_states/Auto_PDO_scRef_centred_state_concordance.R
#   Methodology: analysis/methodology/cell_states/Auto_PDO_scRef_centred_state_concordance_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     - PDOs_outs/PDOs_merged.rds
#     - PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     - PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     - PDOs_outs/centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv
#     - scRef live ref_outs/EAC_Ref_epi.rds
#     - scRef live centred state_definition/intermediate/centred_refined_noreg_states.rds
#     - scRef live centred mp_refinement/intermediate/merged_refined_mp_genes.rds
#     - scRef live centred mp_refinement/tables/centred_refined_mp_state_grouping.csv
#   Outputs:
#     - PDOs_outs/Auto_PDO_scRef_centred_state_concordance/figures/Auto_PDO_scRef_centred_state_concordance.pdf
#     - PDOs_outs/Auto_PDO_scRef_centred_state_concordance/tables/*.csv
#     - PDOs_outs/Auto_PDO_scRef_centred_state_concordance/intermediate/Auto_state_concordance_results.rds
#     - PDOs_outs/Auto_PDO_scRef_centred_state_concordance/Auto_scRef_current_MP_states_on_PDO.rds
#     - PDOs_outs/Auto_PDO_scRef_centred_state_concordance/logs/Auto_PDO_scRef_centred_state_concordance_run_summary.txt
#   Downstream use: none; terminal concordance evidence and auditable inputs
####################

suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
  library("edgeR")
  library("limma")
  library("UCell")
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("ggrepel")
  library("patchwork")
})

source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_helpers.R")

####################
# Paths, state definitions, and expected biological correspondence
####################
script_path <- "analysis/cell_states/Auto_PDO_scRef_centred_state_concordance.R"
out_dir <- file.path(PDO_LIVE_OUTS, "Auto_PDO_scRef_centred_state_concordance")
out_paths <- pdo_ensure_output_tiers(out_dir)
cache_policy <- pdo_cache_policy()

pdo_seurat_file <- file.path(PDO_LIVE_OUTS, "PDOs_merged.rds")
pdo_state_file <- file.path(
  PDO_LIVE_OUTS, "centred_mp_refinement", "centred_refined_noreg_states.rds"
)
pdo_gene_file <- file.path(
  PDO_LIVE_OUTS, "centred_mp_refinement", "merged_refined_mp_genes.rds"
)
pdo_group_file <- file.path(
  PDO_LIVE_OUTS, "centred_mp_refinement", "tables",
  "centred_refined_mp_state_grouping.csv"
)

sc_ref_outs <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"
sc_seurat_file <- file.path(sc_ref_outs, "EAC_Ref_epi.rds")
sc_state_file <- file.path(
  sc_ref_outs, "Metaprogrammes_Results", "centred", "state_definition",
  "intermediate", "centred_refined_noreg_states.rds"
)
sc_gene_file <- file.path(
  sc_ref_outs, "Metaprogrammes_Results", "centred", "mp_refinement",
  "intermediate", "merged_refined_mp_genes.rds"
)
sc_group_file <- file.path(
  sc_ref_outs, "Metaprogrammes_Results", "centred", "mp_refinement",
  "tables", "centred_refined_mp_state_grouping.csv"
)

input_files <- c(
  pdo_seurat_file, pdo_state_file, pdo_gene_file, pdo_group_file,
  sc_seurat_file, sc_state_file, sc_gene_file, sc_group_file
)
pdo_require_files(input_files)

pdo_state_order <- c(
  "Classic proliferation",
  "Columnar-to-intestinal",
  "Glandular differentiation",
  "Stress-adaptive",
  "ECM-remodelling",
  "Motile-cilia differentiation"
)
sc_state_order <- c(
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "Cancer-cell immune mimicry"
)
expected_state_map <- c(
  "Classic proliferation" = "Classic proliferation",
  "Columnar-to-intestinal" = "Squamous-to-intestinal",
  "Glandular differentiation" = "Glandular-to-intestinal",
  "Stress-adaptive" = "Stress-adaptive"
)
shared_state_order <- names(expected_state_map)

state_colors <- c(
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Motile-cilia differentiation" = "#F781BF",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Cancer-cell immune mimicry" = "#377EB8"
)

min_cells_pseudobulk <- 20L
top_marker_n <- 200L
effect_gene_n <- 500L
cache_file <- file.path(
  out_paths[["intermediate"]], "Auto_state_concordance_results.rds"
)
####################

####################
# Input helpers: current step-05 MP genes are grouped by state, never by MP ID
####################
normalise_state_vector <- function(state_vec, dataset_name) {
  state_names <- names(state_vec)
  state_vec <- as.character(state_vec)
  names(state_vec) <- state_names
  if (is.null(names(state_vec)) || any(names(state_vec) == "")) {
    stop(dataset_name, " state vector lacks cell-barcode names.")
  }
  state_vec
}

load_state_mp_signatures <- function(gene_file, grouping_file, allowed_states, prefix) {
  mp_genes <- readRDS(gene_file)
  grouping <- read.csv(grouping_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("state", "mp")
  if (!all(required_cols %in% colnames(grouping))) {
    stop("MP grouping table lacks required columns: ", grouping_file)
  }
  grouping <- grouping[grouping$state %in% allowed_states, , drop = FALSE]
  missing_mps <- setdiff(grouping$mp, names(mp_genes))
  if (length(missing_mps) > 0) {
    stop("Grouped MP(s) absent from current gene lists: ",
         paste(missing_mps, collapse = ", "))
  }
  signatures <- lapply(allowed_states, function(state_name) {
    state_mps <- grouping$mp[grouping$state == state_name]
    unique(unlist(mp_genes[state_mps], use.names = FALSE))
  })
  names(signatures) <- paste0(prefix, "_S", seq_along(allowed_states))
  mapping <- data.frame(
    signature_id = names(signatures),
    source_state = allowed_states,
    source_dataset = prefix,
    n_genes = lengths(signatures),
    stringsAsFactors = FALSE
  )
  if (any(lengths(signatures) == 0)) {
    stop("At least one current MP-derived state signature is empty: ", prefix)
  }
  list(signatures = signatures, mapping = mapping)
}

pdo_mp_signature_info <- load_state_mp_signatures(
  pdo_gene_file, pdo_group_file, pdo_state_order, "PDO"
)
sc_mp_signature_info <- load_state_mp_signatures(
  sc_gene_file, sc_group_file, sc_state_order, "SC"
)
all_mp_signatures <- c(
  pdo_mp_signature_info$signatures,
  sc_mp_signature_info$signatures
)
signature_map <- bind_rows(
  pdo_mp_signature_info$mapping,
  sc_mp_signature_info$mapping
)

sc_current_mp_genes <- readRDS(sc_gene_file)
sc_current_mp_grouping <- read.csv(
  sc_group_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  filter(state %in% sc_state_order)
sc_projection_ids <- paste0("SCMP_", seq_len(nrow(sc_current_mp_grouping)))
sc_projection_signatures <- sc_current_mp_genes[sc_current_mp_grouping$mp]
names(sc_projection_signatures) <- sc_projection_ids
sc_projection_map <- data.frame(
  signature_id = sc_projection_ids,
  mp = sc_current_mp_grouping$mp,
  state = sc_current_mp_grouping$state,
  stringsAsFactors = FALSE
)
####################

####################
# Construct sample-by-state pseudobulks and cross-dataset MP signature scores
####################
z_normalise_projected_scores <- function(mat, sample_var, study_var) {
  score_df <- as.data.frame(mat)
  score_df$.cell <- rownames(mat)
  score_df$.sample <- sample_var[rownames(mat)]
  score_df$.study <- study_var[rownames(mat)]

  study_sd <- score_df %>%
    group_by(.study) %>%
    summarise(
      across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  score_centered <- score_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  score_adj <- as.matrix(score_centered[, colnames(mat), drop = FALSE])
  rownames(score_adj) <- score_centered$.cell
  for (feature in colnames(score_adj)) {
    score_adj[, feature] <- score_adj[, feature] /
      study_sd[score_centered$.study, feature]
  }
  score_adj[!is.finite(score_adj)] <- 0
  score_adj
}

define_projected_states <- function(
    score_mat,
    sample_vec,
    study_vec,
    projection_map,
    target_order,
    threshold = 0.5,
    hybrid_gap = 0.3) {
  score_adj <- z_normalise_projected_scores(score_mat, sample_vec, study_vec)
  group_max <- sapply(target_order, function(state_name) {
    features <- projection_map$signature_id[projection_map$state == state_name]
    features <- intersect(features, colnames(score_adj))
    if (length(features) == 0) return(rep(NA_real_, nrow(score_adj)))
    if (length(features) == 1) return(as.numeric(score_adj[, features]))
    apply(score_adj[, features, drop = FALSE], 1, max)
  })
  group_max <- as.matrix(group_max)
  rownames(group_max) <- rownames(score_adj)
  group_max[!is.finite(group_max)] <- NA_real_

  best_idx <- max.col(group_max, ties.method = "first")
  best_value <- apply(group_max, 1, max, na.rm = TRUE)
  projected <- target_order[best_idx]
  projected[!is.finite(best_value) | best_value < threshold] <- "Unresolved"

  sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
  gap <- sorted_groups[, 1] - sorted_groups[, 2]
  projected[(gap < hybrid_gap) & projected != "Unresolved"] <- "Hybrid"
  names(projected) <- rownames(group_max)
  list(states = projected, group_max = group_max)
}

prepare_dataset <- function(
    seurat_file,
    state_file,
    state_order,
    dataset_name,
    excluded_samples = character(0),
    projection_signatures = NULL,
    projection_map = NULL,
    projection_state_order = NULL) {
  cat("Loading ", dataset_name, " Seurat object...\n", sep = "")
  seurat_obj <- readRDS(seurat_file)
  state_vec <- normalise_state_vector(readRDS(state_file), dataset_name)

  common_cells <- intersect(Cells(seurat_obj), names(state_vec))
  if (length(common_cells) == 0) {
    stop("No overlapping cells between ", dataset_name, " object and state vector.")
  }
  sample_vec <- as.character(seurat_obj$orig.ident)
  names(sample_vec) <- Cells(seurat_obj)
  study_vec <- if ("Batch" %in% colnames(seurat_obj@meta.data)) {
    as.character(seurat_obj$Batch)
  } else if ("study" %in% colnames(seurat_obj@meta.data)) {
    as.character(seurat_obj$study)
  } else {
    sample_vec
  }
  names(study_vec) <- Cells(seurat_obj)
  keep_cells <- common_cells[
    state_vec[common_cells] %in% state_order &
      !sample_vec[common_cells] %in% excluded_samples
  ]
  if (length(keep_cells) == 0) {
    stop("No target-state cells remain for ", dataset_name, ".")
  }

  state_vec <- state_vec[keep_cells]
  sample_vec <- sample_vec[keep_cells]
  study_vec <- study_vec[keep_cells]
  counts <- pdo_get_assay_matrix(seurat_obj, assay = "RNA", layer = "counts")
  counts <- counts[, keep_cells, drop = FALSE]
  rm(seurat_obj)
  gc()

  group_id <- paste(sample_vec, state_vec, sep = "|||")
  group_counts <- table(group_id)
  valid_groups <- names(group_counts)[group_counts >= min_cells_pseudobulk]
  pb_keep <- group_id %in% valid_groups
  if (sum(pb_keep) == 0) {
    stop("No ", dataset_name, " sample-state groups meet the minimum cell count.")
  }

  group_factor <- factor(group_id[pb_keep], levels = valid_groups)
  group_design <- sparse.model.matrix(~ 0 + group_factor)
  colnames(group_design) <- levels(group_factor)
  pseudobulk_counts <- counts[, pb_keep, drop = FALSE] %*% group_design
  pseudobulk_counts <- as.matrix(pseudobulk_counts)

  split_group <- strsplit(colnames(pseudobulk_counts), "\\|\\|\\|", fixed = FALSE)
  pb_meta <- data.frame(
    pseudobulk = colnames(pseudobulk_counts),
    sample = vapply(split_group, `[`, character(1), 1),
    state = vapply(split_group, function(x) paste(x[-1], collapse = "|||"), character(1)),
    cell_count = as.integer(group_counts[colnames(pseudobulk_counts)]),
    stringsAsFactors = FALSE
  )
  rownames(pb_meta) <- pb_meta$pseudobulk

  scoring_features <- all_mp_signatures
  if (!is.null(projection_signatures)) {
    scoring_features <- c(scoring_features, projection_signatures)
  }
  cat("Scoring current PDO and scRef state MP signatures in ", dataset_name, "...\n", sep = "")
  score_df <- as.data.frame(ScoreSignatures_UCell(
    matrix = counts,
    features = scoring_features,
    name = "",
    ncores = 1
  ), check.names = FALSE)

  projected_result <- NULL
  if (!is.null(projection_signatures)) {
    cat("Defining current scRef-signature states in ", dataset_name, "...\n", sep = "")
    projection_ids <- names(projection_signatures)
    projected_result <- define_projected_states(
      as.matrix(score_df[, projection_ids, drop = FALSE]),
      sample_vec,
      study_vec,
      projection_map,
      projection_state_order
    )
  }

  score_df$state <- state_vec[rownames(score_df)]
  score_df$sample <- sample_vec[rownames(score_df)]

  score_state_means <- score_df %>%
    group_by(state) %>%
    summarise(across(all_of(names(all_mp_signatures)), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  score_sample_state <- score_df %>%
    group_by(sample, state) %>%
    summarise(
      cells = n(),
      across(all_of(names(all_mp_signatures)), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  cell_summary <- data.frame(
    dataset = dataset_name,
    sample = sample_vec,
    state = state_vec,
    stringsAsFactors = FALSE
  ) %>%
    count(dataset, sample, state, name = "cells")

  rm(counts, score_df, group_design)
  gc()
  list(
    pseudobulk_counts = pseudobulk_counts,
    pseudobulk_meta = pb_meta,
    score_state_means = score_state_means,
    score_sample_state = score_sample_state,
    cell_summary = cell_summary,
    native_states = state_vec,
    projected_states = if (is.null(projected_result)) NULL else projected_result$states,
    projected_group_max = if (is.null(projected_result)) NULL else projected_result$group_max
  )
}

fit_state_effects <- function(pb_counts, pb_meta, state_order, dataset_name) {
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  pb_meta$state <- factor(pb_meta$state, levels = state_order)
  pb_meta$sample <- factor(pb_meta$sample)

  y <- DGEList(counts = pb_counts)
  keep_genes <- filterByExpr(y, group = pb_meta$state)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  state_design <- model.matrix(~ 0 + state + sample, data = pb_meta)
  state_cols <- grep("^state", colnames(state_design))
  colnames(state_design)[state_cols] <- state_order
  if (qr(state_design)$rank < ncol(state_design)) {
    non_estimable <- nonEstimable(state_design)
    stop(dataset_name, " pseudobulk design is not full rank: ",
         paste(non_estimable, collapse = ", "))
  }

  contrast_mat <- matrix(
    0, nrow = ncol(state_design), ncol = length(state_order),
    dimnames = list(colnames(state_design), state_order)
  )
  for (state_name in state_order) {
    other_states <- setdiff(state_order, state_name)
    contrast_mat[state_name, state_name] <- 1
    contrast_mat[other_states, state_name] <- -1 / length(other_states)
  }

  v <- voom(y, state_design, plot = FALSE)
  fit <- lmFit(v, state_design)
  fit <- contrasts.fit(fit, contrast_mat)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  tables <- lapply(state_order, function(state_name) {
    topTable(fit, coef = state_name, number = Inf, sort.by = "none") %>%
      tibble::rownames_to_column("gene") %>%
      mutate(dataset = dataset_name, state = state_name)
  })
  names(tables) <- state_order
  list(tables = tables, tested_genes = rownames(v$E), design = state_design)
}
####################

####################
# Build or reuse the computational cache
####################
if (cache_policy$replot_only && file.exists(cache_file) && !cache_policy$force_rebuild) {
  cat("PDO_REPLOT_ONLY=1: loading cached concordance results...\n")
  results <- readRDS(cache_file)
} else {
  pdo_data <- prepare_dataset(
    pdo_seurat_file,
    pdo_state_file,
    pdo_state_order,
    "PDO",
    excluded_samples = PDO_EXCLUDED_SAMPLE,
    projection_signatures = sc_projection_signatures,
    projection_map = sc_projection_map,
    projection_state_order = sc_state_order
  )
  sc_data <- prepare_dataset(
    sc_seurat_file,
    sc_state_file,
    sc_state_order,
    "scRef"
  )

  cat("Fitting sample-aware state-vs-rest pseudobulk models...\n")
  pdo_dge <- fit_state_effects(
    pdo_data$pseudobulk_counts,
    pdo_data$pseudobulk_meta,
    pdo_state_order,
    "PDO"
  )
  sc_dge <- fit_state_effects(
    sc_data$pseudobulk_counts,
    sc_data$pseudobulk_meta,
    sc_state_order,
    "scRef"
  )

  results <- list(
    pdo_data = pdo_data,
    sc_data = sc_data,
    pdo_dge = pdo_dge,
    sc_dge = sc_dge,
    signature_map = signature_map,
    pdo_state_order = pdo_state_order,
    sc_state_order = sc_state_order,
    shared_state_order = shared_state_order,
    parameters = list(
      min_cells_pseudobulk = min_cells_pseudobulk,
      top_marker_n = top_marker_n,
      effect_gene_n = effect_gene_n
    )
  )
  saveRDS(results, cache_file)
}

pdo_data <- results$pdo_data
sc_data <- results$sc_data
pdo_dge <- results$pdo_dge
sc_dge <- results$sc_dge
####################

####################
# MP-signature cross-projection matrices
####################
extract_score_matrix <- function(score_means, signature_ids, target_order) {
  score_means <- score_means %>%
    filter(state %in% target_order) %>%
    mutate(state = factor(state, levels = target_order)) %>%
    arrange(state)
  mat <- as.matrix(score_means[, signature_ids, drop = FALSE])
  rownames(mat) <- as.character(score_means$state)
  t(mat)
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

pdo_sig_ids <- signature_map$signature_id[signature_map$source_dataset == "PDO"]
sc_sig_ids <- signature_map$signature_id[signature_map$source_dataset == "SC"]
pdo_sig_state <- setNames(
  signature_map$source_state[match(pdo_sig_ids, signature_map$signature_id)],
  pdo_sig_ids
)
sc_sig_state <- setNames(
  signature_map$source_state[match(sc_sig_ids, signature_map$signature_id)],
  sc_sig_ids
)

pdo_signatures_in_sc <- extract_score_matrix(
  sc_data$score_state_means, pdo_sig_ids, sc_state_order
)
sc_signatures_in_pdo <- extract_score_matrix(
  pdo_data$score_state_means, sc_sig_ids, pdo_state_order
)
rownames(pdo_signatures_in_sc) <- unname(pdo_sig_state[rownames(pdo_signatures_in_sc)])
rownames(sc_signatures_in_pdo) <- unname(sc_sig_state[rownames(sc_signatures_in_pdo)])
pdo_signatures_in_sc_z <- row_zscore(pdo_signatures_in_sc)
sc_signatures_in_pdo_z <- row_zscore(sc_signatures_in_pdo)
####################

####################
# Direct gene overlap of current step-05 MP-derived state signatures
####################
pdo_state_mp_genes <- setNames(pdo_mp_signature_info$signatures, pdo_state_order)
sc_state_mp_genes <- setNames(sc_mp_signature_info$signatures, sc_state_order)
mp_gene_universe <- unique(c(
  unlist(pdo_state_mp_genes, use.names = FALSE),
  unlist(sc_state_mp_genes, use.names = FALSE)
))
mp_gene_overlap_stats <- list()
idx <- 1L
for (pdo_state in pdo_state_order) {
  for (sc_state in sc_state_order) {
    a <- pdo_state_mp_genes[[pdo_state]]
    b <- sc_state_mp_genes[[sc_state]]
    overlap <- intersect(a, b)
    contingency <- matrix(c(
      length(overlap),
      length(setdiff(a, b)),
      length(setdiff(b, a)),
      length(setdiff(mp_gene_universe, union(a, b)))
    ), nrow = 2)
    ft <- fisher.test(contingency, alternative = "greater")
    mp_gene_overlap_stats[[idx]] <- data.frame(
      pdo_state = pdo_state,
      sc_ref_state = sc_state,
      pdo_signature_genes = length(a),
      sc_ref_signature_genes = length(b),
      overlap_n = length(overlap),
      jaccard = length(overlap) / length(union(a, b)),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      shared_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}
mp_gene_overlap_stats <- bind_rows(mp_gene_overlap_stats) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
mp_gene_jaccard <- xtabs(
  jaccard ~ pdo_state + sc_ref_state,
  data = mp_gene_overlap_stats
)
mp_gene_jaccard <- mp_gene_jaccard[pdo_state_order, sc_state_order, drop = FALSE]
####################

####################
# DGE effect correlation and independently ranked marker overlap
####################
pdo_dge_long <- bind_rows(pdo_dge$tables)
sc_dge_long <- bind_rows(sc_dge$tables)
common_tested_genes <- intersect(pdo_dge$tested_genes, sc_dge$tested_genes)

top_abs_genes <- function(dge_tables, n) {
  unique(unlist(lapply(dge_tables, function(tbl) {
    tbl %>%
      filter(gene %in% common_tested_genes) %>%
      arrange(desc(abs(t))) %>%
      slice_head(n = n) %>%
      pull(gene)
  }), use.names = FALSE))
}
informative_genes <- union(
  top_abs_genes(pdo_dge$tables, effect_gene_n),
  top_abs_genes(sc_dge$tables, effect_gene_n)
)
informative_genes <- intersect(informative_genes, common_tested_genes)

make_effect_matrix <- function(dge_tables, state_order, genes) {
  mat <- vapply(state_order, function(state_name) {
    tbl <- dge_tables[[state_name]]
    setNames(tbl$logFC, tbl$gene)[genes]
  }, numeric(length(genes)))
  rownames(mat) <- genes
  colnames(mat) <- state_order
  mat
}

pdo_effects <- make_effect_matrix(pdo_dge$tables, pdo_state_order, informative_genes)
sc_effects <- make_effect_matrix(sc_dge$tables, sc_state_order, informative_genes)
effect_cor <- cor(pdo_effects, sc_effects, method = "spearman", use = "pairwise.complete.obs")

top_positive_markers <- function(dge_tables, state_order, n) {
  result <- lapply(state_order, function(state_name) {
    dge_tables[[state_name]] %>%
      filter(gene %in% common_tested_genes, logFC > 0) %>%
      arrange(desc(t)) %>%
      slice_head(n = n) %>%
      pull(gene)
  })
  names(result) <- state_order
  result
}

pdo_markers <- top_positive_markers(pdo_dge$tables, pdo_state_order, top_marker_n)
sc_markers <- top_positive_markers(sc_dge$tables, sc_state_order, top_marker_n)
marker_stats <- list()
idx <- 1L
for (pdo_state in pdo_state_order) {
  for (sc_state in sc_state_order) {
    a <- pdo_markers[[pdo_state]]
    b <- sc_markers[[sc_state]]
    overlap <- intersect(a, b)
    contingency <- matrix(c(
      length(overlap),
      length(setdiff(a, b)),
      length(setdiff(b, a)),
      length(setdiff(common_tested_genes, union(a, b)))
    ), nrow = 2)
    ft <- fisher.test(contingency, alternative = "greater")
    marker_stats[[idx]] <- data.frame(
      pdo_state = pdo_state,
      sc_ref_state = sc_state,
      overlap_n = length(overlap),
      jaccard = length(overlap) / length(union(a, b)),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      shared_genes = paste(overlap, collapse = ";"),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}
marker_stats <- bind_rows(marker_stats) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
marker_jaccard <- xtabs(jaccard ~ pdo_state + sc_ref_state, data = marker_stats)
marker_jaccard <- marker_jaccard[pdo_state_order, sc_state_order, drop = FALSE]
####################

####################
# Exact state-label permutation tests and best-match specificity summaries
####################
all_permutations <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

diagonal_specificity_test <- function(mat, metric_name) {
  shared_mat <- mat[shared_state_order, expected_state_map[shared_state_order], drop = FALSE]
  observed_diag <- mean(diag(shared_mat), na.rm = TRUE)
  observed_off <- mean(shared_mat[row(shared_mat) != col(shared_mat)], na.rm = TRUE)
  observed_margin <- observed_diag - observed_off
  perms <- all_permutations(seq_along(shared_state_order))
  perm_margin <- apply(perms, 1, function(p) {
    perm_diag <- mean(shared_mat[cbind(seq_along(p), p)], na.rm = TRUE)
    perm_mask <- matrix(TRUE, nrow = nrow(shared_mat), ncol = ncol(shared_mat))
    perm_mask[cbind(seq_along(p), p)] <- FALSE
    perm_off <- mean(shared_mat[perm_mask], na.rm = TRUE)
    perm_diag - perm_off
  })
  data.frame(
    metric = metric_name,
    matched_mean = observed_diag,
    nonmatched_mean = observed_off,
    diagonal_margin = observed_margin,
    exact_permutation_p = mean(perm_margin >= observed_margin),
    permutations = nrow(perms),
    stringsAsFactors = FALSE
  )
}

combined_mp_projection <- (
  pdo_signatures_in_sc_z[pdo_state_order, sc_state_order, drop = FALSE] +
    t(sc_signatures_in_pdo_z[sc_state_order, pdo_state_order, drop = FALSE])
) / 2

specificity_summary <- bind_rows(
  diagonal_specificity_test(effect_cor, "Pseudobulk DGE-effect correlation"),
  diagonal_specificity_test(marker_jaccard, "Top-marker Jaccard overlap"),
  diagonal_specificity_test(mp_gene_jaccard, "Current MP-state gene Jaccard overlap"),
  diagonal_specificity_test(combined_mp_projection, "Reciprocal current-MP projection")
)

state_expected_margins <- bind_rows(lapply(shared_state_order, function(state_name) {
  metric_mats <- list(
    "DGE-effect correlation" = effect_cor,
    "Top-marker Jaccard" = marker_jaccard,
    "Current MP-gene Jaccard" = mp_gene_jaccard,
    "Reciprocal MP projection" = combined_mp_projection
  )
  bind_rows(lapply(names(metric_mats), function(metric_name) {
    values <- metric_mats[[metric_name]][state_name, ]
    expected_sc_ref_name <- expected_state_map[[state_name]]
    expected_value <- unname(values[expected_sc_ref_name])
    alternatives <- values[names(values) != expected_sc_ref_name]
    best_alternative_state <- names(which.max(alternatives))
    data.frame(
      pdo_state = state_name,
      metric = metric_name,
      expected_sc_ref_state = expected_sc_ref_name,
      expected_value = expected_value,
      best_alternative_state = best_alternative_state,
      best_alternative_value = unname(alternatives[best_alternative_state]),
      expected_margin = expected_value - unname(alternatives[best_alternative_state]),
      expected_is_best = expected_value > unname(alternatives[best_alternative_state]),
      stringsAsFactors = FALSE
    )
  }))
}))

best_match_summary <- bind_rows(lapply(pdo_state_order, function(pdo_state) {
  effect_row <- effect_cor[pdo_state, ]
  marker_row <- marker_jaccard[pdo_state, ]
  mp_gene_row <- mp_gene_jaccard[pdo_state, ]
  data.frame(
    pdo_state = pdo_state,
    expected_sc_ref_state = ifelse(
      pdo_state %in% shared_state_order, expected_state_map[pdo_state], NA_character_
    ),
    effect_best_sc_ref_state = names(which.max(effect_row)),
    effect_best_value = max(effect_row),
    effect_second_value = sort(effect_row, decreasing = TRUE)[2],
    marker_best_sc_ref_state = names(which.max(marker_row)),
    marker_best_value = max(marker_row),
    marker_second_value = sort(marker_row, decreasing = TRUE)[2],
    mp_gene_best_sc_ref_state = names(which.max(mp_gene_row)),
    mp_gene_best_value = max(mp_gene_row),
    mp_gene_second_value = sort(mp_gene_row, decreasing = TRUE)[2],
    stringsAsFactors = FALSE
  )
}))
####################

####################
# Current scRef MP signatures applied to PDO cells: matched-proportion table
####################
if (is.null(pdo_data$projected_states) || is.null(pdo_data$native_states)) {
  stop(
    "The cache predates current scRef-signature projection. ",
    "Rerun with PDO_FORCE_REBUILD=1."
  )
}
projection_cells <- intersect(
  names(pdo_data$native_states),
  names(pdo_data$projected_states)
)
projection_cell_df <- data.frame(
  cell = projection_cells,
  pdo_native_state = as.character(pdo_data$native_states[projection_cells]),
  sc_ref_projected_state = as.character(pdo_data$projected_states[projection_cells]),
  stringsAsFactors = FALSE
) %>%
  filter(pdo_native_state %in% pdo_state_order)

projection_state_order <- c(sc_state_order, "Unresolved", "Hybrid")
projection_concordance <- projection_cell_df %>%
  count(pdo_native_state, sc_ref_projected_state, name = "cells") %>%
  tidyr::complete(
    pdo_native_state = pdo_state_order,
    sc_ref_projected_state = projection_state_order,
    fill = list(cells = 0L)
  ) %>%
  group_by(pdo_native_state) %>%
  mutate(
    total_pdo_cells = sum(cells),
    proportion = cells / total_pdo_cells,
    percentage = 100 * proportion,
    expected_match = pdo_native_state %in% shared_state_order &
      pdo_native_state == sc_ref_projected_state
  ) %>%
  ungroup()
####################

####################
# Persistent analysis tables
####################
write_matrix_csv <- function(mat, row_label, filename) {
  write.csv(
    data.frame(setNames(list(rownames(mat)), row_label), mat, check.names = FALSE),
    filename,
    row.names = FALSE
  )
}

table_files <- c(
  cell_counts = file.path(out_paths[["tables"]], "Auto_state_cell_counts_by_sample.csv"),
  mp_signatures = file.path(out_paths[["tables"]], "Auto_current_mp_state_signatures.csv"),
  pdo_mp_in_sc = file.path(out_paths[["tables"]], "Auto_PDO_MP_signatures_in_scRef_states_z.csv"),
  sc_mp_in_pdo = file.path(out_paths[["tables"]], "Auto_scRef_MP_signatures_in_PDO_states_z.csv"),
  dge = file.path(out_paths[["tables"]], "Auto_state_pseudobulk_DGE_all.csv"),
  effect_cor = file.path(out_paths[["tables"]], "Auto_state_DGE_effect_correlation.csv"),
  marker_overlap = file.path(out_paths[["tables"]], "Auto_state_marker_overlap.csv"),
  mp_gene_overlap = file.path(out_paths[["tables"]], "Auto_current_MP_state_gene_overlap.csv"),
  specificity = file.path(out_paths[["tables"]], "Auto_state_diagonal_specificity_tests.csv"),
  state_margins = file.path(out_paths[["tables"]], "Auto_state_expected_match_margins.csv"),
  best_match = file.path(out_paths[["tables"]], "Auto_state_best_match_summary.csv"),
  projected_concordance = file.path(
    out_paths[["tables"]], "Auto_scRef_signature_states_on_PDO_concordance.csv"
  )
)
projected_state_file <- file.path(
  out_dir, "Auto_scRef_current_MP_states_on_PDO.rds"
)

write.csv(
  bind_rows(pdo_data$cell_summary, sc_data$cell_summary),
  table_files[["cell_counts"]], row.names = FALSE
)
write.csv(signature_map, table_files[["mp_signatures"]], row.names = FALSE)
write_matrix_csv(
  pdo_signatures_in_sc_z, "PDO_source_state", table_files[["pdo_mp_in_sc"]]
)
write_matrix_csv(
  sc_signatures_in_pdo_z, "scRef_source_state", table_files[["sc_mp_in_pdo"]]
)
write.csv(
  bind_rows(pdo_dge_long, sc_dge_long),
  table_files[["dge"]], row.names = FALSE
)
write_matrix_csv(effect_cor, "PDO_state", table_files[["effect_cor"]])
write.csv(marker_stats, table_files[["marker_overlap"]], row.names = FALSE)
write.csv(mp_gene_overlap_stats, table_files[["mp_gene_overlap"]], row.names = FALSE)
write.csv(specificity_summary, table_files[["specificity"]], row.names = FALSE)
write.csv(state_expected_margins, table_files[["state_margins"]], row.names = FALSE)
write.csv(best_match_summary, table_files[["best_match"]], row.names = FALSE)
write.csv(
  projection_concordance,
  table_files[["projected_concordance"]],
  row.names = FALSE
)
saveRDS(pdo_data$projected_states, projected_state_file)

results$derived <- list(
  pdo_signatures_in_sc_z = pdo_signatures_in_sc_z,
  sc_signatures_in_pdo_z = sc_signatures_in_pdo_z,
  combined_mp_projection = combined_mp_projection,
  effect_cor = effect_cor,
  marker_stats = marker_stats,
  marker_jaccard = marker_jaccard,
  mp_gene_overlap_stats = mp_gene_overlap_stats,
  mp_gene_jaccard = mp_gene_jaccard,
  informative_genes = informative_genes,
  pdo_effects = pdo_effects,
  sc_effects = sc_effects,
  specificity_summary = specificity_summary,
  state_expected_margins = state_expected_margins,
  best_match_summary = best_match_summary,
  projection_concordance = projection_concordance
)
saveRDS(results, cache_file)
####################

####################
# Plot helpers: all figures intentionally omit subtitles
####################
matrix_to_long <- function(mat, row_name, col_name, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c(row_name, col_name, value_name)
  df
}

make_heatmap_plot <- function(
    mat, title_text, row_title, col_title,
    limits = NULL, midpoint = 0, low = "#2166AC", high = "#B2182B",
    digits = 2, expected_outline = TRUE, sequential = FALSE) {
  plot_df <- matrix_to_long(mat, "row_state", "col_state", "value")
  plot_df$row_state <- factor(plot_df$row_state, levels = rev(rownames(mat)))
  plot_df$col_state <- factor(plot_df$col_state, levels = colnames(mat))
  plot_df$label <- sprintf(paste0("%.", digits, "f"), plot_df$value)
  plot_df$expected <- expected_outline &
    as.character(plot_df$row_state) %in% names(expected_state_map) &
    expected_state_map[as.character(plot_df$row_state)] == as.character(plot_df$col_state)

  p <- ggplot(plot_df, aes(x = col_state, y = row_state, fill = value)) +
    geom_tile(aes(color = expected), linewidth = 1.2) +
    geom_text(aes(label = label), size = 3.8, fontface = "bold") +
    scale_color_manual(values = c("FALSE" = "white", "TRUE" = "black"), guide = "none") +
    labs(title = title_text, x = col_title, y = row_title, fill = NULL) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right"
    )
  if (sequential) {
    p + scale_fill_gradient(low = low, high = high)
  } else if (is.null(limits)) {
    p + scale_fill_gradient2(low = low, mid = "white", high = high, midpoint = midpoint)
  } else {
    p + scale_fill_gradient2(
      low = low, mid = "white", high = high, midpoint = midpoint, limits = limits
    )
  }
}

p_effect <- make_heatmap_plot(
  effect_cor,
  "Independent state-vs-rest transcriptional effect correlation",
  "PDO state",
  "scRef state",
  limits = c(-1, 1)
)
p_marker <- make_heatmap_plot(
  marker_jaccard,
  paste0("Overlap of independently ranked top ", top_marker_n, " positive markers"),
  "PDO state",
  "scRef state",
  low = "white",
  high = "#006D2C",
  sequential = TRUE
)
p_pdo_mp_sc <- make_heatmap_plot(
  pdo_signatures_in_sc_z,
  "Current PDO MP-state signatures projected into scRef",
  "PDO MP source state",
  "Independently assigned scRef state"
)
p_sc_mp_pdo <- make_heatmap_plot(
  sc_signatures_in_pdo_z,
  "Current scRef MP-state signatures projected into PDO",
  "scRef MP source state",
  "Independently assigned PDO state"
)
p_mp_gene <- make_heatmap_plot(
  mp_gene_jaccard,
  "Direct overlap of current MP-state genes",
  "PDO MP source state",
  "scRef MP source state",
  low = "white",
  high = "#6A3D9A",
  sequential = TRUE
)

specificity_long <- specificity_summary %>%
  select(metric, matched_mean, nonmatched_mean, exact_permutation_p) %>%
  pivot_longer(
    cols = c(matched_mean, nonmatched_mean),
    names_to = "comparison",
    values_to = "value"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      matched_mean = "Expected same-state pairs",
      nonmatched_mean = "Other cross-state pairs"
    ),
    metric = factor(metric, levels = rev(specificity_summary$metric))
  )

p_specificity <- ggplot(
  specificity_long,
  aes(x = value, y = metric, color = comparison, group = metric)
) +
  geom_line(color = "grey60", linewidth = 0.8) +
  geom_point(size = 4) +
  scale_color_manual(values = c(
    "Expected same-state pairs" = "#B2182B",
    "Other cross-state pairs" = "#2166AC"
  )) +
  geom_text(
    data = specificity_summary %>%
      mutate(metric = factor(metric, levels = rev(specificity_summary$metric))),
    aes(
      x = pmax(matched_mean, nonmatched_mean),
      y = metric,
      label = paste0("exact p=", format.pval(exact_permutation_p, digits = 2))
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 3.5
  ) +
  labs(
    title = "Expected same-state pairs versus all alternative pairings",
    x = "Metric value",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    legend.position = "bottom",
    plot.margin = margin(8, 80, 8, 8)
  ) +
  coord_cartesian(clip = "off")

margin_plot_df <- state_expected_margins %>%
  mutate(
    pdo_state = factor(pdo_state, levels = rev(shared_state_order)),
    metric = factor(metric, levels = unique(state_expected_margins$metric)),
    label = sprintf("%+.2f", expected_margin)
  )
margin_limit <- max(abs(margin_plot_df$expected_margin), na.rm = TRUE)
p_state_margins <- ggplot(
  margin_plot_df,
  aes(x = metric, y = pdo_state, fill = expected_margin)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-margin_limit, margin_limit),
    name = "Expected minus\nbest alternative"
  ) +
  labs(
    title = "State-specific evidence for the prespecified match",
    x = NULL,
    y = "PDO state"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )

projection_colors <- c(
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Motile-cilia differentiation" = "#F781BF",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)
projection_plot_df <- projection_concordance %>%
  mutate(
    pdo_native_state = factor(pdo_native_state, levels = pdo_state_order),
    sc_ref_projected_state = factor(
      sc_ref_projected_state,
      levels = rev(projection_state_order)
    )
  )
projection_match_labels <- projection_concordance %>%
  filter(expected_match) %>%
  mutate(
    pdo_native_state = factor(pdo_native_state, levels = pdo_state_order),
    label = paste0("Matched: ", sprintf("%.1f%%", percentage))
  )
p_projection_bar <- ggplot(
  projection_plot_df,
  aes(x = pdo_native_state, y = percentage, fill = sc_ref_projected_state)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.78) +
  geom_label(
    data = projection_match_labels,
    aes(x = pdo_native_state, y = 97, label = label),
    inherit.aes = FALSE,
    color = "black",
    fill = "white",
    label.size = 0.2,
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_manual(
    values = projection_colors,
    breaks = projection_state_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Current scRef MP signatures applied to independently defined PDO states",
    x = "PDO centred noreg state",
    y = "PDO cells (%)",
    fill = "scRef-signature projected state"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

effect_axis_limit <- max(
  abs(c(pdo_effects, sc_effects)),
  na.rm = TRUE
) * 1.03

make_effect_scatter <- function(pdo_state, sc_state) {
  plot_df <- data.frame(
    gene = informative_genes,
    pdo_logfc = pdo_effects[informative_genes, pdo_state],
    sc_ref_logfc = sc_effects[informative_genes, sc_state],
    stringsAsFactors = FALSE
  )
  pdo_tbl <- pdo_dge$tables[[pdo_state]]
  sc_tbl <- sc_dge$tables[[sc_state]]
  pdo_t <- setNames(abs(pdo_tbl$t), pdo_tbl$gene)
  sc_t <- setNames(abs(sc_tbl$t), sc_tbl$gene)
  plot_df$label_rank <- pdo_t[plot_df$gene] + sc_t[plot_df$gene]
  plot_df$label <- NA_character_
  concordant_direction <- sign(plot_df$pdo_logfc) == sign(plot_df$sc_ref_logfc) &
    plot_df$pdo_logfc != 0
  label_idx <- which(concordant_direction)[
    order(plot_df$label_rank[concordant_direction], decreasing = TRUE)
  ]
  label_idx <- head(label_idx, 8)
  plot_df$label[label_idx] <- plot_df$gene[label_idx]
  rho <- cor(plot_df$pdo_logfc, plot_df$sc_ref_logfc, method = "spearman")

  ggplot(plot_df, aes(x = sc_ref_logfc, y = pdo_logfc)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point(color = unname(state_colors[pdo_state]), size = 0.7, alpha = 0.25) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      color = "black", linetype = "dashed", linewidth = 0.5
    ) +
    geom_text_repel(
      aes(label = label), size = 2.2, max.overlaps = Inf,
      min.segment.length = 0, seed = 1, na.rm = TRUE
    ) +
    annotate(
      "text",
      x = -0.96 * effect_axis_limit,
      y = 0.96 * effect_axis_limit,
      hjust = 0,
      vjust = 1,
      label = sprintf("Spearman rho = %.2f", rho),
      size = 3, fontface = "bold"
    ) +
    coord_fixed(
      xlim = c(-effect_axis_limit, effect_axis_limit),
      ylim = c(-effect_axis_limit, effect_axis_limit),
      expand = FALSE
    ) +
    labs(
      title = sc_state,
      x = "scRef state-vs-rest logFC",
      y = "PDO state-vs-rest logFC"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(
        color = unname(state_colors[sc_state]),
        face = "bold", hjust = 0.5, size = 9
      ),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      plot.margin = margin(5, 4, 5, 4)
    )
}
####################

####################
# Single multi-page concordance PDF
####################
pdf_file <- file.path(
  out_paths[["figures"]], "Auto_PDO_scRef_centred_state_concordance.pdf"
)
grDevices::pdf(pdf_file, width = 20, height = 9, useDingbats = FALSE)
print(
  (p_effect | p_marker) +
    plot_annotation(
      title = "Independent PDO and scRef state concordance",
      theme = theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))
    )
)
print(
  (p_pdo_mp_sc | p_sc_mp_pdo | p_mp_gene) +
    plot_annotation(
      title = "Reciprocal projection of current centred refined MP state signatures",
      theme = theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5))
    )
)
print(p_state_margins | p_specificity)
print(p_projection_bar)
for (pdo_state in pdo_state_order) {
  effect_panels <- lapply(
    sc_state_order,
    function(sc_state) make_effect_scatter(pdo_state, sc_state)
  )
  print(
    wrap_plots(effect_panels, nrow = 1) +
      plot_annotation(
        title = paste0(
          "PDO state: ", pdo_state,
          " — transcriptional effects versus all scRef states"
        ),
        theme = theme(
          plot.title = element_text(
            color = unname(state_colors[pdo_state]),
            face = "bold",
            size = 19,
            hjust = 0.5
          )
        )
      )
  )
}
grDevices::dev.off()
####################

####################
# Persistent run summary
####################
output_files <- c(
  pdf_file,
  cache_file,
  projected_state_file,
  unname(table_files)
)
log_file <- pdo_write_run_summary(
  script = script_path,
  out_dir = out_dir,
  inputs = input_files,
  outputs = output_files,
  parameters = list(
    state_source = "independent centred refined noreg step-06 vectors",
    mp_source = "current centred refined step-05 state groupings and gene lists",
    minimum_cells_per_sample_state_pseudobulk = min_cells_pseudobulk,
    top_positive_markers_per_state = top_marker_n,
    informative_effect_genes_per_state_dataset = effect_gene_n,
    pdo_excluded_sample = PDO_EXCLUDED_SAMPLE,
    unmatched_states = "ECM-remodelling; Motile-cilia differentiation; Cancer-cell immune mimicry"
  ),
  cache = cache_policy,
  status = "completed"
)
cat("Completed current PDO/scRef centred-state concordance analysis.\n")
cat("PDF: ", pdf_file, "\n", sep = "")
cat("Run summary: ", log_file, "\n", sep = "")
####################
