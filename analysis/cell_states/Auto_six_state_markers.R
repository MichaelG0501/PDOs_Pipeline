####################
# Analysis registry:
#   Status: active canonical centred-state marker workflow
#   Script: analysis/cell_states/Auto_six_state_markers.R
#   Methodology: analysis/methodology/cell_states/Auto_centred_six_state_marker_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     PDOs_outs/PDOs_merged.rds
#     PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#   Outputs:
#     PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_marker_summary.csv
#     PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_markers_ranked.csv
#     PDOs_outs/Auto_six_state_markers/intermediate/*.rds
#     PDOs_outs/Auto_six_state_markers/figures/*.pdf
#   Downstream:
#     Surface-marker prioritization and marker-comparison workflows.
#   Cache/replot behavior:
#     Set PDO_FORCE_REBUILD=1 to ignore all cached embedding/DGE/specificity
#     objects. Cache files remain persistent live downstream inputs.
#   Run command: qsub Auto_run_centred_six_state_markers.sh
####################

####################
# Auto_six_state_markers.R
#
# Rebuild a six-state PDO embedding from the finalized states,
# then derive robust state markers using a sample-aware recurrent DGE workflow.
#
# Methodology adapted from scRef_Pipeline's Auto_six_state_markers.R
#
# Inputs:
#   PDOs_outs/PDOs_merged.rds
#   PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#
# Outputs:
#   PDOs_outs/Auto_six_state_markers/figures/Auto_six_state_umap.pdf
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_umap_embeddings.csv
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_cluster_state_table.csv
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_global_marker_screen.csv.gz
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_sample_state_eligibility.csv
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_per_sample_dge.csv.gz
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_marker_summary.csv
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_markers_ranked.csv
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_markers_final.csv
#   PDOs_outs/Auto_six_state_markers/tables/Auto_six_state_markers_top5_*.csv*
#   PDOs_outs/Auto_six_state_markers/figures/Auto_six_state_marker_heatmap.pdf
#   PDOs_outs/Auto_six_state_markers/reports/Auto_six_state_marker_methodology.md
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(parallel)
  library(data.table)
  library(grid)
})

####################
# setup
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(file.path(project_dir, "PDOs_outs"))

source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_config.R"))

out_dir <- "Auto_six_state_markers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
cache_dir <- file.path(out_dir, "intermediate")
reports_dir <- file.path(out_dir, "reports")
for (path in c(tables_dir, figures_dir, cache_dir, reports_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

set.seed(1234)
force_rebuild <- tolower(Sys.getenv("PDO_FORCE_REBUILD", unset = "0")) %in%
  c("1", "true", "yes")

####################
# constants
####################
state_order <- PDO_STATE_ORDER
state_cols <- PDO_STATE_COLORS[state_order]

params <- list(
  n_variable_features = 3000,
  n_pcs = 30,
  cluster_resolution = 0.5,
  min_cells_feature = 10,  # Adjusted for PDO (fewer cells than ref)
  min_cells_state = PDO_THRESHOLDS$marker_min_cells_state,
  min_cells_rest = PDO_THRESHOLDS$marker_min_cells_rest,
  candidate_pool_per_state = 1000,
  top_markers_per_state = 5,
  mc_cores = 1L
)

####################
# helpers
####################
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

pick_logfc_col <- function(df) {
  candidates <- c("avg_log2FC", "avg_logFC")
  out <- candidates[candidates %in% colnames(df)][1]
  if (is.na(out)) stop("Could not find a logFC column in the marker table.")
  out
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

run_sample_state_markers <- function(sample_id, sample_cells, obj, state_levels, candidate_map, min_cells_state, min_cells_rest) {
  sample_obj <- obj[, sample_cells]
  sample_meta <- sample_obj@meta.data
  sample_batch <- unique(as.character(sample_meta$batch))
  if (length(sample_batch) == 0 || all(is.na(sample_batch))) sample_batch <- "unknown"
  sample_batch <- sample_batch[1]
  Idents(sample_obj) <- "state"

  eligibility_rows <- vector("list", length(state_levels))
  marker_rows <- list()

  for (i in seq_along(state_levels)) {
    state_name <- state_levels[i]
    state_cells <- colnames(sample_obj)[which(sample_obj$state == state_name)]
    other_cells <- colnames(sample_obj)[which(sample_obj$state != state_name)]
    state_n <- length(state_cells)
    other_n <- length(other_cells)
    eligible <- state_n >= min_cells_state && other_n >= min_cells_rest

    eligibility_rows[[i]] <- data.frame(
      sample = sample_id,
      batch = sample_batch,
      state = state_name,
      state_cell_n = state_n,
      other_cell_n = other_n,
      eligible = eligible,
      stringsAsFactors = FALSE
    )

    if (!eligible) next

    features_use <- unique(candidate_map[[state_name]] %||% character())
    if (length(features_use) == 0) next
    features_use <- intersect(features_use, rownames(sample_obj))
    if (length(features_use) == 0) next
    other_states <- intersect(
      setdiff(state_levels, state_name),
      unique(as.character(sample_obj$state))
    )
    if (length(other_states) == 0) next

    res <- tryCatch(
      FindMarkers(
        object = sample_obj,
        ident.1 = state_name,
        ident.2 = other_states,
        features = features_use,
        logfc.threshold = 0,
        min.pct = 0,
        test.use = "wilcox",
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0) next

    logfc_col <- pick_logfc_col(res)
    res$gene <- rownames(res)
    res$avg_log2FC <- res[[logfc_col]]
    res$sample <- sample_id
    res$batch <- sample_batch
    res$state <- state_name
    res$state_cell_n <- state_n
    res$other_cell_n <- other_n
    res$pct_state <- res$pct.1 %||% NA_real_
    res$pct_other <- res$pct.2 %||% NA_real_
    res$pct_delta <- res$pct_state - res$pct_other

    marker_rows[[state_name]] <- res %>%
      select(
        gene,
        sample,
        batch,
        state,
        p_val,
        p_val_adj,
        avg_log2FC,
        pct_state,
        pct_other,
        pct_delta,
        state_cell_n,
        other_cell_n
      )
  }

  list(
    eligibility = bind_rows(eligibility_rows),
    markers = bind_rows(marker_rows)
  )
}

####################
# data loading
####################
message("Loading Seurat object and finalized state labels.")

# Using the canonical centred refined state vector.
pdos_all <- readRDS("PDOs_merged.rds")
state_labels <- readRDS(
  file.path("centred_mp_refinement", "centred_refined_noreg_states.rds")
)

DefaultAssay(pdos_all) <- "RNA"

# Re-derive batch column (Consistent with sample_abundance_pdo.R)
pdos_all$batch <- ifelse(grepl("_Treated_|_Untreated_", pdos_all$orig.ident), "New_batch", "Cynthia_batch")

common_cells <- intersect(colnames(pdos_all), names(state_labels))
state_labels <- state_labels[common_cells] # Retain names
keep_cells <- common_cells[state_labels %in% state_order]

message("Building a lean six-state count matrix.")

meta_state5 <- pdos_all@meta.data[keep_cells, c("orig.ident", "batch"), drop = FALSE]
meta_state5$state <- as.character(state_labels[keep_cells])

counts_mat <- GetAssayData(pdos_all, assay = "RNA", layer = "counts")[, keep_cells, drop = FALSE]
gene_detect_n <- Matrix::rowSums(counts_mat > 0)
keep_features <- rownames(counts_mat)[gene_detect_n >= params$min_cells_feature]
counts_mat <- counts_mat[keep_features, , drop = FALSE]

rm(pdos_all, state_labels, common_cells, gene_detect_n, keep_features)
invisible(gc())

pdos_state5 <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta_state5,
  assay = "RNA"
)

rm(counts_mat, meta_state5)
invisible(gc())

pdos_state5$state <- factor(as.character(pdos_state5$state), levels = state_order)
Idents(pdos_state5) <- "state"

state_counts <- table(pdos_state5$state)
sample_counts <- pdos_state5@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  count(state, orig.ident, batch, name = "n_cells")

####################
# six-state embedding
####################
pdos_state5_cache <- file.path(cache_dir, "pdos_state5_embedded.rds")

if (!force_rebuild && file.exists(pdos_state5_cache)) {
  message("Loading cached six-state embedded Seurat object.")
  pdos_state5 <- readRDS(pdos_state5_cache)
} else {
  message("Rebuilding six-state embedding and clustering.")

  pdos_state5 <- NormalizeData(pdos_state5, verbose = FALSE)
  pdos_state5 <- FindVariableFeatures(
    pdos_state5,
    selection.method = "vst",
    nfeatures = params$n_variable_features,
    verbose = FALSE
  )
  pdos_state5 <- ScaleData(
    pdos_state5,
    features = VariableFeatures(pdos_state5),
    verbose = FALSE
  )
  pdos_state5 <- RunPCA(
    pdos_state5,
    features = VariableFeatures(pdos_state5),
    npcs = params$n_pcs,
    verbose = FALSE
  )
  pdos_state5 <- FindNeighbors(
    pdos_state5,
    dims = seq_len(params$n_pcs),
    verbose = FALSE
  )
  pdos_state5 <- FindClusters(
    pdos_state5,
    resolution = params$cluster_resolution,
    verbose = FALSE
  )
  pdos_state5 <- RunUMAP(
    pdos_state5,
    dims = seq_len(params$n_pcs),
    verbose = FALSE
  )
  
  message("Caching six-state embedded Seurat object.")
  saveRDS(pdos_state5, pdos_state5_cache)
}

umap_df <- Embeddings(pdos_state5, reduction = "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    state = as.character(pdos_state5$state[match(cell, colnames(pdos_state5))]),
    sample = as.character(pdos_state5$orig.ident[match(cell, colnames(pdos_state5))]),
    batch = as.character(pdos_state5$batch[match(cell, colnames(pdos_state5))]),
    cluster = as.character(pdos_state5$seurat_clusters[match(cell, colnames(pdos_state5))])
  )

fwrite(
  umap_df,
  file.path(tables_dir, "Auto_six_state_umap_embeddings.csv")
)

cluster_state_table <- pdos_state5@meta.data %>%
  count(seurat_clusters, state, name = "n_cells") %>%
  group_by(seurat_clusters) %>%
  mutate(cluster_pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup()

fwrite(
  cluster_state_table,
  file.path(tables_dir, "Auto_six_state_cluster_state_table.csv")
)

p_state <- DimPlot(
  pdos_state5,
  reduction = "umap",
  group.by = "state",
  cols = state_cols,
  pt.size = 0.50,
  raster = TRUE
) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Finalized six-state PDO UMAP",
    subtitle = paste0(
      ncol(pdos_state5),
      " cells across ",
      length(unique(pdos_state5$orig.ident)),
      " samples / ",
      length(unique(pdos_state5$batch)),
      " batches"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_cluster <- DimPlot(
  pdos_state5,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.50,
  raster = TRUE
) +
  labs(
    title = "Unsupervised clusters on the six-state subset",
    subtitle = paste0("resolution = ", params$cluster_resolution)
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  filename = file.path(figures_dir, "Auto_six_state_umap.pdf"),
  plot = p_state + p_cluster + plot_layout(widths = c(1.1, 1)),
  width = 16,
  height = 7,
  useDingbats = FALSE
)

####################
# global marker screen
####################
global_markers_cache <- file.path(cache_dir, "global_marker_screen.rds")

if (!force_rebuild && file.exists(global_markers_cache)) {
  message("Loading cached global marker screen.")
  global_markers <- readRDS(global_markers_cache)
} else {
  message("Running pooled six-state descriptive enrichment screen.")

  expr_global <- GetAssayData(pdos_state5, assay = "RNA", layer = "data")

  global_markers <- bind_rows(lapply(state_order, function(state_name) {
    state_cells <- colnames(pdos_state5)[which(pdos_state5$state == state_name)]
    other_cells <- colnames(pdos_state5)[which(pdos_state5$state != state_name)]

    if (length(state_cells) == 0) return(NULL)

    mean_state <- Matrix::rowMeans(expr_global[, state_cells, drop = FALSE])
    mean_other <- Matrix::rowMeans(expr_global[, other_cells, drop = FALSE])
    pct_state <- Matrix::rowMeans(expr_global[, state_cells, drop = FALSE] > 0)
    pct_other <- Matrix::rowMeans(expr_global[, other_cells, drop = FALSE] > 0)

    data.frame(
      gene = rownames(expr_global),
      state = state_name,
      global_mean_state = mean_state,
      global_mean_other = mean_other,
      global_mean_diff = mean_state - mean_other,
      global_pct_state = pct_state,
      global_pct_other = pct_other,
      global_pct_delta = pct_state - pct_other,
      stringsAsFactors = FALSE
    )
  })) %>%
    arrange(state, desc(global_mean_diff), desc(global_pct_delta), desc(global_pct_state))
  
  message("Caching global marker screen.")
  saveRDS(global_markers, global_markers_cache)
}

fwrite(
  global_markers,
  file.path(tables_dir, "Auto_six_state_global_marker_screen.csv.gz")
)

candidate_map <- global_markers %>%
  filter(global_mean_diff > 0) %>%
  group_by(state) %>%
  arrange(desc(global_mean_diff), desc(global_pct_delta), desc(global_pct_state), .by_group = TRUE) %>%
  slice_head(n = params$candidate_pool_per_state) %>%
  summarise(candidate_genes = list(unique(gene)), .groups = "drop")

candidate_map <- setNames(candidate_map$candidate_genes, candidate_map$state)
candidate_map <- candidate_map[state_order]
candidate_map[sapply(candidate_map, is.null)] <- list(character())

if (exists("expr_global")) rm(expr_global)
invisible(gc())

####################
# per-sample DGE
####################
sample_dge_cache <- file.path(cache_dir, "per_sample_dge_ranked_top5.rds")

if (!force_rebuild && file.exists(sample_dge_cache)) {
  message("Loading cached per-sample DGE results.")
  cached_sample_res <- readRDS(sample_dge_cache)
  eligibility_df <- cached_sample_res$eligibility
  per_sample_dge <- cached_sample_res$markers
} else {
  message("Running per-sample recurrent DGE validation.")

  sample_cells_map <- split(colnames(pdos_state5), pdos_state5$orig.ident)
  sample_ids <- names(sample_cells_map)

  sample_res <- mclapply(
    sample_ids,
    function(sample_id) {
      run_sample_state_markers(
        sample_id = sample_id,
        sample_cells = sample_cells_map[[sample_id]],
        obj = pdos_state5,
        state_levels = state_order,
        candidate_map = candidate_map,
        min_cells_state = params$min_cells_state,
        min_cells_rest = params$min_cells_rest
      )
    },
    mc.cores = params$mc_cores
  )

  eligibility_df <- bind_rows(lapply(sample_res, `[[`, "eligibility")) %>%
    arrange(state, sample)
  
  per_sample_dge <- bind_rows(lapply(sample_res, `[[`, "markers"))

  if (nrow(per_sample_dge) == 0) {
    stop("Per-sample DGE produced no marker rows. Check sample/state eligibility and FindMarkers inputs.")
  }
  
  message("Caching per-sample DGE results.")
  saveRDS(list(eligibility = eligibility_df, markers = per_sample_dge), sample_dge_cache)
}

# Always write the eligibility CSV for downstream reference
fwrite(
  eligibility_df,
  file.path(tables_dir, "Auto_six_state_sample_state_eligibility.csv")
)

per_sample_dge <- per_sample_dge %>%
  mutate(
    hit_flag = !is.na(p_val_adj) &
      p_val_adj < 0.05 &
      avg_log2FC > 0
  ) %>%
  arrange(state, gene, sample)

fwrite(
  per_sample_dge,
  file.path(tables_dir, "Auto_six_state_per_sample_dge.csv.gz")
)

state_coverage <- eligibility_df %>%
  filter(eligible) %>%
  group_by(state) %>%
  summarise(
    eligible_sample_n = n_distinct(sample),
    eligible_batch_n = n_distinct(batch),
    .groups = "drop"
  )

marker_summary <- per_sample_dge %>%
  group_by(state, gene) %>%
  summarise(
    tested_sample_n = n_distinct(sample),
    tested_batch_n = n_distinct(batch),
    hit_sample_n = n_distinct(sample[hit_flag]),
    hit_batch_n = n_distinct(batch[hit_flag]),
    hit_sample_pct = hit_sample_n / tested_sample_n,
    median_log2FC_hit = safe_median(avg_log2FC[hit_flag]),
    median_pct_state_hit = safe_median(pct_state[hit_flag]),
    median_pct_other_hit = safe_median(pct_other[hit_flag]),
    median_pct_delta_hit = safe_median(pct_delta[hit_flag]),
    max_log2FC_hit = safe_max(avg_log2FC[hit_flag]),
    min_p_adj_hit = safe_min(p_val_adj[hit_flag]),
    .groups = "drop"
  ) %>%
  left_join(
    global_markers %>%
      select(state, gene, global_mean_state, global_mean_other, global_mean_diff,
             global_pct_state, global_pct_other, global_pct_delta),
    by = c("state", "gene")
  ) %>%
  left_join(state_coverage, by = "state") %>%
  mutate(
    sample_recurrence = ifelse(eligible_sample_n > 0, hit_sample_n / eligible_sample_n, 0),
    batch_recurrence = ifelse(eligible_batch_n > 0, hit_batch_n / eligible_batch_n, 0),
    reproducibility_score = 0.5 * sample_recurrence + 0.5 * batch_recurrence
  )

####################
# state-level specificity check
####################
specificity_cache <- file.path(cache_dir, "state_specificity.rds")

if (!force_rebuild && file.exists(specificity_cache)) {
  message("Loading cached state specificity results.")
  specificity_long <- readRDS(specificity_cache)
} else {
  message("Computing sample-aware state specificity summaries for candidate genes.")

  candidate_union <- sort(unique(marker_summary$gene))
  expr_data <- GetAssayData(pdos_state5, assay = "RNA", layer = "data")[candidate_union, , drop = FALSE]

  state_expr_list <- lapply(state_order, function(state_name) {
    eligible_samples <- eligibility_df %>%
      filter(state == state_name, eligible) %>%
      pull(sample) %>%
      unique()

    if (length(eligible_samples) == 0) {
      return(
        data.frame(
          gene = candidate_union,
          state = state_name,
          median_sample_state_expr = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }

    state_cells_by_sample <- split(
      colnames(pdos_state5)[pdos_state5$state == state_name & pdos_state5$orig.ident %in% eligible_samples],
      pdos_state5$orig.ident[pdos_state5$state == state_name & pdos_state5$orig.ident %in% eligible_samples]
    )

    state_cells_by_sample <- state_cells_by_sample[lengths(state_cells_by_sample) >= params$min_cells_state]

    if (length(state_cells_by_sample) == 0) {
      return(
        data.frame(
          gene = candidate_union,
          state = state_name,
          median_sample_state_expr = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }

    sample_means <- vapply(
      state_cells_by_sample,
      function(cells_use) Matrix::rowMeans(expr_data[, cells_use, drop = FALSE]),
      FUN.VALUE = numeric(length(candidate_union))
    )

    if (is.null(dim(sample_means))) {
      sample_means <- matrix(sample_means, ncol = 1)
    }

    data.frame(
      gene = candidate_union,
      state = state_name,
      median_sample_state_expr = apply(sample_means, 1, median, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  state_expr_wide <- bind_rows(state_expr_list) %>%
    tidyr::pivot_wider(
      names_from = state,
      values_from = median_sample_state_expr
    )

  state_expr_mat <- state_expr_wide %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  best_state_idx <- max.col(state_expr_mat, ties.method = "first")
  best_state <- colnames(state_expr_mat)[best_state_idx]

  specificity_long <- bind_rows(lapply(state_order, function(state_name) {
    state_vals <- state_expr_mat[, state_name]
    other_vals <- state_expr_mat[, setdiff(state_order, state_name), drop = FALSE]
    other_max <- apply(other_vals, 1, max, na.rm = TRUE)

    data.frame(
      gene = rownames(state_expr_mat),
      state = state_name,
      state_median_expr = state_vals,
      off_state_max_median_expr = other_max,
      specificity_gap = state_vals - other_max,
      best_state = best_state,
      stringsAsFactors = FALSE
    )
  }))
  
  message("Caching state specificity results.")
  saveRDS(specificity_long, specificity_cache)
}

# Re-derive state_expr_mat from specificity_long for downstream use (heatmap_expr relies on it)
state_expr_mat <- specificity_long %>%
  select(gene, state, state_median_expr) %>%
  tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

marker_summary <- marker_summary %>%
  left_join(specificity_long, by = c("state", "gene")) %>%
  mutate(
    best_state_match = best_state == state
  )

fwrite(
  marker_summary,
  file.path(tables_dir, "Auto_six_state_marker_summary.csv")
)

####################
# final marker definition
####################
message("Selecting top publication-facing markers by reproducibility, effect size, and specificity.")

ranked_markers <- marker_summary %>%
  filter(
    hit_sample_n > 0,
    best_state_match,
    specificity_gap > 0
  ) %>%
  group_by(state) %>%
  mutate(
    reproducibility_rank = dplyr::percent_rank(reproducibility_score),
    effect_rank = dplyr::percent_rank(median_log2FC_hit),
    specificity_rank = dplyr::percent_rank(specificity_gap),
    ranking_score = reproducibility_rank + effect_rank + specificity_rank
  ) %>%
  ungroup() %>%
  arrange(
    state,
    desc(ranking_score),
    desc(reproducibility_score),
    desc(median_log2FC_hit),
    desc(specificity_gap),
    gene
  )

fwrite(
  ranked_markers,
  file.path(tables_dir, "Auto_six_state_markers_ranked.csv")
)

final_markers <- ranked_markers %>%
  group_by(state) %>%
  slice_head(n = params$top_markers_per_state) %>%
  ungroup()

fwrite(
  final_markers,
  file.path(tables_dir, "Auto_six_state_markers_final.csv")
)

####################
# recurrence reporting
####################
message("Writing sample/batch recurrence tables for the top markers.")

top_marker_recurrence_summary <- final_markers %>%
  mutate(
    # Adjusted thresholds for PDO (smaller dataset)
    legacy_required_sample_n = pmax(2, ceiling(0.20 * eligible_sample_n)),
    legacy_required_batch_n = pmax(1, ceiling(0.35 * eligible_batch_n)),
    is_multi_sample = hit_sample_n >= 2,
    is_multi_batch = hit_batch_n >= 2,
    support_class = case_when(
      is_multi_batch ~ "multi-batch",
      is_multi_sample ~ "multi-sample_single-batch",
      hit_sample_n == 1 ~ "single-sample_single-batch",
      TRUE ~ "no-positive-hit"
    ),
    passes_legacy_strict_recurrence =
      hit_sample_n >= legacy_required_sample_n &
      hit_batch_n >= legacy_required_batch_n
  ) %>%
  select(
    state,
    gene,
    hit_sample_n,
    eligible_sample_n,
    sample_recurrence,
    hit_batch_n,
    eligible_batch_n,
    batch_recurrence,
    reproducibility_score,
    median_log2FC_hit,
    specificity_gap,
    ranking_score,
    is_multi_sample,
    is_multi_batch,
    support_class,
    legacy_required_sample_n,
    legacy_required_batch_n,
    passes_legacy_strict_recurrence
  )

fwrite(
  top_marker_recurrence_summary,
  file.path(tables_dir, "Auto_six_state_markers_top5_recurrence_summary.csv")
)

top_marker_sample_support <- per_sample_dge %>%
  semi_join(
    final_markers %>% select(state, gene),
    by = c("state", "gene")
  ) %>%
  select(
    state,
    gene,
    sample,
    batch,
    hit_flag,
    p_val_adj,
    avg_log2FC,
    pct_state,
    pct_other,
    pct_delta,
    state_cell_n,
    other_cell_n
  ) %>%
  arrange(state, gene, desc(hit_flag), batch, sample)

fwrite(
  top_marker_sample_support,
  file.path(tables_dir, "Auto_six_state_markers_top5_sample_support.csv.gz")
)

top_marker_batch_support <- top_marker_sample_support %>%
  group_by(state, gene, batch) %>%
  summarise(
    eligible_sample_n_in_batch = n_distinct(sample),
    hit_sample_n_in_batch = n_distinct(sample[hit_flag]),
    batch_hit_flag = hit_sample_n_in_batch > 0,
    best_log2FC_in_batch = safe_max(avg_log2FC[hit_flag]),
    .groups = "drop"
  )

fwrite(
  top_marker_batch_support,
  file.path(tables_dir, "Auto_six_state_markers_top5_batch_support.csv")
)

####################
# heatmap construction (exact scRef style)
####################
message("Building marker heatmap with scRef style.")

if (nrow(final_markers) == 0) {
  stop("No final markers passed the current recurrent marker thresholds.")
}

heatmap_genes <- final_markers$gene
heatmap_expr <- state_expr_mat[heatmap_genes, state_order, drop = FALSE]

dup_counter <- ave(seq_along(heatmap_genes), heatmap_genes, FUN = seq_along)
heatmap_row_names <- ifelse(
  dup_counter == 1,
  heatmap_genes,
  paste0(heatmap_genes, "__", dup_counter)
)

rownames(heatmap_expr) <- heatmap_row_names
heatmap_z <- row_zscore(heatmap_expr)

heatmap_state_factor <- factor(final_markers$state, levels = state_order)
names(heatmap_state_factor) <- heatmap_row_names

# Publication-facing labels and annotations
# FIX: Use explicitly named vectors to prevent ComplexHeatmap from getting confused
row_ann <- rowAnnotation(
  State = heatmap_state_factor,
  Sample_support = anno_barplot(
    final_markers$hit_sample_pct * 100,
    gp = gpar(fill = "grey35", col = NA),
    border = FALSE,
    axis_param = list(side = "bottom", at = c(0, 50), labels = c("0", "50")),
    width = unit(18, "mm")
  ),
  Batch_support = anno_barplot(
    final_markers$hit_batch_n,
    gp = gpar(fill = "grey55", col = NA),
    border = FALSE,
    axis_param = list(side = "bottom"),
    width = unit(15, "mm")
  ),
  col = list(State = state_cols),
  # Use named vectors so the package knows exactly which rule applies to which column
  show_annotation_name = c(State = FALSE, Sample_support = TRUE, Batch_support = TRUE), 
  annotation_label = c(State = "State", Sample_support = "Sample\nSupport", Batch_support = "Batch\nSupport"), 
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 9, lineheight = 0.9),
  gap = unit(c(2, 6), "mm"), 
  simple_anno_size = unit(4, "mm")
)

# FIX: Renamed to 'Top_State' to completely prevent the legend-merging bug
top_ann <- HeatmapAnnotation(
  Top_State = factor(state_order, levels = state_order),
  col = list(Top_State = state_cols),
  show_annotation_name = FALSE, 
  show_legend = FALSE, # We only need the legend from row_ann
  simple_anno_size = unit(4, "mm")
)

col_fun <- colorRamp2(
  c(-2.0, 0, 2.5),
  c("#1D4E89", "#F8F4EC", "#B22222") # scRef specific colors
)

ht <- Heatmap(
  heatmap_z,
  name = "Row Z-score",
  top_annotation = top_ann,
  left_annotation = row_ann,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = heatmap_state_factor,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, 
  border = TRUE,
  heatmap_legend_param = list(
    title = "Relative\nexpression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

pdf(
  file.path(figures_dir, "Auto_six_state_marker_heatmap.pdf"),
  width = 11.5,
  height = 11, 
  useDingbats = FALSE
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(2.5, 1), c("cm", "null")))))

# Title block
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text(
  "Recurrent state markers across six finalized PDO states",
  x = unit(0.5, "npc"),
  y = unit(0.70, "npc"), 
  gp = gpar(fontsize = 16, fontface = "bold")
)
grid.text(
  "Rows are per-state recurrent markers ranked by recurrence, effect size, and state specificity. Columns are median sample-level expression per state, scaled by row.",
  x = unit(0.5, "npc"),
  y = unit(0.30, "npc"), 
  gp = gpar(fontsize = 10)
)
popViewport()

# Heatmap block
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht, newpage = FALSE)
popViewport(2)
dev.off()

####################
# ============ EXCEL OUTPUT ============
####################
message("Building Excel output for top markers and complete marker lists...")

suppressPackageStartupMessages(library(openxlsx))

wb <- createWorkbook()

# Define common styles
header_style <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                            halign = "center", valign = "center", wrapText = TRUE,
                            border = "TopBottomLeftRight", borderStyle = "thin",
                            borderColour = "#2C3E50", bgFill = "#ECF0F1")
sc_header_style <- createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                               halign = "center", valign = "center", fontColour = "#FFFFFF",
                               fgFill = "#34495E", border = "TopBottomLeftRight", borderColour = "#2C3E50")
expr_sc_header_style <- createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold",
                                    halign = "center", valign = "center", fontColour = "#FFFFFF",
                                    fgFill = "#2980B9", border = "TopBottomLeftRight", borderColour = "#2C3E50")
sep_style <- createStyle(bgFill = "#ECF0F1")
pct_fmt <- createStyle(numFmt = "0.0%")
num3 <- createStyle(numFmt = "0.000")
num1 <- createStyle(numFmt = "0.0")

# Helper function to build a formatted worksheet
build_marker_sheet <- function(wb, sheet_name, df, state_cols_map, state_expr_mat, state_order) {
  addWorksheet(wb, sheet_name)
  
  excel_df <- data.frame()
  state_block_rows <- list()
  current_excel_row <- 3
  row_colors <- character()
  is_data_row <- logical()
  
  states_present <- intersect(state_order, unique(df$state))
  
  for (i in seq_along(states_present)) {
    st <- states_present[i]
    st_df <- df %>% filter(state == st) %>% arrange(desc(ranking_score))
    if (nrow(st_df) == 0) next
    
    st_genes <- st_df$gene
    
    st_data <- data.frame(
      StateLabel = rep(as.character(st), length(st_genes)),
      Gene = st_genes, 
      stringsAsFactors = FALSE
    )
    
    # Recurrence & Specificity & Pct
    rec_data <- st_df %>% 
      select(gene, 
             sample_recurrence = hit_sample_pct, 
             batch_recurrence = hit_batch_n, 
             median_log2FC = median_log2FC_hit,
             specificity_gap,
             ranking_score,
             pct_target = global_pct_state,
             pct_other = global_pct_other)
    
    st_data <- left_join(st_data, rec_data, by = c("Gene" = "gene"))
    
    # Separator
    st_data$Sep <- ""
    
    # Expression Group
    sc_e <- state_expr_mat[match(st_genes, rownames(state_expr_mat)), state_order, drop = FALSE]
    colnames(sc_e) <- state_order
    st_data <- cbind(st_data, sc_e)
    
    block_start <- current_excel_row
    block_end <- current_excel_row + length(st_genes) - 1
    state_block_rows[[as.character(st)]] <- c(block_start, block_end)
    
    excel_df <- rbind(excel_df, st_data)
    
    st_col <- state_cols_map[as.character(st)]
    if (is.na(st_col)) st_col <- "#000000"
    row_colors <- c(row_colors, rep(st_col, nrow(st_data)))
    is_data_row <- c(is_data_row, rep(TRUE, nrow(st_data)))
    
    current_excel_row <- block_end + 1
    
    # Empty Row separator
    if (i < length(states_present)) {
      empty_row <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(st_data)))
      colnames(empty_row) <- colnames(st_data)
      empty_row$StateLabel <- ""
      empty_row$Gene <- ""
      empty_row$Sep <- ""
      excel_df <- rbind(excel_df, empty_row)
      row_colors <- c(row_colors, NA)
      is_data_row <- c(is_data_row, FALSE)
      current_excel_row <- current_excel_row + 1
    }
  }
  
  if (nrow(excel_df) == 0) {
    writeData(wb, sheet_name, "No markers found", startCol=1, startRow=1)
    return()
  }
  
  # Column indices
  gene_col <- 2
  recur_cols <- 3:9
  sep_col <- 10
  expr_start <- 11
  expr_end <- 10 + length(state_order)
  
  # Row 1 Headers
  writeData(wb, sheet_name, "State", startCol=1, startRow=1)
  mergeCells(wb, sheet_name, cols=1, rows=1:2)
  writeData(wb, sheet_name, "Gene", startCol=2, startRow=1)
  mergeCells(wb, sheet_name, cols=2, rows=1:2)
  
  writeData(wb, sheet_name, "Recurrence & Specificity", startCol=recur_cols[1], startRow=1)
  mergeCells(wb, sheet_name, cols=recur_cols, rows=1)
  addStyle(wb, sheet_name, sc_header_style, rows=1, cols=recur_cols, gridExpand=TRUE, stack=TRUE)
  
  writeData(wb, sheet_name, "Mean Expression (Sample Median)", startCol=expr_start, startRow=1)
  mergeCells(wb, sheet_name, cols=expr_start:expr_end, rows=1)
  addStyle(wb, sheet_name, expr_sc_header_style, rows=1, cols=expr_start:expr_end, gridExpand=TRUE, stack=TRUE)
  
  # Row 2 Headers
  headers2 <- c("State", "Gene", 
                "Sample\nRecurrence", "Batch\nRecurrence", "Median\nlog2FC", 
                "Specificity\nGap", "Ranking\nScore", "% Target\nExpr", "% Other\nExpr", 
                "")
  headers2 <- c(headers2, state_order)
  
  if (length(states_present) == 1) {
    target_idx <- which(state_order == states_present[1])
    if (length(target_idx) == 1) {
      headers2[expr_start - 1 + target_idx] <- paste0(headers2[expr_start - 1 + target_idx], "\n(TARGET)")
    }
  }

  writeData(wb, sheet_name, t(headers2), startCol=1, startRow=2, colNames=FALSE)
  addStyle(wb, sheet_name, header_style, rows=2, cols=1:ncol(excel_df), gridExpand=TRUE, stack=TRUE)
  addStyle(wb, sheet_name, sc_header_style, rows=2, cols=recur_cols, gridExpand=TRUE, stack=TRUE)
  addStyle(wb, sheet_name, expr_sc_header_style, rows=2, cols=expr_start:expr_end, gridExpand=TRUE, stack=TRUE)
  
  # Data
  writeData(wb, sheet_name, excel_df, startCol=1, startRow=3, colNames=FALSE)
  
  # Styling for Column 1 and 2
  for (st_name in names(state_block_rows)) {
    rows <- state_block_rows[[st_name]]
    mergeCells(wb, sheet_name, cols=1, rows=rows[1]:rows[2])
    
    st_col_hex <- state_cols_map[st_name]
    if (is.na(st_col_hex)) st_col_hex <- "#000000"
    
    vert_style <- createStyle(textRotation = 90, halign = "center", valign = "center", textDecoration = "bold",
                              fontColour = "#FFFFFF", fgFill = st_col_hex, border = "TopBottomLeftRight")
    addStyle(wb, sheet_name, vert_style, rows=rows[1]:rows[2], cols=1, stack=TRUE)
  }
  
  for (r in seq_len(nrow(excel_df))) {
    if (is_data_row[r]) {
      st_col_hex <- row_colors[r]
      gene_st_style <- createStyle(textDecoration = "bold", fontName = "Consolas", fontColour = st_col_hex)
      addStyle(wb, sheet_name, gene_st_style, rows=r+2, cols=2, stack=TRUE)
    }
  }
  
  # Metric formatting
  addStyle(wb, sheet_name, pct_fmt, rows=3:(nrow(excel_df)+2), cols=c(3, 8, 9), gridExpand=TRUE, stack=TRUE)
  addStyle(wb, sheet_name, num1, rows=3:(nrow(excel_df)+2), cols=c(5, 6, 7), gridExpand=TRUE, stack=TRUE)
  addStyle(wb, sheet_name, createStyle(numFmt = "0"), rows=3:(nrow(excel_df)+2), cols=4, gridExpand=TRUE, stack=TRUE)
  addStyle(wb, sheet_name, num3, rows=3:(nrow(excel_df)+2), cols=expr_start:expr_end, gridExpand=TRUE, stack=TRUE)
  
  # Separator and borders
  addStyle(wb, sheet_name, sep_style, rows=1:(nrow(excel_df)+2), cols=sep_col, stack=TRUE)
  border_style_med <- createStyle(border = "Left", borderStyle = "medium", borderColour = "#2C3E50")
  addStyle(wb, sheet_name, border_style_med, rows=1:(nrow(excel_df)+2), cols=c(2, 3, sep_col, expr_start), gridExpand=TRUE, stack=TRUE)
  
  # Conditional formatting for expression
  expr_vals <- unlist(excel_df[is_data_row, expr_start:expr_end])
  expr_vals <- expr_vals[is.finite(expr_vals)]
  if (length(expr_vals) > 0) {
    q_low <- quantile(expr_vals, 0.05, na.rm=TRUE)
    q_mid <- quantile(expr_vals, 0.40, na.rm=TRUE)
    q_high <- quantile(expr_vals, 0.85, na.rm=TRUE)
    conditionalFormatting(wb, sheet_name, cols=expr_start:expr_end, rows=3:(nrow(excel_df)+2),
                          style=c("#FFFFFF", "#FB8A8A", "#B22222"), rule=c(q_low, q_mid, q_high), type="colourScale")
  }
  
  # Layout
  setColWidths(wb, sheet_name, cols=1, widths=5)
  setColWidths(wb, sheet_name, cols=2, widths=20)
  setColWidths(wb, sheet_name, cols=recur_cols, widths=13)
  setColWidths(wb, sheet_name, cols=sep_col, widths=3)
  setColWidths(wb, sheet_name, cols=expr_start:expr_end, widths=16)
  freezePane(wb, sheet_name, firstActiveRow=3, firstActiveCol=3)
}

# 1. Build Top 5 Markers sheet
top5_df <- final_markers
build_marker_sheet(wb, "Top 5 Markers", top5_df, state_cols, state_expr_mat, state_order)

# 2. Build individual sheets for all markers per state
# Note: Excel worksheet names are limited to 31 chars, so we truncate if necessary
for (st in state_order) {
  st_df <- ranked_markers %>% filter(state == st)
  if (nrow(st_df) > 0) {
    sheet_name_st <- substr(as.character(st), 1, 31)
    build_marker_sheet(wb, sheet_name_st, st_df, state_cols, state_expr_mat, state_order)
  }
}

out_xlsx_final <- file.path(tables_dir, "Auto_six_state_markers_all.xlsx")
saveWorkbook(wb, out_xlsx_final, overwrite = TRUE)
message("Saved Excel to: ", out_xlsx_final)

####################
# methodology markdown
####################
message("Generating methodology documentation.")

methodology_text <- c(
  "# Six-State PDO Marker Methodology",
  "",
  "This document describes the sample-aware recurrent DGE workflow used to derive robust markers for the six finalized PDO states.",
  "",
  "## 1. Subset & Re-embedding",
  "- Epithelial cells were subsetted to include only those assigned to one of the six finalized states.",
  "- A lean Seurat object was built, excluding 'Unresolved' and 'Hybrid' cells.",
  "- The subset was re-normalized and re-embedded (PCA, UMAP) to resolve state-defining features without interference from excluded populations.",
  "",
  "## 2. Global Marker Screen",
  "- A pooled enrichment screen was performed to identify the top 1000 candidate genes per state based on global mean difference and delta percent expressed.",
  "- These candidates served as the search space for per-sample validation to reduce computational overhead.",
  "",
  "## 3. Per-Sample Recurrent DGE",
  "- Standard Wilcoxon tests were performed for each state within each eligible sample (minimum 10 cells in state, 10 cells in rest).",
  "- Markers were required to be significant (FDR < 0.05) and upregulated (log2FC > 0) in the specific sample.",
  "- **Sample Recurrence**: Fraction of eligible samples where the gene was a hit.",
  "- **Batch Recurrence**: Fraction of eligible batches where the gene was a hit.",
  "",
  "## 4. State Specificity Gap",
  "- For each gene, the median expression across samples was calculated for every state.",
  "- **Specificity Gap**: The difference between the expression in the target state and the expression in the next-best state.",
  "- Only genes where the target state was the primary expression site (gap > 0) were considered.",
  "",
  "## 5. Ranking & Final Selection",
  "- Genes were ranked by a composite **Ranking Score**: `percent_rank(reproducibility) + percent_rank(effect_size) + percent_rank(specificity_gap)`.",
  "- The top 5 genes per state were selected as finalized markers for publication-quality visualization and downstream analysis.",
  ""
)

writeLines(methodology_text, file.path(reports_dir, "Auto_six_state_marker_methodology.md"))

message("=== Pipeline Complete ===")
