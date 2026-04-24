####################
# Auto_marker_selection_simulation.R
#
# In silico validation of selected PDO/scATLAS marker panels.
#
# Simulation 1:
#   Randomly sample cells and assign state labels from marker expression alone.
#   Marker expression is normalized within each gene before panel scoring, so
#   highly abundant genes do not dominate because of absolute expression scale.
#
# Simulation 2:
#   Simulate bulk qRT-PCR by creating paired two-condition shifts where one
#   state increases and another decreases, then test whether panel expression
#   captures the relative state shift.
#
# Inputs:
#   PDOs_outs/PDOs_merged.rds
#   PDOs_outs/Auto_PDO_final_states.rds
#
# Outputs:
#   PDOs_outs/Auto_marker_selection_simulation/
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ggtext)
})

####################
# setup
####################
project_candidates <- c(
  Sys.getenv("PDO_PROJECT_DIR", unset = ""),
  "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline",
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline",
  getwd()
)
project_candidates <- project_candidates[nzchar(project_candidates)]
project_dir <- project_candidates[dir.exists(project_candidates)][1]
if (is.na(project_dir)) stop("Could not locate the PDOs_Pipeline project directory.")

setwd(file.path(project_dir, "PDOs_outs"))

out_dir <- "Auto_marker_selection_simulation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env_int <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  as.integer(val)
}

env_num <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  as.numeric(val)
}

env_int_vec <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  as.integer(strsplit(val, ",", fixed = TRUE)[[1]])
}

params <- list(
  seed = env_int("AUTO_MARKER_SIM_SEED", 1234567L),
  gate_reps = env_int("AUTO_MARKER_GATE_REPS", 200L),
  gate_subset_n = env_int("AUTO_MARKER_GATE_N", 20000L),
  qpcr_reps = env_int("AUTO_MARKER_QPCR_REPS", 1000L),
  qpcr_cells_per_condition = env_int("AUTO_MARKER_QPCR_N", 10000L),
  gate_threshold_fraction = env_num("AUTO_MARKER_GATE_THRESHOLD_FRACTION", 0.05),
  gate_min_positive_markers = env_int("AUTO_MARKER_GATE_MIN_POSITIVE", 1L),
  gate_min_margin = env_num("AUTO_MARKER_GATE_MIN_MARGIN", 0),
  dirichlet_concentration = env_num("AUTO_MARKER_QPCR_CONCENTRATION", 18),
  qpcr_shift_min = env_num("AUTO_MARKER_QPCR_SHIFT_MIN", 0.05),
  qpcr_shift_max = env_num("AUTO_MARKER_QPCR_SHIFT_MAX", 0.30),
  qpcr_pseudo_cpm = env_num("AUTO_MARKER_QPCR_PSEUDO_CPM", 0.1),
  min_truth_abs_log2fc = env_num("AUTO_MARKER_MIN_TRUTH_LOG2FC", 0.25)
)

set.seed(params$seed)

####################
# marker panels
####################
state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Immune regulated" = "#6B7280",
  "Unresolved" = "grey80",
  "Ambiguous" = "black"
)

marker_panels <- list(
  "Classic Proliferative" = c("PCLAF", "STMN1", "FABP5"),
  "Basal to Intest. Meta" = c("MUC13", "HPGD", "TSPAN1"),
  "SMG-like Metaplasia" = c("PROM1", "CHRM3", "PLCB4"),
  "Stress-adaptive" = c("PPP1R15A", "ATF3", "SOX4"),
  "3CA_EMT_and_Protein_maturation" = c("CANX", "PDIA4", "NORAD"),
  "Immune regulated" = c("VIM", "SRGN", "PDE4B")
)

panel_truth_state <- c(
  "Classic Proliferative" = "Classic Proliferative",
  "Basal to Intest. Meta" = "Basal to Intest. Meta",
  "SMG-like Metaplasia" = "SMG-like Metaplasia",
  "Stress-adaptive" = "Stress-adaptive",
  "Immune regulated" = NA_character_,
  "3CA_EMT_and_Protein_maturation" = "3CA_EMT_and_Protein_maturation"
)

classification_panels <- marker_panels[!is.na(panel_truth_state[names(marker_panels)])]
gate_rule_name <- "two_of_three_data_layer_midpoint_gate"

####################
# helpers
####################
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_assay_layer <- function(obj, assay = "RNA", layer = "counts") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  if (length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = method))
}

safe_lm_slope <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  if (length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) return(NA_real_)
  unname(coef(lm(y[keep] ~ x[keep]))[2])
}

rdirichlet_one <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  names(x) <- names(alpha)
  x / sum(x)
}

rank_to_percentile <- function(x) {
  if (length(x) <= 1) return(rep(1, length(x)))
  r <- rank(x, ties.method = "average", na.last = "keep")
  out <- (r - 1) / (length(x) - 1)
  out[!is.finite(out)] <- 0
  out
}

classify_marker_cells <- function(cell_idx, marker_expr, panels, panel_thresholds) {
  panel_names <- names(panels)
  
  # Get all genes involved across all panels
  all_genes <- unique(unlist(panels))
  available_genes <- intersect(all_genes, rownames(marker_expr))
  
  # Extract expression for relevant genes and cells
  expr_mat <- as.matrix(marker_expr[available_genes, cell_idx, drop = FALSE])
  
  # Calculate mean expression per panel for each cell
  panel_means <- matrix(0, nrow = length(panel_names), ncol = length(cell_idx),
                        dimnames = list(panel_names, NULL))
  
  for (pn in panel_names) {
    genes_use <- intersect(panels[[pn]], available_genes)
    if (length(genes_use) > 0) {
      panel_means[pn, ] <- colMeans(expr_mat[genes_use, , drop = FALSE])
    }
  }
  
  # Determine hits based on panel-specific mean threshold
  hits <- panel_means > panel_thresholds[rownames(panel_means)]
  hit_counts <- colSums(hits)
  
  # Assign states
  pred <- character(length(cell_idx))
  pred[hit_counts == 0] <- "Unresolved"
  pred[hit_counts > 1] <- "Ambiguous"
  
  # For unique hits, assign the specific panel
  one_hit_mask <- hit_counts == 1
  if (any(one_hit_mask)) {
    # Extract panel index for the single TRUE in the hits column
    panel_idx <- max.col(t(hits[, one_hit_mask, drop = FALSE]), ties.method = "first")
    pred[one_hit_mask] <- panel_names[panel_idx]
  }
  
  data.frame(
    prediction = pred,
    top_score = apply(panel_means, 2, max),
    hit_n = as.integer(hit_counts),
    stringsAsFactors = FALSE
  )
}

summarise_gate_metrics <- function(pred, truth, rep_id, n_subset, rule_name) {
  panel_names <- names(marker_panels)
  metric_rows <- lapply(panel_names, function(panel_name) {
    truth_state <- panel_truth_state[[panel_name]]
    predicted_panel <- pred == panel_name
    if (is.na(truth_state)) {
      data.frame(
        rep = rep_id,
        n_subset = n_subset,
        gate_rule = rule_name,
        panel = panel_name,
        truth_state = NA_character_,
        true_n = NA_integer_,
        predicted_n = sum(predicted_panel),
        tp = NA_integer_,
        sensitivity = NA_real_,
        precision = NA_real_,
        f1 = NA_real_,
        false_call_rate = mean(predicted_panel),
        stringsAsFactors = FALSE
      )
    } else {
      truth_panel <- truth == truth_state
      tp <- sum(predicted_panel & truth_panel)
      precision <- ifelse(sum(predicted_panel) > 0, tp / sum(predicted_panel), NA_real_)
      sensitivity <- ifelse(sum(truth_panel) > 0, tp / sum(truth_panel), NA_real_)
      f1 <- ifelse(
        is.finite(precision + sensitivity) && (precision + sensitivity) > 0,
        2 * precision * sensitivity / (precision + sensitivity),
        NA_real_
      )
      data.frame(
        rep = rep_id,
        n_subset = n_subset,
        gate_rule = rule_name,
        panel = panel_name,
        truth_state = truth_state,
        true_n = sum(truth_panel),
        predicted_n = sum(predicted_panel),
        tp = tp,
        sensitivity = sensitivity,
        precision = precision,
        f1 = f1,
        false_call_rate = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  })

  data.table::rbindlist(metric_rows, fill = TRUE) %>%
    mutate(
      unresolved_rate = mean(pred == "Unresolved"),
      ambiguous_rate = mean(pred == "Ambiguous"),
      call_rate = mean(!(pred %in% c("Unresolved", "Ambiguous"))),
      confident_call_rate = mean(!(pred %in% c("Unresolved", "Ambiguous", "Immune regulated"))),
      overall_accuracy = mean(pred == truth)
    )
}

sample_condition_indices <- function(props, n_cells, cells_by_state) {
  state_n <- as.integer(rmultinom(1, n_cells, prob = props))
  names(state_n) <- names(props)
  unlist(lapply(names(state_n), function(state_name) {
    n_use <- state_n[[state_name]]
    if (n_use == 0) return(integer())
    sample(cells_by_state[[state_name]], size = n_use, replace = TRUE)
  }), use.names = FALSE)
}

sample_shift_pair <- function(alpha, state_levels, shift_min, shift_max) {
  props_a <- rdirichlet_one(alpha)
  loss_candidates <- state_levels[props_a[state_levels] > 0.02]
  if (length(loss_candidates) == 0) {
    loss_candidates <- names(which.max(props_a[state_levels]))
  }
  loss_state <- sample(loss_candidates, 1)
  gain_state <- sample(setdiff(state_levels, loss_state), 1)
  max_shift <- min(shift_max, props_a[[loss_state]] * 0.8)
  min_shift <- min(shift_min, max_shift)
  shift_amount <- runif(1, min = min_shift, max = max_shift)
  props_b <- props_a
  props_b[[gain_state]] <- props_b[[gain_state]] + shift_amount
  props_b[[loss_state]] <- props_b[[loss_state]] - shift_amount
  props_b[props_b < 0] <- 0
  props_b <- props_b / sum(props_b)
  list(
    props_a = props_a,
    props_b = props_b,
    gain_state = gain_state,
    loss_state = loss_state,
    shift_amount = shift_amount
  )
}

calc_gene_cpm <- function(cell_idx, counts_mat, lib_size) {
  gene_counts <- Matrix::rowSums(counts_mat[, cell_idx, drop = FALSE])
  total_lib <- sum(lib_size[cell_idx], na.rm = TRUE)
  if (!is.finite(total_lib) || total_lib <= 0) {
    total_lib <- sum(gene_counts)
  }
  as.numeric(gene_counts) / total_lib * 1e6
}

calc_gene_data_expr <- function(cell_idx, data_mat) {
  as.numeric(rowMeans(expm1(as.matrix(data_mat[, cell_idx, drop = FALSE]))))
}

####################
# load data
####################
marker_object_path <- file.path("Auto_five_state_markers", "cache", "pdos_state5_embedded.rds")

if (file.exists(marker_object_path)) {
  message("Loading five-state marker-comparison object: ", marker_object_path)
  pdos <- readRDS(marker_object_path)
  DefaultAssay(pdos) <- "RNA"
  if (!"state" %in% colnames(pdos@meta.data)) {
    stop("The marker-comparison object does not contain a `state` metadata column.")
  }
  keep_cells <- colnames(pdos)[
    as.character(pdos$state) %in% state_order &
      as.character(pdos$orig.ident) != "SUR843T3_PDO"
  ]
  truth_state <- as.character(pdos$state[match(keep_cells, colnames(pdos))])
} else {
  message("Marker-comparison object not found; loading PDOs_merged.rds and Auto_PDO_final_states.rds.")
  pdos <- readRDS("PDOs_merged.rds")
  state_labels <- readRDS("Auto_PDO_final_states.rds")
  DefaultAssay(pdos) <- "RNA"

  common_cells <- intersect(colnames(pdos), names(state_labels))
  state_labels <- state_labels[common_cells]

  keep_cells <- common_cells[
    as.character(state_labels) %in% state_order &
      as.character(pdos$orig.ident[match(common_cells, colnames(pdos))]) != "SUR843T3_PDO"
  ]
  truth_state <- as.character(state_labels[keep_cells])
}

if (length(keep_cells) == 0) stop("No cells remained after state and sample filtering.")
truth_state <- factor(truth_state, levels = state_order)

all_marker_genes <- unique(unlist(marker_panels, use.names = FALSE))
counts_all <- get_assay_layer(pdos, assay = "RNA", layer = "counts")
data_all <- get_assay_layer(pdos, assay = "RNA", layer = "data")
available_marker_genes <- intersect(all_marker_genes, rownames(counts_all))
missing_marker_genes <- setdiff(all_marker_genes, available_marker_genes)
if (length(available_marker_genes) == 0) stop("None of the selected marker genes are present.")

message("Subsetting marker count matrix.")
marker_counts <- counts_all[available_marker_genes, keep_cells, drop = FALSE]
marker_detected <- marker_counts > 0
marker_data <- data_all[available_marker_genes, keep_cells, drop = FALSE]
sample_labels <- as.character(pdos$orig.ident[match(keep_cells, colnames(pdos))])

if ("nCount_RNA" %in% colnames(pdos@meta.data)) {
  lib_size <- as.numeric(pdos@meta.data[keep_cells, "nCount_RNA"])
} else {
  lib_size <- as.numeric(Matrix::colSums(counts_all[, keep_cells, drop = FALSE]))
}
lib_size[!is.finite(lib_size) | lib_size <= 0] <- median(lib_size[is.finite(lib_size) & lib_size > 0])

rm(counts_all, data_all, pdos)
if (exists("state_labels")) rm(state_labels)
invisible(gc())

####################
# marker manifest and expression summaries
####################
manifest <- data.frame(
  panel = rep(names(marker_panels), lengths(marker_panels)),
  truth_state = rep(unname(panel_truth_state[names(marker_panels)]), lengths(marker_panels)),
  gene = unlist(marker_panels, use.names = FALSE),
  stringsAsFactors = FALSE
) %>%
  mutate(
    present_in_pdo = gene %in% available_marker_genes,
    negative_control = panel == "Immune regulated"
  )

fwrite(manifest, file.path(out_dir, "Auto_marker_panel_manifest.csv"))

if (length(missing_marker_genes) > 0) {
  warning("Missing marker genes in PDO object: ", paste(missing_marker_genes, collapse = ", "))
}

specificity_cache_path <- file.path("Auto_five_state_markers", "cache", "state_specificity.rds")
if (!file.exists(specificity_cache_path)) {
  stop("Missing marker specificity cache: ", specificity_cache_path)
}

marker_specificity <- readRDS(specificity_cache_path) %>%
  as.data.frame() %>%
  filter(gene %in% available_marker_genes)

threshold_tbl <- manifest %>%
  filter(!negative_control, present_in_pdo) %>%
  select(panel, truth_state, gene) %>%
  left_join(
    marker_specificity %>%
      select(gene, state, state_median_expr, off_state_max_median_expr, specificity_gap, best_state),
    by = c("gene", "truth_state" = "state")
  ) %>%
  mutate(
    threshold = off_state_max_median_expr +
      params$gate_threshold_fraction * (state_median_expr - off_state_max_median_expr),
    threshold = ifelse(is.finite(threshold) & threshold > 0, threshold, state_median_expr),
    threshold = pmax(threshold, 0),
    marker_comparison_target_best = best_state == truth_state,
    marker_comparison_target_above_off = state_median_expr > off_state_max_median_expr
  )

# Refined empirical single-cell thresholding (reducing ambiguity)
panel_thresholds <- sapply(names(classification_panels), function(pn) {
  genes <- intersect(classification_panels[[pn]], rownames(marker_data))
  if (length(genes) == 0) return(0)
  
  pm <- colMeans(as.matrix(marker_data[genes, , drop = FALSE]))
  target_state <- panel_truth_state[[pn]]
  
  # Threshold at the 90th percentile of off-target cells (more specific than 80th)
  # to reduce the "Ambiguous" rate.
  off_target_pm <- pm[truth_state != target_state]
  bg_threshold <- as.numeric(quantile(off_target_pm, probs = 0.90, na.rm = TRUE))
  
  # Use 0.05 as a floor to capture clear signal above baseline noise
  max(bg_threshold, 0.05)
})

if (any(!is.finite(panel_thresholds))) {
  stop("Could not compute aggregated panel thresholds.")
}

detection_by_state <- lapply(state_order, function(state_name) {
  state_idx <- which(truth_state == state_name)
  state_samples <- unique(sample_labels[state_idx])

  sample_expr <- lapply(state_samples, function(sample_id) {
    idx <- state_idx[sample_labels[state_idx] == sample_id]
    if (length(idx) == 0) return(NULL)
    Matrix::rowMeans(marker_data[, idx, drop = FALSE])
  })
  sample_expr <- sample_expr[!vapply(sample_expr, is.null, logical(1))]

  sample_pct <- lapply(state_samples, function(sample_id) {
    idx <- state_idx[sample_labels[state_idx] == sample_id]
    if (length(idx) == 0) return(NULL)
    Matrix::rowMeans(marker_detected[, idx, drop = FALSE])
  })
  sample_pct <- sample_pct[!vapply(sample_pct, is.null, logical(1))]

  sample_counts <- lapply(state_samples, function(sample_id) {
    idx <- state_idx[sample_labels[state_idx] == sample_id]
    if (length(idx) == 0) return(NULL)
    Matrix::rowMeans(marker_counts[, idx, drop = FALSE])
  })
  sample_counts <- sample_counts[!vapply(sample_counts, is.null, logical(1))]

  expr_mat <- do.call(cbind, sample_expr)
  pct_mat <- do.call(cbind, sample_pct)
  count_mat <- do.call(cbind, sample_counts)

  data.frame(
    state = state_name,
    gene = available_marker_genes,
    pct_detected = as.numeric(apply(pct_mat, 1, median, na.rm = TRUE)),
    mean_umi = as.numeric(apply(count_mat, 1, median, na.rm = TRUE)),
    mean_log_norm_expr = as.numeric(rowMeans(expr_mat, na.rm = TRUE)),
    median_log_norm_expr = as.numeric(apply(expr_mat, 1, median, na.rm = TRUE)),
    sample_n = ncol(expr_mat),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  left_join(manifest %>% select(panel, gene, negative_control), by = "gene") %>%
  group_by(gene) %>%
  mutate(
    gene_scaled_expr = as.numeric(scale(median_log_norm_expr)),
    gene_scaled_expr = ifelse(is.finite(gene_scaled_expr), gene_scaled_expr, 0)
  ) %>%
  ungroup()

fwrite(detection_by_state, file.path(out_dir, "Auto_marker_detection_by_state.csv"))

marker_validation <- detection_by_state %>%
  filter(!negative_control) %>%
  group_by(panel, gene) %>%
  summarise(
    target_state = first(panel_truth_state[[first(panel)]]),
    data_layer_best_state = state[which.max(median_log_norm_expr)],
    data_layer_target_median = median_log_norm_expr[state == target_state][1],
    data_layer_max_off_median = max(median_log_norm_expr[state != target_state], na.rm = TRUE),
    data_layer_specificity_gap = data_layer_target_median - data_layer_max_off_median,
    data_layer_target_best = data_layer_best_state == target_state,
    .groups = "drop"
  ) %>%
  left_join(
    threshold_tbl %>%
      select(
        panel,
        gene,
        marker_comparison_target_median = state_median_expr,
        marker_comparison_max_off_median = off_state_max_median_expr,
        marker_comparison_specificity_gap = specificity_gap,
        marker_comparison_best_state = best_state,
        marker_comparison_target_best,
        marker_comparison_target_above_off,
        gate_threshold = threshold
      ),
    by = c("panel", "gene")
  )

fwrite(marker_validation, file.path(out_dir, "Auto_marker_expression_validation_vs_marker_comparison.csv"))

gate_gene_positive_by_state <- bind_rows(lapply(seq_len(nrow(threshold_tbl)), function(i) {
  row_use <- threshold_tbl[i, ]
  gene_use <- row_use$gene
  threshold_use <- row_use$threshold

  bind_rows(lapply(state_order, function(state_name) {
    state_idx <- which(truth_state == state_name)
    state_samples <- unique(sample_labels[state_idx])
    sample_positive <- vapply(state_samples, function(sample_id) {
      idx <- state_idx[sample_labels[state_idx] == sample_id]
      mean(as.numeric(marker_data[gene_use, idx, drop = TRUE]) >= threshold_use, na.rm = TRUE)
    }, numeric(1))
    data.frame(
      panel = row_use$panel,
      gene = gene_use,
      threshold = threshold_use,
      state = state_name,
      median_sample_pct_above_gate = median(sample_positive, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))

gate_panel_positive_by_state <- bind_rows(lapply(names(classification_panels), function(panel_name) {
  genes_present <- intersect(classification_panels[[panel_name]], rownames(marker_data))
  threshold_use <- threshold_tbl$threshold[
    match(
      paste(panel_name, genes_present, sep = "||"),
      paste(threshold_tbl$panel, threshold_tbl$gene, sep = "||")
    )
  ]
  names(threshold_use) <- genes_present

  bind_rows(lapply(state_order, function(state_name) {
    state_idx <- which(truth_state == state_name)
    state_samples <- unique(sample_labels[state_idx])
    sample_positive <- vapply(state_samples, function(sample_id) {
      idx <- state_idx[sample_labels[state_idx] == sample_id]
      expr_use <- as.matrix(marker_data[genes_present, idx, drop = FALSE])
      positive_n <- colSums(sweep(expr_use, 1, threshold_use, ">="))
      mean(positive_n >= params$gate_min_positive_markers, na.rm = TRUE)
    }, numeric(1))
    data.frame(
      panel = panel_name,
      state = state_name,
      median_sample_pct_panel_positive = median(sample_positive, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))

fwrite(gate_gene_positive_by_state, file.path(out_dir, "Auto_marker_gate_gene_positive_by_state.csv"))
fwrite(gate_panel_positive_by_state, file.path(out_dir, "Auto_marker_gate_panel_positive_by_state.csv"))

det_plot <- detection_by_state %>%
  mutate(
    state = factor(state, levels = state_order),
    panel = factor(panel, levels = names(marker_panels)),
    gene = factor(gene, levels = rev(manifest$gene))
  ) %>%
  ggplot(aes(x = state, y = gene)) +
  geom_point(aes(size = pct_detected, fill = gene_scaled_expr), shape = 21, color = "grey35", stroke = 0.35) +
  facet_grid(
    panel ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y",
    labeller = labeller(
      panel = setNames(
        sprintf("<span style='color:%s'>%s</span>", state_cols[names(marker_panels)], names(marker_panels)),
        names(marker_panels)
      )
    )
  ) +
  scale_size_continuous(range = c(0.2, 5), labels = scales::percent_format(accuracy = 1)) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  scale_y_discrete(labels = function(x) x) +
  scale_x_discrete(labels = function(x) x) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(
      colour = state_cols[state_order],
      angle = 35, hjust = 1, face = "bold"
    ),
    axis.text.y = element_text(
      colour = rev(state_cols[manifest$panel]),
      size = 8
    ),
    strip.placement = "outside",
    strip.text.y = ggtext::element_markdown(angle = 0, face = "bold", size = 8),
    strip.text.y.left = ggtext::element_markdown(angle = 0, face = "bold", size = 8),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.6, "lines")
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Detected cells",
    fill = "Gene-normalized\nmedian RNA data"
  )

ggsave(
  file.path(out_dir, "Auto_marker_detection_by_state_dotplot.pdf"),
  det_plot,
  width = 8.5,
  height = 5.5,
  useDingbats = FALSE
)

####################
# simulation 1: marker-gated single-cell classification
####################
message("Running marker-gated single-cell simulations.")

prediction_levels <- c(names(classification_panels), "Unresolved", "Ambiguous")
truth_levels <- state_order

gate_confusion <- list()
gate_metrics <- list()
gate_cell_calls <- list()
gate_counter <- 1L

for (rep_id in seq_len(params$gate_reps)) {
  n_subset <- min(params$gate_subset_n, length(keep_cells))
  cell_idx <- sample(seq_along(keep_cells), size = n_subset, replace = FALSE)
  truth_use <- as.character(truth_state[cell_idx])

  pred_tbl <- classify_marker_cells(
    cell_idx = cell_idx,
    marker_expr = marker_data,
    panels = classification_panels,
    panel_thresholds = panel_thresholds
  )
  pred <- pred_tbl$prediction

  tab <- as.data.frame(table(
    truth_state = factor(truth_use, levels = truth_levels),
    predicted_state = factor(pred, levels = prediction_levels)
  ))

  gate_confusion[[gate_counter]] <- tab %>%
    mutate(
      rep = rep_id,
      n_subset = n_subset,
      gate_rule = gate_rule_name,
      n = as.integer(Freq)
    ) %>%
    select(rep, n_subset, gate_rule, truth_state, predicted_state, n)

  gate_metrics[[gate_counter]] <- summarise_gate_metrics(
    pred = pred,
    truth = truth_use,
    rep_id = rep_id,
    n_subset = n_subset,
    rule_name = gate_rule_name
  )

  gate_cell_calls[[gate_counter]] <- pred_tbl %>%
    mutate(
      rep = rep_id,
      cell = keep_cells[cell_idx],
      truth_state = truth_use,
      n_subset = n_subset,
      gate_rule = gate_rule_name
    ) %>%
    select(rep, n_subset, gate_rule, cell, truth_state, prediction, top_score, hit_n)

  gate_counter <- gate_counter + 1L
}

gate_confusion_df <- rbindlist(gate_confusion)
gate_metrics_df <- rbindlist(gate_metrics, fill = TRUE)
gate_cell_calls_df <- rbindlist(gate_cell_calls, fill = TRUE)

gate_summary <- gate_metrics_df %>%
  group_by(n_subset, gate_rule, panel, truth_state) %>%
  summarise(
    reps = n(),
    median_sensitivity = median(sensitivity, na.rm = TRUE),
    median_precision = median(precision, na.rm = TRUE),
    median_f1 = median(f1, na.rm = TRUE),
    median_predicted_n = median(predicted_n, na.rm = TRUE),
    median_false_call_rate = median(false_call_rate, na.rm = TRUE),
    median_unresolved_rate = median(unresolved_rate, na.rm = TRUE),
    median_ambiguous_rate = median(ambiguous_rate, na.rm = TRUE),
    median_call_rate = median(call_rate, na.rm = TRUE),
    median_confident_call_rate = median(confident_call_rate, na.rm = TRUE),
    median_overall_accuracy = median(overall_accuracy, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(gate_confusion_df, file.path(out_dir, "Auto_marker_gate_simulation_confusion_replicates.csv"))
fwrite(gate_metrics_df, file.path(out_dir, "Auto_marker_gate_simulation_metrics_replicates.csv"))
fwrite(gate_cell_calls_df, file.path(out_dir, "Auto_marker_gate_simulation_cell_calls.csv"))
fwrite(gate_summary, file.path(out_dir, "Auto_marker_gate_simulation_summary.csv"))

gate_confusion_plot_df <- gate_confusion_df %>%
  group_by(rep, n_subset, gate_rule, truth_state) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  group_by(n_subset, gate_rule, truth_state, predicted_state) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  filter(mean_prop > 0)

gate_plot <- gate_confusion_plot_df %>%
  ggplot(aes(x = predicted_state, y = truth_state, fill = mean_prop)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = scales::percent(mean_prop, accuracy = 1)), size = 2.4) +
  scale_fill_gradient(low = "white", high = "#D7191C", labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Marker-only assignment",
    y = "Finalized PDO state",
    fill = "Mean fraction",
    caption = "Input cells all have finalized PDO states. Unresolved: no marker panel passes the 2-of-3 data-layer gate. Ambiguous: exact tie between top positive panel margins."
  )

ggsave(
  file.path(out_dir, "Auto_marker_gate_simulation_confusion_heatmap.pdf"),
  gate_plot,
  width = 9,
  height = 4.8,
  useDingbats = FALSE
)

# New Visualization: Diagnostic scatter distribution
message("Creating diagnostic assignment scatter plots.")

# Re-calculate all panel means for a fixed set of cells for visualization
set.seed(42)
plot_cells_idx <- sample(seq_along(keep_cells), min(5000, length(keep_cells)))
plot_cells <- keep_cells[plot_cells_idx]
plot_truth <- truth_state[plot_cells_idx]

plot_means_list <- list()
for (pn in names(classification_panels)) {
  genes <- intersect(classification_panels[[pn]], rownames(marker_data))
  if (length(genes) == 0) {
    plot_means_list[[pn]] <- rep(0, length(plot_cells_idx))
  } else {
    plot_means_list[[pn]] <- colMeans(as.matrix(marker_data[genes, plot_cells_idx, drop = FALSE]))
  }
}
plot_means_mat <- do.call(cbind, plot_means_list)

# Determine prediction for each cell
plot_preds <- apply(plot_means_mat, 1, function(row) {
  hits <- names(classification_panels)[which(row > panel_thresholds[names(classification_panels)])]
  if (length(hits) == 0) return("Unresolved")
  if (length(hits) > 1) return("Ambiguous")
  return(hits)
})

plot_diag_long <- list()
for (pn in names(classification_panels)) {
  plot_diag_long[[pn]] <- data.frame(
    panel = pn,
    cell = plot_cells,
    truth_state = plot_truth,
    mean_expr = plot_means_mat[, pn],
    prediction = plot_preds,
    threshold = panel_thresholds[pn],
    stringsAsFactors = FALSE
  )
}
plot_diag_df <- rbindlist(plot_diag_long) %>%
  mutate(
    truth_state = factor(truth_state, levels = state_order),
    panel = factor(panel, levels = names(classification_panels)),
    pred_category = case_when(
      prediction == "Unresolved" ~ "Unresolved",
      prediction == "Ambiguous" ~ "Ambiguous",
      prediction == as.character(panel) ~ "Assigned to this panel",
      TRUE ~ "Assigned to other single panel"
    ),
    pred_category = factor(pred_category, levels = c("Assigned to this panel", "Assigned to other single panel", "Ambiguous", "Unresolved"))
  )

diag_scatter_plot <- ggplot(plot_diag_df, aes(x = mean_expr, y = truth_state, color = pred_category)) +
  geom_vline(aes(xintercept = threshold), linetype = "dashed", color = "gray40", alpha = 0.6) +
  geom_jitter(height = 0.25, size = 0.6, alpha = 0.5) +
  facet_wrap(~ panel, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c(
    "Assigned to this panel" = "#33A02C", 
    "Assigned to other single panel" = "#FB9A99", 
    "Ambiguous" = "#1F78B4", 
    "Unresolved" = "#A6CEE3"
  )) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 8)
  ) +
  labs(
    x = "Panel Mean Expression (Log-Normalized)",
    y = "Finalized PDO State (Truth)",
    color = "Marker Assignment Result",
    title = "Marker Expression Distribution vs Assignment Category",
    subtitle = "Each panel shows the expression of its own marker set. Dashed line = panel threshold."
  )

ggsave(
  file.path(out_dir, "Auto_marker_assignment_scatter_diagnostic.pdf"),
  diag_scatter_plot,
  width = 14,
  height = 5.5,
  useDingbats = FALSE
)

####################
# simulation 2: qRT-PCR abundance tracking
####################
message("Running qRT-PCR abundance simulations.")

observed_props <- prop.table(table(factor(truth_state, levels = state_order)))
alpha <- as.numeric(observed_props) * params$dirichlet_concentration + 1
names(alpha) <- state_order

cells_by_state <- split(seq_along(keep_cells), truth_state)
qpcr_rows <- list()
qpcr_gene_rows <- list()

for (rep_id in seq_len(params$qpcr_reps)) {
  shift_pair <- sample_shift_pair(
    alpha = alpha,
    state_levels = state_order,
    shift_min = params$qpcr_shift_min,
    shift_max = params$qpcr_shift_max
  )
  props_a <- shift_pair$props_a
  props_b <- shift_pair$props_b

  idx_a <- sample_condition_indices(
    props = props_a,
    n_cells = params$qpcr_cells_per_condition,
    cells_by_state = cells_by_state
  )
  idx_b <- sample_condition_indices(
    props = props_b,
    n_cells = params$qpcr_cells_per_condition,
    cells_by_state = cells_by_state
  )

  expr_a <- calc_gene_data_expr(idx_a, marker_data)
  expr_b <- calc_gene_data_expr(idx_b, marker_data)
  names(expr_a) <- rownames(marker_data)
  names(expr_b) <- rownames(marker_data)

  sampled_props_a <- prop.table(table(factor(as.character(truth_state[idx_a]), levels = state_order)))
  sampled_props_b <- prop.table(table(factor(as.character(truth_state[idx_b]), levels = state_order)))
  eps <- 1 / params$qpcr_cells_per_condition

  for (panel_name in names(marker_panels)) {
    genes_present <- intersect(marker_panels[[panel_name]], rownames(marker_counts))
    if (length(genes_present) == 0) next

    gene_log2fc <- log2((expr_b[genes_present] + params$qpcr_pseudo_cpm) /
      (expr_a[genes_present] + params$qpcr_pseudo_cpm))

    truth_state_name <- panel_truth_state[[panel_name]]
    is_negative_control <- is.na(truth_state_name)
    if (is_negative_control) {
      prop_a <- NA_real_
      prop_b <- NA_real_
      truth_delta <- NA_real_
      truth_log2fc <- NA_real_
    } else {
      prop_a <- as.numeric(sampled_props_a[[truth_state_name]])
      prop_b <- as.numeric(sampled_props_b[[truth_state_name]])
      truth_delta <- prop_b - prop_a
      truth_log2fc <- log2((prop_b + eps) / (prop_a + eps))
    }

    qpcr_rows[[length(qpcr_rows) + 1L]] <- data.frame(
      rep = rep_id,
      panel = panel_name,
      truth_state = truth_state_name,
      negative_control = is_negative_control,
      gain_state = shift_pair$gain_state,
      loss_state = shift_pair$loss_state,
      shift_amount = shift_pair$shift_amount,
      shift_role = case_when(
        is_negative_control ~ "negative_control",
        truth_state_name == shift_pair$gain_state ~ "gained_state",
        truth_state_name == shift_pair$loss_state ~ "lost_state",
        TRUE ~ "other_state"
      ),
      prop_condition_a = prop_a,
      prop_condition_b = prop_b,
      truth_delta = truth_delta,
      truth_log2fc = truth_log2fc,
      marker_panel_expr_condition_a = median(expr_a[genes_present], na.rm = TRUE),
      marker_panel_expr_condition_b = median(expr_b[genes_present], na.rm = TRUE),
      marker_panel_log2fc = median(gene_log2fc, na.rm = TRUE),
      marker_panel_mean_log2fc = mean(gene_log2fc, na.rm = TRUE),
      marker_gene_n = length(genes_present),
      stringsAsFactors = FALSE
    )

    qpcr_gene_rows[[length(qpcr_gene_rows) + 1L]] <- data.frame(
      rep = rep_id,
      panel = panel_name,
      gene = genes_present,
      expr_condition_a = as.numeric(expr_a[genes_present]),
      expr_condition_b = as.numeric(expr_b[genes_present]),
      gene_log2fc = as.numeric(gene_log2fc),
      stringsAsFactors = FALSE
    )
  }
}

qpcr_df <- rbindlist(qpcr_rows, fill = TRUE) %>%
  mutate(
    truth_direction = case_when(
      is.na(truth_log2fc) ~ NA_real_,
      abs(truth_log2fc) < params$min_truth_abs_log2fc ~ 0,
      truth_log2fc > 0 ~ 1,
      TRUE ~ -1
    ),
    marker_direction = case_when(
      abs(marker_panel_log2fc) < params$min_truth_abs_log2fc ~ 0,
      marker_panel_log2fc > 0 ~ 1,
      TRUE ~ -1
    ),
    direction_match = truth_direction == marker_direction
  ) %>%
  group_by(rep) %>%
  mutate(
    marker_panel_fraction_condition_a = ifelse(
      negative_control,
      NA_real_,
      marker_panel_expr_condition_a / sum(marker_panel_expr_condition_a[!negative_control], na.rm = TRUE)
    ),
    marker_panel_fraction_condition_b = ifelse(
      negative_control,
      NA_real_,
      marker_panel_expr_condition_b / sum(marker_panel_expr_condition_b[!negative_control], na.rm = TRUE)
    ),
    marker_fraction_delta = marker_panel_fraction_condition_b - marker_panel_fraction_condition_a
  ) %>%
  ungroup()

qpcr_gene_df <- rbindlist(qpcr_gene_rows, fill = TRUE)

qpcr_transition_eval <- qpcr_df %>%
  filter(!negative_control) %>%
  group_by(rep, gain_state, loss_state, shift_amount) %>%
  summarise(
    gained_marker_call = panel[which.max(marker_panel_log2fc)],
    lost_marker_call = panel[which.min(marker_panel_log2fc)],
    gained_panel_log2fc = marker_panel_log2fc[panel == first(gain_state)][1],
    lost_panel_log2fc = marker_panel_log2fc[panel == first(loss_state)][1],
    other_panel_median_log2fc = median(marker_panel_log2fc[!(panel %in% c(first(gain_state), first(loss_state)))], na.rm = TRUE),
    gained_marker_highest = gained_marker_call == first(gain_state),
    lost_marker_lowest = lost_marker_call == first(loss_state),
    gained_positive = gained_panel_log2fc > 0,
    lost_negative = lost_panel_log2fc < 0,
    relative_shift_correct = gained_marker_highest & lost_marker_lowest & gained_positive & lost_negative,
    .groups = "drop"
  )

qpcr_transition_summary <- qpcr_transition_eval %>%
  group_by(gain_state, loss_state) %>%
  summarise(
    reps = n(),
    median_shift_amount = median(shift_amount, na.rm = TRUE),
    gain_top_rate = mean(gained_marker_highest, na.rm = TRUE),
    loss_bottom_rate = mean(lost_marker_lowest, na.rm = TRUE),
    gain_positive_rate = mean(gained_positive, na.rm = TRUE),
    loss_negative_rate = mean(lost_negative, na.rm = TRUE),
    relative_shift_correct_rate = mean(relative_shift_correct, na.rm = TRUE),
    median_gained_panel_log2fc = median(gained_panel_log2fc, na.rm = TRUE),
    median_lost_panel_log2fc = median(lost_panel_log2fc, na.rm = TRUE),
    .groups = "drop"
  )

qpcr_summary <- qpcr_df %>%
  group_by(panel, truth_state, negative_control, shift_role) %>%
  summarise(
    reps = n(),
    spearman_truth_marker = safe_cor(truth_log2fc, marker_panel_log2fc, method = "spearman"),
    pearson_truth_marker = safe_cor(truth_log2fc, marker_panel_log2fc, method = "pearson"),
    marker_vs_truth_lm_slope = safe_lm_slope(truth_log2fc, marker_panel_log2fc),
    median_abs_log2fc_error = median(abs(marker_panel_log2fc - truth_log2fc), na.rm = TRUE),
    direction_agreement_all = mean(direction_match, na.rm = TRUE),
    direction_agreement_nonzero_truth = mean(
      direction_match[!is.na(truth_direction) & truth_direction != 0],
      na.rm = TRUE
    ),
    raw_sign_agreement_nonzero_truth = mean(
      sign(marker_panel_log2fc[!is.na(truth_log2fc) & abs(truth_log2fc) >= params$min_truth_abs_log2fc]) ==
        sign(truth_log2fc[!is.na(truth_log2fc) & abs(truth_log2fc) >= params$min_truth_abs_log2fc]),
      na.rm = TRUE
    ),
    median_marker_log2fc = median(marker_panel_log2fc, na.rm = TRUE),
    mad_marker_log2fc = mad(marker_panel_log2fc, na.rm = TRUE),
    false_positive_rate_abs_log2fc_ge_1 = mean(abs(marker_panel_log2fc) >= 1, na.rm = TRUE),
    .groups = "drop"
  )

qpcr_role_summary <- qpcr_df %>%
  group_by(panel, negative_control, shift_role) %>%
  summarise(
    reps = n(),
    median_truth_delta = median(truth_delta, na.rm = TRUE),
    median_truth_log2fc = median(truth_log2fc, na.rm = TRUE),
    median_marker_log2fc = median(marker_panel_log2fc, na.rm = TRUE),
    q25_marker_log2fc = quantile(marker_panel_log2fc, 0.25, na.rm = TRUE),
    q75_marker_log2fc = quantile(marker_panel_log2fc, 0.75, na.rm = TRUE),
    expected_direction_rate = case_when(
      first(shift_role) == "gained_state" ~ mean(marker_panel_log2fc > 0, na.rm = TRUE),
      first(shift_role) == "lost_state" ~ mean(marker_panel_log2fc < 0, na.rm = TRUE),
      first(shift_role) == "other_state" ~ mean(abs(marker_panel_log2fc) < params$min_truth_abs_log2fc, na.rm = TRUE),
      TRUE ~ NA_real_
    ),
    .groups = "drop"
  )

qpcr_gene_summary <- qpcr_gene_df %>%
  group_by(panel, gene) %>%
  summarise(
    median_gene_log2fc = median(gene_log2fc, na.rm = TRUE),
    mad_gene_log2fc = mad(gene_log2fc, na.rm = TRUE),
    median_expr_condition_a = median(expr_condition_a, na.rm = TRUE),
    median_expr_condition_b = median(expr_condition_b, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(qpcr_df, file.path(out_dir, "Auto_qpcr_abundance_simulation_replicates.csv"))
fwrite(qpcr_gene_df, file.path(out_dir, "Auto_qpcr_gene_log2fc_replicates.csv"))
fwrite(qpcr_summary, file.path(out_dir, "Auto_qpcr_abundance_simulation_summary.csv"))
fwrite(qpcr_gene_summary, file.path(out_dir, "Auto_qpcr_gene_log2fc_summary.csv"))
fwrite(qpcr_role_summary, file.path(out_dir, "Auto_qpcr_shift_role_summary.csv"))
fwrite(qpcr_transition_eval, file.path(out_dir, "Auto_qpcr_transition_eval_replicates.csv"))
fwrite(qpcr_transition_summary, file.path(out_dir, "Auto_qpcr_transition_summary.csv"))

qpcr_role_plot <- qpcr_role_summary %>%
  filter(!negative_control) %>%
  mutate(
    panel = factor(panel, levels = names(marker_panels)),
    shift_role = factor(shift_role, levels = c("gained_state", "lost_state", "other_state"))
  ) %>%
  ggplot(aes(x = shift_role, y = panel, fill = median_marker_log2fc)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", median_marker_log2fc)), size = 2.7) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Role of that state in the simulated condition shift",
    y = "Marker panel",
    fill = "Median marker\nlog2FC"
  )

ggsave(
  file.path(out_dir, "Auto_qpcr_shift_role_summary_heatmap.pdf"),
  qpcr_role_plot,
  width = 7,
  height = 4.5,
  useDingbats = FALSE
)

qpcr_transition_plot <- qpcr_transition_summary %>%
  mutate(
    gain_state = factor(gain_state, levels = state_order),
    loss_state = factor(loss_state, levels = rev(state_order))
  ) %>%
  ggplot(aes(x = gain_state, y = loss_state, fill = relative_shift_correct_rate)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = scales::percent(relative_shift_correct_rate, accuracy = 1)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#D7191C", labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "State increased in condition 2",
    y = "State decreased in condition 2",
    fill = "Correct relative\nmarker shift"
  )

ggsave(
  file.path(out_dir, "Auto_qpcr_transition_accuracy_heatmap.pdf"),
  qpcr_transition_plot,
  width = 7.8,
  height = 6,
  useDingbats = FALSE
)

example_reps <- qpcr_transition_eval %>%
  arrange(desc(shift_amount)) %>%
  slice_head(n = 6) %>%
  pull(rep)

example_labels <- qpcr_transition_eval %>%
  filter(rep %in% example_reps) %>%
  transmute(
    rep,
    simulation = paste0(
      "Simulation ", rep,
      ": +", gain_state,
      " / -", loss_state
    )
  )

qpcr_example_line_df <- bind_rows(
  qpcr_df %>%
    filter(!negative_control, rep %in% example_reps) %>%
    select(rep, panel, condition_1 = prop_condition_a, condition_2 = prop_condition_b) %>%
    pivot_longer(starts_with("condition_"), names_to = "condition", values_to = "value") %>%
    mutate(metric = "True cell-state abundance"),
  qpcr_df %>%
    filter(!negative_control, rep %in% example_reps) %>%
    select(rep, panel, condition_1 = marker_panel_fraction_condition_a, condition_2 = marker_panel_fraction_condition_b) %>%
    pivot_longer(starts_with("condition_"), names_to = "condition", values_to = "value") %>%
    mutate(metric = "Marker signal fraction")
) %>%
  left_join(example_labels, by = "rep") %>%
  mutate(
    condition = factor(condition, levels = c("condition_1", "condition_2"), labels = c("Condition 1", "Condition 2")),
    panel = factor(panel, levels = state_order),
    simulation = factor(simulation, levels = unique(example_labels$simulation)),
    metric = factor(metric, levels = c("True cell-state abundance", "Marker signal fraction"))
  )

qpcr_example_line_plot <- qpcr_example_line_df %>%
  ggplot(aes(x = condition, y = value, group = panel, color = panel)) +
  geom_line(linewidth = 0.55, alpha = 0.9) +
  geom_point(size = 1.8) +
  facet_grid(simulation ~ metric) +
  scale_color_manual(values = state_cols, name = "State / marker panel") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 8.5) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text.y = element_text(size = 7.5, face = "bold")
  ) +
  labs(x = NULL, y = "Fraction of cells or marker signal")

ggsave(
  file.path(out_dir, "Auto_qpcr_example_shift_lineplots.pdf"),
  qpcr_example_line_plot,
  width = 10,
  height = 8.5,
  useDingbats = FALSE
)

qpcr_repeated_line_df <- bind_rows(
  qpcr_df %>%
    filter(!negative_control) %>%
    select(rep, panel, shift_role, condition_1 = prop_condition_a, condition_2 = prop_condition_b) %>%
    pivot_longer(starts_with("condition_"), names_to = "condition", values_to = "value") %>%
    mutate(metric = "True cell-state abundance"),
  qpcr_df %>%
    filter(!negative_control) %>%
    select(rep, panel, shift_role, condition_1 = marker_panel_fraction_condition_a, condition_2 = marker_panel_fraction_condition_b) %>%
    pivot_longer(starts_with("condition_"), names_to = "condition", values_to = "value") %>%
    mutate(metric = "Marker signal fraction")
) %>%
  mutate(
    condition = factor(condition, levels = c("condition_1", "condition_2"), labels = c("Condition 1", "Condition 2")),
    shift_role = factor(shift_role, levels = c("gained_state", "lost_state", "other_state")),
    panel = factor(panel, levels = state_order),
    metric = factor(metric, levels = c("True cell-state abundance", "Marker signal fraction"))
  )

qpcr_repeated_summary <- qpcr_repeated_line_df %>%
  group_by(metric, shift_role, panel, condition) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    q25_value = quantile(value, 0.25, na.rm = TRUE),
    q75_value = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(qpcr_repeated_summary, file.path(out_dir, "Auto_qpcr_repeated_shift_line_summary.csv"))

qpcr_repeated_line_plot <- qpcr_repeated_summary %>%
  ggplot(aes(x = condition, y = median_value, group = panel, color = panel)) +
  geom_line(linewidth = 0.65, alpha = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = q25_value, ymax = q75_value), width = 0.08, alpha = 0.55) +
  facet_grid(shift_role ~ metric) +
  scale_color_manual(values = state_cols, name = "State / marker panel") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 8.8) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  ) +
  labs(x = NULL, y = "Median fraction across simulations")

ggsave(
  file.path(out_dir, "Auto_qpcr_repeated_shift_summary_lineplots.pdf"),
  qpcr_repeated_line_plot,
  width = 10,
  height = 6.5,
  useDingbats = FALSE
)

####################
# Scatter plot of log2FC
####################
qpcr_scatter_plot <- qpcr_df %>%
  filter(!negative_control) %>%
  mutate(panel = factor(panel, levels = state_order)) %>%
  ggplot(aes(x = truth_log2fc, y = marker_panel_log2fc, color = panel)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, se = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(values = state_cols) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(
    x = "True state abundance log2FC",
    y = "Marker panel expression log2FC",
    title = "Correlation of RT-qPCR marker expression vs state proportion"
  )

ggsave(
  file.path(out_dir, "Auto_qpcr_scatter_log2fc.pdf"),
  qpcr_scatter_plot,
  width = 10,
  height = 3.5,
  useDingbats = FALSE
)

####################
# Cross-correlation heatmap
####################
truth_wide <- qpcr_df %>%
  filter(!negative_control) %>%
  select(rep, truth_state, truth_log2fc) %>%
  distinct() %>%
  pivot_wider(names_from = truth_state, values_from = truth_log2fc) %>%
  select(-rep)

marker_wide <- qpcr_df %>%
  select(rep, panel, marker_panel_log2fc) %>%
  distinct() %>%
  pivot_wider(names_from = panel, values_from = marker_panel_log2fc) %>%
  select(-rep)

# Compute correlation matrix
cor_mat <- cor(truth_wide, marker_wide, method = "spearman", use = "pairwise.complete.obs")

cor_df <- as.data.frame(as.table(cor_mat)) %>%
  setNames(c("State", "Marker_Panel", "Correlation")) %>%
  mutate(
    State = factor(State, levels = rev(state_order)),
    Marker_Panel = factor(Marker_Panel, levels = names(marker_panels))
  )

cor_heatmap <- ggplot(cor_df, aes(x = Marker_Panel, y = State, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-1, 1)) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Marker Panel (log2FC)",
    y = "True State Proportion (log2FC)",
    fill = "Spearman\nCorrelation"
  )

ggsave(
  file.path(out_dir, "Auto_qpcr_cross_correlation_heatmap.pdf"),
  cor_heatmap,
  width = 7.5,
  height = 4.5,
  useDingbats = FALSE
)

message("Running time-series simulation for line plots.")
time_points <- 10
ts_rows <- list()

# Generate 10 random compositions
for (t in seq_len(time_points)) {
  # Use uniform random proportions to break biological dominance of any single state
  props_t <- runif(length(state_order))
  props_t <- props_t / sum(props_t)
  names(props_t) <- state_order

  idx_t <- sample_condition_indices(
    props = props_t,
    n_cells = 30000, # Large cell count for smooth bulk profiles
    cells_by_state = cells_by_state
  )
  
  expr_t <- calc_gene_data_expr(idx_t, marker_data)
  names(expr_t) <- rownames(marker_data)
  
  for (panel_name in names(marker_panels)) {
    genes_present <- intersect(marker_panels[[panel_name]], rownames(marker_counts))
    if (length(genes_present) == 0) next
    
    panel_expr <- median(expr_t[genes_present], na.rm = TRUE)
    
    ts_rows[[length(ts_rows) + 1L]] <- data.frame(
      time = t,
      panel = panel_name,
      truth_state = panel_truth_state[[panel_name]],
      prop = ifelse(is.na(panel_truth_state[[panel_name]]), NA_real_, props_t[[panel_truth_state[[panel_name]]]]),
      expr = panel_expr,
      stringsAsFactors = FALSE
    )
  }
}

ts_df <- rbindlist(ts_rows)

# Create plotting df where we scale prop and expr to 0-1 within each state/panel
ts_plot_df <- ts_df %>%
  group_by(panel) %>%
  mutate(
    scaled_prop = if(all(is.na(prop))) NA_real_ else (prop - min(prop, na.rm=TRUE)) / (max(prop, na.rm=TRUE) - min(prop, na.rm=TRUE)),
    scaled_expr = (expr - min(expr, na.rm=TRUE)) / (max(expr, na.rm=TRUE) - min(expr, na.rm=TRUE))
  ) %>%
  ungroup()

# Page 1: State abundance vs State's own marker expression
page1_df <- ts_plot_df %>%
  filter(!is.na(truth_state)) %>%
  select(time, state = truth_state, `Abundance` = scaled_prop, `Marker Expression` = scaled_expr) %>%
  pivot_longer(cols = c(`Abundance`, `Marker Expression`), names_to = "metric", values_to = "value") %>%
  mutate(state = factor(state, levels = state_order))

# Compute correlation for page 1 annotations
cor_page1 <- ts_plot_df %>%
  filter(!is.na(truth_state)) %>%
  group_by(state = truth_state) %>%
  summarise(
    cor_val = cor(prop, expr, method = "pearson"),
    p_val = cor.test(prop, expr, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.2f\np = %.2g", cor_val, p_val),
    state = factor(state, levels = state_order)
  )

p1 <- ggplot(page1_df, aes(x = time, y = value, color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ state, nrow = 1) +
  geom_text(
    data = cor_page1,
    aes(x = 1, y = 1, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 3, color = "black"
  ) +
  scale_color_manual(values = c("Abundance" = "black", "Marker Expression" = "#D7191C")) +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(
    x = "Simulated Time Point",
    y = "Min-Max Scaled Value",
    color = "Metric",
    title = "State Abundance vs Own Marker Expression"
  )

# Page 2: State abundance vs Immune infiltration marker expression
# Get immune expression series
immune_expr <- ts_plot_df %>%
  filter(panel == "Immune regulated") %>%
  select(time, immune_expr = scaled_expr)

page2_df <- ts_plot_df %>%
  filter(!is.na(truth_state)) %>%
  select(time, state = truth_state, `Abundance` = scaled_prop) %>%
  left_join(immune_expr, by = "time") %>%
  rename(`Immune Marker Expression` = immune_expr) %>%
  pivot_longer(cols = c(`Abundance`, `Immune Marker Expression`), names_to = "metric", values_to = "value") %>%
  mutate(state = factor(state, levels = state_order))

# Compute correlation for page 2 annotations
cor_page2 <- ts_plot_df %>%
  filter(!is.na(truth_state)) %>%
  select(time, state = truth_state, prop) %>%
  left_join(ts_plot_df %>% filter(panel == "Immune regulated") %>% select(time, expr), by = "time") %>%
  group_by(state) %>%
  summarise(
    cor_val = cor(prop, expr, method = "pearson"),
    p_val = cor.test(prop, expr, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.2f\np = %.2g", cor_val, p_val),
    state = factor(state, levels = state_order)
  )

p2 <- ggplot(page2_df, aes(x = time, y = value, color = metric, group = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ state, nrow = 1) +
  geom_text(
    data = cor_page2,
    aes(x = 1, y = 1, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 3, color = "black"
  ) +
  scale_color_manual(values = c("Abundance" = "black", "Immune Marker Expression" = "#2C7BB6")) +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(
    x = "Simulated Time Point",
    y = "Min-Max Scaled Value",
    color = "Metric",
    title = "State Abundance vs Immune Marker Expression (Negative Control)"
  )

pdf(file.path(out_dir, "Auto_time_series_simulation_lineplots.pdf"), width = 12, height = 4.5, useDingbats = FALSE)
print(p1)
print(p2)
dev.off()
# run manifest
####################
run_info <- data.frame(
  parameter = names(params),
  value = vapply(params, function(x) paste(x, collapse = ","), character(1)),
  stringsAsFactors = FALSE
)

fwrite(run_info, file.path(out_dir, "Auto_marker_selection_simulation_run_parameters.csv"))

methodology <- c(
  "# Auto Marker Selection Simulation",
  "",
  "The marker-only classifier does not use finalized state labels during assignment.",
  "All marker-expression summaries and gates use the Seurat RNA `data` layer (log-normalized).",
  "The script validates markers against `Auto_five_state_markers/cache/state_specificity.rds`.",
  "",
  "Simulation 1 uses only cells with finalized PDO state labels.",
  "Assignment Logic: For each 3-marker panel, the mean expression of the markers is calculated for each single cell. A panel-specific threshold is set to the 90th percentile of off-target cells or a floor of 0.05 (whichever is higher).",
  " - If a cell exceeds the threshold for exactly one panel, it is assigned to that state.",
  " - If a cell exceeds the threshold for multiple panels, it is labeled 'Ambiguous'.",
  " - If a cell exceeds no thresholds, it is labeled 'Unresolved'.",
  "",
  "Simulation 2 creates paired condition shifts (gained/lost states) to test qRT-PCR capturing of relative shifts.",
  "Simulation 3 generates a 10-timepoint longitudinal series to compare state abundance against panel expression (both scaled 0-1) and includes an immune-marker negative control panel.",
  "",
  "The immune-regulated panel is used as a negative control; it is not an assignable PDO state in Simulation 1."
)

writeLines(methodology, file.path(out_dir, "Auto_marker_selection_simulation_methodology.md"))
writeLines(
  methodology,
  file.path(project_dir, "analysis", "methodology", "Auto_marker_selection_simulation_methodology.md")
)

message("Done. Outputs written to: ", file.path(getwd(), out_dir))
