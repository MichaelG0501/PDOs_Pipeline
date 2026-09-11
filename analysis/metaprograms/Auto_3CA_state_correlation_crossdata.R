####################
# Analysis registry:
#   Status: terminal comparison
#   Script: analysis/metaprograms/Auto_3CA_state_correlation_crossdata.R
#   Methodology: analysis/methodology/metaprograms/Auto_3CA_state_correlation_crossdata_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     - PDOs_outs/UCell_3CA_MPs.rds
#     - PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     - scRef_Pipeline/ref_outs/UCell_3CA_MPs.rds
#     - scRef centred state_definition/intermediate/centred_refined_noreg_states.rds
#   Outputs:
#     - PDOs_outs/Auto_3CA_state_correlation_crossdata/figures/Auto_3CA_state_correlation_crossdata.pdf
#     - PDOs_outs/Auto_3CA_state_correlation_crossdata/tables/*.csv
#     - PDOs_outs/Auto_3CA_state_correlation_crossdata/logs/Auto_3CA_state_correlation_crossdata_run_summary.txt
#   Downstream use: none; terminal figure and auditable plotting tables
####################

suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("ggrepel")
  library("patchwork")
  library("tidyr")
})

source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_helpers.R")

####################
# Paths and state orders
####################
script_path <- "analysis/metaprograms/Auto_3CA_state_correlation_crossdata.R"
out_dir <- file.path(PDO_LIVE_OUTS, "Auto_3CA_state_correlation_crossdata")
out_paths <- pdo_ensure_output_tiers(out_dir)

pdo_score_file <- file.path(PDO_LIVE_OUTS, "UCell_3CA_MPs.rds")
pdo_state_file <- file.path(PDO_LIVE_OUTS, "centred_mp_refinement", "centred_refined_noreg_states.rds")
sc_ref_outs <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"
sc_score_file <- file.path(sc_ref_outs, "UCell_3CA_MPs.rds")
sc_state_file <- file.path(
  sc_ref_outs, "Metaprogrammes_Results", "centred", "state_definition",
  "intermediate", "centred_refined_noreg_states.rds"
)
input_files <- c(pdo_score_file, pdo_state_file, sc_score_file, sc_state_file)
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
label_threshold <- 0.1
####################

####################
# Normalize and validate named cell-level inputs
####################
normalise_state_vector <- function(state_vec, dataset_name) {
  state_names <- names(state_vec)
  state_vec <- as.character(state_vec)
  names(state_vec) <- state_names
  if (is.null(names(state_vec)) || any(names(state_vec) == "")) {
    stop(dataset_name, " state vector must have cell-barcode names.")
  }
  state_vec
}

normalise_score_matrix <- function(score_obj, dataset_name) {
  score_mat <- as.data.frame(score_obj, check.names = FALSE)
  if (is.null(rownames(score_mat)) || any(rownames(score_mat) == "")) {
    stop(dataset_name, " 3CA score matrix must have cell-barcode row names.")
  }
  numeric_ok <- vapply(score_mat, is.numeric, logical(1))
  if (!all(numeric_ok)) {
    stop(dataset_name, " score matrix contains non-numeric columns: ",
         paste(names(numeric_ok)[!numeric_ok], collapse = ", "))
  }
  score_mat
}

prepare_dataset <- function(score_mat, state_vec, target_states, dataset_name) {
  common_cells <- intersect(rownames(score_mat), names(state_vec))
  if (length(common_cells) == 0) {
    stop("No shared cell barcodes between ", dataset_name, " scores and states.")
  }
  score_mat <- score_mat[common_cells, , drop = FALSE]
  state_vec <- state_vec[common_cells]
  keep <- !is.na(state_vec) & state_vec %in% target_states
  score_mat <- score_mat[keep, , drop = FALSE]
  state_vec <- state_vec[keep]
  state_counts <- table(factor(state_vec, levels = target_states))
  missing_states <- names(state_counts)[state_counts == 0]
  if (length(missing_states) > 0) {
    stop(dataset_name, " has no cells for target state(s): ",
         paste(missing_states, collapse = ", "))
  }
  list(scores = score_mat, states = state_vec, counts = state_counts)
}

clean_3ca_label <- function(x) {
  x <- gsub("^X3CA_mp_", "3CA_", x)
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("^XMP", "3CA_MP", x)
  x <- gsub("^MP", "3CA_MP", x)
  x <- gsub("^3CA_3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x <- gsub(" +", " ", x)
  x
}

cat("Loading state-resolved 3CA score inputs...\n")
pdo_scores <- normalise_score_matrix(readRDS(pdo_score_file), "PDO")
sc_scores <- normalise_score_matrix(readRDS(sc_score_file), "scRef")
pdo_states <- normalise_state_vector(readRDS(pdo_state_file), "PDO")
sc_states <- normalise_state_vector(readRDS(sc_state_file), "scRef")

colnames(pdo_scores) <- clean_3ca_label(colnames(pdo_scores))
colnames(sc_scores) <- clean_3ca_label(colnames(sc_scores))

common_mps <- intersect(colnames(pdo_scores), colnames(sc_scores))
if (length(common_mps) < 3) {
  stop("Fewer than three shared 3CA MPs were found between PDO and scRef scores.")
}
pdo_data <- prepare_dataset(
  pdo_scores[, common_mps, drop = FALSE], pdo_states, pdo_state_order, "PDO"
)
sc_data <- prepare_dataset(
  sc_scores[, common_mps, drop = FALSE], sc_states, sc_state_order, "scRef"
)
####################

####################
# State means and all 25 cross-dataset Spearman comparisons
####################
state_mean_matrix <- function(score_mat, state_vec, state_order) {
  result <- vapply(
    state_order,
    function(state_name) {
      colMeans(score_mat[state_vec == state_name, , drop = FALSE], na.rm = TRUE)
    },
    numeric(ncol(score_mat))
  )
  rownames(result) <- colnames(score_mat)
  colnames(result) <- state_order
  result
}

pdo_means <- state_mean_matrix(pdo_data$scores, pdo_data$states, pdo_state_order)
sc_means <- state_mean_matrix(sc_data$scores, sc_data$states, sc_state_order)
cor_mat <- matrix(
  NA_real_, nrow = length(pdo_state_order), ncol = length(sc_state_order),
  dimnames = list(pdo_state_order, sc_state_order)
)
cor_stats <- list()
scatter_data <- list()
result_idx <- 1L

for (pdo_state in pdo_state_order) {
  for (sc_state in sc_state_order) {
    pair_df <- data.frame(
      mp = common_mps,
      sc_ref_score = sc_means[common_mps, sc_state],
      pdo_score = pdo_means[common_mps, pdo_state],
      pdo_state = pdo_state,
      sc_ref_state = sc_state,
      stringsAsFactors = FALSE
    )
    pair_df <- pair_df[
      is.finite(pair_df$sc_ref_score) & is.finite(pair_df$pdo_score), , drop = FALSE
    ]
    cor_result <- suppressWarnings(cor.test(
      pair_df$sc_ref_score, pair_df$pdo_score, method = "spearman", exact = FALSE
    ))
    cor_mat[pdo_state, sc_state] <- unname(cor_result$estimate)
    cor_stats[[result_idx]] <- data.frame(
      pdo_state = pdo_state,
      sc_ref_state = sc_state,
      spearman_rho = unname(cor_result$estimate),
      p_value = cor_result$p.value,
      n_3ca_mps = nrow(pair_df),
      pdo_cells = unname(pdo_data$counts[pdo_state]),
      sc_ref_cells = unname(sc_data$counts[sc_state]),
      stringsAsFactors = FALSE
    )
    scatter_data[[result_idx]] <- pair_df
    result_idx <- result_idx + 1L
  }
}
cor_stats <- bind_rows(cor_stats)
scatter_df <- bind_rows(scatter_data)

mean_score_df <- bind_rows(
  as.data.frame(pdo_means, check.names = FALSE) %>%
    tibble::rownames_to_column("mp") %>%
    pivot_longer(-mp, names_to = "state", values_to = "mean_score") %>%
    mutate(dataset = "PDO"),
  as.data.frame(sc_means, check.names = FALSE) %>%
    tibble::rownames_to_column("mp") %>%
    pivot_longer(-mp, names_to = "state", values_to = "mean_score") %>%
    mutate(dataset = "scRef")
) %>%
  select(dataset, state, mp, mean_score)

cor_matrix_file <- file.path(out_paths[["tables"]], "Auto_3CA_state_correlation_matrix.csv")
cor_stats_file <- file.path(out_paths[["tables"]], "Auto_3CA_state_correlation_statistics.csv")
mean_scores_file <- file.path(out_paths[["tables"]], "Auto_3CA_state_mean_scores.csv")
scatter_data_file <- file.path(out_paths[["tables"]], "Auto_3CA_state_pair_scatter_data.csv")
write.csv(
  data.frame(pdo_state = rownames(cor_mat), cor_mat, check.names = FALSE),
  cor_matrix_file, row.names = FALSE
)
write.csv(cor_stats, cor_stats_file, row.names = FALSE)
write.csv(mean_score_df, mean_scores_file, row.names = FALSE)
write.csv(scatter_df, scatter_data_file, row.names = FALSE)
####################

####################
# Page 1: cross-state correlation matrix
####################
heatmap_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
colnames(heatmap_df) <- c("pdo_state", "sc_ref_state", "spearman_rho")
heatmap_df$pdo_state <- factor(heatmap_df$pdo_state, levels = rev(pdo_state_order))
heatmap_df$sc_ref_state <- factor(heatmap_df$sc_ref_state, levels = sc_state_order)
heatmap_df$label <- sprintf("%.2f", heatmap_df$spearman_rho)

p_heatmap <- ggplot(
  heatmap_df, aes(x = sc_ref_state, y = pdo_state, fill = spearman_rho)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-1, 1), name = "Spearman\nrho"
  ) +
  coord_fixed() +
  labs(
    title = "Cross-dataset state correlation across 3CA metaprograms",
    x = "scRef / scATLAS centred noreg state",
    y = "PDO centred noreg state"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 18)
  )
####################

####################
# Pages 2-6: five scRef panels for each successive PDO state
####################
axis_max <- max(c(pdo_means, sc_means), na.rm = TRUE) * 1.06
if (!is.finite(axis_max) || axis_max <= 0) {
  stop("Unable to derive a positive common scatter-plot axis limit.")
}

make_scatter_panel <- function(pdo_state, sc_state) {
  plot_df <- scatter_df %>%
    filter(pdo_state == .env$pdo_state, sc_ref_state == .env$sc_state) %>%
    mutate(
      label = ifelse(
        pdo_score >= label_threshold | sc_ref_score >= label_threshold,
        clean_3ca_label(mp), NA_character_
      ),
      status = ifelse(
        pdo_score >= label_threshold | sc_ref_score >= label_threshold,
        "At least one score >= 0.1", "Both scores < 0.1"
      )
    )
  stat_row <- cor_stats %>%
    filter(pdo_state == .env$pdo_state, sc_ref_state == .env$sc_state)

  ggplot(plot_df, aes(x = sc_ref_score, y = pdo_score)) +
    geom_vline(
      xintercept = label_threshold, linetype = "dotted",
      color = "grey35", linewidth = 0.35
    ) +
    geom_hline(
      yintercept = label_threshold, linetype = "dotted",
      color = "grey35", linewidth = 0.35
    ) +
    geom_point(aes(color = status), size = 2.2, alpha = 0.75) +
    scale_color_manual(values = c(
      "At least one score >= 0.1" = "black",
      "Both scores < 0.1" = "grey65"
    )) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = TRUE, color = "red3", fill = "red",
      linetype = "dashed", linewidth = 0.55, alpha = 0.08
    ) +
    geom_text_repel(
      aes(label = label), size = 2.2, max.overlaps = 20,
      min.segment.length = 0, box.padding = 0.2, point.padding = 0.15,
      na.rm = TRUE, seed = 1
    ) +
    annotate(
      "text", x = 0.015 * axis_max, y = 0.985 * axis_max,
      hjust = 0, vjust = 1, size = 3.1, fontface = "bold",
      label = sprintf(
        "Spearman rho = %.2f\np = %.2g",
        stat_row$spearman_rho, stat_row$p_value
      )
    ) +
    coord_fixed(
      xlim = c(0, axis_max), ylim = c(0, axis_max), expand = FALSE, clip = "on"
    ) +
    labs(
      title = sc_state,
      x = "scRef mean 3CA UCell score",
      y = "PDO mean 3CA UCell score"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        color = unname(state_colors[sc_state]), face = "bold", size = 10, hjust = 0.5
      ),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      plot.margin = margin(5, 5, 5, 5)
    )
}

pdf_file <- file.path(
  out_paths[["figures"]], "Auto_3CA_state_correlation_crossdata.pdf"
)
grDevices::pdf(pdf_file, width = 20, height = 7, useDingbats = FALSE)
print(p_heatmap)
for (pdo_state in pdo_state_order) {
  state_panels <- lapply(
    sc_state_order, function(sc_state) make_scatter_panel(pdo_state, sc_state)
  )
  state_page <- wrap_plots(state_panels, nrow = 1) +
    plot_annotation(
      title = paste0("PDO state: ", pdo_state, " — comparison with all scRef states"),
      theme = theme(
        plot.title = element_text(
          color = unname(state_colors[pdo_state]),
          face = "bold", size = 18, hjust = 0.5
        )
      )
    )
  print(state_page)
}
grDevices::dev.off()
####################

####################
# Persistent run record
####################
output_files <- c(
  pdf_file, cor_matrix_file, cor_stats_file, mean_scores_file, scatter_data_file
)
log_file <- pdo_write_run_summary(
  script = script_path,
  out_dir = out_dir,
  inputs = input_files,
  outputs = output_files,
  parameters = list(
    shared_3ca_mps = length(common_mps),
    pdo_target_states = length(pdo_state_order),
    sc_ref_target_states = length(sc_state_order),
    scatter_label_threshold = label_threshold,
    excluded_state_labels = "Unresolved; Hybrid"
  ),
  status = "completed"
)
cat("Completed state-resolved cross-data 3CA correlation.\n")
cat("PDF: ", pdf_file, "\n", sep = "")
cat("Run summary: ", log_file, "\n", sep = "")
####################
