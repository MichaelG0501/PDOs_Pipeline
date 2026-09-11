####################
# Analysis registry:
#   Status: active terminal personalised-regulon analysis
#   Script: analysis/cell_states/Auto_PDO_personalised_state_regulon_rss.R
#   Methodology: analysis/methodology/cell_states/Auto_PDO_personalised_state_regulon_rss_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Projects four prespecified scATLAS SCENIC regulons into all canonical-state
#     PDO cells, then quantifies expected-state activity and RSS separately in
#     each PDO specimen. Treated and untreated specimens remain separate.
#   Inputs:
#     - PDOs_outs/PDOs_merged.rds
#     - PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     - fixed scATLAS-derived regulon target sets embedded below
#   Outputs:
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/intermediate/Auto_projected_regulon_auc.rds
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/tables/Auto_regulon_signature_genes.csv
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/tables/Auto_sample_state_cell_counts.csv
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/tables/Auto_per_sample_regulon_state_activity.csv
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/tables/Auto_per_sample_regulon_specificity.csv
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/tables/Auto_regulon_personalisation_summary.csv
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/figures/Auto_per_sample_regulon_evidence_matrix.{pdf,png}
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/figures/Auto_per_sample_regulon_profiles.pdf
#     - PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/logs/Auto_PDO_personalised_state_regulon_rss_run_summary.txt
#   Downstream use: none; terminal personalised-regulon evidence and source data.
#   Cache/replot behavior:
#     PDO_FORCE_REBUILD=1 recomputes projected per-cell AUCell scores.
#     PDO_REPLOT_ONLY=1 requires and reuses the live projected-score cache.
#   Run command:
#     qsub analysis/cell_states/Auto_PDO_personalised_state_regulon_rss.sh
#   Conda env: dmtcp
####################

library(Seurat)
library(Matrix)
library(AUCell)
library(SCENIC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

source("analysis/shared/Auto_pdo_analysis_config.R")
source("analysis/shared/Auto_pdo_analysis_helpers.R")

script_path <- "analysis/cell_states/Auto_PDO_personalised_state_regulon_rss.R"
merged_path <- pdo_output_path("PDOs_merged.rds")
state_path <- pdo_output_path(PDO_PREFERRED_STATE_VECTOR)
out_dir <- pdo_output_path(
  "final_mp_scenic",
  "Auto_personalised_state_regulons"
)
out_tiers <- pdo_ensure_output_tiers(out_dir)
cache_path <- file.path(
  out_tiers[["intermediate"]],
  "Auto_projected_regulon_auc.rds"
)

pdo_require_files(
  c(merged_path, state_path),
  c("merged PDO Seurat object", "centred PDO state vector")
)

cache_policy <- pdo_cache_policy()
if (cache_policy$force_rebuild && cache_policy$replot_only) {
  stop("PDO_FORCE_REBUILD and PDO_REPLOT_ONLY cannot both be enabled.")
}

n_cores <- pdo_get_env_integer("PDO_REGULON_NCORES", 8L)
min_expected_cells <- as.integer(PDO_THRESHOLDS$marker_min_cells_state)
min_rest_cells <- as.integer(PDO_THRESHOLDS$marker_min_cells_rest)

####################
# Fixed scATLAS regulon definitions
####################
# These are frozen from:
# /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/
# ref_outs/final_mp_scenic/int/2.6_regulons_asGeneSet.Rds
# SNAI1_extended is used because it is the SNAI1 regulon present in the
# scATLAS AUCell result; the four-gene non-extended set was not scored there.
regulon_gene_sets <- list(
  "SNAI1" = c(
    "DLL4", "ID2", "LAPTM5", "MAFA", "MYCL", "MYLIP", "POU3F1",
    "RASD1", "RGS16", "SNAI1", "SOX7"
  ),
  "RXRB" = c(
    "FOXO4", "HOXA1", "IER5L", "LMX1B", "MNX1", "NUAK2", "OTX1",
    "RXRB", "SP6", "STAT5A", "UBTF", "ZNF579"
  ),
  "ZBTB14" = c(
    "FABP5", "FOXD1", "GINS2", "LEF1", "PDX1", "S100A2", "TFAP2C",
    "ZBTB14", "ZFP82", "ZNF324B", "ZNF396"
  ),
  "TP73" = c(
    "ASF1B", "ATAD2", "ATAD5", "BRCA1", "BRCA2", "BRIP1", "CCNE2",
    "CDC25A", "CDC45", "CDT1", "CHAF1A", "CLSPN", "DMC1", "DNMT1",
    "DSCC1", "DTL", "E2F1", "E2F2", "E2F3", "E2F7", "E2F8", "EXO1",
    "EZH2", "FANCI", "GINS1", "HELLS", "MCM2", "MCM5", "MCM6", "MCM7",
    "MMS22L", "MYBL2", "NR2C2", "PAX6", "PBX3", "PCNA", "PKMYT1",
    "PTPRC", "RAD51", "RBBP5", "RBL1", "RFC2", "RFC3", "RGS16",
    "RUNX3", "SIX3", "TP73", "TYMS", "UHRF1", "ZNF367", "ZNF653"
  )
)

regulon_definitions <- data.frame(
  regulon = c("SNAI1", "RXRB", "ZBTB14", "TP73"),
  source_regulon_id = c("SNAI1_extended", "RXRB", "ZBTB14", "TP73"),
  expected_state = c(
    "Stress-adaptive", "Stress-adaptive",
    "Classic proliferation", "Classic proliferation"
  ),
  stringsAsFactors = FALSE
)

if (!identical(names(regulon_gene_sets), regulon_definitions$regulon)) {
  stop("Regulon definitions and gene-set order do not agree.")
}
if (!all(regulon_definitions$expected_state %in% PDO_STATE_ORDER)) {
  stop("A requested expected state is not in the canonical PDO state order.")
}

make_sample_display <- function(x) {
  x <- gsub("_new4samples_PDO$", "", x)
  x <- gsub("_Untreated_PDO$", " (untreated)", x)
  x <- gsub("_Treated_PDO$", " (treated)", x)
  x <- gsub("_PDO$", "", x)
  x
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

calc_one_rss <- function(values, labels, label_order) {
  if (length(values) == 0 || length(unique(labels)) < 2 ||
      !all(label_order %in% labels) || sum(values, na.rm = TRUE) <= 0) {
    return(setNames(rep(NA_real_, length(label_order)), label_order))
  }
  auc_one <- matrix(values, nrow = 1L)
  rownames(auc_one) <- "selected_regulon"
  rss <- SCENIC::calcRSS(
    AUC = auc_one,
    cellAnnotation = labels,
    cellTypes = label_order
  )
  setNames(as.numeric(rss[1, label_order]), label_order)
}

####################
# Project the fixed regulons into all canonical-state PDO cells
####################
if (cache_policy$replot_only && !file.exists(cache_path)) {
  stop("PDO_REPLOT_ONLY=1 but projected AUCell cache is missing: ", cache_path)
}

if (file.exists(cache_path) && !cache_policy$force_rebuild) {
  message("Reusing projected regulon AUCell cache: ", cache_path)
  projected <- readRDS(cache_path)
} else {
  message("Loading merged PDO object and canonical state vector.")
  pdo_obj <- readRDS(merged_path)
  final_states <- pdo_normalise_final_state_vector(readRDS(state_path))

  common_cells <- intersect(Seurat::Cells(pdo_obj), names(final_states))
  if (length(common_cells) == 0) {
    stop("No cells overlap between the merged PDO object and state vector.")
  }

  sample_col <- PDO_METADATA_COLUMNS$sample
  if (!sample_col %in% colnames(pdo_obj@meta.data)) {
    stop("Merged PDO metadata lacks required sample column: ", sample_col)
  }

  cell_metadata <- pdo_obj@meta.data[common_cells, , drop = FALSE]
  cell_metadata$cell <- rownames(cell_metadata)
  cell_metadata$sample <- as.character(cell_metadata[[sample_col]])
  cell_metadata$final_state <- as.character(final_states[cell_metadata$cell])
  cell_metadata <- cell_metadata %>%
    filter(
      final_state %in% PDO_STATE_ORDER,
      !is.na(sample),
      nzchar(sample),
      sample != PDO_EXCLUDED_SAMPLE
    ) %>%
    select(cell, sample, final_state, everything())

  if (nrow(cell_metadata) == 0) {
    stop("No canonical-state PDO cells remain after required filtering.")
  }
  if (anyDuplicated(cell_metadata$cell)) {
    stop("Cell identifiers are duplicated in the projected-score metadata.")
  }

  message("Reading RNA counts for ", nrow(cell_metadata), " canonical-state cells.")
  counts_mat <- pdo_get_assay_matrix(pdo_obj, assay = "RNA", layer = "counts")
  missing_count_cells <- setdiff(cell_metadata$cell, colnames(counts_mat))
  if (length(missing_count_cells) > 0) {
    stop(
      length(missing_count_cells),
      " canonical-state cells are absent from the RNA counts matrix."
    )
  }
  counts_mat <- counts_mat[, cell_metadata$cell, drop = FALSE]
  detected_genes <- Matrix::rowSums(counts_mat > 0) > 0
  counts_mat <- counts_mat[detected_genes, , drop = FALSE]

  matched_gene_sets <- lapply(
    regulon_gene_sets,
    function(genes) intersect(unique(genes), rownames(counts_mat))
  )
  if (any(lengths(matched_gene_sets) < 5L)) {
    stop(
      "Fewer than five reference targets are present for: ",
      paste(names(matched_gene_sets)[lengths(matched_gene_sets) < 5L], collapse = ", ")
    )
  }

  message(
    "Building AUCell rankings for ", nrow(counts_mat), " genes by ",
    ncol(counts_mat), " cells."
  )
  rankings <- AUCell::AUCell_buildRankings(
    counts_mat,
    plotStats = FALSE,
    splitByBlocks = TRUE,
    nCores = n_cores,
    verbose = TRUE
  )
  auc_max_rank <- ceiling(0.05 * nrow(rankings))
  projected_auc_obj <- AUCell::AUCell_calcAUC(
    geneSets = matched_gene_sets,
    rankings = rankings,
    nCores = n_cores,
    normAUC = TRUE,
    aucMaxRank = auc_max_rank,
    verbose = TRUE
  )
  projected_auc <- as.matrix(AUCell::getAUC(projected_auc_obj))
  projected_auc <- projected_auc[names(regulon_gene_sets), cell_metadata$cell, drop = FALSE]

  signature_table <- bind_rows(lapply(names(regulon_gene_sets), function(regulon) {
    genes <- regulon_gene_sets[[regulon]]
    data.frame(
      regulon = regulon,
      source_regulon_id = regulon_definitions$source_regulon_id[
        match(regulon, regulon_definitions$regulon)
      ],
      expected_state = regulon_definitions$expected_state[
        match(regulon, regulon_definitions$regulon)
      ],
      target_gene = genes,
      present_in_pdo_counts = genes %in% rownames(counts_mat),
      stringsAsFactors = FALSE
    )
  }))

  projected <- list(
    auc = projected_auc,
    cell_metadata = cell_metadata %>% select(cell, sample, final_state),
    signature_table = signature_table,
    auc_max_rank = auc_max_rank,
    n_ranked_genes = nrow(rankings),
    score_method = "AUCell normalized AUC; top 5% expression rank",
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  saveRDS(projected, cache_path)
  message("Saved projected regulon AUCell cache: ", cache_path)
  rm(pdo_obj, counts_mat, rankings, projected_auc_obj)
  invisible(gc())
}

required_cache_fields <- c(
  "auc", "cell_metadata", "signature_table", "auc_max_rank",
  "n_ranked_genes", "score_method"
)
if (!all(required_cache_fields %in% names(projected))) {
  stop("Projected-score cache is missing required fields.")
}

projected_auc <- as.matrix(projected$auc)
cell_metadata <- projected$cell_metadata
if (!identical(colnames(projected_auc), cell_metadata$cell)) {
  stop("Projected AUCell columns do not exactly match cached cell metadata.")
}
if (!all(regulon_definitions$regulon %in% rownames(projected_auc))) {
  stop("Projected AUCell cache lacks one or more requested regulons.")
}

write.csv(
  projected$signature_table,
  file.path(out_tiers[["tables"]], "Auto_regulon_signature_genes.csv"),
  row.names = FALSE
)

####################
# Per-state activity and per-specimen RSS
####################
sample_ids <- sort(unique(cell_metadata$sample))
sample_state_counts <- tidyr::expand_grid(
  sample = sample_ids,
  final_state = PDO_STATE_ORDER
) %>%
  left_join(
    cell_metadata %>% count(sample, final_state, name = "n_cells"),
    by = c("sample", "final_state")
  ) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L),
    sample_display = make_sample_display(sample),
    final_state = factor(final_state, levels = PDO_STATE_ORDER)
  )

write.csv(
  sample_state_counts,
  file.path(out_tiers[["tables"]], "Auto_sample_state_cell_counts.csv"),
  row.names = FALSE
)

activity_long <- as.data.frame(t(projected_auc), check.names = FALSE) %>%
  mutate(cell = rownames(.), .before = 1) %>%
  left_join(cell_metadata, by = "cell") %>%
  pivot_longer(
    cols = all_of(regulon_definitions$regulon),
    names_to = "regulon",
    values_to = "regulon_auc"
  ) %>%
  left_join(regulon_definitions, by = "regulon")

state_activity_observed <- activity_long %>%
  group_by(sample, regulon, source_regulon_id, expected_state, final_state) %>%
  summarise(
    n_cells = n(),
    mean_auc = mean(regulon_auc, na.rm = TRUE),
    median_auc = median(regulon_auc, na.rm = TRUE),
    sd_auc = sd(regulon_auc, na.rm = TRUE),
    .groups = "drop"
  )

state_activity <- tidyr::expand_grid(
  sample = sample_ids,
  regulon = regulon_definitions$regulon,
  final_state = PDO_STATE_ORDER
) %>%
  left_join(regulon_definitions, by = "regulon") %>%
  left_join(
    state_activity_observed,
    by = c(
      "sample", "regulon", "source_regulon_id", "expected_state",
      "final_state"
    )
  ) %>%
  mutate(
    n_cells = replace_na(n_cells, 0L),
    sample_display = make_sample_display(sample),
    is_expected_state = final_state == expected_state
  )

write.csv(
  state_activity,
  file.path(out_tiers[["tables"]], "Auto_per_sample_regulon_state_activity.csv"),
  row.names = FALSE
)

specificity_rows <- vector(
  "list",
  length(sample_ids) * nrow(regulon_definitions)
)
row_index <- 0L

for (sample_id in sample_ids) {
  sample_cells <- cell_metadata$cell[cell_metadata$sample == sample_id]
  sample_states <- cell_metadata$final_state[match(sample_cells, cell_metadata$cell)]

  for (definition_index in seq_len(nrow(regulon_definitions))) {
    row_index <- row_index + 1L
    regulon <- regulon_definitions$regulon[definition_index]
    source_regulon_id <- regulon_definitions$source_regulon_id[definition_index]
    expected_state <- regulon_definitions$expected_state[definition_index]
    values <- as.numeric(projected_auc[regulon, sample_cells])
    expected_index <- sample_states == expected_state
    rest_index <- !expected_index
    n_expected <- sum(expected_index)
    n_rest <- sum(rest_index)

    binary_labels <- ifelse(expected_index, "Expected state", "Other canonical states")
    binary_rss <- calc_one_rss(
      values,
      binary_labels,
      c("Expected state", "Other canonical states")
    )

    present_states <- PDO_STATE_ORDER[PDO_STATE_ORDER %in% unique(sample_states)]
    multiclass_rss <- calc_one_rss(values, sample_states, present_states)
    expected_multiclass_rss <- unname(multiclass_rss[expected_state])
    other_multiclass <- multiclass_rss[names(multiclass_rss) != expected_state]
    if (length(other_multiclass) > 0 && any(is.finite(other_multiclass))) {
      best_other_index <- which.max(other_multiclass)
      max_other_multiclass_rss <- unname(other_multiclass[best_other_index])
      best_other_state <- names(other_multiclass)[best_other_index]
    } else {
      max_other_multiclass_rss <- NA_real_
      best_other_state <- NA_character_
    }

    state_means <- vapply(PDO_STATE_ORDER, function(state_name) {
      safe_mean(values[sample_states == state_name])
    }, numeric(1))
    other_state_means <- state_means[names(state_means) != expected_state]
    if (any(is.finite(other_state_means))) {
      highest_other_activity_state <- names(which.max(other_state_means))
      highest_other_state_mean_auc <- max(other_state_means, na.rm = TRUE)
    } else {
      highest_other_activity_state <- NA_character_
      highest_other_state_mean_auc <- NA_real_
    }

    support_status <- case_when(
      n_expected == 0L ~ "expected state absent",
      n_rest == 0L ~ "canonical-state rest absent",
      n_expected < min_expected_cells & n_rest < min_rest_cells ~ "low expected and rest cell counts",
      n_expected < min_expected_cells ~ "low expected-state cell count",
      n_rest < min_rest_cells ~ "low rest cell count",
      TRUE ~ "adequate cell counts"
    )

    specificity_rows[[row_index]] <- data.frame(
      sample = sample_id,
      sample_display = make_sample_display(sample_id),
      regulon = regulon,
      source_regulon_id = source_regulon_id,
      expected_state = expected_state,
      n_canonical_cells = length(sample_cells),
      n_expected_cells = n_expected,
      n_rest_cells = n_rest,
      expected_mean_auc = safe_mean(values[expected_index]),
      expected_median_auc = safe_median(values[expected_index]),
      rest_mean_auc = safe_mean(values[rest_index]),
      rest_median_auc = safe_median(values[rest_index]),
      mean_auc_difference = safe_mean(values[expected_index]) - safe_mean(values[rest_index]),
      binary_rss_expected = unname(binary_rss["Expected state"]),
      binary_rss_rest = unname(binary_rss["Other canonical states"]),
      binary_rss_gap = unname(binary_rss["Expected state"] - binary_rss["Other canonical states"]),
      multiclass_rss_expected = expected_multiclass_rss,
      multiclass_max_other_rss = max_other_multiclass_rss,
      multiclass_rss_gap = expected_multiclass_rss - max_other_multiclass_rss,
      multiclass_best_other_state = best_other_state,
      highest_other_activity_state = highest_other_activity_state,
      highest_other_state_mean_auc = highest_other_state_mean_auc,
      adequate_cell_counts = support_status == "adequate cell counts",
      support_status = support_status,
      stringsAsFactors = FALSE
    )
  }
}

specificity <- bind_rows(specificity_rows) %>%
  group_by(regulon) %>%
  mutate(
    expected_activity_percentile = if_else(
      is.na(expected_mean_auc),
      NA_real_,
      rank(expected_mean_auc, ties.method = "average", na.last = "keep") /
        sum(!is.na(expected_mean_auc))
    )
  ) %>%
  ungroup() %>%
  arrange(match(regulon, regulon_definitions$regulon), desc(multiclass_rss_gap), sample)

write.csv(
  specificity,
  file.path(out_tiers[["tables"]], "Auto_per_sample_regulon_specificity.csv"),
  row.names = FALSE
)

personalisation_summary <- specificity %>%
  group_by(regulon, source_regulon_id, expected_state) %>%
  summarise(
    n_specimens = n(),
    n_adequate = sum(adequate_cell_counts),
    n_positive_multiclass_gap_adequate = sum(
      adequate_cell_counts & multiclass_rss_gap > 0,
      na.rm = TRUE
    ),
    n_positive_binary_gap_adequate = sum(
      adequate_cell_counts & binary_rss_gap > 0,
      na.rm = TRUE
    ),
    top_next_best_gap_sample = if (any(adequate_cell_counts & is.finite(multiclass_rss_gap))) {
      sample[which.max(ifelse(adequate_cell_counts, multiclass_rss_gap, -Inf))]
    } else {
      NA_character_
    },
    top_next_best_rss_gap = if (any(adequate_cell_counts & is.finite(multiclass_rss_gap))) {
      max(multiclass_rss_gap[adequate_cell_counts], na.rm = TRUE)
    } else {
      NA_real_
    },
    top_activity_sample = if (any(adequate_cell_counts & is.finite(expected_mean_auc))) {
      sample[which.max(ifelse(adequate_cell_counts, expected_mean_auc, -Inf))]
    } else {
      NA_character_
    },
    top_expected_mean_auc = if (any(adequate_cell_counts & is.finite(expected_mean_auc))) {
      max(expected_mean_auc[adequate_cell_counts], na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

write.csv(
  personalisation_summary,
  file.path(out_tiers[["tables"]], "Auto_regulon_personalisation_summary.csv"),
  row.names = FALSE
)

####################
# Overview evidence matrix
####################
sample_order <- specificity %>%
  group_by(sample, sample_display) %>%
  summarise(
    best_supported_gap = if (any(adequate_cell_counts & is.finite(multiclass_rss_gap))) {
      max(multiclass_rss_gap[adequate_cell_counts], na.rm = TRUE)
    } else {
      -Inf
    },
    .groups = "drop"
  ) %>%
  arrange(desc(best_supported_gap), sample) %>%
  pull(sample_display)

matrix_df <- specificity %>%
  mutate(
    sample_display = factor(sample_display, levels = rev(sample_order)),
    regulon = factor(regulon, levels = regulon_definitions$regulon),
    displayed_rss_gap = ifelse(adequate_cell_counts, multiclass_rss_gap, NA_real_),
    gap_label = ifelse(
      adequate_cell_counts & is.finite(multiclass_rss_gap),
      sprintf("%.2f", multiclass_rss_gap),
      ""
    )
  )

p_matrix <- ggplot(matrix_df, aes(x = regulon, y = sample_display)) +
  geom_tile(aes(fill = displayed_rss_gap), color = "white", linewidth = 0.5) +
  geom_text(aes(label = gap_label), size = 3.0, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    na.value = "white"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Personalised activity of scATLAS regulons in PDO states",
    subtitle = paste0(
      "Gap = RSS(expected state) - maximum RSS(other individual canonical states).\n",
      "Positive values support the expected state.\n",
      "Blank: <", min_expected_cells, " expected-state or <", min_rest_cells,
      " other canonical-state cells."
    ),
    x = NULL,
    y = NULL,
    fill = "Next-best RSS gap"
  ) +
  coord_fixed(ratio = 0.42) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8.5),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

matrix_pdf <- file.path(
  out_tiers[["figures"]],
  "Auto_per_sample_regulon_evidence_matrix.pdf"
)
matrix_png <- file.path(
  out_tiers[["figures"]],
  "Auto_per_sample_regulon_evidence_matrix.png"
)
ggsave(matrix_pdf, p_matrix, width = 10.5, height = 9.5, useDingbats = FALSE)
ggsave(matrix_png, p_matrix, width = 10.5, height = 9.5, dpi = 300)

####################
# Detailed per-regulon profiles
####################
profile_pdf <- file.path(
  out_tiers[["figures"]],
  "Auto_per_sample_regulon_profiles.pdf"
)
grDevices::pdf(profile_pdf, width = 13.333, height = 7.5, useDingbats = FALSE)

for (regulon_name in regulon_definitions$regulon) {
  profile_df <- specificity %>%
    filter(regulon == regulon_name) %>%
    arrange(multiclass_rss_gap, sample) %>%
    mutate(
      sample_display = factor(sample_display, levels = sample_display),
      displayed_expected_rss = ifelse(adequate_cell_counts, multiclass_rss_expected, NA_real_),
      displayed_next_best_rss = ifelse(adequate_cell_counts, multiclass_max_other_rss, NA_real_),
      displayed_rss_gap = ifelse(adequate_cell_counts, multiclass_rss_gap, NA_real_),
      displayed_expected_auc = ifelse(adequate_cell_counts, expected_mean_auc, NA_real_),
      displayed_expected_cells = ifelse(adequate_cell_counts, n_expected_cells, NA_real_)
    )

  expected_state <- unique(profile_df$expected_state)
  state_color <- unname(PDO_STATE_COLORS[expected_state])

  p_rss <- ggplot(profile_df, aes(y = sample_display)) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_segment(
      aes(x = displayed_next_best_rss, xend = displayed_expected_rss, yend = sample_display),
      color = "grey65",
      linewidth = 0.7
    ) +
    geom_point(aes(x = displayed_next_best_rss), color = "grey45", size = 2.5) +
    geom_point(
      aes(x = displayed_expected_rss),
      color = state_color,
      fill = state_color,
      size = 3.1,
      stroke = 1
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = paste0(regulon_name, ": expected-state RSS versus next-best state"),
      subtitle = paste0(
        expected_state, " point in colour; grey point = highest-RSS alternative canonical state. ",
        "Low-cell-support specimen rows are blank."
      ),
      x = "Regulon specificity score (RSS)",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 7.2),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5)
    )

  label_df <- profile_df %>%
    filter(adequate_cell_counts, is.finite(multiclass_rss_gap)) %>%
    slice_max(multiclass_rss_gap, n = 5, with_ties = FALSE)

  p_activity <- ggplot(
    profile_df,
    aes(x = displayed_expected_auc, y = displayed_rss_gap)
  ) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.5) +
    geom_point(
      aes(size = displayed_expected_cells),
      shape = 21,
      fill = state_color,
      color = "black",
      stroke = 0.35,
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = sample_display),
      size = 2.6,
      min.segment.length = 0,
      seed = 19,
      max.overlaps = Inf
    ) +
    scale_size_area(max_size = 8) +
    labs(
      title = "Activity-specificity map",
      subtitle = "Upper-right specimens combine stronger expected-state activity and positive specificity",
      x = paste0(regulon_name, " mean AUCell activity in ", expected_state),
      y = "Expected-state minus next-best-state RSS",
      size = "Expected-state\ncells"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5),
      legend.position = "bottom"
    )

  combined_profile <- p_rss + p_activity +
    patchwork::plot_layout(widths = c(1.35, 1)) +
    patchwork::plot_annotation(
      title = paste0(
        "Personalised PDO regulon evidence: ", regulon_name,
        " (expected state: ", expected_state, ")"
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
    )
  print(combined_profile)
}
grDevices::dev.off()

output_paths <- c(
  cache_path,
  file.path(out_tiers[["tables"]], "Auto_regulon_signature_genes.csv"),
  file.path(out_tiers[["tables"]], "Auto_sample_state_cell_counts.csv"),
  file.path(out_tiers[["tables"]], "Auto_per_sample_regulon_state_activity.csv"),
  file.path(out_tiers[["tables"]], "Auto_per_sample_regulon_specificity.csv"),
  file.path(out_tiers[["tables"]], "Auto_regulon_personalisation_summary.csv"),
  matrix_pdf,
  matrix_png,
  profile_pdf
)

missing_outputs <- output_paths[!file.exists(output_paths)]
if (length(missing_outputs) > 0) {
  stop("Expected output(s) were not created: ", paste(missing_outputs, collapse = ", "))
}

pdo_write_run_summary(
  script = script_path,
  out_dir = out_dir,
  inputs = c(merged_path, state_path),
  outputs = output_paths,
  parameters = list(
    specimen_unit = "orig.ident; treated and untreated kept separate",
    requested_regulons = paste(regulon_definitions$regulon, collapse = ", "),
    included_states = paste(PDO_STATE_ORDER, collapse = ", "),
    excluded_sample = PDO_EXCLUDED_SAMPLE,
    auc_max_rank = projected$auc_max_rank,
    n_ranked_genes = projected$n_ranked_genes,
    min_expected_cells_reliability = min_expected_cells,
    min_rest_cells_reliability = min_rest_cells,
    n_cores = n_cores
  ),
  cache = cache_policy,
  status = "completed",
  include_session = TRUE
)

message("Completed personalised PDO regulon analysis: ", out_dir)
