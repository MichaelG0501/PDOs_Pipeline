####################
# Analysis registry:
#   Status: active terminal; scRef finalised MP expression trends across matched FLOT PDO pairs
#   Script: analysis/metaprograms/Auto_scref_mps_flot_matched_trends.R
#   Methodology: analysis/methodology/metaprograms/Auto_scref_mps_flot_matched_trends_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Scores malignant cells from 4 matched FLOT-treated PDO pairs (8 samples total)
#     using the 17 finalized centred-refined scRef metaprograms (MPs) via UCell.
#     Visualizes expression mean and median change across treatment pairs, boxplots,
#     and heatmaps, computing paired Wilcoxon statistics and direction concordance.
#   Inputs:
#     - scRef MP genes: /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#     - scRef MP grouping: /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#     - PDO post-QC list: PDOs_outs/PDOs_list_PDOs.rds
#   Outputs:
#     - live: PDOs_outs/Auto_scref_mps_flot_matched_trends/
#       UCell scores, sample summaries, paired delta tables, trend statistics,
#       boxplots, mean/median pair trend PDFs, and heatmaps.
#   Downstream use:
#     - Informs treatment response dynamics of reference OAC metaprograms in patient organoids.
#   Cache/replot behavior:
#     - --force or PDO_FORCE_REBUILD=1 rebuilds UCell scores and analytical caches.
#   Run command: qsub analysis/metaprograms/Auto_scref_mps_flot_matched_trends.sh
#   Conda env: dmtcp
####################

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
ncore_arg <- grep("^--ncores=", args, value = TRUE)
score_ncores <- if (length(ncore_arg) > 0) as.integer(sub("^--ncores=", "", ncore_arg[1])) else 1
score_ncores <- max(1, min(score_ncores, 2))

suppressPackageStartupMessages({
  library(UCell)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

####################
# Paths and setup
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

force <- force || pdo_get_env_flag("PDO_FORCE_REBUILD", FALSE)
out_dir <- file.path(PDO_LIVE_OUTS, "Auto_scref_mps_flot_matched_trends")
tier_dirs <- pdo_ensure_output_tiers(out_dir)

# scRef inputs
scref_root <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline"
scref_gene_file <- file.path(
  scref_root, "ref_outs", "Metaprogrammes_Results", "centred", "mp_refinement",
  "intermediate", "merged_refined_mp_genes.rds"
)
scref_group_file <- file.path(
  scref_root, "ref_outs", "Metaprogrammes_Results", "centred", "mp_refinement",
  "tables", "centred_refined_mp_state_grouping.csv"
)
pdo_list_path <- file.path(PDO_LIVE_OUTS, "PDOs_list_PDOs.rds")

pdo_require_files(c(scref_gene_file, scref_group_file, pdo_list_path))

patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
sample_order <- as.vector(rbind(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
))
sample_plot_levels <- c(
  "SUR1070_Untreated_PDO", "SUR1070_Treated_PDO", "gap_1070",
  "SUR1072_Untreated_PDO", "SUR1072_Treated_PDO", "gap_1072",
  "SUR1090_Untreated_PDO", "SUR1090_Treated_PDO", "gap_1090",
  "SUR1181_Untreated_PDO", "SUR1181_Treated_PDO"
)
sample_plot_labels <- c(
  "SUR1070_Untreated_PDO" = "SUR1070_Untreated_PDO",
  "SUR1070_Treated_PDO" = "SUR1070_Treated_PDO",
  "gap_1070" = "",
  "SUR1072_Untreated_PDO" = "SUR1072_Untreated_PDO",
  "SUR1072_Treated_PDO" = "SUR1072_Treated_PDO",
  "gap_1072" = "",
  "SUR1090_Untreated_PDO" = "SUR1090_Untreated_PDO",
  "SUR1090_Treated_PDO" = "SUR1090_Treated_PDO",
  "gap_1090" = "",
  "SUR1181_Untreated_PDO" = "SUR1181_Untreated_PDO",
  "SUR1181_Treated_PDO" = "SUR1181_Treated_PDO"
)

patient_cols <- c(
  SUR1070 = "#4C78A8",
  SUR1072 = "#59A14F",
  SUR1090 = "#B07AA1",
  SUR1181 = "#F28E2B"
)
treatment_cols <- c(Untreated = "#D58B2D", Treated = "#374151")
sample_cols <- unlist(lapply(patient_order, function(patient) {
  base_col <- patient_cols[[patient]]
  c(
    setNames(grDevices::adjustcolor(base_col, alpha.f = 1), paste0(patient, "_Untreated_PDO")),
    setNames(grDevices::adjustcolor(base_col, alpha.f = 0.62), paste0(patient, "_Treated_PDO"))
  )
}), use.names = TRUE)

type_levels <- c("increase", "decrease", "mixed")
type_labels <- c(
  increase = "Increase in treated pairs",
  decrease = "Decrease in treated pairs",
  mixed = "Mixed / no consistent trend"
)

state_order <- c(
  "Cell cycle",
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "Cancer-cell immune mimicry"
)

state_cols <- c(
  "Cell cycle" = "#7F7F7F",
  "Classic proliferation" = "#E41A1C",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)

####################
# Helper functions
####################
get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

chunk_vector <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

format_p_label <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
}

paired_trend_call <- function(mean_vec, median_vec, min_support = 3) {
  if (!all(sample_order %in% names(mean_vec)) ||
      !all(sample_order %in% names(median_vec)) ||
      any(!is.finite(mean_vec[sample_order])) ||
      any(!is.finite(median_vec[sample_order]))) {
    return(list(direction = "mixed", support = 0L))
  }

  pair_calls <- vapply(patient_order, function(patient) {
    untreated <- paste0(patient, "_Untreated_PDO")
    treated <- paste0(patient, "_Treated_PDO")
    mean_delta <- mean_vec[[treated]] - mean_vec[[untreated]]
    median_delta <- median_vec[[treated]] - median_vec[[untreated]]
    if (mean_delta > 0 && median_delta > 0) {
      "increase"
    } else if (mean_delta < 0 && median_delta < 0) {
      "decrease"
    } else {
      "mixed"
    }
  }, character(1))

  increase_n <- sum(pair_calls == "increase")
  decrease_n <- sum(pair_calls == "decrease")
  if (increase_n >= min_support && increase_n > decrease_n) {
    list(direction = "increase", support = increase_n)
  } else if (decrease_n >= min_support && decrease_n > increase_n) {
    list(direction = "decrease", support = decrease_n)
  } else {
    list(direction = "mixed", support = max(increase_n, decrease_n))
  }
}

pairwise_trend_stats <- function(mean_vec, median_vec) {
  mean_deltas <- vapply(patient_order, function(patient) {
    mean_vec[[paste0(patient, "_Treated_PDO")]] - mean_vec[[paste0(patient, "_Untreated_PDO")]]
  }, numeric(1))
  median_deltas <- vapply(patient_order, function(patient) {
    median_vec[[paste0(patient, "_Treated_PDO")]] - median_vec[[paste0(patient, "_Untreated_PDO")]]
  }, numeric(1))
  data.frame(
    mean_wilcox_p = tryCatch(stats::wilcox.test(mean_deltas, mu = 0, paired = FALSE, exact = FALSE)$p.value, error = function(e) NA_real_),
    median_wilcox_p = tryCatch(stats::wilcox.test(median_deltas, mu = 0, paired = FALSE, exact = FALSE)$p.value, error = function(e) NA_real_),
    mean_signed_delta = mean(mean_deltas, na.rm = TRUE),
    median_signed_delta = mean(median_deltas, na.rm = TRUE),
    mean_abs_delta_min = min(abs(mean_deltas), na.rm = TRUE),
    median_abs_delta_min = min(abs(median_deltas), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

add_delta_columns <- function(df, prefix) {
  for (patient in patient_order) {
    untreated_col <- paste0(prefix, "_", patient, "_Untreated_PDO")
    treated_col <- paste0(prefix, "_", patient, "_Treated_PDO")
    delta_col <- paste0(prefix, "_", patient, "_treated_minus_untreated")
    df[[delta_col]] <- df[[treated_col]] - df[[untreated_col]]
  }
  delta_cols <- paste0(prefix, "_", patient_order, "_treated_minus_untreated")
  df[[paste0(prefix, "_delta_mean")]] <- rowMeans(df[, delta_cols, drop = FALSE], na.rm = TRUE)
  df[[paste0(prefix, "_delta_min_abs")]] <- apply(abs(df[, delta_cols, drop = FALSE]), 1, min, na.rm = TRUE)
  df
}

paired_stats <- function(sample_summary, mp_order) {
  stats_list <- lapply(mp_order, function(mp) {
    mp_df <- sample_summary |>
      dplyr::filter(MP == mp) |>
      dplyr::select(sample, mean_score, median_score) |>
      tidyr::pivot_wider(names_from = sample, values_from = c(mean_score, median_score))
    untreated_mean <- as.numeric(mp_df[1, paste0("mean_score_", patient_order, "_Untreated_PDO")])
    treated_mean <- as.numeric(mp_df[1, paste0("mean_score_", patient_order, "_Treated_PDO")])
    untreated_median <- as.numeric(mp_df[1, paste0("median_score_", patient_order, "_Untreated_PDO")])
    treated_median <- as.numeric(mp_df[1, paste0("median_score_", patient_order, "_Treated_PDO")])
    mean_delta <- treated_mean - untreated_mean
    median_delta <- treated_median - untreated_median
    mean_p <- tryCatch(stats::wilcox.test(treated_mean, untreated_mean, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
    median_p <- tryCatch(stats::wilcox.test(treated_median, untreated_median, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
    trend_p <- if (all(is.na(c(mean_p, median_p)))) NA_real_ else max(mean_p, median_p, na.rm = TRUE)
    data.frame(
      MP = mp,
      statistical_unit = "sample_level_pseudobulk",
      test = "paired Wilcoxon signed-rank",
      n_patient_pairs = sum(is.finite(mean_delta) & is.finite(median_delta)),
      mean_wilcox_p = mean_p,
      median_wilcox_p = median_p,
      p_value = trend_p,
      mean_delta = mean(mean_delta, na.rm = TRUE),
      median_delta = mean(median_delta, na.rm = TRUE),
      min_abs_delta = min(abs(c(mean_delta, median_delta)), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(stats_list) |>
    dplyr::mutate(
      mean_wilcox_p_adj = stats::p.adjust(mean_wilcox_p, method = "BH"),
      median_wilcox_p_adj = stats::p.adjust(median_wilcox_p, method = "BH"),
      p_adj = stats::p.adjust(p_value, method = "BH"),
      significance = dplyr::case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      stat_label = paste0(
        "mean p = ", format_p_label(mean_wilcox_p),
        "\nmedian p = ", format_p_label(median_wilcox_p)
      )
    )
}

make_activity_plot <- function(ucell_scores, cell_meta, sample_summary, mp_order, title_text, label_map) {
  activity_long <- as.data.frame(ucell_scores[, mp_order, drop = FALSE]) |>
    tibble::rownames_to_column("cell") |>
    dplyr::left_join(cell_meta, by = "cell") |>
    tidyr::pivot_longer(cols = dplyr::all_of(mp_order), names_to = "MP", values_to = "score") |>
    dplyr::mutate(
      sample = factor(sample, levels = sample_order),
      sample_plot = factor(as.character(sample), levels = sample_plot_levels),
      MP = factor(MP, levels = mp_order),
      display_label = factor(label_map[as.character(MP)], levels = label_map[mp_order])
    )

  activity_stats <- paired_stats(sample_summary, mp_order)
  annot_df <- activity_long |>
    dplyr::group_by(MP, display_label) |>
    dplyr::summarise(y_pos = max(score, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(activity_stats, by = "MP") |>
    dplyr::mutate(y_pos = y_pos + 0.012, label = stat_label)

  p <- ggplot2::ggplot(activity_long, ggplot2::aes(x = sample_plot, y = score, fill = sample)) +
    ggplot2::geom_boxplot(
      width = 0.78,
      outlier.shape = NA,
      alpha = 0.78,
      linewidth = 0.28,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = annot_df,
      ggplot2::aes(x = "SUR1090_Untreated_PDO", y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 2.65,
      lineheight = 0.92,
      fontface = "bold"
    ) +
    ggplot2::facet_wrap(~display_label, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = sample_cols, drop = FALSE, name = "Sample") +
    ggplot2::scale_x_discrete(labels = sample_plot_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.10))) +
    ggplot2::labs(title = title_text, x = NULL, y = "UCell score") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 22),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      axis.line.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 18, 12, 12)
    )

  list(plot = p, stats = activity_stats)
}

####################
# 1. Load scRef MP definitions
####################
message("Loading finalized scRef MP genes: ", scref_gene_file)
scref_genes <- readRDS(scref_gene_file)
scref_grouping <- read.csv(scref_group_file, stringsAsFactors = FALSE, check.names = FALSE)

# Ensure canonical order based on grouping
scref_grouping$state <- factor(scref_grouping$state, levels = state_order)
scref_grouping <- scref_grouping[order(scref_grouping$state), ]
canonical_mp_order <- scref_grouping$mp[scref_grouping$mp %in% names(scref_genes)]

# Create clean two-line display labels
display_label_map <- setNames(
  paste0(scref_grouping$mp, "\n", scref_grouping$description),
  scref_grouping$mp
)

# Export scRef MP gene table to local output
scref_gene_table <- do.call(rbind, lapply(canonical_mp_order, function(mp) {
  data.frame(
    MP = mp,
    State = scref_grouping$state[scref_grouping$mp == mp][1],
    Description = scref_grouping$description[scref_grouping$mp == mp][1],
    Rank = seq_along(scref_genes[[mp]]),
    Gene = scref_genes[[mp]],
    stringsAsFactors = FALSE
  )
}))
write.csv(scref_gene_table, file.path(tier_dirs[["tables"]], "Auto_scref_finalised_mp_genes.csv"), row.names = FALSE)

####################
# 2. Score matched PDO cells with UCell
####################
ucell_path <- file.path(tier_dirs[["intermediate"]], "Auto_scref_mps_UCell_scores_matched_pdos.rds")
cell_meta_path <- file.path(tier_dirs[["intermediate"]], "Auto_scref_mps_cell_metadata_matched_pdos.rds")

if (file.exists(ucell_path) && file.exists(cell_meta_path) && !force) {
  message("Loading existing UCell scores: ", ucell_path)
  ucell_scores <- readRDS(ucell_path)
  cell_meta <- readRDS(cell_meta_path)
} else {
  message("Loading persistent post-QC PDO list for matched-sample UCell scoring.")
  pdos_list <- readRDS(pdo_list_path)
  pdos_list[[PDO_EXCLUDED_SAMPLE]] <- NULL
  missing_samples <- setdiff(sample_order, names(pdos_list))
  if (length(missing_samples) > 0L) {
    stop("Matched PDO sample(s) missing: ", paste(missing_samples, collapse = ", "))
  }
  pdos_list <- pdos_list[sample_order]
  sample_genes <- lapply(pdos_list, rownames)
  common_genes <- Reduce(intersect, sample_genes)
  if (length(common_genes) == 0) {
    stop("No common genes found across matched PDO samples.")
  }

  counts_list <- list()
  meta_list <- list()
  for (sample in sample_order) {
    message("Loading counts for ", sample)
    obj <- pdos_list[[sample]]
    old_cells <- colnames(obj)
    new_cells <- paste(sample, old_cells, sep = "_")
    counts <- get_counts(obj)[common_genes, , drop = FALSE]
    colnames(counts) <- new_cells
    counts_list[[sample]] <- counts
    meta_list[[sample]] <- data.frame(
      cell = new_cells,
      original_cell = old_cells,
      sample = sample,
      patient = sub("_(Untreated|Treated)_PDO$", "", sample),
      treatment = ifelse(grepl("_Treated_", sample), "Treated", "Untreated"),
      stringsAsFactors = FALSE
    )
    rm(obj, counts)
    gc()
  }
  cell_meta <- dplyr::bind_rows(meta_list)
  cell_meta$sample <- factor(cell_meta$sample, levels = sample_order)
  cell_meta$patient <- factor(cell_meta$patient, levels = patient_order)
  cell_meta$treatment <- factor(cell_meta$treatment, levels = c("Untreated", "Treated"))
  rownames(cell_meta) <- cell_meta$cell
  counts_all <- do.call(cbind, counts_list)

  score_features <- lapply(scref_genes[canonical_mp_order], intersect, rownames(counts_all))
  feature_counts <- lengths(score_features)
  write.csv(
    data.frame(MP = names(feature_counts), scored_gene_count = as.integer(feature_counts), stringsAsFactors = FALSE),
    file.path(tier_dirs[["tables"]], "Auto_scref_mps_scored_gene_counts.csv"),
    row.names = FALSE
  )
  score_features <- score_features[feature_counts > 0]

  message("Scoring ", length(score_features), " scRef MP signatures on ", ncol(counts_all), " cells with UCell.")
  ucell_scores <- UCell::ScoreSignatures_UCell(
    matrix = counts_all,
    features = score_features,
    maxRank = 1500,
    chunk.size = 1000,
    ncores = score_ncores,
    force.gc = TRUE
  )
  ucell_scores <- as.data.frame(ucell_scores)
  colnames(ucell_scores) <- sub("_UCell$", "", colnames(ucell_scores))
  ucell_scores <- as.matrix(ucell_scores)

  saveRDS(ucell_scores, ucell_path, compress = FALSE)
  saveRDS(cell_meta, cell_meta_path, compress = FALSE)
  rm(counts_all, counts_list, pdos_list)
  gc()
}

####################
# 3. Calculate sample-level summaries & paired trends
####################
mp_order <- intersect(canonical_mp_order, colnames(ucell_scores))
if (length(mp_order) == 0) {
  stop("UCell score matrix does not contain any expected scRef MP columns.")
}

score_long <- as.data.frame(ucell_scores[, mp_order, drop = FALSE]) |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(cell_meta[, c("cell", "sample", "patient", "treatment")], by = "cell") |>
  tidyr::pivot_longer(cols = dplyr::all_of(mp_order), names_to = "MP", values_to = "score")

sample_summary <- score_long |>
  dplyr::group_by(MP, sample, patient, treatment) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = stats::median(score, na.rm = TRUE),
    q1 = stats::quantile(score, 0.25, na.rm = TRUE),
    q3 = stats::quantile(score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    sample = factor(as.character(sample), levels = sample_order),
    patient = factor(as.character(patient), levels = patient_order),
    treatment = factor(as.character(treatment), levels = c("Untreated", "Treated"))
  )
write.csv(sample_summary, file.path(tier_dirs[["tables"]], "Auto_scref_mps_sample_ucell_summary.csv"), row.names = FALSE)

pair_delta_summary <- sample_summary |>
  dplyr::select(MP, patient, treatment, mean_score, median_score) |>
  tidyr::pivot_wider(names_from = treatment, values_from = c(mean_score, median_score)) |>
  dplyr::mutate(
    mean_delta = mean_score_Treated - mean_score_Untreated,
    median_delta = median_score_Treated - median_score_Untreated
  )
write.csv(pair_delta_summary, file.path(tier_dirs[["tables"]], "Auto_scref_mps_paired_delta_summary.csv"), row.names = FALSE)

mean_wide <- sample_summary |>
  dplyr::select(MP, sample, mean_score) |>
  tidyr::pivot_wider(names_from = sample, values_from = mean_score, names_prefix = "mean_")
median_wide <- sample_summary |>
  dplyr::select(MP, sample, median_score) |>
  tidyr::pivot_wider(names_from = sample, values_from = median_score, names_prefix = "median_")
trend_summary <- dplyr::left_join(mean_wide, median_wide, by = "MP")

trend_calls <- lapply(seq_len(nrow(trend_summary)), function(i) {
  mean_vec <- setNames(as.numeric(trend_summary[i, paste0("mean_", sample_order)]), sample_order)
  median_vec <- setNames(as.numeric(trend_summary[i, paste0("median_", sample_order)]), sample_order)
  call <- paired_trend_call(mean_vec, median_vec, min_support = 3)
  stats <- pairwise_trend_stats(mean_vec, median_vec)
  dplyr::bind_cols(
    data.frame(
      treatment_direction = call$direction,
      pair_support_n = as.integer(call$support),
      stringsAsFactors = FALSE
    ),
    stats
  )
})
trend_calls <- dplyr::bind_rows(trend_calls)
trend_summary <- dplyr::bind_cols(trend_summary, trend_calls)
trend_summary$concordant_trend <- trend_summary$treatment_direction %in% c("increase", "decrease")
trend_summary$trend_p_value <- pmax(trend_summary$mean_wilcox_p, trend_summary$median_wilcox_p, na.rm = TRUE)
trend_summary$trend_p_value[!is.finite(trend_summary$trend_p_value)] <- NA_real_
trend_summary$mean_wilcox_p_adj <- stats::p.adjust(trend_summary$mean_wilcox_p, method = "BH")
trend_summary$median_wilcox_p_adj <- stats::p.adjust(trend_summary$median_wilcox_p, method = "BH")
trend_summary$trend_p_adj <- pmax(trend_summary$mean_wilcox_p_adj, trend_summary$median_wilcox_p_adj, na.rm = TRUE)
trend_summary$trend_p_adj[!is.finite(trend_summary$trend_p_adj)] <- NA_real_
trend_summary <- add_delta_columns(trend_summary, "mean")
trend_summary <- add_delta_columns(trend_summary, "median")

# Annotate with State and scRef MP metadata
trend_summary <- trend_summary |>
  dplyr::left_join(scref_grouping[, c("mp", "state", "description", "plot_label")], by = c("MP" = "mp")) |>
  dplyr::arrange(factor(state, levels = state_order), MP)
write.csv(trend_summary, file.path(tier_dirs[["tables"]], "Auto_scref_mps_trend_summary.csv"), row.names = FALSE)

####################
# 4. Generate Visualizations
####################

# Compute activity stats for annotation
paired_stat_df <- paired_stats(sample_summary, mp_order)

# Prepare Long summary for trend plotting
trend_long <- sample_summary |>
  dplyr::select(MP, sample, patient, treatment, mean_score, median_score) |>
  dplyr::mutate(
    sample_plot = factor(as.character(sample), levels = sample_plot_levels),
    MP = factor(MP, levels = mp_order)
  ) |>
  tidyr::pivot_longer(cols = c(mean_score, median_score), names_to = "summary_stat", values_to = "score") |>
  dplyr::mutate(
    summary_stat = dplyr::recode(summary_stat, mean_score = "Mean", median_score = "Median"),
    display_label = factor(display_label_map[as.character(MP)], levels = display_label_map[mp_order])
  )

trend_y_limits <- range(trend_long$score, na.rm = TRUE)
trend_y_pad <- diff(trend_y_limits) * 0.03
if (!is.finite(trend_y_pad) || trend_y_pad == 0) trend_y_pad <- 0.01
trend_y_limits <- trend_y_limits + c(-trend_y_pad, trend_y_pad)

trend_annot_df <- paired_stat_df |>
  dplyr::distinct(MP, stat_label) |>
  dplyr::mutate(
    sample_plot = factor("SUR1090_Untreated_PDO", levels = sample_plot_levels),
    score = trend_y_limits[2] - diff(trend_y_limits) * 0.015,
    display_label = factor(display_label_map[MP], levels = display_label_map[mp_order])
  )

# --- Plot 1: Mean & Median Pair Trends PDF (All 17 scRef MPs grouped by state & canonical order) ---
trend_pdf_path <- file.path(tier_dirs[["figures"]], "Auto_scref_mps_mean_median_pair_trends.pdf")
pdf(trend_pdf_path, width = 15, height = 11, useDingbats = FALSE)

# Page 1: All 17 scRef MPs in canonical order (3 pages max, 6 per page or 12 per page)
chunks_all <- chunk_vector(mp_order, 12)
for (i in seq_along(chunks_all)) {
  chunk_mps <- chunks_all[[i]]
  chunk_labels <- display_label_map[chunk_mps]
  p_chunk <- trend_long |>
    dplyr::filter(MP %in% chunk_mps) |>
    dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)) |>
    ggplot2::ggplot(ggplot2::aes(x = sample_plot, y = score, group = interaction(patient, summary_stat), linetype = summary_stat)) +
    ggplot2::geom_line(color = "grey25", linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(fill = sample), shape = 21, size = 2.9, color = "black") +
    ggplot2::geom_text(
      data = dplyr::filter(trend_annot_df, MP %in% chunk_mps) |> dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)),
      ggplot2::aes(x = sample_plot, y = score, label = stat_label),
      inherit.aes = FALSE, size = 2.45, lineheight = 0.92, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = sample_cols, guide = "none") +
    ggplot2::scale_linetype_manual(values = c("Mean" = "solid", "Median" = "dashed"), name = NULL) +
    ggplot2::scale_x_discrete(labels = sample_plot_labels, drop = FALSE) +
    ggplot2::facet_wrap(~display_label, scales = "fixed", ncol = 4) +
    ggplot2::coord_cartesian(ylim = trend_y_limits) +
    ggplot2::labs(
      title = paste0("scRef Finalised Metaprograms — Mean & Median UCell Pair Trends (Part ", i, " of ", length(chunks_all), ")"),
      subtitle = "4 Matched Treated vs Untreated FLOT Organoid Pairs",
      x = NULL, y = "UCell score"
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 13, color = "grey30"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 9),
      axis.text.y = ggplot2::element_text(colour = "black", size = 10),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 16)
    )
  print(p_chunk)
}

# Subsequent Pages: Stratified by Trend Direction (Increase, Decrease, Mixed)
for (direction_name in type_levels) {
  dir_mps <- trend_summary$MP[trend_summary$treatment_direction == direction_name]
  dir_mps <- intersect(mp_order, dir_mps)
  if (length(dir_mps) == 0) next
  dir_chunks <- chunk_vector(dir_mps, 12)
  for (j in seq_along(dir_chunks)) {
    chunk_mps <- dir_chunks[[j]]
    chunk_labels <- display_label_map[chunk_mps]
    p_dir <- trend_long |>
      dplyr::filter(MP %in% chunk_mps) |>
      dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)) |>
      ggplot2::ggplot(ggplot2::aes(x = sample_plot, y = score, group = interaction(patient, summary_stat), linetype = summary_stat)) +
      ggplot2::geom_line(color = "grey25", linewidth = 0.55, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(fill = sample), shape = 21, size = 2.9, color = "black") +
      ggplot2::geom_text(
        data = dplyr::filter(trend_annot_df, MP %in% chunk_mps) |> dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)),
        ggplot2::aes(x = sample_plot, y = score, label = stat_label),
        inherit.aes = FALSE, size = 2.45, lineheight = 0.92, fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = sample_cols, guide = "none") +
      ggplot2::scale_linetype_manual(values = c("Mean" = "solid", "Median" = "dashed"), name = NULL) +
      ggplot2::scale_x_discrete(labels = sample_plot_labels, drop = FALSE) +
      ggplot2::facet_wrap(~display_label, scales = "fixed", ncol = 4) +
      ggplot2::coord_cartesian(ylim = trend_y_limits) +
      ggplot2::labs(
        title = paste0("scRef MP Pair Trends: ", type_labels[[direction_name]]),
        x = NULL, y = "UCell score"
      ) +
      ggplot2::theme_classic(base_size = 18) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 20),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 9),
        axis.text.y = ggplot2::element_text(colour = "black", size = 10),
        strip.text = ggplot2::element_text(face = "bold", size = 10),
        legend.position = "top",
        legend.text = ggplot2::element_text(size = 16)
      )
    print(p_dir)
  }
}
dev.off()
message("Saved mean/median trend PDF: ", trend_pdf_path)

# --- Plot 2: Boxplots PDF (All 17 scRef MPs per sample) ---
box_pdf_path <- file.path(tier_dirs[["figures"]], "Auto_scref_mps_activity_boxplots.pdf")
pdf(box_pdf_path, width = 20, height = 9.5, useDingbats = FALSE)
box_chunks <- chunk_vector(mp_order, 8)
for (k in seq_along(box_chunks)) {
  act_plot <- make_activity_plot(
    ucell_scores,
    cell_meta,
    sample_summary,
    box_chunks[[k]],
    paste0("scRef Finalised MP UCell Activity across Matched PDO Pairs (Part ", k, " of ", length(box_chunks), ")"),
    display_label_map
  )
  print(act_plot$plot)
}
dev.off()
message("Saved boxplot PDF: ", box_pdf_path)

# --- Plot 3: Heatmap of Mean UCell Activity across Matched PDO Samples ---
mat_df <- sample_summary |>
  dplyr::select(MP, sample, mean_score) |>
  tidyr::pivot_wider(names_from = sample, values_from = mean_score) |>
  as.data.frame()
row.names(mat_df) <- mat_df$MP
mat <- as.matrix(mat_df[mp_order, sample_order, drop = FALSE])

annotation_row <- trend_summary |>
  dplyr::filter(MP %in% rownames(mat)) |>
  dplyr::mutate(
    Pair_support = paste0(pair_support_n, "/4"),
    Direction = treatment_direction,
    State = state
  ) |>
  dplyr::select(MP, State, Direction, Pair_support) |>
  as.data.frame()
row.names(annotation_row) <- annotation_row$MP
annotation_row <- annotation_row[rownames(mat), c("State", "Direction", "Pair_support"), drop = FALSE]

anno_colors <- list(
  State = state_cols,
  Direction = c(increase = "#de2d26", decrease = "#3182bd", mixed = "#969696"),
  Pair_support = c("0/4" = "#f7f7f7", "1/4" = "#d9d9d9", "2/4" = "#bdbdbd", "3/4" = "#c2855a", "4/4" = "#2b6a8e")
)

heat_cols <- grDevices::colorRampPalette(c("white", "#fee0d2", "#fc9272", "#de2d26", "#67000d"))(100)
heat_breaks <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = length(heat_cols) + 1)
labels_row <- display_label_map[rownames(mat)]
labels_col <- sample_plot_labels[sample_order]

heatmap_pdf_path <- file.path(tier_dirs[["figures"]], "Auto_scref_mps_mean_activity_heatmap.pdf")
pdf(heatmap_pdf_path, width = 11, height = max(7, 0.38 * nrow(mat) + 2.5), useDingbats = FALSE)
pheatmap::pheatmap(
  mat,
  color = heat_cols,
  breaks = heat_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(2, 4, 6),
  annotation_row = annotation_row,
  annotation_colors = anno_colors,
  labels_row = labels_row,
  labels_col = labels_col,
  border_color = "white",
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 11,
  main = "scRef Finalised Metaprograms — Mean UCell Activity in Matched PDO Pairs"
)
dev.off()

heatmap_png_path <- file.path(tier_dirs[["figures"]], "Auto_scref_mps_mean_activity_heatmap.png")
png(heatmap_png_path, width = 3200, height = max(2000, 100 * nrow(mat) + 800), res = 300)
pheatmap::pheatmap(
  mat,
  color = heat_cols,
  breaks = heat_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(2, 4, 6),
  annotation_row = annotation_row,
  annotation_colors = anno_colors,
  labels_row = labels_row,
  labels_col = labels_col,
  border_color = "white",
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 11,
  main = "scRef Finalised Metaprograms — Mean UCell Activity in Matched PDO Pairs"
)
dev.off()

####################
# 5. Export Summary Excel Workbook
####################
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  
  # Sheet 1: Trend Summary
  openxlsx::addWorksheet(wb, "Trend_Summary")
  openxlsx::writeData(wb, sheet = "Trend_Summary", x = trend_summary)
  
  # Sheet 2: Sample Summary
  openxlsx::addWorksheet(wb, "Sample_Summary")
  openxlsx::writeData(wb, sheet = "Sample_Summary", x = sample_summary)
  
  # Sheet 3: Paired Deltas
  openxlsx::addWorksheet(wb, "Paired_Deltas")
  openxlsx::writeData(wb, sheet = "Paired_Deltas", x = pair_delta_summary)
  
  # Sheet 4: MP Genes
  openxlsx::addWorksheet(wb, "MP_Genes")
  openxlsx::writeData(wb, sheet = "MP_Genes", x = scref_gene_table)
  
  wb_path <- file.path(tier_dirs[["tables"]], "Auto_scref_mps_matched_flot_summary.xlsx")
  openxlsx::saveWorkbook(wb, wb_path, overwrite = TRUE)
  message("Saved Excel summary workbook: ", wb_path)
}

# Write summary run log
run_log_path <- file.path(tier_dirs[["logs"]], "Auto_scref_mps_flot_matched_trends_summary.txt")
writeLines(
  c(
    paste0("scRef Finalised MP Matched FLOT PDO Analysis Summary"),
    paste0("Timestamp: ", Sys.time()),
    paste0("Total scRef MPs analyzed: ", length(mp_order)),
    paste0("Matched PDO Samples: ", paste(sample_order, collapse = ", ")),
    paste0("Increased in FLOT (>=3/4 pairs): ", sum(trend_summary$treatment_direction == "increase")),
    paste0("Decreased in FLOT (>=3/4 pairs): ", sum(trend_summary$treatment_direction == "decrease")),
    paste0("Mixed/Inconclusive: ", sum(trend_summary$treatment_direction == "mixed")),
    "",
    "Outputs generated:",
    paste0(" - Figures: ", tier_dirs[["figures"]]),
    paste0(" - Tables: ", tier_dirs[["tables"]]),
    paste0(" - Intermediates: ", tier_dirs[["intermediate"]])
  ),
  run_log_path
)

message("Auto_scref_mps_flot_matched_trends.R finished successfully.")
