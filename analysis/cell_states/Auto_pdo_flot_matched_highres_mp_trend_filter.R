####################
# Auto_pdo_flot_matched_highres_mp_trend_filter.R
#
# High-resolution matched-FLOT PDO metaprogram trend filtering.
# Starts from Auto_pdo_flot_matched_geneNMF.R outputs, sets nMP to half of all
# NMF programmes, scores all eight matched PDO samples with UCell, and retains
# MPs whose mean and median UCell activity are both higher in treated samples
# than untreated samples, or both lower, in at least three of four matched pairs.
#
# Env: gnmf
####################

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
ncore_arg <- grep("^--ncores=", args, value = TRUE)
score_ncores <- if (length(ncore_arg) > 0) as.integer(sub("^--ncores=", "", ncore_arg[1])) else 2
score_ncores <- min(score_ncores, 2)

suppressPackageStartupMessages({
  library(GeneNMF)
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
# setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_pdo_flot_highres_metaprogram_trends"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

max_mps_per_boxplot_page <- 8
max_mps_per_trend_page <- 12
type_levels <- c("increase", "decrease")
type_labels <- c(
  increase = "Increase in treated pairs",
  decrease = "Decrease in treated pairs"
)

geneNMF_program_path <- file.path(out_dir, "Auto_pdo_flot_matched_geneNMF_outs.rds")

####################
# helpers
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

program_count <- function(x) {
  if (is.list(x) && !is.null(x$w)) {
    return(ncol(x$w))
  }
  if (is.matrix(x)) {
    return(ncol(x))
  }
  NA_integer_
}

write_mp_gene_table <- function(mp_genes, path) {
  if (length(mp_genes) == 0) {
    gene_table <- data.frame(MP = character(), rank = integer(), gene = character(), stringsAsFactors = FALSE)
  } else {
    gene_table <- do.call(rbind, lapply(names(mp_genes), function(mp) {
      data.frame(
        MP = mp,
        rank = seq_along(mp_genes[[mp]]),
        gene = mp_genes[[mp]],
        stringsAsFactors = FALSE
      )
    }))
  }
  write.csv(gene_table, path, row.names = FALSE)
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

paired_trend_call <- function(mean_vec, median_vec, min_support = 3) {
  if (!all(sample_order %in% names(mean_vec)) ||
      !all(sample_order %in% names(median_vec)) ||
      any(!is.finite(mean_vec[sample_order])) ||
      any(!is.finite(median_vec[sample_order]))) {
    return(list(direction = "not_retained", support = 0L))
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
    list(direction = "not_retained", support = max(increase_n, decrease_n))
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

order_mps_by_metrics <- function(mp_vec, trend_summary) {
  if (length(mp_vec) == 0) return(character())
  trend_summary |>
    dplyr::filter(MP %in% mp_vec) |>
    dplyr::mutate(
      numberPrograms = suppressWarnings(as.integer(numberPrograms)),
      numberPrograms = tidyr::replace_na(numberPrograms, -Inf),
      silhouette = tidyr::replace_na(silhouette, -Inf),
      pair_support_n = tidyr::replace_na(pair_support_n, -Inf),
      trend_p_value = tidyr::replace_na(trend_p_value, Inf),
      mean_delta_min_abs = tidyr::replace_na(mean_delta_min_abs, -Inf)
    ) |>
    dplyr::arrange(trend_p_value, dplyr::desc(pair_support_n), dplyr::desc(mean_delta_min_abs), dplyr::desc(numberPrograms), dplyr::desc(silhouette), MP) |>
    dplyr::pull(MP)
}

order_retained_by_type <- function(trend_summary) {
  unlist(lapply(type_levels, function(type_name) {
    order_mps_by_metrics(trend_summary$MP[trend_summary$retained & trend_summary$treatment_direction == type_name], trend_summary)
  }), use.names = FALSE)
}

read_cell_cycle_genes <- function() {
  cc_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"
  if (!file.exists(cc_path)) return(character())
  cc <- read.csv(cc_path, check.names = FALSE, stringsAsFactors = FALSE)
  unique(stats::na.omit(as.character(unlist(cc, use.names = FALSE))))
}

load_3ca_gene_sets <- function() {
  mp_csv <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
  if (!file.exists(mp_csv)) return(list())
  mp_list <- read.csv(mp_csv, check.names = FALSE, stringsAsFactors = FALSE)
  mp_list <- as.list(mp_list)
  mp_list <- lapply(mp_list, function(x) unique(x[x != "" & !is.na(x)]))
  names(mp_list) <- sub("^MP", "3CA_mp", names(mp_list))
  mp_list[lengths(mp_list) > 0]
}

make_3ca_label_table <- function(mp_genes, nMP) {
  sets_3ca <- load_3ca_gene_sets()
  out_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_top_3CA_noncellcycle_nMP", nMP, ".csv"))
  if (length(sets_3ca) == 0 || length(mp_genes) == 0) {
    empty <- data.frame(MP = names(mp_genes), top_3ca_noncc = NA_character_, top_3ca_noncc_p_adj = NA_real_)
    write.csv(empty, out_path, row.names = FALSE)
    return(empty)
  }

  cc_genes <- read_cell_cycle_genes()
  universe <- unique(c(unlist(mp_genes, use.names = FALSE), unlist(sets_3ca, use.names = FALSE)))
  cc_terms <- vapply(sets_3ca, function(genes) {
    cc_overlap <- length(intersect(genes, cc_genes))
    grepl("cell.?cycle|g1s|g2m|mitotic|prolifer|replication", paste(genes, collapse = " "), ignore.case = TRUE) ||
      (cc_overlap >= 5 && cc_overlap / length(genes) >= 0.05)
  }, logical(1))

  enrich_rows <- do.call(rbind, lapply(names(mp_genes), function(mp) {
    genes <- unique(mp_genes[[mp]])
    do.call(rbind, lapply(names(sets_3ca), function(term) {
      term_genes <- unique(sets_3ca[[term]])
      overlap <- length(intersect(genes, term_genes))
      mat <- matrix(
        c(
          overlap,
          length(genes) - overlap,
          length(term_genes) - overlap,
          length(universe) - length(genes) - length(term_genes) + overlap
        ),
        nrow = 2
      )
      data.frame(
        MP = mp,
        term = term,
        overlap = overlap,
        p_value = tryCatch(stats::fisher.test(mat, alternative = "greater")$p.value, error = function(e) NA_real_),
        is_cell_cycle_term = isTRUE(cc_terms[[term]]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  enrich_rows$p_adj <- stats::p.adjust(enrich_rows$p_value, method = "BH")
  write.csv(enrich_rows, file.path(out_dir, paste0("Auto_pdo_flot_highres_3CA_base_enrichment_nMP", nMP, ".csv")), row.names = FALSE)

  label_table <- enrich_rows |>
    dplyr::filter(!is_cell_cycle_term, overlap > 0, is.finite(p_adj)) |>
    dplyr::arrange(MP, p_adj, dplyr::desc(overlap)) |>
    dplyr::group_by(MP) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(MP, top_3ca_noncc = term, top_3ca_noncc_p_adj = p_adj)

  label_table <- data.frame(MP = names(mp_genes), stringsAsFactors = FALSE) |>
    dplyr::left_join(label_table, by = "MP")
  write.csv(label_table, out_path, row.names = FALSE)
  label_table
}

make_display_labels <- function(trend_summary, nMP) {
  label_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_top_3CA_noncellcycle_nMP", nMP, ".csv"))
  label_table <- if (file.exists(label_path)) {
    read.csv(label_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame(MP = trend_summary$MP, top_3ca_noncc = NA_character_, stringsAsFactors = FALSE)
  }
  label_df <- trend_summary |>
    dplyr::left_join(label_table, by = "MP") |>
    dplyr::mutate(
      top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
      single_programme_label = ifelse(numberPrograms == 1, "\n[1 programme]", ""),
      display_label = paste0(MP, "\n", top_3ca_noncc, single_programme_label)
    )
  setNames(label_df$display_label, label_df$MP)
}

####################
format_p_label <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
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
####################

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
      plot.title = ggplot2::element_text(face = "bold", size = 24),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      axis.line.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "none",
      plot.margin = ggplot2::margin(12, 18, 12, 12)
    )

  list(plot = p, stats = activity_stats)
}

plot_retained_outputs <- function(ucell_scores, cell_meta, sample_summary, trend_summary, nMP) {
  retained_mps <- order_retained_by_type(trend_summary)
  if (length(retained_mps) == 0) {
    writeLines(
      "No metaprograms passed the strict all-four-pairs mean UCell trend filter.",
      file.path(out_dir, paste0("Auto_pdo_flot_highres_no_retained_nMP", nMP, ".txt"))
    )
    return(invisible(NULL))
  }

  label_map <- make_display_labels(trend_summary, nMP)
  boxplot_stats <- list()
  pdf(
    file.path(out_dir, paste0("Auto_pdo_flot_highres_activity_boxplots_nMP", nMP, "_selected.pdf")),
    width = 20,
    height = 9,
    useDingbats = FALSE
  )
  for (direction_name in type_levels) {
    direction_mps <- order_mps_by_metrics(
      trend_summary$MP[trend_summary$retained & trend_summary$treatment_direction == direction_name],
      trend_summary
    )
    if (length(direction_mps) == 0) next
    chunks <- chunk_vector(direction_mps, max_mps_per_boxplot_page)
    for (i in seq_along(chunks)) {
      activity <- make_activity_plot(
        ucell_scores,
        cell_meta,
        sample_summary,
        chunks[[i]],
        paste0("MP Activity, ", type_labels[[direction_name]]),
        label_map
      )
      print(activity$plot)
      boxplot_stats[[paste(direction_name, i, sep = "_")]] <- activity$stats |>
        dplyr::mutate(treatment_direction = direction_name, treatment_direction_label = type_labels[[direction_name]])
    }
  }
  dev.off()
  write.csv(
    dplyr::bind_rows(boxplot_stats),
    file.path(out_dir, paste0("Auto_pdo_flot_highres_activity_boxplots_nMP", nMP, "_selected_stats.csv")),
    row.names = FALSE
  )

  retained_summary <- sample_summary |>
    dplyr::filter(MP %in% retained_mps) |>
    dplyr::mutate(
      sample = factor(sample, levels = sample_order),
      sample_plot = factor(as.character(sample), levels = sample_plot_levels),
      MP = factor(MP, levels = retained_mps)
    )

  trend_long <- retained_summary |>
    dplyr::select(MP, sample, sample_plot, patient, treatment, mean_score, median_score) |>
    tidyr::pivot_longer(cols = c(mean_score, median_score), names_to = "summary_stat", values_to = "score") |>
    dplyr::mutate(
      summary_stat = dplyr::recode(summary_stat, mean_score = "Mean", median_score = "Median"),
      display_label = label_map[as.character(MP)]
    )
  trend_y_limits <- range(trend_long$score, na.rm = TRUE)
  trend_y_pad <- diff(trend_y_limits) * 0.03
  if (!is.finite(trend_y_pad) || trend_y_pad == 0) trend_y_pad <- 0.01
  trend_y_limits <- trend_y_limits + c(-trend_y_pad, trend_y_pad)
  trend_annot_df <- if (length(boxplot_stats) > 0) {
      dplyr::bind_rows(boxplot_stats) |>
      dplyr::distinct(MP, stat_label) |>
      dplyr::mutate(sample_plot = factor("SUR1090_Untreated_PDO", levels = sample_plot_levels), score = trend_y_limits[2] - diff(trend_y_limits) * 0.015, display_label = label_map[MP])
  } else {
    data.frame(MP = character(), stat_label = character(), sample_plot = factor(character(), levels = sample_plot_levels), score = numeric(), display_label = character())
  }

  pdf(
    file.path(out_dir, paste0("Auto_pdo_flot_highres_mean_median_pair_trends_nMP", nMP, "_selected.pdf")),
    width = 14,
    height = 10,
    useDingbats = FALSE
  )
  for (direction_name in type_levels) {
    direction_mps <- order_mps_by_metrics(
      trend_summary$MP[trend_summary$retained & trend_summary$treatment_direction == direction_name],
      trend_summary
    )
    if (length(direction_mps) == 0) next
    chunks <- chunk_vector(direction_mps, max_mps_per_trend_page)
    for (i in seq_along(chunks)) {
      chunk_labels <- label_map[chunks[[i]]]
      p <- trend_long |>
        dplyr::filter(MP %in% chunks[[i]]) |>
        dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)) |>
        ggplot2::ggplot(ggplot2::aes(x = sample_plot, y = score, group = interaction(patient, summary_stat), linetype = summary_stat)) +
        ggplot2::geom_line(color = "grey25", linewidth = 0.55, na.rm = TRUE) +
        ggplot2::geom_point(ggplot2::aes(fill = sample), shape = 21, size = 2.9, color = "black") +
        ggplot2::geom_text(data = dplyr::filter(trend_annot_df, MP %in% chunks[[i]]) |> dplyr::mutate(display_label = factor(display_label, levels = chunk_labels)), ggplot2::aes(x = sample_plot, y = score, label = stat_label), inherit.aes = FALSE, size = 2.45, lineheight = 0.92, fontface = "bold") +
        ggplot2::scale_fill_manual(values = sample_cols, guide = "none") +
        ggplot2::scale_linetype_manual(values = c("Mean" = "solid", "Median" = "dashed"), name = NULL) +
        ggplot2::scale_x_discrete(labels = sample_plot_labels, drop = FALSE) +
        ggplot2::facet_wrap(~display_label, scales = "fixed", ncol = 4) +
        ggplot2::coord_cartesian(ylim = trend_y_limits) +
        ggplot2::labs(title = paste0("Mean and Median UCell Pair Trend, ", type_labels[[direction_name]]), x = NULL, y = "UCell score") +
        ggplot2::theme_classic(base_size = 18) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 22),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 9),
          axis.text.y = ggplot2::element_text(colour = "black", size = 10),
          strip.text = ggplot2::element_text(face = "bold", size = 10),
          legend.position = "top",
          legend.text = ggplot2::element_text(size = 18)
        )
      print(p)
    }
  }
  dev.off()

  mat_df <- retained_summary |>
    dplyr::select(MP, sample, mean_score) |>
    tidyr::pivot_wider(names_from = sample, values_from = mean_score) |>
    as.data.frame()
  row.names(mat_df) <- mat_df$MP
  mat <- as.matrix(mat_df[, sample_order, drop = FALSE])
  ordered_rows <- unlist(lapply(type_levels, function(direction_name) {
    direction_mps <- trend_summary$MP[trend_summary$retained & trend_summary$treatment_direction == direction_name]
    direction_mps <- intersect(direction_mps, rownames(mat))
    if (length(direction_mps) <= 2) return(direction_mps)
    hc <- hclust(dist(mat[direction_mps, , drop = FALSE]), method = "ward.D2")
    direction_mps[hc$order]
  }), use.names = FALSE)
  mat <- mat[ordered_rows, , drop = FALSE]

  annotation_row <- trend_summary |>
    dplyr::filter(MP %in% rownames(mat)) |>
    dplyr::mutate(Pair_support = paste0(pair_support_n, "/4")) |>
    dplyr::select(MP, Direction = treatment_direction, Pair_support) |>
    as.data.frame()
  row.names(annotation_row) <- annotation_row$MP
  annotation_row <- annotation_row[rownames(mat), c("Direction", "Pair_support"), drop = FALSE]
  row_gaps <- which(annotation_row$Direction[-nrow(annotation_row)] != annotation_row$Direction[-1])
  anno_colors <- list(
    Direction = c(increase = "#4a4a4a", decrease = "#a0a0a0"),
    Pair_support = c("3/4" = "#c2855a", "4/4" = "#2b6a8e")
  )
  heat_cols <- grDevices::colorRampPalette(c("white", "#fee0d2", "#fc9272", "#de2d26", "#67000d"))(100)
  heat_breaks <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = length(heat_cols) + 1)
  if (length(unique(as.numeric(mat))) == 1) {
    heat_breaks <- seq(0, max(mat, na.rm = TRUE) + 1e-6, length.out = length(heat_cols) + 1)
  }
  labels_row <- label_map[rownames(mat)]
  labels_col <- sample_plot_labels[sample_order]

  pdf(
    file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mean_activity_heatmap_nMP", nMP, ".pdf")),
    width = 9.5,
    height = max(6, 0.32 * nrow(mat) + 2.5),
    useDingbats = FALSE
  )
  pheatmap::pheatmap(
    mat,
    color = heat_cols,
    breaks = heat_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = row_gaps,
    gaps_col = c(2, 4, 6),
    annotation_row = annotation_row,
    annotation_colors = anno_colors,
    labels_row = labels_row,
    labels_col = labels_col,
    border_color = "white",
    fontsize = 18,
    fontsize_row = 10,
    fontsize_col = 11,
    main = "Selected MP Mean UCell Activity"
  )
  dev.off()

  png(
    file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mean_activity_heatmap_nMP", nMP, ".png")),
    width = 2850,
    height = max(1800, 96 * nrow(mat) + 750),
    res = 300
  )
  pheatmap::pheatmap(
    mat,
    color = heat_cols,
    breaks = heat_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = row_gaps,
    gaps_col = c(2, 4, 6),
    annotation_row = annotation_row,
    annotation_colors = anno_colors,
    labels_row = labels_row,
    labels_col = labels_col,
    border_color = "white",
    fontsize = 16,
    fontsize_row = 10,
    fontsize_col = 11,
    main = "Selected MP Mean UCell Activity"
  )
  dev.off()
}

####################
# metaprograms
####################
if (!file.exists(geneNMF_program_path)) {
  stop("Missing GeneNMF programmes: ", geneNMF_program_path, ". Run Auto_pdo_flot_matched_geneNMF.R first.")
}

message("Loading GeneNMF programmes: ", geneNMF_program_path)
geneNMF.programs <- readRDS(geneNMF_program_path)
program_counts <- vapply(geneNMF.programs, program_count, numeric(1))
total_nmf_programs <- sum(program_counts, na.rm = TRUE)
nMP <- as.integer(round(total_nmf_programs / 2))
if (nMP < 2 || nMP >= total_nmf_programs) {
  stop("Calculated invalid nMP=", nMP, " from total NMF programmes=", total_nmf_programs)
}

write.csv(
  data.frame(
    total_nmf_programs = total_nmf_programs,
    nMP = nMP,
    programmes_per_MP_target = total_nmf_programs / nMP,
    number_samples = length(sample_order),
    samples = paste(sample_order, collapse = ","),
    nMP_rule = "round(total_nmf_programs / 2)",
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, paste0("Auto_pdo_flot_highres_nMP", nMP, "_config.csv")),
  row.names = FALSE
)
message("Using nMP=", nMP, " from ", total_nmf_programs, " total NMF programmes.")

metaprogram_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_geneNMF_metaprograms_nMP", nMP, ".rds"))
if (file.exists(metaprogram_path) && !force) {
  message("Loading existing high-resolution metaprograms: ", metaprogram_path)
  geneNMF.metaprograms <- readRDS(metaprogram_path)
} else {
  message("Running getMetaPrograms for high-resolution nMP=", nMP)
  geneNMF.metaprograms <- GeneNMF::getMetaPrograms(
    geneNMF.programs,
    metric = "cosine",
    specificity.weight = 5,
    weight.explained = 0.5,
    nMP = nMP,
    min.confidence = 0.5
  )
  saveRDS(geneNMF.metaprograms, metaprogram_path, compress = FALSE)
}

mp_genes <- geneNMF.metaprograms$metaprograms.genes
write_mp_gene_table(mp_genes, file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_genes_nMP", nMP, ".csv")))
make_3ca_label_table(mp_genes, nMP)

program_membership <- data.frame(
  nmf_programme = names(geneNMF.metaprograms$programs.clusters),
  MP = paste0("MP", as.integer(geneNMF.metaprograms$programs.clusters)),
  sample = sub("\\.k[0-9]+.*$", "", names(geneNMF.metaprograms$programs.clusters)),
  k = suppressWarnings(as.integer(sub("^.*\\.k([0-9]+).*$", "\\1", names(geneNMF.metaprograms$programs.clusters)))),
  stringsAsFactors = FALSE
)
write.csv(program_membership, file.path(out_dir, paste0("Auto_pdo_flot_highres_program_membership_nMP", nMP, ".csv")), row.names = FALSE)

heatmap_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_metaprograms_heatmap_nMP", nMP, ".png"))
if (!file.exists(heatmap_path) || force) {
  base_colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
  anno_colors <- grDevices::colorRampPalette(base_colors)(length(mp_genes))
  names(anno_colors) <- names(mp_genes)
  png(heatmap_path, width = 3000, height = 2500, res = 300)
  GeneNMF::plotMetaPrograms(
    geneNMF.metaprograms,
    annotation_colors = anno_colors,
    similarity.cutoff = c(0, 1)
  )
  dev.off()
}

####################
# UCell scoring
####################
sample_files <- setNames(
  file.path("by_samples", sample_order, paste0(sample_order, ".rds")),
  sample_order
)
missing_files <- names(sample_files)[!file.exists(sample_files)]
if (length(missing_files) > 0) {
  stop("Missing sample RDS file(s): ", paste(missing_files, collapse = ", "))
}

ucell_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_UCell_scores_nMP", nMP, ".rds"))
cell_meta_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_cell_metadata_nMP", nMP, ".rds"))
if (file.exists(ucell_path) && file.exists(cell_meta_path) && !force) {
  message("Loading existing UCell scores: ", ucell_path)
  ucell_scores <- readRDS(ucell_path)
  cell_meta <- readRDS(cell_meta_path)
} else {
  message("Finding common genes across matched PDO samples.")
  sample_genes <- lapply(sample_files, function(path) {
    obj <- readRDS(path)
    genes <- rownames(obj)
    rm(obj)
    gc()
    genes
  })
  common_genes <- Reduce(intersect, sample_genes)
  if (length(common_genes) == 0) {
    stop("No common genes found across matched PDO samples.")
  }

  counts_list <- list()
  meta_list <- list()
  for (sample in sample_order) {
    message("Loading counts for ", sample)
    obj <- readRDS(sample_files[[sample]])
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

  score_features <- lapply(mp_genes, intersect, rownames(counts_all))
  feature_counts <- lengths(score_features)
  write.csv(
    data.frame(MP = names(feature_counts), scored_gene_count = as.integer(feature_counts), stringsAsFactors = FALSE),
    file.path(out_dir, paste0("Auto_pdo_flot_highres_scored_gene_counts_nMP", nMP, ".csv")),
    row.names = FALSE
  )
  score_features <- score_features[feature_counts > 0]
  if (length(score_features) == 0) {
    stop("No high-resolution MP gene sets had genes present in the count matrix.")
  }

  message("Scoring ", length(score_features), " high-resolution MP signatures with UCell.")
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
  rm(counts_all, counts_list)
  gc()
}

####################
# paired trend filter
####################
mp_order <- intersect(names(mp_genes), colnames(ucell_scores))
if (length(mp_order) == 0) {
  stop("UCell score matrix does not contain any expected MP columns.")
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
write.csv(sample_summary, file.path(out_dir, paste0("Auto_pdo_flot_highres_sample_ucell_summary_nMP", nMP, ".csv")), row.names = FALSE)

pair_delta_summary <- sample_summary |>
  dplyr::select(MP, patient, treatment, mean_score, median_score) |>
  tidyr::pivot_wider(names_from = treatment, values_from = c(mean_score, median_score)) |>
  dplyr::mutate(
    mean_delta = mean_score_Treated - mean_score_Untreated,
    median_delta = median_score_Treated - median_score_Untreated
  )
write.csv(pair_delta_summary, file.path(out_dir, paste0("Auto_pdo_flot_highres_paired_delta_summary_nMP", nMP, ".csv")), row.names = FALSE)

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
trend_summary$retained <- trend_summary$treatment_direction %in% c("increase", "decrease")
trend_summary$trend_p_value <- pmax(trend_summary$mean_wilcox_p, trend_summary$median_wilcox_p, na.rm = TRUE)
trend_summary$trend_p_value[!is.finite(trend_summary$trend_p_value)] <- NA_real_
trend_summary$mean_wilcox_p_adj <- stats::p.adjust(trend_summary$mean_wilcox_p, method = "BH")
trend_summary$median_wilcox_p_adj <- stats::p.adjust(trend_summary$median_wilcox_p, method = "BH")
trend_summary$trend_p_adj <- pmax(trend_summary$mean_wilcox_p_adj, trend_summary$median_wilcox_p_adj, na.rm = TRUE)
trend_summary$trend_p_adj[!is.finite(trend_summary$trend_p_adj)] <- NA_real_
trend_summary <- add_delta_columns(trend_summary, "mean")
trend_summary <- add_delta_columns(trend_summary, "median")

metrics <- as.data.frame(geneNMF.metaprograms$metaprograms.metrics)
metrics$MP <- rownames(metrics)
if (is.null(metrics$MP) || any(!grepl("^MP", metrics$MP))) {
  metrics$MP <- paste0("MP", seq_len(nrow(metrics)))
}
metrics$numberPrograms <- suppressWarnings(as.integer(metrics$numberPrograms))
scored_counts_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_scored_gene_counts_nMP", nMP, ".csv"))
scored_counts <- if (file.exists(scored_counts_path)) {
  read.csv(scored_counts_path, check.names = FALSE)
} else {
  data.frame(MP = names(mp_genes), scored_gene_count = NA_integer_, stringsAsFactors = FALSE)
}
scored_counts$original_gene_count <- lengths(mp_genes[scored_counts$MP])

trend_summary <- trend_summary |>
  dplyr::left_join(metrics, by = "MP") |>
  dplyr::left_join(scored_counts[, c("MP", "original_gene_count", "scored_gene_count")], by = "MP") |>
  dplyr::arrange(
    dplyr::desc(retained),
    factor(treatment_direction, levels = type_levels),
    trend_p_value,
    dplyr::desc(pair_support_n),
    dplyr::desc(mean_delta_min_abs),
    dplyr::desc(numberPrograms),
    dplyr::desc(silhouette),
    MP
  )
write.csv(trend_summary, file.path(out_dir, paste0("Auto_pdo_flot_highres_trend_summary_nMP", nMP, ".csv")), row.names = FALSE)

retained_mps <- order_retained_by_type(trend_summary)
retained_genes <- mp_genes[retained_mps]
saveRDS(retained_genes, file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mp_genes_nMP", nMP, ".rds")), compress = FALSE)
write_mp_gene_table(retained_genes, file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mp_genes_nMP", nMP, ".csv")))
write.csv(
  dplyr::filter(trend_summary, retained),
  file.path(out_dir, paste0("Auto_pdo_flot_highres_trend_retained_nMP", nMP, ".csv")),
  row.names = FALSE
)

message("Retained ", length(retained_mps), " MPs after strict mean+median paired trend filtering with at least 3/4 supporting pairs.")
plot_retained_outputs(ucell_scores, cell_meta, sample_summary, trend_summary, nMP)

run_excel <- function() {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("Skipping excel generation because openxlsx is missing.")
    return(invisible(NULL))
  }
  
  trend_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_trend_summary_nMP", nMP, ".csv"))
  genes_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mp_genes_nMP", nMP, ".rds"))
  
  if (!file.exists(trend_path) || !file.exists(genes_path)) stop("Missing selected MP outputs.")
  
  trend_summary <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
  mp_genes <- readRDS(genes_path)
  
  retained_df <- trend_summary[trend_summary$retained == TRUE | trend_summary$retained == "TRUE", ]
  unique_labels <- unique(retained_df$treatment_direction)
  unique_labels <- unique_labels[!is.na(unique_labels) & unique_labels %in% names(type_labels)]
  
  if (length(unique_labels) == 0) {
    message("No retained categories found for excel generation.")
    return(invisible(NULL))
  }
  
  build_mp_matrix <- function(mp_names_vec) {
    if (length(mp_names_vec) == 0) return(NULL)
    max_g <- max(sapply(mp_names_vec, function(x) length(mp_genes[[x]])))
    n_mp <- length(mp_names_vec)
    
    n_rows <- max_g + 2
    mat <- matrix(NA_character_, nrow = n_rows, ncol = n_mp)
    for (i in seq_along(mp_names_vec)) {
      mp <- mp_names_vec[i]
      n_prog <- retained_df$numberPrograms[retained_df$MP == mp]
      if (length(n_prog) > 0 && !is.na(n_prog[1])) {
        prog_label <- ifelse(n_prog[1] == 1, "1 programme", paste0(n_prog[1], " programmes"))
        mat[1, i] <- paste0(mp, " (", prog_label, ")")
      } else {
        mat[1, i] <- mp
      }
      mat[2, i] <- "" 
      genes <- mp_genes[[mp]]
      if (length(genes) > 0) {
        mat[3:(length(genes)+2), i] <- genes
      }
    }
    return(as.data.frame(mat, stringsAsFactors = FALSE))
  }
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "MP_Genes")
  
  current_col <- 1
  section_style <- openxlsx::createStyle(fontSize = 14, textDecoration = "bold", fgFill = "#FFC000")
  mp_name_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  desc_style <- openxlsx::createStyle(fgFill = "#F2F2F2")
  
  unique_labels <- unique_labels[order(match(unique_labels, type_levels))]
  
  for (label in unique_labels) {
    mps_in_group <- order_mps_by_metrics(
      trend_summary$MP[trend_summary$retained & trend_summary$treatment_direction == label],
      trend_summary
    )
    mps_in_group <- mps_in_group[mps_in_group %in% names(mp_genes)]
    
    if (length(mps_in_group) == 0) next
    
    df_group <- build_mp_matrix(mps_in_group)
    
    display_label <- type_labels[[label]]
    openxlsx::writeData(wb, sheet = 1, x = toupper(display_label), startCol = current_col, startRow = 1)
    openxlsx::writeData(wb, sheet = 1, x = df_group, startCol = current_col, startRow = 2, colNames = FALSE)
    
    openxlsx::addStyle(wb, sheet = 1, section_style, rows = 1, cols = current_col, gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet = 1, mp_name_style, rows = 2, cols = current_col:(current_col + ncol(df_group) - 1), gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet = 1, desc_style, rows = 3, cols = current_col:(current_col + ncol(df_group) - 1), gridExpand = TRUE)
    
    for (i in current_col:(current_col + ncol(df_group) - 1)) {
      openxlsx::setColWidths(wb, 1, cols = i, widths = 25)
    }
    
    current_col <- current_col + ncol(df_group) + 1
  }
  
  output_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_mp_genes_summary_nMP", nMP, ".xlsx"))
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Excel summary written to: ", output_path)
}

run_excel()

message("Auto_pdo_flot_matched_highres_mp_trend_filter.R completed successfully.")
