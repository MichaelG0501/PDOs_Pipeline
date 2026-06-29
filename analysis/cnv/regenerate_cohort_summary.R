# Read all required libs
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
  library(ggrepel)
  library(patchwork)
})
root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)
out_dir <- "Auto_PDO_numbat_subclone_mp_conservative"

# Read intermediates
cell_df <- fread(file.path(out_dir, "Auto_PDO_numbat_subclone_cells.csv"))
sample_df <- fread(file.path(out_dir, "Auto_PDO_numbat_subclone_summary.csv"))
mp_tests_df <- fread(file.path(out_dir, "Auto_PDO_numbat_subclone_mp_tests.csv"))
sub_tests_df <- fread(file.path(out_dir, "Auto_PDO_numbat_subclone_mp_subclone_tests.csv"))

# Force a strict 1.00 threshold as requested
sub_tests_df <- sub_tests_df %>%
  mutate(
    score_threshold = 1.00,
    score_threshold_pass = abs_difference_score >= 1.00,
    significant = statistically_significant & score_threshold_pass,
    threshold_method = "Fixed 1.00"
  )

state_tests_df <- fread(file.path(out_dir, "Auto_PDO_numbat_subclone_state_tests.csv"))
clone_comp_df <- fread(file.path(out_dir, "Auto_PDO_numbat_subclone_compartment_summary.csv"))

# We also need ucell_scores and target_mps from the main script
ucell_scores <- readRDS("UCell_scores_filtered.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")
mp_descriptions <- c("MP6" = "G2M Cell Cycle", "MP7" = "DNA repair", "MP5" = "MYC-related Proliferation", "MP1" = "G2M checkpoint", "MP3" = "G1S Cell Cycle", "MP8" = "Columnar Progenitor", "MP10" = "Inflammatory Stress Epi.", "MP9" = "ECM Remodeling Epi.", "MP4" = "Intestinal Metaplasia")
mp_cols <- c(
  "MP6_G2M Cell Cycle" = "#B0B0B0", "MP7_DNA repair" = "#999999", "MP1_G2M checkpoint" = "#808080",
  "MP3_G1S Cell Cycle" = "#C0C0C0", "MP5_MYC-related Proliferation" = "#E41A1C",
  "MP4_Intestinal Metaplasia" = "#4DAF4A", "MP8_Columnar Progenitor" = "#FF7F00",
  "MP10_Inflammatory Stress Epi." = "#984EA3", "MP9_ECM Remodeling Epi." = "#C77CFF", "Unassigned" = "grey85"
)
state_cols <- c(
  "Classic Proliferative" = "#E41A1C", "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00", "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8", "Unresolved" = "grey80",
  "Hybrid" = "black", "Unassigned" = "grey85"
)
label_mp <- function(mps) {
  desc <- mp_descriptions[mps]
  desc[is.na(desc)] <- mps[is.na(desc)]
  paste0(mps, "_", desc)
}
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
if (!is.null(coverage_tbl)) {
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
  mp.genes <- mp.genes[!names(mp.genes) %in% low_coverage_mps]
}
retained_mps <- names(mp.genes)
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- paste0("MP", rev(unique(ordered_clusters)))
mp_tree_order <- mp_tree_order[mp_tree_order %in% retained_mps]
state_groups <- list("Classic Proliferative" = c("MP5"), "Basal to Intest. Meta" = c("MP4"), "SMG-like Metaplasia" = c("MP8"), "Stress-adaptive" = c("MP10", "MP9"))
ordered_state_mps <- unlist(lapply(state_groups, function(mps) { x <- mp_tree_order[mp_tree_order %in% mps]; c(x, setdiff(mps, x)) }), use.names = FALSE)
cc_mps <- c("MP6", "MP7", "MP1", "MP3")
mp_names <- unique(c(cc_mps[cc_mps %in% retained_mps], ordered_state_mps, mp_tree_order))
mp_names <- mp_names[mp_names %in% names(mp_descriptions)]
mp_labels <- setNames(label_mp(mp_names), mp_names)

make_palette <- function(values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(character(0))
  base <- suppressWarnings(brewer.pal(max(3, min(12, length(values))), palette))
  setNames(colorRampPalette(base)(length(values)), values)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values[is.na(values) | !nzchar(values)] <- "Unassigned"
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(c(Unassigned = "grey85"))
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  missing_values <- setdiff(values, names(cols))
  if (length(missing_values) > 0) cols <- c(cols, make_palette(missing_values, palette))
  out <- cols[values]
  names(out) <- values
  out
}

nullGrob <- function() {
  grid::nullGrob()
}
if (nrow(mp_tests_df) > 0) {
  multi_clone_samples <- sample_df$sample[sample_df$n_display_clones >= 2]
  sig_counts_pair <- sub_tests_df %>%
    filter(.data$sample %in% multi_clone_samples) %>%
    group_by(.data$sample, .data$pair_label) %>%
    summarise(n_sig_mps = sum(.data$significant, na.rm = TRUE), .groups = "drop") %>%
    mutate(category = case_when(
      n_sig_mps == 0 ~ "None",
      n_sig_mps == 1 ~ "One significant",
      TRUE ~ "More than one"
    ))
  mp_summary_sample <- sub_tests_df %>%
    group_by(.data$mp, .data$mp_label) %>%
    summarise(
      n_pairwise_tests = n(),
      n_statistically_significant_pairwise = sum(.data$statistically_significant, na.rm = TRUE),
      n_score_threshold_pairwise = sum(.data$score_threshold_pass, na.rm = TRUE),
      n_significant_pairwise = sum(.data$significant, na.rm = TRUE),
      pct_significant_pairwise = 100 * mean(.data$significant, na.rm = TRUE),
      median_abs_difference_score = median(.data$abs_difference_score, na.rm = TRUE),
      max_abs_difference_score = max(.data$abs_difference_score, na.rm = TRUE),
      score_threshold = dplyr::first(.data$score_threshold),
      .groups = "drop"
    )
  write.csv(sig_counts_pair, file.path(out_dir, "Auto_PDO_numbat_subclone_sig_count_summary.csv"), row.names = FALSE)
  write.csv(mp_summary_sample, file.path(out_dir, "Auto_PDO_numbat_subclone_mp_cohort_summary.csv"), row.names = FALSE)

  p_counts <- sig_counts_pair %>%
    count(.data$category, name = "n") %>%
    mutate(category = factor(.data$category, levels = c("None", "One significant", "More than one")),
           pct = 100 * .data$n / sum(.data$n)) %>%
    ggplot(aes(.data$category, .data$pct, fill = .data$category)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = paste0(.data$n, " (", round(.data$pct, 1), "%)")), vjust = -0.4, size = 4) +
    scale_fill_manual(values = c("None" = "grey70", "One significant" = "#FDB863", "More than one" = "#B2182B")) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = "Significant MP differences per clone pair", x = NULL, y = "Percentage of clone pairs") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  target_mps <- mp_names[mp_names %in% names(mp_labels)]
  mp_plot_df <- mp_tests_df %>%
    filter(.data$sample %in% multi_clone_samples, .data$mp %in% target_mps) %>%
    mutate(mp_label = factor(.data$mp_label, levels = mp_labels[target_mps]),
           val = -log10(pmax(.data$p_adj, .Machine$double.xmin)),
           val_plot = pmax(.data$val, 1e-3))
  mp_pcts <- mp_plot_df %>%
    group_by(.data$mp_label) %>%
    summarise(pct = 100 * mean(.data$significant, na.rm = TRUE), .groups = "drop")
  p_mp <- ggplot(mp_plot_df, aes(.data$mp_label, .data$val_plot, fill = .data$mp_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = .data$significant), width = 0.2, size = 1, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = mp_pcts, aes(x = .data$mp_label, y = max(mp_plot_df$val_plot, na.rm = TRUE) * 1.15, label = sprintf("%.1f%%", .data$pct)), inherit.aes = FALSE, size = 3) +
    scale_y_log10(expand = expansion(mult = c(0.1, 0.3))) +
    scale_fill_manual(values = complete_palette(mp_cols, mp_labels[target_mps], "Paired")) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey35")) +
    labs(title = "MP association with Numbat clones", subtitle = "Points are final calls requiring FDR < 0.05 and MP-specific difference-score threshold", x = NULL, y = "-log10(BH p)", color = "Final significant") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

  p_state <- state_tests_df %>%
    filter(.data$sample %in% multi_clone_samples) %>%
    mutate(sample = factor(.data$sample, levels = .data$sample[order(.data$cramers_v, decreasing = TRUE)])) %>%
    ggplot(aes(.data$cramers_v, .data$sample, color = .data$p_adj < 0.05)) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey50"), na.value = "grey70") +
    labs(title = "State abundance association with Numbat clones", x = "Cramer's V", y = NULL, color = "BH p < 0.05") +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_text(size = 7), legend.position = "top")

  score_plot_df <- sub_tests_df %>%
    filter(.data$sample %in% multi_clone_samples, .data$mp %in% target_mps) %>%
    mutate(mp_label = factor(.data$mp_label, levels = rev(mp_labels[target_mps])))
  score_label_df <- score_plot_df %>%
    group_by(.data$mp, .data$mp_label) %>%
    summarise(
      label_x = max(abs(.data$difference_score), .data$score_threshold, na.rm = TRUE) + 0.25,
      score_threshold = dplyr::first(.data$score_threshold),
      n_sig = sum(.data$significant, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("thr=", sprintf("%.2f", .data$score_threshold), "; sig=", .data$n_sig))
  p_score <- ggplot(score_plot_df, aes(.data$difference_score, .data$mp_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
    geom_boxplot(width = 0.52, outlier.shape = NA, fill = "grey86", color = "black", linewidth = 0.55) +
    geom_point(aes(color = .data$significant, shape = .data$statistically_significant),
               position = position_jitter(height = 0.12, width = 0), alpha = 0.65, size = 1.8) +
    geom_text(data = score_label_df, aes(x = .data$label_x, y = .data$mp_label, label = .data$label),
              inherit.aes = FALSE, hjust = 0, size = 3.1) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey45")) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
    labs(title = "Pairwise subclone MP expression difference scores",
         subtitle = "Score = pairwise median difference / sample MP MAD; final calls also require BH FDR < 0.05",
         x = "Median difference / MAD",
         y = NULL,
         color = "Final significant",
         shape = "FDR < 0.05") +
    theme_classic(base_size = 14) +
    theme(axis.text.y = element_text(size = 10), legend.position = "top",
          plot.title = element_text(face = "bold"))

  # ---------------------------------------------------------
  # Multi-page Summary PDF (Mirrors scRef Logic)
  # ---------------------------------------------------------

  analysed_df <- sample_df
  subclone_count_df <- analysed_df %>%
    mutate(n_subclones_cat = case_when(
      n_display_clones == 1 ~ "1",
      n_display_clones == 2 ~ "2",
      n_display_clones == 3 ~ "3",
      n_display_clones >= 4 ~ "4+"
    )) %>%
    count(n_subclones_cat, name = "n_samples") %>%
    mutate(n_subclones_cat = factor(n_subclones_cat, levels = c("1", "2", "3", "4+")),
           pct = 100 * n_samples / sum(n_samples))

  neutral_cols <- c("1" = "#C6DBEF", "2" = "#9ECAE1", "3" = "#4292C6", "4+" = "#084594")

  p_subclone_dist <- ggplot(subclone_count_df, aes(n_subclones_cat, n_samples, fill = n_subclones_cat)) +
    geom_col(color = "black", linewidth = 0.4, width = 0.6) +
    geom_text(aes(label = paste0(n_samples, " (", round(pct, 1), "%)")), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = neutral_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = "Numbat clones per sample", subtitle = paste0("n = ", sum(subclone_count_df$n_samples), " analysed samples"), x = "Number of clones", y = "Samples") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9), plot.margin = margin(4, 6, 4, 6))

  main_states <- c("Classic Proliferative", "Basal to Intestinal Metaplasia", "SMG-like Metaplasia", "Stress-adaptive", "Immune Infiltrating")
  mp_ucell_thresh <- 0.10

  valid_cells_ms <- intersect(cell_df$cell, rownames(ucell_scores))
  cell_df_ms <- cell_df[match(valid_cells_ms, cell_df$cell), ]
  ucell_mps <- intersect(target_mps, colnames(ucell_scores))
  ucell_ms <- as.matrix(ucell_scores[valid_cells_ms, ucell_mps, drop = FALSE])

  excl_res <- list()
  for (samp in multi_clone_samples) {
    samp_mask <- cell_df_ms$sample == samp
    samp_cells <- cell_df_ms$cell[samp_mask]
    samp_subs <- cell_df_ms$subclone[samp_mask]
    if (length(samp_cells) < 10) next
    
    any_excl_mp <- FALSE
    for (mp in ucell_mps) {
      sc_means <- tapply(ucell_ms[samp_cells, mp], samp_subs, mean, na.rm = TRUE)
      if (sum(sc_means > mp_ucell_thresh, na.rm = TRUE) == 1) { any_excl_mp <- TRUE; break }
    }
    
    any_excl_state <- FALSE
    samp_states <- cell_df_ms$state_label[samp_mask]
    for (st in main_states) {
      st_presence <- tapply(samp_states == st, samp_subs, any, na.rm = TRUE)
      if (sum(st_presence, na.rm = TRUE) == 1) { any_excl_state <- TRUE; break }
    }
    excl_res[[samp]] <- data.frame(sample = samp, has_excl_mp = any_excl_mp, has_excl_state = any_excl_state)
  }
  excl_df <- bind_rows(excl_res)

  if(nrow(excl_df) > 0) {
    mp_excl_summary <- excl_df %>%
      count(category = ifelse(has_excl_mp, "Has exclusive MP(s)", "All MPs shared")) %>%
      mutate(category = factor(category, levels = c("Has exclusive MP(s)", "All MPs shared")), pct = 100 * n / sum(n))
    p_mp_excl <- ggplot(mp_excl_summary, aes(category, pct, fill = category)) +
      geom_col(color = "black", linewidth = 0.4, width = 0.5) +
      geom_text(aes(label = paste0(n, " (", round(pct, 1), "%)")), vjust = -0.4, size = 3.5) +
      scale_fill_manual(values = c("Has exclusive MP(s)" = "#08519C", "All MPs shared" = "#9ECAE1")) +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.12))) +
      labs(title = "MP exclusivity", subtitle = paste0("n = ", nrow(excl_df), " multi-clone; UCell > ", mp_ucell_thresh), x = NULL, y = "% samples") +
      theme_classic(base_size = 11) +
      theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9), plot.margin = margin(4, 6, 4, 6))

    st_excl_summary <- excl_df %>%
      count(category = ifelse(has_excl_state, "Has exclusive State(s)", "All States shared")) %>%
      mutate(category = factor(category, levels = c("Has exclusive State(s)", "All States shared")), pct = 100 * n / sum(n))
    p_st_excl <- ggplot(st_excl_summary, aes(category, pct, fill = category)) +
      geom_col(color = "black", linewidth = 0.4, width = 0.5) +
      geom_text(aes(label = paste0(n, " (", round(pct, 1), "%)")), vjust = -0.4, size = 3.5) +
      scale_fill_manual(values = c("Has exclusive State(s)" = "#08519C", "All States shared" = "#9ECAE1")) +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.12))) +
      labs(title = "State exclusivity", subtitle = paste0("n = ", nrow(excl_df), " multi-clone; 5 main states"), x = NULL, y = "% samples") +
      theme_classic(base_size = 11) +
      theme(legend.position = "none", plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9), plot.margin = margin(4, 6, 4, 6))
  } else {
    p_mp_excl <- nullGrob()
    p_st_excl <- nullGrob()
  }

  sig_counts_pair_no_cc <- sub_tests_df %>%
    filter(.data$sample %in% multi_clone_samples, !(.data$mp %in% cc_mps)) %>%
    group_by(.data$sample, .data$pair_label) %>%
    summarise(n_sig_mps = sum(.data$significant, na.rm = TRUE), .groups = "drop") %>%
    mutate(category = case_when(
      n_sig_mps == 0 ~ "None",
      n_sig_mps == 1 ~ "One significant",
      TRUE ~ "More than one"
    ))
  p_counts_no_cc <- sig_counts_pair_no_cc %>%
    count(.data$category, name = "n") %>%
    mutate(category = factor(.data$category, levels = c("None", "One significant", "More than one")),
           pct = 100 * .data$n / sum(.data$n)) %>%
    ggplot(aes(.data$category, .data$pct, fill = .data$category)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_text(aes(label = paste0(.data$n, " (", round(.data$pct, 1), "%)")), vjust = -0.4, size = 4) +
    scale_fill_manual(values = c("None" = "grey70", "One significant" = "#FDB863", "More than one" = "#B2182B")) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = "Significant MP diff. per pair (excl. CC MPs)", subtitle = "Excludes G2M, G1S, DNA Repair", x = NULL, y = "Percentage of clone pairs") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face="bold"))

  plot_list_qc <- list(nullGrob())
  plot_list_mp <- list(nullGrob())
  plot_list_state <- list(nullGrob())

  qc_summary_rows <- list()
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "cna_signal")
  for (samp in multi_clone_samples) {
    samp_df <- cell_df %>% filter(sample == samp)
    top_subs <- samp_df %>% count(subclone) %>% arrange(desc(n)) %>% head(2) %>% pull(subclone)
    if (length(top_subs) < 2) next
    row_data <- data.frame(sample = samp, clone1 = top_subs[1], clone2 = top_subs[2])
    for (q in qc_metrics) {
      if (q %in% colnames(samp_df)) {
        cl_means <- samp_df %>% group_by(subclone) %>% summarise(m = mean(.data[[q]], na.rm = TRUE), .groups = "drop")
        v1 <- max(cl_means$m, na.rm = TRUE)
        v2 <- min(cl_means$m, na.rm = TRUE)
        row_data[[paste0("X_", q)]] <- v1
        row_data[[paste0("Y_", q)]] <- v2
      }
    }
    qc_summary_rows[[length(qc_summary_rows) + 1]] <- row_data
  }
  if(length(qc_summary_rows) > 0) {
    qc_summary_df <- bind_rows(qc_summary_rows)
    plot_list_qc <- list()
    for (q in qc_metrics) {
      if (paste0("X_", q) %in% colnames(qc_summary_df)) {
        x_col <- paste0("X_", q)
        y_col <- paste0("Y_", q)
        wt <- tryCatch(wilcox.test(qc_summary_df[[x_col]], qc_summary_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
        p_val <- wt$p.value
        diff_stat <- mean(qc_summary_df[[x_col]] - qc_summary_df[[y_col]], na.rm = TRUE)
        p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val)
        subtitle <- sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat)
        qc_all_vals <- unlist(qc_summary_df[, c(x_col, y_col)])
        qc_lims <- quantile(qc_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
        p <- ggplot(qc_summary_df, aes(.data[[x_col]], .data[[y_col]])) +
          geom_point(alpha = 0.5, size = 1.2, color = "black") +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
          coord_cartesian(xlim = qc_lims, ylim = qc_lims) +
          labs(title = q, subtitle = subtitle, x = "Highest subclone", y = "Lowest subclone") +
          theme_classic(base_size = 9) + theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 7.5))
        plot_list_qc[[q]] <- p
      }
    }
  }

  valid_cells <- intersect(cell_df$cell, rownames(ucell_scores))
  cell_df_valid <- cell_df[match(valid_cells, cell_df$cell), ]
  mp_scores_valid <- ucell_scores[valid_cells, target_mps, drop = FALSE]
  
  subclone_means <- cell_df_valid %>%
    select(cell, sample, subclone) %>%
    bind_cols(as.data.frame(mp_scores_valid)) %>%
    filter(sample %in% multi_clone_samples) %>%
    group_by(sample, subclone) %>%
    summarise(across(all_of(target_mps), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    
  pair_rows <- list()
  for (samp in unique(subclone_means$sample)) {
    samp_df <- subclone_means %>% filter(sample == samp) %>% arrange(subclone)
    if (nrow(samp_df) >= 2) {
      combos <- combn(nrow(samp_df), 2)
      for (i in 1:ncol(combos)) {
        idx1 <- combos[1, i]
        idx2 <- combos[2, i]
        row_data <- data.frame(sample = samp, clone1 = samp_df$subclone[idx1], clone2 = samp_df$subclone[idx2])
        for (mp in target_mps) {
          val1 <- samp_df[[mp]][idx1]
          val2 <- samp_df[[mp]][idx2]
          row_data[[paste0("X_", mp)]] <- max(val1, val2, na.rm = TRUE)
          row_data[[paste0("Y_", mp)]] <- min(val1, val2, na.rm = TRUE)
        }
        pair_rows[[length(pair_rows) + 1]] <- row_data
      }
    }
  }
  if (length(pair_rows) > 0) {
    pairs_df <- bind_rows(pair_rows)
    mp_x_cols <- paste0("X_", target_mps)
    mp_y_cols <- paste0("Y_", target_mps)
    mp_all_vals <- unlist(pairs_df[, c(mp_x_cols, mp_y_cols)])
    mp_global_lims <- quantile(mp_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
    plot_list_mp <- list()
    for (mp in target_mps) {
      x_col <- paste0("X_", mp)
      y_col <- paste0("Y_", mp)
      wt <- tryCatch(wilcox.test(pairs_df[[x_col]], pairs_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
      p_val <- wt$p.value
      diff_stat <- mean(pairs_df[[x_col]] - pairs_df[[y_col]], na.rm = TRUE)
      title_label <- mp_labels[mp]
      p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val)
      subtitle <- sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat)
      mp_color <- if (!is.na(mp_cols[mp_labels[mp]])) mp_cols[mp_labels[mp]] else "grey50"
      p <- ggplot(pairs_df, aes(.data[[x_col]], .data[[y_col]])) +
        geom_point(alpha = 0.5, color = mp_color, size = 1.2) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
        labs(title = title_label, subtitle = subtitle, x = "Higher-expressed clone", y = "Lower-expressed clone") +
        coord_cartesian(xlim = mp_global_lims, ylim = mp_global_lims) +
        theme_classic(base_size = 9) + theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 7.5))
      plot_list_mp[[mp]] <- p
    }
  }

  target_states <- c("Classic Proliferative", "Basal to Intestinal Metaplasia", "SMG-like Metaplasia", "Stress-adaptive", "Immune Infiltrating", "3CA_EMT_and_Protein_maturation", "Unresolved", "Hybrid")
  state_counts <- cell_df %>%
    filter(sample %in% multi_clone_samples) %>%
    count(sample, subclone, state_label) %>%
    group_by(sample, subclone) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  pair_rows_state <- list()
  for (samp in unique(state_counts$sample)) {
    samp_df_st <- state_counts %>% filter(sample == samp)
    subs <- sort(unique(samp_df_st$subclone))
    if (length(subs) >= 2) {
      combos <- combn(length(subs), 2)
      for (i in 1:ncol(combos)) {
        sub1 <- subs[combos[1, i]]
        sub2 <- subs[combos[2, i]]
        row_data <- data.frame(sample = samp, clone1 = sub1, clone2 = sub2)
        for (st in target_states) {
          p1 <- samp_df_st$prop[samp_df_st$subclone == sub1 & samp_df_st$state_label == st]
          p2 <- samp_df_st$prop[samp_df_st$subclone == sub2 & samp_df_st$state_label == st]
          val1 <- if (length(p1) > 0) p1[1] else 0
          val2 <- if (length(p2) > 0) p2[1] else 0
          row_data[[paste0("X_", st)]] <- max(val1, val2, na.rm = TRUE)
          row_data[[paste0("Y_", st)]] <- min(val1, val2, na.rm = TRUE)
        }
        pair_rows_state[[length(pair_rows_state) + 1]] <- row_data
      }
    }
  }
  if (length(pair_rows_state) > 0) {
    pairs_state_df <- bind_rows(pair_rows_state)
    state_x_cols <- paste0("X_", target_states)
    state_y_cols <- paste0("Y_", target_states)
    state_all_vals <- unlist(pairs_state_df[, c(state_x_cols, state_y_cols)])
    state_global_lims <- quantile(state_all_vals, probs = c(0.01, 0.99), na.rm = TRUE)
    plot_list_state <- list()
    for (st in target_states) {
      x_col <- paste0("X_", st)
      y_col <- paste0("Y_", st)
      wt <- tryCatch(wilcox.test(pairs_state_df[[x_col]], pairs_state_df[[y_col]], paired = TRUE), error = function(e) list(p.value = NA_real_))
      p_val <- wt$p.value
      diff_stat <- mean(pairs_state_df[[x_col]] - pairs_state_df[[y_col]], na.rm = TRUE)
      st_col <- if (st %in% names(state_cols)) state_cols[[st]] else "grey50"
      p_val_display <- if (is.na(p_val)) "NA" else if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.3f", p_val)
      subtitle <- sprintf("p = %s | Diff = %.3f", p_val_display, diff_stat)
      p <- ggplot(pairs_state_df, aes(.data[[x_col]], .data[[y_col]])) +
        geom_point(alpha = 0.5, color = st_col, size = 1.2) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
        labs(title = st, subtitle = subtitle, x = "Higher-abundance clone", y = "Lower-abundance clone") +
        coord_cartesian(xlim = state_global_lims, ylim = state_global_lims) +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        theme_classic(base_size = 9) + theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 7.5))
      plot_list_state[[st]] <- p
    }
  }

  tmp_page0 <- tempfile(fileext = ".pdf")
  pdf(tmp_page0, width = 15, height = 4.2, useDingbats = FALSE)
  grid.arrange(p_subclone_dist, p_counts, p_mp_excl, p_st_excl, ncol = 4)
  dev.off()

  tmp_rest <- tempfile(fileext = ".pdf")
  pdf(tmp_rest, width = 15, height = 9, useDingbats = FALSE)
  grid.arrange(p_score, p_counts_no_cc, ncol = 2, widths = c(2, 1))
  grid.arrange(grobs = plot_list_qc, ncol = 3)
  grid.arrange(grobs = plot_list_mp, ncol = 4)
  grid.arrange(grobs = plot_list_state, ncol = 4)
  dev.off()

  final_summary_pdf <- file.path(out_dir, "Auto_PDO_numbat_subclone_mp_cohort_summary.pdf")
  if (requireNamespace("qpdf", quietly = TRUE)) {
    qpdf::pdf_combine(c(tmp_page0, tmp_rest), output = final_summary_pdf)
  } else {
    py_cmd <- sprintf("python -c 'import PyPDF2; m = PyPDF2.PdfMerger(); m.append(\"%s\"); m.append(\"%s\"); m.write(\"%s\"); m.close()'", tmp_page0, tmp_rest, final_summary_pdf)
    system(py_cmd)
  }
  unlink(c(tmp_page0, tmp_rest))

  proportion_pdf <- file.path(out_dir, "Auto_PDO_numbat_subclone_mp_proportion_plots.pdf")
  pdf(proportion_pdf, width = 12, height = 6, useDingbats = FALSE)
  grid.arrange(p_counts, p_counts_no_cc, ncol = 2)
  dev.off()

  ####################
  # Diagnostic PDF: one MP per page
  diagnostic_pdf <- file.path(out_dir, "Auto_PDO_numbat_subclone_mp_difference_diagnostics.pdf")
  pdf(diagnostic_pdf, width = 14, height = 8, useDingbats = FALSE)
  
  p_diag_all <- sub_tests_df %>%
    filter(.data$sample %in% multi_clone_samples, .data$mp %in% target_mps) %>%
    mutate(mp_label = factor(.data$mp_label, levels = mp_labels[target_mps])) %>%
    ggplot(aes(x = mp_label, y = abs_difference_score, fill = mp_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, aes(color = significant)) +
    scale_fill_manual(values = complete_palette(mp_cols, mp_labels[target_mps], "Paired")) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey55")) +
    labs(title = "Subclone MP Difference Scores", x = NULL, y = "Abs. Difference Score", color = "Significant") +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", plot.title = element_text(face = "bold"))
  print(p_diag_all)

  for (mp_id in target_mps) {
    diag_df <- sub_tests_df %>%
      filter(.data$sample %in% multi_clone_samples, .data$mp == mp_id) %>%
      arrange(.data$abs_difference_score) %>%
      mutate(comparison_rank = row_number(),
             comparison_label = paste(.data$sample, .data$pair_label, sep = ": "),
             call_class = case_when(
               .data$significant ~ "Final significant",
               .data$statistically_significant ~ "FDR only",
               .data$score_threshold_pass ~ "Score only",
               TRUE ~ "Not significant"
             ))
    if (nrow(diag_df) == 0) next
    threshold_value <- unique(diag_df$score_threshold)[1]
    
    set.seed(42)
    label_subset <- diag_df %>%
      group_by(call_class) %>%
      slice_sample(n = 1) %>%
      ungroup()
    
    p_diag <- ggplot(diag_df, aes(.data$comparison_rank, .data$abs_difference_score)) +
      geom_segment(aes(xend = .data$comparison_rank, y = 0, yend = .data$abs_difference_score, color = .data$call_class), linewidth = 0.6, alpha = 0.82) +
      geom_point(aes(color = .data$call_class), size = 2.5, alpha = 0.9) +
      geom_hline(yintercept = threshold_value, linetype = "dashed", color = "#B2182B", linewidth = 1) +
      geom_text_repel(data = label_subset, aes(label = .data$comparison_label, color = .data$call_class), size = 3.5, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.4, segment.alpha = 0.6, box.padding = 0.8, point.padding = 0.5, nudge_y = 0.4, direction = "both") +
      scale_color_manual(values = c("Final significant" = "#B2182B", "FDR only" = "#2166AC", "Score only" = "#F4A582", "Not significant" = "grey55")) +
      labs(title = mp_labels[mp_id],
           x = "Pairwise comparison rank", y = "Abs. Difference Score", color = "Call class") +
      theme_classic(base_size = 16) +
      theme(legend.position = "top", plot.title = element_text(face = "bold"))
    print(p_diag)
  }
  dev.off()
}
