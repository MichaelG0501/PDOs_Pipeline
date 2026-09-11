####################
# Analysis registry:
#   Status: active terminal; selected centred high-resolution MP TCGA survival
#   Script: analysis/cell_states/Auto_pdo_flot_matched_survival_and_state_plots.R
#   Methodology: analysis/methodology/cell_states/Auto_pdo_flot_centred_highres_metaprogram_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Applies the TCGA EAC GSVA/Cox sensitivity analysis from centred
#     Auto_08_tcga_mp_survival_volcano_centred.R to the diverse high-resolution
#     MPs retained by the matched untreated-versus-FLOT trend filter. Features
#     are shown with their best non-cell-cycle 3CA enrichment labels. No state
#     definition or manual MP grouping is performed.
#   Inputs:
#     - selected high-resolution MP genes and best non-cell-cycle labels
#     - scRef reconstructed TCGA ESCA metadata and TPM matrix
#   Outputs:
#     - live survival CSV, GSVA score RDS and three-panel volcano PDF
#   Downstream use: none; terminal association analysis.
#   Cache/replot behavior: rebuilds deterministically from persistent inputs.
#   Volcano labels:
#     Labels all nominally significant MPs (P_value < 0.05) and five
#     reproducibly sampled non-significant MPs from the 15 closest to the
#     nominal significance boundary in each split-method panel.
#   Run command: use PBS in dmtcp after the trend filter.
#   Conda env: dmtcp
####################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(survival)
library(GSVA)

####################
# Persistent paths and selected high-resolution MPs
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

out_dir <- file.path(PDO_LIVE_OUTS, "Auto_pdo_flot_centred_highres_metaprogram_trends")
survival_dir <- file.path(out_dir, "survival")
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)

config_path <- file.path(out_dir, "Auto_pdo_flot_highres_current_config.csv")
pdo_require_files(config_path)
config <- read.csv(config_path, check.names = FALSE, stringsAsFactors = FALSE)
n_mp <- as.integer(config$nMP[[1]])

meta_path <- file.path(
  PDO_EXTERNAL_PATHS$sc_ref_pipeline,
  "ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds"
)
tpm_path <- file.path(
  PDO_EXTERNAL_PATHS$sc_ref_pipeline,
  "ref_outs/TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
)
gene_path <- file.path(
  out_dir,
  paste0("Auto_pdo_flot_highres_selected_mp_genes_nMP", n_mp, ".rds")
)
label_path <- file.path(
  out_dir,
  paste0("Auto_pdo_flot_highres_top_3CA_noncellcycle_nMP", n_mp, ".csv")
)
trend_path <- file.path(
  out_dir,
  paste0("Auto_pdo_flot_highres_trend_summary_nMP", n_mp, ".csv")
)
input_paths <- c(meta = meta_path, tpm = tpm_path, genes = gene_path, labels = label_path, trend = trend_path)
pdo_require_files(input_paths, names(input_paths))

selected_genes <- readRDS(gene_path)
label_table <- read.csv(label_path, check.names = FALSE, stringsAsFactors = FALSE)
trend_summary <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
trend_summary$retained <- trend_summary$retained %in% c(TRUE, "TRUE")
retained_order <- trend_summary |>
  filter(retained, MP %in% names(selected_genes)) |>
  mutate(direction_order = match(treatment_direction, c("increase", "decrease"))) |>
  arrange(direction_order, desc(pair_support_n), trend_p_value, MP) |>
  pull(MP)
selected_genes <- selected_genes[retained_order]
selected_genes <- selected_genes[lengths(selected_genes) >= 5L]
if (length(selected_genes) == 0L) stop("No retained MP contains at least five genes.")

label_map <- setNames(label_table$top_3ca_noncc, label_table$MP)
label_map <- label_map[names(selected_genes)]
label_map[is.na(label_map) | !nzchar(label_map)] <- "no non-cell-cycle 3CA match"
display_map <- setNames(
  paste0(names(selected_genes), " | ", unname(label_map)),
  names(selected_genes)
)

####################
# TCGA EAC GSVA scores
####################
infer_histology <- function(type_vector) {
  type_vector <- tolower(as.character(type_vector))
  ifelse(grepl("adeno", type_vector), "EAC", "Other")
}

meta_tcga <- readRDS(meta_path)
if (!"HistologyGroup" %in% colnames(meta_tcga)) {
  if (!"type" %in% colnames(meta_tcga)) stop("TCGA metadata lacks type and HistologyGroup.")
  meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)
}

tpm_df <- data.table::fread(tpm_path)
if (ncol(tpm_df) < 2L || !"GeneSymbol" %in% colnames(tpm_df)) {
  stop("TCGA TPM input lacks GeneSymbol or sample columns.")
}
tpm_matrix <- as.matrix(tpm_df[, -1])
storage.mode(tpm_matrix) <- "numeric"
rownames(tpm_matrix) <- tpm_df$GeneSymbol
tpm_matrix <- tpm_matrix[!duplicated(rownames(tpm_matrix)) & nzchar(rownames(tpm_matrix)), , drop = FALSE]

score_sets <- lapply(selected_genes, function(genes) {
  intersect(unique(genes), rownames(tpm_matrix))
})
score_sets <- score_sets[lengths(score_sets) >= 5L]
if (length(score_sets) == 0L) stop("No selected MP has at least five genes in the TCGA TPM matrix.")

message("Scoring ", length(score_sets), " retained high-resolution MPs in TCGA ESCA.")
gsva_scores <- GSVA::gsva(
  tpm_matrix,
  score_sets,
  method = "gsva",
  kcdf = "Gaussian"
)
saveRDS(
  gsva_scores,
  file.path(survival_dir, paste0("Auto_pdo_flot_highres_tcga_gsva_scores_nMP", n_mp, ".rds"))
)

score_df <- as.data.frame(t(gsva_scores), check.names = FALSE)
score_df$sample_barcode <- rownames(score_df)
survival_data <- meta_tcga |>
  filter(sample_type_code == "01", HistologyGroup == "EAC") |>
  inner_join(score_df, by = "sample_barcode")

if (nrow(survival_data) < 20L) {
  stop("Fewer than 20 primary EAC cases remained after joining TCGA metadata and scores.")
}

####################
# Cox continuous and split sensitivity analyses
####################
run_cox <- function(data, feature, split_method) {
  model_data <- data |>
    filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feature]]))
  if (nrow(model_data) < 20L || var(model_data[[feature]], na.rm = TRUE) == 0) return(NULL)

  if (split_method == "continuous") {
    model_data$split_value <- model_data[[feature]]
  } else if (split_method == "median") {
    cut_value <- median(model_data[[feature]], na.rm = TRUE)
    model_data$split_value <- factor(
      ifelse(model_data[[feature]] > cut_value, "High", "Low"),
      levels = c("Low", "High")
    )
  } else if (split_method == "q1q4") {
    cut_values <- quantile(model_data[[feature]], c(0.25, 0.75), na.rm = TRUE)
    model_data <- model_data |>
      filter(.data[[feature]] <= cut_values[[1]] | .data[[feature]] >= cut_values[[2]])
    if (nrow(model_data) < 20L) return(NULL)
    model_data$split_value <- factor(
      ifelse(model_data[[feature]] >= cut_values[[2]], "High", "Low"),
      levels = c("Low", "High")
    )
  } else {
    stop("Unknown split method: ", split_method)
  }

  if (length(unique(model_data$split_value)) < 2L) return(NULL)
  fit <- try(
    survival::coxph(survival::Surv(OS_time, OS_event) ~ split_value, data = model_data),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  fit_summary <- summary(fit)
  data.frame(
    cohort = "TCGA EAC",
    method = "whole_tcga_gsva",
    feature = feature,
    display_label = unname(display_map[feature]),
    treatment_direction = trend_summary$treatment_direction[match(feature, trend_summary$MP)],
    pair_support_n = trend_summary$pair_support_n[match(feature, trend_summary$MP)],
    split_method = split_method,
    HR = fit_summary$coefficients[1, "exp(coef)"],
    P_value = fit_summary$coefficients[1, "Pr(>|z|)"],
    n = fit$n,
    events = fit$nevent,
    stringsAsFactors = FALSE
  )
}

split_methods <- c("continuous", "median", "q1q4")
cox_results <- bind_rows(lapply(split_methods, function(split_method) {
  bind_rows(lapply(names(score_sets), function(feature) {
    run_cox(survival_data, feature, split_method)
  }))
}))
if (nrow(cox_results) == 0L) stop("No high-resolution MP Cox model completed.")

cox_results <- cox_results |>
  group_by(split_method) |>
  mutate(padj = p.adjust(P_value, method = "BH")) |>
  ungroup()

cox_csv <- file.path(
  survival_dir,
  paste0("Auto_pdo_flot_highres_tcga_survival_cox_splits_nMP", n_mp, ".csv")
)
write.csv(cox_results, cox_csv, row.names = FALSE)

####################
# Auto_08-style survival volcano figure
####################
plot_volcano <- function(plot_data, split_method) {
  plot_data <- plot_data |>
    mutate(
      significant = P_value < 0.05,
      log2_hr = log2(HR),
      neg_log10_p = -log10(P_value),
      display_label = factor(display_label, levels = display_map[names(score_sets)])
    )

  significant_labels <- plot_data |>
    filter(significant)
  borderline_pool <- plot_data |>
    filter(!significant, is.finite(P_value)) |>
    arrange(P_value) |>
    slice_head(n = 15L)
  set.seed(20260813L + match(split_method, split_methods))
  borderline_labels <- if (nrow(borderline_pool) > 5L) {
    borderline_pool[sample(seq_len(nrow(borderline_pool)), 5L), , drop = FALSE]
  } else {
    borderline_pool
  }
  label_data <- bind_rows(significant_labels, borderline_labels) |>
    distinct(feature, .keep_all = TRUE)

  ggplot(plot_data, aes(log2_hr, neg_log10_p)) +
    geom_point(aes(colour = significant, shape = treatment_direction), size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, colour = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey45") +
    geom_text_repel(
      data = label_data,
      aes(label = display_label),
      size = 2.6,
      max.overlaps = Inf
    ) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    scale_shape_manual(values = c(increase = 16, decrease = 17), name = "FLOT trend") +
    theme_minimal(base_size = 11) +
    labs(
      title = split_method,
      x = "log2(HR)",
      y = "-log10(p)"
    )
}

volcano_plots <- lapply(split_methods, function(split_method) {
  plot_volcano(
    filter(cox_results, .data$split_method == .env$split_method),
    split_method
  )
})
volcano_grob <- gridExtra::arrangeGrob(
  grobs = volcano_plots,
  ncol = 3,
  top = grid::textGrob(
    "TCGA EAC survival: FLOT-selected centred high-resolution MPs",
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
)
volcano_pdf <- file.path(
  survival_dir,
  paste0("Auto_pdo_flot_highres_tcga_survival_volcano_nMP", n_mp, ".pdf")
)
grDevices::cairo_pdf(volcano_pdf, width = 18, height = 8, onefile = TRUE)
grid::grid.draw(volcano_grob)
dev.off()

pdo_write_run_summary(
  script = "analysis/cell_states/Auto_pdo_flot_matched_survival_and_state_plots.R",
  out_dir = out_dir,
  inputs = unname(input_paths),
  outputs = c(cox_csv, volcano_pdf),
  parameters = list(
    n_mp = n_mp,
    tested_mp_n = length(score_sets),
    primary_eac_n = nrow(survival_data),
    split_methods = paste(split_methods, collapse = ","),
    manual_grouping = FALSE
  )
)

message("Selected high-resolution MP TCGA survival analysis completed.")
