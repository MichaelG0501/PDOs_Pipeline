####################
# Analysis registry:
#   Status: active terminal clinical figure script
#   Script: analysis/clinical/legacy_clinical_association_final_figures.R
#   Methodology: analysis/methodology/clinical/Auto_centred_clinical_auc_response_associations_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Visualises baseline PDO state proportions and MP activities across clinical
#     variables using the updated V3 clinical metadata across all 20 available
#     untreated PDO samples (21 untreated samples excluding SUR843T3_PDO). Produces
#     stacked barplots and scRef-style boxplots with Wilcoxon / Kruskal-Wallis stats.
#   Inputs:
#     PDOs_outs/PDOs_merged.rds
#     PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds
#     PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     PDOs_outs/centred_mp_refinement/centred_refined_mp_strict_order.rds
#     /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDO_ClinicalMetadata_V3.xlsx
#   Outputs:
#     PDOs_outs/Auto_clinical_assoc_stacked_final/Auto_PDO_clinical_assoc_stacked_final_combined.pdf
#     PDOs_outs/Auto_clinical_assoc_stacked_final/Auto_PDO_clinical_assoc_stacked_final_per_batch.pdf
#     PDOs_outs/Auto_clinical_assoc_stacked_final/Auto_PDO_clinical_assoc_stacked_final_summary.csv
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_mp_boxplots_scref_style.pdf
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_mp_boxplots_scref_style_stats.csv
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_mp_boxplots_scref_style_per_batch.pdf
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_mp_boxplots_scref_style_per_batch_stats.csv
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_state_boxplots_scref_style.pdf
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_state_boxplots_scref_style_stats.csv
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_state_boxplots_scref_style_per_batch.pdf
#     PDOs_outs/Auto_clinical_assoc_boxplots_scref_style/Auto_PDO_clinical_assoc_state_boxplots_scref_style_per_batch_stats.csv
#   Downstream: Terminal
#   Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################

####################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
  library(readxl)
  library(data.table)
})
####################

####################
# Paths, constants, and guards
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

model_dir <- file.path(PDO_OUTPUT_DIR, "centred_mp_refinement")
stacked_out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_clinical_assoc_stacked_final")
boxplot_out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_clinical_assoc_boxplots_scref_style")

dir.create(stacked_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(boxplot_out_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  pdo = file.path(PDO_OUTPUT_DIR, "PDOs_merged.rds"),
  states = file.path(model_dir, "centred_refined_noreg_states.rds"),
  ucell = file.path(model_dir, "merged_refined_ucell_scores.rds"),
  mp_genes = file.path(model_dir, "merged_refined_mp_genes.rds"),
  strict_order = file.path(model_dir, "centred_refined_mp_strict_order.rds"),
  clinical = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDO_ClinicalMetadata_V3.xlsx"
)
pdo_require_files(input_paths)

clean_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "N/A", "NA", "n/a", "Unknown")] <- NA_character_
  x
}

normalise_response <- function(x) {
  x <- clean_missing(x)
  case_when(
    str_detect(x, regex("^non[- ]?responder$|^NR$", ignore_case = TRUE)) ~ "Non-responder",
    str_detect(x, regex("^responder$|^R$", ignore_case = TRUE)) ~ "Responder",
    TRUE ~ x
  )
}

message("=== Loading Data ===")
tmdata_all <- readRDS(input_paths[["pdo"]])
state_vec <- pdo_normalise_final_state_vector(readRDS(input_paths[["states"]]))
ucell_scores <- readRDS(input_paths[["ucell"]])
strict_object <- readRDS(input_paths[["strict_order"]])
mp_order <- strict_object$strict_refined_mp_order
if (is.null(mp_order)) stop("Strict-order object lacks strict_refined_mp_order.")

missing_mps <- setdiff(mp_order, colnames(ucell_scores))
if (length(missing_mps) > 0L) stop("UCell matrix misses current MPs: ", paste(missing_mps, collapse = "; "))

####################
# Load and standardize V3 clinical metadata
####################
clinical <- readxl::read_excel(input_paths[["clinical"]], sheet = "Sheet1") %>%
  transmute(
    sur = paste0("SUR", as.integer(SUR)),
    gender = clean_missing(Gender),
    age = suppressWarnings(as.numeric(Age)),
    age_group = case_when(is.na(age) ~ NA_character_, age > 60 ~ ">60", TRUE ~ "<=60"),
    tumour_location = clean_missing(`Tumour location`),
    tumour_type = clean_missing(`Tumour type`),
    histology = clean_missing(Histology),
    ajcc = clean_missing(AJCC),
    clinical_response = normalise_response(`Clinical response at OG MDT: responder/non-responder`),
    mandard_score = clean_missing(`Mandard tumour regression score`),
    mandard_response = normalise_response(`Response based on Mandard`),
    pdo_origin = case_when(
      clean_missing(`Origin of PDO (primary tumour vs lymph node)`) == "T" ~ "Primary tumour",
      clean_missing(`Origin of PDO (primary tumour vs lymph node)`) == "LN" ~ "Lymph node",
      TRUE ~ clean_missing(`Origin of PDO (primary tumour vs lymph node)`)
    ),
    ctx_timepoint = clean_missing(`Pre CTX or Post CTX`),
    auc_rel = suppressWarnings(as.numeric(clean_missing(`AUCrel (FLOT)`))),
    t_stage = clean_missing(str_extract(`cTNM stage at diagnosis`, "T[0-4][a-b]?")),
    histology_at_surgery = case_when(
      str_detect(`Tumour histology at surgery`, regex("adenocarcinoma", ignore_case = TRUE)) ~ "Adenocarcinoma",
      str_detect(`Tumour histology at surgery`, regex("squamous", ignore_case = TRUE)) ~ "Adenosquamous",
      TRUE ~ clean_missing(`Tumour histology at surgery`)
    )
  )

meta_df <- tmdata_all@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    orig.ident = as.character(orig.ident),
    sample = as.character(orig.ident),
    batch = ifelse(is.na(Batch) | !nzchar(as.character(Batch)), "Unknown", as.character(Batch)),
    Batch = batch,
    sur = str_extract(sample, "^SUR[0-9]+"),
    state = as.character(state_vec[cell])
  ) %>%
  filter(
    sample != PDO_EXCLUDED_SAMPLE,
    !str_detect(sample, "_Treated_PDO$")
  ) %>%
  left_join(clinical, by = "sur")

sample_count <- n_distinct(meta_df$sample)
if (sample_count != 20) {
  stop("Expected 20 untreated PDO samples (21 untreated excluding SUR843T3_PDO); observed ", sample_count)
}
message("PDO samples retained after untreated filtering: ", sample_count)
message("Retained samples: ", paste(sort(unique(meta_df$sample)), collapse = ", "))

plot_configs <- list(
  list(var = "gender", title = "Gender", filter = NULL),
  list(var = "age_group", title = "Age (>60)", filter = NULL),
  list(var = "batch", title = "Batch", filter = NULL),
  list(var = "clinical_response", title = "Clinical Response", filter = NULL),
  list(var = "tumour_location", title = "Tumour Location", filter = NULL),
  list(var = "tumour_type", title = "Tumour Type", filter = NULL),
  list(var = "histology", title = "Histology", filter = NULL),
  list(var = "ajcc", title = "AJCC Stage", filter = NULL),
  list(var = "mandard_score", title = "Mandard Score", filter = NULL),
  list(var = "mandard_response", title = "Mandard Response", filter = NULL),
  list(var = "pdo_origin", title = "Origin of PDO", filter = NULL),
  list(var = "t_stage", title = "cTNM T-Stage", filter = NULL),
  list(var = "ctx_timepoint", title = "Collection Timepoint", filter = NULL)
)

apply_optional_filter <- function(data, filter_expr = NULL) {
  if (is.null(filter_expr) || !nzchar(filter_expr)) {
    return(data)
  }
  data %>% filter(!!rlang::parse_expr(filter_expr))
}

####################
# STACKED BAR PLOTS
####################
state_levels <- PDO_STATE_ORDER_WITH_OPTIONAL
present_states <- unique(as.character(meta_df$state))
state_levels_present <- c(
  PDO_STATE_ORDER[PDO_STATE_ORDER %in% present_states],
  c("Unresolved", "Hybrid")[c("Unresolved", "Hybrid") %in% present_states],
  setdiff(present_states, PDO_STATE_ORDER_WITH_OPTIONAL)
)

state_colors <- PDO_STATE_COLORS
missing_color_states <- setdiff(state_levels_present, names(state_colors))
if (length(missing_color_states) > 0) {
  state_colors <- c(state_colors, setNames(scales::hue_pal()(length(missing_color_states)), missing_color_states))
}

compute_plot_data <- function(data, group_var, filter_expr = NULL, facet_var = NULL) {
  data <- apply_optional_filter(data, filter_expr) %>%
    filter(!is.na(.data[[group_var]]), !is.na(sample), !is.na(state), !is.na(batch))
  if (nrow(data) == 0) return(list(plot_data = NULL, label_data = NULL))

  sample_group_cols <- unique(c(facet_var, "sample", group_var, "state"))
  sample_total_cols <- unique(c(facet_var, "sample", group_var))
  mean_group_cols <- unique(c(facet_var, group_var, "state"))
  total_group_cols <- unique(c(facet_var, group_var))

  sample_prop <- data %>%
    group_by(across(all_of(sample_group_cols))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(sample_total_cols))) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()

  plot_data <- sample_prop %>%
    group_by(across(all_of(mean_group_cols))) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(total_group_cols))) %>%
    mutate(group_total = sum(mean_prop), pct = ifelse(group_total > 0, 100 * mean_prop / group_total, 0)) %>%
    ungroup() %>%
    mutate(state = factor(state, levels = state_levels_present))

  label_data <- data %>%
    distinct(across(all_of(unique(c(facet_var, "sample", group_var, "batch"))))) %>%
    group_by(across(all_of(total_group_cols))) %>%
    summarise(
      n_samples = n_distinct(sample),
      n_batches = n_distinct(batch),
      batch_list = paste(sort(unique(batch)), collapse = "\n"),
      .groups = "drop"
    ) %>%
    left_join(
      data %>%
        group_by(across(all_of(total_group_cols))) %>%
        summarise(n_cells = n(), .groups = "drop"),
      by = total_group_cols
    )

  if (is.null(facet_var)) {
    total_samples <- n_distinct(data$sample)
    total_batches <- n_distinct(data$batch)
    label_data <- label_data %>%
      mutate(label = paste0("N=", comma(n_cells), "\n", n_samples, "/", total_samples, "\n(", n_batches, "/", total_batches, ")\n\n", batch_list))
  } else {
    facet_stats <- data %>%
      group_by(across(all_of(facet_var))) %>%
      summarise(facet_cells = n(), facet_samples = n_distinct(sample), .groups = "drop") %>%
      mutate(facet_label = paste0(.data[[facet_var]], "\n(N=", comma(facet_cells), "; samples=", facet_samples, ")"))
    plot_data <- plot_data %>% left_join(facet_stats[, c(facet_var, "facet_label")], by = facet_var)
    label_data <- label_data %>%
      left_join(facet_stats, by = facet_var) %>%
      mutate(label = paste0(n_samples, "/", facet_samples), facet_label = paste0(.data[[facet_var]], "\n(N=", comma(facet_cells), "; samples=", facet_samples, ")"))
  }

  list(plot_data = plot_data, label_data = label_data)
}

plot_stacked_assoc <- function(data, group_var, title, filter_expr = NULL, facet_var = NULL) {
  prep <- compute_plot_data(data, group_var, filter_expr, facet_var)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)

  x_levels <- unique(as.character(prep$plot_data[[group_var]]))
  detail_lines <- if (is.null(facet_var)) max(stringr::str_count(prep$label_data$label, "\\n") + 1) else 2
  expansion_mult <- if (is.null(facet_var)) 0.05 + 0.035 * detail_lines else 0.16

  p <- ggplot(prep$plot_data, aes(x = .data[[group_var]], y = pct / 100, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.72) +
    geom_text(
      data = prep$plot_data %>% mutate(label_text = ifelse(pct >= 8, sprintf("%.0f%%", pct), "")),
      aes(label = label_text, group = state),
      position = position_stack(vjust = 0.5),
      size = 3.8,
      fontface = "bold",
      colour = "black",
      show.legend = FALSE
    ) +
    geom_text(
      data = prep$label_data,
      aes(x = .data[[group_var]], y = 1.02, label = label),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4,
      lineheight = 1.0,
      fontface = "bold"
    ) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, expansion_mult))) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"), width = 16)) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = ifelse(length(x_levels) > 3, 45, 25), hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 15),
      plot.margin = margin(10, 20, 10, 10)
    )

  if (!is.null(facet_var)) {
    p <- p +
      facet_wrap(~ facet_label, scales = "free_x", ncol = 3) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 10)
      )
  }

  p
}

safe_plot <- function(cfg, facet_var = NULL) {
  tryCatch(
    plot_stacked_assoc(meta_df, cfg$var, cfg$title, cfg$filter, facet_var),
    error = function(e) {
      message("Skipping ", cfg$title, ": ", e$message)
      NULL
    }
  )
}

pdf(file.path(stacked_out_dir, "Auto_PDO_clinical_assoc_stacked_final_combined.pdf"), width = 17, height = 12, useDingbats = FALSE)
for (cfg in plot_configs) {
  p <- safe_plot(cfg)
  if (!is.null(p)) print(p)
}
dev.off()

pdf(file.path(stacked_out_dir, "Auto_PDO_clinical_assoc_stacked_final_per_batch.pdf"), width = 17, height = 11, useDingbats = FALSE)
for (cfg in plot_configs) {
  p <- safe_plot(cfg, facet_var = "batch")
  if (!is.null(p)) print(p)
}
dev.off()

summary_df <- bind_rows(lapply(plot_configs, function(cfg) {
  prep <- compute_plot_data(meta_df, cfg$var, cfg$filter)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
  prep$plot_data %>%
    left_join(prep$label_data, by = cfg$var) %>%
    transmute(
      clinical_variable = cfg$var,
      plot_title = cfg$title,
      filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter),
      level = .data[[cfg$var]],
      state = as.character(state),
      mean_sample_prop_pct = pct,
      samples_in_level = n_samples,
      cells_in_level = n_cells
    )
}))

write.csv(
  summary_df,
  file.path(stacked_out_dir, "Auto_PDO_clinical_assoc_stacked_final_summary.csv"),
  row.names = FALSE
)
message("Saved PDO final stacked clinical association outputs.")

####################
# BOX PLOTS (scRef style)
####################
clinical_level_orders <- list(
  gender = c("Female", "Male"),
  age_group = c("<=60", ">60"),
  batch = c("PDO", "Untreated_PDO", "new4samples"),
  clinical_response = c("Responder", "Non-responder"),
  mandard_response = c("Responder", "Non-responder"),
  mandard_score = c("2", "3", "4", "5"),
  ajcc = c("2", "3", "4"),
  histology = c("MD", "PD"),
  pdo_origin = c("Primary tumour", "Lymph node"),
  ctx_timepoint = c("Pre", "Post"),
  t_stage = c("T1", "T2", "T3", "T4")
)

mp_descriptions <- PDO_MP_DESCRIPTIONS
state_groups <- PDO_MP_STATE_GROUPS

state_axis_labels <- c(
  "Classic proliferation" = "Classic\nproliferation",
  "Columnar-to-intestinal" = "Columnar-to-\nintestinal",
  "Glandular differentiation" = "Glandular\ndifferentiation",
  "Stress-adaptive" = "Stress-\nadaptive",
  "ECM-remodelling" = "ECM-\nremodelling",
  "Motile-cilia differentiation" = "Motile-cilia\ndifferentiation",
  "Unresolved" = "Unresolved",
  "Hybrid" = "Hybrid"
)

make_group_order <- function(values, group_var) {
  present_levels <- unique(as.character(values))
  present_levels <- present_levels[!is.na(present_levels) & nzchar(present_levels)]
  preferred_levels <- clinical_level_orders[[group_var]]
  ordered_levels <- character(0)
  if (!is.null(preferred_levels)) {
    ordered_levels <- preferred_levels[preferred_levels %in% present_levels]
  }
  remaining_levels <- setdiff(present_levels, ordered_levels)
  c(ordered_levels, sort(remaining_levels))
}

make_group_palette <- function(levels_use) {
  base_cols <- c(
    "Female" = "#CC79A7",
    "Male" = "#0072B2",
    "<=60" = "#00A087",
    ">60" = "#E69F00",
    "Untreated_PDO" = "#F39B7F",
    "PDO" = "#8491B4",
    "new4samples" = "#2CA02C",
    "Responder" = "#009E73",
    "Non-responder" = "#D55E00",
    "Pre" = "#56B4E9",
    "Post" = "#E69F00",
    "Primary tumour" = "#4DAF4A",
    "Lymph node" = "#984EA3",
    "MD" = "#4DBBD5",
    "PD" = "#E64B35",
    "2" = "#8DD3C7",
    "3" = "#BEBADA",
    "4" = "#FB8072",
    "5" = "#80B1D3",
    "T1" = "#A6CEE3",
    "T2" = "#1F78B4",
    "T3" = "#B2DF8A",
    "T4" = "#33A02C"
  )
  palette <- base_cols[levels_use]
  missing_levels <- setdiff(levels_use, names(base_cols))
  if (length(missing_levels) > 0) {
    palette[missing_levels] <- scales::hue_pal(l = 65, c = 100)(length(missing_levels))
  }
  palette[levels_use]
}

significance_label <- function(p_val) {
  dplyr::case_when(
    is.na(p_val) ~ "",
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

compute_feature_stats <- function(sample_df, clinical_var, feature_type, feature_label_map, digits = 3) {
  feature_splits <- split(sample_df, as.character(sample_df$feature))
  stats_rows <- lapply(feature_splits, function(df_feature) {
    df_feature <- df_feature %>%
      filter(!is.na(group), !is.na(value))

    present_groups <- unique(as.character(df_feature$group))
    feature_name <- as.character(df_feature$feature[1])
    if (length(present_groups) < 2 || n_distinct(df_feature$sample) < 3) {
      return(data.frame(
        clinical_variable = clinical_var,
        feature_type = feature_type,
        feature = feature_name,
        feature_label = unname(feature_label_map[feature_name]),
        test = NA_character_,
        n_groups = length(present_groups),
        n_samples = n_distinct(df_feature$sample),
        p_value = NA_real_,
        group_summary = NA_character_,
        median_range = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    test_name <- if (length(present_groups) == 2) "wilcoxon" else "kruskal"
    test_res <- tryCatch(
      suppressWarnings(
        if (test_name == "wilcoxon") {
          wilcox.test(value ~ group, data = df_feature, exact = FALSE)
        } else {
          kruskal.test(value ~ group, data = df_feature)
        }
      ),
      error = function(e) NULL
    )

    median_df <- df_feature %>%
      group_by(group) %>%
      summarise(
        n_samples = n_distinct(sample),
        median_value = median(value, na.rm = TRUE),
        mean_value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

    data.frame(
      clinical_variable = clinical_var,
      feature_type = feature_type,
      feature = feature_name,
      feature_label = unname(feature_label_map[feature_name]),
      test = test_name,
      n_groups = length(present_groups),
      n_samples = n_distinct(df_feature$sample),
      p_value = if (is.null(test_res)) NA_real_ else test_res$p.value,
      group_summary = paste0(
        as.character(median_df$group),
        " (n=",
        median_df$n_samples,
        ", median=",
        formatC(median_df$median_value, format = "f", digits = digits),
        ")",
        collapse = " | "
      ),
      median_range = diff(range(median_df$median_value, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })

  stats_df <- bind_rows(stats_rows)
  if (nrow(stats_df) == 0) return(stats_df)

  stats_df %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significance = significance_label(p_value)
    )
}

plot_feature_boxplot <- function(sample_df, stats_df, title_text, subtitle_text, y_label, feature_labels, feature_type) {
  sample_df <- sample_df %>%
    filter(!is.na(group), !is.na(feature), !is.na(value))

  if (nrow(sample_df) == 0 || n_distinct(sample_df$group) < 2) return(NULL)

  group_levels <- levels(sample_df$group)
  if (is.null(group_levels)) {
    group_levels <- unique(as.character(sample_df$group))
  }

  legend_counts <- sample_df %>%
    distinct(sample, group) %>%
    count(group, name = "n_samples")

  legend_labels <- setNames(
    paste0(
      iconv(as.character(legend_counts$group), from = "UTF-8", to = "ASCII//TRANSLIT"),
      " (n=",
      legend_counts$n_samples,
      ")"
    ),
    as.character(legend_counts$group)
  )

  palette <- make_group_palette(group_levels)
  y_range <- range(sample_df$value, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- if (feature_type == "state") 5 else 0.02
  }
  y_offset <- if (feature_type == "state") max(2, 0.08 * y_span) else max(0.01, 0.08 * y_span)

  annot_df <- sample_df %>%
    group_by(feature) %>%
    summarise(y_pos = max(value, na.rm = TRUE) + y_offset, .groups = "drop") %>%
    left_join(stats_df %>% select(feature, significance, p_value), by = "feature") %>%
    mutate(label = ifelse(!is.na(p_value) & p_value < 0.05, significance, ""))

  p <- ggplot(sample_df, aes(x = feature, y = value, fill = group, color = group)) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.8,
      linewidth = 0.4,
      color = "black"
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      alpha = 0.7,
      size = 1.0,
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = annot_df %>% filter(label != ""),
      aes(x = feature, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = palette, labels = legend_labels, drop = FALSE) +
    scale_color_manual(values = palette, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = feature_labels) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = y_label,
      fill = "Clinical group"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "grey35"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line.x = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(12, 16, 12, 12)
    )

  if (feature_type == "state") {
    p <- p +
      scale_y_continuous(
        labels = function(x) paste0(formatC(x, format = "f", digits = 0), "%"),
        expand = expansion(mult = c(0.02, 0.18))
      )
  } else {
    p <- p +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)))
  }

  p
}

####################
# MP labels and grouping
####################
mp_feature_labels <- setNames(
  ifelse(!is.na(mp_descriptions[mp_order]), paste0(mp_order, " ", mp_descriptions[mp_order]), mp_order),
  mp_order
)

mp_axis_labels <- setNames(
  ifelse(!is.na(mp_descriptions[mp_order]), paste0(mp_order, "\n", mp_descriptions[mp_order]), mp_order),
  mp_order
)

mp_to_group <- setNames(rep("Other", length(mp_order)), mp_order)
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], mp_order)
  mp_to_group[grp_mps] <- grp
}

cell_long <- as.data.frame(ucell_scores[meta_df$cell, mp_order, drop = FALSE]) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(mp_order), names_to = "MP", values_to = "score") %>%
  left_join(meta_df, by = "cell")

state_cell_df <- meta_df

build_mp_sample_df <- function(data, group_var, filter_expr = NULL) {
  data_use <- apply_optional_filter(data, filter_expr) %>%
    filter(!is.na(.data[[group_var]]), !is.na(sample), !is.na(MP), !is.na(score), !is.na(batch))

  if (nrow(data_use) == 0) return(NULL)

  group_levels <- make_group_order(data_use[[group_var]], group_var)
  if (length(group_levels) < 2) return(NULL)

  data_use %>%
    mutate(group = factor(as.character(.data[[group_var]]), levels = group_levels)) %>%
    group_by(sample, batch, group, MP) %>%
    summarise(value = mean(score, na.rm = TRUE), .groups = "drop") %>%
    transmute(
      orig.ident = sample,
      sample = sample,
      Batch = batch,
      batch = batch,
      group = group,
      feature = factor(MP, levels = mp_order),
      value = value
    )
}

build_state_sample_df <- function(data, group_var, filter_expr = NULL) {
  data_use <- apply_optional_filter(data, filter_expr) %>%
    filter(!is.na(.data[[group_var]]), !is.na(sample), !is.na(state), !is.na(batch), !(state %in% c("Unresolved", "Hybrid")))

  if (nrow(data_use) == 0) return(NULL)

  group_levels <- make_group_order(data_use[[group_var]], group_var)
  if (length(group_levels) < 2) return(NULL)

  data_use <- data_use %>%
    mutate(group = factor(as.character(.data[[group_var]]), levels = group_levels))

  state_levels_use <- PDO_STATE_ORDER

  sample_meta <- data_use %>%
    transmute(
      sample = sample,
      batch = batch,
      group = group
    ) %>%
    distinct()

  sample_totals <- data_use %>%
    group_by(sample, group) %>%
    summarise(total_cells = n(), .groups = "drop")

  state_counts <- data_use %>%
    group_by(sample, group, state) %>%
    summarise(n_cells = n(), .groups = "drop")

  sample_meta %>%
    tidyr::crossing(state = state_levels_use) %>%
    left_join(state_counts, by = c("sample", "group", "state")) %>%
    mutate(n_cells = replace_na(n_cells, 0L)) %>%
    left_join(sample_totals, by = c("sample", "group")) %>%
    mutate(
      group = factor(group, levels = group_levels),
      feature = factor(state, levels = state_levels_use),
      value = 100 * n_cells / pmax(total_cells, 1)
    ) %>%
    transmute(
      orig.ident = sample,
      sample = sample,
      Batch = batch,
      batch = batch,
      group = group,
      feature = feature,
      value = value
    )
}

####################
# MP clinical boxplots
####################
mp_stats_all <- list()
pdf(file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_mp_boxplots_scref_style.pdf"), width = 18, height = 9, useDingbats = FALSE)
for (cfg in plot_configs) {
  mp_sample_df <- build_mp_sample_df(cell_long, cfg$var, cfg$filter)
  if (is.null(mp_sample_df) || nrow(mp_sample_df) == 0 || n_distinct(mp_sample_df$group) < 2) next

  stats_df <- compute_feature_stats(
    mp_sample_df,
    clinical_var = cfg$var,
    feature_type = "MP",
    feature_label_map = mp_feature_labels,
    digits = 3
  ) %>%
    mutate(
      plot_title = cfg$title,
      filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter),
      feature_group = mp_to_group[as.character(feature)]
    )

  mp_stats_all[[cfg$title]] <- stats_df

  p <- plot_feature_boxplot(
    mp_sample_df,
    stats_df,
    title_text = paste0(cfg$title, " - PDO MP activity"),
    subtitle_text = "Sample-level mean UCell scores; stars mark BH-adjusted p < 0.05 across clinical groups.",
    y_label = "Mean sample UCell score",
    feature_labels = mp_axis_labels,
    feature_type = "mp"
  )
  if (!is.null(p)) print(p)
}
dev.off()

write.csv(
  bind_rows(mp_stats_all),
  file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_mp_boxplots_scref_style_stats.csv"),
  row.names = FALSE
)

####################
# Final-state clinical boxplots
####################
state_stats_all <- list()
pdf(file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_state_boxplots_scref_style.pdf"), width = 16, height = 9, useDingbats = FALSE)
for (cfg in plot_configs) {
  state_sample_df <- build_state_sample_df(state_cell_df, cfg$var, cfg$filter)
  if (is.null(state_sample_df) || nrow(state_sample_df) == 0 || n_distinct(state_sample_df$group) < 2) next

  stats_df <- compute_feature_stats(
    state_sample_df,
    clinical_var = cfg$var,
    feature_type = "state",
    feature_label_map = setNames(PDO_STATE_ORDER, PDO_STATE_ORDER),
    digits = 1
  ) %>%
    mutate(
      plot_title = cfg$title,
      filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter)
    )

  state_stats_all[[cfg$title]] <- stats_df

  p <- plot_feature_boxplot(
    state_sample_df,
    stats_df,
    title_text = paste0(cfg$title, " - PDO final state proportions"),
    subtitle_text = "Sample-level state proportions; stars mark BH-adjusted p < 0.05 across clinical groups.",
    y_label = "Sample state proportion",
    feature_labels = state_axis_labels,
    feature_type = "state"
  )
  if (!is.null(p)) print(p)
}
dev.off()

write.csv(
  bind_rows(state_stats_all),
  file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_state_boxplots_scref_style_stats.csv"),
  row.names = FALSE
)

####################
# Per-Batch MP clinical boxplots
####################
all_batches <- sort(unique(cell_long$batch))
all_batches <- all_batches[!is.na(all_batches) & nzchar(all_batches)]

mp_stats_per_batch_all <- list()
pdf(file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_mp_boxplots_scref_style_per_batch.pdf"), width = 18, height = 9, useDingbats = FALSE)
for (cfg in plot_configs) {
  for (batch_name in all_batches) {
    batch_data <- cell_long %>% filter(batch == batch_name)

    mp_sample_df <- build_mp_sample_df(batch_data, cfg$var, cfg$filter)
    if (is.null(mp_sample_df) || nrow(mp_sample_df) == 0 || n_distinct(mp_sample_df$group) < 2) next

    stats_df <- compute_feature_stats(
      mp_sample_df,
      clinical_var = cfg$var,
      feature_type = "MP",
      feature_label_map = mp_feature_labels,
      digits = 3
    ) %>%
      mutate(
        batch = batch_name,
        Batch = batch_name,
        plot_title = cfg$title,
        filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter),
        feature_group = mp_to_group[as.character(feature)]
      )

    mp_stats_per_batch_all[[paste0(cfg$title, "_", batch_name)]] <- stats_df

    p <- plot_feature_boxplot(
      mp_sample_df,
      stats_df,
      title_text = paste0("[", batch_name, "] ", cfg$title, " - PDO MP activity"),
      subtitle_text = "Sample-level mean UCell scores; stars mark BH-adjusted p < 0.05 across clinical groups.",
      y_label = "Mean sample UCell score",
      feature_labels = mp_axis_labels,
      feature_type = "mp"
    )
    if (!is.null(p)) print(p)
  }
}
dev.off()

write.csv(
  bind_rows(mp_stats_per_batch_all),
  file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_mp_boxplots_scref_style_per_batch_stats.csv"),
  row.names = FALSE
)

####################
# Per-Batch final-state clinical boxplots
####################
state_stats_per_batch_all <- list()
pdf(file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_state_boxplots_scref_style_per_batch.pdf"), width = 16, height = 9, useDingbats = FALSE)
for (cfg in plot_configs) {
  for (batch_name in all_batches) {
    batch_state_data <- state_cell_df %>% filter(batch == batch_name)

    state_sample_df <- build_state_sample_df(batch_state_data, cfg$var, cfg$filter)
    if (is.null(state_sample_df) || nrow(state_sample_df) == 0 || n_distinct(state_sample_df$group) < 2) next

    stats_df <- compute_feature_stats(
      state_sample_df,
      clinical_var = cfg$var,
      feature_type = "state",
      feature_label_map = setNames(PDO_STATE_ORDER, PDO_STATE_ORDER),
      digits = 1
    ) %>%
      mutate(
        batch = batch_name,
        Batch = batch_name,
        plot_title = cfg$title,
        filter_expr = ifelse(is.null(cfg$filter), "", cfg$filter)
      )

    state_stats_per_batch_all[[paste0(cfg$title, "_", batch_name)]] <- stats_df

    p <- plot_feature_boxplot(
      state_sample_df,
      stats_df,
      title_text = paste0("[", batch_name, "] ", cfg$title, " - PDO final state proportions"),
      subtitle_text = "Sample-level state proportions; stars mark BH-adjusted p < 0.05 across clinical groups.",
      y_label = "Sample state proportion",
      feature_labels = state_axis_labels,
      feature_type = "state"
    )
    if (!is.null(p)) print(p)
  }
}
dev.off()

write.csv(
  bind_rows(state_stats_per_batch_all),
  file.path(boxplot_out_dir, "Auto_PDO_clinical_assoc_state_boxplots_scref_style_per_batch_stats.csv"),
  row.names = FALSE
)

message("All PDO scRef-style clinical boxplot and stacked association outputs complete.")
