####################
# Adapted from scRef_Pipeline for PDOs
# Auto_clinical_variable_plots.R
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(readxl)

####################
# 1) Setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# 2) Load metadata + placeholder state definitions
####################
tmdata_all <- readRDS("PDOs_merged.rds")
meta_df <- tmdata_all@meta.data
meta_df$cell <- rownames(meta_df)

####################
# UPDATE AFTER geneNMF: Fill in MP descriptions and state groups
####################
retained_mps <- character(0)
mp_descriptions <- c()  # Fill after enrichment
state_groups <- list()  # Fill after state definition
# Auto-populate if empty
if (length(mp_descriptions) == 0) mp_descriptions <- setNames(retained_mps, retained_mps)

####################
# UPDATE AFTER geneNMF: Cell-cycle MP list
####################
cc_mps <- c()  # UPDATE: identify cell cycle MPs

####################
# 3) Load PDO clinical xlsx
####################
clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx"
)
clinical_df <- as.data.frame(clinical_sheet)
if (!"SUR" %in% colnames(clinical_df)) {
  rn <- rownames(clinical_df)
  if (!is.null(rn) && all(nchar(rn) > 0)) {
    clinical_df$SUR <- ifelse(grepl("^SUR", rn), rn, paste0("SUR", rn))
  }
}
clinical_df$SUR <- as.character(clinical_df$SUR)
clinical_df <- clinical_df %>% rename_with(~ paste0("clinical_", .x), -SUR)

####################
# 4) Build cell-level table
####################
cell_df <- meta_df %>%
  mutate(
    orig.ident = as.character(orig.ident),
    SUR = ifelse("SUR" %in% colnames(.), as.character(SUR), sub("^(SUR[0-9]+).*", "\\1", as.character(orig.ident))),
    Batch = as.character(Batch)
  ) %>%
  left_join(clinical_df, by = "SUR") %>%
  filter(orig.ident != "SUR843T3_PDO")

age_vec <- if ("Age" %in% colnames(cell_df)) {
  suppressWarnings(as.numeric(cell_df$Age))
} else if ("clinical_Age" %in% colnames(cell_df)) {
  suppressWarnings(as.numeric(cell_df$clinical_Age))
} else {
  rep(NA_real_, nrow(cell_df))
}

gender_vec <- if ("Gender" %in% colnames(cell_df)) {
  as.character(cell_df$Gender)
} else if ("clinical_Gender" %in% colnames(cell_df)) {
  as.character(cell_df$clinical_Gender)
} else {
  rep(NA_character_, nrow(cell_df))
}

response_vec <- if ("Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(cell_df)) {
  as.character(cell_df$Clinical.response.at.OG.MDT..responder.non.responder)
} else if ("clinical_Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(cell_df)) {
  as.character(cell_df$clinical_Clinical.response.at.OG.MDT..responder.non.responder)
} else {
  rep(NA_character_, nrow(cell_df))
}

cell_df <- cell_df %>%
  mutate(
    state = ifelse(!is.na(Batch), Batch, "Unassigned"),
    Age = age_vec,
    Age_Group = ifelse(!is.na(Age) & Age > 60, ">60", "<=60"),
    Gender = ifelse(gender_vec %in% c("F", "female", "Female"), "Female", ifelse(!is.na(gender_vec), "Male", NA_character_)),
    Clinical.response.at.OG.MDT..responder.non.responder = response_vec,
    Clinical.response.at.OG.MDT..responder.non.responder = ifelse(
      Clinical.response.at.OG.MDT..responder.non.responder %in% c("R", "Responder", "responder"),
      "Responder",
      ifelse(
        Clinical.response.at.OG.MDT..responder.non.responder %in% c("NR", "Nonresponder", "non-responder", "nonresponder"),
        "Nonresponder",
        Clinical.response.at.OG.MDT..responder.non.responder
      )
    )
  )

####################
# 5) State order + colours
####################
state_levels <- sort(unique(cell_df$state))
state_levels <- state_levels[!is.na(state_levels)]
state_colors <- setNames(scales::hue_pal()(length(state_levels)), state_levels)

####################
# 6) Per-sample proportions -> group mean
####################
compute_plot_data <- function(data, group_var, filter_expr = NULL, batchname = NULL) {
  if (!is.null(batchname)) {
    data <- data %>% filter(Batch == batchname)
  }
  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  }

  data <- data %>% filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(state), !is.na(Batch))

  if (nrow(data) == 0) {
    return(list(plot_data = NULL, label_data = NULL, total_samples = 0L))
  }

  sample_prop <- data %>%
    group_by(orig.ident, .data[[group_var]], state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident, .data[[group_var]]) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()

  mean_prop <- sample_prop %>%
    group_by(.data[[group_var]], state) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop")

  group_sum <- mean_prop %>%
    group_by(.data[[group_var]]) %>%
    summarise(group_total = sum(mean_prop), .groups = "drop")

  plot_data <- mean_prop %>%
    left_join(group_sum, by = group_var) %>%
    mutate(pct = ifelse(group_total > 0, 100 * mean_prop / group_total, 0), state = factor(state, levels = state_levels))

  total_samples <- data %>% distinct(orig.ident) %>% nrow()
  total_batches <- data %>% distinct(Batch) %>% nrow()

  n_sample_group <- data %>%
    distinct(orig.ident, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_samples = n_distinct(orig.ident), .groups = "drop")

  n_cell_group <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_cells = n(), .groups = "drop")

  batch_group <- data %>%
    distinct(orig.ident, Batch, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n_batches = n_distinct(Batch),
      batch_list = paste(sort(unique(Batch)), collapse = "\n"),
      .groups = "drop"
    )

  label_data <- n_sample_group %>%
    left_join(n_cell_group, by = group_var) %>%
    left_join(batch_group, by = group_var) %>%
    mutate(
      label = if (!is.null(batchname)) {
        paste0("N=", comma(n_cells), "\n", n_samples, "/", total_samples)
      } else {
        paste0("N=", comma(n_cells), "\n", n_samples, "/", total_samples, "\n(", n_batches, "/", total_batches, ")", "\n\n", batch_list)
      }
    )

  list(plot_data = plot_data, label_data = label_data, total_samples = total_samples)
}

####################
# 7) Combined-style plot helper
####################
plot_clinical_assoc <- function(data, group_var, title, filter_expr = NULL, batchname = NULL) {
  prep <- compute_plot_data(data, group_var, filter_expr, batchname)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)

  show_details <- is.null(batchname)
  max_lines <- max(stringr::str_count(prep$label_data$label, "\\n") + 1)
  expansion_mult <- if (show_details) {
    0.04 + (max_lines * 0.04)
  } else {
    0.15
  }

  ggplot(prep$plot_data, aes(x = .data[[group_var]], y = pct / 100, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    geom_text(
      data = prep$label_data,
      aes(x = .data[[group_var]], y = 1.02, label = label),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3,
      lineheight = 0.9,
      color = "black"
    ) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, expansion_mult))) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")) +
    theme_minimal(base_size = 16) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 18)
    )
}

####################
# 8) Per-Batch faceted helper
####################
plot_variable_facet <- function(data, group_var, title, filter_expr = NULL) {
  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  }

  data <- data %>% filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(Batch), !is.na(state))
  if (nrow(data) == 0) return(NULL)

  sample_prop <- data %>%
    group_by(Batch, orig.ident, .data[[group_var]], state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Batch, orig.ident, .data[[group_var]]) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()

  mean_prop <- sample_prop %>%
    group_by(Batch, .data[[group_var]], state) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop")

  group_sum <- mean_prop %>%
    group_by(Batch, .data[[group_var]]) %>%
    summarise(group_total = sum(mean_prop), .groups = "drop")

  plot_data <- mean_prop %>%
    left_join(group_sum, by = c("Batch", group_var)) %>%
    mutate(freq = ifelse(group_total > 0, mean_prop / group_total, 0), state = factor(state, levels = state_levels))

  batch_stats <- data %>%
    group_by(Batch) %>%
    summarise(n_cells = n(), n_samples = n_distinct(orig.ident), .groups = "drop") %>%
    mutate(batch_label = paste0(Batch, "\n(N=", comma(n_cells), "; samples=", n_samples, ")"))

  label_data <- data %>%
    distinct(Batch, orig.ident, .data[[group_var]]) %>%
    group_by(Batch, .data[[group_var]]) %>%
    summarise(n_samples_level = n_distinct(orig.ident), .groups = "drop") %>%
    left_join(batch_stats, by = "Batch") %>%
    mutate(label = paste0(n_samples_level, "/", n_samples), batch_label = paste0(Batch, "\n(N=", comma(n_cells), "; samples=", n_samples, ")"))

  plot_data <- plot_data %>% left_join(batch_stats[, c("Batch", "batch_label")], by = "Batch")

  ggplot(plot_data, aes(x = .data[[group_var]], y = freq, fill = state)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    geom_text(data = label_data, aes(x = .data[[group_var]], y = 1.02, label = label), inherit.aes = FALSE, size = 3) +
    facet_wrap(~ batch_label, scales = "free_x", ncol = 4) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = state_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")) +
    theme_bw(base_size = 14) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 11)
    )
}

####################
# 9) Plot configuration (PDO clinical variables)
####################
plot_configs <- list(
  list(var = "Gender", title = "Gender Distribution", filter = NULL),
  list(var = "Age_Group", title = "Age (>60) Distribution", filter = NULL),
  list(var = "Batch", title = "Batch Distribution", filter = NULL),
  list(var = "Clinical.response.at.OG.MDT..responder.non.responder", title = "Response Distribution", filter = NULL)
)

####################
# 10) Combined cohort report
####################
p1 <- plot_clinical_assoc(cell_df, "Gender", "Gender")
p2 <- plot_clinical_assoc(cell_df, "Age_Group", "Age (>60)")
p3 <- plot_clinical_assoc(cell_df, "Batch", "Batch")
p4 <- plot_clinical_assoc(cell_df, "Clinical.response.at.OG.MDT..responder.non.responder", "Clinical Response")

pdf("Auto_clinical_assoc_topmp_v2B_combined.pdf", width = 14, height = 16)
if (!is.null(p1) && !is.null(p2) && !is.null(p3) && !is.null(p4)) {
  print(((p1 + p2) / (p3 + p4)) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = "PDO Clinical Variable Overview"))
} else {
  for (cfg in plot_configs) {
    p <- plot_clinical_assoc(cell_df, cfg$var, cfg$title, cfg$filter)
    if (!is.null(p)) print(p)
  }
}
dev.off()

####################
# 11) Per-Batch report
####################
pdf("Auto_clinical_assoc_topmp_v2B_per_batch.pdf", width = 16, height = 10)
for (cfg in plot_configs) {
  p <- plot_variable_facet(
    data = cell_df,
    group_var = cfg$var,
    filter_expr = cfg$filter,
    title = cfg$title
  )
  if (!is.null(p)) print(p)
}
dev.off()

####################
# 12) Save machine-readable summary
####################
summary_rows <- lapply(plot_configs, function(cfg) {
  prep <- compute_plot_data(
    data = cell_df,
    group_var = cfg$var,
    filter_expr = cfg$filter,
    batchname = NULL
  )

  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) {
    return(NULL)
  }

  prep$plot_data %>%
    left_join(prep$label_data, by = cfg$var) %>%
    transmute(
      clinical_variable = cfg$var,
      level = .data[[cfg$var]],
      state = as.character(state),
      mean_sample_prop_pct = pct,
      samples_in_level = n_samples,
      total_samples = prep$total_samples
    )
})

summary_df <- bind_rows(summary_rows)
summary_path <- "Auto_clinical_assoc_topmp_v2B_summary.csv"
write.csv(summary_df, summary_path, row.names = FALSE)

message("Saved: Auto_clinical_assoc_topmp_v2B_combined.pdf")
message("Saved: Auto_clinical_assoc_topmp_v2B_per_batch.pdf")
message(sprintf("Saved: %s", summary_path))
