####################
# Adapted from scRef_Pipeline for PDOs
# Auto_clinical_mp_ucell_plots.R
####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(scales)
library(patchwork)
library(ggnewscale)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# 1) Load data
####################
tmdata_all <- readRDS("PDOs_merged.rds")
ucell_scores <- readRDS("UCell_scores_filtered.rds") # UPDATE: Verify filename
geneNMF.metaprograms <- readRDS("MP_outs_default.rds") # UPDATE: Change to optimal nMP result

clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx"
)

####################
# 2) MP order and placeholder annotations
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
}
retained_mps <- names(mp.genes)

####################
# UPDATE AFTER geneNMF: Fill in MP descriptions and state groups
####################
mp_descriptions <- c(
  "MP6"  = "G2M Cell Cycle",
  "MP7"  = "DNA repair",
  "MP5"  = "MYC-related Proliferation",
  "MP1"  = "G2M checkpoint",
  "MP3"  = "G1S Cell Cycle",
  "MP8"  = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9"  = "ECM Remodeling Epi.",
  "MP4"  = "Intestinal Metaplasia"
)
state_groups <- list(
  "Classic_Proliferative" = c("MP5"),
  "Columnar_Progenitor"   = c("MP8"),
  "EMT_related"           = c("MP10", "MP9"),
  "Intestinal_Metaplasia" = c("MP4"),
  "3CA_EMT_and_Protein_maturation" = character(0)
)
# Auto-populate if empty
if (length(mp_descriptions) == 0) mp_descriptions <- setNames(retained_mps, retained_mps)

####################
# UPDATE AFTER geneNMF: Cell-cycle MP list
####################
cc_mps <- c()  # UPDATE: identify cell cycle MPs

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)

mp_cols <- intersect(mp_tree_order_names, colnames(ucell_scores))
mp_labels <- mp_descriptions
mp_labels[setdiff(mp_cols, names(mp_labels))] <- setdiff(mp_cols, names(mp_labels))

mp_display <- setNames(
  paste0(names(mp_labels[mp_cols]), " ", mp_labels[mp_cols]),
  mp_cols
)

if (length(state_groups) == 0) {
  state_groups <- list("All MPs" = mp_cols)
}

mp_to_group <- setNames(rep("Other", length(mp_cols)), mp_cols)
for (grp in names(state_groups)) {
  for (mp in intersect(state_groups[[grp]], mp_cols)) {
    mp_to_group[mp] <- grp
  }
}

group_order <- names(state_groups)
mp_ordered <- c()
for (grp in group_order) {
  grp_mps <- mp_cols[mp_cols %in% state_groups[[grp]]]
  mp_ordered <- c(mp_ordered, grp_mps)
}
remaining <- setdiff(mp_cols, mp_ordered)
mp_ordered <- c(mp_ordered, remaining)

####################
# 3) Cell alignment and metadata (PDO-specific)
####################
meta_df <- tmdata_all@meta.data
meta_df$cell <- rownames(meta_df)
meta_df <- meta_df %>%
  mutate(
    SUR = ifelse("SUR" %in% colnames(meta_df), as.character(SUR), sub("^(SUR[0-9]+).*", "\\1", as.character(orig.ident))),
    Batch = as.character(Batch),
    orig.ident = as.character(orig.ident)
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

cell_ids <- Reduce(intersect, list(rownames(ucell_scores), rownames(meta_df)))

cell_meta <- meta_df[cell_ids, , drop = FALSE] %>%
  left_join(clinical_df, by = "SUR") %>%
  filter(orig.ident != "SUR843T3_PDO")

age_vec <- if ("Age" %in% colnames(cell_meta)) {
  suppressWarnings(as.numeric(cell_meta$Age))
} else if ("clinical_Age" %in% colnames(cell_meta)) {
  suppressWarnings(as.numeric(cell_meta$clinical_Age))
} else {
  rep(NA_real_, nrow(cell_meta))
}

gender_vec <- if ("Gender" %in% colnames(cell_meta)) {
  as.character(cell_meta$Gender)
} else if ("clinical_Gender" %in% colnames(cell_meta)) {
  as.character(cell_meta$clinical_Gender)
} else {
  rep(NA_character_, nrow(cell_meta))
}

response_vec <- if ("Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(cell_meta)) {
  as.character(cell_meta$Clinical.response.at.OG.MDT..responder.non.responder)
} else if ("clinical_Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(cell_meta)) {
  as.character(cell_meta$clinical_Clinical.response.at.OG.MDT..responder.non.responder)
} else {
  rep(NA_character_, nrow(cell_meta))
}

cell_meta <- cell_meta %>%
  mutate(
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
    ),
    Batch = as.character(Batch)
  )

cell_long <- as.data.frame(ucell_scores[cell_meta$cell, mp_cols, drop = FALSE]) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(mp_cols), names_to = "MP", values_to = "score") %>%
  left_join(cell_meta, by = "cell")

####################
# 4) Core aggregation helper
####################
compute_mp_group_means <- function(data, group_var, filter_expr = NULL, batchname = NULL) {
  if (!is.null(batchname)) data <- data %>% filter(Batch == batchname)
  if (!is.null(filter_expr)) data <- data %>% filter(!!rlang::parse_expr(filter_expr))

  data <- data %>%
    filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(MP), !is.na(score))

  if (nrow(data) == 0) return(list(plot_data = NULL, label_data = NULL, total_samples = 0L))

  sample_mean <- data %>%
    group_by(orig.ident, .data[[group_var]], MP) %>%
    summarise(sample_mean = mean(score, na.rm = TRUE), .groups = "drop")

  group_mean <- sample_mean %>%
    group_by(.data[[group_var]], MP) %>%
    summarise(group_mean = mean(sample_mean, na.rm = TRUE), .groups = "drop")

  total_samples <- data %>% distinct(orig.ident) %>% nrow()
  n_sample_group <- data %>%
    distinct(orig.ident, .data[[group_var]]) %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_samples = n_distinct(orig.ident), .groups = "drop")

  n_cell_group <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(n_cells = n_distinct(cell), .groups = "drop")

  label_data <- n_sample_group %>%
    left_join(n_cell_group, by = group_var) %>%
    mutate(label = paste0(.data[[group_var]], "\n(n=", n_samples, " samples; ", comma(n_cells), " cells)"))

  list(plot_data = group_mean, label_data = label_data, total_samples = total_samples)
}

####################
# 5) Dot-plot heatmap
####################
plot_dot_heatmap <- function(data, group_var, title, filter_expr = NULL, batchname = NULL) {
  prep <- compute_mp_group_means(data, group_var, filter_expr, batchname)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)

  plot_df <- prep$plot_data %>%
    group_by(MP) %>%
    mutate(
      mp_mean = mean(group_mean, na.rm = TRUE),
      mp_sd = sd(group_mean, na.rm = TRUE),
      z_score = ifelse(mp_sd > 0, (group_mean - mp_mean) / mp_sd, 0)
    ) %>%
    ungroup()

  z_lim <- max(abs(plot_df$z_score), na.rm = TRUE)
  z_lim <- max(z_lim, 0.5)

  plot_df$mp_display <- factor(mp_display[plot_df$MP], levels = rev(mp_display[mp_ordered]))
  plot_df$mp_group <- factor(mp_to_group[plot_df$MP], levels = names(state_groups))

  level_labels <- prep$label_data$label
  names(level_labels) <- as.character(prep$label_data[[group_var]])
  plot_df$x_label <- level_labels[as.character(plot_df[[group_var]])]
  plot_df$x_label <- factor(plot_df$x_label, levels = unique(level_labels))

  n_x <- length(unique(plot_df$x_label))
  plot_df$mp_short <- sub("^MP\\d+\\s+", "", as.character(plot_df$mp_display))
  mp_short_levels <- sub("^MP\\d+\\s+", "", rev(mp_display[mp_ordered]))
  plot_df$mp_short <- factor(plot_df$mp_short, levels = unique(mp_short_levels))

  ggplot(plot_df, aes(x = x_label, y = mp_short)) +
    geom_point(aes(fill = z_score), size = 7, shape = 21, colour = "grey30", stroke = 0.4) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-z_lim, z_lim), name = "Relative\nactivity\n(Z-score)"
    ) +
    {if (n_x <= 5) {
      geom_text(aes(label = sprintf("%.3f", group_mean)), size = 2.0, colour = "grey30", nudge_y = -0.3)
    }} +
    facet_grid(mp_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL, title = title, subtitle = paste0("n=", prep$total_samples, " samples")) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16, margin = margin(b = 0)),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 10, colour = "grey45", margin = margin(b = 8)),
      axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 10, lineheight = 0.9),
      axis.text.y = element_text(size = 9.5),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(0.3, "lines"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8, hjust = 1, margin = margin(r = 0)),
      strip.background = element_rect(fill = "grey97", colour = NA),
      legend.position = "right"
    )
}

####################
# 6) Lollipop difference plot
####################
plot_lollipop_diff <- function(data, group_var, title, filter_expr = NULL) {
  prep <- compute_mp_group_means(data, group_var, filter_expr)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)
  if (length(unique(prep$plot_data[[group_var]])) != 2) return(NULL)

  levels_vec <- sort(unique(as.character(prep$plot_data[[group_var]])))
  ref_level <- levels_vec[1]
  alt_level <- levels_vec[2]

  wide_df <- prep$plot_data %>%
    mutate(level = as.character(.data[[group_var]])) %>%
    dplyr::select(MP, level, group_mean) %>%
    pivot_wider(names_from = level, values_from = group_mean)

  wide_df$diff <- wide_df[[alt_level]] - wide_df[[ref_level]]
  wide_df$direction <- ifelse(wide_df$diff > 0, paste0("Higher in ", alt_level), paste0("Higher in ", ref_level))
  wide_df$mp_display <- factor(mp_display[wide_df$MP], levels = rev(mp_display[mp_ordered]))

  ggplot(wide_df, aes(x = diff, y = mp_display, colour = direction)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.5, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = diff, y = mp_display, yend = mp_display), linewidth = 0.8) +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("%+.4f", diff)), hjust = ifelse(wide_df$diff >= 0, -0.3, 1.3), size = 2.8, show.legend = FALSE) +
    scale_colour_manual(
      values = setNames(c("#B2182B", "#2166AC"), c(paste0("Higher in ", alt_level), paste0("Higher in ", ref_level))),
      name = "Direction"
    ) +
    labs(
      x = paste0("Δ Mean UCell score (", alt_level, " − ", ref_level, ")"),
      y = NULL,
      title = paste0(title, " — Difference"),
      subtitle = paste0("n=", prep$total_samples, " samples")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 10, colour = "grey40"),
      axis.text.y = element_text(size = 9, face = "italic"),
      panel.grid.major.y = element_line(colour = "grey95"),
      panel.grid.major.x = element_line(colour = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

####################
# 7) Faceted per-Batch dot plot
####################
plot_dot_faceted <- function(data, group_var, title, filter_expr = NULL) {
  if (!is.null(filter_expr)) data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  data <- data %>% filter(!is.na(.data[[group_var]]), !is.na(orig.ident), !is.na(Batch), !is.na(MP), !is.na(score))
  if (nrow(data) == 0) return(NULL)

  total_samples <- n_distinct(data$orig.ident)

  sample_mean <- data %>%
    group_by(Batch, orig.ident, .data[[group_var]], MP) %>%
    summarise(sample_mean = mean(score, na.rm = TRUE), .groups = "drop")

  plot_data <- sample_mean %>%
    group_by(Batch, .data[[group_var]], MP) %>%
    summarise(group_mean = mean(sample_mean, na.rm = TRUE), .groups = "drop") %>%
    group_by(Batch, MP) %>%
    mutate(mp_mean = mean(group_mean, na.rm = TRUE), mp_sd = sd(group_mean, na.rm = TRUE), z_score = ifelse(mp_sd > 0, (group_mean - mp_mean) / mp_sd, 0)) %>%
    ungroup()

  z_lim <- max(abs(plot_data$z_score), na.rm = TRUE)
  z_lim <- max(z_lim, 0.5)

  batch_stats <- data %>%
    group_by(Batch) %>%
    summarise(n_cells = n_distinct(cell), n_samples = n_distinct(orig.ident), .groups = "drop") %>%
    mutate(batch_label = paste0(Batch, " (n=", n_samples, "; ", comma(n_cells), " cells)"))

  plot_data <- plot_data %>% left_join(batch_stats[, c("Batch", "batch_label")], by = "Batch")
  plot_data$mp_display <- factor(mp_display[plot_data$MP], levels = rev(mp_display[mp_ordered]))

  ggplot(plot_data, aes(x = .data[[group_var]], y = mp_display)) +
    geom_point(aes(fill = z_score), size = 5.5, shape = 21, colour = "grey30", stroke = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-z_lim, z_lim), name = "Z") +
    geom_text(aes(label = sprintf("%.3f", group_mean)), size = 1.3, colour = "grey20") +
    facet_wrap(~ batch_label, scales = "free_x", ncol = 3) +
    labs(x = NULL, y = NULL, title = title, subtitle = paste0("n=", total_samples, " samples")) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 10, colour = "grey45", margin = margin(b = 8)),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7, face = "italic"),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey96", colour = NA),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right"
    )
}

####################
# 8) Variables (PDO clinical columns)
####################
plot_configs <- list(
  list(var = "Gender", title = "Gender", filter = NULL),
  list(var = "Age_Group", title = "Age (>60)", filter = NULL),
  list(var = "Batch", title = "Batch", filter = NULL),
  list(var = "Clinical.response.at.OG.MDT..responder.non.responder", title = "Clinical Response", filter = NULL)
)

####################
# 9) Combined report
####################
cat("Generating combined report...\n")
pdf("Auto_clinical_assoc_mp_ucell_combined.pdf", width = 14, height = 10)
for (cfg in plot_configs) {
  p_dot <- tryCatch(
    plot_dot_heatmap(cell_long, cfg$var, cfg$title, cfg$filter),
    error = function(e) { message("  Skipping dot plot for ", cfg$title, ": ", e$message); NULL }
  )
  if (!is.null(p_dot)) print(p_dot)

  p_lol <- tryCatch(
    plot_lollipop_diff(cell_long, cfg$var, cfg$title, cfg$filter),
    error = function(e) { message("  Skipping lollipop for ", cfg$title, ": ", e$message); NULL }
  )
  if (!is.null(p_lol)) print(p_lol)
}
dev.off()

####################
# 10) Per-Batch report
####################
cat("Generating per-Batch report...\n")
pdf("Auto_clinical_assoc_mp_ucell_per_batch.pdf", width = 18, height = 12)
for (cfg in plot_configs) {
  p <- tryCatch(
    plot_dot_faceted(cell_long, cfg$var, paste0(cfg$title, " — per Batch"), cfg$filter),
    error = function(e) { message("  Skipping faceted plot for ", cfg$title, ": ", e$message); NULL }
  )
  if (!is.null(p)) print(p)
}
dev.off()

####################
# 11) Summary CSV
####################
summary_rows <- lapply(plot_configs, function(cfg) {
  prep <- compute_mp_group_means(cell_long, cfg$var, cfg$filter)
  if (is.null(prep$plot_data) || nrow(prep$plot_data) == 0) return(NULL)

  prep$plot_data %>%
    mutate(MP_label = mp_labels[MP]) %>%
    left_join(prep$label_data, by = cfg$var) %>%
    transmute(
      clinical_variable = cfg$var,
      level = .data[[cfg$var]],
      MP = MP,
      MP_label = MP_label,
      MP_group = mp_to_group[MP],
      mean_of_sample_means = group_mean,
      samples_in_level = n_samples,
      cells_in_level = n_cells,
      total_samples = prep$total_samples
    )
})

summary_df <- bind_rows(summary_rows)
summary_path <- "Auto_clinical_assoc_mp_ucell_summary.csv"
write.csv(summary_df, summary_path, row.names = FALSE)

message("Saved: Auto_clinical_assoc_mp_ucell_combined.pdf")
message("Saved: Auto_clinical_assoc_mp_ucell_per_batch.pdf")
message(sprintf("Saved: %s", summary_path))
