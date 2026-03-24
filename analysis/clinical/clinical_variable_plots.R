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
library(stringr)

####################
# 1) Setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# UPDATE: Standardised MP and State info
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
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "SMG-like Metaplasia"   = c("MP8")
)
# Cell-cycle MP list
cc_mps <- c("MP6", "MP7", "MP1", "MP3")

####################
# 3) Load PDO clinical xlsx
####################
clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx"
)

# Transpose clinical data carefully
# clinical_sheet has columns = patients, first column = variable names
clinical_vars <- as.character(clinical_sheet[[1]])
clinical_vars[is.na(clinical_vars)] <- paste0("VAR_", which(is.na(clinical_vars)))
# Clean names for columns
clean_vars <- gsub("[^[:alnum:]]", "_", clinical_vars)
clean_vars <- gsub("_+", "_", clean_vars)
clean_vars <- gsub("^_|_$", "", clean_vars)

clinical_data <- as.data.frame(t(clinical_sheet[, -1]))
colnames(clinical_data) <- clean_vars

# Set patient ID from rownames (headers in original sheet)
clinical_data$SUR <- rownames(clinical_data)
# Ensure SUR is character and has prefix if it only contains numbers
clinical_data$SUR <- ifelse(grepl("^[0-9]+$", clinical_data$SUR), paste0("SUR", clinical_data$SUR), clinical_data$SUR)

# Prefix all clinical columns with clinical_
colnames(clinical_data) <- ifelse(colnames(clinical_data) == "SUR", "SUR", paste0("clinical_", colnames(clinical_data)))
clinical_df <- clinical_data

message("Loading finalized states...")
tmdata_all <- readRDS("PDOs_final.rds")
state_final <- readRDS("Auto_PDO_final_states.rds")
tmdata_all$state <- as.character(state_final[Cells(tmdata_all)])
message("Unique states in tmdata_all:")
print(unique(tmdata_all$state))
  
meta_df <- tmdata_all@meta.data
meta_df$cell <- rownames(meta_df)

####################
# 4) Build cell-level table
####################
cell_df <- meta_df %>%
  mutate(
    orig.ident = as.character(orig.ident),
    SUR = sub("^(SUR[0-9]+).*", "\\1", as.character(orig.ident)),
    Batch = as.character(Batch)
  ) %>%
  left_join(clinical_df, by = "SUR") %>%
  filter(orig.ident != "SUR843T3_PDO")

# Force state column to stay correct character strings
cell_df$state <- as.character(cell_df$state)
message("States in cell_df after join:")
print(table(cell_df$state, useNA="always"))

# ... rest of mutate (Age, Gender, etc.) remains ...
cell_df <- cell_df %>%
  mutate(
    Age = suppressWarnings(as.numeric(clinical_Age)),
    # fallback
    Age = ifelse(is.na(Age), suppressWarnings(as.numeric(meta_df$Age[match(cell, meta_df$cell)])), Age),
    Age_Group = ifelse(!is.na(Age) & Age > 60, ">60", "<=60"),
    Gender = ifelse(clinical_Gender %in% c("F", "female", "Female"), "Female", 
                    ifelse(!is.na(clinical_Gender), "Male", NA_character_)),
    Clinical.response.at.OG.MDT..responder.non.responder = case_when(
      clinical_Clinical_response_at_OG_MDT_responder_non_responder == "Responder" ~ "Responder",
      clinical_Clinical_response_at_OG_MDT_responder_non_responder == "Non-responder" ~ "Nonresponder",
      TRUE ~ as.character(clinical_Clinical_response_at_OG_MDT_responder_non_responder)
    ),
    Gender = factor(Gender, levels = c("Male", "Female")),
    Clinical.response.at.OG.MDT..responder.non.responder = factor(Clinical.response.at.OG.MDT..responder.non.responder, levels = c("Responder", "Nonresponder")),
    Batch = factor(Batch, levels = c("Untreated_PDO", "Treated_PDO", "PDO")),
    T_Stage = str_extract(clinical_cTNM_stage_at_diagnosis, "T[0-4][a-b]?"),
    T_Stage = ifelse(is.na(T_Stage), "Unknown", T_Stage),
    Collection_Timepoint = factor(clinical_Collection_timepoint, levels = c("Pre", "Post")),
    Histology_at_Surgery = ifelse(grepl("adenocarcinoma", clinical_Tumour_histology_at_surgery, ignore.case=T), "Adenocarcinoma",
                                  ifelse(grepl("squamous", clinical_Tumour_histology_at_surgery, ignore.case=T), "Adenosquamous",
                                         as.character(clinical_Tumour_histology_at_surgery))),
    Tumour_Type_Grouped = case_when(
      clinical_Tumour_type_Oesophageal_GOJ_type_I_III_gastric %in% c("Distal", "multiple lesions - mid and distal") ~ "Distal",
      clinical_Tumour_type_Oesophageal_GOJ_type_I_III_gastric %in% c("GOJ 1", "GOJ 2", "GOJ 3", "GOJ3", "Oesophageal") ~ "Esophageal/GOJ",
      clinical_Tumour_type_Oesophageal_GOJ_type_I_III_gastric == "Gastric" ~ "Gastric",
      TRUE ~ as.character(clinical_Tumour_type_Oesophageal_GOJ_type_I_III_gastric)
    )
  )
  
state_levels <- c(
  "Classic Proliferative", "Basal to Intest. Meta", "Stress-adaptive", "SMG-like Metaplasia",
  "3CA_EMT_and_Protein_maturation", "Unresolved", "Hybrid"
)
state_colors <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

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

  group_vars_sample <- unique(c("orig.ident", group_var, "state"))
  sample_prop <- data %>%
    group_by(across(all_of(group_vars_sample))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(unique(c("orig.ident", group_var))))) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()

  mean_prop <- sample_prop %>%
    group_by(across(all_of(unique(c(group_var, "state"))))) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop")

  group_sum <- mean_prop %>%
    group_by(across(all_of(group_var))) %>%
    summarise(group_total = sum(mean_prop), .groups = "drop")

  plot_data <- mean_prop %>%
    left_join(group_sum, by = unique(group_var)) %>%
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
    scale_x_discrete(labels = function(x) stringr::str_wrap(iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"), width = 15)) +
    theme_minimal(base_size = 14) +
    labs(x = NULL, y = "Proportion", title = title, fill = "Cell State") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 16)
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

  group_vars_sample <- unique(c("Batch", "orig.ident", group_var, "state"))
  sample_prop <- data %>%
    group_by(across(all_of(group_vars_sample))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(unique(c("Batch", "orig.ident", group_var))))) %>%
    mutate(sample_prop = n / sum(n)) %>%
    ungroup()

  mean_prop <- sample_prop %>%
    group_by(across(all_of(unique(c("Batch", group_var, "state"))))) %>%
    summarise(mean_prop = mean(sample_prop, na.rm = TRUE), .groups = "drop")

  group_sum <- mean_prop %>%
    group_by(across(all_of(unique(c("Batch", group_var))))) %>%
    summarise(group_total = sum(mean_prop), .groups = "drop")

  plot_data <- mean_prop %>%
    left_join(group_sum, by = unique(c("Batch", group_var))) %>%
    mutate(freq = ifelse(group_total > 0, mean_prop / group_total, 0), state = factor(state, levels = state_levels))

  batch_stats <- data %>%
    group_by(Batch) %>%
    summarise(n_cells = n(), n_samples = n_distinct(orig.ident), .groups = "drop") %>%
    mutate(batch_label = paste0(Batch, "\n(N=", comma(n_cells), "; samples=", n_samples, ")"))

  label_data <- data %>%
    distinct(across(all_of(unique(c("Batch", "orig.ident", group_var))))) %>%
    group_by(across(all_of(unique(c("Batch", group_var))))) %>%
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
  list(var = "Clinical.response.at.OG.MDT..responder.non.responder", title = "Response Distribution", filter = NULL),
  list(var = "Tumour_Type_Grouped", title = "Tumour Type Distribution", filter = NULL),
  list(var = "Histology_at_Surgery", title = "Histology at Surgery", filter = "Histology_at_Surgery != 'N/A'"),
  list(var = "T_Stage", title = "cTNM T-Stage Distribution", filter = NULL),
  list(var = "Collection_Timepoint", title = "Collection Timepoint", filter = NULL)
)

# Multi-page PDF layout - wider for Better visibility
pdf("Auto_clinical_assoc_topmp_v2B_combined.pdf", width = 18, height = 14)

# Page 1: Demographics & Batch
p1 <- plot_clinical_assoc(cell_df, "Gender", "Gender")
p2 <- plot_clinical_assoc(cell_df, "Age_Group", "Age (>60)")
p3 <- plot_clinical_assoc(cell_df, "Batch", "Batch")
p4 <- plot_clinical_assoc(cell_df, "Clinical.response.at.OG.MDT..responder.non.responder", "Clinical Response")

if (!is.null(p1) && !is.null(p2) && !is.null(p3) && !is.null(p4)) {
  print(((p1 + p2) / (p3 + p4)) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = "PDO Clinical Variable Overview - Page 1"))
}

# Page 2: Tumour Characteristics
p5 <- plot_clinical_assoc(cell_df, "Tumour_Type_Grouped", "Tumour Type")
p6 <- plot_clinical_assoc(cell_df, "Histology_at_Surgery", "Histology at Surgery", "Histology_at_Surgery != 'N/A'")
p7 <- plot_clinical_assoc(cell_df, "T_Stage", "cTNM T-Stage")
p8 <- plot_clinical_assoc(cell_df, "Collection_Timepoint", "Collection Timepoint")

if (!is.null(p5) && !is.null(p6) && !is.null(p7) && !is.null(p8)) {
  print(((p5 + p6) / (p7 + p8)) + plot_layout(heights = c(1, 1.2)) + plot_annotation(title = "PDO Clinical Variable Overview - Page 2"))
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
