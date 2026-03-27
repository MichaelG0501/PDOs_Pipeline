####################
# Auto_sample_abundance_pdo.R
#
# PDO Sample Abundance with clinical annotations
#
# Requested fixes:
# 1) add gap between annotation block and main plot
# 2) keep ONLY state / MP legend at bottom
# 3) keep clinical annotation legends on the right
# 4) bottom state legend max 3 entries per row
# 5) add small separation between samples in annotation bars
# 6) slightly larger font overall relative to plot
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(readxl)
  library(stringr)
  library(RColorBrewer)
  library(grid)
  library(rlang)
})

####################
# 1) Setup & Load
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "sample_abundance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading data...")
pdos <- readRDS("PDOs_merged.rds")

state_path <- "Auto_PDO_final_states.rds"
if (!file.exists(state_path)) {
  state_path <- "Auto_PDO_states_noreg.rds"
}
state_B <- readRDS(state_path)

ucell_scores <- readRDS("UCell_scores_filtered.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")

####################
# 2) Standardize Metadata
####################
message("Processing clinical data...")

pdos$Batch <- case_when(
  grepl("_Treated_|_Untreated_", pdos$orig.ident) ~ "New_batch",
  TRUE ~ "Cynthia_batch"
)

clinical_sheet <- read_excel(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx"
)

vars_raw <- as.character(clinical_sheet[[1]])
na_idx <- is.na(vars_raw)
vars_raw[na_idx] <- paste0("VAR", seq_len(sum(na_idx)))

clean_vars <- gsub("[^[:alnum:]]", "_", vars_raw)
clean_vars <- gsub("_+", "_", clean_vars)
clean_vars <- gsub("^_|_$", "", clean_vars)

clinical_data <- as.data.frame(t(clinical_sheet[, -1]), stringsAsFactors = FALSE)
colnames(clinical_data) <- clean_vars
clinical_data$SUR <- rownames(clinical_data)
clinical_data$SUR <- ifelse(
  grepl("^[0-9]+$", clinical_data$SUR),
  paste0("SUR", clinical_data$SUR),
  clinical_data$SUR
)
colnames(clinical_data) <- ifelse(
  colnames(clinical_data) == "SUR",
  "SUR",
  paste0("clinical_", colnames(clinical_data))
)

meta_df <- pdos@meta.data
meta_df$cell <- rownames(meta_df)

meta_df <- meta_df %>%
  mutate(
    SUR = sub("^(SUR[0-9]+).*", "\\1", as.character(orig.ident))
  ) %>%
  left_join(clinical_data, by = "SUR")

meta_df <- meta_df %>%
  mutate(
    Age_raw = suppressWarnings(as.numeric(clinical_Age)),
    Age_Group = case_when(
      is.na(Age_raw) ~ "Unknown",
      Age_raw > 60 ~ ">60",
      TRUE ~ "<=60"
    ),
    Gender = case_when(
      clinical_Gender %in% c("F", "female", "Female") ~ "Female",
      clinical_Gender %in% c("M", "male", "Male") ~ "Male",
      TRUE ~ "Unknown"
    ),
    Response = case_when(
      grepl(
        "Non-responder|Nonresponder",
        clinical_Clinical_response_at_OG_MDT_responder_non_responder,
        ignore.case = TRUE
      ) ~ "Nonresponder",
      grepl(
        "Responder",
        clinical_Clinical_response_at_OG_MDT_responder_non_responder,
        ignore.case = TRUE
      ) ~ "Responder",
      TRUE ~ "Unknown"
    ),
    T_Stage = str_extract(clinical_cTNM_stage_at_diagnosis, "T[0-4]"),
    T_Stage = ifelse(is.na(T_Stage), "Unknown", T_Stage),
    Timepoint = case_when(
      clinical_Collection_timepoint %in% c("Pre", "Post") ~ clinical_Collection_timepoint,
      TRUE ~ "Unknown"
    ),
    Histology = case_when(
      grepl("adenocarcinoma", clinical_Tumour_histology_at_surgery, ignore.case = TRUE) ~ "Adeno",
      grepl("squamous", clinical_Tumour_histology_at_surgery, ignore.case = TRUE) ~ "Squamous",
      TRUE ~ "Unknown"
    ),
    Tumour_Type = case_when(
      clinical_Tumour_type_Oesophageal_GOJ_type_I_III_gastric %in% c("Distal", "multiple lesions - mid and distal") ~ "Distal",
      clinical_Tumour_type_Oesophageal_GOJ_type_I_III_gastric %in% c("Gastric") ~ "Gastric",
      TRUE ~ "E/GOJ"
    )
  )

####################
# 3) Constants
####################
mp_desc <- c(
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

group_cols <- c(
  "Classic Proliferative"          = "#E41A1C",
  "Basal to Intest. Meta"          = "#4DAF4A",
  "Stress-adaptive"                = "#984EA3",
  "SMG-like Metaplasia"            = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved"                     = "grey80",
  "Hybrid"                         = "black"
)

mp_cols <- c(
  "MP6_G2M Cell Cycle"             = "#E78AC3",
  "MP7_DNA repair"                 = "#999999",
  "MP5_MYC-related Proliferation"  = "#E41A1C",
  "MP1_G2M checkpoint"             = "#B3B3B3",
  "MP3_G1S Cell Cycle"             = "#8DA0CB",
  "MP8_Columnar Progenitor"        = "#FF7F00",
  "MP10_Inflammatory Stress Epi."  = "#984EA3",
  "MP9_ECM Remodeling Epi."        = "#A6D854",
  "MP4_Intestinal Metaplasia"      = "#4DAF4A"
)

cl_cols <- list(
  Batch = c(
    "Cynthia_batch" = "brown",
    "New_batch" = "darkgreen"
  ),
  Gender = c(
    "Male" = "#0072B2",
    "Female" = "#CC79A7",
    "Unknown" = "grey85"
  ),
  Age_Group = c(
    "<=60" = "#F0E442",
    ">60" = "#D55E00",
    "Unknown" = "grey85"
  ),
  Response = c(
    "Responder" = "#009E73",
    "Nonresponder" = "red",
    "Unknown" = "grey85"
  ),
  Tumour_Type = c(
    "Distal" = "#66C2A5",
    "E/GOJ" = "#FC8D62",
    "Gastric" = "#8DA0CB"
  ),
  Histology = c(
    "Adeno" = "#E78AC3",
    "Squamous" = "#A6D854",
    "Unknown" = "grey85"
  ),
  T_Stage = c(
    "T1" = "#FC8D62",
    "T2" = "#8DA0CB",
    "T3" = "#E78AC3",
    "T4" = "#A6D854",
    "Unknown" = "grey85"
  ),
  Timepoint = c(
    "Pre" = "#66C2A5",
    "Post" = "#FC8D62",
    "Unknown" = "grey85"
  )
)

####################
# 4) Filter / Align
####################
mp_genes <- geneNMF.metaprograms$metaprograms.genes
coverage <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage) <- paste0("MP", seq_along(coverage))
silhouette <- geneNMF.metaprograms$metaprograms.metrics$silhouette

retained_mps <- setdiff(
  names(mp_genes),
  c(
    paste0("MP", which(silhouette < 0)),
    names(coverage)[coverage < 0.25]
  )
)

cc_list <- intersect(c("MP6", "MP7", "MP1", "MP3"), retained_mps)
non_cc_list <- setdiff(retained_mps, cc_list)

common <- Reduce(
  intersect,
  list(rownames(ucell_scores), Cells(pdos), names(state_B))
)

pdos <- pdos[, common]
ucell_scores <- ucell_scores[common, , drop = FALSE]
state_B <- state_B[common]
meta_df <- meta_df %>% filter(cell %in% common)

####################
# 5) Label helpers
####################
label_f <- function(mps) {
  d <- mp_desc[mps]
  d[is.na(d)] <- mps[is.na(d)]
  out <- paste0(mps, "_", d)
  names(out) <- names(mps)
  out
}

topmp_noncc_vec <- colnames(ucell_scores[, non_cc_list, drop = FALSE])[
  max.col(ucell_scores[, non_cc_list, drop = FALSE], ties.method = "first")
]
names(topmp_noncc_vec) <- rownames(ucell_scores)
topmp_noncc <- label_f(topmp_noncc_vec)

topmp_all_vec <- colnames(ucell_scores[, retained_mps, drop = FALSE])[
  max.col(ucell_scores[, retained_mps, drop = FALSE], ties.method = "first")
]
names(topmp_all_vec) <- rownames(ucell_scores)
topmp_all <- label_f(topmp_all_vec)

states <- as.character(state_B)
names(states) <- names(state_B)

####################
# 6) Sort orders
####################
sample_info <- meta_df %>%
  distinct(orig.ident, Batch, Response, Tumour_Type, Timepoint)

sort_list <- list(
  "Sort: Batch" = sample_info %>%
    arrange(Batch, Response, Tumour_Type, Timepoint, orig.ident) %>%
    pull(orig.ident),
  
  "Sort: Response" = sample_info %>%
    arrange(Response, Batch, Tumour_Type, Timepoint, orig.ident) %>%
    pull(orig.ident),
  
  "Sort: Tumour Type" = sample_info %>%
    arrange(Tumour_Type, Batch, Response, Timepoint, orig.ident) %>%
    pull(orig.ident),
  
  "Sort: Timepoint" = sample_info %>%
    arrange(Timepoint, Batch, Response, Tumour_Type, orig.ident) %>%
    pull(orig.ident)
)

####################
# 7) Orders / colours
####################
state_mps_full <- c("MP5", "MP4", "MP10", "MP9", "MP8")
state_mps <- label_f(state_mps_full)

cc_lbls <- label_f(cc_list)

all_nocc <- sort(unique(topmp_noncc_vec))
all_nocc <- label_f(all_nocc)
ord1 <- c(
  state_mps[state_mps %in% all_nocc],
  setdiff(all_nocc, state_mps)
)

all_all <- sort(unique(topmp_all_vec))
all_all <- label_f(all_all)
ord2 <- c(
  cc_lbls[cc_lbls %in% all_all],
  state_mps[state_mps %in% all_all],
  setdiff(all_all, c(cc_lbls, state_mps))
)

ord3 <- intersect(names(group_cols), unique(states))

c1 <- setNames(scales::hue_pal()(length(all_nocc)), all_nocc)
for (n in names(mp_cols)) {
  if (n %in% names(c1)) c1[n] <- mp_cols[n]
}

c2 <- setNames(scales::hue_pal()(length(all_all)), all_all)
for (n in names(mp_cols)) {
  if (n %in% names(c2)) c2[n] <- mp_cols[n]
}

####################
# 8) Theme helpers
####################
anno_base_theme <- theme_void(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 19, face = "bold", color = "black"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(0, 4, 0, 4)
  )

right_legend_theme <- theme(
  legend.position = "right",
  legend.direction = "vertical",
  legend.box = "vertical",
  legend.key.size = unit(0.80, "cm"),
  legend.text = element_text(size = 22),
  legend.title = element_text(size = 26, face = "bold"),
  legend.margin = margin(t = 20, b = 20, l = 10, r = 10),
  legend.box.margin = margin(t = 5, b = 5),
  legend.spacing.y = unit(0.85, "cm")
)

bottom_state_legend_theme <- theme(
  legend.position = "right",
  legend.direction = "vertical",
  legend.box = "vertical",
  legend.key.size = unit(0.80, "cm"),
  legend.text = element_text(size = 22),
  legend.title = element_text(size = 26, face = "bold"),
  legend.margin = margin(20, 10, 20, 10),
  legend.box.margin = margin(2, 2, 2, 2),
  legend.spacing.y = unit(0.85, "cm")
)

####################
# 9) Annotation panel
####################
make_anno <- function(var, label, s_ord, cols) {
  d <- meta_df %>%
    filter(orig.ident %in% s_ord) %>%
    group_by(orig.ident) %>%
    summarise(!!sym(var) := first(!!sym(var)), .groups = "drop") %>%
    mutate(
      orig.ident = factor(orig.ident, levels = s_ord),
      ann_label = factor(label, levels = label)
    )
  
  ggplot(d, aes(x = orig.ident, y = ann_label, fill = !!sym(var))) +
    geom_tile(
      width = 0.96,     # slight separation between samples
      height = 1,
      colour = NA
    ) +
    scale_fill_manual(
      values = cols,
      na.value = "grey95",
      drop = FALSE,
      name = label
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    anno_base_theme +
    right_legend_theme
}

####################
# 10) Main plotting function
####################
plot_final <- function(lbl_vec, s_ord, col_map, title, legend_title, lbl_ord = NULL) {
  dat <- data.frame(
    cell = names(lbl_vec),
    label = as.character(lbl_vec),
    stringsAsFactors = FALSE
  ) %>%
    left_join(meta_df %>% distinct(cell, orig.ident), by = "cell")
  
  if (is.null(lbl_ord)) {
    lbl_ord <- names(col_map)
  }
  
  props <- dat %>%
    count(orig.ident, label, name = "n") %>%
    group_by(orig.ident) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    filter(orig.ident %in% s_ord) %>%
    mutate(
      orig.ident = factor(orig.ident, levels = s_ord),
      label = factor(label, levels = rev(lbl_ord))
    )
  
  tot <- meta_df %>%
    count(orig.ident, name = "n") %>%
    filter(orig.ident %in% s_ord) %>%
    mutate(orig.ident = factor(orig.ident, levels = s_ord))
  
  sf <- max(tot$n, na.rm = TRUE) / 100
  
  p_main <- ggplot(props, aes(x = orig.ident, y = pct, fill = label)) +
    geom_col(width = 0.90, colour = NA) +
    geom_point(
      data = tot,
      aes(x = orig.ident, y = n / sf),
      inherit.aes = FALSE,
      colour = "black",
      size = 2.4,
      shape = 16
    ) +
    geom_line(
      data = tot,
      aes(x = orig.ident, y = n / sf, group = 1),
      inherit.aes = FALSE,
      colour = "black",
      alpha = 0.45,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    scale_fill_manual(
      values = col_map,
      breaks = lbl_ord,
      drop = FALSE,
      name = legend_title
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(
      name = "Proportion (%)",
      limits = c(0, 100),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * sf, name = "Total Cells")
    ) +
    labs(x = NULL) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 17,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(size = 17, color = "black"),
      axis.title.y = element_text(size = 19, face = "bold"),
      axis.title.y.right = element_text(size = 19, face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(12, 4, 2, 4)
    ) +
    guides(
      fill = guide_legend(
        ncol = 1,   # 1 column for sidebar
        byrow = TRUE
      )
    ) +
    bottom_state_legend_theme
  
  a1 <- make_anno("Batch",       "Batch",       s_ord, cl_cols$Batch)
  a2 <- make_anno("Gender",      "Gender",      s_ord, cl_cols$Gender)
  a3 <- make_anno("Age_Group",   "Age Group",   s_ord, cl_cols$Age_Group)
  a4 <- make_anno("Response",    "Response",    s_ord, cl_cols$Response)
  a5 <- make_anno("Tumour_Type", "Tumour Type", s_ord, cl_cols$Tumour_Type)
  a6 <- make_anno("Histology",   "Histology",   s_ord, cl_cols$Histology)
  a7 <- make_anno("T_Stage",     "T Stage",     s_ord, cl_cols$T_Stage)
  a8 <- make_anno("Timepoint",   "Timepoint",   s_ord, cl_cols$Timepoint)
  spacer_between <- patchwork::plot_spacer()
  inner_p <- (a1 / a2 / a3 / a4 / a5 / a6 / a7 / a8 / spacer_between / p_main) +
    plot_layout(heights = c(rep(0.35, 8), 0.2, 7))

  final_plot <- (inner_p | patchwork::guide_area()) +
    plot_layout(
      widths = c(4.2, 1), # 4.2:1 ratio for plot vs sidebar
      guides = "collect"
    ) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 40, face = "bold", hjust = 0.5), # larger title
        plot.margin = margin(15, 40, 15, 15) # extra right space
      )
    )
  
  final_plot
}

####################
# 11) Generate output
####################
message("Writing PDF...")

pdf(
  file.path(out_dir, "Auto_sample_abundance_pdo.pdf"),
  width = 34,
  height = 24,
  onefile = TRUE
)

views <- list(
  "Top nCC MP" = list(
    l = topmp_noncc,
    o = ord1,
    c = c1,
    legend_title = "Top nCC MP"
  ),
  "Top All MP" = list(
    l = topmp_all,
    o = ord2,
    c = c2,
    legend_title = "Top All MP"
  ),
  "Cell States" = list(
    l = states,
    o = ord3,
    c = group_cols,
    legend_title = "Cell State"
  )
)

for (v in names(views)) {
  for (s in names(sort_list)) {
    p <- plot_final(
      lbl_vec = views[[v]]$l,
      s_ord = sort_list[[s]],
      col_map = views[[v]]$c,
      title = paste(v, "|", s),
      legend_title = views[[v]]$legend_title,
      lbl_ord = views[[v]]$o
    )
    print(p)
  }
}

dev.off()

message("Success! PDF saved to: sample_abundance/Auto_sample_abundance_pdo.pdf")