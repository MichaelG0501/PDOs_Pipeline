####################
# Auto_pdo_sn_matched_pair_comparison.R
# Compare matched PDO organoid samples against matched snRNA-seq biopsies
# using finalized cell-state proportions and top-metaprogram proportions.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

####################
# setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_pdo_sn_matched_pair_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sn_base_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs"
sn_state_path <- file.path(sn_base_dir, "Auto_final_states.rds")
if (!file.exists(sn_state_path)) {
  sn_state_path <- file.path(sn_base_dir, "Auto_topmp_v2_noreg_states_B.rds")
}

####################
# path helpers
####################
resolve_first_existing <- function(candidates, desc) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(desc, " not found in expected locations: ", paste(candidates, collapse = "; "))
  }
  existing[1]
}

sn_base_dir <- resolve_first_existing(
  candidates = c(
    "/rds/general/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs"
  ),
  desc = "snSeq sn_outs directory"
)

sn_state_path <- resolve_first_existing(
  candidates = c(
    file.path(sn_base_dir, "Auto_final_states.rds"),
    file.path(sn_base_dir, "Auto_topmp_v2_noreg_states_B.rds")
  ),
  desc = "snRNA-seq state RDS"
)

pdo_merged_path <- resolve_first_existing(
  candidates = c(
    "PDOs_merged.rds",
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/PDOs_merged.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/PDOs_merged.rds",
    "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDOs_merged.rds"
  ),
  desc = "PDOs_merged.rds"
)

pdo_states_path <- resolve_first_existing(
  candidates = c(
    "Auto_PDO_final_states.rds",
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_PDO_final_states.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_PDO_final_states.rds"
  ),
  desc = "Auto_PDO_final_states.rds"
)

pdo_ucell_path <- resolve_first_existing(
  candidates = c(
    "UCell_scores_filtered.rds",
    "Metaprogrammes_Results/UCell_scores_filtered.rds",
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/UCell_scores_filtered.rds",
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/UCell_scores_filtered.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/UCell_scores_filtered.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/UCell_scores_filtered.rds"
  ),
  desc = "UCell_scores_filtered.rds"
)

pdo_gnmf_path <- resolve_first_existing(
  candidates = c(
    "Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds",
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds"
  ),
  desc = "PDO geneNMF_metaprograms_nMP_13.rds"
)

sn_malignant_path <- resolve_first_existing(
  candidates = c(
    file.path(sn_base_dir, "snSeq_malignant_epi.rds")
  ),
  desc = "snSeq_malignant_epi.rds"
)

sn_ucell_path <- resolve_first_existing(
  candidates = c(
    file.path(sn_base_dir, "Metaprogrammes_Results/UCell_nMP19_filtered.rds")
  ),
  desc = "sn UCell_nMP19_filtered.rds"
)

sn_gnmf_path <- resolve_first_existing(
  candidates = c(
    file.path(sn_base_dir, "Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
  ),
  desc = "sn geneNMF_metaprograms_nMP_19.rds"
)

####################
# constants
####################
pair_map <- data.frame(
  patient_id = c("SUR680", "SUR791"),
  pdo_sample = c("SUR680T3_PDO", "SUR791T3_PDO"),
  sn_sample = c("H_post_T1_biopsy", "L_post_T1_biopsy"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    pair_id = paste0(patient_id, "_matched_pair"),
    pair_title = paste0(patient_id, ": PDO organoid vs matched snRNA-seq biopsy")
  )

state_levels <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation",
  "Unresolved",
  "Hybrid"
)

fill_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#1B9E77",
  "Unresolved" = "grey80",
  "Hybrid" = "black",
  "Cell cycle" = "grey45",
  "Other" = "grey70"
)

pdo_mp_desc <- c(
  "MP6" = "G2M Cell Cycle",
  "MP7" = "DNA repair",
  "MP5" = "MYC-related Proliferation",
  "MP1" = "G2M checkpoint",
  "MP3" = "G1S Cell Cycle",
  "MP8" = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9" = "ECM Remodeling Epi.",
  "MP4" = "Intestinal Metaplasia"
)

pdo_mp_group <- c(
  "MP6" = "Cell cycle",
  "MP7" = "Cell cycle",
  "MP5" = "Classic Proliferative",
  "MP1" = "Cell cycle",
  "MP3" = "Cell cycle",
  "MP8" = "SMG-like Metaplasia",
  "MP10" = "Stress-adaptive",
  "MP9" = "Stress-adaptive",
  "MP4" = "Basal to Intestinal Metaplasia"
)

sn_mp_desc <- c(
  "MP1" = "G2M Cell Cycle",
  "MP9" = "G1S Cell Cycle",
  "MP2" = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5" = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8" = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7" = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

sn_mp_group <- c(
  "MP1" = "Cell cycle",
  "MP9" = "Cell cycle",
  "MP2" = "Classic Proliferative",
  "MP17" = "Basal to Intestinal Metaplasia",
  "MP14" = "Basal to Intestinal Metaplasia",
  "MP5" = "Basal to Intestinal Metaplasia",
  "MP10" = "Basal to Intestinal Metaplasia",
  "MP8" = "Basal to Intestinal Metaplasia",
  "MP13" = "Stress-adaptive",
  "MP7" = "Cell cycle",
  "MP18" = "SMG-like Metaplasia",
  "MP16" = "SMG-like Metaplasia",
  "MP15" = "Immune Infiltrating",
  "MP12" = "Stress-adaptive"
)

####################
# helpers
####################
harmonise_state_names <- function(x) {
  x <- as.character(x)
  x[x == "Basal to Intest. Meta"] <- "Basal to Intestinal Metaplasia"
  x
}

label_mp <- function(mp_ids, desc_map) {
  desc <- unname(desc_map[mp_ids])
  desc[is.na(desc)] <- mp_ids[is.na(desc)]
  paste0(mp_ids, " ", desc)
}

get_pdo_retained_mps <- function(geneNMF.metaprograms) {
  mp_genes <- geneNMF.metaprograms$metaprograms.genes
  coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
  bad_mps <- paste0("MP", which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0))
  setdiff(names(mp_genes), c(bad_mps, low_coverage_mps))
}

get_sn_retained_mps <- function(geneNMF.metaprograms) {
  mp_genes <- geneNMF.metaprograms$metaprograms.genes
  bad_mps <- paste0("MP", which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0))
  setdiff(names(mp_genes), bad_mps)
}

####################
# z-normalisation (must match state assignment scripts exactly)
# sample-mean center, then batch-SD scale
####################
z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

####################
# keep only state-defining MPs (non-cell-cycle)
####################
get_state_defining_mps <- function(retained_mps, group_map) {
  mp_group <- unname(group_map[retained_mps])
  kept <- retained_mps[!is.na(mp_group) & mp_group != "Cell cycle"]
  if (length(kept) == 0) {
    stop("No non-cell-cycle retained MPs are available for top-MP scoring.")
  }
  kept
}

check_samples_present <- function(all_samples, required_samples, dataset_name) {
  missing_samples <- setdiff(required_samples, all_samples)
  if (length(missing_samples) > 0) {
    stop(
      dataset_name,
      " is missing required samples: ",
      paste(missing_samples, collapse = ", ")
    )
  }
}

assign_top_mp <- function(ucell_scores, retained_mps, desc_map, group_map) {
  keep_mps <- intersect(retained_mps, colnames(ucell_scores))
  if (length(keep_mps) == 0) {
    stop("No retained MPs found in the supplied UCell matrix.")
  }

  top_mp_id <- colnames(ucell_scores[, keep_mps, drop = FALSE])[
    max.col(ucell_scores[, keep_mps, drop = FALSE], ties.method = "first")
  ]

  data.frame(
    cell = rownames(ucell_scores),
    mp_id = top_mp_id,
    label = label_mp(top_mp_id, desc_map),
    mp_group = unname(group_map[top_mp_id]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      mp_group = ifelse(is.na(mp_group), "Other", mp_group)
    )
}

make_prop_table <- function(sample_by_cell, label_df, sample_meta, label_levels = NULL) {
  prop_df <- label_df %>%
    mutate(
      sample = unname(sample_by_cell[cell])
    ) %>%
    filter(!is.na(sample)) %>%
    inner_join(sample_meta, by = "sample")

  if (is.null(label_levels)) {
    label_levels <- sort(unique(prop_df$label))
  }

  prop_df %>%
    count(pair_id, pair_title, patient_id, modality, sample, label, .drop = FALSE, name = "n") %>%
    complete(
      nesting(pair_id, pair_title, patient_id, modality, sample),
      label = label_levels,
      fill = list(n = 0L)
    ) %>%
    group_by(pair_id, pair_title, patient_id, modality, sample) %>%
    mutate(
      total_n = sum(n),
      pct = ifelse(total_n > 0, 100 * n / total_n, 0)
    ) %>%
    ungroup()
}

####################
# recompute state percentages on selected labels only
####################
focus_state_prop <- function(state_prop_df, keep_labels) {
  state_prop_df %>%
    filter(label %in% keep_labels) %>%
    group_by(pair_id, pair_title, patient_id, modality, sample) %>%
    mutate(
      total_focus_n = sum(n),
      pct = ifelse(total_focus_n > 0, 100 * n / total_focus_n, 0),
      total_n = total_focus_n
    ) %>%
    ungroup() %>%
    select(-total_focus_n)
}

####################
# derive state-like proportions from top MP groups
####################
state_from_mp_prop <- function(mp_prop_df, keep_labels) {
  mp_prop_df %>%
    mutate(label = mp_group) %>%
    filter(label %in% keep_labels) %>%
    select(pair_id, pair_title, patient_id, modality, sample, label, n) %>%
    group_by(pair_id, pair_title, patient_id, modality, sample, label) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(pair_id, pair_title, patient_id, modality, sample) %>%
    mutate(
      total_n = sum(n),
      pct = ifelse(total_n > 0, 100 * n / total_n, 0)
    ) %>%
    ungroup() %>%
    arrange(patient_id, modality, label)
}

make_state_plot <- function(
  state_df,
  pair_info,
  state_order = state_levels,
  plot_title = "Finalized Cell States"
) {
  legend_rows <- ifelse(length(state_order) > 6, 3, 2)

  plot_df <- state_df %>%
    mutate(
      modality_label = factor(
        paste0(modality, "\n", sample, "\nN=", comma(total_n)),
        levels = paste0(
          c("PDO organoid", "snRNA-seq biopsy"),
          "\n",
          c(pair_info$pdo_sample, pair_info$sn_sample),
          "\nN=",
          comma(c(
            state_df$total_n[state_df$modality == "PDO organoid"][1],
            state_df$total_n[state_df$modality == "snRNA-seq biopsy"][1]
          ))
        )
      ),
      label = factor(label, levels = rev(state_order))
    )

  ggplot(plot_df, aes(x = modality_label, y = pct, fill = label)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(
      data = subset(plot_df, pct >= 7),
      aes(label = sprintf("%.1f%%", pct)),
      position = position_stack(vjust = 0.5),
      size = 3.2,
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(values = fill_cols, breaks = state_order, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = c(0, 0),
      labels = label_number(accuracy = 1)
    ) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Proportion (%)",
      fill = NULL
    ) +
    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 9, face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.justification = "center",
      legend.box.just = "center",
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.x = grid::unit(0.2, "cm"),
      legend.box.margin = margin(0, 12, 0, 12),
      legend.text = element_text(size = 8.3),
      plot.margin = margin(6, 16, 14, 16)
    )
}

make_mp_plot <- function(mp_df, plot_title) {
  plot_df <- mp_df %>%
    filter(pct > 0) %>%
    arrange(pct) %>%
    mutate(
      label = factor(label, levels = label)
    )

  ggplot(plot_df, aes(x = pct, y = label, fill = mp_group)) +
    geom_col(width = 0.75, colour = NA) +
    geom_text(
      aes(label = sprintf("%.1f%%", pct)),
      hjust = -0.1,
      size = 3.1
    ) +
    scale_fill_manual(values = fill_cols, breaks = names(fill_cols), drop = FALSE) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.12)),
      labels = label_number(accuracy = 1)
    ) +
    labs(
      title = plot_title,
      x = "Proportion (%)",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    ) +
    guides(fill = "none")
}

make_pair_summary <- function(state_props, mp_props, pair_info) {
  pdo_state_vec <- state_props %>%
    filter(modality == "PDO organoid") %>%
    group_by(label) %>%
    summarise(pdo_pct = sum(pct), .groups = "drop")

  sn_state_vec <- state_props %>%
    filter(modality == "snRNA-seq biopsy") %>%
    group_by(label) %>%
    summarise(sn_pct = sum(pct), .groups = "drop")

  state_wide <- data.frame(
    label = state_levels,
    stringsAsFactors = FALSE
  ) %>%
    left_join(pdo_state_vec, by = "label") %>%
    left_join(sn_state_vec, by = "label") %>%
    mutate(
      pdo_pct = ifelse(is.na(pdo_pct), 0, pdo_pct),
      sn_pct = ifelse(is.na(sn_pct), 0, sn_pct)
    )

  pdo_state <- subset(state_props, modality == "PDO organoid")
  sn_state <- subset(state_props, modality == "snRNA-seq biopsy")
  pdo_mp <- subset(mp_props, modality == "PDO organoid")
  sn_mp <- subset(mp_props, modality == "snRNA-seq biopsy")

  data.frame(
    pair_id = pair_info$pair_id,
    patient_id = pair_info$patient_id,
    pdo_sample = pair_info$pdo_sample,
    sn_sample = pair_info$sn_sample,
    pdo_cells = unique(pdo_state$total_n),
    sn_cells = unique(sn_state$total_n),
    dominant_pdo_state = pdo_state$label[which.max(pdo_state$pct)],
    dominant_sn_state = sn_state$label[which.max(sn_state$pct)],
    state_overlap_pct = sum(pmin(state_wide$pdo_pct, state_wide$sn_pct)),
    state_pearson = suppressWarnings(cor(
      state_wide$pdo_pct,
      state_wide$sn_pct,
      method = "pearson"
    )),
    state_spearman = suppressWarnings(cor(
      state_wide$pdo_pct,
      state_wide$sn_pct,
      method = "spearman"
    )),
    dominant_pdo_mp = pdo_mp$label[which.max(pdo_mp$pct)],
    dominant_sn_mp = sn_mp$label[which.max(sn_mp$pct)],
    stringsAsFactors = FALSE
  )
}

####################
# load data
####################
message("Loading PDO inputs ...")
pdos <- readRDS(pdo_merged_path)
pdo_states <- readRDS(pdo_states_path)
pdo_ucell <- readRDS(pdo_ucell_path)
pdo_gnmf <- readRDS(pdo_gnmf_path)

message("Loading snRNA-seq inputs ...")
sn <- readRDS(sn_malignant_path)
sn_states <- readRDS(sn_state_path)
sn_ucell <- readRDS(sn_ucell_path)
sn_gnmf <- readRDS(sn_gnmf_path)

####################
# validate samples
####################
check_samples_present(
  all_samples = unique(as.character(pdos$orig.ident)),
  required_samples = pair_map$pdo_sample,
  dataset_name = "PDO object"
)
check_samples_present(
  all_samples = unique(as.character(sn$orig.ident)),
  required_samples = pair_map$sn_sample,
  dataset_name = "snRNA-seq object"
)

####################
# harmonise matched metadata
####################
sample_meta <- bind_rows(
  pair_map %>%
    transmute(
      pair_id,
      pair_title,
      patient_id,
      modality = "PDO organoid",
      sample = pdo_sample
    ),
  pair_map %>%
    transmute(
      pair_id,
      pair_title,
      patient_id,
      modality = "snRNA-seq biopsy",
      sample = sn_sample
    )
)

pdo_cells_keep <- Cells(pdos)[pdos$orig.ident %in% pair_map$pdo_sample]
sn_cells_keep <- Cells(sn)[sn$orig.ident %in% pair_map$sn_sample]

pdo_state_df <- data.frame(
  cell = intersect(names(pdo_states), pdo_cells_keep),
  label = harmonise_state_names(pdo_states[intersect(names(pdo_states), pdo_cells_keep)]),
  stringsAsFactors = FALSE
)

sn_state_df <- data.frame(
  cell = intersect(names(sn_states), sn_cells_keep),
  label = harmonise_state_names(sn_states[intersect(names(sn_states), sn_cells_keep)]),
  stringsAsFactors = FALSE
)

pdo_sample_by_cell <- setNames(as.character(pdos$orig.ident), Cells(pdos))
sn_sample_by_cell <- setNames(as.character(sn$orig.ident), Cells(sn))

pdo_retained_mps <- get_pdo_retained_mps(pdo_gnmf)
sn_retained_mps <- get_sn_retained_mps(sn_gnmf)

pdo_state_defining_mps <- get_state_defining_mps(
  retained_mps = pdo_retained_mps,
  group_map = pdo_mp_group
)
sn_state_defining_mps <- get_state_defining_mps(
  retained_mps = sn_retained_mps,
  group_map = sn_mp_group
)

# Identify cells that are NOT Unresolved or Hybrid for focused MP analysis
pdo_cells_focus <- pdo_state_df$cell[!(pdo_state_df$label %in% c("Unresolved", "Hybrid"))]
sn_cells_focus <- sn_state_df$cell[!(sn_state_df$label %in% c("Unresolved", "Hybrid"))]

####################
# Z-normalise UCell scores BEFORE top-MP assignment
####################
message("Z-normalising PDO UCell scores for top-MP assignment ...")
pdo_common_cells <- intersect(rownames(pdo_ucell), Cells(pdos))
pdo_sample_var_full <- setNames(as.character(pdos$orig.ident), Cells(pdos))
pdo_batch_var_full <- setNames(
  ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch"),
  Cells(pdos)
)
pdo_ucell_znorm <- z_normalise(
  mat = as.matrix(pdo_ucell[pdo_common_cells, , drop = FALSE]),
  sample_var = pdo_sample_var_full,
  study_var = pdo_batch_var_full
)

message("Z-normalising snRNA-seq UCell scores for top-MP assignment ...")
sn_common_cells <- intersect(rownames(sn_ucell), Cells(sn))
sn_sample_var_full <- setNames(as.character(sn$orig.ident), Cells(sn))
sn_study_var_full <- setNames(as.character(sn$reference_batch), Cells(sn))
sn_ucell_znorm <- z_normalise(
  mat = as.matrix(sn_ucell[sn_common_cells, , drop = FALSE]),
  sample_var = sn_sample_var_full,
  study_var = sn_study_var_full
)

pdo_mp_df <- assign_top_mp(
  ucell_scores = pdo_ucell_znorm[intersect(rownames(pdo_ucell_znorm), pdo_cells_keep), , drop = FALSE],
  retained_mps = pdo_retained_mps,
  desc_map = pdo_mp_desc,
  group_map = pdo_mp_group
)

pdo_mp_state_df <- assign_top_mp(
  ucell_scores = pdo_ucell_znorm[intersect(rownames(pdo_ucell_znorm), pdo_cells_focus), , drop = FALSE],
  retained_mps = pdo_state_defining_mps,
  desc_map = pdo_mp_desc,
  group_map = pdo_mp_group
)

sn_mp_df <- assign_top_mp(
  ucell_scores = sn_ucell_znorm[intersect(rownames(sn_ucell_znorm), sn_cells_keep), , drop = FALSE],
  retained_mps = sn_retained_mps,
  desc_map = sn_mp_desc,
  group_map = sn_mp_group
)

sn_mp_state_df <- assign_top_mp(
  ucell_scores = sn_ucell_znorm[intersect(rownames(sn_ucell_znorm), sn_cells_focus), , drop = FALSE],
  retained_mps = sn_state_defining_mps,
  desc_map = sn_mp_desc,
  group_map = sn_mp_group
)

####################
# proportions
####################
state_prop <- bind_rows(
  make_prop_table(
    sample_by_cell = pdo_sample_by_cell,
    label_df = pdo_state_df,
    sample_meta = subset(sample_meta, modality == "PDO organoid"),
    label_levels = state_levels
  ),
  make_prop_table(
    sample_by_cell = sn_sample_by_cell,
    label_df = sn_state_df,
    sample_meta = subset(sample_meta, modality == "snRNA-seq biopsy"),
    label_levels = state_levels
  )
)

pdo_mp_levels <- pdo_mp_df %>%
  distinct(label, mp_group, mp_id) %>%
  arrange(match(mp_group, names(fill_cols)), mp_id) %>%
  pull(label)

sn_mp_levels <- sn_mp_df %>%
  distinct(label, mp_group, mp_id) %>%
  arrange(match(mp_group, names(fill_cols)), mp_id) %>%
  pull(label)

pdo_mp_prop <- make_prop_table(
  sample_by_cell = pdo_sample_by_cell,
  label_df = pdo_mp_df[, c("cell", "label")],
  sample_meta = subset(sample_meta, modality == "PDO organoid"),
  label_levels = pdo_mp_levels
) %>%
  left_join(
    distinct(pdo_mp_df[, c("label", "mp_id", "mp_group")]),
    by = "label"
  )

sn_mp_prop <- make_prop_table(
  sample_by_cell = sn_sample_by_cell,
  label_df = sn_mp_df[, c("cell", "label")],
  sample_meta = subset(sample_meta, modality == "snRNA-seq biopsy"),
  label_levels = sn_mp_levels
) %>%
  left_join(
    distinct(sn_mp_df[, c("label", "mp_id", "mp_group")]),
    by = "label"
  )

mp_prop <- bind_rows(pdo_mp_prop, sn_mp_prop)

####################
# focused first-page inputs
####################
state_focus_levels <- setdiff(state_levels, c("Unresolved", "Hybrid"))
state_prop_focus <- focus_state_prop(
  state_prop_df = state_prop,
  keep_labels = state_focus_levels
)

state_focus_from_mp_levels <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia"
)

pdo_mp_state_levels <- pdo_mp_state_df %>%
  distinct(label, mp_group, mp_id) %>%
  arrange(match(mp_group, names(fill_cols)), mp_id) %>%
  pull(label)

sn_mp_state_levels <- sn_mp_state_df %>%
  distinct(label, mp_group, mp_id) %>%
  arrange(match(mp_group, names(fill_cols)), mp_id) %>%
  pull(label)

pdo_mp_state_prop <- make_prop_table(
  sample_by_cell = pdo_sample_by_cell,
  label_df = pdo_mp_state_df[, c("cell", "label")],
  sample_meta = subset(sample_meta, modality == "PDO organoid"),
  label_levels = pdo_mp_state_levels
) %>%
  left_join(
    distinct(pdo_mp_state_df[, c("label", "mp_id", "mp_group")]),
    by = "label"
  )

sn_mp_state_prop <- make_prop_table(
  sample_by_cell = sn_sample_by_cell,
  label_df = sn_mp_state_df[, c("cell", "label")],
  sample_meta = subset(sample_meta, modality == "snRNA-seq biopsy"),
  label_levels = sn_mp_state_levels
) %>%
  left_join(
    distinct(sn_mp_state_df[, c("label", "mp_id", "mp_group")]),
    by = "label"
  )

mp_prop_focus <- bind_rows(pdo_mp_state_prop, sn_mp_state_prop)

state_prop_focus_from_mp <- state_from_mp_prop(
  mp_prop_df = mp_prop_focus,
  keep_labels = state_focus_from_mp_levels
)

####################
# write summary tables
####################
state_csv <- state_prop %>%
  arrange(patient_id, modality, factor(label, levels = state_levels))
write.csv(
  state_csv,
  file.path(out_dir, "Auto_pdo_sn_matched_pair_state_proportions.csv"),
  row.names = FALSE
)

mp_csv <- mp_prop %>%
  arrange(patient_id, modality, desc(pct), label)
write.csv(
  mp_csv,
  file.path(out_dir, "Auto_pdo_sn_matched_pair_mp_proportions.csv"),
  row.names = FALSE
)

state_from_topmp_csv <- state_prop_focus_from_mp %>%
  arrange(patient_id, modality, factor(label, levels = state_focus_from_mp_levels))
write.csv(
  state_from_topmp_csv,
  file.path(out_dir, "Auto_pdo_sn_matched_pair_state_from_topmp_groups.csv"),
  row.names = FALSE
)

pair_summaries <- lapply(seq_len(nrow(pair_map)), function(i) {
  pair_info <- pair_map[i, ]
  make_pair_summary(
    state_props = subset(state_prop, pair_id == pair_info$pair_id),
    mp_props = subset(mp_prop, pair_id == pair_info$pair_id),
    pair_info = pair_info
  )
})
pair_summary_df <- bind_rows(pair_summaries) %>%
  mutate(
    dominant_state_match = dominant_pdo_state == dominant_sn_state
  )

write.csv(
  pair_summary_df,
  file.path(out_dir, "Auto_pdo_sn_matched_pair_summary.csv"),
  row.names = FALSE
)

####################
# build plots
####################
message("Building matched-pair comparison plots ...")

####################
# page builder
####################
build_pair_page <- function(
  state_prop_input,
  mp_prop_input,
  state_order,
  state_title,
  mp_title_prefix,
  page_title
) {
  pair_plots <- lapply(seq_len(nrow(pair_map)), function(i) {
    pair_info <- pair_map[i, ]
    pair_state <- subset(state_prop_input, pair_id == pair_info$pair_id)
    pair_mp <- subset(mp_prop_input, pair_id == pair_info$pair_id)

    p_state <- make_state_plot(
      state_df = pair_state,
      pair_info = pair_info,
      state_order = state_order,
      plot_title = state_title
    )

    p_pdo_mp <- make_mp_plot(
      mp_df = subset(pair_mp, modality == "PDO organoid"),
      plot_title = paste0(mp_title_prefix, "\n", pair_info$pdo_sample)
    )

    p_sn_mp <- make_mp_plot(
      mp_df = subset(pair_mp, modality == "snRNA-seq biopsy"),
      plot_title = paste0(mp_title_prefix, "\n", pair_info$sn_sample)
    )

    (p_state | p_pdo_mp | p_sn_mp) +
      plot_layout(widths = c(1.25, 1.35, 1.35)) +
      plot_annotation(title = pair_info$pair_title)
  })

  wrap_plots(pair_plots, ncol = 1, guides = "collect") +
    plot_annotation(title = page_title)
}

####################
# page 1: state-defining focus
####################
overview_plot_page1 <- build_pair_page(
  state_prop_input = state_prop_focus,
  mp_prop_input = mp_prop_focus,
  state_order = state_focus_levels,
  state_title = "Finalized States (Filtered: No Unresolved/Hybrid)",
  mp_title_prefix = "Top State-defining MP (Filtered Cells)",
  page_title = "Matched PDO vs snRNA-seq Pair Comparison - Page 1"
)

####################
# page 2: original full view (no subtitles)
####################
overview_plot_page2 <- build_pair_page(
  state_prop_input = state_prop,
  mp_prop_input = mp_prop,
  state_order = state_levels,
  state_title = "Finalized Cell States",
  mp_title_prefix = "Top MP Proportions",
  page_title = "Matched PDO vs snRNA-seq Pair Comparison - Page 2"
)

####################
# write two-page pdf (no png)
####################
pdf_path <- file.path(out_dir, "Auto_pdo_sn_matched_pair_comparison.pdf")
pdf(pdf_path, width = 18, height = 14, onefile = TRUE)
print(overview_plot_page1)
print(overview_plot_page2)
dev.off()

stale_pngs <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
if (length(stale_pngs) > 0) {
  file.remove(stale_pngs)
}

message("Finished. Outputs written to: ", file.path(getwd(), out_dir))
