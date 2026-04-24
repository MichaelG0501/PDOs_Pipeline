####################
# Auto_pdo_flot_matched_response.R
# Paired FLOT response analysis for the four matched treated/untreated PDO pairs.
# Produces presentation-ready abundance, pathway-response, and pseudobulk DEG
# summaries in PDOs_outs/Auto_pdo_flot_matched_response/.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(scales)
  library(Matrix)
  library(msigdbr)
  library(edgeR)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
})

####################
# setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_pdo_flot_matched_response"
deg_dir <- file.path(out_dir, "pseudobulk_deg")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

matched_samples <- c(
  "SUR1070_Treated_PDO", "SUR1070_Untreated_PDO",
  "SUR1090_Treated_PDO", "SUR1090_Untreated_PDO",
  "SUR1072_Treated_PDO", "SUR1072_Untreated_PDO",
  "SUR1181_Treated_PDO", "SUR1181_Untreated_PDO"
)

patient_order <- c("SUR1070", "SUR1090", "SUR1072", "SUR1181")
patient_labels <- c(
  SUR1070 = "SUR1070 (Pre-responder)",
  SUR1090 = "SUR1090 (Post-responder)",
  SUR1072 = "SUR1072 (Post-nonresponder)",
  SUR1181 = "SUR1181 (Pre-nonresponder)"
)
patient_short <- c(
  SUR1070 = "1070",
  SUR1090 = "1090",
  SUR1072 = "1072",
  SUR1181 = "1181"
)

state_levels_main <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "3CA_EMT_and_Protein_maturation"
)

state_levels_all <- c(
  state_levels_main,
  "Unresolved",
  "Hybrid"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

patient_cols <- c(
  "SUR1070 (Pre-responder)" = "#4C78A8",
  "SUR1090 (Post-responder)" = "#B07AA1",
  "SUR1072 (Post-nonresponder)" = "#59A14F",
  "SUR1181 (Pre-nonresponder)" = "#F28E2B"
)

treatment_cols <- c(
  "Untreated" = "#D58B2D",
  "Treated" = "#374151"
)

pathway_ids <- c(
  "HALLMARK_E2F_TARGETS" = "E2F targets",
  "HALLMARK_G2M_CHECKPOINT" = "G2M checkpoint",
  "HALLMARK_APOPTOSIS" = "Apoptosis",
  "HALLMARK_P53_PATHWAY" = "p53 pathway",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFa/NFkB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_HYPOXIA" = "Hypoxia",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = "Unfolded protein response",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation",
  "HALLMARK_DNA_REPAIR" = "DNA repair",
  "HALLMARK_XENOBIOTIC_METABOLISM" = "Xenobiotic metabolism"
)

pathway_order <- c(
  unname(pathway_ids),
  "Interferon response"
)

composite_pathways <- list(
  "Cell-cycle / proliferation" = c("E2F targets", "G2M checkpoint"),
  "Injury / checkpoint" = c("Apoptosis", "p53 pathway", "DNA repair"),
  "Adaptive / persistence" = c("TNFa/NFkB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response"),
  "Proteostasis / transition" = c("Unfolded protein response", "Oxidative phosphorylation")
)

min_cells_per_pseudobulk <- 20
min_pairs_for_deg <- 3
top_deg_per_direction <- 8

####################
# helpers
####################
flot_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

state_key <- function(x) {
  gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

mean_score <- function(mat, genes) {
  genes_use <- intersect(unique(genes), rownames(mat))
  if (length(genes_use) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  Matrix::colMeans(mat[genes_use, , drop = FALSE], na.rm = TRUE)
}

zscore_rows <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z
}

build_discrete_umap <- function(plot_df, value_col, color_map, title_text) {
  ggplot(plot_df, aes(UMAP_1, UMAP_2, color = .data[[value_col]])) +
    geom_point(size = 0.22, alpha = 0.82) +
    scale_color_manual(values = color_map, drop = FALSE) +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    coord_equal() +
    flot_theme(10) +
    theme(legend.position = "right")
}

build_continuous_umap <- function(plot_df, value_col, title_text) {
  ggplot(plot_df, aes(UMAP_1, UMAP_2, color = .data[[value_col]])) +
    geom_point(size = 0.22, alpha = 0.82) +
    scale_color_gradientn(
      colors = c("#1B365D", "#4F9DD9", "#F2E6B6", "#B33C2F"),
      na.value = "grey85"
    ) +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    coord_equal() +
    flot_theme(10) +
    theme(legend.position = "right")
}

score_gene_sets <- function(expr_mat, gene_sets) {
  score_mat <- sapply(gene_sets, function(genes) mean_score(expr_mat, genes))
  score_mat <- t(score_mat)
  rownames(score_mat) <- names(gene_sets)
  score_mat
}

classify_state_fate <- function(abundance_log2fc, proliferation_loss, injury, adaptation, proteostasis) {
  if (!is.na(abundance_log2fc) && abundance_log2fc <= -0.35 &&
      (( !is.na(injury) && injury >= 0.15) || (!is.na(proliferation_loss) && proliferation_loss >= 0.15))) {
    return("Likely sensitive")
  }
  if (!is.na(abundance_log2fc) && abundance_log2fc >= 0.10 &&
      !is.na(adaptation) && adaptation >= 0.15) {
    return("Adaptive / persistent")
  }
  if (!is.na(proteostasis) && proteostasis >= 0.15 &&
      !is.na(adaptation) && adaptation >= 0.05) {
    return("Transition-like")
  }
  if (!is.na(abundance_log2fc) && abs(abundance_log2fc) < 0.15 &&
      !is.na(injury) && abs(injury) < 0.10 &&
      !is.na(adaptation) && abs(adaptation) < 0.10 &&
      !is.na(proteostasis) && abs(proteostasis) < 0.10) {
    return("Buffered / tolerant")
  }
  "Mixed / context-dependent"
}

select_recurrent_genes <- function(res_tbl, patient_logfc, top_per_direction = 8) {
  common_genes <- intersect(res_tbl$gene, rownames(patient_logfc))
  if (length(common_genes) == 0) {
    return(character(0))
  }

  patient_logfc <- patient_logfc[common_genes, , drop = FALSE]
  res_sub <- res_tbl %>%
    filter(gene %in% common_genes) %>%
    mutate(
      consistency_n = rowSums(sign(patient_logfc[gene, , drop = FALSE]) == sign(logFC), na.rm = TRUE),
      patient_median_logFC = apply(patient_logfc[gene, , drop = FALSE], 1, median, na.rm = TRUE)
    )

  min_consistency <- max(2, ceiling(0.75 * ncol(patient_logfc)))

  up_tbl <- res_sub %>%
    filter(logFC > 0, consistency_n >= min_consistency) %>%
    arrange(FDR, desc(abs(logFC)), desc(consistency_n))
  down_tbl <- res_sub %>%
    filter(logFC < 0, consistency_n >= min_consistency) %>%
    arrange(FDR, desc(abs(logFC)), desc(consistency_n))

  if (nrow(up_tbl) < top_per_direction) {
    up_tbl <- res_sub %>%
      filter(logFC > 0) %>%
      arrange(FDR, desc(abs(logFC)), desc(consistency_n))
  }
  if (nrow(down_tbl) < top_per_direction) {
    down_tbl <- res_sub %>%
      filter(logFC < 0) %>%
      arrange(FDR, desc(abs(logFC)), desc(consistency_n))
  }

  unique(c(
    head(up_tbl$gene, top_per_direction),
    head(down_tbl$gene, top_per_direction)
  ))
}

####################
# load data
####################
message("Loading PDO object and matched metadata ...")
pdos <- readRDS("PDOs_merged.rds")
state_vec <- readRDS("Auto_PDO_final_states.rds")
pdos$state <- state_vec[Cells(pdos)]

pdos_matched <- subset(pdos, subset = orig.ident %in% matched_samples)
pdos_matched$Patient <- str_extract(pdos_matched$orig.ident, "^SUR\\d+")
pdos_matched$Patient <- factor(pdos_matched$Patient, levels = patient_order)
pdos_matched$Patient_label <- factor(
  unname(patient_labels[as.character(pdos_matched$Patient)]),
  levels = unname(patient_labels[patient_order])
)
pdos_matched$Treatment <- ifelse(grepl("_Treated_", pdos_matched$orig.ident), "Treated", "Untreated")
pdos_matched$Treatment <- factor(pdos_matched$Treatment, levels = c("Untreated", "Treated"))
pdos_matched$state <- factor(as.character(pdos_matched$state), levels = state_levels_all)

meta_matched <- pdos_matched@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    Patient = factor(as.character(Patient), levels = patient_order),
    Patient_label = factor(as.character(Patient_label), levels = unname(patient_labels[patient_order])),
    Treatment = factor(as.character(Treatment), levels = c("Untreated", "Treated")),
    state = factor(as.character(state), levels = state_levels_all)
  )

####################
# gene sets and support scores
####################
message("Preparing Hallmark pathway gene sets ...")
hallmark_tbl <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- lapply(names(pathway_ids), function(gs_name) {
  unique(hallmark_tbl$gene_symbol[hallmark_tbl$gs_name == gs_name])
})
names(gene_sets) <- unname(pathway_ids)
gene_sets[["Interferon response"]] <- unique(
  hallmark_tbl$gene_symbol[
    hallmark_tbl$gs_name %in% c(
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE"
    )
  ]
)

message("Scoring cell-level support signatures for UMAP overlays ...")
data_mat <- GetAssayData(pdos_matched, assay = "RNA", layer = "data")
support_scores <- data.frame(
  cell = colnames(data_mat),
  injury_score = mean_score(data_mat, unique(c(gene_sets[["Apoptosis"]], gene_sets[["p53 pathway"]]))),
  adaptive_score = mean_score(
    data_mat,
    unique(c(gene_sets[["TNFa/NFkB"]], gene_sets[["EMT"]], gene_sets[["Unfolded protein response"]]))
  ),
  stringsAsFactors = FALSE
)

umap_embed <- as.data.frame(Embeddings(pdos_matched, reduction = "umap")) %>%
  rownames_to_column("cell") %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2)

umap_plot_df <- umap_embed %>%
  left_join(
    meta_matched %>%
      select(cell, state, Treatment, Patient_label),
    by = "cell"
  ) %>%
  left_join(support_scores, by = "cell")

p_umap_state <- build_discrete_umap(
  plot_df = umap_plot_df,
  value_col = "state",
  color_map = state_cols,
  title_text = "Matched cells: finalized states"
)

p_umap_treatment <- build_discrete_umap(
  plot_df = umap_plot_df,
  value_col = "Treatment",
  color_map = treatment_cols,
  title_text = "Matched cells: untreated vs FLOT-treated"
)

p_umap_injury <- build_continuous_umap(
  plot_df = umap_plot_df,
  value_col = "injury_score",
  title_text = "Matched cells: injury score"
)

p_umap_adaptive <- build_continuous_umap(
  plot_df = umap_plot_df,
  value_col = "adaptive_score",
  title_text = "Matched cells: adaptation score"
)

umap_panel <- (p_umap_state | p_umap_treatment) / (p_umap_injury | p_umap_adaptive) +
  plot_annotation(
    title = "PDO matched-pair support UMAPs",
    subtitle = "UMAP is supportive only; response interpretation below is based on paired state abundance and paired pseudobulk summaries."
  )

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_umap_support.pdf"),
  plot = umap_panel,
  width = 14,
  height = 11,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_umap_support.png"),
  plot = umap_panel,
  width = 14,
  height = 11,
  units = "in",
  dpi = 300
)

####################
# state abundance
####################
message("Computing paired state abundance shifts ...")
cell_counts <- meta_matched %>%
  count(Patient, Patient_label, Treatment, orig.ident, state, name = "cell_n") %>%
  complete(
    nesting(Patient, Patient_label, Treatment, orig.ident),
    state = state_levels_all,
    fill = list(cell_n = 0L)
  ) %>%
  group_by(Patient, Patient_label, Treatment, orig.ident) %>%
  mutate(total_cells = sum(cell_n)) %>%
  ungroup() %>%
  mutate(
    pct = 100 * cell_n / total_cells,
    state = factor(state, levels = state_levels_all)
  )

write.csv(
  cell_counts,
  file.path(out_dir, "Auto_pdo_flot_matched_state_cell_counts.csv"),
  row.names = FALSE
)

bar_totals <- cell_counts %>%
  distinct(Patient, Patient_label, Treatment, total_cells)

p_state_bar <- ggplot(cell_counts, aes(Treatment, pct, fill = state)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
  geom_text(
    data = subset(cell_counts, pct >= 8),
    aes(label = sprintf("%.1f%%", pct)),
    position = position_stack(vjust = 0.5),
    size = 3.0,
    color = "black"
  ) +
  geom_text(
    data = bar_totals,
    aes(x = Treatment, y = 104, label = paste0("N=", comma(total_cells))),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  facet_wrap(~ Patient_label, nrow = 1) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  scale_x_discrete(labels = c(Untreated = "Untreated", Treated = "FLOT-treated")) +
  scale_y_continuous(
    limits = c(0, 108),
    expand = expansion(mult = c(0, 0)),
    labels = label_number(accuracy = 1)
  ) +
  labs(
    title = "State composition per matched PDO pair",
    subtitle = "All finalized malignant-state calls are shown; grey/black capture unresolved and hybrid cells.",
    x = NULL,
    y = "Fraction of malignant cells (%)",
    fill = NULL
  ) +
  flot_theme(10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(face = "bold")
  )

abundance_main <- cell_counts %>%
  filter(state %in% state_levels_main) %>%
  mutate(state = factor(as.character(state), levels = state_levels_main))

abundance_pairs <- abundance_main %>%
  select(Patient, Patient_label, Treatment, state, cell_n, total_cells, pct) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(cell_n, total_cells, pct)
  ) %>%
  mutate(
    untreated_prop_adj = (cell_n_Untreated + 0.5) / (total_cells_Untreated + 0.5 * length(state_levels_all)),
    treated_prop_adj = (cell_n_Treated + 0.5) / (total_cells_Treated + 0.5 * length(state_levels_all)),
    log2FC = log2(treated_prop_adj / untreated_prop_adj),
    delta_pct = pct_Treated - pct_Untreated
  )

abundance_summary <- abundance_pairs %>%
  group_by(state) %>%
  summarise(
    median_log2FC = median(log2FC, na.rm = TRUE),
    median_delta_pct = median(delta_pct, na.rm = TRUE),
    .groups = "drop"
  )

state_plot_order <- rev(state_levels_main)
abundance_pairs$state_plot <- match(as.character(abundance_pairs$state), state_plot_order)
abundance_summary$state_plot <- match(as.character(abundance_summary$state), state_plot_order)

write.csv(
  abundance_pairs,
  file.path(out_dir, "Auto_pdo_flot_matched_state_abundance_changes.csv"),
  row.names = FALSE
)

p_abundance_lfc <- ggplot(
  abundance_pairs,
  aes(log2FC, state_plot)
) +
  geom_vline(xintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_segment(
    data = abundance_summary,
    aes(
      x = median_log2FC,
      xend = median_log2FC,
      y = state_plot - 0.28,
      yend = state_plot + 0.28
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.1
  ) +
  geom_point(
    aes(color = Patient_label),
    size = 2.9,
    alpha = 0.95,
    position = position_jitter(height = 0.09, width = 0)
  ) +
  scale_color_manual(values = patient_cols) +
  scale_y_continuous(
    breaks = seq_along(state_plot_order),
    labels = state_plot_order,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    title = "State abundance change after FLOT",
    subtitle = "Points are patient-level log2 fold-changes; black ticks show state medians.",
    x = "log2 fold-change in state fraction (FLOT-treated / untreated)",
    y = NULL,
    color = NULL
  ) +
  flot_theme(10) +
  theme(legend.position = "bottom")

state_abundance_panel <- p_state_bar / p_abundance_lfc + plot_layout(heights = c(1.25, 0.9))

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_state_abundance.pdf"),
  plot = state_abundance_panel,
  width = 15,
  height = 12,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_state_abundance.png"),
  plot = state_abundance_panel,
  width = 15,
  height = 12,
  units = "in",
  dpi = 300
)

####################
# pseudobulk aggregation
####################
message("Aggregating pseudobulk counts by sample and state ...")
main_cells <- meta_matched %>%
  filter(state %in% state_levels_main) %>%
  mutate(sample_state = paste(orig.ident, state, sep = "__"))

counts_mat <- GetAssayData(pdos_matched, assay = "RNA", layer = "counts")
group_factor <- factor(main_cells$sample_state, levels = unique(main_cells$sample_state))
group_mm <- sparse.model.matrix(~ 0 + group_factor)
colnames(group_mm) <- levels(group_factor)

pseudobulk_counts <- counts_mat[, main_cells$cell, drop = FALSE] %*% group_mm

pseudobulk_meta <- tibble(sample_state = colnames(pseudobulk_counts)) %>%
  tidyr::separate(sample_state, into = c("orig.ident", "state"), sep = "__", extra = "merge", remove = FALSE) %>%
  mutate(
    state = factor(state, levels = state_levels_main),
    Patient = factor(str_extract(orig.ident, "^SUR\\d+"), levels = patient_order),
    Patient_label = factor(unname(patient_labels[as.character(Patient)]), levels = unname(patient_labels[patient_order])),
    Treatment = factor(ifelse(grepl("_Treated_", orig.ident), "Treated", "Untreated"),
      levels = c("Untreated", "Treated")
    ),
    patient_short = unname(patient_short[as.character(Patient)]),
    cell_n = as.integer(colSums(group_mm))
  )

write.csv(
  pseudobulk_meta,
  file.path(out_dir, "Auto_pdo_flot_matched_pseudobulk_samples.csv"),
  row.names = FALSE
)

pair_eligibility <- pseudobulk_meta %>%
  mutate(valid_for_pseudobulk = cell_n >= min_cells_per_pseudobulk) %>%
  count(state, Patient, valid_for_pseudobulk, name = "sample_n") %>%
  pivot_wider(
    names_from = valid_for_pseudobulk,
    values_from = sample_n,
    values_fill = 0
  ) %>%
  rename(
    invalid_sample_n = `FALSE`,
    valid_sample_n = `TRUE`
  ) %>%
  mutate(valid_pair = valid_sample_n == 2)

write.csv(
  pair_eligibility,
  file.path(out_dir, "Auto_pdo_flot_matched_pseudobulk_pair_eligibility.csv"),
  row.names = FALSE
)

pb_dge_all <- DGEList(counts = pseudobulk_counts)
pb_dge_all <- calcNormFactors(pb_dge_all)
pb_logcpm_all <- cpm(pb_dge_all, log = TRUE, prior.count = 2)
pb_logcpm_z <- zscore_rows(pb_logcpm_all)

pathway_score_mat <- score_gene_sets(pb_logcpm_z, gene_sets)
pathway_scores_long <- as.data.frame(t(pathway_score_mat)) %>%
  rownames_to_column("sample_state") %>%
  left_join(
    pseudobulk_meta %>% select(sample_state, orig.ident, state, Patient, Patient_label, Treatment, cell_n),
    by = "sample_state"
  ) %>%
  pivot_longer(
    cols = all_of(pathway_order),
    names_to = "pathway",
    values_to = "score"
  )

write.csv(
  pathway_scores_long,
  file.path(out_dir, "Auto_pdo_flot_matched_pseudobulk_pathway_scores.csv"),
  row.names = FALSE
)

valid_pathway_scores <- pathway_scores_long %>%
  filter(cell_n >= min_cells_per_pseudobulk)

pathway_delta <- valid_pathway_scores %>%
  select(Patient, Patient_label, state, pathway, Treatment, score, cell_n) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(score, cell_n)
  ) %>%
  filter(!is.na(score_Untreated), !is.na(score_Treated)) %>%
  mutate(
    delta_score = score_Treated - score_Untreated,
    pair_id = paste(state, Patient, sep = "__")
  ) %>%
  arrange(factor(state, levels = state_levels_main), factor(Patient, levels = patient_order), factor(pathway, levels = pathway_order))

write.csv(
  pathway_delta,
  file.path(out_dir, "Auto_pdo_flot_matched_pathway_deltas.csv"),
  row.names = FALSE
)

pathway_delta_matrix <- pathway_delta %>%
  mutate(pair_label = as.character(Patient)) %>%
  select(pathway, pair_id, delta_score) %>%
  pivot_wider(names_from = pair_id, values_from = delta_score) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

pair_meta <- pathway_delta %>%
  distinct(pair_id, state, Patient, Patient_label) %>%
  mutate(
    state = factor(state, levels = state_levels_main),
    Patient = factor(Patient, levels = patient_order),
    Patient_label = factor(Patient_label, levels = unname(patient_labels[patient_order]))
  ) %>%
  arrange(state, Patient)

pathway_delta_matrix <- pathway_delta_matrix[pathway_order, pair_meta$pair_id, drop = FALSE]
colnames(pathway_delta_matrix) <- as.character(pair_meta$Patient)

top_pathway_annotation <- HeatmapAnnotation(
  Patient = as.character(pair_meta$Patient_label),
  col = list(Patient = patient_cols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontface = "bold")
)

pathway_heatmap <- Heatmap(
  pathway_delta_matrix,
  name = "Treated - untreated",
  col = colorRamp2(c(-0.8, 0, 0.8), c("#245F7B", "white", "#B63E2F")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_split = pair_meta$state,
  top_annotation = top_pathway_annotation,
  column_title = "State-resolved pathway response matrix",
  column_title_gp = gpar(fontface = "bold", fontsize = 12),
  heatmap_legend_param = list(title_position = "topcenter")
)

pdf(file.path(out_dir, "Auto_pdo_flot_matched_pathway_heatmap.pdf"), width = 15, height = 7)
draw(pathway_heatmap, merge_legend = TRUE)
dev.off()

png(file.path(out_dir, "Auto_pdo_flot_matched_pathway_heatmap.png"), width = 4200, height = 2100, res = 300)
draw(pathway_heatmap, merge_legend = TRUE)
dev.off()

####################
# composite response deltas
####################
message("Summarising composite pathway responses ...")
composite_delta <- bind_rows(lapply(names(composite_pathways), function(metric_name) {
  pathway_delta %>%
    filter(pathway %in% composite_pathways[[metric_name]]) %>%
    group_by(Patient, Patient_label, state) %>%
    summarise(
      metric = metric_name,
      delta_score = mean(delta_score, na.rm = TRUE),
      pathway_n = n_distinct(pathway),
      .groups = "drop"
    )
})) %>%
  mutate(
    state = factor(state, levels = state_levels_main),
    metric = factor(metric, levels = names(composite_pathways)),
    Patient = factor(Patient, levels = patient_order),
    Patient_label = factor(Patient_label, levels = unname(patient_labels[patient_order]))
  )

write.csv(
  composite_delta,
  file.path(out_dir, "Auto_pdo_flot_matched_composite_response_deltas.csv"),
  row.names = FALSE
)

composite_summary <- composite_delta %>%
  group_by(state, metric) %>%
  summarise(median_delta = median(delta_score, na.rm = TRUE), .groups = "drop")

p_composite <- ggplot(
  composite_delta,
  aes(state, delta_score, color = Patient_label, group = Patient_label)
) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_line(alpha = 0.45, linewidth = 0.45) +
  geom_point(size = 2.4) +
  geom_point(
    data = composite_summary,
    aes(x = state, y = median_delta),
    inherit.aes = FALSE,
    color = "black",
    shape = 18,
    size = 3.0
  ) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = patient_cols) +
  labs(
    title = "Within-state FLOT response scores",
    subtitle = "Patient-level paired pseudobulk deltas; black diamonds mark state medians.",
    x = NULL,
    y = "Composite pathway delta (treated - untreated)",
    color = NULL
  ) +
  flot_theme(10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 40, hjust = 1)
  )

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_composite_response.pdf"),
  plot = p_composite,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_composite_response.png"),
  plot = p_composite,
  width = 14,
  height = 10,
  units = "in",
  dpi = 300
)

summary_panel <- (p_state_bar / p_abundance_lfc / p_composite) +
  plot_layout(heights = c(1.25, 0.9, 1.1)) +
  plot_annotation(
    title = "Matched PDO FLOT response summary",
    subtitle = "Primary evidence is paired state abundance plus paired pseudobulk response scoring."
  )

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_summary_panels.pdf"),
  plot = summary_panel,
  width = 16,
  height = 18,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_summary_panels.png"),
  plot = summary_panel,
  width = 16,
  height = 18,
  units = "in",
  dpi = 300
)

####################
# heuristic state-fate summary
####################
message("Building heuristic state-fate summary ...")
fate_tbl <- abundance_summary %>%
  transmute(
    state = state,
    abundance_log2FC = median_log2FC,
    median_delta_pct = median_delta_pct
  ) %>%
  left_join(
    composite_summary %>%
      select(state, metric, median_delta) %>%
      pivot_wider(names_from = metric, values_from = median_delta),
    by = "state"
  ) %>%
  mutate(
    `Proliferation loss` = -`Cell-cycle / proliferation`,
    `Injury` = `Injury / checkpoint`,
    `Adaptation` = `Adaptive / persistence`,
    `Proteostasis` = `Proteostasis / transition`,
    contributing_pairs = sapply(state, function(st) {
      sum(pathway_delta$state == st & pathway_delta$pathway == "Apoptosis")
    }),
    heuristic_call = mapply(
      classify_state_fate,
      abundance_log2fc = abundance_log2FC,
      proliferation_loss = `Proliferation loss`,
      injury = `Injury`,
      adaptation = `Adaptation`,
      proteostasis = `Proteostasis`
    )
  ) %>%
  mutate(state = factor(state, levels = state_levels_main))

write.csv(
  fate_tbl,
  file.path(out_dir, "Auto_pdo_flot_matched_state_fate_summary.csv"),
  row.names = FALSE
)

fate_metric_long <- fate_tbl %>%
  select(
    state,
    abundance_log2FC,
    `Proliferation loss`,
    `Injury`,
    `Adaptation`,
    `Proteostasis`
  ) %>%
  rename(`Abundance log2FC` = abundance_log2FC) %>%
  pivot_longer(
    cols = -state,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("Abundance log2FC", "Proliferation loss", "Injury", "Adaptation", "Proteostasis")
    ),
    state = factor(state, levels = rev(state_levels_main))
  )

p_fate_metrics <- ggplot(fate_metric_long, aes(metric, state, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.1) +
  scale_fill_gradient2(
    low = "#245F7B",
    mid = "white",
    high = "#B63E2F",
    midpoint = 0
  ) +
  labs(
    title = "Heuristic state-fate summary",
    subtitle = "Positive proliferation loss means stronger treatment-associated cell-cycle suppression.",
    x = NULL,
    y = NULL,
    fill = "Median delta"
  ) +
  flot_theme(10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

fate_call_df <- bind_rows(
  fate_tbl %>%
    transmute(
      state = factor(state, levels = rev(state_levels_main)),
      panel = "Interpretation",
      label = heuristic_call
    ),
  fate_tbl %>%
    transmute(
      state = factor(state, levels = rev(state_levels_main)),
      panel = "Pairs",
      label = paste0("n=", contributing_pairs)
    )
)

p_fate_calls <- ggplot(fate_call_df, aes(panel, state)) +
  geom_tile(fill = "#F8F5EF", color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = label),
    size = 3.1,
    fontface = "bold",
    lineheight = 0.95
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(t = 36, r = 18, b = 18, l = 0)
  )

fate_panel <- p_fate_metrics | p_fate_calls

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_state_fate_summary.pdf"),
  plot = fate_panel,
  width = 14,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Auto_pdo_flot_matched_state_fate_summary.png"),
  plot = fate_panel,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300
)

####################
# paired state-specific pseudobulk DEG
####################
message("Running paired state-specific pseudobulk edgeR ...")
deg_summary <- list()
deg_results <- list()
recurrent_gene_tables <- list()

pdf(file.path(out_dir, "Auto_pdo_flot_matched_recurrent_deg_heatmaps.pdf"), width = 9, height = 8)

for (state_name in state_levels_main) {
  message("  DEG: ", state_name)

  state_meta <- pseudobulk_meta %>%
    filter(state == state_name, cell_n >= min_cells_per_pseudobulk) %>%
    group_by(Patient) %>%
    filter(n_distinct(Treatment) == 2) %>%
    ungroup() %>%
    mutate(
      Patient = factor(as.character(Patient), levels = intersect(patient_order, unique(as.character(Patient)))),
      Treatment = factor(Treatment, levels = c("Untreated", "Treated"))
    ) %>%
    arrange(Patient, Treatment)

  valid_pairs <- n_distinct(state_meta$Patient)
  state_stub <- state_key(state_name)

  if (valid_pairs < min_pairs_for_deg) {
    deg_summary[[state_name]] <- data.frame(
      state = state_name,
      valid_pairs = valid_pairs,
      tested_genes = 0,
      sig_fdr_0_05 = 0,
      sig_fdr_0_10 = 0,
      note = "Skipped: fewer than 3 eligible matched pairs",
      stringsAsFactors = FALSE
    )
    next
  }

  y <- DGEList(counts = pseudobulk_counts[, state_meta$sample_state, drop = FALSE])
  design <- model.matrix(~ Patient + Treatment, data = state_meta)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "TreatmentTreated")

  res_tbl <- topTags(qlf, n = Inf, sort.by = "PValue")$table %>%
    rownames_to_column("gene") %>%
    mutate(
      state = state_name,
      valid_pairs = valid_pairs,
      pair_ids = paste(as.character(unique(state_meta$Patient)), collapse = ";")
    )

  logcpm_state <- cpm(y, log = TRUE, prior.count = 2)
  patient_logfc <- sapply(as.character(unique(state_meta$Patient)), function(patient_id) {
    treated_col <- state_meta$sample_state[state_meta$Patient == patient_id & state_meta$Treatment == "Treated"]
    untreated_col <- state_meta$sample_state[state_meta$Patient == patient_id & state_meta$Treatment == "Untreated"]
    logcpm_state[, treated_col, drop = FALSE][, 1] - logcpm_state[, untreated_col, drop = FALSE][, 1]
  })
  if (is.null(dim(patient_logfc))) {
    patient_logfc <- matrix(patient_logfc, ncol = 1)
    rownames(patient_logfc) <- rownames(logcpm_state)
    colnames(patient_logfc) <- as.character(unique(state_meta$Patient))
  }

  consistent_n <- rowSums(sign(patient_logfc[res_tbl$gene, , drop = FALSE]) == sign(res_tbl$logFC), na.rm = TRUE)
  res_tbl$consistent_direction_n <- consistent_n
  res_tbl$patient_median_logFC <- apply(patient_logfc[res_tbl$gene, , drop = FALSE], 1, median, na.rm = TRUE)

  deg_results[[state_name]] <- res_tbl
  write.csv(
    res_tbl,
    file.path(deg_dir, paste0("Auto_pdo_flot_matched_deg_", state_stub, ".csv")),
    row.names = FALSE
  )

  patient_logfc_df <- as.data.frame(patient_logfc) %>%
    rownames_to_column("gene")
  write.csv(
    patient_logfc_df,
    file.path(deg_dir, paste0("Auto_pdo_flot_matched_patient_logFC_", state_stub, ".csv")),
    row.names = FALSE
  )

  top_genes <- select_recurrent_genes(res_tbl, patient_logfc, top_per_direction = top_deg_per_direction)
  recurrent_gene_tables[[state_name]] <- res_tbl %>% filter(gene %in% top_genes)

  write.csv(
    recurrent_gene_tables[[state_name]],
    file.path(deg_dir, paste0("Auto_pdo_flot_matched_recurrent_genes_", state_stub, ".csv")),
    row.names = FALSE
  )

  deg_summary[[state_name]] <- data.frame(
    state = state_name,
    valid_pairs = valid_pairs,
    tested_genes = nrow(res_tbl),
    sig_fdr_0_05 = sum(res_tbl$FDR < 0.05, na.rm = TRUE),
    sig_fdr_0_10 = sum(res_tbl$FDR < 0.10, na.rm = TRUE),
    note = "Paired edgeR completed",
    stringsAsFactors = FALSE
  )

  if (length(top_genes) == 0) {
    next
  }

  heat_mat <- patient_logfc[top_genes, as.character(unique(state_meta$Patient)), drop = FALSE]
  heat_mat <- heat_mat[match(top_genes, rownames(heat_mat)), , drop = FALSE]
  colnames(heat_mat) <- unname(patient_short[colnames(heat_mat)])

  direction_vec <- ifelse(res_tbl$logFC[match(rownames(heat_mat), res_tbl$gene)] > 0, "Up in treated", "Down in treated")
  direction_vec <- factor(direction_vec, levels = c("Up in treated", "Down in treated"))

  row_ha <- rowAnnotation(
    Direction = direction_vec,
    col = list(Direction = c("Up in treated" = "#B63E2F", "Down in treated" = "#245F7B")),
    show_annotation_name = TRUE
  )

  top_ha <- HeatmapAnnotation(
    Patient = as.character(factor(unname(patient_labels[match(colnames(heat_mat), patient_short)]), levels = unname(patient_labels[patient_order]))),
    col = list(Patient = patient_cols),
    show_annotation_name = TRUE
  )

  ht <- Heatmap(
    heat_mat,
    name = "Per-patient\nlogFC",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#245F7B", "white", "#B63E2F")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),
    top_annotation = top_ha,
    right_annotation = row_ha,
    column_title = paste0(state_name, " recurrent paired pseudobulk DEGs"),
    column_title_gp = gpar(fontface = "bold", fontsize = 12)
  )

  draw(ht, merge_legend = TRUE)

  png(file.path(deg_dir, paste0("Auto_pdo_flot_matched_recurrent_deg_", state_stub, ".png")), width = 2400, height = 2200, res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
}

dev.off()

deg_summary_df <- bind_rows(deg_summary)
write.csv(
  deg_summary_df,
  file.path(out_dir, "Auto_pdo_flot_matched_deg_summary.csv"),
  row.names = FALSE
)

####################
# output notes and cache
####################
presentation_notes <- c(
  "# PDO matched FLOT presentation notes",
  "",
  "- Main panels: `Auto_pdo_flot_matched_summary_panels.*`, `Auto_pdo_flot_matched_pathway_heatmap.*`, and `Auto_pdo_flot_matched_recurrent_deg_heatmaps.pdf`.",
  "- State abundance is shown with paired composition plots plus patient-level log2 fold-change points.",
  "- Within-state response is summarised with paired pseudobulk composite deltas rather than cell-level violin plots to avoid pseudoreplication.",
  "- The pathway heatmap shows paired treated-minus-untreated deltas for each patient-state combination.",
  "- Alluvial plots were intentionally not used because these are matched samples, not lineage-traced cells; an alluvial would overstate redistribution.",
  "- Recurrent DEG heatmaps are based on paired edgeR pseudobulk analysis within each state and highlight genes with consistent direction across patients.",
  "- `Auto_pdo_flot_matched_state_fate_summary.*` is heuristic and intended as an interpretation panel, not a standalone statistical test."
)

writeLines(
  presentation_notes,
  con = file.path(out_dir, "Auto_pdo_flot_matched_presentation_notes.md")
)

saveRDS(
  list(
    cell_counts = cell_counts,
    abundance_pairs = abundance_pairs,
    pseudobulk_meta = pseudobulk_meta,
    pathway_delta = pathway_delta,
    composite_delta = composite_delta,
    fate_summary = fate_tbl,
    deg_summary = deg_summary_df,
    deg_results = deg_results,
    recurrent_genes = recurrent_gene_tables
  ),
  file = file.path(out_dir, "Auto_pdo_flot_matched_response_results.rds")
)

message("Finished. Outputs written to: ", file.path(getwd(), out_dir))
