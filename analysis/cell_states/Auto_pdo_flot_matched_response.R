####################
# Auto_pdo_flot_matched_response.R
# Merged FLOT matched-pair response analysis for PDOs.
# Combines original analysis and presentation PDF generation.
#
# Produces:
#   1. All original individual plot files (abundance, MP, hybrid, pathway, DEG, etc.)
#   2. Presentation PDF (pathway heatmap, composite, abundance, MP, hybrid, DEG, fate, UMAP)
#   3. NEW: Node plot PDF — state/hybrid network for untreated vs treated (5 pages, scaled)
#   4. NEW: Paired boxplot PDF — state abundance, hybrid abundance, MP expression (3 pages)
#   5. NEW: Improved pathway heatmap PDF — with CCSIG, Lineage block, unified color scale, sig labels
#
# Inputs: PDOs_merged.rds, Auto_PDO_final_states.rds, UCell_scores_filtered.rds,
#         Auto_PDO_mp_adj_noreg.rds, geneNMF_metaprograms_nMP_13.rds,
#         Cell_Cycle_Genes.csv (for CC signature)
# Env: dmtcp
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
  library(grid)
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
patient_short <- c(SUR1070 = "1070", SUR1090 = "1090", SUR1072 = "1072", SUR1181 = "1181")

state_levels_main <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)
state_levels_all <- c(state_levels_main, "Unresolved", "Hybrid")

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
  "SUR1070 (Pre-responder)" = "#0072B2",
  "SUR1090 (Post-responder)" = "#56B4E9",
  "SUR1072 (Post-nonresponder)" = "#D55E00",
  "SUR1181 (Pre-nonresponder)" = "#E69F00"
)

treatment_cols <- c(Untreated = "#D58B2D", Treated = "#374151")

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

pathway_order <- c(unname(pathway_ids), "Interferon response")

composite_pathways <- list(
  "Cell-cycle / proliferation" = c("E2F targets", "G2M checkpoint"),
  "Injury / checkpoint" = c("Apoptosis", "p53 pathway", "DNA repair"),
  "Adaptive / persistence" = c("TNFa/NFkB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response"),
  "Proteostasis / transition" = c("Unfolded protein response", "Oxidative phosphorylation")
)

pathway_block_df <- data.frame(
  pathway = c(
    "E2F targets", "G2M checkpoint",
    "Apoptosis", "p53 pathway", "DNA repair",
    "TNFa/NFkB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response",
    "Unfolded protein response", "Oxidative phosphorylation"
  ),
  block = c(
    rep("Cell-cycle", 2), rep("Injury", 3),
    rep("Adaptive", 5), rep("Proteostasis", 2)
  ),
  stringsAsFactors = FALSE
)
pathway_order_grouped <- pathway_block_df$pathway

block_cols <- c(
  "Cell-cycle" = "#F0C75E", "Injury" = "#E07B54",
  "Adaptive" = "#5DAA68", "Proteostasis" = "#5B9BD5"
)

mp_descriptions <- c(
  "MP6"  = "MP6_G2M Cell Cycle",
  "MP7"  = "MP7_DNA repair",
  "MP5"  = "MP5_MYC-related Proliferation",
  "MP1"  = "MP1_G2M checkpoint",
  "MP3"  = "MP3_G1S Cell Cycle",
  "MP8"  = "MP8_Columnar Progenitor",
  "MP10" = "MP10_Inflammatory Stress Epi.",
  "MP9"  = "MP9_ECM Remodeling Epi.",
  "MP4"  = "MP4_Intestinal Metaplasia"
)

state_mp_map <- list(
  "Classic Proliferative" = c("MP5", "MP6", "MP7", "MP1", "MP3"),
  "Basal to Intest. Meta" = c("MP4"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9")
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

pres_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title    = element_text(face = "bold", size = base_size + 3),
      strip.text    = element_text(face = "bold"),
      legend.title  = element_text(face = "bold")
    )
}

state_key <- function(x) gsub("_+$", "", gsub("[^A-Za-z0-9]+", "_", x))

normalise_final_state_labels <- function(state_vec) {
  state_names <- names(state_vec)
  state_vec <- as.character(state_vec)
  names(state_vec) <- state_names
  state_vec[state_vec %in% c("Basal to Intestinal Metaplasia")] <- "Basal to Intest. Meta"
  state_vec[state_vec %in% c("3CA_mp_30 Respiration 1", "3CA_mp_3 Cell Cylce HMG-rich", "3CA_mp_3 Cell Cycle HMG-rich")] <- "Classic Proliferative"
  state_vec[state_vec %in% c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III", "3CA_mp_17 EMT-III")] <- "3CA_EMT_and_Protein_maturation"
  state_vec
}

mean_score <- function(mat, genes) {
  genes_use <- intersect(unique(genes), rownames(mat))
  if (length(genes_use) == 0) return(rep(NA_real_, ncol(mat)))
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
    coord_equal() + flot_theme(10) + theme(legend.position = "right")
}

build_continuous_umap <- function(plot_df, value_col, title_text) {
  ggplot(plot_df, aes(UMAP_1, UMAP_2, color = .data[[value_col]])) +
    geom_point(size = 0.22, alpha = 0.82) +
    scale_color_gradientn(colors = c("#1B365D", "#4F9DD9", "#F2E6B6", "#B33C2F"), na.value = "grey85") +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    coord_equal() + flot_theme(10) + theme(legend.position = "right")
}

score_gene_sets <- function(expr_mat, gene_sets) {
  score_mat <- sapply(gene_sets, function(genes) mean_score(expr_mat, genes))
  score_mat <- t(score_mat)
  rownames(score_mat) <- names(gene_sets)
  score_mat
}

classify_state_fate <- function(abundance_log2fc, proliferation_loss, injury, adaptation, proteostasis) {
  if (!is.na(abundance_log2fc) && abundance_log2fc <= -0.35 &&
      ((!is.na(injury) && injury >= 0.15) || (!is.na(proliferation_loss) && proliferation_loss >= 0.15))) {
    return("Likely sensitive")
  }
  if (!is.na(abundance_log2fc) && abundance_log2fc >= 0.10 && !is.na(adaptation) && adaptation >= 0.15) {
    return("Adaptive / persistent")
  }
  if (!is.na(proteostasis) && proteostasis >= 0.15 && !is.na(adaptation) && adaptation >= 0.05) {
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
  if (length(common_genes) == 0) return(character(0))
  patient_logfc <- patient_logfc[common_genes, , drop = FALSE]
  res_sub <- res_tbl %>%
    filter(gene %in% common_genes) %>%
    mutate(
      consistency_n = rowSums(sign(patient_logfc[gene, , drop = FALSE]) == sign(logFC), na.rm = TRUE),
      patient_median_logFC = apply(patient_logfc[gene, , drop = FALSE], 1, median, na.rm = TRUE)
    )
  min_consistency <- max(2, ceiling(0.75 * ncol(patient_logfc)))
  up_tbl <- res_sub %>% filter(logFC > 0, consistency_n >= min_consistency) %>% arrange(FDR, desc(abs(logFC)))
  down_tbl <- res_sub %>% filter(logFC < 0, consistency_n >= min_consistency) %>% arrange(FDR, desc(abs(logFC)))
  if (nrow(up_tbl) < top_per_direction) up_tbl <- res_sub %>% filter(logFC > 0) %>% arrange(FDR, desc(abs(logFC)))
  if (nrow(down_tbl) < top_per_direction) down_tbl <- res_sub %>% filter(logFC < 0) %>% arrange(FDR, desc(abs(logFC)))
  unique(c(head(up_tbl$gene, top_per_direction), head(down_tbl$gene, top_per_direction)))
}


####################
# load data
####################
message("Loading PDO object and matched metadata ...")
pdos <- readRDS("PDOs_merged.rds")
state_vec <- readRDS("Auto_PDO_final_states.rds")
state_vec <- normalise_final_state_labels(state_vec)
pdos$state <- state_vec[Cells(pdos)]

pdos_matched <- subset(pdos, subset = orig.ident %in% matched_samples)
matched_state_vec <- normalise_final_state_labels(state_vec[Cells(pdos_matched)])
pdos_matched$state <- factor(matched_state_vec, levels = state_levels_all)
pdos_matched$Patient <- str_extract(pdos_matched$orig.ident, "^SUR\\d+")
pdos_matched$Patient <- factor(pdos_matched$Patient, levels = patient_order)
pdos_matched$Patient_label <- factor(
  unname(patient_labels[as.character(pdos_matched$Patient)]),
  levels = unname(patient_labels[patient_order])
)
pdos_matched$Treatment <- ifelse(grepl("_Treated_", pdos_matched$orig.ident), "Treated", "Untreated")
pdos_matched$Treatment <- factor(pdos_matched$Treatment, levels = c("Untreated", "Treated"))

meta_matched <- pdos_matched@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    Patient = factor(as.character(Patient), levels = patient_order),
    Patient_label = factor(as.character(Patient_label), levels = unname(patient_labels[patient_order])),
    Treatment = factor(as.character(Treatment), levels = c("Untreated", "Treated")),
    state = factor(matched_state_vec[cell], levels = state_levels_all)
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
    hallmark_tbl$gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
  ]
)

message("Scoring cell-level support signatures ...")
data_mat <- GetAssayData(pdos_matched, assay = "RNA", layer = "data")
support_scores <- data.frame(
  cell = colnames(data_mat),
  injury_score = mean_score(data_mat, unique(c(gene_sets[["Apoptosis"]], gene_sets[["p53 pathway"]]))),
  adaptive_score = mean_score(data_mat, unique(c(gene_sets[["TNFa/NFkB"]], gene_sets[["EMT"]], gene_sets[["Unfolded protein response"]]))),
  stringsAsFactors = FALSE
)

umap_embed <- as.data.frame(Embeddings(pdos_matched, reduction = "umap")) %>%
  rownames_to_column("cell") %>% rename(UMAP_1 = umap_1, UMAP_2 = umap_2)
umap_plot_df <- umap_embed %>%
  left_join(meta_matched %>% select(cell, state, Treatment, Patient_label), by = "cell") %>%
  left_join(support_scores, by = "cell")

# UMAP plots
p_umap_state <- build_discrete_umap(umap_plot_df, "state", state_cols, "Matched cells: finalized states")
p_umap_treatment <- build_discrete_umap(umap_plot_df, "Treatment", treatment_cols, "Matched cells: untreated vs FLOT-treated")
p_umap_injury <- build_continuous_umap(umap_plot_df, "injury_score", "Matched cells: injury score")
p_umap_adaptive <- build_continuous_umap(umap_plot_df, "adaptive_score", "Matched cells: adaptation score")
umap_panel <- (p_umap_state | p_umap_treatment) / (p_umap_injury | p_umap_adaptive) +
  plot_annotation(title = "PDO matched-pair support UMAPs")
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_umap_support.pdf"), umap_panel, width = 14, height = 11)
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_umap_support.png"), umap_panel, width = 14, height = 11, dpi = 300)

####################
# state abundance
####################
message("Computing paired state abundance shifts ...")
cell_counts <- meta_matched %>%
  count(Patient, Patient_label, Treatment, orig.ident, state, name = "cell_n") %>%
  complete(nesting(Patient, Patient_label, Treatment, orig.ident), state = state_levels_all, fill = list(cell_n = 0L)) %>%
  group_by(Patient, Patient_label, Treatment, orig.ident) %>% mutate(total_cells = sum(cell_n)) %>% ungroup() %>%
  mutate(pct = 100 * cell_n / total_cells, state = factor(state, levels = state_levels_all))
write.csv(cell_counts, file.path(out_dir, "Auto_pdo_flot_matched_state_cell_counts.csv"), row.names = FALSE)

bar_totals <- cell_counts %>% distinct(Patient, Patient_label, Treatment, total_cells)
p_state_bar <- ggplot(cell_counts, aes(Treatment, pct, fill = state)) +
  geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
  geom_text(data = subset(cell_counts, pct >= 8), aes(label = sprintf("%.1f%%", pct)),
    position = position_stack(vjust = 0.5), size = 3.0, color = "black") +
  geom_text(data = bar_totals, aes(x = Treatment, y = 104, label = paste0("N=", comma(total_cells))),
    inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  facet_wrap(~ Patient_label, nrow = 1) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  scale_x_discrete(labels = c(Untreated = "Untreated", Treated = "FLOT-treated")) +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
  labs(title = "State composition per matched PDO pair", x = NULL, y = "Fraction of malignant cells (%)", fill = NULL) +
  flot_theme(10) + theme(legend.position = "bottom", axis.text.x = element_text(face = "bold"))

abundance_main <- cell_counts %>% filter(state %in% state_levels_main) %>%
  mutate(state = factor(as.character(state), levels = state_levels_main))
abundance_pairs <- abundance_main %>%
  select(Patient, Patient_label, Treatment, state, cell_n, total_cells, pct) %>%
  pivot_wider(names_from = Treatment, values_from = c(cell_n, total_cells, pct)) %>%
  mutate(
    untreated_prop_adj = (cell_n_Untreated + 0.5) / (total_cells_Untreated + 0.5 * length(state_levels_all)),
    treated_prop_adj = (cell_n_Treated + 0.5) / (total_cells_Treated + 0.5 * length(state_levels_all)),
    log2FC = log2(treated_prop_adj / untreated_prop_adj),
    delta_pct = pct_Treated - pct_Untreated
  )
abundance_summary <- abundance_pairs %>% group_by(state) %>%
  summarise(median_log2FC = median(log2FC, na.rm = TRUE), median_delta_pct = median(delta_pct, na.rm = TRUE), .groups = "drop")

state_plot_order <- state_levels_main
abundance_pairs$state <- factor(as.character(abundance_pairs$state), levels = state_plot_order)
abundance_summary$state <- factor(as.character(abundance_summary$state), levels = state_plot_order)
abundance_pairs$state_x <- match(as.character(abundance_pairs$state), state_plot_order)
abundance_summary$state_x <- match(as.character(abundance_summary$state), state_plot_order)
write.csv(abundance_pairs, file.path(out_dir, "Auto_pdo_flot_matched_state_abundance_changes.csv"), row.names = FALSE)

p_abundance_lfc <- ggplot(abundance_pairs, aes(state_x, log2FC)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = abundance_summary, aes(x = as.numeric(state) - 0.28, xend = as.numeric(state) + 0.28, y = median_log2FC, yend = median_log2FC), inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_point(aes(color = Patient_label), size = 2.9, alpha = 0.95, position = position_jitter(width = 0.09, height = 0)) +
  geom_text(data = abundance_summary, aes(x = state_x, y = median_log2FC, label = sprintf("%.2f", median_log2FC)), inherit.aes = FALSE, vjust = -0.9, size = 3.1, fontface = "bold") +
  scale_color_manual(values = patient_cols) +
  scale_x_continuous(breaks = seq_along(state_plot_order), labels = state_plot_order) +
  labs(title = "State abundance change after FLOT", x = NULL, y = "log2 FC (FLOT / untreated)", color = NULL) +
  flot_theme(10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))

state_abundance_panel <- p_state_bar / p_abundance_lfc + plot_layout(heights = c(1.25, 0.9))
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_state_abundance.pdf"), state_abundance_panel, width = 15, height = 12)
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_state_abundance.png"), state_abundance_panel, width = 15, height = 12, dpi = 300)

####################
# per MP expression change
####################
message("Computing per-MP expression changes ...")
ucell_scores <- readRDS("UCell_scores_filtered.rds")
common_cells_mp <- intersect(colnames(pdos_matched), rownames(ucell_scores))

geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- rev(unique(ordered_clusters))
mp_tree_order_names <- paste0("MP", mp_tree_order)

ordered_mp_list <- character()
for (grp in c("Classic Proliferative", "Basal to Intest. Meta", "SMG-like Metaplasia", "Stress-adaptive")) {
  mps_in_grp <- state_mp_map[[grp]]
  mps_sorted <- mps_in_grp[order(match(mps_in_grp, mp_tree_order_names))]
  ordered_mp_list <- c(ordered_mp_list, mps_sorted)
}
ordered_mp_list <- intersect(ordered_mp_list, colnames(ucell_scores))
mp_plot_order <- ordered_mp_list

mp_sample_scores <- as.data.frame(ucell_scores[common_cells_mp, , drop = FALSE]) %>%
  rownames_to_column("cell") %>%
  left_join(meta_matched %>% select(cell, Patient, Patient_label, Treatment), by = "cell") %>%
  pivot_longer(cols = starts_with("MP"), names_to = "MP", values_to = "score") %>%
  group_by(Patient, Patient_label, Treatment, MP) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  filter(MP %in% ordered_mp_list)

mp_pairs <- mp_sample_scores %>%
  pivot_wider(names_from = Treatment, values_from = mean_score) %>%
  mutate(log2FC = log2((Treated + 1e-4) / (Untreated + 1e-4)))
mp_summary <- mp_pairs %>% group_by(MP) %>% summarise(median_log2FC = median(log2FC, na.rm = TRUE), .groups = "drop")
mp_pairs$MP <- factor(as.character(mp_pairs$MP), levels = mp_plot_order)
mp_summary$MP <- factor(as.character(mp_summary$MP), levels = mp_plot_order)
mp_pairs$MP_x <- match(as.character(mp_pairs$MP), mp_plot_order)
mp_summary$MP_x <- match(as.character(mp_summary$MP), mp_plot_order)
mp_x_labels <- mp_descriptions[mp_plot_order]

p_mp_lfc <- ggplot(mp_pairs, aes(MP_x, log2FC)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = mp_summary, aes(x = MP_x - 0.28, xend = MP_x + 0.28, y = median_log2FC, yend = median_log2FC), inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_point(aes(color = Patient_label), size = 2.9, alpha = 0.95, position = position_jitter(width = 0.09, height = 0)) +
  geom_text(data = mp_summary, aes(x = MP_x, y = median_log2FC, label = sprintf("%.2f", median_log2FC)), inherit.aes = FALSE, vjust = -0.9, size = 3.0, fontface = "bold") +
  scale_color_manual(values = patient_cols) +
  scale_x_continuous(breaks = seq_along(mp_plot_order), labels = mp_x_labels) +
  labs(title = "MP expression change after FLOT", x = NULL, y = "log2 FC in mean MP UCell score", color = NULL) +
  flot_theme(10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"))

ggsave(file.path(out_dir, "Auto_pdo_flot_matched_mp_expression.pdf"), p_mp_lfc, width = 10, height = 7)
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_mp_expression.png"), p_mp_lfc, width = 10, height = 7, dpi = 300)
write.csv(mp_pairs, file.path(out_dir, "Auto_pdo_flot_matched_mp_expression_changes.csv"), row.names = FALSE)


####################
# hybrid pairwise state abundance change
####################
message("Computing hybrid pairwise state abundance changes ...")
mp_adj_noncc <- readRDS("Auto_PDO_mp_adj_noreg.rds")
hybrid_state_vec <- normalise_final_state_labels(readRDS("Auto_PDO_final_states.rds"))
matched_hybrid_cells <- base::intersect(names(hybrid_state_vec)[hybrid_state_vec == "Hybrid"], Cells(pdos_matched))

hybrid_meta_source <- pdos_matched@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    Patient = factor(as.character(Patient), levels = patient_order),
    Patient_label = factor(as.character(Patient_label), levels = unname(patient_labels[patient_order])),
    Treatment = factor(as.character(Treatment), levels = c("Untreated", "Treated")),
    state = factor(hybrid_state_vec[cell], levels = state_levels_all)
  )

common_hybrid <- base::intersect(matched_hybrid_cells, rownames(mp_adj_noncc))
hybrid_base_states <- c("Classic Proliferative", "Basal to Intest. Meta", "SMG-like Metaplasia", "Stress-adaptive")

state_groups_hybrid <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9")
)

if (length(common_hybrid) == 0) stop("No matched hybrid cells overlapped.")

group_max_list <- lapply(state_groups_hybrid, function(mps) {
  mps_avail <- base::intersect(mps, colnames(mp_adj_noncc))
  if (length(mps_avail) == 0) return(rep(NA_real_, length(common_hybrid)))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc[common_hybrid, mps_avail, drop = TRUE]))
  apply(mp_adj_noncc[common_hybrid, mps_avail, drop = FALSE], 1, max)
})
group_max <- do.call(cbind, group_max_list)
group_max <- matrix(as.numeric(group_max), nrow = length(common_hybrid), ncol = length(state_groups_hybrid),
                    dimnames = list(common_hybrid, names(state_groups_hybrid)))

hybrid_subtypes <- apply(group_max, 1, function(scores) {
  top2 <- names(sort(scores, decreasing = TRUE))[1:2]
  top2_ordered <- top2[order(match(top2, hybrid_base_states))]
  paste(top2_ordered, collapse = "
-
")
})

hybrid_pairs_list <- combn(hybrid_base_states, 2, simplify = FALSE)
hybrid_levels_ordered <- vapply(hybrid_pairs_list, function(x) paste(x, collapse = "
-
"), character(1))

hybrid_meta <- hybrid_meta_source %>% filter(cell %in% common_hybrid)
hybrid_meta$hybrid_subtype <- factor(hybrid_subtypes[hybrid_meta$cell], levels = hybrid_levels_ordered)

hybrid_counts <- hybrid_meta %>%
  count(Patient, Patient_label, Treatment, orig.ident, hybrid_subtype, name = "cell_n") %>%
  complete(nesting(Patient, Patient_label, Treatment, orig.ident), hybrid_subtype = hybrid_levels_ordered, fill = list(cell_n = 0L)) %>%
  group_by(Patient, Patient_label, Treatment, orig.ident) %>% mutate(total_hybrid_cells = sum(cell_n)) %>% ungroup() %>%
  mutate(pct = 100 * cell_n / total_hybrid_cells)

hybrid_pairs_abund <- hybrid_counts %>%
  select(Patient, Patient_label, Treatment, hybrid_subtype, cell_n, total_hybrid_cells, pct) %>%
  pivot_wider(names_from = Treatment, values_from = c(cell_n, total_hybrid_cells, pct)) %>%
  mutate(
    untreated_prop_adj = (cell_n_Untreated + 0.5) / (total_hybrid_cells_Untreated + 0.5 * length(hybrid_levels_ordered)),
    treated_prop_adj = (cell_n_Treated + 0.5) / (total_hybrid_cells_Treated + 0.5 * length(hybrid_levels_ordered)),
    log2FC = log2(treated_prop_adj / untreated_prop_adj)
  )

hybrid_summary <- hybrid_pairs_abund %>% group_by(hybrid_subtype) %>%
  summarise(median_log2FC = median(log2FC, na.rm = TRUE), .groups = "drop")
hybrid_pairs_abund$hybrid_subtype <- factor(as.character(hybrid_pairs_abund$hybrid_subtype), levels = hybrid_levels_ordered)
hybrid_summary$hybrid_subtype <- factor(as.character(hybrid_summary$hybrid_subtype), levels = hybrid_levels_ordered)
hybrid_pairs_abund$hybrid_x <- match(as.character(hybrid_pairs_abund$hybrid_subtype), hybrid_levels_ordered)
hybrid_summary$hybrid_x <- match(as.character(hybrid_summary$hybrid_subtype), hybrid_levels_ordered)

p_hybrid_lfc <- ggplot(hybrid_pairs_abund, aes(hybrid_x, log2FC)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = hybrid_summary, aes(x = hybrid_x - 0.28, xend = hybrid_x + 0.28, y = median_log2FC, yend = median_log2FC), inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_point(aes(color = Patient_label), size = 2.9, alpha = 0.95, position = position_jitter(width = 0.09, height = 0)) +
  geom_text(data = hybrid_summary, aes(hybrid_x, median_log2FC, label = sprintf("%.2f", median_log2FC)), inherit.aes = FALSE, vjust = -0.9, size = 3.0, fontface = "bold") +
  scale_color_manual(values = patient_cols) +
  scale_x_continuous(breaks = seq_along(hybrid_levels_ordered), labels = hybrid_levels_ordered) +
  labs(title = "Hybrid pairwise abundance change after FLOT", x = NULL, y = "log2 FC in hybrid fraction", color = NULL) +
  flot_theme(10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", lineheight = 0.8))

ggsave(file.path(out_dir, "Auto_pdo_flot_matched_hybrid_abundance.pdf"), p_hybrid_lfc, width = 12, height = 7)
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_hybrid_abundance.png"), p_hybrid_lfc, width = 12, height = 7, dpi = 300)
write.csv(hybrid_pairs_abund, file.path(out_dir, "Auto_pdo_flot_matched_hybrid_abundance_changes.csv"), row.names = FALSE)

####################
# pseudobulk aggregation
####################
message("Aggregating pseudobulk counts ...")
main_cells <- meta_matched %>% filter(state %in% state_levels_main) %>%
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
    Treatment = factor(ifelse(grepl("_Treated_", orig.ident), "Treated", "Untreated"), levels = c("Untreated", "Treated")),
    patient_short = unname(patient_short[as.character(Patient)]),
    cell_n = as.integer(colSums(group_mm))
  )
write.csv(pseudobulk_meta, file.path(out_dir, "Auto_pdo_flot_matched_pseudobulk_samples.csv"), row.names = FALSE)

pb_dge_all <- DGEList(counts = pseudobulk_counts)
pb_dge_all <- calcNormFactors(pb_dge_all)
pb_logcpm_all <- cpm(pb_dge_all, log = TRUE, prior.count = 2)
pb_logcpm_z <- zscore_rows(pb_logcpm_all)

pathway_score_mat <- score_gene_sets(pb_logcpm_z, gene_sets)
pathway_scores_long <- as.data.frame(t(pathway_score_mat)) %>%
  rownames_to_column("sample_state") %>%
  left_join(pseudobulk_meta %>% select(sample_state, orig.ident, state, Patient, Patient_label, Treatment, cell_n), by = "sample_state") %>%
  pivot_longer(cols = all_of(pathway_order), names_to = "pathway", values_to = "score")
write.csv(pathway_scores_long, file.path(out_dir, "Auto_pdo_flot_matched_pseudobulk_pathway_scores.csv"), row.names = FALSE)

valid_pathway_scores <- pathway_scores_long %>% filter(cell_n >= min_cells_per_pseudobulk)
pathway_delta <- valid_pathway_scores %>%
  select(Patient, Patient_label, state, pathway, Treatment, score, cell_n) %>%
  pivot_wider(names_from = Treatment, values_from = c(score, cell_n)) %>%
  filter(!is.na(score_Untreated), !is.na(score_Treated)) %>%
  mutate(delta_score = score_Treated - score_Untreated, pair_id = paste(state, Patient, sep = "__")) %>%
  arrange(factor(state, levels = state_levels_main), factor(Patient, levels = patient_order), factor(pathway, levels = pathway_order))
write.csv(pathway_delta, file.path(out_dir, "Auto_pdo_flot_matched_pathway_deltas.csv"), row.names = FALSE)

# Pathway heatmap (original)
pathway_delta_matrix <- pathway_delta %>%
  select(pathway, pair_id, delta_score) %>%
  pivot_wider(names_from = pair_id, values_from = delta_score) %>%
  column_to_rownames("pathway") %>% as.matrix()
pair_meta <- pathway_delta %>%
  distinct(pair_id, state, Patient, Patient_label) %>%
  mutate(state = factor(state, levels = state_levels_main), Patient = factor(Patient, levels = patient_order)) %>%
  arrange(state, Patient)
pathway_delta_matrix <- pathway_delta_matrix[pathway_order, pair_meta$pair_id, drop = FALSE]
colnames(pathway_delta_matrix) <- as.character(pair_meta$Patient)

top_pathway_annotation <- HeatmapAnnotation(
  Patient = as.character(pair_meta$Patient_label),
  col = list(Patient = patient_cols),
  show_annotation_name = TRUE, annotation_name_gp = gpar(fontface = "bold")
)
pathway_heatmap <- Heatmap(
  pathway_delta_matrix, name = "Treated - untreated",
  col = colorRamp2(c(-0.8, 0, 0.8), c("#245F7B", "white", "#B63E2F")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
  column_split = pair_meta$state, top_annotation = top_pathway_annotation,
  column_title = "State-resolved pathway response matrix",
  column_title_gp = gpar(fontface = "bold", fontsize = 12)
)
pdf(file.path(out_dir, "Auto_pdo_flot_matched_pathway_heatmap.pdf"), width = 15, height = 7)
draw(pathway_heatmap, merge_legend = TRUE); dev.off()

####################
# composite response deltas
####################
message("Summarising composite pathway responses ...")
composite_delta <- bind_rows(lapply(names(composite_pathways), function(metric_name) {
  pathway_delta %>% filter(pathway %in% composite_pathways[[metric_name]]) %>%
    group_by(Patient, Patient_label, state) %>%
    summarise(metric = metric_name, delta_score = mean(delta_score, na.rm = TRUE), pathway_n = n_distinct(pathway), .groups = "drop")
})) %>% mutate(
  state = factor(state, levels = state_levels_main),
  metric = factor(metric, levels = names(composite_pathways)),
  Patient = factor(Patient, levels = patient_order),
  Patient_label = factor(Patient_label, levels = unname(patient_labels[patient_order]))
)
write.csv(composite_delta, file.path(out_dir, "Auto_pdo_flot_matched_composite_response_deltas.csv"), row.names = FALSE)
composite_summary <- composite_delta %>% group_by(state, metric) %>%
  summarise(median_delta = median(delta_score, na.rm = TRUE), .groups = "drop")

p_composite <- ggplot(composite_delta, aes(state, delta_score, color = Patient_label, group = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_line(alpha = 0.45, linewidth = 0.45) + geom_point(size = 2.4) +
  geom_point(data = composite_summary, aes(x = state, y = median_delta), inherit.aes = FALSE, color = "black", shape = 18, size = 3.0) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") + scale_color_manual(values = patient_cols) +
  labs(title = "Within-state FLOT response scores", x = NULL, y = "Composite pathway delta", color = NULL) +
  flot_theme(10) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_composite_response.pdf"), p_composite, width = 14, height = 10)

# Summary panel
summary_panel <- (p_state_bar / p_abundance_lfc / p_composite / p_hybrid_lfc / p_mp_lfc) +
  plot_layout(heights = c(1.25, 0.9, 1.1, 1.0, 1.0))
ggsave(file.path(out_dir, "Auto_pdo_flot_matched_summary_panels.pdf"), summary_panel, width = 16, height = 18)

####################
# heuristic state-fate summary
####################
message("Building state-fate summary ...")
fate_tbl <- abundance_summary %>%
  transmute(state, abundance_log2FC = median_log2FC, median_delta_pct) %>%
  left_join(composite_summary %>% select(state, metric, median_delta) %>%
    pivot_wider(names_from = metric, values_from = median_delta), by = "state") %>%
  mutate(
    `Proliferation loss` = -`Cell-cycle / proliferation`,
    Injury = `Injury / checkpoint`,
    Adaptation = `Adaptive / persistence`,
    Proteostasis = `Proteostasis / transition`,
    contributing_pairs = sapply(state, function(st) sum(pathway_delta$state == st & pathway_delta$pathway == "Apoptosis")),
    heuristic_call = mapply(classify_state_fate, abundance_log2FC, `Proliferation loss`, Injury, Adaptation, Proteostasis)
  ) %>% mutate(state = factor(state, levels = state_levels_main))
write.csv(fate_tbl, file.path(out_dir, "Auto_pdo_flot_matched_state_fate_summary.csv"), row.names = FALSE)

####################
# paired pseudobulk DEG
####################
message("Running paired pseudobulk edgeR ...")
deg_summary <- list(); deg_results <- list(); recurrent_gene_tables <- list()
pdf(file.path(out_dir, "Auto_pdo_flot_matched_recurrent_deg_heatmaps.pdf"), width = 9, height = 8)
for (state_name in state_levels_main) {
  message("  DEG: ", state_name)
  state_meta <- pseudobulk_meta %>%
    filter(state == state_name, cell_n >= min_cells_per_pseudobulk) %>%
    group_by(Patient) %>% filter(n_distinct(Treatment) == 2) %>% ungroup() %>%
    mutate(Patient = factor(as.character(Patient), levels = intersect(patient_order, unique(as.character(Patient)))),
           Treatment = factor(Treatment, levels = c("Untreated", "Treated"))) %>% arrange(Patient, Treatment)
  valid_pairs <- n_distinct(state_meta$Patient)
  state_stub <- state_key(state_name)
  if (valid_pairs < min_pairs_for_deg) {
    deg_summary[[state_name]] <- data.frame(state = state_name, valid_pairs = valid_pairs, tested_genes = 0, sig_fdr_0_05 = 0, sig_fdr_0_10 = 0, note = "Skipped", stringsAsFactors = FALSE)
    next
  }
  y <- DGEList(counts = pseudobulk_counts[, state_meta$sample_state, drop = FALSE])
  design <- model.matrix(~ Patient + Treatment, data = state_meta)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y); y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE); qlf <- glmQLFTest(fit, coef = "TreatmentTreated")
  res_tbl <- topTags(qlf, n = Inf, sort.by = "PValue")$table %>% rownames_to_column("gene") %>%
    mutate(state = state_name, valid_pairs = valid_pairs, pair_ids = paste(unique(state_meta$Patient), collapse = ";"))
  logcpm_state <- cpm(y, log = TRUE, prior.count = 2)
  patient_logfc <- sapply(as.character(unique(state_meta$Patient)), function(pid) {
    tc <- state_meta$sample_state[state_meta$Patient == pid & state_meta$Treatment == "Treated"]
    uc <- state_meta$sample_state[state_meta$Patient == pid & state_meta$Treatment == "Untreated"]
    logcpm_state[, tc] - logcpm_state[, uc]
  })
  if (is.null(dim(patient_logfc))) {
    patient_logfc <- matrix(patient_logfc, ncol = 1); rownames(patient_logfc) <- rownames(logcpm_state)
    colnames(patient_logfc) <- as.character(unique(state_meta$Patient))
  }
  res_tbl$consistent_direction_n <- rowSums(sign(patient_logfc[res_tbl$gene, , drop = FALSE]) == sign(res_tbl$logFC), na.rm = TRUE)
  deg_results[[state_name]] <- res_tbl
  write.csv(res_tbl, file.path(deg_dir, paste0("Auto_pdo_flot_matched_deg_", state_stub, ".csv")), row.names = FALSE)
  write.csv(as.data.frame(patient_logfc) %>% rownames_to_column("gene"), file.path(deg_dir, paste0("Auto_pdo_flot_matched_patient_logFC_", state_stub, ".csv")), row.names = FALSE)
  top_genes <- select_recurrent_genes(res_tbl, patient_logfc, top_deg_per_direction)
  recurrent_gene_tables[[state_name]] <- res_tbl %>% filter(gene %in% top_genes)
  write.csv(recurrent_gene_tables[[state_name]], file.path(deg_dir, paste0("Auto_pdo_flot_matched_recurrent_genes_", state_stub, ".csv")), row.names = FALSE)
  deg_summary[[state_name]] <- data.frame(state = state_name, valid_pairs = valid_pairs, tested_genes = nrow(res_tbl), sig_fdr_0_05 = sum(res_tbl$FDR < 0.05), sig_fdr_0_10 = sum(res_tbl$FDR < 0.10), note = "OK", stringsAsFactors = FALSE)
  if (length(top_genes) == 0) next
  heat_mat <- patient_logfc[top_genes, , drop = FALSE]; colnames(heat_mat) <- unname(patient_short[colnames(heat_mat)])
  dir_vec <- factor(ifelse(res_tbl$logFC[match(rownames(heat_mat), res_tbl$gene)] > 0, "Up", "Down"), levels = c("Up", "Down"))
  ht <- Heatmap(heat_mat, name = "logFC", col = colorRamp2(c(-1.5, 0, 1.5), c("#245F7B", "white", "#B63E2F")),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8),
    column_title = paste0(state_name, " recurrent DEGs"), column_title_gp = gpar(fontface = "bold", fontsize = 12))
  draw(ht, merge_legend = TRUE)
}
dev.off()
deg_summary_df <- bind_rows(deg_summary)
write.csv(deg_summary_df, file.path(out_dir, "Auto_pdo_flot_matched_deg_summary.csv"), row.names = FALSE)

####################
# cache
####################
saveRDS(list(
  cell_counts = cell_counts, abundance_pairs = abundance_pairs, mp_pairs = mp_pairs,
  hybrid_pairs_abund = hybrid_pairs_abund, pseudobulk_meta = pseudobulk_meta,
  pathway_delta = pathway_delta, composite_delta = composite_delta, fate_summary = fate_tbl,
  deg_summary = deg_summary_df, deg_results = deg_results, recurrent_genes = recurrent_gene_tables
), file = file.path(out_dir, "Auto_pdo_flot_matched_response_results.rds"))

# ============================================================
# PRESENTATION PDF (same 8 pages as before)
# ============================================================
message("=== Writing presentation PDF ===")
pdf_path <- file.path(out_dir, "Auto_pdo_flot_presentation_final.pdf")
pdf(pdf_path, width = 16, height = 9.5)

# PAGE 1 — Pathway heatmap
message("  Page 1: Pathway heatmap")
pdm <- pathway_delta %>% mutate(pair_id = paste(state, Patient, sep = "__")) %>%
  select(pathway, pair_id, delta_score) %>%
  pivot_wider(names_from = pair_id, values_from = delta_score) %>%
  column_to_rownames("pathway") %>% as.matrix()
pdm <- pdm[pathway_order_grouped, pair_meta$pair_id, drop = FALSE]
colnames(pdm) <- as.character(pair_meta$Patient)
ha_top <- HeatmapAnnotation(State = as.character(pair_meta$state), Patient = as.character(pair_meta$Patient_label),
  col = list(State = state_cols[state_levels_main], Patient = patient_cols),
  show_annotation_name = TRUE, annotation_name_gp = gpar(fontface = "bold", fontsize = 9))
row_blocks <- pathway_block_df$block[match(rownames(pdm), pathway_block_df$pathway)]
ha_row <- rowAnnotation(Block = row_blocks, col = list(Block = block_cols),
  show_annotation_name = TRUE, annotation_name_gp = gpar(fontface = "bold", fontsize = 9))
clip_val <- max(0.4, quantile(abs(pdm), 0.95, na.rm = TRUE))
ht1 <- Heatmap(pdm, name = "Treated-Untreated",
  col = colorRamp2(c(-clip_val, 0, clip_val), c("#245F7B", "white", "#B63E2F")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = factor(row_blocks, levels = names(block_cols)), row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(3, "mm"), row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
  column_split = factor(pair_meta$state, levels = state_levels_main), column_gap = unit(3, "mm"),
  top_annotation = ha_top, left_annotation = ha_row,
  column_title = "Per State Pathway Response Matrix", column_title_gp = gpar(fontface = "bold", fontsize = 14),
  cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.2f", pdm[i, j]), x, y, gp = gpar(fontsize = 6.5)))
draw(ht1, merge_legend = TRUE)

# PAGE 2 — Composite response
message("  Page 2: Composite response")
cd2 <- composite_delta %>% mutate(state = factor(state, levels = state_levels_main),
  metric = factor(metric, levels = names(composite_pathways)))
cs2 <- composite_summary %>% mutate(state = factor(state, levels = state_levels_main),
  metric = factor(metric, levels = names(composite_pathways)))
p2 <- ggplot(cd2, aes(x = state, y = delta_score, color = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.55), size = 4) +
  geom_point(data = cs2, aes(x = state, y = median_delta), inherit.aes = FALSE, shape = 18, size = 2.5, color = "black") +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") + scale_color_manual(values = patient_cols) +
  labs(title = "Per State Pathway Response Score", x = NULL, y = "Delta (Treated vs Untreated)", color = NULL) +
  pres_theme(11) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1, size = 12))
print(p2)

# PAGE 3 — Abundance change
message("  Page 3: Abundance change")
ap3 <- abundance_pairs %>% mutate(state = factor(as.character(state), levels = state_plot_order), state_x = match(as.character(state), state_plot_order))
as3 <- abundance_summary %>% mutate(state = factor(as.character(state), levels = state_plot_order), state_x = match(as.character(state), state_plot_order))
p3 <- ggplot(ap3, aes(state_x, log2FC, color = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.7) +
  geom_segment(data = as3, aes(x = state_x - 0.28, xend = state_x + 0.28, y = median_log2FC, yend = median_log2FC), inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 4) +
  geom_text(data = as3, aes(x = state_x, y = median_log2FC, label = sprintf("%.2f", median_log2FC)), inherit.aes = FALSE, vjust = -0.9, size = 4, fontface = "bold") +
  scale_color_manual(values = patient_cols) +
  scale_x_continuous(breaks = seq_along(state_plot_order), labels = state_plot_order) +
  labs(title = "State Abundance Change After FLOT", x = NULL, y = expression(log[2]~"FC"), color = NULL) +
  pres_theme(12) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", size = 12))
print(p3)

# PAGE 4 — MP expression change
message("  Page 4: MP expression change")
mp4 <- mp_pairs %>% mutate(MP = factor(as.character(MP), levels = mp_plot_order), MP_x = match(as.character(MP), mp_plot_order))
ms4 <- mp_summary %>% mutate(MP = factor(as.character(MP), levels = mp_plot_order), MP_x = match(as.character(MP), mp_plot_order))
p4 <- ggplot(mp4, aes(MP_x, log2FC, color = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.7) +
  geom_segment(data = ms4, aes(x = MP_x - 0.28, xend = MP_x + 0.28, y = median_log2FC, yend = median_log2FC), inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 4) +
  scale_x_continuous(breaks = seq_along(mp_plot_order), labels = mp_descriptions[mp_plot_order]) +
  scale_color_manual(values = patient_cols) +
  labs(title = "MP Expression Change After FLOT", x = NULL, y = "log2 FC UCell score", color = NULL) +
  pres_theme(12) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 10))
print(p4)

# PAGE 5 — Hybrid pairwise
message("  Page 5: Hybrid pairwise")
hp5 <- hybrid_pairs_abund %>% mutate(hybrid_subtype = factor(as.character(hybrid_subtype), levels = hybrid_levels_ordered), hybrid_x = match(as.character(hybrid_subtype), hybrid_levels_ordered))
hs5 <- hybrid_summary %>% mutate(hybrid_x = match(as.character(hybrid_subtype), hybrid_levels_ordered))
p5 <- ggplot(hp5, aes(hybrid_x, log2FC, color = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.7) +
  geom_segment(data = hs5, aes(x = hybrid_x - 0.28, xend = hybrid_x + 0.28, y = median_log2FC, yend = median_log2FC), inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 4) +
  scale_color_manual(values = patient_cols) +
  scale_x_continuous(breaks = seq_along(hybrid_levels_ordered), labels = hybrid_levels_ordered) +
  labs(title = "Hybrid Pairwise Abundance Change After FLOT", x = NULL, y = "log2 FC hybrid fraction", color = NULL) +
  pres_theme(12) + theme(legend.position = "bottom", axis.text.x = element_text(face = "bold", size = 10, angle = 45, hjust = 1, lineheight = 0.8))
print(p5)

# PAGE 6 — DEG heatmaps
message("  Page 6: DEG heatmaps")
build_deg_heatmap <- function(sn) {
  sk <- state_key(sn)
  f1 <- file.path(deg_dir, paste0("Auto_pdo_flot_matched_recurrent_genes_", sk, ".csv"))
  f2 <- file.path(deg_dir, paste0("Auto_pdo_flot_matched_patient_logFC_", sk, ".csv"))
  if (!file.exists(f1) || !file.exists(f2)) return(NULL)
  recur <- read.csv(f1); plfc <- read.csv(f2) %>% column_to_rownames("gene") %>% as.matrix()
  tg <- recur %>% filter(gene %in% rownames(plfc)) %>% arrange(desc(logFC)) %>% pull(gene)
  if (length(tg) == 0) return(NULL)
  hm <- plfc[tg, , drop = FALSE]; colnames(hm) <- paste0("SUR", unname(patient_short[colnames(hm)]))
  Heatmap(hm, name = "logFC", col = colorRamp2(c(-1.5, 0, 1.5), c("#245F7B", "white", "#B63E2F")),
    cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8),
    column_title = sn, column_title_gp = gpar(fontface = "bold", fontsize = 11),
    show_heatmap_legend = FALSE, width = unit(ncol(hm)*8, "mm"), height = unit(nrow(hm)*4, "mm"))
}
ht_list <- lapply(state_levels_main, build_deg_heatmap)
grid.newpage(); grid.text("Recurrent Paired Pseudobulk DEGs", y = unit(0.97, "npc"), gp = gpar(fontface = "bold", fontsize = 14))
for (i in seq_along(ht_list)) {
  if (is.null(ht_list[[i]])) next
  ri <- (i - 1) %/% 3; ci <- (i - 1) %% 3
  pushViewport(viewport(x = ci/3, y = 0.93 - (ri+1)*0.45, width = 1/3, height = 0.45, just = c("left", "bottom")))
  draw(ht_list[[i]], newpage = FALSE); popViewport()
}

# PAGE 7 — State fate summary
message("  Page 7: Fate summary")
fate_mat_df <- composite_summary %>% pivot_wider(names_from = metric, values_from = median_delta) %>%
  rename(Proliferation = `Cell-cycle / proliferation`, Injury = `Injury / checkpoint`,
         Adaptation = `Adaptive / persistence`, Proteostasis = `Proteostasis / transition`) %>%
  left_join(abundance_summary %>% select(state, abundance_log2FC = median_log2FC), by = "state") %>%
  column_to_rownames("state")
abund_col <- as.matrix(fate_mat_df[, "abundance_log2FC", drop = FALSE]); colnames(abund_col) <- "Abundance
log2FC"
path_mat <- as.matrix(fate_mat_df[, c("Proliferation", "Injury", "Adaptation", "Proteostasis")])
row_order_fate <- rev(state_levels_main)
abund_col <- abund_col[row_order_fate, , drop = FALSE]; path_mat <- path_mat[row_order_fate, , drop = FALSE]
ht_a <- Heatmap(abund_col, name = "Abund", col = colorRamp2(c(-max(abs(abund_col)), 0, max(abs(abund_col))), c("#1B4965", "white", "#9B2226")),
  cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  show_column_names = FALSE, width = unit(2.5, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.2f", abund_col[i, j]), x, y, gp = gpar(fontsize = 14)),
  column_title = "Abundance", column_title_gp = gpar(fontface = "bold", fontsize = 14))
ht_p <- Heatmap(path_mat, name = "Delta", col = colorRamp2(c(-max(abs(path_mat)), 0, max(abs(path_mat))), c("#245F7B", "white", "#B63E2F")),
  cluster_rows = FALSE, cluster_columns = FALSE, column_names_gp = gpar(fontsize = 16, fontface = "bold"), column_names_rot = 35,
  width = unit(12, "cm"), show_row_names = FALSE,
  cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.2f", path_mat[i, j]), x, y, gp = gpar(fontsize = 14)),
  column_title = "Pathway composite deltas", column_title_gp = gpar(fontface = "bold", fontsize = 14))
draw(ht_a + ht_p, column_title = "Per State Summary", column_title_gp = gpar(fontface = "bold", fontsize = 18),
     merge_legend = FALSE, show_heatmap_legend = FALSE)

# PAGE 8 — UMAP support
message("  Page 8: UMAP support")
umap_base <- function(df, ...) {
  ggplot(df, aes(UMAP_1, UMAP_2, ...)) + geom_point(size = 0.18, alpha = 0.8) + coord_equal() + pres_theme(9) + theme(legend.position = "right")
}
pu1 <- umap_base(umap_plot_df, color = state) + scale_color_manual(values = state_cols, drop = FALSE) + labs(title = "Finalized states", color = NULL) +
  theme(legend.text = element_text(size = 14)) + guides(color = guide_legend(override.aes = list(size = 6)))
pu2 <- umap_base(umap_plot_df, color = Treatment) + scale_color_manual(values = treatment_cols) + labs(title = "Untreated vs FLOT", color = NULL) +
  theme(legend.text = element_text(size = 14)) + guides(color = guide_legend(override.aes = list(size = 6)))
pu3 <- umap_base(umap_plot_df, color = injury_score) + scale_color_gradientn(colors = c("#1B365D", "#4F9DD9", "#F2E6B6", "#B33C2F"), na.value = "grey85") +
  labs(title = "Injury score", color = NULL)
pu4 <- umap_base(umap_plot_df, color = adaptive_score) + scale_color_gradientn(colors = c("#1B365D", "#4F9DD9", "#F2E6B6", "#B33C2F"), na.value = "grey85") +
  labs(title = "Adaptation score", color = NULL)
print((pu1 | pu2) / (pu3 | pu4) + plot_annotation(title = "PDO Matched-Pair Support UMAPs",
  theme = theme(plot.title = element_text(size = 16, face = "bold"))))
dev.off()
message("Presentation PDF done: ", pdf_path)

####################
# NEW PDF 1: Node plots (5 pages)
####################
message("=== NEW PDF 1: Node plots ===")

real_states <- c("Classic Proliferative","Basal to Intest. Meta","SMG-like Metaplasia","Stress-adaptive")
n_rs <- length(real_states)
theta <- seq(0, 2*pi, length.out = n_rs + 1)[1:n_rs]
layout_df <- data.frame(state = real_states, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)

# Build per-sample node data
build_node_data <- function(patient_id, trt) {
  samp <- paste0(patient_id, "_", trt, "_PDO")
  cells_in <- names(hybrid_state_vec)[hybrid_state_vec %in% c(real_states, "Hybrid")]
  cells_in <- intersect(cells_in, Cells(pdos_matched)[pdos_matched$orig.ident == samp])
  if (length(cells_in) == 0) return(list(node_df = NULL, edge_df = NULL))
  st <- hybrid_state_vec[cells_in]
  tot <- length(cells_in)
  # State nodes
  sdf <- data.frame(state = st[st %in% real_states], stringsAsFactors = FALSE) %>%
    count(state, name = "cells") %>% mutate(pct = 100 * cells / tot)
  # Hybrid edges
  hyb_cells <- intersect(names(st)[st == "Hybrid"], rownames(mp_adj_noncc))
  if (length(hyb_cells) > 0) {
    gm <- do.call(cbind, lapply(state_groups_hybrid, function(mps) {
      ma <- intersect(mps, colnames(mp_adj_noncc))
      if (length(ma) == 1) return(as.numeric(mp_adj_noncc[hyb_cells, ma]))
      if (length(ma) == 0) return(rep(0, length(hyb_cells)))
      apply(mp_adj_noncc[hyb_cells, ma, drop = FALSE], 1, max)
    }))
    colnames(gm) <- names(state_groups_hybrid); rownames(gm) <- hyb_cells
    pair_lab <- vapply(hyb_cells, function(cl) {
      t2 <- names(sort(gm[cl,], decreasing = TRUE))[1:2]
      paste(sort(t2), collapse = "__")
    }, character(1))
    edf <- data.frame(pair = pair_lab, stringsAsFactors = FALSE) %>%
      count(pair, name = "hybrid_cells") %>%
      tidyr::separate(pair, into = c("from","to"), sep = "__", remove = FALSE) %>%
      mutate(pct = 100 * hybrid_cells / tot)
  } else {
    edf <- data.frame(pair = character(), from = character(), to = character(), hybrid_cells = integer(), pct = numeric())
  }
  ndf <- left_join(layout_df, sdf, by = "state")
  ndf$cells[is.na(ndf$cells)] <- 0; ndf$pct[is.na(ndf$pct)] <- 0
  ndf$label_x <- ndf$x * 1.25; ndf$label_y <- ndf$y * 1.25
  edf2 <- edf %>% left_join(layout_df, by = c("from" = "state")) %>%
    left_join(layout_df, by = c("to" = "state"), suffix = c("","_to")) %>%
    rename(xend = x_to, yend = y_to)
  list(node_df = ndf, edge_df = edf2)
}

draw_nodeplot <- function(ndf, edf, title_text, max_node, max_edge) {
  ggplot() +
    { if (nrow(edf) > 0) geom_segment(data = edf, aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct), color = "grey35", alpha = 0.8) } +
    geom_point(data = ndf, aes(x = x, y = y, size = pct, color = state)) +
    geom_text(data = ndf, aes(x = label_x, y = label_y, label = paste0(state, "
", sprintf("%.1f%%", pct))), size = 3.5, fontface = "bold") +
    { if (nrow(edf) > 0) geom_label(data = edf, aes(x = (x+xend)/2, y = (y+yend)/2, label = sprintf("%.1f%%", pct)), size = 2.6, fill = "white", label.size = 0, fontface = "bold") } +
    scale_color_manual(values = state_cols) +
    scale_size(limits = c(0, max_node), range = c(8, 30), guide = "none") +
    scale_linewidth(limits = c(0, max_edge), range = c(0.6, 12), guide = "none") +
    coord_equal() + expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
    theme_void(base_size = 14) +
    labs(title = title_text) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
}

# Collect per-patient per-treatment data
all_node_data <- list()
for (pid in patient_order) {
  for (trt in c("Untreated", "Treated")) {
    key <- paste(pid, trt, sep = "_")
    all_node_data[[key]] <- build_node_data(pid, trt)
  }
}

# Combined: median node pct and edge pct
combine_node_edge <- function(trt_label) {
  keys <- paste(patient_order, trt_label, sep = "_")
  node_list <- lapply(keys, function(k) all_node_data[[k]]$node_df)
  node_all <- bind_rows(node_list) %>% group_by(state, x, y, label_x, label_y) %>%
    summarise(pct = median(pct, na.rm = TRUE), cells = median(cells, na.rm = TRUE), .groups = "drop")
  edge_list <- lapply(keys, function(k) all_node_data[[k]]$edge_df)
  edge_all <- bind_rows(edge_list)
  if (nrow(edge_all) > 0) {
    edge_all <- edge_all %>% group_by(from, to, x, y, xend, yend) %>%
      summarise(pct = median(pct, na.rm = TRUE), hybrid_cells = median(hybrid_cells, na.rm = TRUE), .groups = "drop")
  }
  list(node_df = node_all, edge_df = edge_all)
}

node_pdf_path <- file.path(out_dir, "Auto_pdo_flot_nodeplot_untreated_vs_treated.pdf")
message("Writing: ", node_pdf_path)
pdf(node_pdf_path, width = 14, height = 7)

# Calculate global max for scale matching
comb_ut <- combine_node_edge("Untreated"); comb_tr <- combine_node_edge("Treated")
max_node_comb <- max(c(comb_ut$node_df$pct, comb_tr$node_df$pct), na.rm = TRUE)
max_edge_comb <- max(c(comb_ut$edge_df$pct, comb_tr$edge_df$pct, 0.1), na.rm = TRUE)

# Page 1: Combined
p_comb <- draw_nodeplot(comb_ut$node_df, comb_ut$edge_df, "Untreated (median)", max_node_comb, max_edge_comb) |
  draw_nodeplot(comb_tr$node_df, comb_tr$edge_df, "Treated (median)", max_node_comb, max_edge_comb)
print(p_comb + plot_annotation(title = "Combined Across 4 Patients", theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))

# Pages 2-5: Per patient
for (pid in patient_order) {
  ut <- all_node_data[[paste0(pid, "_Untreated")]]; tr <- all_node_data[[paste0(pid, "_Treated")]]
  if (is.null(ut$node_df) || is.null(tr$node_df)) next
  max_node_pt <- max(c(ut$node_df$pct, tr$node_df$pct), na.rm = TRUE)
  max_edge_pt <- max(c(ut$edge_df$pct, tr$edge_df$pct, 0.1), na.rm = TRUE)
  
  pp <- draw_nodeplot(ut$node_df, ut$edge_df, "Untreated", max_node_pt, max_edge_pt) |
    draw_nodeplot(tr$node_df, tr$edge_df, "FLOT-treated", max_node_pt, max_edge_pt)
  print(pp + plot_annotation(title = paste0(pid, " — ", patient_labels[pid]),
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))
}
dev.off()
message("Node plot PDF done.")

####################
# NEW PDF 2: Paired boxplots (3 pages)
####################
message("=== NEW PDF 2: Paired boxplots ===")

sig_label <- function(p) {
  if (is.na(p)) return("")
  # Do not use multiple testing correction, report raw pairs as requested
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

# State abundance per patient per treatment (raw pct)
state_abund_long <- cell_counts %>%
  filter(state %in% state_levels_main) %>%
  mutate(state = factor(as.character(state), levels = state_levels_main),
         Treatment = factor(as.character(Treatment), levels = c("Untreated", "Treated")))

state_sig <- state_abund_long %>%
  group_by(state) %>%
  summarise(
    p = tryCatch(wilcox.test(pct[Treatment == "Untreated"], pct[Treatment == "Treated"], paired = TRUE)$p.value, error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(label = sapply(p, function(pval) {
    if (is.na(pval)) return("")
    paste0("p=", format(round(pval, 3), nsmall=3))
  }))

state_ymax <- state_abund_long %>% group_by(state) %>% summarise(ymax = max(pct, na.rm = TRUE), .groups = "drop")
state_sig <- left_join(state_sig, state_ymax, by = "state") %>% mutate(y_pos = ymax * 1.08)

# Hybrid abundance per patient per treatment
hybrid_abund_long <- hybrid_counts %>%
  mutate(hybrid_subtype = factor(as.character(hybrid_subtype), levels = hybrid_levels_ordered),
         Treatment = factor(as.character(Treatment), levels = c("Untreated", "Treated")))

hybrid_sig <- hybrid_abund_long %>%
  group_by(hybrid_subtype) %>%
  summarise(
    p = tryCatch(wilcox.test(pct[Treatment == "Untreated"], pct[Treatment == "Treated"], paired = TRUE)$p.value, error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(label = sapply(p, function(pval) {
    if (is.na(pval)) return("")
    paste0("p=", format(round(pval, 3), nsmall=3))
  }))
hybrid_ymax <- hybrid_abund_long %>% group_by(hybrid_subtype) %>% summarise(ymax = max(pct, na.rm = TRUE), .groups = "drop")
hybrid_sig <- left_join(hybrid_sig, hybrid_ymax, by = "hybrid_subtype") %>% mutate(y_pos = ymax * 1.08)

# MP expression per patient per treatment (mean UCell per sample)
mp_expr_long <- mp_sample_scores %>%
  filter(MP %in% mp_plot_order) %>%
  mutate(MP = factor(MP, levels = mp_plot_order),
         Treatment = factor(as.character(Treatment), levels = c("Untreated", "Treated")))

mp_sig <- mp_expr_long %>%
  group_by(MP) %>%
  summarise(
    p = tryCatch(wilcox.test(mean_score[Treatment == "Untreated"], mean_score[Treatment == "Treated"], paired = TRUE)$p.value, error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(label = sapply(p, function(pval) {
    if (is.na(pval)) return("")
    paste0("p=", format(round(pval, 3), nsmall=3))
  }))
mp_ymax <- mp_expr_long %>% group_by(MP) %>% summarise(ymax = max(mean_score, na.rm = TRUE), .groups = "drop")
mp_sig <- left_join(mp_sig, mp_ymax, by = "MP") %>% mutate(y_pos = ymax * 1.08)

box_pdf_path <- file.path(out_dir, "Auto_pdo_flot_paired_boxplots.pdf")
message("Writing: ", box_pdf_path)
pdf(box_pdf_path, width = 16, height = 9)

# Page 1: State abundance
p_box_state <- ggplot(state_abund_long, aes(x = state, y = pct, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_point(aes(color = Patient_label, group = Treatment), position = position_dodge(width = 0.75), size = 2.5, alpha = 0.9) +
  geom_text(data = state_sig %>% filter(label != ""), aes(x = state, y = y_pos, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = treatment_cols) +
  scale_color_manual(values = patient_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "State Abundance: Untreated vs FLOT-treated", x = NULL, y = "% of malignant cells", fill = NULL, color = "Patient") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 35, hjust = 1, size = 12, face = "bold"),
        legend.position = "top", legend.text = element_text(size = 12))
print(p_box_state)

# Page 2: Hybrid abundance
p_box_hybrid <- ggplot(hybrid_abund_long, aes(x = hybrid_subtype, y = pct, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_point(aes(color = Patient_label, group = Treatment), position = position_dodge(width = 0.75), size = 2.5, alpha = 0.9) +
  geom_text(data = hybrid_sig %>% filter(label != ""), aes(x = hybrid_subtype, y = y_pos, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = treatment_cols) +
  scale_color_manual(values = patient_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "Hybrid Abundance: Untreated vs FLOT-treated", x = NULL, y = "% of hybrid cells", fill = NULL, color = "Patient") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", lineheight = 0.8),
        legend.position = "top", legend.text = element_text(size = 12))
print(p_box_hybrid)

# Page 3: MP expression
mp_axis_labels_box <- setNames(ifelse(!is.na(mp_descriptions[mp_plot_order]), paste0(mp_plot_order, "
", sub("^MP[0-9]+_", "", mp_descriptions[mp_plot_order])), mp_plot_order), mp_plot_order)
p_box_mp <- ggplot(mp_expr_long, aes(x = MP, y = mean_score, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_point(aes(color = Patient_label, group = Treatment), position = position_dodge(width = 0.75), size = 2.5, alpha = 0.9) +
  geom_text(data = mp_sig %>% filter(label != ""), aes(x = MP, y = y_pos, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = treatment_cols) +
  scale_color_manual(values = patient_cols) +
  scale_x_discrete(labels = mp_axis_labels_box) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "MP Expression: Untreated vs FLOT-treated", x = NULL, y = "Mean UCell score per sample", fill = NULL, color = "Patient") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        legend.position = "top", legend.text = element_text(size = 12))
print(p_box_mp)
dev.off()
message("Paired boxplot PDF done.")

####################
# NEW PDF 3: Improved pathway heatmap
####################
message("=== NEW PDF 3: Improved pathway heatmap ===")

# --- CC signature score ---
cc_genes_df <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE, stringsAsFactors = FALSE
)[, 1:3]
cc_consensus <- cc_genes_df$Gene[cc_genes_df$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(pdos_matched))
cc_top50 <- names(sort(rowMeans(as.matrix(
  GetAssayData(pdos_matched, assay = "RNA", layer = "data")[cc_consensus, , drop = FALSE]
), na.rm = TRUE), decreasing = TRUE))[1:50]

# Score CC per pseudobulk sample using logCPM z-scored matrix
cc_pb_score <- mean_score(pb_logcpm_z, cc_top50)
names(cc_pb_score) <- colnames(pb_logcpm_z)

# Score MP4 and MP8 per pseudobulk sample
mp_genes_all <- geneNMF.metaprograms$metaprograms.genes
mp4_genes <- mp_genes_all[["MP4"]]
mp8_genes <- mp_genes_all[["MP8"]]
mp4_pb_score <- mean_score(pb_logcpm_z, mp4_genes)
mp8_pb_score <- mean_score(pb_logcpm_z, mp8_genes)
names(mp4_pb_score) <- colnames(pb_logcpm_z)
names(mp8_pb_score) <- colnames(pb_logcpm_z)

# Add to extra_scores
extra_scores <- data.frame(
  sample_state = colnames(pb_logcpm_z),
  CCSIG = cc_pb_score,
  `Intestinal Metaplasia` = mp4_pb_score,
  `Columnar Progenitor` = mp8_pb_score,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  left_join(pseudobulk_meta %>% select(sample_state, state, Patient, Patient_label, Treatment, cell_n), by = "sample_state") %>%
  filter(cell_n >= min_cells_per_pseudobulk)

extra_delta <- extra_scores %>%
  select(Patient, Patient_label, state, Treatment, CCSIG, `Intestinal Metaplasia`, `Columnar Progenitor`, cell_n) %>%
  pivot_longer(cols = c(CCSIG, `Intestinal Metaplasia`, `Columnar Progenitor`), names_to = "pathway", values_to = "score") %>%
  pivot_wider(names_from = Treatment, values_from = c(score, cell_n)) %>%
  filter(!is.na(score_Untreated), !is.na(score_Treated)) %>%
  mutate(delta_score = score_Treated - score_Untreated) %>%
  select(Patient, Patient_label, state, pathway, delta_score)

# Combine with pathway_delta
comb_delta <- bind_rows(
  pathway_delta %>% select(Patient, Patient_label, state, pathway, delta_score),
  extra_delta
)

# Mean delta
comb_mean <- comb_delta %>%
  group_by(state, pathway) %>%
  summarise(mean_delta = mean(delta_score, na.rm = TRUE), .groups = "drop")

# Pairwise sig (paired test of Treated - Untreated == 0)
comb_sig <- comb_delta %>%
  group_by(state, pathway) %>%
  summarise(
    p = tryCatch(wilcox.test(delta_score, mu = 0)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = sapply(p, sig_label))

# Build matrices (rows = features, cols = states)
states_for_hm <- state_levels_main

pw_mat <- comb_mean %>%
  pivot_wider(names_from = state, values_from = mean_delta) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

pw_order <- c(pathway_order_grouped, "CCSIG", "Intestinal Metaplasia", "Columnar Progenitor")
pw_mat <- pw_mat[pw_order, states_for_hm, drop = FALSE]

pw_sig_mat <- comb_sig %>%
  select(state, pathway, label) %>%
  pivot_wider(names_from = state, values_from = label) %>%
  column_to_rownames("pathway") %>% as.matrix()
pw_sig_mat <- pw_sig_mat[pw_order, states_for_hm, drop = FALSE]

row_blocks_pw <- c(pathway_block_df$block[match(pathway_order_grouped, pathway_block_df$pathway)], "CCSIG", "Lineage", "Lineage")
block_cols_ext <- c(block_cols, "CCSIG" = "#6A0572", "Lineage" = "#1A535C")

# Color scale (unified for all)
pw_clip <- max(0.3, quantile(abs(pw_mat), 0.95, na.rm = TRUE))
col_pw <- colorRamp2(c(-pw_clip, 0, pw_clip), c("#245F7B", "white", "#B63E2F"))

# Row block annotation for pathways + CCSIG + Lineage
ha_row_pw <- rowAnnotation(Block = row_blocks_pw, col = list(Block = block_cols_ext),
  show_annotation_name = FALSE, show_legend = TRUE,
  annotation_legend_param = list(Block = list(title = "Signature Type")))

# State annotation on top
ha_top_imp <- HeatmapAnnotation(
  State = states_for_hm,
  col = list(State = state_cols[states_for_hm]),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 10)
)

# Build heatmaps
cell_fun_sig <- function(j, i, x, y, w, h, fill) {
  grid.text(sprintf("%.2f", pw_mat[i, j]), x, y - unit(1, "mm"), gp = gpar(fontsize = 7))
  lbl <- pw_sig_mat[i, j]
  if (!is.na(lbl) && lbl != "" && lbl != "ns") {
    grid.text(lbl, x, y + unit(2, "mm"), gp = gpar(fontsize = 9, fontface = "bold", col = "black"))
  }
}

ht_pw <- Heatmap(pw_mat, name = "Mean Delta (Treated - Untreated)",
  col = col_pw, cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  show_column_names = FALSE,
  row_split = factor(row_blocks_pw, levels = names(block_cols_ext)),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(3, "mm"),
  left_annotation = ha_row_pw,
  top_annotation = ha_top_imp,
  cell_fun = cell_fun_sig
)

improved_hm_path <- file.path(out_dir, "Auto_pdo_flot_improved_pathway_heatmap.pdf")
message("Writing: ", improved_hm_path)
pdf(improved_hm_path, width = 12, height = 11)
draw(ht_pw,
  column_title = "Mean Pathway Response Matrix (Treated - Untreated)",
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  merge_legend = TRUE
)
dev.off()
message("Improved pathway heatmap PDF done.")

####################
# notes
####################
writeLines(c(
  "# PDO matched FLOT presentation notes",
  "",
  "- Merged script: Auto_pdo_flot_matched_response.R",
  "- All original outputs preserved.",
  "- NEW: Auto_pdo_flot_nodeplot_untreated_vs_treated.pdf (5 pages: combined + 4 per-patient).",
  "- NEW: Auto_pdo_flot_paired_boxplots.pdf (3 pages: state abundance, hybrid abundance, MP expression).",
  "- NEW: Auto_pdo_flot_improved_pathway_heatmap.pdf (mean across patients, CCSIG + lineage + significance, unified scale)."
), con = file.path(out_dir, "Auto_pdo_flot_matched_presentation_notes.md"))

message("=== All outputs complete. ===")
