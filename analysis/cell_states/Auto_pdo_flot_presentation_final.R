####################
# Auto_pdo_flot_presentation_final.R
# Polished 6-page "Presentation Final" PDF for the matched FLOT response analysis.
# Pages: 1-Pathway heatmap, 2-Composite response, 3-Abundance change,
#         4-Recurrent DEG heatmaps (2 states), 5-State fate summary, 6-UMAP support.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(grid)
})

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
out_dir <- "Auto_pdo_flot_matched_response"

####################
# shared constants (mirrored from original script)
####################
patient_order <- c("SUR1070", "SUR1090", "SUR1072", "SUR1181")
patient_labels <- c(
  SUR1070 = "SUR1070 (Pre-responder)",
  SUR1090 = "SUR1090 (Post-responder)",
  SUR1072 = "SUR1072 (Post-nonresponder)",
  SUR1181 = "SUR1181 (Pre-nonresponder)"
)
patient_short <- c(SUR1070 = "1070", SUR1090 = "1090", SUR1072 = "1072", SUR1181 = "1181")

state_levels_main <- c(
  "Classic Proliferative", "Basal to Intest. Meta",
  "Stress-adaptive", "SMG-like Metaplasia",
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

# pathway ordering grouped by functional block
pathway_block_df <- data.frame(
  pathway = c(
    "E2F targets", "G2M checkpoint",
    "Apoptosis", "p53 pathway", "DNA repair",
    "TNFa/NFkB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response",
    "Unfolded protein response", "Oxidative phosphorylation"
  ),
  block = c(
    rep("Cell-cycle", 2),
    rep("Injury", 3),
    rep("Adaptive", 5),
    rep("Proteostasis", 2)
  ),
  stringsAsFactors = FALSE
)
pathway_order_grouped <- pathway_block_df$pathway

block_cols <- c(
  "Cell-cycle"    = "#F0C75E",
  "Injury"        = "#E07B54",
  "Adaptive"      = "#5DAA68",
  "Proteostasis"  = "#5B9BD5"
)

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

####################
# load cached results
####################
message("Loading cached results ...")
cache <- readRDS(file.path(out_dir, "Auto_pdo_flot_matched_response_results.rds"))

pathway_delta     <- cache$pathway_delta
composite_delta   <- cache$composite_delta
abundance_pairs   <- cache$abundance_pairs
fate_tbl          <- cache$fate_summary
pseudobulk_meta   <- cache$pseudobulk_meta

# rebuild small derived objects
abundance_summary <- abundance_pairs %>%
  group_by(state) %>%
  summarise(
    median_log2FC   = median(log2FC, na.rm = TRUE),
    median_delta_pct = median(delta_pct, na.rm = TRUE),
    .groups = "drop"
  )

composite_summary <- composite_delta %>%
  group_by(state, metric) %>%
  summarise(median_delta = median(delta_score, na.rm = TRUE), .groups = "drop")

####################
# open PDF
####################
pdf_path <- file.path(out_dir, "Auto_pdo_flot_presentation_final.pdf")
message("Writing: ", pdf_path)
pdf(pdf_path, width = 16, height = 9.5)

# ================================================================
# PAGE 1 — Pathway heatmap (ComplexHeatmap)
# ================================================================
message("  Page 1: Pathway heatmap ...")

# rebuild delta matrix in grouped pathway order
pair_meta <- pathway_delta %>%
  distinct(pair_id = paste(state, Patient, sep = "__"), state, Patient, Patient_label) %>%
  mutate(
    state   = factor(state, levels = state_levels_main),
    Patient = factor(Patient, levels = patient_order)
  ) %>%
  arrange(state, Patient)

pdm <- pathway_delta %>%
  mutate(pair_id = paste(state, Patient, sep = "__")) %>%
  select(pathway, pair_id, delta_score) %>%
  pivot_wider(names_from = pair_id, values_from = delta_score) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

pdm <- pdm[pathway_order_grouped, pair_meta$pair_id, drop = FALSE]
colnames(pdm) <- as.character(pair_meta$Patient)

# state and patient annotations
ha_top <- HeatmapAnnotation(
  State   = as.character(pair_meta$state),
  Patient = as.character(pair_meta$Patient_label),
  col = list(
    State   = state_cols[state_levels_main],
    Patient = patient_cols
  ),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 9)
)

# row (pathway) block annotation
row_blocks <- pathway_block_df$block[match(rownames(pdm), pathway_block_df$pathway)]
ha_row <- rowAnnotation(
  Block = row_blocks,
  col   = list(Block = block_cols),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 9)
)

# determine colour threshold — clip at 95th percentile of abs values (min 0.4)
clip_val <- max(0.4, quantile(abs(pdm), 0.95, na.rm = TRUE))
clip_val <- round(clip_val, 2)

ht_pathway <- Heatmap(
  pdm,
  name = "Treated − Untreated",
  col  = colorRamp2(c(-clip_val, 0, clip_val), c("#245F7B", "white", "#B63E2F")),
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  row_split        = factor(row_blocks, levels = names(block_cols)),
  row_title_gp     = gpar(fontsize = 9, fontface = "bold"),
  row_gap          = unit(3, "mm"),
  row_names_gp     = gpar(fontsize = 9),
  column_names_gp  = gpar(fontsize = 9),
  column_split     = factor(pair_meta$state, levels = state_levels_main),
  column_gap       = unit(3, "mm"),
  top_annotation   = ha_top,
  left_annotation  = ha_row,
  column_title     = "Per State Pathway Response Matrix",
  column_title_gp  = gpar(fontface = "bold", fontsize = 14),
  heatmap_legend_param = list(title_position = "topcenter"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", pdm[i, j]), x, y, gp = gpar(fontsize = 6.5))
  }
)

draw(ht_pathway, merge_legend = TRUE)

# ================================================================
# PAGE 2 — Composite response (ggplot)
# ================================================================
message("  Page 2: Composite response ...")

composite_delta2 <- composite_delta %>%
  mutate(
    state  = factor(state, levels = state_levels_main),
    metric = factor(metric, levels = c(
      "Cell-cycle / proliferation", "Injury / checkpoint",
      "Adaptive / persistence", "Proteostasis / transition"
    ))
  )
composite_summary2 <- composite_summary %>%
  mutate(
    state  = factor(state, levels = state_levels_main),
    metric = factor(metric, levels = levels(composite_delta2$metric))
  )

p2 <- ggplot(composite_delta2, aes(x = state, y = delta_score, color = Patient_label)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#B63E2F", alpha = 0.04) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#245F7B", alpha = 0.04) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.55), size = 4) +
  geom_point(
    data = composite_summary2,
    aes(x = state, y = median_delta),
    inherit.aes = FALSE,
    shape = 18, size = 2.5, color = "black"
  ) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  scale_color_manual(values = patient_cols) +
  labs(
    title = "Per State Pathway Response Score",
    x     = NULL,
    y     = "PathwayScoreDelta (Treated vs Untreated)",
    color = NULL
  ) +
  pres_theme(11) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 12),
    strip.text      = element_text(size = 13),
    panel.spacing.y = unit(3, "lines")
  )

print(p2)

# ================================================================
# PAGE 3 — State abundance change (ggplot)
# ================================================================
message("  Page 3: Abundance change ...")

abund_plot <- abundance_pairs %>%
  mutate(state = factor(state, levels = rev(state_levels_main)))
abund_sum_plot <- abundance_summary %>%
  mutate(state = factor(state, levels = rev(state_levels_main)))

p3 <- ggplot(abund_plot, aes(x = factor(state, levels = state_levels_main), y = log2FC, color = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.55), size = 6) +
  geom_point(
    data = abund_sum_plot,
    aes(x = factor(state, levels = state_levels_main), y = median_log2FC),
    inherit.aes = FALSE,
    shape = 18, size = 4, color = "black"
  ) +
  scale_color_manual(values = patient_cols) +
  labs(
    title = "State Abundance Change After FLOT",
    x     = NULL,
    y     = expression(log[2]~"fold-change (FLOT-treated / untreated)"),
    color = NULL
  ) +
  pres_theme(12) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(face = "bold", size = 12, angle = 35, hjust = 1)
  )

print(p3)

# ================================================================
# PAGE 4 — Recurrent DEG heatmaps (All Malignant States)
# ================================================================
message("  Page 4: Recurrent DEG heatmaps ...")

deg_dir <- file.path(out_dir, "pseudobulk_deg")

build_deg_heatmap <- function(state_name, show_legend = FALSE) {
  sk <- state_key(state_name)
  recur <- read.csv(file.path(deg_dir, paste0("Auto_pdo_flot_matched_recurrent_genes_", sk, ".csv")))
  plfc  <- read.csv(file.path(deg_dir, paste0("Auto_pdo_flot_matched_patient_logFC_", sk, ".csv")))

  plfc_mat <- plfc %>% column_to_rownames("gene") %>% as.matrix()

  # Order genes by average logFC (descending) so Up-regulated genes are at top
  top_genes <- recur %>%
    filter(gene %in% rownames(plfc_mat)) %>%
    arrange(desc(logFC)) %>%
    pull(gene)

  if (length(top_genes) == 0) return(NULL)

  heat_mat <- plfc_mat[top_genes, , drop = FALSE]
  colnames(heat_mat) <- paste0("SUR", unname(patient_short[colnames(heat_mat)]))

  direction_vec <- ifelse(recur$logFC[match(rownames(heat_mat), recur$gene)] > 0,
                          "Up in treated", "Down in treated")
  direction_vec <- factor(direction_vec, levels = c("Up in treated", "Down in treated"))

  row_ha <- rowAnnotation(
    Direction = direction_vec,
    col = list(Direction = c("Up in treated" = "#B63E2F", "Down in treated" = "#245F7B")),
    show_annotation_name = show_legend,
    show_legend = show_legend
  )

  col_labels <- colnames(heat_mat)
  
  # Map colors from patient_cols to SUR prefixed IDs
  patient_ids_from_cols <- str_extract(names(patient_cols), "SUR\\d+")
  patient_cols_map <- patient_cols
  names(patient_cols_map) <- patient_ids_from_cols
  
  top_ha <- HeatmapAnnotation(
    Patient = col_labels,
    col = list(Patient = patient_cols_map),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  Heatmap(
    heat_mat,
    name = "Per-patient\nlogFC",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#245F7B", "white", "#B63E2F")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    top_annotation = top_ha,
    right_annotation = row_ha,
    column_title = state_name,
    column_title_gp = gpar(fontface = "bold", fontsize = 11),
    show_heatmap_legend = FALSE,
    width = unit(ncol(heat_mat)*8, "mm"),
    height = unit(nrow(heat_mat)*4, "mm")
  )
}

ht_list <- lapply(state_levels_main, build_deg_heatmap)
names(ht_list) <- state_levels_main

# Grid layout for 5 heatmaps (2 rows, 3 columns)
grid.newpage()
grid.text("Recurrent Paired Pseudobulk DEGs", y = unit(0.97, "npc"),
          gp = gpar(fontface = "bold", fontsize = 14))

v_width  <- 1/3
v_height <- 0.45 # 2 rows fit in ~0.9npc

for (i in seq_along(ht_list)) {
  if (is.null(ht_list[[i]])) next
  
  row_idx <- (i - 1) %/% 3 # 0 for top row, 1 for bottom row
  col_idx <- (i - 1) %% 3
  
  # y starts from bottom. Row 0 (top) is at y=0.48. Row 1 (bottom) is at y=0.03.
  y_start <- 0.93 - (row_idx + 1) * v_height
  x_start <- col_idx * v_width
  
  pushViewport(viewport(x = x_start, y = y_start, 
                        width = v_width, height = v_height,
                        just = c("left", "bottom")))
  draw(ht_list[[i]], newpage = FALSE)
  popViewport()
}

# ================================================================
# PAGE 5 — State fate summary heatmap (ComplexHeatmap)
# ================================================================
message("  Page 5: State fate summary ...")

# Re-derive fate matrix directly from summaries to ensure absolute consistency with Pages 2 and 3
fate_mat_df <- composite_summary %>%
  pivot_wider(names_from = metric, values_from = median_delta) %>%
  rename(
    `Proliferation` = `Cell-cycle / proliferation`,
    `Injury`        = `Injury / checkpoint`,
    `Adaptation`    = `Adaptive / persistence`,
    `Proteostasis`  = `Proteostasis / transition`
  ) %>%
  left_join(abundance_summary %>% select(state, abundance_log2FC = median_log2FC), by = "state") %>%
  column_to_rownames("state")

# separate abundance column from pathway columns
abund_col <- as.matrix(fate_mat_df[, "abundance_log2FC", drop = FALSE])
colnames(abund_col) <- "Abundance\nlog2FC"
path_mat  <- as.matrix(fate_mat_df[, c("Proliferation", "Injury", "Adaptation", "Proteostasis")])

# reorder rows
row_order_fate <- rev(state_levels_main)
abund_col <- abund_col[row_order_fate, , drop = FALSE]
path_mat  <- path_mat[row_order_fate, , drop = FALSE]

# colour scales
abund_max <- max(abs(abund_col), na.rm = TRUE)
path_max  <- max(abs(path_mat), na.rm = TRUE)

ht_abund <- Heatmap(
  abund_col,
  name = "Abundance\nlog2FC",
  col = colorRamp2(
    c(-abund_max, 0, abund_max),
    c("#1B4965", "white", "#9B2226")
  ),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "left",
  show_column_names = FALSE,
  width = unit(2.5, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", abund_col[i, j]), x, y, gp = gpar(fontsize = 14))
  },
  rect_gp = gpar(col = "grey40", lwd = 1.5),
  column_title = "Abundance",
  column_title_gp = gpar(fontface = "bold", fontsize = 14)
)

ht_path <- Heatmap(
  path_mat,
  name = "Pathway\ndelta",
  col = colorRamp2(c(-path_max, 0, path_max), c("#245F7B", "white", "#B63E2F")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  column_names_rot = 35,
  width = unit(12, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", path_mat[i, j]), x, y, gp = gpar(fontsize = 14))
  },
  column_title = "Pathway composite deltas",
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  show_row_names = FALSE
)

ht_fate <- ht_abund + ht_path
draw(ht_fate,
     column_title = "Per State Summary",
     column_title_gp = gpar(fontface = "bold", fontsize = 18),
     merge_legend = FALSE,
     show_heatmap_legend = FALSE)

# ================================================================
# PAGE 6 — UMAP support (ggplot + patchwork)
# ================================================================
message("  Page 6: UMAP support ...")

matched_samples <- c(
  "SUR1070_Treated_PDO", "SUR1070_Untreated_PDO",
  "SUR1090_Treated_PDO", "SUR1090_Untreated_PDO",
  "SUR1072_Treated_PDO", "SUR1072_Untreated_PDO",
  "SUR1181_Treated_PDO", "SUR1181_Untreated_PDO"
)

message("    Loading Seurat for UMAP ...")
pdos <- readRDS("PDOs_merged.rds")
state_vec <- readRDS("Auto_PDO_final_states.rds")
pdos$state <- state_vec[Cells(pdos)]
pdos_matched <- subset(pdos, subset = orig.ident %in% matched_samples)
pdos_matched$Treatment <- ifelse(grepl("_Treated_", pdos_matched$orig.ident), "Treated", "Untreated")
pdos_matched$Treatment <- factor(pdos_matched$Treatment, levels = c("Untreated", "Treated"))
pdos_matched$state <- factor(as.character(pdos_matched$state), levels = state_levels_all)
rm(pdos); invisible(gc())

# gene-set scores
message("    Computing support scores ...")
library(msigdbr)
hallmark_tbl <- msigdbr(species = "Homo sapiens", category = "H")
mean_score <- function(mat, genes) {
  g <- intersect(unique(genes), rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
}
data_mat <- GetAssayData(pdos_matched, assay = "RNA", layer = "data")
injury_genes <- unique(hallmark_tbl$gene_symbol[hallmark_tbl$gs_name %in%
  c("HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY")])
adaptive_genes <- unique(hallmark_tbl$gene_symbol[hallmark_tbl$gs_name %in%
  c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")])

umap_df <- as.data.frame(Embeddings(pdos_matched, reduction = "umap")) %>%
  rownames_to_column("cell") %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2) %>%
  mutate(
    state = pdos_matched$state[match(cell, Cells(pdos_matched))],
    Treatment = pdos_matched$Treatment[match(cell, Cells(pdos_matched))],
    injury_score = mean_score(data_mat, injury_genes)[match(cell, colnames(data_mat))],
    adaptive_score = mean_score(data_mat, adaptive_genes)[match(cell, colnames(data_mat))]
  )

rm(pdos_matched, data_mat); invisible(gc())

umap_base <- function(df, ...) {
  ggplot(df, aes(UMAP_1, UMAP_2, ...)) +
    geom_point(size = 0.18, alpha = 0.8) +
    coord_equal() +
    pres_theme(9) +
    theme(legend.position = "right")
}

pu1 <- umap_base(umap_df, color = state) +
  scale_color_manual(values = state_cols, drop = FALSE) +
  labs(title = "Finalized states", color = NULL)

pu2 <- umap_base(umap_df, color = Treatment) +
  scale_color_manual(values = treatment_cols) +
  labs(title = "Untreated vs FLOT-treated", color = NULL)

pu3 <- umap_base(umap_df, color = injury_score) +
  scale_color_gradientn(colors = c("#1B365D", "#4F9DD9", "#F2E6B6", "#B33C2F"),
                        na.value = "grey85") +
  labs(title = "Injury score", color = NULL)

pu4 <- umap_base(umap_df, color = adaptive_score) +
  scale_color_gradientn(colors = c("#1B365D", "#4F9DD9", "#F2E6B6", "#B33C2F"),
                        na.value = "grey85") +
  labs(title = "Adaptation score", color = NULL)

pu1 <- pu1 + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15), legend.key.size = unit(1, "cm")) + guides(color = guide_legend(override.aes = list(size = 6)))
pu2 <- pu2 + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15), legend.key.size = unit(1, "cm")) + guides(color = guide_legend(override.aes = list(size = 6)))
pu3 <- pu3 + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15), legend.key.size = unit(1, "cm")) + guides(color = guide_colorbar(barwidth = 1.5, barheight = 10, title.theme = element_text(size = 15), label.theme = element_text(size = 14)))
pu4 <- pu4 + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 15), legend.key.size = unit(1, "cm")) + guides(color = guide_colorbar(barwidth = 1.5, barheight = 10, title.theme = element_text(size = 15), label.theme = element_text(size = 14)))

p_umap_panel <- (pu1 | pu2) / (pu3 | pu4) +
  plot_annotation(title = "PDO Matched-Pair Support UMAPs", theme = theme(plot.title = element_text(size = 16, face = "bold")))

print(p_umap_panel)

####################
# close PDF
####################
dev.off()
message("Done. Presentation PDF: ", file.path(getwd(), out_dir,
        "Auto_pdo_flot_presentation_final.pdf"))
