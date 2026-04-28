####################
# Auto_PDO_finalize_states.R
#
# 1. Finalize PDO cell states (merges) - Standardize naming to scRef conventions
# 2. Regenerate per-cell heatmap (Exact scRef Style: %v% partitions)
# 3. Regenerate TCGA survival volcano (scRef names, GSVA aggregation)
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(grid)
library(survival)
library(survminer)
library(GSVA)
library(data.table)
library(stringr)
library(ggrepel)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

message("=== Loading data ===")

# Load core data - use PDOs_merged.rds as PDOs_final.rds may be corrupted
pdos <- readRDS("PDOs_merged.rds")
ucell_scores <- readRDS("UCell_scores_filtered.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")

# Load relabeled states (These should now have the correct scRef space-separated names)
states_path <- "unresolved_states/Auto_PDO_unresolved_relabel_states.rds"
if (!file.exists(states_path)) {
  stop("Relabeled states file not found at: ", states_path)
}
states <- readRDS(states_path)

message("=== Performing finalized merges (Upstream states already corrected) ===")
# Convert to character vector
states_final <- as.character(states)
names(states_final) <- names(states)

# Perform requested merges using correct space-separated names
# 1) Respiration 1 (from 3CA) -> Classic Proliferative
states_final[states_final == "3CA_mp_30 Respiration 1"] <- "Classic Proliferative"
# 1b) Cell-cycle/HMG-rich 3CA state also folds into the classic proliferative compartment
####################
states_final[states_final %in% c("3CA_mp_3 Cell Cylce HMG-rich", "3CA_mp_3 Cell Cycle HMG-rich")] <- "Classic Proliferative"
####################
# 2) EMT III + Protein maturation -> combined state
states_final[states_final %in% c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III", "3CA_mp_17 EMT-III")] <- "3CA_EMT_and_Protein_maturation"
# 3) Ensure Basal name matches visualization (avoid premature renaming)
states_final[states_final == "Basal to Intestinal Metaplasia"] <- "Basal to Intest. Meta"

# Save finalized states vector
saveRDS(states_final, "Auto_PDO_final_states.rds")

# Update Seurat (and save the 2GB object as per user request to resave)
pdos$state <- states_final[Cells(pdos)]
pdos$Batch_fixed <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")
#saveRDS(pdos, "PDOs_final.rds")

####################
# Visualization Setup
####################

# Helper: Z-normalise function
z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]
  study_sd <- clust_df %>%
    dplyr::group_by(.study) %>%
    dplyr::summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1
  clust_centered <- clust_df %>%
    dplyr::group_by(.sample) %>%
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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

# Retained 3CA IDs (matching Slide 1 exactly)
retained_3ca_targets <- c("3CA_mp_30 Respiration 1", "3CA_mp_17 EMT III", "3CA_mp_12 Protein maturation")
retained_3ca_ids <- colnames(ucell_3ca)[clean_3ca_name(colnames(ucell_3ca)) %in% retained_3ca_targets]

# MP filtering and tree order
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]

mp.genes <- geneNMF.metaprograms$metaprograms.genes
if (length(bad_mps) > 0 || length(low_coverage_mps) > 0) {
  bad_targets <- unique(c(paste0("MP", bad_mps), low_coverage_mps))
  mp.genes <- mp.genes[!names(mp.genes) %in% bad_targets]
}
retained_mps <- names(mp.genes)

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order_names <- paste0("MP", rev(unique(ordered_clusters)))
mp_tree_order_names <- mp_tree_order_names[mp_tree_order_names %in% retained_mps]

# MP descriptions (matching PDO resolution)
mp_descriptions <- c(
  "MP6"  = "G2M Cell Cycle", "MP7"  = "DNA repair",
  "MP5"  = "MYC-related Proliferation", "MP1"  = "G2M checkpoint",
  "MP3"  = "G1S Cell Cycle", "MP8"  = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.", "MP9"  = "ECM Remodeling Epi.",
  "MP4"  = "Intestinal Metaplasia"
)
# Add MP ID prefixes for vibrant style matching
mp_descriptions <- setNames(paste0(names(mp_descriptions), "_", mp_descriptions), names(mp_descriptions))
cc_mps <- c("MP6", "MP7", "MP1", "MP3")

# State definition
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "SMG-like Metaplasia"   = c("MP8")
)

####################
# Canonical finalized PDO state order:
# Classic Proliferative -> Basal to Intest. Meta -> SMG-like Metaplasia ->
# Stress-adaptive -> 3CA_EMT_and_Protein_maturation.
####################
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9")
)
####################

state_level_order <- c("Classic Proliferative", "Basal to Intest. Meta", "Stress-adaptive", "SMG-like Metaplasia", "3CA_EMT_and_Protein_maturation", "Unresolved", "Hybrid")
####################
state_level_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation",
  "Unresolved",
  "Hybrid"
)
####################

group_cols <- c(
  "Classic Proliferative" = "#E41A1C", # Red
  "Basal to Intest. Meta" = "#4DAF4A", # Green
  "Stress-adaptive"       = "#984EA3", # Purple
  "SMG-like Metaplasia"   = "#FF7F00", # Orange
  "3CA_EMT_and_Protein_maturation" = "#377EB8", # Blue
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

####################
# Part 1: Per-cell Heatmap
####################
message("=== Generating per-cell heatmap (Exact ScRef Slide 1 Style) ===")

common_cells <- intersect(rownames(ucell_scores), Cells(pdos))
sample_var <- setNames(pdos$orig.ident, Cells(pdos))
study_var <- setNames(pdos$Batch_fixed, Cells(pdos))

# Normalize matrices
mp_adj_p <- z_normalise(as.matrix(ucell_scores[common_cells, retained_mps]), sample_var, study_var)
mp_adj_3 <- z_normalise(as.matrix(ucell_3ca[common_cells, retained_3ca_ids]), sample_var, study_var)
mp_raw_3 <- as.matrix(ucell_3ca[common_cells, retained_3ca_ids])

# CC scores and Diversity
cell_cycle_genes <- read.csv("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv")[,1:3]
cc_consensus <- intersect(cell_cycle_genes$Gene[cell_cycle_genes$Consensus==1], rownames(pdos))
cc_score <- colMeans(as.matrix(pdos@assays$RNA$data[intersect(cc_consensus, rownames(pdos)), , drop=FALSE]))
state_div_map <- pdos@meta.data %>% group_by(state) %>%
  summarise(score = 0.5*(n_distinct(orig.ident)/length(unique(pdos$orig.ident))) + 0.5*(n_distinct(Batch_fixed)/length(unique(pdos$Batch_fixed)))) %>% tibble::deframe()

# Sampling
set.seed(42)
MAX_CELLS <- 8000
state_assigned <- pdos$state[common_cells]
state_cells <- split(common_cells, state_assigned)
cells_to_plot <- unlist(mapply(function(cells, n) sample(cells, min(length(cells), n)),
                               state_cells, floor(MAX_CELLS * table(state_assigned) / sum(table(state_assigned))), SIMPLIFY=FALSE))

# Row blocks - NEW ORDER: Cell cycle -> State defining (by state_groups order) -> Other 3CA MPs
# 1. Cell cycle MPs first
cc_block_order <- cc_mps[cc_mps %in% retained_mps]

# 2. State defining MPs in state_groups order (Classic Prolif -> Basal to Intest -> SMG-like -> Stress-adaptive)
state_mp_order <- unlist(state_groups)  # MP5, MP4, MP8, MP10, MP9
state_mp_order <- state_mp_order[state_mp_order %in% retained_mps]

# 3. Remaining non-CC, non-state MPs (other 3CA MPs in tree order)
other_mps <- setdiff(retained_mps, c(cc_block_order, state_mp_order))
other_mps <- other_mps[order(match(other_mps, mp_tree_order_names))]

adj_rows_pipeline <- c(cc_block_order, state_mp_order, other_mps)
# Manual order for 3CA as in Slide 1: Resp1, EMT III, ProtMat
adj_rows_3ca <- c(
  grep("mp_30", colnames(ucell_3ca), value=T),
  grep("mp_17", colnames(ucell_3ca), value=T),
  grep("mp_12", colnames(ucell_3ca), value=T)
)
adj_rows_3ca <- intersect(adj_rows_3ca, retained_3ca_ids)

sub_p <- t(mp_adj_p[cells_to_plot, adj_rows_pipeline])
sub_3 <- t(mp_adj_3[cells_to_plot, adj_rows_3ca])
rownames(sub_3) <- clean_3ca_name(rownames(sub_3))
sub_scores_raw <- t(mp_raw_3[cells_to_plot, adj_rows_3ca])
rownames(sub_scores_raw) <- paste0(clean_3ca_name(rownames(sub_scores_raw)), " (raw)")

# Full Adj score matrix
sub_scores_adj <- rbind(sub_p, sub_3)
row_ids_adj <- c(adj_rows_pipeline, rownames(sub_3))

# Rename display labels
rownames(sub_scores_adj) <- ifelse(!is.na(mp_descriptions[rownames(sub_scores_adj)]), 
                                   mp_descriptions[rownames(sub_scores_adj)], rownames(sub_scores_adj))

# MP Group left annotation
mp_group_names <- c(
  "Cell_cycle" = "Cell cycle",
  "Classic Proliferative" = "Classic\nProlif",
  "Basal to Intest. Meta" = "Basal-IM",
  "Stress-adaptive" = "Stress\nadaptive",
  "SMG-like Metaplasia" = "SMG-like\nMeta",
  "PanCancer" = "Pan-Cancer\n(norm)",
  "PanCancerRaw" = "Pan-Cancer\n(raw)",
  "Other" = "Other"
)

group_colors_v2 <- c(
  "Cell cycle" = "gold",
  "Classic\nProlif" = "#E41A1C",
  "Basal-IM" = "#4DAF4A",
  "Stress\nadaptive" = "#984EA3",
  "SMG-like\nMeta" = "#FF7F00",
  "Pan-Cancer\n(norm)" = "#4B0082",
  "Pan-Cancer\n(raw)" = "#7F0000",
  "Other" = "grey70"
)

get_mp_group <- function(rn_raw) {
  if (rn_raw %in% cc_mps) return("Cell_cycle")
  if (rn_raw %in% state_groups[["Classic Proliferative"]]) return("Classic Proliferative")
  if (rn_raw %in% state_groups[["Basal to Intest. Meta"]]) return("Basal to Intest. Meta")
  if (rn_raw %in% state_groups[["Stress-adaptive"]]) return("Stress-adaptive")
  if (rn_raw %in% state_groups[["SMG-like Metaplasia"]]) return("SMG-like Metaplasia")
  if (clean_3ca_name(rn_raw) %in% retained_3ca_targets) return("PanCancer")
  return("Other")
}

mp_group_adj_vec <- factor(mp_group_names[sapply(row_ids_adj, get_mp_group)], levels=mp_group_names)
mp_group_raw_vec <- factor(rep("Pan-Cancer\n(raw)", nrow(sub_scores_raw)), levels=mp_group_names)

row_ann_adj <- rowAnnotation(MP_group = mp_group_adj_vec, col = list(MP_group = group_colors_v2), show_annotation_name = FALSE,
                             annotation_legend_param = list(MP_group = list(title = "MP Group", title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))))
row_ann_raw <- rowAnnotation(MP_group = mp_group_raw_vec, col = list(MP_group = group_colors_v2), show_annotation_name = FALSE, show_legend = FALSE)

# Top Annotation
split_cols <- factor(as.character(pdos$state[cells_to_plot]), levels=state_level_order)
col_ann <- HeatmapAnnotation(
  State = split_cols,
  CNA = rep("cna_malignant", length(cells_to_plot)),
  CC_score = cc_score[cells_to_plot],
  Diversity = state_div_map[as.character(split_cols)],
  Study = pdos$Batch_fixed[cells_to_plot],
  col = list(
    State = group_cols,
    CNA = c(cna_malignant="black", cna_unresolved="grey70"),
    CC_score = colorRamp2(c(0, max(cc_score)), c("white", "darkgreen")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Study = c(Cynthia_batch="brown", New_batch="darkgreen")
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    State = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    CNA = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    CC_score = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    Diversity = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    Study = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))
  )
)

# Row split factor for Adjusted part - match new order
row_split_adj_factor <- factor(
  case_when(
    row_ids_adj %in% cc_block_order ~ "Cell cycle\nMPs",
    row_ids_adj %in% state_mp_order ~ "State defining\nMPs",
    clean_3ca_name(row_ids_adj) %in% retained_3ca_targets ~ "Pan-Cancer\n(norm)",
    TRUE ~ "Other\nMPs"
  ),
  levels = c("Cell cycle\nMPs", "State defining\nMPs", "Pan-Cancer\n(norm)", "Other\nMPs")
)

# Colors (dynamically scaled to 98th percentile)
lim_adj <- as.numeric(quantile(abs(sub_scores_adj), 0.98, na.rm = TRUE))
col_fun_adj <- colorRamp2(c(-lim_adj, 0, lim_adj), c("navy", "white", "firebrick3"))

raw_lim <- as.numeric(quantile(sub_scores_raw, 0.98, na.rm = TRUE))
if (is.na(raw_lim) || raw_lim == 0) raw_lim <- 0.15 # Fallback if data is too sparse
col_fun_raw <- colorRamp2(c(0, raw_lim), c("white", "firebrick3"))

# Column order logic: Cluster columns within each State split
col_order <- (function() {
  col_order_list <- lapply(levels(split_cols), function(lvl) {
    idx <- which(as.character(split_cols) == lvl)
    if (length(idx) <= 1) return(idx)
    mat_lvl <- sub_scores_adj[, idx, drop = FALSE]
    dcols <- dist(t(mat_lvl))
    hc <- hclust(dcols, method = "ward.D2")
    idx[hc$order]
  })
  full_ord <- unlist(col_order_list, use.names = FALSE)
  if (length(full_ord) != ncol(sub_scores_adj) || !setequal(full_ord, seq_len(ncol(sub_scores_adj)))) {
    return(seq_len(ncol(sub_scores_adj)))
  }
  full_ord
})()

# Heatmap Draw
ht_adj <- Heatmap(sub_scores_adj, name = "Adj score", col = col_fun_adj, top_annotation = col_ann, left_annotation = row_ann_adj,
                  column_split = split_cols, column_order = col_order, row_split = row_split_adj_factor, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
                  row_names_side = "left", row_names_gp = gpar(fontsize=8, fontface="italic"), column_title_rot = 30, column_title_gp = gpar(fontsize=10, fontface="bold"),
                  row_gap = unit(2, "mm"), column_gap = unit(1.5, "mm"), border = FALSE, rect_gp = gpar(col = NA),
                  use_raster = TRUE, raster_quality = 3,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)))

ht_raw <- Heatmap(sub_scores_raw, name = "Raw score", col = col_fun_raw, left_annotation = row_ann_raw,
                  column_split = split_cols, column_order = col_order,cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
                  row_names_side = "left", row_names_gp = gpar(fontsize=8, fontface="italic"), column_title = NULL,
                  column_gap = unit(1.5, "mm"), border = FALSE, rect_gp = gpar(col = NA),
                  use_raster = TRUE, raster_quality = 3,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)))

pdf("Auto_PDO_final_states_heatmap.pdf", width = 18, height = 11)
draw(ht_adj %v% ht_raw, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
grid.text("Unresolved cells subclass: 3CA-based relabeling (noreg)", x = unit(5, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

####################
# Part 1.5: Proportions Plot
####################
message("=== Generating state proportions plot ===")
prop_df <- data.frame(
  state = as.character(pdos$state),
  Batch = as.character(pdos$Batch_fixed),
  stringsAsFactors = FALSE
)
overall <- prop_df %>% count(state) %>% mutate(Batch = "Overall", pct = 100 * n / sum(n))
per_batch <- prop_df %>%
  count(Batch, state) %>%
  group_by(Batch) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
plot_df <- bind_rows(overall, per_batch)
plot_df$state <- factor(plot_df$state, levels = state_level_order)

p_bar <- ggplot(plot_df, aes(Batch, pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(title = "Finalized State Proportions", x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pie_df <- overall %>% 
  mutate(label = ifelse(pct < 5, "", paste0(state, "\n", sprintf("%.1f%%", pct))))

p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(title = "Overall Pie Chart", fill = "State") +
  theme_void(base_size = 11)

pdf("Auto_PDO_final_states_proportion.pdf", width = 12, height = 6)
print(p_bar + p_pie + plot_layout(widths = c(2, 1)))
dev.off()

####################
# Part 2: TCGA Survival Volcano
####################
# Aggregates OS hazard by state scores (scRef names)
####################
message("=== Generating TCGA survival volcano ===")
ref_3ca <- fread("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv", header=TRUE)
sets_3ca <- list(
  "3CA_mp_30" = na.omit(ref_3ca[[grep("MP30 Respiration 1", colnames(ref_3ca), value=TRUE)]]),
  "3CA_mp_12" = na.omit(ref_3ca[[grep("MP12 Protein maturation", colnames(ref_3ca), value=TRUE)]]),
  "3CA_mp_17" = na.omit(ref_3ca[[grep("MP17 EMT-III", colnames(ref_3ca), value=TRUE)]])
)
gsva_sets <- c(lapply(mp.genes, unique), sets_3ca)

tpm_df <- fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1]); rownames(tpm_mat) <- tpm_df$GeneSymbol
meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")

gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat))); gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]
gsva_res <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
gsva_df <- as.data.frame(t(gsva_res)) %>% mutate(sample_barcode = rownames(.))

surv_data <- meta_tcga %>% inner_join(gsva_df, by = "sample_barcode") %>% 
  filter(sample_type_code == "01", grepl("adeno", tolower(as.character(type))))

# Aggregate
surv_data[["Classic Proliferative"]] <- apply(surv_data[, intersect(c("MP5", "3CA_mp_30"), colnames(surv_data)), drop=FALSE], 1, max)
surv_data[["Basal to Intest. Meta"]] <- apply(surv_data[, intersect(c("MP4"), colnames(surv_data)), drop=FALSE], 1, max)
surv_data[["SMG-like Metaplasia"]]   <- apply(surv_data[, intersect(c("MP8"), colnames(surv_data)), drop=FALSE], 1, max)
surv_data[["Stress-adaptive"]]       <- apply(surv_data[, intersect(c("MP10", "MP9"), colnames(surv_data)), drop=FALSE], 1, max)
surv_data[["3CA_EMT_and_Protein_maturation"]] <- apply(surv_data[, intersect(c("3CA_mp_12", "3CA_mp_17"), colnames(surv_data)), drop=FALSE], 1, max)

features <- c("Classic Proliferative", "Basal to Intest. Meta", "Stress-adaptive", "SMG-like Metaplasia", "3CA_EMT_and_Protein_maturation")
####################
features <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)
####################

cox_res <- list()
for (f in features) {
  d <- surv_data %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[f]]))
  if (nrow(d) < 20) next
  fit <- try(coxph(Surv(OS_time, OS_event) ~ .data[[f]], data = d), silent=TRUE)
  if (!inherits(fit, "try-error")) {
    ss <- summary(fit)
    cox_res[[f]] <- data.frame(feature=f, HR=ss$coefficients[1,"exp(coef)"], P_value=ss$coefficients[1,"Pr(>|z|)"])
  }
}
cox_res <- bind_rows(cox_res)
if (nrow(cox_res) > 0) {
  p_volc <- ggplot(cox_res %>% mutate(log2HR = log2(HR), neglog10 = -log10(P_value), sig = P_value < 0.05),
                   aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 4) +
    scale_color_manual(values = c("grey70", "firebrick3"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(aes(label = feature), size = 3, box.padding = 0.5) +
    theme_minimal() + labs(title = "Finalized States: Survival Volcano (EAC)", x = "log2(HR)", y = "-log10(p)")

  pdf("Auto_PDO_final_states_volcano.pdf", width = 8, height = 6)
  print(p_volc)
  dev.off()
} else {
  message("No survival results found (EAC). Skipping volcano.")
}

message("=== DONE ===")
