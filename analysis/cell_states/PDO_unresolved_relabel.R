####################
# Auto_PDO_unresolved_relabel.R
#
# Relabel unresolved Approach-B (noreg) cells by top pan-cancer (3CA) MP for PDO,
# retain 3-5 pan-cancer MPs by sample coverage thresholds,
# regenerate: proportions, heatmap, TCGA survival volcano
#
# Input:
#   PDOs_outs/PDOs_final.rds
#   PDOs_outs/Auto_PDO_states_noreg.rds
#   PDOs_outs/Auto_PDO_mp_adj_noreg.rds
#   PDOs_outs/UCell_scores_filtered.rds
#   PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#   /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#   /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds
#   /rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
#
# Output:
#   PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_states.rds
#   PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_proportion.pdf
#   PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_heatmap.pdf
#   PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_volcano.pdf
#   PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_cox_results.csv
#   PDOs_outs/unresolved_states/Auto_PDO_unresolved_relabel_mp_coverage.csv
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)
library(survival)
library(survminer)
library(GSVA)
library(patchwork)
library(data.table)
library(stringr)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
dir.create("unresolved_states/", recursive = TRUE, showWarnings = FALSE)

####################
# constants
####################
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "Basal to Intest. Meta" = c("MP4")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "SMG-like Metaplasia"   = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "Basal to Intest. Meta" = "#FF7F00",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

# MP descriptions for PDO
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

cc_mps <- c("MP6", "MP7", "MP1", "MP3")

batch_cols <- c(Cynthia_batch = "brown", New_batch = "darkgreen")

####################
# helper functions
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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  x <- gsub("\\.", " ", x)
  x
}

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}

run_cox_for_group <- function(df, features, cohort_name) {
  out <- list()
  for (feat in features) {
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]))
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    fit <- try(coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`")), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      cohort = cohort_name,
      feature = feat,
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

plot_volcano <- function(df, title_text) {
  if (nrow(df) == 0) return(NULL)
  df <- df %>% mutate(sig = P_value < 0.05, neglog10 = -log10(P_value), log2HR = log2(HR))
  ggplot(df, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 3, alpha = 0.85) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text_repel(aes(label = feature), size = 3.1, max.overlaps = 100) +
    theme_minimal(base_size = 13) +
    labs(title = title_text, x = "log2(HR)", y = "-log10(p)")
}

####################
# data loading
####################
message("Loading Seurat object...")
pdos <- readRDS("PDOs_final.rds")

# Fix batch
pdos$Batch <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")

message("Loading state assignments...")
state_B <- readRDS("Auto_PDO_states_noreg.rds")
pdos$state <- state_B[Cells(pdos)]

message("Loading MP adjustments...")
mp_adj_noncc <- readRDS("Auto_PDO_mp_adj_noreg.rds")

message("Loading UCell scores...")
ucell_scores <- readRDS("UCell_scores_filtered.rds")

####################
# STEP 1: Load existing 3CA MP scores
# (If not available, uncomment the following to compute from scratch)
####################
message("Loading existing 3CA UCell scores...")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")

# Alternative: Compute 3CA MP scores from scratch (slow - ~15 min)
# message("Reading New_NMFs.csv for 3CA MPs...")
# MP_df <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv", check.names = FALSE)
# MP_list <- as.list(MP_df)
# MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
# MP_list <- MP_list[sapply(MP_list, length) > 0]
# names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
# message("Adding 3CA MP modules to Seurat object (this may take ~15 min)...")
# pdos <- AddModuleScore_UCell(pdos, features = MP_list, ncores = 1, name = "")
# score_cols_3ca <- grep("^X3CA_mp|^3CA_mp", colnames(pdos@meta.data), value = TRUE)
# ucell_3ca <- pdos@meta.data[, score_cols_3ca, drop = FALSE]
# saveRDS(ucell_3ca, file = "UCell_3CA_MPs.rds")

message("Loading GeneNMF metaprograms...")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")

# Get MP genes and filtering
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
mp.genes <- mp.genes[!names(mp.genes) %in% c(bad_mp_names, low_coverage_mps)]
retained_mps <- names(mp.genes)

# MP tree order
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- rev(unique(ordered_clusters))
mp_tree_order_names <- paste0("MP", mp_tree_order)
mp_tree_order_names <- mp_tree_order_names[mp_tree_order_names %in% retained_mps]

# Filter state groups
state_groups <- lapply(state_groups, function(mps) mps[mps %in% retained_mps])
state_groups <- state_groups[sapply(state_groups, length) > 0]
ordered_group_names <- names(state_groups)
state_level_order <- c(ordered_group_names, "Unresolved", "Hybrid")

# Common cells - include 3CA scores
common_cells <- Reduce(
  function(x, y) intersect(x, y),
  list(names(state_B), Cells(pdos), rownames(mp_adj_noncc), rownames(ucell_scores), rownames(ucell_3ca))
)

pdos <- pdos[, common_cells]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
ucell_3ca <- ucell_3ca[common_cells, , drop = FALSE]

sample_var <- pdos$orig.ident
study_var <- pdos$Batch
names(sample_var) <- Cells(pdos)
names(study_var) <- Cells(pdos)

####################
# STEP 1: relabel unresolved cells by top 3CA MP
####################
unresolved_cells <- names(state_B)[state_B == "Unresolved"]
unresolved_3ca <- ucell_3ca[unresolved_cells, , drop = FALSE]

CC_FIXED <- c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
)

# remove cell-cycle MPs
unresolved_3ca_nocc <- unresolved_3ca[, !colnames(unresolved_3ca) %in% CC_FIXED, drop = FALSE]

# Calculate max column using raw UCell scores
top_3ca_mp <- colnames(unresolved_3ca_nocc)[max.col(unresolved_3ca_nocc, ties.method = "first")]
names(top_3ca_mp) <- unresolved_cells

####################
# STEP 2: coverage-threshold selection of 3CA MPs
# Adjust threshold for PDO (only 20 samples)
####################
unresolved_meta <- data.frame(
  cell = unresolved_cells,
  mp_label = top_3ca_mp,
  orig.ident = as.character(pdos$orig.ident[unresolved_cells]),
  Batch = as.character(pdos$Batch[unresolved_cells]),
  stringsAsFactors = FALSE
)

total_samples <- length(unique(pdos$orig.ident[common_cells]))
total_batches <- length(unique(pdos$Batch[common_cells]))

mp_coverage <- unresolved_meta %>%
  group_by(mp_label) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(orig.ident),
    n_batches = n_distinct(Batch),
    pct_samples = 100 * n_distinct(orig.ident) / total_samples,
    pct_batches = 100 * n_distinct(Batch) / total_batches,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

message("Pan-cancer MP coverage in unresolved cells:")
print(mp_coverage)

# Adjust threshold: for PDO with ~20 samples, use n_samples >= 10 (50%) or n_batches >= 2
# This is more lenient than scRef due to smaller sample size
candidates <- mp_coverage %>%
  filter(n_samples >= 10 | n_batches >= 2) %>%
  arrange(desc(n_cells))

# Keep top 3-5
retained_3ca <- head(candidates %>% pull(mp_label), 5)

message(paste("Retained pan-cancer MPs:", paste(retained_3ca, collapse = ", ")))
write.csv(mp_coverage, "unresolved_states/Auto_PDO_unresolved_relabel_mp_coverage.csv", row.names = FALSE)

####################
# STEP 3: update state labels
####################
state_updated <- state_B

for (cell in unresolved_cells) {
  mp <- top_3ca_mp[cell]
  if (mp %in% retained_3ca) {
    state_updated[cell] <- clean_3ca_name(mp)
  }
}

new_state_names <- unique(clean_3ca_name(retained_3ca))
state_level_order_updated <- c(
  ordered_group_names,
  sort(new_state_names),
  "Unresolved", "Hybrid"
)

new_state_cols <- setNames(
  scales::hue_pal()(length(new_state_names)),
  new_state_names
)
group_cols_updated <- c(group_cols, new_state_cols)

saveRDS(state_updated, "unresolved_states/Auto_PDO_unresolved_relabel_states.rds")

####################
# STEP 4: updated proportion plot
####################
prop_df <- data.frame(
  state = as.character(state_updated[common_cells]),
  Batch = as.character(pdos$Batch[common_cells]),
  stringsAsFactors = FALSE
)

overall <- prop_df %>% count(state) %>% mutate(Batch = "Overall", pct = 100 * n / sum(n))
per_batch <- prop_df %>%
  count(Batch, state) %>%
  group_by(Batch) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
plot_df <- bind_rows(overall, per_batch)
plot_df$state <- factor(plot_df$state, levels = state_level_order_updated)

p_bar <- ggplot(plot_df, aes(Batch, pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Updated state proportions (4 original + 3CA relabeled)", x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pie_df <- overall %>% 
  mutate(label = ifelse(pct < 5, "", paste0(state, "\n", sprintf("%.1f%%", pct))))

p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Updated state proportions: overall pie", fill = "State") +
  theme_void(base_size = 11)

ggsave(
  "unresolved_states/Auto_PDO_unresolved_relabel_proportion.pdf",
  p_bar + p_pie + plot_layout(widths = c(2, 1)),
  width = 18,
  height = 8
)

####################
# STEP 5: per-cell heatmap (all cells)
####################
cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc)

# Cell cycle score
cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]
cc_consensus <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(pdos))
cc_top50 <- names(sort(rowMeans(pdos@assays$RNA$data[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:50]
cc_score <- colMeans(as.matrix(pdos@assays$RNA$data[cc_top50, , drop = FALSE]))

set.seed(42)
MAX_CELLS_TOTAL <- 8000
state_counts <- table(state_updated)
state_fracs <- state_counts / sum(state_counts)
cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)
state_cells <- split(names(state_updated), state_updated)
cells_to_plot <- unlist(
  mapply(
    function(cells, n) sample(cells, min(length(cells), n)),
    state_cells,
    cells_per_state[names(state_cells)],
    SIMPLIFY = FALSE
  ),
  use.names = FALSE
)

sub_scores_orig <- t(mp_adj_all[cells_to_plot, , drop = FALSE])
cc_block_order <- cc_mps[cc_mps %in% rownames(sub_scores_orig)]
non_cc_block_order <- mp_tree_order_names[
  mp_tree_order_names %in% rownames(sub_scores_orig) & !(mp_tree_order_names %in% cc_mps)
]
mp_row_order <- c(cc_block_order, non_cc_block_order)
sub_scores <- sub_scores_orig[mp_row_order, , drop = FALSE]

# Apply MP descriptions
mp_label_map <- mp_descriptions
missing_mps <- setdiff(rownames(sub_scores), names(mp_label_map))
if (length(missing_mps) > 0) mp_label_map[missing_mps] <- missing_mps
rownames(sub_scores) <- mp_label_map[rownames(sub_scores)]

present_states <- intersect(state_level_order_updated, unique(as.character(state_updated[cells_to_plot])))
if (length(present_states) == 0) present_states <- unique(as.character(state_updated[cells_to_plot]))
split_vec <- factor(as.character(state_updated[cells_to_plot]), levels = present_states)

# Batch values
batch_vals <- pdos@meta.data[cells_to_plot, "Batch"]

max_cc <- max(cc_score[cells_to_plot], na.rm = TRUE)

col_ann <- HeatmapAnnotation(
  State = split_vec,
  CC_score = cc_score[cells_to_plot],
  Batch = batch_vals,
  col = list(
    State = group_cols_updated[present_states],
    CC_score = colorRamp2(c(0, max_cc), c("white", "darkgreen")),
    Batch = batch_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE,
  na_col = "white"
)

mp_to_group <- rep("Other", length(mp_row_order))
names(mp_to_group) <- mp_row_order
mp_to_group[cc_mps[cc_mps %in% names(mp_to_group)]] <- "Cell_cycle"
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], names(mp_to_group))
  mp_to_group[grp_mps] <- grp
}
group_colors_row <- c(group_cols[ordered_group_names], Cell_cycle = "gold", Other = "grey70")
mp_group_label <- mp_to_group
names(mp_group_label) <- rownames(sub_scores)

row_ann <- rowAnnotation(
  MP_group = factor(mp_group_label, levels = c("Cell_cycle", ordered_group_names, "Other")),
  col = list(MP_group = group_colors_row),
  show_annotation_name = FALSE
)

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

ht <- Heatmap(
  sub_scores,
  name = "Adj score",
  col = col_fun_sc,
  top_annotation = col_ann,
  left_annotation = row_ann,
  column_split = split_vec,
  column_order = (function() {
    col_order_list <- lapply(levels(split_vec), function(lvl) {
      idx <- which(as.character(split_vec) == lvl)
      if (length(idx) <= 1) return(idx)
      mat_lvl <- sub_scores[, idx, drop = FALSE]
      dcols <- dist(t(mat_lvl))
      hc <- hclust(dcols, method = "ward.D2")
      idx[hc$order]
    })
    full_ord <- unlist(col_order_list, use.names = FALSE)
    if (length(full_ord) != ncol(sub_scores) || !setequal(full_ord, seq_len(ncol(sub_scores)))) {
      return(seq_len(ncol(sub_scores)))
    }
    full_ord
  })(),
  column_gap = unit(1.5, "mm"),
  row_split = factor(
    ifelse(mp_row_order %in% cc_mps, "Cell_cycle_MPs", "Other_MPs"),
    levels = c("Cell_cycle_MPs", "Other_MPs")
  ),
  row_gap = unit(2.5, "mm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  show_column_names = FALSE,
  column_title_rot = 30,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 5,
  border = FALSE,
  rect_gp = gpar(col = NA)
)

pdf("unresolved_states/Auto_PDO_unresolved_relabel_heatmap.pdf", width = 18, height = 8, useDingbats = FALSE)
draw(ht, merge_legend = TRUE)
dev.off()

####################
# STEP 6: survival volcano using GSVA
####################
message("Running TCGA survival analysis...")

# Load TCGA data
meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

# Original MP genes
gsva_sets <- lapply(mp.genes, unique)

# Load 3CA genes
MP_df <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv", check.names = FALSE)
MP_list <- as.list(MP_df)
MP_list <- lapply(MP_list, function(x) x[x != "" & !is.na(x)])
names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))

# Filter to retained 3CA MPs
new_state_sigs <- MP_list[retained_3ca]
names(new_state_sigs) <- clean_3ca_name(names(new_state_sigs))

# Append 3CA gene sets
gsva_sets <- c(gsva_sets, new_state_sigs)

# Filter and run GSVA
gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]

gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

# Merge with TCGA metadata
surv_data <- meta_tcga %>%
  inner_join(gsva_df, by = "sample_barcode") %>%
  filter(sample_type_code == "01")

surv_data$HistologyGroup <- infer_histology(surv_data$type)

# Aggregate MP scores into State scores
for (nm in names(state_groups)) {
  mps <- intersect(state_groups[[nm]], colnames(surv_data))
  if (length(mps) == 0) next
  surv_data[[nm]] <- apply(as.matrix(surv_data[, mps, drop = FALSE]), 1, max)
}

# Combine original states and 3CA states
state_cols <- intersect(c(names(state_groups), new_state_names), colnames(surv_data))

all_cox <- list()

pdf("unresolved_states/Auto_PDO_unresolved_relabel_volcano.pdf", width = 9, height = 7)
for (coh in c("EAC")) {
  cox_df <- run_cox_for_group(
    surv_data %>% filter(HistologyGroup == coh),
    state_cols,
    cohort_name = coh
  )
  all_cox[[coh]] <- cox_df
  p <- plot_volcano(cox_df, paste0("Updated states: TCGA survival volcano (", coh, ")"))
  if (!is.null(p)) print(p)
}
dev.off()

cox_res <- bind_rows(all_cox)
write.csv(cox_res, "unresolved_states/Auto_PDO_unresolved_relabel_cox_results.csv", row.names = FALSE)

message("=== DONE ===")
message("All outputs saved to PDOs_outs/unresolved_states/")
