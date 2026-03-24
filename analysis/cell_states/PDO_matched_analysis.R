####################
# Auto_PDO_matched_analysis_v2.R
# Fixes for matched sample analysis
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(survival)
library(survminer)
library(GSVA)
library(data.table)
library(stringr)
library(ggrepel)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

message("=== Loading data ===")

# Load data
message("Loading PDOs_final.rds...")
pdos <- readRDS("PDOs_final.rds")
message("Loading UCell_scores_filtered.rds...")
ucell_scores <- readRDS("UCell_scores_filtered.rds")
message("Loading geneNMF_metaprograms_nMP_13.rds...")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")

# Fix batch
pdos$Batch <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")

# Get MP info
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

# MP descriptions - map MP numbers to descriptions (filter to only retained MPs)
mp_desc_map <- c(
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
# Filter to only retained MPs
mp_desc_map <- mp_desc_map[names(mp_desc_map) %in% retained_mps]

# State groups (scRef nomenclature)
state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "Basal to Intest. Meta" = c("MP4"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "SMG-like Metaplasia"   = c("MP8")
)
state_groups <- lapply(state_groups, function(mps) mps[mps %in% retained_mps])
state_groups <- state_groups[sapply(state_groups, length) > 0]

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))
state_level_order <- c(ordered_group_names, "3CA_EMT_and_Protein_maturation", "Unresolved", "Hybrid")

# State colors
group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)
group_cols <- group_cols[names(group_cols) %in% state_level_order]

# Common cells
common_cells <- intersect(rownames(ucell_scores), Cells(pdos))
pdos <- pdos[, common_cells]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]

message("=== Loading state assignments ===")

# Load finalized states
state_B <- readRDS("Auto_PDO_final_states.rds")
pdos$state <- state_B[Cells(pdos)]

message("=== TCGA Survival Analysis ===")

# Load TCGA data
geneNMF.metaprograms_sc <- readRDS("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
meta_tcga <- readRDS("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/tcga_esca_meta.rds")
tpm_df <- data.table::fread("/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt")
tpm_mat <- as.matrix(tpm_df[, -1])
rownames(tpm_mat) <- tpm_df$GeneSymbol

# Use PDO MP genes for scoring
gsva_sets <- lapply(mp.genes, unique)
gsva_sets <- lapply(gsva_sets, function(g) intersect(g, rownames(tpm_mat)))
gsva_sets <- gsva_sets[sapply(gsva_sets, length) >= 5]
message("Calculating GSVA scores (this may take time)...")
gsva_scores <- gsva(tpm_mat, gsva_sets, method = "gsva", kcdf = "Gaussian")
message("Filtering surv_data...")
gsva_df <- as.data.frame(t(gsva_scores))
gsva_df$sample_barcode <- rownames(gsva_df)

# Filter: only EAC, only sample_type_code == "01"
surv_data <- meta_tcga %>% inner_join(gsva_df, by = "sample_barcode") %>% filter(sample_type_code == "01")

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out[grepl("squamous", t)] <- "ESCC"
  out
}
surv_data$HistologyGroup <- infer_histology(surv_data$type)
surv_data <- surv_data %>% filter(HistologyGroup == "EAC")

# Compute state scores
for (nm in names(state_groups)) {
  mps <- intersect(state_groups[[nm]], colnames(surv_data))
  if (length(mps) == 0) next
  surv_data[[nm]] <- apply(as.matrix(surv_data[, mps, drop = FALSE]), 1, max)
}
state_cols <- intersect(names(state_groups), colnames(surv_data))

# Run Cox
run_cox_for_group <- function(df, features, cohort_name, mode_name) {
  out <- list()
  for (feat in features) {
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]))
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    fit <- try(coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`")), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      mode = mode_name,
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

# Volcano plot - NO rotation
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

# State volcano
all_cox_state <- run_cox_for_group(surv_data, state_cols, "EAC", "noreg")

pdf("Auto_PDO_survival_state_volcano_EAC_noreg.pdf", width = 9, height = 7, useDingbats = FALSE)
p <- plot_volcano(all_cox_state, "[noreg] PDO State survival volcano (EAC)")
if (!is.null(p)) print(p)
dev.off()

# MP columns present in surv_data
mp_cols_in_surv <- colnames(surv_data)[grepl("^MP\\d+$", colnames(surv_data))]

# keep only MPs that have a description
mp_cols_in_surv <- mp_cols_in_surv[mp_cols_in_surv %in% names(mp_desc_map)]

# sort by MP number
mp_cols_sorted <- mp_cols_in_surv[order(as.numeric(sub("^MP", "", mp_cols_in_surv)))]

# rename those columns in surv_data to descriptions
colnames(surv_data)[match(mp_cols_sorted, colnames(surv_data))] <- mp_desc_map[mp_cols_sorted]

# now the MP column names are descriptions
mp_cols_sorted <- unname(mp_desc_map[mp_cols_sorted])

# Run cox 
all_cox_mp <- run_cox_for_group(surv_data, mp_cols_sorted, "EAC", "noreg")

pdf("Auto_PDO_survival_mp_volcano_EAC_noreg.pdf", width = 9, height = 7, useDingbats = FALSE)
p <- plot_volcano(all_cox_mp, "[noreg] PDO MP survival volcano (EAC)")
if (!is.null(p)) print(p)
dev.off()

# MP Kaplan-Meier - use descriptions (columns are already descriptions in surv_data)
pdf("Auto_PDO_survival_mp_km_EAC_noreg.pdf", width = 8, height = 7, useDingbats = FALSE)
for (mp_col in mp_cols_sorted) {
  d <- surv_data %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[mp_col]]))
  if (nrow(d) < 20) next
  med <- median(d[[mp_col]], na.rm = TRUE)
  d$Group <- factor(ifelse(d[[mp_col]] >= med, "High", "Low"), levels = c("Low", "High"))
  if (length(unique(d$Group)) < 2) next
  fit <- survfit(Surv(OS_time, OS_event) ~ Group, data = d)
  p <- ggsurvplot(
    fit,
    data = d,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    title = paste0("[noreg] EAC | ", mp_col, " (median split)"),
    xlab = "Days",
    ylab = "Overall survival"
  )
  print(p)
}
dev.off()

message("=== Matched Sample Analysis ===")

# Define matched samples
matched_samples <- c(
  "SUR1070_Treated_PDO", "SUR1070_Untreated_PDO",
  "SUR1072_Treated_PDO", "SUR1072_Untreated_PDO",
  "SUR1090_Treated_PDO", "SUR1090_Untreated_PDO",
  "SUR1181_Treated_PDO", "SUR1181_Untreated_PDO"
)

# Filter to matched samples
matched_cells <- which(pdos$orig.ident %in% matched_samples)
pdos_matched <- pdos[, matched_cells]
ucell_matched <- ucell_scores[matched_cells, ]

# Create treatment and patient variables
pdos_matched$Treatment <- ifelse(grepl("Treated", pdos_matched$orig.ident), "Treated", "Untreated")
pdos_matched$Patient <- str_extract(pdos_matched$orig.ident, "^SUR\\d+")

# Patient colors - distinct colors for each patient
patient_cols_base <- c(
  SUR1070 = "#4C78A8",
  SUR1072 = "#59A14F",
  SUR1090 = "#B07AA1",
  SUR1181 = "#F28E2B"
)

# Sample order: untreated before treated for each patient
patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
sample_order <- c()
for (pat in patient_order) {
  sample_order <- c(sample_order, paste0(pat, "_Untreated_PDO"), paste0(pat, "_Treated_PDO"))
}

# Build sample color mapping - patient color with treatment modifier
sample_colors <- c()
for (pat in patient_order) {
  base <- patient_cols_base[pat]
  untreated_col <- grDevices::adjustcolor(base, alpha.f = 1)
  treated_col <- grDevices::adjustcolor(base, alpha.f = 0.7)
  sample_colors[paste0(pat, "_Untreated_PDO")] <- untreated_col
  sample_colors[paste0(pat, "_Treated_PDO")] <- treated_col
}

# State colors
state_cols_plot <- group_cols

# Filter out unresolved and hybrid cells for UMAP
# UMAP plots will use all finalized states (no filtering)

# Helper for colors
lighten_col <- function(col, amount = 0.55) {
  if (is.na(col)) return(NA)
  rgb_col <- col2rgb(col) / 255
  new_col <- rgb_col + (1 - rgb_col) * amount
  grDevices::rgb(new_col[1], new_col[2], new_col[3])
}

# UMAP plot function
make_umap_plot <- function(seurat_obj, col_var, col_map, title, reduction = "umap") {
  if (reduction == "umap" && is.null(seurat_obj@reductions$umap)) {
    message("UMAP not found, skipping...")
    return(NULL)
  }
  df_plot <- data.frame(
    UMAP1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
    UMAP2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
    val = seurat_obj@meta.data[[col_var]]
  )
  ggplot(df_plot, aes(UMAP1, UMAP2, color = val)) +
    geom_point(size = 0.3, alpha = 0.8) +
    scale_color_manual(values = col_map, name = col_var, breaks = names(col_map)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_minimal() +
    labs(title = title) +
    theme(aspect.ratio = 1, legend.position = "right", plot.title = element_text(hjust = 0.5))
}

# Generate UMAP and Proportion Summary Report
pdf("Auto_PDO_matched_summary_report_v2.pdf", width = 18, height = 6, useDingbats = FALSE)

# helper for state barplots (stacked)
make_prop_plot <- function(df, title_text, fill_map) {
  df_pct <- df %>%
    group_by(orig.ident, state) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(pct = 100 * n / sum(n))
  
  df_pct$state <- factor(df_pct$state, levels = names(fill_map))
  df_pct$orig.ident <- factor(df_pct$orig.ident, levels = sample_order)
  
  ggplot(df_pct, aes(orig.ident, pct, fill = state)) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_manual(values = fill_map, breaks = names(fill_map)) +
    labs(title = title_text, x = NULL, y = "% cells", fill = "State") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 1. Page 1: ALL Patients Merged (using global coordinates but all finalized states)
pdos_all <- pdos_matched # No filtering, keep finalized states

p1_tr <- make_umap_plot(pdos_all, "Treatment", c(Untreated = "#E69F00", Treated = "#4D4D4D"), "Merged Treatment")
p1_st <- make_umap_plot(pdos_all, "state", group_cols, "Merged State")
p1_pb <- make_prop_plot(pdos_all@meta.data, "All Patients Proportions", group_cols)

print(p1_tr + p1_st + p1_pb + plot_layout(ncol = 3))

# 2. Pages 2-5: Per Patient (RECALCULATED UMAP)
for (pat in patient_order) {
  message("Processing patient (UMAP recalculation): ", pat)
  # Subset cells for this patient
  cells_pat <- which(pdos_matched$Patient == pat)
  p_obj <- pdos_matched[, cells_pat]
  
  if (ncol(p_obj) < 30) {
    message("Too few cells for PCA/UMAP: ", pat)
    next
  }
  
  # Recalculate UMAP per patient
  p_obj <- FindVariableFeatures(p_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  p_obj <- ScaleData(p_obj, verbose = FALSE)
  p_obj <- RunPCA(p_obj, verbose = FALSE)
  
  # Check if PCA succeeded, then run UMAP
  p_obj <- RunUMAP(p_obj, dims = 1:15, verbose = FALSE)
  
  # Individual Treatment UMAP (new coords)
  p_p_tr <- make_umap_plot(p_obj, "Treatment", c(Untreated = "#E69F00", Treated = "#4D4D4D"), paste0(pat, " - Treatment (Re-calc)"))
  # Individual State UMAP (new coords)
  p_p_st <- make_umap_plot(p_obj, "state", group_cols, paste0(pat, " - State (Re-calc)"))
  # Individual Proportion Barplot
  p_p_pb <- make_prop_plot(p_obj@meta.data, paste0(pat, " - Proportions"), group_cols)
  
  print(p_p_tr + p_p_st + p_p_pb + plot_layout(ncol = 3))
}
dev.off()

# Keep high-res PNG of merged state UMAP (global)
p3_highres <- make_umap_plot(pdos_matched, "state", group_cols, "Finalized States (Matched)") +
  geom_point(size = 0.15, alpha = 0.8) 
png("Auto_PDO_matched_UMAP_state_highres.png", width = 2400, height = 1800, res = 300)
print(p3_highres)
dev.off()

# Existing MP expression plot can be updated to use finalized colors or kept as is
# I will keep the MP expression plot as a separate file for reference
message("=== MP expression bar plots (patient colors, gaps between patients) ===")
# ... (rest of the MP expression logic from previous script)
sample_treatment <- pdos_matched@meta.data %>% select(orig.ident, Treatment, Patient) %>% distinct()
mp_cols_plot <- colnames(ucell_matched)
mp_cols_plot <- mp_cols_plot[mp_cols_plot %in% names(mp_desc_map)]
mp_cols_plot <- mp_cols_plot[order(match(gsub(".*_", "", mp_cols_plot), mp_tree_order_names))]

mp_by_sample <- data.frame(
  sample = rep(pdos_matched$orig.ident, length(mp_cols_plot)),
  MP = rep(mp_cols_plot, each = length(pdos_matched$orig.ident)),
  score = as.vector(as.matrix(ucell_matched[, mp_cols_plot]))
) %>%
  group_by(sample, MP) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  left_join(sample_treatment, by = c("sample" = "orig.ident"))

mp_by_sample$MP_desc <- factor(mp_desc_map[mp_by_sample$MP], levels = mp_desc_map[mp_cols_plot])
mp_by_sample$color_col <- NA_character_
for (pat in patient_order) {
  base <- unname(patient_cols_base[pat])
  for (trt in c("Untreated", "Treated")) {
    idx <- mp_by_sample$Patient == pat & mp_by_sample$Treatment == trt
    mp_by_sample$color_col[idx] <- if (trt == "Untreated") base else lighten_col(base, amount = 0.55)
  }
}

sample_order_gap <- c(
  "SUR1070_Untreated_PDO", "SUR1070_Treated_PDO", "gap_1070",
  "SUR1072_Untreated_PDO", "SUR1072_Treated_PDO", "gap_1072",
  "SUR1090_Untreated_PDO", "SUR1090_Treated_PDO", "gap_1090",
  "SUR1181_Untreated_PDO", "SUR1181_Treated_PDO"
)
mp_by_sample$sample_plot <- factor(mp_by_sample$sample, levels = sample_order_gap)
fill_palette <- mp_by_sample %>% filter(!is.na(color_col)) %>% select(sample_plot, color_col) %>% distinct() %>% tibble::deframe()

# Add gaps
gap_df <- expand.grid(sample_plot = factor(c("gap_1070", "gap_1072", "gap_1090"), levels = sample_order_gap), MP_desc = levels(mp_by_sample$MP_desc))
gap_df$mean_score <- NA_real_

p_mp_expr <- ggplot(bind_rows(mp_by_sample, gap_df), aes(sample_plot, mean_score, fill = sample_plot)) +
  geom_col(color = "black", linewidth = 0.2, na.rm = TRUE) +
  facet_wrap(~ MP_desc, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = fill_palette, guide = "none", drop = FALSE) +
  scale_x_discrete(labels = function(x) ifelse(grepl("gap_", x), "", x), drop = FALSE) +
  labs(title = "Mean MP expression (Matched Samples)", x = NULL, y = "Mean UCell") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Auto_PDO_matched_MP_expression.pdf", width = 14, height = 10, useDingbats = FALSE)
print(p_mp_expr)
dev.off()

message("=== DONE ===")
