####################
# Auto_PDO_state_concordance.R
# Computes scRef (nMP=19) states on PDO cells using noreg approach B,
# then compares with PDO's native states computed with the same approach.
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(UCell)
library(reshape2)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

message("=== Loading PDO data ===")
pdos <- readRDS("PDOs_final.rds")

# Fix batch logic as in PDO_states_analysis.R
pdos$Batch <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")
sample_var <- pdos$orig.ident
study_var <- pdos$Batch

message("=== Loading scRef MPs ===")
sc_mp_obj <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")

# Filter scRef MPs (Silhouette < 0)
mp.genes.sc <- sc_mp_obj$metaprograms.genes
bad_mps.sc <- which(sc_mp_obj$metaprograms.metrics$silhouette < 0)
if (length(bad_mps.sc) > 0) {
  mp.genes.sc <- mp.genes.sc[!names(mp.genes.sc) %in% paste0("MP", bad_mps.sc)]
}
retained_mps.sc <- names(mp.genes.sc)

message("Retained scRef MPs: ", paste(retained_mps.sc, collapse=", "))

message("=== Computing UCell scores for scRef MPs on PDOs ===")
pdos <- AddModuleScore_UCell(pdos, features = mp.genes.sc, ncores = 1, name = "")

# Extract UCell scores
sc_mp_cols <- retained_mps.sc
ucell_scores_scref <- pdos@meta.data[, sc_mp_cols, drop=FALSE]
common_cells <- rownames(ucell_scores_scref)

message("=== Defining scRef State logic ===")
# Exact scRef logic
cc_mps.sc <- c("MP1", "MP7", "MP9")
non_cc_mps.sc <- setdiff(retained_mps.sc, cc_mps.sc)

state_groups.sc <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intest. Meta" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
)
state_groups.sc <- lapply(state_groups.sc, function(x) intersect(x, retained_mps.sc))
state_groups.sc <- state_groups.sc[sapply(state_groups.sc, length) > 0]

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

# noreg mode setup
Y_use <- ucell_scores_scref[, intersect(non_cc_mps.sc, colnames(ucell_scores_scref)), drop=FALSE]
mp_adj_noncc.sc <- z_normalise(Y_use, sample_var, study_var)

group_max.sc <- sapply(state_groups.sc, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_noncc.sc))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc.sc[, mps_avail]))
  if (length(mps_avail) == 0) return(rep(0, nrow(mp_adj_noncc.sc)))
  apply(mp_adj_noncc.sc[, mps_avail, drop = FALSE], 1, max)
})
group_max.sc <- as.matrix(group_max.sc)
rownames(group_max.sc) <- rownames(mp_adj_noncc.sc)

# thresholds
THRESHOLD <- 0.5
HYBRID_GAP_B <- 0.3

best_group_idx <- max.col(group_max.sc, ties.method = "first")
best_group_val <- apply(group_max.sc, 1, max)
base_state <- colnames(group_max.sc)[best_group_idx]
base_state[best_group_val < THRESHOLD] <- "Unresolved"

sorted_groups <- t(apply(group_max.sc, 1, sort, decreasing = TRUE))
gap <- sorted_groups[, 1] - sorted_groups[, 2]
state_scref_on_pdo <- base_state
state_scref_on_pdo[(gap < HYBRID_GAP_B) & (base_state != "Unresolved")] <- "Hybrid"
names(state_scref_on_pdo) <- rownames(group_max.sc)

# rename "Basal to Intestinal Metaplasia" to "Basal to Intest. Meta" to match PDO for plotting brevity
state_scref_on_pdo[state_scref_on_pdo == "Basal to Intestinal Metaplasia"] <- "Basal to Intest. Meta"

# Save scref states
saveRDS(state_scref_on_pdo, "Auto_PDO_scref_states.rds")

message("=== Loading PDO States ===")
state_pdo_on_pdo <- readRDS("Auto_PDO_states_noreg.rds")

# Ensure common cells
cells <- intersect(names(state_scref_on_pdo), names(state_pdo_on_pdo))
s_sc <- state_scref_on_pdo[cells]
s_pdo <- state_pdo_on_pdo[cells]

message("=== Concordance Analysis ===")

df <- data.frame(
  cell = cells,
  PDO_State = s_pdo,
  scRef_State = s_sc,
  stringsAsFactors = FALSE
)

# Confusion Matrix
tbl <- table(PDO_State = df$PDO_State, scRef_State = df$scRef_State)
tbl_pct <- t(apply(tbl, 1, function(x) x / sum(x) * 100))

# Custom order for state levels
state_levels <- c("Classic Proliferative", "Basal to Intest. Meta", "Stress-adaptive", "SMG-like Metaplasia", "Immune Infiltrating", "Unresolved", "Hybrid")

ord_pdo <- intersect(state_levels, rownames(tbl))
ord_sc <- intersect(state_levels, colnames(tbl))

tbl <- tbl[ord_pdo, ord_sc]
tbl_pct <- tbl_pct[ord_pdo, ord_sc]

pdf("Auto_PDO_concordance_heatmap.pdf", width = 9, height = 7, useDingbats = FALSE)
col_fun <- colorRamp2(
  c(0, 50, 100),
  c("#FFFFFF", "#E31A1C", "#7F0000")
)

ht <- Heatmap(tbl_pct, name = "Overlap (% of PDO state)",
        cluster_rows = FALSE, cluster_columns = FALSE,
        col = col_fun,
        rect_gp = gpar(col = "white", lwd=1),
        row_names_side = "left",
        row_title = "PDO native state (assigned using PDO MPs)",
        column_title = "scRef state (assigned using scRef MPs)",
        column_title_side = "bottom",
        cell_fun = function(j, i, x, y, width, height, fill) {
          lab <- sprintf("%.1f%%\n(N=%d)", tbl_pct[i, j], tbl[i, j])
          grid.text(lab, x, y, gp = gpar(fontsize = 10, col = ifelse(tbl_pct[i, j] > 50, "white", "black")))
        })
draw(ht)
dev.off()

# Stacked Barplot (Normalized)
df_plot <- as.data.frame(tbl)
colnames(df_plot) <- c("PDO_State", "scRef_State", "Count")
df_plot <- df_plot %>%
  group_by(PDO_State) %>%
  mutate(Pct = Count / sum(Count) * 100) %>%
  ungroup()

df_plot$PDO_State <- factor(df_plot$PDO_State, levels = ord_pdo)
df_plot$scRef_State <- factor(df_plot$scRef_State, levels = rev(ord_sc))

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

p <- ggplot(df_plot, aes(x = PDO_State, y = Pct, fill = scRef_State)) +
  geom_bar(stat = "identity", color="black", linewidth=0.3) +
  scale_fill_manual(values = group_cols, drop=FALSE) +
  labs(title = "Concordance of State Definition Approaches (noreg B) on PDOs",
       x = "PDO Native State", y = "Percentage of Cells", fill = "Assigned scRef State") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Auto_PDO_concordance_barplot.pdf", p, width = 9, height = 7, useDingbats = FALSE)

message("=== Checking if ggalluvial is available for alluvial plot ===")
has_alluvial <- require(ggalluvial, quietly = TRUE)
if (has_alluvial) {
  df$PDO_State <- factor(df$PDO_State, levels=ord_pdo)
  df$scRef_State <- factor(df$scRef_State, levels=ord_sc)
  
  # Prepare for ggalluvial
  df_aggr <- df %>% count(PDO_State, scRef_State)
  
  p_alluv <- ggplot(df_aggr,
                    aes(axis1 = PDO_State, axis2 = scRef_State, y = n)) +
    geom_alluvium(aes(fill = PDO_State), width = 1/12) +
    geom_stratum(width = 1/4, fill = "grey90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_fill_manual(values = group_cols) +
    scale_x_discrete(limits = c("PDO Native State", "scRef State"), expand = c(.1, .1)) +
    labs(title = "State transitions (PDO MPs vs scRef MPs on PDO cells)", y = "Number of Cells") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", panel.grid = element_blank())
    
  ggsave("Auto_PDO_concordance_alluvial.pdf", p_alluv, width = 10, height = 7, useDingbats = FALSE)
  message("Saved alluvial plot.")
} else {
  message("ggalluvial not installed, skipping alluvial plot.")
}

message("=== SUCCESS ===")
