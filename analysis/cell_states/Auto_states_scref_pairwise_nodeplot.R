####################
# Auto_states_scref_pairwise_nodeplot.R
# Node plot for PDO cells using scRef-derived metaprogram state definitions.
# Style matches scRef_Pipeline's node plot exactly.
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(UCell)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
out_dir <- "hybrid_pairwise_scref"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# scRef Definition - matching scRef_Pipeline names exactly
real_states <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating"
)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive"       = c("MP13", "MP12"),
  "SMG-like Metaplasia"   = c("MP18", "MP16"),
  "Immune Infiltrating"   = c("MP15")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "Immune Infiltrating"   = "#377EB8",
  Unresolved = "grey80",
  Hybrid = "black"
)

message("=== Loading PDO data ===")
pdos <- readRDS("PDOs_final.rds")

# Batch/Study logic for normalization
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

message("=== Computing UCell scores for scRef MPs on PDOs ===")
# ncores=1 to avoid "subassignment" error on some login nodes
pdos <- AddModuleScore_UCell(pdos, features = mp.genes.sc, ncores = 1, name = "")

# scRef logic groups
cc_mps.sc <- c("MP1", "MP7", "MP9")
non_cc_mps.sc <- setdiff(retained_mps.sc, cc_mps.sc)

# Normalise
ucell_scores_scref <- pdos@meta.data[, retained_mps.sc, drop=FALSE]
Y_use <- ucell_scores_scref[, intersect(non_cc_mps.sc, colnames(ucell_scores_scref)), drop=FALSE]
mp_adj_noncc.sc <- z_normalise(Y_use, sample_var, study_var)

# Compute group max
group_max.sc <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_noncc.sc))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_noncc.sc[, mps_avail]))
  if (length(mps_avail) == 0) return(rep(0, nrow(mp_adj_noncc.sc)))
  apply(mp_adj_noncc.sc[, mps_avail, drop = FALSE], 1, max)
})
group_max.sc <- as.matrix(group_max.sc)
rownames(group_max.sc) <- rownames(mp_adj_noncc.sc)

# State assignment logic (matching scRef)
THRESHOLD <- 0.5
HYBRID_GAP_B <- 0.3

best_group_idx <- max.col(group_max.sc, ties.method = "first")
best_group_val <- apply(group_max.sc, 1, max)
base_state <- colnames(group_max.sc)[best_group_idx]
base_state[best_group_val < THRESHOLD] <- "Unresolved"

state_scref <- base_state
sorted_groups <- t(apply(group_max.sc, 1, sort, decreasing = TRUE))
gap <- sorted_groups[, 1] - sorted_groups[, 2]
state_scref[(gap < HYBRID_GAP_B) & (base_state != "Unresolved")] <- "Hybrid"
names(state_scref) <- rownames(group_max.sc)

# Hybrid Pairwise Assignment
hybrid_cells <- names(state_scref)[state_scref == "Hybrid"]
assign_pair <- function(x, names_vec) {
  ord <- names(sort(x, decreasing = TRUE))[1:2]
  # Maintain consistent order based on real_states vector
  ord <- names_vec[names_vec %in% ord]
  paste(ord, collapse = "__")
}

if (length(hybrid_cells) > 0) {
  pair_labels <- vapply(hybrid_cells, function(cl) assign_pair(group_max.sc[cl, real_states], real_states), character(1))
  pair_df <- data.frame(pair = pair_labels, stringsAsFactors = FALSE) %>%
    count(pair, name = "hybrid_cells") %>%
    separate(pair, into = c("from", "to"), sep = "__", remove = FALSE)
} else {
  pair_df <- data.frame(pair=character(), from=character(), to=character(), hybrid_cells=numeric(), stringsAsFactors=FALSE)
}

state_df <- data.frame(state = state_scref, stringsAsFactors = FALSE) %>%
  filter(state %in% real_states) %>%
  count(state, name = "cells")

tot_cells <- length(state_scref)
state_df <- state_df %>% mutate(pct = 100 * cells / tot_cells)
pair_df <- pair_df %>% mutate(pct = 100 * hybrid_cells / tot_cells)

# Node Layout
n <- length(real_states)
theta <- seq(0, 2 * pi, length.out = n + 1)[1:n]
layout_df <- data.frame(
  state = real_states,
  x = cos(theta),
  y = sin(theta),
  stringsAsFactors = FALSE
)

node_df <- left_join(layout_df, state_df, by = c("state" = "state"))
node_df$cells[is.na(node_df$cells)] <- 0
node_df$pct[is.na(node_df$pct)] <- 0
node_df$label_x <- node_df$x * 1.3
node_df$label_y <- node_df$y * 1.3

edge_df <- pair_df %>%
  left_join(layout_df, by = c("from" = "state")) %>%
  rename(x = x, y = y) %>%
  left_join(layout_df, by = c("to" = "state"), suffix = c("", "_to")) %>%
  rename(xend = x_to, yend = y_to)

p <- ggplot() +
  geom_segment(
    data = edge_df,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct),
    color = "grey35",
    alpha = 0.8
  ) +
  geom_point(
    data = node_df,
    aes(x = x, y = y, size = pct, color = state)
  ) +
  geom_text(
    data = node_df,
    aes(x = label_x, y = label_y, label = paste0(state, "\n", sprintf("%.1f%%", pct))),
    size = 3.2,
    fontface = "bold"
  ) +
  geom_label(
    data = edge_df,
    aes(
      x = (x + xend) / 2,
      y = (y + yend) / 2,
      label = sprintf("%.1f%%", pct)
    ),
    size = 2.6,
    fill = "white",
    label.size = 0,
    fontface = "bold"
  ) +
  scale_color_manual(values = group_cols) +
  scale_size(range = c(8, 22), guide = "none") +
  scale_linewidth(range = c(0.6, 6), guide = "none") +
  coord_equal() +
  expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
  theme_void(base_size = 14) +
  labs(title = "PDO Node Plot (scRef-derived states)") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(t = 10, b = 10)))

pdf(file.path(out_dir, "Auto_PDO_scref_hybrid_nodeplot.pdf"), width = 7, height = 7)
print(p)
dev.off()

# Heatmap
mat <- matrix(0, nrow = length(real_states), ncol = length(real_states),
              dimnames = list(real_states, real_states))
if (nrow(pair_df) > 0) {
  for (k in seq_len(nrow(pair_df))) {
    a <- pair_df$from[k]
    b <- pair_df$to[k]
    v <- pair_df$pct[k]
    if (a %in% real_states && b %in% real_states) {
      mat[a, b] <- v
      mat[b, a] <- v
    }
  }
}

hm_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
colnames(hm_df) <- c("StateA", "StateB", "Pct")
p_hm <- ggplot(hm_df, aes(StateB, StateA, fill = Pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Pct)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick3") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "PDO Pairwise hybrid heatmap (scRef-derived)",
       x = NULL, y = NULL, fill = "% of all cells")

pdf(file.path(out_dir, "Auto_PDO_scref_hybrid_heatmap.pdf"), width = 8, height = 6)
print(p_hm)
dev.off()

# Summary
summary_rows <- bind_rows(
  state_df %>% transmute(type = "state", label = state, cells = cells, pct = pct),
  pair_df %>% transmute(type = "pairwise_hybrid", label = pair, cells = hybrid_cells, pct = pct)
)
write.csv(summary_rows, file.path(out_dir, "Auto_PDO_scref_hybrid_summary.csv"), row.names = FALSE)

message("Saved scRef-derived node plot and summary to ", out_dir)
