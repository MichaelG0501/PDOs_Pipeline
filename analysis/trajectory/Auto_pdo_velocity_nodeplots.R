####################
# Auto_pdo_velocity_nodeplots.R
#
# Directed node plots for PDO RNA velocity state transitions.
# Nodes are four-state abundances before unresolved relabeling; arrows are
# positive velocity alignment from source-state mean velocity toward target
# state centroids.
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(data.table)
})

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_velocity_PDO"
node_path <- file.path(out_dir, "tables", "Auto_pdo_velocity_state_nodes.csv")
edge_path <- file.path(out_dir, "tables", "Auto_pdo_velocity_state_direction_edges.csv")
if (!file.exists(node_path) || !file.exists(edge_path)) {
  stop("Missing scVelo state transition tables.")
}

nodes <- fread(node_path)
edges <- fread(edge_path)
if (!"pct_total_pre_relabel" %in% names(nodes)) {
  nodes$pct_total_pre_relabel <- nodes$pct_major
}

state_levels <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive"
)
state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3"
)
layout_df <- data.frame(
  state = state_levels,
  x = c(-1, 1, 1, -1),
  y = c(0.72, 0.72, -0.72, -0.72),
  stringsAsFactors = FALSE
)

nodes <- nodes %>%
  mutate(
    state = as.character(state),
    patient = sub("^(SUR[0-9]+).*", "\\1", sample),
    pct_plot = pct_total_pre_relabel
  ) %>%
  left_join(layout_df, by = c("state" = "state")) %>%
  mutate(state = factor(state, levels = state_levels))

edges <- edges %>%
  mutate(
    source = as.character(source),
    target = as.character(target),
    patient = sub("^(SUR[0-9]+).*", "\\1", sample),
    plotted = velocity_alignment > 0.10
  ) %>%
  left_join(layout_df, by = c("source" = "state")) %>%
  left_join(layout_df, by = c("target" = "state"), suffix = c("", "_to")) %>%
  rename(xend = x_to, yend = y_to) %>%
  mutate(
    source = factor(source, levels = state_levels),
    target = factor(target, levels = state_levels)
  )

combine_nodes <- function(df) {
  df %>%
    group_by(state, x, y) %>%
    summarise(
      cells = median(cells, na.rm = TRUE),
      pct_major = median(pct_major, na.rm = TRUE),
      pct_plot = median(pct_plot, na.rm = TRUE),
      .groups = "drop"
    )
}

combine_edges <- function(df) {
  df %>%
    group_by(source, target, x, y, xend, yend) %>%
    summarise(
      velocity_alignment = mean(velocity_alignment, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(velocity_alignment > 0.10)
}

draw_nodeplot <- function(ndf, edf, title_text, max_node, max_edge) {
  ndf <- ndf %>%
    mutate(
      label = dplyr::recode(
        as.character(state),
        "Classic Proliferative" = "Classic\nProlif.",
        "Basal to Intest. Meta" = "Basal to\nIntest. Meta",
        "SMG-like Metaplasia" = "SMG-like\nMetaplasia",
        "Stress-adaptive" = "Stress-\nadaptive"
      ),
      label = paste0(label, "\n", sprintf("%.1f%%", pct_plot))
    )
  edf <- edf %>%
    filter(velocity_alignment > 0.10) %>%
    mutate(
      edge_dx = xend - x,
      edge_dy = yend - y,
      edge_len = sqrt(edge_dx^2 + edge_dy^2),
      edge_ux = ifelse(edge_len > 0, edge_dx / edge_len, 0),
      edge_uy = ifelse(edge_len > 0, edge_dy / edge_len, 0),
      source_idx = match(as.character(source), state_levels),
      target_idx = match(as.character(target), state_levels),
      curve_side = ifelse(source_idx < target_idx, 1, -1),
      perp_x = -edge_uy,
      perp_y = edge_ux,
      x_plot = x + 0.28 * edge_ux,
      y_plot = y + 0.28 * edge_uy,
      xend_plot = xend - 0.34 * edge_ux,
      yend_plot = yend - 0.34 * edge_uy,
      label_x = (x + xend) / 2 + 0.22 * curve_side * perp_x,
      label_y = (y + yend) / 2 + 0.22 * curve_side * perp_y
    )

  ggplot() +
    {
      if (nrow(edf) > 0) {
        geom_curve(
          data = edf,
          aes(x = x_plot, y = y_plot, xend = xend_plot, yend = yend_plot, linewidth = velocity_alignment),
          curvature = 0.28,
          color = "grey10",
          alpha = 0.95,
          arrow = arrow(length = unit(0.28, "inches"), type = "closed"),
          lineend = "round"
        )
      }
    } +
    {
      if (nrow(edf) > 0) {
        geom_label(
          data = edf,
          aes(x = label_x, y = label_y, label = sprintf("%.2f", velocity_alignment)),
          size = 3.8,
          fill = "white",
          linewidth = 0.2,
          fontface = "bold"
        )
      }
    } +
    geom_point(data = ndf, aes(x = x, y = y, size = pct_plot, color = state)) +
    geom_text(data = ndf, aes(x = x, y = y, label = label), size = 3.8, fontface = "bold", lineheight = 0.86, color = "black") +
    scale_color_manual(values = state_cols, drop = FALSE) +
    scale_size(limits = c(0, max_node), range = c(12, 34), guide = "none") +
    scale_linewidth(limits = c(0, max_edge), range = c(1.0, 6.5), guide = "none") +
    coord_equal() +
    expand_limits(x = c(-1.55, 1.55), y = c(-1.22, 1.22)) +
    theme_void(base_size = 14) +
    labs(title = title_text) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
}

safe_max <- function(x, fallback = 1) {
  out <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(out) || out <= 0) fallback else out
}

pdf_path <- file.path(out_dir, "figures", "Auto_pdo_velocity_nodeplot_untreated_vs_treated.pdf")
dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
message("Writing: ", pdf_path)
pdf(pdf_path, width = 14, height = 7)

untreated_nodes <- combine_nodes(nodes %>% filter(treatment == "Untreated"))
treated_nodes <- combine_nodes(nodes %>% filter(treatment == "Treated"))
untreated_edges <- combine_edges(edges %>% filter(treatment == "Untreated"))
treated_edges <- combine_edges(edges %>% filter(treatment == "Treated"))
max_node_comb <- safe_max(c(untreated_nodes$pct_plot, treated_nodes$pct_plot))
max_edge_comb <- safe_max(c(untreated_edges$velocity_alignment, treated_edges$velocity_alignment), 0.2)

p_comb <- draw_nodeplot(untreated_nodes, untreated_edges, "Untreated (median nodes, mean arrows)", max_node_comb, max_edge_comb) |
  draw_nodeplot(treated_nodes, treated_edges, "Treated (median nodes, mean arrows)", max_node_comb, max_edge_comb)
print(p_comb + plot_annotation(
  title = "New Batch PDO Velocity State Direction"
))

paired_patients <- intersect(
  unique(nodes$patient[nodes$treatment == "Untreated"]),
  unique(nodes$patient[nodes$treatment == "Treated"])
)
paired_patients <- sort(paired_patients)
for (pid in paired_patients) {
  n_ut <- nodes %>% filter(patient == pid, treatment == "Untreated")
  n_tr <- nodes %>% filter(patient == pid, treatment == "Treated")
  e_ut <- edges %>% filter(patient == pid, treatment == "Untreated", plotted)
  e_tr <- edges %>% filter(patient == pid, treatment == "Treated", plotted)
  max_node_pt <- safe_max(c(n_ut$pct_plot, n_tr$pct_plot))
  max_edge_pt <- safe_max(c(e_ut$velocity_alignment, e_tr$velocity_alignment), 0.2)
  pp <- draw_nodeplot(n_ut, e_ut, "Untreated", max_node_pt, max_edge_pt) |
    draw_nodeplot(n_tr, e_tr, "Treated", max_node_pt, max_edge_pt)
  print(pp + plot_annotation(
    title = paste0(pid, " PDO velocity state direction")
  ))
}

cynthia_samples <- sort(unique(nodes$sample[nodes$batch_type == "Cynthia_batch"]))
for (sample_id in cynthia_samples) {
  n_sm <- nodes %>% filter(sample == sample_id)
  e_sm <- edges %>% filter(sample == sample_id, plotted)
  max_node_sm <- safe_max(n_sm$pct_plot)
  max_edge_sm <- safe_max(e_sm$velocity_alignment, 0.2)
  print(draw_nodeplot(n_sm, e_sm, sample_id, max_node_sm, max_edge_sm))
}

dev.off()

message("Node plot PDF done.")
