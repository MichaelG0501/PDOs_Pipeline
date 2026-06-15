####################
# Auto_PDO_pseudotime_linear_plot.R
#
# PDO per-sample trajectory reports matching the Parse/scRef pseudotime report
# style. Uses cached per-sample Monocle3 trajectories or builds them if missing.
####################

suppressPackageStartupMessages({
  library(ggridges)
})

helper_candidates <- c(
  "analysis/trajectory/Auto_PDO_pseudotime_helpers.R",
  "analysis/cell_states/Auto_PDO_pseudotime_helpers.R",
  "Auto_PDO_pseudotime_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Could not find Auto_PDO_pseudotime_helpers.R")
source(helper_path)

message("=== PDO per-sample pseudotime trajectory reports ===")
results <- get_or_build_pdo_sample_trajectories(rebuild = FALSE)
if (length(results) == 0) stop("No sample trajectories available.")

pdf_path <- file.path(pseudotime_dir, "Auto_PDO_pseudotime_linear_reports.pdf")
pdf(pdf_path, width = 16, height = 12, onefile = TRUE)
for (sample_id in names(results)) {
  message("Processing report for sample: ", sample_id)
  cds <- results[[sample_id]]$cds
  proj_df <- results[[sample_id]]$projections %>% mutate(state = factor(as.character(state), levels = primary_states))
  graph_bits <- extract_graph_structure(cds)
  root_cells <- proj_df$cell[as.character(proj_df$state) == root_state]
  root_label_df <- get_root_label_position(cds, graph_bits, root_cells)
  ridge_bits <- build_weighted_ridge_df(proj_df, primary_states)
  state_counts <- table(factor(proj_df$state, levels = primary_states))
  legend_labels <- setNames(paste0(primary_states, " (", as.integer(state_counts[primary_states]), ")"), primary_states)
  x_span <- diff(range(c(graph_bits$nodes$x, proj_df$graph_x), na.rm = TRUE))
  y_span <- diff(range(c(graph_bits$nodes$y, proj_df$graph_y), na.rm = TRUE))
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  root_n <- sum(as.character(proj_df$state) == root_state, na.rm = TRUE)

  p_proj_state <- ggplot() +
    geom_segment(data = graph_bits$edges, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.9, color = "grey55", alpha = 0.65) +
    geom_point(data = proj_df, aes(x = graph_x, y = graph_y, color = state), size = 1.0, alpha = 0.85) +
    { if (!is.null(root_label_df)) geom_point(data = root_label_df, aes(x = x, y = y), shape = 8, size = 4.2, stroke = 1.0, color = "black") } +
    { if (!is.null(root_label_df)) geom_label(data = root_label_df, aes(x = x + 0.03 * x_span, y = y + 0.04 * y_span, label = paste0("ROOT\n", root_state, "\n(n=", root_n, ")")), size = 3.2, label.size = 0.25, hjust = 0, vjust = 0, fill = "white", color = "black", alpha = 0.95) } +
    scale_color_manual(values = state_cols, breaks = primary_states, labels = legend_labels, drop = FALSE) +
    coord_equal() +
    theme_classic() +
    labs(title = "Principal Graph: state-projected cells", x = "Dim 1", y = "Dim 2", color = "State")

  p_proj_pt <- ggplot() +
    geom_segment(data = graph_bits$edges, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.9, color = "grey55", alpha = 0.65) +
    geom_point(data = proj_df, aes(x = graph_x, y = graph_y, color = pseudotime), size = 1.0, alpha = 0.85) +
    { if (!is.null(root_label_df)) geom_point(data = root_label_df, aes(x = x, y = y), shape = 8, size = 4.2, stroke = 1.0, color = "black") } +
    { if (!is.null(root_label_df)) geom_label(data = root_label_df, aes(x = x + 0.03 * x_span, y = y + 0.04 * y_span, label = "ROOT"), size = 3.2, label.size = 0.25, hjust = 0, vjust = 0, fill = "white", color = "black", alpha = 0.95) } +
    scale_color_viridis_c(option = "D", na.value = "grey85") +
    coord_equal() +
    theme_classic() +
    labs(title = "Principal Graph: pseudotime-projected cells", x = "Dim 1", y = "Dim 2", color = "Pseudotime")

  p_umap_pt <- ggplot(proj_df %>% filter(is.finite(umap_x), is.finite(umap_y)), aes(x = umap_x, y = umap_y, color = pseudotime)) +
    geom_point(size = 0.9, alpha = 0.85) +
    scale_color_viridis_c(option = "D", na.value = "grey85") +
    coord_equal() +
    theme_classic() +
    labs(title = "UMAP: pseudotime", x = "UMAP 1", y = "UMAP 2", color = "Pseudotime")

  if (!is.null(ridge_bits) && nrow(ridge_bits$ridge_df) > 0) {
    p_ridges <- ggplot(ridge_bits$ridge_df, aes(x = x, y = y, height = height, group = state, fill = state)) +
      geom_ridgeline(stat = "identity", scale = 1, alpha = 0.85, color = "white", linewidth = 0.35) +
      scale_fill_manual(values = state_cols, drop = FALSE) +
      coord_cartesian(xlim = ridge_bits$focus_range) +
      scale_y_continuous(breaks = ridge_bits$state_labels$y, labels = ridge_bits$state_labels$label, expand = expansion(mult = c(0.03, 0.1))) +
      theme_classic() +
      theme(legend.position = "none", axis.title.y = element_blank()) +
      labs(title = "Pseudotime density ridges by state (1st-99th percentile focus)", x = "Pseudotime")
  } else {
    p_ridges <- ggplot() + theme_void() + labs(title = "Pseudotime density ridges") + annotate("text", x = 0.5, y = 0.5, label = "Insufficient finite pseudotime values", size = 5)
  }

  print((p_proj_state + p_proj_pt) / (p_umap_pt + p_ridges) +
    plot_annotation(
      title = paste0("PDO Trajectory Report: ", sample_id),
      subtitle = paste0("Total cells: ", nrow(proj_df), " | Root state: ", root_state, " | Root cells: ", root_n)
    ))
}
dev.off()
message("Saved PDO trajectory reports to: ", pdf_path)
