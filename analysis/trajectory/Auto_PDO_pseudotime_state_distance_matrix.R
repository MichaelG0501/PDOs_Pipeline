####################
# Auto_PDO_pseudotime_state_distance_matrix.R
#
# Monocle3-based state-to-state distance comparison for PDO four-state,
# pre-unresolved-relabel trajectories. Distances are computed per sample and
# summarised across valid samples.
####################

suppressPackageStartupMessages({
  library(ggrepel)
})

helper_candidates <- c(
  "analysis/trajectory/Auto_PDO_pseudotime_helpers.R",
  "analysis/cell_states/Auto_PDO_pseudotime_helpers.R",
  "Auto_PDO_pseudotime_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Could not find Auto_PDO_pseudotime_helpers.R")
source(helper_path)

distance_dir <- file.path(pseudotime_dir, "state_distance_pseudotime")
dir.create(distance_dir, recursive = TRUE, showWarnings = FALSE)

build_heatmap_pdf <- function(long_df, file_path) {
  if (nrow(long_df) == 0) return(invisible(NULL))
  plot_df <- long_df %>%
    group_by(method) %>%
    mutate(
      method_range = max(distance, na.rm = TRUE) - min(distance, na.rm = TRUE),
      distance_scaled = ifelse(is.finite(method_range) & method_range > 0, (distance - min(distance, na.rm = TRUE)) / method_range, 0)
    ) %>%
    ungroup() %>%
    mutate(state_a = factor(state_a, levels = primary_states), state_b = factor(state_b, levels = primary_states))
  p <- ggplot(plot_df, aes(x = state_a, y = state_b, fill = distance_scaled)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = ifelse(is.na(distance), "NA", sprintf("%.2f", distance))), size = 3.5, fontface = "bold") +
    scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "grey90", name = "Within-method\nscaled distance") +
    facet_wrap(~method) +
    theme_minimal(base_size = 13) +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.text.y = element_text(face = "bold"), panel.grid = element_blank(), strip.text = element_text(face = "bold")) +
    labs(title = "PDO Four-State Distance Method Comparison", subtitle = "Tile colour is scaled within each method; labels show raw mean distances.")
  ggsave(file_path, p, width = 15, height = 10)
}

build_interconnected_network_pdf <- function(matrix_list, file_path) {
  plot_list <- list()
  for (method_name in names(matrix_list)) {
    mat <- matrix_list[[method_name]]
    mat_sub <- mat[primary_states, primary_states, drop = FALSE]
    mat_sub[is.na(mat_sub)] <- 0
    if (nrow(mat_sub) < 2) next
    set.seed(12345)
    mds_res <- tryCatch(cmdscale(as.dist(mat_sub), k = 2, eig = TRUE), error = function(e) NULL)
    if (is.null(mds_res) || is.null(mds_res$points)) next
    layout_df <- as.data.frame(mds_res$points)
    colnames(layout_df) <- c("x", "y")
    layout_df$state <- rownames(layout_df)
    layout_df$state <- factor(layout_df$state, levels = primary_states)
    edge_df <- as.data.frame(as.table(mat_sub), stringsAsFactors = FALSE) %>%
      rename(state_a = Var1, state_b = Var2, distance = Freq) %>%
      mutate(state_a = as.character(state_a), state_b = as.character(state_b)) %>%
      filter(state_a < state_b, distance > 0) %>%
      left_join(layout_df %>% select(state, x, y), by = c("state_a" = "state")) %>%
      left_join(layout_df %>% select(state, xend = x, yend = y), by = c("state_b" = "state")) %>%
      mutate(inv_dist = 1 / (distance + 1e-6), rel_inv_dist = (inv_dist - min(inv_dist)) / (max(inv_dist) - min(inv_dist) + 1e-6))
    p <- ggplot() +
      geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend, alpha = rel_inv_dist, linewidth = rel_inv_dist), colour = "grey70", lineend = "round") +
      geom_point(data = layout_df, aes(x = x, y = y, fill = state), shape = 21, size = 10, colour = "black", stroke = 1.1) +
      ggrepel::geom_text_repel(data = layout_df, aes(x = x, y = y, label = state), fontface = "bold", size = 5.2, box.padding = 0.6, point.padding = 0.4) +
      scale_fill_manual(values = state_cols) +
      scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
      scale_linewidth_continuous(range = c(1.2, 4.0), guide = "none") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), legend.position = "none") +
      labs(title = gsub("_", " ", method_name))
    plot_list[[method_name]] <- p
  }
  if (length(plot_list) > 0) {
    ggsave(file_path, wrap_plots(plot_list, ncol = 3) + plot_annotation(title = "PDO Four-State Distance Network"), width = 24, height = 14)
  }
}

message("=== PDO four-state distance matrix ===")
results <- get_or_build_pdo_sample_trajectories(rebuild = FALSE)
if (length(results) == 0) stop("No sample trajectories available.")

sample_summary_rows <- list()
directed_rows <- list()
geodesic_medoid_rows <- list()
geodesic_centroid_rows <- list()
umap_centroid_rows <- list()

for (sample_id in names(results)) {
  message("Computing distances for sample: ", sample_id)
  cds <- results[[sample_id]]$cds
  meta_out <- results[[sample_id]]$metadata %>% filter(state %in% primary_states)
  graph_bits <- tryCatch(extract_graph_structure(cds), error = function(e) NULL)
  if (is.null(graph_bits)) next
  graph_obj <- compute_graph_weights(graph_bits$graph, graph_bits$graph_coords)
  umap_mat <- reducedDims(cds)$UMAP
  umap_df <- data.frame(cell = rownames(umap_mat), UMAP_1 = umap_mat[, 1], UMAP_2 = umap_mat[, 2]) %>%
    left_join(meta_out %>% select(cell, state), by = "cell") %>%
    filter(state %in% primary_states)
  state_summary <- umap_df %>%
    group_by(state) %>%
    summarise(n_cells = n(), medoid_cell = get_state_medoid_cell(pick(cell, UMAP_1, UMAP_2, state)), centroid_x = mean(UMAP_1), centroid_y = mean(UMAP_2), .groups = "drop") %>%
    rowwise() %>%
    mutate(medoid_vertex = coerce_graph_vertex_name(graph_bits$closest_vertex[medoid_cell], graph_obj), centroid_vertex = nearest_graph_vertex(c(centroid_x, centroid_y), graph_bits$graph_coords)) %>%
    ungroup()
  sample_summary_rows[[length(sample_summary_rows) + 1]] <- state_summary %>% mutate(sample = sample_id, .before = 1)
  present_states <- intersect(primary_states, state_summary$state)
  for (state_a in present_states) {
    for (state_b in present_states) {
      row_a <- state_summary %>% filter(state == state_a)
      row_b <- state_summary %>% filter(state == state_b)
      medoid_distance <- suppressWarnings(igraph::distances(graph_obj, v = row_a$medoid_vertex[1], to = row_b$medoid_vertex[1], weights = igraph::E(graph_obj)$weight)[1, 1])
      centroid_distance <- suppressWarnings(igraph::distances(graph_obj, v = row_a$centroid_vertex[1], to = row_b$centroid_vertex[1], weights = igraph::E(graph_obj)$weight)[1, 1])
      euclidean_distance <- sqrt((row_a$centroid_x[1] - row_b$centroid_x[1]) ^ 2 + (row_a$centroid_y[1] - row_b$centroid_y[1]) ^ 2)
      geodesic_medoid_rows[[length(geodesic_medoid_rows) + 1]] <- data.frame(sample = sample_id, state_a = state_a, state_b = state_b, distance = if (is.finite(medoid_distance)) as.numeric(medoid_distance) else NA_real_)
      geodesic_centroid_rows[[length(geodesic_centroid_rows) + 1]] <- data.frame(sample = sample_id, state_a = state_a, state_b = state_b, distance = if (is.finite(centroid_distance)) as.numeric(centroid_distance) else NA_real_)
      umap_centroid_rows[[length(umap_centroid_rows) + 1]] <- data.frame(sample = sample_id, state_a = state_a, state_b = state_b, distance = euclidean_distance)
    }
  }
  root_candidates <- present_states[state_summary$n_cells[match(present_states, state_summary$state)] >= ROOT_MIN_CELLS]
  for (root_state_i in root_candidates) {
    root_cells <- meta_out$cell[meta_out$state == root_state_i]
    root_node <- get_root_pr_node(cds, root_cells)
    if (is.null(root_node)) next
    ordered_cds <- tryCatch(order_cells(cds, root_pr_nodes = root_node), error = function(e) NULL)
    if (is.null(ordered_cds)) next
    pt_df <- data.frame(cell = names(pseudotime(ordered_cds)), pseudotime = as.numeric(pseudotime(ordered_cds))) %>%
      left_join(meta_out %>% select(cell, state), by = "cell")
    state_pt <- pt_df %>%
      group_by(state) %>%
      summarise(median_pseudotime = safe_median(pseudotime), mean_pseudotime = safe_mean(pseudotime), n_cells_with_pt = sum(is.finite(pseudotime)), .groups = "drop")
    root_median <- state_pt$median_pseudotime[state_pt$state == root_state_i]
    root_mean <- state_pt$mean_pseudotime[state_pt$state == root_state_i]
    for (target_state in primary_states) {
      target_row <- state_pt %>% filter(state == target_state)
      if (nrow(target_row) == 0) next
      directed_rows[[length(directed_rows) + 1]] <- data.frame(sample = sample_id, root_state = root_state_i, target_state = target_state, root_state_cells = sum(meta_out$state == root_state_i), target_state_cells = target_row$n_cells_with_pt[1], median_distance = target_row$median_pseudotime[1] - root_median[1], mean_distance = target_row$mean_pseudotime[1] - root_mean[1])
    }
  }
}

sample_summary_df <- bind_rows(sample_summary_rows)
directed_df <- bind_rows(directed_rows)
geodesic_medoid_df <- bind_rows(geodesic_medoid_rows)
geodesic_centroid_df <- bind_rows(geodesic_centroid_rows)
umap_centroid_df <- bind_rows(umap_centroid_rows)
write_csv(sample_summary_df, file.path(distance_dir, "Auto_PDO_state_distance_sample_summary.csv"))
write_csv(directed_df, file.path(distance_dir, "Auto_PDO_state_distance_directed_pseudotime.csv"))
write_csv(geodesic_medoid_df, file.path(distance_dir, "Auto_PDO_state_distance_geodesic_medoid.csv"))
write_csv(geodesic_centroid_df, file.path(distance_dir, "Auto_PDO_state_distance_geodesic_centroid.csv"))
write_csv(umap_centroid_df, file.path(distance_dir, "Auto_PDO_state_distance_umap_centroid.csv"))

directed_pairwise_median <- directed_df %>% filter(root_state != target_state) %>% mutate(state_a = pmin(root_state, target_state), state_b = pmax(root_state, target_state), distance = abs(median_distance)) %>% group_by(sample, state_a, state_b) %>% summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop")
directed_pairwise_mean <- directed_df %>% filter(root_state != target_state) %>% mutate(state_a = pmin(root_state, target_state), state_b = pmax(root_state, target_state), distance = abs(mean_distance)) %>% group_by(sample, state_a, state_b) %>% summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop")
geodesic_medoid_sym <- geodesic_medoid_df %>% mutate(state_min = pmin(state_a, state_b), state_max = pmax(state_a, state_b)) %>% group_by(sample, state_min, state_max) %>% summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>% rename(state_a = state_min, state_b = state_max)
geodesic_centroid_sym <- geodesic_centroid_df %>% mutate(state_min = pmin(state_a, state_b), state_max = pmax(state_a, state_b)) %>% group_by(sample, state_min, state_max) %>% summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>% rename(state_a = state_min, state_b = state_max)
umap_centroid_sym <- umap_centroid_df %>% mutate(state_min = pmin(state_a, state_b), state_max = pmax(state_a, state_b)) %>% group_by(sample, state_min, state_max) %>% summarise(distance = mean(distance, na.rm = TRUE), .groups = "drop") %>% rename(state_a = state_min, state_b = state_max)

summary_df <- bind_rows(
  directed_pairwise_median %>% group_by(state_a, state_b) %>% summarise(method = "directed_pseudotime_median", mean_distance = safe_mean(distance), median_distance = safe_median(distance), n_samples_used = sum(is.finite(distance)), .groups = "drop"),
  directed_pairwise_mean %>% group_by(state_a, state_b) %>% summarise(method = "directed_pseudotime_mean", mean_distance = safe_mean(distance), median_distance = safe_median(distance), n_samples_used = sum(is.finite(distance)), .groups = "drop"),
  geodesic_medoid_sym %>% group_by(state_a, state_b) %>% summarise(method = "principal_graph_geodesic_medoid", mean_distance = safe_mean(distance), median_distance = safe_median(distance), n_samples_used = sum(is.finite(distance)), .groups = "drop"),
  geodesic_centroid_sym %>% group_by(state_a, state_b) %>% summarise(method = "principal_graph_geodesic_centroid", mean_distance = safe_mean(distance), median_distance = safe_median(distance), n_samples_used = sum(is.finite(distance)), .groups = "drop"),
  umap_centroid_sym %>% group_by(state_a, state_b) %>% summarise(method = "umap_centroid_euclidean", mean_distance = safe_mean(distance), median_distance = safe_median(distance), n_samples_used = sum(is.finite(distance)), .groups = "drop")
) %>% arrange(method, state_a, state_b)
write_csv(summary_df, file.path(distance_dir, "Auto_PDO_state_distance_summary.csv"))
write_csv(summary_df, file.path(summary_dir, "Auto_PDO_state_distance_summary.csv"))

matrix_list <- list(
  directed_pseudotime_median = make_symmetric_matrix(summary_df, "directed_pseudotime_median", primary_states),
  directed_pseudotime_mean = make_symmetric_matrix(summary_df, "directed_pseudotime_mean", primary_states),
  principal_graph_geodesic_medoid = make_symmetric_matrix(summary_df, "principal_graph_geodesic_medoid", primary_states),
  principal_graph_geodesic_centroid = make_symmetric_matrix(summary_df, "principal_graph_geodesic_centroid", primary_states),
  umap_centroid_euclidean = make_symmetric_matrix(summary_df, "umap_centroid_euclidean", primary_states)
)
saveRDS(matrix_list, file.path(distance_dir, "Auto_PDO_state_distance_matrices.rds"))
for (method_name in names(matrix_list)) {
  write.csv(matrix_list[[method_name]], file.path(distance_dir, paste0("Auto_PDO_", method_name, "_state_matrix.csv")), quote = FALSE)
}
heatmap_df <- bind_rows(lapply(names(matrix_list), function(method_name) matrix_to_long(matrix_list[[method_name]], method_name)))
write_csv(heatmap_df, file.path(distance_dir, "Auto_PDO_state_distance_long.csv"))
build_heatmap_pdf(heatmap_df, file.path(distance_dir, "Auto_PDO_state_distance_method_comparison_heatmap.pdf"))
build_interconnected_network_pdf(matrix_list, file.path(distance_dir, "Auto_PDO_state_distance_nodeplot.pdf"))
message("Saved PDO state distance outputs to: ", distance_dir)
