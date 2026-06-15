####################
# Shared helpers for PDO pre-relabel four-state Monocle3 pseudotime workflows.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(igraph)
  library(patchwork)
  library(scales)
})

set.seed(12345)

get_pdo_root_dir <- function() {
  env_root <- Sys.getenv("AUTO_PDO_ROOT_DIR", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = FALSE))
  candidates <- c(
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
  )
  candidates[file.exists(candidates)][1]
}

root_dir <- get_pdo_root_dir()
setwd(root_dir)
qc_dir <- file.path(root_dir, "PDOs_outs")
pseudotime_dir <- file.path(qc_dir, "Auto_PDO_pseudotime_pre_relabel")
asset_dir <- file.path(pseudotime_dir, "sample_trajectory_assets")
summary_dir <- file.path(qc_dir, "summary")
dir.create(pseudotime_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

primary_states <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive"
)

root_state <- Sys.getenv("AUTO_PDO_ROOT_STATE", unset = "Basal to Intest. Meta")

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3"
)

ROOT_MIN_CELLS <- as.integer(Sys.getenv("AUTO_PDO_PSEUDOTIME_ROOT_MIN_CELLS", unset = "20"))
OTHER_MIN_CELLS <- as.integer(Sys.getenv("AUTO_PDO_PSEUDOTIME_OTHER_MIN_CELLS", unset = "20"))
TOTAL_MIN_CELLS <- as.integer(Sys.getenv("AUTO_PDO_PSEUDOTIME_TOTAL_MIN_CELLS", unset = "80"))
MIN_STATES_OVER_THRESHOLD <- as.integer(Sys.getenv("AUTO_PDO_PSEUDOTIME_MIN_STATES", unset = "2"))

load_pdo_object <- function() {
  obj_path <- file.path(qc_dir, "PDOs_merged.rds")
  if (!file.exists(obj_path)) stop("Missing PDO Seurat object: ", obj_path)
  readRDS(obj_path)
}

load_pre_relabel_states <- function() {
  state_path <- file.path(qc_dir, "Auto_PDO_states_noreg.rds")
  if (!file.exists(state_path)) stop("Missing pre-relabel PDO state vector: ", state_path)
  state_vec <- readRDS(state_path)
  state_names <- names(state_vec)
  state_vec <- as.character(state_vec)
  names(state_vec) <- state_names
  state_vec
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

sanitize_sample_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

make_empty_matrix <- function(labels) {
  mat <- matrix(NA_real_, nrow = length(labels), ncol = length(labels))
  rownames(mat) <- labels
  colnames(mat) <- labels
  diag(mat) <- 0
  mat
}

matrix_to_long <- function(mat, method_name) {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    rename(state_a = Var1, state_b = Var2, distance = Freq) %>%
    mutate(method = method_name)
}

make_symmetric_matrix <- function(summary_df, method_name, states) {
  out <- make_empty_matrix(states)
  df_use <- summary_df %>%
    filter(method == method_name) %>%
    mutate(state_min = pmin(state_a, state_b), state_max = pmax(state_a, state_b)) %>%
    group_by(state_min, state_max) %>%
    summarise(distance = mean(mean_distance, na.rm = TRUE), .groups = "drop")
  if (nrow(df_use) == 0) return(out)
  for (i in seq_len(nrow(df_use))) {
    out[df_use$state_min[i], df_use$state_max[i]] <- df_use$distance[i]
    out[df_use$state_max[i], df_use$state_min[i]] <- df_use$distance[i]
  }
  diag(out) <- 0
  out
}

prepare_pdo_sample_trajectory <- function(seurat_obj) {
  if (inherits(seurat_obj[["RNA"]], "Assay5")) {
    seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
  }
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  n_pcs <- min(30, ncol(seurat_obj) - 1)
  if (n_pcs < 2) stop("Too few cells for PCA.")
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = n_pcs, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:min(15, n_pcs), verbose = FALSE)
  cds <- as.cell_data_set(seurat_obj)
  cds <- cluster_cells(cds, verbose = FALSE)
  learn_graph(cds, verbose = FALSE, use_partition = FALSE)
}

flatten_closest_vertex <- function(closest_vertex) {
  if (is.null(closest_vertex)) return(character())
  if (is.matrix(closest_vertex) || is.data.frame(closest_vertex)) {
    out <- as.character(closest_vertex[, 1, drop = TRUE])
    names(out) <- rownames(closest_vertex)
    return(out)
  }
  as.character(closest_vertex)
}

coerce_graph_vertex_name <- function(raw_vertex, graph_obj) {
  graph_nodes <- igraph::V(graph_obj)$name
  raw_vertex <- as.character(raw_vertex)[1]
  if (length(raw_vertex) == 0 || is.na(raw_vertex) || raw_vertex == "") return(NA_character_)
  if (raw_vertex %in% graph_nodes) return(raw_vertex)
  raw_num <- suppressWarnings(as.numeric(raw_vertex))
  if (!is.na(raw_num)) {
    raw_idx <- as.integer(raw_num)
    if (raw_idx >= 1 && raw_idx <= length(graph_nodes)) return(graph_nodes[raw_idx])
  }
  NA_character_
}

extract_graph_structure <- function(cds) {
  graph_obj <- principal_graph(cds)[["UMAP"]]
  graph_coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(graph_obj) || is.null(graph_coords)) stop("Missing Monocle3 principal graph.")
  edge_df <- igraph::as_data_frame(graph_obj, what = "edges") %>%
    mutate(
      x = graph_coords[1, from],
      y = graph_coords[2, from],
      xend = graph_coords[1, to],
      yend = graph_coords[2, to]
    )
  node_df <- data.frame(node = colnames(graph_coords), x = graph_coords[1, ], y = graph_coords[2, ])
  list(
    graph = graph_obj,
    graph_coords = graph_coords,
    closest_vertex = flatten_closest_vertex(closest_vertex),
    edges = edge_df,
    nodes = node_df
  )
}

compute_graph_weights <- function(graph_obj, graph_coords) {
  edge_df <- igraph::as_data_frame(graph_obj, what = "edges")
  if (nrow(edge_df) == 0) return(graph_obj)
  igraph::E(graph_obj)$weight <- purrr::map2_dbl(edge_df$from, edge_df$to, function(from_node, to_node) {
    sqrt(sum((graph_coords[, from_node] - graph_coords[, to_node]) ^ 2))
  })
  graph_obj
}

nearest_graph_vertex <- function(point_xy, graph_coords) {
  vertex_names <- colnames(graph_coords)
  dists <- apply(graph_coords, 2, function(node_xy) sqrt(sum((point_xy - node_xy) ^ 2)))
  vertex_names[which.min(dists)]
}

get_state_medoid_cell <- function(umap_df) {
  if (nrow(umap_df) == 1) return(umap_df$cell[1])
  dist_mat <- as.matrix(dist(umap_df[, c("UMAP_1", "UMAP_2"), drop = FALSE]))
  umap_df$cell[which.min(rowSums(dist_mat))]
}

get_root_pr_node <- function(cds, root_cells) {
  graph_bits <- extract_graph_structure(cds)
  root_cells <- intersect(root_cells, names(graph_bits$closest_vertex))
  if (length(root_cells) == 0) return(NULL)
  root_vertex_raw <- names(sort(table(as.character(graph_bits$closest_vertex[root_cells])), decreasing = TRUE))[1]
  root_vertex <- coerce_graph_vertex_name(root_vertex_raw, graph_bits$graph)
  if (!is.na(root_vertex)) return(root_vertex)
  NULL
}

get_root_label_position <- function(cds, graph_bits, root_cells) {
  root_node <- NA_character_
  common_root <- intersect(root_cells, names(graph_bits$closest_vertex))
  if (length(common_root) > 0) {
    vals <- as.character(graph_bits$closest_vertex[common_root])
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      root_node <- coerce_graph_vertex_name(names(sort(table(vals), decreasing = TRUE))[1], graph_bits$graph)
    }
  }
  if (!is.na(root_node) && root_node %in% graph_bits$nodes$node) {
    return(graph_bits$nodes %>% filter(node == root_node) %>% slice(1) %>% mutate(label = "ROOT"))
  }
  umap_mat <- reducedDims(cds)$UMAP
  common_root <- intersect(root_cells, rownames(umap_mat))
  if (length(common_root) > 0) {
    xy <- umap_mat[common_root, , drop = FALSE]
    return(data.frame(node = "root_centroid", x = median(xy[, 1]), y = median(xy[, 2]), label = "ROOT"))
  }
  NULL
}

project_point_to_segment <- function(px, py, ax, ay, bx, by) {
  abx <- bx - ax
  aby <- by - ay
  ab2 <- abx * abx + aby * aby
  if (!is.finite(ab2) || ab2 == 0) {
    return(list(x = ax, y = ay, t = 0, dist2 = (px - ax)^2 + (py - ay)^2))
  }
  t <- max(0, min(1, ((px - ax) * abx + (py - ay) * aby) / ab2))
  proj_x <- ax + t * abx
  proj_y <- ay + t * aby
  list(x = proj_x, y = proj_y, t = t, dist2 = (px - proj_x)^2 + (py - proj_y)^2)
}

project_cells_to_graph <- function(cds, edges_df, state_map) {
  umap_mat <- reducedDims(cds)$UMAP
  cell_names <- rownames(umap_mat)
  pt_vec <- pseudotime(cds)
  bind_rows(lapply(seq_along(cell_names), function(i) {
    px <- umap_mat[i, 1]
    py <- umap_mat[i, 2]
    best_dist2 <- Inf
    best_p <- NULL
    best_e <- 1L
    for (e in seq_len(nrow(edges_df))) {
      pr <- project_point_to_segment(px, py, edges_df$x[e], edges_df$y[e], edges_df$xend[e], edges_df$yend[e])
      if (pr$dist2 < best_dist2) {
        best_dist2 <- pr$dist2
        best_p <- pr
        best_e <- e
      }
    }
    data.frame(
      cell = cell_names[i],
      umap_x = px,
      umap_y = py,
      graph_x = best_p$x,
      graph_y = best_p$y,
      projection_distance = sqrt(best_dist2),
      edge_index = best_e,
      state = factor(state_map[cell_names[i]], levels = primary_states),
      pseudotime = as.numeric(pt_vec[cell_names[i]])
    )
  }))
}

build_weighted_ridge_df <- function(proj_df, state_order, n_points = 512) {
  df <- proj_df %>%
    filter(is.finite(pseudotime), !is.na(state)) %>%
    mutate(state = factor(as.character(state), levels = state_order))
  if (nrow(df) == 0) return(NULL)
  counts <- df %>%
    count(state, name = "n") %>%
    complete(state = factor(state_order, levels = state_order), fill = list(n = 0)) %>%
    mutate(weight = ifelse(max(n) > 0, n / max(n), 0), y = rev(seq_along(state_order)))
  pt_range <- range(df$pseudotime, na.rm = TRUE)
  if (!all(is.finite(pt_range))) return(NULL)
  if (diff(pt_range) == 0) pt_range <- pt_range + c(-1e-6, 1e-6)
  focus_range <- as.numeric(quantile(df$pseudotime, c(0.01, 0.99), na.rm = TRUE, names = FALSE))
  if (!all(is.finite(focus_range)) || diff(focus_range) <= 0) focus_range <- pt_range
  ridge_df <- bind_rows(lapply(seq_len(nrow(counts)), function(i) {
    state_i <- counts$state[i]
    n_i <- counts$n[i]
    sub <- df %>% filter(state == state_i)
    if (n_i < 2) {
      x_grid <- seq(pt_range[1], pt_range[2], length.out = n_points)
      dens_y <- rep(0, length(x_grid))
    } else {
      dens <- density(sub$pseudotime, from = pt_range[1], to = pt_range[2], n = n_points, na.rm = TRUE, bw = "nrd0")
      x_grid <- dens$x
      dens_y <- dens$y
    }
    data.frame(state = factor(as.character(state_i), levels = state_order), x = x_grid, height = dens_y * counts$weight[i], y = counts$y[i], n_cells = n_i)
  }))
  max_h <- max(ridge_df$height, na.rm = TRUE)
  if (is.finite(max_h) && max_h > 0) ridge_df <- ridge_df %>% mutate(height = height / max_h * 0.9)
  list(
    ridge_df = ridge_df,
    state_labels = counts %>% mutate(label = paste0(as.character(state), " (n=", n, ")")) %>% select(state, y, label),
    counts = counts,
    focus_range = focus_range
  )
}

build_pdo_inputs <- function() {
  message("Loading PDO Seurat object and pre-relabel state vector ...")
  pdos <- load_pdo_object()
  state_vec <- load_pre_relabel_states()
  common_cells <- intersect(Cells(pdos), names(state_vec))
  if (length(common_cells) == 0) stop("No overlapping cells between PDO object and state vector.")
  pdos <- pdos[, common_cells]
  state_vec <- state_vec[common_cells]
  pdos$state_pre_relabel <- state_vec[Cells(pdos)]
  meta_df <- data.frame(
    cell = Cells(pdos),
    sample = as.character(pdos$orig.ident),
    state = as.character(pdos$state_pre_relabel)
  ) %>% filter(state %in% primary_states)
  list(pdos = pdos, state_vec = state_vec, meta_df = meta_df)
}

summarise_pdo_sample_inclusion <- function(meta_df) {
  sample_state_counts <- meta_df %>%
    count(sample, state, name = "n_cells") %>%
    complete(sample, state = primary_states, fill = list(n_cells = 0))
  sample_summary <- sample_state_counts %>%
    group_by(sample) %>%
    summarise(
      total_primary_cells = sum(n_cells),
      n_states_over_threshold = sum(n_cells >= OTHER_MIN_CELLS),
      root_state_cells = n_cells[state == root_state][1],
      qualifies_total = total_primary_cells >= TOTAL_MIN_CELLS,
      qualifies_multistate = n_states_over_threshold >= MIN_STATES_OVER_THRESHOLD,
      qualifies_root = root_state_cells >= ROOT_MIN_CELLS,
      valid_sample = qualifies_total & qualifies_multistate & qualifies_root,
      .groups = "drop"
    ) %>% arrange(sample)
  root_summary <- sample_state_counts %>%
    left_join(sample_summary %>% select(sample, valid_sample), by = "sample") %>%
    group_by(sample) %>%
    mutate(other_state_max = vapply(state, function(this_state) max(n_cells[state != this_state], na.rm = TRUE), numeric(1))) %>%
    ungroup() %>%
    mutate(valid_root = valid_sample & n_cells >= ROOT_MIN_CELLS & other_state_max >= OTHER_MIN_CELLS) %>%
    rename(root_state = state, root_state_cells = n_cells)
  list(sample_summary = sample_summary, root_summary = root_summary)
}

build_pdo_sample_trajectories <- function(rebuild = FALSE) {
  inputs <- build_pdo_inputs()
  pdos <- inputs$pdos
  state_vec <- inputs$state_vec
  meta_df <- inputs$meta_df
  inclusion <- summarise_pdo_sample_inclusion(meta_df)
  write_csv(inclusion$sample_summary, file.path(pseudotime_dir, "Auto_PDO_pseudotime_sample_inclusion_summary.csv"))
  write_csv(inclusion$root_summary, file.path(pseudotime_dir, "Auto_PDO_pseudotime_root_state_summary.csv"))
  valid_samples <- inclusion$sample_summary %>% filter(valid_sample) %>% pull(sample)
  if (length(valid_samples) == 0) stop("No PDO samples passed pseudotime inclusion thresholds.")
  message("Valid samples detected: ", length(valid_samples))
  results <- list()
  for (sample_id in valid_samples) {
    sample_tag <- sanitize_sample_id(sample_id)
    cds_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_cds.rds"))
    pt_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_pseudotime.rds"))
    meta_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_metadata.csv"))
    proj_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_projections.csv"))
    if (!rebuild && all(file.exists(cds_path, pt_path, meta_path, proj_path))) {
      message("Loading cached trajectory for sample: ", sample_id)
      results[[sample_id]] <- list(sample = sample_id, cds = readRDS(cds_path), pseudotime = readRDS(pt_path), metadata = read_csv(meta_path, show_col_types = FALSE), projections = read_csv(proj_path, show_col_types = FALSE))
      next
    }
    message("Building trajectory for sample: ", sample_id)
    sample_cells <- meta_df %>% filter(sample == sample_id) %>% pull(cell)
    sample_obj <- pdos[, sample_cells]
    sample_obj$state_pre_relabel <- factor(state_vec[Cells(sample_obj)], levels = primary_states)
    root_cells <- colnames(sample_obj)[as.character(sample_obj$state_pre_relabel) == root_state]
    sample_result <- tryCatch({
      cds <- prepare_pdo_sample_trajectory(sample_obj)
      cds <- order_cells(cds, root_cells = intersect(root_cells, colnames(cds)))
      pt <- pseudotime(cds)
      pt[is.infinite(pt)] <- NA_real_
      graph_bits <- extract_graph_structure(cds)
      proj_df <- project_cells_to_graph(cds, graph_bits$edges, state_vec)
      meta_out <- data.frame(cell = names(pt), sample = sample_id, state = as.character(state_vec[names(pt)]), pseudotime = as.numeric(pt))
      saveRDS(cds, cds_path)
      saveRDS(pt, pt_path)
      write_csv(meta_out, meta_path)
      write_csv(proj_df, proj_path)
      list(sample = sample_id, cds = cds, pseudotime = pt, metadata = meta_out, projections = proj_df)
    }, error = function(e) {
      message("Skipping sample after pseudotime failure: ", sample_id, " | ", e$message)
      NULL
    })
    if (!is.null(sample_result)) results[[sample_id]] <- sample_result
  }
  results
}

get_or_build_pdo_sample_trajectories <- function(rebuild = FALSE) {
  if (rebuild) return(build_pdo_sample_trajectories(rebuild = TRUE))
  inclusion_path <- file.path(pseudotime_dir, "Auto_PDO_pseudotime_sample_inclusion_summary.csv")
  if (!file.exists(inclusion_path)) return(build_pdo_sample_trajectories(rebuild = FALSE))
  valid_samples <- read_csv(inclusion_path, show_col_types = FALSE) %>% filter(valid_sample) %>% pull(sample)
  cached <- list()
  missing_cache <- character()
  for (sample_id in valid_samples) {
    sample_tag <- sanitize_sample_id(sample_id)
    cds_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_cds.rds"))
    pt_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_pseudotime.rds"))
    meta_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_metadata.csv"))
    proj_path <- file.path(asset_dir, paste0("Auto_PDO_", sample_tag, "_projections.csv"))
    if (!all(file.exists(cds_path, pt_path, meta_path, proj_path))) {
      missing_cache <- c(missing_cache, sample_id)
      next
    }
    cached[[sample_id]] <- list(sample = sample_id, cds = readRDS(cds_path), pseudotime = readRDS(pt_path), metadata = read_csv(meta_path, show_col_types = FALSE), projections = read_csv(proj_path, show_col_types = FALSE))
  }
  if (length(missing_cache) > 0 || length(cached) == 0) {
    message("Missing cached trajectories for: ", paste(missing_cache, collapse = ", "))
    return(build_pdo_sample_trajectories(rebuild = FALSE))
  }
  cached
}
