####################
# Auto_PDO_numbat_phylogeny_visualisation.R
#
# Analysis registry
# Status: active terminal figure workflow.
# Script: analysis/cnv/Auto_PDO_numbat_phylogeny_visualisation.R
# Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
# - PDOs_outs/Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv
# - PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/tree_final_<iter>.rds
# - PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/clones_<iter>.rds
# Outputs:
# - PDOs_outs/Auto_PDO_numbat/phylogeny/Auto_PDO_numbat_phylogenetic_trees.pdf
# - PDOs_outs/Auto_PDO_numbat/phylogeny/Auto_PDO_numbat_phylogenetic_tree_summary.csv
# Downstream use: terminal visualization and audit table; not a required input
# for clone calling, concordance, or MP/state tests.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(RColorBrewer)
  library(scales)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

args <- commandArgs(trailingOnly = TRUE)
sample_arg <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "all"

manifest_path <- "Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv"
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path)

manifest <- fread(manifest_path)
manifest <- manifest[sample != "SUR843T3_PDO"]
if (!identical(sample_arg, "all")) {
  requested <- trimws(unlist(strsplit(sample_arg, ",")))
  manifest <- manifest[sample %in% requested]
}
if (nrow(manifest) == 0) stop("No samples found for argument: ", sample_arg)

out_dir <- "Auto_PDO_numbat/phylogeny"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_palette <- function(values) {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(character(0))
  base <- brewer.pal(max(3, min(12, length(values))), if (length(values) <= 8) "Dark2" else "Paired")
  setNames(colorRampPalette(base)(length(values)), values)
}

final_iter_from <- function(numbat_dir, prefix = "tree_final") {
  files <- Sys.glob(file.path(numbat_dir, paste0(prefix, "_*.rds")))
  if (length(files) == 0) return(NA_integer_)
  iter <- suppressWarnings(as.integer(sub(paste0("^", prefix, "_([0-9]+)\\.rds$"), "\\1", basename(files))))
  iter <- iter[is.finite(iter)]
  if (length(iter) == 0) NA_integer_ else max(iter)
}

as_vertex_df <- function(g) {
  attrs <- vertex_attr(g)
  as.data.frame(lapply(attrs, function(x) {
    if (is.list(x)) {
      vapply(x, function(y) paste(y, collapse = ","), character(1))
    } else {
      x
    }
  }), stringsAsFactors = FALSE)
}

load_clone_members <- function(numbat_dir, iter, clone_counts) {
  clone_file <- file.path(numbat_dir, paste0("clones_", iter, ".rds"))
  if (!file.exists(clone_file)) {
    return(data.frame(
      clone = names(clone_counts),
      n_cells = as.integer(clone_counts),
      frac = as.numeric(clone_counts) / sum(clone_counts),
      events = "",
      stringsAsFactors = FALSE
    ))
  }
  clones <- readRDS(clone_file)
  out <- bind_rows(lapply(names(clones), function(cl) {
    x <- clones[[cl]]
    data.frame(
      clone = as.character(cl),
      n_cells = if (!is.null(x$size)) as.integer(x$size) else length(x$cells),
      events = if (!is.null(x$members)) paste(x$members, collapse = ",") else "",
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(out) == 0) {
    out <- data.frame(clone = names(clone_counts), n_cells = as.integer(clone_counts), events = "", stringsAsFactors = FALSE)
  }
  out <- out %>%
    mutate(n_cells = ifelse(is.na(.data$n_cells) | .data$n_cells <= 0,
                            as.integer(clone_counts[.data$clone]),
                            .data$n_cells),
           frac = .data$n_cells / sum(.data$n_cells, na.rm = TRUE)) %>%
    arrange(desc(.data$n_cells), .data$clone)
  out
}

plot_full_tree <- function(g, vertex_df, clone_cols, sample_id, iter) {
  root_idx <- which(as.logical(vertex_df$root))
  if (length(root_idx) == 0) root_idx <- 1L
  root_idx <- root_idx[1]
  leaf <- as.logical(vertex_df$leaf)
  clone <- as.character(vertex_df$clone)
  clone[is.na(clone) | !nzchar(clone)] <- "unknown"
  v_cols <- ifelse(leaf, clone_cols[clone], adjustcolor("grey55", alpha.f = 0.35))
  v_cols[is.na(v_cols)] <- "grey55"
  leaf_n <- sum(leaf, na.rm = TRUE)
  v_size <- ifelse(leaf, ifelse(leaf_n > 6000, 0.35, ifelse(leaf_n > 2500, 0.55, 0.85)), 0.05)
  edge_cols <- adjustcolor("grey45", alpha.f = ifelse(leaf_n > 6000, 0.22, 0.35))
  edge_w <- ifelse(leaf_n > 6000, 0.18, 0.25)
  layout <- layout_as_tree(g, root = root_idx, circular = TRUE, mode = "out")
  plot(
    g,
    layout = layout,
    vertex.label = NA,
    vertex.size = v_size,
    vertex.color = v_cols,
    vertex.frame.color = NA,
    edge.arrow.size = 0,
    edge.color = edge_cols,
    edge.width = edge_w,
    margin = -0.06,
    asp = 1
  )
  title(main = paste0(sample_id, " Numbat phylogeny"), sub = paste0("tree_final_", iter, " | tips=", comma(leaf_n)), line = 0.2, cex.main = 1.1, cex.sub = 0.75)
}

plot_clone_bar <- function(clone_df, clone_cols) {
  clone_df <- clone_df %>% arrange(.data$n_cells)
  cols <- clone_cols[clone_df$clone]
  cols[is.na(cols)] <- "grey70"
  bp <- barplot(
    clone_df$frac,
    names.arg = paste0("C", clone_df$clone),
    horiz = TRUE,
    las = 1,
    col = cols,
    border = NA,
    xlim = c(0, max(clone_df$frac, na.rm = TRUE) * 1.18),
    xlab = "Fraction of cells",
    main = "Clone sizes",
    cex.names = 0.75,
    cex.axis = 0.75,
    cex.lab = 0.85,
    cex.main = 0.95
  )
  text(clone_df$frac, bp, labels = paste0(comma(clone_df$n_cells), " (", percent(clone_df$frac, accuracy = 0.1), ")"), pos = 4, cex = 0.62, xpd = TRUE)
}

clone_transition_graph <- function(g, vertex_df) {
  edges <- as_data_frame(g, what = "edges")
  if (nrow(edges) == 0) return(make_empty_graph())
  clone <- as.character(vertex_df$clone)
  clone[is.na(clone) | !nzchar(clone)] <- "unknown"
  edge_df <- data.frame(
    from_clone = clone[edges$from],
    to_clone = clone[edges$to],
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(.data$from_clone), !is.na(.data$to_clone), .data$from_clone != .data$to_clone) %>%
    count(.data$from_clone, .data$to_clone, name = "n")
  if (nrow(edge_df) == 0) return(make_empty_graph())
  graph_from_data_frame(edge_df, directed = TRUE)
}

plot_clone_tree <- function(g, vertex_df, clone_df, clone_cols) {
  cg <- clone_transition_graph(g, vertex_df)
  if (vcount(cg) == 0) {
    plot.new()
    title("Clone-level transitions", cex.main = 0.95)
    text(0.5, 0.55, "No between-clone\ntransition edges", cex = 0.9)
    return(invisible(NULL))
  }
  counts <- setNames(clone_df$n_cells, clone_df$clone)
  clone_names <- V(cg)$name
  v_cols <- clone_cols[clone_names]
  v_cols[is.na(v_cols)] <- "grey75"
  v_size <- scales::rescale(sqrt(pmax(counts[clone_names], 1)), to = c(10, 28), from = range(sqrt(pmax(counts, 1)), na.rm = TRUE))
  v_size[!is.finite(v_size)] <- 16
  edge_w <- scales::rescale(E(cg)$n, to = c(0.8, 4), from = range(E(cg)$n, na.rm = TRUE))
  edge_w[!is.finite(edge_w)] <- 1.5
  lay <- tryCatch(layout_as_tree(cg, circular = FALSE), error = function(e) layout_nicely(cg))
  plot(
    cg,
    layout = lay,
    vertex.label = paste0("C", clone_names),
    vertex.label.cex = 0.8,
    vertex.label.color = "black",
    vertex.color = v_cols,
    vertex.frame.color = "white",
    vertex.size = v_size,
    edge.arrow.size = 0.35,
    edge.color = adjustcolor("grey35", alpha.f = 0.75),
    edge.width = edge_w,
    margin = 0.1
  )
  title("Clone-level transitions", cex.main = 0.95)
}

plot_event_text <- function(clone_df) {
  plot.new()
  title("Clone CNV event sets", cex.main = 0.95)
  y <- 0.96
  text(0, y, "Clone", adj = c(0, 1), font = 2, cex = 0.72)
  text(0.18, y, "Cells", adj = c(0, 1), font = 2, cex = 0.72)
  text(0.38, y, "Numbat events", adj = c(0, 1), font = 2, cex = 0.72)
  y <- y - 0.07
  max_rows <- min(nrow(clone_df), 12)
  for (i in seq_len(max_rows)) {
    events <- clone_df$events[i]
    if (is.na(events) || !nzchar(events)) events <- "(root/no additional event set)"
    if (nchar(events) > 70) events <- paste0(substr(events, 1, 67), "...")
    text(0, y, paste0("C", clone_df$clone[i]), adj = c(0, 1), cex = 0.68)
    text(0.18, y, paste0(comma(clone_df$n_cells[i]), " / ", percent(clone_df$frac[i], accuracy = 0.1)), adj = c(0, 1), cex = 0.68)
    text(0.38, y, events, adj = c(0, 1), cex = 0.62)
    y <- y - 0.065
  }
  if (nrow(clone_df) > max_rows) {
    text(0, y, paste0("... ", nrow(clone_df) - max_rows, " additional clones omitted"), adj = c(0, 1), cex = 0.65, col = "grey40")
  }
}

summary_rows <- list()
pdf_path <- file.path(out_dir, "Auto_PDO_numbat_phylogenetic_trees.pdf")
pdf(pdf_path, width = 16, height = 10, useDingbats = FALSE)

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample[i]
  numbat_dir <- manifest$numbat_dir[i]
  iter <- final_iter_from(numbat_dir, "tree_final")
  if (!is.finite(iter)) {
    plot.new()
    title(sample_id)
    text(0.5, 0.5, "No tree_final_*.rds found")
    summary_rows[[sample_id]] <- data.frame(sample = sample_id, status = "missing_tree", stringsAsFactors = FALSE)
    next
  }
  tree_file <- file.path(numbat_dir, paste0("tree_final_", iter, ".rds"))
  message("Plotting ", sample_id, " from ", basename(tree_file))
  g <- readRDS(tree_file)
  vertex_df <- as_vertex_df(g)
  leaf <- as.logical(vertex_df$leaf)
  clone <- as.character(vertex_df$clone)
  clone[is.na(clone) | !nzchar(clone)] <- "unknown"
  clone_counts <- sort(table(clone[leaf]), decreasing = TRUE)
  clone_cols <- make_palette(names(clone_counts))
  clone_df <- load_clone_members(numbat_dir, iter, clone_counts)
  clone_df$clone <- as.character(clone_df$clone)
  clone_df <- clone_df[clone_df$clone %in% names(clone_counts), , drop = FALSE]
  if (nrow(clone_df) == 0) {
    clone_df <- data.frame(clone = names(clone_counts), n_cells = as.integer(clone_counts), frac = as.numeric(clone_counts) / sum(clone_counts), events = "", stringsAsFactors = FALSE)
  }

  layout(matrix(c(1, 2, 1, 3, 1, 4), nrow = 3, byrow = TRUE), widths = c(2.25, 1), heights = c(1, 1, 1))
  par(mar = c(0.3, 0.3, 2.0, 0.3))
  plot_full_tree(g, vertex_df, clone_cols, sample_id, iter)
  par(mar = c(4.0, 6.5, 2.0, 1.0))
  plot_clone_bar(clone_df, clone_cols)
  par(mar = c(0.8, 0.8, 2.0, 0.8))
  plot_clone_tree(g, vertex_df, clone_df, clone_cols)
  par(mar = c(0.5, 0.5, 2.0, 0.5))
  plot_event_text(clone_df)

  summary_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    status = "plotted",
    iter = iter,
    tree_file = tree_file,
    n_vertices = vcount(g),
    n_edges = ecount(g),
    n_tips = sum(leaf, na.rm = TRUE),
    n_clones = length(clone_counts),
    largest_clone = names(clone_counts)[1],
    largest_clone_cells = as.integer(clone_counts[1]),
    stringsAsFactors = FALSE
  )
}

dev.off()

summary_df <- bind_rows(summary_rows)
fwrite(summary_df, file.path(out_dir, "Auto_PDO_numbat_phylogenetic_tree_summary.csv"))

message("Wrote: ", file.path(out_root, pdf_path))
message("Wrote: ", file.path(out_root, out_dir, "Auto_PDO_numbat_phylogenetic_tree_summary.csv"))
