####################
# Analysis registry:
#   Status: active shared helper library
#   Script: analysis/shared/Auto_pdo_analysis_helpers.R
#   Methodology: analysis/methodology/shared/shared_config_and_logging_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs: analysis/shared/Auto_pdo_analysis_config.R
#   Outputs: reusable helper functions for config, cache policy, file checks,
#            state-vector normalization, output tiers, plotting, and run logs
####################

####################
# Shared PDO analysis helpers
####################

if (!exists("PDO_PROJECT_DIR")) {
  source("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")
}

pdo_get_env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

pdo_get_env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) default else parsed
}

pdo_get_env_integer <- function(name, default) {
  as.integer(round(pdo_get_env_numeric(name, default)))
}

pdo_cache_policy <- function() {
  list(
    force_rebuild = pdo_get_env_flag(PDO_CACHE_ENV$force_rebuild, FALSE),
    replot_only = pdo_get_env_flag(PDO_CACHE_ENV$replot_only, FALSE)
  )
}

pdo_output_path <- function(..., base_dir = PDO_OUTPUT_DIR) {
  file.path(base_dir, ...)
}

pdo_ensure_output_tiers <- function(out_dir, tiers = PDO_OUTPUT_TIERS) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tier_paths <- file.path(out_dir, tiers)
  invisible(vapply(tier_paths, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
  setNames(tier_paths, tiers)
}

pdo_require_files <- function(paths, labels = basename(paths)) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    msg <- paste0(labels[match(missing, paths)], ": ", missing, collapse = "\n")
    stop("Required input file(s) not found:\n", msg)
  }
  invisible(paths)
}

pdo_normalise_final_state_vector <- function(state_vec) {
  state_names <- names(state_vec)
  state_vec <- as.character(state_vec)
  names(state_vec) <- state_names
  state_vec
}

pdo_filter_state_levels <- function(states, include_unresolved = TRUE, include_hybrid = TRUE) {
  levels <- PDO_STATE_ORDER
  if (include_unresolved) {
    levels <- c(levels, "Unresolved")
  }
  if (include_hybrid) {
    levels <- c(levels, "Hybrid")
  }
  factor(states, levels = levels)
}

pdo_get_assay_matrix <- function(seurat_obj, assay = "RNA", layer = c("data", "counts")) {
  layer <- match.arg(layer)
  mat <- tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer),
      error = function(e) NULL
    )
  }
  if (is.null(mat)) {
    mat <- tryCatch(
      SeuratObject::LayerData(seurat_obj, assay = assay, layer = layer),
      error = function(e) NULL
    )
  }
  if (is.null(mat)) {
    stop("Could not read assay matrix: assay=", assay, ", layer=", layer)
  }
  mat[, Seurat::Cells(seurat_obj), drop = FALSE]
}

pdo_theme_slide <- function(base_size = PDO_PLOT_DEFAULTS$base_size) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = PDO_PLOT_DEFAULTS$axis_text_size),
      axis.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = PDO_PLOT_DEFAULTS$legend_text_size),
      legend.title = ggplot2::element_text(size = PDO_PLOT_DEFAULTS$legend_title_size),
      strip.text = ggplot2::element_text(size = PDO_PLOT_DEFAULTS$strip_text_size),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold", hjust = 0.5)
    )
}

pdo_save_slide_pdf <- function(plot, filename, width = PDO_PLOT_DEFAULTS$slide_pdf_width, height = PDO_PLOT_DEFAULTS$slide_pdf_height) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(filename)
}

pdo_write_run_summary <- function(script,
                                  out_dir,
                                  inputs = character(),
                                  outputs = character(),
                                  parameters = list(),
                                  cache = pdo_cache_policy(),
                                  status = "completed",
                                  include_session = TRUE) {
  tier_paths <- pdo_ensure_output_tiers(out_dir)
  log_file <- file.path(tier_paths[["logs"]], paste0(tools::file_path_sans_ext(basename(script)), "_run_summary.txt"))
  lines <- c(
    paste0("script: ", script),
    paste0("status: ", status),
    paste0("start_or_write_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "inputs:",
    if (length(inputs) == 0) "  none recorded" else paste0("  - ", inputs),
    "",
    "outputs:",
    if (length(outputs) == 0) "  none recorded" else paste0("  - ", outputs),
    "",
    "parameters:",
    if (length(parameters) == 0) "  none recorded" else paste0("  - ", names(parameters), ": ", unlist(parameters, use.names = FALSE)),
    "",
    "cache:",
    paste0("  - ", names(cache), ": ", unlist(cache, use.names = FALSE))
  )
  if (include_session) {
    lines <- c(lines, "", "sessionInfo:", paste0("  ", capture.output(utils::sessionInfo())))
  }
  writeLines(lines, log_file)
  invisible(log_file)
}
