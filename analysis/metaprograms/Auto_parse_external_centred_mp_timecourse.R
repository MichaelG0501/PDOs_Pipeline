####################
# Analysis registry:
#   Status: active terminal cross-dataset comparison
#   Script: analysis/metaprograms/Auto_parse_external_centred_mp_timecourse.R
#   Methodology: analysis/methodology/metaprograms/Auto_parse_external_centred_mp_timecourse_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Method:
#     Scores the current retained centred-refined PDO and scRef MP gene sets in
#     the six Parse timepoints with UCell. MP numbers are source-specific and
#     are never matched across PDO and scRef. Rows are ordered by the biological
#     groupings defined by each current ordered-heatmap workflow, with source
#     shown explicitly in every display label and heatmap annotation.
#   Exact inputs:
#     - PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     - PDOs_outs/centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv
#     - /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#     - /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv
#     - /rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/parse_outs/by_samples/{T0,T1,T2,T4,R4,eR4}/Auto_{sample}_final.rds
#   Exact output directory:
#     - PDOs_outs/Auto_parse_external_centred_mp_timecourse/
#   Key persistent downstream-capable outputs:
#     - intermediate/Auto_parse_external_centred_mp_signatures.rds
#     - intermediate/Auto_parse_external_centred_mp_ucell_scores.rds
#     - intermediate/Auto_parse_external_centred_mp_cell_metadata.rds
#     - tables/Auto_parse_external_centred_mp_signature_manifest.csv
#     - tables/Auto_parse_external_centred_mp_timepoint_summary.csv
#     - tables/Auto_parse_external_centred_mp_timepoint_tests.csv
#   Terminal figure outputs:
#     - figures/Auto_parse_external_centred_mp_timecourse_report.pdf
#     - figures/Auto_parse_external_centred_mp_timecourse_raw_mean_heatmap.pdf
#     - figures/Auto_parse_external_centred_mp_timecourse_zscore_heatmap.pdf
#     - figures/Auto_parse_external_centred_mp_timecourse_zscore_heatmap.png
#     - figures/Auto_parse_external_centred_mp_timecourse_abundance.pdf
#   Cache/replot:
#     - PDO_FORCE_REBUILD=1 recomputes every per-sample UCell cache.
#     - PDO_REPLOT_ONLY=1 requires the combined live caches and only rebuilds
#       tables and figures.
#   Conda env: dmtcp
#   Downstream status: score/signature tables are reusable; figures are terminal.
####################

####################
# Libraries, shared configuration, and fixed paths
####################
suppressPackageStartupMessages({
  library("SeuratObject")
  library("UCell")
  library("Matrix")
  library("dplyr")
  library("tidyr")
  library("tibble")
  library("ggplot2")
  library("ComplexHeatmap")
  library("circlize")
  library("grid")
})

script_start <- Sys.time()
script_path <- "analysis/metaprograms/Auto_parse_external_centred_mp_timecourse.R"
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(project_dir)

source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_helpers.R"))

out_dir <- file.path(PDO_LIVE_OUTS, "Auto_parse_external_centred_mp_timecourse")
out_tiers <- pdo_ensure_output_tiers(out_dir)
cache_policy <- pdo_cache_policy()

if (cache_policy$force_rebuild && cache_policy$replot_only) {
  stop("PDO_FORCE_REBUILD=1 and PDO_REPLOT_ONLY=1 cannot be used together.")
}

parse_samples <- c("T0", "T1", "T2", "T4", "R4", "eR4")
parse_colours <- c(
  "T0" = "#0072B2",
  "T1" = "#E69F00",
  "T2" = "#009E73",
  "T4" = "#D55E00",
  "R4" = "#CC79A7",
  "eR4" = "#56B4E9"
)
source_colours <- c("PDO" = "#D55E00", "scRef" = "#0072B2")
group_levels <- c(
  "Cell cycle",
  "Classic proliferation",
  "Columnar-to-intestinal",
  "Squamous-to-intestinal",
  "Glandular differentiation",
  "Glandular-to-intestinal",
  "Stress-adaptive",
  "ECM-remodelling",
  "Cancer-cell immune mimicry",
  "Motile-cilia differentiation"
)
group_colours <- c(
  "Cell cycle" = "#6B7280",
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Cancer-cell immune mimicry" = "#377EB8",
  "Motile-cilia differentiation" = "#F781BF"
)

pdo_gene_file <- file.path(
  PDO_LIVE_OUTS,
  "centred_mp_refinement",
  "merged_refined_mp_genes.rds"
)
pdo_group_file <- file.path(
  PDO_LIVE_OUTS,
  "centred_mp_refinement",
  "tables",
  "centred_refined_mp_state_grouping.csv"
)
scref_base <- paste0(
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/",
  "ref_outs/Metaprogrammes_Results/centred/mp_refinement"
)
scref_gene_file <- file.path(
  scref_base,
  "intermediate",
  "merged_refined_mp_genes.rds"
)
scref_group_file <- file.path(
  scref_base,
  "tables",
  "centred_refined_mp_state_grouping.csv"
)
parse_base <- paste0(
  "/rds/general/project/spatialtranscriptomics/live/",
  "Parse_Pipeline/parse_outs/by_samples"
)
parse_files <- setNames(
  file.path(
    parse_base,
    parse_samples,
    paste0("Auto_", parse_samples, "_final.rds")
  ),
  parse_samples
)

input_files <- c(
  pdo_gene_file,
  pdo_group_file,
  scref_gene_file,
  scref_group_file,
  unname(parse_files)
)
pdo_require_files(input_files)
####################

####################
# Build the source-safe, biologically grouped signature manifest
####################
pdo_genes <- readRDS(pdo_gene_file)
scref_genes <- readRDS(scref_gene_file)
pdo_groups <- read.csv(
  pdo_group_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
scref_groups <- read.csv(
  scref_group_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_group_cols <- c("state", "mp", "description")
for (tbl_name in c("pdo_groups", "scref_groups")) {
  tbl <- get(tbl_name)
  missing_cols <- setdiff(required_group_cols, colnames(tbl))
  if (length(missing_cols) > 0) {
    stop(
      tbl_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

# The current scRef ordered-heatmap definition explicitly labels these as
# Excluded. They remain documented in the input table but are not active MP
# signatures in this cross-dataset score comparison.
scref_groups <- scref_groups |>
  filter(state != "Excluded")

make_manifest <- function(group_table, gene_lists, source_name) {
  missing_gene_lists <- setdiff(group_table$mp, names(gene_lists))
  if (length(missing_gene_lists) > 0) {
    stop(
      source_name,
      " current grouping contains MP(s) absent from its gene-list object: ",
      paste(missing_gene_lists, collapse = ", ")
    )
  }

  group_table |>
    transmute(
      source = source_name,
      biological_group = as.character(state),
      mp = as.character(mp),
      description = as.character(description),
      signature_key = paste0(source_name, "__", mp),
      n_genes_reference = lengths(gene_lists[mp])
    )
}

signature_manifest <- bind_rows(
  make_manifest(pdo_groups, pdo_genes, "PDO"),
  make_manifest(scref_groups, scref_genes, "scRef")
) |>
  mutate(
    biological_group = factor(biological_group, levels = group_levels),
    source = factor(source, levels = c("PDO", "scRef")),
    description = ifelse(
      is.na(description) | description == "" | description == mp,
      "",
      description
    )
  ) |>
  arrange(biological_group, source) |>
  mutate(
    row_order = row_number(),
    display_label = ifelse(
      description == "",
      paste0("[", source, "] ", mp),
      paste0("[", source, "] ", description)
    )
  )

if (any(is.na(signature_manifest$biological_group))) {
  stop(
    "Unrecognised biological grouping(s): ",
    paste(
      unique(as.character(signature_manifest$biological_group[
        is.na(signature_manifest$biological_group)
      ])),
      collapse = ", "
    )
  )
}
if (anyDuplicated(signature_manifest$signature_key) > 0) {
  stop("Source-prefixed signature keys are not unique.")
}

signature_genes <- c(
  setNames(
    lapply(pdo_groups$mp, function(mp) unique(as.character(pdo_genes[[mp]]))),
    paste0("PDO__", pdo_groups$mp)
  ),
  setNames(
    lapply(scref_groups$mp, function(mp) unique(as.character(scref_genes[[mp]]))),
    paste0("scRef__", scref_groups$mp)
  )
)
signature_genes <- signature_genes[signature_manifest$signature_key]
signature_genes <- lapply(
  signature_genes,
  function(genes) genes[!is.na(genes) & nzchar(genes)]
)

manifest_path <- file.path(
  out_tiers[["tables"]],
  "Auto_parse_external_centred_mp_signature_manifest.csv"
)
signature_path <- file.path(
  out_tiers[["intermediate"]],
  "Auto_parse_external_centred_mp_signatures.rds"
)
write.csv(signature_manifest, manifest_path, row.names = FALSE)
saveRDS(signature_genes, signature_path)
####################

####################
# Per-sample UCell scoring with persistent live caches
####################
score_path <- file.path(
  out_tiers[["intermediate"]],
  "Auto_parse_external_centred_mp_ucell_scores.rds"
)
metadata_path <- file.path(
  out_tiers[["intermediate"]],
  "Auto_parse_external_centred_mp_cell_metadata.rds"
)
coverage_path <- file.path(
  out_tiers[["tables"]],
  "Auto_parse_external_centred_mp_gene_coverage.csv"
)

get_counts <- function(obj) {
  counts <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts)) {
    counts <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(counts)) {
    stop("Could not read the RNA counts layer from a Parse Seurat object.")
  }
  counts
}

score_one_sample <- function(sample_name, sample_file, force_rebuild = FALSE) {
  sample_cache <- file.path(
    out_tiers[["intermediate"]],
    paste0("Auto_", sample_name, "_external_centred_mp_ucell_scores.rds")
  )
  sample_metadata_cache <- file.path(
    out_tiers[["intermediate"]],
    paste0("Auto_", sample_name, "_external_centred_mp_cell_metadata.rds")
  )
  sample_coverage_cache <- file.path(
    out_tiers[["intermediate"]],
    paste0("Auto_", sample_name, "_external_centred_mp_gene_coverage.rds")
  )

  cache_files <- c(sample_cache, sample_metadata_cache, sample_coverage_cache)
  if (!force_rebuild && all(file.exists(cache_files))) {
    message("Reusing live UCell cache for ", sample_name)
    return(list(
      scores = readRDS(sample_cache),
      metadata = readRDS(sample_metadata_cache),
      coverage = readRDS(sample_coverage_cache)
    ))
  }

  message("Loading Parse counts for ", sample_name, ": ", sample_file)
  obj <- readRDS(sample_file)
  counts <- get_counts(obj)
  old_cells <- colnames(counts)
  new_cells <- paste(sample_name, old_cells, sep = "_")
  colnames(counts) <- new_cells

  sample_features <- lapply(
    signature_genes,
    function(genes) intersect(genes, rownames(counts))
  )
  coverage <- data.frame(
    sample = sample_name,
    signature_key = names(sample_features),
    n_genes_reference = lengths(signature_genes),
    n_genes_detectable = lengths(sample_features),
    detectable_fraction = lengths(sample_features) / lengths(signature_genes),
    stringsAsFactors = FALSE
  )
  if (any(coverage$n_genes_detectable < 5)) {
    failed <- coverage$signature_key[coverage$n_genes_detectable < 5]
    stop(
      "Fewer than five detectable genes for ",
      sample_name,
      " signature(s): ",
      paste(failed, collapse = ", ")
    )
  }

  message(
    "Scoring ",
    length(sample_features),
    " external centred-refined MPs in ",
    sample_name,
    " (",
    ncol(counts),
    " cells)"
  )
  scored <- UCell::ScoreSignatures_UCell(
    matrix = counts,
    features = sample_features,
    maxRank = 1500,
    chunk.size = 1000,
    ncores = 6,
    force.gc = TRUE
  )
  scored <- as.matrix(scored)
  colnames(scored) <- sub("_UCell$", "", colnames(scored))
  missing_scores <- setdiff(names(signature_genes), colnames(scored))
  if (length(missing_scores) > 0) {
    stop(
      "UCell output is missing signature(s): ",
      paste(missing_scores, collapse = ", ")
    )
  }
  scored <- scored[, names(signature_genes), drop = FALSE]
  sample_metadata <- data.frame(
    cell = rownames(scored),
    original_cell = old_cells,
    sample = sample_name,
    stringsAsFactors = FALSE,
    row.names = rownames(scored)
  )

  saveRDS(scored, sample_cache, compress = FALSE)
  saveRDS(sample_metadata, sample_metadata_cache, compress = FALSE)
  saveRDS(coverage, sample_coverage_cache)
  rm(obj, counts)
  invisible(gc())

  list(scores = scored, metadata = sample_metadata, coverage = coverage)
}

if (cache_policy$replot_only) {
  pdo_require_files(c(score_path, metadata_path, coverage_path))
  message("PDO_REPLOT_ONLY=1: loading combined live UCell caches.")
  combined_scores <- readRDS(score_path)
  cell_metadata <- readRDS(metadata_path)
  coverage_table <- read.csv(
    coverage_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  sample_results <- lapply(parse_samples, function(sample_name) {
    score_one_sample(
      sample_name,
      parse_files[[sample_name]],
      force_rebuild = cache_policy$force_rebuild
    )
  })
  names(sample_results) <- parse_samples

  combined_scores <- do.call(
    rbind,
    lapply(sample_results, `[[`, "scores")
  )
  cell_metadata <- bind_rows(lapply(sample_results, `[[`, "metadata"))
  rownames(cell_metadata) <- cell_metadata$cell
  coverage_table <- bind_rows(lapply(sample_results, `[[`, "coverage"))

  if (!identical(rownames(combined_scores), cell_metadata$cell)) {
    stop("Combined UCell score rows and cell metadata are not aligned.")
  }
  saveRDS(combined_scores, score_path, compress = FALSE)
  saveRDS(cell_metadata, metadata_path, compress = FALSE)
  write.csv(coverage_table, coverage_path, row.names = FALSE)
  rm(sample_results)
  invisible(gc())
}

if (!all(signature_manifest$signature_key %in% colnames(combined_scores))) {
  stop("Combined UCell score columns do not contain all signatures in the manifest.")
}
combined_scores <- combined_scores[, signature_manifest$signature_key, drop = FALSE]
sample_vector <- factor(cell_metadata[rownames(combined_scores), "sample"], levels = parse_samples)
if (any(is.na(sample_vector))) {
  stop("Cell metadata contains missing or unexpected Parse sample labels.")
}
####################

####################
# Timepoint summaries and across-timepoint tests
####################
summarise_signature <- function(signature_key) {
  score <- combined_scores[, signature_key]
  bind_rows(lapply(parse_samples, function(sample_name) {
    values <- score[sample_vector == sample_name]
    data.frame(
      signature_key = signature_key,
      sample = sample_name,
      n_cells = sum(is.finite(values)),
      mean_score = mean(values, na.rm = TRUE),
      median_score = median(values, na.rm = TRUE),
      q05 = quantile(values, 0.05, na.rm = TRUE, names = FALSE),
      q1 = quantile(values, 0.25, na.rm = TRUE, names = FALSE),
      q3 = quantile(values, 0.75, na.rm = TRUE, names = FALSE),
      q95 = quantile(values, 0.95, na.rm = TRUE, names = FALSE),
      stringsAsFactors = FALSE
    )
  }))
}

timepoint_summary <- bind_rows(
  lapply(signature_manifest$signature_key, summarise_signature)
) |>
  left_join(
    signature_manifest |>
      select(
        signature_key,
        source,
        biological_group,
        mp,
        description,
        display_label,
        row_order
      ),
    by = "signature_key"
  ) |>
  mutate(
    sample = factor(sample, levels = parse_samples),
    biological_group = factor(biological_group, levels = group_levels)
  ) |>
  arrange(row_order, sample)

downsample_balanced <- function(values, target_n, seed_key) {
  values <- values[is.finite(values)]
  if (length(values) <= target_n) {
    return(values)
  }
  set.seed(sum(utf8ToInt(seed_key)) %% .Machine$integer.max)
  sample(values, target_n, replace = FALSE)
}

test_signature <- function(signature_key) {
  score <- combined_scores[, signature_key]
  response <- score[sample_vector %in% c("T2", "T4")]
  reference <- score[sample_vector %in% c("T0", "eR4")]
  balanced_n <- min(sum(is.finite(response)), sum(is.finite(reference)))
  response_balanced <- downsample_balanced(
    response,
    balanced_n,
    paste0(signature_key, "_T2T4")
  )
  reference_balanced <- downsample_balanced(
    reference,
    balanced_n,
    paste0(signature_key, "_T0eR4")
  )

  data.frame(
    signature_key = signature_key,
    kruskal_p = tryCatch(
      kruskal.test(score ~ sample_vector)$p.value,
      error = function(e) NA_real_
    ),
    comparison = "T2+T4 versus T0+eR4",
    balanced_n_per_group = balanced_n,
    response_mean = mean(response_balanced, na.rm = TRUE),
    reference_mean = mean(reference_balanced, na.rm = TRUE),
    mean_difference = mean(response_balanced, na.rm = TRUE) -
      mean(reference_balanced, na.rm = TRUE),
    response_median = median(response_balanced, na.rm = TRUE),
    reference_median = median(reference_balanced, na.rm = TRUE),
    median_difference = median(response_balanced, na.rm = TRUE) -
      median(reference_balanced, na.rm = TRUE),
    wilcox_p = tryCatch(
      wilcox.test(
        response_balanced,
        reference_balanced,
        alternative = "two.sided",
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    stringsAsFactors = FALSE
  )
}

timepoint_tests <- bind_rows(
  lapply(signature_manifest$signature_key, test_signature)
) |>
  mutate(
    kruskal_p_adj = p.adjust(kruskal_p, method = "BH"),
    wilcox_p_adj = p.adjust(wilcox_p, method = "BH")
  ) |>
  left_join(
    signature_manifest |>
      select(
        signature_key,
        source,
        biological_group,
        mp,
        description,
        display_label,
        row_order
      ),
    by = "signature_key"
  ) |>
  arrange(row_order)

summary_path <- file.path(
  out_tiers[["tables"]],
  "Auto_parse_external_centred_mp_timepoint_summary.csv"
)
test_path <- file.path(
  out_tiers[["tables"]],
  "Auto_parse_external_centred_mp_timepoint_tests.csv"
)
write.csv(timepoint_summary, summary_path, row.names = FALSE)
write.csv(timepoint_tests, test_path, row.names = FALSE)
####################

####################
# Ordered heatmaps: raw mean UCell enrichment and within-MP temporal z-score
####################
mean_wide <- timepoint_summary |>
  select(signature_key, sample, mean_score) |>
  pivot_wider(names_from = sample, values_from = mean_score)
mean_wide <- as.data.frame(mean_wide)
rownames(mean_wide) <- mean_wide$signature_key
mean_matrix <- as.matrix(
  mean_wide[signature_manifest$signature_key, parse_samples, drop = FALSE]
)
z_matrix <- t(scale(t(mean_matrix)))
z_matrix[!is.finite(z_matrix)] <- 0

row_annotation <- signature_manifest |>
  transmute(
    Source = as.character(source),
    States = as.character(biological_group)
  )
rownames(row_annotation) <- signature_manifest$signature_key

annotation_colours <- list(
  Source = source_colours,
  States = group_colours
)

left_anno <- rowAnnotation(
  States = row_annotation$States,
  col = list(States = annotation_colours$States),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

right_anno <- rowAnnotation(
  Source = row_annotation$Source,
  col = list(Source = annotation_colours$Source),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

make_heatmap <- function(mat, title, col_fun, breaks = NA) {
  Heatmap(
    mat,
    name = ifelse(grepl("z-score", title), "z-score", "UCell"),
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(row_annotation$States, levels = group_levels),
    row_title = NULL,
    row_gap = unit(4, "mm"),
    left_annotation = left_anno,
    right_annotation = right_anno,
    row_labels = signature_manifest$display_label,
    column_title = title,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_names_gp = gpar(fontsize = 10.5),
    column_names_gp = gpar(fontsize = 14),
    column_names_rot = 0,
    column_names_centered = TRUE,
    rect_gp = gpar(col = "white", lwd = 1)
  )
}

raw_col_fun <- colorRamp2(
  seq(min(mean_matrix, na.rm=TRUE), max(mean_matrix, na.rm=TRUE), length.out = 9),
  c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
)

z_limit <- max(2, min(3, max(abs(z_matrix), na.rm = TRUE)))
z_col_fun <- colorRamp2(
  seq(-z_limit, z_limit, length.out = 11),
  rev(RColorBrewer::brewer.pal(11, "RdBu"))
)

raw_heatmap <- make_heatmap(
  mean_matrix,
  "Raw Ucell score across Parse timepoints",
  raw_col_fun
)
z_heatmap <- make_heatmap(
  z_matrix,
  "Normalised z-score across Parse timepoints",
  z_col_fun
)

raw_heatmap_pdf <- file.path(
  out_tiers[["figures"]],
  "Auto_parse_external_centred_mp_timecourse_raw_mean_heatmap.pdf"
)
z_heatmap_pdf <- file.path(
  out_tiers[["figures"]],
  "Auto_parse_external_centred_mp_timecourse_zscore_heatmap.pdf"
)
z_heatmap_png <- file.path(
  out_tiers[["figures"]],
  "Auto_parse_external_centred_mp_timecourse_zscore_heatmap.png"
)

pdf(raw_heatmap_pdf, width = 18, height = 20, useDingbats = FALSE)
draw(raw_heatmap, merge_legend = TRUE)
dev.off()

pdf(z_heatmap_pdf, width = 18, height = 20, useDingbats = FALSE)
draw(z_heatmap, merge_legend = TRUE)
dev.off()

png(z_heatmap_png, width = 5400, height = 6000, res = 300)
draw(z_heatmap, merge_legend = TRUE)
dev.off()
####################

####################
# Group-specific timepoint trend and distribution-summary pages
####################
make_group_trend <- function(group_name) {
  plot_df <- timepoint_summary |>
    filter(biological_group == group_name) |>
    mutate(
      display_label = factor(
        display_label,
        levels = signature_manifest$display_label[
          signature_manifest$biological_group == group_name
        ]
      )
    )

  ggplot(
    plot_df,
    aes(
      x = sample,
      y = mean_score,
      group = signature_key,
      colour = source
    )
  ) +
    geom_ribbon(
      aes(ymin = q1, ymax = q3, fill = source),
      alpha = 0.13,
      colour = NA
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(
      aes(fill = sample),
      shape = 21,
      colour = "black",
      size = 2.8,
      stroke = 0.35
    ) +
    facet_wrap(
      ~display_label,
      scales = "free_y",
      ncol = 3,
      labeller = labeller(display_label = label_wrap_gen(width = 42))
    ) +
    scale_colour_manual(values = source_colours, drop = FALSE) +
    scale_fill_manual(
      values = c(source_colours, parse_colours),
      breaks = parse_samples,
      name = "Timepoint",
      drop = FALSE
    ) +
    labs(
      title = paste0(group_name, ": mean UCell enrichment across Parse timepoints"),
      x = NULL,
      y = "UCell score"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        colour = "black",
        size = 13
      ),
      axis.text.y = element_text(colour = "black", size = 12),
      strip.text = element_text(face = "bold", size = 12, lineheight = 1.05),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13),
      panel.spacing = unit(1.35, "lines")
    ) +
    guides(colour = guide_legend(title = "MP source", order = 1))
}

make_group_distribution <- function(group_name) {
  plot_df <- timepoint_summary |>
    filter(biological_group == group_name) |>
    mutate(
      display_label = factor(
        display_label,
        levels = signature_manifest$display_label[
          signature_manifest$biological_group == group_name
        ]
      )
    )

  ggplot(plot_df, aes(x = sample, colour = source)) +
    geom_linerange(aes(ymin = q05, ymax = q95), linewidth = 0.55) +
    geom_crossbar(
      aes(y = median_score, ymin = q1, ymax = q3, fill = sample),
      width = 0.62,
      alpha = 0.8,
      colour = "black",
      linewidth = 0.35
    ) +
    facet_wrap(
      ~display_label,
      scales = "free_y",
      ncol = 3,
      labeller = labeller(display_label = label_wrap_gen(width = 42))
    ) +
    scale_colour_manual(values = source_colours, drop = FALSE) +
    scale_fill_manual(values = parse_colours, drop = FALSE, name = "Timepoint") +
    labs(
      title = paste0(group_name, ": Parse-cell UCell score distributions"),
      x = NULL,
      y = "UCell score"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", size = 20),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        colour = "black",
        size = 13
      ),
      axis.text.y = element_text(colour = "black", size = 12),
      strip.text = element_text(face = "bold", size = 12, lineheight = 1.05),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13),
      panel.spacing = unit(1.35, "lines")
    ) +
    guides(colour = "none")
}

####################
# State and MP abundance visualization (Approach B noreg)
####################
mp_desc <- PDO_MP_DESCRIPTIONS

label_f <- function(mps) {
  mps_clean <- sub("^PDO__", "", mps)
  d <- mp_desc[mps_clean]
  d[is.na(d)] <- mps_clean[is.na(d)]
  out <- paste0(mps_clean, "_", d)
  names(out) <- names(mps)
  out
}

mp_cols <- c(
  "MP11_Single-nucleus-associated cell cycle"  = "#E78AC3",
  "MP1_G2/M cell cycle"                        = "#B3B3B3",
  "MP2_G1/S cell cycle"                        = "#8DA0CB",
  "MP3_Replication-dependent histones"         = "#999999",
  "MP19+_MYC-associated proliferation"         = "#E41A1C",
  "MP15_Intestinal metaplasia"                 = "#4DAF4A",
  "MP5+_Inflammatory-reactive columnar epithelium" = "#984EA3",
  "MP12_KRAS-active columnar epithelium"       = "#74C476",
  "MP13b_Metabolic-detox columnar epithelium"  = "#BAE4B3",
  "MP14b_Proliferative epithelial plasticity"  = "#A1D99B",
  "MP16b_EMT/KRAS adaptive plasticity"         = "#377EB8",
  "MP17+_Ciliated progenitor epithelium"       = "#FDBF6F",
  "MP8+_Secretory-transport glandular epithelium" = "#FF7F00",
  "MP9_ECM-remodelling epithelium"             = "#A6D854",
  "MP18_Motile-cilia differentiation"          = "#1F78B4"
)

state_groups <- list(
  "Classic proliferation" = c("PDO__MP19+"),
  "Columnar-to-intestinal" = c("PDO__MP14b", "PDO__MP13b", "PDO__MP5+", "PDO__MP12", "PDO__MP15"),
  "Glandular differentiation" = c("PDO__MP17+", "PDO__MP8+"),
  "Stress-adaptive" = c("PDO__MP16b"),
  "ECM-remodelling" = c("PDO__MP9"),
  "Motile-cilia differentiation" = c("PDO__MP18")
)
state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Motile-cilia differentiation" = "#F781BF",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)
state_level_order <- c(names(state_groups), "Unresolved", "Hybrid")

pdo_all_keys <- paste0("PDO__", names(mp_desc))
pdo_cc_keys <- paste0("PDO__", PDO_CELL_CYCLE_MPS)
pdo_noncc_keys <- unlist(state_groups, use.names = FALSE)

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

make_prop_data <- function(label_vec, cell_meta, label_order, sample_order) {
  data.frame(
    sample = cell_meta[names(label_vec), "sample"],
    label = as.character(label_vec),
    stringsAsFactors = FALSE
  ) %>%
    filter(sample %in% sample_order, label %in% label_order) %>%
    count(sample, label, name = "n") %>%
    right_join(tidyr::expand_grid(sample = sample_order, label = label_order), by = c("sample", "label")) %>%
    mutate(n = replace_na(n, 0L)) %>%
    group_by(sample) %>%
    mutate(pct = 100 * n / pmax(sum(n), 1)) %>%
    ungroup()
}

plot_abundance <- function(label_vec, cell_meta, label_order, col_map, sample_order, title_text) {
  totals_df <- cell_meta %>%
    filter(sample %in% sample_order) %>%
    count(sample, name = "total_cells") %>%
    mutate(sample = factor(sample, levels = sample_order))
  scale_factor <- max(totals_df$total_cells, na.rm = TRUE) / 100
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1

  prop_df <- make_prop_data(label_vec, cell_meta, label_order, sample_order) %>%
    mutate(
      sample = factor(sample, levels = sample_order),
      label = factor(label, levels = rev(label_order))
    )

  ggplot(prop_df, aes(x = sample, y = pct, fill = label)) +
    geom_col(width = 0.78, colour = NA) +
    geom_point(
      data = totals_df,
      aes(x = sample, y = total_cells / scale_factor),
      inherit.aes = FALSE,
      colour = "black",
      size = 2.4
    ) +
    geom_line(
      data = totals_df,
      aes(x = sample, y = total_cells / scale_factor, group = 1),
      inherit.aes = FALSE,
      colour = "black",
      alpha = 0.45,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    scale_fill_manual(values = col_map, breaks = label_order, drop = FALSE, name = NULL) +
    scale_y_continuous(
      name = "Proportion (%)",
      limits = c(0, 100),
      expand = c(0, 0),
      sec.axis = sec_axis(~ . * scale_factor, name = "Total Cells", labels = scales::comma)
    ) +
    labs(x = NULL, title = title_text) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 17, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 15, colour = "black"),
      axis.title.y = element_text(size = 17, face = "bold"),
      axis.title.y.right = element_text(size = 17, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 13),
      legend.key.size = unit(0.55, "cm")
    )
}

sample_var <- cell_metadata$sample
names(sample_var) <- rownames(cell_metadata)
study_var <- sample_var

# Z-normalise all PDO MPs
pdo_adj_all <- z_normalise(as.matrix(combined_scores[, pdo_all_keys, drop = FALSE]), sample_var, study_var)
mp_adj_all <- pdo_adj_all[, pdo_noncc_keys, drop = FALSE]

# 1) State assignment
group_max <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj_all))
  if (length(mps_avail) == 0) return(rep(NA_real_, nrow(mp_adj_all)))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj_all[, mps_avail]))
  apply(mp_adj_all[, mps_avail, drop = FALSE], 1, max)
})
group_max <- as.matrix(group_max)
rownames(group_max) <- rownames(mp_adj_all)
group_max[!is.finite(group_max)] <- NA_real_

threshold <- 0.5
hybrid_gap <- 0.3

best_group_idx <- max.col(group_max, ties.method = "first")
best_group_val <- apply(group_max, 1, max, na.rm = TRUE)
state_vec <- names(state_groups)[best_group_idx]
state_vec[!is.finite(best_group_val) | best_group_val < threshold] <- "Unresolved"

sorted_groups <- t(apply(group_max, 1, sort, decreasing = TRUE))
gap <- sorted_groups[, 1] - sorted_groups[, 2]
state_vec[(gap < hybrid_gap) & (state_vec != "Unresolved")] <- "Hybrid"
names(state_vec) <- rownames(group_max)

abundance_state_plot <- plot_abundance(
  label_vec = state_vec,
  cell_meta = cell_metadata,
  label_order = state_level_order,
  col_map = state_cols,
  sample_order = parse_samples,
  title_text = "PDO Centred Refined State Abundance Across Parse Timepoints"
)

# 2) Top Non-CC MP abundance
topmp_noncc_keys <- colnames(mp_adj_all)[max.col(mp_adj_all, ties.method = "first")]
topmp_noncc_vec <- label_f(topmp_noncc_keys)
names(topmp_noncc_vec) <- rownames(mp_adj_all)
ord_noncc <- label_f(pdo_noncc_keys)
col_noncc <- mp_cols[ord_noncc]

abundance_noncc_plot <- plot_abundance(
  label_vec = topmp_noncc_vec,
  cell_meta = cell_metadata,
  label_order = ord_noncc,
  col_map = col_noncc,
  sample_order = parse_samples,
  title_text = "PDO Centred Refined Non-CC Top MP Abundance Across Parse Timepoints"
)

# 3) Top All MP abundance
topmp_all_keys <- colnames(pdo_adj_all)[max.col(pdo_adj_all, ties.method = "first")]
topmp_all_vec <- label_f(topmp_all_keys)
names(topmp_all_vec) <- rownames(pdo_adj_all)
ord_all <- c(label_f(pdo_cc_keys), ord_noncc)
col_all <- mp_cols[ord_all]

abundance_all_plot <- plot_abundance(
  label_vec = topmp_all_vec,
  cell_meta = cell_metadata,
  label_order = ord_all,
  col_map = col_all,
  sample_order = parse_samples,
  title_text = "PDO Centred Refined All Top MP Abundance Across Parse Timepoints"
)
####################

####################
# Multi-page Report and Standalone Figures
####################
groups_present <- group_levels[
  group_levels %in% as.character(signature_manifest$biological_group)
]
report_pdf <- file.path(
  out_tiers[["figures"]],
  "Auto_parse_external_centred_mp_timecourse_report.pdf"
)
pdf(report_pdf, width = 18, height = 14, onefile = TRUE, useDingbats = FALSE)
draw(z_heatmap, merge_legend = TRUE)
draw(raw_heatmap, merge_legend = TRUE)
for (group_name in groups_present) {
  print(make_group_trend(group_name))
  print(make_group_distribution(group_name))
}
print(abundance_state_plot)
print(abundance_noncc_plot)
print(abundance_all_plot)
dev.off()

abundance_pdf <- file.path(
  out_tiers[["figures"]],
  "Auto_parse_external_centred_mp_timecourse_abundance.pdf"
)
pdf(abundance_pdf, width = 18, height = 14, onefile = TRUE, useDingbats = FALSE)
print(abundance_state_plot)
print(abundance_noncc_plot)
print(abundance_all_plot)
dev.off()
####################

####################
# Persistent plot matrices and run summary
####################
matrix_path <- file.path(
  out_tiers[["intermediate"]],
  "Auto_parse_external_centred_mp_timepoint_matrices.rds"
)
saveRDS(
  list(
    mean_ucell = mean_matrix,
    temporal_zscore = z_matrix,
    signature_manifest = signature_manifest,
    sample_order = parse_samples,
    biological_group_order = group_levels,
    source_colours = source_colours,
    biological_group_colours = group_colours
  ),
  matrix_path
)

output_files <- c(
  signature_path,
  score_path,
  metadata_path,
  matrix_path,
  manifest_path,
  coverage_path,
  summary_path,
  test_path,
  raw_heatmap_pdf,
  z_heatmap_pdf,
  z_heatmap_png,
  report_pdf,
  abundance_pdf
)
missing_outputs <- output_files[!file.exists(output_files)]
if (length(missing_outputs) > 0) {
  stop(
    "Expected output file(s) were not produced: ",
    paste(missing_outputs, collapse = ", ")
  )
}

run_log <- pdo_write_run_summary(
  script = script_path,
  out_dir = out_dir,
  inputs = input_files,
  outputs = output_files,
  parameters = list(
    parse_samples = paste(parse_samples, collapse = ","),
    pdo_signature_n = sum(signature_manifest$source == "PDO"),
    scref_signature_n = sum(signature_manifest$source == "scRef"),
    scref_excluded_omitted = "MP11c,MP18a",
    ucell_max_rank = 1500,
    ucell_ncores = 6,
    start_time = format(script_start, "%Y-%m-%d %H:%M:%S %Z"),
    end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  cache = cache_policy,
  status = "completed",
  include_session = TRUE
)

message("Completed external centred MP scoring across Parse timepoints.")
message("Report: ", report_pdf)
message("Run log: ", run_log)
####################
