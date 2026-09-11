####################
# Analysis registry:
#   Status: active terminal centred MP/state selective-dependency workflow
#   Script: analysis/cell_states/Auto_crispr_selective_dependency_analysis.R
#   Methodology: analysis/methodology/cell_states/Auto_crispr_selective_dependency_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Tests canonical centred PDO and scATLAS MP genes, combined state-component
#     MP genes, and state DGE/marker signatures against the significant
#     oesophageal differential-dependency genes in Herranz-Ors et al. (2026).
#   Inputs:
#     /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/CRISPR_screen_results_PDO.xlsx
#     PDOs_outs/PDOs_merged.rds
#     PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     PDOs_outs/centred_mp_refinement/merged_refined_mp_gene_weights.rds
#     /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds
#     /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_gene_weights.rds
#     /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds
#     /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds
#     /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv
#     /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R
#     analysis/shared/Auto_pdo_analysis_config.R
#   Outputs:
#     PDOs_outs/Auto_crispr_selective_dependency/intermediate/*.rds
#     PDOs_outs/Auto_crispr_selective_dependency/tables/*.csv*
#     PDOs_outs/Auto_crispr_selective_dependency/tables/Auto_crispr_selective_dependency_results.xlsx
#     PDOs_outs/Auto_crispr_selective_dependency/figures/*.pdf
#     PDOs_outs/Auto_crispr_selective_dependency/logs/*.txt
#     PDOs_outs/Auto_crispr_selective_dependency/reports/Auto_crispr_selective_dependency_methodology.md
#   Downstream:
#     Terminal prioritization of selectively targetable centred MPs/states.
#     Does not redefine canonical MPs or states.
#   Cache/replot behavior:
#     Set PDO_FORCE_REBUILD=1 to rebuild edgeR results and scoring cache.
#     Set PDO_REPLOT_ONLY=1 to require and reuse the existing scoring cache.
#   Run command: qsub Auto_run_crispr_selective_dependency.sh
#   Conda environment: dmtcp
####################

####################
# Libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(readxl)
  library(openxlsx)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

####################
# Setup and fixed inputs
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis", "shared", "Auto_pdo_analysis_helpers.R"))

scref_config_file <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline",
  "analysis/shared/scRef_config.R"
)
source(scref_config_file)

out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_crispr_selective_dependency")
out_paths <- pdo_ensure_output_tiers(out_dir)
set.seed(10830)

workbook_file <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix",
  "CRISPR_screen_results_PDO.xlsx"
)
seurat_file <- file.path(PDO_OUTPUT_DIR, "PDOs_merged.rds")
state_file <- file.path(PDO_OUTPUT_DIR, PDO_PREFERRED_STATE_VECTOR)
mp_gene_file <- file.path(PDO_OUTPUT_DIR, PDO_PREFERRED_MP_GENES)
mp_weight_file <- file.path(
  PDO_OUTPUT_DIR,
  "centred_mp_refinement",
  "merged_refined_mp_gene_weights.rds"
)
scatlas_mp_gene_file <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline",
  "ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate",
  "merged_refined_mp_genes.rds"
)
scatlas_mp_weight_file <- file.path(
  SCREF_CENTRED_MP_DIR,
  "intermediate",
  "merged_refined_mp_gene_weights.rds"
)
scatlas_seurat_file <- SCREF_EPI_RDS
scatlas_state_file <- SCREF_STATE_NOREG_RDS
scatlas_state_marker_file <- file.path(
  "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline",
  "ref_outs/Metaprogrammes_Results/centred/state_markers",
  "Auto_five_state_markers_ranked.csv"
)
methodology_file <- file.path(
  project_dir,
  "analysis",
  "methodology",
  "cell_states",
  "Auto_crispr_selective_dependency_methodology.md"
)

pdo_require_files(
  c(
    workbook_file, seurat_file, state_file, mp_gene_file, mp_weight_file,
    scatlas_mp_gene_file, scatlas_mp_weight_file, scatlas_seurat_file,
    scatlas_state_file, scatlas_state_marker_file, scref_config_file,
    methodology_file
  )
)

force_rebuild <- pdo_get_env_flag("PDO_FORCE_REBUILD", FALSE)
replot_only <- pdo_get_env_flag("PDO_REPLOT_ONLY", FALSE)
params <- list(
  reference_sheet = "diff_dep_analysis_oesophageal",
  reference_fdr_guard = 0.05,
  dge_top_n = 150L,
  scatlas_dge_top_n = 150L,
  min_detected_cells = 10L,
  min_cells_per_pseudobulk = 20L,
  min_paired_samples = 3L,
  expression_bins = 10L,
  matched_permutations = 10000L,
  empirical_fdr = 0.10
)

dge_cache_file <- file.path(
  out_paths[["intermediate"]],
  "Auto_state_edger_dge_all.rds"
)
score_cache_file <- file.path(
  out_paths[["intermediate"]],
  "Auto_crispr_selective_dependency_scoring_cache.rds"
)
scatlas_score_cache_file <- file.path(
  out_paths[["intermediate"]],
  "Auto_scatlas_crispr_selective_dependency_scoring_cache.rds"
)

####################
# Helpers
####################
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

identifier_keys <- function(x) {
  x <- trimws(as.character(x))
  upper <- toupper(x)
  data.frame(
    symbol_key = upper,
    ensembl_key = ifelse(
      grepl("^ENSG[0-9]+", upper),
      sub("\\..*$", "", upper),
      NA_character_
    ),
    stringsAsFactors = FALSE
  )
}

match_identifiers <- function(x, target_symbol_key, target_ensembl_key) {
  keys <- identifier_keys(x)
  symbol_idx <- match(keys$symbol_key, target_symbol_key)
  ensembl_idx <- match(keys$ensembl_key, target_ensembl_key)
  ifelse(!is.na(symbol_idx), symbol_idx, ensembl_idx)
}

####################
# Classic-proliferation three-way overlap helpers
####################
canonicalize_gene_set <- function(x, reference_df) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(trimws(x))]
  keys <- identifier_keys(x)
  reference_idx <- match_identifiers(
    x,
    reference_df$gene_key,
    reference_df$ensembl_key
  )
  canonical <- ifelse(
    !is.na(reference_idx),
    reference_df$gene_key[reference_idx],
    ifelse(!is.na(keys$ensembl_key), keys$ensembl_key, keys$symbol_key)
  )
  sort(unique(canonical[!is.na(canonical) & nzchar(canonical)]))
}

build_threeway_membership <- function(panel_name, scatlas_genes, crispr_genes, pdo_genes) {
  gene_union <- sort(unique(c(scatlas_genes, crispr_genes, pdo_genes)))
  data.frame(
    panel = panel_name,
    gene = gene_union,
    in_scatlas = gene_union %in% scatlas_genes,
    in_crispr = gene_union %in% crispr_genes,
    in_pdo = gene_union %in% pdo_genes,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      region = case_when(
        in_scatlas & in_crispr & in_pdo ~ "scATLAS & CRISPR & PDO",
        in_scatlas & in_crispr ~ "scATLAS & CRISPR only",
        in_scatlas & in_pdo ~ "scATLAS & PDO only",
        in_crispr & in_pdo ~ "CRISPR & PDO only",
        in_scatlas ~ "scATLAS only",
        in_crispr ~ "CRISPR only",
        in_pdo ~ "PDO only",
        TRUE ~ NA_character_
      )
    )
}

safe_fisher_enrichment <- function(hit_n, set_n, reference_n, universe_n) {
  if (set_n <= 0 || universe_n <= 0 || reference_n <= 0) {
    return(c(odds_ratio = NA_real_, p_value = NA_real_))
  }
  outside_hit <- reference_n - hit_n
  outside_nonhit <- universe_n - set_n - outside_hit
  if (any(c(hit_n, set_n - hit_n, outside_hit, outside_nonhit) < 0)) {
    return(c(odds_ratio = NA_real_, p_value = NA_real_))
  }
  result <- fisher.test(
    matrix(
      c(hit_n, set_n - hit_n, outside_hit, outside_nonhit),
      nrow = 2,
      byrow = TRUE
    ),
    alternative = "greater"
  )
  c(odds_ratio = unname(result$estimate), p_value = result$p.value)
}

run_pseudobulk_state_dge <- function(counts_mat, meta, state_name, all_states, params) {
  state_meta <- meta %>%
    filter(state %in% all_states) %>%
    mutate(
      pb_group = ifelse(state == state_name, "target", "rest"),
      pb_id = paste(orig.ident, pb_group, sep = "___")
    )

  cells_use <- intersect(state_meta$cell, colnames(counts_mat))
  state_meta <- state_meta[match(cells_use, state_meta$cell), , drop = FALSE]
  counts_use <- counts_mat[, cells_use, drop = FALSE]

  pb_design <- Matrix::sparse.model.matrix(~ 0 + factor(state_meta$pb_id))
  colnames(pb_design) <- sub("^factor\\(state_meta\\$pb_id\\)", "", colnames(pb_design))
  pb_counts <- counts_use %*% pb_design

  pb_meta <- data.frame(
    pb_id = colnames(pb_counts),
    sample = sub("___(target|rest)$", "", colnames(pb_counts)),
    group = sub("^.*___", "", colnames(pb_counts)),
    n_cells = as.integer(table(state_meta$pb_id)[colnames(pb_counts)]),
    stringsAsFactors = FALSE
  ) %>%
    filter(n_cells >= params$min_cells_per_pseudobulk)

  paired_samples <- pb_meta %>%
    count(sample, name = "n_groups") %>%
    filter(n_groups == 2L) %>%
    pull(sample)
  pb_meta <- pb_meta %>%
    filter(sample %in% paired_samples) %>%
    arrange(sample, group)

  if (length(unique(pb_meta$sample)) < params$min_paired_samples) {
    stop(
      "State ", state_name, " has only ", length(unique(pb_meta$sample)),
      " eligible paired pseudobulk samples; minimum is ", params$min_paired_samples, "."
    )
  }

  pb_counts <- as.matrix(pb_counts[, pb_meta$pb_id, drop = FALSE])
  colnames(pb_counts) <- pb_meta$pb_id
  pb_meta$group <- factor(pb_meta$group, levels = c("rest", "target"))
  pb_meta$sample <- factor(pb_meta$sample)
  design <- model.matrix(~ sample + group, data = pb_meta)

  if (!"grouptarget" %in% colnames(design) || qr(design)$rank < ncol(design)) {
    stop("Sample-blocked pseudobulk design is not estimable for ", state_name, ".")
  }

  dge <- DGEList(counts = pb_counts, group = pb_meta$group)
  keep <- filterByExpr(dge, design = design)
  if (!any(keep)) stop("edgeR filterByExpr retained no genes for ", state_name, ".")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = design, robust = TRUE)
  fit <- glmQLFit(dge, design = design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "grouptarget")

  result <- topTags(qlf, n = Inf, sort.by = "none")$table %>%
    rownames_to_column("gene") %>%
    as_tibble()

  paired_cells <- state_meta %>% filter(orig.ident %in% paired_samples)
  target_cells <- paired_cells$cell[paired_cells$state == state_name]
  rest_cells <- paired_cells$cell[paired_cells$state != state_name]
  pct_target <- Matrix::rowMeans(counts_mat[result$gene, target_cells, drop = FALSE] > 0)
  pct_rest <- Matrix::rowMeans(counts_mat[result$gene, rest_cells, drop = FALSE] > 0)

  result %>%
    transmute(
      state = state_name,
      gene,
      avg_log2FC = logFC,
      p_val = PValue,
      p_val_adj = FDR,
      pct_target = as.numeric(pct_target[gene]),
      pct_rest = as.numeric(pct_rest[gene]),
      logCPM,
      F,
      n_paired_samples = n_distinct(pb_meta$sample),
      n_target_cells = length(target_cells),
      n_rest_cells = length(rest_cells),
      dge_method = "pseudobulk_edgeR_QLF_sample_blocked"
    ) %>%
    arrange(p_val_adj, desc(avg_log2FC), desc(pct_target), gene)
}

collapse_signature_to_universe <- function(signature_df, universe_df) {
  universe_idx <- match_identifiers(
    signature_df$gene,
    universe_df$symbol_key,
    universe_df$ensembl_key
  )
  signature_df$universe_idx <- universe_idx
  mapped <- signature_df %>%
    filter(!is.na(universe_idx)) %>%
    group_by(universe_idx) %>%
    summarise(
      gene = paste(unique(gene), collapse = ";"),
      definition_weight = sum(definition_weight, na.rm = TRUE),
      component_mps = paste(unique(component_mps[!is.na(component_mps)]), collapse = ";"),
      .groups = "drop"
    )
  mapped$definition_weight[!is.finite(mapped$definition_weight) | mapped$definition_weight <= 0] <- 1
  mapped$definition_weight <- mapped$definition_weight / sum(mapped$definition_weight)
  mapped
}

score_signature <- function(signature_name, signature_type, signature_df, universe_df, reference_df, n_perm) {
  reference_idx_all <- match_identifiers(
    signature_df$gene,
    reference_df$gene_key,
    reference_df$ensembl_key
  )
  signature_df$reference_idx <- reference_idx_all
  signature_df$is_selective_dependency <- !is.na(reference_idx_all)
  signature_df$definition_weight[!is.finite(signature_df$definition_weight) |
                                   signature_df$definition_weight <= 0] <- 1
  signature_df$definition_weight_normalized <-
    signature_df$definition_weight / sum(signature_df$definition_weight)

  ref_columns <- c(
    "gene", "ensembl_id", "mean_lfc", "median_lfc",
    "mean_lfc_not_depleted_group", "mean_lfc_depleted_group",
    "median_lfc_not_depleted_group", "median_lfc_depleted_group",
    "delta_mean", "delta_median", "n_not_depleted", "n_depleted",
    "p_value", "p_value_adj", "depletion_prevalence",
    "conditional_depletion_strength", "prevalence_weighted_strength",
    "overall_depletion_strength"
  )
  for (column_name in ref_columns) {
    signature_df[[paste0("crispr_", column_name)]] <-
      reference_df[[column_name]][reference_idx_all]
  }
  signature_df$weighted_gene_burden <- ifelse(
    signature_df$is_selective_dependency,
    signature_df$definition_weight_normalized *
      signature_df$crispr_prevalence_weighted_strength,
    0
  )

  mapped <- collapse_signature_to_universe(signature_df, universe_df)
  mapped_universe <- universe_df[mapped$universe_idx, , drop = FALSE]
  weights <- mapped$definition_weight
  observed_overlap <- sum(weights * mapped_universe$is_selective_dependency)
  observed_burden <- sum(weights * mapped_universe$prevalence_weighted_strength)
  observed_overall_lfc_burden <- sum(weights * mapped_universe$overall_depletion_strength)

  null_overlap <- numeric(n_perm)
  null_burden <- numeric(n_perm)
  null_overall_lfc_burden <- numeric(n_perm)
  bins <- mapped_universe$expression_bin

  for (iteration in seq_len(n_perm)) {
    sampled_idx <- integer(nrow(mapped))
    for (bin_value in unique(bins)) {
      positions <- which(bins == bin_value)
      pool <- which(universe_df$expression_bin == bin_value)
      sampled_idx[positions] <- sample(
        pool,
        size = length(positions),
        replace = length(pool) < length(positions)
      )
    }
    null_overlap[iteration] <- sum(
      weights * universe_df$is_selective_dependency[sampled_idx]
    )
    null_burden[iteration] <- sum(
      weights * universe_df$prevalence_weighted_strength[sampled_idx]
    )
    null_overall_lfc_burden[iteration] <- sum(
      weights * universe_df$overall_depletion_strength[sampled_idx]
    )
  }

  overlap_sd <- sd(null_overlap)
  burden_sd <- sd(null_burden)
  overlap_z <- ifelse(overlap_sd > 0, (observed_overlap - mean(null_overlap)) / overlap_sd, NA_real_)
  burden_z <- ifelse(burden_sd > 0, (observed_burden - mean(null_burden)) / burden_sd, NA_real_)
  overlap_p <- (1 + sum(null_overlap >= observed_overlap)) / (n_perm + 1)
  burden_p <- (1 + sum(null_burden >= observed_burden)) / (n_perm + 1)
  dependency_index <- 100 * (sum(null_burden <= observed_burden) + 0.5) / (n_perm + 1)

  set_n <- nrow(mapped)
  hit_n <- sum(mapped_universe$is_selective_dependency)
  fisher_result <- safe_fisher_enrichment(
    hit_n = hit_n,
    set_n = set_n,
    reference_n = sum(universe_df$is_selective_dependency),
    universe_n = nrow(universe_df)
  )

  summary_row <- data.frame(
    signature_type,
    signature = signature_name,
    signature_gene_n = nrow(signature_df),
    expression_universe_gene_n = set_n,
    expression_universe_coverage = set_n / nrow(signature_df),
    overlap_gene_n = sum(signature_df$is_selective_dependency),
    overlap_fraction = mean(signature_df$is_selective_dependency),
    definition_weighted_overlap_fraction = sum(
      signature_df$definition_weight_normalized * signature_df$is_selective_dependency
    ),
    median_depletion_prevalence_among_overlap = ifelse(
      any(signature_df$is_selective_dependency),
      median(signature_df$crispr_depletion_prevalence[signature_df$is_selective_dependency]),
      NA_real_
    ),
    median_depleted_group_strength_among_overlap = ifelse(
      any(signature_df$is_selective_dependency),
      median(signature_df$crispr_conditional_depletion_strength[signature_df$is_selective_dependency]),
      NA_real_
    ),
    prevalence_weighted_depletion_burden = sum(signature_df$weighted_gene_burden),
    expression_matched_overlap_z = overlap_z,
    expression_matched_overlap_p = overlap_p,
    expression_matched_burden_z = burden_z,
    expression_matched_burden_p = burden_p,
    dependency_index_0_100 = dependency_index,
    overall_lfc_burden = observed_overall_lfc_burden,
    fisher_odds_ratio = fisher_result[["odds_ratio"]],
    fisher_p_value = fisher_result[["p_value"]],
    stringsAsFactors = FALSE
  )

  detail <- signature_df %>%
    mutate(
      signature_type = signature_type,
      signature = signature_name,
      in_expression_universe = !is.na(match_identifiers(
        gene,
        universe_df$symbol_key,
        universe_df$ensembl_key
      ))
    ) %>%
    select(signature_type, signature, everything())

  list(
    summary = summary_row,
    detail = detail,
    null = list(
      overlap = null_overlap,
      prevalence_weighted_burden = null_burden,
      overall_lfc_burden = null_overall_lfc_burden
    )
  )
}

add_empirical_classification <- function(summary_df, fdr_threshold) {
  summary_df %>%
    group_by(signature_type) %>%
    mutate(
      expression_matched_overlap_q = p.adjust(expression_matched_overlap_p, method = "BH"),
      expression_matched_burden_q = p.adjust(expression_matched_burden_p, method = "BH"),
      selective_vulnerability_class = case_when(
        expression_matched_overlap_q <= fdr_threshold &
          expression_matched_burden_q <= fdr_threshold ~
          "Strong: overlap and depletion burden enriched",
        expression_matched_burden_q <= fdr_threshold ~
          "Depletion-strength enriched",
        expression_matched_overlap_q <= fdr_threshold ~
          "Overlap enriched",
        expression_matched_overlap_p < 0.05 &
          expression_matched_burden_p < 0.05 ~
          "Nominal overlap and burden; not FDR-significant",
        expression_matched_burden_p < 0.05 ~
          "Nominal burden enrichment; not FDR-significant",
        expression_matched_overlap_p < 0.05 ~
          "Nominal overlap enrichment; not FDR-significant",
        expression_matched_burden_z >= 0 ~
          "Above-background trend; not significant",
        expression_matched_burden_z < 0 ~
          "Below expression-matched background",
        TRUE ~ "Indeterminate"
      )
    ) %>%
    ungroup()
}

####################
# Replot-only cache branch
####################
if (replot_only) {
  if (!file.exists(score_cache_file)) {
    stop("PDO_REPLOT_ONLY=1 but scoring cache is missing: ", score_cache_file)
  }
  message("Replot-only mode: loading existing scoring cache.")
  scoring_cache <- readRDS(score_cache_file)
  reference_df <- scoring_cache$reference
  dge_all <- scoring_cache$dge_all
  dge_signature <- scoring_cache$dge_signature
  summary_all <- scoring_cache$summary
  detail_all <- scoring_cache$detail
  null_results <- scoring_cache$null
} else {
  ####################
  # Oesophageal selective-dependency reference
  ####################
  message("Reading the oesophageal differential-dependency sheet only.")
  reference_df <- read_excel(workbook_file, sheet = params$reference_sheet) %>%
    as.data.frame(stringsAsFactors = FALSE)
  expected_columns <- c(
    "gene", "ENSEMBL_ID", "MEAN_LFC", "MEDIAN_LFC",
    "MEAN_LFC_not_depleted_group", "MEAN_LFC_depleted_group",
    "MEDIAN_LFC_not_depleted_group", "MEDIAN_LFC_depleted_group",
    "delta_MEAN", "delta_MEDIAN", "N_not_depleted", "N_depleted",
    "p.value", "p.value.adj"
  )
  missing_reference_columns <- setdiff(expected_columns, colnames(reference_df))
  if (length(missing_reference_columns) > 0) {
    stop("CRISPR workbook is missing expected column(s): ", paste(missing_reference_columns, collapse = ", "))
  }

  reference_df <- reference_df %>%
    transmute(
      gene = as.character(gene),
      ensembl_id = sub("\\..*$", "", toupper(as.character(ENSEMBL_ID))),
      mean_lfc = as.numeric(MEAN_LFC),
      median_lfc = as.numeric(MEDIAN_LFC),
      mean_lfc_not_depleted_group = as.numeric(MEAN_LFC_not_depleted_group),
      mean_lfc_depleted_group = as.numeric(MEAN_LFC_depleted_group),
      median_lfc_not_depleted_group = as.numeric(MEDIAN_LFC_not_depleted_group),
      median_lfc_depleted_group = as.numeric(MEDIAN_LFC_depleted_group),
      delta_mean = as.numeric(delta_MEAN),
      delta_median = as.numeric(delta_MEDIAN),
      n_not_depleted = as.integer(N_not_depleted),
      n_depleted = as.integer(N_depleted),
      p_value = as.numeric(p.value),
      p_value_adj = as.numeric(p.value.adj),
      depletion_prevalence = n_depleted / (n_depleted + n_not_depleted),
      conditional_depletion_strength = pmax(0, -mean_lfc_depleted_group),
      prevalence_weighted_strength = depletion_prevalence * conditional_depletion_strength,
      overall_depletion_strength = pmax(0, -mean_lfc),
      gene_key = toupper(gene),
      ensembl_key = ensembl_id
    ) %>%
    arrange(p_value_adj, gene)

  if (nrow(reference_df) != 4082L) {
    warning("Expected 4,082 oesophageal differential-dependency genes; found ", nrow(reference_df), ".")
  }
  if (any(reference_df$p_value_adj >= params$reference_fdr_guard, na.rm = TRUE)) {
    stop("The selected oesophageal sheet contains rows with adjusted P >= 0.05.")
  }
  if (anyDuplicated(reference_df$gene_key)) {
    stop("Duplicated gene symbols found in the oesophageal dependency reference.")
  }

  fwrite(
    reference_df %>% select(-gene_key, -ensembl_key),
    file.path(out_paths[["tables"]], "Auto_oesophageal_selective_dependency_reference.csv")
  )

  ####################
  # Canonical six-state cells and expression universe
  ####################
  message("Loading canonical centred states and building the state-restricted count matrix.")
  pdos_all <- readRDS(seurat_file)
  state_vector <- pdo_normalise_final_state_vector(readRDS(state_file))
  common_cells <- intersect(colnames(pdos_all), names(state_vector))
  if (length(common_cells) == 0) stop("No cells overlap between merged PDO and state vector.")

  meta <- pdos_all@meta.data[common_cells, , drop = FALSE] %>%
    rownames_to_column("cell") %>%
    transmute(
      cell,
      orig.ident = as.character(orig.ident),
      state = as.character(state_vector[cell])
    ) %>%
    filter(
      state %in% PDO_STATE_ORDER,
      orig.ident != PDO_EXCLUDED_SAMPLE
    )
  keep_cells <- meta$cell
  if (length(keep_cells) == 0) stop("No canonical six-state cells remain after filtering.")

  counts_mat <- pdo_get_assay_matrix(pdos_all, assay = "RNA", layer = "counts")
  counts_mat <- counts_mat[, keep_cells, drop = FALSE]
  keep_features <- Matrix::rowSums(counts_mat > 0) >= params$min_detected_cells
  counts_mat <- counts_mat[keep_features, , drop = FALSE]
  rm(pdos_all, state_vector)
  invisible(gc())

  total_gene_counts <- Matrix::rowSums(counts_mat)
  total_library_counts <- sum(total_gene_counts)
  universe_keys <- identifier_keys(rownames(counts_mat))
  universe_df <- data.frame(
    feature = rownames(counts_mat),
    symbol_key = universe_keys$symbol_key,
    ensembl_key = universe_keys$ensembl_key,
    detected_cell_n = as.integer(Matrix::rowSums(counts_mat > 0)),
    log_cpm = log2((total_gene_counts + 0.5) / (total_library_counts + 1) * 1e6),
    stringsAsFactors = FALSE
  ) %>%
    mutate(expression_bin = ntile(log_cpm, params$expression_bins))

  universe_reference_idx <- match_identifiers(
    universe_df$feature,
    reference_df$gene_key,
    reference_df$ensembl_key
  )
  universe_df$is_selective_dependency <- !is.na(universe_reference_idx)
  universe_df$depletion_prevalence <- ifelse(
    universe_df$is_selective_dependency,
    reference_df$depletion_prevalence[universe_reference_idx],
    0
  )
  universe_df$conditional_depletion_strength <- ifelse(
    universe_df$is_selective_dependency,
    reference_df$conditional_depletion_strength[universe_reference_idx],
    0
  )
  universe_df$prevalence_weighted_strength <- ifelse(
    universe_df$is_selective_dependency,
    reference_df$prevalence_weighted_strength[universe_reference_idx],
    0
  )
  universe_df$overall_depletion_strength <- ifelse(
    universe_df$is_selective_dependency,
    reference_df$overall_depletion_strength[universe_reference_idx],
    0
  )

  fwrite(
    universe_df,
    file.path(out_paths[["tables"]], "Auto_pdo_expression_dependency_universe.csv.gz")
  )

  ####################
  # Fresh sample-blocked edgeR state-vs-rest DGE
  ####################
  if (file.exists(dge_cache_file) && !force_rebuild) {
    message("Reusing cached six-state edgeR DGE.")
    dge_all <- readRDS(dge_cache_file)
  } else {
    message("Running fresh sample-blocked edgeR QL state-vs-rest DGE.")
    dge_results <- lapply(PDO_STATE_ORDER, function(state_name) {
      message("  edgeR: ", state_name)
      state_result <- run_pseudobulk_state_dge(
        counts_mat = counts_mat,
        meta = meta,
        state_name = state_name,
        all_states = PDO_STATE_ORDER,
        params = params
      )
      fwrite(
        state_result,
        file.path(
          out_paths[["tables"]],
          paste0("Auto_state_edger_", safe_name(state_name), ".csv.gz")
        )
      )
      state_result
    })
    dge_all <- bind_rows(dge_results)
    saveRDS(dge_all, dge_cache_file)
  }
  if (nrow(dge_all) == 0) stop("The six-state edgeR workflow returned no genes.")
  if (!all(dge_all$dge_method == "pseudobulk_edgeR_QLF_sample_blocked")) {
    stop("Unexpected DGE method in the cached or newly generated table.")
  }
  fwrite(
    dge_all,
    file.path(out_paths[["tables"]], "Auto_state_edger_dge_all.csv.gz")
  )

  dge_signature <- bind_rows(lapply(PDO_STATE_ORDER, function(state_name) {
    dge_all %>%
      filter(state == state_name, avg_log2FC > 0) %>%
      arrange(p_val_adj, desc(avg_log2FC), desc(pct_target), gene) %>%
      slice_head(n = params$dge_top_n) %>%
      mutate(
        signature_rank = row_number(),
        definition_weight = avg_log2FC,
        component_mps = NA_character_
      )
  }))
  dge_signature_counts <- dge_signature %>% count(state, name = "signature_gene_n")
  if (any(dge_signature_counts$signature_gene_n < params$dge_top_n)) {
    warning("At least one state has fewer than ", params$dge_top_n, " positive edgeR genes.")
  }
  fwrite(
    dge_signature,
    file.path(out_paths[["tables"]], "Auto_state_edger_signature_top150.csv")
  )

  ####################
  # MP and state signature construction
  ####################
  mp_genes <- readRDS(mp_gene_file)
  mp_weights <- readRDS(mp_weight_file)
  expected_mps <- c(PDO_CELL_CYCLE_MPS, unlist(PDO_MP_STATE_GROUPS, use.names = FALSE))
  missing_mps <- setdiff(expected_mps, names(mp_genes))
  if (length(missing_mps) > 0) stop("Canonical MP genes are missing MP(s): ", paste(missing_mps, collapse = ", "))

  mp_signatures <- lapply(expected_mps, function(mp_name) {
    genes <- unique(as.character(mp_genes[[mp_name]]))
    weights <- as.numeric(mp_weights[[mp_name]][genes])
    weights[!is.finite(weights) | weights <= 0] <- 1
    data.frame(
      gene = genes,
      definition_weight = weights,
      component_mps = mp_name,
      stringsAsFactors = FALSE
    )
  })
  names(mp_signatures) <- expected_mps

  state_mp_signatures <- lapply(names(PDO_MP_STATE_GROUPS), function(state_name) {
    component_mps <- PDO_MP_STATE_GROUPS[[state_name]]
    membership <- bind_rows(lapply(component_mps, function(mp_name) {
      data.frame(gene = unique(as.character(mp_genes[[mp_name]])), component_mp = mp_name)
    }))
    membership %>%
      group_by(gene) %>%
      summarise(
        definition_weight = 1,
        component_mps = paste(component_mp, collapse = ";"),
        .groups = "drop"
      )
  })
  names(state_mp_signatures) <- names(PDO_MP_STATE_GROUPS)

  state_dge_signatures <- lapply(PDO_STATE_ORDER, function(state_name) {
    dge_signature %>%
      filter(state == state_name) %>%
      transmute(
        gene,
        definition_weight,
        component_mps = NA_character_
      )
  })
  names(state_dge_signatures) <- PDO_STATE_ORDER

  ####################
  # Expression-matched overlap and depletion-strength scoring
  ####################
  message("Scoring MP and state signatures with expression-matched permutations.")
  scoring_results <- list()
  for (mp_name in names(mp_signatures)) {
    scoring_results[[paste0("mp__", mp_name)]] <- score_signature(
      signature_name = mp_name,
      signature_type = "MP gene list",
      signature_df = mp_signatures[[mp_name]],
      universe_df = universe_df,
      reference_df = reference_df,
      n_perm = params$matched_permutations
    )
  }
  for (state_name in names(state_mp_signatures)) {
    scoring_results[[paste0("state_mp__", state_name)]] <- score_signature(
      signature_name = state_name,
      signature_type = "State combined MP genes",
      signature_df = state_mp_signatures[[state_name]],
      universe_df = universe_df,
      reference_df = reference_df,
      n_perm = params$matched_permutations
    )
  }
  for (state_name in names(state_dge_signatures)) {
    scoring_results[[paste0("state_dge__", state_name)]] <- score_signature(
      signature_name = state_name,
      signature_type = "State edgeR top-150 genes",
      signature_df = state_dge_signatures[[state_name]],
      universe_df = universe_df,
      reference_df = reference_df,
      n_perm = params$matched_permutations
    )
  }

  summary_all <- bind_rows(lapply(scoring_results, `[[`, "summary")) %>%
    add_empirical_classification(params$empirical_fdr)
  detail_all <- bind_rows(lapply(scoring_results, `[[`, "detail"))
  null_results <- lapply(scoring_results, `[[`, "null")

  scoring_cache <- list(
    reference = reference_df,
    dge_all = dge_all,
    dge_signature = dge_signature,
    universe = universe_df,
    summary = summary_all,
    detail = detail_all,
    null = null_results,
    parameters = params
  )
  saveRDS(scoring_cache, score_cache_file)
  rm(counts_mat)
  invisible(gc())
}

####################
# Independent scATLAS MP/state selective-dependency analysis
####################
if (file.exists(scatlas_score_cache_file) && !force_rebuild) {
  message("Reusing scATLAS selective-dependency scoring cache.")
  scatlas_scoring_cache <- readRDS(scatlas_score_cache_file)
  scatlas_summary_all <- scatlas_scoring_cache$summary
  scatlas_detail_all <- scatlas_scoring_cache$detail
  scatlas_null_results <- scatlas_scoring_cache$null
  scatlas_marker_signature <- scatlas_scoring_cache$marker_signature
  scatlas_universe_df <- scatlas_scoring_cache$universe
} else {
  if (replot_only) {
    stop(
      "PDO_REPLOT_ONLY=1 but the scATLAS scoring cache is missing: ",
      scatlas_score_cache_file
    )
  }
  message("Building the independent scATLAS expression universe and signatures.")
  scatlas_obj <- readRDS(scatlas_seurat_file)
  scatlas_states <- readRDS(scatlas_state_file)
  if (is.null(names(scatlas_states))) {
    stop("The current scATLAS state vector is not named by cell barcode.")
  }
  scatlas_common_cells <- intersect(colnames(scatlas_obj), names(scatlas_states))
  scatlas_keep_cells <- scatlas_common_cells[
    as.character(scatlas_states[scatlas_common_cells]) %in% SCREF_PRIMARY_STATE_ORDER
  ]
  if (length(scatlas_keep_cells) == 0) {
    stop("No current canonical scATLAS state cells overlap EAC_Ref_epi.rds.")
  }
  scatlas_counts <- pdo_get_assay_matrix(scatlas_obj, assay = "RNA", layer = "counts")
  scatlas_counts <- scatlas_counts[, scatlas_keep_cells, drop = FALSE]
  scatlas_keep_features <- Matrix::rowSums(scatlas_counts > 0) >= params$min_detected_cells
  scatlas_counts <- scatlas_counts[scatlas_keep_features, , drop = FALSE]
  rm(scatlas_obj, scatlas_states)
  invisible(gc())

  scatlas_total_gene_counts <- Matrix::rowSums(scatlas_counts)
  scatlas_total_library_counts <- sum(scatlas_total_gene_counts)
  scatlas_universe_keys <- identifier_keys(rownames(scatlas_counts))
  scatlas_universe_df <- data.frame(
    feature = rownames(scatlas_counts),
    symbol_key = scatlas_universe_keys$symbol_key,
    ensembl_key = scatlas_universe_keys$ensembl_key,
    detected_cell_n = as.integer(Matrix::rowSums(scatlas_counts > 0)),
    log_cpm = log2(
      (scatlas_total_gene_counts + 0.5) /
        (scatlas_total_library_counts + 1) * 1e6
    ),
    stringsAsFactors = FALSE
  ) %>%
    mutate(expression_bin = ntile(log_cpm, params$expression_bins))

  scatlas_universe_reference_idx <- match_identifiers(
    scatlas_universe_df$feature,
    reference_df$gene_key,
    reference_df$ensembl_key
  )
  scatlas_universe_df$is_selective_dependency <-
    !is.na(scatlas_universe_reference_idx)
  scatlas_universe_df$depletion_prevalence <- ifelse(
    scatlas_universe_df$is_selective_dependency,
    reference_df$depletion_prevalence[scatlas_universe_reference_idx],
    0
  )
  scatlas_universe_df$conditional_depletion_strength <- ifelse(
    scatlas_universe_df$is_selective_dependency,
    reference_df$conditional_depletion_strength[scatlas_universe_reference_idx],
    0
  )
  scatlas_universe_df$prevalence_weighted_strength <- ifelse(
    scatlas_universe_df$is_selective_dependency,
    reference_df$prevalence_weighted_strength[scatlas_universe_reference_idx],
    0
  )
  scatlas_universe_df$overall_depletion_strength <- ifelse(
    scatlas_universe_df$is_selective_dependency,
    reference_df$overall_depletion_strength[scatlas_universe_reference_idx],
    0
  )
  fwrite(
    scatlas_universe_df,
    file.path(out_paths[["tables"]], "Auto_scatlas_expression_dependency_universe.csv.gz")
  )
  rm(scatlas_counts)
  invisible(gc())

  scatlas_mp_genes_all <- readRDS(scatlas_mp_gene_file)
  scatlas_mp_weights_all <- readRDS(scatlas_mp_weight_file)
  scatlas_expected_mps <- c(
    SCREF_CC_MPS,
    unlist(SCREF_STATE_GROUPS, use.names = FALSE)
  )
  scatlas_missing_mps <- setdiff(scatlas_expected_mps, names(scatlas_mp_genes_all))
  if (length(scatlas_missing_mps) > 0) {
    stop("Current scATLAS MP genes are missing: ", paste(scatlas_missing_mps, collapse = ", "))
  }
  scatlas_missing_weights <- setdiff(scatlas_expected_mps, names(scatlas_mp_weights_all))
  if (length(scatlas_missing_weights) > 0) {
    stop("Current scATLAS MP weights are missing: ", paste(scatlas_missing_weights, collapse = ", "))
  }

  scatlas_mp_signatures <- lapply(scatlas_expected_mps, function(mp_name) {
    genes <- unique(as.character(scatlas_mp_genes_all[[mp_name]]))
    weights <- as.numeric(scatlas_mp_weights_all[[mp_name]][genes])
    weights[!is.finite(weights) | weights <= 0] <- 1
    data.frame(
      gene = genes,
      definition_weight = weights,
      component_mps = mp_name,
      stringsAsFactors = FALSE
    )
  })
  names(scatlas_mp_signatures) <- scatlas_expected_mps

  scatlas_state_mp_signatures <- lapply(names(SCREF_STATE_GROUPS), function(state_name) {
    component_mps <- SCREF_STATE_GROUPS[[state_name]]
    bind_rows(lapply(component_mps, function(mp_name) {
      data.frame(
        gene = unique(as.character(scatlas_mp_genes_all[[mp_name]])),
        component_mp = mp_name,
        stringsAsFactors = FALSE
      )
    })) %>%
      group_by(gene) %>%
      summarise(
        definition_weight = 1,
        component_mps = paste(component_mp, collapse = ";"),
        .groups = "drop"
      )
  })
  names(scatlas_state_mp_signatures) <- names(SCREF_STATE_GROUPS)

  scatlas_ranked_markers_all <- fread(scatlas_state_marker_file)
  required_scatlas_columns <- c(
    "state", "gene", "ranking_score", "reproducibility_score",
    "median_log2FC_hit", "specificity_gap"
  )
  missing_scatlas_columns <- setdiff(
    required_scatlas_columns,
    colnames(scatlas_ranked_markers_all)
  )
  if (length(missing_scatlas_columns) > 0) {
    stop(
      "The current scATLAS marker table is missing: ",
      paste(missing_scatlas_columns, collapse = ", ")
    )
  }
  scatlas_marker_signature <- bind_rows(lapply(SCREF_PRIMARY_STATE_ORDER, function(state_name) {
    scatlas_ranked_markers_all %>%
      filter(state == state_name) %>%
      arrange(
        desc(ranking_score),
        desc(reproducibility_score),
        desc(median_log2FC_hit),
        desc(specificity_gap),
        gene
      ) %>%
      distinct(gene, .keep_all = TRUE) %>%
      slice_head(n = params$scatlas_dge_top_n) %>%
      mutate(
        signature_rank = row_number(),
        definition_weight = median_log2FC_hit,
        component_mps = NA_character_
      )
  }))
  scatlas_marker_counts <- scatlas_marker_signature %>%
    count(state, name = "signature_gene_n")
  if (
    nrow(scatlas_marker_counts) != length(SCREF_PRIMARY_STATE_ORDER) ||
      any(scatlas_marker_counts$signature_gene_n != params$scatlas_dge_top_n)
  ) {
    stop("At least one current scATLAS state lacks 150 ranked positive markers.")
  }
  fwrite(
    scatlas_marker_signature,
    file.path(out_paths[["tables"]], "Auto_scatlas_state_marker_signature_top150.csv")
  )

  scatlas_state_marker_signatures <- lapply(SCREF_PRIMARY_STATE_ORDER, function(state_name) {
    scatlas_marker_signature %>%
      filter(state == state_name) %>%
      transmute(gene, definition_weight, component_mps)
  })
  names(scatlas_state_marker_signatures) <- SCREF_PRIMARY_STATE_ORDER

  message("Scoring all current scATLAS MP and state signatures.")
  set.seed(10831)
  scatlas_scoring_results <- list()
  for (mp_name in names(scatlas_mp_signatures)) {
    scatlas_scoring_results[[paste0("mp__", mp_name)]] <- score_signature(
      signature_name = mp_name,
      signature_type = "scATLAS MP gene list",
      signature_df = scatlas_mp_signatures[[mp_name]],
      universe_df = scatlas_universe_df,
      reference_df = reference_df,
      n_perm = params$matched_permutations
    )
  }
  for (state_name in names(scatlas_state_mp_signatures)) {
    scatlas_scoring_results[[paste0("state_mp__", state_name)]] <- score_signature(
      signature_name = state_name,
      signature_type = "scATLAS state combined MP genes",
      signature_df = scatlas_state_mp_signatures[[state_name]],
      universe_df = scatlas_universe_df,
      reference_df = reference_df,
      n_perm = params$matched_permutations
    )
  }
  for (state_name in names(scatlas_state_marker_signatures)) {
    scatlas_scoring_results[[paste0("state_marker__", state_name)]] <- score_signature(
      signature_name = state_name,
      signature_type = "scATLAS state marker top-150 genes",
      signature_df = scatlas_state_marker_signatures[[state_name]],
      universe_df = scatlas_universe_df,
      reference_df = reference_df,
      n_perm = params$matched_permutations
    )
  }

  scatlas_summary_all <- bind_rows(lapply(scatlas_scoring_results, `[[`, "summary")) %>%
    add_empirical_classification(params$empirical_fdr)
  scatlas_detail_all <- bind_rows(lapply(scatlas_scoring_results, `[[`, "detail"))
  scatlas_null_results <- lapply(scatlas_scoring_results, `[[`, "null")
  scatlas_scoring_cache <- list(
    summary = scatlas_summary_all,
    detail = scatlas_detail_all,
    null = scatlas_null_results,
    marker_signature = scatlas_marker_signature,
    universe = scatlas_universe_df,
    parameters = params
  )
  saveRDS(scatlas_scoring_cache, scatlas_score_cache_file)
}

####################
# Classic proliferation: scATLAS / CRISPR / PDO three-way overlaps
####################
message("Building Classic-proliferation cross-dataset three-way overlaps.")
if (!identical(unname(PDO_MP_STATE_GROUPS[["Classic proliferation"]]), "MP19+")) {
  stop("PDO Classic proliferation is no longer defined by MP19+ alone; update the Venn specification.")
}

scatlas_mp_genes <- readRDS(scatlas_mp_gene_file)
pdo_mp_genes <- readRDS(mp_gene_file)
if (!"MP2+" %in% names(scatlas_mp_genes)) {
  stop("The current scATLAS MP object does not contain Classic-proliferation MP2+.")
}
if (!"MP19+" %in% names(pdo_mp_genes)) {
  stop("The current PDO MP object does not contain Classic-proliferation MP19+.")
}

scatlas_ranked_markers <- fread(scatlas_state_marker_file)
required_scatlas_marker_columns <- c(
  "state", "gene", "ranking_score", "reproducibility_score",
  "median_log2FC_hit", "specificity_gap"
)
missing_scatlas_marker_columns <- setdiff(
  required_scatlas_marker_columns,
  colnames(scatlas_ranked_markers)
)
if (length(missing_scatlas_marker_columns) > 0) {
  stop(
    "The current scATLAS ranked-marker table is missing column(s): ",
    paste(missing_scatlas_marker_columns, collapse = ", ")
  )
}
scatlas_classic_dge <- scatlas_ranked_markers %>%
  filter(state == "Classic proliferation") %>%
  arrange(
    desc(ranking_score),
    desc(reproducibility_score),
    desc(median_log2FC_hit),
    desc(specificity_gap),
    gene
  ) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_head(n = params$scatlas_dge_top_n) %>%
  mutate(signature_rank = row_number())
if (nrow(scatlas_classic_dge) != params$scatlas_dge_top_n) {
  stop(
    "Expected ", params$scatlas_dge_top_n,
    " current scATLAS Classic-proliferation DGE genes; found ",
    nrow(scatlas_classic_dge), "."
  )
}

pdo_classic_dge <- dge_signature %>%
  filter(state == "Classic proliferation") %>%
  arrange(signature_rank) %>%
  distinct(gene, .keep_all = TRUE)
if (nrow(pdo_classic_dge) != params$dge_top_n) {
  stop(
    "Expected ", params$dge_top_n,
    " current PDO Classic-proliferation edgeR genes; found ",
    nrow(pdo_classic_dge), "."
  )
}

classic_signature_sets <- list(
  scatlas_mp = canonicalize_gene_set(scatlas_mp_genes[["MP2+"]], reference_df),
  scatlas_dge = canonicalize_gene_set(scatlas_classic_dge$gene, reference_df),
  crispr = sort(unique(reference_df$gene_key)),
  pdo_mp = canonicalize_gene_set(pdo_mp_genes[["MP19+"]], reference_df),
  pdo_dge = canonicalize_gene_set(pdo_classic_dge$gene, reference_df)
)
if (any(lengths(classic_signature_sets) == 0)) {
  stop("At least one Classic-proliferation Venn set is empty.")
}

classic_panel_spec <- tibble(
  panel = c(
    "scATLAS MP2+ | PDO MP19+",
    "scATLAS MP2+ | PDO Classic DGE top 150",
    "scATLAS Classic DGE top 150 | PDO MP19+",
    "scATLAS Classic DGE top 150 | PDO Classic DGE top 150"
  ),
  scatlas_key = c("scatlas_mp", "scatlas_mp", "scatlas_dge", "scatlas_dge"),
  scatlas_label = c("MP2+", "MP2+", "Classic DGE top 150", "Classic DGE top 150"),
  pdo_key = c("pdo_mp", "pdo_dge", "pdo_mp", "pdo_dge"),
  pdo_label = c("MP19+", "Classic DGE top 150", "MP19+", "Classic DGE top 150")
)

classic_venn_membership <- bind_rows(lapply(seq_len(nrow(classic_panel_spec)), function(i) {
  build_threeway_membership(
    panel_name = classic_panel_spec$panel[i],
    scatlas_genes = classic_signature_sets[[classic_panel_spec$scatlas_key[i]]],
    crispr_genes = classic_signature_sets$crispr,
    pdo_genes = classic_signature_sets[[classic_panel_spec$pdo_key[i]]]
  )
}))

classic_region_order <- c(
  "scATLAS only", "CRISPR only", "PDO only",
  "scATLAS & CRISPR only", "scATLAS & PDO only", "CRISPR & PDO only",
  "scATLAS & CRISPR & PDO"
)
classic_venn_counts <- classic_venn_membership %>%
  count(panel, region, name = "gene_n") %>%
  right_join(
    tidyr::crossing(
      panel = classic_panel_spec$panel,
      region = classic_region_order
    ),
    by = c("panel", "region")
  ) %>%
  mutate(gene_n = replace_na(gene_n, 0L)) %>%
  left_join(classic_panel_spec, by = "panel") %>%
  mutate(
    scatlas_set_n = lengths(classic_signature_sets[scatlas_key]),
    crispr_set_n = length(classic_signature_sets$crispr),
    pdo_set_n = lengths(classic_signature_sets[pdo_key]),
    panel = factor(panel, levels = classic_panel_spec$panel),
    region = factor(region, levels = classic_region_order)
  ) %>%
  arrange(panel, region) %>%
  mutate(panel = as.character(panel), region = as.character(region))

classic_venn_membership <- classic_venn_membership %>%
  mutate(panel = factor(panel, levels = classic_panel_spec$panel)) %>%
  arrange(panel, gene) %>%
  mutate(panel = as.character(panel))

classic_threeway_genes <- classic_venn_membership %>%
  filter(in_scatlas, in_crispr, in_pdo) %>%
  select(panel, gene) %>%
  arrange(match(panel, classic_panel_spec$panel), gene)

classic_venn_counts_file <- file.path(
  out_paths[["tables"]],
  "Auto_classic_proliferation_threeway_venn_counts.csv"
)
classic_venn_membership_file <- file.path(
  out_paths[["tables"]],
  "Auto_classic_proliferation_threeway_venn_gene_membership.csv.gz"
)
classic_threeway_genes_file <- file.path(
  out_paths[["tables"]],
  "Auto_classic_proliferation_threeway_all_overlap_genes.csv"
)
classic_venn_sets_file <- file.path(
  out_paths[["intermediate"]],
  "Auto_classic_proliferation_threeway_venn_sets.rds"
)
fwrite(classic_venn_counts, classic_venn_counts_file)
fwrite(classic_venn_membership, classic_venn_membership_file)
fwrite(classic_threeway_genes, classic_threeway_genes_file)
saveRDS(
  list(
    panel_specification = classic_panel_spec,
    gene_sets = classic_signature_sets,
    scatlas_classic_dge_top150 = scatlas_classic_dge,
    pdo_classic_dge_top150 = pdo_classic_dge,
    counts = classic_venn_counts
  ),
  classic_venn_sets_file
)

####################
# Persistent result tables and workbook
####################
summary_all <- summary_all %>%
  select(
    -any_of(c(
      "expression_matched_overlap_q",
      "expression_matched_burden_q",
      "selective_vulnerability_class"
    ))
  ) %>%
  add_empirical_classification(params$empirical_fdr)

if (replot_only) {
  scoring_cache$summary <- summary_all
  saveRDS(scoring_cache, score_cache_file)
}

mp_summary <- summary_all %>%
  filter(signature_type == "MP gene list") %>%
  mutate(
    state_group = case_when(
      signature %in% PDO_CELL_CYCLE_MPS ~ "Cell cycle/QC only",
      TRUE ~ vapply(signature, function(mp_name) {
        groups <- names(PDO_MP_STATE_GROUPS)[vapply(PDO_MP_STATE_GROUPS, function(x) mp_name %in% x, logical(1))]
        ifelse(length(groups) == 0, "Other", groups[1])
      }, character(1))
    ),
    description = unname(PDO_MP_DESCRIPTIONS[signature])
  ) %>%
  arrange(match(signature, c(PDO_CELL_CYCLE_MPS, unlist(PDO_MP_STATE_GROUPS, use.names = FALSE))))

state_summary <- summary_all %>%
  filter(signature_type != "MP gene list") %>%
  arrange(match(signature, PDO_STATE_ORDER), signature_type)

mp_detail <- detail_all %>% filter(signature_type == "MP gene list")
state_mp_detail <- detail_all %>% filter(signature_type == "State combined MP genes")
state_dge_detail <- detail_all %>% filter(signature_type == "State edgeR top-150 genes")

scatlas_summary_all <- scatlas_summary_all %>%
  select(
    -any_of(c(
      "expression_matched_overlap_q",
      "expression_matched_burden_q",
      "selective_vulnerability_class"
    ))
  ) %>%
  add_empirical_classification(params$empirical_fdr)
scatlas_mp_summary <- scatlas_summary_all %>%
  filter(signature_type == "scATLAS MP gene list") %>%
  mutate(
    state_group = case_when(
      signature %in% SCREF_CC_MPS ~ "Cell cycle/QC only",
      TRUE ~ vapply(signature, function(mp_name) {
        groups <- names(SCREF_STATE_GROUPS)[vapply(
          SCREF_STATE_GROUPS,
          function(x) mp_name %in% x,
          logical(1)
        )]
        ifelse(length(groups) == 0, "Other", groups[1])
      }, character(1))
    ),
    description = unname(SCREF_MP_DESCRIPTIONS[signature])
  ) %>%
  arrange(match(signature, c(SCREF_CC_MPS, unlist(SCREF_STATE_GROUPS, use.names = FALSE))))
scatlas_state_summary <- scatlas_summary_all %>%
  filter(signature_type != "scATLAS MP gene list") %>%
  arrange(match(signature, SCREF_PRIMARY_STATE_ORDER), signature_type)
scatlas_mp_detail <- scatlas_detail_all %>%
  filter(signature_type == "scATLAS MP gene list")
scatlas_state_mp_detail <- scatlas_detail_all %>%
  filter(signature_type == "scATLAS state combined MP genes")
scatlas_state_marker_detail <- scatlas_detail_all %>%
  filter(signature_type == "scATLAS state marker top-150 genes")

fwrite(mp_summary, file.path(out_paths[["tables"]], "Auto_mp_selective_dependency_summary.csv"))
fwrite(state_summary, file.path(out_paths[["tables"]], "Auto_state_selective_dependency_summary.csv"))
fwrite(mp_detail, file.path(out_paths[["tables"]], "Auto_mp_selective_dependency_gene_details.csv"))
fwrite(state_mp_detail, file.path(out_paths[["tables"]], "Auto_state_combined_mp_dependency_gene_details.csv"))
fwrite(state_dge_detail, file.path(out_paths[["tables"]], "Auto_state_edger_dependency_gene_details.csv"))
fwrite(
  scatlas_mp_summary,
  file.path(out_paths[["tables"]], "Auto_scatlas_mp_selective_dependency_summary.csv")
)
fwrite(
  scatlas_state_summary,
  file.path(out_paths[["tables"]], "Auto_scatlas_state_selective_dependency_summary.csv")
)
fwrite(
  scatlas_mp_detail,
  file.path(out_paths[["tables"]], "Auto_scatlas_mp_dependency_gene_details.csv")
)
fwrite(
  scatlas_state_mp_detail,
  file.path(out_paths[["tables"]], "Auto_scatlas_state_combined_mp_dependency_gene_details.csv")
)
fwrite(
  scatlas_state_marker_detail,
  file.path(out_paths[["tables"]], "Auto_scatlas_state_marker_dependency_gene_details.csv")
)

workbook_out <- file.path(
  out_paths[["tables"]],
  "Auto_crispr_selective_dependency_results.xlsx"
)
write.xlsx(
  list(
    MP_summary = mp_summary,
    state_summary = state_summary,
    MP_gene_details = mp_detail,
    state_MP_gene_details = state_mp_detail,
    state_edgeR_gene_details = state_dge_detail,
    state_edgeR_top150 = dge_signature,
    scATLAS_MP_summary = scatlas_mp_summary,
    scATLAS_state_summary = scatlas_state_summary,
    scATLAS_MP_gene_details = scatlas_mp_detail,
    scATLAS_state_MP_details = scatlas_state_mp_detail,
    scATLAS_marker_details = scatlas_state_marker_detail,
    scATLAS_marker_top150 = scatlas_marker_signature,
    classic_venn_counts = classic_venn_counts,
    classic_venn_membership = classic_venn_membership,
    classic_threeway_genes = classic_threeway_genes,
    oesophageal_reference = reference_df %>% select(-gene_key, -ensembl_key)
  ),
  workbook_out,
  overwrite = TRUE,
  asTable = TRUE
)

####################
# Presentation figures
####################
class_colors <- c(
  "Strong: overlap and depletion burden enriched" = "#B2182B",
  "Depletion-strength enriched" = "#EF8A62",
  "Overlap enriched" = "#FDDDBC",
  "Nominal overlap and burden; not FDR-significant" = "#B35806",
  "Nominal burden enrichment; not FDR-significant" = "#F1A340",
  "Nominal overlap enrichment; not FDR-significant" = "#998EC3",
  "Above-background trend; not significant" = "#7F7F7F",
  "Below expression-matched background" = "#2166AC",
  "Indeterminate" = "#CCCCCC"
)

####################
# Independent scATLAS summary figures
####################
scatlas_mp_plot_df <- scatlas_mp_summary %>%
  mutate(
    plot_label = paste0(signature, " - ", description),
    plot_label = factor(plot_label, levels = rev(plot_label))
  )
scatlas_mp_plot <- ggplot(
  scatlas_mp_plot_df,
  aes(
    x = dependency_index_0_100,
    y = plot_label,
    size = overlap_fraction,
    color = selective_vulnerability_class
  )
) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = class_colors, drop = TRUE) +
  scale_size_continuous(labels = percent_format(accuracy = 1), range = c(3, 10)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Selective oesophageal dependency of current scATLAS MPs",
    subtitle = "Point size = raw overlap; x-axis = scATLAS expression-matched depletion burden",
    x = "Selective-dependency index (0-100)",
    y = NULL,
    size = "Gene overlap",
    color = "Evidence class"
  ) +
  pdo_theme_slide(12) +
  theme(legend.position = "right")
scatlas_mp_figure_file <- file.path(
  out_paths[["figures"]],
  "Auto_scatlas_mp_selective_dependency_summary.pdf"
)
pdo_save_slide_pdf(
  scatlas_mp_plot,
  scatlas_mp_figure_file,
  width = 14,
  height = 9
)

scatlas_state_plot_df <- scatlas_state_summary %>%
  mutate(
    signature = factor(signature, levels = rev(SCREF_PRIMARY_STATE_ORDER)),
    signature_type = factor(
      signature_type,
      levels = c(
        "scATLAS state combined MP genes",
        "scATLAS state marker top-150 genes"
      ),
      labels = c("Combined component-MP genes", "Ranked top-150 state markers")
    )
  )
scatlas_state_plot <- ggplot(
  scatlas_state_plot_df,
  aes(
    x = dependency_index_0_100,
    y = signature,
    shape = signature_type,
    color = selective_vulnerability_class,
    size = overlap_fraction
  )
) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.9, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = class_colors, drop = TRUE) +
  scale_size_continuous(labels = percent_format(accuracy = 1), range = c(3, 10)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Selective oesophageal dependency of current scATLAS states",
    subtitle = "Independent scATLAS expression-matched background",
    x = "Selective-dependency index (0-100)",
    y = NULL,
    size = "Gene overlap",
    shape = "State definition",
    color = "Evidence class"
  ) +
  pdo_theme_slide(12) +
  theme(legend.position = "right")
scatlas_state_figure_file <- file.path(
  out_paths[["figures"]],
  "Auto_scatlas_state_selective_dependency_summary.pdf"
)
pdo_save_slide_pdf(
  scatlas_state_plot,
  scatlas_state_figure_file,
  width = 14,
  height = 8
)

mp_plot_df <- mp_summary %>%
  mutate(
    plot_label = paste0(signature, " - ", description),
    plot_label = factor(plot_label, levels = rev(plot_label))
  )
mp_plot <- ggplot(
  mp_plot_df,
  aes(
    x = dependency_index_0_100,
    y = plot_label,
    size = overlap_fraction,
    color = selective_vulnerability_class
  )
) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = class_colors, drop = TRUE) +
  scale_size_continuous(labels = percent_format(accuracy = 1), range = c(3, 10)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Selective oesophageal dependency of canonical PDO MPs",
    subtitle = "Point size = raw overlap; x-axis = expression-matched prevalence x depletion-strength burden",
    x = "Selective-dependency index (0-100)",
    y = NULL,
    size = "Gene overlap",
    color = "Evidence class"
  ) +
  pdo_theme_slide(12) +
  theme(legend.position = "right")
pdo_save_slide_pdf(
  mp_plot,
  file.path(out_paths[["figures"]], "Auto_mp_selective_dependency_summary.pdf"),
  width = 14,
  height = 9
)

state_plot_df <- state_summary %>%
  mutate(
    signature = factor(signature, levels = rev(PDO_STATE_ORDER)),
    signature_type = factor(
      signature_type,
      levels = c("State combined MP genes", "State edgeR top-150 genes"),
      labels = c("Combined component-MP genes", "edgeR top-150 state genes")
    )
  )
state_plot <- ggplot(
  state_plot_df,
  aes(
    x = dependency_index_0_100,
    y = signature,
    shape = signature_type,
    color = selective_vulnerability_class,
    size = overlap_fraction
  )
) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.9, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = class_colors, drop = TRUE) +
  scale_size_continuous(labels = percent_format(accuracy = 1), range = c(3, 10)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = "Selective oesophageal dependency of canonical PDO states",
    subtitle = "Point size = raw overlap; x-axis = expression-matched burden; colour = enrichment evidence",
    x = "Selective-dependency index (0-100)",
    y = NULL,
    size = "Gene overlap",
    shape = "State definition",
    color = "Evidence class"
  ) +
  pdo_theme_slide(12) +
  theme(legend.position = "right")
pdo_save_slide_pdf(
  state_plot,
  file.path(out_paths[["figures"]], "Auto_state_selective_dependency_summary.pdf"),
  width = 14,
  height = 8
)

metric_plot_df <- summary_all %>%
  mutate(
    mp_description = unname(PDO_MP_DESCRIPTIONS[signature]),
    label = case_when(
      signature_type == "MP gene list" ~ paste0(signature, " - ", mp_description),
      signature_type == "State combined MP genes" ~ paste0(signature, " - combined component-MP genes"),
      signature_type == "State edgeR top-150 genes" ~ paste0(signature, " - edgeR top-150 state genes"),
      TRUE ~ paste(signature, signature_type, sep = " - ")
    ),
    label = reorder(label, overlap_fraction)
  )
overlap_plot <- ggplot(metric_plot_df, aes(x = overlap_fraction, y = label, color = signature_type)) +
  geom_point(size = 3) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Direct overlap with significant oesophageal differential dependencies",
    x = "Fraction of definition genes overlapping the filtered CRISPR list",
    y = NULL,
    color = "Signature type"
  ) +
  pdo_theme_slide(11) +
  theme(legend.position = "bottom")
pdo_save_slide_pdf(
  overlap_plot,
  file.path(out_paths[["figures"]], "Auto_all_signature_crispr_overlap_fraction.pdf"),
  width = 14,
  height = 10
)

####################
# Classic-proliferation three-way Venn figure
####################
venn_centres <- tibble(
  set = c("scATLAS", "CRISPR", "PDO"),
  centre_x = c(-0.58, 0.58, 0),
  centre_y = c(0.30, 0.30, -0.55),
  colour = c("#377EB8", "#4D4D4D", "#E41A1C")
)
venn_theta <- seq(0, 2 * pi, length.out = 361)
venn_circle_df <- bind_rows(lapply(seq_len(nrow(classic_panel_spec)), function(i) {
  bind_rows(lapply(seq_len(nrow(venn_centres)), function(j) {
    tibble(
      panel = classic_panel_spec$panel[i],
      set = venn_centres$set[j],
      x = venn_centres$centre_x[j] + cos(venn_theta),
      y = venn_centres$centre_y[j] + sin(venn_theta)
    )
  }))
})) %>%
  mutate(panel = factor(panel, levels = classic_panel_spec$panel))

venn_region_positions <- tibble(
  region = classic_region_order,
  x = c(-1.08, 1.08, 0, 0, -0.52, 0.52, 0),
  y = c(0.58, 0.58, -1.20, 0.82, -0.42, -0.42, -0.05)
)
venn_count_labels <- classic_venn_counts %>%
  left_join(venn_region_positions, by = "region") %>%
  mutate(
    panel = factor(panel, levels = classic_panel_spec$panel),
    label = comma(gene_n)
  )

venn_set_labels <- bind_rows(lapply(seq_len(nrow(classic_panel_spec)), function(i) {
  sc_key <- classic_panel_spec$scatlas_key[i]
  pdo_key <- classic_panel_spec$pdo_key[i]
  tibble(
    panel = classic_panel_spec$panel[i],
    set = c("scATLAS", "CRISPR", "PDO"),
    x = c(-1.48, 1.48, 0),
    y = c(1.48, 1.48, -1.72),
    label = c(
      paste0(
        "scATLAS\n", classic_panel_spec$scatlas_label[i],
        "\n(n = ", comma(length(classic_signature_sets[[sc_key]])), ")"
      ),
      paste0("Oesophageal CRISPR\n(n = ", comma(length(classic_signature_sets$crispr)), ")"),
      paste0(
        "PDO\n", classic_panel_spec$pdo_label[i],
        "\n(n = ", comma(length(classic_signature_sets[[pdo_key]])), ")"
      )
    )
  )
})) %>%
  mutate(panel = factor(panel, levels = classic_panel_spec$panel))

classic_venn_plot <- ggplot() +
  geom_polygon(
    data = venn_circle_df,
    aes(x = x, y = y, group = interaction(panel, set), fill = set),
    alpha = 0.10,
    colour = NA
  ) +
  geom_path(
    data = venn_circle_df,
    aes(x = x, y = y, group = interaction(panel, set), colour = set),
    linewidth = 1.05
  ) +
  geom_text(
    data = venn_count_labels,
    aes(x = x, y = y, label = label),
    size = 4.2,
    fontface = "bold"
  ) +
  geom_text(
    data = venn_set_labels,
    aes(x = x, y = y, label = label, colour = set),
    size = 3.6,
    fontface = "bold",
    lineheight = 0.92
  ) +
  facet_wrap(~panel, ncol = 2) +
  scale_colour_manual(values = setNames(venn_centres$colour, venn_centres$set)) +
  scale_fill_manual(values = setNames(venn_centres$colour, venn_centres$set)) +
  coord_fixed(xlim = c(-1.85, 1.85), ylim = c(-1.85, 1.65), clip = "off") +
  labs(
    title = "Classic proliferation: three-way gene overlap"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey94", colour = "grey75"),
    panel.spacing = grid::unit(1.1, "lines"),
    plot.margin = margin(10, 18, 10, 18)
  )
classic_venn_figure_file <- file.path(
  out_paths[["figures"]],
  "Auto_classic_proliferation_scatlas_crispr_pdo_threeway_venn.pdf"
)
pdo_save_slide_pdf(
  classic_venn_plot,
  classic_venn_figure_file,
  width = 16,
  height = 12
)

####################
# Simplified summary figures
####################
build_simple_plot <- function(df, title, y_var, shape_var = NULL) {
  df <- df %>%
    mutate(
      log10_p = -log10(expression_matched_burden_p),
      is_sig = expression_matched_burden_p < 0.05
    )
  
  p <- ggplot(df, aes(x = log10_p, y = .data[[y_var]], color = is_sig)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey50")
    
  if (!is.null(shape_var)) {
    p <- p + geom_point(aes(shape = .data[[shape_var]]), size = 5, alpha = 0.9) +
             labs(shape = "Signature Type")
  } else {
    p <- p + geom_point(size = 5, alpha = 0.9)
  }
  
  p <- p +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"), guide = "none") +
    labs(
      title = title,
      x = "-log10(P-value) [Depletion burden]",
      y = NULL
    ) +
    pdo_theme_slide(12) +
    theme(legend.position = "bottom")
  
  return(p)
}

p_scatlas_mp <- build_simple_plot(scatlas_mp_plot_df, "scATLAS MPs", "plot_label")
p_scatlas_state <- build_simple_plot(scatlas_state_plot_df, "scATLAS States", "signature", "signature_type")

p_pdo_mp <- build_simple_plot(mp_plot_df, "PDO MPs", "plot_label")
p_pdo_state <- build_simple_plot(state_plot_df, "PDO States", "signature", "signature_type")

p_scatlas_page <- p_scatlas_mp + p_scatlas_state + 
  plot_layout(ncol = 2) + 
  plot_annotation(title = "scATLAS Selective Dependency Summary", 
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

p_pdo_page <- p_pdo_mp + p_pdo_state + 
  plot_layout(ncol = 2) + 
  plot_annotation(title = "PDO Selective Dependency Summary", 
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

simple_summary_file <- file.path(
  out_paths[["figures"]],
  "Auto_simplified_selective_dependency_summary.pdf"
)

pdf(simple_summary_file, width = 20, height = 9)
print(p_scatlas_page)
print(p_pdo_page)
invisible(dev.off())

####################
# Methodology copy and auditable run summary
####################
invisible(file.copy(
  methodology_file,
  file.path(out_paths[["reports"]], basename(methodology_file)),
  overwrite = TRUE
))

output_files <- c(
  workbook_out,
  file.path(out_paths[["tables"]], "Auto_mp_selective_dependency_summary.csv"),
  file.path(out_paths[["tables"]], "Auto_state_selective_dependency_summary.csv"),
  file.path(out_paths[["tables"]], "Auto_state_edger_dge_all.csv.gz"),
  file.path(out_paths[["tables"]], "Auto_state_edger_signature_top150.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_mp_selective_dependency_summary.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_state_selective_dependency_summary.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_mp_dependency_gene_details.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_state_combined_mp_dependency_gene_details.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_state_marker_dependency_gene_details.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_state_marker_signature_top150.csv"),
  file.path(out_paths[["tables"]], "Auto_scatlas_expression_dependency_universe.csv.gz"),
  classic_venn_counts_file,
  classic_venn_membership_file,
  classic_threeway_genes_file,
  classic_venn_sets_file,
  file.path(out_paths[["figures"]], "Auto_mp_selective_dependency_summary.pdf"),
  file.path(out_paths[["figures"]], "Auto_state_selective_dependency_summary.pdf"),
  scatlas_mp_figure_file,
  scatlas_state_figure_file,
  classic_venn_figure_file,
  simple_summary_file,
  score_cache_file,
  scatlas_score_cache_file
)
missing_outputs <- output_files[!file.exists(output_files)]
if (length(missing_outputs) > 0) {
  stop("Required output(s) were not produced: ", paste(missing_outputs, collapse = ", "))
}

pdo_write_run_summary(
  script = "analysis/cell_states/Auto_crispr_selective_dependency_analysis.R",
  out_dir = out_dir,
  inputs = c(
    workbook_file, seurat_file, state_file, mp_gene_file, mp_weight_file,
    scatlas_mp_gene_file, scatlas_mp_weight_file, scatlas_seurat_file,
    scatlas_state_file, scatlas_state_marker_file, scref_config_file
  ),
  outputs = output_files,
  parameters = params,
  cache = list(force_rebuild = force_rebuild, replot_only = replot_only),
  status = "completed"
)

message("CRISPR selective-dependency analysis completed successfully.")
