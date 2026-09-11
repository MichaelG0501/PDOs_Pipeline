####################
# Analysis registry:
#   Status: active terminal clinical and FLOT-sensitivity association workflow
#   Script: analysis/clinical/Auto_centred_clinical_auc_response_associations.R
#   Methodology: analysis/methodology/clinical/Auto_centred_clinical_auc_response_associations_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Tests and visualises baseline PDO state proportions and
#     finalized centred-refined MP activity across clinical variables, models
#     their linear relationships with FLOT AUCrel, and compares the three most
#     sensitive with the three most resistant PDOs using state proportions,
#     overall/state-resolved MP activity, and selected state-resolved pathway
#     scores.
#   Inputs:
#     PDOs_outs/PDOs_merged.rds
#     PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#     PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds
#     PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     PDOs_outs/centred_mp_refinement/centred_refined_mp_strict_order.rds
#     /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDO_ClinicalMetadata_V3.xlsx
#     /rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv
#   Outputs:
#     PDOs_outs/Auto_centred_clinical_auc_response_associations/intermediate/*
#     PDOs_outs/Auto_centred_clinical_auc_response_associations/tables/*
#     PDOs_outs/Auto_centred_clinical_auc_response_associations/figures/*
#     PDOs_outs/Auto_centred_clinical_auc_response_associations/logs/*
#     PDOs_outs/Auto_centred_clinical_auc_response_associations/reports/*
#   Downstream: terminal presentation and source-data outputs only
#   Cache/replot: PDO_FORCE_REBUILD=1 recomputes the persistent live cache;
#     PDO_REPLOT_ONLY=1 requires and replots that cache.
#   Run: qsub analysis/clinical/Auto_centred_clinical_auc_response_associations.sh
#   Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################

####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(data.table)
  library(msigdbr)
})
####################

####################
# Paths, constants, and guards
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

model_dir <- file.path(PDO_OUTPUT_DIR, "centred_mp_refinement")
out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_centred_clinical_auc_response_associations")
out_paths <- pdo_ensure_output_tiers(out_dir)
intermediate_dir <- out_paths[["intermediate"]]
tables_dir <- out_paths[["tables"]]
figures_dir <- out_paths[["figures"]]
reports_dir <- out_paths[["reports"]]

input_paths <- c(
  pdo = file.path(PDO_OUTPUT_DIR, "PDOs_merged.rds"),
  states = file.path(model_dir, "centred_refined_noreg_states.rds"),
  ucell = file.path(model_dir, "merged_refined_ucell_scores.rds"),
  mp_genes = file.path(model_dir, "merged_refined_mp_genes.rds"),
  strict_order = file.path(model_dir, "centred_refined_mp_strict_order.rds"),
  clinical = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDO_ClinicalMetadata_V3.xlsx",
  cell_cycle = PDO_EXTERNAL_PATHS$cell_cycle_genes
)
pdo_require_files(input_paths)

state_order <- PDO_STATE_ORDER
state_order_all <- PDO_STATE_ORDER_WITH_OPTIONAL
state_colours <- PDO_STATE_COLORS
sensitive_sur <- c("SUR680", "SUR629", "SUR1090")
resistant_sur <- c("SUR727", "SUR1181", "SUR1363")
comparison_sur <- c(sensitive_sur, resistant_sur)
comparison_colours <- c(Sensitive = "#0072B2", Resistant = "#D55E00")
min_cells_per_pseudobulk <- pdo_get_env_integer("PDO_CLINICAL_MIN_CELLS", 20L)
cache_policy <- pdo_cache_policy()
cache_path <- file.path(intermediate_dir, "Auto_centred_clinical_auc_response_results.rds")
start_time <- Sys.time()

clinical_specs <- tibble::tribble(
  ~variable, ~label,
  "gender", "Gender",
  "age_group", "Age (>60)",
  "tumour_location", "Tumour location",
  "tumour_type", "Tumour type",
  "histology", "Histology",
  "ajcc", "AJCC",
  "clinical_response", "Clinical response at OG MDT",
  "mandard_score", "Mandard tumour regression score",
  "mandard_response", "Response based on Mandard",
  "pdo_origin", "Origin of PDO",
  "ctx_timepoint", "Pre CTX or Post CTX"
)

state_axis_labels <- c(
  "Classic proliferation" = "Classic\nproliferation",
  "Columnar-to-intestinal" = "Columnar-to-\nintestinal",
  "Glandular differentiation" = "Glandular\ndifferentiation",
  "Stress-adaptive" = "Stress-\nadaptive",
  "ECM-remodelling" = "ECM-\nremodelling",
  "Motile-cilia differentiation" = "Motile-cilia\ndifferentiation"
)

state_axis_labels_all <- c(
  state_axis_labels,
  "Unresolved" = "Unresolved",
  "Hybrid" = "Hybrid"
)

clean_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "N/A", "NA", "n/a", "Unknown")] <- NA_character_
  x
}

normalise_response <- function(x) {
  x <- clean_missing(x)
  case_when(
    str_detect(x, regex("^non[- ]?responder$|^NR$", ignore_case = TRUE)) ~ "Non-responder",
    str_detect(x, regex("^responder$|^R$", ignore_case = TRUE)) ~ "Responder",
    TRUE ~ x
  )
}

safe_wilcox_or_kruskal <- function(data) {
  data <- data %>% filter(is.finite(value), !is.na(group))
  group_n <- n_distinct(data$group)
  if (group_n < 2L || n_distinct(data$sample) < 3L) {
    return(tibble(test = NA_character_, p_value = NA_real_))
  }
  if (group_n == 2L) {
    result <- tryCatch(wilcox.test(value ~ group, data = data, exact = FALSE), error = function(e) NULL)
    return(tibble(test = "Wilcoxon rank-sum", p_value = if (is.null(result)) NA_real_ else result$p.value))
  }
  result <- tryCatch(kruskal.test(value ~ group, data = data), error = function(e) NULL)
  tibble(test = "Kruskal-Wallis", p_value = if (is.null(result)) NA_real_ else result$p.value)
}

significance_label <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
}

compute_categorical_stats <- function(data, family_name) {
  data %>%
    group_by(clinical_variable, clinical_label, feature) %>%
    group_modify(~safe_wilcox_or_kruskal(.x)) %>%
    ungroup() %>%
    group_by(clinical_variable) %>%
    mutate(fdr_within_variable = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(feature_family = family_name)
}

compute_auc_lm <- function(data, family_name) {
  data %>%
    filter(is.finite(auc_rel), is.finite(value)) %>%
    group_by(feature) %>%
    group_modify(~{
      fit <- lm(value ~ auc_rel, data = .x)
      coef_table <- summary(fit)$coefficients
      ci <- confint(fit, "auc_rel", level = 0.95)
      tibble(
        n_samples = nrow(.x),
        slope = unname(coef(fit)[["auc_rel"]]),
        slope_ci_low = ci[1],
        slope_ci_high = ci[2],
        p_value = coef_table["auc_rel", "Pr(>|t|)"],
        r_squared = summary(fit)$r.squared
      )
    }) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "BH"), feature_family = family_name)
}

compute_two_group_stats <- function(data, family_name, strata = character()) {
  grouping <- c(strata, "feature")
  data %>%
    filter(is.finite(value), !is.na(sensitivity_group)) %>%
    group_by(across(all_of(grouping))) %>%
    group_modify(~{
      group_summary <- .x %>%
        group_by(sensitivity_group) %>%
        summarise(n = n_distinct(sample), median = median(value), mean = mean(value), .groups = "drop")
      p_value <- if (n_distinct(.x$sensitivity_group) == 2L && all(group_summary$n >= 2L)) {
        tryCatch(wilcox.test(value ~ sensitivity_group, data = .x, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
      sensitive_median <- group_summary$median[group_summary$sensitivity_group == "Sensitive"]
      resistant_median <- group_summary$median[group_summary$sensitivity_group == "Resistant"]
      sensitive_mean <- group_summary$mean[group_summary$sensitivity_group == "Sensitive"]
      resistant_mean <- group_summary$mean[group_summary$sensitivity_group == "Resistant"]
      tibble(
        sensitive_n = group_summary$n[group_summary$sensitivity_group == "Sensitive"] %||% NA_integer_,
        resistant_n = group_summary$n[group_summary$sensitivity_group == "Resistant"] %||% NA_integer_,
        sensitive_mean = sensitive_mean %||% NA_real_,
        resistant_mean = resistant_mean %||% NA_real_,
        mean_difference_resistant_minus_sensitive = (resistant_mean %||% NA_real_) - (sensitive_mean %||% NA_real_),
        sensitive_median = sensitive_median %||% NA_real_,
        resistant_median = resistant_median %||% NA_real_,
        median_difference_resistant_minus_sensitive = (resistant_median %||% NA_real_) - (sensitive_median %||% NA_real_),
        p_value = p_value
      )
    }) %>%
    ungroup() %>%
    group_by(across(all_of(strata))) %>%
    mutate(fdr_within_stratum = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(feature_family = family_name)
}

`%||%` <- function(x, y) if (length(x) == 0L || is.na(x[1])) y else x[1]

mean_signature <- function(expression_matrix, genes) {
  genes_use <- intersect(unique(genes), rownames(expression_matrix))
  if (length(genes_use) == 0L) return(rep(NA_real_, ncol(expression_matrix)))
  colMeans(expression_matrix[genes_use, , drop = FALSE])
}

write_table <- function(x, filename) {
  data.table::fwrite(as.data.frame(x), file.path(tables_dir, filename))
}
####################

####################
# Compute and cache auditable source data
####################
if (!file.exists(cache_path) || cache_policy$force_rebuild) {
  if (cache_policy$replot_only) {
    stop("PDO_REPLOT_ONLY=1 but cache is missing: ", cache_path)
  }

  message("Loading current centred PDO inputs ...")
  pdos <- readRDS(input_paths[["pdo"]])
  state_vector <- pdo_normalise_final_state_vector(readRDS(input_paths[["states"]]))
  ucell <- readRDS(input_paths[["ucell"]])
  strict_object <- readRDS(input_paths[["strict_order"]])
  refined_mp_genes <- readRDS(input_paths[["mp_genes"]])
  mp_order <- strict_object$strict_refined_mp_order
  if (is.null(mp_order)) stop("Strict-order object lacks strict_refined_mp_order.")
  missing_mps <- setdiff(mp_order, colnames(ucell))
  if (length(missing_mps) > 0L) stop("UCell input lacks current MPs: ", paste(missing_mps, collapse = "; "))
  if (is.null(names(state_vector))) stop("Current state vector is not cell-named.")

  clinical <- readxl::read_excel(input_paths[["clinical"]], sheet = "Sheet1") %>%
    transmute(
      sur = paste0("SUR", as.integer(SUR)),
      gender = clean_missing(Gender),
      age = suppressWarnings(as.numeric(Age)),
      age_group = case_when(is.na(age) ~ NA_character_, age > 60 ~ ">60", TRUE ~ "<=60"),
      tumour_location = clean_missing(`Tumour location`),
      tumour_type = clean_missing(`Tumour type`),
      histology = clean_missing(Histology),
      ajcc = clean_missing(AJCC),
      clinical_response = normalise_response(`Clinical response at OG MDT: responder/non-responder`),
      mandard_score = clean_missing(`Mandard tumour regression score`),
      mandard_response = normalise_response(`Response based on Mandard`),
      pdo_origin = case_when(
        clean_missing(`Origin of PDO (primary tumour vs lymph node)`) == "T" ~ "Primary tumour",
        clean_missing(`Origin of PDO (primary tumour vs lymph node)`) == "LN" ~ "Lymph node",
        TRUE ~ clean_missing(`Origin of PDO (primary tumour vs lymph node)`)
      ),
      ctx_timepoint = clean_missing(`Pre CTX or Post CTX`),
      auc_rel = suppressWarnings(as.numeric(clean_missing(`AUCrel (FLOT)`)))
    )

  meta <- pdos@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
      cell,
      sample = as.character(orig.ident),
      batch = ifelse(is.na(Batch) | !nzchar(as.character(Batch)), "Unknown", as.character(Batch)),
      sur = str_extract(sample, "^SUR[0-9]+"),
      state = state_vector[cell]
    ) %>%
    filter(
      sample != PDO_EXCLUDED_SAMPLE,
      !str_detect(sample, "_Treated_PDO$")
    ) %>%
    left_join(clinical, by = "sur")

  if (anyNA(meta$state)) stop("Baseline cells are missing current state assignments.")
  duplicated_sur <- meta %>% distinct(sample, sur) %>% count(sur) %>% filter(n > 1L)
  if (nrow(duplicated_sur) > 0L) {
    stop("More than one baseline PDO sample was found for: ", paste(duplicated_sur$sur, collapse = "; "))
  }
  missing_clinical_sur <- meta %>% distinct(sur) %>% filter(!sur %in% clinical$sur)
  if (nrow(missing_clinical_sur) > 0L) {
    stop("Baseline PDOs lack V3 clinical rows: ", paste(missing_clinical_sur$sur, collapse = "; "))
  }

  baseline_cells <- meta$cell
  missing_ucell_cells <- setdiff(baseline_cells, rownames(ucell))
  if (length(missing_ucell_cells) > 0L) stop("UCell matrix misses baseline cells: ", length(missing_ucell_cells))

  sample_meta <- meta %>%
    distinct(sample, sur, batch, gender, age, age_group, tumour_location, tumour_type,
             histology, ajcc, clinical_response, mandard_score, mandard_response,
             pdo_origin, ctx_timepoint, auc_rel) %>%
    mutate(
      sensitivity_group = case_when(
        sur %in% sensitive_sur ~ "Sensitive",
        sur %in% resistant_sur ~ "Resistant",
        TRUE ~ NA_character_
      ),
      sensitivity_group = factor(sensitivity_group, levels = c("Sensitive", "Resistant"))
    )

  auc_samples <- sample_meta %>% filter(is.finite(auc_rel))
  if (nrow(auc_samples) != 12L) stop("Expected all 12 numeric AUCrel PDOs; observed ", nrow(auc_samples))
  if (!setequal(comparison_sur, sample_meta$sur[!is.na(sample_meta$sensitivity_group)])) {
    stop("Sensitive/resistant PDO membership did not resolve exactly to the requested six PDOs.")
  }

  state_sample_all <- meta %>%
    count(sample, sur, state, name = "cell_n") %>%
    complete(
      nesting(sample, sur),
      state = state_order_all,
      fill = list(cell_n = 0L)
    ) %>%
    group_by(sample, sur) %>%
    mutate(total_cells = sum(cell_n), value = 100 * cell_n / total_cells) %>%
    ungroup() %>%
    mutate(feature = factor(state, levels = state_order_all)) %>%
    left_join(sample_meta, by = c("sample", "sur"))
  state_sample <- state_sample_all %>%
    filter(state %in% state_order) %>%
    mutate(feature = factor(as.character(state), levels = state_order))

  mp_sample <- as.data.frame(ucell[baseline_cells, mp_order, drop = FALSE]) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = all_of(mp_order), names_to = "feature", values_to = "ucell_score") %>%
    left_join(meta %>% select(cell, sample, sur), by = "cell") %>%
    group_by(sample, sur, feature) %>%
    summarise(value = mean(ucell_score, na.rm = TRUE), cell_n = n(), .groups = "drop") %>%
    mutate(feature = factor(feature, levels = mp_order)) %>%
    left_join(sample_meta, by = c("sample", "sur"))

  categorical_state <- bind_rows(lapply(seq_len(nrow(clinical_specs)), function(i) {
    variable_name <- clinical_specs$variable[i]
    state_sample %>%
      transmute(
        sample, sur, feature, value,
        clinical_variable = variable_name,
        clinical_label = clinical_specs$label[i],
        group = clean_missing(.data[[variable_name]])
      ) %>%
      filter(!is.na(group))
  }))
  categorical_mp <- bind_rows(lapply(seq_len(nrow(clinical_specs)), function(i) {
    variable_name <- clinical_specs$variable[i]
    mp_sample %>%
      transmute(
        sample, sur, feature, value,
        clinical_variable = variable_name,
        clinical_label = clinical_specs$label[i],
        group = clean_missing(.data[[variable_name]])
      ) %>%
      filter(!is.na(group))
  }))
  categorical_state_stats <- compute_categorical_stats(categorical_state, "State proportion")
  categorical_mp_stats <- compute_categorical_stats(categorical_mp, "MP UCell")

  auc_state <- state_sample %>% filter(is.finite(auc_rel)) %>% select(sample, sur, feature, value, auc_rel)
  auc_mp <- mp_sample %>% filter(is.finite(auc_rel)) %>% select(sample, sur, feature, value, auc_rel)
  auc_state_stats <- compute_auc_lm(auc_state, "State proportion")
  auc_mp_stats <- compute_auc_lm(auc_mp, "MP UCell")

  comparison_state <- state_sample %>% filter(!is.na(sensitivity_group))
  comparison_state_all <- state_sample_all %>% filter(!is.na(sensitivity_group))
  comparison_mp <- mp_sample %>% filter(!is.na(sensitivity_group))
  comparison_state_stats <- compute_two_group_stats(comparison_state, "State proportion")
  comparison_mp_stats <- compute_two_group_stats(comparison_mp, "MP UCell")

  comparison_meta <- meta %>%
    filter(sur %in% comparison_sur, state %in% state_order) %>%
    mutate(
      state = factor(state, levels = state_order),
      sensitivity_group = factor(ifelse(sur %in% sensitive_sur, "Sensitive", "Resistant"), levels = c("Sensitive", "Resistant")),
      sample_state = paste(sample, state, sep = "__")
    )
  comparison_cells <- comparison_meta$cell

  state_mp <- as.data.frame(ucell[comparison_cells, mp_order, drop = FALSE]) %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = all_of(mp_order), names_to = "feature", values_to = "ucell_score") %>%
    left_join(comparison_meta %>% select(cell, sample, sur, state, sensitivity_group), by = "cell") %>%
    group_by(sample, sur, sensitivity_group, state, feature) %>%
    summarise(value = mean(ucell_score, na.rm = TRUE), cell_n = n(), .groups = "drop") %>%
    filter(cell_n >= min_cells_per_pseudobulk) %>%
    mutate(feature = factor(feature, levels = mp_order))
  state_mp_stats <- compute_two_group_stats(state_mp, "State-resolved MP UCell", strata = "state")

  message("Computing state-resolved selected pathway scores for the requested six PDOs ...")
  comparison_object <- subset(pdos, cells = comparison_cells)
  counts <- pdo_get_assay_matrix(comparison_object, assay = "RNA", layer = "counts")
  grouping <- factor(comparison_meta$sample_state, levels = unique(comparison_meta$sample_state))
  design_sparse <- Matrix::sparse.model.matrix(~0 + grouping)
  colnames(design_sparse) <- levels(grouping)
  pseudobulk_counts <- counts[, comparison_meta$cell, drop = FALSE] %*% design_sparse
  pseudobulk_meta <- comparison_meta %>%
    count(sample_state, sample, sur, sensitivity_group, state, name = "cell_n") %>%
    arrange(state, sensitivity_group, sur)
  pseudobulk_counts <- pseudobulk_counts[, pseudobulk_meta$sample_state, drop = FALSE]
  pseudobulk_dge <- edgeR::DGEList(counts = pseudobulk_counts)
  pseudobulk_dge <- edgeR::calcNormFactors(pseudobulk_dge)
  pseudobulk_logcpm <- edgeR::cpm(pseudobulk_dge, log = TRUE, prior.count = 2)
  pseudobulk_z <- t(scale(t(pseudobulk_logcpm)))
  pseudobulk_z[!is.finite(pseudobulk_z)] <- 0

  hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  hallmark_ids <- c(
    HALLMARK_E2F_TARGETS = "E2F targets",
    HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
    HALLMARK_APOPTOSIS = "Apoptosis",
    HALLMARK_P53_PATHWAY = "p53 pathway",
    HALLMARK_DNA_REPAIR = "DNA repair",
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF/NF-kB",
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
    HALLMARK_HYPOXIA = "Hypoxia",
    HALLMARK_XENOBIOTIC_METABOLISM = "Xenobiotic metabolism",
    HALLMARK_UNFOLDED_PROTEIN_RESPONSE = "Unfolded protein response",
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxidative phosphorylation"
  )
  pathway_sets <- lapply(names(hallmark_ids), function(id) unique(hallmark$gene_symbol[hallmark$gs_name == id]))
  names(pathway_sets) <- unname(hallmark_ids)
  pathway_sets[["Interferon response"]] <- unique(hallmark$gene_symbol[
    hallmark$gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
  ])

  cell_cycle <- read.csv(input_paths[["cell_cycle"]], stringsAsFactors = FALSE)
  if (!all(c("Gene", "Consensus") %in% colnames(cell_cycle))) {
    stop("Cell-cycle reference lacks Gene and Consensus columns.")
  }
  cc_consensus <- intersect(cell_cycle$Gene[cell_cycle$Consensus == 1], rownames(pseudobulk_logcpm))
  cc_top50 <- names(sort(rowMeans(pseudobulk_logcpm[cc_consensus, , drop = FALSE]), decreasing = TRUE))[
    seq_len(min(50L, length(cc_consensus)))
  ]
  pathway_sets[["CCSIG"]] <- cc_top50
  pathway_sets[["Intestinal metaplasia"]] <- refined_mp_genes[["MP15"]]
  pathway_sets[["Ciliated progenitor epithelium"]] <- refined_mp_genes[["MP17+"]]
  if (any(vapply(pathway_sets, is.null, logical(1)))) stop("One or more selected pathway signatures are unavailable.")

  pathway_matrix <- vapply(pathway_sets, function(genes) mean_signature(pseudobulk_z, genes), numeric(ncol(pseudobulk_z)))
  rownames(pathway_matrix) <- colnames(pseudobulk_z)
  pathway_scores <- as.data.frame(pathway_matrix, check.names = FALSE) %>%
    rownames_to_column("sample_state") %>%
    left_join(pseudobulk_meta, by = "sample_state") %>%
    pivot_longer(cols = all_of(names(pathway_sets)), names_to = "feature", values_to = "value") %>%
    filter(cell_n >= min_cells_per_pseudobulk, is.finite(value))
  pathway_stats <- compute_two_group_stats(pathway_scores, "State-resolved pathway score", strata = "state")

  result_cache <- list(
    sample_meta = sample_meta,
    state_sample_all = state_sample_all,
    state_sample = state_sample,
    mp_sample = mp_sample,
    categorical_state = categorical_state,
    categorical_mp = categorical_mp,
    categorical_state_stats = categorical_state_stats,
    categorical_mp_stats = categorical_mp_stats,
    auc_state = auc_state,
    auc_mp = auc_mp,
    auc_state_stats = auc_state_stats,
    auc_mp_stats = auc_mp_stats,
    comparison_state = comparison_state,
    comparison_state_all = comparison_state_all,
    comparison_mp = comparison_mp,
    comparison_state_stats = comparison_state_stats,
    comparison_mp_stats = comparison_mp_stats,
    state_mp = state_mp,
    state_mp_stats = state_mp_stats,
    pseudobulk_meta = pseudobulk_meta,
    pseudobulk_logcpm = pseudobulk_logcpm,
    pathway_sets = pathway_sets,
    pathway_scores = pathway_scores,
    pathway_stats = pathway_stats,
    clinical_specs = clinical_specs,
    state_order = state_order,
    mp_order = mp_order
  )
  saveRDS(result_cache, cache_path)
  saveRDS(pathway_sets, file.path(intermediate_dir, "Auto_selected_pathway_gene_sets.rds"))
  saveRDS(pseudobulk_logcpm, file.path(intermediate_dir, "Auto_sensitive_resistant_state_pseudobulk_logcpm.rds"))
} else {
  message("Reusing persistent clinical/AUC response cache ...")
  result_cache <- readRDS(cache_path)
}

list2env(result_cache, envir = environment())
####################

####################
# Persistent source-data tables
####################
write_table(sample_meta, "Auto_sample_clinical_metadata.csv")
write_table(state_sample_all, "Auto_state_proportion_all_states_source_data.csv")
write_table(state_sample, "Auto_state_proportion_source_data.csv")
write_table(mp_sample, "Auto_mp_ucell_source_data.csv")
write_table(categorical_state_stats, "Auto_categorical_state_association_stats.csv")
write_table(categorical_mp_stats, "Auto_categorical_mp_association_stats.csv")
write_table(auc_state_stats, "Auto_auc_state_linear_model_stats.csv")
write_table(auc_mp_stats, "Auto_auc_mp_linear_model_stats.csv")
write_table(comparison_state_stats, "Auto_sensitive_resistant_state_stats.csv")
write_table(comparison_mp_stats, "Auto_sensitive_resistant_mp_stats.csv")
write_table(state_mp, "Auto_sensitive_resistant_state_resolved_mp_source_data.csv")
write_table(state_mp_stats, "Auto_sensitive_resistant_state_resolved_mp_stats.csv")
write_table(pseudobulk_meta, "Auto_sensitive_resistant_pseudobulk_metadata.csv")
write_table(pathway_scores, "Auto_sensitive_resistant_state_pathway_scores.csv")
write_table(pathway_stats, "Auto_sensitive_resistant_state_pathway_stats.csv")
####################

####################
# Plot helpers and categorical clinical figures
####################
mp_labels <- setNames(
  ifelse(is.na(PDO_MP_DESCRIPTIONS[mp_order]), mp_order, paste0(mp_order, "\n", PDO_MP_DESCRIPTIONS[mp_order])),
  mp_order
)

clinical_palette <- function(groups) {
  base <- c(
    Female = "#CC79A7", Male = "#0072B2", `<=60` = "#00A087", `>60` = "#E69F00",
    Responder = "#009E73", `Non-responder` = "#D55E00", R = "#009E73", NR = "#D55E00",
    Pre = "#56B4E9", Post = "#E69F00", `Primary tumour` = "#4DAF4A", `Lymph node` = "#984EA3"
  )
  missing <- setdiff(groups, names(base))
  if (length(missing) > 0L) base[missing] <- hue_pal()(length(missing))
  base[groups]
}

plot_categorical_features <- function(data, stats, label_map, y_label, title_suffix, percent_axis = FALSE) {
  groups <- sort(unique(as.character(data$group)))
  legend_counts <- data %>% distinct(sample, group) %>% count(group, name = "sample_n")
  legend_labels <- setNames(
    paste0(legend_counts$group, " (n=", legend_counts$sample_n, ")"),
    legend_counts$group
  )
  stats_label <- stats %>%
    mutate(label = significance_label(p_value))
  y_top <- data %>%
    group_by(feature) %>%
    summarise(
      y_max = max(value, na.rm = TRUE),
      y_min = min(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_span = pmax(y_max - y_min, ifelse(percent_axis, 5, 0.02)),
      y = y_max + pmax(ifelse(percent_axis, 2, 0.01), 0.08 * y_span)
    )
  stats_label <- left_join(stats_label, y_top, by = "feature")
  p <- ggplot(data, aes(feature, value, fill = group, colour = group)) +
    geom_boxplot(position = position_dodge(width = 0.78), width = 0.65, outlier.shape = NA, alpha = 0.78, colour = "black") +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.78), size = 1.5, alpha = 0.85, show.legend = FALSE) +
    geom_text(data = stats_label %>% filter(label != ""), aes(feature, y, label = label), inherit.aes = FALSE, size = 5.5, fontface = "bold") +
    scale_fill_manual(values = clinical_palette(groups), labels = legend_labels, drop = FALSE) +
    scale_colour_manual(values = clinical_palette(groups), guide = "none", drop = FALSE) +
    scale_x_discrete(labels = label_map) +
    labs(title = paste0(unique(data$clinical_label), " - ", title_suffix), x = NULL, y = y_label, fill = "Clinical group") +
    pdo_theme_slide(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
    coord_cartesian(clip = "off")
  if (percent_axis) p <- p + scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0.02, 0.18)))
  else p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)))
  p
}

stacked_pdf <- file.path(figures_dir, "Auto_categorical_state_composition_stacked.pdf")
pdf(stacked_pdf, width = 12, height = 8, useDingbats = FALSE)
for (i in seq_len(nrow(clinical_specs))) {
  variable_name <- clinical_specs$variable[i]
  stacked_source <- state_sample_all %>%
    mutate(group = clean_missing(.data[[variable_name]])) %>%
    filter(!is.na(group))
  data_use <- stacked_source %>%
    group_by(group, feature) %>%
    summarise(mean_percent = mean(value), .groups = "drop")
  total_samples <- n_distinct(stacked_source$sample)
  total_batches <- n_distinct(stacked_source$batch)
  group_counts <- stacked_source %>%
    distinct(sample, group, batch, total_cells) %>%
    group_by(group) %>%
    summarise(
      cell_n = sum(total_cells),
      sample_n = n_distinct(sample),
      batch_n = n_distinct(batch),
      batch_list = paste(sort(unique(batch)), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(
      "Cells=", comma(cell_n), "\n",
      "Samples=", sample_n, "/", total_samples, "\n",
      "Batches=", batch_n, "/", total_batches, "\n",
      batch_list
    ))
  p <- ggplot(data_use, aes(group, mean_percent / 100, fill = feature)) +
    geom_col(width = 0.72) +
    geom_text(
      aes(label = ifelse(mean_percent >= 8, sprintf("%.0f%%", mean_percent), "")),
      position = position_stack(vjust = 0.5), size = 3.8, fontface = "bold"
    ) +
    geom_text(data = group_counts, aes(group, 1.02, label = label), inherit.aes = FALSE, vjust = 0, size = 3.8, lineheight = 1, fontface = "bold") +
    scale_fill_manual(values = state_colours[state_order_all], labels = state_axis_labels_all[state_order_all], drop = FALSE) +
    scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.34))) +
    labs(title = clinical_specs$label[i], x = NULL, y = "Mean sample state proportion", fill = "Cell state") +
    pdo_theme_slide(13) +
    coord_cartesian(clip = "off") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right", plot.margin = margin(10, 20, 10, 10))
  print(p)
}
dev.off()

pdf(file.path(figures_dir, "Auto_categorical_state_association_boxplots.pdf"), width = 15, height = 8.5, useDingbats = FALSE)
for (i in seq_len(nrow(clinical_specs))) {
  variable_name <- clinical_specs$variable[i]
  data_use <- categorical_state %>% filter(clinical_variable == variable_name)
  stats_use <- categorical_state_stats %>% filter(clinical_variable == variable_name)
  if (n_distinct(data_use$group) >= 2L) print(plot_categorical_features(
    data_use, stats_use, state_axis_labels, "Sample state proportion", "PDO states", TRUE
  ))
}
dev.off()

pdf(file.path(figures_dir, "Auto_categorical_mp_association_boxplots.pdf"), width = 18, height = 9, useDingbats = FALSE)
for (i in seq_len(nrow(clinical_specs))) {
  variable_name <- clinical_specs$variable[i]
  data_use <- categorical_mp %>% filter(clinical_variable == variable_name)
  stats_use <- categorical_mp_stats %>% filter(clinical_variable == variable_name)
  if (n_distinct(data_use$group) >= 2L) print(plot_categorical_features(
    data_use, stats_use, mp_labels, "Mean sample UCell score", "centred-refined MPs", FALSE
  ))
}
dev.off()
####################

####################
# AUCrel linear-pattern figures (all 12 numeric PDOs)
####################
plot_auc_linear <- function(data, stats, label_map, y_label, title_text) {
  annotation <- stats %>%
    mutate(
      feature = factor(feature, levels = levels(data$feature)),
      significant = ifelse(is.finite(p_value) & p_value < 0.05, "P < 0.05", "P >= 0.05"),
      label = paste0("R²=", sprintf("%.2f", r_squared), "\nP=", format.pval(p_value, digits = 2), ifelse(significant == "P < 0.05", " *", ""))
    )
  plot_data <- data %>%
    left_join(annotation %>% select(feature, significant), by = "feature")
  ggplot(plot_data, aes(auc_rel, value)) +
    geom_smooth(aes(colour = significant, fill = significant), method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
    geom_point(size = 2.2, alpha = 0.9, colour = "#333333") +
    geom_text(aes(label = sur), nudge_y = 0.015 * diff(range(plot_data$value, na.rm = TRUE)), size = 2.4, check_overlap = TRUE, colour = "#333333") +
    geom_text(data = annotation, aes(x = -Inf, y = Inf, label = label, colour = significant), inherit.aes = FALSE, hjust = -0.08, vjust = 1.15, size = 3, fontface = "bold") +
    facet_wrap(~feature, scales = "free_y", labeller = as_labeller(label_map)) +
    scale_colour_manual(values = c("P < 0.05" = "#D73027", "P >= 0.05" = "#333333"), guide = "none") +
    scale_fill_manual(values = c("P < 0.05" = "#F4A6A1", "P >= 0.05" = "grey75"), guide = "none") +
    labs(title = title_text, x = "AUCrel (FLOT; higher = more resistant)", y = y_label) +
    pdo_theme_slide(10) +
    theme(strip.text = element_text(size = 9, face = "bold"))
}

p_auc_state <- plot_auc_linear(auc_state, auc_state_stats, state_axis_labels, "Sample state proportion (%)", "State proportions across the FLOT AUCrel gradient")
p_auc_mp <- plot_auc_linear(auc_mp, auc_mp_stats, mp_labels, "Mean sample UCell score", "Centred-refined MP activity across the FLOT AUCrel gradient")
ggsave(file.path(figures_dir, "Auto_auc_state_linear_patterns.pdf"), p_auc_state, width = 13.33, height = 8, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Auto_auc_mp_linear_patterns.pdf"), p_auc_mp, width = 16, height = 10, useDingbats = FALSE)
####################

####################
# Requested sensitive-versus-resistant figures
####################
plot_two_group <- function(data, stats, label_map, y_label, title_text, percent_axis = FALSE) {
  annotation <- data %>%
    group_by(feature) %>%
    summarise(
      y_max = max(value, na.rm = TRUE),
      y_min = min(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y = y_max + pmax(ifelse(percent_axis, 2, 0.01), 0.1 * pmax(y_max - y_min, ifelse(percent_axis, 5, 0.02)))) %>%
    left_join(stats %>% select(feature, p_value), by = "feature") %>%
    mutate(label = ifelse(is.finite(p_value) & p_value < 0.05, significance_label(p_value), "ns"))
  p <- ggplot(data, aes(feature, value, fill = sensitivity_group)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.62, outlier.shape = NA, alpha = 0.8, colour = "black") +
    stat_summary(aes(group = sensitivity_group), fun = mean, geom = "point", position = position_dodge(width = 0.75), shape = 23, size = 3.4, stroke = 0.8, fill = "white", colour = "black") +
    geom_point(aes(group = sensitivity_group, shape = sur), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75), size = 2.4, colour = "#222222", stroke = 0.8) +
    geom_text(data = annotation, aes(feature, y, label = label), inherit.aes = FALSE, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = comparison_colours, drop = FALSE) +
    scale_x_discrete(labels = label_map) +
    labs(title = title_text, x = NULL, y = y_label, fill = NULL, shape = "PDO") +
    pdo_theme_slide(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
    coord_cartesian(clip = "off")
  if (percent_axis) p <- p + scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0.02, 0.18)))
  else p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)))
  p
}

p_group_state <- plot_two_group(comparison_state, comparison_state_stats, state_axis_labels, "Sample state proportion (%)", "State proportions: FLOT-sensitive versus FLOT-resistant PDOs", TRUE)
p_group_mp <- plot_two_group(comparison_mp, comparison_mp_stats, mp_labels, "Mean sample UCell score", "MP activity: FLOT-sensitive versus FLOT-resistant PDOs")
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_state_proportions.pdf"), p_group_state, width = 13.33, height = 7.5, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_mp_activity.pdf"), p_group_mp, width = 16, height = 8.5, useDingbats = FALSE)

group_stacked <- comparison_state_all %>%
  group_by(sensitivity_group, feature) %>%
  summarise(mean_percent = mean(value), .groups = "drop")
p_group_stacked <- ggplot(group_stacked, aes(sensitivity_group, mean_percent, fill = feature)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = ifelse(mean_percent >= 8, sprintf("%.0f%%", mean_percent), "")),
    position = position_stack(vjust = 0.5), size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = state_colours[state_order_all], labels = state_axis_labels_all[state_order_all], drop = FALSE) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = c(0, 0)) +
  labs(title = "Mean PDO state composition", subtitle = "All states sum to 100%; equal weight per PDO (n=3 per response extreme)", x = NULL, y = "Mean sample proportion", fill = "Cell state") +
  pdo_theme_slide(13)
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_state_composition_stacked.pdf"), p_group_stacked, width = 10, height = 7.5, useDingbats = FALSE)

state_mp_heat <- state_mp_stats %>%
  mutate(
    state = factor(state, levels = state_order),
    feature = factor(feature, levels = rev(mp_order))
  )
mp_limit <- max(abs(state_mp_heat$mean_difference_resistant_minus_sensitive), na.rm = TRUE)
if (!is.finite(mp_limit) || mp_limit == 0) mp_limit <- 0.05
p_state_mp_heat <- ggplot(state_mp_heat, aes(state, feature, fill = mean_difference_resistant_minus_sensitive)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(is.finite(mean_difference_resistant_minus_sensitive), sprintf("%.3f", mean_difference_resistant_minus_sensitive), "")), size = 2.3) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-mp_limit, mp_limit), oob = squish, name = "Mean difference\nresistant - sensitive") +
  scale_x_discrete(labels = state_axis_labels[state_order]) +
  scale_y_discrete(labels = mp_labels[rev(mp_order)]) +
  labs(title = "State-resolved MP activity", x = NULL, y = NULL) +
  pdo_theme_slide(9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_state_resolved_mp_heatmap.pdf"), p_state_mp_heat, width = 13.33, height = 9, useDingbats = FALSE)

pathway_order <- names(pathway_sets)
pathway_heat <- pathway_stats %>%
  mutate(
    state = factor(state, levels = state_order),
    feature = factor(feature, levels = rev(pathway_order))
  )
pathway_limit <- max(abs(pathway_heat$mean_difference_resistant_minus_sensitive), na.rm = TRUE)
if (!is.finite(pathway_limit) || pathway_limit == 0) pathway_limit <- 0.25
p_pathway_heat <- ggplot(pathway_heat, aes(state, feature, fill = mean_difference_resistant_minus_sensitive)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(aes(size = pmin(sensitive_n, resistant_n)), shape = 21, fill = NA, colour = "#333333", stroke = 0.35) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-pathway_limit, pathway_limit), oob = squish, name = "Mean difference\nresistant - sensitive") +
  scale_size_continuous(range = c(0.5, 2.8), breaks = 1:3, name = "PDOs/group") +
  scale_x_discrete(labels = state_axis_labels[state_order]) +
  labs(title = "State-resolved selected pathway activity", x = NULL, y = NULL) +
  pdo_theme_slide(9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_state_pathway_heatmap.pdf"), p_pathway_heat, width = 13.33, height = 9, useDingbats = FALSE)

composite <- (p_group_state / p_group_mp / p_state_mp_heat / p_pathway_heat) +
  plot_layout(heights = c(0.8, 0.9, 1.05, 1.05)) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_composite.pdf"), composite, width = 16, height = 30, limitsize = FALSE, useDingbats = FALSE)
ggsave(file.path(figures_dir, "Auto_sensitive_resistant_composite.png"), composite, width = 16, height = 30, dpi = 300, bg = "white", limitsize = FALSE)
####################

####################
# Report and auditable run summary
####################
report <- tibble(
  item = c(
    "baseline PDO samples", "numeric AUCrel PDOs", "sensitive PDOs", "resistant PDOs",
    "canonical states", "finalized MPs", "selected pathway signatures", "minimum cells per sample-state"
  ),
  value = c(
    nrow(sample_meta), nrow(sample_meta %>% filter(is.finite(auc_rel))),
    paste(sensitive_sur, collapse = "; "), paste(resistant_sur, collapse = "; "),
    length(state_order), length(mp_order), length(pathway_sets), min_cells_per_pseudobulk
  )
)
data.table::fwrite(report, file.path(reports_dir, "Auto_centred_clinical_auc_response_summary.csv"))

####################
# UMAP for Sensitive vs Resistant PDOs
####################
message("Generating UMAP for sensitive vs resistant PDOs...")
if (!exists("pdos")) {
  message("Loading merged PDO object for UMAP.")
  pdos <- readRDS(input_paths[["pdo"]])
}

umap_df <- as.data.frame(Seurat::Embeddings(pdos, "umap"))
umap_df$cell <- rownames(umap_df)
umap_col1 <- colnames(umap_df)[1]
umap_col2 <- colnames(umap_df)[2]

meta_df <- pdos@meta.data
meta_df$cell <- rownames(meta_df)
meta_df$sample <- as.character(meta_df[[PDO_METADATA_COLUMNS$sample]])
meta_df$sur <- stringr::str_extract(meta_df$sample, "^SUR[0-9]+")

umap_meta <- merge(umap_df, meta_df, by = "cell")
umap_filtered <- umap_meta %>%
  filter(
    !grepl("_Treated_PDO$", sample),
    sur %in% c(sensitive_sur, resistant_sur)
  ) %>%
  mutate(
    sensitivity_group = factor(
      ifelse(sur %in% sensitive_sur, "Sensitive", "Resistant"),
      levels = c("Sensitive", "Resistant")
    ),
    sample_display = sur
  )

umap_pdf <- file.path(figures_dir, "Auto_sensitive_resistant_umap.pdf")

if (nrow(umap_filtered) > 0) {
  message("Computing rederived UMAP on the 6 selected samples...")
  pdos_subset <- subset(pdos, cells = umap_filtered$cell)
  pdos_subset <- Seurat::FindVariableFeatures(pdos_subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  pdos_subset <- Seurat::ScaleData(pdos_subset, verbose = FALSE)
  pdos_subset <- Seurat::RunPCA(pdos_subset, npcs = 30, verbose = FALSE)
  pdos_subset <- Seurat::RunUMAP(pdos_subset, dims = 1:30, verbose = FALSE)

  new_umap_df <- as.data.frame(Seurat::Embeddings(pdos_subset, "umap"))
  new_umap_df$cell <- rownames(new_umap_df)
  new_umap_col1 <- colnames(new_umap_df)[1]
  new_umap_col2 <- colnames(new_umap_df)[2]
  new_umap_meta <- merge(new_umap_df, umap_filtered[, c("cell", "sensitivity_group", "sample_display")], by = "cell")

  p_new_sample <- ggplot(new_umap_meta, aes(x = .data[[new_umap_col1]], y = .data[[new_umap_col2]], color = sample_display)) +
    geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
    theme_classic(base_size = 11) +
    labs(title = "By Sample (Rederived UMAP)", color = "Sample") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(legend.position = "bottom")

  p_new_group <- ggplot(new_umap_meta, aes(x = .data[[new_umap_col1]], y = .data[[new_umap_col2]], color = sensitivity_group)) +
    geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
    scale_color_manual(values = c("Sensitive" = "#0072B2", "Resistant" = "#D55E00")) +
    theme_classic(base_size = 11) +
    labs(title = "By Clinical Response (Rederived UMAP)", color = "Response") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(legend.position = "bottom")

  p_new_combined <- p_new_sample + p_new_group + 
    patchwork::plot_annotation(
      title = "UMAP of Sensitive and Resistant PDOs (Rederived)",
      theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
    )

  p_umap_sample <- ggplot(umap_filtered, aes(x = .data[[umap_col1]], y = .data[[umap_col2]], color = sample_display)) +
    geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
    theme_classic(base_size = 11) +
    labs(title = "By Sample (Stored UMAP)", color = "Sample") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(legend.position = "bottom")

  p_umap_group <- ggplot(umap_filtered, aes(x = .data[[umap_col1]], y = .data[[umap_col2]], color = sensitivity_group)) +
    geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
    scale_color_manual(values = c("Sensitive" = "#0072B2", "Resistant" = "#D55E00")) +
    theme_classic(base_size = 11) +
    labs(title = "By Clinical Response (Stored UMAP)", color = "Response") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(legend.position = "bottom")

  p_umap_combined <- p_umap_sample + p_umap_group + 
    patchwork::plot_annotation(
      title = "UMAP of Sensitive and Resistant PDOs (Stored)",
      theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
    )

  grDevices::pdf(umap_pdf, width = 12, height = 6, useDingbats = FALSE)
  print(p_new_combined)
  print(p_umap_combined)
  grDevices::dev.off()
  message("Saved 2-page UMAP plot: ", umap_pdf)
}

expected_outputs <- c(
  cache_path,
  file.path(figures_dir, "Auto_categorical_state_composition_stacked.pdf"),
  file.path(figures_dir, "Auto_categorical_state_association_boxplots.pdf"),
  file.path(figures_dir, "Auto_categorical_mp_association_boxplots.pdf"),
  file.path(figures_dir, "Auto_auc_state_linear_patterns.pdf"),
  file.path(figures_dir, "Auto_auc_mp_linear_patterns.pdf"),
  file.path(figures_dir, "Auto_sensitive_resistant_composite.pdf"),
  file.path(figures_dir, "Auto_sensitive_resistant_composite.png")
)

if (nrow(umap_filtered) > 0) {
  expected_outputs <- c(expected_outputs, umap_pdf)
}

pdo_require_files(expected_outputs)

pdo_write_run_summary(
  script = "analysis/clinical/Auto_centred_clinical_auc_response_associations.R",
  out_dir = out_dir,
  inputs = input_paths,
  outputs = expected_outputs,
  parameters = list(
    baseline_samples = nrow(sample_meta),
    auc_samples = sum(is.finite(sample_meta$auc_rel)),
    sensitive_sur = paste(sensitive_sur, collapse = ","),
    resistant_sur = paste(resistant_sur, collapse = ","),
    min_cells_per_pseudobulk = min_cells_per_pseudobulk,
    elapsed_minutes = round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)
  ),
  cache = cache_policy
)
message("Completed current-centred clinical/AUC response workflow: ", out_dir)
####################
