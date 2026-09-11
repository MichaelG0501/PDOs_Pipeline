####################
# Analysis registry:
#   Status: active terminal successor
#   Script: analysis/cell_states/centred_refined_pdo_flot_matched_response.R
#   Methodology: analysis/methodology/cell_states/centred_refined_pdo_flot_matched_response_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Matched untreated/FLOT PDO response using the finalized
#     centred-refined MPs and noreg states.
# Inputs (all persistent live files):
#   PDOs_outs/PDOs_merged.rds
#   PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds
#   PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds
#   PDOs_outs/centred_mp_refinement/centred_refined_noreg_mp_adj.rds
#   PDOs_outs/centred_mp_refinement/centred_refined_noreg_group_max.rds
#   PDOs_outs/centred_mp_refinement/centred_refined_mp_strict_order.rds
# Outputs:
#   intermediate/: compact matched-cell, pseudobulk, gene-set, and result RDS files
#   tables/: figure source data, paired effects, pathway deltas, and state-resolved DEG tables
#   figures/: Nature-width vector PDF and 600-dpi PNG panels
#   logs/: run summary and session information
#   reports/: compact workflow summary table
# Cache/replot: PDO_FORCE_REBUILD=TRUE recomputes; PDO_REPLOT_ONLY=TRUE requires the live result cache.
# Run: qsub analysis/cell_states/centred_refined_pdo_flot_matched_response.sh
# Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################

####################
suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
  library("edgeR")
  library("dplyr")
  library("tidyr")
  library("tibble")
  library("stringr")
  library("ggplot2")
  library("patchwork")
  library("ComplexHeatmap")
  library("circlize")
  library("data.table")
  library("msigdbr")
  library("grid")
})
####################

####################
# Persistent paths and finalized model constants
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
input_dir <- file.path(project_dir, "PDOs_outs")
model_dir <- file.path(input_dir, "centred_mp_refinement")
out_dir <- file.path(input_dir, "centred_refined_flot_matched_response")
intermediate_dir <- file.path(out_dir, "intermediate")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
logs_dir <- file.path(out_dir, "logs")
summary_dir <- file.path(out_dir, "reports")
invisible(lapply(
  c(intermediate_dir, tables_dir, figures_dir, logs_dir, summary_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

input_paths <- c(
  pdo = file.path(input_dir, "PDOs_merged.rds"),
  ucell = file.path(model_dir, "merged_refined_ucell_scores.rds"),
  states = file.path(model_dir, "centred_refined_noreg_states.rds"),
  mp_adj = file.path(model_dir, "centred_refined_noreg_mp_adj.rds"),
  group_max = file.path(model_dir, "centred_refined_noreg_group_max.rds"),
  strict_order = file.path(model_dir, "centred_refined_mp_strict_order.rds")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Missing required persistent live inputs: ", paste(missing_inputs, collapse = "; "))
}

matched_samples <- c(
  "SUR1070_Treated_PDO", "SUR1070_Untreated_PDO",
  "SUR1090_Treated_PDO", "SUR1090_Untreated_PDO",
  "SUR1072_Treated_PDO", "SUR1072_Untreated_PDO",
  "SUR1181_Treated_PDO", "SUR1181_Untreated_PDO"
)
patient_order <- c("SUR1070", "SUR1090", "SUR1072", "SUR1181")
patient_labels <- c(
  SUR1070 = "SUR1070 (pre-responder)",
  SUR1090 = "SUR1090 (post-responder)",
  SUR1072 = "SUR1072 (post-non-responder)",
  SUR1181 = "SUR1181 (pre-non-responder)"
)
patient_colours <- c(
  "SUR1070 (pre-responder)" = "#0072B2",
  "SUR1090 (post-responder)" = "#56B4E9",
  "SUR1072 (post-non-responder)" = "#D55E00",
  "SUR1181 (pre-non-responder)" = "#E69F00"
)
state_order <- PDO_STATE_ORDER
state_order_all <- PDO_STATE_ORDER_WITH_OPTIONAL
state_colours <- PDO_STATE_COLORS
treatment_colours <- c(Untreated = "#B8B8B8", Treated = "#B43C3C")

strict_object <- readRDS(input_paths[["strict_order"]])
mp_order <- strict_object$strict_refined_mp_order
mp_to_state <- strict_object$mp_to_state[mp_order]
mp_colours <- unname(strict_object$state_cols[mp_to_state])
names(mp_colours) <- mp_order

mp_desc_map <- PDO_MP_DESCRIPTIONS
mp_axis_labels <- setNames(
  ifelse(!is.na(mp_desc_map[mp_order]),
         paste0(mp_order, "\n", mp_desc_map[mp_order]),
         mp_order),
  mp_order
)

composite_pathways <- list(
  "Cell-cycle / proliferation" = c("E2F targets", "G2M checkpoint"),
  "Injury / checkpoint" = c("Apoptosis", "p53 pathway", "DNA repair"),
  "Adaptive / persistence" = c("TNF/NF-kB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response"),
  "Proteostasis / transition" = c("Unfolded protein response", "Oxidative phosphorylation"),
  "Lineage" = c("Intestinal metaplasia", "Ciliated progenitor epithelium"),
  "CCSIG" = c("CCSIG")
)

min_cells_per_pseudobulk <- as.integer(Sys.getenv("PDO_FLOT_MIN_CELLS", "20"))
min_pairs_for_deg <- as.integer(Sys.getenv("PDO_FLOT_MIN_PAIRS", "3"))
force_rebuild <- identical(Sys.getenv("PDO_FORCE_REBUILD"), "TRUE")
replot_only <- identical(Sys.getenv("PDO_REPLOT_ONLY"), "TRUE")
paper_font_family <- Sys.getenv("OAC_PAPER_FONT", "sans")
start_time <- Sys.time()
cache_path <- file.path(intermediate_dir, "centred_refined_pdo_flot_matched_results.rds")

####################
# Helpers
####################
nature_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = paper_font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.3),
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      panel.grid = element_blank()
    )
}

mean_gene_set <- function(expression_matrix, genes) {
  genes_use <- intersect(unique(genes), rownames(expression_matrix))
  if (length(genes_use) == 0L) return(rep(NA_real_, ncol(expression_matrix)))
  Matrix::colMeans(expression_matrix[genes_use, , drop = FALSE])
}

paired_wilcox_summary <- function(data, feature_col, value_col) {
  feature_name <- rlang::ensym(feature_col)
  value_name <- rlang::ensym(value_col)
  data %>%
    group_by(!!feature_name) %>%
    summarise(
      paired_n = sum(is.finite(!!value_name)),
      median_delta = median(!!value_name, na.rm = TRUE),
      min_delta = min(!!value_name, na.rm = TRUE),
      max_delta = max(!!value_name, na.rm = TRUE),
      p_value = if (sum(is.finite(!!value_name)) >= 3L) {
        wilcox.test((!!value_name)[is.finite(!!value_name)], mu = 0, exact = TRUE)$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(fdr = p.adjust(p_value, method = "BH"))
}

state_axis_labels <- c(
  "Classic proliferation" = "Classic\nproliferation",
  "Columnar-to-intestinal" = "Columnar-to-\nintestinal",
  "Glandular differentiation" = "Glandular\ndifferentiation",
  "Stress-adaptive" = "Stress-\nadaptive",
  "ECM-remodelling" = "ECM-\nremodelling",
  "Motile-cilia differentiation" = "Motile-cilia\ndifferentiation"
)

####################
# Current-centred computation and persistent source-data cache
####################
if (!file.exists(cache_path) || force_rebuild) {
  if (replot_only) {
    stop("PDO_REPLOT_ONLY=TRUE but the persistent live cache is missing: ", cache_path)
  }

  message("Loading PDO object and finalized centred-refined states ...")
  pdos <- readRDS(input_paths[["pdo"]])
  state_vector <- readRDS(input_paths[["states"]])
  if (is.null(names(state_vector))) stop("Finalized state vector is not cell-named.")
  missing_state_cells <- setdiff(Cells(pdos), names(state_vector))
  if (length(missing_state_cells) > 0L) {
    stop("Finalized state vector is missing ", length(missing_state_cells), " PDO cells.")
  }
  pdos$centred_refined_state <- state_vector[Cells(pdos)]
  sample_present <- sort(unique(as.character(pdos$orig.ident)))
  absent_samples <- setdiff(matched_samples, sample_present)
  if (length(absent_samples) > 0L) {
    stop("Matched-FLOT samples absent from PDO object: ", paste(absent_samples, collapse = "; "))
  }

  matched <- subset(pdos, subset = orig.ident %in% matched_samples)
  matched$Patient <- str_extract(as.character(matched$orig.ident), "^SUR[0-9]+")
  matched$Treatment <- ifelse(grepl("_Treated_", matched$orig.ident), "Treated", "Untreated")
  matched$Patient_label <- unname(patient_labels[matched$Patient])

  matched_meta <- matched@meta.data %>%
    rownames_to_column("cell") %>%
    transmute(
      cell,
      orig.ident = as.character(orig.ident),
      Patient = factor(Patient, levels = patient_order),
      Patient_label = factor(Patient_label, levels = unname(patient_labels[patient_order])),
      Treatment = factor(Treatment, levels = c("Untreated", "Treated")),
      state = factor(centred_refined_state, levels = state_order_all)
    )
  if (anyNA(matched_meta$state)) stop("Matched cells contain state labels outside the finalized ontology.")

  reduction_names <- Reductions(matched)
  umap_name <- if ("umap" %in% reduction_names) "umap" else reduction_names[grepl("umap", reduction_names, ignore.case = TRUE)][1]
  if (is.na(umap_name) || length(umap_name) == 0L) stop("No UMAP reduction found in PDOs_merged.rds.")
  umap_matrix <- Embeddings(matched, reduction = umap_name)
  umap_source <- data.frame(
    cell = rownames(umap_matrix), UMAP_1 = umap_matrix[, 1], UMAP_2 = umap_matrix[, 2],
    stringsAsFactors = FALSE
  ) %>% left_join(matched_meta, by = "cell")

  message("Computing paired state-composition effects ...")
  state_counts <- matched_meta %>%
    count(Patient, Patient_label, Treatment, orig.ident, state, name = "cell_n") %>%
    complete(
      nesting(Patient, Patient_label, Treatment, orig.ident),
      state = factor(state_order_all, levels = state_order_all),
      fill = list(cell_n = 0L)
    ) %>%
    group_by(Patient, Patient_label, Treatment, orig.ident) %>%
    mutate(total_cells = sum(cell_n), proportion = cell_n / total_cells, percentage = 100 * proportion) %>%
    ungroup()

  state_pairs <- state_counts %>%
    filter(as.character(state) %in% state_order) %>%
    mutate(state = factor(as.character(state), levels = state_order)) %>%
    select(Patient, Patient_label, Treatment, state, cell_n, total_cells, proportion, percentage) %>%
    pivot_wider(names_from = Treatment, values_from = c(cell_n, total_cells, proportion, percentage)) %>%
    mutate(
      delta_percentage_points = percentage_Treated - percentage_Untreated,
      log2_odds_ratio = log2(
        ((cell_n_Treated + 0.5) / (total_cells_Treated - cell_n_Treated + 0.5)) /
          ((cell_n_Untreated + 0.5) / (total_cells_Untreated - cell_n_Untreated + 0.5))
      )
    )
  state_effect_summary <- paired_wilcox_summary(state_pairs, state, delta_percentage_points)

  message("Computing current MP-score effects ...")
  ucell <- readRDS(input_paths[["ucell"]])
  matched_cells <- intersect(matched_meta$cell, rownames(ucell))
  if (length(matched_cells) != nrow(matched_meta)) {
    stop("Current UCell matrix does not cover every matched PDO cell.")
  }
  absent_mps <- setdiff(mp_order, colnames(ucell))
  if (length(absent_mps) > 0L) stop("Current UCell matrix lacks MPs: ", paste(absent_mps, collapse = "; "))
  mp_cell_source <- as.data.frame(ucell[matched_cells, mp_order, drop = FALSE]) %>%
    rownames_to_column("cell") %>%
    left_join(matched_meta %>% select(cell, Patient, Patient_label, Treatment), by = "cell")
  mp_means <- mp_cell_source %>%
    pivot_longer(cols = all_of(mp_order), names_to = "mp", values_to = "ucell_score") %>%
    group_by(Patient, Patient_label, Treatment, mp) %>%
    summarise(mean_ucell = mean(ucell_score, na.rm = TRUE), cell_n = n(), .groups = "drop")
  mp_pairs <- mp_means %>%
    pivot_wider(names_from = Treatment, values_from = c(mean_ucell, cell_n)) %>%
    mutate(
      delta_ucell = mean_ucell_Treated - mean_ucell_Untreated,
      mp = factor(mp, levels = mp_order),
      mp_group = factor(unname(mp_to_state[as.character(mp)]), levels = names(strict_object$state_groups))
    )
  mp_effect_summary <- paired_wilcox_summary(mp_pairs, mp, delta_ucell)

  message("Computing state-resolved pseudobulk Hallmark response ...")
  hallmark <- msigdbr(species = "Homo sapiens", category = "H")
  pathway_ids <- c(
    HALLMARK_E2F_TARGETS = "E2F targets",
    HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
    HALLMARK_APOPTOSIS = "Apoptosis",
    HALLMARK_P53_PATHWAY = "p53 pathway",
    HALLMARK_DNA_REPAIR = "DNA repair",
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF/NF-kB",
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
    HALLMARK_HYPOXIA = "Hypoxia",
    HALLMARK_UNFOLDED_PROTEIN_RESPONSE = "Unfolded protein response",
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxidative phosphorylation",
    HALLMARK_XENOBIOTIC_METABOLISM = "Xenobiotic metabolism"
  )
  pathway_sets <- lapply(names(pathway_ids), function(pathway_id) {
    unique(hallmark$gene_symbol[hallmark$gs_name == pathway_id])
  })
  names(pathway_sets) <- unname(pathway_ids)
  pathway_sets[["Interferon response"]] <- unique(hallmark$gene_symbol[
    hallmark$gs_name %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
  ])
  pathway_order <- names(pathway_sets)

  main_meta <- matched_meta %>%
    filter(as.character(state) %in% state_order) %>%
    mutate(
      state = factor(as.character(state), levels = state_order),
      sample_state = paste(orig.ident, as.character(state), sep = "__")
    )
  counts_matrix <- GetAssayData(matched, assay = "RNA", layer = "counts")
  grouping <- factor(main_meta$sample_state, levels = unique(main_meta$sample_state))
  design_sparse <- Matrix::sparse.model.matrix(~0 + grouping)
  colnames(design_sparse) <- levels(grouping)
  pseudobulk_counts <- counts_matrix[, main_meta$cell, drop = FALSE] %*% design_sparse
  pseudobulk_meta <- main_meta %>%
    count(sample_state, orig.ident, Patient, Patient_label, Treatment, state, name = "cell_n") %>%
    arrange(state, Patient, Treatment)
  pseudobulk_counts <- pseudobulk_counts[, pseudobulk_meta$sample_state, drop = FALSE]

  pseudobulk_dge <- DGEList(counts = pseudobulk_counts)
  pseudobulk_dge <- calcNormFactors(pseudobulk_dge)
  pseudobulk_logcpm <- cpm(pseudobulk_dge, log = TRUE, prior.count = 2)
  pseudobulk_z <- t(scale(t(pseudobulk_logcpm)))
  pseudobulk_z[!is.finite(pseudobulk_z)] <- 0
  pathway_matrix <- vapply(pathway_sets, function(genes) {
    genes_use <- intersect(genes, rownames(pseudobulk_z))
    if (length(genes_use) == 0L) return(rep(NA_real_, ncol(pseudobulk_z)))
    colMeans(pseudobulk_z[genes_use, , drop = FALSE])
  }, numeric(ncol(pseudobulk_z)))
  rownames(pathway_matrix) <- colnames(pseudobulk_z)
  pathway_scores <- as.data.frame(pathway_matrix) %>%
    mutate(sample_state = rownames(.)) %>%
    left_join(pseudobulk_meta, by = "sample_state") %>%
    pivot_longer(cols = all_of(pathway_order), names_to = "pathway", values_to = "score")
  pathway_pairs <- pathway_scores %>%
    filter(cell_n >= min_cells_per_pseudobulk) %>%
    select(Patient, Patient_label, state, pathway, Treatment, score, cell_n) %>%
    pivot_wider(names_from = Treatment, values_from = c(score, cell_n)) %>%
    filter(is.finite(score_Untreated), is.finite(score_Treated)) %>%
    mutate(delta_score = score_Treated - score_Untreated)
  pathway_summary <- pathway_pairs %>%
    group_by(state, pathway) %>%
    summarise(
      paired_n = n(), median_delta = median(delta_score),
      p_value = if (n() >= 3L) wilcox.test(delta_score, mu = 0, exact = TRUE)$p.value else NA_real_,
      .groups = "drop"
    ) %>%
    group_by(state) %>% mutate(fdr_within_state = p.adjust(p_value, method = "BH")) %>% ungroup()

  message("Running state-resolved paired edgeR models ...")
  deg_results <- list()
  deg_summary <- list()
  for (state_name in state_order) {
    state_meta <- pseudobulk_meta %>%
      filter(as.character(state) == state_name, cell_n >= min_cells_per_pseudobulk) %>%
      group_by(Patient) %>% filter(n_distinct(Treatment) == 2L) %>% ungroup() %>%
      mutate(
        Patient = droplevels(factor(as.character(Patient), levels = patient_order)),
        Treatment = factor(Treatment, levels = c("Untreated", "Treated"))
      ) %>%
      arrange(Patient, Treatment)
    pair_n <- n_distinct(state_meta$Patient)
    if (pair_n < min_pairs_for_deg) {
      deg_summary[[state_name]] <- data.frame(
        state = state_name, paired_n = pair_n, tested_genes = 0L,
        fdr_0_05 = 0L, fdr_0_10 = 0L, status = "insufficient paired pseudobulks"
      )
      next
    }
    y <- DGEList(pseudobulk_counts[, state_meta$sample_state, drop = FALSE])
    design <- model.matrix(~Patient + Treatment, data = state_meta)
    keep <- filterByExpr(y, design = design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    test <- glmQLFTest(fit, coef = "TreatmentTreated")
    result <- topTags(test, n = Inf, sort.by = "PValue")$table %>%
      rownames_to_column("gene") %>%
      mutate(state = state_name, paired_n = pair_n)
    deg_results[[state_name]] <- result
    deg_summary[[state_name]] <- data.frame(
      state = state_name, paired_n = pair_n, tested_genes = nrow(result),
      fdr_0_05 = sum(result$FDR < 0.05), fdr_0_10 = sum(result$FDR < 0.10), status = "completed"
    )
  }
  deg_all <- bind_rows(deg_results)
  deg_summary <- bind_rows(deg_summary)

  result_cache <- list(
    model_version = "centred refined MPs; noreg state definition",
    input_paths = input_paths,
    matched_meta = matched_meta,
    umap_source = umap_source,
    state_counts = state_counts,
    state_pairs = state_pairs,
    state_effect_summary = state_effect_summary,
    mp_means = mp_means,
    mp_pairs = mp_pairs,
    mp_effect_summary = mp_effect_summary,
    pathway_sets = pathway_sets,
    pseudobulk_meta = pseudobulk_meta,
    pseudobulk_logcpm = pseudobulk_logcpm,
    pathway_scores = pathway_scores,
    pathway_pairs = pathway_pairs,
    pathway_summary = pathway_summary,
    deg_all = deg_all,
    deg_summary = deg_summary,
    state_order = state_order,
    state_order_all = state_order_all,
    mp_order = mp_order
  )
  saveRDS(result_cache, cache_path)
  saveRDS(pathway_sets, file.path(intermediate_dir, "hallmark_gene_sets_used.rds"))
  saveRDS(pseudobulk_logcpm, file.path(intermediate_dir, "matched_state_pseudobulk_logcpm.rds"))

  data.table::fwrite(umap_source, file.path(tables_dir, "figure5b_matched_umap_source_data.csv"))
  data.table::fwrite(state_counts, file.path(tables_dir, "figure5c_state_composition_source_data.csv"))
  data.table::fwrite(state_pairs, file.path(tables_dir, "figure5c_state_paired_effects.csv"))
  data.table::fwrite(state_effect_summary, file.path(tables_dir, "figure5c_state_effect_summary.csv"))
  data.table::fwrite(mp_pairs, file.path(tables_dir, "figure5c_mp_paired_effects.csv"))
  data.table::fwrite(mp_effect_summary, file.path(tables_dir, "figure5c_mp_effect_summary.csv"))
  data.table::fwrite(pathway_pairs, file.path(tables_dir, "figure5d_pathway_paired_effects.csv"))
  data.table::fwrite(pathway_summary, file.path(tables_dir, "figure5d_pathway_effect_summary.csv"))
  data.table::fwrite(pseudobulk_meta, file.path(tables_dir, "matched_state_pseudobulk_metadata.csv"))
  data.table::fwrite(deg_summary, file.path(tables_dir, "matched_state_deg_summary.csv"))
  if (nrow(deg_all) > 0L) {
    data.table::fwrite(deg_all, file.path(tables_dir, "matched_state_deg_all.csv"))
  }
} else {
  message("Reusing persistent live matched-FLOT result cache ...")
  result_cache <- readRDS(cache_path)
}

list2env(result_cache, envir = environment())

####################
# Nature-width vector figure layer (replot-only from live cache)
####################
umap_plot <- umap_source %>%
  mutate(
    state = factor(as.character(state), levels = state_order_all),
    Treatment = factor(Treatment, levels = c("Untreated", "Treated"))
  )
p_a <- ggplot(umap_plot, aes(UMAP_1, UMAP_2, colour = state)) +
  geom_point(size = 0.05, alpha = 0.6, stroke = 0) +
  facet_wrap(~Treatment, nrow = 1) +
  scale_colour_manual(values = state_colours, drop = FALSE) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  nature_theme(6.5) +
  theme(legend.position = "bottom", legend.key.width = unit(6, "pt")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

composition_plot_data <- state_counts %>%
  mutate(
    state = factor(as.character(state), levels = rev(state_order_all)),
    pair_label = factor(
      paste0(as.character(Patient), "\n", ifelse(Treatment == "Treated", "FLOT", "Untreated")),
      levels = unlist(lapply(patient_order, function(patient) c(
        paste0(patient, "\nUntreated"), paste0(patient, "\nFLOT")
      )))
    )
  )
p_b <- ggplot(composition_plot_data, aes(pair_label, percentage, fill = state)) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = state_colours, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cells (%)", fill = NULL) +
  nature_theme(6.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

state_pair_plot <- state_pairs %>%
  mutate(
    state = factor(as.character(state), levels = state_order),
    state_x = as.numeric(state)
  )
state_summary_plot <- state_effect_summary %>%
  mutate(state = factor(as.character(state), levels = state_order), state_x = as.numeric(state))
p_c <- ggplot(state_pair_plot, aes(state_x, delta_percentage_points, colour = Patient_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_segment(
    data = state_summary_plot,
    aes(x = state_x - 0.27, xend = state_x + 0.27, y = median_delta, yend = median_delta),
    inherit.aes = FALSE, linewidth = 0.8, colour = "black"
  ) +
  geom_point(size = 1.8, position = position_jitter(width = 0.06, height = 0)) +
  scale_colour_manual(values = patient_colours) +
  scale_x_continuous(breaks = seq_along(state_order), labels = state_axis_labels[state_order]) +
  labs(x = NULL, y = expression(Delta*" cells (percentage points)"), colour = NULL) +
  nature_theme(6.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")

mp_heatmap_data <- mp_pairs %>%
  mutate(Patient = factor(Patient, levels = patient_order), mp = factor(as.character(mp), levels = mp_order)) %>%
  select(Patient, mp, delta_ucell) %>%
  pivot_wider(names_from = mp, values_from = delta_ucell) %>%
  column_to_rownames("Patient") %>%
  as.matrix()
mp_heatmap_data <- mp_heatmap_data[, intersect(mp_order, colnames(mp_heatmap_data))]
mp_limit <- max(abs(mp_heatmap_data), na.rm = TRUE)
if (!is.finite(mp_limit) || mp_limit == 0) mp_limit <- 0.05
mp_heatmap_col_labels <- ifelse(
  !is.na(mp_desc_map[colnames(mp_heatmap_data)]),
  paste0(colnames(mp_heatmap_data), ": ", mp_desc_map[colnames(mp_heatmap_data)]),
  colnames(mp_heatmap_data)
)
mp_heatmap <- Heatmap(
  mp_heatmap_data,
  name = "Delta UCell",
  col = colorRamp2(c(-mp_limit, 0, mp_limit), c("#2C6AA0", "white", "#B33C3C")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontfamily = paper_font_family, fontsize = 6.5),
  column_labels = mp_heatmap_col_labels,
  column_names_gp = gpar(fontfamily = paper_font_family, fontsize = 6.2),
  column_names_rot = 45,
  top_annotation = HeatmapAnnotation(
    MP_group = unname(mp_to_state[colnames(mp_heatmap_data)]),
    col = list(MP_group = strict_object$state_cols),
    show_annotation_name = FALSE,
    annotation_legend_param = list(labels_gp = gpar(fontfamily = paper_font_family, fontsize = 6.2))
  ),
  heatmap_legend_param = list(
    title_gp = gpar(fontfamily = paper_font_family, fontsize = 6.5),
    labels_gp = gpar(fontfamily = paper_font_family, fontsize = 6.2)
  )
)
mp_grob <- grid::grid.grabExpr(draw(mp_heatmap, merge_legend = TRUE))

pathway_plot_data <- pathway_summary %>%
  mutate(
    state = factor(as.character(state), levels = state_order),
    pathway = factor(pathway, levels = rev(unique(pathway)))
  )
pathway_limit <- max(abs(pathway_plot_data$median_delta), na.rm = TRUE)
if (!is.finite(pathway_limit) || pathway_limit == 0) pathway_limit <- 0.25
p_e <- ggplot(pathway_plot_data, aes(state, pathway, fill = median_delta)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(aes(size = paired_n), shape = 21, fill = NA, colour = "#333333", stroke = 0.35) +
  scale_fill_gradient2(
    low = "#2C6AA0", mid = "white", high = "#B33C3C", midpoint = 0,
    limits = c(-pathway_limit, pathway_limit), oob = scales::squish,
    name = expression("Median "*Delta*" score")
  ) +
  scale_size_continuous(range = c(0.4, 2.3), breaks = 3:4, name = "Paired PDOs") +
  scale_x_discrete(labels = state_axis_labels[state_order]) +
  labs(x = NULL, y = NULL) +
  nature_theme(6.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

pdf_files <- c(
  figure5b = file.path(figures_dir, "figure5b_matched_state_umap.pdf"),
  figure5c_composition = file.path(figures_dir, "figure5c_state_composition.pdf"),
  figure5c_state = file.path(figures_dir, "figure5c_state_paired_change.pdf"),
  figure5c_mp = file.path(figures_dir, "figure5c_mp_delta_heatmap.pdf"),
  figure5d = file.path(figures_dir, "figure5d_pathway_delta_heatmap.pdf"),
  composite = file.path(figures_dir, "figure5_centred_refined_flot_response.pdf")
)
ggsave(pdf_files[["figure5b"]], p_a, device = cairo_pdf, width = 183, height = 82, units = "mm")
ggsave(pdf_files[["figure5c_composition"]], p_b, device = cairo_pdf, width = 120, height = 82, units = "mm")
ggsave(pdf_files[["figure5c_state"]], p_c, device = cairo_pdf, width = 120, height = 82, units = "mm")
cairo_pdf(pdf_files[["figure5c_mp"]], width = 183 / 25.4, height = 72 / 25.4, family = paper_font_family)
draw(mp_heatmap, merge_legend = TRUE)
dev.off()
ggsave(pdf_files[["figure5d"]], p_e, device = cairo_pdf, width = 183, height = 90, units = "mm")

composite <- (p_a / (p_b | p_c) / patchwork::wrap_elements(full = mp_grob) / p_e) +
  plot_layout(heights = c(1.05, 1, 0.9, 1.25)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(family = paper_font_family, face = "bold", size = 8)))
ggsave(pdf_files[["composite"]], composite, device = cairo_pdf, width = 183, height = 285, units = "mm", limitsize = FALSE)
ggsave(
  sub("pdf$", "png", pdf_files[["composite"]]), composite,
  width = 183, height = 285, units = "mm", dpi = 600, bg = "white", limitsize = FALSE
)

####################
# Canonical Figure 5 panel source-data names (the first compute run retained
# earlier provisional names; these persistent copies match the manuscript).
####################
data.table::fwrite(umap_source, file.path(tables_dir, "figure5b_matched_umap_source_data.csv"))
data.table::fwrite(state_counts, file.path(tables_dir, "figure5c_state_composition_source_data.csv"))
data.table::fwrite(state_pairs, file.path(tables_dir, "figure5c_state_paired_effects.csv"))
data.table::fwrite(state_effect_summary, file.path(tables_dir, "figure5c_state_effect_summary.csv"))
data.table::fwrite(mp_pairs, file.path(tables_dir, "figure5c_mp_paired_effects.csv"))
data.table::fwrite(mp_effect_summary, file.path(tables_dir, "figure5c_mp_effect_summary.csv"))
data.table::fwrite(pathway_pairs, file.path(tables_dir, "figure5d_pathway_paired_effects.csv"))
data.table::fwrite(pathway_summary, file.path(tables_dir, "figure5d_pathway_effect_summary.csv"))

compact_summary <- state_effect_summary %>%
  transmute(
    analysis = "matched state abundance",
    feature = as.character(state), paired_n, median_delta,
    p_value, fdr,
    model_version = "centred refined noreg"
  ) %>%
  bind_rows(mp_effect_summary %>% transmute(
    analysis = "matched MP UCell",
    feature = as.character(mp), paired_n, median_delta,
    p_value, fdr,
    model_version = "centred refined noreg"
  ))
data.table::fwrite(compact_summary, file.path(summary_dir, "centred_refined_pdo_flot_matched_summary.csv"))

####################
# Extended / Legacy Plots (Boxplots and Node Plots)
####################
message("=== Extended / Legacy Plots ===")

# --- Helper functions ---
sig_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

draw_nodeplot <- function(ndf, edf, title_text, max_node, max_edge) {
  p <- ggplot() +
    { if (nrow(edf) > 0) geom_segment(data = edf, aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct), color = "grey35", alpha = 0.8) } +
    geom_point(data = ndf, aes(x = x, y = y, size = pct, color = state)) +
    geom_text(data = ndf, aes(x = label_x, y = label_y, label = paste0(state, "\n", sprintf("%.1f%%", pct))), size = 3.5, fontface = "bold") +
    { if (nrow(edf) > 0) geom_label(data = edf, aes(x = (x+xend)/2, y = (y+yend)/2, label = sprintf("%.1f%%", pct)), size = 2.6, fill = "white", label.size = 0, fontface = "bold") } +
    scale_color_manual(values = state_colours) +
    scale_size(limits = c(0, max_node), range = c(8, 30), guide = "none") +
    scale_linewidth(limits = c(0, max_edge), range = c(0.6, 12), guide = "none") +
    coord_equal() + expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5)) +
    theme_void(base_size = 14) +
    labs(title = title_text) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  p
}

# --- Node Plots ---
message("Generating Node Plots...")
real_states <- state_order
n_rs <- length(real_states)
theta <- seq(0, 2*pi, length.out = n_rs + 1)[1:n_rs]
layout_df <- data.frame(state = real_states, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)

# Load group_max for hybrid assignment
group_max_mat <- readRDS(input_paths[["group_max"]])
matched_hybrid_cells <- intersect(matched_meta$cell[matched_meta$state == "Hybrid"], rownames(group_max_mat))

hybrid_subtypes <- apply(group_max_mat[matched_hybrid_cells, real_states, drop = FALSE], 1, function(scores) {
  top2 <- names(sort(scores, decreasing = TRUE))[1:2]
  top2_ordered <- top2[order(match(top2, real_states))]
  paste(top2_ordered, collapse = "\n-\n")
})
hybrid_pairs_list <- combn(real_states, 2, simplify = FALSE)
hybrid_levels_ordered <- vapply(hybrid_pairs_list, function(x) paste(x, collapse = "\n-\n"), character(1))

build_node_data <- function(patient_id, trt) {
  samp <- paste0(patient_id, "_", trt, "_PDO")
  cells_in <- matched_meta$cell[matched_meta$orig.ident == samp & matched_meta$state %in% c(real_states, "Hybrid")]
  if (length(cells_in) == 0) return(list(node_df = NULL, edge_df = NULL))
  
  st <- as.character(matched_meta$state[match(cells_in, matched_meta$cell)])
  names(st) <- cells_in
  tot <- length(cells_in)
  
  sdf <- data.frame(state = st[st %in% real_states], stringsAsFactors = FALSE) %>%
    count(state, name = "cells") %>% mutate(pct = 100 * cells / tot)
  
  hyb_cells <- intersect(names(st)[st == "Hybrid"], matched_hybrid_cells)
  if (length(hyb_cells) > 0) {
    pair_lab <- hybrid_subtypes[hyb_cells]
    edf <- data.frame(pair = pair_lab, stringsAsFactors = FALSE) %>%
      count(pair, name = "hybrid_cells") %>%
      tidyr::separate(pair, into = c("from","to"), sep = "\n-\n", remove = FALSE) %>%
      mutate(pct = 100 * hybrid_cells / tot)
  } else {
    edf <- data.frame(pair = character(), from = character(), to = character(), hybrid_cells = integer(), pct = numeric())
  }
  
  ndf <- left_join(layout_df, sdf, by = "state")
  ndf$cells[is.na(ndf$cells)] <- 0; ndf$pct[is.na(ndf$pct)] <- 0
  ndf$label_x <- ndf$x * 1.25; ndf$label_y <- ndf$y * 1.25
  edf2 <- edf %>% left_join(layout_df, by = c("from" = "state")) %>%
    left_join(layout_df, by = c("to" = "state"), suffix = c("","_to")) %>%
    rename(xend = x_to, yend = y_to)
  list(node_df = ndf, edge_df = edf2)
}

all_node_data <- list()
for (pid in patient_order) {
  for (trt in c("Untreated", "Treated")) {
    key <- paste(pid, trt, sep = "_")
    all_node_data[[key]] <- build_node_data(pid, trt)
  }
}

combine_node_edge <- function(trt_label) {
  keys <- paste(patient_order, trt_label, sep = "_")
  node_list <- lapply(keys, function(k) all_node_data[[k]]$node_df)
  node_all <- bind_rows(node_list) %>% group_by(state, x, y, label_x, label_y) %>%
    summarise(pct = median(pct, na.rm = TRUE), cells = median(cells, na.rm = TRUE), .groups = "drop")
  edge_list <- lapply(keys, function(k) all_node_data[[k]]$edge_df)
  edge_all <- bind_rows(edge_list)
  if (nrow(edge_all) > 0) {
    edge_all <- edge_all %>% group_by(from, to, x, y, xend, yend) %>%
      summarise(pct = median(pct, na.rm = TRUE), hybrid_cells = median(hybrid_cells, na.rm = TRUE), .groups = "drop")
  }
  list(node_df = node_all, edge_df = edge_all)
}

node_pdf_path <- file.path(figures_dir, "Auto_pdo_flot_nodeplot_untreated_vs_treated.pdf")
pdf(node_pdf_path, width = 14, height = 7)
comb_ut <- combine_node_edge("Untreated"); comb_tr <- combine_node_edge("Treated")
max_node_comb <- max(c(comb_ut$node_df$pct, comb_tr$node_df$pct), na.rm = TRUE)
max_edge_comb <- max(c(comb_ut$edge_df$pct, comb_tr$edge_df$pct, 0.1), na.rm = TRUE)
p_comb <- draw_nodeplot(comb_ut$node_df, comb_ut$edge_df, "Untreated (median)", max_node_comb, max_edge_comb) |
  draw_nodeplot(comb_tr$node_df, comb_tr$edge_df, "FLOT-treated (median)", max_node_comb, max_edge_comb)
print(p_comb + plot_annotation(title = "Combined Across 4 Patients", theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))

for (pid in patient_order) {
  ut <- all_node_data[[paste0(pid, "_Untreated")]]; tr <- all_node_data[[paste0(pid, "_Treated")]]
  if (is.null(ut$node_df) || is.null(tr$node_df)) next
  max_node_pt <- max(c(ut$node_df$pct, tr$node_df$pct), na.rm = TRUE)
  max_edge_pt <- max(c(ut$edge_df$pct, tr$edge_df$pct, 0.1), na.rm = TRUE)
  pp <- draw_nodeplot(ut$node_df, ut$edge_df, "Untreated", max_node_pt, max_edge_pt) |
    draw_nodeplot(tr$node_df, tr$edge_df, "FLOT-treated", max_node_pt, max_edge_pt)
  print(pp + plot_annotation(title = paste0(pid, " - ", patient_labels[pid]),
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))
}
dev.off()

# --- Boxplots ---
message("Generating Boxplots...")

state_abund_long <- state_counts %>% filter(state %in% state_order) %>%
  mutate(state = factor(as.character(state), levels = state_order))
state_sig <- state_abund_long %>% group_by(state) %>%
  summarise(p = tryCatch(wilcox.test(percentage[Treatment == "Untreated"], percentage[Treatment == "Treated"], paired = TRUE)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  mutate(label = sapply(p, function(pval) if(is.na(pval)) "" else paste0("p=", format(round(pval, 3), nsmall=3))))
state_ymax <- state_abund_long %>% group_by(state) %>% summarise(ymax = max(percentage, na.rm = TRUE), .groups = "drop")
state_sig <- left_join(state_sig, state_ymax, by = "state") %>% mutate(y_pos = ymax * 1.08)

hybrid_meta <- matched_meta %>% filter(cell %in% matched_hybrid_cells)
hybrid_meta$hybrid_subtype <- factor(hybrid_subtypes[hybrid_meta$cell], levels = hybrid_levels_ordered)
hybrid_counts <- hybrid_meta %>% count(Patient, Patient_label, Treatment, orig.ident, hybrid_subtype, name = "cell_n") %>%
  complete(nesting(Patient, Patient_label, Treatment, orig.ident), hybrid_subtype = hybrid_levels_ordered, fill = list(cell_n = 0L)) %>%
  group_by(Patient, Patient_label, Treatment, orig.ident) %>% mutate(total_hybrid_cells = sum(cell_n)) %>% ungroup() %>%
  mutate(pct = 100 * cell_n / total_hybrid_cells)

hybrid_sig <- hybrid_counts %>% group_by(hybrid_subtype) %>%
  summarise(p = tryCatch(wilcox.test(pct[Treatment == "Untreated"], pct[Treatment == "Treated"], paired = TRUE)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  mutate(label = sapply(p, function(pval) if(is.na(pval)) "" else paste0("p=", format(round(pval, 3), nsmall=3))))
hybrid_ymax <- hybrid_counts %>% group_by(hybrid_subtype) %>% summarise(ymax = max(pct, na.rm = TRUE), .groups = "drop")
hybrid_sig <- left_join(hybrid_sig, hybrid_ymax, by = "hybrid_subtype") %>% mutate(y_pos = ymax * 1.08)

mp_expr_long <- mp_means %>% mutate(MP = factor(as.character(mp), levels = mp_order))
mp_sig <- mp_expr_long %>% group_by(MP) %>%
  summarise(p = tryCatch(wilcox.test(mean_ucell[Treatment == "Untreated"], mean_ucell[Treatment == "Treated"], paired = TRUE)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  mutate(label = sapply(p, function(pval) if(is.na(pval)) "" else paste0("p=", format(round(pval, 3), nsmall=3))))
mp_ymax <- mp_expr_long %>% group_by(MP) %>% summarise(ymax = max(mean_ucell, na.rm = TRUE), .groups = "drop")
mp_sig <- left_join(mp_sig, mp_ymax, by = "MP") %>% mutate(y_pos = ymax * 1.08)

box_pdf_path <- file.path(figures_dir, "Auto_pdo_flot_paired_boxplots.pdf")
pdf(box_pdf_path, width = 16, height = 9)

p_box_state <- ggplot(state_abund_long, aes(x = state, y = percentage, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_point(aes(color = Patient_label, group = Treatment), position = position_dodge(width = 0.75), size = 2.5, alpha = 0.9) +
  geom_text(data = state_sig %>% filter(label != ""), aes(x = state, y = y_pos, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = treatment_colours) + scale_color_manual(values = patient_colours) +
  scale_x_discrete(labels = state_axis_labels[state_order]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "State Abundance: Untreated vs FLOT-treated", x = NULL, y = "% of malignant cells", fill = NULL, color = "Patient") +
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 35, hjust = 1, size = 12, face = "bold"), legend.position = "top")
print(p_box_state)

p_box_hybrid <- ggplot(hybrid_counts, aes(x = hybrid_subtype, y = pct, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_point(aes(color = Patient_label, group = Treatment), position = position_dodge(width = 0.75), size = 2.5, alpha = 0.9) +
  geom_text(data = hybrid_sig %>% filter(label != ""), aes(x = hybrid_subtype, y = y_pos, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = treatment_colours) + scale_color_manual(values = patient_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "Hybrid Abundance: Untreated vs FLOT-treated", x = NULL, y = "% of hybrid cells", fill = NULL, color = "Patient") +
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", lineheight = 0.8), legend.position = "top")
print(p_box_hybrid)

p_box_mp <- ggplot(mp_expr_long, aes(x = MP, y = mean_ucell, fill = Treatment)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_point(aes(color = Patient_label, group = Treatment), position = position_dodge(width = 0.75), size = 2.5, alpha = 0.9) +
  geom_text(data = mp_sig %>% filter(label != ""), aes(x = MP, y = y_pos, label = label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = treatment_colours) + scale_color_manual(values = patient_colours) +
  scale_x_discrete(labels = mp_axis_labels[mp_order]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(title = "MP Expression: Untreated vs FLOT-treated", x = NULL, y = "Mean UCell score per sample", fill = NULL, color = "Patient") +
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold", size = 18), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), legend.position = "top")
print(p_box_mp)
dev.off()

# --- Improved Pathway Heatmap ---
message("Generating Improved Pathway Heatmap...")

# Recompute pseudobulk z-scores from cached pseudobulk_logcpm
pseudobulk_z <- t(scale(t(pseudobulk_logcpm)))
pseudobulk_z[!is.finite(pseudobulk_z)] <- 0

# CC signature: select top-50 consensus cell-cycle genes by mean pseudobulk logCPM
cc_genes_df <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE, stringsAsFactors = FALSE
)[, 1:3]
cc_consensus <- cc_genes_df$Gene[cc_genes_df$Consensus == 1]
cc_consensus <- intersect(cc_consensus, rownames(pseudobulk_z))
cc_top50 <- names(sort(rowMeans(pseudobulk_logcpm[cc_consensus, , drop = FALSE], na.rm = TRUE), decreasing = TRUE))[1:min(50, length(cc_consensus))]
cc_pb_score <- mean_gene_set(pseudobulk_z, cc_top50)
names(cc_pb_score) <- colnames(pseudobulk_z)

# Lineage MP gene scores from the merged refined MP gene lists
refined_mp_genes <- readRDS(file.path(model_dir, "merged_refined_mp_genes.rds"))
mp15_genes <- refined_mp_genes[["MP15"]]
mp17_genes <- refined_mp_genes[["MP17+"]]

mp15_pb_score <- mean_gene_set(pseudobulk_z, mp15_genes)
mp17_pb_score <- mean_gene_set(pseudobulk_z, mp17_genes)
names(mp15_pb_score) <- colnames(pseudobulk_z)
names(mp17_pb_score) <- colnames(pseudobulk_z)

extra_scores <- data.frame(
  sample_state = colnames(pseudobulk_z),
  CCSIG = cc_pb_score,
  `Intestinal metaplasia` = mp15_pb_score,
  `Ciliated progenitor epithelium` = mp17_pb_score,
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  left_join(pseudobulk_meta %>% select(sample_state, state, Patient, Patient_label, Treatment, cell_n), by = "sample_state") %>%
  filter(cell_n >= min_cells_per_pseudobulk)

extra_delta <- extra_scores %>%
  select(Patient, Patient_label, state, Treatment, CCSIG, `Intestinal metaplasia`, `Ciliated progenitor epithelium`, cell_n) %>%
  pivot_longer(cols = c(CCSIG, `Intestinal metaplasia`, `Ciliated progenitor epithelium`), names_to = "pathway", values_to = "score") %>%
  pivot_wider(names_from = Treatment, values_from = c(score, cell_n)) %>%
  filter(!is.na(score_Untreated), !is.na(score_Treated)) %>%
  mutate(delta_score = score_Treated - score_Untreated) %>%
  select(Patient, Patient_label, state, pathway, delta_score)

comb_delta <- bind_rows(
  pathway_pairs %>% select(Patient, Patient_label, state, pathway, delta_score),
  extra_delta
)

data.table::fwrite(comb_delta, file.path(tables_dir, "improved_pathway_heatmap_comb_delta.csv"))

# --- Composite Pathway Response Plot ---
message("Generating Composite Pathway Response Plot...")
composite_delta <- bind_rows(lapply(names(composite_pathways), function(metric_name) {
  comb_delta %>% filter(pathway %in% composite_pathways[[metric_name]]) %>%
    group_by(Patient, Patient_label, state) %>%
    summarise(metric = metric_name, delta_score = mean(delta_score, na.rm = TRUE), pathway_n = n_distinct(pathway), .groups = "drop")
})) %>% mutate(
  state = factor(as.character(state), levels = state_order),
  metric = factor(metric, levels = names(composite_pathways)),
  Patient = factor(Patient, levels = patient_order),
  Patient_label = factor(Patient_label, levels = unname(patient_labels[patient_order]))
)
data.table::fwrite(composite_delta, file.path(tables_dir, "composite_response_deltas.csv"))
composite_summary <- composite_delta %>% group_by(state, metric) %>%
  summarise(median_delta = median(delta_score, na.rm = TRUE), .groups = "drop")

p_composite <- ggplot(composite_delta, aes(state, delta_score, color = Patient_label, group = Patient_label)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_line(alpha = 0.45, linewidth = 0.45) + geom_point(size = 2.4) +
  geom_point(data = composite_summary, aes(x = state, y = median_delta), inherit.aes = FALSE, color = "black", shape = 18, size = 3.0) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") + scale_color_manual(values = patient_colours) +
  scale_x_discrete(labels = state_axis_labels[state_order]) +
  labs(title = "Per State Pathway Response Score", x = NULL, y = "Delta (Treated vs Untreated)", color = NULL) +
  nature_theme(7) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(figures_dir, "Auto_pdo_flot_composite_pathway_response.pdf"), p_composite, device = cairo_pdf, width = 183, height = 140, units = "mm")
ggsave(file.path(figures_dir, "Auto_pdo_flot_composite_pathway_response.png"), p_composite, width = 183, height = 140, units = "mm", dpi = 600, bg = "white")

# Mean delta
comb_mean <- comb_delta %>%
  group_by(state, pathway) %>%
  summarise(mean_delta = mean(delta_score, na.rm = TRUE), .groups = "drop")

# Pairwise sig (paired test of Treated - Untreated == 0)
comb_sig <- comb_delta %>%
  group_by(state, pathway) %>%
  summarise(
    p = tryCatch(wilcox.test(delta_score, mu = 0)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = sapply(p, sig_label))

# Build matrices (rows = features, cols = states)
# Pathway order grouped by functional block (matching legacy style)
pathway_order_grouped <- c(
  "E2F targets", "G2M checkpoint",
  "Apoptosis", "p53 pathway", "DNA repair",
  "TNF/NF-kB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response",
  "Unfolded protein response", "Oxidative phosphorylation"
)
pw_order <- c(pathway_order_grouped, "CCSIG", "Intestinal metaplasia", "Ciliated progenitor epithelium")

pw_mat <- comb_mean %>%
  pivot_wider(names_from = state, values_from = mean_delta) %>%
  column_to_rownames("pathway") %>%
  as.matrix()
pw_mat <- pw_mat[pw_order, state_order, drop = FALSE]

pw_sig_mat <- comb_sig %>%
  select(state, pathway, label) %>%
  pivot_wider(names_from = state, values_from = label) %>%
  column_to_rownames("pathway") %>% as.matrix()
pw_sig_mat <- pw_sig_mat[pw_order, state_order, drop = FALSE]

pathway_block_df_imp <- data.frame(
  pathway = c(
    "E2F targets", "G2M checkpoint",
    "Apoptosis", "p53 pathway", "DNA repair",
    "TNF/NF-kB", "EMT", "Hypoxia", "Xenobiotic metabolism", "Interferon response",
    "Unfolded protein response", "Oxidative phosphorylation"
  ),
  block = c(
    rep("Cell-cycle", 2), rep("Injury", 3),
    rep("Adaptive", 5), rep("Proteostasis", 2)
  ),
  stringsAsFactors = FALSE
)

row_blocks_pw <- c(pathway_block_df_imp$block[match(pathway_order_grouped, pathway_block_df_imp$pathway)], "CCSIG", "Lineage", "Lineage")
block_cols_ext <- c(
  "Cell-cycle" = "#F0C75E", "Injury" = "#E07B54",
  "Adaptive" = "#5DAA68", "Proteostasis" = "#5B9BD5",
  "CCSIG" = "#6A0572", "Lineage" = "#1A535C"
)

# Color scale (unified for all)
pw_clip <- max(0.3, quantile(abs(pw_mat), 0.95, na.rm = TRUE))
col_pw <- colorRamp2(c(-pw_clip, 0, pw_clip), c("#245F7B", "white", "#B63E2F"))

# Row block annotation for pathways + CCSIG + Lineage
ha_row_pw <- rowAnnotation(Block = row_blocks_pw, col = list(Block = block_cols_ext),
  show_annotation_name = FALSE, show_legend = TRUE,
  annotation_legend_param = list(Block = list(title = "Signature Type")))

# State annotation on top
ha_top_imp <- HeatmapAnnotation(
  State = state_order,
  col = list(State = state_colours[state_order]),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 10)
)

# Build heatmap with significance labels and values
cell_fun_sig <- function(j, i, x, y, w, h, fill) {
  grid.text(sprintf("%.2f", pw_mat[i, j]), x, y - unit(1, "mm"), gp = gpar(fontsize = 7))
  lbl <- pw_sig_mat[i, j]
  if (!is.na(lbl) && lbl != "" && lbl != "ns") {
    grid.text(lbl, x, y + unit(2, "mm"), gp = gpar(fontsize = 9, fontface = "bold", col = "black"))
  }
}

ht_pw <- Heatmap(pw_mat, name = "Mean Delta (Treated - Untreated)",
  col = col_pw, cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  row_split = factor(row_blocks_pw, levels = names(block_cols_ext)),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(3, "mm"),
  left_annotation = ha_row_pw,
  top_annotation = ha_top_imp,
  cell_fun = cell_fun_sig
)

improved_hm_path <- file.path(figures_dir, "Auto_pdo_flot_improved_pathway_heatmap.pdf")
message("Writing: ", improved_hm_path)
pdf(improved_hm_path, width = 12, height = 11)
draw(ht_pw,
  column_title = "Mean Pathway Response Matrix (Treated - Untreated)",
  column_title_gp = gpar(fontface = "bold", fontsize = 14),
  merge_legend = TRUE
)
dev.off()
message("Improved pathway heatmap PDF done.")


end_time <- Sys.time()
writeLines(
  c(
    "centred_refined_pdo_flot_matched_response.R",
    "Model: finalized centred-refined MPs; noreg state definition",
    paste("Start:", format(start_time, tz = "Europe/London")),
    paste("End:", format(end_time, tz = "Europe/London")),
    paste("Force rebuild:", force_rebuild),
    paste("Replot only:", replot_only),
    paste("Matched patients:", paste(patient_order, collapse = "; ")),
    paste("Inputs:", paste(input_paths, collapse = "; ")),
    paste("Cache:", cache_path),
    paste("Figures:", paste(pdf_files, collapse = "; ")),
    capture.output(sessionInfo())
  ),
  file.path(logs_dir, "centred_refined_pdo_flot_matched_run_summary.txt")
)
message("Current-centred matched-FLOT workflow completed: ", out_dir)
####################
