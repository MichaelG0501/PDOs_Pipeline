####################
# Analysis registry:
#   Status: active terminal extension of the clinical/AUC workflow
#   Script: analysis/clinical/Auto_centred_clinical_pseudobulk_de.R
#   Methodology: analysis/methodology/clinical/Auto_centred_clinical_auc_response_associations_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Batch-adjusted resistant-minus-sensitive pseudobulk edgeR
#     quasi-likelihood differential expression and all-Hallmark CAMERA testing;
#     overall all-cell comparison only; two figure-only PDF pages.
#   Inputs: PDOs_outs/PDOs_merged.rds;
#     /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDO_ClinicalMetadata_V3.xlsx;
#     human Hallmark gene sets from installed msigdbr 7.5.1 (MSigDB 7.5.1).
#   Outputs: under PDOs_outs/Auto_centred_clinical_auc_response_associations/Auto_sensitive_resistant_DE/:
#     intermediate/Auto_pseudobulk_de_results.rds
#     tables/Auto_{sample_metadata,cell_counts,model_audit,differential_expression,hallmark_camera,design_matrices,gene_filtering,report_hallmark_plot}.csv
#     figures/Auto_sensitive_resistant_pseudobulk_DE_pathways.pdf (only figure)
#     reports/Auto_summary.txt; logs/Auto_centred_clinical_pseudobulk_de_run_summary.txt
#   Downstream: terminal report and reproducible source data only.
#   Cache/replot: PDO_FORCE_REBUILD=1 rebuilds; PDO_REPLOT_ONLY=1 requires
#     matching inputs, parameters, package versions and existing cache.
#   Run: qsub analysis/clinical/Auto_centred_clinical_pseudobulk_de.sh
#     Called by Auto_centred_clinical_auc_response_associations.R.
#   Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(msigdbr)
  library(ggplot2)
  library(data.table)
})
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")
source(file.path(PDO_PROJECT_DIR, "analysis/shared/Auto_pdo_analysis_helpers.R"))

run_clinical_pseudobulk_de <- function(input_paths, sensitive_sur, resistant_sur,
                                      min_cells = 20L, pdos = NULL) {
  started <- Sys.time()
  if (!is.finite(min_cells) || min_cells < 20L) stop("DE requires at least 20 cells per pseudobulk.")
  out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_centred_clinical_auc_response_associations",
                       "Auto_sensitive_resistant_DE")
  paths <- pdo_ensure_output_tiers(out_dir)
  cache_path <- file.path(paths[["intermediate"]], "Auto_pseudobulk_de_results.rds")
  pdf_path <- file.path(paths[["figures"]], "Auto_sensitive_resistant_pseudobulk_DE_pathways.pdf")
  inputs <- input_paths[c("pdo", "clinical")]
  pdo_require_files(inputs)
  input_info <- file.info(inputs)
  signature <- list(inputs = inputs, size = input_info$size, mtime = as.numeric(input_info$mtime),
                    sensitive = sensitive_sur, resistant = resistant_sur, min_cells = min_cells,
                    method_version = 2L,
                    packages = vapply(c("edgeR", "limma", "msigdbr", "Seurat"),
                                      function(x) as.character(packageVersion(x)), character(1)))
  policy <- pdo_cache_policy()
  if (policy$force_rebuild && policy$replot_only) stop("Force-rebuild and replot-only cannot both be set.")
  reuse <- FALSE
  if (file.exists(cache_path) && !policy$force_rebuild) {
    result <- readRDS(cache_path)
    reuse <- identical(result$signature, signature)
  }
  if (!reuse && policy$replot_only) stop("Replot requires a current matching DE cache.")
  if (!reuse) {
    message("Loading raw counts for six independent baseline PDOs ...")
    if (is.null(pdos)) pdos <- readRDS(inputs[["pdo"]])
    required <- c("orig.ident", "Batch")
    if (!all(required %in% names(pdos@meta.data))) stop("Missing sample or Batch metadata.")
    meta <- pdos@meta.data %>% rownames_to_column("cell") %>%
      transmute(cell, sample = as.character(orig.ident), batch = as.character(Batch),
                sur = stringr::str_extract(sample, "^SUR[0-9]+")) %>%
      filter(sample != PDO_EXCLUDED_SAMPLE, !grepl("_Treated_PDO$", sample))
    clinical <- readxl::read_excel(inputs[["clinical"]], sheet = "Sheet1") %>%
      transmute(sur = paste0("SUR", as.integer(SUR)),
                auc_rel = suppressWarnings(as.numeric(`AUCrel (FLOT)`)))
    if (anyDuplicated(clinical$sur)) stop("Clinical SUR identifiers are duplicated.")
    baseline <- meta %>% distinct(sample, sur) %>% left_join(clinical, by = "sur")
    auc_order <- baseline %>% filter(is.finite(auc_rel)) %>% arrange(auc_rel)
    if (nrow(auc_order) != 12L || anyDuplicated(auc_order$sur)) stop("Expected 12 independent numeric-AUC baseline PDOs.")
    if (!setequal(head(auc_order$sur, 3L), sensitive_sur) ||
        !setequal(tail(auc_order$sur, 3L), resistant_sur)) stop("Fixed groups no longer match the three AUC extremes.")
    meta <- meta %>% filter(sur %in% c(sensitive_sur, resistant_sur)) %>%
      mutate(sensitivity_group = factor(ifelse(sur %in% sensitive_sur, "Sensitive", "Resistant"),
                                        levels = c("Sensitive", "Resistant")))
    if (anyNA(meta$batch) || any(!nzchar(meta$batch)) ||
        any(meta$batch == "Unknown")) stop("Selected cells lack known Batch metadata.")
    sample_meta <- meta %>% distinct(sample, sur, batch, sensitivity_group) %>%
      left_join(clinical, by = "sur") %>% arrange(sensitivity_group, auc_rel)
    if (nrow(sample_meta) != 6L || anyDuplicated(sample_meta$sur) ||
        anyDuplicated(sample_meta$sample)) stop("Require exactly six independent PDOs, with one Batch per PDO.")
    sample_meta <- sample_meta %>% left_join(meta %>% count(sample, name = "cell_n"), by = "sample")
    if (any(sample_meta$cell_n < min_cells)) stop("Overall PDO has too few cells.")
    counts <- Seurat::GetAssayData(pdos, assay = "RNA", layer = "counts")
    if (!all(meta$cell %in% colnames(counts))) stop("Raw RNA counts lack selected cells.")
    counts <- counts[, meta$cell, drop = FALSE]
    if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != round(counts@x))) {
      stop("RNA counts must be finite nonnegative raw integers.")
    }
    if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) stop("Gene identifiers must be unique.")
    overall_group <- factor(meta$sample, levels = sample_meta$sample)
    overall_design <- Matrix::sparse.model.matrix(~0 + overall_group)
    colnames(overall_design) <- levels(overall_group)
    overall_counts <- as.matrix(counts %*% overall_design)
    coverage <- sample_meta
    hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    gene_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
    gene_sets <- lapply(gene_sets, unique)
    if (length(gene_sets) != 50L) stop("Expected exactly 50 human Hallmark pathways.")
    models <- list()
    audit <- list()
    for (scope in "Overall") {
      message("Pseudobulk model: ", scope)
      sm <- sample_meta
      cm <- overall_counts[, sm$sample, drop = FALSE]
      group_sizes <- table(factor(sm$sensitivity_group, levels = c("Sensitive", "Resistant")))
      reason <- ""
      design <- NULL
      if (any(group_sizes < 2L)) {
        reason <- "Fewer than two independent PDOs per group"
      } else {
        sm$batch <- droplevels(factor(sm$batch))
        sm$sensitivity_group <- factor(sm$sensitivity_group, levels = c("Sensitive", "Resistant"))
        model_formula <- if (nlevels(sm$batch) > 1L) ~ batch + sensitivity_group else ~ sensitivity_group
        design <- model.matrix(model_formula, data = sm)
        rownames(design) <- sm$sample
        if (qr(design)$rank < ncol(design)) reason <- "Batch and response are confounded; model not estimable"
        else if (nrow(design) - ncol(design) < 2L) reason <- "Fewer than two residual degrees of freedom"
      }
      audit[[scope]] <- data.frame(scope = scope, sensitive_n = unname(group_sizes[1]),
        resistant_n = unname(group_sizes[2]), batch_n = n_distinct(sm$batch),
        design = if (is.null(design)) NA_character_ else paste(colnames(design), collapse = "; "),
        residual_df = if (is.null(design)) NA_integer_ else nrow(design) - qr(design)$rank,
        status = if (nzchar(reason)) "Not tested" else "Tested", reason = reason)
      if (nzchar(reason)) {
        if (scope == "Overall") stop("Primary model cannot be tested: ", reason)
        next
      }
      dge <- edgeR::DGEList(counts = cm)
      keep <- edgeR::filterByExpr(dge, design = design)
      if (sum(keep) < 100L) stop("Too few expressed genes in ", scope)
      filtering <- data.frame(scope = scope, gene = rownames(dge), retained = keep)
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
      fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
      contrast <- rep(0, ncol(design))
      names(contrast) <- colnames(design)
      contrast["sensitivity_groupResistant"] <- 1
      test <- edgeR::glmQLFTest(fit, contrast = contrast)
      de <- edgeR::topTags(test, n = Inf, sort.by = "PValue")$table %>%
        rownames_to_column("gene") %>% mutate(scope = scope,
          direction = ifelse(logFC > 0, "Higher in resistant", "Higher in sensitive"))
      indices <- lapply(gene_sets, function(g) which(rownames(dge) %in% g))
      indices <- indices[lengths(indices) >= 10L]
      if (!length(indices)) stop("No Hallmark pathway has at least 10 tested genes.")
      camera <- limma::camera(dge, index = indices, design = design, contrast = contrast,
                             inter.gene.cor = 0.01, sort = FALSE) %>%
        rownames_to_column("pathway") %>% mutate(scope = scope,
          interpretation = ifelse(Direction == "Up", "Higher in resistant", "Higher in sensitive"))
      models[[scope]] <- list(metadata = sm, counts = cm, dge = dge, design = design,
                             contrast = contrast, de = de, camera = camera, filtering = filtering,
                             logcpm = edgeR::cpm(dge, log = TRUE, prior.count = 2),
                             mds = limma::plotMDS(dge, plot = FALSE),
                             fit = fit, test = test)
    }
    result <- list(signature = signature, sample_meta = sample_meta, cell_counts = coverage,
                   raw_overall_counts = overall_counts,
                   gene_sets = gene_sets, models = models, audit = bind_rows(audit))
    saveRDS(result, cache_path)
  } else message("Reusing matching persistent DE cache.")
  de_all <- bind_rows(lapply(result$models, `[[`, "de"))
  camera_all <- bind_rows(lapply(result$models, `[[`, "camera"))
  de_all$FDR_across_tested_scopes <- p.adjust(de_all$PValue, "BH")
  camera_all$FDR_across_tested_scopes <- p.adjust(camera_all$PValue, "BH")
  write_csv <- function(x, name) data.table::fwrite(as.data.frame(x), file.path(paths[["tables"]], name))
  write_csv(result$sample_meta, "Auto_sample_metadata.csv")
  write_csv(result$cell_counts, "Auto_cell_counts.csv")
  write_csv(result$audit, "Auto_model_audit.csv")
  write_csv(de_all, "Auto_differential_expression.csv")
  write_csv(camera_all, "Auto_hallmark_camera.csv")
  write_csv(bind_rows(lapply(result$models, `[[`, "filtering")), "Auto_gene_filtering.csv")
  designs <- bind_rows(lapply(names(result$models), function(scope) {
    as.data.frame(result$models[[scope]]$design, check.names = FALSE) %>%
      rownames_to_column("sample") %>% mutate(scope = scope)
  }))
  write_csv(designs, "Auto_design_matrices.csv")
  summary_lines <- c(
    "FLOT response extremes: batch-adjusted pseudobulk differential expression",
    "Sensitive: SUR680, SUR629, SUR1090; resistant: SUR727, SUR1181, SUR1363.",
    "Contrast: resistant minus sensitive. Raw RNA counts summed per independent PDO.",
    "Overall model includes every selected baseline cell; no state-stratified analyses.",
    "Design: batch + response group, except a single batch requires response group only.",
    "edgeR: filterByExpr, TMM, robust dispersion and quasi-likelihood F-test; BH FDR < 0.05.",
    "Hallmark: CAMERA competitive directional test; inter-gene correlation 0.01; >=10 tested genes/set.",
    "Up means higher in resistant; Down means higher in sensitive. Enrichment is not proof of activation.",
    paste0("Eligibility: >=", min_cells, " cells/PDO; six independent PDOs; full-rank batch-adjusted design."),
    "Primary overall results may reflect different proportions of states as well as within-state expression.",
    "One batch contains only SUR1363; that sample does not supply a within-batch response comparison.",
    "Only six PDOs and two primary residual df: exploratory results, sensitive to individual PDOs.",
    "Batch denotes recorded acquisition groups; unrecorded technical effects cannot be adjusted.",
    "BH FDR is calculated across all tested genes and, separately, all 50 Hallmark pathways.",
    paste0("Cache reused: ", reuse), paste0("msigdbr version: ", packageVersion("msigdbr")),
    "", capture.output(print(result$audit, row.names = FALSE)), "",
    unlist(lapply(names(result$models), function(scope) {
      d <- result$models[[scope]]$de
      g <- result$models[[scope]]$camera
      paste0(scope, ": ", nrow(d), " genes tested; ", sum(d$FDR < 0.05),
             " genes and ", sum(g$FDR < 0.05), " Hallmarks at within-scope FDR < 0.05.")
    }))
  )
  writeLines(summary_lines, file.path(paths[["reports"]], "Auto_summary.txt"))
  ####################
  # Two figure-only pages: paired QC panels, then paired DE/pathway panels.
  # Full methods and non-significant statistics remain in reports/tables.
  model <- result$models[["Overall"]]
  md <- model$metadata
  colours <- c(Sensitive = "#0072B2", Resistant = "#D55E00")
  plot_theme <- theme_classic(base_size = 15) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 11), legend.position = "bottom",
          legend.box = "vertical", legend.text = element_text(size = 11),
          legend.title = element_text(size = 11), plot.margin = margin(14, 16, 14, 16))
  mds <- data.frame(MDS1 = model$mds$x, MDS2 = model$mds$y, md)
  p_mds <- ggplot(mds, aes(MDS1, MDS2, colour = sensitivity_group, shape = batch)) +
    geom_point(size = 4) +
    geom_text(aes(label = sur), vjust = -0.9, size = 4, show.legend = FALSE) +
    scale_colour_manual(values = colours) +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    scale_y_continuous(expand = expansion(mult = 0.18)) +
    labs(title = "Sample expression structure", subtitle = "TMM-normalized pseudobulks; batch remains visible",
         x = "Leading logFC dimension 1", y = "Leading logFC dimension 2",
         colour = "Response", shape = "Batch") + plot_theme +
    guides(colour = guide_legend(order = 1, nrow = 1),
           shape = guide_legend(order = 2, nrow = 1))
  md$sur <- factor(md$sur, levels = md$sur)
  p_cells <- ggplot(md, aes(sur, cell_n, fill = sensitivity_group)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = scales::comma(cell_n)), vjust = -0.5, size = 4, show.legend = FALSE) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Cells contributing per PDO", subtitle = "Each PDO is one biological replicate",
         x = NULL, y = "Cells", fill = "Response") + plot_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  d <- model$de
  d$significance <- factor(ifelse(d$FDR < 0.05, d$direction, "FDR >= 0.05"),
                          levels = c("Higher in sensitive", "FDR >= 0.05", "Higher in resistant"))
  d$minus_log10_fdr <- -log10(pmax(d$FDR, .Machine$double.xmin))
  p_volcano <- ggplot(d, aes(logFC, minus_log10_fdr, colour = significance)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey35") +
    geom_vline(xintercept = 0, colour = "grey65") +
    scale_colour_manual(values = c("Higher in sensitive" = colours[["Sensitive"]],
                                   "FDR >= 0.05" = "grey60",
                                   "Higher in resistant" = colours[["Resistant"]]), drop = FALSE) +
    labs(title = "Overall differential gene expression",
         subtitle = paste0("Batch-adjusted edgeR QL; ", sum(d$FDR < 0.05), "/", nrow(d),
                           " genes at BH FDR < 0.05; dashed line marks FDR 0.05"),
         x = "log2 fold change: resistant minus sensitive", y = "-log10 BH FDR", colour = NULL) + plot_theme
  # Display only significant pathways; the full 50-pathway family determines FDR.
  cam <- model$camera %>% filter(FDR < 0.05) %>% arrange(FDR, PValue)
  cam$label <- stringr::str_to_sentence(gsub("_", " ", sub("^HALLMARK_", "", cam$pathway)))
  cam$label <- factor(cam$label, levels = rev(cam$label))
  cam$signed_evidence <- ifelse(cam$Direction == "Up", 1, -1) * -log10(pmax(cam$FDR, .Machine$double.xmin))
  p_pathways <- ggplot(cam, aes(signed_evidence, label, colour = interpretation)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_vline(xintercept = c(-1, 1) * -log10(0.05), linetype = 2, colour = "grey65") +
    geom_segment(aes(x = 0, xend = signed_evidence, yend = label), linewidth = 1, show.legend = FALSE) +
    geom_point(aes(size = NGenes)) +
    scale_colour_manual(values = c("Higher in sensitive" = colours[["Sensitive"]],
                                   "Higher in resistant" = colours[["Resistant"]])) +
    scale_size_continuous(range = c(4, 8), name = "Tested genes in set") +
    scale_x_continuous(limits = c(-1, 1) * max(2, abs(cam$signed_evidence)) * 1.15) +
    labs(title = "Hallmark pathways associated with FLOT sensitivity",
         subtitle = paste0(nrow(cam), " of 50 pathways at BH FDR < 0.05; batch-adjusted CAMERA"),
         x = "Signed -log10 BH FDR: higher in sensitive (left), higher in resistant (right)",
         y = NULL, colour = NULL) + plot_theme
  if (!nrow(cam)) {
    p_pathways <- ggplot() + annotate("text", x = 0, y = 0,
                                    label = "No Hallmark pathways at BH FDR < 0.05", size = 7) +
      theme_void() + labs(title = "Overall Hallmark enrichment")
  }
  grDevices::pdf(pdf_path, width = 13.33, height = 7.5, useDingbats = FALSE)
  pdf_device <- dev.cur()
  on.exit(if (pdf_device %in% dev.list()) dev.off(pdf_device), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_mds, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_cells, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
  # Fit titles, axes and separate legend rows to half-slide panel widths.
  p_volcano <- p_volcano +
    labs(title = "Overall differential expression",
         subtitle = paste0("Batch-adjusted edgeR QL; ", sum(d$FDR < 0.05), "/", nrow(d),
                           " genes at FDR < 0.05\nDashed line: BH FDR 0.05"),
         x = "log2 fold change\nResistant minus sensitive") +
    guides(colour = guide_legend(ncol = 1))
  p_pathways <- p_pathways +
    labs(title = "Significant Hallmark pathways",
         subtitle = paste0(nrow(cam), "/50 pathways at BH FDR < 0.05\nBatch-adjusted CAMERA"),
         x = "Signed -log10 BH FDR\nSensitive (left); resistant (right)") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 22)) +
    guides(colour = guide_legend(order = 1, ncol = 1),
           size = guide_legend(order = 2, nrow = 1))
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_volcano, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_pathways, newpage = FALSE,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
  dev.off(pdf_device)
  write_csv(cam, "Auto_report_hallmark_plot.csv")
  ####################
  outputs <- c(cache_path, pdf_path, list.files(paths[["tables"]], full.names = TRUE),
               file.path(paths[["reports"]], "Auto_summary.txt"))
  pdo_require_files(outputs)
  if (file.info(pdf_path)$size < 1000L) stop("Report PDF is unexpectedly small.")
  pdo_write_run_summary("analysis/clinical/Auto_centred_clinical_pseudobulk_de.R", out_dir,
    inputs = inputs, outputs = outputs,
    parameters = list(contrast = "resistant minus sensitive", design = "batch + sensitivity_group",
      min_cells = min_cells, min_group_n = 2L, min_residual_df = 2L, gene_fdr = 0.05,
      camera_inter_gene_cor = 0.01, cache_reused = reuse,
      elapsed_minutes = round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2)), cache = policy)
  message("Completed pseudobulk DE and pathway PDF: ", pdf_path)
  invisible(result)
}
